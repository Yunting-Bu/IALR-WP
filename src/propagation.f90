module propMod
    use machina_basic, only : f8, c8
    use gPara
!    use potentMod
!    use basisMod
!    use initWPMod
    implicit none

    public
    private :: getCPMat, getRotMat, getZKinMat, getVintMat, lambdaMinus, lambdaPlus

contains 

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine preWF0AndWF1()
        implicit none
        integer :: iEtot
       
        allocate(ALR_TDWF(IALR%nZ_IALR,nChannels))
        allocate(ALR_auxWFm(IALR%nZ_IALR,nChannels))
        allocate(ALR_auxWFp(IALR%nZ_IALR,nChannels))
        allocate(ALR_TIDWF(IALR%nZ_IALR,nChannels,nEtot))

        !> |WFp> = D(2*\hat{H_s} |TDWF> - D* |WFm>)
        !> For k = 1, |WFm> = |phi0>, |TDWF> = |phi1>
        !> |phi0> = |initWP>, |phi1> = \hat{H_s}|phi0>

        ALR_auxWFm = initTotWP
        ALR_TDWF = initTotWP
        
        call lrHamAction(ALR_TDWF)
        call asyHamAction(ALR_TDWF)
        !> Add interaction range later

        do iEtot = 1, nEtot
            ALR_TIDWF(:,:,iEtot) = cmplx(initTotWP(:,:),0.0_f8, kind=c8)
        end do 

    end subroutine preWF0AndWF1
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine ChebyshevRecursion()
        implicit none
        real(f8) :: auxWF(IALR%nZ_IALR,nChannels)
        real(f8) :: auxWF1(IALR%nZ_IALR,nChannels)

        !> |WFp> = D(2*\hat{H_s} |TDWF> - D* |WFm>)
        !> Long-range 
        call lrHamAction(ALR_TDWF)
        !> Asymptotic range
        call asyHamAction(ALR_TDWF)
        !> interaction range
        !> ...

        auxWF = ALR_auxWFm
        call dampAction(auxWF)
        auxWF1 = ALR_TDWF
        auxWF = 2.0_f8*ALR_TDWF - auxWF
        call dampAction(auxWF)

    end subroutine ChebyshevRecursion
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine dampAction(WF)
        implicit none
        real(f8), intent(inout) :: WF(:,:)
        integer :: ichnl, v, j, iPES 
                
        !> Damp in lr
        do ichnl = 1, nChannels
            WF(:,ichnl) = Flr(:)*WF(:,ichnl)
        end do 

        !> Damp in asy
        !> only action when v /= v0, j /= j0, iPES /= initPES
        do ichnl = 1, nChannels
            v = qn_channel(ichnl,1)
            j = qn_channel(ichnl,2)
            iPES = qn_channel(ichnl,4)
            if (v == initWP%v0 .or. j == initWP%j0 .or. iPES == initWP%initPES) cycle 
            WF(:,ichnl) = Fasy(:)*WF(:,ichnl)
        end do 

        !> Damp in int
        !> ...

    end subroutine dampAction
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine HamScale()
        implicit none
        real(f8) :: TZmax, TZmin, Trmax, Trmin 
        real(f8) :: Vmax, Vmin, Umax, Umin, Rotmax, Rotmin 
        real(f8) :: Hmax, Hmin
        real(f8) :: dtMin, dtMax
        integer :: iEtot

        call getZKinMat()
        call getCPMat()
        call getRotMat()
        call getVintMat()

        TZmin = 0.0_f8
        Trmin = 0.0_f8
        Rotmin = 0.0_f8
        Umin = minval(CPMat)
        Vmin = minval(INT_Vadia)
        Hmin = TZmin+Trmin+Rotmin+Umin+Vmin

        TZmax = (pi/(IALR%Z_range(2)-IALR%Z_range(1)))**2/(2.0_f8*massTot)
        Trmax = (pi/(IALR%r_range(2)-IALR%r_range(1)))**2/(2.0_f8*massBC)
        Rotmax = maxval(rotMat)
        Umax = maxval(CPMat)
        Vmax = maxval(INT_Vadia)
        Hmax = TZmax+Trmax+Rotmax+Umax+Vmax

        Hplus = (Hmax+Hmin)/2.0_f8
        Hminus = (Hmax-Hmin)/2.0_f8

        Z_KinMat(:,:) = (Z_KinMat(:,:)-Hplus)/Hminus
        CPMat(:,:,:) = (CPMat(:,:,:)-Hplus)/Hminus
        rotMat(:,:) = (rotMat(:,:)-Hplus)/Hminus
        asyBC_Evj(:,:,:) = (asyBC_Evj(:,:,:)-Hplus)/Hminus
        adiaVBC(:,:) = (adiaVBC(:,:)-Hplus)/Hminus
        INT_Vadia(:,:,:,:) = (INT_Vadia(:,:,:,:)-Hplus)/Hminus
        VintMat(:,:,:) = (VintMat(:,:,:)-Hplus)/Hminus

        allocate(VeffMat(IALR%nZ_IALR,nChannels,nChannels))
        VeffMat(:,:,:) = CPMat(:,:,:)+VintMat(:,:,:)

        do iEtot = 1, nEtot
            ChebyAngle(iEtot) = acos((Etot(iEtot)-Hplus)/Hminus)
        end do
        
        dtMin = 1.0_f8/(Hminus*sin(maxval(ChebyAngle)))
        dtMax = 1.0_f8/(Hminus*sin(minval(ChebyAngle)))
        timeStep = (dtMax+dtMin)/2.0_f8

        write(outFileUnit,'(1x,a)') '=====> Chebyshev scale information <====='
        write(outFileUnit,'(1x,a,f12.5,2x,a,f12.5)') 'Hmax = ', Hmax, 'Hmin = ', Hmin 
        write(outFileUnit,'(1x,a,f12.5)') 'TimeStep Max = ', dtMax
        write(outFileUnit,'(1x,a,f12.5)') 'TimeStep Min = ', dtMin 
        write(outFileUnit,'(1x,a,f12.5)') 'TimeStep Avg = ', timeStep
        write(outFileUnit,*) ''


    end subroutine HamScale
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine asyHamAction(WF)
        implicit none
        real(f8), intent(inout) :: WF(:,:)
        real(f8), allocatable :: Wtemp(:,:)
        real(f8) :: alpha, beta 
        integer :: iZ, ichnl, v, j, iPES

        allocate(Wtemp(IALR%nZ_IALR,nChannels))

        alpha = 1.0_f8
        beta = 0.0_f8
        !> Action of \hat{T_R}
        call dgemm('N', 'N', IALR%nZ_IALR, nChannels, IALR%nZ_IALR, alpha, & 
                    Z_KinMat, IALR%nZ_IALR, WF, IALR%nZ_IALR, beta, Wtemp, IALR%nZ_IALR)
        WF = Wtemp

        !> Action of \hat{h(r)}
        do ichnl = 1, nChannels
            v = qn_channel(ichnl,1)
            j = qn_channel(ichnl,2)
            iPES = qn_channel(ichnl,4)
            WF(:,ichnl) = WF(:,ichnl) * asyBC_Evj(iPES,v,j)
        end do 

        !> Action of \hat{U}+\hat{Vint}
        do iZ = 1, IALR%nZ_IALR 
            call dgemv('N', nChannels, nChannels, alpha, VeffMat(iZ,:,:), &
                        nChannels, WF(iZ,:), 1, beta, Wtemp(iZ,:), 1)
        end do 
        WF = Wtemp
        deallocate(Wtemp)

    end subroutine asyHamAction
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine lrHamAction(WF)
        implicit none
        real(f8), intent(inout) :: WF(:,:)
        integer :: iZ, i, j
        real(f8), allocatable :: Wtemp(:,:), auxWF(:,:), Veff(:,:,:)
        real(f8) :: alpha, beta

        allocate(Wtemp(IALR%nZ_IALR,nLrChnl))
        allocate(auxWF(IALR%nZ_IALR,nLrChnl))
        allocate(Veff(nZ,nLrChnl,nLrChnl))

        !> Construction of Veff and auxWF
        do iZ = 1, IALR%nZ_IALR
            do i = 1, nLrChnl
                auxWF(iZ,i) = WF(iZ, lrChnlNo(i))
                do j = 1, nLrChnl
                    Veff(iZ,i,j) = VeffMat(iZ , lrChnlNo(i) , lrChnlNo(j))
                end do
            end do
        end do

        alpha = 1.0_f8
        beta  = 0.0_f8

        !> \hat{T_R}
        call dgemm('N','N', IALR%nZ_IALR, nLrChnl, IALR%nZ_IALR, alpha, &
                    Z_KinMat, nZ_IALR, auxWF, nZ_IALR, beta, Wtemp, nZ_IALR)
        auxWF = Wtemp

        !> \hat{h(r)}
        auxWF(:,:) = auxWF(:,:) * asyBC_Evj(initWP%initPES,initWP%v0,initWP%j0)

        !> \hat{U}+\hat{Vint}
        do iZ = 1, IALR%nZ_IALR
            call dgemv('N', nLrChnl, nLrChnl, alpha, Veff(iZ,:,:), &
                    nLrChnl, auxWF(iZ,:), 1, beta, Wtemp(iZ,:), 1)
        end do
        auxWF = Wtemp

        do i = nLrChnl
            WF(:,lrChnlNo(i)) = auxWF(:,i)
        end do

        deallocate(Wtemp, Veff)

    end subroutine lrHamAction
!> ------------------------------------------------------------------------------------------------------------------ <!
        
!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine getZKinMat()
        implicit none
        real(f8) :: range_Z, fact
        integer :: i, j 

        range_Z = IALR%Z_range(2)-IALR%Z_range(1)

        allocate(Z_KinMat(IALR%nZ_IALR,IALR%nZ_IALR))

!> FFT could be better
        fact = pi*pi/(4.0_f8*massTot*range_Z*range_Z)
        do i = 1, IALR%nZ_IALR 
            do j = 1, i-1
                Z_KinMat(i,j) = fact  * &
                                ((dsin(pi * (i-j) / (2.0_f8 * (IALR%nZ_IALR+1))))**(-2) - &
                                (dsin(pi * (i+j) / (2.0_f8 * (IALR%nZ_IALR+1))))**(-2)) * (-1.0_f8) **(i-j)
            end do
            Z_KinMat(i,i) = fact * &
                            ((2.0_f8 * (IALR%nZ_IALR+1)**2 +1.0_f8) / 3.0_f8 - &
                            (dsin(pi * i / (IALR%nZ_IALR+1)))**(-2))
        end do 
        do i = 1, IALR%nZ_IALR
            do j = i+1, IALR%nZ_IALR
                Z_KinMat(i,j) = Z_KinMat(j,i)
            end do
        end do
    
    end subroutine getZKinMat
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine getRotMat()
        implicit none
        real(f8) :: fact, jeigen
        integer :: ir, j
        
        allocate(rotMat(IALR%vint,IALR%jint))
        do ir = 1, IALR%vint 
            fact = 1.0_f8/(2.0_f8*massTot*r_Int(ir)**2)
            do j = 1, IALR%jint
                jeigen = real(j*(j+1),f8)
                rotMat(ir,j) = min(fact*jeigen, TMaxCut)
            end do 
        end do 

    end subroutine getRotMat
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine getCPMat()
        implicit none
        real(f8) ::  delta, term1, term2, term3
        real(f8), allocatable :: fact(:)
        real(f8), allocatable :: CPE(:), CPM(:,:)
        real(f8), allocatable :: work(:)
        integer :: v, qn_j, vp, qn_jp, K, Kp
        integer :: iZ, ichnl, jchnl, lwork, info

        allocate(CPMat(IALR%nZ_IALR, nChannels, nChannels))
        allocate(fact(IALR%nZ_IALR))
        allocate(CPE(nChannels), CPM(nChannels,nChannels))

        do iZ = 1, IALR%nZ_IALR
            fact(iZ) = 1.0_f8/(2.0_f8*massTot*Z_IALR(iZ)**2)
        end do

        do ichnl = 1, nChannels
            v = qn_channel(ichnl,1)
            qn_j = qn_channel(ichnl,2)
            K = qn_channel(ichnl,3)
            do jchnl = 1, ichnl
                vp = qn_channel(jchnl,1)
                qn_jp = qn_channel(jchnl,2)
                Kp = qn_channel(jchnl,3)
                if (v /= vp .or. qn_j /= qn_jp) cycle

                term1 = merge((initWP%Jtot*(initWP%Jtot+1.0_f8)+qn_j*(qn_j+1.0_f8)-2.0_f8*K*K), &
                               0.0_f8, K == Kp)
            
                delta = merge(dsqrt(2.0_f8), 1.0_f8, K == 0)
                term2 = merge(delta*lambdaPlus(initWP%Jtot,K)*lambdaPlus(qn_j,K), &
                              0.0_f8, Kp == K+1)
                
                delta = merge(dsqrt(2.0_f8), 1.0_f8, K == 1)
                term3 = merge(delta*lambdaMinus(initWP%Jtot,K)*lambdaMinus(qn_j,K), &
                              0.0_f8, Kp == K-1)

                CPM(ichnl, jchnl) = term1 - term2 - term3
                CPM(jchnl, ichnl) = CPM(ichnl, jchnl)
            end do 
        end do
                

        allocate(work(1))
        call dsyev('V', 'U', nChannels, CPM, nChannels, CPE, work, -1, info)
        lwork = int(work(1))
        deallocate(work)

        allocate(work(lwork))
        call dsyev('V', 'U', nChannels, CPM, nChannels, CPE, work, lwork, info)
        if (info /= 0) then
            write(outFileUnit,*) "Error in DSYEV of CPE, info =", info
            write(outFileUnit,*) "POSITION: propagation.f90, subroutine getCPMat()"
            stop
        end if
        deallocate(work)

        do iZ = 1, IALR%nZ_IALR
            do ichnl = 1, nChannels
                do jchnl = 1, nChannels
                    delta = 0.0_f8
                    do K = 1, nChannels
                        CPE(K) = min(CPE(K)*fact(iZ), TMaxCut)
                        delta = delta + CPM(jchnl,K)*CPE(K)*CPM(ichnl,K)
                    end do 
                    CPMat(iZ,jchnl,ichnl) = delta
                end do 
            end do 
        end do 

        deallocate(fact, CPM, CPE)
    end subroutine getCPMat
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine getVintMat()
        implicit none
        real(f8) :: YjK, WFvjK, YjKp, WFvjKp
        real(f8) :: w, cth
        real(f8), external :: spgndr
        real(f8), allocatable :: YjK_table(:,:,:)
        real(f8), allocatable :: WFvjK_PO(:,:)
        real(f8), allocatable :: WFvjK_DVR(:,:)
        integer :: v, vp, j, jp, K, Kp, iPES, jPES, ichnl, jchnl
        integer :: iZ, irD, irP, ith
        integer :: nPtot, nDtot 

        allocate(VintMat(IALR%nZ_IALR,nChannels,nChannels))
        allocate(YjK_table(IALR%jasy,0:IALR%jasy,0:max(IALR%jasy,initWP%Jtot)))
        allocate(WFvjK_PO(IALR%nPODVR,nChannels))
        allocate(WFvjK_DVR(IALR%vint,nChannels))

!> Vint(R) = <v'j'K'i'|V_i'i(R,r,theta)|vjKi>, i for PES
!> See J. Chem. Phys. 162, 074301 (2025)

        do ith = 1, IALR%jasy 
            do j = 0, IALR%jasy 
                do K = 0, max(IALR%jasy,initWP%Jtot)
                    if (K > j) cycle
                    w = dsqrt(asyAWeight(ith))
                    YjK_table(ith,j,K) = w*spgndr(j,K,cth)
                end do 
            end do 
        end do
        
        do ichnl = 1, nChannels
            v = qn_channel(ichnl,1)
            j = qn_channel(ichnl,2)
            K = qn_channel(ichnl,3)
            iPES = qn_channel(ichnl,4)
            do ith = 1, IALR%jasy 
                WFvjK_PO(:,ichnl) = YjK_table(ith,j,K)*asyBC_POWF(iPES,:,v,j)
                if (nPES == 1) cycle 
                WFvjK_DVR(:,ichnl) = YjK_table(ith,j,K)*asyBC_DVRWF(iPES,:,v,j)
            end do 
        end do 

        VintMat = 0.0_f8
        do ichnl = 1, nChannels
            do jchnl = 1, nChannels
                K = qn_channel(ichnl,3)
                iPES = qn_channel(ichnl,4)

                Kp = qn_channel(jchnl,3)
                jPES = qn_channel(jchnl,4)

                if (K /= Kp) cycle 
                if (iPES == jPES) then 
                    !> Use PODVR WFvj 
                    do iZ = 1, IALR%nZ_IALR
                        do irP = 1, IALR%nPODVR
                            do ith = 1, IALR%jasy 
                                VintMat(iZ,ichnl,jchnl) = VintMat(iZ,ichnl,jchnl) &
                                                          + WFvjK_PO(irP,ichnl)*ALR_Vdiag(iPES,iZloc,irP,ith) &
                                                          * WFvjK_PO(irP,jchnl)
                            end do 
                        end do 
                        VintMat(iZ,ichnl,jchnl) = min(VintMat(iZ,ichnl,jchnl), VMaxCut)
                    end do 
                else
                    !> Use DVR WFvj
                    do iZ = 1, IALR%nZ_IALR
                        do irD = 1, IALR%vint 
                            do ith = 1, IALR%jasy 
                                VintMat(iZ,ichnl,jchnl) = VintMat(iZ,ichnl,jchnl) &
                                                          + WFvjK_DVR(irD,ichnl)*ALR_Voff(iPES,jPES,iZloc,irD,ith) &
                                                          * WFvjK_DVR(irD,jchnl)
                            end do 
                        end do 
                        VintMat(iZ,ichnl,jchnl) = min(VintMat(iZ,ichnl,jchnl), VMaxCut)
                    end do 
                end if 
            end do 
        end do 
        deallocate(YjK_table,WFvjK_DVR,WFvjK_PO)

    end subroutine getVintMat
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    real(f8) function lambdaPlus(J,K) result(res)
        implicit none
        integer, intent(in) :: J, K

        res = dsqrt(J*(J+1.0_f8) - K*(K+1.0_f8))
        return

    end function lambdaPlus
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    real(f8) function lambdaMinus(J,K) result(res)
        implicit none
        integer, intent(in) :: J, K

        res = dsqrt(J*(J+1.0_f8) - K*(K-1.0_f8))
        return

    end function lambdaMinus
!> ------------------------------------------------------------------------------------------------------------------ <!

end module propMod
    