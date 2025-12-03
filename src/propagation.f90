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
    subroutine HamScale()
        implicit none
        real(f8) :: TZmax, TZmin, Trmax, Trmin 
        real(f8) :: Vmax, Vmin, Umax, Umin, Rotmax, Rotmin 
        real(f8) :: Hmax, Hmin

        call getZKinMat()
        call getCPMat()
        call getRotMat()
        call getVintMat()

        TZmin = 0.0_f8
        Trmin = 0.0_f8
        Rotmin = 0.0_f8
        Umin = minval(CPMat)
        Vmin = minval(VintMat)
        Hmin = TZmin+Trmin+Rotmin+Umin+Vmin

        TZmax = (pi/(IALR%Z_range(2)-IALR%Z_range(1)))**2/(2.0_f8*massTot)
        Trmax = (pi/(IALR%r_range(2)-IALR%r_range(1)))**2/(2.0_f8*massBC)
        Rotmax = maxval(rotMat)
        Umax = maxval(CPMat)
        Vmax = maxval(VintMat)
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


    end subroutine HamScale
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
        integer :: v, vp, j, jp, K, Kp, iPES, jPES, ichnl, jchnl
        integer :: iZ, irD, irP, ith, iZloc

        allocate(VintMat(IALR%nZ_IALR,nChannels,nChannels))

!> Vint(R) = <v'j'K'i'|V_i'i(R,r,theta)|vjKi>, i for PES
!> See J. Chem. Phys. 162, 074301 (2025)
        VintMat = 0.0_f8
        do ichnl = 1, nChannels
            do jchnl = 1, nChannels
                v = qn_channel(ichnl,1)
                j = qn_channel(ichnl,2)
                K = qn_channel(ichnl,3)
                iPES = qn_channel(ichnl,4)

                vp = qn_channel(jchnl,1)
                jp = qn_channel(jchnl,2)
                Kp = qn_channel(jchnl,3)
                jPES = qn_channel(jchnl,4)

                if (K /= Kp) cycle 
                if (iPES == jPES) then 
                    !> Use PODVR WFvj 
                    do iZ = IALR%nZ_I+1, IALR%nZ_IALR
                        do irP = 1, IALR%nPODVR
                            do ith = 1, IALR%jasy 
                                iZloc = iZ - IALR%nZ_I
                                cth = asyANode(ith)
                                w = dsqrt(asyAWeight(ith))
                                YjK = spgndr(j,K,cth)
                                YjKp = spgndr(jp,K,cth)
                                WFvjK = w*YjK*asyBC_POWF(iPES,irP,v,j)
                                WFvjKp = w*YjKp*asyBC_POWF(jPES,irP,vp,jp)
                                VintMat(iZ,ichnl,jchnl) = VintMat(iZ,ichnl,jchnl) &
                                                          + WFvjK*ALR_Vdiag(iPES,iZloc,irP,ith)*WFvjKp
                                VintMat(iZ,ichnl,jchnl) = min(VintMat(iZ,ichnl,jchnl), VMaxCut)
                            end do 
                        end do 
                    end do 
                else
                    !> Use DVR WFvj
                    do iZ = IALR%nZ_I+1, IALR%nZ_IALR
                        do irD = 1, IALR%vint 
                            do ith = 1, IALR%jasy 
                                iZloc = iZ - IALR%nZ_I
                                cth = asyANode(ith)
                                w = dsqrt(asyAWeight(ith))
                                YjK = spgndr(j,K,cth)
                                YjKp = spgndr(jp,K,cth)
                                WFvjK = w*YjK*asyBC_DVRWF(iPES,irD,v,j)
                                WFvjKp = w*YjKp*asyBC_DVRWF(jPES,irD,vp,jp)
                                VintMat(iZ,ichnl,jchnl) = VintMat(iZ,ichnl,jchnl) &
                                                          + WFvjK*ALR_Voff(iPES,jPES,iZloc,irD,ith)*WFvjKp
                                VintMat(iZ,ichnl,jchnl) = min(VintMat(iZ,ichnl,jchnl), VMaxCut)
                            end do 
                        end do 
                    end do 
                end if 
            end do 
        end do 

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
    