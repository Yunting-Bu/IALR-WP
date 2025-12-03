module potentMod
    use machina_basic, only : f8, BinReadWrite 
    use gPara
    implicit none

    public
    private :: getVabs
    
contains

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine Jacobi2Bond(Z, r, theta, massB, massC, bond)
        implicit none
        real(f8), intent(in) :: Z, r, theta 
        real(f8), intent(in) :: massB, massC 
        real(f8), intent(inout) :: bond(3)
        real(f8) :: rOB, rOC

!> bond(1) = rAB, bond(2) = rBC, bond(3) = rAC
!> theta in radian
        bond(2) = r
        rOB = massC / (massB+massC) * r 
        rOC = massB / (massB+massC) * r 
        bond(1) = dsqrt(Z*Z+rOB*rOB-2.0_f8*Z*rOB*dcos(theta))
        bond(3) = dsqrt(Z*Z+rOC*rOC+2.0_f8*Z*rOC*dcos(theta))

    end subroutine Jacobi2Bond
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine bond2Jacobi(bond, Z, r, theta, massB, massC)
        implicit none
        real(f8), intent(in) :: bond(3)
        real(f8), intent(in) :: massB, massC 
        real(f8), intent(out) :: Z, r, theta 
        real(f8) :: r1, r2, a1, a2

!> theta in radian
        r = bond(2)
        r1 = massB * r / (massB+massC)
        r2 = massC * r / (massB+massC)

        a1 = (bond(1)**2+r**2-bond(3)**2) / (2.0_f8*r*bond(1))
        Z = dsqrt(bond(1)**2+r1**2-2.0_f8*r1*bond(1)*a1)
        a2 = (Z**2+r2**2-bond(3)**2) / (2.0_f8*r2*Z)

        if (a2 > 1.0_f8) a2 = 1.0_f8
        if (a2 < -1.0_f8) a2 = -1.0_f8
        theta = acos(a2)
        if (Z == 0.0_f8) theta = pi / 2.0_f8
    
    end subroutine bond2Jacobi
!> ------------------------------------------------------------------------------------------------------------------ <!


!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine diagDiaVmat(bond, AtDMat, adiaV)
        implicit none
        real(f8), intent(in) :: bond(3)
        real(f8), intent(inout) :: AtDMat(nPES,nPES), adiaV(nPES)
        real(f8) :: diaV(nPES,nPES)
        real(f8), allocatable :: work(:)
        integer :: lwork, info

!> bond(3) is the three bonds coordinate
!> POT0 is the common PES interface
        call POT0(nPES,diaV,bond)
        if (nPES > 1) then 
            allocate(work(1))
            call dsyev('V', 'U', nPES, diaV, nPES, adiaV, work, -1, info)
            lwork = int(work(1))
            deallocate(work)

            allocate(work(lwork))
            call dsyev('V', 'U', nPES, diaV, nPES, adiaV, work, lwork, info)
            if (info /= 0) then
                write(outFileUnit,*) "Error in DSYEV of diaV, info =", info
                write(outFileUnit,*) "POSITION: potent.f90, subroutine diagDiaVmat()"
                stop
            end if
            deallocate(work)
        else 
            AtDMat = 1.0_f8
            adiaV(1) = diaV(1,1)
        end if 
    
    end subroutine diagDiaVmat
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine getVabs(range, Cabs, nGird, grid, dt, nabs, VabsMat)
        implicit none
        real(f8), intent(in) :: range, Cabs 
        integer, intent(in) :: nGird 
        real(f8), intent(in) :: grid(nGird)
        integer, intent(in) :: dt
        integer, intent(out) :: nabs
        real(f8), allocatable, intent(inout) :: VabsMat(:)
        real(f8) :: rangeAll, rangeNoAbs
        integer :: i

!> Since DVR grid don't have the boundary point
        rangeAll = grid(nGird) + grid(2) - grid(1)
        rangeNoAbs = rangeAll - range
        do i = 1, nGird
            if (grid(i) >= rangeNoAbs) then
                VabsMat(i) = dexp(-Cabs * ((grid(i)-rangeNoAbs)/(rangeAll - rangeNoAbs))**2 * dt)
            else
                VabsMat(i) = 1.0_f8
                nabs = i
            end if
        end do

        write(outFileUnit,'(1x,a,f15.9,a,f15.9,a)') 'Range of Vabs: [', rangeNoAbs, ', ', rangeAll, ' ]' 
        write(outFileUnit,'(1x,a,i5,a)') 'When nDVR >=', nabs, ', the absorbing potential starts to work.'
        write(outFileUnit,'(1x,a)') 'Please ensure that the intermediate coordinate lies outside the absorbing region !!'
        write(outFileUnit,'(1x,a)') ''

    end subroutine getVabs
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine initAllVabs()
        implicit none

        write(outFileUnit,'(1x,a)') '=====> Vasb information <====='
        write(outFileUnit,'(1x,a)') ''
!> Vabs in r
        allocate(Fabs(IALR%vint))
        write(outFileUnit,'(1x,a)') 'Absorbing potential in r:'
        call getVabs(Vabs%rabs_range,Vabs%Cr,IALR%vint,r_Int,timeStep,Vabs%nrabs,Fabs)
!> Vabs in long-range
        allocate(Flr(IALR%nZ_IALR))
        write(outFileUnit,'(1x,a)') 'Absorbing potential in Z_lr:'
        call getVabs(Vabs%Zlr_range,Vabs%Clr,IALR%nZ_IALR,Z_IALR,timeStep,Vabs%nZlr,Flr)
!> Vabs in asymptotic
        allocate(Fasy(IALR%nZ_IA))
        write(outFileUnit,'(1x,a)') 'Absorbing potential in Z_asy:'
        write(outFileUnit,'(1x,a)') 'Note that the Vabs only works for the channel with (v, j, iPES) /= (v0, j0, initPES)!'
        call getVabs(Vabs%Zasy_range,Vabs%Casy,IALR%nZ_IA,Z_IA,timeStep,Vabs%nZasy,Fasy)
        write(outFileUnit,'(1x,a)') ''

    end subroutine initAllVabs
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine interactionPot(nZ, nr, nth, ZGrid, rGrid, thGrid, type, Vdiag, Voff, VINT, AtD)
        implicit none
        integer, intent(in) :: nZ, nr, nth 
        real(f8), intent(in) :: ZGrid(nZ), rGrid(nr), thGrid(nth) 
        character(len=*), intent(in) :: type 
        real(f8), optional, intent(inout) :: VINT(:,:,:,:)
        real(f8), optional, intent(inout) :: Vdiag(:,:,:,:)
        real(f8), optional, intent(inout) :: Voff(:,:,:,:,:)
        real(f8), optional, intent(inout) :: AtD(:,:,:,:,:)
        real(f8) :: bond(3)
        real(f8) :: adiaV(nPES), diaV(nPES,nPES), AtDMat(nPES,nPES)
        real(f8) :: Z, r, th
        integer :: iZ, ir, ith, iPES, jPES
        character(len=256) :: fileName

        do iZ = 1, nZ
            Z = ZGrid(iZ)
            do ir = 1, nr
                r = rGrid(ir)
                do ith = 1, nth
                !> The grids of Gauss-Legendre are in [-1, 1], need to convert to [0, pi]
                    th = acos(thGrid(ith))
                    call Jacobi2Bond(Z, r, th, atomMass(2), atomMass(3), bond)
                    if (type == 'INT') then 
                        call diagDiaVmat(bond, AtDMat, adiaV) 
                        VINT(:,iZ,ir,ith) = adiaV(:)
                        AtD(:,:,iZ,ir,ith) = AtDMat(:,:)
                    else if (type == 'ALR' .and. present(Vdiag)) then 
                        call POT0(nPES,diaV,bond)
                        do iPES = 1, nPES 
                            Vdiag(iPES,iZ,ir,ith) = diaV(iPES,iPES)-adiaVBC(iPES,ir)
                        end do 
                    end if
                end do
            end do
        end do

        if (type == 'ALR' .and. present(Voff)) then 
            Voff = 0.0_f8
            do iZ = 1, nZ
                Z = ZGrid(iZ)
                do ir = 1, nr
                    r = rGrid(ir)
                    do ith = 1, nth
                        th = acos(thGrid(ith))
                        call Jacobi2Bond(Z, r, th, atomMass(2), atomMass(3), bond)
                        call POT0(nPES,diaV,bond)
                        do iPES = 1, nPES 
                            do jPES = 1, nPES
                                if (iPES == jPES) cycle
                                Voff(iPES,jPES,iZ,ir,ith) = diaV(iPES,jPES)
                            end do
                        end do 
                    end do
                end do
            end do
        end if

!> Type = 'INT', for the vint * nZ_I range
!> Type = 'ALR', for the nPODVR/nDVR * (nZ_asy+nZ_lr) range 
        if (present(VINT)) then
            fileName = 'Vint_'//trim(outfile)//'_'//trim(type)//'.bin'
            call BinReadWrite(fileName, VINT, 'write')
        end if

        if (present(Voff)) then
            fileName = 'Voff_'//trim(outfile)//'_'//trim(type)//'.bin'
            call BinReadWrite(fileName, Voff, 'write')
        end if

        if (present(Vdiag)) then
            fileName = 'Vdiag_'//trim(outfile)//'_'//trim(type)//'.bin'
            call BinReadWrite(fileName, Vdiag, 'write')
        end if

        if (present(AtD)) then
            fileName = 'AtD_'//trim(outfile)//'_'//trim(type)//'.bin'
            call BinReadWrite(fileName, AtD, 'write')
        end if


    end subroutine interactionPot
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine getIntPot()
        implicit none
        integer :: nZ_ALR
        integer :: streamUnit 
        real(f8), allocatable :: Z_ALR(:)
        character(len=3) :: type
        character(len=256) :: fileName
        

!> Interaction potential in the interaction region
        allocate(INT_Vadia(nPES,IALR%nZ_I,IALR%vint,IALR%jint))
        allocate(INT_AtD(nPES,nPES,IALR%nZ_I,IALR%vint,IALR%jint))
        write(outFileUnit,'(1x,a)') '=====> Interaction potential information <====='
        write(outFileUnit,'(1x,a)') ''
        type = 'INT'
        if (trim(potentialType) == 'New') then 
            write(outFileUnit,'(1x,a)') 'Calculating interaction potential in the interaction region...'
            call interactionPot(nZ=IALR%nZ_I, nr=IALR%vint, nth=IALR%jint, ZGrid=Z_I, &
                                rGrid=r_Int, thGrid=intANode, type=type, VINT=INT_Vadia, AtD=INT_AtD)
        else if (trim(potentialType) == 'Read') then 

            fileName = 'Vint_'//trim(outfile)//'_'//trim(type)//'.bin'
            call BinReadWrite(fileName, INT_Vadia, 'read')
            write(outFileUnit,'(1x,a,a)') 'Read INT_Vadia in file: ', trim(fileName)

            fileName = 'AtD_'//trim(outfile)//'_'//trim(type)//'.bin'
            call BinReadWrite(fileName, INT_AtD, 'read')
            write(outFileUnit,'(1x,a,a)') 'Read INT_AtD in file: ', trim(fileName)
        else 
            write(outFileUnit,'(1x,a)') 'Error: unknown potentialType, must Read or New !'
            write(outFileUnit,'(1x,a)') 'POSITION: potent.f90, subroutine getIntPot()'
            stop
        end if

!> Interaction potential in the asymptotic and long-range region
        nZ_ALR = IALR%nZ_IALR-IALR%nZ_I
        allocate(Z_ALR(nZ_ALR))
        Z_ALR(1:nZ_ALR) = Z_IALR(IALR%nZ_I+1:IALR%nZ_IALR)

        allocate(ALR_Vdiag(nPES,nZ_ALR,IALR%nPODVR,IALR%jasy))
        allocate(ALR_Voff(nPES,nPES,nZ_ALR,IALR%vint,IALR%jasy))
        type = 'ALR'
        if (trim(potentialType) == 'New') then 
            write(outFileUnit,'(1x,a)') 'Calculating interaction potential in the asymptotic and long-range region...'
            call interactionPot(nZ=nZ_ALR, nr=IALR%nPODVR, nth=IALR%jasy, ZGrid=Z_ALR, &
                                rGrid=r_Asy, thGrid=asyANode, type=type, Vdiag=ALR_Vdiag)
            call interactionPot(nZ=nZ_ALR, nr=IALR%vint, nth=IALR%jasy, ZGrid=Z_ALR, &
                                rGrid=r_Int, thGrid=asyANode, type=type, Voff=ALR_Voff)
        else if (trim(potentialType) == 'Read') then 

            fileName = 'Vdiag_'//trim(outfile)//'_'//trim(type)//'.bin'
            call BinReadWrite(fileName, ALR_Vdiag, 'read')
            write(outFileUnit,'(1x,a,a)') 'Read ALR_Vdia in file: ', trim(fileName)

            fileName = 'Voff_'//trim(outfile)//'_'//trim(type)//'.bin'
            call BinReadWrite(fileName, ALR_Voff, 'read')
            write(outFileUnit,'(1x,a,a)') 'Read ALR_Voff in file: ', trim(fileName)

        else 
            write(outFileUnit,'(1x,a)') 'Error: unknown potentialType, must Read or New !'
            write(outFileUnit,'(1x,a)') 'POSITION: potent.f90, subroutine getIntPot()'
            stop
        end if
        write(outFileUnit,'(1x,a)') ''

        deallocate(Z_ALR)

    end subroutine getIntPot
!> ------------------------------------------------------------------------------------------------------------------ <!


end module potentMod