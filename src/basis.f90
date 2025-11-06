module basisMod
    use machina_basic , only : f8
    use gPara
    implicit none
   
contains

    subroutine DVR_Grid(nGrid, rMin, rMax, Grid)
        implicit none
        integer, intent(in) :: nGrid
        real(f8), intent(in) :: rMin, rMax
        real(f8), intent(inout) :: Grid(:)
        integer :: i
       
        do i = 1, nGrid
            Grid(i) = rMin + i*(rMax - rMin)/(nGrid + 1)
        end do
    end subroutine DVR_Grid

    subroutine DVR_TransMatAndKinetic(range, nGrid, mass, Kinetic, TransMat)
        implicit none
        real(f8), intent(in) :: range
        integer, intent(in) :: nGrid
        real(f8), intent(in) :: mass
        real(f8), intent(inout) :: Kinetic, TransMat(:,:)
        integer :: i, j
        real(f8) :: fact 

        Kinetic = pi*pi/(4.0_f8*mass*range*range)

        fact = dsqrt(2.0_f8/(nGrid + 1))
        do i = 1, nGrid
            do j = 1, nGrid
                TransMat(j,i) = fact * dsin(pi*i*j/(nGrid + 1))
            end do
        end do
    end subroutine DVR_TransMatAndKinetic

    subroutine DVR_IALR()
        implicit none
        real(f8) :: range_IALR, range_IA, range_I, range_r
        real(f8) :: dR

        range_IALR = IALR%Z_range(2) - IALR%Z_range(1)
        dR = range_IALR / real(IALR%nZ_IALR + 1, f8)
        call DVR_Grid(IALR%nZ_IALR, IALR%Z_range(1), IALR%Z_range(2), Z_IALR)
        call DVR_TransMatAndKinetic(range_IALR, IALR%nZ_IALR, massTot, kinZ_IALR, BZ_IALR)

        range_IA = dR * (IALR%nZ_IA + 1)
        call DVR_TransMatAndKinetic(range_IA, IALR%nZ_IA, massTot, kinZ_IA, BZ_IA)

        range_I = dR * (IALR%nZ_I + 1)
        call DVR_TransMatAndKinetic(range_I, IALR%nZ_I, massTot, kinZ_I, BZ_I)

        range_r = IALR%r_range(2) - IALR%r_range(1)
        call DVR_Grid(IALR%nr_DVR, IALR%r_range(1), IALR%r_range(2), r_DVR)
        call DVR_TransMatAndKinetic(range_r, IALR%nr_DVR, massBC, kin_r, B_r)

    end subroutine DVR_IALR

    subroutine setChannel()
        implicit none
        integer :: ichnl, v, j, K, Kmax

        ichnl = 0
        do v = 0, IALR%vint
            do j = initWP%jmin, IALR%jint, initWP%jinc
                Kmax = min(j, initWP%Jtot)
                do K = initWP%Kmin, Kmax
                    ichnl = ichnl + 1
                end do
            end do
        end do
        nChannels = ichnl

        allocate(qn_channel(nChannels,3))
        ichnl = 0
        do v = 0, IALR%vint
            do j = initWP%jmin, IALR%jint, initWP%jinc
                Kmax = min(j, initWP%Jtot)
                do K = initWP%Kmin, Kmax
                    ichnl = ichnl + 1
                    qn_channel(ichnl,1) = v
                    qn_channel(ichnl,2) = j
                    qn_channel(ichnl,3) = K
                end do
            end do
        end do

    end subroutine setChannel

    subroutine DVR_calc(nDVR, kinetic, Vdiatom, eigVal, eigVec)
!> Calculate the DVR eigenvalues and eigenvectors for a given diatomic potential
!> See J. Chem. Phys. Vol. 96 (3), 1 Feb 1992, pp 1982-1991
        implicit none
        integer, intent(in) :: nDVR
        real(f8), intent(in) :: kinetic
!> Vdiatom: nPES x nDVR
        real(f8), intent(in) :: Vdiatom(:,:)
!> Only calculate initPES vibrational states
        real(f8), intent(inout) :: eigVal(:)
        real(f8), intent(inout) :: eigVec(:,:) 
        integer :: i, j, info, lwork
        real(f8), allocatable :: work(:)

        associate(n => nDVR, C => eigVec, V => Vdiatom)
!> Since diagonal T matrix is C matrix
            do i = 1, n 
                do j = 1, i-1
                    C(i,j) = kinetic * &
                             (dsin(pi * (i-j) / (2.0_f8 * (n+1)))**(-2) - &
                             dsin(pi * (i+j) / (2.0_f8 * (n+1)))**(-2)) * (-1.0_f8) **(i-j)
                end do
                C(i,i) = kinetic * &
                         ((2.0_f8 * (n+1)**2 +1.0_f8) / 3.0_f8 - &
                         dsin(pi * i / (n+1))**(-2)) + V(initWP%initPES,i)
            end do 
            do i = 1, n
                do j = i+1, n
                    C(i,j) = C(j,i)
                end do
            end do
!> Test for the need space of work 
            allocate(work(1))
            call dsyev('V', 'U', n, C, n, eigVal, work, -1, info)
            lwork = int(work(1))
            deallocate(work)

            allocate(work(lwork))
            call dsyev('V', 'U', n, C, n, eigVal, work, lwork, info)
            if (info /= 0) then
                write(outFileUnit,*) "Error in DSYEV of DVR eigen, info =", info
                write(outFileUnit,*) "POSITION: basis.f90, subroutine DVR_calc()"
                stop
            end if
            deallocate(work)


        end associate

    end subroutine DVR_calc

    subroutine PODVR(nPODVR, nDVR, DVRCoeff, DVRGrid, DVREig, vmax, jmax, mass, POGrid, POWF, Evj)
!> Calculate the PODVR basis and eigenvalues/eigenfunctions for given DVR basis
!> See Chemical Physics Letters 1992, 190 (3–4), 225–230.
        implicit none
        integer, intent(in) :: nPODVR, nDVR
        real(f8), intent(in) :: DVRCoeff(:,:)
        real(f8), intent(in) :: DVRGrid(:)
        real(f8), intent(in) :: DVREig(:)
        integer, intent(in) :: vmax, jmax
        real(f8), intent(in) :: mass
        real(f8), intent(inout) :: POGrid(:)
        real(f8), intent(inout) :: POWF(:,:,:)
        real(f8), intent(inout) :: Evj(:,:)
        real(f8) :: POEig(nPODVR)
        real(f8) :: POTransMat(nDVR,nPODVR)
        real(f8) :: Xmat(nPODVR,nPODVR), HRefMat(nPODVR,nPODVR), EMat(nPODVR,nPODVR)
        real(f8), allocatable :: work(:)
        integer :: i, j, l, info, lwork

!> DVR to PODVR transformation matrix
        do i = 1, nDVR
            do j = 1, nPODVR
                POTransMat(i,j) = DVRCoeff(i,j)
            end do
        end do

        Evj = 0.0_f8
!> Phase, the sign is choosen such that the first non-zero element is positive 
!> See J. Chem. Phys. 158, 054801 (2023), Sec. II H
        do i = 1, nPODVR
            call phaseTrans(nDVR,POTransMat(:,i))
        end do

        block 
            real(f8) :: tmp 
!> X_PO = TransMat * X_DVR * TransMat^T
            do i = 1, nPODVR
                do j = i, nPODVR
                    tmp = 0.0_f8
                    do l = 1, nDVR
                        tmp = tmp + POTransMat(l,i) * DVRGrid(l) * POTransMat(l,j)
                    end do
                    Xmat(i,j) = tmp
                    Xmat(j,i) = tmp
                end do
            end do
        end block

        allocate(work(1))
        call dsyev('V', 'U', nPODVR, Xmat, nPODVR, POGrid, work, -1, info)
        lwork = int(work(1))
        deallocate(work)
!> Diagonalize Xmat, Eigenvalues are POGrid, eigenvectors are TransMat named Xmat
        allocate(work(lwork))
        call dsyev('V', 'U', nPODVR, Xmat, nPODVR, POGrid, work, lwork, info)
        if (info /= 0) then
            write(outFileUnit,*) "Error in DSYEV of PODVR, info =", info
            write(outFileUnit,*) "POSITION: basis.f90, subroutine PODVR()"
            stop
        end if
        deallocate(work)

        do i = 1, nPODVR
            call phaseTrans(nPODVR, Xmat(:,i))
        end do

        block 
            real(f8) :: tmp 
            integer :: qn_j, qn_v
!> HRef(i,j) = Sum_l Xmat(l,i) * DvrEig(l) * Xmat(l,j)
            HRefMat = 0.0_f8
            do i = 1, nPODVR
                do j = 1, i 
                    tmp = 0.0_f8
                    do l = 1, nPODVR
                        tmp = tmp + Xmat(l,i) * DVREig(l) * Xmat(l,j)
                    end do
                    HRefMat(i,j) = tmp
                    HRefMat(j,i) = tmp
                end do
            end do

            do qn_j = 0, jmax 
                EMat = HRefMat
                do i = 1, nPODVR
                    EMat(i,i) = EMat(i,i) + real(qn_j*(qn_j+1),f8)/(2.0_f8*mass*POGrid(i)**2)
                end do

                allocate(work(1))
                call dsyev('V', 'U', nPODVR, EMat, nPODVR, POEig, work, -1, info)
                lwork = int(work(1))
                deallocate(work)

                allocate(work(lwork))
                call dsyev('V', 'U', nPODVR, EMat, nPODVR, POEig, work, lwork, info)
                if (info /= 0) then
                    write(outFileUnit,*) "Error in DSYEV of PODVR EMat, info =", info
                    write(outFileUnit,*) "POSITION: basis.f90, subroutine PODVR()"
                    stop
                end if
                deallocate(work)

                do i = 1, nPODVR
                    call phaseTrans(nPODVR, EMat(:,i))
                end do 

                do qn_v = 0, vmax 
                    Evj(qn_v,qn_j) = POEig(qn_v+1)
                end do
                do i = 1, nPODVR
                    do qn_v = 0, vmax 
                        POWF(i,qn_v,qn_j) = EMat(i,qn_v+1)
                    end do
                end do
            end do
        end block
    end subroutine PODVR

    subroutine phaseTrans(n, vec)
        implicit none
        integer, intent(in) :: n
        real(f8), intent(inout) :: vec(:)
        integer :: i
        logical :: found = .false.

        do i = 1, n 
            if (abs(vec(i)) > 1.0e-4_f8) then 
                found = .true.
                exit
            end if
        end do

        if (.not. found) then 
            write(outFileUnit,*) "Error: All vector elements near 0."
            write(outFileUnit,*) "POSITION: basis.f90, subroutine phaseTrans()"
            stop
        end if

        if (vec(i) < 0.0_f8) then
            vec = -vec
        end if

    end subroutine phaseTrans

    
end module basisMod