module basisMod
    use machina_basic, only : f8
    use gPara
    use potentMod
    implicit none
    public

    private :: setChannel, asyBC_vibRotThetaWF
   
contains

!> ------------------------------------------------------------------------------------------------------------------ <!
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
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
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
            do j = 1, i
                TransMat(i,j) = fact * dsin(pi*i*j/(nGrid + 1))
                TransMat(j,i) = TransMat(i,j)
            end do
        end do
    end subroutine DVR_TransMatAndKinetic
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine DVR_IALR()
        implicit none
        real(f8) :: range_IALR, range_IA, range_I, range_r, range_rA
        real(f8) :: dR

!> Interaction-Asympotic-Lonng-Range (IALR) range
!> r3   |----|
!> r2   |----|----|
!>      |----|----|----|
!> r1   |----|----|----|
!>      Z1   Z2   Z3   Z4
!>  nDVR in range [r1, r3] = vint, nDVR in range [r1, r2] = vasy 
!>  VAbs: Z3 = Zabs_asy_end, Z4 = Zabs_lr_end

!> Allocate IALR DVR grids, kinetic energy and transMatrix
        allocate(Z_IALR(IALR%nZ_IALR))
        allocate(Z_IA(IALR%nZ_IA))
        allocate(Z_I(IALR%nZ_I))
        allocate(BZ_IALR(IALR%nZ_IALR, IALR%nZ_IALR))
        allocate(BZ_IA(IALR%nZ_IA, IALR%nZ_IA))
        allocate(BZ_I(IALR%nZ_I, IALR%nZ_I))
        allocate(r_Asy(IALR%vasy), r_All(IALR%vint))
        allocate(B_rAll(IALR%vint, IALR%vint))
        allocate(B_rAsy(IALR%vasy, IALR%vasy))

        range_IALR = IALR%Z_range(2) - IALR%Z_range(1)
        dR = range_IALR / real(IALR%nZ_IALR + 1, f8)
        call DVR_Grid(IALR%nZ_IALR, IALR%Z_range(1), IALR%Z_range(2), Z_IALR)
        call DVR_TransMatAndKinetic(range_IALR, IALR%nZ_IALR, massTot, kinZ_IALR, BZ_IALR)

!> ======== DVR calculate =========
        range_IA = dR * (IALR%nZ_IA + 1)
        call DVR_TransMatAndKinetic(range_IA, IALR%nZ_IA, massTot, kinZ_IA, BZ_IA)

        range_I = dR * (IALR%nZ_I + 1)
        call DVR_TransMatAndKinetic(range_I, IALR%nZ_I, massTot, kinZ_I, BZ_I)

        range_r = IALR%r_range(2) - IALR%r_range(1)
        dR = range_r / real(IALR%vint + 1, f8)
        call DVR_Grid(IALR%vasy, IALR%r_range(1), IALR%r_range(2), r_All)
        call DVR_TransMatAndKinetic(range_r, IALR%vint, massBC, kin_rAll, B_rAll)

        range_rA = dR * (IALR%vasy + 1)
        call DVR_TransMatAndKinetic(range_rA, IALR%vasy, massBC, kin_rAsy, B_rAsy)
        r_Asy(:) = r_All(1:IALR%vasy)

    end subroutine DVR_IALR
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine asyBC_vibRotThetaWF()
        implicit none
        real(f8), parameter :: longDistance = 30.0_f8
        real(f8) :: bond(3), AtDMat(nPES,nPES), Vadia(nPES)
        real(f8), allocatable :: adiaV(:,:)
        real(f8), allocatable :: DVREig(:)
        real(f8), allocatable :: DVRWF(:,:)
        real(f8), allocatable :: asyBC_AtDMat(:,:)
        real(f8) :: wA, cth, YjK, normWF
        real(f8), external :: spgndr
        integer :: ir, v, j, K, i, ith
        integer :: normWFunit

!> Allocate asympotic range angular quadrature grids and weights
        allocate(asyANode(IALR%jasy))
        allocate(asyAWeight(IALR%jasy))
!> Allocate wave function
        allocate(asyWFvjK(nChannels,IALR%nr_PODVR,IALR%jasy))
!> Channel set as (v,j,K)
        call setChannel()
        call getANodeAndWeight(initWP%jpar, jasy, asyANode, asyAWeight)

!> Allocate PODVR array
        allocate(adiaV(nPES,IALR%vasy))
        allocate(DVREig(IALR%vasy), DVRWF(IALR%vasy,IALR%vasy))
        allocate(asyBC_AtDMat(nPES,IALR%vasy))
        allocate(asyBC_Evj(0:IALR%vasy,0:IALR%jasy))
        allocate(r_PODVR(IALR%nr_PODVR))
        allocate(asyBC_POWF(IALR%nr_PODVR,0:IALR%vasy,0:IALR%jasy))

!> In this situation, bond(2) is the length of BC
!> You should check the PES interface!
        bond(1) = longDistance
        bond(3) = longDistance
        do ir = 1, IALR%vasy
            bond(2) = r_Asy(ir)
            call diagDiaVmat(bond,AtDMat,Vadia)
            asyBC_AtDMat(:,ir) = AtDMat(:)
            adiaV(:,ir) = Vadia
        end do 
        call DVR_calc(IALR%vasy,kin_rAsy,adiaV,DVREig,DVRWF)
        call PODVR(IALR%nr_PODVR,IALR%vasy,DVRWF,r_Asy,DVREig,IALR%vasy,IALR%jasy,massBC,r_PODVR,asyBC_POWF,asyBC_Evj)

        normWF = 0.0_f8
        do ir = 1, IALR%vasy 
            normWF = normWF + asyBC_POWF(i,initWP%v0,initWP%j0)**2
        end do 

        write(outFileUnit,'(1x,a)') '==================================================================================='
        write(outFileUnit,'(1x,a,2i,a)') 'Initial ro-vibrational energy of (v0, j0) = ', initWP%v0, initWP%j0, ' state.'
        write(outFileUnit,'(1x,a,f15.9,a)') 'Evj of BC = ', asyBC_Evj(initWP%v0,initWP%j0)*au2ev, ' eV.'
        write(outFileUnit,'(1x,a,f15.9,a)') 'Evj of BC = ', asyBC_Evj(initWP%v0,initWP%j0)*au2cm, ' cm-1.'
        write(outFileUnit,'(1x,a)') 'Please check the energy!'
        write(outFileUnit,*) ''
        write(outFileUnit,'(1x,a,f15.0)') 'Initial ro-vibrational wave function normalization check: ', normWF
        write(outFileUnit,'(1x,a)') '==================================================================================='

        do i = 1, nChannels
            v = qn_channel(i,1)
            j = qn_channel(i,2)
            K = qn_channel(i,3)
            do ith = 1, IALR%jasy
                cth = asyANode(ith)
                wA = dsqrt(asyAWeight(ith))
                YjK = spgndr(j,K,cth)
                asyWFvjK(i,:,ith) = wA*YjK*asyBC_POWF(:,v,j)
            end do 
        end do

    end subroutine asyBC_vibRotThetaWF
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine lrBC_vibRotThetaWF()
        implicit none
        integer :: ichnl, K, Kmax, nchnl 

!> The long range WF is based on asyWF
        call asyBC_vibRotThetaWF()
        
        allocate(lrBC_POWF(IALR%nr_PODVR))

        lrBC_Evj = asyBC_Evj(initWP%v0,initWP%j0)
        lrBC_POWF(:) = asyBC_POWF(:,initWP%v0,initWP%j0)

        Kmax = min(initWP%Jtot,initWP%j0)
        nchnl = Kmax - initWP%Kmin + 1
        allocate(lrWFvjK(nchnl,IALR%nr_PODVR,IALR%jasy))
        do K = initWP%Kmin, Kmax 
            ichnl = seq_channel(initWP%v0,initWP%j0,K)
            lrWFvjK(ichnl,:,:) = asyWFvjK(ichnl,:,:)
        end do 
    
    end subroutine lrBC_vibRotThetaWF
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine getANodeAndWeight(jpar, nth, nodes, weights)
        implicit none
        integer, intent(in) :: jpar, nth
        real(f8), intent(inout) :: nodes(:), weights(:)
        integer :: i, ntmp 

!> Take advantage of the symmetry of potential the grid used for angle could be used half.
        if (jpar /= 0) then 
            ntmp = 2*nth 
        else
            ntmp = nth 
        end if 

!> First call ntmp = 2*nth for the system has symmetry.
        block 
            real(f8) :: tmpNode(ntmp), tmpWeight(ntmp)
            
            call GAULEG(-1.0_f8, 1.0_f8, tmpNode, tmpWeight, ntmp)
!> Take half of grid
            do i = 1, nth 
                nodes(i) = tmpNode(i)
                if (jpar /= 0) then
                    weights(i) = 2.0_f8*tmpWeight(i)
                else 
                    weights(i) = tmpWeight(i)
                end if 
            end do 
        end block

    end subroutine getANodeAndWeight
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setChannel()
        implicit none
        integer :: ichnl, v, j, K, Kmax

        ichnl = 0
        do v = 0, IALR%vasy
            do j = initWP%jmin, IALR%jasy, initWP%jinc
                Kmax = min(j, initWP%Jtot)
                do K = initWP%Kmin, Kmax
                    ichnl = ichnl + 1
                end do
            end do
        end do
        nChannels = ichnl

        allocate(qn_channel(nChannels,3))
        ichnl = 0
        do v = 0, IALR%vasy
            do j = initWP%jmin, IALR%jasy, initWP%jinc
                Kmax = min(j, initWP%Jtot)
                do K = initWP%Kmin, Kmax
                    ichnl = ichnl + 1
                    qn_channel(ichnl,1) = v
                    qn_channel(ichnl,2) = j
                    qn_channel(ichnl,3) = K
                end do
            end do
        end do

        Kmax = min(initWP%Jtot,IALR%jasy)
        allocate(seq_channel(0:IALR%vasy,initWP%jmin:IALR%jasy,initWP%Kmin:Kmax))
        seq_channel = -1
        do ichnl = 1, nChannels
            v = qn_channel(ichnl,1)
            j = qn_channel(ichnl,2)
            K = qn_channel(ichnl,3)
            seq_channel(v,j,K) = ichnl
        end do

    end subroutine setChannel
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
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
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
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
!> HRef(i,j) = Sum_l Xmat(l,i) * DVREig(l) * Xmat(l,j)
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
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
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
!> ------------------------------------------------------------------------------------------------------------------ <!

    
end module basisMod