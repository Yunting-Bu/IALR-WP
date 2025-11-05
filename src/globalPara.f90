module gPara
    use machina_basic , only : f8, c8
    implicit none
    public

    private :: getMass, diatomParity

!> ========== Constants ==========

    real(f8), parameter :: au2ev = 27.211386245988_f8
    real(f8), parameter :: ev2au = 1.0_f8 / au2ev
    real(f8), parameter :: au2cm = 219474.6313705_f8
    real(f8), parameter :: cm2au = 1.0_f8 / au2cm 
    real(f8), parameter :: au2K = 315777.663_f8
    real(f8), parameter :: K2au = 1.0_f8 / au2K
    real(f8), parameter :: amu2au = 1822.888486209_f8
    real(f8), parameter :: pi = 4.0_f8*atan(1.0_f8)
    real(f8), parameter :: one = 1.0_f8
    real(c8), parameter :: img = (0.0_f8, 1.0_f8)
   
!> ========== Global Parameters ==========

    type :: initWP_class
        integer :: j0, v0, Jtot, tpar, jpar
        real(f8) :: Zc, delta, Ec
        integer :: initPES
        integer :: Kmin, jmin, jinc 
    end type initWP_class

    type :: IALR_class
       integer :: nZ_IALR, nZ_IA, nZ_I, nr_DVR 
       integer :: vint, jint, vasy, jasy 
       real(f8) :: Z_range(2), r_range(2)
!> I - interaction region
!> A - asymptotic region
!> LR - long range region
!> kin for Kinetic energy of 1D box 
!> B for DVR-FBR transformation matrix
       real(f8), allocatable :: Z_IALR(:), Z_IA(:), Z_I(:), r_DVR(:)
       real(f8), allocatable :: kinZ_IALR(:), kinZ_IA(:), kinZ_I(:), kin_r(:)
       real(f8), allocatable :: BZ_IALR(:,:), BZ_IA(:,:), BZ_I(:,:), B_r(:,:)
    end type IALR_class

    type :: Vabs_class
        real(f8) :: Casy, Clr, Cr 
        real(f8) :: Zasy_range(2), Zlr_range(2), rabs_range(2)
        real(f8), allocatable :: Zasy(:), ZLr(:), rabs(:)
        real(f8), allocatable :: Fasy(:), Flr(:), Fabs(:) 
    end type Vabs_class

    type :: ine_class
        integer :: vmax, jmax
        real(f8) :: project_Zine
        integer :: T_ine
    end type ine_class

    type :: channel1_class
        integer :: nrp, vpmax, jpmax, jpar, midcoor 
        real(f8) :: rp_range(2), rinf 
        integer :: jmin, jinc
        real(f8), allocatable :: rp(:)
    end type channel1_class

    type :: channel2_class
        integer :: nrp, vpmax, jpmax, jpar, midcoor 
        real(f8) :: rp_range(2), rinf 
        integer :: jmin, jinc
        real(f8), allocatable :: rp(:)
    end type channel2_class

!> ========== Parameters from input file ==========

    integer :: reactChannel, IF_inelastic, nPES, energyUnit
    integer :: T_tot, timeStep, timePrint
    integer :: nEtot  
    real(f8) :: E_range(2), dE, Etot
    real(f8) :: atomMass(3), massBC, massTot
    character(len=2) :: Atoms(3), potentialType
    character(len=50) :: outfile
    real(f8) :: energyUnitTrans(4) = [cm2au, ev2au, K2au, 1.0_f8]

    type(initWP_class) :: initWP
    type(IALR_class) :: IALR 
    type(Vabs_class) :: Vabs 
    type(ine_class) :: ine 
    type(channel1_class) :: channel1
    type(channel2_class) :: channel2 

!> ========== Namelists ==========

    namelist /task/ reactChannel, IF_inelastic, Atoms, nPES, energyUnit, potentialType, outfile
    namelist /energy/ E_range, dE
    namelist /initWP/ initWP
    namelist /IALR/ IALR
    namelist /Vabs/ Vabs
    namelist /propagation/ T_tot, timeStep, timePrint
    namelist /inelastic/ ine
    namelist /channel1/ channel1
    namelist /channel2/ channel2


contains

    subroutine initPara()
        implicit none
        integer :: iEtot, outFileUnit
        real(f8) :: temp

!> ========== Read parameters from input file ==========

        open(unit=233, file="ABC.inf", status='old')
        read(233, nml=task) 
        read(233, nml=energy) 
        read(233, nml=initWP)
        read(233, nml=IALR)
        read(233, nml=Vabs)
        read(233, nml=propagation)
        read(233, nml=channel1)

        if (reactChannel == 2) read(233, nml=channel2)
        if (IF_inelastic == 0) read(233, nml=inelastic) 

        close(233)

!> ========== Construct collision energy ==========

        nEtot = nint((E_range(2) - E_range(1))/dE)
        allocate(Etot(nEtot))
        temp = E_range(1)
        do iEtot = 1, nEtot 
            Etot(iEtot) = temp + dE 
            if (Etot(iEtot) > E_range(2)) then
                Etot(iEtot) = E_range(2)
            end if 
            !> convert to au
            Etot(iEtot) = Etot(iEtot) * energyUnitTrans(energyUnit)
        end do

        if (initWP%tpar == 1) then 
            initWP%Kmin = 0
        else if (initWP%tpar == -1) then 
            initWP%Kmin = 1
        else 
            write(*,*) "Error: total parity tpar must be 1 or -1."
            write(*,*) "POSITION: globalPara.f90, subroutine initPara()"
            stop
        end if

        call getMass()
        call diatomParity()

!> Interaction-Asympotic-Lonng-Range (IALR) range
!>  |----|
!>  |----|----|
!>  |----|----|----|
!>  |----|----|----|
!>  Z1   Z2   Z3   Z4
!>  VAbs: Z3 = Zabs_asy_end, Z4 = Zabs_lr_end

!> Allocate IALR DVR grids, kinetic energy and transMatrix
        allocate(IALR%Z_IALR(IALR%nZ_IALR))
        allocate(IALR%Z_IA(IALR%nZ_IA))
        allocate(IALR%Z_I(IALR%nZ_I))
        allocate(IALR%kinZ_IALR(IALR%nZ_IALR))
        allocate(IALR%kinZ_IA(IALR%nZ_IA))
        allocate(IALR%kinZ_I(IALR%nZ_I))
        allocate(IALR%BZ_IALR(IALR%nZ_IALR, IALR%nZ_IALR))
        allocate(IALR%BZ_IA(IALR%nZ_IA, IALR%nZ_IA))
        allocate(IALR%BZ_I(IALR%nZ_I, IALR%nZ_I))
        allocate(IALR%r_DVR(IALR%nr_DVR))
        allocate(IALR%kin_r(IALR%nr_DVR))
        allocate(IALR%B_r(IALR%nr_DVR, IALR%nr_DVR))
!> Allocate Vabs grids and Vabs
        allocate(Vabs%Zasy(Vabs%nZasy), Vabs%ZLr(Vabs%nZlr), Vabs%rabs(Vabs%nZint))
        allocate(Vabs%Fasy(Vabs%nZasy), Vabs%Flr(Vabs%nZlr), Vabs%Fabs(Vabs%nZint))
!> Allocate channel rp grids
        allocate(channel1%rp(channel1%nrp))
        if (reactChannel == 2) allocate(channel2%rp(channel2%nrp))

!> ========== Output ==========

        open(newunit=outFileUnit, file="output.inf", status='replace')
        write(outFileUnit,'(1x,a)') " ========== Input Parameters =========="
        write(outFileUnit,'(1x,a,a,a)') "Atoms A-B-C: ", Atoms(1), Atoms(2), Atoms(3)
        write(outFileUnit,'(1x,a,2i4)') "v0, j0 = ", initWP%v0, initWP%j0
        write(outFileUnit,'(1x,a,i4)') "Total angular momentum Jtot = ", initWP%Jtot
        write(outFileUnit,'(1x,a,i4)') "Total parity = ", initWP%tpar
        write(outFileUnit,'(1x,a,i4)') "Number of PESs = ", nPES
        


    end subroutine initPara

    subroutine getMass()
        implicit none
        integer, parameter :: nSupportedElements = 10
        character(len=2), parameter :: elemSymbol(nSupportedElements) = &
            ['H ', 'D', 'He', 'Li', 'N ', 'O ', 'F ', 'S', 'Cl', 'Ar']
        real(f8), parameter :: elemMass(nSupportedElements) = &
            [1.00782503223_f8, 2.01410177812_f8, 4.002602_f8, 6.938_f8, 14.00307400443_f8, &
             15.99491461957_f8, 18.99840316273_f8, 32.06_f8, 34.968852682_f8, 39.9623831237_f8]
        
        integer :: iatm, ipos

        do iatm = 1, 3 
            ipos = findloc(elemSymbol, Atoms(iatm))
            if (ipos == 0) then
                write(*,*) "Error: Unsupported element ", Atoms(iatm)
                write(*,*) "POSITION: globalPara.f90, subroutine getMass()"
                stop
            else
                atomMass(iatm) = elemMass(ipos)*amu2au
            end if
        end do

        !> Reduce mass of BC : mB*mC/(mB+mC)
        !> Reduce mass of ABC : mA*(mB+mC)/(mA+mB+mC)
        massBC = atomMass(2)*atomMass(3)/(atomMass(2)+atomMass(3))
        massTot = atomMass(1)*(atomMass(2)+atomMass(3))/(atomMass(1)+atomMass(2)+atomMass(3))

    end subroutine getMass
    
    subroutine diatomParity()
        implicit none
        !> jpar = -1, 0, 1, jmin = 1, 0, 0, jinc = 2, 0, 2
        integer :: jmin(-1:1) = [1, 0, 0]
        integer :: jinc(-1:1) = [2, 0, 2]
        real(f8) :: eps 
        
        eps = 100.0_f8 * epsilon(atomMass(1))

        if (abs(atomMass(2)-atomMass(3)) < eps) then
            !> B and C are identical atoms
            if (initWP%jpar /= 0) then
                write(*,*) "Error: For identical atoms, BC parity jpar must be 0."
                write(*,*) "POSITION: globalPara.f90, subroutine diatomParity()"
                stop
            end if
        end if
        initWP%jmin = jmin(initWP%jpar)
        initWP%jinc = jinc(initWP%jpar)

        if (abs(atomMass(1)-atomMass(2)) < eps) then
            !> A and B are identical atoms
            if (channel1%jpar /= 0) then
                write(*,*) "Error: For identical atoms, AB parity jpar must be 0."
                write(*,*) "POSITION: globalPara.f90, subroutine diatomParity()"
                stop
            end if
        end if
        channel1%jmin = jmin(channel1%jpar)
        channel1%jinc = jinc(channel1%jpar)

        if (reactChannel == 2) then
            if (abs(atomMass(1)-atomMass(3)) < eps) then
                !> A and C are identical atoms
                if (channel2%jpar /= 0) then
                    write(*,*) "Error: For identical atoms, AC parity jpar must be 0."
                    write(*,*) "POSITION: globalPara.f90, subroutine diatomParity()"
                    stop
                end if
            end if
            channel2%jmin = jmin(channel2%jpar)
            channel2%jinc = jinc(channel2%jpar)
        end if

    end subroutine diatomParity

end module gPara
