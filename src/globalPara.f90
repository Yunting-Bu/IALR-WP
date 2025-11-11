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
        integer :: j0, v0, l0, Jtot, tpar, jpar
        real(f8) :: Zc, delta, Ec
        integer :: initPES
        integer :: Kmin, jmin, jinc 
    end type initWP_class

    type :: IALR_class
       integer :: nZ_IALR, nZ_IA, nZ_I, nr_PODVR 
       integer :: vint, jint, vasy, jasy 
       real(f8) :: Z_range(2), r_range(2)
    end type IALR_class

    type :: Vabs_class
        real(f8) :: Casy, Clr, Cr 
        real(f8) :: Zasy_range, Zlr_range, rabs_range
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
    end type channel1_class

    type :: channel2_class
        integer :: nrp, vpmax, jpmax, jpar, midcoor 
        real(f8) :: rp_range(2), rinf 
        integer :: jmin, jinc
    end type channel2_class

!> ========== Parameters from input file ==========

!> Task types and propagation settings
    integer :: reactChannel, nPES, energyUnit
    integer :: timeTot, timeStep, timePrint
!> Energy input
    integer :: nEtot, outFileUnit
    real(f8) :: E_range(2), dE
    real(f8), allocatable :: Etot(:)
    real(f8) :: energyUnitTrans(4) = [cm2au, ev2au, K2au, 1.0_f8]
!> Atoms and masses
    real(f8) :: atomMass(3), massBC, massTot
    character(len=2) :: Atoms(3)
!> Output and other flags
    character(len=3) :: potentialType
    character(len=50) :: outfile
    logical :: IF_inelastic
!> Flux
    character(len=1) :: IDflux 
    real(f8) :: fluxPos
!> Channels
    integer :: nChannels
    integer, allocatable :: qn_channel(:,:)
    integer, allocatable :: seq_channel(:,:,:)
!> I - interaction region
!> A - asymptotic region
!> LR - long range region
!> B for DVR-FBR transformation matrix
    real(f8), allocatable :: Z_IALR(:), Z_IA(:), Z_I(:), r_All(:), r_Asy(:)
    real(f8), allocatable :: BZ_IALR(:,:), BZ_IA(:,:), BZ_I(:,:), B_rAll(:,:), B_rAsy(:,:)
    real(f8), allocatable :: r_PODVR(:)
!> Grids and weights for K independent Gauss-Legendre quadrature
    real(f8), allocatable :: asyANode(:), asyAWeight(:)
!> Vib-rotational basis and K-independent Guass-Legrendre basis of BC in asymptotic range
    real(f8), allocatable :: asyWFvjK(:,:,:)
    real(f8), allocatable :: asyBC_Evj(:,:)
    real(f8), allocatable :: asyBC_POWF(:,:,:)
!> Vib-rotational basis and K-independent Guass-Legrendre basis of BC in long-range
!> Only have one state (v0, j0)
    real(f8) :: lrBC_Evj
    real(f8), allocatable :: lrWFvjK(:,:,:)
    real(f8), allocatable :: lrBC_POWF(:)
!> Initial Gaussian-shape WP and initial total WP
    complex(c8), allocatable :: initGaussWP(:)
    complex(c8), allocatable :: initTotWP(:,:,:,:,:)
    real(f8), allocatable :: initWP_BLK(:,:)
!> Vabs
    real(f8), allocatable :: Fasy(:), Flr(:), Fabs(:) 
!> Product channel grids
    real(f8), allocatable :: rp1(:)
    real(f8), allocatable :: rp2(:)
!> Tpyes declarations
    type(initWP_class) :: initWP
    type(IALR_class) :: IALR 
    type(Vabs_class) :: abs 
    type(ine_class) :: ine 
    type(channel1_class) :: channel1
    type(channel2_class) :: channel2 

!> ========== Namelists ==========

    namelist /task/ reactChannel, IF_inelastic, IDflux, fluxPos, Atoms, nPES, energyUnit, potentialType, outfile
    namelist /energy/ E_range, dE
    namelist /initWavePacket/ initWP
    namelist /IALRset/ IALR
    namelist /VabsAndDump/ abs
    namelist /propagation/ timeTot, timeStep, timePrint
    namelist /inelastic/ ine
    namelist /productChannel1/ channel1
    namelist /productChannel2/ channel2


contains

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine initPara()
        implicit none
        integer :: iEtot, inpFileUnit
        real(f8) :: temp

!> ========== Read parameters from input file ==========

        open(unit=inpFileUnit, file="ABC.inf", status='old')
        read(inpFileUnit, nml=task) 
        read(inpFileUnit, nml=initWavePacket)
        read(inpFileUnit, nml=IALRset)
        read(inpFileUnit, nml=propagation)
        read(inpFileUnit, nml=energy) 
        read(inpFileUnit, nml=VabsAndDump)
        read(inpFileUnit, nml=productChannel1)

        rewind(inpFileUnit)

        if (reactChannel == 2) read(inpFileUnit, nml=productChannel2)
        rewind(inpFileUnit)
        if (IF_inelastic) read(inpFileUnit, nml=inelastic) 

        close(inpFileUnit)
        
        open(newunit=outFileUnit, file=trim(outfile)//".out", status='replace')
        !> Should add more output files later

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
            write(outFileUnit,*) "Error: total parity tpar must be 1 or -1."
            write(outFileUnit,*) "POSITION: globalPara.f90, subroutine initPara()"
            stop
        end if

        if (initWP%tpar /= (-1)**(initWP%Jtot+initWP%j0+initWP%l0)) then 
            write(outFileUnit,*) "Error: total pairty incoorent with l0."
            write(outFileUnit,*) "POSITION: globalPara.f90, subroutine initPara()"
            stop
        end if

        if (initWP%l0 > initWP%j0+initWP%Jtot) then 
            write(outFileUnit,*) "Error: l0 > Jtot + j0."
            write(outFileUnit,*) "POSITION: globalPara.f90, subroutine initPara()"
            stop
        end if

        if (initWP%l0 < abs(initWP%j0-initWP%Jtot)) then 
            write(outFileUnit,*) "Error: l0 < |Jtot - j0|."
            write(outFileUnit,*) "POSITION: globalPara.f90, subroutine initPara()"
            stop
        end if

        call getMass()
        call diatomParity()

!> Allocate channel rp grids
        allocate(rp1(channel1%nrp))
        if (reactChannel == 2) allocate(rp2(channel2%nrp))

!> ========== Output ==========

        write(outFileUnit,'(1x,a)') " ============= Input Parameters ============="
        write(outFileUnit,'(1x,a,a,a)') "Reaction Channel: ", &
            merge("A + BC -> AB + C           ", "A + BC -> AB + C and AC + B", reactChannel==1)
        write(outFileUnit,'(1x,6a)') "Atoms A, B, C: ", Atoms(1),', ', Atoms(2),', ', Atoms(3)
        write(outFileUnit,'(1x,a,f15.9)') "Reduced Masses of BC (a.u.):", massBC
        write(outFileUnit,'(1x,a,f15.9)') "Reduced Masses of ABC (a.u.):", massTot
        write(outFileUnit,'(1x,a,3i4)') "v0, j0, l0 : ", initWP%v0, initWP%j0, initWP%l0
        write(outFileUnit,'(1x,a,i4)') "Total angular momentum Jtot: ", initWP%Jtot
        write(outFileUnit,'(1x,a,i4)') "Total parity: ", initWP%tpar
        write(outFileUnit,'(1x,a,i4)') "Number of PESs: ", nPES
        


    end subroutine initPara
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine getMass()
        implicit none
        integer, parameter :: nSupportedElements = 10
        character(len=2), parameter :: elemSymbol(nSupportedElements) = &
            ['H ', 'D ', 'He', 'Li', 'N ', 'O ', 'F ', 'S ', 'Cl', 'Ar']
        real(f8), parameter :: elemMass(nSupportedElements) = &
            [1.00782503223_f8, 2.01410177812_f8, 4.002602_f8, 6.938_f8, 14.00307400443_f8, &
             15.99491461957_f8, 18.99840316273_f8, 32.06_f8, 34.968852682_f8, 39.9623831237_f8]
        
        integer :: iatm, ipos

        do iatm = 1, 3 
            ipos = findloc(elemSymbol, Atoms(iatm), dim=1)
            if (ipos == 0) then
                write(outFileUnit,*) "Error: Unsupported element ", Atoms(iatm)
                write(outFileUnit,*) "POSITION: globalPara.f90, subroutine getMass()"
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
!> ------------------------------------------------------------------------------------------------------------------ <!
    
!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine diatomParity()
        implicit none
        !> jpar = -1, 0, 1, jmin = 1, 0, 0, jinc = 2, 1, 2
        integer :: jmin(-1:1) = [1, 0, 0]
        integer :: jinc(-1:1) = [2, 1, 2]
        real(f8) :: eps 
        
        eps = 100.0_f8 * epsilon(atomMass(1))

        if (abs(atomMass(2)-atomMass(3)) < eps) then
            !> B and C are identical atoms
            if (initWP%jpar /= 0) then
                write(outFileUnit,*) "Error: For identical atoms, BC parity jpar must be 0."
                write(outFileUnit,*) "POSITION: globalPara.f90, subroutine diatomParity()"
                stop
            end if
        end if
        initWP%jmin = jmin(initWP%jpar)
        initWP%jinc = jinc(initWP%jpar)

        if (abs(atomMass(1)-atomMass(2)) < eps) then
            !> A and B are identical atoms
            if (channel1%jpar /= 0) then
                write(outFileUnit,*) "Error: For identical atoms, AB parity jpar must be 0."
                write(outFileUnit,*) "POSITION: globalPara.f90, subroutine diatomParity()"
                stop
            end if
        end if
        channel1%jmin = jmin(channel1%jpar)
        channel1%jinc = jinc(channel1%jpar)

        if (reactChannel == 2) then
            if (abs(atomMass(1)-atomMass(3)) < eps) then
                !> A and C are identical atoms
                if (channel2%jpar /= 0) then
                    write(outFileUnit,*) "Error: For identical atoms, AC parity jpar must be 0."
                    write(outFileUnit,*) "POSITION: globalPara.f90, subroutine diatomParity()"
                    stop
                end if
            end if
            channel2%jmin = jmin(channel2%jpar)
            channel2%jinc = jinc(channel2%jpar)
        end if

    end subroutine diatomParity
!> ------------------------------------------------------------------------------------------------------------------ <!

end module gPara
