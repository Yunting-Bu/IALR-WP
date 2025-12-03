module machina_basic
    use iso_fortran_env
!    use MKL_DFTI
    implicit none
    private
    public :: i4, i8
    public :: f4, f8, f16
    public :: c4, c8, c16
    public :: BinReadWrite

    ! kind specifier of 4 byte integer
    integer, parameter :: i4 = int32
    ! kind specifier of 8 byte integer
    integer, parameter :: i8 = int64
    ! kind specifier of 4 byte real
    integer, parameter :: f4 = real32
    ! kind specifier of 8 byte real
    integer, parameter :: f8 = real64
    ! kind specifier of 16 byte real
    integer, parameter :: f16 = real128
    ! kind specifier of 4 byte complex
    integer, parameter :: c4 = real32
    ! kind specifier of 8 byte complex
    integer, parameter :: c8 = real64
    ! kind specifier of 16 byte complex
    integer, parameter :: c16 = real128

!    type, public :: DSTClass
!        type(DFTI_DESCRIPTOR), pointer :: handle => null()
!        integer :: err 
!    contains 
!        procedure :: create
!        procedure :: forward 
!        procedure :: backward 
!        procedure :: destroy
!    end type DSTClass

    interface BinReadWrite
        module procedure realBinReadWrite4D
        module procedure realBinReadWrite5D
    end interface BinReadWrite

contains 


!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine realBinReadWrite4D(file, data, action)
        implicit none
        
        character(len=*), intent(in) :: file
        character(len=*), intent(in) :: action
        real(f8), intent(inout) :: data(:,:,:,:)
        integer :: streamUnit

        if (action == 'read') then
            open(unit=streamUnit, file=trim(file), access='stream', form='unformatted', status='old')
            read(streamUnit) data
            close(streamUnit)
        else if (action == 'write') then
            open(unit=streamUnit, file=trim(file), access='stream', form='unformatted', status='replace')
            write(streamUnit) data
            close(streamUnit)
        else
            write(*,*) 'Error: action must be either "read" or "write".'
        end if  

    end subroutine realBinReadWrite4D
!> ------------------------------------------------------------------------------------------------------------------ <!
    
!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine realBinReadWrite5D(file, data, action)
        implicit none
        
        character(len=*), intent(in) :: file
        character(len=*), intent(in) :: action
        real(f8), intent(inout) :: data(:,:,:,:,:)
        integer :: streamUnit

        if (action == 'read') then
            open(unit=streamUnit, file=trim(file), access='stream', form='unformatted', status='old')
            read(streamUnit) data
            close(streamUnit)
        else if (action == 'write') then
            open(unit=streamUnit, file=trim(file), access='stream', form='unformatted', status='replace')
            write(streamUnit) data
            close(streamUnit)
        else
            write(*,*) 'Error: action must be either "read" or "write".'
        end if  

    end subroutine realBinReadWrite5D
!> ------------------------------------------------------------------------------------------------------------------ <!

!> FFTClass
!> ------------------------------------------------------------------------------------------------------------------ <!
!    subroutine create(this, n)
!        implicit none
!        class(CSTClass) :: this
!        integer, intent(in) :: n

!        this%err = DftiCreateDescriptor(this%handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, n)
!        this%err = DftiSetValue(this%handle, DFTI_PLACEMENT, DFTI_INPLACE)
!        this%err = DftiCommitDescriptor(this%handle)
!    end subroutine create
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
!    subroutine forward(this, DataF)
!        implicit none
!        class(CSTClass) :: this 
!        complex(c8), intent(inout) :: DataF(:)

!        this%err = DftiComputeForward(this%handle, DataF)
!    end subroutine forward
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
!    subroutine backward(this, DataI)
!        implicit none
!        class(CSTClass) :: this 
!        complex(c8), intent(inout) :: DataI(:)

!        this%err = DftiComputeBackward(this%handle, DataI)
!    end subroutine backward
!> ------------------------------------------------------------------------------------------------------------------ <!
    
!> ------------------------------------------------------------------------------------------------------------------ <!
!    subroutine destroy(this)
!        implicit none
!        class(CSTClass) :: this 

!        this%err = DftiFreeDescriptor(this%handle)
!    end subroutine destroy
!> ------------------------------------------------------------------------------------------------------------------ <!

end module machina_basic