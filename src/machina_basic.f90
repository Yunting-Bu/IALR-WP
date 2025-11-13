module machina_basic
use iso_fortran_env
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

    interface BinReadWrite
        module procedure realBinReadWrite4D
        module procedure realBinReadWrite5D
    end interface BinReadWrite

contains 

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



end module machina_basic