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
    end subroutine DVRGrid

    subroutine DVR_TransMatAndKinetic(range, nGrid, mass, Kinetic, TransMat)
        implicit none
        real(f8), intent(in) :: range
        integer, intent(in) :: nGrid
        real(f8), intent(in) :: mass
        real(f8), intent(inout) :: Kinetic(:), TransMat(:,:)
        integer :: i, j
        real(f8) :: fact 

        do i = 1, nGrid
            Kinetic(i) = pi*pi/(4.0_f8*mass*range*range)
        end do 

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
        call DVR_Grid(IALR%nZ_IALR, IALR%Z_range(1), IALR%Z_range(2), IALR%Z_IALR)
        call DVR_TransMatAndKinetic(range_IALR, IALR%nZ_IALR, massTot, IALR%kinZ_IALR, IALR%BZ_IALR)

        range_IA = dR * (IALR%nZ_IA + 1)
        call DVR_TransMatAndKinetic(range_IA, IALR%nZ_IA, massTot, IALR%kinZ_IA, IALR%BZ_IA)

        range_I = dR * (IALR%nZ_I + 1)
        call DVR_TransMatAndKinetic(range_I, IALR%nZ_I, massTot, IALR%kinZ_I, IALR%BZ_I)

        range_r = IALR%r_range(2) - IALR%r_range(1)
        call DVR_Grid(IALR%nr_DVR, IALR%r_range(1), IALR%r_range(2), IALR%r_DVR)
        call DVR_TransMatAndKinetic(range_r, IALR%nr_DVR, massBC, IALR%kin_r, IALR%B_r)

    end subroutine DVR_IALR

    subroutine setChannel()
        implicit none
        integer :: ichnl, v, j, K, Kmax

        ichnl = 0
        do v = 0, ine%vmax
            do j = initWP%jmin, initWP%jmax, initWP%jinc
                Kmax = min(j, initWP%Jtot)
                do K = initWP%Kmin, Kmax
                    ichnl = ichnl + 1
                end do
            end do
        end do
        initWP%nChannels = ichnl

    end subroutine setChannel
    
end module basisMod