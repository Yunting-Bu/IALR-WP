module propMod
    use machina_basic, only : f8, c8, FFTClass
    use MKL_DFTI
    use gPara
!    use potentMod
!    use basisMod
!    use initWPMod
    implicit none

    public
    private :: lambdaMinus, lambdaPlus

contains 

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine getKinMat()
        implicit none
        real(f8) :: range_Z, range_r 
        integer :: iZ, ir 

        range_Z = IALR%Z_range(2)-IALR%Z_range(1)
        range_r = IALR%r_range(2)-IALR%r_range(1)

        allocate(Z_KinMat(IALR%nZ_IALR))
        allocate(r_KinMat(IALR%nr_int))

        do iZ = 1, IALR%nZ_IALR
            Z_KinMat(iZ) = (iZ*pi/range_Z)**2/(2.0_f8*massTot)
        end do

        do ir = 1, IALR%nr_int
            r_KinMat(ir) = (ir*pi/range_r)**2/(2.0_f8*massBC)
        end do
    
    end subroutine getKinMat
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine getRotMat()
        implicit none
        real(f8) :: fact, jeigen
        integer :: ir, j
        
        allocate(rotMat(IALR%nr_int,IALR%jint))
        do ir = 1, IALR%nr_int 
            fact = 1.0_f8/(2.0_f8*massTot*r_Int(ir)**2)
            do j = 1, IALR%jint
                jeigen = real(j*(j+1),f8)
                rotMat(ir,j) = fact*jeigen
            end do 
        end do 

    end subroutine getRotMat
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine getCPMat()
        implicit none
        real(f8) ::  delta 
        real(f8), allocatable :: fact(:)
        integer :: iZ, KVeryMax, Kmaxj, K, j 

        KVeryMax = min(initWP%Jtot, IALR%jint)
        allocate(CPDiag(IALR%nZ_IALR, initWP%jmin:IALR%jint,initWP%Kmin:KVeryMax))
        allocate(CPKMinus(IALR%nZ_IALR, initWP%jmin:IALR%jint,initWP%Kmin:KVeryMax))
        allocate(CPKPlus(IALR%nZ_IALR, initWP%jmin:IALR%jint,initWP%Kmin:KVeryMax))
        allocate(fact(IALR%nZ_IALR))

        do iZ = 1, IALR%nZ_IALR
            fact(iZ) = 1.0_f8/(2.0_f8*massTot*Z_IALR(iZ)**2)
        end do

        do j = initWP%jmin, IALR%jint, initWP%jinc
            Kmaxj = min(j, initWP%Jtot, KVeryMax) 
            do K = initWP%Kmin, Kmaxj 
                CPDiag(:,j,K) = (initWP%Jtot*(initWP%Jtot+1.0_f8)+j*(j+1.0_f8)-K**2)*fact(:)
                delta = merge(dsqrt(2.0_f8), 1.0_f8, K == 0)
                if ((K+1 > Kmaxj) .or. (K-1 < initWP%Kmin)) cycle
                CPKPlus(:,j,K+1) = delta*lambdaPlus(j,K)*lambdaPlus(initWP%Jtot,K)*fact(:)
                CPKMinus(:,j,K-1) = delta*lambdaMinus(j,K)*lambdaMinus(initWP%Jtot,K)*fact(:)
            end do 
        end do 
        deallocate(fact)
    end subroutine getCPMat
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
    