module propMod
    use machina_basic, only : f8, c8 
    use gPara
!    use potentMod
!    use basisMod
!    use initWPMod
    implicit none

contains 

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
        allocate(CPDiag(IALR%nZ_IALR, IALR%jint:,initWP%Kmin:KVeryMax))
        allocate(CPKMinus(IALR%nZ_IALR, IALR%jint:,initWP%Kmin:KVeryMax))
        allocate(CPKPlus(IALR%nZ_IALR, IALR%jint:,initWP%Kmin:KVeryMax))
        allocate(fact(IALR%nZ_IALR))

        do iZ = 1, IALR%nZ_IALR
            fact(iZ) = 1.0_f8/(2.0_f8*massTot*Z_IALR(iZ)**2)
        end do

        do j = initWP%jmin, IALR%jint, initWP%jinc
            Kmaxj = min(j, init%Jtot, KVeryMax) 
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

    end function lambdaPlus
!> ------------------------------------------------------------------------------------------------------------------ <!
    