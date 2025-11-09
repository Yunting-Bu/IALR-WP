module potentMod
    use machina_basic, only : f8 
    use gPara
    implicit none
    
contains

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine Jacobi2Bond(Z, r, theta, massB, massC, bond)
        implicit none
        real(f8), intent(in) :: Z, r, theta 
        real(f8), intent(in) :: massB, massC 
        real(f8), intent(out) :: bond(3)
        real(f8) :: rOB, rOC, th

!> bond(1) = rAB, bond(2) = rBC, bond(3) = rAC
!> theta in degree
        bond(2) = r
        th = theta * pi / 180.0_f8
        rOB = massC / (massB+massC) * r 
        rOC = massB / (massB+massC) * r 
        bond(1) = dsqrt(Z*Z+rOB*rOB-2.0_f8*Z*rOB*dcos(th))
        bond(3) = dsqrt(Z*Z+rOC*rOC+2.0_f8*Z*rOC*dcos(th))

    end subroutine Jacobi2Bond
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine bond2Jacobi(bond, Z, r, theta, massB, massC)
        implicit none
        real(f8), intent(in) :: bond(3)
        real(f8), intent(in) :: massB, massC 
        real(f8), intent(out) :: Z, r, theta 
        real(f8) :: r1, r2, a1, a2

!> theta in degree
        r = bond(2)
        r1 = massB * r / (massB+massC)
        r2 = massC * r / (massB+massC)

        a1 = (bond(1)**2+r**2-bond(3)**2) / (2.0_f8*r*bond(1))
        Z = dsqrt(bond(1)**2+r1**2-2.0_f8*r1*bond(1)*a1)
        a2 = (Z**2+r2**2-bond(3)**2) / (2.0_f8*r2*Z)

        if (a2 > 1.0_f8) a2 = 1.0_f8
        if (a2 < -1.0_f8) a2 = -1.0_f8
        theta = dacos(a2) * 180.0_f8 / pi 
        if (Z == 0.0_f8) theta = 90.0_f8
    
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
    
end module potentMod