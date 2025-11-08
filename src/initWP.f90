module initWPMod
    use machina_basic, only :: f8 
    use gPara
    use basisMod
    implicit none
    
contains

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine SF2BF(Lmin, Lmax, Kmin, Kmax, j0, Jtot, BLK)
        implicit none
        integer, intent(in) :: Lmin, Lmax, Kmin, Kmax, j0, Jtot
        real(f8), intent(inout) :: BLK(:,:)
        real(f8), external :: CG 
        real(f8) :: delta, fact, CGtmp
        integer :: iL, iK

        do iL = Lmin, Lmax
            do iK = Kmin, Kmax
                delta = merge(1.0_f8, dsqrt(2.0_f8), iK == 0)
                fact = dsqrt((2.0_f8*iL+1)/(2.0_f8*Jtot+1))
                CGtmp = CG(j0,iK,iL,0,Jtot)
                BLK(iL,iK) = delta*fact*CGtmp
            end do  
        end do
    end subroutine SF2BF
!> ------------------------------------------------------------------------------------------------------------------ <!
    
end module initWPMod