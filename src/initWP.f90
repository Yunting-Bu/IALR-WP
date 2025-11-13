module initWPMod
    use machina_basic, only : f8, c8 
    use gPara
    use basisMod
    implicit none

    public
    private :: GaussianWavePacket, Z_initGaussWP
    
contains

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine SF2BFMat(Lmin, Lmax, Kmin, Kmax, j0, Jtot)
        implicit none
        integer, intent(in) :: Lmin, Lmax, Kmin, Kmax, j0, Jtot
        real(f8), external :: CG 
        real(f8) :: delta, fact, CGtmp
        integer :: iL, iK

!> SF to BF transMat
!> BLK = sqrt(2-delta_K0)*sqrt((2L+1)/(2J+1))*<j,K,L,0|J,K>
        do iL = Lmin, Lmax
            do iK = Kmin, Kmax
                delta = merge(1.0_f8, dsqrt(2.0_f8), iK == 0)
                fact = dsqrt((2.0_f8*iL+1)/(2.0_f8*Jtot+1))
                CGtmp = CG(j0,iK,iL,0,Jtot)
                initWP_BLK(iL,iK) = delta*fact*CGtmp
            end do  
        end do
    end subroutine SF2BFMat
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine Z_initGaussWP()
        implicit none
        real(f8) :: normWP
        integer :: iZ 

        allocate(initGaussWP(IALR%nZ_IALR))

!> Gaussian wave-packet in DVR representation
        normWP = 0.0_f8
        do iZ = 1, IALR%nZ_IALR
            initGaussWP(iZ) = GaussianWavePacket(initWP%Ec,initWP%delta,initWP%Zc,massTot,Z_IALR(iZ)) &
                              * dsqrt(Z_IALR(2)-Z_IALR(1))
            normWP = normWP + abs(initGaussWP(iZ))**2
        end do 

        write(outFileUnit,'(1x,a,f15.9)') 'Initial Gaussian wave-packet normalization check: ', normWP

    end subroutine Z_initGaussWP
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine getAdiaInitTotWP()
        implicit none
        integer :: ir, iZ, iPES, ith 
        integer :: Kmax, K, ichnl, nchnl

!> WF at Z in DVR
        call Z_initGaussWP()
!> WF at r in PODVR and theta in FBR
        call lrBC_vibRotThetaWF()

        Kmax = min(initWP%j0, initWP%Jtot)
        nchnl = Kmax - initWP%Kmin + 1
        allocate(initWP_BLK(0:0,initWP%Kmin:Kmax))
        allocate(initAdiaTotWP(nPES,IALR%nZ_IALR,IALR%nr_PODVR,IALR%jasy,nchnl))

        call SF2BFMat(initWP%l0,initWP%l0,initWP%Kmin,Kmax,initWP%j0,initWP%Jtot)
!> Construct initial WP in adiabatic representation
        do iZ = 1, IALR%nZ_IALR
            do ir = 1, IALR%nr_PODVR
                do ith = 1, IALR%jasy 
                    do K = initWP%Kmin, Kmax 
                        ichnl = seqAsy_channel(initWP%v0,initWP%j0,K)
                        initAdiaTotWP(:,iZ,ir,ith,ichnl) = initGaussWP(iZ)*lrWFvjK(ichnl,ir,ith)*initWP_BLK(initWP%l0,K)
                    end do 
                end do 
            end do 
        end do

    end subroutine getAdiaInitTotWP
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
!> Has some problems
    subroutine getDiaInitTotWP()
        implicit none
        integer :: iPES, iZ, ir, ith 
        integer :: Kmax, K, nchnl

        Kmax = min(initWP%j0, initWP%Jtot)
        nchnl = Kmax - initWP%Kmin + 1
        allocate(initDiaTotWP(nPES,IALR%nZ_IALR,IALR%nr_PODVR,IALR%jasy,nchnl))

        do iPES = 1, nPES
            do iZ = 1, IALR%nZ_IALR
                do ir = 1, IALR%nr_PODVR
                    do ith = 1, IALR%jasy 
                        initDiaTotWP(iPES,iZ,ir,ith,:) = initAdiaTotWP(iPES,iZ,ir,ith,:) * asyBC_AtDMat(iPES,ir)
                    end do 
                end do 
            end do 
        end do
    end subroutine getDiaInitTotWP
!> ------------------------------------------------------------------------------------------------------------------ <!
        

!> ------------------------------------------------------------------------------------------------------------------ <!
    complex(c8) function GaussianWavePacket(E0, delta, Z0, mass, Z) result(WP)
        implicit none
        real(f8), intent(in) :: E0, delta, Z0, mass, Z
        real(f8) :: fact, realPart, imgPart 

!> Since we use SO method
!> The WP isn't real packet
        fact = (2.0_f8/(pi*delta*delta))**(0.25)
        realPart = -(Z-Z0)**2 / delta**2
        imgPart = Z * dsqrt(2.0_f8*E0*mass)
        WP = fact*dexp(realPart-img*imgPart)

        return
    end function GaussianWavePacket
!> ------------------------------------------------------------------------------------------------------------------ <!
    
end module initWPMod