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
        integer :: iZ 

        allocate(initGaussWP(IALR%nZ_IALR))

!> Gaussian wave-packet in DVR representation
        do iZ = 1, IALR%nZ_IALR
            initGaussWP(iZ) = GaussianWavePacket(initWP%Ec,initWP%delta,initWP%Zc,massTot,Z_IALR(iZ)) &
                              * dsqrt(Z_IALR(2)-Z_IALR(1))
        end do 


    end subroutine Z_initGaussWP
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine getAdiaInitTotWP()
        implicit none
        real(f8) :: normWPTot, normWPZ
        integer :: ir, iZ, iPES, ith 
        integer :: Kmax, K, ichnl, nchnl

!> WF at Z in DVR
        call Z_initGaussWP()
!> WF at r in PODVR and theta in FBR
        call lrBC_vibRotThetaWF()

        Kmax = min(initWP%j0, initWP%Jtot)
        nchnl = Kmax - initWP%Kmin + 1
        allocate(initWP_BLK(0:0,initWP%Kmin:Kmax))
        allocate(initAdiaTotWP(nPES,IALR%nZ_IALR,IALR%nr_PODVR,nchnl))

        call SF2BFMat(initWP%l0,initWP%l0,initWP%Kmin,Kmax,initWP%j0,initWP%Jtot)
!> Construct initial WP in adiabatic representation
        do iZ = 1, IALR%nZ_IALR
            do ir = 1, IALR%nr_PODVR
                do K = initWP%Kmin, Kmax 
                    ichnl = seq_channel(initWP%v0,initWP%j0,K)
                    initAdiaTotWP(:,iZ,ir,ichnl) = initGaussWP(iZ)*lrWFvjK(ichnl,ir)*initWP_BLK(initWP%l0,K)
                end do 
            end do 
        end do

        write(outFileUnit,'(1x,a)') "=====> Initial WP infromation <====="
        write(outFileUnit,'(1x,a)') ''
        write(outFileUnit,'(1x,a,f15.9,a,f15.9,a,f15.9)') 'Z0 = ', initWP%Zc, ' , E0 = ', initWP%Ec, ' , delta = ', initWP%delta
        normWPTot = sum( abs(initAdiaTotWP)**2 )
        normWPZ = sum( abs(initGaussWP)**2 )
        write(outFileUnit,'(1x,a,f15.9)') 'Initial Gaussian wave-packet normalization check: ', normWPZ
        write(outFileUnit,'(1x,a,f15.9)') 'Initial adiabatic total wave-packet normalization check: ', normWPTot
        write(outFileUnit,'(1x,a)') ''

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
                        initDiaTotWP(iPES,iZ,ir,:) = initAdiaTotWP(iPES,iZ,ir,:) * asyBC_AtDMat(iPES,ir)
                end do 
            end do 
        end do
    end subroutine getDiaInitTotWP
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine getEnergyAM()
        implicit none
        !> For Ricatti-Bessel function 
        real(f8) :: rbZ, rbZU, rbZP, rbZUP 
        real(f8) :: fact, wZ
        !> Phase for trans the Bessel to Hankel
        complex(c8) :: phase
        integer :: iZ, iEtot

        fact = dsqrt(massTot/(2.0_f8*pi))
        phase = exp(-img*pi*0.5_f8*initWP%l0)
        wZ = dsqrt(Z_IALR(2)-Z_IALR(1))

        do iEtot = 1, nEtot
            energyAM(iEtot) = imgZore
            do iZ = 1, IALR%nZ_IALR
                call rbesjy(initWP%l0*1.0_f8, kReact(iEtot)*Z_IALR(iZ), rbZ, rbZP, rbZU, rbZUP)
                energyAM(iEtot) = energyAM(iEtot) + fact/dsqrt(kReact(iEtot)) * &
                                  phase * (-rbZU+img*rbZ) * wZ * initGaussWP(iZ)
            end do 
        end do

        write(outFileUnit,*) ''
        write(outFileUnit,'(1x,a)') '=====> Energy information <====='
        write(outFileUnit,'(1x,a)') ''
        write(outFileUnit,'(1x,a5,3x,a20,3x,a20,3x,a20,3x,a20)') 'nEtot','Etot', 'Wave Number', 'Amplitude (real)', 'Amplitude (img)'
        write(outFileUnit,'(1x,5("-"),3x,20("-"),3x,20("-"),3x,20("-"),3x,20("-"))')
        do iEtot = 1, nEtot
            write(outFileUnit,'(1x,i5,3x,f20.10,3x,f20.10,3x,f20.10,3x,f20.10)') iEtot, Etot(iEtot)*au2ev, kReact(iEtot)*au2ev, &
                                                             energyAM(iEtot)%re*au2ev, energyAM(iEtot)%im*au2ev
        end do 
        write(outFileUnit,*) ''

    end subroutine getEnergyAM
!> ------------------------------------------------------------------------------------------------------------------ <!
                                  
!> ------------------------------------------------------------------------------------------------------------------ <!
    real(f8) function GaussianWavePacket(E0, delta, Z0, mass, Z) result(WP)
        implicit none
        real(f8), intent(in) :: E0, delta, Z0, mass, Z
        real(f8) :: fact, expPart, cosPart 

!> Since we use SO method
!> The WP isn't real packet
        fact = (2.0_f8/(pi*delta*delta))**(0.25)
        expPart = -(Z-Z0)**2 / delta**2
        cosPart = Z * dsqrt(2.0_f8*E0*mass)
        WP = fact * exp(expPart)*cos(cosPart)

        return
    end function GaussianWavePacket
!> ------------------------------------------------------------------------------------------------------------------ <!
    
end module initWPMod