program IALR_WP_SO
    use machina_basic, only : f8 
    use gPara
    use potentMod
    use basisMod
    use initWPMod
    use propMod
    implicit none

!> Test
    call initPara()
    call DVR_IALR()
    call getAdiaInitTotWP()
    call initAllVabs()
    call getIntPot()
    call getEnergyAM()
    call getKinMat()
    call getRotMat()
    call getCPMat()
    
end program IALR_WP_SO
