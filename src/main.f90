program IALR_WP_SO
    use machina_basic, only : f8 
    use gPara
    use basisMod
    use initWPMod
    implicit none

!> Test
    call initPara()
    call DVR_IALR()
    call getAdiaInitTotWP()
    
end program IALR_WP_SO
