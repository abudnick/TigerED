! Soil biogeochemistry module for ED2.

Module soil_bgc
  implicit none

  !-------------------------------------------------
  ! 1. Define all soil biogeochemistry variables.
  !-------------------------------------------------
  type sbgc_vars

     !====================
     ! Soil Carbon
     !====================
     
     ! Soil carbon, fast pool (kgC/m2)
     real, allocatable, dimension(:) :: fast_soil_C 
     
     ! Soil carbon, slow pool (kgC/m2)
     real, allocatable, dimension(:)  :: slow_soil_C 
     
     ! Soil carbon, structural pool (kgC/m2)
     real, allocatable, dimension(:)  :: struct_soil_C
     
     !---------------
     ! This one is a legacy variable.  Probably can delete this.
     ! Soil lignin , structural pool (kgC/m2)
     real, allocatable, dimension(:)  :: struct_soil_L
     !---------------
     
     !====================
     ! Soil Nitrogen
     !====================
     
     ! Soil nitrogen, fast pool (kgN/m2)
     real, allocatable, dimension(:) :: fast_soil_N  
     
     ! Soil nitrogen, structural pool (kgN/m2)
     real, allocatable, dimension(:) :: struct_soil_N  
     
     ! Soil nitrogen, mineralized pool (kgN/m2)
     real, allocatable, dimension(:) :: miner_soil_N  
     
     ! Soil nitrogen, slow pool (kgN/m2)
     real, allocatable, dimension(:) :: slow_soil_N  
     
     
     !====================
     ! Soil Phosphorus
     !====================
     
     ! Soil phosphorus, fast pool (kgP/m2)
     real, allocatable, dimension(:)  :: fast_soil_P  
     ! Soil phosphorus, structural pool (kgP/m2)
     real, allocatable, dimension(:)  :: struct_soil_P  
     ! Soil phosphorus, mineralized pool (kgP/m2)
     real, allocatable, dimension(:)  :: miner_soil_P 
     ! Soil phosphorus, slow pool (kgP/m2)
     real, allocatable, dimension(:)  :: slow_soil_P 
  
  end type sbgc_vars

Contains

  !-------------------------------------------------
  ! 2. Allocate variables.  Required in all cases.
  !-------------------------------------------------
  subroutine allocate_sitetype_sbgc(sbgc, npatches)
    implicit none

    type(sbgc_vars) :: sbgc
    integer, intent(in) :: npatches

    allocate(sbgc%fast_soil_C(npatches))
    allocate(sbgc%struct_soil_C(npatches))
    allocate(sbgc%slow_soil_C(npatches))
    allocate(sbgc%struct_soil_L(npatches))
    allocate(sbgc%fast_soil_N(npatches))
    allocate(sbgc%struct_soil_N(npatches))
    allocate(sbgc%miner_soil_N(npatches))
    allocate(sbgc%slow_soil_N(npatches))
    allocate(sbgc%fast_soil_P(npatches))
    allocate(sbgc%struct_soil_P(npatches))
    allocate(sbgc%miner_soil_P(npatches))
    allocate(sbgc%slow_soil_P(npatches))

    return
  end subroutine allocate_sitetype_sbgc

  !===================================================================

  !-------------------------------------------------
  ! 3. Deallocate variables.  Required in all cases.
  !-------------------------------------------------
  subroutine deallocate_sitetype_sbgc(sbgc)
    implicit none

    type(sbgc_vars) :: sbgc

    deallocate(sbgc%fast_soil_C)
    deallocate(sbgc%struct_soil_C)
    deallocate(sbgc%slow_soil_C)
    deallocate(sbgc%struct_soil_L)
    deallocate(sbgc%fast_soil_N)
    deallocate(sbgc%struct_soil_N)
    deallocate(sbgc%miner_soil_N)
    deallocate(sbgc%slow_soil_N)
    deallocate(sbgc%fast_soil_P)
    deallocate(sbgc%struct_soil_P)
    deallocate(sbgc%miner_soil_P)
    deallocate(sbgc%slow_soil_P)

    return
  end subroutine deallocate_sitetype_sbgc

  !===================================================================

  !-------------------------------------------------
  ! 4. Enable copying of patches.  Required in all cases.
  !-------------------------------------------------
  subroutine copy_sitetype_sbgc(osbgc, isbgc, opa, ipa)
    implicit none

    type(sbgc_vars) :: osbgc, isbgc
    integer, intent(in) :: opa,ipa

    osbgc%fast_soil_C(opa) = isbgc%fast_soil_C(ipa)
    osbgc%slow_soil_C(opa) = isbgc%slow_soil_C(ipa)
    osbgc%struct_soil_C(opa) = isbgc%struct_soil_C(ipa)
    osbgc%struct_soil_L(opa) = isbgc%struct_soil_L(ipa)
    osbgc%fast_soil_N(opa) = isbgc%fast_soil_N(ipa)
    osbgc%slow_soil_N(opa) = isbgc%slow_soil_N(ipa)
    osbgc%struct_soil_N(opa) = isbgc%struct_soil_N(ipa)
    osbgc%miner_soil_N(opa) = isbgc%miner_soil_N(ipa)
    osbgc%fast_soil_P(opa) = isbgc%fast_soil_P(ipa)
    osbgc%slow_soil_P(opa) = isbgc%slow_soil_P(ipa)
    osbgc%struct_soil_P(opa) = isbgc%struct_soil_P(ipa)
    osbgc%miner_soil_P(opa) = isbgc%miner_soil_P(ipa)

    return
  end subroutine copy_sitetype_sbgc

  !===================================================================

  !-------------------------------------------------
  ! 5. Enable selective copying of patches.  Very good to have.  Required for disturbance.
  !-------------------------------------------------
  subroutine copy_sitetype_mask_sbgc(masksz,logmask,sbgc_out,sbgc_in)
    implicit none

    integer, intent(in) :: masksz
    logical, dimension(masksz), intent(in) :: logmask
    integer :: inc, ipa
    type(sbgc_vars) :: sbgc_in
    type(sbgc_vars) :: sbgc_out

    inc = 0
    do ipa = 1,masksz
       if(logmask(ipa))then
          inc = inc + 1
          sbgc_out%fast_soil_C(inc) = sbgc_in%fast_soil_C(ipa)
          sbgc_out%struct_soil_C(inc) = sbgc_in%struct_soil_C(ipa)
          sbgc_out%slow_soil_C(inc) = sbgc_in%slow_soil_C(ipa)
          sbgc_out%struct_soil_L(inc) = sbgc_in%struct_soil_L(ipa)

          sbgc_out%fast_soil_N(inc) = sbgc_in%fast_soil_N(ipa)
          sbgc_out%struct_soil_N(inc) = sbgc_in%struct_soil_N(ipa)
          sbgc_out%slow_soil_N(inc) = sbgc_in%slow_soil_N(ipa)
          sbgc_out%miner_soil_N(inc) = sbgc_in%miner_soil_N(ipa)

          sbgc_out%fast_soil_P(inc) = sbgc_in%fast_soil_P(ipa)
          sbgc_out%struct_soil_P(inc) = sbgc_in%struct_soil_P(ipa)
          sbgc_out%slow_soil_P(inc) = sbgc_in%slow_soil_P(ipa)
          sbgc_out%miner_soil_P(inc) = sbgc_in%miner_soil_P(ipa)
       endif
    enddo

    return
  end subroutine copy_sitetype_mask_sbgc

  !===========================================================

  !-------------------------------------------------
  ! 6. Specify which files to write variables to.  Required if you want to write output.
  !-------------------------------------------------
  subroutine filltab_sitetype_sbgc(nvar, npts, sbgc, igr, init, paglob_id, &
       var_len, var_len_global, max_ptrs)

    use ed_var_tables, only: vtable_edio_r, metadata_edio

    implicit none

    integer, intent(in) :: npts, igr, init, paglob_id, var_len,   &
         max_ptrs, var_len_global
    type(sbgc_vars) :: sbgc
    integer, intent(inout) :: nvar

    nvar=nvar+1
    call vtable_edio_r(npts,sbgc%fast_soil_C(1:npts),  &
         nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,'FAST_SOIL_C :31:hist:year') 
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    
    nvar=nvar+1
    call vtable_edio_r(npts,sbgc%struct_soil_C(1:npts),  &
         nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,'STRUCT_SOIL_C :31:hist:year') 
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    
    nvar=nvar+1
    call vtable_edio_r(npts,sbgc%slow_soil_C(1:npts),  &
         nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,'SLOW_SOIL_C :31:hist:year') 
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    
    nvar=nvar+1
    call vtable_edio_r(npts,sbgc%struct_soil_L(1:npts),  &
         nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,'STRUCT_SOIL_L :31:hist:year') 
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    
    nvar=nvar+1
    call vtable_edio_r(npts,sbgc%fast_soil_N(1:npts),  &
         nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,'FAST_SOIL_N :31:hist:year') 
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    
    nvar=nvar+1
    call vtable_edio_r(npts,sbgc%struct_soil_N(1:npts),  &
         nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,'STRUCT_SOIL_N :31:hist:year') 
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    
    nvar=nvar+1
    call vtable_edio_r(npts,sbgc%miner_soil_N(1:npts),  &
         nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,'MINER_SOIL_N :31:hist:year') 
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    
    nvar=nvar+1
    call vtable_edio_r(npts,sbgc%slow_soil_N(1:npts),  &
         nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,'SLOW_SOIL_N :31:hist:year') 
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    
    
    nvar=nvar+1
    call vtable_edio_r(npts,sbgc%fast_soil_P(1:npts),  &
         nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,'FAST_SOIL_P :31:hist:year') 
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    
    nvar=nvar+1
    call vtable_edio_r(npts,sbgc%struct_soil_P(1:npts),  &
         nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,'STRUCT_SOIL_P :31:hist:year') 
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    
    nvar=nvar+1
    call vtable_edio_r(npts,sbgc%miner_soil_P(1:npts),  &
         nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,'MINER_SOIL_P :31:hist:year') 
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    
    nvar=nvar+1
    call vtable_edio_r(npts,sbgc%slow_soil_P(1:npts),  &
         nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,'SLOW_SOIL_P :31:hist:year') 
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    return
  end subroutine filltab_sitetype_sbgc

  !===========================================================

  !-------------------------------------------------
  ! 7. Required if you want disturbance.
  !-------------------------------------------------
  subroutine init_disturb_patch_sbgc(sbgc, np)
    implicit none

    type(sbgc_vars) :: sbgc
    integer, intent(in) :: np

    sbgc%fast_soil_C(np) = 0.0
    sbgc%slow_soil_C(np) = 0.0
    sbgc%struct_soil_C(np) = 0.0
    sbgc%struct_soil_L(np) = 0.0
    sbgc%miner_soil_N(np) = 0.0
    sbgc%fast_soil_N(np) = 0.0
    sbgc%struct_soil_N(np) = 0.0
    sbgc%slow_soil_N(np) = 0.0
    sbgc%fast_soil_P(np) = 0.0
    sbgc%struct_soil_P(np) = 0.0
    sbgc%miner_soil_P(np) = 0.0
    sbgc%slow_soil_P(np) = 0.0

    return
  end subroutine init_disturb_patch_sbgc

  !===========================================================

  !-------------------------------------------------
  ! 8. Required if you want disturbance.
  !-------------------------------------------------
  subroutine inc_patch_vars_sbgc(sbgc, np, cp, area_fac)
    implicit none

    type(sbgc_vars) :: sbgc
    integer, intent(in) :: np, cp
    real, intent(in) :: area_fac

    sbgc%fast_soil_C(np) = sbgc%fast_soil_C(np) + sbgc%fast_soil_C(cp) * area_fac
    sbgc%slow_soil_C(np) = sbgc%slow_soil_C(np) + sbgc%slow_soil_C(cp) * area_fac
    sbgc%struct_soil_C(np) = sbgc%struct_soil_C(np) + sbgc%struct_soil_C(cp) * area_fac
    sbgc%struct_soil_L(np) = sbgc%struct_soil_L(np) + sbgc%struct_soil_L(cp) * area_fac
    sbgc%miner_soil_N(np) = sbgc%miner_soil_N(np) + sbgc%miner_soil_N(cp) * area_fac
    sbgc%fast_soil_N(np) = sbgc%fast_soil_N(np) + sbgc%fast_soil_N(cp) * area_fac
    sbgc%struct_soil_N(np) = sbgc%struct_soil_N(np) + sbgc%struct_soil_N(cp) * area_fac
    sbgc%slow_soil_N(np) = sbgc%slow_soil_N(np) + sbgc%slow_soil_N(cp) * area_fac 
    sbgc%fast_soil_P(np) = sbgc%fast_soil_P(np) + sbgc%fast_soil_P(cp) * area_fac
    sbgc%struct_soil_P(np) = sbgc%struct_soil_P(np) + sbgc%struct_soil_P(cp) * area_fac
    sbgc%miner_soil_P(np) = sbgc%miner_soil_P(np) + sbgc%miner_soil_P(cp) * area_fac
    sbgc%slow_soil_P(np) = sbgc%slow_soil_P(np) + sbgc%slow_soil_P(cp) * area_fac

    return
  end subroutine inc_patch_vars_sbgc


  !===========================================================

  !-------------------------------------------------
  ! 9. Required if you want forest harvesting.
  !-------------------------------------------------
  subroutine norm_harv_patch_sbgc(sbgc, np, area_fac)
    implicit none

    type(sbgc_vars) :: sbgc
    integer, intent(in) :: np
    real, intent(in) :: area_fac

    sbgc%fast_soil_C(np) = sbgc%fast_soil_C(np) * area_fac
    sbgc%slow_soil_C(np) = sbgc%slow_soil_C(np) * area_fac
    sbgc%struct_soil_C(np) = sbgc%struct_soil_C(np) * area_fac
    sbgc%struct_soil_L(np) = sbgc%struct_soil_L(np) * area_fac
    sbgc%miner_soil_N(np) = sbgc%miner_soil_N(np) * area_fac
    sbgc%fast_soil_N(np) = sbgc%fast_soil_N(np) * area_fac
    sbgc%struct_soil_N(np) = sbgc%struct_soil_N(np) * area_fac
    sbgc%slow_soil_N(np) = sbgc%slow_soil_N(np) * area_fac
    sbgc%fast_soil_P(np) = sbgc%fast_soil_P(np) * area_fac
    sbgc%struct_soil_P(np) = sbgc%struct_soil_P(np) * area_fac
    sbgc%miner_soil_P(np) = sbgc%miner_soil_P(np) * area_fac
    sbgc%slow_soil_P(np) = sbgc%slow_soil_P(np) * area_fac

    return
  end subroutine norm_harv_patch_sbgc

  !===========================================================

  !-------------------------------------------------
  ! 10. Required for patch fusion.
  !-------------------------------------------------
  subroutine fuse_2_patches_sbgc(sbgc, recp, donp, newareai, aread, arear)
    implicit none
    
    type(sbgc_vars) :: sbgc
    integer, intent(in) :: recp, donp
    real, intent(in) :: newareai, aread, arear

    sbgc%fast_soil_C(recp) = newareai * ( sbgc%fast_soil_C(donp) * aread   &
         + sbgc%fast_soil_C(recp) * arear )
    sbgc%slow_soil_C(recp)        = newareai * ( sbgc%slow_soil_C(donp) * aread   &
         + sbgc%slow_soil_C(recp)        * arear )
    sbgc%struct_soil_C(recp)  = newareai *                                          &
         ( sbgc%struct_soil_C(donp)  * aread   &
         + sbgc%struct_soil_C(recp)  * arear )
    sbgc%struct_soil_L(recp)  = newareai *                                          &
         ( sbgc%struct_soil_L(donp)  * aread   &
         + sbgc%struct_soil_L(recp)  * arear )
    
    sbgc%fast_soil_P(recp) = newareai *                                          &
         ( sbgc%fast_soil_P(donp) * aread   &
         + sbgc%fast_soil_P(recp) * arear )
    sbgc%struct_soil_P(recp) = newareai *                                          &
         ( sbgc%struct_soil_P(donp) * aread   &
         + sbgc%struct_soil_P(recp) * arear )
    sbgc%miner_soil_P(recp) = newareai *                                          &
         ( sbgc%miner_soil_P(donp) * aread   &
         + sbgc%miner_soil_P(recp) * arear )
    sbgc%slow_soil_P(recp) = newareai *                                          &
         ( sbgc%slow_soil_P(donp) * aread   &
         + sbgc%slow_soil_P(recp) * arear )
    
    sbgc%miner_soil_N(recp) = newareai *                                          &
         ( sbgc%miner_soil_N(donp) * aread   &
         + sbgc%miner_soil_N(recp) * arear )
    sbgc%fast_soil_N(recp) = newareai *                                          &
         ( sbgc%fast_soil_N(donp)        * aread   &
         + sbgc%fast_soil_N(recp)        * arear )
    sbgc%struct_soil_N(recp) = newareai *                                          &
         ( sbgc%struct_soil_N(donp)        * aread   &
         + sbgc%struct_soil_N(recp)        * arear )
    sbgc%slow_soil_N(recp) = newareai *                                          &
         ( sbgc%slow_soil_N(donp)        * aread   &
         + sbgc%slow_soil_N(recp)        * arear )
    

    return
  end subroutine fuse_2_patches_sbgc

  !======================================================

  !-------------------------------------------------
  ! 11. Required for some aspects of history file I/O
  !-------------------------------------------------
  subroutine fill_history_site_sbgc(sbgc, dsetrank, iparallel)
    implicit none
    
    type(sbgc_vars) :: sbgc
    integer :: dsetrank, iparallel

    call hdf_getslab_r(sbgc%fast_soil_C,'FAST_SOIL_C ',dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%slow_soil_C,'SLOW_SOIL_C ',dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%struct_soil_C,'STRUCT_SOIL_C ',dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%struct_soil_L,'STRUCT_SOIL_L ',dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%miner_soil_N,'MINER_SOIL_N ',dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%fast_soil_P,'FAST_SOIL_P ',dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%struct_soil_P,'STRUCT_SOIL_P ',dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%miner_soil_P,'MINER_SOIL_P ',dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%slow_soil_P,'SLOW_SOIL_P ',dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%fast_soil_N,'FAST_SOIL_N ',dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%slow_soil_N,'SLOW_SOIL_N ',dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%struct_soil_N,'STRUCT_SOIL_N ',dsetrank,iparallel,.true.)
    
    
    return
  end subroutine fill_history_site_sbgc

  !======================================================

  !-------------------------------------------------
  ! 12. Required for some aspects of history file I/O
  !-------------------------------------------------
  subroutine read_ed21_history_sbgc(sbgc, dsetrank, iparallel)
    implicit none
    
    type(sbgc_vars) :: sbgc
    integer :: dsetrank, iparallel

    call hdf_getslab_r(sbgc%fast_soil_C       ,'FAST_SOIL_C '                  &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%slow_soil_C       ,'SLOW_SOIL_C '                  &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%fast_soil_N       ,'FAST_SOIL_N '                  &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%struct_soil_C ,'STRUCT_SOIL_C '            &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%struct_soil_L ,'STRUCT_SOIL_L '            &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%miner_soil_N,'MINER_SOIL_N '           &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%slow_soil_N,'SLOW_SOIL_N '           &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%struct_soil_N,'STRUCT_SOIL_N '           &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%fast_soil_P,'FAST_SOIL_P '           &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%struct_soil_P,'STRUCT_SOIL_P '           &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%miner_soil_P,'MINER_SOIL_P '           &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%slow_soil_P,'SLOW_SOIL_P '           &
         ,dsetrank,iparallel,.true.)
    
    
    return
  end subroutine read_ed21_history_sbgc

  !===============================================================
  
  !-------------------------------------------------
  ! 13. Required for some aspects of history file I/O
  !-------------------------------------------------
  subroutine read_ed21_history_unstr_sbgc(sbgc,dsetrank,iparallel)
    implicit none
    type(sbgc_vars) :: sbgc
    integer :: dsetrank,iparallel

    call hdf_getslab_r(sbgc%fast_soil_C       ,'FAST_SOIL_C '               &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%slow_soil_C       ,'SLOW_SOIL_C '               &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%fast_soil_N       ,'FAST_SOIL_N '               &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%struct_soil_C ,'STRUCT_SOIL_C '         &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%struct_soil_N ,'STRUCT_SOIL_N '         &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%slow_soil_N ,'SLOW_SOIL_N '         &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%struct_soil_L ,'STRUCT_SOIL_L '         &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%miner_soil_N,'MINER_SOIL_N '        &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%fast_soil_P,'FAST_SOIL_P '        &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%struct_soil_P,'STRUCT_SOIL_P '        &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%miner_soil_P,'MINER_SOIL_P '        &
         ,dsetrank,iparallel,.true.)
    call hdf_getslab_r(sbgc%slow_soil_P,'SLOW_SOIL_P '        &
         ,dsetrank,iparallel,.true.)
    
    
    return
  end subroutine read_ed21_history_unstr_sbgc


end Module soil_bgc

