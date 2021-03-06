module fuse_fiss_utils

   use ed_state_vars,only :   copy_patchtype        & ! subroutine
                            , deallocate_patchtype  & ! subroutine
                            , allocate_patchtype    & ! subroutine
                            , allocate_sitetype     & ! subroutine
                            , deallocate_sitetype   & ! subroutine
                            , copy_sitetype_mask    & ! subroutine
                            , copy_sitetype         & ! subroutine
                            , copy_patchtype_mask   ! ! subroutine

   contains
   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine will sort the cohorts by size (1st = tallest, last = shortest.)     !
   !---------------------------------------------------------------------------------------!
   subroutine sort_cohorts(cpatch)

      use ed_state_vars,only :  patchtype   ! ! Structure
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype), target  :: cpatch     ! Current patch, to have cohorts sorted.
      !----- Local variables --------------------------------------------------------------!
      type(patchtype), pointer :: temppatch  ! Scratch patch structure
      integer                  :: ico        ! Counters
      integer                  :: tallco     ! Index of tallest cohort
      logical                  :: sorted     ! Flag: this patch is already sorted
      !------------------------------------------------------------------------------------!
      
      !----- No need to sort an empty patch or a patch with a single cohort. --------------!
      if (cpatch%ncohorts < 2) return


      !------------------------------------------------------------------------------------!
      !     Check whether this patch is already sorted.   We don't want to do the entire   !
      ! deallocating/copying/allocating thing if it's not needed as this takes up too much !
      ! time.                                                                              !
      !------------------------------------------------------------------------------------!
      sorted = .true.
      sortcheck: do ico=1,cpatch%ncohorts-1
         sorted = cpatch%costate%hite(ico) >= cpatch%costate%hite(ico+1)
         if (.not. sorted) exit sortcheck
      end do sortcheck
      if (sorted) return
      !------------------------------------------------------------------------------------!



      !----- Assign a scratch patch. ------------------------------------------------------!
      nullify(temppatch)
      allocate(temppatch)
      call allocate_patchtype(temppatch,cpatch%ncohorts)
      
      ico = 0
      !---- Loop until all cohorts were sorted. -------------------------------------------!
      do while(ico < cpatch%ncohorts)
         ico = ico + 1
      
         !----- Find the tallest cohort. --------------------------------------------------!
         tallco = maxloc(cpatch%costate%hite,dim=1)
         
         !----- Copy to the scratch structure. --------------------------------------------!
         call copy_patchtype(cpatch,temppatch,tallco,tallco,ico,ico)
         
         !----- Put a non-sense height so this will never "win" again. --------------------!
         cpatch%costate%hite(tallco) = -huge(1.)

      end do

      !------ Copy the scratch patch to the regular one and deallocate it. ----------------!
      call copy_patchtype(temppatch,cpatch,1,cpatch%ncohorts,1,cpatch%ncohorts)
      call deallocate_patchtype(temppatch)
      deallocate(temppatch)

      return

   end subroutine sort_cohorts
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will eliminate cohorts based on their sizes. This is intended to   !
   ! eliminate cohorts that have little contribution and thus we can speed up the run.     !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_cohorts(csite,ipa,elim_nplant,elim_lai)
      use pft_coms           , only : min_cohort_size  & ! intent(in)
                                    , l2n_stem         & ! intent(in)
                                    , c2n_stem         & ! intent(in)
                                    , c2n_storage      & ! intent(in), lookup table
                                    , c2n_leaf         ! ! intent(in), lookup table
      use pft_coms, only: c2p_alive, c2p_storage, c2p_dead
      use decomp_coms        , only : f_labile         ! ! intent(in), lookup table

      use ed_state_vars      , only : patchtype        & ! structure
                                    , sitetype         ! ! structure
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)       , target      :: csite        ! Current site
      integer              , intent(in)  :: ipa          ! Current patch ID
      real                 , intent(out) :: elim_nplant  ! Nplants eliminated here
      real                 , intent(out) :: elim_lai     ! LAI eliminated here
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)      , pointer     :: cpatch       ! Current patch
      type(patchtype)      , pointer     :: temppatch    ! Scratch patch structure
      logical, dimension(:), allocatable :: remain_table ! Flag: this cohort will remain.
      integer                            :: ico, inew    ! Counters
      integer                            :: ipft         ! PFT size
      real                               :: csize        ! Size of current cohort
      !------------------------------------------------------------------------------------!
      
      cpatch => csite%patch(ipa)
      elim_nplant = 0.
      elim_lai    = 0.

      !----- Initialize the temporary patch structures and the remain/terminate table -----!
      nullify(temppatch)
      allocate(temppatch)
      allocate(remain_table(cpatch%ncohorts))
      remain_table(:) = .true.
     
      !----- Main loop --------------------------------------------------------------------!
      do ico = 1,cpatch%ncohorts

         !----- Save the PFT type in a convenient alias. ----------------------------------!
         ipft = cpatch%costate%pft(ico)

         !----- Checking whether the cohort size is too small -----------------------------!
         csize = cpatch%costate%nplant(ico)                                                        &
               * (cpatch%costate%balive(ico) + cpatch%costate%bdead(ico) + cpatch%costate%bstorage(ico))

         if ( csize < min_cohort_size(ipft) ) then
            !----- Cohort is indeed too small, it won't remain ----------------------------!
            remain_table(ico) = .false.
            elim_nplant = elim_nplant + cpatch%costate%nplant(ico) * csite%area(ipa)
            elim_lai    = elim_lai    + cpatch%costate%lai(ico)    * csite%area(ipa)

            !----- Update litter pools ----------------------------------------------------!
            csite%fsc_in(ipa) = csite%fsc_in(ipa) + cpatch%costate%nplant(ico)                     &
                              * (f_labile(ipft)*cpatch%costate%balive(ico) + cpatch%costate%bstorage(ico))

            csite%fsn_in(ipa) = csite%fsn_in(ipa) + cpatch%costate%nplant(ico)                     &
                              * ( f_labile(ipft) * cpatch%costate%balive(ico) / c2n_leaf(ipft)     &
                                + cpatch%costate%bstorage(ico) / c2n_storage)
            
            csite%ssc_in(ipa) = csite%ssc_in(ipa) + cpatch%costate%nplant(ico)                     &
                              * ( (1.0 - f_labile(ipft)) * cpatch%costate%balive(ico)              &
                                + cpatch%costate%bdead(ico))
            
            csite%ssl_in(ipa) = csite%ssl_in(ipa) + cpatch%costate%nplant(ico)                     &
                              * ( (1.0 - f_labile(ipft)) * cpatch%costate%balive(ico)              &
                                + cpatch%costate%bdead(ico) ) * l2n_stem/c2n_stem(ipft)




            csite%fsp_in(ipa) = csite%fsp_in(ipa) + cpatch%costate%nplant(ico)                     &
                              * ( f_labile(ipft) * cpatch%costate%balive(ico) / c2p_alive(ipft)     &
                                + cpatch%costate%bstorage(ico) / c2p_storage(ipft))
            
            csite%stsp_in(ipa) = csite%stsp_in(ipa) + cpatch%costate%nplant(ico)                     &
                              * ( (1.0 - f_labile(ipft)) * cpatch%costate%balive(ico)              &
                                + cpatch%costate%bdead(ico)) / c2p_dead(ipft)
            
         end if
      end do

      !----- Copy the remaining cohorts to a temporary patch ------------------------------!
      call allocate_patchtype(temppatch,count(remain_table))
      call copy_patchtype_mask(cpatch,temppatch,remain_table,size(remain_table)            &
                              ,count(remain_table))

      !----- Reallocate the new patch and populate with the saved cohorts -----------------!
      call deallocate_patchtype(cpatch)
      call allocate_patchtype(cpatch,count(remain_table))
      call copy_patchtype(temppatch,cpatch,1,cpatch%ncohorts,1,cpatch%ncohorts)
      call sort_cohorts(cpatch)
     
      !----- Deallocate the temporary patch -----------------------------------------------!     
      call deallocate_patchtype(temppatch)
      deallocate(temppatch)
      deallocate(remain_table)

      !----- Update the cohort census at the site level -----------------------------------!
      csite%cohort_count(ipa) = cpatch%ncohorts
      if (elim_lai < 0. .or. elim_nplant < 0.) then
         write (unit=*,fmt='(a,1x,es12.5)') 'TERMINATE: ELIM_LAI=',elim_lai
         write (unit=*,fmt='(a,1x,es12.5)') 'TERMINATE: ELIM_NPLANT=',elim_nplant
      end if
      return
   end subroutine terminate_cohorts
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will eliminate tiny or empty patches. This is intended to          !
   ! eliminate patches that have little contribution and thus we can speed up the run.     !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_patches(csite)

      use ed_state_vars, only : polygontype        & ! Structure
                              , sitetype           & ! Structure
                              , patchtype          ! ! Structure
      use disturb_coms , only : min_new_patch_area ! ! intent(in)
      use ed_misc_coms , only : iqoutput           & ! intent(in)
                              , imoutput           & ! intent(in)
                              , idoutput           ! ! intent(in)
      use cohort_state, only: terminate_patches_state
      use cohort_mort, only: terminate_patches_mort
      use cohort_therm, only: terminate_patches_therm
      use cohort_resp, only: terminate_patches_resp
      use cohort_photo, only: terminate_patches_photo
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)       , target      :: csite        ! Current site
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)       , pointer     :: tempsite     ! Scratch site
      type(patchtype)      , pointer     :: cpatch       ! Pointer to current site
      integer                            :: ipa,ico      ! Counters
      logical, dimension(:), allocatable :: remain_table ! Flag: this patch will remain.
      real                               :: elim_area    ! Area of removed patches
      real                               :: new_area     ! Just to make sure area is 1.
      real                               :: area_scale   ! Scaling area factor.
      !------------------------------------------------------------------------------------!

      allocate (remain_table(csite%npatches))
      remain_table(:) = .true.

      !------------------------------------------------------------------------------------!
      !     Loop through all the patches in this site and determine which of these patches !
      ! is too small in area to be valid. Remove these patches via the mask function.      !
      ! Realocate a new site with only the valid patches, and normalize their areas and    !
      ! plant densities to reflect the area loss.                                          !
      !------------------------------------------------------------------------------------!
      elim_area = 0.0
      do ipa = 1,csite%npatches
         if (csite%area(ipa) < min_new_patch_area) then
            elim_area = elim_area + csite%area(ipa)
            remain_table(ipa) = .false.
         end if
      end do

      !----- Use the mask to resize the patch vectors in the current site. ----------------!
      allocate(tempsite)
      call allocate_sitetype(tempsite,count(remain_table))
      call copy_sitetype_mask(csite,tempsite,remain_table,size(remain_table)               &
                             ,count(remain_table))
      call deallocate_sitetype(csite)
      call allocate_sitetype(csite,count(remain_table))

      remain_table(:)                   = .false.
      remain_table(1:tempsite%npatches) = .true.
      call copy_sitetype_mask(tempsite,csite,remain_table(1:tempsite%npatches)             &
                             ,count(remain_table),count(remain_table))
      call deallocate_sitetype(tempsite)
      deallocate(tempsite)

      !------------------------------------------------------------------------------------!
      !    Renormalize the total area.  We must also rescale all extensive properties from !
      ! cohorts, since they are per unit area and we are effectively changing the area.    !
      ! IMPORTANT: Only cohort-level variables that have units per area (m2) should be     !
      !            rescaled.  Variables whose units are per plant should _NOT_ be included !
      !            here.                                                                   !
      !------------------------------------------------------------------------------------!
      new_area=0.
      area_scale = 1./(1. - elim_area)
      do ipa = 1,csite%npatches
         csite%area(ipa) = csite%area(ipa) * area_scale
         new_area = new_area + csite%area(ipa)

         cpatch => csite%patch(ipa)

         call terminate_patches_state(cpatch%costate, cpatch%ncohorts, area_scale)
         call terminate_patches_photo(cpatch%cophoto, cpatch%ncohorts, area_scale)
         call terminate_patches_resp(cpatch%coresp, cpatch%ncohorts, area_scale)
         call terminate_patches_therm(cpatch%cotherm, cpatch%ncohorts, area_scale)
         call terminate_patches_mort(cpatch%comort, cpatch%ncohorts, area_scale)
      enddo

      if (abs(new_area-1.0) > 1.e-5) then
         write (unit=*,fmt='(a,1x,es12.5)') ' + ELIM_AREA:',elim_area
         write (unit=*,fmt='(a,1x,es12.5)') ' + NEW_AREA: ',new_area
         call fatal_error('New_area should be 1 but it isn''t!!!','terminate_patches'      &
                         ,'fuse_fiss_utils.f90')
      end if 
      
      return
   end subroutine terminate_patches
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine will perform cohort fusion based on various similarity criteria to  !
   ! determine whether they can be fused with no significant loss of information. The user !
   ! is welcome to set up a benchmark, but should be aware that no miracles will happen    !
   ! here. If there are more very distinct cohorts than maxcohort, then the user will need !
   ! to live with that and accept life is not always fair with those with limited          !
   ! computational resources.                                                              !
   !---------------------------------------------------------------------------------------!

   subroutine fuse_cohorts(csite,ipa, green_leaf_factor, lsl)

      use ed_state_vars       , only : sitetype            & ! Structure
                                     , patchtype           ! ! Structure
      use pft_coms            , only : rho                 & ! intent(in)
                                     , b1Ht                & ! intent(in)
                                     , hgt_max             & ! intent(in)
                                     , sla                 & ! intent(in)
                                     , hgt_ref             ! ! intent(in)
      use fusion_fission_coms , only : fusetol_h           & ! intent(in)
                                     , fusetol             & ! intent(in)
                                     , lai_fuse_tol        & ! intent(in)
                                     , fuse_relax          & ! intent(in)
                                     , coh_tolerance_max   ! ! intent(in)
      use ed_max_dims         , only : n_pft               ! ! intent(in)
      use mem_polygons        , only : maxcohort           ! ! intent(in)
      use canopy_layer_coms   , only : crown_mod           ! ! intent(in)
      use allometry           , only : dbh2h               & ! function
                                     , dbh2bl              ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target      :: csite             ! Current site
      integer                , intent(in)  :: ipa               ! Current patch ID
      real, dimension(n_pft) , intent(in)  :: green_leaf_factor ! 
      integer                , intent(in)  :: lsl               ! Lowest soil level
      !----- Local variables --------------------------------------------------------------!
      logical, dimension(:)  , allocatable :: fuse_table     ! Flag, remaining cohorts
      type(patchtype)        , pointer     :: cpatch         ! Current patch
      type(patchtype)        , pointer     :: temppatch      ! Scratch patch
      integer                              :: donc,recc,ico3 ! Counters
      logical                              :: fusion_test    ! Flag: proceed with fusion?
      real                                 :: newn           ! new nplants of merged coh.
      real                                 :: lai_max        ! Maximum LAI the fused 
                                                             !    cohort could have.
      real                                 :: total_size     ! Total size
      real                                 :: tolerance_mult ! Multiplication factor
      integer                              :: ncohorts_old   ! # of coh. before fusion test
      real                                 :: mean_dbh       ! Mean DBH           (???)
      real                                 :: mean_hite      ! Mean height        (???)
      real                                 :: new_size       ! New size
      integer                              :: ntall          ! # of tall cohorts  (???)
      integer                              :: nshort         ! # of short cohorts (???)
      logical                              :: any_fusion     ! Flag: was there any fusion?
      !------------------------------------------------------------------------------------!


      !----- Start with no factor ---------------------------------------------------------!
      tolerance_mult = 1.0

      cpatch => csite%patch(ipa)

      !------------------------------------------------------------------------------------!
      !     Return if maxcohort is 0 (flag for no cohort fusion), or if the patch is empty !
      ! or has a single cohort.                                                            !
      !------------------------------------------------------------------------------------!
      if (maxcohort == 0 .or. cpatch%ncohorts < 2) return

      !------------------------------------------------------------------------------------!
      !    Calculate mean DBH and HITE to help with the normalization of differences mean  !
      ! hite is not being used right now, but can be optioned in the future if it seems    !
      ! advantageous.                                                                      !
      !------------------------------------------------------------------------------------!
      mean_dbh  = 0.0
      mean_hite = 0.0
      nshort    = 0
      ntall     = 0
      do ico3 = 1,cpatch%ncohorts
         !---------------------------------------------------------------------------------!
         !    Get fusion height threshold.  Height is a good predictor when plants are     !
         ! growing in height, but it approaches the maximum height DBH becomes the only    !
         ! possible predictor because height saturates.                                    !
         !---------------------------------------------------------------------------------!
         if (cpatch%costate%hite(ico3) < (0.95 * hgt_max(cpatch%costate%pft(ico3))) ) then
            mean_hite = mean_hite + cpatch%costate%hite(ico3)
            nshort    = nshort + 1
         else
            mean_dbh  = mean_dbh + cpatch%costate%dbh(ico3)
            ntall     = ntall + 1
         end if
      end do 
      !------------------------------------------------------------------------------------!
      if (ntall  > 0) mean_dbh = mean_dbh   / real(ntall)
      if (nshort > 0) mean_hite= mean_hite  / real(nshort)

      !----- Initialize table. In principle, all cohorts stay. ----------------------------!
      allocate(fuse_table(cpatch%ncohorts))
      fuse_table(:) = .true.

      force_fusion: do
         
         ncohorts_old =  count(fuse_table) ! Save current number of cohorts ---------------!
         
         donloop:do donc = 1,cpatch%ncohorts-1
            if (.not. fuse_table(donc)) cycle donloop ! This one is gone, move to next.

            recloop: do recc = donc+1,cpatch%ncohorts
               if (.not. fuse_table(recc)) cycle recloop ! This one is gone, move to next.
                                                         ! Hope it never happens...

               !---------------------------------------------------------------------------!
               !     Test for similarity.  Again, we use height to assess similarity only  !
               ! when the cohort is not approaching the maximum height.  If this is the    !
               ! case, then we use DBH to test.                                            !
               !---------------------------------------------------------------------------!
               if (cpatch%costate%hite(donc) >= (0.95 * hgt_max(cpatch%costate%pft(donc))) ) then
                  mean_dbh=0.5*(cpatch%costate%dbh(donc)+cpatch%costate%dbh(recc))
                  fusion_test = ( abs(cpatch%costate%dbh(donc) - cpatch%costate%dbh(recc)))/mean_dbh       &
                              < fusetol * tolerance_mult
               elseif (fuse_relax) then
                  fusion_test = ( abs(cpatch%costate%hite(donc) - cpatch%costate%hite(recc))               &
                                     / (0.5*(cpatch%costate%hite(donc) + cpatch%costate%hite(recc)))  <    &
                                fusetol * tolerance_mult)  
               else
                  fusion_test = (abs(cpatch%costate%hite(donc) - cpatch%costate%hite(recc))  <             &
                                fusetol_h * tolerance_mult)
               end if

               if (fusion_test) then

                  !----- New cohort has the total number of plants ------------------------!
                  newn = cpatch%costate%nplant(donc) + cpatch%costate%nplant(recc)

                  !------------------------------------------------------------------------!
                  !     We now check the maximum LAI the fused cohorts could have.  We     !
                  ! don't want the cohort to have a very large LAI.  If both cohorts have  !
                  ! leaves fully flushed, this is the same as adding the individual LAIs,  !
                  ! but if they are not, we need to consider that LAI may grow...          !
                  !------------------------------------------------------------------------!
                  lai_max = ( cpatch%costate%nplant(recc)                                          &
                            * dbh2bl(cpatch%costate%dbh(recc),cpatch%costate%pft(recc))                    &
                            + cpatch%costate%nplant(donc)                                          &
                            * dbh2bl(cpatch%costate%dbh(donc),cpatch%costate%pft(donc)))                   &
                          * cpatch%cophen%sla(recc)

                  !----- Checking the total size of this cohort before and after fusion. --!
                  total_size = cpatch%costate%nplant(donc) * ( cpatch%costate%balive(donc)                 &
                                                     + cpatch%costate%bdead(donc)                  &
                                                     + cpatch%costate%bstorage(donc) )             &
                             + cpatch%costate%nplant(recc) * ( cpatch%costate%balive(recc)                 &
                                                     + cpatch%costate%bdead(recc)                  &
                                                     + cpatch%costate%bstorage(recc) )


                  !------------------------------------------------------------------------!
                  !    Five conditions must be met to allow two cohorts to be fused:       !
                  ! 1. Both cohorts must have the same PFT;                                !
                  ! 2. Combined LAI won't be too large.                                    !
                  ! 3. Both cohorts must have the same status with respect to the first    !
                  !    census.                                                             !
                  ! 4. Both cohorts must have the same recruit status with respect to the  !
                  !    first census.                                                       !
                  ! 5. Both cohorts must have the same phenology status.                   !
                  !------------------------------------------------------------------------!
                  if (     cpatch%costate%pft(donc)              == cpatch%costate%pft(recc)               &
                     .and. lai_max                        < lai_fuse_tol*tolerance_mult    &
                     .and. cpatch%costate%first_census(donc)     == cpatch%costate%first_census(recc)      &
                     .and. cpatch%costate%new_recruit_flag(donc) == cpatch%costate%new_recruit_flag(recc)  &
                     .and. cpatch%cophen%phenology_status(donc) == cpatch%cophen%phenology_status(recc)  &
                     ) then

                     !----- Proceed with fusion -------------------------------------------!
                     call fuse_2_cohorts(cpatch,donc,recc,newn                             &
                                        ,green_leaf_factor(cpatch%costate%pft(donc))               &
                                        ,csite%can_prss(ipa),lsl)

                     !----- Flag donating cohort as gone, so it won't be checked again. ---!
                     fuse_table(donc) = .false.
                     
                     !----- Checking whether total size and LAI are conserved. ------------!
                     new_size = cpatch%costate%nplant(recc) * ( cpatch%costate%balive(recc)                &
                                                      + cpatch%costate%bdead(recc)                 &
                                                      + cpatch%costate%bstorage(recc) )
                     if (new_size < 0.99* total_size .or. new_size > 1.01* total_size )    &
                     then
                        write (unit=*,fmt='(a,1x,es14.7)') 'OLD SIZE: ',total_size
                        write (unit=*,fmt='(a,1x,es14.7)') 'NEW SIZE: ',new_size
                        call fatal_error('Cohort fusion didn''t conserve plant size!!!'    &
                                        &,'fuse_2_cohorts','fuse_fiss_utils.f90')
                     end if
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !    Recalculate the means                                            !
                     !---------------------------------------------------------------------!
                     mean_dbh  = 0.0
                     mean_hite = 0.0
                     nshort    = 0
                     ntall     = 0
                     recalcloop: do ico3 = 1,cpatch%ncohorts
                        if (.not. fuse_table(ico3)) cycle recalcloop
                        !----- Get fusion height threshold --------------------------------!
                        if (cpatch%costate%hite(ico3) < (0.95 * hgt_max(cpatch%costate%pft(ico3))) ) then
                           mean_hite = mean_hite + cpatch%costate%hite(ico3)
                           nshort = nshort+1
                        else
                           mean_dbh = mean_dbh + cpatch%costate%dbh(ico3)
                           ntall=ntall+1
                        end if
                     end do recalcloop
                     !---------------------------------------------------------------------!
                     cycle donloop
                  end if
                  !------------------------------------------------------------------------!
               end if
            end do recloop
         end do donloop

         !------ If we are under maxcohort, no need to continue fusing. -------------------!
         if ( count(fuse_table) <= abs(maxcohort)) exit force_fusion
         !------ If no fusion happened and the tolerance exceeded the maximum, I give up. -!
         if ( count(fuse_table) == ncohorts_old .and. tolerance_mult > coh_tolerance_max ) &
            exit force_fusion

         tolerance_mult = tolerance_mult * 1.01
         ncohorts_old = count(fuse_table)
      end do force_fusion

      !----- If any fusion has happened, then we need to rearrange cohorts. ---------------!
      any_fusion = .not. all(fuse_table)
      if (any_fusion) then

         !---------------------------------------------------------------------------------!
         !     Now copy the merged patch to a temporary patch using the fuse_table as a    !
         ! mask.  Then allocate a temporary patch, copy the remaining cohorts there.       !
         !---------------------------------------------------------------------------------!
         nullify (temppatch)
         allocate(temppatch)
         call allocate_patchtype(temppatch,cpatch%ncohorts)
         call copy_patchtype_mask(cpatch,temppatch,fuse_table,size(fuse_table)             &
                                 ,count(fuse_table))

         !----- Now I reallocate the current patch with its new reduced size. -------------!
         call deallocate_patchtype(cpatch)  
         call allocate_patchtype(cpatch,count(fuse_table))
  
         !----- Make fuse_table true to all remaining cohorts. ----------------------------!
         fuse_table(:)                 = .false.
         fuse_table(1:cpatch%ncohorts) = .true.
         call copy_patchtype_mask(temppatch,cpatch,fuse_table,size(fuse_table)             &
                                 ,count(fuse_table))

         !----- Discard the scratch patch. ------------------------------------------------!
         call deallocate_patchtype(temppatch)
         deallocate(temppatch)  

         !----- Sort cohorts by size again, and update the cohort census for this patch. --!
         call sort_cohorts(cpatch)
         csite%cohort_count(ipa) = count(fuse_table)
      end if

      !----- Deallocate the aux. table ----------------------------------------------------!
      deallocate(fuse_table)
     
      return
   end subroutine fuse_cohorts
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine will split two cohorts if its LAI has become too large.  This is    !
   ! only necessary when we solve radiation cohort by cohort rather than layer by layer.   !
   !---------------------------------------------------------------------------------------!
   subroutine split_cohorts(cpatch, green_leaf_factor, lsl)

      use ed_state_vars        , only : patchtype              ! ! structure
      use pft_coms             , only : q                      & ! intent(in), lookup table
                                      , qsw                    ! ! intent(in), lookup table
      use fusion_fission_coms  , only : lai_tol                ! ! intent(in)
      use ed_max_dims          , only : n_pft                  ! ! intent(in)
      use allometry            , only : dbh2h                  & ! function
                                      , bd2dbh                 & ! function
                                      , dbh2bd                 ! ! function
      use ed_misc_coms         , only : iqoutput               & ! intent(in)
                                      , imoutput               & ! intent(in)
                                      , idoutput               ! ! intent(in)
      use canopy_layer_coms    , only : crown_mod              ! ! intent(in)
      use cohort_state, only: split_cohorts_state
      use cohort_photo, only: split_cohorts_photo
      use cohort_resp, only: split_cohorts_resp
      use cohort_therm, only: split_cohorts_therm
      use cohort_mort, only: split_cohorts_mort
      implicit none
      !----- Constants --------------------------------------------------------------------!
      real                   , parameter   :: epsilon=0.0001    ! Tweak factor...
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype)        , target      :: cpatch            ! Current patch
      real, dimension(n_pft) , intent(in)  :: green_leaf_factor !
      integer                , intent(in)  :: lsl               ! Lowest soil level
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)        , pointer     :: temppatch         ! Temporary patch
      logical, dimension(:)  , allocatable :: split_mask        ! Flag: split this cohort
      integer                              :: ipa,ico,inew      ! Counters
      integer                              :: ncohorts_new      ! New # of cohorts
      integer                              :: tobesplit         ! # of cohorts to be split
      integer                              :: ipft              ! PFT type
      real                                 :: stai              ! Potential TAI
      real                                 :: old_nplant        ! Old nplant
      real                                 :: new_nplant        ! New nplant
      real                                 :: old_size          ! Old size
      real                                 :: new_size          ! New size
      !------------------------------------------------------------------------------------!


      !----- Initialize the vector with splitting table -----------------------------------!
      allocate(split_mask(cpatch%ncohorts))
      split_mask(:) = .false.
      old_nplant = 0.
      old_size   = 0.
      !----- Loop through cohorts ---------------------------------------------------------!
      do ico = 1,cpatch%ncohorts
         ipft = cpatch%costate%pft(ico)

         !---------------------------------------------------------------------------------! 
         !     STAI is the potential TAI that this cohort has when its leaves are fully    !
         ! flushed.                                                                        !
         !---------------------------------------------------------------------------------! 
         stai = cpatch%costate%nplant(ico) * cpatch%costate%balive(ico) * green_leaf_factor(ipft)          &
              * q(ipft) / ( 1.0 + q(ipft) + qsw(ipft) * cpatch%costate%hite(ico) )                 &
              * cpatch%cophen%sla(ico) + cpatch%costate%wai(ico)

         !----- If the resulting TAI is too large, split this cohort. ---------------------!
         split_mask(ico) = stai > lai_tol
         
         old_nplant = old_nplant + cpatch%costate%nplant(ico)
         old_size   = old_size   + cpatch%costate%nplant(ico) * ( cpatch%costate%balive(ico)               &
                                                        + cpatch%costate%bdead(ico)                &
                                                        + cpatch%costate%bstorage(ico) )
      end do

      !----- Compute the new number of cohorts. -------------------------------------------!
      tobesplit    = count(split_mask)
      ncohorts_new = cpatch%ncohorts + tobesplit
      
      if (tobesplit > 0) then

         !----- Allocate the temppatch. ---------------------------------------------------!
         nullify(temppatch)
         allocate(temppatch)
         call allocate_patchtype(temppatch,cpatch%ncohorts)

         !----- Fill the temp space with the current patches. -----------------------------!
         call copy_patchtype(cpatch,temppatch,1,cpatch%ncohorts,1,cpatch%ncohorts)

         !----- Deallocate the current patch. ---------------------------------------------!
         call deallocate_patchtype(cpatch)

         !----- Re-allocate the current patch. --------------------------------------------!
         call allocate_patchtype(cpatch,ncohorts_new)

         !----- Transfer the temp values back in. -----------------------------------------!
         call copy_patchtype(temppatch,cpatch,1,temppatch%ncohorts,1,temppatch%ncohorts)

         !----- Remove the temporary patch. -----------------------------------------------!
         call deallocate_patchtype(temppatch)
         deallocate(temppatch)
     
         inew = size(split_mask)
         do ico = 1,size(split_mask)

            if (split_mask(ico)) then

               !---------------------------------------------------------------------------!
               !   Half the densities of the original cohort.  All "extensive" variables   !
               ! need to be rescaled.                                                      !

               !---------------------------------------------------------------------------!
               call split_cohorts_state(cpatch%costate, ico)
               call split_cohorts_photo(cpatch%cophoto, ico)
               call split_cohorts_resp(cpatch%coresp, ico)
               call split_cohorts_therm(cpatch%cotherm, ico)
               call split_cohorts_mort(cpatch%comort, ico)



               !----- Apply these values to the new cohort. -------------------------------!
               inew = inew+1
               call clone_cohort(cpatch,ico,inew)
               !---------------------------------------------------------------------------!

               !----- Tweaking bdead, to ensure carbon is conserved. ----------------------!
               cpatch%costate%bdead(ico)  = cpatch%costate%bdead(ico) * (1.-epsilon)
               cpatch%costate%dbh  (ico)  = bd2dbh(cpatch%costate%pft(ico), cpatch%costate%bdead(ico))
               cpatch%costate%hite (ico)  = dbh2h(cpatch%costate%pft(ico), cpatch%costate%dbh(ico))

               cpatch%costate%bdead(inew) = cpatch%costate%bdead(inew) * (1.+epsilon)
               cpatch%costate%dbh  (inew) = bd2dbh(cpatch%costate%pft(inew), cpatch%costate%bdead(inew))
               cpatch%costate%hite (inew) = dbh2h(cpatch%costate%pft(inew), cpatch%costate%dbh(inew))
               !---------------------------------------------------------------------------!

            end if
         end do

         !----- After splitting, cohorts may need to be sorted again... -------------------!
         call sort_cohorts(cpatch)

         !----- Checking whether the total # of plants is conserved... --------------------!
         new_nplant = 0.
         new_size   = 0.
         do ico=1,cpatch%ncohorts
            new_nplant = new_nplant + cpatch%costate%nplant(ico)
            new_size   = new_size   + cpatch%costate%nplant(ico) * ( cpatch%costate%balive(ico)            &
                                                           + cpatch%costate%bdead(ico)             &
                                                           + cpatch%costate%bstorage(ico) )
         end do
         if (new_nplant < 0.99 * old_nplant .or. new_nplant > 1.01 * old_nplant .or.       &
             new_size   < 0.99 * old_size   .or. new_size   > 1.01 * old_size) then
            write (unit=*,fmt='(a,1x,es14.7)') 'OLD NPLANT: ',old_nplant
            write (unit=*,fmt='(a,1x,es14.7)') 'NEW NPLANT: ',new_nplant
            write (unit=*,fmt='(a,1x,es14.7)') 'OLD SIZE:   ',old_size
            write (unit=*,fmt='(a,1x,es14.7)') 'NEW SIZE:   ',new_size
            call fatal_error('Cohort splitting didn''t conserve plants!!!'                 &
                                        &,'split_cohorts','fuse_fiss_utils.f90')
         end if
         
      end if
      deallocate(split_mask)
      return
   end subroutine split_cohorts
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine will clone one cohort.                                              !
   !---------------------------------------------------------------------------------------!
   subroutine clone_cohort(cpatch,isc,idt)
   
      use ed_max_dims  , only : n_mort     ! ! intent(in)
      use ed_state_vars, only : patchtype  & ! Structure
                              , stoma_data ! ! Structure
      use ed_misc_coms , only : iqoutput   & ! intent(in)
                              , idoutput   & ! intent(in)
                              , imoutput   ! ! intent(in)
      use cohort_state, only: clone_cohort_state
      use cohort_phen, only: clone_cohort_phen
      use cohort_mort, only: clone_cohort_mort
      use cohort_resp, only: clone_cohort_resp
      use cohort_photo, only: clone_cohort_photo
      use cohort_rad, only: clone_cohort_rad
      use cohort_therm, only: clone_cohort_therm
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype) , target     :: cpatch ! Current patch
      integer         , intent(in) :: isc    ! Index of "Source" cohort
      integer         , intent(in) :: idt    ! Index of "Destination" cohort"
      !----- Local variables --------------------------------------------------------------!
      integer                      :: imonth
      type(stoma_data), pointer    :: osdt,ossc
      !------------------------------------------------------------------------------------!

      call clone_cohort_state(cpatch%costate, isc, idt)
      call clone_cohort_phen(cpatch%cophen, isc, idt)
      call clone_cohort_mort(cpatch%comort, isc, idt)
      call clone_cohort_resp(cpatch%coresp, isc, idt)
      call clone_cohort_photo(cpatch%cophoto, isc, idt)
      call clone_cohort_rad(cpatch%corad, isc, idt)
      call clone_cohort_therm(cpatch%cotherm, isc, idt)

      return
   end subroutine clone_cohort
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will merge two cohorts into 1. The donating cohort (donc) is the  !
   ! one that will be deallocated, while the receptor cohort (recc) will contain the       !
   !  information from both cohorts.                                                       !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_cohorts(cpatch,donc,recc, newn,green_leaf_factor, can_prss,lsl)
      use ed_state_vars , only : patchtype              ! ! Structure
      use pft_coms      , only : q                      & ! intent(in), lookup table
                               , qsw                    ! ! intent(in), lookup table
      use therm_lib     , only : qwtk                   & ! subroutine
                               , rslif                  ! ! function
      use allometry     , only : dbh2krdepth            & ! function
                               , bd2dbh                 & ! function
                               , dbh2h                  ! ! function
      use ed_max_dims   , only : n_mort                 ! ! intent(in)
      use ed_misc_coms  , only : imoutput               & ! intent(in)
                               , iqoutput               & ! intent(in)
                               , idoutput               & ! intent(in)
                               , ndcycle                ! ! intent(in)
      use cohort_state, only: fuse_2_cohorts_state
      use cohort_phen, only: fuse_2_cohorts_phen
      use cohort_mort, only: fuse_2_cohorts_mort
      use cohort_resp, only: fuse_2_cohorts_resp
      use cohort_photo, only: fuse_2_cohorts_photo
      use cohort_therm, only: fuse_2_cohorts_therm
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype) , target     :: cpatch            ! Current patch
      integer                      :: donc              ! Donating cohort.
      integer                      :: recc              ! Receptor cohort.
      real            , intent(in) :: newn              ! New nplant
      real            , intent(in) :: green_leaf_factor ! Green leaf factor
      real            , intent(in) :: can_prss          ! Canopy air pressure
      integer         , intent(in) :: lsl               ! Lowest soil level
      !----- Local variables --------------------------------------------------------------!
      integer                      :: imon              ! Month for cb loop
      integer                      :: icyc              ! Time of day for dcycle loop
      integer                      :: imty              ! Mortality type
      real                         :: newni             ! Inverse of new nplants
      real                         :: newlaii           ! Inverse of new LAI
      real                         :: cb_act            !
      real                         :: cb_max            !
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Find the scaling factor for variables that are not "extensive".                 !
      !  - If the unit is X/plant, then we scale by nplant.                                !
      !  - If the unit is X/m2_leaf, then we scale by LAI.                                 !
      !  - If the unit is X/m2_gnd, then we add, since they are "extensive".               !
      !------------------------------------------------------------------------------------!
      newni   = 1.0 / newn
      if (cpatch%costate%lai(recc) + cpatch%costate%lai(donc) > 0.0) then
         newlaii = 1.0 / (cpatch%costate%lai(recc) + cpatch%costate%lai(donc))
      else
         newlaii = 0.0
      end if
      !------------------------------------------------------------------------------------!

      call fuse_2_cohorts_state(cpatch%costate, recc, donc, newni, cpatch%cophen%phenology_status(recc))
      call fuse_2_cohorts_phen(cpatch%cophen, recc, donc, newni, cpatch%costate%nplant(recc),  &
           cpatch%costate%nplant(donc))
      call fuse_2_cohorts_mort(cpatch%comort, recc, donc, newni, cpatch%costate%nplant(recc),  &
           cpatch%costate%nplant(donc), n_mort)
      call fuse_2_cohorts_resp(cpatch%coresp, recc, donc, newni, cpatch%costate%nplant(recc),  &
           cpatch%costate%nplant(donc),ndcycle)
      call fuse_2_cohorts_photo(cpatch%cophoto, recc, donc, newni, cpatch%costate%nplant(recc),  &
           cpatch%costate%nplant(donc),cpatch%costate%lai(recc), cpatch%costate%nplant(donc), &
           ndcycle, newlaii)
      call fuse_2_cohorts_therm(cpatch%cotherm, recc, donc, newni, cpatch%costate%nplant(recc), &
           cpatch%costate%nplant(donc))

      cpatch%costate%krdepth(recc) = dbh2krdepth(cpatch%costate%hite(recc),  &
           cpatch%costate%dbh(recc),cpatch%costate%pft(recc),lsl)

      !     Lastly, we update nplant and LAI.                                              !
      cpatch%costate%nplant(recc) = newn

      !------------------------------------------------------------------------------------!
      !    LAI must be zero if phenology status is 2.  This is probably done correctly     !
      ! throughout the code, but being safe here.                                          !
      !------------------------------------------------------------------------------------!
      cpatch%costate%lai(recc) = cpatch%costate%lai(recc) + cpatch%costate%lai(donc)
      !------------------------------------------------------------------------------------!

      return
   end subroutine fuse_2_cohorts

   !=======================================================================================!
   !   This subroutine will sort the patches by age (1st = oldest, last = youngest.)       !
   !---------------------------------------------------------------------------------------!
   subroutine sort_patches(csite)

      use ed_state_vars, only  :  sitetype   ! ! Structure
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype), target   :: csite      ! Current site, that will have patches sorted.
      !----- Local variables --------------------------------------------------------------!
      type(sitetype), pointer  :: tempsite   ! Structure to temporarily host the sorted site
      integer                  :: ipa        ! Counters
      integer                  :: oldpa      ! Index of oldest patch
      logical                  :: sorted     ! Flag: the site is already sorted
      !------------------------------------------------------------------------------------!
      
      !----- No need to sort a site with a single patch. ----------------------------------!
      if (csite%npatches < 2) return

      !------------------------------------------------------------------------------------!
      !     Check whether this site is already sorted.   We don't want to do the entire    !
      ! deallocating/copying/allocating thing if it's not needed as this takes up too much !
      ! time.                                                                              !
      !------------------------------------------------------------------------------------!
      sorted = .true.
      sortcheck: do ipa=1,csite%npatches-1
         sorted = csite%age(ipa) >= csite%age(ipa+1)
         if (.not. sorted) exit sortcheck
      end do sortcheck
      if (sorted) return
      !------------------------------------------------------------------------------------!



      !----- Assign a scratch patch. ------------------------------------------------------!
      nullify (tempsite)
      allocate(tempsite)
      call allocate_sitetype(tempsite,csite%npatches)
      
      ipa = 0
      !---- Loop until all patches were sorted. -------------------------------------------!
      do while (ipa < csite%npatches)
         ipa = ipa + 1
      
         !----- Find the oldest site. -----------------------------------------------------!
         oldpa = maxloc(csite%age,dim=1)
         
         !----- Copy to patch the scratch structure. --------------------------------------!
         call copy_sitetype(csite,tempsite,oldpa,oldpa,ipa,ipa)
         
         !----- Put a non-sense age so this patch will never "win" again. -----------------!
         csite%age(oldpa) = -huge(1.)
      end do

      !------ Reset the actual patch, and re-allocate it. ---------------------------------!
      call deallocate_sitetype(csite)
      call allocate_sitetype  (csite,tempsite%npatches)

      !------ Copy the scratch patch to the regular one. ----------------------------------!
      call copy_sitetype(tempsite,csite,1,tempsite%npatches,1,tempsite%npatches)
      call deallocate_sitetype(tempsite)
      deallocate(tempsite)

      return

   end subroutine sort_patches
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine will perform patch fusion based on some similarity criteria to      !
   ! determine whether they can be fused with no significant loss of information. The user !
   ! is welcome to set up a benchmark, but they should be aware that no miracles will      !
   ! happen here. If there are more very distinct patches than maxpatch, then the user     !
   ! will need to live with that and accept life is not always fair with those with        !
   ! limited computational resources.                                                      !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_patches(cgrid,ifm)
      use ed_state_vars       , only : edtype              & ! structure
                                     , polygontype         & ! structure
                                     , sitetype            & ! structure
                                     , patchtype           ! ! structure
      use fusion_fission_coms , only : ff_nhgt             & ! intent(in)
                                     , niter_patfus        & ! intent(in)
                                     , dark_cumlai_min     & ! intent(in)
                                     , dark_cumlai_max     & ! intent(in)
                                     , dark_cumlai_mult    & ! intent(in)
                                     , sunny_cumlai_min    & ! intent(in)
                                     , sunny_cumlai_max    & ! intent(in)
                                     , sunny_cumlai_mult   & ! intent(in)
                                     , print_fuse_details  & ! intent(in)
                                     , light_toler_min     & ! intent(in)
                                     , light_toler_max     & ! intent(in)
                                     , light_toler_mult    & ! intent(in)
                                     , fuse_prefix         ! ! intent(in)
      use ed_max_dims         , only : n_pft               & ! intent(in)
                                     , str_len             ! ! intent(in)
      use mem_polygons        , only : maxpatch            & ! intent(in)
                                     , maxcohort           ! ! intent(in)
      use ed_node_coms        , only : mynum               ! ! intent(in)
      use ed_misc_coms        , only : current_time        ! ! intent(in)
      use grid_coms           , only : nzg                 & ! intent(in)
                                     , nzs                 ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)          , target      :: cgrid           ! Current grid
      integer               , intent(in)  :: ifm             ! Current grid index
      !----- Local variables --------------------------------------------------------------!
      type(polygontype)     , pointer     :: cpoly           ! Current polygon
      type(polygontype)     , pointer     :: jpoly           ! Current polygon
      type(sitetype)        , pointer     :: csite           ! Current site
      type(patchtype)       , pointer     :: cpatch          ! Current patch
      type(patchtype)       , pointer     :: donpatch        ! Donor patch
      type(patchtype)       , pointer     :: recpatch        ! Receptor patch
      type(sitetype)        , pointer     :: tempsite        ! Temporary site
      logical, dimension(:) , allocatable :: fuse_table      ! Flag: this will remain.
      character(len=str_len)              :: fuse_fout       ! Filename for detailed output
      real   , dimension(ff_nhgt)         :: laimax_cum      ! Mean # of plants
      integer                             :: ipy             ! Counters
      integer                             :: isi             ! Counters
      integer                             :: jpy             ! Counters
      integer                             :: jsi             ! Counters
      integer                             :: ipa             ! Counters
      integer                             :: ico             ! Counters
      integer                             :: donp            ! Counters
      integer                             :: recp            ! Counters
      integer                             :: ipft            ! Counters
      integer                             :: ihgt            ! Counters
      integer                             :: ifus            ! Counters
      integer                             :: npatches_new    ! New # of patches
      integer                             :: npatches_old    ! Old # of patches
      integer                             :: npatches_orig   ! Original # of patches
      logical                             :: fuse_flag       ! Flag: fusion will happen
      logical                             :: recp_found      ! Found a receptor candidate
      logical                             :: sunny_donp      ! Donor patch bin too sunny
      logical                             :: sunny_recp      ! Receptor patch bin too sunny
      logical                             :: dark_donp       ! Donor patch bin too small
      logical                             :: dark_recp       ! Receptor patch bin too small
      logical                             :: same_age        ! Patches with same age
      real                                :: diff            ! Absolute difference in prof.
      real                                :: refv            ! Reference value of bin
      real                                :: norm            ! Normalised difference
      real                                :: llevel_donp     ! Light level of donor patch
      real                                :: llevel_recp     ! Light level of rec.  patch
      real                                :: sunny_toler     ! Light layer tolerance.
      real                                :: dark_lai80      ! Minimum dark layer.
      real                                :: dark_toler      ! Dark layer tolerance.
      real                                :: light_toler     ! Light level Relative toler.
      real                                :: old_area        ! For area conservation check
      real                                :: new_area        ! For area conservation check
      real                                :: old_lai_tot     ! Old total LAI
      real                                :: old_nplant_tot  ! Old total nplant
      real                                :: new_lai_tot     ! New total LAI
      real                                :: new_nplant_tot  ! New total nplant
      real                                :: elim_nplant     ! Elim. nplant during 1 fusion
      real                                :: elim_lai        ! Elim. LAI during 1 fusion
      real                                :: elim_nplant_tot ! Total eliminated nplant
      real                                :: elim_lai_tot    ! Elim. eliminated LAI
      real                                :: cumlai_recp     ! Cumulative LAI (receptor)
      real                                :: cumlai_donp     ! Cumulative LAI (donor)
      integer                             :: tot_npolygons   ! Total # of polygons
      integer                             :: tot_nsites      ! Total # of sites
      integer                             :: tot_npatches    ! Total # of patches
      integer                             :: tot_ncohorts    ! Total # of cohorts
      !----- Locally saved variables. --------------------------------------------------------!
      logical                   , save    :: first_time = .true.
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     First time here.  Delete all files.                                            !
      !------------------------------------------------------------------------------------!
      if (first_time .and. print_fuse_details) then
         do jpy = 1, cgrid%npolygons
            jpoly => cgrid%polygon(jpy)
            do jsi = 1, jpoly%nsites
               write (fuse_fout,fmt='(a,2(a,i4.4),a)')                                     &
                     trim(fuse_prefix),'polygon_',jpy,'_site_',jsi,'.txt'
               open (unit=72,file=trim(fuse_fout),status='replace',action='write')
               write(unit=72,fmt='(a)')       '----------------------------------------'
               write(unit=72,fmt='(a)')       ' Patch Fusion log for: '
               write(unit=72,fmt='(a,1x,i5)') ' POLYGON: ',jpy 
               write(unit=72,fmt='(a,1x,i5)') ' SITE:    ',jsi 
               write(unit=72,fmt='(a)')       '----------------------------------------'
               write(unit=72,fmt='(a)')       ' '
               close(unit=72,status='keep')
            end do
         end do
         first_time = .false.
      end if
      !---------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Return if maxpatch is 0, this is a flag for no patch fusion.                   !
      !------------------------------------------------------------------------------------!
      if (maxpatch == 0) return
      !------------------------------------------------------------------------------------!

      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         
         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            write(fuse_fout,fmt='(a,2(a,i4.4),a)')                                         &
                     trim(fuse_prefix),'polygon_',ipy,'_site_',isi,'.txt'

            if (print_fuse_details) then
               open (unit=72,file=trim(fuse_fout),status='old',action='write'              &
                                                 ,position='append')

               write(unit=72,fmt='(2(a,i2.2),a,i4.4)')    ' - Date: ',current_time%month   &
                                                                 ,'/',current_time%date    &
                                                                 ,'/',current_time%year
               write(unit=72,fmt='(a,1x,i6)') '   + Initial number of patches: '           &
                                             ,csite%npatches
               write(unit=72,fmt='(a)')       '   + Looking for empty patches: '
               close(unit=72,status='keep')
            end if

            !----- Skip this site if it contains only one patch... ------------------------!
            if (csite%npatches < 2) cycle siteloop

            !----- Save original number of patches. ---------------------------------------!
            npatches_orig = csite%npatches

            !----- Allocate the swapper patches in the site type. -------------------------!
            nullify(tempsite)
            allocate(tempsite)
            call allocate_sitetype(tempsite, csite%npatches)
            allocate(fuse_table(csite%npatches))

            !------------------------------------------------------------------------------!
            !     Allocate the fusion flag vector, and set all elements to .true., which   !
            ! means that every patch can be fused.  As soon as the patch is fused, we will !
            ! switch the flag to false.                                                    !
            !------------------------------------------------------------------------------!
            fuse_table(:) = .true.
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Find the original number of plants, LAI, and area, which will be used    !
            ! for sanity check.  We also compute the pft and size profile for all patches, !
            ! which will be used for the fusion criterion.                                 !
            !------------------------------------------------------------------------------!
            old_nplant_tot = 0.
            old_lai_tot    = 0.
            old_area       = 0.
            do ipa = 1,csite%npatches
               call patch_pft_size_profile(csite,ipa)

               old_area  = old_area + csite%area(ipa)
               cpatch => csite%patch(ipa)
               do ico = 1, cpatch%ncohorts
                  old_nplant_tot = old_nplant_tot + cpatch%costate%nplant(ico) * csite%area(ipa)
                  old_lai_tot    = old_lai_tot    + cpatch%costate%lai(ico)    * csite%area(ipa)
               end do
            end do
            !------------------------------------------------------------------------------!



            !----- Initialise the total eliminated nplant and LAI to zero. ----------------!
            elim_nplant_tot = 0.
            elim_lai_tot    = 0.
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     This loop will check whether there is at least two patches with exact    !
            ! same age.  This is common in initialisation with site-level measurements.    !
            ! In this case we want to merge any patch, regardless on their "age position", !
            ! which doesn't mean anything in this case.                                    !
            !------------------------------------------------------------------------------!
            same_age=.false.
            donloop_check: do donp=csite%npatches,2,-1
               recloop_check: do recp=donp-1,1,-1
                  same_age = csite%age(recp)       == csite%age(donp) .and.                &
                             csite%dist_type(recp) == csite%dist_type(donp)
                  if (print_fuse_details) then
                        open (unit=72,file=trim(fuse_fout),status='old',action='write'     &
                                              ,position='append')
                        write(unit=72,fmt='(a,1x,l1)') '     * same_age is ',same_age
                        close(unit=72,status='keep')
                  end if
                  !----- At least two patches have the same age. --------------------------!
                  if (same_age) exit donloop_check
              end do recloop_check
            end do donloop_check
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            donloope: do donp=csite%npatches,2,-1
               donpatch => csite%patch(donp)
               
               !----- If patch is not empty, or has already been fused, move on. ----------!
               if ( (.not. fuse_table(donp)) .or. donpatch%ncohorts > 0) cycle donloope


               !---------------------------------------------------------------------------!
               !     If we reach this point, it means that the donor patch is empty and    !
               ! hasn't been fused yet: look for an older empty patch and merge them.      !
               !---------------------------------------------------------------------------!
               if (print_fuse_details) then
                  open (unit=72,file=trim(fuse_fout),status='old',action='write'           &
                                                    ,position='append')
                  write(unit=72,fmt='(a,i6,a)') '     * Patch ',donp,' is empty.'
                  close(unit=72,status='keep')
               end if
               recloope: do recp=donp-1,1,-1
                  recpatch => csite%patch(recp)

                  !------------------------------------------------------------------------!
                  !     Skip the patch if it isn't empty, or it has already been fused, or !
                  ! if the donor and receptor have different disturbance types.            !
                  !------------------------------------------------------------------------!
                  if ( (.not. fuse_table(recp))                       .or.                 &
                       recpatch%ncohorts > 0                          .or.                 &
                       csite%dist_type(donp) /= csite%dist_type(recp)     ) then
                     cycle recloope
                  end if
                  !------------------------------------------------------------------------!

                  !----- Skip the patch if they don't have the same disturbance type. -----!
                  if ( csite%dist_type(donp) /= csite%dist_type(recp)) cycle recloope
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !     Take an average of the patch properties of donpatch and recpatch,  !
                  ! and assign the average recpatch.                                       !
                  !------------------------------------------------------------------------!
                  call fuse_2_patches(csite,donp,recp,nzg,nzs,cpoly%met(isi)%prss          &
                                     ,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)               &
                                     ,cpoly%green_leaf_factor(:,isi),elim_nplant,elim_lai)


                  !----- Record the fusion if requested by the user. ----------------------!
                  if (print_fuse_details) then
                     open (unit=72,file=trim(fuse_fout),status='old',action='write'        &
                                                        ,position='append')
                     write(unit=72,fmt='(2(a,i6),a)') '     * Patches ',donp,' and ',recp  &
                                                     ,' were fused.'
                     close(unit=72,status='keep')
                  end if
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    Update total eliminated nplant and LAI.  This is actually not       !
                  ! necessary in this loop as both patches are empty, but we do it anyway  !
                  ! just to be consistent.                                                 !
                  !------------------------------------------------------------------------!
                  elim_nplant_tot = elim_nplant_tot + elim_nplant
                  elim_lai_tot    = elim_lai_tot    + elim_lai
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Recalculate the pft size profile for the averaged patch at recp.   !
                  ! Again, this is not really necessary as the receptor patch is empty,    !
                  ! but just to be consistent...                                           !
                  !------------------------------------------------------------------------!
                  call patch_pft_size_profile(csite,recp)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     The patch at index donp will be eliminated and should not be       !
                  ! checked for fusion again; we switch the fuse_table flag to .false. so  !
                  ! next time we reach this patch we will skip it.                         !
                  !------------------------------------------------------------------------!
                  fuse_table(donp) = .false.
                  !------------------------------------------------------------------------!

                  !------ We are done with donp, so we quit the recp loop. ----------------!
                  exit recloope

               end do recloope
            end do donloope
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !     Second loop.  Here we will fuse all patches with the same age and        !
            ! disturbance type.                                                            !
            !------------------------------------------------------------------------------!

            if (same_age) then
               !----- Start with no multiplication factor. --------------------------------!
               dark_toler   = dark_cumlai_max
               sunny_toler  = sunny_cumlai_min
               light_toler  = light_toler_min

               mainfuseloopa: do ifus=0,niter_patfus
                  npatches_old = count(fuse_table)
                  npatches_new = npatches_old

                  !------------------------------------------------------------------------!
                  !    Inform that the upcoming fusions are going to be with populated     !
                  ! patches, and record the tolerance used.                                !
                  !------------------------------------------------------------------------!
                  if (print_fuse_details) then
                     open (unit=72,file=trim(fuse_fout),status='old',action='write'        &
                                                        ,position='append')
                     write(unit=72,fmt='(a,1x,3(a,1x,es9.2,1x))')                          &
                                  '   + Looking for similar populated patches of same age' &
                                 ,' - Sunny Tolerance =',sunny_toler                       &
                                 ,' - Dark Tolerance  =',dark_toler                        &
                                 ,' - Rel. Tolerance  =',light_toler
                     close(unit=72,status='keep')
                  end if
                  !------------------------------------------------------------------------!

                  donloopa: do donp=csite%npatches,2,-1
                     donpatch => csite%patch(donp)
                     
                     !----- If patch is not empty, or has already been fused, move on. ----!
                     if ( (.not. fuse_table(donp)) .or. donpatch%ncohorts == 0) then
                        cycle donloopa
                     end if


                     !---------------------------------------------------------------------!
                     !     If we reach this point, it means that the donor patch is        !
                     ! populated and hasn't been fused yet: look for other patches with    !
                     ! the same age and merge them.                                        !
                     !---------------------------------------------------------------------!
                     if (print_fuse_details) then
                        open (unit=72,file=trim(fuse_fout),status='old',action='write'     &
                                                          ,position='append')
                        write(unit=72,fmt='(a,i6,a)') '     * Patch ',donp,' is populated.'
                        close(unit=72,status='keep')
                     end if
                     recloopa: do recp=donp-1,1,-1
                        recpatch => csite%patch(recp)

                        !------------------------------------------------------------------!
                        !     Skip the patch if it isn't empty, or it has already been     !
                        ! fused, or if the donor and receptor have different disturbance   !
                        ! types.                                                           !
                        !------------------------------------------------------------------!
                        if ( (.not. fuse_table(recp))                       .or.           &
                             recpatch%ncohorts == 0                         .or.           &
                             csite%dist_type(donp) /= csite%dist_type(recp) .or.           &
                             csite%age(donp)      /=  csite%age(recp)    ) then
                           cycle recloopa
                        end if
                        !------------------------------------------------------------------!

                        !------------------------------------------------------------------!
                        !     Find the LAI that corresponds to 80% of the maximum LAI, to  !
                        ! avoid relaxing too much for forests.                             !
                        !------------------------------------------------------------------!
                        dark_lai80 = 0.40 * ( sum(csite%cumlai_profile(:,1,recp))          &
                                            + sum(csite%cumlai_profile(:,1,donp)) )


                        !------------------------------------------------------------------!
                        !     Compare the size profile for each PFT.  Here we compare the  !
                        ! maximum LAI for each PFT and height bin.  We switched the        !
                        ! classes from DBH to height because different PFTs may have       !
                        ! different heights for a given DBH, so we want to make sure the   !
                        ! light profile is right.                                          !
                        !------------------------------------------------------------------!
                        hgtloopa: do ihgt=1,ff_nhgt
                           cumlai_recp = sum(csite%cumlai_profile(:,ihgt,recp))
                           cumlai_donp = sum(csite%cumlai_profile(:,ihgt,donp))

                           !---------------------------------------------------------------!
                           !    Check whether these bins contain some LAI.  We don't need  !
                           ! to check the cohorts if the understorey is too dark, so once  !
                           ! both patches becomes very dark (very high LAI), we stop       !
                           ! checking the profiles.                                        !
                           !---------------------------------------------------------------!
                           dark_donp = cumlai_donp > dark_toler
                           dark_recp = cumlai_recp > dark_toler
                           !---------------------------------------------------------------!



                           !---------------------------------------------------------------!
                           if (dark_recp .and. dark_donp) then

                              if (print_fuse_details) then
                                 open  (unit=72,file=trim(fuse_fout),status='old'          &
                                               ,action='write',position='append')
                                 write (unit=72,fmt='(1(a,1x,i6,1x),4(a,1x,es9.2,1x)'//    &
                                                    ',2(a,1x,l1,1x))')                     &
                                    '       * IHGT =',ihgt                                 &
                                            ,'CUMLAI_RECP =',cumlai_recp                   &
                                            ,'CUMLAI_DONP =',cumlai_donp                   &
                                            ,'DARK_TOLER =',dark_toler                     &
                                            ,'DARK_LAI80 =',dark_lai80                     &
                                            ,'DARK_RECP =',dark_recp                       &
                                            ,'DARK_DONP =',dark_donp
                                 close (unit=72,status='keep')
                              end if
                              cycle hgtloopa
                           end if
                           !---------------------------------------------------------------!



                           
                           !---------------------------------------------------------------!
                           !    Check whether these bins contain some LAI.  Bins that have !
                           ! tiny cumulative LAI may differ by a lot in the relative       !
                           ! scale, but the actual value is so small that we don't really  !
                           ! care whether they are relatively different.                   !
                           !---------------------------------------------------------------!
                           sunny_donp = cumlai_donp <= sunny_toler
                           sunny_recp = cumlai_recp <= sunny_toler
                           !---------------------------------------------------------------!





                           !---------------------------------------------------------------!
                           !    If both patches have little or no biomass in this bin,     !
                           ! don't even bother checking the difference.                    !
                           !---------------------------------------------------------------!
                           if (sunny_donp .and. sunny_recp) then
                              if (print_fuse_details) then
                                 open  (unit=72,file=trim(fuse_fout),status='old'          &
                                               ,action='write',position='append')
                                 write (unit=72,fmt='(1(a,1x,i6,1x),3(a,1x,es9.2,1x)'//    &
                                                    ',2(a,1x,l1,1x))')                     &
                                    '       * IHGT=',ihgt                                  &
                                            ,'CUMLAI_RECP =',cumlai_recp                   &
                                            ,'CUMLAI_DONP =',cumlai_donp                   &
                                            ,'SUNNY_TOLER =',sunny_toler                   &
                                            ,'SUNNY_RECP =',sunny_recp                     &
                                            ,'SUNNY_RECP =',sunny_donp
                                 close (unit=72,status='keep')
                              end if
                              cycle hgtloopa
                           end if
                           !---------------------------------------------------------------!


                           !---------------------------------------------------------------!
                           !    Find the normalised difference in the density of this PFT  !
                           ! and size.  If one of the patches is missing any member of the !
                           ! profile the norm will be set to 2.0, which is the highest     !
                           ! value that the norm can be.                                   !
                           !---------------------------------------------------------------!
                           llevel_donp = exp(- 0.5 * cumlai_donp)
                           llevel_recp = exp(- 0.5 * cumlai_recp)
                           
                           diff = abs(llevel_donp - llevel_recp )
                           refv =    (llevel_donp + llevel_recp ) * 0.5
                           norm = diff / refv
                           fuse_flag = norm <= light_toler
                           !---------------------------------------------------------------!



                           !---------------------------------------------------------------!
                           if (print_fuse_details) then
                              open  (unit=72,file=trim(fuse_fout),status='old'             &
                                            ,action='write',position='append')
                              write (unit=72,fmt='(1(a,1x,i6,1x),7(a,1x,es9.2,1x)'//       &
                                                 ',1(a,1x,l1,1x))')                        &
                                 '       * IHGT=',ihgt                                     &
                                ,'CLAI_RECP =',cumlai_recp,'CLAI_DONP =',cumlai_donp       &
                                ,'LL_RECP =',llevel_recp,'LL_DONP =',llevel_donp           &
                                ,'DIFF =',diff,'REFV =',refv,'NORM =',norm                 &
                                ,'FUSE_FLAG =',fuse_flag
                              close (unit=72,status='keep')
                           end if
                           !---------------------------------------------------------------!



                           !---------------------------------------------------------------!
                           !     If fuse_flag is false, the patches aren't similar, move   !
                           ! to the next donor patch.                                      !
                           !---------------------------------------------------------------!
                           if (.not. fuse_flag) cycle recloopa
                           !---------------------------------------------------------------!
                        end do hgtloopa
                        !------------------------------------------------------------------!
                        !------------------------------------------------------------------!
                        !     Take an average of the patch properties of donpatch and      !
                        ! recpatch, and assign the average recpatch.                       !
                        !------------------------------------------------------------------!
                        call fuse_2_patches(csite,donp,recp,nzg,nzs,cpoly%met(isi)%prss    &
                                           ,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)         &
                                           ,cpoly%green_leaf_factor(:,isi),elim_nplant     &
                                           ,elim_lai)


                        !----- Record the fusion if requested by the user. ----------------!
                        if (print_fuse_details) then
                           open (unit=72,file=trim(fuse_fout),status='old',action='write'  &
                                                              ,position='append')
                           write(unit=72,fmt='(2(a,i6),a)') '     * Patches ',donp,' and ' &
                                                                   ,recp,' were fused.'
                           close(unit=72,status='keep')
                        end if
                        !------------------------------------------------------------------!


                        !------------------------------------------------------------------!
                        !    Update total eliminated nplant and LAI.  This is actually not !
                        ! necessary in this loop as both patches are empty, but we do it   !
                        ! anyway just to be consistent.                                    !
                        !------------------------------------------------------------------!
                        elim_nplant_tot = elim_nplant_tot + elim_nplant
                        elim_lai_tot    = elim_lai_tot    + elim_lai
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Recalculate the pft size profile for the averaged patch at   !
                        ! recp.  Again, this is not really necessary as the receptor patch !
                        ! is empty, but just to be consistent...                           !
                        !------------------------------------------------------------------!
                        call patch_pft_size_profile(csite,recp)
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     The patch at index donp will be eliminated and should not be !
                        ! checked for fusion again; we switch the fuse_table flag to       !
                        ! .false. so next time we reach this patch we will skip it.        !
                        !------------------------------------------------------------------!
                        fuse_table(donp) = .false.
                        !------------------------------------------------------------------!
                        !------------------------------------------------------------------!
                        !     Update the number of valid patches.                          !
                        !------------------------------------------------------------------!
                        npatches_new = npatches_new - 1
                        !------------------------------------------------------------------!

                        !------ We are done with donp, so we quit the recp loop. ----------!
                        exit recloopa
                     end do recloopa
                  end do donloopa


                  !------------------------------------------------------------------------!
                  !      Check how many patches are valid.  If the total number of patches !
                  ! is less than the target, or if we have reached the maximum tolerance   !
                  ! and the patch fusion still can't find similar patches, we quit the     !
                  ! fusion loop.                                                           !
                  !------------------------------------------------------------------------!
                  if (npatches_new <= abs(maxpatch)) exit mainfuseloopa
                  !------------------------------------------------------------------------!

                  !----- Increment tolerance ----------------------------------------------!
                  sunny_toler =     sunny_toler * sunny_cumlai_mult
                  dark_toler  = max(dark_toler  * dark_cumlai_mult , dark_lai80 )
                  light_toler =     light_toler * light_toler_mult
                  !------------------------------------------------------------------------!
               end do mainfuseloopa
            end if
            !------------------------------------------------------------------------------!





            !------------------------------------------------------------------------------!
            !     Third patch loop. Now that all empty patches have been fused, we will    !
            ! look for populated patches that have similar size and PFT structure, using   !
            ! the following algorithm:                                                     !
            !                                                                              !
            ! 1. Loop from the youngest to oldest patch;                                   !
            ! 2. Find next older patch with same dist_type;                                !
            ! 3. Check whether fusion criterion is met, and If so, fuse them.              !
            ! 4. After all the fusion, check how many patches we still have:               !
            !    If it is less than maxpatch, we quit, otherwise, we relax the tolerance a !
            !    little and try fusing more.  Notice that we will always try fusing        !
            !    patches at least once, even when the original number is less than         !
            !    maxpatch.                                                                 !
            !------------------------------------------------------------------------------!
            !----- Start with no multiplication factor. -----------------------------------!
            dark_toler   = dark_cumlai_max
            sunny_toler  = sunny_cumlai_min
            light_toler  = light_toler_min

            mainfuseloop: do ifus=0,niter_patfus
               npatches_old = count(fuse_table)
               npatches_new = npatches_old

               !---------------------------------------------------------------------------!
               !    Inform that the upcoming fusions are going to be with populated        !
               ! patches, and record the tolerance used.                                   !
               !---------------------------------------------------------------------------!
               if (print_fuse_details) then
                  open (unit=72,file=trim(fuse_fout),status='old',action='write'           &
                                                     ,position='append')
                  write(unit=72,fmt='(a,1x,3(a,1x,es9.2,1x))')                             &
                                              '   + Looking for similar populated patches' &
                                             ,' - Sunny Tolerance =',sunny_toler           &
                                             ,' - Dark Tolerance  =',dark_toler            &
                                             ,' - Rel. Tolerance  =',light_toler
                  close(unit=72,status='keep')
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Loop from youngest to the second oldest patch.                        !
               !---------------------------------------------------------------------------!
               donloopp: do donp = csite%npatches,2,-1
                  donpatch => csite%patch(donp)

                  !------------------------------------------------------------------------!
                  !     If this is an empty patch, or has already been merged, we skip it. !
                  !------------------------------------------------------------------------!
                  if ((.not. fuse_table(donp))) then
                     if (print_fuse_details) then
                        open  (unit=72,file=trim(fuse_fout),status='old',action='write'    &
                                                           ,position='append')
                        write (unit=72,fmt='(a,1x,i6,1x,a)') '     - DONP:',donp           &
                                                            ,'has been already fused...'
                        close (unit=72,status='keep')
                     end if
                     cycle donloopp
                  end if
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      If we have reached this place, the donor patch can be fused.  Now !
                  ! look for the next oldest patch that has the same disturbance type.  In !
                  ! case we can't find such patch, we will move to the next donor          !
                  ! candidate.                                                             !
                  !------------------------------------------------------------------------!
                  recp_found = .false.
                  recloopp: do recp=donp-1,1,-1
                     recp_found = csite%dist_type(donp) == csite%dist_type(recp) .and.     &
                                  fuse_table(recp) .and.                                   &
                                  (csite%dist_type(recp) == 1 .or. csite%age(recp) > 2.)
                     if (recp_found) then
                        recpatch => csite%patch(recp)
                        exit recloopp
                     end if
                  end do recloopp
                  
                  if (.not. recp_found) then
                     if (print_fuse_details) then
                        open  (unit=72,file=trim(fuse_fout),status='old',action='write'    &
                                                           ,position='append')
                        write (unit=72,fmt='(a)') '     - No receptor patch found. '
                        close (unit=72,status='keep')
                     end if
                     cycle donloopp
                  end if
                  !------------------------------------------------------------------------!


                  if (print_fuse_details) then
                     open  (unit=72,file=trim(fuse_fout),status='old',action='write'       &
                                                        ,position='append')
                     write (unit=72,fmt='(2(a,1x,i6,1x))') '     - DONP =',donp            &
                                                                 ,'RECP =',recp
                     close (unit=72,status='keep')
                  end if


                  !------------------------------------------------------------------------!
                  !     This should never happen because we have already fused all empty   !
                  ! patches, but, just in case... If both patches are empty they cannot be !
                  ! fused in this loop.                                                    !
                  !------------------------------------------------------------------------!
                  if (donpatch%ncohorts == 0 .and. recpatch%ncohorts == 0) cycle donloopp
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Find the LAI that corresponds to 80% of the maximum LAI, to avoid  !
                  ! relaxing too much for forests.                                         !
                  !------------------------------------------------------------------------!
                  dark_lai80 = 0.40 * ( sum(csite%cumlai_profile(:,1,recp))                &
                                      + sum(csite%cumlai_profile(:,1,donp)) )


                  !------------------------------------------------------------------------!
                  !     Compare the size profile for each PFT.  Here we compare the        !
                  ! maximum LAI for each PFT and height bin.  We switched the classes from !
                  ! DBH to height because different PFTs may have different heights for a  !
                  ! given DBH, so we want to make sure the light profile is right.         !
                  !------------------------------------------------------------------------!
                  hgtloop: do ihgt=1,ff_nhgt
                     cumlai_recp = sum(csite%cumlai_profile(:,ihgt,recp))
                     cumlai_donp = sum(csite%cumlai_profile(:,ihgt,donp))

                     !---------------------------------------------------------------------!
                     !    Check whether these bins contain some LAI.  We don't need to     !
                     ! check the cohorts if the understorey is too dark, so once both      !
                     ! patches becomes very dark (very high LAI), we stop checking the     !
                     ! profiles.                                                           !
                     !---------------------------------------------------------------------!
                     dark_donp = cumlai_donp > dark_toler
                     dark_recp = cumlai_recp > dark_toler
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     if (dark_recp .and. dark_donp) then

                        if (print_fuse_details) then
                           open  (unit=72,file=trim(fuse_fout),status='old',action='write' &
                                                              ,position='append')
                           write (unit=72,fmt='(1(a,1x,i6,1x),4(a,1x,es9.2,1x)'//          &
                                              ',2(a,1x,l1,1x))')                           &
                              '       * IHGT=',ihgt                                        &
                             ,'CUMLAI_RECP =',cumlai_recp,'CUMLAI_DONP =',cumlai_donp      &
                             ,'DARK_TOLER =',dark_toler,'DARK_LAI80 =',dark_lai80          &
                             ,'DARK_RECP =',dark_recp,'DARK_DONP =',dark_donp
                           close (unit=72,status='keep')
                        end if
                        cycle hgtloop
                     end if
                     !---------------------------------------------------------------------!



                     
                     !---------------------------------------------------------------------!
                     !    Check whether these bins contain some LAI.  Bins that have       !
                     ! tiny cumulative LAI may differ by a lot in the relative scale,      !
                     ! but the actual value is so small that we don't really care          !
                     ! whether they are relatively different.                              !
                     !---------------------------------------------------------------------!
                     sunny_donp = cumlai_donp <= sunny_toler
                     sunny_recp = cumlai_recp <= sunny_toler
                     !---------------------------------------------------------------------!





                     !---------------------------------------------------------------------!
                     !    If both patches have little or no biomass in this bin, don't     !
                     ! even bother checking the difference.                                !
                     !---------------------------------------------------------------------!
                     if (sunny_donp .and. sunny_recp) then
                        if (print_fuse_details) then
                           open  (unit=72,file=trim(fuse_fout),status='old',action='write' &
                                                              ,position='append')
                           write (unit=72,fmt='(1(a,1x,i6,1x),3(a,1x,es9.2,1x)'//          &
                                              ',2(a,1x,l1,1x))')                           &
                              '       * IHGT=',ihgt                                        &
                             ,'CUMLAI_RECP =',cumlai_recp,'CUMLAI_DONP =',cumlai_donp      &
                             ,'SUNNY_TOLER =',sunny_toler,'SUNNY_RECP =',sunny_recp        &
                             ,'SUNNY_RECP =',sunny_donp
                           close (unit=72,status='keep')
                        end if
                        cycle hgtloop
                     end if
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !    Find the normalised difference in the density of this PFT and    !
                     ! size.  If one of the patches is missing any member of the           !
                     ! profile the norm will be set to 2.0, which is the highest value     !
                     ! that the norm can be.                                               !
                     !---------------------------------------------------------------------!
                     llevel_donp = exp(- 0.5 * cumlai_donp)
                     llevel_recp = exp(- 0.5 * cumlai_recp)
                     
                     diff = abs(llevel_donp - llevel_recp )
                     refv =    (llevel_donp + llevel_recp ) * 0.5
                     norm = diff / refv
                     fuse_flag = norm <= light_toler
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     if (print_fuse_details) then
                        open  (unit=72,file=trim(fuse_fout),status='old',action='write'    &
                                                           ,position='append')
                        write (unit=72,fmt='(1(a,1x,i6,1x),7(a,1x,es9.2,1x)'//             &
                                           ',1(a,1x,l1,1x))')                              &
                           '       * IHGT=',ihgt                                           &
                          ,'CLAI_RECP =',cumlai_recp,'CLAI_DONP =',cumlai_donp             &
                          ,'LL_RECP =',llevel_recp,'LL_DONP =',llevel_donp                 &
                          ,'DIFF =',diff,'REFV =',refv,'NORM =',norm                       &
                          ,'FUSE_FLAG =',fuse_flag
                        close (unit=72,status='keep')
                     end if
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !     If fuse_flag is false, the patches aren't similar, move to      !
                     ! the next donor patch.                                               !
                     !---------------------------------------------------------------------!
                     if (.not. fuse_flag) cycle donloopp
                     !---------------------------------------------------------------------!
                  end do hgtloop
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Reaching this point means that the patches are sufficiently       !
                  ! similar so they will be fused.   We take the average of the patch      !
                  ! properties of donpatch and recpatch, and leave the averaged values at  !
                  ! recpatch.                                                              !
                  !------------------------------------------------------------------------!
                  call fuse_2_patches(csite,donp,recp,nzg,nzs,cpoly%met(isi)%prss          &
                                     ,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)               &
                                     ,cpoly%green_leaf_factor(:,isi),elim_nplant,elim_lai)
                  !------------------------------------------------------------------------!


                  !----- Record the fusion if requested by the user. ----------------------!
                  if (print_fuse_details) then
                     open (unit=72,file=trim(fuse_fout),status='old',action='write'        &
                                                        ,position='append')
                     write(unit=72,fmt='(2(a,i6),a)') '     * Patches ',donp,' and ',recp  &
                                                     ,' were fused.'
                     close(unit=72,status='keep')
                  end if
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !    Some cohorts may have been eliminated during the fusion process,    !
                  ! because they were way too small.  Add the eliminated plant density and !
                  ! LAI because we want to make sure that the fusion routine conserves the !
                  ! total plant density and LAI that remained in the polygon.              !
                  !------------------------------------------------------------------------!
                  elim_nplant_tot = elim_nplant_tot + elim_nplant
                  elim_lai_tot    = elim_lai_tot    + elim_lai
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Recalculate the pft size profile for the updated receptor patch.   !
                  !------------------------------------------------------------------------!
                  call patch_pft_size_profile(csite,recp)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     From now on donpatch should not be used: set fuse_table flag as    !
                  ! .false. so we won't check it again.                                    !
                  !------------------------------------------------------------------------!
                  fuse_table(donp) = .false.
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Update the number of valid patches.                                !
                  !------------------------------------------------------------------------!
                  npatches_new = npatches_new - 1
                  !------------------------------------------------------------------------!
               end do donloopp         ! do donp = csite%npatches,2,-1
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Check how many patches are valid.  If the total number of patches is !
               ! less than the target, or if we have reached the maximum tolerance and the !
               ! patch fusion still can't find similar patches, we quit the fusion loop.   !
               !---------------------------------------------------------------------------!
               if (npatches_new <= abs(maxpatch)) exit mainfuseloop
               !---------------------------------------------------------------------------!

               !----- Increment tolerance -------------------------------------------------!
               sunny_toler =     sunny_toler * sunny_cumlai_mult
               dark_toler  = max(dark_toler  * dark_cumlai_mult , dark_lai80 )
               light_toler =     light_toler * light_toler_mult
               !---------------------------------------------------------------------------!
            end do mainfuseloop

            !------------------------------------------------------------------------------!



            !----- Set the number of patches in the site to "npatches_new" ----------------!
            tempsite%npatches = npatches_new
            !------------------------------------------------------------------------------!



            !----- If there was any patch fusion, need to shrink csite --------------------!
            if (npatches_new < csite%npatches) then
               !---------------------------------------------------------------------------!
               !    Copy the selected data into the temporary space, args 1 and 3 must be  !
               ! dimension of arg 4. Argument 2 must be the dimension of the sum of the    !
               ! 3rd argument.                                                             !
               !---------------------------------------------------------------------------!
               call copy_sitetype_mask(csite,tempsite,fuse_table,size(fuse_table)          &
                                      ,npatches_new)
               call deallocate_sitetype(csite)

               !----- Reallocate the current site. ----------------------------------------!
               call allocate_sitetype(csite,npatches_new)

               !----- Copy the selected temporary data into the orignal site vectors. -----!
               call copy_sitetype(tempsite,csite,1,npatches_new,1,npatches_new)

               !---------------------------------------------------------------------------!
               !     The new and fused csite is now complete, clean up the temporary       !
               ! data. Deallocate it afterwards.                                           !
               !---------------------------------------------------------------------------!
               call deallocate_sitetype(tempsite)
            end if
            !------------------------------------------------------------------------------!



            !----- Deallocation should happen outside the "if" statement ------------------!
            deallocate(tempsite)
            deallocate(fuse_table)
            !------------------------------------------------------------------------------!


            !----- Make sure that patches are sorted from oldest to youngest. -------------!
            call sort_patches(csite)

            if (print_fuse_details) then
               open (unit=72,file=trim(fuse_fout),status='old',action='write'              &
                                                  ,position='append')
               write(unit=72,fmt='(a)')             '   + Patches were sorted. '
               write(unit=72,fmt='(2(a,1x,i6,1x))')                                        &
                                       '   + Number of patches.  Original =',npatches_orig &
                                                                ,'Current =',csite%npatches
               write(unit=72,fmt='(a)')       ' '
               close(unit=72,status='keep')
            end if


            !----- This is for mass conservation check ------------------------------------!
            new_nplant_tot = 0.
            new_lai_tot    = 0.
            new_area       = 0.
            do ipa = 1,csite%npatches
               new_area = new_area + csite%area(ipa)
               cpatch => csite%patch(ipa)
               do ico = 1, cpatch%ncohorts
                  new_nplant_tot = new_nplant_tot + cpatch%costate%nplant(ico)*csite%area(ipa)
                  new_lai_tot    = new_lai_tot    + cpatch%costate%lai(ico)*csite%area(ipa)
               end do
            end do
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Sanity check.  Except for the cohorts that were eliminated because they  !
            ! have become too small after fusion, the total plant count and LAI should be  !
            ! preserved.  In case something went wrong, we stop the run, it is likely to   !
            ! be a bug.                                                                    !
            !------------------------------------------------------------------------------!
            if (new_area       < 0.99 * old_area       .or.                                &
                new_area       > 1.01 * old_area       .or.                                &
                new_nplant_tot < 0.99 * (old_nplant_tot - elim_nplant_tot) .or.            &
                new_nplant_tot > 1.01 * (old_nplant_tot - elim_nplant_tot) .or.            &
                new_lai_tot    < 0.99 * (old_lai_tot    - elim_lai_tot   ) .or.            &
                new_lai_tot    > 1.01 * (old_lai_tot    - elim_lai_tot   )     ) then
               write (unit=*,fmt='(a,1x,es12.5)') 'OLD_AREA:       ',old_area
               write (unit=*,fmt='(a,1x,es12.5)') 'NEW_AREA:       ',new_area
               write (unit=*,fmt='(a,1x,es12.5)') 'NEW_LAI_TOT:    ',new_lai_tot
               write (unit=*,fmt='(a,1x,es12.5)') 'OLD_LAI_TOT:    ',old_lai_tot
               write (unit=*,fmt='(a,1x,es12.5)') 'ELIM_LAI_TOT:   ',elim_lai_tot
               write (unit=*,fmt='(a,1x,es12.5)') 'NEW_NPLANT_TOT: ',new_nplant_tot
               write (unit=*,fmt='(a,1x,es12.5)') 'OLD_NPLANT_TOT: ',old_nplant_tot
               write (unit=*,fmt='(a,1x,es12.5)') 'ELIM_NPLANT_TOT:',elim_nplant_tot
               call fatal_error('Conservation failed while fusing patches'                 &
                              &,'fuse_patches','fuse_fiss_utils.f90')
            end if
            !------------------------------------------------------------------------------!
            
         end do siteloop
      end do polyloop

      !------------------------------------------------------------------------------------!
      !     Print a banner to inform the user how many patches and cohorts exist.          !
      !------------------------------------------------------------------------------------!
      tot_npolygons = cgrid%npolygons
      tot_ncohorts  = 0
      tot_npatches  = 0
      tot_nsites    = 0
      do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         tot_nsites = tot_nsites + cpoly%nsites 
         do isi=1,cpoly%nsites
            csite => cpoly%site(isi)
            tot_npatches = tot_npatches + csite%npatches
            do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)
               tot_ncohorts = tot_ncohorts + cpatch%ncohorts
            end do
         end do
      end do
      write (unit=*,fmt='(6(a,1x,i8,1x))')                                                 &
        'Total count in node',mynum,'for grid',ifm,': POLYGONS=',tot_npolygons             &
       ,'SITES=',tot_nsites,'PATCHES=',tot_npatches,'COHORTS=',tot_ncohorts
      !------------------------------------------------------------------------------------!

      return
   end subroutine fuse_patches
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine will merge two patches into 1.                                      !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_patches(csite,donp,recp,mzg,mzs,prss,lsl,ntext_soil,green_leaf_factor &
                               ,elim_nplant,elim_lai)
      use ed_state_vars      , only : sitetype              & ! Structure 
                                    , patchtype             ! ! Structure
      use soil_coms          , only : soil                  ! ! intent(in), lookup table
      use ed_max_dims        , only : n_pft                 & ! intent(in)
                                    , n_dbh                 ! ! intent(in)
      use mem_polygons       , only : maxcohort             ! ! intent(in)
      use consts_coms        , only : cpi                   & ! intent(in)
                                    , cpor                  & ! intent(in)
                                    , p00                   ! ! intent(in)
      use therm_lib          , only : qwtk                  ! ! function
      use ed_misc_coms       , only : iqoutput              & ! intent(in)
                                    , idoutput              & ! intent(in)
                                    , imoutput              & ! intent(in)
                                    , ndcycle               ! ! intent(in)
      use soil_bgc, only: fuse_2_patches_sbgc
      use cohort_state, only: fuse_2_patches_costate
      use cohort_therm, only: fuse_2_patches_cotherm
      use cohort_photo, only: fuse_2_patches_cophoto
      use cohort_mort, only: fuse_2_patches_comort
      use cohort_resp, only: fuse_2_patches_coresp
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target      :: csite             ! Current site
      integer                , intent(in)  :: donp              ! Donating patch
      integer                , intent(in)  :: recp              ! Receptor patch
      integer                , intent(in)  :: lsl               ! Lowest soil level
      integer                , intent(in)  :: mzg               ! # of soil layers
      integer                , intent(in)  :: mzs               ! # of sfc. water layers
      integer, dimension(mzg), intent(in)  :: ntext_soil        ! Soil type
      real, dimension(n_pft) , intent(in)  :: green_leaf_factor ! Green leaf factor...
      real                   , intent(in)  :: prss              ! Sfc. air density
      real                   , intent(out) :: elim_nplant       ! Eliminated nplant 
      real                   , intent(out) :: elim_lai          ! Eliminated lai
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)        , pointer     :: cpatch            ! Current patch
      type(patchtype)        , pointer     :: temppatch         ! Temporary patch
      integer                              :: ico,iii,icyc      ! Counters
      integer                              :: ndc               ! # of cohorts - donp patch
      integer                              :: nrc               ! # of cohorts - recp patch
      real                                 :: newarea           ! new patch area
      real                                 :: newareai          ! 1./(new patch area)
      real                                 :: area_scale        ! Cohort rescaling factor.
      integer :: iday
      !------------------------------------------------------------------------------------!
     
      !------------------------------------------------------------------------------------!
      !     This function fuses the two patches specified in the argument. It fuses the    !
      ! first patch in the argument (the "donor" = donp ) into the second patch in the     !
      ! argument (the "recipient" = recp ), and frees the memory associated with the donor !
      ! patch.                                                                             !
      !------------------------------------------------------------------------------------!
    
      !----- The new area is simply the sum of each patch area. ---------------------------!
      newarea  = csite%area(donp) + csite%area(recp)
      newareai = 1.0/newarea

      !----- Assign eliminated LAI and nplant to zero (everything stays) ------------------!
      elim_nplant = 0.
      elim_lai    = 0.

      !----- We now take the weighted average, scale by the individual patch area. --------!
      csite%age(recp)                = newareai *                                          &
                                     ( csite%age(donp)                * csite%area(donp)   &
                                     + csite%age(recp)                * csite%area(recp) )

      call fuse_2_patches_sbgc(csite%sbgc, recp, donp, newareai, csite%area(donp), &
           csite%area(recp))


      csite%sum_dgd(recp)            = newareai *                                          &
                                     ( csite%sum_dgd(donp)            * csite%area(donp)   &
                                     + csite%sum_dgd(recp)            * csite%area(recp) )

      csite%sum_chd(recp)            = newareai *                                          &
                                     ( csite%sum_chd(donp)            * csite%area(donp)   &
                                     + csite%sum_chd(recp)            * csite%area(recp) )

      csite%can_co2(recp)            = newareai *                                          &
                                     ( csite%can_co2(donp)            * csite%area(donp)   &
                                     + csite%can_co2(recp)            * csite%area(recp) )

      csite%can_theta(recp)          = newareai *                                          &
                                     ( csite%can_theta(donp)          * csite%area(donp)   &
                                     + csite%can_theta(recp)          * csite%area(recp) )

      csite%can_theiv(recp)          = newareai *                                          &
                                     ( csite%can_theiv(donp)          * csite%area(donp)   &
                                     + csite%can_theiv(recp)          * csite%area(recp) )

      csite%can_prss(recp)           = newareai *                                          &
                                     ( csite%can_prss(donp)           * csite%area(donp)   &
                                     + csite%can_prss(recp)           * csite%area(recp) )

      csite%can_shv(recp)            = newareai *                                          &
                                     ( csite%can_shv(donp)            * csite%area(donp)   &
                                     + csite%can_shv(recp)            * csite%area(recp) )

      csite%can_depth(recp)          = newareai *                                          &
                                     ( csite%can_depth(donp)          * csite%area(donp)   &
                                     + csite%can_depth(recp)          * csite%area(recp) )

      csite%ggbare(recp)             = newareai *                                          &
                                     ( csite%ggbare(donp)             * csite%area(donp)   &
                                     + csite%ggbare(recp)             * csite%area(recp) )

      csite%ggnet(recp)              = newareai *                                          &
                                     ( csite%ggnet(donp)              * csite%area(donp)   &
                                     + csite%ggnet(recp)              * csite%area(recp) )

      csite%ggsoil(recp)             = newareai *                                          &
                                     ( csite%ggsoil(donp)             * csite%area(donp)   &
                                     + csite%ggsoil(recp)             * csite%area(recp) )


      !------------------------------------------------------------------------------------!
      !    There is no guarantee that there will be a minimum amount of mass in the tempo- !
      ! rary layer, nor is there any reason for both patches to have the same number of    !
      ! layers. In order to be safe, the fusion must happen in 5 stages.                   !
      !------------------------------------------------------------------------------------!
      !----- 1. Find the "extensive" sfcwater_energy (convert from J/kg to J/m2); ---------!
      do iii=1,csite%nlev_sfcwater(recp)
         csite%sfcwater_energy(iii,recp) = csite%sfcwater_energy(iii,recp)                 &
                                         * csite%sfcwater_mass(iii,recp)
      end do
      do iii=1,csite%nlev_sfcwater(donp)
         csite%sfcwater_energy(iii,donp) = csite%sfcwater_energy(iii,donp)                 &
                                         * csite%sfcwater_mass(iii,donp)
      end do
      !------------------------------------------------------------------------------------!
      ! 2. Squeeze all layers into one.  If needed, the layer will be split again next     !
      !    time the Runge-Kutta integrator is called.  After adding the value to the first !
      !    layer, discard the value.                                                       !
      !------------------------------------------------------------------------------------!
      do iii=2,csite%nlev_sfcwater(recp)
         csite%sfcwater_energy(1,recp) = csite%sfcwater_energy(1,recp)                     &
                                       + csite%sfcwater_energy(iii,recp)
         csite%sfcwater_depth(1,recp)  = csite%sfcwater_depth(1,recp)                      &
                                       + csite%sfcwater_depth(iii,recp)
         csite%sfcwater_mass(1,recp)   = csite%sfcwater_mass(1,recp)                       &
                                       + csite%sfcwater_mass(iii,recp)
         csite%sfcwater_energy(iii,recp) = 0.
         csite%sfcwater_depth(iii,recp)  = 0.
         csite%sfcwater_mass(iii,recp)   = 0.
      end do
      do iii=2,csite%nlev_sfcwater(donp)
         csite%sfcwater_energy(1,donp) = csite%sfcwater_energy(1,donp)                     &
                                       + csite%sfcwater_energy(iii,donp)
         csite%sfcwater_depth(1,donp)  = csite%sfcwater_depth(1,donp)                      &
                                       + csite%sfcwater_depth(iii,donp)
         csite%sfcwater_mass(1,donp)   = csite%sfcwater_mass(1,donp)                       &
                                       + csite%sfcwater_mass(iii,donp)
         csite%sfcwater_energy(iii,donp) = 0.
         csite%sfcwater_depth(iii,donp)  = 0.
         csite%sfcwater_mass(iii,donp)   = 0.
      end do
      !----- 3. Merge the patches; --------------------------------------------------------!
      if (csite%nlev_sfcwater(recp) > 0 .or. csite%nlev_sfcwater(donp) > 0) then
         csite%sfcwater_mass(1,recp)   = newareai *                                        &
                                         (csite%sfcwater_mass(1,recp)  * csite%area(recp)  &
                                         +csite%sfcwater_mass(1,donp)  * csite%area(donp)  )
         csite%sfcwater_depth(1,recp)  = newareai *                                        &
                                         (csite%sfcwater_depth(1,recp) * csite%area(recp)  &
                                         +csite%sfcwater_depth(1,donp) * csite%area(donp)  )
         csite%sfcwater_energy(1,recp) = newareai *                                        &
                                         (csite%sfcwater_energy(1,recp) * csite%area(recp) &
                                         +csite%sfcwater_energy(1,donp) * csite%area(donp) )
      else
         csite%sfcwater_mass(1,recp)   = 0.
         csite%sfcwater_depth(1,recp)  = 0.
         csite%sfcwater_energy(1,recp) = 0.
      end if
      !------------------------------------------------------------------------------------!
      ! 4. Converting energy back to J/kg;                                                 !
      ! 5. Finding temperature and liquid water fraction;                                  !
      !    (Both are done in new_patch_sfc_props).                                         !
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!


      !----- Merge soil energy and water. -------------------------------------------------!
      do iii=1,mzg
         csite%soil_energy(iii,recp)     = newareai *                                      &
                                         ( csite%soil_energy(iii,donp) * csite%area(donp)  &
                                         + csite%soil_energy(iii,recp) * csite%area(recp))

         csite%soil_water(iii,recp)      = newareai *                                      &
                                         ( csite%soil_water(iii,recp)  * csite%area(recp)  &
                                         + csite%soil_water(iii,donp)  * csite%area(donp))
         csite%current_paw(iii,recp)      = newareai *                                      &
                                         ( csite%current_paw(iii,recp)  * csite%area(recp)  &
                                         + csite%current_paw(iii,donp)  * csite%area(donp))
         do iday = 1, 10
            csite%past_paw(iii,iday,recp)      = newareai *                                      &
                                         ( csite%past_paw(iii,iday,recp)  * csite%area(recp)  &
                                         + csite%past_paw(iii,iday,donp)  * csite%area(donp))
         enddo
      end do

      !------------------------------------------------------------------------------------!
      !    This subroutine takes care of filling:                                          !
      !                                                                                    !
      ! + csite%ground_shv(recp)                                                           !
      ! + csite%ground_ssh(recp)                                                           !
      ! + csite%ground_temp(recp)                                                          !
      ! + csite%ground_fliq(recp)                                                          !
      ! + csite%soil_tempk(k,recp)                                                         !
      ! + csite%soil_fracliq(k,recp)                                                       !
      ! + csite%nlev_sfcwater(recp)                                                        !
      ! + csite%sfcwater_energy(k,recp) (Just converting back to J/kg)                     !
      ! + csite%csite%sfcwater_tempk(k,recp)                                               !
      ! + csite%sfcwater_fracliq(k,recp)                                                   !
      !------------------------------------------------------------------------------------!
      call new_patch_sfc_props(csite,recp,mzg,mzs,ntext_soil)
      !------------------------------------------------------------------------------------!

      csite%mean_rh(recp)                = newareai *                                      &
                                         ( csite%mean_rh(donp)        * csite%area(donp)   &
                                         + csite%mean_rh(recp)        * csite%area(recp) )

      csite%today_A_decomp(recp)         = newareai *                                      &
                                         ( csite%today_A_decomp(donp) * csite%area(donp)   &
                                         + csite%today_A_decomp(recp) * csite%area(recp) )

      csite%today_Af_decomp(recp)        = newareai *                                      &
                                         ( csite%today_Af_decomp(donp)* csite%area(donp)   &
                                         + csite%today_Af_decomp(recp)* csite%area(recp) )

      do iii = 1,n_pft
         csite%repro(iii,recp)           = newareai *                                      &
                                         ( csite%repro(iii,donp)      * csite%area(donp)   &
                                         + csite%repro(iii,recp)      * csite%area(recp) )
      end do
  
      !------------------------------------------------------------------------------------!
      !    Even though these variables are not prognostic, they need to be copied so the   !
      ! output will have the values.  Other variables will probably be scaled here as      !
      ! well.                                                                              !
      !------------------------------------------------------------------------------------!
      csite%avg_rshort_gnd(recp)      = newareai *                                         &
                                      ( csite%avg_rshort_gnd(donp)    * csite%area(donp)   &
                                      + csite%avg_rshort_gnd(recp)    * csite%area(recp) )

      csite%avg_rlong_gnd(recp)       = newareai *                                         &
                                      ( csite%avg_rlong_gnd(donp)     * csite%area(donp)   &
                                      + csite%avg_rlong_gnd(recp)     * csite%area(recp) )

      csite%avg_carbon_ac(recp)       = newareai *                                         &
                                      ( csite%avg_carbon_ac(donp)     * csite%area(donp)   &
                                      + csite%avg_carbon_ac(recp)     * csite%area(recp) )

      csite%avg_vapor_lc(recp)        = newareai *                                         &
                                      ( csite%avg_vapor_lc(donp)      * csite%area(donp)   &
                                      + csite%avg_vapor_lc(recp)      * csite%area(recp) )  

      csite%avg_vapor_wc(recp)        = newareai *                                         &
                                      ( csite%avg_vapor_wc(donp)      * csite%area(donp)   &
                                      + csite%avg_vapor_wc(recp)      * csite%area(recp) )  

      csite%avg_dew_cg(recp)          = newareai *                                         &
                                      ( csite%avg_dew_cg(donp)        * csite%area(donp)   &
                                      + csite%avg_dew_cg(recp)        * csite%area(recp) )  

      csite%avg_vapor_gc(recp)        = newareai *                                         &
                                      ( csite%avg_vapor_gc(donp)      * csite%area(donp)   &
                                      + csite%avg_vapor_gc(recp)      * csite%area(recp) )  

      csite%avg_wshed_vg(recp)        = newareai *                                         &
                                      ( csite%avg_wshed_vg(donp)      * csite%area(donp)   &
                                      + csite%avg_wshed_vg(recp)      * csite%area(recp) )  

      csite%avg_intercepted(recp)     = newareai *                                         &
                                      ( csite%avg_intercepted(donp)   * csite%area(donp)   &
                                      + csite%avg_intercepted(recp)   * csite%area(donp) )

      csite%avg_throughfall(recp)     = newareai *                                         &
                                      ( csite%avg_throughfall(donp)   * csite%area(donp)   &
                                      + csite%avg_throughfall(recp)   * csite%area(donp) )

      csite%avg_vapor_ac(recp)        = newareai *                                         &
                                      ( csite%avg_vapor_ac(donp)      * csite%area(donp)   &
                                      + csite%avg_vapor_ac(recp)      * csite%area(recp) )  

      csite%avg_transp(recp)          = newareai *                                         &
                                      ( csite%avg_transp(donp)        * csite%area(donp)   &
                                      + csite%avg_transp(recp)        * csite%area(recp) )  

      csite%avg_evap(recp)            = newareai *                                         &
                                      ( csite%avg_evap(donp)          * csite%area(donp)   &
                                      + csite%avg_evap(recp)          * csite%area(recp) )  

      csite%avg_runoff(recp)          = newareai *                                         &
                                      ( csite%avg_runoff(donp)        * csite%area(donp)   &
                                      + csite%avg_runoff(recp)        * csite%area(recp) )  

      csite%avg_drainage(recp)        = newareai *                                         &
                                      ( csite%avg_drainage(donp)      * csite%area(donp)   &
                                      + csite%avg_drainage(recp)      * csite%area(recp) )  

      csite%aux(recp)                 = newareai *                                         &
                                      ( csite%aux(donp)               * csite%area(donp)   &
                                      + csite%aux(recp)               * csite%area(recp) )  

      csite%avg_sensible_lc(recp)     = newareai *                                         &
                                      ( csite%avg_sensible_lc(donp)   * csite%area(donp)   &
                                      + csite%avg_sensible_lc(recp)   * csite%area(recp) )  

      csite%avg_sensible_wc(recp)     = newareai *                                         &
                                      ( csite%avg_sensible_wc(donp)   * csite%area(donp)   &
                                      + csite%avg_sensible_wc(recp)   * csite%area(recp) )  

      csite%avg_qwshed_vg(recp)       = newareai *                                         &
                                      ( csite%avg_qwshed_vg(donp)     * csite%area(donp)   &
                                      + csite%avg_qwshed_vg(recp)     * csite%area(recp) )  

      csite%avg_qintercepted(recp)    = newareai *                                         &
                                      ( csite%avg_qintercepted(donp)  * csite%area(donp)   &
                                      + csite%avg_qintercepted(recp)  * csite%area(donp) )

      csite%avg_qthroughfall(recp)    = newareai *                                         &
                                      ( csite%avg_qthroughfall(donp)  * csite%area(donp)   &
                                      + csite%avg_qthroughfall(recp)  * csite%area(donp) )

      csite%avg_sensible_gc(recp)     = newareai *                                         &
                                      ( csite%avg_sensible_gc(donp)   * csite%area(donp)   &
                                      + csite%avg_sensible_gc(recp)   * csite%area(recp) )  

      csite%avg_sensible_ac(recp)     = newareai *                                         &
                                      ( csite%avg_sensible_ac(donp)   * csite%area(donp)   &
                                      + csite%avg_sensible_ac(recp)   * csite%area(recp) )  

      csite%avg_runoff_heat(recp)     = newareai *                                         &
                                      ( csite%avg_runoff_heat(donp)   * csite%area(donp)   &
                                      + csite%avg_runoff_heat(recp)   * csite%area(recp) )  

      csite%avg_drainage_heat(recp)   = newareai *                                         &
                                      ( csite%avg_drainage_heat(donp) * csite%area(donp)   &
                                      + csite%avg_drainage_heat(recp) * csite%area(recp) )  

      csite%avg_leaf_energy(recp)     = newareai *                                         &
                                      ( csite%avg_leaf_energy(donp)   * csite%area(donp)   &
                                      + csite%avg_leaf_energy(recp)   * csite%area(recp) )

      csite%avg_leaf_water(recp)      = newareai *                                         &
                                      ( csite%avg_leaf_water(donp)    * csite%area(donp)   &
                                      + csite%avg_leaf_water(recp)    * csite%area(recp) )

      csite%avg_leaf_hcap(recp)       = newareai *                                         &
                                      ( csite%avg_leaf_hcap(donp)     * csite%area(donp)   &
                                      + csite%avg_leaf_hcap(recp)     * csite%area(recp) )

      csite%avg_wood_energy(recp)     = newareai *                                         &
                                      ( csite%avg_wood_energy(donp)   * csite%area(donp)   &
                                      + csite%avg_wood_energy(recp)   * csite%area(recp) )  

      csite%avg_wood_water(recp)      = newareai *                                         &
                                      ( csite%avg_wood_water(donp)    * csite%area(donp)   &
                                      + csite%avg_wood_water(recp)    * csite%area(recp) )  

      csite%avg_wood_hcap(recp)       = newareai *                                         &
                                      ( csite%avg_wood_hcap(donp)     * csite%area(donp)   &
                                      + csite%avg_wood_hcap(recp)     * csite%area(recp) )

      csite%co2budget_residual(recp)  = newareai *                                         &
                                      ( csite%co2budget_residual(donp)* csite%area(donp)   &
                                      + csite%co2budget_residual(recp)* csite%area(recp) )  

      csite%co2budget_loss2atm(recp)  = newareai *                                         &
                                      ( csite%co2budget_loss2atm(donp)* csite%area(donp)   &
                                      + csite%co2budget_loss2atm(recp)* csite%area(recp) )  

      csite%co2budget_denseffect(recp)= newareai *                                         &
                      ( csite%co2budget_denseffect(donp) * csite%area(donp)                &
                      + csite%co2budget_denseffect(recp) * csite%area(recp) )

      csite%co2budget_gpp(recp)       = newareai *                                         &
                                      ( csite%co2budget_gpp(donp)     * csite%area(donp)   &
                                      + csite%co2budget_gpp(recp)     * csite%area(recp) )  

      csite%co2budget_plresp(recp)    = newareai *                                         &
                                      ( csite%co2budget_plresp(donp)  * csite%area(donp)   &
                                      + csite%co2budget_plresp(recp)  * csite%area(recp) )  

      csite%co2budget_rh(recp)        = newareai *                                         &
                                      ( csite%co2budget_rh(donp)      * csite%area(donp)   &
                                      + csite%co2budget_rh(recp)      * csite%area(recp) )  

      csite%ebudget_residual(recp)    = newareai *                                         &
                                      ( csite%ebudget_residual(donp)  * csite%area(donp)   &
                                      + csite%ebudget_residual(recp)  * csite%area(recp) )

      csite%ebudget_loss2atm(recp)    = newareai *                                         &
                                      ( csite%ebudget_loss2atm(donp)  * csite%area(donp)   &
                                      + csite%ebudget_loss2atm(recp)  * csite%area(recp) )

      csite%ebudget_denseffect(recp)  = newareai *                                         &
                                      ( csite%ebudget_denseffect(donp) * csite%area(donp)  &
                                      + csite%ebudget_denseffect(recp) * csite%area(recp) )

      csite%ebudget_loss2runoff(recp) = newareai *                                         &
                                     ( csite%ebudget_loss2runoff(donp) * csite%area(donp)  &
                                     + csite%ebudget_loss2runoff(recp) * csite%area(recp) )

      csite%ebudget_loss2drainage(recp) = newareai *                                       &
                                   ( csite%ebudget_loss2drainage(donp) * csite%area(donp)  &
                                   + csite%ebudget_loss2drainage(recp) * csite%area(recp) )

      csite%ebudget_netrad(recp)      = newareai *                                         &
                                      ( csite%ebudget_netrad(donp) * csite%area(donp)      &
                                      + csite%ebudget_netrad(recp) * csite%area(recp) )

      csite%ebudget_precipgain(recp)  = newareai *                                         &
                                  ( csite%ebudget_precipgain(donp) * csite%area(donp)      &
                                  + csite%ebudget_precipgain(recp) * csite%area(recp) )

      csite%wbudget_residual(recp)    = newareai *                                         &
                                      ( csite%wbudget_residual(donp)  * csite%area(donp)   &
                                      + csite%wbudget_residual(recp)  * csite%area(recp) )

      csite%wbudget_loss2atm(recp)    = newareai *                                         &
                                      ( csite%wbudget_loss2atm(donp)  * csite%area(donp)   &
                                      + csite%wbudget_loss2atm(recp)  * csite%area(recp) )

      csite%wbudget_denseffect(recp)  = newareai *                                         &
                                      ( csite%wbudget_denseffect(donp) * csite%area(donp)  &
                                      + csite%wbudget_denseffect(recp) * csite%area(recp) )

      csite%wbudget_loss2runoff(recp) = newareai *                                         &
                                     ( csite%wbudget_loss2runoff(donp) * csite%area(donp)  &
                                     + csite%wbudget_loss2runoff(recp) * csite%area(recp) )

      csite%wbudget_loss2drainage(recp) = newareai *                                       &
                                   ( csite%wbudget_loss2drainage(donp) * csite%area(donp)  &
                                   + csite%wbudget_loss2drainage(recp) * csite%area(recp) )

      csite%wbudget_precipgain(recp)  = newareai *                                         &
                                  ( csite%wbudget_precipgain(donp) * csite%area(donp)      &
                                  + csite%wbudget_precipgain(recp) * csite%area(recp) )


      do iii=1,mzg
         csite%avg_smoist_gg(iii,recp)   = newareai *                                      &
              ( csite%avg_smoist_gg(iii,donp)       * csite%area(donp)                     &
              + csite%avg_smoist_gg(iii,recp)       * csite%area(recp) )

         csite%avg_transloss(iii,recp)   = newareai *                                      &
              ( csite%avg_transloss(iii,donp)       * csite%area(donp)                     &
              + csite%avg_transloss(iii,recp)       * csite%area(recp) )

         csite%aux_s(iii,recp)           = newareai *                                      &
              ( csite%aux_s(iii,donp)               * csite%area(donp)                     &
              + csite%aux_s(iii,recp)               * csite%area(recp) )

         csite%avg_sensible_gg(iii,recp) = newareai *                                      &
              ( csite%avg_sensible_gg(iii,donp)     * csite%area(donp)                     &
              + csite%avg_sensible_gg(iii,recp)     * csite%area(recp) )
      end do

      do iii=1,n_dbh
         csite%co2budget_gpp_dbh(iii,recp) = newareai *                                    &
              ( csite%co2budget_gpp_dbh(iii,donp)   * csite%area(donp)                     &
              + csite%co2budget_gpp_dbh(iii,recp)   * csite%area(recp) )
      end do
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------! 
      !    We must also check whether daily/monthly output variables exist.  If they do,   !
      ! then we must fuse them too.                                                        !
      !------------------------------------------------------------------------------------! 
      if (idoutput > 0 .or. imoutput > 0 .or. iqoutput > 0) then
         csite%dmean_rh(recp)           = newareai                                         &
                                        * ( csite%dmean_rh(donp) * csite%area(donp)        &
                                          + csite%dmean_rh(recp) * csite%area(recp) )
         csite%dmean_co2_residual(recp) = newareai                                         &
                                        * ( csite%dmean_co2_residual(donp)                 &
                                          * csite%area(donp)                               &
                                          + csite%dmean_co2_residual(recp)                 &
                                          * csite%area(recp) )
         csite%dmean_energy_residual(recp) = newareai                                      &
                                           * ( csite%dmean_energy_residual(donp)           &
                                             * csite%area(donp)                            &
                                             + csite%dmean_energy_residual(recp)           &
                                             * csite%area(recp) )
         csite%dmean_water_residual(recp)  = newareai                                      &
                                           * ( csite%dmean_water_residual(donp)            &
                                             * csite%area(donp)                            &
                                             + csite%dmean_water_residual(recp)            &
                                             * csite%area(recp) )
         csite%dmean_lambda_light(recp)    = newareai                                      &
                                           * ( csite%dmean_lambda_light(donp)              &
                                             * csite%area(donp)                            &
                                             + csite%dmean_lambda_light(recp)              &
                                             * csite%area(recp) )
         csite%dmean_A_decomp(recp)        = newareai                                      &
                                           * ( csite%dmean_A_decomp(donp)                  &
                                             * csite%area(donp)                            &
                                             + csite%dmean_A_decomp(recp)                  &
                                             * csite%area(recp) )

         csite%dmean_Af_decomp(recp)       = newareai                                      &
                                           * ( csite%dmean_Af_decomp(donp)                 &
                                             * csite%area(donp)                            &
                                             + csite%dmean_Af_decomp(recp)                 &
                                             * csite%area(recp) )
         csite%dmean_rk4step(recp)         = newareai                                      &
                                           * ( csite%dmean_rk4step(donp)                   &
                                             * csite%area(donp)                            &
                                             + csite%dmean_rk4step(recp)                   &
                                             * csite%area(recp) )
      end if
      if (imoutput > 0 .or. iqoutput > 0) then
         csite%mmean_rh(recp)           = newareai                                         &
                                        * ( csite%mmean_rh(donp) * csite%area(donp)        &
                                          + csite%mmean_rh(recp) * csite%area(recp) )
         csite%mmean_co2_residual(recp) = newareai                                         &
                                        * ( csite%mmean_co2_residual(donp)                 &
                                          * csite%area(donp)                               &
                                          + csite%mmean_co2_residual(recp)                 &
                                          * csite%area(recp) )
         csite%mmean_energy_residual(recp) = newareai                                      &
                                           * ( csite%mmean_energy_residual(donp)           &
                                             * csite%area(donp)                            &
                                             + csite%mmean_energy_residual(recp)           &
                                             * csite%area(recp) )
         csite%mmean_water_residual(recp)  = newareai                                      &
                                           * ( csite%mmean_water_residual(donp)            &
                                             * csite%area(donp)                            &
                                             + csite%mmean_water_residual(recp)            &
                                             * csite%area(recp) )
         csite%mmean_lambda_light(recp)    = newareai                                      &
                                           * ( csite%mmean_lambda_light(donp)              &
                                             * csite%area(donp)                            &
                                             + csite%mmean_lambda_light(recp)              &
                                             * csite%area(recp) )
         csite%mmean_A_decomp(recp)        = newareai                                      &
                                           * ( csite%mmean_A_decomp(donp)                  &
                                             * csite%area(donp)                            &
                                             + csite%mmean_A_decomp(recp)                  &
                                             * csite%area(recp) )
         csite%mmean_Af_decomp(recp)       = newareai                                      &
                                           * ( csite%mmean_Af_decomp(donp)                 &
                                             * csite%area(donp)                            &
                                             + csite%mmean_Af_decomp(recp)                 &
                                             * csite%area(recp) )
         csite%mmean_rk4step(recp)         = newareai                                      &
                                           * ( csite%mmean_rk4step(donp)                   &
                                             * csite%area(donp)                            &
                                             + csite%mmean_rk4step(recp)                   &
                                             * csite%area(recp) )
      end if

      if (iqoutput > 0) then
         do icyc=1,ndcycle
            csite%qmean_rh     (icyc,recp) = newareai                                      &
                                           * ( csite%qmean_rh               (icyc,donp)    &
                                             * csite%area                        (donp)    &
                                             + csite%qmean_rh               (icyc,recp)    &
                                             * csite%area                        (recp))
         end do
      end if

      !------------------------------------------------------------------------------------!
      !      Update the leaf and wood temperature and liquid fraction.  We must check      !
      ! whether we can solve the average or not, because it may be an empty patch or the   !
      ! user may have disabled branchwood thermodynamics.                                  !
      !------------------------------------------------------------------------------------!
      if (csite%avg_leaf_hcap(recp) > 0.) then
         call qwtk(csite%avg_leaf_energy(recp),csite%avg_leaf_water(recp)                  &
                  ,csite%avg_leaf_hcap(recp),csite%avg_leaf_temp(recp)                     &
                  ,csite%avg_leaf_fliq(recp))
      else
         csite%avg_leaf_temp(recp) = newareai                                              &
                                   * ( csite%avg_leaf_temp(donp) * csite%area(donp)        &
                                     + csite%avg_leaf_temp(recp) * csite%area(recp) )
         csite%avg_leaf_fliq(recp) = newareai                                              &
                                   * ( csite%avg_leaf_fliq(donp) * csite%area(donp)        &
                                     + csite%avg_leaf_fliq(recp) * csite%area(recp) )
      end if
      if (csite%avg_wood_hcap(recp) > 0.) then
         call qwtk(csite%avg_wood_energy(recp),csite%avg_wood_water(recp)                  &
                  ,csite%avg_wood_hcap(recp),csite%avg_wood_temp(recp)                     &
                  ,csite%avg_wood_fliq(recp))
      else
         csite%avg_wood_temp(recp) = newareai                                              &
                                   * ( csite%avg_wood_temp(donp) * csite%area(donp)        &
                                     + csite%avg_wood_temp(recp) * csite%area(recp) )
         csite%avg_wood_fliq(recp) = newareai                                              &
                                   * ( csite%avg_wood_fliq(donp) * csite%area(donp)        &
                                     + csite%avg_wood_fliq(recp) * csite%area(recp) )
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    We now update the canopy thermodynamic propeties:                               !
      ! + csite%can_temp(recp)                                                             !
      ! + csite%can_rhos(recp)                                                             !
      !------------------------------------------------------------------------------------!
      call update_patch_thermo_props(csite,recp,recp,mzg,mzs,ntext_soil)

      !------------------------------------------------------------------------------------!
      !     Now we need to adjust the densities of cohorts. Because the patch area         !
      ! increased we want to retain the same total amount of mass and energy.              !
      !------------------------------------------------------------------------------------!
      !----- 1. Adjust densities of cohorts in recipient patch ----------------------------!
      cpatch => csite%patch(recp)
      nrc = cpatch%ncohorts
      area_scale = csite%area(recp) * newareai
      !------------------------------------------------------------------------------------!
      ! IMPORTANT: Only cohort-level variables that have units per area (m2) should be     !
      !            rescaled.  Variables whose units are per plant should _NOT_ be included !
      !            here.                                                                   !
      !------------------------------------------------------------------------------------!
      call fuse_2_patches_costate(nrc, cpatch%costate, area_scale)
      call fuse_2_patches_comort(nrc, cpatch%comort, area_scale)
      call fuse_2_patches_coresp(nrc, cpatch%coresp, area_scale)
      call fuse_2_patches_cophoto(nrc, cpatch%cophoto, area_scale)
      call fuse_2_patches_cotherm(nrc, cpatch%cotherm, area_scale)

      !----- 2. Adjust densities of cohorts in donor patch --------------------------------!
      cpatch => csite%patch(donp)
      ndc = cpatch%ncohorts
      area_scale = csite%area(donp) * newareai
      !------------------------------------------------------------------------------------!
      ! IMPORTANT: Only cohort-level variables that have units per area (m2) should be     !
      !            rescaled.  Variables whose units are per plant should _NOT_ be included !
      !            here.                                                                   !
      !------------------------------------------------------------------------------------!

      call fuse_2_patches_costate(nrc, cpatch%costate, area_scale)
      call fuse_2_patches_comort(nrc, cpatch%comort, area_scale)
      call fuse_2_patches_coresp(nrc, cpatch%coresp, area_scale)
      call fuse_2_patches_cophoto(nrc, cpatch%cophoto, area_scale)
      call fuse_2_patches_cotherm(nrc, cpatch%cotherm, area_scale)

      !------------------------------------------------------------------------------------!
      !    Fill a new patch with the donor and recipient cohort vectors.                   !
      !------------------------------------------------------------------------------------!
      !----- Allocate the temporary patch with room for all cohorts. ----------------------!
      if (ndc + nrc > 0 ) then
         nullify(temppatch)
         allocate(temppatch)
         call allocate_patchtype(temppatch,ndc + nrc )
         !----- Copy the recipient and donor cohorts to the temporary patch. --------------!
         call copy_patchtype(csite%patch(recp),temppatch,1,nrc,1,nrc)
         call copy_patchtype(csite%patch(donp),temppatch,1,ndc,nrc+1,nrc+ndc)
         !----- Reallocate the current recipient patch with all cohorts -------------------!
         call deallocate_patchtype(csite%patch(recp))
         call allocate_patchtype(csite%patch(recp),ndc+nrc)
         !----- Copy the temporary patch back to the recipient patch. ---------------------!
         call copy_patchtype(temppatch,csite%patch(recp),1,nrc+ndc,1,nrc+ndc)
         !----- Get rid of the temporary patch --------------------------------------------!
         call deallocate_patchtype(temppatch)
         deallocate(temppatch)
         !----- Sort cohorts in the new patch ---------------------------------------------!
         cpatch => csite%patch(recp)
         call sort_cohorts(cpatch)
         !---------------------------------------------------------------------------------!
         !    We just combined two patches, so we may be able to fuse some cohorts and/or  !
         ! eliminate others.                                                               !
         !---------------------------------------------------------------------------------!
         if (cpatch%ncohorts > 0 .and. maxcohort >= 0) then
            call fuse_cohorts(csite,recp,green_leaf_factor,lsl)
            call terminate_cohorts(csite,recp,elim_nplant,elim_lai)
            call split_cohorts(cpatch,green_leaf_factor,lsl)
         end if
         !---------------------------------------------------------------------------------!
      end if

      !------------------------------------------------------------------------------------!
      !    Now we update some variables that depend on cohort statistics, namely:          !
      ! + csite%veg_height(recp)                                                           !
      ! + csite%veg_displace(recp)                                                         !
      ! + csite%disp_height(recp)                                                          !
      ! + csite%lai(recp)                                                                  !
      ! + csite%veg_rough(recp)                                                            !
      ! + csite%total_sfcw_depth(recp)                                                     !
      ! + csite%snowfac(recp)                                                              !
      ! + csite%opencan_frac(recp)                                                         !
      !------------------------------------------------------------------------------------!
      call update_patch_derived_props(csite,lsl, prss,recp)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Now we update the budget variables:                                             !
      ! + csite%wbudget_initialstorage(recp)                                               !
      ! + csite%ebudget_initialstorage(recp)                                               !
      ! + csite%co2budget_initialstorage(recp)                                             !
      !------------------------------------------------------------------------------------!
      call update_budget(csite,lsl,recp,recp)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    This subroutine will update the size profile within patch.                      !
      ! + csite%cumlai_profile(:,:,recp)                                                   !
      !------------------------------------------------------------------------------------!
      call patch_pft_size_profile(csite,recp)
      !------------------------------------------------------------------------------------!

      !----- Last, but not the least, we update the patch area ----------------------------!
      csite%area(recp) = newarea

      return

   end subroutine fuse_2_patches
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine patch_pft_size_profile(csite,ipa)
      use ed_state_vars       , only : sitetype   & ! structure
                                     , patchtype  ! ! structure
      use fusion_fission_coms , only : ff_nhgt    & ! intent(in)
                                     , hgt_class  ! ! intent(in)
      use allometry           , only : dbh2bl     ! ! intent(in)
      use ed_max_dims         , only : n_pft      ! ! intent(in)
      use pft_coms            , only : hgt_min    ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target     :: csite     ! Current site
      integer                , intent(in) :: ipa       ! Current patch index
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)        , pointer    :: cpatch    ! Current patch
      integer                             :: ipft      ! PFT index
      integer                             :: ihgt      ! Height class index
      integer                             :: ico       ! Counters
      real                                :: lai_pot   ! Potential LAI
      !------------------------------------------------------------------------------------!


      !----- Reset all bins to zero. ------------------------------------------------------!
      do ipft=1,n_pft
         do ihgt=1,ff_nhgt
            csite%cumlai_profile(ipft,ihgt,ipa)=0.0
         end do
      end do
      !------------------------------------------------------------------------------------!



      !----- Update bins ------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      cohortloop: do ico = 1,cpatch%ncohorts

         !----- Find the PFT class. -------------------------------------------------------!
         ipft    = cpatch%costate%pft(ico)
         ihgt    = min(ff_nhgt,1 + count(hgt_class < cpatch%costate%hite(ico)))
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Check whether this cohort is almost at the minimum height given its PFT.    !
         ! If it is, then we will skip it.                                                 !
         !---------------------------------------------------------------------------------!
         if (cpatch%costate%hite(ico) < hgt_min(ipft) + 0.2) cycle cohortloop
         !---------------------------------------------------------------------------------!


         !----- Find the height class. ----------------------------------------------------!
         ihgt    = min(ff_nhgt,1 + count(hgt_class < cpatch%costate%hite(ico)))
         !---------------------------------------------------------------------------------!


         !----- Find the potential (on-allometry) leaf area index. ------------------------!
         lai_pot = cpatch%costate%nplant(ico) * cpatch%cophen%sla(ico)                                    &
                 * dbh2bl(cpatch%costate%dbh(ico),ipft)
         !---------------------------------------------------------------------------------!


         !----- Add the potential LAI to the bin. -----------------------------------------!
         csite%cumlai_profile(ipft,ihgt,ipa) = lai_pot                                     &
                                             + csite%cumlai_profile(ipft,ihgt,ipa)
         !---------------------------------------------------------------------------------!
      end do cohortloop
      !------------------------------------------------------------------------------------!



      !----- Integrate the leaf area index from top to bottom. ----------------------------!
      do ihgt=ff_nhgt-1,1,-1
         do ipft=1,n_pft
            csite%cumlai_profile(ipft,ihgt,ipa) = csite%cumlai_profile(ipft,ihgt  ,ipa)    &
                                                + csite%cumlai_profile(ipft,ihgt+1,ipa)
         end do
      end do
      !------------------------------------------------------------------------------------!

      return
   end subroutine patch_pft_size_profile
   !=======================================================================================!
   !=======================================================================================!
end module fuse_fiss_utils
!==========================================================================================!
!==========================================================================================!
