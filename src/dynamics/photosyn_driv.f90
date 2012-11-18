!==========================================================================================!
!==========================================================================================!
!     This subroutine will control the photosynthesis scheme (Farquar and Leuning).  This  !
! is called every step, but not every sub-step.                                            !
!------------------------------------------------------------------------------------------!
subroutine canopy_photosynthesis(csite,cmet,mzg,ipa,lsl,ntext_soil                         &
                                ,leaf_aging_factor,green_leaf_factor)
   use ed_state_vars  , only : sitetype          & ! structure
                             , patchtype         ! ! structure
   use ed_max_dims    , only : n_pft             ! ! intent(in)
   use pft_coms       , only : leaf_width        & ! intent(in)
                             , water_conductance & ! intent(in)
                             , include_pft       ! ! intent(in)
   use soil_coms      , only : soil              & ! intent(in)
                             , slz               & ! intent(in)
                             , dslz              & ! intent(in)
                             , freezecoef        ! ! intent(in)
   use consts_coms    , only : t00               & ! intent(in)
                             , epi               & ! intent(in)
                             , wdnsi             & ! intent(in)
                             , wdns              & ! intent(in)
                             , kgCday_2_umols, umol_2_kgC    & ! intent(in)
                             , lnexp_min         ! ! intent(in)
   use ed_misc_coms   , only : current_time, dtlsm      ! ! intent(in)
   use met_driver_coms, only : met_driv_state    ! ! structure
   use physiology_coms, only : print_photo_debug & ! intent(in)
                             , h2o_plant_lim     ! ! intent(in)
   use farq_leuning   , only : lphysiol_full     ! ! sub-routine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)            , target      :: csite             ! Current site
   type(met_driv_state)      , target      :: cmet              ! Current met. conditions.
   integer                   , intent(in)  :: ipa               ! Current patch #
   integer                   , intent(in)  :: lsl               ! Lowest soil level
   integer                   , intent(in)  :: mzg               ! Number of soil layers
   integer, dimension(mzg)   , intent(in)  :: ntext_soil        ! Soil class
   real   , dimension(n_pft) , intent(in)  :: leaf_aging_factor ! 
   real   , dimension(n_pft) , intent(in)  :: green_leaf_factor ! 
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)           , pointer     :: cpatch             ! Current site
   integer                                 :: ico                ! Current cohort #
   integer                                 :: tuco               ! Tallest used cohort
   integer                                 :: ipft
   integer                                 :: k1
   integer                                 :: k2
   integer                                 :: kroot
   integer                                 :: nsoil
   integer                                 :: limit_flag
   logical, dimension(mzg+1)               :: root_depth_indices ! 
   logical                                 :: las
   real   , dimension(mzg+1)               :: avg_frozen_water
   real   , dimension(mzg+1)               :: avg_liquid_water
   real   , dimension(mzg+1)               :: available_liquid_water
   real   , dimension(mzg+1)               :: wilting_factor
   real                                    :: leaf_par
   real                                    :: leaf_resp
   real                                    :: d_gsw_open
   real                                    :: d_gsw_closed
   real                                    :: slpotv
   real                                    :: swp
   real                                    :: vm
   real                                    :: compp
   real                                    :: broot_tot
   real                                    :: broot_loc
   real                                    :: wgpfrac
   real                                    :: avg_fracliq
   real                                    :: psi_wilting
   real                                    :: psi_layer
   real                                    :: freezecor
   real                                    :: pss_available_water
   !---------------------------------------------------------------------------------------!


   !----- Point to the cohort structures --------------------------------------------------!
   cpatch => csite%patch(ipa)

   !----- Find the patch-level Total Leaf and Wood Area Index. ----------------------------!
   csite%lai(ipa) = 0.0
   csite%wpa(ipa) = 0.0
   csite%wai(ipa) = 0.0
   do ico=1,cpatch%ncohorts
      csite%lai(ipa)  = csite%lai(ipa)  + cpatch%costate%lai(ico)
      csite%wpa(ipa)  = csite%wpa(ipa)  + cpatch%costate%wpa(ico)
      csite%wai(ipa)  = csite%wai(ipa)  + cpatch%costate%wai(ico)
   end do


   !----- Calculate liquid water available for transpiration. -----------------------------!
   available_liquid_water(mzg+1) = 0.
   do k1 = mzg, lsl, -1
      nsoil = ntext_soil(k1)
      available_liquid_water(k1) = available_liquid_water(k1+1)                            &
                                 + wdns * dslz(k1) * csite%soil_fracliq(k1,ipa)            &
                                 * max(0.0, csite%soil_water(k1,ipa) - soil(nsoil)%soilwp )
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     If we are solving H2O_PLANT_LIM = 2, then we must account for the water potential !
   ! as in CLM.                                                                            !
   !---------------------------------------------------------------------------------------!
   select case (h2o_plant_lim)
   case (2)

      !----- Find the average liquid and frozen water. ------------------------------------!
      avg_liquid_water(:) = 0.0
      avg_frozen_water(:) = 0.0
      do k1=mzg,lsl,-1
         avg_frozen_water(k1) = avg_frozen_water(k1+1)                                     &
                              + dslz(k1) * (1. - csite%soil_fracliq(k1,ipa))               &
                              * csite%soil_water(k1,ipa)
         avg_liquid_water(k1) = avg_liquid_water(k1+1)                                     &
                              + dslz(k1) * csite%soil_fracliq(k1,ipa)                      &
                              * csite%soil_water(k1,ipa)
      end do
      avg_frozen_water(lsl:mzg) = avg_frozen_water(lsl:mzg) / (- slz(lsl:mzg))
      avg_liquid_water(lsl:mzg) = avg_liquid_water(lsl:mzg) / (- slz(lsl:mzg))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Compute the soil potential for transpiration for each layer as in CLM.          !
      !------------------------------------------------------------------------------------!
      wilting_factor(:) = 0.0
      do k1=lsl,mzg
         nsoil = ntext_soil(k1)
         if (avg_liquid_water(k1) > soil(nsoil)%soilwp) then
            !----- Find the average soil wetness. -----------------------------------------!
            wgpfrac = avg_liquid_water(k1) / soil(nsoil)%slmsts
            avg_fracliq  = avg_liquid_water(k1)                                            &
                         / (avg_liquid_water(k1) + avg_frozen_water(k1))

            !----- Apply correction in case the soil is partially frozen. -----------------!
            freezecor    = 10. ** max(lnexp_min,- freezecoef * (1.0 - avg_fracliq))
            psi_wilting  = soil(nsoil)%slpots                                              &
                         / (soil(nsoil)%soilwp / soil(nsoil)%slmsts) ** soil(nsoil)%slbs
            psi_layer    = soil(nsoil)%slpots / wgpfrac ** soil(nsoil)%slbs

            !----- Find the wilting factor which will control the dry soil correction. ----!
            wilting_factor(k1) = max( 0., min(1., (psi_wilting - psi_layer         )       &
                                                / (psi_wilting - soil(nsoil)%slpots) ) )
         end if
      end do
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Initialize the array of maximum photosynthesis rates used in the mortality        !
   ! function.                                                                             !
   !---------------------------------------------------------------------------------------!
   csite%A_o_max(1:n_pft,ipa) = 0.0
   csite%A_c_max(1:n_pft,ipa) = 0.0

   !---------------------------------------------------------------------------------------!
   !     Find the tallest cohort with TAI above minimum, sufficient heat capacity, and not !
   ! buried in snow.  The first two conditions are redundant, but we will keep them for    !
   ! the time being, so it is going to be safer.                                           !
   !---------------------------------------------------------------------------------------!
   las = .false.
   do ico = 1,cpatch%ncohorts
      !----- If this is the tallest cohort to be used, we save its index. -----------------!
      if (.not. las .and. cpatch%costate%leaf_resolvable(ico)) then
         las  = .true.
         tuco = ico
      end if
   end do

   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !    There is at least one cohort that meet requirements.  And this is tallest one, so  !
   ! we can use it to compute the maximum photosynthetic rates, i.e., the rate the cohort  !
   ! would have if it were at the top of the canopy.  This is used for the mortality       !
   ! function.                                                                             !
   !---------------------------------------------------------------------------------------!
   if (las) then
      !----- We now loop over PFTs, not cohorts, skipping those we are not using. ---------!
      do ipft = 1, n_pft
         if (include_pft(ipft)) then

            !------------------------------------------------------------------------------!
            !    Scale photosynthetically active radiation per unit of leaf.               !
            !------------------------------------------------------------------------------!
            leaf_par = csite%par_l_max(ipa) / cpatch%costate%lai(tuco)

            !------------------------------------------------------------------------------!
            !    Call the photosynthesis for maximum photosynthetic rates.  The units      !
            ! of the input and output are the standard in most of ED modules, but many of  !
            ! them are converted inside the photosynthesis model.                          !
            !    Notice that the units that are per unit area are per m² of leaf, not the  !
            ! patch area.                                                                  !
            !------------------------------------------------------------------------------!
            call lphysiol_full(            & !
               csite%can_prss(ipa)         & ! Canopy air pressure              [       Pa]
             , csite%can_rhos(ipa)         & ! Canopy air density               [    kg/m³]
             , csite%can_shv(ipa)          & ! Canopy air sp. humidity          [    kg/kg]
             , csite%can_co2(ipa)          & ! Canopy air CO2 mixing ratio      [ µmol/mol]
             , ipft                        & ! Plant functional type            [      ---]
             , leaf_par                    & ! Absorbed photos. active rad.     [     W/m²]
             , cpatch%cotherm%leaf_temp(tuco)      & ! Leaf temperature                 [        K]
             , green_leaf_factor(ipft)     & ! Greenness rel. to on-allometry   [      ---]
             , leaf_aging_factor(ipft)     & ! Ageing parameter to scale VM     [      ---]
             , cpatch%cophen%llspan(tuco)         & ! Leaf life span                   [       yr]
             , cpatch%cophen%vm_bar(tuco)         & ! Average Vm function              [µmol/m²/s]
             , cpatch%cotherm%leaf_gbw(tuco)       & ! Aerodyn. condct. of water vapour [  kg/m²/s]
             , csite%A_o_max(ipft,ipa)     & ! Photosynthesis rate     (open)   [µmol/m²/s]
             , csite%A_c_max(ipft,ipa)     & ! Photosynthesis rate     (closed) [µmol/m²/s]
             , d_gsw_open                  & ! Stom. condct. of water  (open)   [  kg/m²/s]
             , d_gsw_closed                & ! Stom. condct. of water  (closed) [  kg/m²/s]
             , leaf_resp                   & ! Leaf respiration rate            [µmol/m²/s]
             , vm                          & ! Max. capacity of Rubisco         [µmol/m²/s]
             , compp                       & ! Gross photo. compensation point  [ µmol/mol]
             , limit_flag                  & ! Photosynthesis limitation flag   [      ---]
             )
         end if
      end do
         
   else
      !---- There is no "active" cohort. --------------------------------------------------!
      csite%A_o_max(1:n_pft,ipa) = 0.0
      csite%A_c_max(1:n_pft,ipa) = 0.0
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Initialize some variables.                                                         !
   !---------------------------------------------------------------------------------------!
   !----- Total root biomass (in kgC/m2) and patch sum available water. -------------------!
   pss_available_water = 0.0
   broot_tot           = 0.0
   !----- Initialize variables for transpiration calculation. -----------------------------!
   root_depth_indices(:) = .false.
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Loop over all cohorts, from tallest to shortest.                                   !
   !---------------------------------------------------------------------------------------!
   cohortloop: do ico = 1,cpatch%ncohorts
         
      !------------------------------------------------------------------------------------!
      !     Only need to worry about photosyn if radiative transfer has been  done for     !
      ! this cohort.                                                                       !
      !------------------------------------------------------------------------------------!
      if (cpatch%costate%leaf_resolvable(ico)) then

            !----- Alias for PFT ----------------------------------------------------------!
            ipft = cpatch%costate%pft(ico)

            !------------------------------------------------------------------------------!
            !    Scale photosynthetically active radiation per unit of leaf.               !
            !------------------------------------------------------------------------------!
            leaf_par = cpatch%corad%par_l(ico) / cpatch%costate%lai(ico) 


            !------------------------------------------------------------------------------!
            !    Call the photosynthesis for actual photosynthetic rates.  The units       !
            ! of the input and output are the standard in most of ED modules, but many of  !
            ! them are converted inside the photosynthesis model.                          !
            !    Notice that the units that are per unit area are per m² of leaf, not the  !
            ! patch area.                                                                  !
            !------------------------------------------------------------------------------!
            call lphysiol_full(            & !
               csite%can_prss(ipa)         & ! Canopy air pressure              [       Pa]
             , csite%can_rhos(ipa)         & ! Canopy air density               [    kg/m³]
             , csite%can_shv(ipa)          & ! Canopy air sp. humidity          [    kg/kg]
             , csite%can_co2(ipa)          & ! Canopy air CO2 mixing ratio      [ µmol/mol]
             , ipft                        & ! Plant functional type            [      ---]
             , leaf_par                    & ! Absorbed photos. active rad.     [     W/m²]
             , cpatch%cotherm%leaf_temp(ico)       & ! Leaf temperature                 [        K]
             , green_leaf_factor(ipft)     & ! Greenness rel. to on-allometry   [      ---]
             , leaf_aging_factor(ipft)     & ! Ageing parameter to scale VM     [      ---]
             , cpatch%cophen%llspan(ico)          & ! Leaf life span                   [       yr]
             , cpatch%cophen%vm_bar(ico)          & ! Average Vm function              [µmol/m²/s]
             , cpatch%cotherm%leaf_gbw(ico)        & ! Aerodyn. condct. of water vapour [  kg/m²/s]
             , cpatch%cophoto%A_open(ico)          & ! Photosynthesis rate     (open)   [µmol/m²/s]
             , cpatch%cophoto%A_closed(ico)        & ! Photosynthesis rate     (closed) [µmol/m²/s]
             , cpatch%cophoto%gsw_open(ico)        & ! Stom. condct. of water  (open)   [  kg/m²/s]
             , cpatch%cophoto%gsw_closed(ico)      & ! Stom. condct. of water  (closed) [  kg/m²/s]
             , leaf_resp                   & ! Leaf respiration rate            [µmol/m²/s]
             , vm                          & ! Max. capacity of Rubisco         [µmol/m²/s]
             , compp                       & ! Gross photo. compensation point  [ µmol/mol]
             , limit_flag                  & ! Photosynthesis limitation flag   [      ---]
             )

            !----- Convert leaf respiration to [µmol/m²ground/s] --------------------------!
            cpatch%coresp%leaf_respiration(ico) = leaf_resp * cpatch%costate%lai(ico)
            cpatch%coresp%mean_leaf_resp(ico)   = cpatch%coresp%mean_leaf_resp(ico)                      &
                                         + cpatch%coresp%leaf_respiration(ico)
            cpatch%coresp%dmean_leaf_resp(ico)  = cpatch%coresp%dmean_leaf_resp(ico)                     &
                                         + cpatch%coresp%leaf_respiration(ico) / cpatch%costate%nplant(ico) * umol_2_kgC * dtlsm

            !----- Root biomass [kg/m2]. --------------------------------------------------!
            broot_loc = cpatch%costate%broot(ico)  * cpatch%costate%nplant(ico)

            !----- Supply of water. -------------------------------------------------------!
            cpatch%cophoto%water_supply(ico) = water_conductance(ipft)                             &
                                     * available_liquid_water(cpatch%costate%krdepth(ico))         &
                                     * broot_loc

            root_depth_indices(cpatch%costate%krdepth(ico)) = .true.
            broot_tot = broot_tot + broot_loc
            pss_available_water = pss_available_water                                      &
                                + available_liquid_water(cpatch%costate%krdepth(ico)) * broot_loc

            !------------------------------------------------------------------------------!
            !     Determine the fraction of open stomata due to water limitation.          !
            ! This is a function of the ratio between the potential water demand           !
            ! (cpatch%psi_open, which is the average over the last time step), and the     !
            ! supply (cpatch%water_supply).                                                !
            !------------------------------------------------------------------------------!
            select case (h2o_plant_lim)
            case (0)
               !---- No water limitation, fsw is always 1.0. ------------------------------!
               cpatch%cophoto%fsw(ico) = 1.0

            case (1)
               !---- Original ED-1.0 scheme. ----------------------------------------------!
               cpatch%cophoto%fsw(ico) = cpatch%cophoto%water_supply(ico)                                  &
                               / max( 1.0e-20                                              &
                                    , cpatch%cophoto%water_supply(ico) + cpatch%cophoto%psi_open(ico))
            case (2)
               !---------------------------------------------------------------------------!
               !     Somewhat based on CLM, but we reduce the total amount of available    !
               ! water by the fraction of root biomass belonging to this cohort.  We don't !
               ! have the root profile up to now, assume they are evenly distributed       !
               ! through all layers that have roots.                                       !
               !---------------------------------------------------------------------------!
               cpatch%cophoto%fsw(ico) = wilting_factor(cpatch%costate%krdepth(ico))

            end select
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Photorespiration can become important at high temperatures.  If so,     !
            ! close down the stomata.                                                      !
            !------------------------------------------------------------------------------!
            if (cpatch%cophoto%A_open(ico) < cpatch%cophoto%A_closed(ico)) then
               cpatch%cophoto%fs_open(ico) = 0.0
            else
               cpatch%cophoto%fs_open(ico) = cpatch%cophoto%fsw(ico) * cpatch%cophoto%fsn(ico)
            end if

            !----- Net stomatal conductance. ----------------------------------------------!

            cpatch%cophoto%stomatal_conductance(ico) =  cpatch%cophoto%fs_open(ico) *cpatch%cophoto%gsw_open(ico)  &
                                             + (1.0 - cpatch%cophoto%fs_open(ico))                 &
                                             * cpatch%cophoto%gsw_closed(ico)

            !----- GPP, averaged over frqstate. -------------------------------------------!
            cpatch%cophoto%gpp(ico)       = cpatch%costate%lai(ico)                                        &
                                  * ( cpatch%cophoto%fs_open(ico) * cpatch%cophoto%A_open(ico)             &
                                    + (1.0 - cpatch%cophoto%fs_open(ico)) * cpatch%cophoto%A_closed(ico) ) &
                                  + cpatch%coresp%leaf_respiration(ico)
            cpatch%cophoto%mean_gpp(ico)  = cpatch%cophoto%mean_gpp(ico) + cpatch%cophoto%gpp(ico)

            !----- GPP, summed over 1 day. [kgC/plant/day] --------------------------------!
            cpatch%cophoto%today_gpp(ico) = cpatch%cophoto%today_gpp(ico) + cpatch%cophoto%gpp(ico) / cpatch%costate%nplant(ico) * umol_2_kgC * dtlsm

            !----- Potential GPP if no N limitation. [µmol/m²ground] ----------------------!
!            cpatch%cophoto%today_gpp_pot(ico) = cpatch%cophoto%today_gpp_pot(ico)                          &
!                                      + cpatch%lai(ico)                                    &
!                                      * ( cpatch%fsw(ico) * cpatch%A_open(ico)             &
!                                        + (1.0 - cpatch%fsw(ico)) * cpatch%A_closed(ico))  &
!                                      + cpatch%coresp%leaf_respiration(ico)
            cpatch%cophoto%today_gpp_pot(ico) = cpatch%cophoto%today_gpp_pot(ico)                          &
                                      + (cpatch%costate%lai(ico)                                    &
                                      * ( cpatch%cophoto%fsw(ico) * cpatch%cophoto%A_open(ico)             &
                                        + (1.0 - cpatch%cophoto%fsw(ico)) * cpatch%cophoto%A_closed(ico))  &
                                      + cpatch%coresp%leaf_respiration(ico)) / cpatch%costate%nplant(ico) * umol_2_kgC * dtlsm

            !----- Maximum GPP if at the top of the canopy [µmol/m²ground] ----------------!
!            cpatch%cophoto%today_gpp_max(ico) = cpatch%cophoto%today_gpp_max(ico)                          &
!                                      + cpatch%lai(ico)                                    &
!                                      * ( cpatch%fs_open(ico) * csite%A_o_max(ipft,ipa)    &
!                                        + (1.0 - cpatch%fs_open(ico))                      &
!                                          * csite%A_c_max(ipft,ipa))                       &
!                                      + cpatch%coresp%leaf_respiration(ico)
            cpatch%cophoto%today_gpp_max(ico) = cpatch%cophoto%today_gpp_max(ico)                          &
                                      + (cpatch%costate%lai(ico)                                    &
                                      * ( cpatch%cophoto%fs_open(ico) * csite%A_o_max(ipft,ipa)    &
                                        + (1.0 - cpatch%cophoto%fs_open(ico))                      &
                                          * csite%A_c_max(ipft,ipa))                       &
                                      + cpatch%coresp%leaf_respiration(ico)) / cpatch%costate%nplant(ico) * umol_2_kgC * dtlsm

      else
         !----- If the cohort wasn't solved, we must assign some zeroes. ------------------!
         cpatch%cophoto%A_open(ico)               = 0.0
         cpatch%cophoto%A_closed(ico)             = 0.0
         cpatch%cophoto%psi_open(ico)             = 0.0
         cpatch%cophoto%psi_closed(ico)           = 0.0
         cpatch%cophoto%water_supply(ico)         = 0.0
         cpatch%cophoto%gsw_open(ico)             = 0.0
         cpatch%cophoto%gsw_closed(ico)           = 0.0
         cpatch%cotherm%leaf_gbh(ico)             = 0.0
         cpatch%cotherm%leaf_gbw(ico)             = 0.0
         cpatch%cophoto%stomatal_conductance(ico) = 0.0
         cpatch%cophoto%gpp(ico)                  = 0.0
         cpatch%coresp%leaf_respiration(ico)     = 0.0
         vm                               = 0.0
         limit_flag                       = 0
      end if
      
      !------------------------------------------------------------------------------------!
      !    Not really a part of the photosynthesis scheme, but this will do it.  We must   !
      ! integrate the "mean" of the remaining respiration terms, except for the root one.  !
      ! This is done regardless on whether the cohort is doing photosynthesis.  Also, we   !
      ! convert units so all fast respiration terms are in [µmol/m²ground/s].              !
      !------------------------------------------------------------------------------------!
      cpatch%coresp%mean_growth_resp (ico) = cpatch%coresp%mean_growth_resp (ico)                        &
                                    + cpatch%coresp%growth_respiration (ico) * kgCday_2_umols     &
                                    * cpatch%costate%nplant(ico)
      cpatch%coresp%mean_storage_resp(ico) = cpatch%coresp%mean_storage_resp(ico)                        &
                                    + cpatch%coresp%storage_respiration(ico) * kgCday_2_umols     &
                                    * cpatch%costate%nplant(ico)
      cpatch%coresp%mean_vleaf_resp  (ico) = cpatch%coresp%mean_vleaf_resp  (ico)                        &
                                    + cpatch%coresp%vleaf_respiration  (ico) * kgCday_2_umols     &
                                    * cpatch%costate%nplant(ico)                                    
      !------------------------------------------------------------------------------------!

   end do cohortloop

   !---------------------------------------------------------------------------------------!
   !     Add the contribution of this time step to the average available water.            !
   !---------------------------------------------------------------------------------------!
   if (broot_tot > 1.e-20) then
      csite%avg_available_water(ipa) = csite%avg_available_water(ipa)                      &
                                     + pss_available_water / broot_tot
   !else
   !  Add nothing, the contribution of this time is zero since no cohort can transpire... 
   end if

   return
end subroutine canopy_photosynthesis
!==========================================================================================!
!==========================================================================================!




