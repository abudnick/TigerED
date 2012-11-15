!==========================================================================================!
!==========================================================================================!
module growth_balive
   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will update the alive biomass, and compute the respiration terms  !
   ! other than leaf respiration.                                                          !
   ! IMPORTANT: The order of the operations here affect the C/N budgets, so don't change   !
   !            the order of the operations unless you really know what you are doing.     !
   !---------------------------------------------------------------------------------------!
   subroutine dbalive_dt(cgrid, tfact)
      use ed_state_vars   , only : edtype                 & ! structure
                                 , polygontype            & ! structure
                                 , sitetype               & ! structure
                                 , patchtype              ! ! structure
      use pft_coms        , only : q                      & ! intent(in)
                                 , qsw                    & ! intent(in)
                                 , plant_N_supply_scale   & ! intent(in)
                                 , c2n_storage            & ! intent(in)
                                 , growth_resp_factor     & ! intent(in)
                                 , storage_turnover_rate  & ! intent(in)
                                 , phenology              ! ! intent(in)
      use physiology_coms , only : N_plant_lim            ! ! intent(in)
      use grid_coms       , only : nzg                    ! ! intent(in)
      use ed_therm_lib    , only : calc_veg_hcap          & ! function
                                 , update_veg_energy_cweh ! ! function
      use allometry       , only : area_indices           & ! subroutine
                                 , ed_biomass             ! ! function
      use mortality       , only : mortality_rates        ! ! subroutine
      use phenology_coms  , only : theta_crit             ! ! intent(in)
      use pft_coms, only: c2p_storage, c2p_alive, c2p_dead,c2p_recruit
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid
      real             , intent(in) :: tfact
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer    :: cpoly
      type(sitetype)   , pointer    :: csite
      type(patchtype)  , pointer    :: cpatch
      integer                       :: ipy
      integer                       :: isi
      integer                       :: ipa
      integer                       :: ico
      integer                       :: ipft
      real                          :: salloc
      real                          :: salloci
      real                          :: bl
      real                          :: br
      real                          :: daily_C_gain
      real                          :: carbon_balance
      real                          :: carbon_balance_pot
      real                          :: carbon_balance_max
      real                          :: balive_in
      real                          :: nitrogen_supply
      real                          :: dndt
      real                          :: old_leaf_hcap
      real                          :: old_wood_hcap
      real                          :: nitrogen_uptake
      real :: P_uptake
      real                          :: N_uptake_pot
      real                          :: temp_dep
real :: new_soil_P, new_plant_P, old_soil_P, old_plant_P
      !------------------------------------------------------------------------------------!


   old_soil_P = 0.
   old_plant_P = 0.
   new_soil_P = 0.
   new_plant_P = 0.

      do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               !----- Reset averaged variables. -------------------------------------------!
               csite%total_plant_nitrogen_uptake(ipa) = 0.0
               csite%total_plant_P_uptake(ipa) = 0.0

! Phosphorus budgets
old_soil_P = old_soil_P + (csite%sbgc%slow_soil_P(ipa)+csite%sbgc%fast_soil_P(ipa) + csite%sbgc%struct_soil_P(ipa) + csite%sbgc%miner_soil_P(ipa)) * csite%area(ipa)
old_plant_P = old_plant_P + csite%area(ipa) * (csite%repro(1,ipa)/c2p_recruit(1)+csite%repro(2,ipa)/c2p_recruit(2)+csite%repro(3,ipa)/c2p_recruit(3)+csite%repro(4,ipa)/c2p_recruit(4))

               !----- Loop over cohorts. --------------------------------------------------!
               do ico = 1,cpatch%ncohorts

                  !----- Alias for current PFT. -------------------------------------------!
                  ipft = cpatch%pft(ico)

! Phosphorus budgets
old_plant_P = old_plant_P + csite%area(ipa) * cpatch%nplant(ico) * (cpatch%balive(ico)/c2p_alive(cpatch%pft(ico)) + cpatch%bdead(ico)/c2p_dead(cpatch%pft(ico))+cpatch%bstorage(ico)/c2p_storage(cpatch%pft(ico)))

                  !----- Update the elongation factor. ------------------------------------!
                  select case (phenology(ipft))
                  case (4)
                     cpatch%elongf(ico) = max(0.0, min(1.0, cpatch%paw_avg(ico)            &
                                                          / theta_crit))
                  case default
                     cpatch%elongf(ico) = 1.0

                  end select


                  !----- Set allocation factors. ------------------------------------------!
                  salloc  = 1.0 + qsw(ipft) * cpatch%hite(ico) + q(ipft)
                  salloci = 1.0 / salloc

                  !------------------------------------------------------------------------!
                  !     Compute maintenance costs using actual pools.                      !
                  !------------------------------------------------------------------------!
                  call plant_maintenance(cpatch,ico,cpatch%broot(ico),cpatch%bleaf(ico)    &
                                        ,tfact,daily_C_gain,csite%avg_daily_temp(ipa))
!cpatch%leaf_maintenance(ico) = 0.                  
!cpatch%root_maintenance(ico) = 0.                  
!daily_C_gain = 0.
                  !----- Subtract maintenance costs from pools. ---------------------------!
                  cpatch%balive(ico)    = cpatch%balive(ico)                               &
                                        - cpatch%leaf_maintenance(ico)                     &
                                        - cpatch%root_maintenance(ico)
                  cpatch%bleaf(ico)     = cpatch%bleaf(ico)                                &
                                        - cpatch%leaf_maintenance(ico)       
                  cpatch%broot(ico)     = cpatch%broot(ico)                                &
                                        - cpatch%root_maintenance(ico)
                  cpatch%cb(13,ico)     = cpatch%cb(13,ico)                                &
                                        - cpatch%leaf_maintenance(ico)                     &
                                        - cpatch%root_maintenance(ico)
                  cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico)                            &
                                        - cpatch%leaf_maintenance(ico)                     &
                                        - cpatch%root_maintenance(ico)

                  !------------------------------------------------------------------------!
                  !    Storage respiration/turnover_rate.                                  !
                  !    Calculate in same way as leaf and root turnover in kgC/plant/year.  !
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     The commented line is an experimental and arbitrary test, borrowed !
                  ! from maintainence temperature dependency. [[MCD]]                      !
                  !------------------------------------------------------------------------!
                  ! temp_dep = 1.0                                                         &
                  !          / ( 1.0  + exp( 0.4 * (278.15 - csite%avg_daily_temp(ipa))))
                  temp_dep = 1.0
                  !------------------------------------------------------------------------!

                  cpatch%storage_respiration(ico) = cpatch%bstorage(ico)                   &
                                                  * storage_turnover_rate(ipft)            &
                                                  * tfact * temp_dep

                  cpatch%bstorage(ico) = cpatch%bstorage(ico)                              &
                                         - cpatch%storage_respiration(ico)

                  !------------------------------------------------------------------------!
                  !     When storage carbon is lost, allow the associated nitrogen to go   !
                  ! to litter in order to maintain prescribed C2N ratio.                   !
                  !------------------------------------------------------------------------!
                  csite%fsn_in(ipa) = csite%fsn_in(ipa) + cpatch%storage_respiration(ico)  &
                                                        / c2n_storage * cpatch%nplant(ico)

                  csite%fsp_in(ipa) = csite%fsp_in(ipa) + cpatch%storage_respiration(ico) &
                       / c2p_storage(ipft) * cpatch%nplant(ico)

                  !------------------------------------------------------------------------!
                  !      Calculate actual, potential and maximum carbon balances.          !
                  !------------------------------------------------------------------------!
                  call plant_carbon_balances(cpatch,ipa,ico,daily_C_gain,carbon_balance    &
                                            ,carbon_balance_pot,carbon_balance_max)
!carbon_balance = 0.
!carbon_balance_pot = 0.
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Compute respiration rates for coming day [kgC/plant/day].         !
                  !------------------------------------------------------------------------!
                  cpatch%growth_respiration(ico) = max(0.0, daily_C_gain                   &
                                                          * growth_resp_factor(ipft))
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Find the "virtual" leaf respiration.                               !
                  !------------------------------------------------------------------------!
                  cpatch%vleaf_respiration(ico) = (1.0-cpoly%green_leaf_factor(ipft,isi))  &
                                                * salloci * cpatch%balive(ico)             &
                                                * storage_turnover_rate(ipft)              &
                                                * tfact * temp_dep
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Allocate plant carbon balance to balive and bstorage.             !
                  !------------------------------------------------------------------------!
                  balive_in = cpatch%balive(ico)
                  call alloc_plant_c_balance(csite,ipa,ico,salloc,salloci,carbon_balance   &
                                            ,nitrogen_uptake, P_uptake                 &
                                            ,cpoly%green_leaf_factor(ipft,isi))
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Do a shadow calculation to see what would have happened if stomata !
                  ! were open.  This is used to calculate potential nitrogen uptake,       !
                  ! N_uptake_pot.                                                          !
                  !------------------------------------------------------------------------!
                  if (N_plant_lim == 1) then
                     call potential_N_uptake(cpatch,ico,salloc,salloci,balive_in           &
                                            ,carbon_balance_pot,N_uptake_pot               &
                                            ,cpoly%green_leaf_factor(ipft,isi))
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !  Increment the [kgN/m2] taken up during previous day.                  !
                  !------------------------------------------------------------------------!
                  csite%total_plant_nitrogen_uptake(ipa) =                                 &
                                       csite%total_plant_nitrogen_uptake(ipa)              &
                                     + nitrogen_uptake * cpatch%nplant(ico)

                  !------------------------------------------------------------------------!
                  !------------------------------------------------------------------------!
                  !  Increment the [kgP/m2] taken up during previous day.                  !
                  !------------------------------------------------------------------------!
                  csite%total_plant_P_uptake(ipa) =                                 &
                                       csite%total_plant_P_uptake(ipa)              &
                                     + P_uptake * cpatch%nplant(ico)
                  !------------------------------------------------------------------------!

                  !----- Calculate plant N limitation factor. -----------------------------!
                  if (n_plant_lim == 0 .or. N_uptake_pot <= 0.0) then
                     cpatch%fsn(ico) = 1.0
                  else
                     nitrogen_supply = plant_N_supply_scale * cpatch%broot(ico)            &
                                     * csite%sbgc%miner_soil_N(ipa)
                     cpatch%fsn(ico) = nitrogen_supply                                     &
                                     / (nitrogen_supply + N_uptake_pot)
                  end if
                  
                  !------------------------------------------------------------------------!
                  !      Do mortality --- note that only frost mortality changes daily.    !
                  !------------------------------------------------------------------------!
                  call mortality_rates(cpatch,ipa,ico,csite%avg_daily_temp(ipa)            &
                                      ,csite%age(ipa))
                  dndt = - sum(cpatch%mort_rate(:,ico)) * cpatch%nplant(ico) * tfact

                  !----- Update monthly mortality rate [plants/m2/month]. -----------------!
                  cpatch%monthly_dndt(ico) = cpatch%monthly_dndt(ico) + dndt


                  !----- Updating LAI, WPA, and WAI. --------------------------------------!
                  call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico)                   &
                                   ,cpatch%bdead(ico),cpatch%balive(ico),cpatch%dbh(ico)   &
                                   ,cpatch%hite(ico) ,cpatch%pft(ico),cpatch%sla(ico)      &
                                   ,cpatch%lai(ico),cpatch%wpa(ico),cpatch%wai(ico)        &
                                   ,cpatch%crown_area(ico),cpatch%bsapwood(ico))

                  !----- Update above-ground biomass. -------------------------------------!
                  cpatch%agb(ico) = ed_biomass(cpatch%bdead(ico),cpatch%balive(ico)        &
                                              ,cpatch%bleaf(ico),cpatch%pft(ico)           &
                                              ,cpatch%hite(ico),cpatch%bstorage(ico)       &
                                              ,cpatch%bsapwood(ico))

                  !------------------------------------------------------------------------!
                  !     It is likely that biomass has changed, therefore, update           !
                  ! vegetation energy and heat capacity.                                   !
                  !------------------------------------------------------------------------!
                  old_leaf_hcap         = cpatch%leaf_hcap(ico)
                  old_wood_hcap         = cpatch%wood_hcap(ico)
                  call calc_veg_hcap(cpatch%bleaf(ico) ,cpatch%bdead(ico)                  &
                                    ,cpatch%bsapwood(ico),cpatch%nplant(ico)               &
                                    ,cpatch%pft(ico)                                       &
                                    ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico))
                  call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap)
                  !----- Update the stability status. -------------------------------------!
                  call is_resolvable(csite,ipa,ico,cpoly%green_leaf_factor(:,isi))
                  !------------------------------------------------------------------------!
! Phosphorus budgets
new_plant_P = new_plant_P + csite%area(ipa) * cpatch%nplant(ico) * (cpatch%balive(ico)/c2p_alive(cpatch%pft(ico)) + cpatch%bdead(ico)/c2p_dead(cpatch%pft(ico))+cpatch%bstorage(ico)/c2p_storage(cpatch%pft(ico)))
               end do
               
               !----- Update litter. ------------------------------------------------------!
               call litter(csite,ipa)
               
               !----- Update patch LAI, WAI, height, roughness... -------------------------!
               call update_patch_derived_props(csite,cpoly%lsl(isi),cpoly%met(isi)%prss,ipa)

               !----- Recalculate storage terms (for budget assessment). ------------------!
               call update_budget(csite,cpoly%lsl(isi),ipa,ipa)

               !----- It's a new day, reset average daily temperature. --------------------!
               csite%avg_daily_temp(ipa) = 0.0 
! Phosphorus budgets
new_soil_P = new_soil_P + (csite%sbgc%slow_soil_P(ipa)+csite%sbgc%fast_soil_P(ipa) + csite%sbgc%struct_soil_P(ipa) + csite%sbgc%miner_soil_P(ipa) + csite%fsp_in(ipa) + csite%stsp_in(ipa) - csite%total_plant_P_uptake(ipa)) * csite%area(ipa)
new_plant_P = new_plant_P + csite%area(ipa) * (csite%repro(1,ipa)/c2p_recruit(1)+csite%repro(2,ipa)/c2p_recruit(2)+csite%repro(3,ipa)/c2p_recruit(3)+csite%repro(4,ipa)/c2p_recruit(4))
            end do
         end do
      end do
!print*,old_soil_P, old_plant_P, old_soil_P+old_plant_P
!print*,new_soil_P, new_plant_P, new_soil_P+new_plant_P

      return
   end subroutine dbalive_dt
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will compute the respiration terms other than leaf                !
   ! respiration, plus the carbon balance and maintenance costs but without                !
   ! updating the pools.                                                                   !
   !---------------------------------------------------------------------------------------!
   subroutine dbalive_dt_eq_0(cgrid, tfact)
      use ed_state_vars   , only : edtype                 & ! structure
                                 , polygontype            & ! structure
                                 , sitetype               & ! structure
                                 , patchtype              ! ! structure
      use pft_coms        , only : q                      & ! intent(in)
                                 , qsw                    & ! intent(in)
                                 , plant_N_supply_scale   & ! intent(in)
                                 , c2n_storage            & ! intent(in)
                                 , growth_resp_factor     & ! intent(in)
                                 , storage_turnover_rate  & ! intent(in)
                                 , phenology              ! ! intent(in)
      use physiology_coms , only : N_plant_lim            ! ! intent(in)
      use grid_coms       , only : nzg                    ! ! intent(in)
      use ed_therm_lib    , only : calc_veg_hcap          & ! function
                                 , update_veg_energy_cweh ! ! function
      use allometry       , only : area_indices           & ! subroutine
                                 , ed_biomass             ! ! function
      use mortality       , only : mortality_rates        ! ! subroutine
      use phenology_coms  , only : theta_crit             ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid
      real             , intent(in) :: tfact
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer    :: cpoly
      type(sitetype)   , pointer    :: csite
      type(patchtype)  , pointer    :: cpatch
      integer                       :: ipy
      integer                       :: isi
      integer                       :: ipa
      integer                       :: ico
      integer                       :: ipft
      real                          :: salloc
      real                          :: salloci
      real                          :: bl
      real                          :: br
      real                          :: daily_C_gain
      real                          :: carbon_balance
      real                          :: carbon_balance_pot
      real                          :: carbon_balance_max
      real                          :: balive_in
      real                          :: nitrogen_supply
      real                          :: dndt
      real                          :: old_leaf_hcap
      real                          :: old_wood_hcap
      real                          :: nitrogen_uptake, P_uptake
      real                          :: N_uptake_pot
      !------------------------------------------------------------------------------------!


      do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               !----- Reset averaged variables. -------------------------------------------!
               csite%total_plant_nitrogen_uptake(ipa) = 0.0
               csite%total_plant_P_uptake(ipa) = 0.0

               !----- Loop over cohorts. --------------------------------------------------!
               do ico = 1,cpatch%ncohorts

                  !----- Alias for current PFT. -------------------------------------------!
                  ipft = cpatch%pft(ico)

                  !----- Update the elongation factor. ------------------------------------!
                  select case (phenology(ipft))
                  case (4)
                     cpatch%elongf(ico) = max(0.0, min(1.0, cpatch%paw_avg(ico)/theta_crit))
                  case default
                     cpatch%elongf(ico) = 1.0
                  end select

                  !----- Initialize cohort nitrogen uptake. -------------------------------!
                  nitrogen_uptake = 0.0
                  N_uptake_pot    = 0.0
                  P_uptake = 0.

                  !----- Set allocation factors. ------------------------------------------!
                  salloc  = 1.0 + qsw(ipft) * cpatch%hite(ico) + q(ipft)
                  salloci = 1.0 / salloc
                  
                  !----- Leaf and root biomass. -------------------------------------------!
                  bl = cpatch%bleaf(ico)
                  br = cpatch%broot(ico)

                  !------------------------------------------------------------------------!
                  !     Compute maintenance costs.                                         !
                  !------------------------------------------------------------------------!
                  call plant_maintenance(cpatch,ico,br,bl,tfact,daily_C_gain               &
                                        ,csite%avg_daily_temp(ipa))

                  !----- Subtract maintenance costs from balive. --------------------------!
                  cpatch%cb(13,ico)     = cpatch%cb(13,ico)                                &
                                        - cpatch%leaf_maintenance(ico)                     &
                                        - cpatch%root_maintenance(ico)
                  cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico)                            &
                                        - cpatch%leaf_maintenance(ico)                     &
                                        - cpatch%root_maintenance(ico)

                  !------------------------------------------------------------------------!
                  !      Calculate actual, potential and maximum carbon balances.          !
                  !------------------------------------------------------------------------!
                  call plant_carbon_balances(cpatch,ipa,ico,daily_C_gain,carbon_balance    &
                                            ,carbon_balance_pot,carbon_balance_max)

                  !------------------------------------------------------------------------!
                  !      Compute respiration rates for coming day [kgC/plant/day].         !
                  !------------------------------------------------------------------------!
                  cpatch%growth_respiration(ico)  = max(0.0, daily_C_gain                  &
                                                           * growth_resp_factor(ipft))
                  cpatch%storage_respiration(ico) = cpatch%bstorage(ico)                   &
                                                  * storage_turnover_rate(ipft) * tfact
                  cpatch%vleaf_respiration(ico) =                                          &
                                        (1.0 - cpoly%green_leaf_factor(ipft,isi))          &
                                      / (1.0 + q(ipft) + qsw(ipft) * cpatch%hite(ico))     &
                                      * cpatch%balive(ico) * storage_turnover_rate(ipft)   &
                                      * tfact

                  !------------------------------------------------------------------------!
                  !     Do a shadow calculation to see what would have happened if stomata !
                  ! were open.  This is used to calculate potential nitrogen uptake,       !
                  ! N_uptake_pot.                                                          !
                  !------------------------------------------------------------------------!
                  if (N_plant_lim == 1) then
                     call potential_N_uptake(cpatch,ico,salloc,salloci,balive_in           &
                                            ,carbon_balance_pot,N_uptake_pot               &
                                            ,cpoly%green_leaf_factor(ipft,isi))
                  end if

                  !------------------------------------------------------------------------!
                  !  Increment the [kgN/m2] taken up during previous day.                  !
                  !------------------------------------------------------------------------!
                  csite%total_plant_nitrogen_uptake(ipa) =                                 &
                                                   csite%total_plant_nitrogen_uptake(ipa)  &
                                                 + nitrogen_uptake * cpatch%nplant(ico)

                  !------------------------------------------------------------------------!
                  !  Increment the [kgP/m2] taken up during previous day.                  !
                  !------------------------------------------------------------------------!
                  csite%total_plant_P_uptake(ipa) =                                 &
                                                   csite%total_plant_P_uptake(ipa)  &
                                                 + P_uptake * cpatch%nplant(ico)

                  !----- Calculate plant N limitation factor. -----------------------------!
                  if (n_plant_lim == 0 .or. N_uptake_pot <= 0.0) then
                     cpatch%fsn(ico) = 1.0
                  else
                     nitrogen_supply = plant_N_supply_scale * br                           &
                                     * csite%sbgc%miner_soil_N(ipa)
                     cpatch%fsn(ico) = nitrogen_supply / (nitrogen_supply + N_uptake_pot)
                  end if
                  
                  !------------------------------------------------------------------------!
                  !      Do mortality --- note that only frost mortality changes daily.    !
                  !------------------------------------------------------------------------!
                  call mortality_rates(cpatch,ipa,ico,csite%avg_daily_temp(ipa)            &
                                      ,csite%age(ipa))
               end do

               !----- It's a new day, reset average daily temperature. --------------------!
               csite%avg_daily_temp(ipa) = 0.0 
            end do
         end do
      end do

      return
   end subroutine dbalive_dt_eq_0
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will transfer some of the stored carbon to balive in order to put  !
   ! the plant back on allometry.                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine transfer_C_from_storage(cpatch,ico,salloc,salloci,nitrogen_uptake            &
                                     ,N_uptake_pot)
      use ed_state_vars , only : patchtype
      use pft_coms      , only : c2n_leaf    & ! intent(in)
                               , c2n_storage & ! intent(in)
                               , c2n_stem    & ! intent(in)
                               , q           & ! intent(in)
                               , qsw         ! ! intent(in)
      use decomp_coms   , only : f_labile    ! ! intent(in)
      use allometry     , only : dbh2bl      ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: salloc
      real           , intent(in)    :: salloci
      real           , intent(inout) :: nitrogen_uptake
      real           , intent(inout) :: N_uptake_pot
      !----- Local variables. -------------------------------------------------------------!
      integer                        :: ipft
      real                           :: off_allometry_cb
      real                           :: increment
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Only do the transfer if leaves exist.                                          !
      !------------------------------------------------------------------------------------!
      if (cpatch%phenology_status(ico) == 2) return
     
      !----- Alias for pft type. ----------------------------------------------------------!
      ipft = cpatch%pft(ico)
     
      !----- Determine how much biomass we need to go back to allometry. ------------------!
      off_allometry_cb = dbh2bl(cpatch%dbh(ico),ipft) * salloc - cpatch%balive(ico)

      !----- If plants have storage, transfer it to balive. -------------------------------!
      increment            = max(0.0,min(max(0.0, off_allometry_cb),cpatch%bstorage(ico)))
      cpatch%balive(ico)   = cpatch%balive(ico)   + increment
      cpatch%bstorage(ico) = cpatch%bstorage(ico) - increment

      !----- Compute sapwood and fine root biomass. ---------------------------------------!
      cpatch%broot(ico)    = q(ipft) * cpatch%balive(ico) * salloci
      cpatch%bsapwood(ico) = qsw(ipft) * cpatch%hite(ico) * cpatch%balive(ico) * salloci

      !------------------------------------------------------------------------------------!
      !      N uptake is required since c2n_leaf < c2n_storage.  Units are kgN/plant/day.  !
      !------------------------------------------------------------------------------------!
      nitrogen_uptake = increment * (        f_labile(ipft)  / c2n_leaf(ipft)              &
                                    + (1.0 - f_labile(ipft)) / c2n_stem(ipft)              &
                                    -  1.0 / c2n_storage)
      N_uptake_pot    = nitrogen_uptake

      return
   end subroutine transfer_C_from_storage
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine plant_maintenance(cpatch,ico,br,bl,tfact,daily_C_gain,tempk)
      use ed_state_vars, only : patchtype          ! ! structure
      use pft_coms     , only : phenology          & ! intent(in)
                              , root_turnover_rate & ! intent(in)
                              , leaf_turnover_rate ! ! intent(in)
      use consts_coms  , only : umol_2_kgC         & ! intent(in)
                              , day_sec            ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: br
      real           , intent(in)    :: bl
      real           , intent(in)    :: tfact
      real           , intent(in)    :: tempk
      real           , intent(out)   :: daily_C_gain
      !----- Local variables. -------------------------------------------------------------!
      integer                        :: ipft
      real                           :: maintenance_temp_dep
      !------------------------------------------------------------------------------------!

      !------ Alias for plant functional type. --------------------------------------------!
      ipft = cpatch%pft(ico)

      !------ Get the temperature dependence. ---------------------------------------------!
      if (phenology(ipft) == 0) then
         maintenance_temp_dep = 1.0 / (1.0 + exp(0.4 * (278.15 - tempk)))
      else
         maintenance_temp_dep = 1.0
      end if

      !----- Calculate maintenance demand (kgC/plant/year). -------------------------------!
      cpatch%root_maintenance(ico) = root_turnover_rate(ipft) * br * maintenance_temp_dep
      if (phenology(ipft) /= 3)then
         cpatch%leaf_maintenance(ico) = leaf_turnover_rate(ipft) * bl * maintenance_temp_dep
      else
         cpatch%leaf_maintenance(ico) = leaf_turnover_rate(ipft) * bl                      &
                                      * cpatch%turnover_amp(ico) * maintenance_temp_dep
      end if


      !----- Convert units of maintenance to [kgC/plant/day]. -----------------------------!
      cpatch%leaf_maintenance(ico) = cpatch%leaf_maintenance(ico) * tfact
      cpatch%root_maintenance(ico) = cpatch%root_maintenance(ico) * tfact


      !----- Compute daily C uptake [kgC/plant/day]. --------------------------------------!
      if(cpatch%nplant(ico) > tiny(1.0)) then
!         daily_C_gain = umol_2_kgC * day_sec * ( cpatch%today_gpp(ico)                     &
!                                               - cpatch%today_leaf_resp(ico)               &
!                                               - cpatch%today_root_resp(ico))              &
!                                             / cpatch%nplant(ico)
         daily_C_gain = cpatch%today_gpp(ico) - cpatch%today_leaf_resp(ico)               &
              - cpatch%today_root_resp(ico)
      else
         daily_C_gain = 0.0
      end if

      return
   end subroutine plant_maintenance
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine plant_carbon_balances(cpatch,ipa,ico,daily_C_gain,carbon_balance             &
                                   ,carbon_balance_pot,carbon_balance_max)
      use ed_state_vars, only : patchtype          ! ! structure
      use pft_coms     , only : growth_resp_factor ! ! intent(in)
      use consts_coms  , only : umol_2_kgC         & ! intent(in)
                              , day_sec            ! ! intent(in)
      use ed_misc_coms , only : current_time       ! ! intent(in)
      use ed_max_dims  , only : n_pft              ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype)          , target      :: cpatch
      integer                  , intent(in)  :: ipa
      integer                  , intent(in)  :: ico
      real                     , intent(in)  :: daily_C_gain
      real                     , intent(out) :: carbon_balance
      real                     , intent(out) :: carbon_balance_pot
      real                     , intent(out) :: carbon_balance_max
      !----- Local variables. -------------------------------------------------------------!
      real                                   :: daily_C_gain_pot
      real                                   :: daily_C_gain_max
      real                                   :: growth_respiration_pot
      real                                   :: growth_respiration_max
      integer                                :: ipft
      !----- Local constants. -------------------------------------------------------------!
      logical                  , parameter   :: print_debug = .false.
      !----- Locally saved variables. -----------------------------------------------------!
      logical, dimension(n_pft), save        :: first_time  = .true.
      !------------------------------------------------------------------------------------!

      !----- Alias for PFT type. ----------------------------------------------------------!
      ipft = cpatch%pft(ico)

      !------ Calculate actual daily carbon balance: kgC/plant/day. -----------------------!
      carbon_balance = daily_C_gain - cpatch%growth_respiration(ico)                       &
                                    - cpatch%vleaf_respiration(ico)

      if (cpatch%nplant(ico) > tiny(1.0)) then

         !---------------------------------------------------------------------------------!
         !      Calculate potential carbon balance (used for nitrogen demand function).    !
         ! [kgC/plant/day].                                                                !
         !---------------------------------------------------------------------------------!
!         daily_C_gain_pot       = umol_2_kgC * day_sec * ( cpatch%today_gpp_pot(ico)       &
!                                                         - cpatch%today_leaf_resp(ico)     &
!                                                         - cpatch%today_root_resp(ico))    &
!                                                       / cpatch%nplant(ico)
         daily_C_gain_pot       = cpatch%today_gpp_pot(ico)       &
              - cpatch%today_leaf_resp(ico)     &
              - cpatch%today_root_resp(ico)
         growth_respiration_pot = max(0.0, daily_C_gain_pot * growth_resp_factor(ipft))
         carbon_balance_pot     = daily_C_gain_pot - growth_respiration_pot                &
                                - cpatch%vleaf_respiration(ico)

         !----- Calculate maximum carbon balance (used for mortality). --------------------!
!         daily_C_gain_max       = umol_2_kgC * day_sec * ( cpatch%today_gpp_max(ico)       &
!                                                         - cpatch%today_leaf_resp(ico)     &
!                                                         - cpatch%today_root_resp(ico) )   &
!                                                       / cpatch%nplant(ico)
         daily_C_gain_max       = cpatch%today_gpp_max(ico)       &
              - cpatch%today_leaf_resp(ico)     &
              - cpatch%today_root_resp(ico)
         growth_respiration_max = max(0.0, daily_C_gain_max * growth_resp_factor(ipft))
         carbon_balance_max     = daily_C_gain_max - growth_respiration_max                &
                                                   - cpatch%vleaf_respiration(ico)
      else
         carbon_balance_max = 0.0
         carbon_balance_pot = 0.0
      end if

      !----- Carbon balances for mortality. -----------------------------------------------!
      cpatch%cb(13,ico)     = cpatch%cb(13,ico)     + carbon_balance
      cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico) + carbon_balance_max

!      if (print_debug) then

!         if (first_time(ipft)) then
!            first_time(ipft) = .false.
!            write (unit=30+ipft,fmt='(a10,15(1x,a12))')                                    &
!                '      TIME','       PATCH','      COHORT','      NPLANT','    CB_TODAY'   &
!                            ,' GROWTH_RESP','  VLEAF_RESP','   TODAY_GPP','TODAY_GPPMAX'   &
!                            ,'  TODAY_LEAF','  TODAY_ROOT',' CBMAX_TODAY','          CB'   &
!                            ,'       CBMAX','  LEAF_MAINT','  ROOT_MAINT'
!         end if

!         write(unit=30+ipft,fmt='(2(i2.2,a1),i4.4,2(1x,i12),13(1x,es12.5))')               &
!              current_time%month,'/',current_time%date,'/',current_time%year               &
!             ,ipa,ico,cpatch%nplant(ico),carbon_balance,cpatch%growth_respiration(ico)     &
!             ,cpatch%vleaf_respiration(ico),cpatch%today_gpp(ico)                          &
!             ,cpatch%today_gpp_max(ico),cpatch%today_leaf_resp(ico)                        &
!             ,cpatch%today_root_resp(ico),carbon_balance_max,cpatch%cb(13,ico)             &
!             ,cpatch%cb_max(13,ico),cpatch%leaf_maintenance(ico)                           &
!             ,cpatch%root_maintenance(ico)
!      end if

      return
   end subroutine plant_carbon_balances
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   subroutine alloc_plant_c_balance(csite,ipa,ico,salloc,salloci,carbon_balance            &
                                   ,nitrogen_uptake, P_uptake,green_leaf_factor)
      use ed_state_vars , only : sitetype     & ! structure
                               , patchtype    ! ! structure
      use pft_coms      , only : c2n_storage  & ! intent(in)
                               , c2n_leaf     & ! intent(in)
                               , sla          & ! intent(in)
                               , q            & ! intent(in)
                               , qsw          & ! intent(in)
                               , c2n_stem     ! ! intent(in)
      use decomp_coms   , only : f_labile     ! ! intent(in)
      use allometry     , only : dbh2bl       ! ! function
      use pft_coms, only: c2p_storage, c2p_alive, c2p_dead
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype) , target        :: csite
      integer        , intent(in)    :: ipa
      integer        , intent(in)    :: ico
      real           , intent(in)    :: salloc
      real           , intent(in)    :: salloci
      real           , intent(in)    :: carbon_balance
      real           , intent(inout) :: nitrogen_uptake, P_uptake
      real           , intent(in)    :: green_leaf_factor
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype), pointer       :: cpatch
      integer                        :: ipft
      real                           :: bl_max
      real                           :: balive_max
      real                           :: bl_pot
      real                           :: increment
      real                           :: old_status
      real                           :: delta_bleaf
      real                           :: delta_broot
      real                           :: delta_bsapwood
      real                           :: available_carbon
      real                           :: f_total
      real                           :: f_bleaf
      real                           :: f_broot
      real                           :: f_bsapwood
      real                           :: f_resp
      real                           :: tr_bleaf
      real                           :: tr_broot
      real                           :: tr_bsapwood
      real                           :: bl
      logical                        :: on_allometry
      logical                        :: time_to_flush
      real :: balive_target, storage2alive, desired_carbon
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)
      
      ipft = cpatch%pft(ico) 

      ! First, have the daily carbon balance transit through storage.  Adjust N, P as necessary.
      cpatch%bstorage(ico) = cpatch%bstorage(ico) + carbon_balance
      nitrogen_uptake = carbon_balance / c2n_storage
      P_uptake = carbon_balance / c2p_storage(ipft)

      ! Null assumption is that everything that is in storage stays in storage.
      storage2alive = 0.0

      ! However, f there is carbon and leaves want to grow, let them grow.
      if(cpatch%bstorage(ico) > 0. .and. cpatch%phenology_status(ico) == 1)then

         !  Maximum bleaf that the allometric relationship would allow.  
         bl_max     = dbh2bl(cpatch%dbh(ico),ipft) * green_leaf_factor
         balive_target = bl_max * (green_leaf_factor + q(ipft) + qsw(ipft) * cpatch%hite(ico))
         
         desired_carbon = balive_target - cpatch%balive(ico)

         if(cpatch%bstorage(ico) >= desired_carbon)then
            storage2alive = desired_carbon
            cpatch%phenology_status(ico) = 0
         else
            storage2alive = cpatch%bstorage(ico)
         endif

         ! This process is going to require some additional nutrients.
         nitrogen_uptake = nitrogen_uptake + storage2alive * (1./c2n_leaf(ipft) - 1./c2n_storage)
         P_uptake = P_uptake + storage2alive * (1./c2p_alive(ipft) - 1./c2p_storage(ipft))
      endif

      ! If plants have negative storage, take from balive.
      if(cpatch%bstorage(ico) < 0.)then
         storage2alive = cpatch%bstorage(ico)
         if(cpatch%phenology_status(ico) == 0)cpatch%phenology_status(ico) = 1

         ! This process implies the loss of some additional nutrients.
         csite%fsp_in(ipa) = csite%fsp_in(ipa) - cpatch%nplant(ico) *   &
              storage2alive * (1./c2p_alive(ipft) - 1./c2p_storage(ipft))
         csite%fsn_in(ipa) = csite%fsn_in(ipa) - cpatch%nplant(ico) *   &
              storage2alive * (1./c2n_leaf(ipft) - 1./c2n_storage)

      endif

      cpatch%bstorage(ico) = cpatch%bstorage(ico) - storage2alive
      cpatch%balive(ico) = cpatch%balive(ico) + storage2alive
      
      cpatch%bleaf(ico)    = green_leaf_factor * cpatch%balive(ico) / (1.0 + q(ipft)   &
           + qsw(ipft) * cpatch%hite(ico))
      cpatch%broot(ico)    = q(ipft) * cpatch%balive(ico) / (1.0 + q(ipft)   &
           + qsw(ipft) * cpatch%hite(ico))
      cpatch%bsapwood(ico) = qsw(ipft) * cpatch%hite(ico) * cpatch%balive(ico) /  &
           (1.0 + q(ipft) + qsw(ipft) * cpatch%hite(ico))
      
      return
   end subroutine alloc_plant_c_balance
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   subroutine potential_N_uptake(cpatch,ico,salloc,salloci,balive_in,carbon_balance_pot    &
                                ,N_uptake_pot,green_leaf_factor)
      use ed_state_vars , only : patchtype    ! ! structure
      use pft_coms      , only : c2n_storage  & ! intent(in)
                               , c2n_leaf     & ! intent(in)
                               , c2n_stem     ! ! intent(in)
      use decomp_coms   , only : f_labile     ! ! intent(in)
      use allometry     , only : dbh2bl       ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: salloc
      real           , intent(in)    :: salloci
      real           , intent(in)    :: balive_in
      real           , intent(in)    :: carbon_balance_pot
      real           , intent(inout) :: N_uptake_pot
      real           , intent(in)    :: green_leaf_factor
      !----- Local variables. -------------------------------------------------------------!
      integer                        :: ipft
      real                           :: bl_max
      real                           :: bl_pot
      real                           :: increment
      !------------------------------------------------------------------------------------!

      ipft = cpatch%pft(ico) 

      if ( cpatch%phenology_status(ico) == 0 .and. carbon_balance_pot > 0.0 ) then

         !----- Positive carbon balance with plants fully flushed. ------------------------!
         N_uptake_pot = N_uptake_pot + carbon_balance_pot / c2n_storage

      elseif (cpatch%phenology_status(ico) == 1) then

         bl_max = dbh2bl(cpatch%dbh(ico),ipft) * green_leaf_factor * cpatch%elongf(ico)
         bl_pot = cpatch%bleaf(ico) + carbon_balance_pot

         if (bl_pot > bl_max) then
            !------------------------------------------------------------------------------!
            !     This increment would take us over the limit, so we assign all that can   !
            ! go for leaves to them, and put the remainder in storage.                     !
            !------------------------------------------------------------------------------!
            increment    = carbon_balance_pot - (bl_max-cpatch%bleaf(ico))
            N_uptake_pot = N_uptake_pot + increment / c2n_storage
            increment    = bl_max-cpatch%bleaf(ico)
            N_uptake_pot = N_uptake_pot + increment                                        &
                         * ( f_labile(ipft) / c2n_leaf(ipft)                               &
                           + (1.0 - f_labile(ipft)) / c2n_stem(ipft))
         elseif (carbon_balance_pot > 0.0) then

            !------------------------------------------------------------------------------!
            !      This increment did not exceed the limit, put everything in leaves.  We  !
            ! don't compute the uptake if carbon balance is negative, just because there   !
            ! will be no uptake...                                                         !
            !------------------------------------------------------------------------------!
            N_uptake_pot = N_uptake_pot + carbon_balance_pot                               &
                         * ( f_labile(ipft) / c2n_leaf(ipft)                               &
                           + (1.0 - f_labile(ipft)) / c2n_stem(ipft))
         end if
      end if

      return
   end subroutine potential_N_uptake
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine litter(csite,ipa)

      use ed_state_vars, only : patchtype & ! structure
                              , sitetype  ! ! structure
      use pft_coms     , only : c2n_leaf  & ! intent(in)
                              , c2n_stem  & ! intent(in)
                              , l2n_stem  ! ! intent(in)
      use pft_coms, only: c2p_alive, c2p_dead
      use decomp_coms  , only : f_labile  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)  , target     :: csite
      integer         , intent(in) :: ipa
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype) , pointer    :: cpatch
      integer                      :: ico
      integer                      :: ipft
      real                         :: plant_litter
      real                         :: plant_litter_f
      real                         :: plant_litter_s
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)

      !------------------------------------------------------------------------------------!
      !      Add fine root and leaf turnover to the litter.                                !
      !------------------------------------------------------------------------------------!
      do ico=1,cpatch%ncohorts
         ipft = cpatch%pft(ico)

         plant_litter   = ( cpatch%leaf_maintenance(ico) + cpatch%root_maintenance(ico) )  &
                        * cpatch%nplant(ico)
         plant_litter_f = plant_litter * f_labile(ipft)
         plant_litter_s = plant_litter - plant_litter_f

         csite%fsc_in(ipa) = csite%fsc_in(ipa) + plant_litter_f
         csite%fsn_in(ipa) = csite%fsn_in(ipa) + plant_litter_f / c2n_leaf(ipft)

         csite%fsp_in(ipa) = csite%fsp_in(ipa) + plant_litter_f / c2p_alive(ipft)
         csite%stsp_in(ipa) = csite%stsp_in(ipa) + plant_litter_s / c2p_dead(ipft)

         csite%ssc_in(ipa) = csite%ssc_in(ipa) + plant_litter_s
         csite%ssl_in(ipa) = csite%ssl_in(ipa) + plant_litter_s * l2n_stem / c2n_stem(ipft)
      end do

      return
   end subroutine litter
   !=======================================================================================!
   !=======================================================================================!
end module growth_balive
!==========================================================================================!
!==========================================================================================!
