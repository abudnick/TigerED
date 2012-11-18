!==========================================================================================!
!==========================================================================================!
!     This module contains the main variable structures in ED, and several important sub-  !
! routines for allocation, de-allocation.                                                  !
!------------------------------------------------------------------------------------------!
module ed_state_vars

  use grid_coms, only: nzg,nzs,ngrids
  use c34constants, only : stoma_data,n_stoma_atts
  use ed_max_dims, only: max_site,n_pft,n_dbh,n_age,n_mort,n_dist_types      &
                       ,maxmach,maxgrds, str_len
  use disturb_coms, only : lutime,num_lu_trans,max_lu_years
  use met_driver_coms, only: met_driv_data,met_driv_state
  use fusion_fission_coms, only: ff_nhgt,hgt_class
  use phenology_coms, only: prescribed_phen
  use ed_misc_coms, only: idoutput, imoutput, iqoutput, ndcycle
  use soil_bgc, only: sbgc_vars
  use cohort_state, only: costate_vars
  use cohort_phen, only: cophen_vars
  use cohort_mort, only: comort_vars
  use cohort_therm, only: cotherm_vars
  use cohort_photo, only: cophoto_vars
  use cohort_rad, only: corad_vars
  use cohort_resp, only: coresp_vars

  implicit none
!============================================================================!
!============================================================================!


!============================================================================!
!============================================================================!

  !--------------------------------------------------------------------------!  
  ! Patch type:
  ! The following are arrays of cohorts that populate the current patch.
  ! This is the lowest (most infantile) part of the memory structure.
  !--------------------------------------------------------------------------!  

  type patchtype

     !-------------------------
     ! Book-keeping variables.
     !-------------------------

     ! Number of cohorts in the patch.
     integer :: ncohorts

     ! Global index of the first cohort across all cohorts.
     integer :: coglob_id

     !-------------------------
     ! Module 1: State variables
     !-------------------------

     type(costate_vars) :: costate

     !--------------------------
     ! Module 2: Phenology
     !--------------------------

     type(cophen_vars) :: cophen

     !--------------------------
     ! Module 3: Mortality
     !--------------------------

     type(comort_vars) :: comort

     !--------------------------
     ! Module 4: Thermodynamics
     !--------------------------

     type(cotherm_vars) :: cotherm

     !--------------------------
     ! Module 5: Photosynthesis
     !--------------------------

     type(cophoto_vars) :: cophoto

     !----------------------------------------------
     ! Module 6: Radiation
     !----------------------------------------------

     type(corad_vars) :: corad

     !---------------------------------
     ! Module 7: Respiration and Maintenance
     !---------------------------------

     type(coresp_vars) :: coresp

  end type patchtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  !---------------------------------------------------------------------------!  
  ! Site type:
  ! The following are the patch level arrays that populate the current site.
  !---------------------------------------------------------------------------!  
  
  type sitetype
     
     ! About 226 variables here.


     integer :: npatches

     !  The global index of the first cohort in all patches
     integer,pointer,dimension(:) :: paco_id

     ! The number of cohorts in each patch
     integer,pointer,dimension(:) :: paco_n

     ! Global index of the first patch in this vector, across all patches
     ! on the grid

     integer :: paglob_id

     ! The patches containing the cohort arrays
     type(patchtype),pointer,dimension(:) :: patch


     ! Leaf area index (m2 leaf / m2 ground)
     real, pointer,dimension(:) :: lai

     ! Wood projected area (m2 wood / m2 ground)
     real, pointer,dimension(:) :: wpa

     ! Wood area index (m2 wood / m2 ground)
     real, pointer,dimension(:) :: wai

     ! Patch type index:
     !       Agriculture = 1
     !       Secondary Forest = 2
     !       Primary Forest = 3
     integer , pointer,dimension(:) :: dist_type

     ! Time since last disturbance (years)
     real , pointer,dimension(:) :: age

     ! Fractional area of the patch
     real , pointer,dimension(:) :: area

     ! Fractional area of the patch (considers lai weighting)
     real, pointer,dimension(:) :: laiarea

     ! Number of degree days
     ! Degree days --- sum of daily average temperatures above 278.15 K 
     real , pointer,dimension(:) :: sum_dgd

     ! Chill days --- number of days with average temperatures below 278.15 K
     real , pointer,dimension(:) :: sum_chd  

     ! Flag specifying whether (1) or not (0) this patch is a plantation
     integer , pointer,dimension(:) :: plantation

     ! Ice-vapour equivalent potential temperature of canopy air space [K]
     real , pointer,dimension(:) :: can_theiv

     ! Temperature (K) of canopy air
     real , pointer,dimension(:) :: can_temp

     ! Water vapor specific humidity (kg/kg) of canopy air
     real , pointer,dimension(:) :: can_shv

     ! CO2 mixing ratio [umol/mol] of canopy air
     real , pointer,dimension(:) :: can_co2

     ! Density [kg/m3] of canopy air
     real , pointer,dimension(:) :: can_rhos

     ! Canopy air pressure.
     real , pointer,dimension(:) :: can_prss

     ! Canopy air potential temperature
     real , pointer,dimension(:) :: can_theta

     ! Canopy depth (m)
     real , pointer,dimension(:) :: can_depth

     ! Fraction of canopy that is open (---)
     real , pointer,dimension(:) :: opencan_frac
     
     ! Ground to canopy conductances
     real , pointer, dimension(:) :: ggbare
     real , pointer, dimension(:) :: ggveg
     real , pointer, dimension(:) :: ggnet
     real , pointer, dimension(:) :: ggsoil

     ! Mean light extinction coefficient, its diurnal and monthly cycle.
     real , pointer, dimension(:) :: lambda_light
     real , pointer, dimension(:) :: dmean_lambda_light
     real , pointer, dimension(:) :: mmean_lambda_light


     ! Number of cohorts in the patch
     integer , pointer,dimension(:) :: cohort_count

     ! Patch name (only used when restarting from ED1)
     character(len=str_len), pointer,dimension(:) :: pname
     
     ! Number of surface water layers
     integer,pointer,dimension(:) :: nlev_sfcwater

     ! Surface water mass (kg/m2)
     real, pointer,dimension(:,:) :: sfcwater_mass

     ! Surface water internal energy (J/kg)
     real, pointer,dimension(:,:) :: sfcwater_energy

     ! Depth of surface water (m)
     real, pointer,dimension(:,:) :: sfcwater_depth

     ! Temperature of surface water (K)
     real, pointer,dimension(:,:) :: sfcwater_tempk
    
     ! Liquid fraction of surface water
     real, pointer,dimension(:,:) :: sfcwater_fracliq

     ! Short wave radiation absorbed by the surface water (W/m2)
     real, pointer,dimension(:,:) :: rshort_s

     ! Short wave radiation absorbed by the surface water, 
     ! beam component (W/m2)
     real, pointer,dimension(:,:) :: rshort_s_beam

     ! Short wave radiation absorbed by the surface water, 
     ! diffuse component (W/m2)
     real, pointer,dimension(:,:) :: rshort_s_diffuse
     
     ! Soil internal energy (J/m3)
     real,    pointer,dimension(:,:) :: soil_energy
     
     ! Soil water (m3 water / m3 soil)
     real,    pointer,dimension(:,:) :: soil_water
          
     ! Temperature of soil (K)
     real,    pointer,dimension(:,:) :: soil_tempk
     
     ! Liquid fraction of soil
     real,    pointer,dimension(:,:) :: soil_fracliq

     ! Root density [kg/m3] over the profile
     ! Typically diagnosed for history or monthly output
     real, pointer, dimension (:,:) :: rootdense
     
     ! Effective specific humidity (kg/kg) just above soil
     real,    pointer,dimension(:) :: ground_shv
     
     ! Surface saturation specific humidity (kg/kg)
     real,    pointer,dimension(:) :: ground_ssh

     ! Ground temperature (either the top soil or top surface snow/water) (K)
     real, pointer, dimension (:) :: ground_temp

     ! Ground temperature (either the top soil or top surface snow/water) (K)
     real, pointer, dimension (:) :: ground_fliq

     ! Net roughness length (m)
     real,    pointer,dimension(:)  :: rough

     ! Maximum PAR possibly attained in this patch.  This is found by calling the radiation
     ! scheme for a single cohort (the tallest that can be solved).
     real,    pointer, dimension(:) :: par_l_max
     real,    pointer, dimension(:) :: par_l_beam_max
     real,    pointer, dimension(:) :: par_l_diffuse_max

     ! Photosynthetic rates for different PFTs, if they were at the top 
     ! of the canopy (umol/m2 leaf/s).  Used by mortality function.
     real, pointer,dimension(:,:) :: A_o_max  ! open stomata
     real, pointer,dimension(:,:) :: A_c_max  ! closed stomata
     
     ! This will hold the stomatal conductance data from the previous 
     ! time step corresponding to A_o_max
     type(stoma_data), pointer,dimension(:,:)  :: old_stoma_data_max
     real, pointer,dimension(:,:,:)            :: old_stoma_vector_max

     ! Average daily temperature [K]
     real , pointer,dimension(:) :: avg_daily_temp  

     ! average of rh [umol/m2/s] over FRQSTATE
     real , pointer,dimension(:) :: mean_rh

     ! average of rh [umol/m2/s] over a day and a month, respectively
     real , pointer,dimension(:  ) :: dmean_rh
     real , pointer,dimension(:,:) :: qmean_rh
     real , pointer,dimension(:  ) :: mmean_rh

     ! average of net ecosystem productivity (NEP) [umol/m2/s] over FRQSTATE
     real , pointer,dimension(:) :: mean_nep

     !=======================================================================!
     !      These variables are terms of the water, carbon, and energy       !
     ! budget, and will compute what enters, what stays, and what leaves ED. !
     ! They are averaged over FRQSTATE, so the units are now per unit of     !
     ! time for all fluxes, and the total storage at the beginning of the    !
     ! next time step (or the final storage between the previous analysis    !
     ! and this one).                                                        !
     !=======================================================================!
     ! Mean moisture transfer from the canopy air to the atmosphere 
     ! (kg_H2O/m2/s)
     real , pointer,dimension(:) :: wbudget_loss2atm

     ! Contribution of density change on change in the final storage
     ! (kg_H2O/m2/s)
     real , pointer,dimension(:) :: wbudget_denseffect

     ! Precipitation [kg_H2O/m2/s]
     real , pointer,dimension(:) :: wbudget_precipgain

     ! Mean runoff (kg_H2O/m2/s)
     real , pointer,dimension(:) :: wbudget_loss2runoff

     ! Mean drainage (kg_H2O/m2/s)
     real , pointer,dimension(:) :: wbudget_loss2drainage

     ! Total water (soil, sfcwater, can_shv, veg_water) at the beginning
     ! of budget-averaging time. [kg_H2O/m2]
     real , pointer,dimension(:) :: wbudget_initialstorage

     ! Residual of the water budget (kg_H2O/m2/s)
     real , pointer,dimension(:) :: wbudget_residual

     ! Mean sensible heat transfer from the canopy air to the atmosphere
     ! (J/m2/s)
     real , pointer,dimension(:) :: ebudget_loss2atm

     ! Mean change in storage due to change in density
     ! (J/m2/s)
     real , pointer,dimension(:) :: ebudget_denseffect

     ! Energy associated with runoff (J/m2/s)
     real , pointer,dimension(:) :: ebudget_loss2runoff

     ! Energy associated with drainage (J/m2/s)
     real , pointer,dimension(:) :: ebudget_loss2drainage

     ! Net absorbed radiation by soil, sfcwater, vegetation [J/m2/s]
     real , pointer,dimension(:) :: ebudget_netrad

     ! Energy associated with precipitation (J/m2/s)
     real , pointer,dimension(:) :: ebudget_precipgain

     ! Total energy (soil, sfcwater, can_shv, veg_water) at the beginning
     ! of budget-averaging time. [J/m2]
     real , pointer,dimension(:) :: ebudget_initialstorage

     ! Residual of the energy budget [J/m2/s]
     real, pointer, dimension(:) :: ebudget_residual

     ! Total CO2 (can_shv) at the beginning of budget-averaging 
     ! time. [umol_CO2/m2]
     real , pointer,dimension(:) :: co2budget_initialstorage

     ! Residual of the CO2 budget [umol_CO2/m2/s]
     real , pointer,dimension(:) :: co2budget_residual

     ! Flux of CO2 from the canopy air to the atmosphere [umol_CO2/m2/s]
     real , pointer,dimension(:) :: co2budget_loss2atm

     ! Change in CO2 total storage due to change in density [umol_CO2/m2/s]
     real , pointer,dimension(:) :: co2budget_denseffect

     ! Average GPP [umol_CO2/m2/s]
     real , pointer,dimension(:) :: co2budget_gpp

     ! Average GPP by GPP class[umol_CO2/m2/s]
     real , pointer,dimension(:,:) :: co2budget_gpp_dbh

     ! Average plant respiration [umol_CO2/m2/s]
     real , pointer,dimension(:) :: co2budget_plresp

     ! Average heterotrophic respiration [umol_CO2/m2/s]
     real , pointer,dimension(:) :: co2budget_rh

     ! Daily average residuals
     real, pointer, dimension(:) :: dmean_co2_residual    ! [umol_CO2/m2]
     real, pointer, dimension(:) :: dmean_energy_residual ! [J/m2]
     real, pointer, dimension(:) :: dmean_water_residual  ! [kg_H20/m2]
     ! Monthly average residuals
     real, pointer, dimension(:) :: mmean_co2_residual    ! [umol_CO2/m2]
     real, pointer, dimension(:) :: mmean_energy_residual ! [J/m2]
     real, pointer, dimension(:) :: mmean_water_residual  ! [kg_H20/m2]

     ! Daily average of A_decomp, the temperature and moisture dependence
     ! of heterotrophic respiration.  The "today" variable is used in the
     ! model, whereas dmean and mmean are for output only.
     real , pointer,dimension(:) :: today_A_decomp
     real , pointer,dimension(:) :: dmean_A_decomp
     real , pointer,dimension(:) :: mmean_A_decomp

     ! Daily average of the product A_decomp * f_decomp, which incorporates
     ! temperature, moisture, and N dependence of decomposition.  The "today" 
     ! variable is used in the model, whereas dmean and mmean are for output only.
     real , pointer,dimension(:) :: today_Af_decomp
     real , pointer,dimension(:) :: dmean_Af_decomp
     real , pointer,dimension(:) :: mmean_Af_decomp

     ! Carbon available to establish recruits [kgC/m2]
     real, pointer,dimension(:,:) :: repro    !(n_pft,npatches)  

     ! Vegetation roughness length (m)
     real , pointer,dimension(:) :: veg_rough

     ! Vegetation height (m)
     real , pointer,dimension(:) :: veg_height

     ! 0-plane displacement height (m)
     real , pointer,dimension(:) :: veg_displace

     ! Input to fast soil carbon pool [kgC/m2/day]
     real , pointer,dimension(:) :: fsc_in

     ! Input to structural soil carbon pool [kgC/m2/day]
     real , pointer,dimension(:) :: ssc_in

     ! Input to soil lignin pool [kgC/m2/day]
     real , pointer,dimension(:) :: ssl_in  

     ! Input to fast soil nitrogen pool [kgN/m2/day]
     real , pointer,dimension(:) :: fsn_in  

     ! Input to fast soil phosphorus pool [kgP/m2/day]
     real , pointer,dimension(:) :: fsp_in  
     ! Input to structural soil phosphorus pool [kgP/m2/day]
     real , pointer,dimension(:) :: stsp_in  

     ! Plant nitrogen update summed over all cohorts [kgN/m2/day]
     real , pointer,dimension(:) :: total_plant_nitrogen_uptake

     type(sbgc_vars) :: sbgc

     ! Plant phosphorus update summed over all cohorts [kgP/m2/day]
     real , pointer,dimension(:) :: total_plant_P_uptake

     !real, pointer,dimension(:) :: Nnet_min
     !real, pointer,dimension(:) :: Ngross_min
     real, pointer,dimension(:) :: mineralized_N_input
     real, pointer,dimension(:) :: mineralized_N_loss

     ! Short wave radiation absorbed by the ground (W/m2)
     real , pointer,dimension(:) :: rshort_g
     
     ! Short wave radiation absorbed by the ground, beam component (W/m2)
     real , pointer,dimension(:) :: rshort_g_beam
     
     ! Short wave radiation absorbed by the ground, diffuse component (W/m2)
     real , pointer,dimension(:) :: rshort_g_diffuse
     
     ! Long wave radiation absorbed by the ground (W/m2)
     real , pointer,dimension(:) :: rlong_g

     ! Long wave radiation absorbed by the ground (W/m2), due to the 
     ! surface and vegetation alone
     real , pointer,dimension(:) :: rlong_g_surf

     ! Long wave radiation absorbed by the ground (W/m2), due to the 
     ! incident long wave alone
     real , pointer,dimension(:) :: rlong_g_incid

     ! Long wave radiation absorbed by the surface water (W/m2)
     real , pointer,dimension(:) :: rlong_s

     ! Long wave radiation absorbed by the surface water (W/m2), due to 
     ! the surface and vegetation alone
     real , pointer,dimension(:) :: rlong_s_surf

     ! Long wave radiation absorbed by the surface water (W/m2), due 
     ! to the incident atmospheric long wave alone
     real , pointer,dimension(:) :: rlong_s_incid

     ! Patch albedo
     real , pointer,dimension(:) :: albedo

     ! Patch albedo, beam component
     real , pointer,dimension(:) :: albedo_beam

     ! Patch albedo, diffuse component
     real , pointer,dimension(:) :: albedo_diffuse

     ! Upward long wave radiation at the top of the canopy (W/m2)
     real , pointer,dimension(:) :: rlongup

     ! Albedo for long wave radiation
     real , pointer,dimension(:) :: rlong_albedo

     ! Total snow depth as calculated in the radiation scheme.  Used for 
     ! checking if cohorts are buried.  [m]
     real , pointer,dimension(:) :: total_sfcw_depth

     ! Fraction of vegetation covered with snow.  Used for computing 
     ! surface roughness.
     real , pointer,dimension(:) :: snowfac

     ! limitation of heterotrophic respiration due to physical 
     ! environmental factors (0-1 coefficient)
     real , pointer,dimension(:) :: A_decomp

     ! damping of decomposition due to nitrogen immobilization 
     ! (0-1 coefficient)
     real , pointer,dimension(:) :: f_decomp

     ! total heterotrophic respiration (umol/m2/s)
     real , pointer,dimension(:) :: rh

     ! coarse woody debris contribution to rh (umol/m2/s)
     real , pointer,dimension(:) :: cwd_rh

     ! Integer flag specifying whether this patch is to be fused
     integer , pointer,dimension(:) :: fuse_flag

     ! Plant density broken down into size and PFT bins.  Used in patch fusion
     real, pointer,dimension(:,:,:) :: cumlai_profile !(n_pft,ff_nhgt,npatches)

     ! Above ground biomass in this patch [kgC/m2]
     real , pointer,dimension(:) :: plant_ag_biomass

     ! Mean water flux from the canopy air to the atmosphere (kg_H2O/mÂ²/s).
     real, pointer,dimension(:) :: mean_wflux

     ! Mean latent heat flux from the canopy air to the atmosphere (W/mÂ²).
     real, pointer,dimension(:) :: mean_latflux

     ! Mean sensible heat flux from the canopy air to the atmosphere (W/mÂ²).
     real, pointer,dimension(:) :: mean_hflux
     
     ! Mean runoff (kg_H2O/mÂ²/s).
     real, pointer,dimension(:) :: mean_runoff

     ! Energy associated with runoff (W/mÂ²), averaged over FRQSTATE
     real, pointer,dimension(:) :: mean_qrunoff

     ! Last time step successfully completed by integrator.
     real, pointer,dimension(:)  :: htry

     ! Average time step used by the integrator, its daily and monthly mean
     real, pointer, dimension(:) :: avg_rk4step
     real, pointer, dimension(:) :: dmean_rk4step
     real, pointer, dimension(:) :: mmean_rk4step

     ! Average available water for transpiration.
     real, pointer,dimension(:)  :: avg_available_water

     real, pointer,dimension(:)  :: ustar ! Friction velocity [m/s]

     real, pointer,dimension(:)  :: tstar ! Characteristic temperature fluctuation scale [K]

     real, pointer,dimension(:)  :: qstar ! Characteristic specific humidity fluct. scale [kg/kg]

     real, pointer,dimension(:)  :: cstar ! Characteristic CO2 mix. ratio fluct. scale [ppm]

     real, pointer,dimension(:)  :: zeta  ! Height / Obukhov length

     real, pointer,dimension(:)  :: ribulk ! Bulk Richardson number

     real, pointer,dimension(:)  :: upwp !eddy mom. flux u-prime w-prime

     real, pointer,dimension(:)  :: tpwp !eddy heat flux t-prime w-prime

     real, pointer,dimension(:)  :: qpwp !eddy moist. flux q-prime w-prime

     real, pointer,dimension(:)  :: cpwp !eddy CO2 flux flux c-prime w-prime

     real, pointer,dimension(:)  :: wpwp !eddy v. mom. flux w-prime w-prime

     ! -----------------------------------------------------------------
     ! Fast time flux diagnostic variables
     ! The following variables are all averaged over 
     ! the fast time flux period frqfast
     !------------------------------------------------------------------
     
     !----- Radiation --------------------------------------------------
     real, pointer, dimension(:)   :: avg_rshort_gnd     ! Absorbed SW rad. of soil + top water/snow layer   [ W/m2]
     real, pointer, dimension(:)   :: avg_rlong_gnd      ! Emitted LW rad. from soil + top water/snow layer  [ W/m2]
     real, pointer, dimension(:)   :: avg_rlongup        ! Emitted LW rad. from soil + top water/snow layer  [ W/m2]
     real, pointer, dimension(:)   :: avg_albedo         ! Albedo                                            [ ----]
     real, pointer, dimension(:)   :: avg_albedo_beam    ! Direct Albedo                                     [ ----]
     real, pointer, dimension(:)   :: avg_albedo_diffuse ! Diffuse Albedo                                    [ ----]
     real, pointer, dimension(:)   :: avg_rlong_albedo   ! Longwave Albedo                                   [ ----]
     real, pointer, dimension(:)   :: dmean_albedo         ! Albedo                                          [ ----]
     real, pointer, dimension(:)   :: dmean_albedo_beam    ! Direct Albedo                                   [ ----]
     real, pointer, dimension(:)   :: dmean_albedo_diffuse ! Diffuse Albedo                                  [ ----]
     real, pointer, dimension(:)   :: mmean_albedo         ! Albedo                                          [ ----]
     real, pointer, dimension(:)   :: mmean_albedo_beam    ! Direct Albedo                                   [ ----]
     real, pointer, dimension(:)   :: mmean_albedo_diffuse ! Diffuse Albedo                                  [ ----]
     real, pointer, dimension(:,:) :: qmean_albedo         ! Albedo                                          [ ----]
     real, pointer, dimension(:,:) :: qmean_albedo_beam    ! Direct Albedo                                   [ ----]
     real, pointer, dimension(:,:) :: qmean_albedo_diffuse ! Diffuse Albedo                                  [ ----]


     ! Notation
     ! ATM - atmosphere, CAS - canopy air space

     !----- Moisture ---------------------------------------------------

     real,pointer,dimension(:)   :: avg_carbon_ac   ! Average carbon flux, ATM -> CAS,       [umol/m2/s]
     real,pointer,dimension(:)   :: avg_vapor_lc    ! Average water vapor flux, Leaf -> CAS     [kg/m2/s]
     real,pointer,dimension(:)   :: avg_vapor_wc    ! Average water vapor flux, wood -> CAS     [kg/m2/s]
     real,pointer,dimension(:)   :: avg_dew_cg      ! Average dew flux, CAS -> GND             [kg/m2/s]
     real,pointer,dimension(:)   :: avg_vapor_gc    ! Average water vapor flux, GND -> CAS     [kg/m2/s]
     real,pointer,dimension(:)   :: avg_wshed_vg    ! Average precip that falls from leaves    [kg/m2/s]
     real,pointer,dimension(:)   :: avg_intercepted ! Average precip that is intercepted       [kg/m2/s]
     real,pointer,dimension(:)   :: avg_throughfall ! Average precip that is never intercepted [kg/m2/s]
     real,pointer,dimension(:)   :: avg_vapor_ac    ! Average water vapor flux, ATM-CAS        [kg/m2/s]
     real,pointer,dimension(:)   :: avg_transp      ! Average transpiration of water vapor     [kg/m2/s]
     real,pointer,dimension(:)   :: avg_evap        ! Average evaporation from leaf and
                                                    !  and ground surfaces -> CAS              [kg/m2/s]
     real,pointer,dimension(:,:) :: avg_smoist_gg   ! Moisture flux between soil layers
                                                    !  where layer (nzg) is the flux from the
                                                    !  surface water into the top layer        [kg/m2/s]
     real,pointer,dimension(:,:) :: avg_transloss   ! Transpired soil moisture sink in each layer
                                                    !                                          [kg/m2/s]
     real,pointer,dimension(:)   :: avg_runoff      ! Average surface water runoff             [kg/m2/s]
     real,pointer,dimension(:)   :: avg_drainage    ! Average water drainage through the lower 
                                                    !  soil layer                              [kg/m2/s]
     real,pointer,dimension(:)   :: avg_drainage_heat! Average internal energy loss due to water 
                                                     ! drainage through the lower soil layer   [kg/m2/s]

     !----- Auxillary Variables (user can modify to view any variable ----------------------!
     real,pointer,dimension(:)   :: aux             ! Auxillary surface variable
     real,pointer,dimension(:,:) :: aux_s           ! Auxillary soil variable
     
     !----- Sensible heat ------------------------------------------------------------------!

     real,pointer,dimension(:) :: avg_sensible_lc   ! Sensible heat flux, Leaf -> CAS           [  W/m2]
     real,pointer,dimension(:) :: avg_sensible_wc   ! Sensible heat flux, wood -> CAS           [  W/m2]
     real,pointer,dimension(:) :: avg_qwshed_vg     ! Average precip that falls from leaves    [  W/m2]
     real,pointer,dimension(:) :: avg_qintercepted  ! Average precip that is intercepted       [  W/m2]
     real,pointer,dimension(:) :: avg_qthroughfall  ! Average precip that is never intercepted [  W/m2]
     real,pointer,dimension(:) :: avg_sensible_gc   ! Sensible heat flux, GND -> CAS           [  W/m2]
     real,pointer,dimension(:) :: avg_sensible_ac   ! Sensible heat flux, ATM -> CAS           [  W/m2]
     real,pointer,dimension(:,:) :: avg_sensible_gg ! Net soil heat flux between layers        [  W/m2]
     real,pointer,dimension(:) :: avg_runoff_heat   ! Surface runoff internal energy flux      [  W/m2]

     !----- Mass and Energy ---------------------------------------------------------------!
     real,pointer,dimension(:) :: avg_leaf_energy ! Average leaf internal energy   [  J/m2]
     real,pointer,dimension(:) :: avg_leaf_temp   ! Average leaf temperature       [     K]
     real,pointer,dimension(:) :: avg_leaf_fliq   ! Avg. liq. frac. of leaf water  [   ---]
     real,pointer,dimension(:) :: avg_leaf_water  ! Avg. water on top of leaves    [ kg/m2]
     real,pointer,dimension(:) :: avg_leaf_hcap   ! Avg. leaf heat capacity        [J/m2/K]
     real,pointer,dimension(:) :: avg_wood_energy ! Average wood internal energy   [  J/m2]
     real,pointer,dimension(:) :: avg_wood_temp   ! Average wood temperature       [     K]
     real,pointer,dimension(:) :: avg_wood_fliq   ! Avg. liq. frac. of wood water  [   ---]
     real,pointer,dimension(:) :: avg_wood_water  ! Avg. water on top of woods     [ kg/m2]
     real,pointer,dimension(:) :: avg_wood_hcap   ! Avg. wood heat capacity        [J/m2/K]

     !----- Hydrology variables -----------------------------------------------------------!

     real,pointer,dimension(:)   :: watertable
     real,pointer,dimension(:)   :: moist_dz
     real,pointer,dimension(:)   :: ksat
     real,pointer,dimension(:)   :: soil_sat_energy
     real,pointer,dimension(:)   :: soil_sat_water
     real,pointer,dimension(:)   :: soil_sat_heat
     real,pointer,dimension(:,:) :: runoff_A ! Runoff parameters, 1st dimension is 3.
     real,pointer,dimension(:)   :: runoff_rate
     real,pointer,dimension(:)   :: runoff

     real,pointer,dimension(:,:,:) :: past_paw  ! k,day,patch
     real,pointer,dimension(:,:) :: current_paw  ! k, patch
          
  end type sitetype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  !----------------------------------------------------------------------!
  ! Polygon Type:
  ! The following are the arrays of site level variables
  ! that populate the current polygon.
  ! ----------------------------------------------------

  ! About 159 variables here.
  type polygontype



     integer :: nsites

     !  The global index of the first patch in all each site
     integer,pointer,dimension(:) :: sipa_id

     ! The number of patches in each site
     integer,pointer,dimension(:) :: sipa_n

     ! Global index of the first site in this vector across all sites
     ! on the grid
     integer :: siglob_id     ! Global index of the first patch across all cohorts

     integer,pointer,dimension(:) :: patch_count    ! number of patches per site
     
     ! Pointer to the patch vectors
     type(sitetype),pointer,dimension(:) :: site

     integer,pointer,dimension(:) :: sitenum

     ! The vectorized landuse matrix is allocated in landuse_init
     type(lutime), pointer,dimension(:,:) :: clutimes !(luyears,nsites)

     !-------------------------------------------------------------------------------------!
     !-------------------------------------------------------------------------------------!
     !     Additional land use-related variables.                                          !
     !-------------------------------------------------------------------------------------!

     !-------------------------------------------------------------------------------------!
     !     The number of years of landuse data calculated at read-in of data during        !
     ! landuse_init.                                                                       !
     !-------------------------------------------------------------------------------------!
     integer,pointer,dimension(:) :: num_landuse_years
     
     !----- Minimum DBH and probability of being harvest when selective logging happens. --!
     real   , pointer, dimension(:,:) :: mindbh_primary
     real   , pointer, dimension(:,:) :: probharv_primary
     real   , pointer, dimension(:,:) :: mindbh_secondary
     real   , pointer, dimension(:,:) :: probharv_secondary
     !-------------------------------------------------------------------------------------!
     !-------------------------------------------------------------------------------------!

     real,   pointer,dimension(:,:) :: lai_pft    ! Site level mean LAI, grouped by cohort PFT
                                                  ! [m2/m2] (n_pft,nsites)
     real,   pointer,dimension(:,:) :: wpa_pft    ! Woody Projected Area, grouped by cohort PFT
                                                  ! [m2/m2] (n_pft,nsites)
     real,   pointer,dimension(:,:) :: wai_pft    ! Woody area index, grouped by cohort PFT
                                                  ! [m2/m2] (n_pft,nsites)

     !-----------------------------------
     ! BASIC INFO
     !-----------------------------------

     real,pointer,dimension(:) :: area            ! Proportion of a grid cell occupied by the site
     real,pointer,dimension(:) :: patch_area      ! Un-normalized sum of patch areas
     real,pointer,dimension(:) :: elevation       ! mean site elevation (meters)
     real,pointer,dimension(:) :: slope           ! mean site slope (degrees)
     real,pointer,dimension(:) :: aspect          ! mean site aspect (degrees)

     
     !--------------------------
     !  hydrologic routing info
     !--------------------------
     real,pointer,dimension(:) :: TCI             ! topographic convergence index

     integer,pointer,dimension(:) :: lsl          ! Index of the lowest soil layer

     real,pointer,dimension(:) :: pptweight       ! lapse weight for precip

     ! The Following variables used to point to structures
     ! Now they point to the array index of the site
     !   type(site), pointer :: hydro_next  !site run-off goes to
     !   type(site), pointer :: hydro_prev  !site run-on comes from
     integer, pointer,dimension(:) :: hydro_next  !site run-off goes to
     integer, pointer,dimension(:) :: hydro_prev  !site run-on comes from
     
     real,pointer,dimension(:)  :: moist_W   !Wetness index
     real,pointer,dimension(:)  :: moist_f   ! decay of soil conductance w/ depth
     real,pointer,dimension(:)  :: moist_tau ! tau characteristic time scale for water redistribution
     real,pointer,dimension(:)  :: moist_zi  ! TOPMODEL "equilibrium" water table depth 
     real,pointer,dimension(:)  :: baseflow  ! loss of water from site to watershed discharge (kg/m2/s)

     integer,pointer,dimension(:,:) :: ntext_soil  ! The soil classifications index of the soil layer
                                                   ! refer to soil_coms to view the soil parameters
                                                   ! associated with each of these classes

     !-----------------------------------
     ! TEMPERATURE VARIABLES
     !-----------------------------------
     ! minimum daily-averaged temperature in this month -- determines whether
     ! or not recruits come up this month.
     real,pointer,dimension(:)  :: min_monthly_temp

     !-----------------------------------
     ! FORESTRY
     !-----------------------------------

     ! is the site under plantation management? (1=yes, 0=no)
     integer,pointer,dimension(:) :: plantation  ! initialized to zero in site creation

     ! Upon creating an agriculture patch in this site, stock it with this 
     ! PFT.  Set, along with the other stocking parameters, in 
     ! init_ed_site_vars().
     integer,pointer,dimension(:) :: agri_stocking_pft
     
     ! Upon creating an agriculture patch in this site, stock it with 
     ! this density of plants [plants/m2]
     real,pointer,dimension(:) :: agri_stocking_density
     
     ! Upon creating a plantation patch in this site, stock it with this PFT
     integer,pointer,dimension(:) :: plantation_stocking_pft
     
     ! Upon creating an plantation patch in this site, stock it with 
     ! this density of plants [plants/m2]
     real,pointer,dimension(:) :: plantation_stocking_density

     ! Unapplied primary forest harvest from previous years (save until 
     ! harvest is above minimum threshold.) [kgC/m2].  Initialized 
     ! together with secondary memory in init_ed_site_vars().
     real,pointer,dimension(:) :: primary_harvest_memory
     
     ! Unapplied secondary forest harvest from previous years (save until 
     ! harvest is above minimum threshold.) [kgC/m2]
     real ,pointer,dimension(:):: secondary_harvest_memory


     !-----------------------------------
     ! FIRE
     !-----------------------------------

     ! site average fire disturbance rate
     real,pointer,dimension(:) :: fire_disturbance_rate
     ! total fuel in the dry patches
     real,pointer,dimension(:) :: ignition_rate

     real,pointer, dimension(:,:) :: lambda_fire ! initialized in create_site !(12,nsites)

     type(prescribed_phen),pointer, dimension(:) :: phen_pars

     !-----------------------------------
     ! DISTURBANCE
     !-----------------------------------
     ! rate of natural disturbance
     real,pointer, dimension(:) :: nat_disturbance_rate

     ! disturbance dist_type id: 
     !        dist_type = 0 if treefall was last disturbance
     !        dist_type = 1 if fire was last disturbance
     integer,pointer, dimension(:) :: nat_dist_type

     ! if new patch is less than min size, store information in the memory
     real,pointer, dimension(:,:,:) :: disturbance_memory !(n_dist_types,n_dist_types,nsites)

     ! the disturbance matrix (to,from)
     real,pointer,dimension(:,:,:) :: disturbance_rates !(n_dist_types,n_dist_types,nsites)

     ! fraction of non-surviving a.g. material removed from patch
     real,pointer,dimension(:,:) :: loss_fraction !(n_dist_types,nsites)

     real,pointer,dimension(:,:):: green_leaf_factor !(n_pft,nsites)
     real,pointer,dimension(:,:) :: leaf_aging_factor !(n_pft,nsites)

     type(met_driv_state),pointer,dimension(:) :: met

     real,pointer, dimension(:,:,:) :: basal_area !(n_pft,n_dbh,nsites), cm2/m2
     real,pointer, dimension(:,:,:) :: agb        !(n_pft,n_dbh,nsites), kgC/m2
     real,pointer, dimension(:,:,:) :: pldens     !(n_pft,n_dbh,nsites)  plant/m2
     real,pointer, dimension(:,:,:) :: bseeds     !(n_pft,n_dbh,nsites)  kgC/m2

     real,pointer, dimension(:,:,:) :: basal_area_growth ! cm2/m2/yr
     real,pointer, dimension(:,:,:) :: basal_area_mort   ! c2/m2/yr
     real,pointer, dimension(:,:,:) :: basal_area_cut    ! c2/m2/yr

     real,pointer, dimension(:,:,:) :: agb_growth        ! kgC/m2/yr
     real,pointer, dimension(:,:,:) :: agb_mort          ! kgC/m2/yr
     real,pointer, dimension(:,:,:) :: agb_cut           ! kgC/m2/yr


     real,pointer,dimension(:) :: cosaoi
     real,pointer,dimension(:) :: avg_albedo_beam
     real,pointer,dimension(:) :: avg_albedo_diffuse
     real,pointer,dimension(:) :: avg_rlong_albedo
     real,pointer,dimension(:) :: avg_albedo
     real,pointer,dimension(:) :: avg_rlongup

     !----- Length of day light, used to average light levels properly. -------------------!
     real, pointer, dimension(:) :: daylight

     ! ------------------------------------------
     ! Fast time flux diagnostic variables == Polygon Level
     !-------------------------------------------
     real,pointer,dimension(:)   :: avg_rshort_gnd ! Total absorbed radiation at the ground
     real,pointer,dimension(:)   :: avg_rlong_gnd  ! Total absorbed radiation at the ground
     real,pointer,dimension(:)   :: avg_carbon_ac  ! Vegetation to Canopy carbon dioxide flux
     
     !----- Moisture ----------------------------
     !                                              | Description
     real,pointer,dimension(:)   :: avg_vapor_lc    ! Leaf to canopy air water flux [kg/m2/s]
     real,pointer,dimension(:)   :: avg_vapor_wc    ! wood to canopy air water flux [kg/m2/s]
     real,pointer,dimension(:)   :: avg_dew_cg      ! Dew to ground flux
     real,pointer,dimension(:)   :: avg_vapor_gc    ! Ground to canopy air latent heat flux [kg/m2/s]
     real,pointer,dimension(:)   :: avg_wshed_vg    ! Water shedding
     real,pointer,dimension(:)   :: avg_intercepted ! Intercepted
     real,pointer,dimension(:)   :: avg_throughfall ! Throughfall
     real,pointer,dimension(:)   :: avg_vapor_ac    ! Canopy to atmosphere water flux [kg/m2/s]
     real,pointer,dimension(:)   :: avg_transp      ! Transpiration
     real,pointer,dimension(:)   :: avg_evap        ! Evaporation
     real,pointer,dimension(:,:) :: avg_smoist_gg   ! Moisture flux between layers
     real,pointer,dimension(:,:) :: avg_transloss   ! Trabspired soil moisture sink
     real,pointer,dimension(:)   :: avg_runoff      ! Total runoff
     real,pointer,dimension(:)   :: avg_drainage    ! Total drainage
     real,pointer,dimension(:)   :: avg_drainage_heat! Total drainage heat flux

     !----- Auxiliary Variables (user can modify to view any variable -----------------------------------------------!
     real,pointer,dimension(:)   :: aux             ! Auxillary surface variable
     real,pointer,dimension(:,:) :: aux_s           ! Auxillary soil variable

     !---- Polygon LAI, WPA, and WAI ------------------------------------------------------!
     real, pointer, dimension(:) :: lai
     real, pointer, dimension(:) :: avg_lma
     real, pointer, dimension(:) :: wpa
     real, pointer, dimension(:) :: wai

     !----- Sensible heat -------------------------------------------------------------------------------------------!
     !                                              | Description
     real,pointer,dimension(:) :: avg_sensible_lc   ! Leaf to Canopy sensible heat flux
     real,pointer,dimension(:) :: avg_sensible_wc   ! wood to Canopy sensible heat flux
     real,pointer,dimension(:) :: avg_qwshed_vg     ! Internal energy of water shedding from leaves
     real,pointer,dimension(:) :: avg_qintercepted  ! Internal energy of intercepted precipitation
     real,pointer,dimension(:) :: avg_qthroughfall  ! Internal energy of throughfall precipitation
     real,pointer,dimension(:) :: avg_sensible_gc   ! Ground to canopy air sensible heat flux
     real,pointer,dimension(:) :: avg_sensible_ac   ! Canopy to atmosphere sensible heat flux
     real,pointer,dimension(:,:) :: avg_sensible_gg ! Net soil heat flux between layers
     real,pointer,dimension(:) :: avg_runoff_heat   ! Total runoff internal energy flux

     !----- Mass and Energy ---------------------------------------------------------------!
     real,pointer,dimension(:) :: avg_leaf_energy ! Average leaf internal energy   [  J/m2]
     real,pointer,dimension(:) :: avg_leaf_temp   ! Average leaf temperature       [     K]
     real,pointer,dimension(:) :: avg_leaf_fliq   ! Avg. liq. frac. of leaf water  [   ---]
     real,pointer,dimension(:) :: avg_leaf_water  ! Avg. water on top of leaves    [ kg/m2]
     real,pointer,dimension(:) :: avg_leaf_hcap   ! Avg. leaf heat capacity        [J/m2/K]
     real,pointer,dimension(:) :: avg_wood_energy ! Average wood internal energy   [  J/m2]
     real,pointer,dimension(:) :: avg_wood_temp   ! Average wood temperature       [     K]
     real,pointer,dimension(:) :: avg_wood_fliq   ! Avg. liq. frac. of wood water  [   ---]
     real,pointer,dimension(:) :: avg_wood_water  ! Avg. water on top of woods     [ kg/m2]
     real,pointer,dimension(:) :: avg_wood_hcap   ! Avg. wood heat capacity        [J/m2/K]
     
     real,pointer,dimension(:) :: avg_can_temp
     real,pointer,dimension(:) :: avg_can_shv
     real,pointer,dimension(:) :: avg_can_co2
     real,pointer,dimension(:) :: avg_can_rhos
     real,pointer,dimension(:) :: avg_can_prss
     real,pointer,dimension(:) :: avg_can_theta
     real,pointer,dimension(:) :: avg_can_theiv
     real,pointer,dimension(:) :: avg_can_depth

     real,pointer,dimension(:,:) :: avg_soil_energy
     real,pointer,dimension(:,:) :: avg_soil_mstpot
     real,pointer,dimension(:,:) :: avg_soil_water
     real,pointer,dimension(:,:) :: avg_soil_temp
     real,pointer,dimension(:,:) :: avg_soil_fracliq
     real,pointer,dimension(:)   :: avg_soil_wetness	!relative to wilting point
     real,pointer,dimension(:,:) :: avg_soil_rootfrac
     
     real,pointer,dimension(:) :: avg_skin_temp
     real,pointer,dimension(:) :: avg_available_water

     real,pointer,dimension(:) :: runoff

     !-----  Phenology
     real, pointer, dimension(:) :: rad_avg

     !----- Meteorological forcing
     real,pointer,dimension(:) :: avg_atm_tmp
     real,pointer,dimension(:) :: avg_atm_shv
     real,pointer,dimension(:) :: avg_atm_prss


     !----- NACP intercomparison ---------------------------------------------!
     real,pointer,dimension(:) :: avg_sfcw_depth
     real,pointer,dimension(:) :: avg_sfcw_energy
     real,pointer,dimension(:) :: avg_sfcw_mass
     real,pointer,dimension(:) :: avg_sfcw_tempk
     real,pointer,dimension(:) :: avg_sfcw_fracliq
     real,pointer,dimension(:) :: avg_fsc
     real,pointer,dimension(:) :: avg_ssc
     real,pointer,dimension(:) :: avg_stsc
     real,pointer,dimension(:) :: avg_balive
     real,pointer,dimension(:) :: avg_bleaf
     real,pointer,dimension(:) :: avg_broot
     real,pointer,dimension(:) :: avg_bsapwood
     real,pointer,dimension(:) :: avg_bdead
     real,pointer,dimension(:) :: avg_fsn
     real,pointer,dimension(:) :: avg_msn
     real,pointer,dimension(:) :: avg_bstorage
     real,pointer,dimension(:) :: avg_bseeds

     ! Daily average residuals
     real, pointer, dimension(:) :: dmean_co2_residual    ! [umol_CO2/m2]
     real, pointer, dimension(:) :: dmean_energy_residual ! [J/m2]
     real, pointer, dimension(:) :: dmean_water_residual  ! [kg_H20/m2]
     ! Monthly average residuals
     real, pointer, dimension(:) :: mmean_co2_residual    ! [umol_CO2/m2]
     real, pointer, dimension(:) :: mmean_energy_residual ! [J/m2]
     real, pointer, dimension(:) :: mmean_water_residual  ! [kg_H20/m2]

  end type polygontype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  !---------------------------------------------------
  ! Ed Type:
  ! The following are the arrays of polygon level variables
  ! that populate the current grid.
  !---------------------------------------------------

  ! About 450 of these.
  type edtype

     !---- Data space indexing variables ---------------
     !
     ! These variables are the total number of each
     ! respective hierarchical type, summed across
     ! all compute nodes if it is a parallel process
     ! These values are mostly needed for output
     ! dataspace dimensioning.

     integer :: npolygons_global
     integer :: nsites_global
     integer :: npatches_global
     integer :: ncohorts_global

     ! Index offsets for total counts of cohorts, patches
     ! and sites.  If this is the nth machine writing to
     ! a file, it needs to know how many of these types
     ! had come before, so it writes data contiguously
     ! into the same file with the others

     integer :: mach_cohort_offset_index
     integer :: mach_patch_offset_index
     integer :: mach_site_offset_index
     integer :: mach_polygon_offset_index

     !
     !---------------------------------------------------


     integer :: npolygons

     integer :: pyglob_id


     !  The global index of the first site in each polygon
     integer,pointer,dimension(:) :: pysi_id

     ! The number of sites in each polygon
     integer,pointer,dimension(:) :: pysi_n

     real(kind=8),pointer,dimension(:) :: walltime_py

     real,pointer,dimension(:) :: lat

     real,pointer,dimension(:) :: lon

     integer,pointer,dimension(:) :: xatm
     
     integer,pointer,dimension(:) :: yatm

     ! Sensible heat flux from the canopy to atmosphere, averaged across sites
   
     real,pointer,dimension(:) :: sensflux_py       ! of dimension npolys

     integer,pointer,dimension(:,:) :: ntext_soil   ! Soil texture classification
                                                    ! (nzg,polygon)
     
     integer,pointer,dimension(:) :: lsl            ! Layer of lowest soil
     
     ! matrix of site hydrologic adjacency
     real,pointer,dimension(:,:,:) :: site_adjacency

     real,pointer,dimension(:) :: wbar
     real,pointer,dimension(:) :: Te
     real,pointer,dimension(:) :: zbar
     real,pointer,dimension(:) :: tau
     real,pointer,dimension(:) :: sheat
     real,pointer,dimension(:) :: baseflow
     real,pointer,dimension(:) :: runoff
     real,pointer,dimension(:) :: swliq

     type(polygontype),pointer,dimension(:) :: polygon

     integer,pointer,dimension(:) :: ilon   ! index for matching met. data
     integer,pointer,dimension(:) :: ilat   ! index for matching met. data

     ! Polygon AGB (kgC/m2)
     real,pointer,dimension(:) :: total_agb

     ! Polygon basal area (cm2/m2)
     real,pointer,dimension(:) :: total_basal_area

     ! AGB accruing due to growth (kgC/m2/yr)
     real,pointer,dimension(:) :: total_agb_growth

     ! AGB lost due to mortality (kgC/m2/yr)
     real,pointer,dimension(:) :: total_agb_mort

     ! AGB used to generate recruits (kgC/m2/yr)
     real,pointer,dimension(:) :: total_agb_recruit



     ! CHANGED THE FOLLOWING UNIT DESCRIPTORS: FROM (cm2/m2/yr) RGK 6-13-08
     !------------
     ! BASAL_AREA accruing due to growth (cm2/m2/yr)
     real,pointer,dimension(:) :: total_basal_area_growth

     ! BASAL_AREA lost due to mortality (cm2/m2/yr)
     real,pointer,dimension(:) :: total_basal_area_mort

     ! BASAL_AREA used to generate recruits (cm2/m2/yr)
     real,pointer,dimension(:) :: total_basal_area_recruit
     !------------


     integer,pointer,dimension(:) :: nsites
     ! list of site numbers
     integer,pointer,dimension(:,:) :: sitenums !(max_site,npolygons)
 
     ! specification if the adjacency table was loaded from a file
     integer,pointer,dimension(:) :: load_adjacency

     !! Meteorological driver data
     type(met_driv_data),pointer,dimension(:) :: metinput
     
     !! Lapse rate transfer data
     type(met_driv_state),pointer,dimension(:) :: met, lapse

     real,pointer,dimension(:) :: cosz
     real,pointer,dimension(:) :: mean_gpp
     real,pointer,dimension(:) :: mean_precip
     real,pointer,dimension(:) :: mean_qprecip
     real,pointer,dimension(:) :: mean_netrad

     ! Total carbon (vegetation plus soil) at the beginning of budget-averaging
     ! time [kgC/m2]
     real,pointer,dimension(:) :: cbudget_initialstorage

     ! Average NEP (GPP - plant respiration - heterotrophic respiration)
     ! [kgC/m2/day], used for evaluating daily carbon budget.
     real,pointer,dimension(:) :: cbudget_nep

     ! Total nitrogen (vegetation plus soil) at the beginning of  
     ! budget-averaging time [kgN/m2]
     real,pointer,dimension(:) :: nbudget_initialstorage

     ! Polygon basal area profile [cm2/m2]
     real,pointer,dimension(:,:,:) :: basal_area !(n_pft,n_dbh,npolygons)

     ! Polygon above-ground biomass [kgC/m2]
     real,pointer,dimension(:,:,:) :: agb  !(n_pft,n_dbh,npolygons)
     ! Polygon plant density [plant/m2]
     real,pointer, dimension(:,:,:) :: pldens     !(n_pft,n_dbh,npolygons)

     ! Seed biomass [kgC/m2]
     real,pointer, dimension(:,:,:) :: bseeds     !(n_pft,n_dbh,npolygons)

     ! ------------------------------------------
     ! Diagnostic variables
     !-------------------------------------------
     
     real,pointer,dimension(:)   :: avg_rshort_gnd ! Total absorbed radiation at the ground
     real,pointer,dimension(:)   :: avg_rlong_gnd  ! Total absorbed radiation at the ground
     real,pointer,dimension(:)   :: avg_carbon_ac  ! Vegetation to Canopy carbon dioxide flux
     !----- Moisture Flux ----------------------------

     real,pointer,dimension(:)   :: avg_vapor_lc    ! Leaf to canopy air water flux
     real,pointer,dimension(:)   :: avg_vapor_wc    ! wood to canopy air water flux
     real,pointer,dimension(:)   :: avg_dew_cg      ! Dew to ground flux
     real,pointer,dimension(:)   :: avg_vapor_gc    ! Ground to canopy air latent heat flux [kg/m2/s]
     real,pointer,dimension(:)   :: avg_wshed_vg    ! Water shedding from the leaves
     real,pointer,dimension(:)   :: avg_intercepted ! Intercepted
     real,pointer,dimension(:)   :: avg_throughfall ! Throughfall
     real,pointer,dimension(:)   :: avg_vapor_ac    ! Canopy to atmosphere water flux [kg/m2/s]
     real,pointer,dimension(:)   :: avg_transp      ! Transpiration
     real,pointer,dimension(:)   :: avg_evap        ! Evaporation
     real,pointer,dimension(:,:) :: avg_smoist_gg   ! Moisture flux between layers
     real,pointer,dimension(:,:) :: avg_transloss   ! Trabspired soil moisture sink
     real,pointer,dimension(:)   :: avg_runoff      ! Total runoff
     real,pointer,dimension(:)   :: avg_drainage      ! Total drainage through the soil bottom
     real,pointer,dimension(:)   :: avg_drainage_heat ! Total drainage internal heat loss

     !----- Auxillary Variables (user can modify to view any variable --------!
     real,pointer,dimension(:)   :: aux             ! Auxillary surface variable
     real,pointer,dimension(:,:) :: aux_s           ! Auxillary soil variable

     !----- LAI, WPA, and WAI ------------------------------------------------!
     real,pointer,dimension(:) :: lai
     real,pointer,dimension(:) :: avg_lma
     real,pointer,dimension(:) :: wpa
     real,pointer,dimension(:) :: wai


     !----- Sensible heat Flux -----------------------------------------------!

     real,pointer,dimension(:)   :: avg_sensible_lc  ! Leaf to Canopy sensible heat flux
     real,pointer,dimension(:)   :: avg_sensible_wc  ! wood to Canopy sensible heat flux
     real,pointer,dimension(:)   :: avg_qwshed_vg    ! Internal energy of leaf water shedding
     real,pointer,dimension(:)   :: avg_qintercepted ! Internal energy of intercepted precipitation
     real,pointer,dimension(:)   :: avg_qthroughfall ! Internal energy of intercepted precipitation
     real,pointer,dimension(:)   :: avg_sensible_gc  ! Ground to canopy air sensible heat flux
     real,pointer,dimension(:)   :: avg_sensible_ac  ! Canopy to atmosphere sensible heat flux
     real,pointer,dimension(:,:) :: avg_sensible_gg  ! Net soil heat flux between layers
     real,pointer,dimension(:)   :: avg_runoff_heat  ! Total runoff internal energy flux

     !----- Mass and Energy --------------------------------------------------!

     ! New variables for testing the model stability
     real,pointer,dimension(:) :: max_leaf_temp
     real,pointer,dimension(:) :: min_leaf_temp
     real,pointer,dimension(:) :: max_wood_temp
     real,pointer,dimension(:) :: min_wood_temp
     real,pointer,dimension(:) :: max_soil_temp
     real,pointer,dimension(:) :: min_soil_temp
     !----- Mass and Energy ---------------------------------------------------------------!
     real,pointer,dimension(:) :: avg_leaf_energy ! Average leaf internal energy   [  J/m2]
     real,pointer,dimension(:) :: avg_leaf_temp   ! Average leaf temperature       [     K]
     real,pointer,dimension(:) :: avg_leaf_fliq   ! Avg. liq. frac. of leaf water  [   ---]
     real,pointer,dimension(:) :: avg_leaf_water  ! Avg. water on top of leaves    [ kg/m2]
     real,pointer,dimension(:) :: avg_leaf_hcap   ! Avg. leaf heat capacity        [J/m2/K]
     real,pointer,dimension(:) :: avg_wood_energy ! Average wood internal energy   [  J/m2]
     real,pointer,dimension(:) :: avg_wood_temp   ! Average wood temperature       [     K]
     real,pointer,dimension(:) :: avg_wood_fliq   ! Avg. liq. frac. of wood water  [   ---]
     real,pointer,dimension(:) :: avg_wood_water  ! Avg. water on top of woods     [ kg/m2]
     real,pointer,dimension(:) :: avg_wood_hcap   ! Avg. wood heat capacity        [J/m2/K]
     
     real,pointer,dimension(:) :: avg_can_temp
     real,pointer,dimension(:) :: avg_can_shv
     real,pointer,dimension(:) :: avg_can_co2
     real,pointer,dimension(:) :: avg_can_rhos
     real,pointer,dimension(:) :: avg_can_prss
     real,pointer,dimension(:) :: avg_can_theta
     real,pointer,dimension(:) :: avg_can_theiv
     real,pointer,dimension(:) :: avg_can_depth

     real,pointer,dimension(:,:) :: avg_soil_energy
     real,pointer,dimension(:,:) :: avg_soil_mstpot
     real,pointer,dimension(:,:) :: avg_soil_water
     real,pointer,dimension(:,:) :: avg_soil_temp
     real,pointer,dimension(:,:) :: avg_soil_fracliq
     real,pointer,dimension(:,:) :: avg_soil_rootfrac
     
     real,pointer,dimension(:) :: avg_soil_wetness
     real,pointer,dimension(:) :: avg_skin_temp
     real,pointer,dimension(:) :: avg_available_water


     real,pointer,dimension(:) :: avg_sfcw_depth     !sfcwater_depth
     real,pointer,dimension(:) :: avg_sfcw_energy    !sfcwater_mass
     real,pointer,dimension(:) :: avg_sfcw_mass      !sfcwater_mass
     real,pointer,dimension(:) :: avg_sfcw_tempk     !sfcwater_tempk
     real,pointer,dimension(:) :: avg_sfcw_fracliq   !sfcwater_fracliq
     real,pointer,dimension(:) :: avg_bdead
     real,pointer,dimension(:) :: avg_balive
     real,pointer,dimension(:) :: avg_bleaf
     real,pointer,dimension(:) :: avg_broot
     real,pointer,dimension(:) :: avg_bsapwood
     real,pointer,dimension(:) :: avg_bseeds
     real,pointer,dimension(:) :: avg_bstorage
     real,pointer,dimension(:) :: avg_fsc
     real,pointer,dimension(:) :: avg_ssc
     real,pointer,dimension(:) :: avg_stsc

     real,pointer,dimension(:) :: avg_fsn
     real,pointer,dimension(:) :: avg_msn

     !-------- TOTAL CARBON AND NITROGEN POOLS  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
     real,pointer,dimension(:) :: Cleaf
     real,pointer,dimension(:) :: Croot
     real,pointer,dimension(:) :: Cstore
     real,pointer,dimension(:) :: Ccwd
     real,pointer,dimension(:) :: Nleaf
     real,pointer,dimension(:) :: Ndead
     real,pointer,dimension(:) :: Nroot
     real,pointer,dimension(:) :: Nstore
     real,pointer,dimension(:) :: Ncwd
     
     
     !-------- TOTAL CARBON AND NITROGEN FLUX  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
     real,pointer,dimension(:) :: Cleaf_grow
     real,pointer,dimension(:) :: Croot_grow
     real,pointer,dimension(:) :: Cdead_grow
     real,pointer,dimension(:) :: Cstore_grow
     real,pointer,dimension(:) :: Cleaf_litter_flux
     real,pointer,dimension(:) :: Croot_litter_flux
     real,pointer,dimension(:) :: Ccwd_flux
     real,pointer,dimension(:) :: Nleaf_grow
     real,pointer,dimension(:) :: Ndead_grow
     real,pointer,dimension(:) :: Nroot_grow
     real,pointer,dimension(:) :: Nstore_grow
     real,pointer,dimension(:) :: Nleaf_litter_flux
     real,pointer,dimension(:) :: Nroot_litter_flux
     real,pointer,dimension(:) :: Ncwd_flux
     real,pointer,dimension(:) :: Nbiomass_uptake
     real,pointer,dimension(:) :: Ngross_min
     real,pointer,dimension(:) :: Nnet_min

     !----- Meteorologic Conditions ----------------------------------------------------!
     real,pointer,dimension(:) :: avg_nir_beam
     real,pointer,dimension(:) :: avg_nir_diffuse
     real,pointer,dimension(:) :: avg_par_beam
     real,pointer,dimension(:) :: avg_par_diffuse
     real,pointer,dimension(:) :: avg_atm_tmp
     real,pointer,dimension(:) :: avg_atm_shv
     real,pointer,dimension(:) :: avg_rshort
     real,pointer,dimension(:) :: avg_rshort_diffuse
     real,pointer,dimension(:) :: avg_rlong
     real,pointer,dimension(:) :: avg_pcpg
     real,pointer,dimension(:) :: avg_qpcpg
     real,pointer,dimension(:) :: avg_dpcpg
     real,pointer,dimension(:) :: avg_vels
     real,pointer,dimension(:) :: avg_atm_prss
     real,pointer,dimension(:) :: avg_exner
     real,pointer,dimension(:) :: avg_geoht
     real,pointer,dimension(:) :: avg_atm_co2
     real,pointer,dimension(:) :: avg_albedo
     real,pointer,dimension(:) :: avg_albedo_beam
     real,pointer,dimension(:) :: avg_albedo_diffuse
     real,pointer,dimension(:) :: avg_rlong_albedo
     real,pointer,dimension(:) :: avg_rlongup
     real,pointer,dimension(:,:,:) :: avg_lai_ebalvars   ! This diagnostic partitions energy flux
                                                         ! variables into LAI regimes.
                                                         ! The matrix is arranged as follows
                                                         ! (LAI,VARIABLE,POLYGON)
                                                         ! Where LAI is (<2.0,2-4,>4)
                                                         ! Where Variable is (RNET,LHF,SHF,CANTEMP)
     real,pointer,dimension(:) :: avg_gpp
     real,pointer,dimension(:) :: avg_leaf_resp
     real,pointer,dimension(:) :: avg_root_resp
     real,pointer,dimension(:) :: avg_growth_resp
     real,pointer,dimension(:) :: avg_storage_resp
     real,pointer,dimension(:) :: avg_vleaf_resp
     real,pointer,dimension(:) :: avg_plant_resp
     real,pointer,dimension(:) :: avg_htroph_resp
     real,pointer,dimension(:) :: avg_leaf_drop
     real,pointer,dimension(:) :: avg_leaf_maintenance
     real,pointer,dimension(:) :: avg_root_maintenance
!     real,pointer,dimension(:) :: avg_nppleaf
!     real,pointer,dimension(:) :: avg_nppfroot
!     real,pointer,dimension(:) :: avg_nppsapwood
!     real,pointer,dimension(:) :: avg_nppcroot
!     real,pointer,dimension(:) :: avg_nppseeds
!     real,pointer,dimension(:) :: avg_nppwood
!     real,pointer,dimension(:) :: avg_nppdaily


     !-------------------------------------------------------------------!
     ! These variables carry the daily mean, and are allocated only when !
     ! daily or monthly means are requested by the user. For now only    !
     ! polygon-level averages are available, and these are (almost) the  !
     ! same as the old structure, now only at the polygon level. This is !
     ! not a requirement at all, if you feel like looking at daily means !
     ! of site/patch/cohort level variables, feel free to include them.  !
     !-------------------------------------------------------------------!

     real, pointer, dimension(:)   :: dmean_pcpg           ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_runoff         ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_drainage       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_evap           ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_transp         ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_vapor_lc       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_vapor_wc       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_vapor_gc       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_vapor_ac       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: dmean_sensible_lc    ! [      W/m2]
     real, pointer, dimension(:)   :: dmean_sensible_wc    ! [      W/m2]
     real, pointer, dimension(:)   :: dmean_sensible_gc    ! [      W/m2]
     real, pointer, dimension(:)   :: dmean_sensible_ac    ! [      W/m2]
     real, pointer, dimension(:)   :: dmean_nee            ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_gpp            ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_nep            ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_plresp         ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_rh             ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_leaf_resp      ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_root_resp      ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_growth_resp    ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_storage_resp   ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: dmean_vleaf_resp     ! [ kgC/m²/yr]
     
!     real, pointer, dimension(:)   :: dmean_nppleaf
!     real, pointer, dimension(:)   :: dmean_nppfroot
!     real, pointer, dimension(:)   :: dmean_nppsapwood
!     real, pointer, dimension(:)   :: dmean_nppcroot
!     real, pointer, dimension(:)   :: dmean_nppseeds
!     real, pointer, dimension(:)   :: dmean_nppwood
!     real, pointer, dimension(:)   :: dmean_nppdaily
     
     real, pointer, dimension(:,:) :: dmean_soil_temp      !(nzg,npolygons)
     real, pointer, dimension(:,:) :: dmean_soil_water     !(nzg,npolygons)
     real, pointer, dimension(:,:) :: dmean_soil_mstpot    !(nzg,npolygons)
     real, pointer, dimension(:,:) :: dmean_transloss     !(nzg,npolygons)
     real, pointer, dimension(:)   :: dmean_fs_open
     real, pointer, dimension(:)   :: dmean_fsw
     real, pointer, dimension(:)   :: dmean_fsn
     real, pointer, dimension(:)   :: dmean_can_temp      ! (npolygons)
     real, pointer, dimension(:)   :: dmean_can_shv       ! (npolygons)
     real, pointer, dimension(:)   :: dmean_can_co2       ! (npolygons)
     real, pointer, dimension(:)   :: dmean_can_rhos      ! (npolygons)
     real, pointer, dimension(:)   :: dmean_can_prss      ! (npolygons)
     real, pointer, dimension(:)   :: dmean_can_theta     ! (npolygons)
     real, pointer, dimension(:)   :: dmean_can_theiv     ! (npolygons)
     real, pointer, dimension(:)   :: dmean_gnd_temp      ! (npolygons)
     real, pointer, dimension(:)   :: dmean_gnd_shv       ! (npolygons)
     real, pointer, dimension(:)   :: dmean_leaf_energy   ! (npolygons)
     real, pointer, dimension(:)   :: dmean_leaf_water    ! (npolygons)
     real, pointer, dimension(:)   :: dmean_leaf_hcap     ! (npolygons)
     real, pointer, dimension(:)   :: dmean_leaf_temp     ! (npolygons)
     real, pointer, dimension(:)   :: dmean_wood_energy   ! (npolygons)
     real, pointer, dimension(:)   :: dmean_wood_water    ! (npolygons)
     real, pointer, dimension(:)   :: dmean_wood_hcap     ! (npolygons)
     real, pointer, dimension(:)   :: dmean_wood_temp     ! (npolygons)
     real, pointer, dimension(:)   :: dmean_atm_temp      ! (npolygons)
     real, pointer, dimension(:)   :: dmean_atm_shv       ! (npolygons)
     real, pointer, dimension(:)   :: dmean_atm_co2       ! (npolygons)
     real, pointer, dimension(:)   :: dmean_atm_prss      ! (npolygons)
     real, pointer, dimension(:)   :: dmean_atm_vels      ! (npolygons)
     real, pointer, dimension(:)   :: dmean_rshort        ! (npolygons)
     real, pointer, dimension(:)   :: dmean_rshort_diff   ! (npolygons)
     real, pointer, dimension(:)   :: dmean_rlong         ! (npolygons)
     real, pointer, dimension(:)   :: dmean_rshort_gnd    ! (npolygons)
     real, pointer, dimension(:)   :: dmean_rlong_gnd     ! (npolygons)
     real, pointer, dimension(:)   :: dmean_albedo        ! (npolygons)
     real, pointer, dimension(:)   :: dmean_albedo_beam   ! (npolygons)
     real, pointer, dimension(:)   :: dmean_albedo_diffuse! (npolygons)
     real, pointer, dimension(:)   :: dmean_rlong_albedo  ! (npolygons)
     real, pointer, dimension(:)   :: dmean_rlongup       ! (npolygons)

     real, pointer, dimension(:,:)   :: qmean_pcpg           ! [   kg/m²/s]
     real, pointer, dimension(:,:)   :: qmean_runoff         ! [   kg/m²/s]
     real, pointer, dimension(:,:)   :: qmean_drainage       ! [   kg/m²/s]
     real, pointer, dimension(:,:)   :: qmean_evap           ! [   kg/m²/s]
     real, pointer, dimension(:,:)   :: qmean_transp         ! [   kg/m²/s]
     real, pointer, dimension(:,:)   :: qmean_vapor_lc       ! [   kg/m²/s]
     real, pointer, dimension(:,:)   :: qmean_vapor_wc       ! [   kg/m²/s]
     real, pointer, dimension(:,:)   :: qmean_vapor_gc       ! [   kg/m²/s]
     real, pointer, dimension(:,:)   :: qmean_vapor_ac       ! [   kg/m²/s]
     real, pointer, dimension(:,:)   :: qmean_sensible_lc    ! [      W/m2]
     real, pointer, dimension(:,:)   :: qmean_sensible_wc    ! [      W/m2]
     real, pointer, dimension(:,:)   :: qmean_sensible_gc    ! [      W/m2]
     real, pointer, dimension(:,:)   :: qmean_sensible_ac    ! [      W/m2]
     real, pointer, dimension(:,:)   :: qmean_nee            ! [ kgC/m²/yr]
     real, pointer, dimension(:,:)   :: qmean_gpp            ! [ kgC/m²/yr]
     real, pointer, dimension(:,:)   :: qmean_nep            ! [ kgC/m²/yr]
     real, pointer, dimension(:,:)   :: qmean_plresp         ! [ kgC/m²/yr]
     real, pointer, dimension(:,:)   :: qmean_rh             ! [ kgC/m²/yr]
     real, pointer, dimension(:,:)   :: qmean_leaf_resp      ! [ kgC/m²/yr]
     real, pointer, dimension(:,:)   :: qmean_root_resp      ! [ kgC/m²/yr]

     real, pointer, dimension(:,:,:) :: qmean_soil_temp      ! (nzg,ndcycle,npolygons)
     real, pointer, dimension(:,:,:) :: qmean_soil_water     ! (nzg,ndcycle,npolygons)
     real, pointer, dimension(:,:,:) :: qmean_soil_mstpot    ! (nzg,ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_fs_open        ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_fsw            ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_fsn            ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_can_temp       ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_can_shv        ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_can_co2        ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_can_rhos       ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_can_prss       ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_can_theta      ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_can_theiv      ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_gnd_temp       ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_gnd_shv        ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_leaf_energy    ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_leaf_water     ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_leaf_hcap      ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_leaf_temp      ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_wood_energy    ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_wood_water     ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_wood_hcap      ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_wood_temp      ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_atm_temp       ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_atm_shv        ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_atm_co2        ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_atm_prss       ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_atm_vels       ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_rshort         ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_rshort_diff    ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_rlong          ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_rshort_gnd     ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_rlong_gnd      ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_albedo         ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_albedo_beam    ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_albedo_diffuse ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_rlong_albedo   ! (ndcycle,npolygons)
     real, pointer, dimension(:,:)   :: qmean_rlongup        ! (ndcycle,npolygons)

     !-------------------------------------------------------------------------------------!
     ! These variables have the monthly mean square,  We save this instead of the standard !
     ! deviation because we can aggregate several months after the output.                 !
     !-------------------------------------------------------------------------------------!
     real, pointer, dimension(:,:) :: qmsqu_gpp
     real, pointer, dimension(:,:) :: qmsqu_leaf_resp
     real, pointer, dimension(:,:) :: qmsqu_root_resp
     real, pointer, dimension(:,:) :: qmsqu_plresp
     real, pointer, dimension(:,:) :: qmsqu_nee
     real, pointer, dimension(:,:) :: qmsqu_nep
     real, pointer, dimension(:,:) :: qmsqu_rh
     real, pointer, dimension(:,:) :: qmsqu_sensible_ac
     real, pointer, dimension(:,:) :: qmsqu_sensible_lc
     real, pointer, dimension(:,:) :: qmsqu_sensible_wc
     real, pointer, dimension(:,:) :: qmsqu_sensible_gc
     real, pointer, dimension(:,:) :: qmsqu_evap
     real, pointer, dimension(:,:) :: qmsqu_transp
     real, pointer, dimension(:,:) :: qmsqu_vapor_ac
     real, pointer, dimension(:,:) :: qmsqu_vapor_lc
     real, pointer, dimension(:,:) :: qmsqu_vapor_wc
     real, pointer, dimension(:,:) :: qmsqu_vapor_gc

     ! Daily average residuals
     real, pointer, dimension(:) :: dmean_co2_residual    ! [umol_CO2/m2]
     real, pointer, dimension(:) :: dmean_energy_residual ! [J/m2]
     real, pointer, dimension(:) :: dmean_water_residual  ! [kg_H20/m2]
     ! Monthly average residuals
     real, pointer, dimension(:) :: mmean_co2_residual    ! [umol_CO2/m2]
     real, pointer, dimension(:) :: mmean_energy_residual ! [J/m2]
     real, pointer, dimension(:) :: mmean_water_residual  ! [kg_H20/m2]

     !-------------------------------------------------------------------!
     !   These are variables updated at a daily basis, so they are not   !
     ! averages but they are written at the daily analysis               !
     !-------------------------------------------------------------------!
     real, pointer, dimension(:,:) :: lai_pft    ! (n_pft       , npolygons)
     real, pointer, dimension(:,:) :: wpa_pft    ! (n_pft       , npolygons)
     real, pointer, dimension(:,:) :: wai_pft    ! (n_pft       , npolygons)
     real, pointer, dimension(:,:) :: dmean_gpp_dbh    !(n_dbh       ,npolygons)


     !-------------------------------------------------------------------!
     !   These are variables updated at a monthly basis, so they are not !
     ! averages but they are written at the monthly analysis             !
     !-------------------------------------------------------------------!
     real, pointer, dimension(:,:) :: agb_pft    !(n_pft       ,npolygons), kgC/m2
     real, pointer, dimension(:,:) :: ba_pft     !(n_pft       ,npolygons), cm2/m2
     real, pointer, dimension(:,:) :: bseeds_pft !(n_pft      ,npolygons), kgC/m2

     real, pointer, dimension(:,:) :: mmean_gpp_dbh  !(n_dbh       ,npolygons)
     real, pointer, dimension(:,:) :: mmean_lai_pft  !(n_pft       ,npolygons)
     real, pointer, dimension(:,:) :: mmean_wpa_pft  !(n_pft       ,npolygons)
     real, pointer, dimension(:,:) :: mmean_wai_pft  !(n_pft       ,npolygons)

     !-------------------------------------------------------------------!
     ! These variables carry the montlhly mean, and are allocated only   !
     ! when monthly means are requested by the user. For now only        !
     ! polygon-level averages are available, and these are (almost) the  !
     ! same as the old structure, now only at the polygon level. This is !
     ! not a requirement at all, if you feel like looking at daily means !
     ! of site/patch/cohort level variables, feel free to include them.  !
     !-------------------------------------------------------------------!
     real, pointer, dimension(:)   :: mmean_pcpg           ! [   kg/m²/s]
     real, pointer, dimension(:)   :: mmean_runoff         ! [   kg/m²/s]
     real, pointer, dimension(:)   :: mmean_drainage       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: mmean_evap           ! [   kg/m²/s]
     real, pointer, dimension(:)   :: mmean_transp         ! [   kg/m²/s]
     real, pointer, dimension(:)   :: mmean_vapor_lc       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: mmean_vapor_wc       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: mmean_vapor_gc       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: mmean_vapor_ac       ! [   kg/m²/s]
     real, pointer, dimension(:)   :: mmean_sensible_lc    ! [      W/m²]
     real, pointer, dimension(:)   :: mmean_sensible_wc    ! [      W/m²]
     real, pointer, dimension(:)   :: mmean_sensible_gc    ! [      W/m²]
     real, pointer, dimension(:)   :: mmean_sensible_ac    ! [      W/m²]
     real, pointer, dimension(:)   :: mmean_nee            ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: mmean_gpp            ! [ kgC/m²/yr]
!     real, pointer, dimension(:)   :: mmean_nppleaf        ! [ kgC/m²/yr]
!     real, pointer, dimension(:)   :: mmean_nppfroot       ! [ kgC/m²/yr]
!     real, pointer, dimension(:)   :: mmean_nppsapwood     ! [ kgC/m²/yr]
!     real, pointer, dimension(:)   :: mmean_nppcroot       ! [ kgC/m²/yr]
!     real, pointer, dimension(:)   :: mmean_nppseeds       ! [ kgC/m²/yr]
!     real, pointer, dimension(:)   :: mmean_nppwood        ! [ kgC/m²/yr]
!     real, pointer, dimension(:)   :: mmean_nppdaily       ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: mmean_nep            ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: mmean_plresp         ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: mmean_rh             ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: mmean_leaf_resp      ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: mmean_root_resp      ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: mmean_growth_resp    ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: mmean_storage_resp   ! [ kgC/m²/yr]
     real, pointer, dimension(:)   :: mmean_vleaf_resp     ! [ kgC/m²/yr]
     real, pointer, dimension(:,:) :: mmean_soil_temp      !(nzg,npolygons)
     real, pointer, dimension(:,:) :: mmean_soil_water     !(nzg,npolygons)
     real, pointer, dimension(:,:) :: mmean_soil_mstpot    !(nzg,npolygons)
     real, pointer, dimension(:,:) :: mmean_transloss     !(nzg,npolygons)
     real, pointer, dimension(:)   :: mmean_fs_open
     real, pointer, dimension(:)   :: mmean_fsw
     real, pointer, dimension(:)   :: mmean_fsn

     real, pointer, dimension(:)   :: mmean_can_temp       ! (npolygons)
     real, pointer, dimension(:)   :: mmean_can_shv        ! (npolygons)
     real, pointer, dimension(:)   :: mmean_can_co2        ! (npolygons)
     real, pointer, dimension(:)   :: mmean_can_rhos       ! (npolygons)
     real, pointer, dimension(:)   :: mmean_can_prss       ! (npolygons)
     real, pointer, dimension(:)   :: mmean_can_theta      ! (npolygons)
     real, pointer, dimension(:)   :: mmean_can_theiv      ! (npolygons)
     real, pointer, dimension(:)   :: mmean_gnd_temp       ! (npolygons)
     real, pointer, dimension(:)   :: mmean_gnd_shv        ! (npolygons)
     real, pointer, dimension(:)   :: mmean_leaf_energy    ! (npolygons)
     real, pointer, dimension(:)   :: mmean_leaf_water     ! (npolygons)
     real, pointer, dimension(:)   :: mmean_leaf_temp      ! (npolygons)
     real, pointer, dimension(:)   :: mmean_leaf_hcap      ! (npolygons)
     real, pointer, dimension(:)   :: mmean_wood_energy    ! (npolygons)
     real, pointer, dimension(:)   :: mmean_wood_water     ! (npolygons)
     real, pointer, dimension(:)   :: mmean_wood_temp      ! (npolygons)
     real, pointer, dimension(:)   :: mmean_wood_hcap      ! (npolygons)
     real, pointer, dimension(:)   :: mmean_atm_temp       ! (npolygons)
     real, pointer, dimension(:)   :: mmean_rshort         ! (npolygons)
     real, pointer, dimension(:)   :: mmean_rshort_diff    ! (npolygons)
     real, pointer, dimension(:)   :: mmean_rlong          ! (npolygons)
     real, pointer, dimension(:)   :: mmean_rshort_gnd     ! (npolygons)
     real, pointer, dimension(:)   :: mmean_rlong_gnd      ! (npolygons)
     real, pointer, dimension(:)   :: mmean_atm_shv        ! (npolygons)
     real, pointer, dimension(:)   :: mmean_atm_co2        ! (npolygons)
     real, pointer, dimension(:)   :: mmean_atm_prss       ! (npolygons)
     real, pointer, dimension(:)   :: mmean_atm_vels       ! (npolygons)
     real, pointer, dimension(:)   :: mmean_albedo         ! (npolygons)
     real, pointer, dimension(:)   :: mmean_albedo_beam    ! (npolygons)
     real, pointer, dimension(:)   :: mmean_albedo_diffuse ! (npolygons)
     real, pointer, dimension(:)   :: mmean_rlong_albedo   ! (npolygons)
     real, pointer, dimension(:)   :: mmean_rlongup        ! (npolygons)

     !-------------------------------------------------------------------------------------!
     ! These variables have the monthly mean square,  We save this instead of the standard !
     ! deviation because we can aggregate several months afterwards.                       !
     !-------------------------------------------------------------------------------------!
     real, pointer, dimension(:) :: mmsqu_gpp
     real, pointer, dimension(:) :: mmsqu_leaf_resp
     real, pointer, dimension(:) :: mmsqu_root_resp
     real, pointer, dimension(:) :: mmsqu_plresp
     real, pointer, dimension(:) :: mmsqu_nee
     real, pointer, dimension(:) :: mmsqu_nep
     real, pointer, dimension(:) :: mmsqu_rh
     real, pointer, dimension(:) :: mmsqu_sensible_ac
     real, pointer, dimension(:) :: mmsqu_sensible_lc
     real, pointer, dimension(:) :: mmsqu_sensible_wc
     real, pointer, dimension(:) :: mmsqu_sensible_gc
     real, pointer, dimension(:) :: mmsqu_evap
     real, pointer, dimension(:) :: mmsqu_transp
     real, pointer, dimension(:) :: mmsqu_vapor_ac
     real, pointer, dimension(:) :: mmsqu_vapor_lc
     real, pointer, dimension(:) :: mmsqu_vapor_wc
     real, pointer, dimension(:) :: mmsqu_vapor_gc


     !----- Disturbance rates. ----------------------------------------------!
     real, pointer, dimension(:,:,:) :: disturbance_rates

     !-----------------------------------------------------------------------!
     !     This variable storages the workload of this polygon, and this is  !
     ! "measured" by the total number of Runge-Kutta time steps this polygon !
     ! takes on average every day.  This is not averaged by patch, instead   ! 
     ! this is added (because having two patches is computationally more     !
     ! expensive than having one).  Since the number of patches changes by   !
     ! season, and this may be dramatically different, we store the previous !
     ! 12 months.                                                            !
     !-----------------------------------------------------------------------!
     real, pointer, dimension(:,:) :: workload

  end type edtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  !----------------------------------------------!  
  !  GLOBAL VARIABLES



  !------------------------------------------------------------------------------------!
  ! These variables are allocated and assigned before the parallel distribution, so we !
  ! don't want to keep it inside the structure. In case of serial runs, we should have !
  ! offset full of zeroes and the polygons numbered from 1 to the number of polygons.  !
  ! In a non-POI parallel run, we want to allocate "mpolygons" in each node, but then  !
  ! we need to keep track of the actual polygon "ID" to write the output correctly.    !
  !------------------------------------------------------------------------------------!
  
  ! Number of polygons in each grid, for each machine. so this has a ngrids size. 
  integer, dimension(maxmach,maxgrds) :: gdpy

  ! Number of sites in each grid, for each machine.
  integer, dimension(maxmach,maxgrds) :: gdsi

  ! Number of patches in each grid, for each machine.
  integer, dimension(maxmach,maxgrds) :: gdpa

  ! Number of cohorts in each grid, for each machine.
  integer, dimension(maxmach,maxgrds) :: gdco
  

  ! Offset for each machine, so this has a nmachs size.
  integer, dimension(maxmach,maxgrds) :: py_off 
  
  ! Offset for each machine, so this has a nmachs size.
  integer, dimension(maxmach,maxgrds) :: si_off 

  ! Offset for each machine, so this has a nmachs size.
  integer, dimension(maxmach,maxgrds) :: pa_off 

  ! Offset for each machine, so this has a nmachs size.
  integer, dimension(maxmach,maxgrds) :: co_off 


  type(edtype),pointer,dimension(:) :: edgrid_g

  ! The following are swap variables used during
  ! deallocation-reallocation

  type(edtype)      :: edswap_g

  type(polygontype) :: polyswap_g

  type(sitetype)    :: siteswap_g

  type(patchtype)   :: patchswap_g
  
!------------------------------------------------------------------------------------------!
!   The following variables are for tracking the number of variables written in the output,!
! this way we avoid having the same ID used twice.                                         !
!------------------------------------------------------------------------------------------!
  integer :: nioglobal, niogrid, niopoly, niosite

!------------------------------------------------------------------------------------------!
! Logical switch that decides if the pointer tables for IO need to be updated
! The number and allocation of cohorts and patches dictates this, and they 
! change at a monthly frequency typically.
!------------------------------------------------------------------------------------------!
  logical :: filltables

  
  
contains
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  ! ===================================================
  !  Allocation subroutines of ed state variables
  ! ===================================================

  subroutine allocate_edglobals(ngrids)

    use ed_var_tables,only:num_var,vt_info,maxvars

    implicit none
    integer :: ngrids
    
    if (associated(edgrid_g))  then
       print*,"SHOULD NOT HAVE ASSOCIATED GLOBALS"
    else
       nullify(edgrid_g)
       allocate(edgrid_g(ngrids))
    end if

    !  Allocate the basic var-table structures

    allocate(num_var(ngrids))
    num_var = 0
    
    allocate(vt_info(maxvars,ngrids))
    vt_info(:,:)%first=.true.


    !  Initialize the global offsets

    edgrid_g(:)%mach_cohort_offset_index = 0
    edgrid_g(:)%mach_patch_offset_index = 0
    edgrid_g(:)%mach_site_offset_index = 0
    edgrid_g(:)%mach_polygon_offset_index = 0


    return
  end subroutine allocate_edglobals
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine allocate_edtype(cgrid,npolygons)
    
    implicit none

    integer :: npolygons
    type(edtype),target :: cgrid

    call nullify_edtype(cgrid)
    cgrid%npolygons = npolygons
    ! This if is needed for coupled runs: some nodes may receive areas exclusively over oceans, and 
    ! npolygons will be 0 in these nodes. Nothing should be allocated in these circumstances. 
    if (npolygons > 0) then
       allocate(cgrid%polygon(npolygons))
       allocate(cgrid%lat(npolygons))
       allocate(cgrid%lon(npolygons))
       allocate(cgrid%xatm(npolygons))
       allocate(cgrid%yatm(npolygons))
       allocate(cgrid%ntext_soil(nzg,npolygons))
       allocate(cgrid%lsl(npolygons))
       
       allocate(cgrid%pysi_id(npolygons))
       allocate(cgrid%pysi_n(npolygons))
       allocate(cgrid%walltime_py(npolygons))
       allocate(cgrid%sensflux_py(npolygons))
       allocate(cgrid%site_adjacency(max_site,(max_site+1),npolygons))
       allocate(cgrid%wbar(npolygons))
       allocate(cgrid%Te(npolygons))
       allocate(cgrid%zbar(npolygons))
!! NOT USED       allocate(cgrid%tau(npolygons))
       allocate(cgrid%sheat(npolygons))
       allocate(cgrid%baseflow(npolygons))
       allocate(cgrid%runoff(npolygons))
       allocate(cgrid%swliq(npolygons))
       
       allocate(cgrid%ilon(npolygons))
       allocate(cgrid%ilat(npolygons))
       allocate(cgrid%total_agb(npolygons))
       allocate(cgrid%total_basal_area(npolygons))
       allocate(cgrid%total_agb_growth(npolygons))
       allocate(cgrid%total_agb_mort(npolygons))
       allocate(cgrid%total_agb_recruit(npolygons))
       allocate(cgrid%total_basal_area_growth(npolygons))
       allocate(cgrid%total_basal_area_mort(npolygons))
       allocate(cgrid%total_basal_area_recruit(npolygons))
       allocate(cgrid%nsites(npolygons))
       allocate(cgrid%sitenums(max_site,npolygons))
       allocate(cgrid%load_adjacency(npolygons))
       allocate(cgrid%cosz(npolygons))
       allocate(cgrid%mean_gpp(npolygons))
       allocate(cgrid%mean_precip(npolygons))
       allocate(cgrid%mean_qprecip(npolygons))
       allocate(cgrid%mean_netrad(npolygons))
       allocate(cgrid%cbudget_initialstorage(npolygons))
       allocate(cgrid%cbudget_nep(npolygons))
       allocate(cgrid%nbudget_initialstorage(npolygons))
       allocate(cgrid%basal_area(n_pft,n_dbh,npolygons))
       allocate(cgrid%agb(n_pft,n_dbh,npolygons))
       allocate(cgrid%pldens(n_pft,n_dbh,npolygons))
       allocate(cgrid%bseeds(n_pft,n_dbh,npolygons))
       
       allocate(cgrid%metinput(npolygons))
       allocate(cgrid%met(npolygons))

       allocate(cgrid%lapse(npolygons))

       allocate(cgrid%lai  (npolygons))
       allocate(cgrid%avg_lma(npolygons))
       allocate(cgrid%wpa  (npolygons))
       allocate(cgrid%wai  (npolygons))

       ! Fast time flux diagnostics
       ! ---------------------------------------------
       allocate(cgrid%avg_vapor_lc  (npolygons))
       allocate(cgrid%avg_vapor_wc  (npolygons))
       allocate(cgrid%avg_dew_cg    (npolygons))
       allocate(cgrid%avg_vapor_gc  (npolygons))
       allocate(cgrid%avg_wshed_vg  (npolygons))
       allocate(cgrid%avg_intercepted (npolygons))
       allocate(cgrid%avg_throughfall (npolygons))
       allocate(cgrid%avg_vapor_ac  (npolygons))
       allocate(cgrid%avg_transp    (npolygons))
       allocate(cgrid%avg_evap      (npolygons))
       allocate(cgrid%avg_smoist_gg (nzg,npolygons))
       allocate(cgrid%avg_transloss (nzg,npolygons))
       allocate(cgrid%avg_runoff        (npolygons))
       allocate(cgrid%avg_drainage       (npolygons))
       allocate(cgrid%avg_drainage_heat  (npolygons))
       allocate(cgrid%aux           (npolygons))
       allocate(cgrid%aux_s         (nzg,npolygons))
       allocate(cgrid%avg_rshort_gnd   (npolygons))
       allocate(cgrid%avg_rlong_gnd    (npolygons))
       allocate(cgrid%avg_carbon_ac    (npolygons))
       allocate(cgrid%avg_qwshed_vg    (npolygons))
       allocate(cgrid%avg_qintercepted (npolygons))
       allocate(cgrid%avg_qthroughfall (npolygons))
       allocate(cgrid%avg_sensible_lc  (npolygons))
       allocate(cgrid%avg_sensible_wc  (npolygons))
       allocate(cgrid%avg_sensible_gc  (npolygons))
       allocate(cgrid%avg_sensible_ac  (npolygons))
       allocate(cgrid%avg_sensible_gg  (nzg,npolygons))
       allocate(cgrid%avg_runoff_heat  (npolygons))

       allocate(cgrid%max_leaf_temp(npolygons))
       allocate(cgrid%min_leaf_temp(npolygons))
       allocate(cgrid%max_wood_temp(npolygons))
       allocate(cgrid%min_wood_temp(npolygons))
       allocate(cgrid%max_soil_temp(npolygons))
       allocate(cgrid%min_soil_temp(npolygons))

       ! Fast time state diagnostics
       allocate(cgrid%avg_leaf_energy  (npolygons))
       allocate(cgrid%avg_leaf_hcap    (npolygons))
       allocate(cgrid%avg_leaf_temp    (npolygons))
       allocate(cgrid%avg_leaf_fliq    (npolygons))
       allocate(cgrid%avg_leaf_water   (npolygons))
       allocate(cgrid%avg_wood_energy  (npolygons))
       allocate(cgrid%avg_wood_hcap    (npolygons))
       allocate(cgrid%avg_wood_temp    (npolygons))
       allocate(cgrid%avg_wood_fliq    (npolygons))
       allocate(cgrid%avg_wood_water   (npolygons))
       allocate(cgrid%avg_can_temp    (npolygons))
       allocate(cgrid%avg_can_shv     (npolygons))
       allocate(cgrid%avg_can_co2     (npolygons))
       allocate(cgrid%avg_can_rhos    (npolygons))
       allocate(cgrid%avg_can_prss    (npolygons))
       allocate(cgrid%avg_can_theta   (npolygons))
       allocate(cgrid%avg_can_theiv(npolygons))
       allocate(cgrid%avg_can_depth   (npolygons))
       allocate(cgrid%avg_soil_energy(nzg,npolygons))
       allocate(cgrid%avg_soil_mstpot(nzg,npolygons))
       allocate(cgrid%avg_soil_water(nzg,npolygons))
       allocate(cgrid%avg_soil_temp (nzg,npolygons))
       allocate(cgrid%avg_soil_fracliq (nzg,npolygons))
       allocate(cgrid%avg_soil_rootfrac(nzg,npolygons))

       allocate(cgrid%avg_soil_wetness   (npolygons))
       allocate(cgrid%avg_skin_temp      (npolygons))
       allocate(cgrid%avg_available_water(npolygons))

       allocate(cgrid%avg_lai_ebalvars (3,4,npolygons))

       allocate(cgrid%avg_gpp  (npolygons))
       allocate(cgrid%avg_leaf_resp  (npolygons))
       allocate(cgrid%avg_root_resp  (npolygons))
       allocate(cgrid%avg_growth_resp  (npolygons))
       allocate(cgrid%avg_storage_resp  (npolygons))
       allocate(cgrid%avg_vleaf_resp  (npolygons))
       allocate(cgrid%avg_plant_resp  (npolygons))
       allocate(cgrid%avg_growth_resp  (npolygons))
       allocate(cgrid%avg_storage_resp  (npolygons))
       allocate(cgrid%avg_vleaf_resp  (npolygons))
       allocate(cgrid%avg_htroph_resp  (npolygons))
       allocate(cgrid%avg_leaf_drop       (npolygons))
       allocate(cgrid%avg_leaf_maintenance(npolygons))
       allocate(cgrid%avg_root_maintenance(npolygons))

       !! added MCD for NACP intercomparison
       allocate(cgrid%avg_sfcw_depth   (npolygons))
       allocate(cgrid%avg_sfcw_energy  (npolygons))
       allocate(cgrid%avg_sfcw_mass    (npolygons))
       allocate(cgrid%avg_sfcw_tempk   (npolygons))
       allocate(cgrid%avg_sfcw_fracliq (npolygons))
       allocate(cgrid%avg_bdead       (npolygons))
       allocate(cgrid%avg_balive      (npolygons))
       allocate(cgrid%avg_bleaf       (npolygons))
       allocate(cgrid%avg_broot       (npolygons))
       allocate(cgrid%avg_bsapwood    (npolygons))
       allocate(cgrid%avg_bstorage    (npolygons))
       allocate(cgrid%avg_bseeds      (npolygons))
       allocate(cgrid%avg_fsc         (npolygons))
       allocate(cgrid%avg_ssc         (npolygons))
       allocate(cgrid%avg_stsc        (npolygons))
       allocate(cgrid%avg_fsn         (npolygons))
       allocate(cgrid%avg_msn         (npolygons))

       !! C & N fluxes and pools for NCEAS/FACE
       allocate(cgrid%Cleaf  (npolygons))
     allocate(cgrid%Croot  (npolygons))
     allocate(cgrid%Cstore (npolygons))
     allocate(cgrid%Ccwd   (npolygons))
     allocate(cgrid%Nleaf  (npolygons))
     allocate(cgrid%Ndead  (npolygons))
     allocate(cgrid%Nroot  (npolygons))
     allocate(cgrid%Nstore (npolygons))
     allocate(cgrid%Ncwd   (npolygons))
     allocate(cgrid%Cleaf_grow       (npolygons))
     allocate(cgrid%Croot_grow       (npolygons))
     allocate(cgrid%Cdead_grow       (npolygons))
     allocate(cgrid%Cstore_grow      (npolygons))
     allocate(cgrid%Cleaf_litter_flux(npolygons))
     allocate(cgrid%Croot_litter_flux(npolygons))
     allocate(cgrid%Ccwd_flux        (npolygons))
     allocate(cgrid%Nleaf_grow       (npolygons))
     allocate(cgrid%Ndead_grow       (npolygons))
     allocate(cgrid%Nroot_grow       (npolygons))
     allocate(cgrid%Nstore_grow      (npolygons))
     allocate(cgrid%Nleaf_litter_flux(npolygons))
     allocate(cgrid%Nroot_litter_flux(npolygons))
     allocate(cgrid%Ncwd_flux        (npolygons))
     allocate(cgrid%Nbiomass_uptake  (npolygons))
     allocate(cgrid%Ngross_min       (npolygons))
     allocate(cgrid%Nnet_min         (npolygons))



       ! Meteorologic conditions (forcing)
       allocate(cgrid%avg_nir_beam       (npolygons))
       allocate(cgrid%avg_nir_diffuse    (npolygons))
       allocate(cgrid%avg_par_beam       (npolygons))
       allocate(cgrid%avg_par_diffuse    (npolygons))
       allocate(cgrid%avg_atm_tmp        (npolygons))
       allocate(cgrid%avg_atm_shv        (npolygons))
       allocate(cgrid%avg_rshort         (npolygons))
       allocate(cgrid%avg_rshort_diffuse (npolygons))
       allocate(cgrid%avg_rlong          (npolygons))
       allocate(cgrid%avg_pcpg           (npolygons))
       allocate(cgrid%avg_qpcpg          (npolygons))
       allocate(cgrid%avg_dpcpg          (npolygons))
       allocate(cgrid%avg_vels           (npolygons))
       allocate(cgrid%avg_atm_prss       (npolygons))
       allocate(cgrid%avg_exner          (npolygons))
       allocate(cgrid%avg_geoht          (npolygons))
       allocate(cgrid%avg_atm_co2        (npolygons))
       allocate(cgrid%avg_albedo         (npolygons))
       allocate(cgrid%avg_albedo_beam    (npolygons))
       allocate(cgrid%avg_albedo_diffuse (npolygons))
       allocate(cgrid%avg_rlong_albedo   (npolygons))
       allocate(cgrid%avg_rlongup        (npolygons))
       
       allocate(cgrid%lai_pft            (n_pft       ,npolygons))
       allocate(cgrid%wpa_pft            (n_pft       ,npolygons))
       allocate(cgrid%wai_pft            (n_pft       ,npolygons))

       allocate(cgrid%workload             (13          ,npolygons))

       !-----------------------------------------------------------------------------------!
       !     Allocate the daily means, only if daily means, monthly means, or mean diurnal !
       ! cycles are requested by the user.                                                  !
       !-----------------------------------------------------------------------------------!
       if (idoutput > 0 .or. imoutput > 0 .or. iqoutput > 0) then
          allocate(cgrid%dmean_pcpg           (             npolygons))
          allocate(cgrid%dmean_runoff         (             npolygons))
          allocate(cgrid%dmean_drainage       (             npolygons))
          allocate(cgrid%dmean_vapor_ac       (             npolygons))
          allocate(cgrid%dmean_vapor_gc       (             npolygons))
          allocate(cgrid%dmean_vapor_lc       (             npolygons))
          allocate(cgrid%dmean_vapor_wc       (             npolygons))
          allocate(cgrid%dmean_nee            (             npolygons))
          allocate(cgrid%dmean_gpp            (             npolygons))
          allocate(cgrid%dmean_evap           (             npolygons))
          allocate(cgrid%dmean_transp         (             npolygons))
          allocate(cgrid%dmean_sensible_ac    (             npolygons))
          allocate(cgrid%dmean_sensible_lc    (             npolygons))
          allocate(cgrid%dmean_sensible_wc    (             npolygons))
          allocate(cgrid%dmean_sensible_gc    (             npolygons))
          allocate(cgrid%dmean_plresp         (             npolygons))
          allocate(cgrid%dmean_rh             (             npolygons))
          allocate(cgrid%dmean_leaf_resp      (             npolygons))
          allocate(cgrid%dmean_root_resp      (             npolygons))
          allocate(cgrid%dmean_growth_resp    (             npolygons))
          allocate(cgrid%dmean_storage_resp   (             npolygons))
          allocate(cgrid%dmean_vleaf_resp     (             npolygons))
          allocate(cgrid%dmean_nep            (             npolygons))
          allocate(cgrid%dmean_soil_temp      (nzg,         npolygons))
          allocate(cgrid%dmean_soil_water     (nzg,         npolygons))
          allocate(cgrid%dmean_soil_mstpot    (nzg,         npolygons))
          allocate(cgrid%dmean_transloss      (nzg,         npolygons))
          allocate(cgrid%dmean_fs_open        (             npolygons))
          allocate(cgrid%dmean_fsw            (             npolygons))
          allocate(cgrid%dmean_fsn            (             npolygons))
          allocate(cgrid%dmean_gpp_dbh        (n_dbh       ,npolygons))
          allocate(cgrid%dmean_can_temp       (             npolygons))
          allocate(cgrid%dmean_can_shv        (             npolygons))
          allocate(cgrid%dmean_can_co2        (             npolygons))
          allocate(cgrid%dmean_can_rhos       (             npolygons))
          allocate(cgrid%dmean_can_prss       (             npolygons))
          allocate(cgrid%dmean_can_theta      (             npolygons))
          allocate(cgrid%dmean_can_theiv      (             npolygons))
          allocate(cgrid%dmean_gnd_temp       (             npolygons))
          allocate(cgrid%dmean_gnd_shv        (             npolygons))
          allocate(cgrid%dmean_leaf_energy    (             npolygons))
          allocate(cgrid%dmean_leaf_water     (             npolygons))
          allocate(cgrid%dmean_leaf_hcap      (             npolygons))
          allocate(cgrid%dmean_leaf_temp      (             npolygons))
          allocate(cgrid%dmean_wood_energy    (             npolygons))
          allocate(cgrid%dmean_wood_water     (             npolygons))
          allocate(cgrid%dmean_wood_hcap      (             npolygons))
          allocate(cgrid%dmean_wood_temp      (             npolygons))
          allocate(cgrid%dmean_atm_temp       (             npolygons))
          allocate(cgrid%dmean_atm_shv        (             npolygons))
          allocate(cgrid%dmean_atm_co2        (             npolygons))
          allocate(cgrid%dmean_atm_prss       (             npolygons))
          allocate(cgrid%dmean_atm_vels       (             npolygons))
          allocate(cgrid%dmean_rshort         (             npolygons))
          allocate(cgrid%dmean_rshort_diff    (             npolygons))
          allocate(cgrid%dmean_rlong          (             npolygons))
          allocate(cgrid%dmean_rshort_gnd     (             npolygons))
          allocate(cgrid%dmean_rlong_gnd      (             npolygons))
          allocate(cgrid%dmean_albedo         (             npolygons))
          allocate(cgrid%dmean_albedo_beam    (             npolygons))
          allocate(cgrid%dmean_albedo_diffuse (             npolygons))
          allocate(cgrid%dmean_rlong_albedo   (             npolygons))
          allocate(cgrid%dmean_rlongup        (             npolygons))

          allocate(cgrid%dmean_co2_residual   (             npolygons))
          allocate(cgrid%dmean_energy_residual(             npolygons))
          allocate(cgrid%dmean_water_residual (             npolygons))

       end if
       !-----------------------------------------------------------------------------------!




       !-----------------------------------------------------------------------------------!
       ! Allocate the monthly means only if monthly means or mean diurnal cycles are       !
       ! requested by the user.                                                            !
       !-----------------------------------------------------------------------------------!
       if (imoutput > 0 .or. iqoutput > 0) then
          allocate(cgrid%mmean_pcpg           (             npolygons))
          allocate(cgrid%mmean_runoff         (             npolygons))
          allocate(cgrid%mmean_drainage       (             npolygons))
          allocate(cgrid%mmean_evap           (             npolygons))
          allocate(cgrid%mmean_transp         (             npolygons))
          allocate(cgrid%mmean_vapor_ac       (             npolygons))
          allocate(cgrid%mmean_vapor_gc       (             npolygons))
          allocate(cgrid%mmean_vapor_lc       (             npolygons))
          allocate(cgrid%mmean_vapor_wc       (             npolygons))
          allocate(cgrid%mmean_sensible_ac    (             npolygons))
          allocate(cgrid%mmean_sensible_gc    (             npolygons))
          allocate(cgrid%mmean_sensible_lc    (             npolygons))
          allocate(cgrid%mmean_sensible_wc    (             npolygons))
          allocate(cgrid%mmean_nee            (             npolygons))
          allocate(cgrid%mmean_gpp            (             npolygons))
          allocate(cgrid%mmean_nep            (             npolygons))
          allocate(cgrid%mmean_plresp         (             npolygons))
          allocate(cgrid%mmean_rh             (             npolygons))
          allocate(cgrid%mmean_leaf_resp      (             npolygons))
          allocate(cgrid%mmean_root_resp      (             npolygons))
          allocate(cgrid%mmean_growth_resp    (             npolygons))
          allocate(cgrid%mmean_storage_resp   (             npolygons))
          allocate(cgrid%mmean_vleaf_resp     (             npolygons))
          allocate(cgrid%mmean_soil_temp      (nzg,         npolygons))
          allocate(cgrid%mmean_soil_water     (nzg,         npolygons))
          allocate(cgrid%mmean_soil_mstpot    (nzg,         npolygons))
          allocate(cgrid%mmean_transloss      (nzg,         npolygons))
          allocate(cgrid%mmean_fs_open        (             npolygons))
          allocate(cgrid%mmean_fsw            (             npolygons))
          allocate(cgrid%mmean_fsn            (             npolygons))
          allocate(cgrid%mmean_gpp_dbh        (n_dbh       ,npolygons))
          allocate(cgrid%mmean_lai_pft        (n_pft       ,npolygons))
          allocate(cgrid%mmean_wpa_pft        (n_pft       ,npolygons))
          allocate(cgrid%mmean_wai_pft        (n_pft       ,npolygons))
          allocate(cgrid%mmean_can_temp       (             npolygons))
          allocate(cgrid%mmean_can_shv        (             npolygons))
          allocate(cgrid%mmean_can_co2        (             npolygons))
          allocate(cgrid%mmean_can_rhos       (             npolygons))
          allocate(cgrid%mmean_can_prss       (             npolygons))
          allocate(cgrid%mmean_can_theta      (             npolygons))
          allocate(cgrid%mmean_can_theiv      (             npolygons))
          allocate(cgrid%mmean_gnd_temp       (             npolygons))
          allocate(cgrid%mmean_gnd_shv        (             npolygons))
          allocate(cgrid%mmean_leaf_energy    (             npolygons))
          allocate(cgrid%mmean_leaf_water     (             npolygons))
          allocate(cgrid%mmean_leaf_temp      (             npolygons))
          allocate(cgrid%mmean_leaf_hcap      (             npolygons))
          allocate(cgrid%mmean_wood_energy    (             npolygons))
          allocate(cgrid%mmean_wood_water     (             npolygons))
          allocate(cgrid%mmean_wood_temp      (             npolygons))
          allocate(cgrid%mmean_wood_hcap      (             npolygons))
          allocate(cgrid%mmean_atm_temp       (             npolygons))
          allocate(cgrid%mmean_rshort         (             npolygons))
          allocate(cgrid%mmean_rshort_diff    (             npolygons))
          allocate(cgrid%mmean_rlong          (             npolygons))
          allocate(cgrid%mmean_rshort_gnd     (             npolygons))
          allocate(cgrid%mmean_rlong_gnd      (             npolygons))
          allocate(cgrid%mmean_albedo         (             npolygons))
          allocate(cgrid%mmean_albedo_beam    (             npolygons))
          allocate(cgrid%mmean_albedo_diffuse (             npolygons))
          allocate(cgrid%mmean_rlong_albedo   (             npolygons))
          allocate(cgrid%mmean_rlongup        (             npolygons))
          allocate(cgrid%mmean_atm_shv        (             npolygons))
          allocate(cgrid%mmean_atm_co2        (             npolygons))
          allocate(cgrid%mmean_atm_prss       (             npolygons))
          allocate(cgrid%mmean_atm_vels       (             npolygons))
          allocate(cgrid%mmean_co2_residual   (             npolygons))
          allocate(cgrid%mmean_energy_residual(             npolygons))
          allocate(cgrid%mmean_water_residual (             npolygons))
          allocate(cgrid%bseeds_pft           (n_pft       ,npolygons))
          allocate(cgrid%agb_pft              (n_pft       ,npolygons))
          allocate(cgrid%ba_pft               (n_pft       ,npolygons))
          allocate(cgrid%mmsqu_gpp            (             npolygons))
          allocate(cgrid%mmsqu_leaf_resp      (             npolygons))
          allocate(cgrid%mmsqu_root_resp      (             npolygons))
          allocate(cgrid%mmsqu_plresp         (             npolygons))
          allocate(cgrid%mmsqu_nee            (             npolygons))
          allocate(cgrid%mmsqu_nep            (             npolygons))
          allocate(cgrid%mmsqu_rh             (             npolygons))
          allocate(cgrid%mmsqu_sensible_ac    (             npolygons))
          allocate(cgrid%mmsqu_sensible_lc    (             npolygons))
          allocate(cgrid%mmsqu_sensible_wc    (             npolygons))
          allocate(cgrid%mmsqu_sensible_gc    (             npolygons))
          allocate(cgrid%mmsqu_evap           (             npolygons))
          allocate(cgrid%mmsqu_transp         (             npolygons))
          allocate(cgrid%mmsqu_vapor_ac       (             npolygons))
          allocate(cgrid%mmsqu_vapor_lc       (             npolygons))
          allocate(cgrid%mmsqu_vapor_wc       (             npolygons))
          allocate(cgrid%mmsqu_vapor_gc       (             npolygons))

          allocate(cgrid%disturbance_rates    (n_dist_types,n_dist_types,npolygons))

       end if
       !-----------------------------------------------------------------------------------!




       !-----------------------------------------------------------------------------------!
       !      Allocate the mean diurnal cycles only if mean diurnal cycles are requested   !
       ! by the user.                                                                      !
       !-----------------------------------------------------------------------------------!
       if (iqoutput > 0) then
          allocate(cgrid%qmean_pcpg           (     ndcycle, npolygons))
          allocate(cgrid%qmean_runoff         (     ndcycle, npolygons))
          allocate(cgrid%qmean_drainage       (     ndcycle, npolygons))
          allocate(cgrid%qmean_evap           (     ndcycle, npolygons))
          allocate(cgrid%qmean_transp         (     ndcycle, npolygons))
          allocate(cgrid%qmean_vapor_ac       (     ndcycle, npolygons))
          allocate(cgrid%qmean_vapor_gc       (     ndcycle, npolygons))
          allocate(cgrid%qmean_vapor_lc       (     ndcycle, npolygons))
          allocate(cgrid%qmean_vapor_wc       (     ndcycle, npolygons))
          allocate(cgrid%qmean_sensible_ac    (     ndcycle, npolygons))
          allocate(cgrid%qmean_sensible_gc    (     ndcycle, npolygons))
          allocate(cgrid%qmean_sensible_lc    (     ndcycle, npolygons))
          allocate(cgrid%qmean_sensible_wc    (     ndcycle, npolygons))
          allocate(cgrid%qmean_nee            (     ndcycle, npolygons))
          allocate(cgrid%qmean_gpp            (     ndcycle, npolygons))
          allocate(cgrid%qmean_nep            (     ndcycle, npolygons))
          allocate(cgrid%qmean_plresp         (     ndcycle, npolygons))
          allocate(cgrid%qmean_rh             (     ndcycle, npolygons))
          allocate(cgrid%qmean_leaf_resp      (     ndcycle, npolygons))
          allocate(cgrid%qmean_root_resp      (     ndcycle, npolygons))
          allocate(cgrid%qmean_soil_temp      (nzg, ndcycle, npolygons))
          allocate(cgrid%qmean_soil_water     (nzg, ndcycle, npolygons))
          allocate(cgrid%qmean_soil_mstpot    (nzg, ndcycle, npolygons))
          allocate(cgrid%qmean_fs_open        (     ndcycle, npolygons))
          allocate(cgrid%qmean_fsw            (     ndcycle, npolygons))
          allocate(cgrid%qmean_fsn            (     ndcycle, npolygons))
          allocate(cgrid%qmean_can_temp       (     ndcycle, npolygons))
          allocate(cgrid%qmean_can_shv        (     ndcycle, npolygons))
          allocate(cgrid%qmean_can_co2        (     ndcycle, npolygons))
          allocate(cgrid%qmean_can_rhos       (     ndcycle, npolygons))
          allocate(cgrid%qmean_can_prss       (     ndcycle, npolygons))
          allocate(cgrid%qmean_can_theta      (     ndcycle, npolygons))
          allocate(cgrid%qmean_can_theiv      (     ndcycle, npolygons))
          allocate(cgrid%qmean_gnd_temp       (     ndcycle, npolygons))
          allocate(cgrid%qmean_gnd_shv        (     ndcycle, npolygons))
          allocate(cgrid%qmean_leaf_energy    (     ndcycle, npolygons))
          allocate(cgrid%qmean_leaf_water     (     ndcycle, npolygons))
          allocate(cgrid%qmean_leaf_temp      (     ndcycle, npolygons))
          allocate(cgrid%qmean_leaf_hcap      (     ndcycle, npolygons))
          allocate(cgrid%qmean_wood_energy    (     ndcycle, npolygons))
          allocate(cgrid%qmean_wood_water     (     ndcycle, npolygons))
          allocate(cgrid%qmean_wood_temp      (     ndcycle, npolygons))
          allocate(cgrid%qmean_wood_hcap      (     ndcycle, npolygons))
          allocate(cgrid%qmean_atm_temp       (     ndcycle, npolygons))
          allocate(cgrid%qmean_rshort         (     ndcycle, npolygons))
          allocate(cgrid%qmean_rshort_diff    (     ndcycle, npolygons))
          allocate(cgrid%qmean_rlong          (     ndcycle, npolygons))
          allocate(cgrid%qmean_rshort_gnd     (     ndcycle, npolygons))
          allocate(cgrid%qmean_rlong_gnd      (     ndcycle, npolygons))
          allocate(cgrid%qmean_albedo         (     ndcycle, npolygons))
          allocate(cgrid%qmean_albedo_beam    (     ndcycle, npolygons))
          allocate(cgrid%qmean_albedo_diffuse (     ndcycle, npolygons))
          allocate(cgrid%qmean_rlong_albedo   (     ndcycle, npolygons))
          allocate(cgrid%qmean_rlongup        (     ndcycle, npolygons))
          allocate(cgrid%qmean_atm_shv        (     ndcycle, npolygons))
          allocate(cgrid%qmean_atm_co2        (     ndcycle, npolygons))
          allocate(cgrid%qmean_atm_prss       (     ndcycle, npolygons))
          allocate(cgrid%qmean_atm_vels       (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_gpp            (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_leaf_resp      (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_root_resp      (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_plresp         (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_nee            (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_nep            (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_rh             (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_sensible_ac    (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_sensible_lc    (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_sensible_wc    (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_sensible_gc    (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_evap           (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_transp         (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_vapor_ac       (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_vapor_lc       (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_vapor_wc       (     ndcycle, npolygons))
          allocate(cgrid%qmsqu_vapor_gc       (     ndcycle, npolygons))

       end if
    end if
    return
  end subroutine allocate_edtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine allocate_polygontype(cpoly,nsites)

    implicit none

    integer :: nsites
    type(polygontype),target :: cpoly

    call nullify_polygontype(cpoly)

    cpoly%nsites = nsites

    allocate(cpoly%sipa_id(nsites))
    allocate(cpoly%sipa_n(nsites))

    allocate(cpoly%patch_count(nsites))  
    allocate(cpoly%site(nsites))
    allocate(cpoly%sitenum(nsites))

    allocate(cpoly%lsl(nsites))   

    allocate(cpoly%area(nsites))
    allocate(cpoly%patch_area(nsites))
    allocate(cpoly%elevation(nsites))
    allocate(cpoly%slope(nsites))
    allocate(cpoly%aspect(nsites))
    
    allocate(cpoly%num_landuse_years(nsites))
    allocate(cpoly%mindbh_primary(n_pft,nsites))
    allocate(cpoly%probharv_primary(n_pft,nsites))
    allocate(cpoly%mindbh_secondary(n_pft,nsites))
    allocate(cpoly%probharv_secondary(n_pft,nsites))

    allocate(cpoly%lai_pft(n_pft,nsites))
    allocate(cpoly%wpa_pft(n_pft,nsites))
    allocate(cpoly%wai_pft(n_pft,nsites))

    allocate(cpoly%TCI(nsites))      
    allocate(cpoly%pptweight(nsites))      
    allocate(cpoly%hydro_next(nsites))
    allocate(cpoly%hydro_prev(nsites))
    allocate(cpoly%moist_W(nsites))
    allocate(cpoly%moist_f(nsites))  
    allocate(cpoly%moist_tau(nsites))
    allocate(cpoly%moist_zi(nsites)) 
    allocate(cpoly%baseflow(nsites)) 
    allocate(cpoly%ntext_soil(nzg,nsites))
    allocate(cpoly%min_monthly_temp(nsites))
    allocate(cpoly%plantation(nsites)) 
    allocate(cpoly%agri_stocking_pft(nsites))
    allocate(cpoly%agri_stocking_density(nsites))
    allocate(cpoly%plantation_stocking_pft(nsites))
    allocate(cpoly%plantation_stocking_density(nsites))
    allocate(cpoly%primary_harvest_memory(nsites))
    allocate(cpoly%secondary_harvest_memory(nsites))
    allocate(cpoly%fire_disturbance_rate(nsites))
    allocate(cpoly%ignition_rate(nsites))
    allocate(cpoly%lambda_fire(12,nsites))
    allocate(cpoly%phen_pars(nsites))!THIS PTR IS ALLOCATED IN PHENOLOGY_INIT
    allocate(cpoly%nat_disturbance_rate(nsites))
    allocate(cpoly%nat_dist_type(nsites))
    allocate(cpoly%disturbance_memory(n_dist_types,n_dist_types,nsites))
    allocate(cpoly%disturbance_rates(n_dist_types,n_dist_types,nsites))
    allocate(cpoly%loss_fraction(n_dist_types,nsites))

    allocate(cpoly%green_leaf_factor(n_pft,nsites))
    allocate(cpoly%leaf_aging_factor(n_pft,nsites))
    
    allocate(cpoly%met(nsites))

    allocate(cpoly%basal_area  (n_pft,n_dbh,nsites))
    allocate(cpoly%agb         (n_pft,n_dbh,nsites))
    allocate(cpoly%pldens      (n_pft,n_dbh,nsites))
    allocate(cpoly%bseeds      (n_pft,n_dbh,nsites))

    allocate(cpoly%basal_area_growth (n_pft,n_dbh,nsites))
    allocate(cpoly%agb_growth        (n_pft,n_dbh,nsites))
    allocate(cpoly%basal_area_mort   (n_pft,n_dbh,nsites))
    allocate(cpoly%basal_area_cut   (n_pft,n_dbh,nsites))
    allocate(cpoly%agb_mort          (n_pft,n_dbh,nsites))
    allocate(cpoly%agb_cut          (n_pft,n_dbh,nsites))

    allocate(cpoly%cosaoi(nsites))
    allocate(cpoly%avg_albedo_beam(nsites))
    allocate(cpoly%avg_albedo_diffuse(nsites))
    allocate(cpoly%avg_rlong_albedo(nsites))
    
    allocate(cpoly%avg_albedo(nsites))
    allocate(cpoly%avg_rlongup(nsites))
    allocate(cpoly%daylight(nsites))

    allocate(cpoly%lai  (nsites)) 
    allocate(cpoly%avg_lma(nsites))
    allocate(cpoly%wpa  (nsites))
    allocate(cpoly%wai  (nsites))
    ! Fast time flux diagnostics
    ! ---------------------------------------------
    allocate(cpoly%avg_vapor_lc  (nsites))
    allocate(cpoly%avg_vapor_wc  (nsites))
    allocate(cpoly%avg_dew_cg    (nsites))
    allocate(cpoly%avg_vapor_gc  (nsites))
    allocate(cpoly%avg_wshed_vg  (nsites))
    allocate(cpoly%avg_intercepted (nsites))
    allocate(cpoly%avg_throughfall (nsites))
    allocate(cpoly%avg_vapor_ac  (nsites))
    allocate(cpoly%avg_transp    (nsites))
    allocate(cpoly%avg_evap      (nsites))
    
    allocate(cpoly%avg_smoist_gg (nzg,nsites))
    allocate(cpoly%avg_transloss (nzg,nsites))
    allocate(cpoly%avg_runoff        (nsites))
    allocate(cpoly%avg_drainage  (nsites))
    allocate(cpoly%avg_drainage_heat  (nsites))
    allocate(cpoly%aux           (nsites))
    allocate(cpoly%aux_s         (nzg,nsites))
    allocate(cpoly%avg_rshort_gnd   (nsites))
    allocate(cpoly%avg_rlong_gnd    (nsites))
    allocate(cpoly%avg_carbon_ac    (nsites))
    allocate(cpoly%avg_sensible_lc  (nsites))
    allocate(cpoly%avg_sensible_wc  (nsites))
    allocate(cpoly%avg_qwshed_vg    (nsites))
    allocate(cpoly%avg_qintercepted (nsites))
    allocate(cpoly%avg_qthroughfall (nsites))
    allocate(cpoly%avg_sensible_gc  (nsites))
    allocate(cpoly%avg_sensible_ac  (nsites))
    allocate(cpoly%avg_sensible_gg  (nzg,nsites))
    allocate(cpoly%avg_runoff_heat         (nsites))

    ! Fast time state diagnostics
    allocate(cpoly%avg_leaf_energy          (nsites))
    allocate(cpoly%avg_leaf_hcap            (nsites))
    allocate(cpoly%avg_leaf_temp            (nsites))
    allocate(cpoly%avg_leaf_fliq            (nsites))
    allocate(cpoly%avg_leaf_water           (nsites))
    allocate(cpoly%avg_wood_energy          (nsites))
    allocate(cpoly%avg_wood_hcap            (nsites))
    allocate(cpoly%avg_wood_temp            (nsites))
    allocate(cpoly%avg_wood_fliq            (nsites))
    allocate(cpoly%avg_wood_water           (nsites))
    allocate(cpoly%avg_can_temp            (nsites))
    allocate(cpoly%avg_can_shv             (nsites))
    allocate(cpoly%avg_can_co2             (nsites))
    allocate(cpoly%avg_can_rhos            (nsites))
    allocate(cpoly%avg_can_prss            (nsites))
    allocate(cpoly%avg_can_theta           (nsites))
    allocate(cpoly%avg_can_theiv           (nsites))
    allocate(cpoly%avg_can_depth           (nsites))
    allocate(cpoly%avg_soil_energy     (nzg,nsites))
    allocate(cpoly%avg_soil_mstpot     (nzg,nsites))
    allocate(cpoly%avg_soil_water      (nzg,nsites))
    allocate(cpoly%avg_soil_temp       (nzg,nsites))
    allocate(cpoly%avg_soil_fracliq    (nzg,nsites))
    allocate(cpoly%avg_soil_rootfrac   (nzg,nsites))
    allocate(cpoly%avg_soil_wetness        (nsites))
    allocate(cpoly%avg_skin_temp           (nsites))
    allocate(cpoly%avg_available_water     (nsites))
    allocate(cpoly%runoff                  (nsites))
    ! Phenology-related
    allocate(cpoly%rad_avg                 (nsites))
    ! Meteorological data
    allocate(cpoly%avg_atm_tmp             (nsites))
    allocate(cpoly%avg_atm_shv             (nsites))
    allocate(cpoly%avg_atm_prss            (nsites))

    !!!NACP
    allocate(cpoly%avg_sfcw_depth           (nsites))
    allocate(cpoly%avg_sfcw_energy          (nsites))
    allocate(cpoly%avg_sfcw_mass            (nsites))
    allocate(cpoly%avg_sfcw_fracliq         (nsites))
    allocate(cpoly%avg_sfcw_tempk           (nsites))
    allocate(cpoly%avg_fsc                 (nsites))
    allocate(cpoly%avg_stsc                (nsites))
    allocate(cpoly%avg_ssc                 (nsites))
    allocate(cpoly%avg_balive              (nsites))
    allocate(cpoly%avg_bleaf               (nsites))
    allocate(cpoly%avg_broot               (nsites))
    allocate(cpoly%avg_bsapwood            (nsites))
    allocate(cpoly%avg_bdead               (nsites))
    allocate(cpoly%avg_fsn                 (nsites))
    allocate(cpoly%avg_msn                 (nsites))


    allocate(cpoly%avg_bstorage            (nsites))
    allocate(cpoly%avg_bseeds              (nsites))


    allocate(cpoly%dmean_co2_residual      (nsites))
    allocate(cpoly%dmean_energy_residual   (nsites))
    allocate(cpoly%dmean_water_residual    (nsites))
    allocate(cpoly%mmean_co2_residual      (nsites))
    allocate(cpoly%mmean_energy_residual   (nsites))
    allocate(cpoly%mmean_water_residual    (nsites))
    
    return
  end subroutine allocate_polygontype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine allocate_sitetype(csite,npatches)

    use soil_bgc, only: allocate_sitetype_sbgc
    implicit none

    integer :: npatches,ipa
    type(sitetype),target :: csite
    
    call nullify_sitetype(csite)
    csite%npatches = npatches
    allocate(csite%paco_id(npatches))
    allocate(csite%paco_n(npatches))
    allocate(csite%patch(npatches))

    ! Initialize zero cohorts
    do ipa=1,npatches
       csite%patch(ipa)%ncohorts = 0
    end do

    allocate(csite%lai(npatches))
    allocate(csite%wpa(npatches))
    allocate(csite%wai(npatches))
    allocate(csite%dist_type(npatches))
    allocate(csite%age(npatches))
    allocate(csite%area(npatches))
    allocate(csite%laiarea(npatches))

    call allocate_sitetype_sbgc(csite%sbgc, npatches)

    allocate(csite%sum_dgd(npatches))
    allocate(csite%sum_chd(npatches))
    allocate(csite%plantation(npatches))
    allocate(csite%can_theiv(npatches))
    allocate(csite%can_temp(npatches))
    allocate(csite%can_shv(npatches))
    allocate(csite%can_co2(npatches))
    allocate(csite%can_rhos(npatches))
    allocate(csite%can_prss(npatches))
    allocate(csite%can_theta(npatches))
    allocate(csite%can_depth(npatches))
    allocate(csite%opencan_frac(npatches))
    allocate(csite%ggbare(npatches))
    allocate(csite%ggveg(npatches))
    allocate(csite%ggnet(npatches))
    allocate(csite%ggsoil(npatches))
    allocate(csite%lambda_light(npatches))
    allocate(csite%cohort_count(npatches))
    allocate(csite%pname(npatches))
    
    allocate(csite%sfcwater_mass(nzs,npatches))
    allocate(csite%sfcwater_energy(nzs,npatches))
    allocate(csite%sfcwater_depth(nzs,npatches))
    allocate(csite%rshort_s(nzs,npatches))
    allocate(csite%rshort_s_beam(nzs,npatches))
    allocate(csite%rshort_s_diffuse(nzs,npatches))
    allocate(csite%sfcwater_tempk(nzs,npatches))
    allocate(csite%sfcwater_fracliq(nzs,npatches))
    allocate(csite%nlev_sfcwater(npatches))
    allocate(csite%soil_energy(nzg,npatches))
    allocate(csite%soil_water(nzg,npatches))
    allocate(csite%soil_tempk(nzg,npatches))
    allocate(csite%soil_fracliq(nzg,npatches))
    allocate(csite%rootdense(nzg,npatches))
    allocate(csite%ground_shv(npatches))
    allocate(csite%ground_ssh(npatches))
    allocate(csite%ground_temp(npatches))
    allocate(csite%ground_fliq(npatches))
    allocate(csite%rough(npatches))
    allocate(csite%par_l_max(npatches))
    allocate(csite%par_l_beam_max(npatches))
    allocate(csite%par_l_diffuse_max(npatches))
    allocate(csite%A_o_max(n_pft,npatches)) 
    allocate(csite%A_c_max(n_pft,npatches)) 

    allocate(csite%old_stoma_data_max(n_pft,npatches))
    allocate(csite%old_stoma_vector_max(n_stoma_atts,n_pft,npatches))

    allocate(csite%avg_daily_temp(npatches))  
    allocate(csite%mean_rh(npatches))
    allocate(csite%mean_nep(npatches))
    allocate(csite%wbudget_loss2atm(npatches))
    allocate(csite%wbudget_denseffect(npatches))
    allocate(csite%wbudget_precipgain(npatches))
    allocate(csite%wbudget_loss2runoff(npatches))
    allocate(csite%wbudget_loss2drainage(npatches))
    allocate(csite%wbudget_initialstorage(npatches))
    allocate(csite%wbudget_residual(npatches))
    allocate(csite%ebudget_loss2atm(npatches))
    allocate(csite%ebudget_denseffect(npatches))
    allocate(csite%ebudget_loss2runoff(npatches))
    allocate(csite%ebudget_loss2drainage(npatches))
    allocate(csite%ebudget_netrad(npatches))
    allocate(csite%ebudget_precipgain(npatches))
    allocate(csite%ebudget_initialstorage(npatches))
    allocate(csite%ebudget_residual(npatches))
    allocate(csite%co2budget_initialstorage(npatches))
    allocate(csite%co2budget_residual(npatches))
    allocate(csite%co2budget_loss2atm(npatches))
    allocate(csite%co2budget_denseffect(npatches))
    allocate(csite%co2budget_gpp(npatches))
    allocate(csite%co2budget_gpp_dbh(n_dbh,npatches))
    allocate(csite%co2budget_plresp(npatches))
    allocate(csite%co2budget_rh(npatches))
    allocate(csite%today_A_decomp(npatches))
    allocate(csite%today_Af_decomp(npatches))
    allocate(csite%repro(n_pft,npatches))
    allocate(csite%veg_rough(npatches))
    allocate(csite%veg_height (npatches))
    allocate(csite%veg_displace (npatches))
    allocate(csite%fsc_in(npatches))
    allocate(csite%ssc_in(npatches))
    allocate(csite%ssl_in(npatches))
    allocate(csite%fsn_in(npatches))


    allocate(csite%fsp_in(npatches))
    allocate(csite%stsp_in(npatches))
    allocate(csite%total_plant_P_uptake(npatches))


    allocate(csite%total_plant_nitrogen_uptake(npatches))
    allocate(csite%mineralized_N_loss(npatches))
    allocate(csite%mineralized_N_input(npatches))
    allocate(csite%rshort_g(npatches))
    allocate(csite%rshort_g_beam(npatches))
    allocate(csite%rshort_g_diffuse(npatches))
    allocate(csite%rlong_g(npatches))
    allocate(csite%rlong_g_surf(npatches))
    allocate(csite%rlong_g_incid(npatches))
    allocate(csite%rlong_s(npatches))
    allocate(csite%rlong_s_surf(npatches))
    allocate(csite%rlong_s_incid(npatches))
    allocate(csite%albedo(npatches))
    allocate(csite%albedo_beam(npatches))
    allocate(csite%albedo_diffuse(npatches))
    allocate(csite%rlongup(npatches))
    allocate(csite%rlong_albedo(npatches))
    allocate(csite%total_sfcw_depth(npatches))
    allocate(csite%snowfac(npatches))
    allocate(csite%A_decomp(npatches))
    allocate(csite%f_decomp(npatches))
    allocate(csite%rh(npatches))
    allocate(csite%cwd_rh(npatches))
    allocate(csite%fuse_flag(npatches))
    allocate(csite%cumlai_profile(n_pft,ff_nhgt,npatches))
    allocate(csite%plant_ag_biomass(npatches))

    allocate(csite%mean_wflux(npatches))
    allocate(csite%mean_latflux(npatches))
    allocate(csite%mean_hflux(npatches))
    allocate(csite%mean_runoff(npatches))
    allocate(csite%mean_qrunoff(npatches))

    allocate(csite%htry       (npatches))
    allocate(csite%avg_rk4step(npatches))

    allocate(csite%avg_available_water(npatches))

    allocate(csite%ustar(npatches))
    allocate(csite%tstar(npatches))
    allocate(csite%qstar(npatches))
    allocate(csite%cstar(npatches))
    
    allocate(csite%zeta  (npatches))
    allocate(csite%ribulk(npatches))

    allocate(csite%upwp(npatches))
    allocate(csite%qpwp(npatches))
    allocate(csite%cpwp(npatches))
    allocate(csite%tpwp(npatches))
    allocate(csite%wpwp(npatches))

    allocate(csite%avg_rshort_gnd     (npatches))
    allocate(csite%avg_rlong_gnd      (npatches))
    allocate(csite%avg_carbon_ac      (npatches))
    allocate(csite%avg_rlongup        (npatches))
    allocate(csite%avg_albedo         (npatches))
    allocate(csite%avg_albedo_beam    (npatches))
    allocate(csite%avg_albedo_diffuse (npatches))
    allocate(csite%avg_rlong_albedo   (npatches))

    ! Fast time flux diagnostics
    ! ---------------------------------------------
    allocate(csite%avg_vapor_lc  (npatches))
    allocate(csite%avg_vapor_wc  (npatches))
    allocate(csite%avg_dew_cg    (npatches))
    allocate(csite%avg_vapor_gc  (npatches))
    allocate(csite%avg_wshed_vg  (npatches))
    allocate(csite%avg_intercepted (npatches))
    allocate(csite%avg_throughfall (npatches))
    allocate(csite%avg_vapor_ac  (npatches))
    allocate(csite%avg_transp    (npatches))
    allocate(csite%avg_evap      (npatches))
    allocate(csite%avg_smoist_gg (nzg,npatches))
    allocate(csite%avg_transloss (nzg,npatches))
    allocate(csite%avg_runoff    (npatches))
    allocate(csite%avg_drainage  (npatches))
    allocate(csite%avg_drainage_heat  (npatches))
    allocate(csite%aux           (npatches))
    allocate(csite%aux_s         (nzg,npatches))
    allocate(csite%avg_sensible_lc  (npatches))
    allocate(csite%avg_sensible_wc  (npatches))
    allocate(csite%avg_qwshed_vg    (npatches))
    allocate(csite%avg_qintercepted (npatches))
    allocate(csite%avg_qthroughfall (npatches))
    allocate(csite%avg_sensible_gc  (npatches))
    allocate(csite%avg_sensible_ac  (npatches))
    allocate(csite%avg_sensible_gg  (nzg,npatches))
    allocate(csite%avg_runoff_heat  (npatches))

    ! ----------------------------------------------

    allocate(csite%avg_leaf_energy(npatches))
    allocate(csite%avg_leaf_temp  (npatches))
    allocate(csite%avg_leaf_hcap  (npatches))
    allocate(csite%avg_leaf_fliq  (npatches))
    allocate(csite%avg_leaf_water (npatches))
    allocate(csite%avg_wood_energy(npatches))
    allocate(csite%avg_wood_temp  (npatches))
    allocate(csite%avg_wood_hcap  (npatches))
    allocate(csite%avg_wood_fliq  (npatches))
    allocate(csite%avg_wood_water (npatches))


    allocate(csite%watertable      (npatches))
    allocate(csite%moist_dz        (npatches))
    allocate(csite%ksat            (npatches))
    allocate(csite%soil_sat_energy (npatches))
    allocate(csite%soil_sat_water  (npatches))
    allocate(csite%soil_sat_heat   (npatches))

    allocate(csite%runoff_A        (3,npatches))
    allocate(csite%runoff_rate     (npatches))
    allocate(csite%runoff          (npatches))

    if (imoutput > 0 .or. idoutput > 0 .or. iqoutput > 0) then
       allocate(csite%dmean_rk4step           (npatches))
       allocate(csite%dmean_lambda_light      (npatches))
       allocate(csite%dmean_co2_residual      (npatches))
       allocate(csite%dmean_energy_residual   (npatches))
       allocate(csite%dmean_water_residual    (npatches))
       allocate(csite%dmean_rh                (npatches))
       allocate(csite%dmean_A_decomp          (npatches))
       allocate(csite%dmean_Af_decomp         (npatches))
       allocate(csite%dmean_albedo            (npatches))
       allocate(csite%dmean_albedo_beam       (npatches))
       allocate(csite%dmean_albedo_diffuse    (npatches))
    end if
    if (imoutput > 0 .or. iqoutput > 0) then
       allocate(csite%mmean_rk4step           (npatches))
       allocate(csite%mmean_lambda_light      (npatches))
       allocate(csite%mmean_co2_residual      (npatches))
       allocate(csite%mmean_energy_residual   (npatches))
       allocate(csite%mmean_water_residual    (npatches))
       allocate(csite%mmean_rh                (npatches))
       allocate(csite%mmean_A_decomp          (npatches))
       allocate(csite%mmean_Af_decomp         (npatches))
       allocate(csite%mmean_albedo            (npatches))
       allocate(csite%mmean_albedo_beam       (npatches))
       allocate(csite%mmean_albedo_diffuse    (npatches))
    end if
    if (iqoutput > 0) then
       allocate(csite%qmean_rh            (ndcycle,npatches))
       allocate(csite%qmean_albedo        (ndcycle,npatches))
       allocate(csite%qmean_albedo_beam   (ndcycle,npatches))
       allocate(csite%qmean_albedo_diffuse(ndcycle,npatches))
    end if

    allocate(csite%past_paw(nzg,10,npatches))
    allocate(csite%current_paw(nzg,npatches))

    return
  end subroutine allocate_sitetype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine allocate_patchtype(cpatch,ncohorts)

    use cohort_state, only: allocate_patchtype_state
    use cohort_mort, only: allocate_patchtype_mort
    use cohort_phen, only: allocate_patchtype_phen
    use cohort_photo, only: allocate_patchtype_photo
    use cohort_rad, only: allocate_patchtype_rad
    use cohort_therm, only: allocate_patchtype_therm
    use cohort_resp, only: allocate_patchtype_resp

    implicit none

    integer :: ncohorts
    type(patchtype),target :: cpatch
    
    cpatch%ncohorts = ncohorts

    !----- We need to leave the subroutine in case this patch is empty. -------------------!
    if (ncohorts == 0) return
    
    call allocate_patchtype_state(cpatch%costate, ncohorts)
    call allocate_patchtype_phen(cpatch%cophen, ncohorts, imoutput, iqoutput)
    call allocate_patchtype_mort(cpatch%comort, ncohorts, idoutput, imoutput,  &
         iqoutput, n_mort)
    call allocate_patchtype_therm(cpatch%cotherm, ncohorts)
    call allocate_patchtype_photo(cpatch%cophoto, ncohorts, idoutput, imoutput,  &
         iqoutput, ndcycle)
    call allocate_patchtype_rad(cpatch%corad, ncohorts)
    call allocate_patchtype_resp(cpatch%coresp, ncohorts, idoutput, imoutput,  &
         iqoutput, ndcycle)


    return
  end subroutine allocate_patchtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine nullify_edtype(cgrid)
    
    implicit none
    type(edtype),target :: cgrid


       nullify(cgrid%polygon                 )
       nullify(cgrid%lat                     )
       nullify(cgrid%lon                     )
       nullify(cgrid%xatm                    )
       nullify(cgrid%yatm                    )
       nullify(cgrid%ntext_soil              )
       nullify(cgrid%lsl                     )
       
       nullify(cgrid%pysi_id                 )
       nullify(cgrid%pysi_n                  )
       nullify(cgrid%walltime_py             )
       nullify(cgrid%sensflux_py             )
       nullify(cgrid%site_adjacency          )
       nullify(cgrid%wbar                    )
       nullify(cgrid%Te                      )
       nullify(cgrid%zbar                    )
       nullify(cgrid%sheat                   )
       nullify(cgrid%baseflow                )
       nullify(cgrid%runoff                  )
       nullify(cgrid%swliq                   )
       
       nullify(cgrid%ilon                    )
       nullify(cgrid%ilat                    )
       nullify(cgrid%total_agb               )
       nullify(cgrid%total_basal_area        )
       nullify(cgrid%total_agb_growth        )
       nullify(cgrid%total_agb_mort          )
       nullify(cgrid%total_agb_recruit       )
       nullify(cgrid%total_basal_area_growth )
       nullify(cgrid%total_basal_area_mort   )
       nullify(cgrid%total_basal_area_recruit)
       nullify(cgrid%nsites                  )
       nullify(cgrid%sitenums                )
       nullify(cgrid%load_adjacency          )
       nullify(cgrid%cosz                    )
       nullify(cgrid%mean_gpp                )
       nullify(cgrid%mean_precip             )
       nullify(cgrid%mean_qprecip            )
       nullify(cgrid%mean_netrad             )
       nullify(cgrid%cbudget_initialstorage  )
       nullify(cgrid%cbudget_nep             )
       nullify(cgrid%nbudget_initialstorage  )
       nullify(cgrid%basal_area              )
       nullify(cgrid%agb                     )
       nullify(cgrid%pldens                  )
       nullify(cgrid%bseeds                  )
       
       nullify(cgrid%metinput                )
       nullify(cgrid%met                     )

       nullify(cgrid%lapse                   )

       nullify(cgrid%lai                     )
       nullify(cgrid%avg_lma                 )
       nullify(cgrid%wpa                     )
       nullify(cgrid%wai                     )

       ! Fast time flux diagnostics
       ! ---------------------------------------------
       nullify(cgrid%avg_vapor_lc            )
       nullify(cgrid%avg_vapor_wc            )
       nullify(cgrid%avg_dew_cg              )
       nullify(cgrid%avg_vapor_gc            )
       nullify(cgrid%avg_wshed_vg            )
       nullify(cgrid%avg_intercepted         )
       nullify(cgrid%avg_throughfall         )
       nullify(cgrid%avg_vapor_ac            )
       nullify(cgrid%avg_transp              )
       nullify(cgrid%avg_evap                )
       nullify(cgrid%avg_smoist_gg           )
       nullify(cgrid%avg_transloss           )
       nullify(cgrid%avg_runoff              )
       nullify(cgrid%avg_drainage            )
       nullify(cgrid%avg_drainage_heat       )
       nullify(cgrid%aux                     )
       nullify(cgrid%aux_s                   )
       nullify(cgrid%avg_rshort_gnd          )
       nullify(cgrid%avg_rlong_gnd           )
       nullify(cgrid%avg_carbon_ac           )
       nullify(cgrid%avg_sensible_lc         )
       nullify(cgrid%avg_sensible_wc         )
       nullify(cgrid%avg_qwshed_vg           )
       nullify(cgrid%avg_qintercepted        )
       nullify(cgrid%avg_qthroughfall        )
       nullify(cgrid%avg_sensible_gc         )
       nullify(cgrid%avg_sensible_ac         )
       nullify(cgrid%avg_sensible_gg         )
       nullify(cgrid%avg_runoff_heat         )

       ! Fast time state diagnostics 
       nullify(cgrid%avg_leaf_energy         )
       nullify(cgrid%avg_leaf_hcap           )
       nullify(cgrid%avg_leaf_temp           )
       nullify(cgrid%avg_leaf_fliq           )
       nullify(cgrid%avg_leaf_water          )
       nullify(cgrid%avg_wood_energy         )
       nullify(cgrid%avg_wood_hcap           )
       nullify(cgrid%avg_wood_temp           )
       nullify(cgrid%avg_wood_fliq           )
       nullify(cgrid%avg_wood_water          )
       nullify(cgrid%avg_can_temp            )
       nullify(cgrid%avg_can_shv             )
       nullify(cgrid%avg_can_co2             )
       nullify(cgrid%avg_can_rhos            )
       nullify(cgrid%avg_can_prss            )
       nullify(cgrid%avg_can_theta           )
       nullify(cgrid%avg_can_theiv           )
       nullify(cgrid%avg_can_depth           )
       nullify(cgrid%avg_soil_energy         )
       nullify(cgrid%avg_soil_mstpot         )
       nullify(cgrid%avg_soil_water          )
       nullify(cgrid%avg_soil_temp           )
       nullify(cgrid%avg_soil_fracliq        )
       nullify(cgrid%avg_soil_rootfrac       )
       nullify(cgrid%avg_soil_wetness        )
       nullify(cgrid%avg_skin_temp           )
       nullify(cgrid%avg_available_water     )

       nullify(cgrid%avg_lai_ebalvars)

       nullify(cgrid%avg_gpp             )
       nullify(cgrid%avg_leaf_resp       )
       nullify(cgrid%avg_root_resp       )
       nullify(cgrid%avg_growth_resp     )
       nullify(cgrid%avg_storage_resp    )
       nullify(cgrid%avg_vleaf_resp      )
       nullify(cgrid%avg_plant_resp      )
       nullify(cgrid%avg_htroph_resp     )
       nullify(cgrid%avg_leaf_drop       )
       nullify(cgrid%avg_leaf_maintenance)
       nullify(cgrid%avg_root_maintenance)
       

       !!! added for NACP intercomparison (MCD)
       nullify(cgrid%avg_sfcw_depth     )
       nullify(cgrid%avg_sfcw_energy    )
       nullify(cgrid%avg_sfcw_mass      )
       nullify(cgrid%avg_sfcw_tempk     )
       nullify(cgrid%avg_sfcw_fracliq   )
       nullify(cgrid%avg_bdead         )
       nullify(cgrid%avg_balive        )
       nullify(cgrid%avg_bleaf          )
       nullify(cgrid%avg_broot          )
       nullify(cgrid%avg_bsapwood       )
       nullify(cgrid%avg_bstorage       )
       nullify(cgrid%avg_bseeds         )

       nullify(cgrid%avg_fsc           )
       nullify(cgrid%avg_ssc           )
       nullify(cgrid%avg_stsc          )
       nullify(cgrid%avg_fsn           )
       nullify(cgrid%avg_msn           )


       
     !-------- TOTAL CARBON AND NITROGEN POOLS  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
     nullify(cgrid%Cleaf  )
     nullify(cgrid%Croot  )
     nullify(cgrid%Cstore )
     nullify(cgrid%Ccwd   )
     nullify(cgrid%Nleaf  )
     nullify(cgrid%Ndead  )
     nullify(cgrid%Nroot  )
     nullify(cgrid%Nstore )
     nullify(cgrid%Ncwd   )
     
     
     !-------- TOTAL CARBON AND NITROGEN FLUX  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
     nullify(cgrid%Cleaf_grow       )
     nullify(cgrid%Croot_grow       )
     nullify(cgrid%Cdead_grow       )
     nullify(cgrid%Cstore_grow      )
     nullify(cgrid%Cleaf_litter_flux)
     nullify(cgrid%Croot_litter_flux)
     nullify(cgrid%Ccwd_flux        )
     nullify(cgrid%Nleaf_grow       )
     nullify(cgrid%Ndead_grow       )
     nullify(cgrid%Nroot_grow       )
     nullify(cgrid%Nstore_grow      )
     nullify(cgrid%Nleaf_litter_flux)
     nullify(cgrid%Nroot_litter_flux)
     nullify(cgrid%Ncwd_flux        )
     nullify(cgrid%Nbiomass_uptake  )
     nullify(cgrid%Ngross_min       )
     nullify(cgrid%Nnet_min         )


     ! Meteorologic conditions (forcing)
     nullify(cgrid%avg_nir_beam            )
     nullify(cgrid%avg_nir_diffuse         )
     nullify(cgrid%avg_par_beam            )
     nullify(cgrid%avg_par_diffuse         )
     nullify(cgrid%avg_atm_tmp             )
     nullify(cgrid%avg_atm_shv             )
     nullify(cgrid%avg_rshort              )
     nullify(cgrid%avg_rshort_diffuse      )
     nullify(cgrid%avg_rlong               )
     nullify(cgrid%avg_pcpg                )
     nullify(cgrid%avg_qpcpg               )
     nullify(cgrid%avg_dpcpg               )
     nullify(cgrid%avg_vels                )
     nullify(cgrid%avg_atm_prss            )
     nullify(cgrid%avg_exner               )
     nullify(cgrid%avg_geoht               )
     nullify(cgrid%avg_atm_co2             )
     nullify(cgrid%avg_albedo              )
     nullify(cgrid%avg_albedo_beam         )
     nullify(cgrid%avg_albedo_diffuse      )
     nullify(cgrid%avg_rlong_albedo        )
     nullify(cgrid%avg_rlongup             )
     nullify(cgrid%dmean_pcpg              )
     nullify(cgrid%dmean_runoff            )
     nullify(cgrid%dmean_drainage          )
     nullify(cgrid%dmean_vapor_ac          )
     nullify(cgrid%dmean_vapor_gc          )
     nullify(cgrid%dmean_vapor_lc          )
     nullify(cgrid%dmean_vapor_wc          )
     nullify(cgrid%dmean_nee               )
     nullify(cgrid%dmean_gpp               )
     nullify(cgrid%dmean_evap              )
     nullify(cgrid%dmean_transp            )
     nullify(cgrid%dmean_sensible_lc       )
     nullify(cgrid%dmean_sensible_wc       )
     nullify(cgrid%dmean_sensible_gc       )
     nullify(cgrid%dmean_sensible_ac       )
     nullify(cgrid%dmean_plresp            )
     nullify(cgrid%dmean_rh                )
     nullify(cgrid%dmean_leaf_resp         )
     nullify(cgrid%dmean_root_resp         )
     nullify(cgrid%dmean_growth_resp       )
     nullify(cgrid%dmean_storage_resp      )
     nullify(cgrid%dmean_vleaf_resp        )
     nullify(cgrid%dmean_nep               )
     nullify(cgrid%dmean_soil_temp         )
     nullify(cgrid%dmean_soil_water        )
     nullify(cgrid%dmean_soil_mstpot       )
     nullify(cgrid%dmean_transloss         )
     nullify(cgrid%dmean_fs_open           )
     nullify(cgrid%dmean_fsw               )
     nullify(cgrid%dmean_fsn               )
     nullify(cgrid%dmean_gpp_dbh           )
     nullify(cgrid%dmean_can_temp          )
     nullify(cgrid%dmean_can_shv           )
     nullify(cgrid%dmean_can_co2           )
     nullify(cgrid%dmean_can_rhos          )
     nullify(cgrid%dmean_can_prss          )
     nullify(cgrid%dmean_can_theta         )
     nullify(cgrid%dmean_can_theiv         )
     nullify(cgrid%dmean_gnd_temp          )
     nullify(cgrid%dmean_gnd_shv           )
     nullify(cgrid%dmean_leaf_energy       )
     nullify(cgrid%dmean_leaf_water        )
     nullify(cgrid%dmean_leaf_hcap         )
     nullify(cgrid%dmean_leaf_temp         )
     nullify(cgrid%dmean_wood_energy       )
     nullify(cgrid%dmean_wood_water        )
     nullify(cgrid%dmean_wood_hcap         )
     nullify(cgrid%dmean_wood_temp         )
     nullify(cgrid%dmean_atm_temp          )
     nullify(cgrid%dmean_atm_shv           )
     nullify(cgrid%dmean_atm_co2           )
     nullify(cgrid%dmean_atm_prss          )
     nullify(cgrid%dmean_atm_vels          )
     nullify(cgrid%dmean_rshort            )
     nullify(cgrid%dmean_rshort_diff       )
     nullify(cgrid%dmean_rlong             )
     nullify(cgrid%dmean_rshort_gnd        )
     nullify(cgrid%dmean_rlong_gnd         )
     nullify(cgrid%dmean_albedo            )
     nullify(cgrid%dmean_albedo_beam       )
     nullify(cgrid%dmean_albedo_diffuse    )
     nullify(cgrid%dmean_rlong_albedo      )
     nullify(cgrid%dmean_rlongup           )
     nullify(cgrid%dmean_co2_residual      )
     nullify(cgrid%dmean_energy_residual   )
     nullify(cgrid%dmean_water_residual    )
     nullify(cgrid%lai_pft                 )
     nullify(cgrid%wpa_pft                 )
     nullify(cgrid%wai_pft                 )
     nullify(cgrid%mmean_nee               )
     nullify(cgrid%mmean_gpp               )
     nullify(cgrid%mmean_evap              )
     nullify(cgrid%mmean_transp            )
     nullify(cgrid%mmean_sensible_lc       )
     nullify(cgrid%mmean_sensible_wc       )
     nullify(cgrid%mmean_sensible_gc       )
     nullify(cgrid%mmean_sensible_ac       )
     nullify(cgrid%mmean_vapor_lc          )
     nullify(cgrid%mmean_vapor_wc          )
     nullify(cgrid%mmean_vapor_gc          )
     nullify(cgrid%mmean_vapor_ac          )
     nullify(cgrid%mmean_nep               )
     nullify(cgrid%mmean_soil_temp         )
     nullify(cgrid%mmean_soil_water        )
     nullify(cgrid%mmean_soil_mstpot       )
     nullify(cgrid%mmean_transloss         )
     nullify(cgrid%mmean_fs_open           )
     nullify(cgrid%mmean_fsw               )
     nullify(cgrid%mmean_fsn               )
     nullify(cgrid%mmean_plresp            )
     nullify(cgrid%mmean_rh                )
     nullify(cgrid%mmean_leaf_resp         )
     nullify(cgrid%mmean_root_resp         )
     nullify(cgrid%mmean_growth_resp       )
     nullify(cgrid%mmean_storage_resp      )
     nullify(cgrid%mmean_vleaf_resp        )
     nullify(cgrid%mmean_gpp_dbh           )
     nullify(cgrid%mmean_lai_pft           )
     nullify(cgrid%mmean_wpa_pft           )
     nullify(cgrid%mmean_wai_pft           )
     nullify(cgrid%mmean_can_temp          )
     nullify(cgrid%mmean_can_shv           )
     nullify(cgrid%mmean_can_co2           )
     nullify(cgrid%mmean_can_rhos          )
     nullify(cgrid%mmean_can_prss          )
     nullify(cgrid%mmean_can_theta         )
     nullify(cgrid%mmean_can_theiv         )
     nullify(cgrid%mmean_gnd_temp          )
     nullify(cgrid%mmean_gnd_shv           )
     nullify(cgrid%mmean_leaf_energy       )
     nullify(cgrid%mmean_leaf_water        )
     nullify(cgrid%mmean_leaf_hcap         )
     nullify(cgrid%mmean_leaf_temp         )
     nullify(cgrid%mmean_wood_energy       )
     nullify(cgrid%mmean_wood_water        )
     nullify(cgrid%mmean_wood_hcap         )
     nullify(cgrid%mmean_wood_temp         )
     nullify(cgrid%mmean_atm_temp          )
     nullify(cgrid%mmean_rshort            )
     nullify(cgrid%mmean_rshort_diff       )
     nullify(cgrid%mmean_rlong             )
     nullify(cgrid%mmean_rshort_gnd        )
     nullify(cgrid%mmean_rlong_gnd         )
     nullify(cgrid%mmean_albedo            )
     nullify(cgrid%mmean_albedo_beam       )
     nullify(cgrid%mmean_albedo_diffuse    )
     nullify(cgrid%mmean_rlong_albedo      )
     nullify(cgrid%mmean_rlongup           )
     nullify(cgrid%mmean_atm_shv           )
     nullify(cgrid%mmean_atm_co2           )
     nullify(cgrid%mmean_atm_prss          )
     nullify(cgrid%mmean_atm_vels          )
     nullify(cgrid%mmean_pcpg              )
     nullify(cgrid%mmean_runoff            )
     nullify(cgrid%mmean_drainage          )
     nullify(cgrid%mmean_co2_residual      )
     nullify(cgrid%mmean_energy_residual   )
     nullify(cgrid%mmean_water_residual    )
     nullify(cgrid%bseeds_pft              )
     nullify(cgrid%agb_pft                 )
     nullify(cgrid%ba_pft                  )

     nullify(cgrid%mmsqu_gpp               )
     nullify(cgrid%mmsqu_leaf_resp         )
     nullify(cgrid%mmsqu_root_resp         )
     nullify(cgrid%mmsqu_plresp            )
     nullify(cgrid%mmsqu_nee               )
     nullify(cgrid%mmsqu_nep               )
     nullify(cgrid%mmsqu_rh                )
     nullify(cgrid%mmsqu_sensible_ac       )
     nullify(cgrid%mmsqu_sensible_lc       )
     nullify(cgrid%mmsqu_sensible_wc       )
     nullify(cgrid%mmsqu_sensible_gc       )
     nullify(cgrid%mmsqu_evap              )
     nullify(cgrid%mmsqu_transp            )
     nullify(cgrid%mmsqu_vapor_ac          )
     nullify(cgrid%mmsqu_vapor_lc          )
     nullify(cgrid%mmsqu_vapor_wc          )
     nullify(cgrid%mmsqu_vapor_gc          )

     nullify(cgrid%disturbance_rates       )
     nullify(cgrid%workload                )

     nullify(cgrid%qmean_pcpg           )
     nullify(cgrid%qmean_runoff         )
     nullify(cgrid%qmean_drainage       )
     nullify(cgrid%qmean_evap           )
     nullify(cgrid%qmean_transp         )
     nullify(cgrid%qmean_vapor_ac       )
     nullify(cgrid%qmean_vapor_gc       )
     nullify(cgrid%qmean_vapor_wc       )
     nullify(cgrid%qmean_vapor_lc       )
     nullify(cgrid%qmean_sensible_ac    )
     nullify(cgrid%qmean_sensible_gc    )
     nullify(cgrid%qmean_sensible_wc    )
     nullify(cgrid%qmean_sensible_lc    )
     nullify(cgrid%qmean_nee            )
     nullify(cgrid%qmean_gpp            )
     nullify(cgrid%qmean_nep            )
     nullify(cgrid%qmean_plresp         )
     nullify(cgrid%qmean_rh             )
     nullify(cgrid%qmean_leaf_resp      )
     nullify(cgrid%qmean_root_resp      )
     nullify(cgrid%qmean_soil_temp      )
     nullify(cgrid%qmean_soil_water     )
     nullify(cgrid%qmean_soil_mstpot    )
     nullify(cgrid%qmean_fs_open        )
     nullify(cgrid%qmean_fsw            )
     nullify(cgrid%qmean_fsn            )
     nullify(cgrid%qmean_can_temp       )
     nullify(cgrid%qmean_can_shv        )
     nullify(cgrid%qmean_can_co2        )
     nullify(cgrid%qmean_can_rhos       )
     nullify(cgrid%qmean_can_prss       )
     nullify(cgrid%qmean_can_theta      )
     nullify(cgrid%qmean_can_theiv      )
     nullify(cgrid%qmean_gnd_temp       )
     nullify(cgrid%qmean_gnd_shv        )
     nullify(cgrid%qmean_leaf_energy    )
     nullify(cgrid%qmean_leaf_water     )
     nullify(cgrid%qmean_leaf_temp      )
     nullify(cgrid%qmean_leaf_hcap      )
     nullify(cgrid%qmean_wood_energy    )
     nullify(cgrid%qmean_wood_water     )
     nullify(cgrid%qmean_wood_temp      )
     nullify(cgrid%qmean_wood_hcap      )
     nullify(cgrid%qmean_atm_temp       )
     nullify(cgrid%qmean_rshort         )
     nullify(cgrid%qmean_rshort_diff    )
     nullify(cgrid%qmean_rlong          )
     nullify(cgrid%qmean_rshort_gnd     )
     nullify(cgrid%qmean_rlong_gnd      )
     nullify(cgrid%qmean_albedo         )
     nullify(cgrid%qmean_albedo_beam    )
     nullify(cgrid%qmean_albedo_diffuse )
     nullify(cgrid%qmean_rlong_albedo   )
     nullify(cgrid%qmean_rlongup        )
     nullify(cgrid%qmean_atm_shv        )
     nullify(cgrid%qmean_atm_co2        )
     nullify(cgrid%qmean_atm_prss       )
     nullify(cgrid%qmean_atm_vels       )
     nullify(cgrid%qmsqu_gpp            )
     nullify(cgrid%qmsqu_leaf_resp      )
     nullify(cgrid%qmsqu_root_resp      )
     nullify(cgrid%qmsqu_plresp         )
     nullify(cgrid%qmsqu_nee            )
     nullify(cgrid%qmsqu_nep            )
     nullify(cgrid%qmsqu_rh             )
     nullify(cgrid%qmsqu_sensible_ac    )
     nullify(cgrid%qmsqu_sensible_lc    )
     nullify(cgrid%qmsqu_sensible_wc    )
     nullify(cgrid%qmsqu_sensible_gc    )
     nullify(cgrid%qmsqu_evap           )
     nullify(cgrid%qmsqu_transp         )
     nullify(cgrid%qmsqu_vapor_ac       )
     nullify(cgrid%qmsqu_vapor_lc       )
     nullify(cgrid%qmsqu_vapor_wc       )
     nullify(cgrid%qmsqu_vapor_gc       )

    return
  end subroutine nullify_edtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine nullify_polygontype(cpoly)

    implicit none

    type(polygontype),target :: cpoly

    nullify(cpoly%sipa_id)
    nullify(cpoly%sipa_n)
    nullify(cpoly%patch_count)  
    nullify(cpoly%site)
    nullify(cpoly%sitenum)

    nullify(cpoly%area)
    nullify(cpoly%patch_area)
    nullify(cpoly%elevation)
    nullify(cpoly%slope)
    nullify(cpoly%aspect)

    nullify(cpoly%num_landuse_years)
    nullify(cpoly%mindbh_primary)
    nullify(cpoly%probharv_primary)
    nullify(cpoly%mindbh_secondary)
    nullify(cpoly%probharv_secondary)


    nullify(cpoly%lai_pft)
    nullify(cpoly%wpa_pft)
    nullify(cpoly%wai_pft)

    nullify(cpoly%TCI)   
    nullify(cpoly%pptweight)   
    nullify(cpoly%lsl)   
    nullify(cpoly%hydro_next)
    nullify(cpoly%hydro_prev)
    nullify(cpoly%moist_W)
    nullify(cpoly%moist_f)  
    nullify(cpoly%moist_tau)
    nullify(cpoly%moist_zi) 
    nullify(cpoly%baseflow) 
    nullify(cpoly%ntext_soil)
    nullify(cpoly%min_monthly_temp)
    nullify(cpoly%plantation) 
    nullify(cpoly%agri_stocking_pft)
    nullify(cpoly%agri_stocking_density)
    nullify(cpoly%plantation_stocking_pft)
    nullify(cpoly%plantation_stocking_density)
    nullify(cpoly%primary_harvest_memory)
    nullify(cpoly%secondary_harvest_memory)
    nullify(cpoly%fire_disturbance_rate)
    nullify(cpoly%ignition_rate)
    nullify(cpoly%lambda_fire)
    nullify(cpoly%phen_pars)
    nullify(cpoly%nat_disturbance_rate)
    nullify(cpoly%nat_dist_type)
    nullify(cpoly%disturbance_memory)
    nullify(cpoly%disturbance_rates)
    nullify(cpoly%loss_fraction)
    nullify(cpoly%green_leaf_factor)
    nullify(cpoly%leaf_aging_factor)
    nullify(cpoly%met)
    nullify(cpoly%basal_area)
    nullify(cpoly%agb)    
    nullify(cpoly%pldens)    
    nullify(cpoly%bseeds)    
    
    nullify(cpoly%basal_area_growth)
    nullify(cpoly%agb_growth       )
    nullify(cpoly%basal_area_mort  )
    nullify(cpoly%agb_mort         )
    nullify(cpoly%basal_area_cut   ) !NOT IN REGISTRY
    nullify(cpoly%agb_cut          )

    nullify(cpoly%cosaoi)
    nullify(cpoly%avg_albedo)
    nullify(cpoly%avg_albedo_beam)
    nullify(cpoly%avg_albedo_diffuse)
    nullify(cpoly%avg_rlong_albedo)
    nullify(cpoly%avg_rlongup)

    nullify(cpoly%avg_lma    )
    nullify(cpoly%daylight)
    nullify(cpoly%lai    )
    nullify(cpoly%wpa    )
    nullify(cpoly%wai    )
    
    ! Fast time flux diagnostics
    ! ---------------------------------------------
    nullify(cpoly%avg_vapor_lc  )
    nullify(cpoly%avg_vapor_wc  )
    nullify(cpoly%avg_dew_cg    )
    nullify(cpoly%avg_vapor_gc  )
    nullify(cpoly%avg_wshed_vg  )
    nullify(cpoly%avg_intercepted)
    nullify(cpoly%avg_throughfall)
    nullify(cpoly%avg_vapor_ac  )
    nullify(cpoly%avg_transp    )
    nullify(cpoly%avg_evap      )
    nullify(cpoly%avg_smoist_gg )
    nullify(cpoly%avg_transloss )
    nullify(cpoly%avg_runoff    )
    nullify(cpoly%avg_drainage  )
    nullify(cpoly%avg_drainage_heat  )
    nullify(cpoly%aux           )
    nullify(cpoly%aux_s         )
    nullify(cpoly%avg_rshort_gnd   )
    nullify(cpoly%avg_rlong_gnd    )
    nullify(cpoly%avg_carbon_ac    )
    nullify(cpoly%avg_sensible_lc  )
    nullify(cpoly%avg_sensible_wc  )
    nullify(cpoly%avg_qwshed_vg    )
    nullify(cpoly%avg_qintercepted )
    nullify(cpoly%avg_qthroughfall )
    nullify(cpoly%avg_sensible_gc  )
    nullify(cpoly%avg_sensible_ac  )
    nullify(cpoly%avg_sensible_gg  )
    nullify(cpoly%avg_runoff_heat  )

    ! ----------------------------------------------
    nullify(cpoly%avg_leaf_energy)
    nullify(cpoly%avg_leaf_hcap  )
    nullify(cpoly%avg_leaf_temp  )
    nullify(cpoly%avg_leaf_fliq  )
    nullify(cpoly%avg_leaf_water )
    nullify(cpoly%avg_wood_energy)
    nullify(cpoly%avg_wood_hcap  )
    nullify(cpoly%avg_wood_temp  )
    nullify(cpoly%avg_wood_fliq  )
    nullify(cpoly%avg_wood_water )
    nullify(cpoly%avg_can_temp   )
    nullify(cpoly%avg_can_shv    )
    nullify(cpoly%avg_can_co2    )
    nullify(cpoly%avg_can_rhos   )
    nullify(cpoly%avg_can_prss   )
    nullify(cpoly%avg_can_theta  )
    nullify(cpoly%avg_can_theiv  )
    nullify(cpoly%avg_can_depth  )
    nullify(cpoly%avg_soil_energy)
    nullify(cpoly%avg_soil_mstpot)
    nullify(cpoly%avg_soil_water)
    nullify(cpoly%avg_soil_temp )
    nullify(cpoly%avg_soil_fracliq )
    nullify(cpoly%avg_soil_rootfrac)
    nullify(cpoly%avg_soil_wetness )
    nullify(cpoly%avg_skin_temp  )
    nullify(cpoly%avg_available_water)
    nullify(cpoly%runoff        )
    ! Phenology
    nullify(cpoly%rad_avg          )
    ! Meteorological conditions
    nullify(cpoly%avg_atm_tmp      )
    nullify(cpoly%avg_atm_shv      )
    nullify(cpoly%avg_atm_prss     )
    ! NACP
    nullify(cpoly%avg_sfcw_depth   )
    nullify(cpoly%avg_sfcw_energy  )
    nullify(cpoly%avg_sfcw_mass    )
    nullify(cpoly%avg_sfcw_tempk   )
    nullify(cpoly%avg_sfcw_fracliq )
    nullify(cpoly%avg_fsc          )
    nullify(cpoly%avg_stsc         )
    nullify(cpoly%avg_ssc          )
    nullify(cpoly%avg_balive       )
    nullify(cpoly%avg_bleaf        )
    nullify(cpoly%avg_broot        )
    nullify(cpoly%avg_bsapwood     )
    nullify(cpoly%avg_bdead        )
    nullify(cpoly%avg_fsn          )
    nullify(cpoly%avg_msn          )
    nullify(cpoly%avg_bstorage     )
    nullify(cpoly%avg_bseeds       )

    nullify(cpoly%dmean_co2_residual      )
    nullify(cpoly%dmean_energy_residual   )
    nullify(cpoly%dmean_water_residual    )
    nullify(cpoly%mmean_co2_residual      )
    nullify(cpoly%mmean_energy_residual   )
    nullify(cpoly%mmean_water_residual    )

    return
  end subroutine nullify_polygontype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine nullify_sitetype(csite)

    implicit none

    type(sitetype),target :: csite
    
    nullify(csite%paco_id)
    nullify(csite%paco_n)

    nullify(csite%dist_type)
    nullify(csite%age)
    nullify(csite%area)
    nullify(csite%laiarea)

    nullify(csite%pname)
    nullify(csite%sum_dgd)
    nullify(csite%sum_chd)
    nullify(csite%plantation) 
    nullify(csite%cohort_count)
    nullify(csite%can_theiv)
    nullify(csite%can_temp)
    nullify(csite%can_shv)
    nullify(csite%can_co2)
    nullify(csite%can_rhos)
    nullify(csite%can_prss)
    nullify(csite%can_theta)
    nullify(csite%can_depth)
    nullify(csite%opencan_frac)
    nullify(csite%ggbare)
    nullify(csite%ggveg)
    nullify(csite%ggnet)
    nullify(csite%ggsoil)
    nullify(csite%lambda_light)
    nullify(csite%dmean_lambda_light)
    nullify(csite%mmean_lambda_light)
    nullify(csite%lai)
    nullify(csite%wpa)
    nullify(csite%wai)

    nullify(csite%sfcwater_mass)
    nullify(csite%sfcwater_energy)
    nullify(csite%sfcwater_depth)
    nullify(csite%rshort_s)
    nullify(csite%rshort_s_beam)
    nullify(csite%rshort_s_diffuse)
    nullify(csite%sfcwater_tempk)
    nullify(csite%sfcwater_fracliq)
    nullify(csite%nlev_sfcwater)
    nullify(csite%soil_energy)
    nullify(csite%soil_water)
    nullify(csite%soil_tempk)
    nullify(csite%soil_fracliq)
    nullify(csite%rootdense)
    nullify(csite%ground_shv)
    nullify(csite%ground_ssh)
    nullify(csite%ground_temp)
    nullify(csite%ground_fliq)
    nullify(csite%rough)
    nullify(csite%par_l_max) 
    nullify(csite%par_l_beam_max) 
    nullify(csite%par_l_diffuse_max) 
    nullify(csite%A_o_max) 
    nullify(csite%A_c_max) 
    nullify(csite%old_stoma_data_max)
    nullify(csite%old_stoma_vector_max)
    nullify(csite%avg_daily_temp)  
    nullify(csite%mean_rh)
    nullify(csite%dmean_rh)
    nullify(csite%mmean_rh)
    nullify(csite%qmean_rh)
    nullify(csite%dmean_albedo)
    nullify(csite%dmean_albedo_beam)
    nullify(csite%dmean_albedo_diffuse)
    nullify(csite%mmean_albedo)
    nullify(csite%mmean_albedo_beam)
    nullify(csite%mmean_albedo_diffuse)
    nullify(csite%qmean_albedo)
    nullify(csite%qmean_albedo_beam)
    nullify(csite%qmean_albedo_diffuse)
    nullify(csite%mean_nep)
    nullify(csite%wbudget_loss2atm)
    nullify(csite%wbudget_denseffect)
    nullify(csite%wbudget_precipgain)
    nullify(csite%wbudget_loss2runoff)
    nullify(csite%wbudget_loss2drainage)
    nullify(csite%wbudget_initialstorage)
    nullify(csite%wbudget_residual)
    nullify(csite%ebudget_loss2atm)
    nullify(csite%ebudget_denseffect)
    nullify(csite%ebudget_loss2runoff)
    nullify(csite%ebudget_loss2drainage)
    nullify(csite%ebudget_netrad)
    nullify(csite%ebudget_precipgain)
    nullify(csite%ebudget_initialstorage)
    nullify(csite%ebudget_residual)
    nullify(csite%co2budget_initialstorage)
    nullify(csite%co2budget_residual)
    nullify(csite%co2budget_loss2atm)
    nullify(csite%co2budget_denseffect)
    nullify(csite%co2budget_gpp)
    nullify(csite%co2budget_gpp_dbh)
    nullify(csite%co2budget_plresp)
    nullify(csite%co2budget_rh)
    nullify(csite%today_A_decomp)
    nullify(csite%today_Af_decomp)
    nullify(csite%dmean_A_decomp)
    nullify(csite%dmean_Af_decomp)
    nullify(csite%mmean_A_decomp)
    nullify(csite%mmean_Af_decomp)
    nullify(csite%repro)
    nullify(csite%veg_rough)
    nullify(csite%veg_height)
    nullify(csite%veg_displace)
    nullify(csite%fsc_in)
    nullify(csite%ssc_in)
    nullify(csite%ssl_in)
    nullify(csite%fsn_in)


    nullify(csite%fsp_in)
    nullify(csite%stsp_in)
    nullify(csite%total_plant_P_uptake)


    nullify(csite%total_plant_nitrogen_uptake)
    nullify(csite%mineralized_N_loss)
    nullify(csite%mineralized_N_input)
    nullify(csite%rshort_g)
    nullify(csite%rshort_g_beam)
    nullify(csite%rshort_g_diffuse)
    nullify(csite%rlong_g)
    nullify(csite%rlong_g_surf)
    nullify(csite%rlong_g_incid)
    nullify(csite%rlong_s)
    nullify(csite%rlong_s_surf)
    nullify(csite%rlong_s_incid)
    nullify(csite%albedo)
    nullify(csite%albedo_beam)
    nullify(csite%albedo_diffuse)
    nullify(csite%rlongup)
    nullify(csite%rlong_albedo)
    nullify(csite%total_sfcw_depth)
    nullify(csite%snowfac)
    nullify(csite%A_decomp)
    nullify(csite%f_decomp)
    nullify(csite%rh)
    nullify(csite%cwd_rh)
    nullify(csite%fuse_flag)
    nullify(csite%cumlai_profile)
    nullify(csite%plant_ag_biomass)

    nullify(csite%mean_wflux)
    nullify(csite%mean_latflux)
    nullify(csite%mean_hflux)
    nullify(csite%mean_runoff)
    nullify(csite%mean_qrunoff)

    nullify(csite%htry)
    nullify(csite%avg_rk4step)
    nullify(csite%dmean_rk4step)
    nullify(csite%mmean_rk4step)

    nullify(csite%avg_available_water)

    nullify(csite%ustar)
    nullify(csite%tstar)
    nullify(csite%qstar)
    nullify(csite%cstar)

    nullify(csite%zeta)
    nullify(csite%ribulk)

    nullify(csite%upwp)
    nullify(csite%qpwp)
    nullify(csite%cpwp)
    nullify(csite%tpwp)
    nullify(csite%wpwp)


    nullify(csite%avg_rshort_gnd    )
    nullify(csite%avg_rlong_gnd     )
    nullify(csite%avg_carbon_ac     )
    nullify(csite%avg_rlongup       )
    nullify(csite%avg_albedo        )
    nullify(csite%avg_albedo_beam   )
    nullify(csite%avg_albedo_diffuse)
    nullify(csite%avg_rlong_albedo  )

    ! Fast time flux diagnostics
    ! ---------------------------------------------
    nullify(csite%avg_vapor_lc  )
    nullify(csite%avg_vapor_wc  )
    nullify(csite%avg_dew_cg    )
    nullify(csite%avg_vapor_gc  )
    nullify(csite%avg_wshed_vg  )
    nullify(csite%avg_intercepted)
    nullify(csite%avg_throughfall)
    nullify(csite%avg_vapor_ac  )
    nullify(csite%avg_transp    )
    nullify(csite%avg_evap      )
    nullify(csite%avg_smoist_gg )
    nullify(csite%avg_transloss )
    nullify(csite%avg_runoff    )
    nullify(csite%avg_drainage  )
    nullify(csite%avg_drainage_heat  )
    nullify(csite%aux           )
    nullify(csite%aux_s         )
    nullify(csite%avg_sensible_lc  )
    nullify(csite%avg_sensible_wc  )
    nullify(csite%avg_qwshed_vg    )
    nullify(csite%avg_qintercepted )
    nullify(csite%avg_qthroughfall )
    nullify(csite%avg_sensible_gc  )
    nullify(csite%avg_sensible_ac  )
    nullify(csite%avg_sensible_gg  )
    nullify(csite%avg_runoff_heat  )

    ! ----------------------------------------------
    nullify(csite%avg_leaf_energy)
    nullify(csite%avg_leaf_temp  )
    nullify(csite%avg_leaf_hcap  )
    nullify(csite%avg_leaf_fliq  )
    nullify(csite%avg_leaf_water )
    nullify(csite%avg_wood_energy)
    nullify(csite%avg_wood_temp  )
    nullify(csite%avg_wood_hcap  )
    nullify(csite%avg_wood_fliq  )
    nullify(csite%avg_wood_water )

    nullify(csite%watertable      )
    nullify(csite%moist_dz        )
    nullify(csite%ksat            )
    nullify(csite%soil_sat_energy )
    nullify(csite%soil_sat_water  )
    nullify(csite%soil_sat_heat   )
    nullify(csite%runoff_A        )
    nullify(csite%runoff_rate     )
    nullify(csite%runoff          )

    nullify(csite%patch)

    nullify(csite%dmean_co2_residual      )
    nullify(csite%dmean_energy_residual   )
    nullify(csite%dmean_water_residual    )
    nullify(csite%mmean_co2_residual      )
    nullify(csite%mmean_energy_residual   )
    nullify(csite%mmean_water_residual    )

    nullify(csite%current_paw)
    nullify(csite%past_paw)

    return
  end subroutine nullify_sitetype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine deallocate_edtype(cgrid)
  ! ===================================================
  !  Deallocation subroutines of ed state variables
  ! ===================================================
    
    implicit none
    
    type(edtype),target :: cgrid

       if(associated(cgrid%polygon                 )) deallocate(cgrid%polygon                 )
       if(associated(cgrid%lat                     )) deallocate(cgrid%lat                     )
       if(associated(cgrid%lon                     )) deallocate(cgrid%lon                     )
       if(associated(cgrid%xatm                    )) deallocate(cgrid%xatm                    )
       if(associated(cgrid%yatm                    )) deallocate(cgrid%yatm                    )
       if(associated(cgrid%ntext_soil              )) deallocate(cgrid%ntext_soil              )
       if(associated(cgrid%lsl                     )) deallocate(cgrid%lsl                     )
       
       if(associated(cgrid%pysi_n                  )) deallocate(cgrid%pysi_n                  )
       if(associated(cgrid%pysi_id                 )) deallocate(cgrid%pysi_id                 )
       if(associated(cgrid%walltime_py             )) deallocate(cgrid%walltime_py             )
       if(associated(cgrid%sensflux_py             )) deallocate(cgrid%sensflux_py             )
       if(associated(cgrid%site_adjacency          )) deallocate(cgrid%site_adjacency          )
       if(associated(cgrid%wbar                    )) deallocate(cgrid%wbar                    )
       if(associated(cgrid%Te                      )) deallocate(cgrid%Te                      )
       if(associated(cgrid%zbar                    )) deallocate(cgrid%zbar                    )
       if(associated(cgrid%sheat                   )) deallocate(cgrid%sheat                   )
       if(associated(cgrid%baseflow                )) deallocate(cgrid%baseflow                )
       if(associated(cgrid%runoff                  )) deallocate(cgrid%runoff                  )
       if(associated(cgrid%swliq                   )) deallocate(cgrid%swliq                   )
       
       if(associated(cgrid%ilon                    )) deallocate(cgrid%ilon                    )
       if(associated(cgrid%ilat                    )) deallocate(cgrid%ilat                    )
       if(associated(cgrid%total_agb               )) deallocate(cgrid%total_agb               )
       if(associated(cgrid%total_basal_area        )) deallocate(cgrid%total_basal_area        )
       if(associated(cgrid%total_agb_growth        )) deallocate(cgrid%total_agb_growth        )
       if(associated(cgrid%total_agb_mort          )) deallocate(cgrid%total_agb_mort          )
       if(associated(cgrid%total_agb_recruit       )) deallocate(cgrid%total_agb_recruit       )
       if(associated(cgrid%total_basal_area_growth )) deallocate(cgrid%total_basal_area_growth )
       if(associated(cgrid%total_basal_area_mort   )) deallocate(cgrid%total_basal_area_mort   )
       if(associated(cgrid%total_basal_area_recruit)) deallocate(cgrid%total_basal_area_recruit)
       if(associated(cgrid%nsites                  )) deallocate(cgrid%nsites                  )
       if(associated(cgrid%sitenums                )) deallocate(cgrid%sitenums                )
       if(associated(cgrid%load_adjacency          )) deallocate(cgrid%load_adjacency          )
       if(associated(cgrid%cosz                    )) deallocate(cgrid%cosz                    )
       if(associated(cgrid%mean_gpp                )) deallocate(cgrid%mean_gpp                )
       if(associated(cgrid%mean_precip             )) deallocate(cgrid%mean_precip             )
       if(associated(cgrid%mean_qprecip            )) deallocate(cgrid%mean_qprecip            )
       if(associated(cgrid%mean_netrad             )) deallocate(cgrid%mean_netrad             )
       if(associated(cgrid%cbudget_initialstorage  )) deallocate(cgrid%cbudget_initialstorage  )
       if(associated(cgrid%cbudget_nep             )) deallocate(cgrid%cbudget_nep             )
       if(associated(cgrid%nbudget_initialstorage  )) deallocate(cgrid%nbudget_initialstorage  )
       if(associated(cgrid%basal_area              )) deallocate(cgrid%basal_area              )
       if(associated(cgrid%agb                     )) deallocate(cgrid%agb                     )
       if(associated(cgrid%pldens                  )) deallocate(cgrid%pldens                  )
       if(associated(cgrid%bseeds                  )) deallocate(cgrid%bseeds                  )
       
       if(associated(cgrid%metinput                )) deallocate(cgrid%metinput                )
       if(associated(cgrid%met                     )) deallocate(cgrid%met                     )
       if(associated(cgrid%lapse                   )) deallocate(cgrid%lapse                   )

       if(associated(cgrid%lai                     )) deallocate(cgrid%lai                     ) 
       if(associated(cgrid%avg_lma                 )) deallocate(cgrid%avg_lma                    )
       if(associated(cgrid%wpa                     )) deallocate(cgrid%wpa                     )
       if(associated(cgrid%wai                     )) deallocate(cgrid%wai                     )

       ! Fast time flux diagnostics
       ! ---------------------------------------------
       if(associated(cgrid%avg_vapor_lc            )) deallocate(cgrid%avg_vapor_lc            )
       if(associated(cgrid%avg_vapor_wc            )) deallocate(cgrid%avg_vapor_wc            )
       if(associated(cgrid%avg_dew_cg              )) deallocate(cgrid%avg_dew_cg              )
       if(associated(cgrid%avg_vapor_gc            )) deallocate(cgrid%avg_vapor_gc            )
       if(associated(cgrid%avg_wshed_vg            )) deallocate(cgrid%avg_wshed_vg            )
       if(associated(cgrid%avg_intercepted         )) deallocate(cgrid%avg_intercepted         )
       if(associated(cgrid%avg_throughfall         )) deallocate(cgrid%avg_throughfall         )
       if(associated(cgrid%avg_vapor_ac            )) deallocate(cgrid%avg_vapor_ac            )
       if(associated(cgrid%avg_transp              )) deallocate(cgrid%avg_transp              )
       if(associated(cgrid%avg_evap                )) deallocate(cgrid%avg_evap                )
       if(associated(cgrid%avg_smoist_gg           )) deallocate(cgrid%avg_smoist_gg           )
       if(associated(cgrid%avg_transloss           )) deallocate(cgrid%avg_transloss           )
       if(associated(cgrid%avg_runoff              )) deallocate(cgrid%avg_runoff              )
       if(associated(cgrid%avg_drainage            )) deallocate(cgrid%avg_drainage            )
       if(associated(cgrid%avg_drainage_heat       )) deallocate(cgrid%avg_drainage_heat       )
       if(associated(cgrid%aux                     )) deallocate(cgrid%aux                     )
       if(associated(cgrid%aux_s                   )) deallocate(cgrid%aux_s                   )
       if(associated(cgrid%avg_rshort_gnd          )) deallocate(cgrid%avg_rshort_gnd          )
       if(associated(cgrid%avg_rlong_gnd           )) deallocate(cgrid%avg_rlong_gnd           )
       if(associated(cgrid%avg_carbon_ac           )) deallocate(cgrid%avg_carbon_ac           )
       if(associated(cgrid%avg_sensible_lc         )) deallocate(cgrid%avg_sensible_lc         )
       if(associated(cgrid%avg_sensible_wc         )) deallocate(cgrid%avg_sensible_wc         )
       if(associated(cgrid%avg_qwshed_vg           )) deallocate(cgrid%avg_qwshed_vg           )
       if(associated(cgrid%avg_qintercepted        )) deallocate(cgrid%avg_qintercepted        )
       if(associated(cgrid%avg_qthroughfall        )) deallocate(cgrid%avg_qthroughfall        )
       if(associated(cgrid%avg_sensible_gc         )) deallocate(cgrid%avg_sensible_gc         )
       if(associated(cgrid%avg_sensible_ac         )) deallocate(cgrid%avg_sensible_ac         )
       if(associated(cgrid%avg_sensible_gg         )) deallocate(cgrid%avg_sensible_gg         )
       if(associated(cgrid%avg_runoff_heat         )) deallocate(cgrid%avg_runoff_heat         )

       if(associated(cgrid%max_leaf_temp          )) deallocate(cgrid%max_leaf_temp          )
       if(associated(cgrid%min_leaf_temp          )) deallocate(cgrid%min_leaf_temp          )
       if(associated(cgrid%max_wood_temp          )) deallocate(cgrid%max_wood_temp          )
       if(associated(cgrid%min_wood_temp          )) deallocate(cgrid%min_wood_temp          )
       if(associated(cgrid%max_soil_temp          )) deallocate(cgrid%max_soil_temp          )
       if(associated(cgrid%min_soil_temp          )) deallocate(cgrid%min_soil_temp          )

       ! Fast time state diagnostics
       if(associated(cgrid%avg_leaf_energy         )) deallocate(cgrid%avg_leaf_energy         )
       if(associated(cgrid%avg_leaf_hcap           )) deallocate(cgrid%avg_leaf_hcap           )
       if(associated(cgrid%avg_leaf_temp           )) deallocate(cgrid%avg_leaf_temp           )
       if(associated(cgrid%avg_leaf_fliq           )) deallocate(cgrid%avg_leaf_fliq           )
       if(associated(cgrid%avg_leaf_water          )) deallocate(cgrid%avg_leaf_water          )
       if(associated(cgrid%avg_wood_energy         )) deallocate(cgrid%avg_wood_energy         )
       if(associated(cgrid%avg_wood_hcap           )) deallocate(cgrid%avg_wood_hcap           )
       if(associated(cgrid%avg_wood_temp           )) deallocate(cgrid%avg_wood_temp           )
       if(associated(cgrid%avg_wood_fliq           )) deallocate(cgrid%avg_wood_fliq           )
       if(associated(cgrid%avg_wood_water          )) deallocate(cgrid%avg_wood_water          )
       if(associated(cgrid%avg_can_temp            )) deallocate(cgrid%avg_can_temp            )
       if(associated(cgrid%avg_can_shv             )) deallocate(cgrid%avg_can_shv             )
       if(associated(cgrid%avg_can_co2             )) deallocate(cgrid%avg_can_co2             )
       if(associated(cgrid%avg_can_rhos            )) deallocate(cgrid%avg_can_rhos            )
       if(associated(cgrid%avg_can_prss            )) deallocate(cgrid%avg_can_prss            )
       if(associated(cgrid%avg_can_theta           )) deallocate(cgrid%avg_can_theta           )
       if(associated(cgrid%avg_can_theiv           )) deallocate(cgrid%avg_can_theiv           )
       if(associated(cgrid%avg_can_depth           )) deallocate(cgrid%avg_can_depth           )
       if(associated(cgrid%avg_soil_energy         )) deallocate(cgrid%avg_soil_energy         )
       if(associated(cgrid%avg_soil_mstpot         )) deallocate(cgrid%avg_soil_mstpot         )
       if(associated(cgrid%avg_soil_water          )) deallocate(cgrid%avg_soil_water          )
       if(associated(cgrid%avg_soil_temp           )) deallocate(cgrid%avg_soil_temp           )
       if(associated(cgrid%avg_soil_fracliq        )) deallocate(cgrid%avg_soil_fracliq        )
       if(associated(cgrid%avg_soil_rootfrac       )) deallocate(cgrid%avg_soil_rootfrac       )
       if(associated(cgrid%avg_soil_wetness        )) deallocate(cgrid%avg_soil_wetness        )
       if(associated(cgrid%avg_skin_temp           )) deallocate(cgrid%avg_skin_temp           )
       if(associated(cgrid%avg_available_water     )) deallocate(cgrid%avg_available_water     )
      
       if(associated(cgrid%avg_lai_ebalvars        )) deallocate(cgrid%avg_lai_ebalvars        )

       if(associated(cgrid%avg_gpp                 )) deallocate(cgrid%avg_gpp                 )
       if(associated(cgrid%avg_leaf_resp           )) deallocate(cgrid%avg_leaf_resp           )
       if(associated(cgrid%avg_root_resp           )) deallocate(cgrid%avg_root_resp           )
       if(associated(cgrid%avg_growth_resp         )) deallocate(cgrid%avg_growth_resp         )
       if(associated(cgrid%avg_storage_resp        )) deallocate(cgrid%avg_storage_resp        )
       if(associated(cgrid%avg_vleaf_resp          )) deallocate(cgrid%avg_vleaf_resp          )
       if(associated(cgrid%avg_plant_resp          )) deallocate(cgrid%avg_plant_resp          )
       if(associated(cgrid%avg_htroph_resp         )) deallocate(cgrid%avg_htroph_resp         )
       if(associated(cgrid%avg_leaf_drop           )) deallocate(cgrid%avg_leaf_drop           )
       if(associated(cgrid%avg_leaf_maintenance    )) deallocate(cgrid%avg_leaf_maintenance    )
       if(associated(cgrid%avg_root_maintenance    )) deallocate(cgrid%avg_root_maintenance    )

       !!! added for NACP intercomparison (MCD)
       if(associated(cgrid%avg_sfcw_depth          )) deallocate(cgrid%avg_sfcw_depth      )
       if(associated(cgrid%avg_sfcw_energy         )) deallocate(cgrid%avg_sfcw_energy     )
       if(associated(cgrid%avg_sfcw_mass           )) deallocate(cgrid%avg_sfcw_mass       )
       if(associated(cgrid%avg_sfcw_tempk          )) deallocate(cgrid%avg_sfcw_tempk      )
       if(associated(cgrid%avg_sfcw_fracliq        )) deallocate(cgrid%avg_sfcw_fracliq    )
       if(associated(cgrid%avg_bdead               )) deallocate(cgrid%avg_bdead           )
       if(associated(cgrid%avg_balive              )) deallocate(cgrid%avg_balive          )
       if(associated(cgrid%avg_broot               )) deallocate(cgrid%avg_broot           )
       if(associated(cgrid%avg_bleaf               )) deallocate(cgrid%avg_bleaf           )
       if(associated(cgrid%avg_bsapwood            )) deallocate(cgrid%avg_bsapwood        )
       if(associated(cgrid%avg_bstorage            )) deallocate(cgrid%avg_bstorage        )
       if(associated(cgrid%avg_bseeds              )) deallocate(cgrid%avg_bseeds          )
       if(associated(cgrid%avg_fsc                 )) deallocate(cgrid%avg_fsc             )
       if(associated(cgrid%avg_ssc                 )) deallocate(cgrid%avg_ssc             )
       if(associated(cgrid%avg_stsc                )) deallocate(cgrid%avg_stsc            )

       if(associated(cgrid%avg_fsn                 )) deallocate(cgrid%avg_fsn             )
       if(associated(cgrid%avg_msn                 )) deallocate(cgrid%avg_msn             )


     !-------- TOTAL CARBON AND NITROGEN POOLS  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
     if(associated(cgrid%Cleaf  )) deallocate(cgrid%Cleaf)
     if(associated(cgrid%Croot  )) deallocate(cgrid%Croot)
     if(associated(cgrid%Cstore )) deallocate(cgrid%Cstore)
     if(associated(cgrid%Ccwd   )) deallocate(cgrid%Ccwd)
     if(associated(cgrid%Nleaf  )) deallocate(cgrid%Nleaf)
     if(associated(cgrid%Ndead  )) deallocate(cgrid%Ndead)
     if(associated(cgrid%Nroot  )) deallocate(cgrid%Nroot)
     if(associated(cgrid%Nstore )) deallocate(cgrid%Nstore)
     if(associated(cgrid%Ncwd   )) deallocate(cgrid%Ncwd)
     
     
     !-------- TOTAL CARBON AND NITROGEN FLUX  ---------------
     ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
     if(associated(cgrid%Cleaf_grow       )) deallocate(cgrid%Cleaf_grow)
     if(associated(cgrid%Croot_grow       )) deallocate(cgrid%Croot_grow)
     if(associated(cgrid%Cdead_grow       )) deallocate(cgrid%Cdead_grow)
     if(associated(cgrid%Cstore_grow      )) deallocate(cgrid%Cstore_grow)
     if(associated(cgrid%Cleaf_litter_flux)) deallocate(cgrid%Cleaf_litter_flux)
     if(associated(cgrid%Croot_litter_flux)) deallocate(cgrid%Croot_litter_flux)
     if(associated(cgrid%Ccwd_flux        )) deallocate(cgrid%Ccwd_flux)
     if(associated(cgrid%Nleaf_grow       )) deallocate(cgrid%Nleaf_grow)
     if(associated(cgrid%Ndead_grow       )) deallocate(cgrid%Ndead_grow)
     if(associated(cgrid%Nroot_grow       )) deallocate(cgrid%Nroot_grow)
     if(associated(cgrid%Nstore_grow      )) deallocate(cgrid%Nstore_grow)
     if(associated(cgrid%Nleaf_litter_flux)) deallocate(cgrid%Nleaf_litter_flux)
     if(associated(cgrid%Nroot_litter_flux)) deallocate(cgrid%Nroot_litter_flux)
     if(associated(cgrid%Ncwd_flux        )) deallocate(cgrid%Ncwd_flux)
     if(associated(cgrid%Nbiomass_uptake  )) deallocate(cgrid%Nbiomass_uptake)
     if(associated(cgrid%Ngross_min       )) deallocate(cgrid%Ngross_min)
     if(associated(cgrid%Nnet_min         )) deallocate(cgrid%Nnet_min)


 


       ! ----------------------------------------------

       if(associated(cgrid%avg_nir_beam            )) deallocate(cgrid%avg_nir_beam            )
       if(associated(cgrid%avg_nir_diffuse         )) deallocate(cgrid%avg_nir_diffuse         )
       if(associated(cgrid%avg_par_beam            )) deallocate(cgrid%avg_par_beam            )
       if(associated(cgrid%avg_par_diffuse         )) deallocate(cgrid%avg_par_diffuse         )
       if(associated(cgrid%avg_atm_tmp             )) deallocate(cgrid%avg_atm_tmp             )
       if(associated(cgrid%avg_atm_shv             )) deallocate(cgrid%avg_atm_shv             )
       if(associated(cgrid%avg_rshort              )) deallocate(cgrid%avg_rshort              )
       if(associated(cgrid%avg_rshort_diffuse      )) deallocate(cgrid%avg_rshort_diffuse      )
       if(associated(cgrid%avg_rlong               )) deallocate(cgrid%avg_rlong               )
       if(associated(cgrid%avg_pcpg                )) deallocate(cgrid%avg_pcpg                )
       if(associated(cgrid%avg_qpcpg               )) deallocate(cgrid%avg_qpcpg               )
       if(associated(cgrid%avg_dpcpg               )) deallocate(cgrid%avg_dpcpg               )
       if(associated(cgrid%avg_vels                )) deallocate(cgrid%avg_vels                )
       if(associated(cgrid%avg_atm_prss            )) deallocate(cgrid%avg_atm_prss            )
       if(associated(cgrid%avg_exner               )) deallocate(cgrid%avg_exner               )
       if(associated(cgrid%avg_geoht               )) deallocate(cgrid%avg_geoht               )
       if(associated(cgrid%avg_atm_co2             )) deallocate(cgrid%avg_atm_co2             )
       if(associated(cgrid%avg_albedo              )) deallocate(cgrid%avg_albedo              )
       if(associated(cgrid%avg_albedo_beam         )) deallocate(cgrid%avg_albedo_beam         )
       if(associated(cgrid%avg_albedo_diffuse      )) deallocate(cgrid%avg_albedo_diffuse      )
       if(associated(cgrid%avg_rlong_albedo        )) deallocate(cgrid%avg_rlong_albedo        )
       if(associated(cgrid%avg_rlongup             )) deallocate(cgrid%avg_rlongup             )


       if(associated(cgrid%runoff                  )) deallocate(cgrid%runoff                  )
       
       if(associated(cgrid%dmean_pcpg              )) deallocate(cgrid%dmean_pcpg              )
       if(associated(cgrid%dmean_runoff            )) deallocate(cgrid%dmean_runoff            )
       if(associated(cgrid%dmean_drainage          )) deallocate(cgrid%dmean_drainage          )
       if(associated(cgrid%dmean_vapor_ac          )) deallocate(cgrid%dmean_vapor_ac          )
       if(associated(cgrid%dmean_vapor_gc          )) deallocate(cgrid%dmean_vapor_gc          )
       if(associated(cgrid%dmean_vapor_lc          )) deallocate(cgrid%dmean_vapor_lc          )
       if(associated(cgrid%dmean_vapor_wc          )) deallocate(cgrid%dmean_vapor_wc          )
       if(associated(cgrid%dmean_nee               )) deallocate(cgrid%dmean_nee               )
       if(associated(cgrid%dmean_gpp               )) deallocate(cgrid%dmean_gpp               )
       if(associated(cgrid%dmean_evap              )) deallocate(cgrid%dmean_evap              )
       if(associated(cgrid%dmean_transp            )) deallocate(cgrid%dmean_transp            )
       if(associated(cgrid%dmean_sensible_lc       )) deallocate(cgrid%dmean_sensible_lc       )
       if(associated(cgrid%dmean_sensible_wc       )) deallocate(cgrid%dmean_sensible_wc       )
       if(associated(cgrid%dmean_sensible_gc       )) deallocate(cgrid%dmean_sensible_gc       )
       if(associated(cgrid%dmean_sensible_ac       )) deallocate(cgrid%dmean_sensible_ac       )
       if(associated(cgrid%dmean_plresp            )) deallocate(cgrid%dmean_plresp            )
       if(associated(cgrid%dmean_rh                )) deallocate(cgrid%dmean_rh                )
       if(associated(cgrid%dmean_leaf_resp         )) deallocate(cgrid%dmean_leaf_resp         )
       if(associated(cgrid%dmean_root_resp         )) deallocate(cgrid%dmean_root_resp         )
       if(associated(cgrid%dmean_growth_resp       )) deallocate(cgrid%dmean_growth_resp       )
       if(associated(cgrid%dmean_storage_resp      )) deallocate(cgrid%dmean_storage_resp      )
       if(associated(cgrid%dmean_vleaf_resp        )) deallocate(cgrid%dmean_vleaf_resp        )
       if(associated(cgrid%dmean_nep               )) deallocate(cgrid%dmean_nep               )
       if(associated(cgrid%dmean_soil_temp         )) deallocate(cgrid%dmean_soil_temp         )
       if(associated(cgrid%dmean_soil_water        )) deallocate(cgrid%dmean_soil_water        )
       if(associated(cgrid%dmean_soil_mstpot       )) deallocate(cgrid%dmean_soil_mstpot       )
       if(associated(cgrid%dmean_transloss         )) deallocate(cgrid%dmean_transloss         )
       if(associated(cgrid%dmean_fs_open           )) deallocate(cgrid%dmean_fs_open           )
       if(associated(cgrid%dmean_fsw               )) deallocate(cgrid%dmean_fsw               )
       if(associated(cgrid%dmean_fsn               )) deallocate(cgrid%dmean_fsn               )
       if(associated(cgrid%dmean_gpp_dbh           )) deallocate(cgrid%dmean_gpp_dbh           )
       if(associated(cgrid%dmean_can_temp          )) deallocate(cgrid%dmean_can_temp          )
       if(associated(cgrid%dmean_can_shv           )) deallocate(cgrid%dmean_can_shv           )
       if(associated(cgrid%dmean_can_co2           )) deallocate(cgrid%dmean_can_co2           )
       if(associated(cgrid%dmean_can_rhos          )) deallocate(cgrid%dmean_can_rhos          )
       if(associated(cgrid%dmean_can_prss          )) deallocate(cgrid%dmean_can_prss          )
       if(associated(cgrid%dmean_can_theta         )) deallocate(cgrid%dmean_can_theta         )
       if(associated(cgrid%dmean_can_theiv         )) deallocate(cgrid%dmean_can_theiv         )
       if(associated(cgrid%dmean_gnd_temp          )) deallocate(cgrid%dmean_gnd_temp          )
       if(associated(cgrid%dmean_gnd_shv           )) deallocate(cgrid%dmean_gnd_shv           )
       if(associated(cgrid%dmean_leaf_energy       )) deallocate(cgrid%dmean_leaf_energy       )
       if(associated(cgrid%dmean_leaf_water        )) deallocate(cgrid%dmean_leaf_water        )
       if(associated(cgrid%dmean_leaf_hcap         )) deallocate(cgrid%dmean_leaf_hcap         )
       if(associated(cgrid%dmean_leaf_temp         )) deallocate(cgrid%dmean_leaf_temp         )
       if(associated(cgrid%dmean_wood_energy       )) deallocate(cgrid%dmean_wood_energy       )
       if(associated(cgrid%dmean_wood_water        )) deallocate(cgrid%dmean_wood_water        )
       if(associated(cgrid%dmean_wood_hcap         )) deallocate(cgrid%dmean_wood_hcap         )
       if(associated(cgrid%dmean_wood_temp         )) deallocate(cgrid%dmean_wood_temp         )
       if(associated(cgrid%dmean_atm_temp          )) deallocate(cgrid%dmean_atm_temp          )
       if(associated(cgrid%dmean_atm_shv           )) deallocate(cgrid%dmean_atm_shv           )
       if(associated(cgrid%dmean_atm_co2           )) deallocate(cgrid%dmean_atm_co2           )
       if(associated(cgrid%dmean_atm_prss          )) deallocate(cgrid%dmean_atm_prss          )
       if(associated(cgrid%dmean_atm_vels          )) deallocate(cgrid%dmean_atm_vels          )
       if(associated(cgrid%dmean_rshort            )) deallocate(cgrid%dmean_rshort            )
       if(associated(cgrid%dmean_rshort_diff       )) deallocate(cgrid%dmean_rshort_diff       )
       if(associated(cgrid%dmean_rlong             )) deallocate(cgrid%dmean_rlong             )
       if(associated(cgrid%dmean_rshort_gnd        )) deallocate(cgrid%dmean_rshort_gnd        )
       if(associated(cgrid%dmean_rlong_gnd         )) deallocate(cgrid%dmean_rlong_gnd         )
       if(associated(cgrid%dmean_albedo            )) deallocate(cgrid%dmean_albedo            )
       if(associated(cgrid%dmean_albedo_beam       )) deallocate(cgrid%dmean_albedo_beam       )
       if(associated(cgrid%dmean_albedo_diffuse    )) deallocate(cgrid%dmean_albedo_diffuse    )
       if(associated(cgrid%dmean_rlong_albedo      )) deallocate(cgrid%dmean_rlong_albedo      )
       if(associated(cgrid%dmean_rlongup           )) deallocate(cgrid%dmean_rlongup           )
       if(associated(cgrid%dmean_co2_residual      )) deallocate(cgrid%dmean_co2_residual      )
       if(associated(cgrid%dmean_energy_residual   )) deallocate(cgrid%dmean_energy_residual   )
       if(associated(cgrid%dmean_water_residual    )) deallocate(cgrid%dmean_water_residual    )

       if(associated(cgrid%lai_pft                 )) deallocate(cgrid%lai_pft                 )
       if(associated(cgrid%wpa_pft                 )) deallocate(cgrid%wpa_pft                 )
       if(associated(cgrid%wai_pft                 )) deallocate(cgrid%wai_pft                 )
       if(associated(cgrid%mmean_gpp               )) deallocate(cgrid%mmean_gpp               )
       if(associated(cgrid%mmean_evap              )) deallocate(cgrid%mmean_evap              )
       if(associated(cgrid%mmean_transp            )) deallocate(cgrid%mmean_transp            )
       if(associated(cgrid%mmean_vapor_lc          )) deallocate(cgrid%mmean_vapor_lc          )
       if(associated(cgrid%mmean_vapor_wc          )) deallocate(cgrid%mmean_vapor_wc          )
       if(associated(cgrid%mmean_vapor_gc          )) deallocate(cgrid%mmean_vapor_gc          )
       if(associated(cgrid%mmean_vapor_ac          )) deallocate(cgrid%mmean_vapor_ac          )
       if(associated(cgrid%mmean_sensible_lc       )) deallocate(cgrid%mmean_sensible_lc       )
       if(associated(cgrid%mmean_sensible_wc       )) deallocate(cgrid%mmean_sensible_wc       )
       if(associated(cgrid%mmean_sensible_gc       )) deallocate(cgrid%mmean_sensible_gc       )
       if(associated(cgrid%mmean_sensible_ac       )) deallocate(cgrid%mmean_sensible_ac       )
       if(associated(cgrid%mmean_nee               )) deallocate(cgrid%mmean_nee               )
       if(associated(cgrid%mmean_nep               )) deallocate(cgrid%mmean_nep               )
       if(associated(cgrid%mmean_soil_temp         )) deallocate(cgrid%mmean_soil_temp         )
       if(associated(cgrid%mmean_soil_water        )) deallocate(cgrid%mmean_soil_water        )
       if(associated(cgrid%mmean_soil_mstpot       )) deallocate(cgrid%mmean_soil_mstpot       )
       if(associated(cgrid%dmean_transloss         )) deallocate(cgrid%dmean_transloss         )
       if(associated(cgrid%mmean_fs_open           )) deallocate(cgrid%mmean_fs_open           )
       if(associated(cgrid%mmean_fsw               )) deallocate(cgrid%mmean_fsw               )
       if(associated(cgrid%mmean_fsn               )) deallocate(cgrid%mmean_fsn               )
       if(associated(cgrid%mmean_plresp            )) deallocate(cgrid%mmean_plresp            )
       if(associated(cgrid%mmean_rh                )) deallocate(cgrid%mmean_rh                )
       if(associated(cgrid%mmean_leaf_resp         )) deallocate(cgrid%mmean_leaf_resp         )
       if(associated(cgrid%mmean_root_resp         )) deallocate(cgrid%mmean_root_resp         )
       if(associated(cgrid%mmean_growth_resp       )) deallocate(cgrid%mmean_growth_resp       )
       if(associated(cgrid%mmean_storage_resp      )) deallocate(cgrid%mmean_storage_resp      )
       if(associated(cgrid%mmean_vleaf_resp        )) deallocate(cgrid%mmean_vleaf_resp        )
       if(associated(cgrid%mmean_gpp_dbh           )) deallocate(cgrid%mmean_gpp_dbh           )
       if(associated(cgrid%mmean_lai_pft           )) deallocate(cgrid%mmean_lai_pft           )
       if(associated(cgrid%mmean_wpa_pft           )) deallocate(cgrid%mmean_wpa_pft           )
       if(associated(cgrid%mmean_wai_pft           )) deallocate(cgrid%mmean_wai_pft           )
       if(associated(cgrid%mmean_can_temp          )) deallocate(cgrid%mmean_can_temp          )
       if(associated(cgrid%mmean_can_shv           )) deallocate(cgrid%mmean_can_shv           )
       if(associated(cgrid%mmean_can_co2           )) deallocate(cgrid%mmean_can_co2           )
       if(associated(cgrid%mmean_can_rhos          )) deallocate(cgrid%mmean_can_rhos          )
       if(associated(cgrid%mmean_can_prss          )) deallocate(cgrid%mmean_can_prss          )
       if(associated(cgrid%mmean_can_theta         )) deallocate(cgrid%mmean_can_theta         )
       if(associated(cgrid%mmean_can_theiv         )) deallocate(cgrid%mmean_can_theiv         )
       if(associated(cgrid%mmean_gnd_temp          )) deallocate(cgrid%mmean_gnd_temp          )
       if(associated(cgrid%mmean_gnd_shv           )) deallocate(cgrid%mmean_gnd_shv           )
       if(associated(cgrid%mmean_leaf_energy       )) deallocate(cgrid%mmean_leaf_energy       )
       if(associated(cgrid%mmean_leaf_water        )) deallocate(cgrid%mmean_leaf_water        )
       if(associated(cgrid%mmean_leaf_hcap         )) deallocate(cgrid%mmean_leaf_hcap         )
       if(associated(cgrid%mmean_leaf_temp         )) deallocate(cgrid%mmean_leaf_temp         )
       if(associated(cgrid%mmean_wood_energy       )) deallocate(cgrid%mmean_wood_energy       )
       if(associated(cgrid%mmean_wood_water        )) deallocate(cgrid%mmean_wood_water        )
       if(associated(cgrid%mmean_wood_hcap         )) deallocate(cgrid%mmean_wood_hcap         )
       if(associated(cgrid%mmean_wood_temp         )) deallocate(cgrid%mmean_wood_temp         )
       if(associated(cgrid%mmean_atm_temp          )) deallocate(cgrid%mmean_atm_temp          )
       if(associated(cgrid%mmean_rshort            )) deallocate(cgrid%mmean_rshort            )
       if(associated(cgrid%mmean_rshort_diff       )) deallocate(cgrid%mmean_rshort_diff       )
       if(associated(cgrid%mmean_rlong             )) deallocate(cgrid%mmean_rlong             )
       if(associated(cgrid%mmean_rshort_gnd        )) deallocate(cgrid%mmean_rshort_gnd        )
       if(associated(cgrid%mmean_rlong_gnd         )) deallocate(cgrid%mmean_rlong_gnd         )
       if(associated(cgrid%mmean_albedo            )) deallocate(cgrid%mmean_albedo            )
       if(associated(cgrid%mmean_albedo_beam       )) deallocate(cgrid%mmean_albedo_beam       )
       if(associated(cgrid%mmean_albedo_diffuse    )) deallocate(cgrid%mmean_albedo_diffuse    )
       if(associated(cgrid%mmean_rlong_albedo      )) deallocate(cgrid%mmean_rlong_albedo      )
       if(associated(cgrid%mmean_rlongup           )) deallocate(cgrid%mmean_rlongup           )
       if(associated(cgrid%mmean_atm_shv           )) deallocate(cgrid%mmean_atm_shv           )
       if(associated(cgrid%mmean_atm_co2           )) deallocate(cgrid%mmean_atm_co2           )
       if(associated(cgrid%mmean_atm_prss          )) deallocate(cgrid%mmean_atm_prss          )
       if(associated(cgrid%mmean_atm_vels          )) deallocate(cgrid%mmean_atm_vels          )
       if(associated(cgrid%mmean_pcpg              )) deallocate(cgrid%mmean_pcpg              )
       if(associated(cgrid%mmean_runoff            )) deallocate(cgrid%mmean_runoff            )
       if(associated(cgrid%mmean_drainage          )) deallocate(cgrid%mmean_drainage          )
       if(associated(cgrid%mmean_co2_residual      )) deallocate(cgrid%mmean_co2_residual      )
       if(associated(cgrid%mmean_energy_residual   )) deallocate(cgrid%mmean_energy_residual   )
       if(associated(cgrid%mmean_water_residual    )) deallocate(cgrid%mmean_water_residual    )
       if(associated(cgrid%bseeds_pft              )) deallocate(cgrid%bseeds_pft              )
       if(associated(cgrid%agb_pft                 )) deallocate(cgrid%agb_pft                 )
       if(associated(cgrid%ba_pft                  )) deallocate(cgrid%ba_pft                  )


       if(associated(cgrid%mmsqu_gpp               )) deallocate(cgrid%mmsqu_gpp               )
       if(associated(cgrid%mmsqu_leaf_resp         )) deallocate(cgrid%mmsqu_leaf_resp         )
       if(associated(cgrid%mmsqu_root_resp         )) deallocate(cgrid%mmsqu_root_resp         )
       if(associated(cgrid%mmsqu_plresp            )) deallocate(cgrid%mmsqu_plresp            )
       if(associated(cgrid%mmsqu_nee               )) deallocate(cgrid%mmsqu_nee               )
       if(associated(cgrid%mmsqu_nep               )) deallocate(cgrid%mmsqu_nep               )
       if(associated(cgrid%mmsqu_rh                )) deallocate(cgrid%mmsqu_rh                )
       if(associated(cgrid%mmsqu_sensible_ac       )) deallocate(cgrid%mmsqu_sensible_ac       )
       if(associated(cgrid%mmsqu_sensible_lc       )) deallocate(cgrid%mmsqu_sensible_lc       )
       if(associated(cgrid%mmsqu_sensible_wc       )) deallocate(cgrid%mmsqu_sensible_wc       )
       if(associated(cgrid%mmsqu_sensible_gc       )) deallocate(cgrid%mmsqu_sensible_gc       )
       if(associated(cgrid%mmsqu_evap              )) deallocate(cgrid%mmsqu_evap              )
       if(associated(cgrid%mmsqu_transp            )) deallocate(cgrid%mmsqu_transp            )
       if(associated(cgrid%mmsqu_vapor_ac          )) deallocate(cgrid%mmsqu_vapor_ac          )
       if(associated(cgrid%mmsqu_vapor_lc          )) deallocate(cgrid%mmsqu_vapor_lc          )
       if(associated(cgrid%mmsqu_vapor_wc          )) deallocate(cgrid%mmsqu_vapor_wc          )
       if(associated(cgrid%mmsqu_vapor_gc          )) deallocate(cgrid%mmsqu_vapor_gc          )

       if(associated(cgrid%disturbance_rates       )) deallocate(cgrid%disturbance_rates       )
       if(associated(cgrid%workload                )) deallocate(cgrid%workload                )

    if(associated(cgrid%qmean_pcpg           )) deallocate(cgrid%qmean_pcpg           )
    if(associated(cgrid%qmean_runoff         )) deallocate(cgrid%qmean_runoff         )
    if(associated(cgrid%qmean_drainage       )) deallocate(cgrid%qmean_drainage       )
    if(associated(cgrid%qmean_evap           )) deallocate(cgrid%qmean_evap           )
    if(associated(cgrid%qmean_transp         )) deallocate(cgrid%qmean_transp         )
    if(associated(cgrid%qmean_vapor_ac       )) deallocate(cgrid%qmean_vapor_ac       )
    if(associated(cgrid%qmean_vapor_gc       )) deallocate(cgrid%qmean_vapor_gc       )
    if(associated(cgrid%qmean_vapor_lc       )) deallocate(cgrid%qmean_vapor_lc       )
    if(associated(cgrid%qmean_vapor_wc       )) deallocate(cgrid%qmean_vapor_wc       )
    if(associated(cgrid%qmean_sensible_ac    )) deallocate(cgrid%qmean_sensible_ac    )
    if(associated(cgrid%qmean_sensible_gc    )) deallocate(cgrid%qmean_sensible_gc    )
    if(associated(cgrid%qmean_sensible_lc    )) deallocate(cgrid%qmean_sensible_lc    )
    if(associated(cgrid%qmean_sensible_wc    )) deallocate(cgrid%qmean_sensible_wc    )
    if(associated(cgrid%qmean_nee            )) deallocate(cgrid%qmean_nee            )
    if(associated(cgrid%qmean_gpp            )) deallocate(cgrid%qmean_gpp            )
    if(associated(cgrid%qmean_nep            )) deallocate(cgrid%qmean_nep            )
    if(associated(cgrid%qmean_plresp         )) deallocate(cgrid%qmean_plresp         )
    if(associated(cgrid%qmean_rh             )) deallocate(cgrid%qmean_rh             )
    if(associated(cgrid%qmean_leaf_resp      )) deallocate(cgrid%qmean_leaf_resp      )
    if(associated(cgrid%qmean_root_resp      )) deallocate(cgrid%qmean_root_resp      )
    if(associated(cgrid%qmean_soil_temp      )) deallocate(cgrid%qmean_soil_temp      )
    if(associated(cgrid%qmean_soil_mstpot    )) deallocate(cgrid%qmean_soil_mstpot    )
    if(associated(cgrid%qmean_soil_water     )) deallocate(cgrid%qmean_soil_water     )
    if(associated(cgrid%qmean_fs_open        )) deallocate(cgrid%qmean_fs_open        )
    if(associated(cgrid%qmean_fsw            )) deallocate(cgrid%qmean_fsw            )
    if(associated(cgrid%qmean_fsn            )) deallocate(cgrid%qmean_fsn            )
    if(associated(cgrid%qmean_can_temp       )) deallocate(cgrid%qmean_can_temp       )
    if(associated(cgrid%qmean_can_shv        )) deallocate(cgrid%qmean_can_shv        )
    if(associated(cgrid%qmean_can_co2        )) deallocate(cgrid%qmean_can_co2        )
    if(associated(cgrid%qmean_can_rhos       )) deallocate(cgrid%qmean_can_rhos       )
    if(associated(cgrid%qmean_can_prss       )) deallocate(cgrid%qmean_can_prss       )
    if(associated(cgrid%qmean_can_theta      )) deallocate(cgrid%qmean_can_theta      )
    if(associated(cgrid%qmean_can_theiv      )) deallocate(cgrid%qmean_can_theiv      )
    if(associated(cgrid%qmean_gnd_temp       )) deallocate(cgrid%qmean_gnd_temp       )
    if(associated(cgrid%qmean_gnd_shv        )) deallocate(cgrid%qmean_gnd_shv        )
    if(associated(cgrid%qmean_leaf_energy    )) deallocate(cgrid%qmean_leaf_energy    )
    if(associated(cgrid%qmean_leaf_water     )) deallocate(cgrid%qmean_leaf_water     )
    if(associated(cgrid%qmean_leaf_temp      )) deallocate(cgrid%qmean_leaf_temp      )
    if(associated(cgrid%qmean_leaf_hcap      )) deallocate(cgrid%qmean_leaf_hcap      )
    if(associated(cgrid%qmean_wood_energy    )) deallocate(cgrid%qmean_wood_energy    )
    if(associated(cgrid%qmean_wood_water     )) deallocate(cgrid%qmean_wood_water     )
    if(associated(cgrid%qmean_wood_temp      )) deallocate(cgrid%qmean_wood_temp      )
    if(associated(cgrid%qmean_wood_hcap      )) deallocate(cgrid%qmean_wood_hcap      )
    if(associated(cgrid%qmean_atm_temp       )) deallocate(cgrid%qmean_atm_temp       )
    if(associated(cgrid%qmean_rshort         )) deallocate(cgrid%qmean_rshort         )
    if(associated(cgrid%qmean_rshort_diff    )) deallocate(cgrid%qmean_rshort_diff    )
    if(associated(cgrid%qmean_rlong          )) deallocate(cgrid%qmean_rlong          )
    if(associated(cgrid%qmean_rshort_gnd     )) deallocate(cgrid%qmean_rshort_gnd     )
    if(associated(cgrid%qmean_rlong_gnd      )) deallocate(cgrid%qmean_rlong_gnd      )
    if(associated(cgrid%qmean_albedo         )) deallocate(cgrid%qmean_albedo         )
    if(associated(cgrid%qmean_albedo_beam    )) deallocate(cgrid%qmean_albedo_beam    )
    if(associated(cgrid%qmean_albedo_diffuse )) deallocate(cgrid%qmean_albedo_diffuse )
    if(associated(cgrid%qmean_rlong_albedo   )) deallocate(cgrid%qmean_rlong_albedo   )
    if(associated(cgrid%qmean_rlongup        )) deallocate(cgrid%qmean_rlongup        )
    if(associated(cgrid%qmean_atm_co2        )) deallocate(cgrid%qmean_atm_co2        )
    if(associated(cgrid%qmean_atm_prss       )) deallocate(cgrid%qmean_atm_prss       )
    if(associated(cgrid%qmean_atm_vels       )) deallocate(cgrid%qmean_atm_vels       )
    if(associated(cgrid%qmsqu_gpp            )) deallocate(cgrid%qmsqu_gpp            )
    if(associated(cgrid%qmsqu_leaf_resp      )) deallocate(cgrid%qmsqu_leaf_resp      )
    if(associated(cgrid%qmsqu_root_resp      )) deallocate(cgrid%qmsqu_root_resp      )
    if(associated(cgrid%qmsqu_plresp         )) deallocate(cgrid%qmsqu_plresp         )
    if(associated(cgrid%qmsqu_nee            )) deallocate(cgrid%qmsqu_nee            )
    if(associated(cgrid%qmsqu_nep            )) deallocate(cgrid%qmsqu_nep            )
    if(associated(cgrid%qmsqu_rh             )) deallocate(cgrid%qmsqu_rh             )
    if(associated(cgrid%qmsqu_sensible_ac    )) deallocate(cgrid%qmsqu_sensible_ac    )
    if(associated(cgrid%qmsqu_sensible_lc    )) deallocate(cgrid%qmsqu_sensible_lc    )
    if(associated(cgrid%qmsqu_sensible_wc    )) deallocate(cgrid%qmsqu_sensible_wc    )
    if(associated(cgrid%qmsqu_sensible_gc    )) deallocate(cgrid%qmsqu_sensible_gc    )
    if(associated(cgrid%qmsqu_evap           )) deallocate(cgrid%qmsqu_evap           )
    if(associated(cgrid%qmsqu_transp         )) deallocate(cgrid%qmsqu_transp         )
    if(associated(cgrid%qmsqu_vapor_ac       )) deallocate(cgrid%qmsqu_vapor_ac       )
    if(associated(cgrid%qmsqu_vapor_lc       )) deallocate(cgrid%qmsqu_vapor_lc       )
    if(associated(cgrid%qmsqu_vapor_wc       )) deallocate(cgrid%qmsqu_vapor_wc       )
    if(associated(cgrid%qmsqu_vapor_gc       )) deallocate(cgrid%qmsqu_vapor_gc       )

    return
  end subroutine deallocate_edtype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine deallocate_polygontype(cpoly)

    implicit none

    type(polygontype),target :: cpoly

    if(associated(cpoly%sipa_id                     )) deallocate(cpoly%sipa_id                     )
    if(associated(cpoly%sipa_n                      )) deallocate(cpoly%sipa_n                      )
    if(associated(cpoly%patch_count                 )) deallocate(cpoly%patch_count                 )
    if(associated(cpoly%site                        )) deallocate(cpoly%site                        )
    if(associated(cpoly%sitenum                     )) deallocate(cpoly%sitenum                     )

    if(associated(cpoly%lsl                         )) deallocate(cpoly%lsl                         )

    if(associated(cpoly%area                        )) deallocate(cpoly%area                        )
    if(associated(cpoly%patch_area                  )) deallocate(cpoly%patch_area                  )
    if(associated(cpoly%elevation                   )) deallocate(cpoly%elevation                   )
    if(associated(cpoly%slope                       )) deallocate(cpoly%slope                       )
    if(associated(cpoly%aspect                      )) deallocate(cpoly%aspect                      )

    if(associated(cpoly%num_landuse_years           )) deallocate(cpoly%num_landuse_years           )
    if(associated(cpoly%mindbh_primary              )) deallocate(cpoly%mindbh_primary              )
    if(associated(cpoly%probharv_primary            )) deallocate(cpoly%probharv_primary            )
    if(associated(cpoly%mindbh_secondary            )) deallocate(cpoly%mindbh_secondary            )
    if(associated(cpoly%probharv_secondary          )) deallocate(cpoly%probharv_secondary          )


    if(associated(cpoly%lai_pft                     )) deallocate(cpoly%lai_pft                     )
    if(associated(cpoly%wpa_pft                     )) deallocate(cpoly%wpa_pft                     )
    if(associated(cpoly%wai_pft                     )) deallocate(cpoly%wai_pft                     )

    if(associated(cpoly%TCI                         )) deallocate(cpoly%TCI                         )
    if(associated(cpoly%pptweight                   )) deallocate(cpoly%pptweight                         )
    if(associated(cpoly%lsl                         )) deallocate(cpoly%lsl                         )
    if(associated(cpoly%hydro_next                  )) deallocate(cpoly%hydro_next                  )
    if(associated(cpoly%hydro_prev                  )) deallocate(cpoly%hydro_prev                  )
    if(associated(cpoly%moist_W                     )) deallocate(cpoly%moist_W                     )
    if(associated(cpoly%moist_f                     )) deallocate(cpoly%moist_f                     )
    if(associated(cpoly%moist_tau                   )) deallocate(cpoly%moist_tau                   )
    if(associated(cpoly%moist_zi                    )) deallocate(cpoly%moist_zi                    )
    if(associated(cpoly%baseflow                    )) deallocate(cpoly%baseflow                    )
    if(associated(cpoly%ntext_soil                  )) deallocate(cpoly%ntext_soil                  )
    if(associated(cpoly%min_monthly_temp            )) deallocate(cpoly%min_monthly_temp            )
    if(associated(cpoly%plantation                  )) deallocate(cpoly%plantation                  )
    if(associated(cpoly%agri_stocking_pft           )) deallocate(cpoly%agri_stocking_pft           )
    if(associated(cpoly%agri_stocking_density       )) deallocate(cpoly%agri_stocking_density       )
    if(associated(cpoly%plantation_stocking_pft     )) deallocate(cpoly%plantation_stocking_pft     )
    if(associated(cpoly%plantation_stocking_density )) deallocate(cpoly%plantation_stocking_density )
    if(associated(cpoly%primary_harvest_memory      )) deallocate(cpoly%primary_harvest_memory      )
    if(associated(cpoly%secondary_harvest_memory    )) deallocate(cpoly%secondary_harvest_memory    )
    if(associated(cpoly%fire_disturbance_rate       )) deallocate(cpoly%fire_disturbance_rate       )
    if(associated(cpoly%ignition_rate               )) deallocate(cpoly%ignition_rate               )
    if(associated(cpoly%lambda_fire                 )) deallocate(cpoly%lambda_fire                 )
    if(associated(cpoly%phen_pars                   )) deallocate(cpoly%phen_pars                   )
    if(associated(cpoly%nat_disturbance_rate        )) deallocate(cpoly%nat_disturbance_rate        )
    if(associated(cpoly%nat_dist_type               )) deallocate(cpoly%nat_dist_type               )
    if(associated(cpoly%disturbance_memory          )) deallocate(cpoly%disturbance_memory          )
    if(associated(cpoly%disturbance_rates           )) deallocate(cpoly%disturbance_rates           )
    if(associated(cpoly%loss_fraction               )) deallocate(cpoly%loss_fraction               )
    if(associated(cpoly%green_leaf_factor           )) deallocate(cpoly%green_leaf_factor           )
    if(associated(cpoly%leaf_aging_factor           )) deallocate(cpoly%leaf_aging_factor           )
    if(associated(cpoly%met                         )) deallocate(cpoly%met                         )
    if(associated(cpoly%basal_area                  )) deallocate(cpoly%basal_area                  )
    if(associated(cpoly%agb                         )) deallocate(cpoly%agb                         )
    if(associated(cpoly%pldens                      )) deallocate(cpoly%pldens                      )
    if(associated(cpoly%bseeds                      )) deallocate(cpoly%bseeds                      )
    
    if(associated(cpoly%basal_area_growth           )) deallocate(cpoly%basal_area_growth           )
    if(associated(cpoly%agb_growth                  )) deallocate(cpoly%agb_growth                  )
    if(associated(cpoly%basal_area_mort             )) deallocate(cpoly%basal_area_mort             )
    if(associated(cpoly%agb_mort                    )) deallocate(cpoly%agb_mort                    )
    if(associated(cpoly%basal_area_cut              )) deallocate(cpoly%basal_area_cut              )!NOT IN REGISTRY
    if(associated(cpoly%agb_cut                     )) deallocate(cpoly%agb_cut                     )

    if(associated(cpoly%cosaoi                      )) deallocate(cpoly%cosaoi                      )
    if(associated(cpoly%avg_albedo_beam             )) deallocate(cpoly%avg_albedo_beam             )
    if(associated(cpoly%avg_albedo_diffuse          )) deallocate(cpoly%avg_albedo_diffuse          )
    if(associated(cpoly%avg_rlong_albedo            )) deallocate(cpoly%avg_rlong_albedo            )
    if(associated(cpoly%avg_albedo                  )) deallocate(cpoly%avg_albedo                  )
    if(associated(cpoly%avg_rlongup                 )) deallocate(cpoly%avg_rlongup                 )

    if(associated(cpoly%daylight                    )) deallocate(cpoly%daylight                    )
    
    if(associated(cpoly%lai                         )) deallocate(cpoly%lai                         )
    if(associated(cpoly%avg_lma                     )) deallocate(cpoly%avg_lma                     )
    if(associated(cpoly%wpa                         )) deallocate(cpoly%wpa                         )
    if(associated(cpoly%wai                         )) deallocate(cpoly%wai                         )

    ! Fast time flux diagnostics
    ! ---------------------------------------------
    if(associated(cpoly%avg_vapor_lc                )) deallocate(cpoly%avg_vapor_lc                )
    if(associated(cpoly%avg_vapor_wc                )) deallocate(cpoly%avg_vapor_wc                )
    if(associated(cpoly%avg_dew_cg                  )) deallocate(cpoly%avg_dew_cg                  )
    if(associated(cpoly%avg_vapor_gc                )) deallocate(cpoly%avg_vapor_gc                )
    if(associated(cpoly%avg_wshed_vg                )) deallocate(cpoly%avg_wshed_vg                )
    if(associated(cpoly%avg_intercepted             )) deallocate(cpoly%avg_intercepted             )
    if(associated(cpoly%avg_throughfall             )) deallocate(cpoly%avg_throughfall             )
    if(associated(cpoly%avg_vapor_ac                )) deallocate(cpoly%avg_vapor_ac                )
    if(associated(cpoly%avg_transp                  )) deallocate(cpoly%avg_transp                  )
    if(associated(cpoly%avg_evap                    )) deallocate(cpoly%avg_evap                    )
    if(associated(cpoly%avg_smoist_gg               )) deallocate(cpoly%avg_smoist_gg               )
    if(associated(cpoly%avg_transloss               )) deallocate(cpoly%avg_transloss               )
    if(associated(cpoly%avg_runoff                  )) deallocate(cpoly%avg_runoff                  )
    if(associated(cpoly%avg_drainage                )) deallocate(cpoly%avg_drainage                )
    if(associated(cpoly%avg_drainage_heat           )) deallocate(cpoly%avg_drainage_heat           )
    if(associated(cpoly%aux                         )) deallocate(cpoly%aux                         )
    if(associated(cpoly%aux_s                       )) deallocate(cpoly%aux_s                       )
    if(associated(cpoly%avg_rshort_gnd              )) deallocate(cpoly%avg_rshort_gnd              )
    if(associated(cpoly%avg_rlong_gnd               )) deallocate(cpoly%avg_rlong_gnd               )
    if(associated(cpoly%avg_carbon_ac               )) deallocate(cpoly%avg_carbon_ac               )
    if(associated(cpoly%avg_sensible_lc             )) deallocate(cpoly%avg_sensible_lc             )
    if(associated(cpoly%avg_sensible_wc             )) deallocate(cpoly%avg_sensible_wc             )
    if(associated(cpoly%avg_qwshed_vg               )) deallocate(cpoly%avg_qwshed_vg               )
    if(associated(cpoly%avg_qintercepted            )) deallocate(cpoly%avg_qintercepted            )
    if(associated(cpoly%avg_qthroughfall            )) deallocate(cpoly%avg_qthroughfall            )
    if(associated(cpoly%avg_sensible_gc             )) deallocate(cpoly%avg_sensible_gc             )
    if(associated(cpoly%avg_sensible_ac             )) deallocate(cpoly%avg_sensible_ac             )
    if(associated(cpoly%avg_sensible_gg             )) deallocate(cpoly%avg_sensible_gg             )
    if(associated(cpoly%avg_runoff_heat             )) deallocate(cpoly%avg_runoff_heat             )
    if(associated(cpoly%avg_leaf_energy             )) deallocate(cpoly%avg_leaf_energy             )
    if(associated(cpoly%avg_leaf_hcap               )) deallocate(cpoly%avg_leaf_hcap               )
    if(associated(cpoly%avg_leaf_temp               )) deallocate(cpoly%avg_leaf_temp               )
    if(associated(cpoly%avg_leaf_fliq               )) deallocate(cpoly%avg_leaf_fliq               )
    if(associated(cpoly%avg_leaf_water              )) deallocate(cpoly%avg_leaf_water              )
    if(associated(cpoly%avg_wood_energy             )) deallocate(cpoly%avg_wood_energy             )
    if(associated(cpoly%avg_wood_hcap               )) deallocate(cpoly%avg_wood_hcap               )
    if(associated(cpoly%avg_wood_temp               )) deallocate(cpoly%avg_wood_temp               )
    if(associated(cpoly%avg_wood_fliq               )) deallocate(cpoly%avg_wood_fliq               )
    if(associated(cpoly%avg_wood_water              )) deallocate(cpoly%avg_wood_water              )
    if(associated(cpoly%avg_can_temp                )) deallocate(cpoly%avg_can_temp                )
    if(associated(cpoly%avg_can_shv                 )) deallocate(cpoly%avg_can_shv                 )
    if(associated(cpoly%avg_can_co2                 )) deallocate(cpoly%avg_can_co2                 )
    if(associated(cpoly%avg_can_rhos                )) deallocate(cpoly%avg_can_rhos                )
    if(associated(cpoly%avg_can_prss                )) deallocate(cpoly%avg_can_prss                )
    if(associated(cpoly%avg_can_theta               )) deallocate(cpoly%avg_can_theta               )
    if(associated(cpoly%avg_can_theiv               )) deallocate(cpoly%avg_can_theiv               )
    if(associated(cpoly%avg_can_depth               )) deallocate(cpoly%avg_can_depth               )
    if(associated(cpoly%avg_soil_energy             )) deallocate(cpoly%avg_soil_energy             )
    if(associated(cpoly%avg_soil_mstpot             )) deallocate(cpoly%avg_soil_mstpot             )
    if(associated(cpoly%avg_soil_water              )) deallocate(cpoly%avg_soil_water              )
    if(associated(cpoly%avg_soil_temp               )) deallocate(cpoly%avg_soil_temp               )
    if(associated(cpoly%avg_soil_fracliq            )) deallocate(cpoly%avg_soil_fracliq            )
    if(associated(cpoly%avg_soil_rootfrac           )) deallocate(cpoly%avg_soil_rootfrac           )
    if(associated(cpoly%avg_soil_wetness            )) deallocate(cpoly%avg_soil_wetness            )
    if(associated(cpoly%avg_skin_temp               )) deallocate(cpoly%avg_skin_temp               )
    if(associated(cpoly%avg_available_water         )) deallocate(cpoly%avg_available_water         )
    if(associated(cpoly%runoff                      )) deallocate(cpoly%runoff                      )
    if(associated(cpoly%rad_avg                     )) deallocate(cpoly%rad_avg                     )
    ! Meteorological information
    if(associated(cpoly%avg_atm_tmp                 )) deallocate(cpoly%avg_atm_tmp                 )
    if(associated(cpoly%avg_atm_shv                 )) deallocate(cpoly%avg_atm_shv                 )
    if(associated(cpoly%avg_atm_prss                )) deallocate(cpoly%avg_atm_prss                )
    ! NACP
    if(associated(cpoly%avg_sfcw_depth            )) deallocate(cpoly%avg_sfcw_depth            )
    if(associated(cpoly%avg_sfcw_energy           )) deallocate(cpoly%avg_sfcw_energy           )
    if(associated(cpoly%avg_sfcw_mass             )) deallocate(cpoly%avg_sfcw_mass             )
    if(associated(cpoly%avg_sfcw_fracliq          )) deallocate(cpoly%avg_sfcw_fracliq          )
    if(associated(cpoly%avg_sfcw_tempk            )) deallocate(cpoly%avg_sfcw_tempk            )
    if(associated(cpoly%avg_fsc                   )) deallocate(cpoly%avg_fsc                   )
    if(associated(cpoly%avg_stsc                  )) deallocate(cpoly%avg_stsc                  )
    if(associated(cpoly%avg_ssc                   )) deallocate(cpoly%avg_ssc                   )
    if(associated(cpoly%avg_bdead                 )) deallocate(cpoly%avg_bdead                 )
    if(associated(cpoly%avg_balive                )) deallocate(cpoly%avg_balive                )
    if(associated(cpoly%avg_fsn                   )) deallocate(cpoly%avg_fsn                   )
    if(associated(cpoly%avg_msn                   )) deallocate(cpoly%avg_msn                   )

    if(associated(cpoly%avg_bleaf                 )) deallocate(cpoly%avg_bleaf                 )
    if(associated(cpoly%avg_broot                 )) deallocate(cpoly%avg_broot                 )
    if(associated(cpoly%avg_bsapwood              )) deallocate(cpoly%avg_bsapwood              )
    if(associated(cpoly%avg_bstorage              )) deallocate(cpoly%avg_bstorage              )
    if(associated(cpoly%avg_bseeds                )) deallocate(cpoly%avg_bseeds                )

    if(associated(cpoly%dmean_co2_residual        )) deallocate(cpoly%dmean_co2_residual        )
    if(associated(cpoly%dmean_energy_residual     )) deallocate(cpoly%dmean_energy_residual     )
    if(associated(cpoly%dmean_water_residual      )) deallocate(cpoly%dmean_water_residual      )
    if(associated(cpoly%mmean_co2_residual        )) deallocate(cpoly%mmean_co2_residual        )
    if(associated(cpoly%mmean_energy_residual     )) deallocate(cpoly%mmean_energy_residual     )
    if(associated(cpoly%mmean_water_residual      )) deallocate(cpoly%mmean_water_residual      )

    return
  end subroutine deallocate_polygontype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine deallocate_sitetype(csite)
    use soil_bgc, only: deallocate_sitetype_sbgc
    implicit none

    type(sitetype),target :: csite
    integer :: ipa
    
    if(associated(csite%paco_id                      )) deallocate(csite%paco_id                      )
    if(associated(csite%paco_n                       )) deallocate(csite%paco_n                       )

    if(associated(csite%dist_type                    )) deallocate(csite%dist_type                    )
    if(associated(csite%age                          )) deallocate(csite%age                          )
    if(associated(csite%area                         )) deallocate(csite%area                         )
    if(associated(csite%laiarea                      )) deallocate(csite%laiarea                      )

    if(allocated(csite%sbgc%fast_soil_C))call deallocate_sitetype_sbgc(csite%sbgc)


    if(associated(csite%pname                        )) deallocate(csite%pname                        )
    if(associated(csite%sum_dgd                      )) deallocate(csite%sum_dgd                      )
    if(associated(csite%sum_chd                      )) deallocate(csite%sum_chd                      )
    if(associated(csite%plantation                   )) deallocate(csite%plantation                   )
    if(associated(csite%cohort_count                 )) deallocate(csite%cohort_count                 )
    if(associated(csite%can_theiv                    )) deallocate(csite%can_theiv                    )
    if(associated(csite%can_temp                     )) deallocate(csite%can_temp                     )
    if(associated(csite%can_shv                      )) deallocate(csite%can_shv                      )
    if(associated(csite%can_co2                      )) deallocate(csite%can_co2                      )
    if(associated(csite%can_rhos                     )) deallocate(csite%can_rhos                     )
    if(associated(csite%can_prss                     )) deallocate(csite%can_prss                     )
    if(associated(csite%can_theta                    )) deallocate(csite%can_theta                    )
    if(associated(csite%can_depth                    )) deallocate(csite%can_depth                    )
    if(associated(csite%opencan_frac                 )) deallocate(csite%opencan_frac                 )
    if(associated(csite%ggbare                       )) deallocate(csite%ggbare                       )
    if(associated(csite%ggveg                        )) deallocate(csite%ggveg                        )
    if(associated(csite%ggnet                        )) deallocate(csite%ggnet                        )
    if(associated(csite%ggsoil                       )) deallocate(csite%ggsoil                       )
    if(associated(csite%lambda_light                 )) deallocate(csite%lambda_light                 )
    if(associated(csite%dmean_lambda_light           )) deallocate(csite%dmean_lambda_light           )
    if(associated(csite%mmean_lambda_light           )) deallocate(csite%mmean_lambda_light           )
    if(associated(csite%lai                          )) deallocate(csite%lai                          )
    if(associated(csite%wpa                          )) deallocate(csite%wpa                          )
    if(associated(csite%wai                          )) deallocate(csite%wai                          )

    if(associated(csite%sfcwater_mass                )) deallocate(csite%sfcwater_mass                )
    if(associated(csite%sfcwater_energy              )) deallocate(csite%sfcwater_energy              )
    if(associated(csite%sfcwater_depth               )) deallocate(csite%sfcwater_depth               )
    if(associated(csite%rshort_s                     )) deallocate(csite%rshort_s                     )
    if(associated(csite%rshort_s_beam                )) deallocate(csite%rshort_s_beam                )
    if(associated(csite%rshort_s_diffuse             )) deallocate(csite%rshort_s_diffuse             )
    if(associated(csite%sfcwater_tempk               )) deallocate(csite%sfcwater_tempk               )
    if(associated(csite%sfcwater_fracliq             )) deallocate(csite%sfcwater_fracliq             )
    if(associated(csite%nlev_sfcwater                )) deallocate(csite%nlev_sfcwater                )
    if(associated(csite%soil_energy                  )) deallocate(csite%soil_energy                  )
    if(associated(csite%soil_water                   )) deallocate(csite%soil_water                   )
    if(associated(csite%soil_tempk                   )) deallocate(csite%soil_tempk                   )
    if(associated(csite%soil_fracliq                 )) deallocate(csite%soil_fracliq                 )
    if(associated(csite%rootdense                    )) deallocate(csite%rootdense                    )
    if(associated(csite%ground_shv                   )) deallocate(csite%ground_shv                   )
    if(associated(csite%ground_ssh                   )) deallocate(csite%ground_ssh                   )
    if(associated(csite%ground_temp                  )) deallocate(csite%ground_temp                  )
    if(associated(csite%ground_fliq                  )) deallocate(csite%ground_fliq                  )
    if(associated(csite%rough                        )) deallocate(csite%rough                        )
    if(associated(csite%par_l_max                    )) deallocate(csite%par_l_max                    )
    if(associated(csite%par_l_beam_max               )) deallocate(csite%par_l_beam_max               )
    if(associated(csite%par_l_diffuse_max            )) deallocate(csite%par_l_diffuse_max            )
    if(associated(csite%A_o_max                      )) deallocate(csite%A_o_max                      )
    if(associated(csite%A_c_max                      )) deallocate(csite%A_c_max                      )

    if(associated(csite%old_stoma_data_max           )) deallocate(csite%old_stoma_data_max           )
    if(associated(csite%old_stoma_vector_max         )) deallocate(csite%old_stoma_vector_max         )

    if(associated(csite%avg_daily_temp               )) deallocate(csite%avg_daily_temp               )
    if(associated(csite%mean_rh                      )) deallocate(csite%mean_rh                      )
    if(associated(csite%dmean_rh                     )) deallocate(csite%dmean_rh                     )
    if(associated(csite%qmean_rh                     )) deallocate(csite%qmean_rh                     )
    if(associated(csite%mmean_rh                     )) deallocate(csite%mmean_rh                     )
    if(associated(csite%dmean_albedo                 )) deallocate(csite%dmean_albedo                 )
    if(associated(csite%dmean_albedo_beam            )) deallocate(csite%dmean_albedo_beam            )
    if(associated(csite%dmean_albedo_diffuse         )) deallocate(csite%dmean_albedo_diffuse         )
    if(associated(csite%mmean_albedo                 )) deallocate(csite%mmean_albedo                 )
    if(associated(csite%mmean_albedo_beam            )) deallocate(csite%mmean_albedo_beam            )
    if(associated(csite%mmean_albedo_diffuse         )) deallocate(csite%mmean_albedo_diffuse         )
    if(associated(csite%qmean_albedo                 )) deallocate(csite%qmean_albedo                 )
    if(associated(csite%qmean_albedo_beam            )) deallocate(csite%qmean_albedo_beam            )
    if(associated(csite%qmean_albedo_diffuse         )) deallocate(csite%qmean_albedo_diffuse         )
    if(associated(csite%mean_nep                     )) deallocate(csite%mean_nep                     )
    if(associated(csite%wbudget_loss2atm             )) deallocate(csite%wbudget_loss2atm             )
    if(associated(csite%wbudget_denseffect           )) deallocate(csite%wbudget_denseffect           )
    if(associated(csite%wbudget_precipgain           )) deallocate(csite%wbudget_precipgain           )
    if(associated(csite%wbudget_loss2runoff          )) deallocate(csite%wbudget_loss2runoff          )
    if(associated(csite%wbudget_loss2drainage        )) deallocate(csite%wbudget_loss2drainage        )
    if(associated(csite%wbudget_initialstorage       )) deallocate(csite%wbudget_initialstorage       )
    if(associated(csite%wbudget_residual             )) deallocate(csite%wbudget_residual             )
    if(associated(csite%ebudget_loss2atm             )) deallocate(csite%ebudget_loss2atm             )
    if(associated(csite%ebudget_denseffect           )) deallocate(csite%ebudget_denseffect           )
    if(associated(csite%ebudget_loss2runoff          )) deallocate(csite%ebudget_loss2runoff          )
    if(associated(csite%ebudget_loss2drainage        )) deallocate(csite%ebudget_loss2drainage        )
    if(associated(csite%ebudget_netrad               )) deallocate(csite%ebudget_netrad               )
    if(associated(csite%ebudget_precipgain           )) deallocate(csite%ebudget_precipgain           )
    if(associated(csite%ebudget_initialstorage       )) deallocate(csite%ebudget_initialstorage       )
    if(associated(csite%ebudget_residual             )) deallocate(csite%ebudget_residual             )
    if(associated(csite%co2budget_initialstorage     )) deallocate(csite%co2budget_initialstorage     )
    if(associated(csite%co2budget_residual           )) deallocate(csite%co2budget_residual           )
    if(associated(csite%co2budget_loss2atm           )) deallocate(csite%co2budget_loss2atm           )
    if(associated(csite%co2budget_denseffect         )) deallocate(csite%co2budget_denseffect         )
    if(associated(csite%co2budget_gpp                )) deallocate(csite%co2budget_gpp                )
    if(associated(csite%co2budget_gpp_dbh            )) deallocate(csite%co2budget_gpp_dbh            )
    if(associated(csite%co2budget_plresp             )) deallocate(csite%co2budget_plresp             )
    if(associated(csite%co2budget_rh                 )) deallocate(csite%co2budget_rh                 )
    if(associated(csite%today_A_decomp               )) deallocate(csite%today_A_decomp               )
    if(associated(csite%today_Af_decomp              )) deallocate(csite%today_Af_decomp              )
    if(associated(csite%dmean_A_decomp               )) deallocate(csite%dmean_A_decomp               )
    if(associated(csite%dmean_Af_decomp              )) deallocate(csite%dmean_Af_decomp              )
    if(associated(csite%mmean_A_decomp               )) deallocate(csite%mmean_A_decomp               )
    if(associated(csite%mmean_Af_decomp              )) deallocate(csite%mmean_Af_decomp              )
    if(associated(csite%repro                        )) deallocate(csite%repro                        )
    if(associated(csite%veg_rough                    )) deallocate(csite%veg_rough                    )
    if(associated(csite%veg_height                   )) deallocate(csite%veg_height                   )
    if(associated(csite%veg_displace                 )) deallocate(csite%veg_displace                 )
    if(associated(csite%fsc_in                       )) deallocate(csite%fsc_in                       )
    if(associated(csite%ssc_in                       )) deallocate(csite%ssc_in                       )
    if(associated(csite%ssl_in                       )) deallocate(csite%ssl_in                       )
    if(associated(csite%fsn_in                       )) deallocate(csite%fsn_in                       )



    if(associated(csite%fsp_in)) deallocate(csite%fsp_in)
    if(associated(csite%stsp_in)) deallocate(csite%stsp_in)
    if(associated(csite%total_plant_P_uptake  )) deallocate(csite%total_plant_P_uptake  )



    if(associated(csite%total_plant_nitrogen_uptake  )) deallocate(csite%total_plant_nitrogen_uptake  )
    if(associated(csite%mineralized_N_input )) deallocate(csite%mineralized_N_input  )
    if(associated(csite%mineralized_N_loss  )) deallocate(csite%mineralized_N_loss  )
    if(associated(csite%rshort_g                     )) deallocate(csite%rshort_g                     )
    if(associated(csite%rshort_g_beam                )) deallocate(csite%rshort_g_beam                )
    if(associated(csite%rshort_g_diffuse             )) deallocate(csite%rshort_g_diffuse             )
    if(associated(csite%rlong_g                      )) deallocate(csite%rlong_g                      )
    if(associated(csite%rlong_g_surf                 )) deallocate(csite%rlong_g_surf                 )
    if(associated(csite%rlong_g_incid                )) deallocate(csite%rlong_g_incid                )
    if(associated(csite%rlong_s                      )) deallocate(csite%rlong_s                      )
    if(associated(csite%rlong_s_surf                 )) deallocate(csite%rlong_s_surf                 )
    if(associated(csite%rlong_s_incid                )) deallocate(csite%rlong_s_incid                )
    if(associated(csite%albedo                       )) deallocate(csite%albedo                       )
    if(associated(csite%albedo_beam                  )) deallocate(csite%albedo_beam                  )
    if(associated(csite%albedo_diffuse               )) deallocate(csite%albedo_diffuse               )
    if(associated(csite%rlongup                      )) deallocate(csite%rlongup                      )
    if(associated(csite%rlong_albedo                 )) deallocate(csite%rlong_albedo                 )
    if(associated(csite%total_sfcw_depth             )) deallocate(csite%total_sfcw_depth             )
    if(associated(csite%snowfac                      )) deallocate(csite%snowfac                      )
    if(associated(csite%A_decomp                     )) deallocate(csite%A_decomp                     )
    if(associated(csite%f_decomp                     )) deallocate(csite%f_decomp                     )
    if(associated(csite%rh                           )) deallocate(csite%rh                           )
    if(associated(csite%cwd_rh                       )) deallocate(csite%cwd_rh                       )
    if(associated(csite%fuse_flag                    )) deallocate(csite%fuse_flag                    )
    if(associated(csite%cumlai_profile               )) deallocate(csite%cumlai_profile               )
    if(associated(csite%plant_ag_biomass             )) deallocate(csite%plant_ag_biomass             )

    if(associated(csite%mean_wflux                   )) deallocate(csite%mean_wflux                   )
    if(associated(csite%mean_latflux                 )) deallocate(csite%mean_latflux                 )
    if(associated(csite%mean_hflux                   )) deallocate(csite%mean_hflux                   )
    if(associated(csite%mean_runoff                  )) deallocate(csite%mean_runoff                  )
    if(associated(csite%mean_qrunoff                 )) deallocate(csite%mean_qrunoff                 )

    if(associated(csite%htry                         )) deallocate(csite%htry                         )
    if(associated(csite%avg_rk4step                  )) deallocate(csite%avg_rk4step                  )
    if(associated(csite%dmean_rk4step                )) deallocate(csite%dmean_rk4step                )
    if(associated(csite%mmean_rk4step                )) deallocate(csite%mmean_rk4step                )

    if(associated(csite%avg_available_water          )) deallocate(csite%avg_available_water          )

    if(associated(csite%ustar                        )) deallocate(csite%ustar                        )
    if(associated(csite%tstar                        )) deallocate(csite%tstar                        )
    if(associated(csite%qstar                        )) deallocate(csite%qstar                        )
    if(associated(csite%cstar                        )) deallocate(csite%cstar                        )

    if(associated(csite%zeta                         )) deallocate(csite%zeta                         )
    if(associated(csite%ribulk                       )) deallocate(csite%ribulk                       )

    if(associated(csite%upwp                         )) deallocate(csite%upwp                         )
    if(associated(csite%qpwp                         )) deallocate(csite%qpwp                         )
    if(associated(csite%cpwp                         )) deallocate(csite%cpwp                         )
    if(associated(csite%tpwp                         )) deallocate(csite%tpwp                         )
    if(associated(csite%wpwp                         )) deallocate(csite%wpwp                         )

    if(associated(csite%avg_rshort_gnd               )) deallocate(csite%avg_rshort_gnd               )
    if(associated(csite%avg_rlong_gnd                )) deallocate(csite%avg_rlong_gnd                )
    if(associated(csite%avg_carbon_ac                )) deallocate(csite%avg_carbon_ac                )
    if(associated(csite%avg_rlongup                  )) deallocate(csite%avg_rlongup                  )
    if(associated(csite%avg_albedo                   )) deallocate(csite%avg_albedo                   )
    if(associated(csite%avg_albedo_beam              )) deallocate(csite%avg_albedo_beam              )
    if(associated(csite%avg_albedo_diffuse           )) deallocate(csite%avg_albedo_diffuse           )
    if(associated(csite%avg_rlong_albedo             )) deallocate(csite%avg_rlong_albedo             )

    if(associated(csite%avg_vapor_lc                 )) deallocate(csite%avg_vapor_lc                 )
    if(associated(csite%avg_vapor_wc                 )) deallocate(csite%avg_vapor_wc                 )
    if(associated(csite%avg_dew_cg                   )) deallocate(csite%avg_dew_cg                   )
    if(associated(csite%avg_vapor_gc                 )) deallocate(csite%avg_vapor_gc                 )
    if(associated(csite%avg_wshed_vg                 )) deallocate(csite%avg_wshed_vg                 )
    if(associated(csite%avg_intercepted              )) deallocate(csite%avg_intercepted              )
    if(associated(csite%avg_throughfall              )) deallocate(csite%avg_throughfall              )
    if(associated(csite%avg_vapor_ac                 )) deallocate(csite%avg_vapor_ac                 )
    if(associated(csite%avg_transp                   )) deallocate(csite%avg_transp                   )
    if(associated(csite%avg_evap                     )) deallocate(csite%avg_evap                     )
    if(associated(csite%avg_smoist_gg                )) deallocate(csite%avg_smoist_gg                )
    if(associated(csite%avg_transloss                )) deallocate(csite%avg_transloss                )
    if(associated(csite%avg_runoff                   )) deallocate(csite%avg_runoff                   )
    if(associated(csite%avg_drainage                 )) deallocate(csite%avg_drainage                 )
    if(associated(csite%avg_drainage_heat            )) deallocate(csite%avg_drainage_heat            )
    if(associated(csite%aux                          )) deallocate(csite%aux                          )
    if(associated(csite%aux_s                        )) deallocate(csite%aux_s                        )
    if(associated(csite%avg_sensible_lc              )) deallocate(csite%avg_sensible_lc              )
    if(associated(csite%avg_sensible_wc              )) deallocate(csite%avg_sensible_wc              )
    if(associated(csite%avg_qwshed_vg                )) deallocate(csite%avg_qwshed_vg                )
    if(associated(csite%avg_qintercepted             )) deallocate(csite%avg_qintercepted             )
    if(associated(csite%avg_qthroughfall             )) deallocate(csite%avg_qthroughfall             )
    if(associated(csite%avg_sensible_gc              )) deallocate(csite%avg_sensible_gc              )
    if(associated(csite%avg_sensible_ac              )) deallocate(csite%avg_sensible_ac              )
    if(associated(csite%avg_sensible_gg              )) deallocate(csite%avg_sensible_gg              )
    if(associated(csite%avg_runoff_heat              )) deallocate(csite%avg_runoff_heat              )
    if(associated(csite%avg_leaf_energy              )) deallocate(csite%avg_leaf_energy              )
    if(associated(csite%avg_leaf_temp                )) deallocate(csite%avg_leaf_temp                )
    if(associated(csite%avg_leaf_hcap                )) deallocate(csite%avg_leaf_hcap                )
    if(associated(csite%avg_leaf_fliq                )) deallocate(csite%avg_leaf_fliq                )
    if(associated(csite%avg_leaf_water               )) deallocate(csite%avg_leaf_water               )
    if(associated(csite%avg_wood_energy              )) deallocate(csite%avg_wood_energy              )
    if(associated(csite%avg_wood_temp                )) deallocate(csite%avg_wood_temp                )
    if(associated(csite%avg_wood_hcap                )) deallocate(csite%avg_wood_hcap                )
    if(associated(csite%avg_wood_fliq                )) deallocate(csite%avg_wood_fliq                )
    if(associated(csite%avg_wood_water               )) deallocate(csite%avg_wood_water               )

    if(associated(csite%watertable                   )) deallocate(csite%watertable                   )
    if(associated(csite%moist_dz                     )) deallocate(csite%moist_dz                     )
    if(associated(csite%ksat                         )) deallocate(csite%ksat                         )
    if(associated(csite%soil_sat_energy              )) deallocate(csite%soil_sat_energy              )
    if(associated(csite%soil_sat_water               )) deallocate(csite%soil_sat_water               )
    if(associated(csite%soil_sat_heat                )) deallocate(csite%soil_sat_heat                )
    if(associated(csite%runoff_A                     )) deallocate(csite%runoff_A                     )
    if(associated(csite%runoff_rate                  )) deallocate(csite%runoff_rate                  )
    if(associated(csite%runoff                       )) deallocate(csite%runoff                       )

    if(associated(csite%current_paw))deallocate(csite%current_paw)
    if(associated(csite%past_paw))deallocate(csite%past_paw)

    if(associated(csite%dmean_co2_residual        )) deallocate(csite%dmean_co2_residual        )
    if(associated(csite%dmean_energy_residual     )) deallocate(csite%dmean_energy_residual     )
    if(associated(csite%dmean_water_residual      )) deallocate(csite%dmean_water_residual      )
    if(associated(csite%mmean_co2_residual        )) deallocate(csite%mmean_co2_residual        )
    if(associated(csite%mmean_energy_residual     )) deallocate(csite%mmean_energy_residual     )
    if(associated(csite%mmean_water_residual      )) deallocate(csite%mmean_water_residual      )

    do ipa=1,csite%npatches
       if (csite%patch(ipa)%ncohorts > 0) call deallocate_patchtype(csite%patch(ipa))
    end do

    if(associated(csite%patch                        )) deallocate(csite%patch                        )

    return
  end subroutine deallocate_sitetype
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine deallocate_patchtype(cpatch)

    use cohort_state, only: deallocate_patchtype_state
    use cohort_phen, only: deallocate_patchtype_phen
    use cohort_mort, only: deallocate_patchtype_mort
    use cohort_resp, only: deallocate_patchtype_resp
    use cohort_photo, only: deallocate_patchtype_photo
    use cohort_rad, only: deallocate_patchtype_rad
    use cohort_therm, only: deallocate_patchtype_therm

    implicit none

    type(patchtype),target :: cpatch
    
    if (cpatch%ncohorts == 0) return
    
    call deallocate_patchtype_state(cpatch%costate)
    call deallocate_patchtype_phen(cpatch%cophen)
    call deallocate_patchtype_mort(cpatch%comort)
    call deallocate_patchtype_resp(cpatch%coresp)
    call deallocate_patchtype_photo(cpatch%cophoto)
    call deallocate_patchtype_rad(cpatch%corad)
    call deallocate_patchtype_therm(cpatch%cotherm)

    return
  end subroutine deallocate_patchtype
!============================================================================!
!============================================================================!

   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine copies the patches from the input site to the output.   The       !
   ! number of patches to be copied should be always match between the input and output,   !
   ! and every single variable, including the cohorts, will be copied.                     !
   ! IMPORTANT.  This subroutine assumes that the output patches still don't have cohorts  !
   !             allocated, so this should be never used in a previously allocated patch.  !
   !---------------------------------------------------------------------------------------!
   subroutine copy_sitetype(isite,osite,ipaa,ipaz,opaa,opaz)

     use soil_bgc, only: copy_sitetype_sbgc

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)  , target     :: isite ! Input  (donor) site
      type(sitetype)  , target     :: osite ! Output (receptor) site
      integer         , intent(in) :: ipaa  ! First input  patch index
      integer         , intent(in) :: ipaz  ! Last  input  patch index
      integer         , intent(in) :: opaa  ! First output patch index
      integer         , intent(in) :: opaz  ! Last  output patch index
      !----- Local variables. -------------------------------------------------------------!
      type(stoma_data), pointer    :: osdi  ! Old stomate data - input patch
      type(stoma_data), pointer    :: osdo  ! Old stomate data - output patch
      integer                      :: ipa   ! Counter for the input  site patches
      integer                      :: opa   ! Counter for the output site patches
      integer                      :: k     ! Vertical layer counter
      integer                      :: ipft  ! PFT counter
      integer                      :: isto  ! Stomate attribute counter
      integer                      :: idbh  ! DBH counter
      integer                      :: ihgt  ! Height counter
      integer                      :: icyc  ! Time of day counter
      !------------------------------------------------------------------------------------!


      !------ Check whether this patch copying makes sense. -------------------------------!
      if (ipaz - ipaa /= opaz - opaa) then
         write (unit=*,fmt='(a)')       '------------------------------------------------'
         write (unit=*,fmt='(a)')       ' - Input site:'
         write (unit=*,fmt='(a,1x,i6)') '   * First patch: ',ipaa
         write (unit=*,fmt='(a,1x,i6)') '   * Last  patch: ',ipaz
         write (unit=*,fmt='(a,1x,i6)') '   * Patch count: ',ipaz-ipaa+1
         write (unit=*,fmt='(a)')       ' - Output site:'
         write (unit=*,fmt='(a,1x,i6)') '   * First patch: ',opaa
         write (unit=*,fmt='(a,1x,i6)') '   * Last  patch: ',opaz
         write (unit=*,fmt='(a,1x,i6)') '   * Patch count: ',opaz-opaa+1
         write (unit=*,fmt='(a)'      ) '------------------------------------------------'
         call fatal_error ('Patch count of input and output paths don''t match'            &
                          ,'copy_sitetype','ed_state_vars.f90')
      end if
      !------------------------------------------------------------------------------------!

      ipa = ipaa - 1
      opaloop: do opa = opaa, opaz
         ipa = ipa + 1

         osite%paco_id(opa)                     = isite%paco_id(ipa)
         osite%paco_n(opa)                      = isite%paco_n(ipa)
         osite%dist_type(opa)                   = isite%dist_type(ipa)
         osite%age(opa)                         = isite%age(ipa)
         osite%area(opa)                        = isite%area(ipa)

         call copy_sitetype_sbgc(osite%sbgc, isite%sbgc, opa, ipa)

         osite%sum_dgd(opa)                     = isite%sum_dgd(ipa)
         osite%sum_chd(opa)                     = isite%sum_chd(ipa)
         osite%plantation(opa)                  = isite%plantation(ipa)
         osite%cohort_count(opa)                = isite%cohort_count(ipa)
         osite%can_theiv(opa)                   = isite%can_theiv(ipa)
         osite%can_temp(opa)                    = isite%can_temp(ipa)
         osite%can_shv(opa)                     = isite%can_shv(ipa)
         osite%can_co2(opa)                     = isite%can_co2(ipa)
         osite%can_rhos(opa)                    = isite%can_rhos(ipa)
         osite%can_prss(opa)                    = isite%can_prss(ipa)
         osite%can_theta(opa)                   = isite%can_theta(ipa)
         osite%can_depth(opa)                   = isite%can_depth(ipa)
         osite%opencan_frac(opa)                = isite%opencan_frac(ipa)
         osite%ggbare(opa)                      = isite%ggbare(ipa)
         osite%ggveg(opa)                       = isite%ggveg(ipa)
         osite%ggnet(opa)                       = isite%ggnet(ipa)
         osite%ggsoil(opa)                      = isite%ggsoil(ipa)
         osite%lambda_light(opa)                = isite%lambda_light(ipa)
         osite%lai(opa)                         = isite%lai(ipa)
         osite%wpa(opa)                         = isite%wpa(ipa)
         osite%wai(opa)                         = isite%wai(ipa)
         osite%avg_daily_temp(opa)              = isite%avg_daily_temp(ipa)
         osite%mean_rh(opa)                     = isite%mean_rh(ipa)
         osite%mean_nep(opa)                    = isite%mean_nep(ipa)
         osite%wbudget_loss2atm(opa)            = isite%wbudget_loss2atm(ipa)
         osite%wbudget_denseffect(opa)          = isite%wbudget_denseffect(ipa)
         osite%wbudget_precipgain(opa)          = isite%wbudget_precipgain(ipa)
         osite%wbudget_loss2runoff(opa)         = isite%wbudget_loss2runoff(ipa)
         osite%wbudget_loss2drainage(opa)       = isite%wbudget_loss2drainage(ipa)
         osite%wbudget_initialstorage(opa)      = isite%wbudget_initialstorage(ipa)
         osite%wbudget_residual(opa)            = isite%wbudget_residual(ipa)
         osite%ebudget_loss2atm(opa)            = isite%ebudget_loss2atm(ipa)
         osite%ebudget_denseffect(opa)          = isite%ebudget_denseffect(ipa)
         osite%ebudget_loss2runoff(opa)         = isite%ebudget_loss2runoff(ipa)
         osite%ebudget_loss2drainage(opa)       = isite%ebudget_loss2drainage(ipa)
         osite%ebudget_netrad(opa)              = isite%ebudget_netrad(ipa)
         osite%ebudget_precipgain(opa)          = isite%ebudget_precipgain(ipa)
         osite%ebudget_initialstorage(opa)      = isite%ebudget_initialstorage(ipa)
         osite%ebudget_residual(opa)            = isite%ebudget_residual(ipa)
         osite%co2budget_initialstorage(opa)    = isite%co2budget_initialstorage(ipa)
         osite%co2budget_residual(opa)          = isite%co2budget_residual(ipa)
         osite%co2budget_loss2atm(opa)          = isite%co2budget_loss2atm(ipa)
         osite%co2budget_denseffect(opa)        = isite%co2budget_denseffect(ipa)
         osite%co2budget_gpp(opa)               = isite%co2budget_gpp(ipa)
         osite%co2budget_plresp(opa)            = isite%co2budget_plresp(ipa)
         osite%co2budget_rh(opa)                = isite%co2budget_rh(ipa)
         osite%today_A_decomp(opa)              = isite%today_A_decomp(ipa)
         osite%today_Af_decomp(opa)             = isite%today_Af_decomp(ipa)
         osite%veg_rough(opa)                   = isite%veg_rough(ipa)
         osite%veg_height(opa)                  = isite%veg_height(ipa)
         osite%veg_displace(opa)                = isite%veg_displace(ipa)
         osite%fsc_in(opa)                      = isite%fsc_in(ipa)
         osite%ssc_in(opa)                      = isite%ssc_in(ipa)
         osite%ssl_in(opa)                      = isite%ssl_in(ipa)
         osite%fsn_in(opa)                      = isite%fsn_in(ipa)



         osite%fsp_in(opa)                      = isite%fsp_in(ipa)
         osite%stsp_in(opa)                      = isite%stsp_in(ipa)
         osite%total_plant_P_uptake(opa) = isite%total_plant_P_uptake(ipa)



         osite%total_plant_nitrogen_uptake(opa) = isite%total_plant_nitrogen_uptake(ipa)
         osite%mineralized_N_loss(opa)          = isite%mineralized_N_loss(ipa)
         osite%mineralized_N_input(opa)         = isite%mineralized_N_input(ipa)
         osite%rshort_g(opa)                    = isite%rshort_g(ipa)
         osite%rshort_g_beam(opa)               = isite%rshort_g_beam(ipa)
         osite%rshort_g_diffuse(opa)            = isite%rshort_g_diffuse(ipa)
         osite%rlong_g(opa)                     = isite%rlong_g(ipa)
         osite%rlong_g_surf(opa)                = isite%rlong_g_surf(ipa)
         osite%rlong_g_incid(opa)               = isite%rlong_g_incid(ipa)
         osite%rlong_s(opa)                     = isite%rlong_s(ipa)
         osite%rlong_s_surf(opa)                = isite%rlong_s_surf(ipa)
         osite%rlong_s_incid(opa)               = isite%rlong_s_incid(ipa)
         osite%albedo(opa)                      = isite%albedo(ipa)
         osite%albedo_beam(opa)                 = isite%albedo_beam(ipa)
         osite%albedo_diffuse(opa)              = isite%albedo_diffuse(ipa)
         osite%rlongup(opa)                     = isite%rlongup(ipa)
         osite%rlong_albedo(opa)                = isite%rlong_albedo(ipa)
         osite%total_sfcw_depth(opa)            = isite%total_sfcw_depth(ipa)
         osite%snowfac(opa)                     = isite%snowfac(ipa)
         osite%A_decomp(opa)                    = isite%A_decomp(ipa)
         osite%f_decomp(opa)                    = isite%f_decomp(ipa)
         osite%rh(opa)                          = isite%rh(ipa)
         osite%cwd_rh(opa)                      = isite%cwd_rh(ipa)
         osite%fuse_flag(opa)                   = isite%fuse_flag(ipa)
         osite%plant_ag_biomass(opa)            = isite%plant_ag_biomass(ipa)
         osite%mean_wflux(opa)                  = isite%mean_wflux(ipa)
         osite%mean_latflux(opa)                = isite%mean_latflux(ipa)
         osite%mean_hflux(opa)                  = isite%mean_hflux(ipa)
         osite%runoff(opa)                      = isite%runoff(ipa)
         osite%mean_runoff(opa)                 = isite%mean_runoff(ipa)
         osite%mean_qrunoff(opa)                = isite%mean_qrunoff(ipa)
         osite%htry(opa)                        = isite%htry(ipa)
         osite%avg_rk4step(opa)                 = isite%avg_rk4step(ipa)
         osite%avg_available_water(opa)         = isite%avg_available_water(ipa)
         osite%ustar(opa)                       = isite%ustar(ipa)
         osite%tstar(opa)                       = isite%tstar(ipa)
         osite%qstar(opa)                       = isite%qstar(ipa)
         osite%cstar(opa)                       = isite%cstar(ipa)

         osite%zeta(opa)                        = isite%zeta(ipa)
         osite%ribulk(opa)                      = isite%ribulk(ipa)

         osite%upwp(opa)                        = isite%upwp(ipa)
         osite%tpwp(opa)                        = isite%tpwp(ipa)
         osite%qpwp(opa)                        = isite%qpwp(ipa)
         osite%cpwp(opa)                        = isite%cpwp(ipa)
         osite%wpwp(opa)                        = isite%wpwp(ipa)

         osite%nlev_sfcwater(opa)               = isite%nlev_sfcwater(ipa)
         osite%ground_shv(opa)                  = isite%ground_shv(ipa)
         osite%ground_ssh(opa)                  = isite%ground_ssh(ipa)
         osite%ground_temp(opa)                 = isite%ground_temp(ipa)
         osite%ground_fliq(opa)                 = isite%ground_fliq(ipa)
         osite%rough(opa)                       = isite%rough(ipa)

         osite%avg_rshort_gnd     (opa)         = isite%avg_rshort_gnd     (ipa)
         osite%avg_rlong_gnd      (opa)         = isite%avg_rlong_gnd      (ipa)
         osite%avg_carbon_ac      (opa)         = isite%avg_carbon_ac      (ipa)
         osite%avg_rlongup        (opa)         = isite%avg_rlongup        (ipa)
         osite%avg_albedo         (opa)         = isite%avg_albedo         (ipa)
         osite%avg_albedo_beam    (opa)         = isite%avg_albedo_beam    (ipa)
         osite%avg_albedo_diffuse (opa)         = isite%avg_albedo_diffuse (ipa)
         osite%avg_rlong_albedo   (opa)         = isite%avg_rlong_albedo   (ipa)

         osite%avg_vapor_lc(opa)                = isite%avg_vapor_lc(ipa)
         osite%avg_vapor_wc(opa)                = isite%avg_vapor_wc(ipa)
         osite%avg_dew_cg(opa)                  = isite%avg_dew_cg(ipa)
         osite%avg_vapor_gc(opa)                = isite%avg_vapor_gc(ipa)
         osite%avg_wshed_vg(opa)                = isite%avg_wshed_vg(ipa)
         osite%avg_intercepted(opa)             = isite%avg_intercepted(ipa)
         osite%avg_throughfall(opa)             = isite%avg_throughfall(ipa)
         osite%avg_vapor_ac(opa)                = isite%avg_vapor_ac(ipa)
         osite%avg_transp(opa)                  = isite%avg_transp(ipa)
         osite%avg_evap(opa)                    = isite%avg_evap(ipa)
         osite%avg_runoff(opa)                  = isite%avg_runoff(ipa)
         osite%avg_drainage(opa)                = isite%avg_drainage(ipa)
         osite%avg_drainage_heat(opa)           = isite%avg_drainage_heat(ipa)
         osite%aux(opa)                         = isite%aux(ipa)
         osite%avg_sensible_lc(opa)             = isite%avg_sensible_lc(ipa)
         osite%avg_sensible_wc(opa)             = isite%avg_sensible_wc(ipa)
         osite%avg_qwshed_vg(opa)               = isite%avg_qwshed_vg(ipa)
         osite%avg_qintercepted(opa)            = isite%avg_qintercepted(ipa)
         osite%avg_qthroughfall(opa)            = isite%avg_qthroughfall(ipa)
         osite%avg_sensible_gc(opa)             = isite%avg_sensible_gc(ipa)
         osite%avg_sensible_ac(opa)             = isite%avg_sensible_ac(ipa)
         osite%avg_runoff_heat(opa)             = isite%avg_runoff_heat(ipa)
         osite%avg_leaf_energy(opa)             = isite%avg_leaf_energy(ipa)
         osite%avg_leaf_temp(opa)               = isite%avg_leaf_temp(ipa)
         osite%avg_leaf_hcap(opa)               = isite%avg_leaf_hcap(ipa)
         osite%avg_leaf_fliq(opa)               = isite%avg_leaf_fliq(ipa)
         osite%avg_leaf_water(opa)              = isite%avg_leaf_water(ipa)
         osite%avg_wood_energy(opa)             = isite%avg_wood_energy(ipa)
         osite%avg_wood_temp(opa)               = isite%avg_wood_temp(ipa)
         osite%avg_wood_hcap(opa)               = isite%avg_wood_hcap(ipa)
         osite%avg_wood_fliq(opa)               = isite%avg_wood_fliq(ipa)
         osite%avg_wood_water(opa)              = isite%avg_wood_water(ipa)
         osite%par_l_max(opa)                   = isite%par_l_max(ipa)
         osite%par_l_beam_max(opa)              = isite%par_l_beam_max(ipa)
         osite%par_l_diffuse_max(opa)           = isite%par_l_diffuse_max(ipa)

         !----- Copy the temporary surface water variables. -------------------------------!
         do k=1,nzs
            osite%sfcwater_mass(k,opa)          = isite%sfcwater_mass(k,ipa)
            osite%sfcwater_energy(k,opa)        = isite%sfcwater_energy(k,ipa)
            osite%sfcwater_depth(k,opa)         = isite%sfcwater_depth(k,ipa)
            osite%rshort_s(k,opa)               = isite%rshort_s(k,ipa)
            osite%rshort_s_beam(k,opa)          = isite%rshort_s_beam(k,ipa)
            osite%rshort_s_diffuse(k,opa)       = isite%rshort_s_diffuse(k,ipa)
            osite%sfcwater_tempk(k,opa)         = isite%sfcwater_tempk(k,ipa)
            osite%sfcwater_fracliq(k,opa)       = isite%sfcwater_fracliq(k,ipa)
         end do

         !----- Copy the soil variables. --------------------------------------------------!
         do k=1,nzg
            osite%soil_energy(k,opa)            =  isite%soil_energy(k,ipa)
            osite%soil_water(k,opa)             =  isite%soil_water(k,ipa)
            osite%soil_tempk(k,opa)             =  isite%soil_tempk(k,ipa)
            osite%soil_fracliq(k,opa)           =  isite%soil_fracliq(k,ipa)
            osite%rootdense(k,opa)              =  isite%rootdense(k,ipa)
            osite%avg_smoist_gg(k,opa)          =  isite%avg_smoist_gg(k,ipa)
            osite%avg_transloss(k,opa)          =  isite%avg_transloss(k,ipa)
            osite%avg_sensible_gg(k,opa)        =  isite%avg_sensible_gg(k,ipa)
            osite%aux_s(k,opa)                  =  isite%aux_s(k,ipa)

            osite%current_paw(k,opa) = isite%current_paw(k,ipa)
            do ihgt = 1, 10
               osite%past_paw(k,ihgt,opa) = isite%past_paw(k,ihgt,ipa)
            enddo
         end do

         !----- PFT types. ----------------------------------------------------------------!
         do ipft=1,n_pft
            osite%repro(ipft,opa)                =  isite%repro(ipft,ipa)
            osite%A_o_max(ipft,opa)              =  isite%A_o_max(ipft,ipa)
            osite%A_c_max(ipft,opa)              =  isite%A_c_max(ipft,ipa)
            
            do isto=1,n_stoma_atts
               osite%old_stoma_vector_max(isto,ipft,opa) =                                 &
                                                  isite%old_stoma_vector_max(isto,ipft,ipa)
            end do

            do ihgt=1,ff_nhgt
               osite%cumlai_profile(ipft,ihgt,opa) = isite%cumlai_profile(ipft,ihgt,ipa)
            end do
            
            !----- This is to copy the old_stoma_data_max structure. ----------------------!
            osdo => osite%old_stoma_data_max(ipft,opa)
            osdi => isite%old_stoma_data_max(ipft,ipa)
           
            osdo%recalc           = osdi%recalc
            osdo%T_L              = osdi%T_L
            osdo%e_A              = osdi%e_A
            osdo%PAR              = osdi%PAR
            osdo%rb_factor        = osdi%rb_factor
            osdo%prss             = osdi%prss
            osdo%phenology_factor = osdi%phenology_factor
            osdo%gsw_open         = osdi%gsw_open
            osdo%ilimit           = osdi%ilimit
            osdo%T_L_residual     = osdi%T_L_residual
            osdo%e_a_residual     = osdi%e_a_residual
            osdo%par_residual     = osdi%par_residual
            osdo%rb_residual      = osdi%rb_residual
            osdo%leaf_residual    = osdi%leaf_residual
            osdo%gsw_residual     = osdi%gsw_residual
         end do

         !----- DBH types. ----------------------------------------------------------------!
         do idbh=1,n_dbh
            osite%co2budget_gpp_dbh(idbh,opa) = isite%co2budget_gpp_dbh(idbh,ipa)
         end do

         !----- Daily averages. -----------------------------------------------------------!
         if (idoutput > 0 .or. imoutput > 0 .or. iqoutput > 0) then
            osite%dmean_rh             (opa) = isite%dmean_rh             (ipa)
            osite%dmean_co2_residual   (opa) = isite%dmean_co2_residual   (ipa)
            osite%dmean_energy_residual(opa) = isite%dmean_energy_residual(ipa)
            osite%dmean_water_residual (opa) = isite%dmean_water_residual (ipa)
            osite%dmean_lambda_light   (opa) = isite%dmean_lambda_light   (ipa)
            osite%dmean_A_decomp       (opa) = isite%dmean_A_decomp       (ipa)
            osite%dmean_Af_decomp      (opa) = isite%dmean_Af_decomp      (ipa)
            osite%dmean_rk4step        (opa) = isite%dmean_rk4step        (ipa)
            osite%dmean_albedo         (opa) = isite%dmean_albedo         (ipa)
            osite%dmean_albedo_beam    (opa) = isite%dmean_albedo_beam    (ipa)
            osite%dmean_albedo_diffuse (opa) = isite%dmean_albedo_diffuse (ipa)
         end if

         if (imoutput > 0 .or. iqoutput > 0) then
            osite%mmean_rh             (opa) = isite%mmean_rh             (ipa)
            osite%mmean_co2_residual   (opa) = isite%mmean_co2_residual   (ipa)
            osite%mmean_energy_residual(opa) = isite%mmean_energy_residual(ipa)
            osite%mmean_water_residual (opa) = isite%mmean_water_residual (ipa)
            osite%mmean_lambda_light   (opa) = isite%mmean_lambda_light   (ipa)
            osite%mmean_A_decomp       (opa) = isite%mmean_A_decomp       (ipa)
            osite%mmean_Af_decomp      (opa) = isite%mmean_Af_decomp      (ipa)
            osite%mmean_rk4step        (opa) = isite%mmean_rk4step        (ipa)
            osite%mmean_albedo         (opa) = isite%mmean_albedo         (ipa)
            osite%mmean_albedo_beam    (opa) = isite%mmean_albedo_beam    (ipa)
            osite%mmean_albedo_diffuse (opa) = isite%mmean_albedo_diffuse (ipa)
         end if

         if (iqoutput > 0) then
            do icyc=1,ndcycle
               osite%qmean_rh             (icyc,opa) = isite%qmean_rh             (icyc,ipa)
               osite%qmean_albedo         (icyc,opa) = isite%qmean_albedo         (icyc,ipa)
               osite%qmean_albedo_beam    (icyc,opa) = isite%qmean_albedo_beam    (icyc,ipa)
               osite%qmean_albedo_diffuse (icyc,opa) = isite%qmean_albedo_diffuse (icyc,ipa)
            end do
         end if

         !----- Copy all cohorts. ---------------------------------------------------------!
         call allocate_patchtype(osite%patch(opa),isite%patch(ipa)%ncohorts)
         call copy_patchtype(isite%patch(ipa),osite%patch(opa),1,isite%patch(ipa)%ncohorts &
                                                              ,1,isite%patch(ipa)%ncohorts)
      end do opaloop

      return
   end subroutine copy_sitetype
   !=======================================================================================!
   !=======================================================================================!





!============================================================================!
!============================================================================!
  subroutine copy_sitetype_mask(sitein,siteout,logmask,masksz,newsz)

    ! This subroutine assumes that the size of vectors in siteout
    ! are the number of true elements in mask, while the size of the
    ! vectors in sitein are of the size of the mask itself. If this
    ! is not true, you will get a segmentation violation and the
    ! code will crash.args 1 and 3 must be dimension of arg 4
    ! argument 2 must be the dimension of the sum of the 3rd argument
    ! 
    ! THIS ROUTINE CURRENTLY ASSUMES THAT THE OUTPUT SITE
    ! HAS NOT ALLOCATED IT'S PATCH'S COHORT VECTORS YET, THIS
    ! IS BECAUSE THE LENGTHS OF THESE VECTORS ARE BASED ON THE
    ! DONOR PATH'S VECTOR SIZES.  DO NOT USE PRE-ALLOCATED
    ! RECIPIENTS
    use soil_bgc, only: copy_sitetype_mask_sbgc
    implicit none

    type(sitetype),target :: sitein,siteout
    integer :: masksz,newsz
    integer,dimension(newsz) :: incmask
    logical,dimension(masksz)           :: logmask
    integer :: i,k,m,inc,ipft,icyc
    type(stoma_data),pointer :: osdi,osdo

    inc = 0
    do i=1,masksz
       if (logmask(i)) then
          inc = inc + 1
          incmask(inc) = i
       end if
    end do

    ! First do all of the true vectors
    siteout%paco_id(1:inc)              = pack(sitein%paco_id,logmask)
    siteout%paco_n(1:inc)               = pack(sitein%paco_n,logmask)
    siteout%dist_type(1:inc)            = pack(sitein%dist_type,logmask)
    siteout%age(1:inc)                  = pack(sitein%age,logmask)
    siteout%area(1:inc)                 = pack(sitein%area,logmask)

    call copy_sitetype_mask_sbgc(masksz, logmask, siteout%sbgc, sitein%sbgc)

    siteout%sum_dgd(1:inc)              = pack(sitein%sum_dgd,logmask)
    siteout%sum_chd(1:inc)              = pack(sitein%sum_chd,logmask)
    siteout%plantation(1:inc)           = pack(sitein%plantation,logmask)
    siteout%cohort_count(1:inc)         = pack(sitein%cohort_count,logmask)
    siteout%can_theiv(1:inc)            = pack(sitein%can_theiv,logmask)
    siteout%can_temp(1:inc)             = pack(sitein%can_temp,logmask)
    siteout%can_shv(1:inc)              = pack(sitein%can_shv,logmask)
    siteout%can_co2(1:inc)              = pack(sitein%can_co2,logmask)
    siteout%can_rhos(1:inc)             = pack(sitein%can_rhos,logmask)
    siteout%can_prss(1:inc)             = pack(sitein%can_prss,logmask)
    siteout%can_theta(1:inc)            = pack(sitein%can_theta,logmask)
    siteout%can_depth(1:inc)            = pack(sitein%can_depth,logmask)
    siteout%opencan_frac(1:inc)         = pack(sitein%opencan_frac,logmask)
    siteout%ggbare(1:inc)               = pack(sitein%ggbare,logmask)
    siteout%ggveg(1:inc)                = pack(sitein%ggveg,logmask)
    siteout%ggnet(1:inc)                = pack(sitein%ggnet,logmask)
    siteout%ggsoil(1:inc)               = pack(sitein%ggsoil,logmask)
    siteout%lambda_light(1:inc)         = pack(sitein%lambda_light,logmask)
    siteout%lai(1:inc)                  = pack(sitein%lai,logmask)
    siteout%wpa(1:inc)                  = pack(sitein%wpa,logmask)
    siteout%wai(1:inc)                  = pack(sitein%wai,logmask)
    siteout%avg_daily_temp(1:inc)       = pack(sitein%avg_daily_temp,logmask)
    siteout%mean_rh(1:inc)              = pack(sitein%mean_rh,logmask)
    siteout%mean_nep(1:inc)             = pack(sitein%mean_nep,logmask)
    siteout%wbudget_loss2atm(1:inc)     = pack(sitein%wbudget_loss2atm,logmask)
    siteout%wbudget_denseffect(1:inc)        = pack(sitein%wbudget_denseffect,logmask)
    siteout%wbudget_precipgain(1:inc)        = pack(sitein%wbudget_precipgain,logmask)
    siteout%wbudget_loss2runoff(1:inc)       = pack(sitein%wbudget_loss2runoff,logmask)
    siteout%wbudget_loss2drainage(1:inc)     = pack(sitein%wbudget_loss2drainage,logmask)
    siteout%wbudget_initialstorage(1:inc)    = pack(sitein%wbudget_initialstorage,logmask)
    siteout%wbudget_residual(1:inc)          = pack(sitein%wbudget_residual,logmask)
    siteout%ebudget_loss2atm(1:inc)          = pack(sitein%ebudget_loss2atm,logmask)
    siteout%ebudget_denseffect(1:inc)        = pack(sitein%ebudget_denseffect,logmask)
    siteout%ebudget_loss2runoff(1:inc)       = pack(sitein%ebudget_loss2runoff,logmask)
    siteout%ebudget_loss2drainage(1:inc)     = pack(sitein%ebudget_loss2drainage,logmask)
    siteout%ebudget_netrad(1:inc)            = pack(sitein%ebudget_netrad,logmask)
    siteout%ebudget_precipgain(1:inc)        = pack(sitein%ebudget_precipgain,logmask)
    siteout%ebudget_initialstorage(1:inc)    = pack(sitein%ebudget_initialstorage,logmask)
    siteout%ebudget_residual(1:inc)          = pack(sitein%ebudget_residual,logmask)
    siteout%co2budget_initialstorage(1:inc)  = pack(sitein%co2budget_initialstorage,logmask)
    siteout%co2budget_residual(1:inc)   = pack(sitein%co2budget_residual,logmask)
    siteout%co2budget_loss2atm(1:inc)   = pack(sitein%co2budget_loss2atm,logmask)
    siteout%co2budget_denseffect(1:inc) = pack(sitein%co2budget_denseffect,logmask)
    siteout%co2budget_gpp(1:inc)        = pack(sitein%co2budget_gpp,logmask)
    siteout%co2budget_plresp(1:inc)     = pack(sitein%co2budget_plresp,logmask)
    siteout%co2budget_rh(1:inc)         = pack(sitein%co2budget_rh,logmask)
    siteout%today_A_decomp(1:inc)       = pack(sitein%today_A_decomp,logmask)
    siteout%today_Af_decomp(1:inc)      = pack(sitein%today_Af_decomp,logmask)
    siteout%veg_rough(1:inc)            = pack(sitein%veg_rough,logmask)
    siteout%veg_height(1:inc)           = pack(sitein%veg_height,logmask)
    siteout%veg_displace(1:inc)         = pack(sitein%veg_displace,logmask)
    siteout%fsc_in(1:inc)               = pack(sitein%fsc_in,logmask)
    siteout%ssc_in(1:inc)               = pack(sitein%ssc_in,logmask)
    siteout%ssl_in(1:inc)               = pack(sitein%ssl_in,logmask)
    siteout%fsn_in(1:inc)               = pack(sitein%fsn_in,logmask)



    siteout%fsp_in(1:inc)               = pack(sitein%fsp_in,logmask)
    siteout%stsp_in(1:inc)               = pack(sitein%stsp_in,logmask)
    siteout%total_plant_P_uptake(1:inc)    = pack(sitein%total_plant_P_uptake,logmask)



    siteout%total_plant_nitrogen_uptake(1:inc)    = pack(sitein%total_plant_nitrogen_uptake,logmask)
    siteout%mineralized_N_loss(1:inc)    = pack(sitein%mineralized_N_loss,logmask)
    siteout%mineralized_N_input(1:inc)    = pack(sitein%mineralized_N_input,logmask)
    siteout%rshort_g(1:inc)             = pack(sitein%rshort_g,logmask)
    siteout%rshort_g_beam(1:inc)        = pack(sitein%rshort_g_beam,logmask)
    siteout%rshort_g_diffuse(1:inc)     = pack(sitein%rshort_g_diffuse,logmask)
    siteout%rlong_g(1:inc)              = pack(sitein%rlong_g,logmask)
    siteout%rlong_g_surf(1:inc)         = pack(sitein%rlong_g_surf,logmask)
    siteout%rlong_g_incid(1:inc)        = pack(sitein%rlong_g_incid,logmask)
    siteout%rlong_s(1:inc)              = pack(sitein%rlong_s,logmask)
    siteout%rlong_s_surf(1:inc)         = pack(sitein%rlong_s_surf,logmask)
    siteout%rlong_s_incid(1:inc)        = pack(sitein%rlong_s_incid,logmask)
    siteout%albedo(1:inc)               = pack(sitein%albedo,logmask)
    siteout%albedo_beam(1:inc)          = pack(sitein%albedo_beam,logmask)
    siteout%albedo_diffuse(1:inc)       = pack(sitein%albedo_diffuse,logmask)
    siteout%rlongup(1:inc)              = pack(sitein%rlongup,logmask)
    siteout%rlong_albedo(1:inc)         = pack(sitein%rlong_albedo,logmask)
    siteout%total_sfcw_depth(1:inc)     = pack(sitein%total_sfcw_depth,logmask)
    siteout%snowfac(1:inc)              = pack(sitein%snowfac,logmask)
    siteout%A_decomp(1:inc)             = pack(sitein%A_decomp,logmask)
    siteout%f_decomp(1:inc)             = pack(sitein%f_decomp,logmask)
    siteout%rh(1:inc)                   = pack(sitein%rh,logmask)
    siteout%cwd_rh(1:inc)               = pack(sitein%cwd_rh,logmask)
    siteout%fuse_flag(1:inc)            = pack(sitein%fuse_flag,logmask)
    siteout%plant_ag_biomass(1:inc)     = pack(sitein%plant_ag_biomass,logmask)
    siteout%mean_wflux(1:inc)           = pack(sitein%mean_wflux,logmask)
    siteout%mean_latflux(1:inc)         = pack(sitein%mean_latflux,logmask)
    siteout%mean_hflux(1:inc)           = pack(sitein%mean_hflux,logmask)
    siteout%runoff(1:inc)               = pack(sitein%runoff,logmask)
    siteout%mean_runoff(1:inc)          = pack(sitein%mean_runoff,logmask)
    siteout%mean_qrunoff(1:inc)         = pack(sitein%mean_qrunoff,logmask)
    siteout%htry(1:inc)                 = pack(sitein%htry,logmask)
    siteout%avg_rk4step(1:inc)          = pack(sitein%avg_rk4step,logmask)
    siteout%avg_available_water(1:inc)  = pack(sitein%avg_available_water,logmask)
    siteout%ustar(1:inc)                = pack(sitein%ustar,logmask)
    siteout%tstar(1:inc)                = pack(sitein%tstar,logmask)
    siteout%qstar(1:inc)                = pack(sitein%qstar,logmask)
    siteout%cstar(1:inc)                = pack(sitein%cstar,logmask)
    
    siteout%zeta(1:inc)                 = pack(sitein%zeta,logmask)
    siteout%ribulk(1:inc)               = pack(sitein%ribulk,logmask)

    siteout%upwp(1:inc)                 = pack(sitein%upwp,logmask)
    siteout%tpwp(1:inc)                 = pack(sitein%tpwp,logmask)
    siteout%qpwp(1:inc)                 = pack(sitein%qpwp,logmask)
    siteout%cpwp(1:inc)                 = pack(sitein%cpwp,logmask)
    siteout%wpwp(1:inc)                 = pack(sitein%wpwp,logmask)

    siteout%nlev_sfcwater(1:inc)        = pack(sitein%nlev_sfcwater,logmask)
    siteout%ground_shv(1:inc)           = pack(sitein%ground_shv,logmask)
    siteout%ground_ssh(1:inc)           = pack(sitein%ground_ssh,logmask)
    siteout%ground_temp(1:inc)          = pack(sitein%ground_temp,logmask)
    siteout%ground_fliq(1:inc)          = pack(sitein%ground_fliq,logmask)
    siteout%rough(1:inc)                = pack(sitein%rough,logmask)

    siteout%avg_rshort_gnd     (1:inc)  = pack(sitein%avg_rshort_gnd    ,logmask)
    siteout%avg_rlong_gnd      (1:inc)  = pack(sitein%avg_rlong_gnd     ,logmask)
    siteout%avg_carbon_ac      (1:inc)  = pack(sitein%avg_carbon_ac     ,logmask)
    siteout%avg_rlongup        (1:inc)  = pack(sitein%avg_rlongup       ,logmask)
    siteout%avg_albedo         (1:inc)  = pack(sitein%avg_albedo        ,logmask)
    siteout%avg_albedo_beam    (1:inc)  = pack(sitein%avg_albedo_beam   ,logmask)
    siteout%avg_albedo_diffuse (1:inc)  = pack(sitein%avg_albedo_diffuse,logmask)
    siteout%avg_rlong_albedo   (1:inc)  = pack(sitein%avg_rlong_albedo  ,logmask)

    siteout%avg_vapor_lc(1:inc)         = pack(sitein%avg_vapor_lc,logmask)
    siteout%avg_vapor_wc(1:inc)         = pack(sitein%avg_vapor_wc,logmask)
    siteout%avg_dew_cg(1:inc)           = pack(sitein%avg_dew_cg,logmask)
    siteout%avg_vapor_gc(1:inc)         = pack(sitein%avg_vapor_gc,logmask)
    siteout%avg_wshed_vg(1:inc)         = pack(sitein%avg_wshed_vg,logmask)
    siteout%avg_intercepted(1:inc)      = pack(sitein%avg_intercepted,logmask)
    siteout%avg_throughfall(1:inc)      = pack(sitein%avg_throughfall,logmask)
    siteout%avg_vapor_ac(1:inc)         = pack(sitein%avg_vapor_ac,logmask)
    siteout%avg_transp(1:inc)           = pack(sitein%avg_transp,logmask)
    siteout%avg_evap(1:inc)             = pack(sitein%avg_evap,logmask)
    siteout%avg_runoff(1:inc)           = pack(sitein%avg_runoff,logmask)
    siteout%avg_drainage(1:inc)         = pack(sitein%avg_drainage,logmask)
    siteout%avg_drainage_heat(1:inc)    = pack(sitein%avg_drainage_heat,logmask)
    siteout%aux(1:inc)                  = pack(sitein%aux,logmask)
    siteout%avg_sensible_lc(1:inc)      = pack(sitein%avg_sensible_lc,logmask)
    siteout%avg_sensible_wc(1:inc)      = pack(sitein%avg_sensible_wc,logmask)
    siteout%avg_qwshed_vg(1:inc)        = pack(sitein%avg_qwshed_vg,logmask)
    siteout%avg_qintercepted(1:inc)     = pack(sitein%avg_qintercepted,logmask)
    siteout%avg_qthroughfall(1:inc)     = pack(sitein%avg_qthroughfall,logmask)
    siteout%avg_sensible_gc(1:inc)      = pack(sitein%avg_sensible_gc,logmask)
    siteout%avg_sensible_ac(1:inc)      = pack(sitein%avg_sensible_ac,logmask)
    siteout%avg_runoff_heat(1:inc)      = pack(sitein%avg_runoff_heat,logmask)
    siteout%avg_leaf_energy(1:inc)      = pack(sitein%avg_leaf_energy,logmask)
    siteout%avg_leaf_temp(1:inc)        = pack(sitein%avg_leaf_temp,logmask)
    siteout%avg_leaf_hcap(1:inc)        = pack(sitein%avg_leaf_hcap,logmask)
    siteout%avg_leaf_fliq(1:inc)        = pack(sitein%avg_leaf_fliq,logmask)
    siteout%avg_leaf_water(1:inc)       = pack(sitein%avg_leaf_water,logmask)
    siteout%avg_wood_energy(1:inc)       = pack(sitein%avg_wood_energy,logmask)
    siteout%avg_wood_temp(1:inc)         = pack(sitein%avg_wood_temp,logmask)
    siteout%avg_wood_hcap(1:inc)         = pack(sitein%avg_wood_hcap,logmask)
    siteout%avg_wood_fliq(1:inc)         = pack(sitein%avg_wood_fliq,logmask)
    siteout%avg_wood_water(1:inc)        = pack(sitein%avg_wood_water,logmask)
    siteout%par_l_max(1:inc)            = pack(sitein%par_l_max,logmask)
    siteout%par_l_beam_max(1:inc)       = pack(sitein%par_l_beam_max,logmask)
    siteout%par_l_diffuse_max(1:inc)    = pack(sitein%par_l_diffuse_max,logmask)

    ! Water layers 1:nzs
    
    do k=1,nzs
       siteout%sfcwater_mass(k,1:inc)    = pack(sitein%sfcwater_mass(k,:),logmask)
       siteout%sfcwater_energy(k,1:inc)  = pack(sitein%sfcwater_energy(k,:),logmask)
       siteout%sfcwater_depth(k,1:inc)   = pack(sitein%sfcwater_depth(k,:),logmask)
       siteout%rshort_s(k,1:inc)         = pack(sitein%rshort_s(k,:),logmask)
       siteout%rshort_s_beam(k,1:inc)    = pack(sitein%rshort_s_beam(k,:),logmask)
       siteout%rshort_s_diffuse(k,1:inc) = pack(sitein%rshort_s_diffuse(k,:),logmask)
       siteout%sfcwater_tempk(k,1:inc)   = pack(sitein%sfcwater_tempk(k,:),logmask)
       siteout%sfcwater_fracliq(k,1:inc)  = pack(sitein%sfcwater_fracliq(k,:),logmask)
    end do

    ! Soil layers 1:nzg

    do k=1,nzg
       siteout%soil_energy(k,1:inc)        = pack(sitein%soil_energy(k,:),logmask)
       siteout%soil_water(k,1:inc)         = pack(sitein%soil_water(k,:),logmask)
       siteout%soil_tempk(k,1:inc)         = pack(sitein%soil_tempk(k,:),logmask)
       siteout%soil_fracliq(k,1:inc)       = pack(sitein%soil_fracliq(k,:),logmask)
       siteout%rootdense(k,1:inc)          = pack(sitein%rootdense(k,:),logmask)
       siteout%avg_smoist_gg(k,1:inc)      = pack(sitein%avg_smoist_gg(k,:),logmask)
       siteout%avg_transloss(k,1:inc)      = pack(sitein%avg_transloss(k,:),logmask)
       siteout%avg_sensible_gg(k,1:inc)    = pack(sitein%avg_sensible_gg(k,:),logmask)
       siteout%aux_s(k,1:inc)              = pack(sitein%aux_s(k,:),logmask)

       siteout%current_paw(k,1:inc) = pack(sitein%current_paw(k,:),logmask)
       do m = 1, 10
          siteout%past_paw(k,m,1:inc) = pack(sitein%past_paw(k,m,:),logmask)
       enddo
    end do

    ! pft types

    do k=1,n_pft
       siteout%repro(k,1:inc)            = pack(sitein%repro(k,:),logmask)
       siteout%A_o_max(k,1:inc)          = pack(sitein%A_o_max(k,:),logmask)
       siteout%A_c_max(k,1:inc)          = pack(sitein%A_c_max(k,:),logmask)
       
       do m=1,n_stoma_atts
          siteout%old_stoma_vector_max(m,k,1:inc) = pack(sitein%old_stoma_vector_max(m,k,:),logmask)
       end do

       do m=1,ff_nhgt
          siteout%cumlai_profile(k,m,1:inc)       = pack(sitein%cumlai_profile(k,m,:),logmask)
       end do
    end do
    
    !dbh types
    do k=1,n_dbh
        siteout%co2budget_gpp_dbh(k,1:inc)        = pack(sitein%co2budget_gpp_dbh(k,:),logmask)
    end do

    
    ! Old_stoma_data_max type
    ! Derived type with n_pft x n_patch, so.... the intrinsic "pack" wont work

    do m=1,newsz
       k=incmask(m)
       call allocate_patchtype(siteout%patch(m),sitein%patch(k)%ncohorts)
       
       call copy_patchtype(sitein%patch(k),siteout%patch(m),1,sitein%patch(k)%ncohorts,1,sitein%patch(k)%ncohorts)

       do ipft=1,n_pft

          osdo => siteout%old_stoma_data_max(ipft,m)
          osdi => sitein%old_stoma_data_max(ipft,k)
          
          osdo%recalc           = osdi%recalc
          osdo%T_L              = osdi%T_L
          osdo%e_A              = osdi%e_A
          osdo%PAR              = osdi%PAR
          osdo%rb_factor        = osdi%rb_factor
          osdo%prss             = osdi%prss
          osdo%phenology_factor = osdi%phenology_factor
          osdo%gsw_open         = osdi%gsw_open
          osdo%ilimit           = osdi%ilimit
          osdo%T_L_residual     = osdi%T_L_residual
          osdo%e_a_residual     = osdi%e_a_residual
          osdo%par_residual     = osdi%par_residual
          osdo%rb_residual      = osdi%rb_residual
          osdo%leaf_residual    = osdi%leaf_residual
          osdo%gsw_residual     = osdi%gsw_residual

       end do
       
    end do
       
    if (idoutput > 0 .or. imoutput > 0 .or. iqoutput > 0) then
       siteout%dmean_rh             (1:inc) = pack(sitein%dmean_rh             ,logmask)
       siteout%dmean_co2_residual   (1:inc) = pack(sitein%dmean_co2_residual   ,logmask)
       siteout%dmean_energy_residual(1:inc) = pack(sitein%dmean_energy_residual,logmask)
       siteout%dmean_water_residual (1:inc) = pack(sitein%dmean_water_residual ,logmask)
       siteout%dmean_lambda_light   (1:inc) = pack(sitein%dmean_lambda_light   ,logmask)
       siteout%dmean_A_decomp       (1:inc) = pack(sitein%dmean_A_decomp       ,logmask)
       siteout%dmean_Af_decomp      (1:inc) = pack(sitein%dmean_Af_decomp      ,logmask)
       siteout%dmean_rk4step        (1:inc) = pack(sitein%dmean_rk4step        ,logmask)
       siteout%dmean_albedo         (1:inc) = pack(sitein%dmean_albedo         ,logmask)
       siteout%dmean_albedo_beam    (1:inc) = pack(sitein%dmean_albedo_beam    ,logmask)
       siteout%dmean_albedo_diffuse (1:inc) = pack(sitein%dmean_albedo_diffuse ,logmask)
    end if
    
    if (imoutput > 0 .or. iqoutput > 0) then
       siteout%mmean_rh             (1:inc) = pack(sitein%mmean_rh             ,logmask)
       siteout%mmean_co2_residual   (1:inc) = pack(sitein%mmean_co2_residual   ,logmask)
       siteout%mmean_energy_residual(1:inc) = pack(sitein%mmean_energy_residual,logmask)
       siteout%mmean_water_residual (1:inc) = pack(sitein%mmean_water_residual ,logmask)
       siteout%mmean_lambda_light   (1:inc) = pack(sitein%mmean_lambda_light   ,logmask)
       siteout%mmean_A_decomp       (1:inc) = pack(sitein%mmean_A_decomp       ,logmask)
       siteout%mmean_Af_decomp      (1:inc) = pack(sitein%mmean_Af_decomp      ,logmask)
       siteout%mmean_rk4step        (1:inc) = pack(sitein%mmean_rk4step        ,logmask)
       siteout%mmean_albedo         (1:inc) = pack(sitein%mmean_albedo         ,logmask)
       siteout%mmean_albedo_beam    (1:inc) = pack(sitein%mmean_albedo_beam    ,logmask)
       siteout%mmean_albedo_diffuse (1:inc) = pack(sitein%mmean_albedo_diffuse ,logmask)
    end if
    
    if (iqoutput > 0) then
       do icyc=1,ndcycle
          siteout%qmean_rh             (icyc,1:inc) = pack(sitein%qmean_rh            (icyc,:) ,logmask)
          siteout%qmean_albedo         (icyc,1:inc) = pack(sitein%qmean_albedo        (icyc,:) ,logmask)
          siteout%qmean_albedo_beam    (icyc,1:inc) = pack(sitein%qmean_albedo_beam   (icyc,:) ,logmask)
          siteout%qmean_albedo_diffuse (icyc,1:inc) = pack(sitein%qmean_albedo_diffuse(icyc,:) ,logmask)
       end do
    end if

    return
  end subroutine copy_sitetype_mask
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine copy_patchtype_mask(patchin,patchout,mask,masksz,newsz)

    ! This subroutine assumes that the size of vectors in siteout
    ! are the number of true elements in mask, while the size of the
    ! vectors in sitein are of the size of the mask itself. If this
    ! is not true, you will get a segmentation violation and the
    ! code will crash.args 1 and 3 must be dimension of arg 4
    ! argument 2 must be the dimension of the sum of the 3rd argument
    ! 
    ! THIS ROUTINE CURRENTLY ASSUMES THAT THE OUTPUT SITE
    ! HAS NOT ALLOCATED IT'S PATCH'S COHORT VECTORS YET, THIS
    ! IS BECAUSE THE LENGTHS OF THESE VECTORS ARE BASED ON THE
    ! DONOR PATH'S VECTOR SIZES.  DO NOT USE PRE-ALLOCATED
    ! RECIPIENTS
    
    use cohort_state, only: copy_patchtype_mask_state
    use cohort_phen, only: copy_patchtype_mask_phen
    use cohort_mort, only: copy_patchtype_mask_mort
    use cohort_resp, only: copy_patchtype_mask_resp
    use cohort_photo, only: copy_patchtype_mask_photo
    use cohort_rad, only: copy_patchtype_mask_rad
    use cohort_therm, only: copy_patchtype_mask_therm
    implicit none

    type(patchtype),target :: patchin,patchout
    integer :: masksz,newsz
    integer,dimension(newsz) :: incmask
    integer,dimension(masksz):: imask
    logical,dimension(masksz)  :: mask
    integer :: i,k,m,inc
    type(stoma_data),pointer :: osdi,osdo

    do i=1,masksz
       imask(i) = i
    end do
    inc=count(mask)                   ! Number of true elements

    if (inc == 0) return

    incmask=pack(imask,mask)   ! List of true elements

    call copy_patchtype_mask_state(patchout%costate, inc, patchin%costate, mask, masksz)
    call copy_patchtype_mask_phen(patchout%cophen, inc, patchin%cophen, mask, masksz)
    call copy_patchtype_mask_mort(patchout%comort, inc, patchin%comort, mask, masksz, incmask,   &
         newsz, n_mort)
    call copy_patchtype_mask_resp(patchout%coresp, inc, patchin%coresp, mask, masksz, ndcycle)
    call copy_patchtype_mask_photo(patchout%cophoto, inc, patchin%cophoto, mask, masksz, ndcycle)
    call copy_patchtype_mask_rad(patchout%corad, inc, patchin%corad, mask, masksz)
    call copy_patchtype_mask_therm(patchout%cotherm, inc, patchin%cotherm, mask, masksz)
    
    return
  end subroutine copy_patchtype_mask
!============================================================================!
!============================================================================!





!============================================================================!
!============================================================================!
  subroutine copy_patchtype(patchin,patchout,ipin1,ipin2,ipout1,ipout2)

    use cohort_state, only: copy_patchtype_state
    use cohort_phen, only: copy_patchtype_phen
    use cohort_mort, only: copy_patchtype_mort
    use cohort_photo, only: copy_patchtype_photo
    use cohort_resp, only: copy_patchtype_resp
    use cohort_rad, only: copy_patchtype_rad
    use cohort_therm, only: copy_patchtype_therm

    implicit none
    integer :: ipin1,ipin2,ipout1,ipout2
    type(patchtype),target :: patchin,patchout
    type(stoma_data),pointer :: osdo,osdi
    integer :: iout,iin

    if (ipout2-ipout1.ne.ipin2-ipin1) then
       print*,"In copy_patchtype:"
       print*,"You specified unequal vector lengths"
       print*,"in the input and output targets"
       print*,"This cannot be..stopping"
       call fatal_error('unequal vector lengths','copy_patchtype','ed_state_vars.f90')
    end if

    call copy_patchtype_state(patchin%costate,patchout%costate,ipin1,ipin2,ipout1,ipout2)
    call copy_patchtype_phen(patchin%cophen,patchout%cophen,ipin1,ipin2,ipout1,ipout2)
    call copy_patchtype_mort(patchin%comort,patchout%comort,ipin1,ipin2,ipout1,ipout2)
    call copy_patchtype_resp(patchin%coresp,patchout%coresp,ipin1,ipin2,ipout1,ipout2)
    call copy_patchtype_photo(patchin%cophoto,patchout%cophoto,ipin1,ipin2,ipout1,ipout2)
    call copy_patchtype_rad(patchin%corad,patchout%corad,ipin1,ipin2,ipout1,ipout2)
    call copy_patchtype_therm(patchin%cotherm,patchout%cotherm,ipin1,ipin2,ipout1,ipout2)

    return
  end subroutine copy_patchtype

!============================================================================!
!============================================================================!
  

  ! ===============================================================
  ! Define the vtables of the state/output variables
  !
  ! The various state scalars, vectors and arrays are
  ! now populate the vtable.  The vtable indexes the array
  ! gives it a name, records its dimensions, when it is to be
  ! used as output and how (averaging and such) and most importantly
  ! saves a pointer to its starting position.  If this routine is
  ! being called as a compute node in parallel, the first position
  ! is not necessarily the first position of the whole datavector,
  ! but will only be the first position of that nodes hyperslab chunk
  ! within the given continuous dataset.
  !
  ! The first number correspond to the data level:
  ! 1. Gridtype     (polygon level)
  ! 2. Polygontype  (site level)
  ! 3. Sitetype     (patch level)
  ! 4. Patchtype    (cohort level)

  ! 9. Scalar
  ! The other numbers correspond to the kind of dimension and variable.
  ! 0. Main vector ordinate only, integer.
  ! 1. Main vector ordinate only, real.
  ! 2. Soil layer
  ! 3. Surface water layer
  ! 4. PFT
  ! 5. Disturbance Type
  ! 6. DBH class
  ! 7. FF_DBH class
  ! 8. Mortality
  ! 9. Month/13 months
  !
  ! An extra dimension for diurnal cycle can be denoted by the negative sign
  !
  ! Of these possible dimensions (2-9), they may be used concurrently
  ! to partition the data into multi-dimensional spaces, but all seven
  ! will not be used simultaneously.  The highest ranks in use are 3.
  ! Each unique combination will have a call number associated with it.
  !
  !  10    : rank 1 : polygon, integer
  !  11    : rank 1 : polygon
  ! -11    : rank 2 : polygon, diurnal cycle
  !  12    : rank 2 : polygon, s-layer
  !  120   : rank 2 : polygon, s-layer, integer
  ! -12    : rank 3 : polygon, s-layer, diurnal cycle
  !  13    : rank 2 : polygon, w-layer
  !  14    : rank 2 : polygon, pft
  !  14567 : rank 5 : polygon, pft, disturbance, dbh, age
  !  146   : rank 3 : polygon, pft, dbh
  !  15    : rank 2 : polygon, disturbance
  !  155   : rank 3 : polygon, disturbance, disturbance
  !  157   : rank 3 : polygon, disturbance, age
  !  16    : rank 2 : polygon, dbh
  !  17    : rank 2 : polygon, age
  !  18    : rank 2 : polygon, mort
  !  19    : rank 2 : polygon, month+1
  !
  !  20    : rank 1 : site, integer
  !  21    : rank 1 : site
  !  22    : rank 2 : site, s-layer
  !  23    : rank 2 : site, w-layer
  !  24    : rank 2 : site, pft
  !  246   : rank 3 : site, pft, dbh
  !  25    : rank 2 : site, disturbance
  !  255   : rank 3 : site, disturbance, disturbance
  !  26    : rank 2 : site, dbh
  !  27    : rank 2 : site, age
  !  28    : rank 2 : site, mort
  !  29    : rank 2 : site, month
  !
  !  30    : rank 1 : patch, integer
  !  31    : rank 1 : patch
  ! -31    : rank 2 : patch, diurnal cycle
  !  32    : rank 2 : patch, s-layer
  !  33    : rank 2 : patch, w-layer
  !  34    : rank 2 : patch, pft
  !  346   : rank 3 : patch, pft, ff_dbh
  !  35    : rank 2 : patch, disturbance
  !  36    : rank 2 : patch, dbh
  !  37    : rank 2 : patch, age
  !  38    : rank 2 : patch, mort
  !
  !  40    : rank 1 : cohort, integer
  !  41    : rank 1 : cohort
  ! -41    : rank 2 : cohort, diurnal cycle
  !  44    : rank 2 : cohort, pft
  !  46    : rank 2 : cohort, dbh
  !  47    : rank 2 : cohort, age
  !  48    : rank 2 : cohort, mort
  !  49    : rank 2 : cohort, month+1
  !
  !  90    : rank 0 : integer scalar 
  !  92    : rank 1 : s-layer
  !===================================================================
  
  subroutine filltab_alltypes

    ! =================================================
    !
    ! This subroutine is the main driver for filling
    ! filling the var_table of ED variables.  On a
    ! serial computing environment, this routine should be
    ! called near the end of the initialization process
    ! after the hierarchical tree structure has been
    ! trimmed via fusion/fission.  Similiarly, this
    ! routine should be called after any fusion/fission
    ! process, assuming that the major vtable structures
    ! have been deallocated prior to reallocation.
    !
    ! In a paralell environment, this routine should
    ! operate in a similiar fashion on each of the compute
    ! nodes.  It is designed such that the compute nodes
    ! will write hyperslabs of data in parallel to a
    ! joing HDF5 dataset as "collective-chunked" data
    ! The modifications that must be made after running this
    ! subroutine, are that the indexes should account
    ! for the offset of the current compute node.
    !
    ! =================================================
    
    
    use ed_var_tables,only:num_var,vt_info,var_table,nullify_vt_vector_pointers
    use ed_node_coms,only:mynum,mchnum,machs,nmachs,nnodetot,sendnum,recvnum,master_num
    use ed_max_dims, only: maxgrds, maxmach
    implicit none
    
    include 'mpif.h'

    integer :: ncohorts_g,npatches_g,nsites_g
    integer :: igr,ipy,isi,ipa,nv,ierr,nm,iptr
    integer,       dimension(MPI_STATUS_SIZE) :: status
    integer :: ping,uniqueid
    logical,save :: model_start = .true.
   
    type(edtype),pointer      :: cgrid
    type(polygontype),pointer :: cpoly
    type(sitetype),pointer    :: csite
    type(patchtype),pointer   :: cpatch
    logical :: verbose = .false.    

!    if (mynum == 1) then
!       write(*,"(a)")'--- Re-hashing the IO pointer tables and mapping arrays'
!    end if


    ! The first loop through populates the info tables

    do igr = 1,ngrids
       cgrid => edgrid_g(igr)
       
       if (num_var(igr)>0) then
          do nv=1,num_var(igr)
             if (associated(vt_info(nv,igr)%vt_vector)) then
                do iptr=1,vt_info(nv,igr)%nptrs
                   call nullify_vt_vector_pointers(vt_info(nv,igr)%vt_vector(iptr))
                end do
                deallocate(vt_info(nv,igr)%vt_vector)
             end if
          end do
       end if

       num_var(igr) = 0

       cgrid%npolygons_global = cgrid%npolygons
       cgrid%nsites_global    = get_nsites(cgrid)
       cgrid%npatches_global  = get_npatches(cgrid)
       cgrid%ncohorts_global  = get_ncohorts(cgrid)
       
       cgrid%mach_cohort_offset_index = 0
       cgrid%mach_patch_offset_index  = 0
       cgrid%mach_site_offset_index   = 0
       cgrid%mach_polygon_offset_index= 0

       if (nnodetot /= 1) then

          ! Send all them sizes to root (CHANGED, NODE 1)
          

          if (mynum == 1) then
          
             gdpy(1,igr) = cgrid%npolygons_global
             gdsi(1,igr) = cgrid%nsites_global
             gdpa(1,igr) = cgrid%npatches_global
             gdco(1,igr) = cgrid%ncohorts_global

             call MPI_Send(ping,1,MPI_INTEGER,sendnum,94,MPI_COMM_WORLD,ierr)
             
             ! Have node 1 recieve the info
             do nm=2,nnodetot
                uniqueid=((igr-1)*maxmach)+nm
                call MPI_Recv(gdpy(nm,igr),1,MPI_INTEGER,machs(nm),500000+uniqueid,MPI_COMM_WORLD,status,ierr)
                call MPI_Recv(gdsi(nm,igr),1,MPI_INTEGER,machs(nm),600000+uniqueid,MPI_COMM_WORLD,status,ierr)
                call MPI_Recv(gdpa(nm,igr),1,MPI_INTEGER,machs(nm),700000+uniqueid,MPI_COMM_WORLD,status,ierr)
                call MPI_Recv(gdco(nm,igr),1,MPI_INTEGER,machs(nm),800000+uniqueid,MPI_COMM_WORLD,status,ierr)
             end do

             ! Broadcast all this info to the nodes
             do nm=2,nnodetot
                uniqueid=((igr-1)*maxmach)+nm
                call MPI_Send(gdpy,maxmach*maxgrds,MPI_INTEGER,machs(nm), 900000+uniqueid,MPI_COMM_WORLD,ierr)
                call MPI_Send(gdsi,maxmach*maxgrds,MPI_INTEGER,machs(nm),1000000+uniqueid,MPI_COMM_WORLD,ierr)
                call MPI_Send(gdpa,maxmach*maxgrds,MPI_INTEGER,machs(nm),1100000+uniqueid,MPI_COMM_WORLD,ierr)
                call MPI_Send(gdco,maxmach*maxgrds,MPI_INTEGER,machs(nm),1200000+uniqueid,MPI_COMM_WORLD,ierr)
             end do

         else

            ! Set the blocking recieve to allow ordering, start with machine 1
            call MPI_Recv(ping,1,MPI_INTEGER,recvnum,94,MPI_COMM_WORLD,status,ierr)

            uniqueid=((igr-1)*maxmach)+mynum
            ! Send the information to node (1)
            call MPI_Send(cgrid%npolygons_global, 1,MPI_INTEGER,machs(1),500000+uniqueid,MPI_COMM_WORLD,ierr)
            call MPI_Send(cgrid%nsites_global   , 1,MPI_INTEGER,machs(1),600000+uniqueid,MPI_COMM_WORLD,ierr)
            call MPI_Send(cgrid%npatches_global , 1,MPI_INTEGER,machs(1),700000+uniqueid,MPI_COMM_WORLD,ierr)
            call MPI_Send(cgrid%ncohorts_global , 1,MPI_INTEGER,machs(1),800000+uniqueid,MPI_COMM_WORLD,ierr)
          
            ! When this node is finished, send the blocking MPI_Send to the next machine
            if (mynum /= nnodetot) call MPI_Send(ping,1,MPI_INTEGER,sendnum,94,MPI_COMM_WORLD,ierr)

            uniqueid=((igr-1)*maxmach)+mynum
            call MPI_Recv(gdpy,maxmach*maxgrds,MPI_INTEGER,machs(1), 900000+uniqueid,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(gdsi,maxmach*maxgrds,MPI_INTEGER,machs(1),1000000+uniqueid,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(gdpa,maxmach*maxgrds,MPI_INTEGER,machs(1),1100000+uniqueid,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(gdco,maxmach*maxgrds,MPI_INTEGER,machs(1),1200000+uniqueid,MPI_COMM_WORLD,status,ierr)
 
         end if


         if(mynum == 1 .and. model_start .and. verbose) then
            
            print*,"Global Polygons: ",gdpy(1:nnodetot,igr)
            print*,"Global Site: "    ,gdsi(1:nnodetot,igr)
            print*,"Global Patches: " ,gdpa(1:nnodetot,igr)
            print*,"Global Cohorts: " ,gdco(1:nnodetot,igr)

         end if

         ! Calculate the offsets that each machine has
         
         py_off(1,igr) = 0
         si_off(1,igr) = 0
         pa_off(1,igr) = 0
         co_off(1,igr) = 0
         do nm=2,nnodetot
            py_off(nm,igr) = py_off(nm-1,igr) + gdpy(nm-1,igr)
            si_off(nm,igr) = si_off(nm-1,igr) + gdsi(nm-1,igr)
            pa_off(nm,igr) = pa_off(nm-1,igr) + gdpa(nm-1,igr)
            co_off(nm,igr) = co_off(nm-1,igr) + gdco(nm-1,igr)
         end do
         
         ! Calculate the total sizes of the arrays
         
         cgrid%npolygons_global = sum(gdpy(1:nnodetot,igr))
         cgrid%nsites_global    = sum(gdsi(1:nnodetot,igr))
         cgrid%npatches_global  = sum(gdpa(1:nnodetot,igr))
         cgrid%ncohorts_global  = sum(gdco(1:nnodetot,igr))

         ! Calculate the local offsets
         
         cgrid%mach_polygon_offset_index = py_off(mynum,igr)
         cgrid%mach_site_offset_index    = si_off(mynum,igr)
         cgrid%mach_patch_offset_index   = pa_off(mynum,igr)
         cgrid%mach_cohort_offset_index  = co_off(mynum,igr)

          
       end if

       call filltab_globtype(igr)

       call filltab_edtype(igr,0)
       
       if (gdpy(mynum,igr)>0) then
          call filltab_polygontype(igr,1,0)
          call filltab_sitetype(igr,1,1,0)
          call filltab_patchtype(igr,1,1,1,0)
       end if
       
       
    end do


    do igr = 1,ngrids

       ! Test to see if the var_table has been initialized. If it has
       ! then deallocate its pointers and reset its first flag. These
       ! will be reallocated on the first pass of the filltab_
       ! subroutines.

       cgrid => edgrid_g(igr)

       cgrid%pyglob_id = 0 + cgrid%mach_polygon_offset_index
       
       ! Determine the total number of variables for each grid
       ! These will determine the length of the vt_vector

       call filltab_edtype(igr,1)
       
       ncohorts_g = 0 + cgrid%mach_cohort_offset_index
       npatches_g = 0 + cgrid%mach_patch_offset_index
       nsites_g   = 0 + cgrid%mach_site_offset_index
       
       do ipy = 1,cgrid%npolygons
          
          cpoly => cgrid%polygon(ipy)
          
          cpoly%siglob_id = nsites_g + 0 ! This is the offset for the vtable write
          
          cgrid%pysi_id(ipy) = nsites_g + 1 ! This is the index written in the file
                                            ! for the user to reference

          cgrid%pysi_n(ipy) = cpoly%nsites

          nsites_g = nsites_g + cpoly%nsites
          
          call filltab_polygontype(igr,ipy,1)
          
          do isi = 1,cpoly%nsites
             
             csite => cpoly%site(isi)
             
             csite%paglob_id = npatches_g + 0

             cpoly%sipa_id(isi) = npatches_g + 1

             cpoly%sipa_n(isi) = csite%npatches

             npatches_g = npatches_g + csite%npatches
             
             call filltab_sitetype(igr,ipy,isi,1)
             
             do ipa = 1,csite%npatches

                cpatch => csite%patch(ipa)
                
                cpatch%coglob_id = ncohorts_g + 0

                csite%paco_id(ipa) = ncohorts_g + 1

                csite%paco_n(ipa) = cpatch%ncohorts

                ncohorts_g = ncohorts_g + cpatch%ncohorts

                if (cpatch%ncohorts > 0 ) then
                   
                   call filltab_patchtype(igr,ipy,isi,ipa,1)
                   
                end if

             end do
             
          end do
          
       end do
!       if (mynum.eq.1) then
!          write(*,"(a)")'--- Mapping Completed'
!       end if

       if (mynum.eq.1 .and. model_start .and. verbose) then
          model_start = .false.
          do nv=1,num_var(igr)
!             write(*,"(a,i4,a,i4,a,a)")'Registering: ',nv,' of',num_var(igr),'  ',vt_info(nv,igr)%name
          end do
       end if 

    end do


    return
  end subroutine filltab_alltypes
!==========================================================================================!
!==========================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine filltab_globtype(igr)

      use ed_var_tables, only : vtable_edio_r       & ! sub-routine
                              , vtable_edio_r_sca   & ! sub-rouitne
                              , vtable_edio_i_sca   ! ! sub-rouitne
      use soil_coms    , only : slz                 & ! intent(in)
                              , slxclay             & ! intent(in)
                              , slxsand             & ! intent(in)
                              , isoilflg            ! ! intent(in)


      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: igr
      !----- Local variables. -------------------------------------------------------------!
      integer                  :: var_len
      integer                  :: max_ptrs
      integer                  :: var_len_global
      integer                  :: nvar
      type(edtype), pointer    :: cgrid
      !------------------------------------------------------------------------------------!


      cgrid => edgrid_g(igr)


      !------------------------------------------------------------------------------------!
      !     Single values (scalars).                                                       !
      !------------------------------------------------------------------------------------!
      var_len        = 1
      var_len_global = 1
      max_ptrs       = 1

      nvar=1
      call vtable_edio_i_sca(cgrid%npolygons_global,nvar,igr,0,0                           &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'NPOLYGONS_GLOBAL :90:hist:anal:dail:mont:dcyc:year')
      call vtable_edio_i_sca(cgrid%npolygons_global,nvar,igr,1,0                           &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'NPOLYGONS_GLOBAL :90:hist:anal:dail:mont:dcyc:year')
      
      nvar=nvar+1
      call vtable_edio_i_sca(cgrid%nsites_global,nvar,igr,0,0                              &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'NSITES_GLOBAL :90:hist:anal:dail:mont:dcyc:year')
      call vtable_edio_i_sca(cgrid%nsites_global,nvar,igr,1,0                              &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'NSITES_GLOBAL :90:hist:anal:dail:mont:dcyc:year')
      
      nvar=nvar+1
      call vtable_edio_i_sca(cgrid%npatches_global,nvar,igr,0,0                            &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'NPATCHES_GLOBAL :90:hist:anal:dail:mont:dcyc:year')
      call vtable_edio_i_sca(cgrid%npatches_global,nvar,igr,1,0                            &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'NPATCHES_GLOBAL :90:hist:anal:dail:mont:dcyc:year')

      nvar=nvar+1
      call vtable_edio_i_sca(cgrid%ncohorts_global,nvar,igr,0,0                            &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'NCOHORTS_GLOBAL :90:hist:anal:dail:mont:dcyc:year')
      call vtable_edio_i_sca(cgrid%ncohorts_global,nvar,igr,1,0                            &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'NCOHORTS_GLOBAL :90:hist:anal:dail:mont:dcyc:year')

      nvar=nvar+1
      call vtable_edio_i_sca(nzg,nvar,igr,0,0                                              &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'NZG :90:hist:anal:dail:mont:dcyc:year')
      call vtable_edio_i_sca(nzg,nvar,igr,1,0                                              &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'NZG :90:hist:anal:dail:mont:dcyc:year')

      nvar=nvar+1
      call vtable_edio_i_sca(nzs,nvar,igr,0,0                                              &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'NZS :90:hist:anal:dail:mont:dcyc:year')
      call vtable_edio_i_sca(nzs,nvar,igr,1,0                                              &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'NZS :90:hist:anal:dail:mont:dcyc:year')

      nvar=nvar+1
      call vtable_edio_i_sca(ff_nhgt,nvar,igr,0,0                                              &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'FF_NHGT :90:hist:anal:dail:mont:dcyc:year')
      call vtable_edio_i_sca(ff_nhgt,nvar,igr,1,0                                          &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'FF_NHGT :90:hist:anal:dail:mont:dcyc:year')

      nvar=nvar+1
      call vtable_edio_i_sca(ndcycle,nvar,igr,0,0                                              &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'NDCYCLE :90:hist:anal:dail:mont:dcyc:year')
      call vtable_edio_i_sca(ndcycle,nvar,igr,1,0                                              &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'NDCYCLE :90:hist:anal:dail:mont:dcyc:year')

      nvar=nvar+1
      call vtable_edio_i_sca(isoilflg(igr),nvar,igr,0,0                                    &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'ISOILFLG :90:hist:anal:dail:mont:dcyc:year')
      call vtable_edio_i_sca(isoilflg(igr),nvar,igr,1,0                                    &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'ISOILFLG :90:hist:anal:dail:mont:dcyc:year')

      nvar=nvar+1
      call vtable_edio_r_sca(slxsand,nvar,igr,0,0                                          &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'SLXSAND :90:hist:anal:dail:mont:dcyc:year')
      call vtable_edio_r_sca(slxsand,nvar,igr,1,0                                          &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'SLXSAND :90:hist:anal:dail:mont:dcyc:year')

      nvar=nvar+1
      call vtable_edio_r_sca(slxclay,nvar,igr,0,0                                          &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'SLXCLAY :90:hist:anal:dail:mont:dcyc:year')
      call vtable_edio_r_sca(slxclay,nvar,igr,1,0                                          &
                            ,var_len,var_len_global,max_ptrs                               &
                            ,'SLXCLAY :90:hist:anal:dail:mont:dcyc:year')


      !------------------------------------------------------------------------------------!
      !    1-D variables, soil layers.                                                     !
      !------------------------------------------------------------------------------------!
      var_len        = nzg
      var_len_global = nzg



      nvar=nvar+1
      call vtable_edio_r(nzg,slz,nvar,igr,0,0                                              &
                        ,var_len,var_len_global,max_ptrs                                   &
                        ,'SLZ :92:hist:anal:dail:mont:dcyc:year')
      call vtable_edio_r(nzg,slz,nvar,igr,1,0                                              &
                        ,var_len,var_len_global,max_ptrs                                   &
                        ,'SLZ :92:hist:anal:dail:mont:dcyc:year')
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    1-D variables, height classes.                                                  !
      !------------------------------------------------------------------------------------!
      var_len        = ff_nhgt
      var_len_global = ff_nhgt



      nvar=nvar+1
      call vtable_edio_r(ff_nhgt,hgt_class,nvar,igr,0,0                                    &
                        ,var_len,var_len_global,max_ptrs                                   &
                        ,'HGT_CLASS :96:hist:anal:dail:mont:dcyc:year')
      call vtable_edio_r(ff_nhgt,hgt_class,nvar,igr,1,0                                    &
                        ,var_len,var_len_global,max_ptrs                                   &
                        ,'HGT_CLASS :96:hist:anal:dail:mont:dcyc:year')
      !------------------------------------------------------------------------------------!


   


      !----- Save the number of global-level variables that go to the output. -------------!
      nioglobal=nvar
      !------------------------------------------------------------------------------------!

      return
   end subroutine filltab_globtype
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine will fill the pointer table with the polygon-level variables         !
   ! (edtype).                                                                             !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_edtype(igr,init)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: init
      integer     , intent(in) :: igr
      !----- Local variables. -------------------------------------------------------------!
      type(edtype), pointer    :: cgrid
      integer                  :: var_len
      integer                  :: max_ptrs
      integer                  :: var_len_global
      integer                  :: nvar
      integer                  :: npts
      !------------------------------------------------------------------------------------!

      cgrid => edgrid_g(igr)

      !------ Define the global dimensions. -----------------------------------------------!
      var_len        = cgrid%npolygons
      var_len_global = cgrid%npolygons_global
      max_ptrs       = 1

      !------ Continue the counting. ------------------------------------------------------!
      nvar = nioglobal

      call filltab_edtype_p10 (cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      call filltab_edtype_p11 (cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      call filltab_edtype_m11 (cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      call filltab_edtype_p120(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      call filltab_edtype_p12 (cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      call filltab_edtype_m12 (cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      call filltab_edtype_p14 (cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      call filltab_edtype_p16 (cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      call filltab_edtype_p19 (cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      call filltab_edtype_p146(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      call filltab_edtype_p199(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      call filltab_edtype_p155(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      !----- Save the number of polygon-level (edtype) variables that go to the output. ---!
      if (init == 0) niogrid=nvar-nioglobal
      !------------------------------------------------------------------------------------!

      return
     
   end subroutine filltab_edtype
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine will fill the pointer table with the polygon-level variables         !
   ! (edtype) that have one dimension and are integer (type 10).                           !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_edtype_p10(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      use ed_var_tables, only : vtable_edio_i  & ! sub-routine
                              , metadata_edio  ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype), target        :: cgrid
      integer     , intent(in)    :: init
      integer     , intent(in)    :: igr
      integer     , intent(in)    :: var_len
      integer     , intent(in)    :: max_ptrs
      integer     , intent(in)    :: var_len_global
      integer     , intent(inout) :: nvar
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: npts
      !------------------------------------------------------------------------------------!
      
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 1-D block.  All variables must have the number of points defined   !
      ! by npts.                                                                           !
      !------------------------------------------------------------------------------------!

      npts = cgrid%npolygons

      if (associated(cgrid%pysi_id)) then
         nvar = nvar + 1
         call vtable_edio_i(npts,cgrid%pysi_id,nvar,igr,init,cgrid%pyglob_id               &
                           ,var_len,var_len_global,max_ptrs                                &
                           ,'PYSI_ID :10:hist:anal:dail:mont:dcyc:year')
         call metadata_edio(nvar,igr,'Polygons first site indices','NA','ipoly')
      end if

      
      if (associated(cgrid%pysi_n)) then
         nvar=nvar+1
         call vtable_edio_i(npts,cgrid%pysi_n,nvar,igr,init,cgrid%pyglob_id                &
                           ,var_len,var_len_global,max_ptrs                                &
                           ,'PYSI_N :10:hist:anal:dail:mont:dcyc:year')
         call metadata_edio(nvar,igr,'Number of sites per polygon','NA','ipoly')
      end if

      if (associated(cgrid%xatm)) then
         nvar=nvar+1
         call vtable_edio_i(npts,cgrid%xatm,nvar,igr,init,cgrid%pyglob_id                  &
                           ,var_len,var_len_global,max_ptrs                                &
                           ,'XATM :10:hist:anal:dail:mont:dcyc:year') 
         call metadata_edio(nvar,igr,'Atm. cell x-indices of polygon','NA','ipoly')
      end if
      
      if (associated(cgrid%yatm)) then
         nvar=nvar+1
         call vtable_edio_i(npts,cgrid%yatm,nvar,igr,init,cgrid%pyglob_id                  &
                           ,var_len,var_len_global,max_ptrs                                &
                           ,'YATM :10:hist:anal:dail:mont:dcyc:year') 
         call metadata_edio(nvar,igr,'Atm cell y-indices of polygon','NA','ipoly')
      end if
      
      if (associated(cgrid%lsl)) then
         nvar=nvar+1
         call vtable_edio_i(npts,cgrid%lsl,nvar,igr,init,cgrid%pyglob_id                   &
                           ,var_len,var_len_global,max_ptrs                                &
                           ,'LSL :10:hist:anal:dail:mont:dcyc:year')
         call metadata_edio(nvar,igr,'Index of lowest soil layer','NA','ipoly')
         
      end if
      
      if (associated(cgrid%load_adjacency)) then
         nvar=nvar+1
         call vtable_edio_i(npts,cgrid%load_adjacency,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'LOAD_ADJACENCY :10:hist') 
         call metadata_edio(nvar,igr,'Load Adjacency','[NA]','ipoly')
      end if

      return     
   end subroutine filltab_edtype_p10
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine will fill the pointer table with the polygon-level variables         !
   ! (edtype) that have one dimension and are real (type 11).                              !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_edtype_p11(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      use ed_var_tables, only : vtable_edio_r  & ! sub-routine
                              , metadata_edio  ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype), target        :: cgrid
      integer     , intent(in)    :: init
      integer     , intent(in)    :: igr
      integer     , intent(in)    :: var_len
      integer     , intent(in)    :: max_ptrs
      integer     , intent(in)    :: var_len_global
      integer     , intent(inout) :: nvar
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: npts
      !------------------------------------------------------------------------------------!


      
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 1-D block.  All variables must have the number of points defined   !
      ! by npts.                                                                           !
      !------------------------------------------------------------------------------------!
      npts = cgrid%npolygons




      if (associated(cgrid%lat)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%lat,nvar,igr,init,cgrid%pyglob_id                   &
                           ,var_len,var_len_global,max_ptrs                                &
                           ,'LATITUDE :11:hist:anal:dail:mont:dcyc:year')
         call metadata_edio(nvar,igr,'Latitude of Polygon','decimal degrees','ipoly')
      end if


      if (associated(cgrid%lon)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%lon,nvar,igr,init,cgrid%pyglob_id                   &
                           ,var_len,var_len_global,max_ptrs                                &
                           ,'LONGITUDE :11:hist:anal:dail:mont:dcyc:year') 
         call metadata_edio(nvar,igr,'Longitude of Polygon','decimal degrees','ipoly')
      end if
      
      if (associated(cgrid%wbar)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%wbar,nvar,igr,init,cgrid%pyglob_id                  &
                           ,var_len,var_len_global,max_ptrs                                &
                           ,'WBAR :11:hist') 
         call metadata_edio(nvar,igr,'Polygon average topographic moisture index'          &
                           ,'NA','ipoly')
      end if
      
      if (associated(cgrid%Te)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Te,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'TE :11:hist') 
         call metadata_edio(nvar,igr,'NA','NA','ipoly')
      end if
      
      if (associated(cgrid%zbar)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%zbar,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'ZBAR :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon average water table depth','[m]','ipoly')
      end if
      
      if (associated(cgrid%sheat)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%sheat,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'SHEAT :11:hist') 
         call metadata_edio(nvar,igr,'soil heat pool for lateral hydrology','NA','ipoly')
      end if
      
      if (associated(cgrid%baseflow)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%baseflow,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'BASEFLOW :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'loss of water from site to watershed discharge','kg/m2/s','ipoly')
      end if
      
      if (associated(cgrid%runoff)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%runoff,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'RUNOFF :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'NA','NA','ipoly')
      end if
      
      if (associated(cgrid%swliq)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%swliq,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'SWLIQ :11:hist:anal') 
         call metadata_edio(nvar,igr,'NA','NA','ipoly')
      end if

      if (associated(cgrid%total_agb)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%total_agb,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'TOTAL_AGB :11:hist:anal:year') 
         call metadata_edio(nvar,igr,'Polygon Total Above Ground Biomass','[kgC/m2]','ipoly')
      end if
      
      if (associated(cgrid%total_basal_area)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%total_basal_area,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'TOTAL_BASAL_AREA :11:hist:anal:year') 
         call metadata_edio(nvar,igr,'Polygon Total Basal Area','[cm2/m2]','ipoly')
         
      end if

      if (associated(cgrid%total_agb_growth)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%total_agb_growth,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'TOTAL_AGB_GROWTH:11:hist:anal:year') 
          call metadata_edio(nvar,igr,'Polygon AGB gain through growth','[kgC/m2/yr]','ipoly')
      end if
      
      if (associated(cgrid%total_agb_mort)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%total_agb_mort,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'TOTAL_AGB_MORT :11:hist:anal:year') 
         call metadata_edio(nvar,igr,'Polygon AGB lost due to mortality','[kgC/m2/yr]','ipoly')
      end if
      
      if (associated(cgrid%total_agb_recruit)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%total_agb_recruit,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'TOTAL_AGB_RECRUIT :11:hist:anal:year') 
         call metadata_edio(nvar,igr,'Polygon AGB used to generate recruits','[kgC/m2/yr]','ipoly')
      end if
      
      if (associated(cgrid%total_basal_area_growth)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%total_basal_area_growth,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'TOTAL_BASAL_AREA_GROWTH :11:hist:anal:year') 
         call metadata_edio(nvar,igr,'Polygon basal area gained through growth ','[kgC/m2/yr]','ipoly')
      end if
      
      if (associated(cgrid%total_basal_area_mort)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%total_basal_area_mort,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'TOTAL_BASAL_AREA_MORT :11:hist:anal:year') 
         call metadata_edio(nvar,igr,'Polygon basal area lost through growth ','[cm2/m2/yr]','ipoly')
      end if
      
      if (associated(cgrid%total_basal_area_recruit)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%total_basal_area_recruit,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'TOTAL_BASAL_AREA_RECRUIT :11:hist:anal:year') 
         call metadata_edio(nvar,igr,'Polygon basal area gained by recruits','[cm2/m2/yr]','ipoly')
      end if
      
      if (associated(cgrid%cosz)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%cosz,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'COSZ :11:hist') 
         call metadata_edio(nvar,igr,'Cosine of the zenith angle','[a/h]','ipoly')
      end if
      
      if (associated(cgrid%cbudget_initialstorage)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%cbudget_initialstorage,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'CBUDGET_INITIALSTORAGE :11:hist') 
         call metadata_edio(nvar,igr,'Vegetation and soil carbon,at start of budget-averaging','[kgC/m2]','ipoly')
      end if
          
      if (associated(cgrid%cbudget_nep)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%cbudget_nep,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'CBUDGET_NEP :11:hist') 
         call metadata_edio(nvar,igr,'Polygon average net ecosystem production','[kgC/m2/day]','ipoly')
      end if
      
      if (associated(cgrid%nbudget_initialstorage)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%nbudget_initialstorage,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NBUDGET_INITIALSTORAGE :11:hist') 
         call metadata_edio(nvar,igr,'Veg and soil nitrogen, at start of budget-averaging','[kgN/m2]','ipoly')
      end if

      if (associated(cgrid%avg_vapor_lc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_vapor_lc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_VAPOR_LC :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'polygon leaf to canopy air vapor flux','[kg/m2/s]','ipoly') 
      end if

      if (associated(cgrid%avg_vapor_wc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_vapor_wc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_VAPOR_WC :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'polygon wood to canopy air vapor flux','[kg/m2/s]','ipoly') 
      end if
      
      if (associated(cgrid%avg_dew_cg)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_dew_cg,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_DEW_CG :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon averaged dew to ground flux','[kg/m2/s]','ipoly') 
      end if
      
      if (associated(cgrid%avg_vapor_gc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_vapor_gc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_VAPOR_GC :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'polygon moisture flux ground to canopy air','[kg/m2/s]','ipoly') 
      end if
      
      if (associated(cgrid%avg_wshed_vg)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_wshed_vg,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_WSHED_VG :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon averaged water shed from vegetation to ground','[kg/m2/s]','ipoly') 
      end if
      
      if (associated(cgrid%avg_intercepted)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_intercepted,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_INTERCEPTED :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon averaged intercepted precipitation by vegetation','[kg/m2/s]','ipoly') 
      end if
      
      if (associated(cgrid%avg_throughfall)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_throughfall,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_THROUGHFALL :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon averaged throughfall precipitation','[kg/m2/s]','ipoly') 
      end if
      
      if (associated(cgrid%avg_vapor_ac)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_vapor_ac,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_VAPOR_AC :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'polygon vapor flux atmosphere to canopy air','[kg/m2/s]','ipoly') 
      end if
      
      if (associated(cgrid%avg_transp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_transp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_TRANSP :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'polygon transpiration from stomata to canopy air space','[kg/m2/s]','ipoly') 
      end if
      
      if (associated(cgrid%avg_evap)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_evap,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_EVAP :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon averaged evap/dew from ground and leaves to CAS','[kg/m2/s]','ipoly') 
      end if

      if (associated(cgrid%avg_runoff)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_runoff,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_RUNOFF :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon average surface runoff','[kg/m2/s]','NA') 
      end if
      
      if (associated(cgrid%avg_drainage)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_drainage,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_DRAINAGE :11:hist:anal') 
         call metadata_edio(nvar,igr,'polygon average water flux through lower soil layer','[kg/m2/s]','ipoly') 
      end if
       
      if (associated(cgrid%avg_drainage_heat)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_drainage_heat,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_DRAINAGE_HEAT :11:hist:anal') 
         call metadata_edio(nvar,igr,'polygon average internal energy loss through lower soil layer','[W/m2]','ipoly') 
      end if
     
      if (associated(cgrid%aux)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%aux,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AUX :11:hist:anal') 
         call metadata_edio(nvar,igr,'Auxillary variable - user discretion,see rk4_derivs.f90','[user-defined]','ipoly') 
      end if
      

      if (associated(cgrid%avg_rshort_gnd)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_rshort_gnd,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_RSHORT_GND :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon averaged ground absorbed SW radiation','[W/m2]','ipoly') 
      end if

      if (associated(cgrid%avg_rlong_gnd)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_rlong_gnd,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_RLONG_GND :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon averaged ground absorbed LW radiation','[W/m2]','ipoly') 
      end if

      if (associated(cgrid%avg_carbon_ac)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_carbon_ac,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_CARBON_AC :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon averaged vegetation to canopy air  CO2 flux','[umol/m2/s]','ipoly') 
      end if

      if (associated(cgrid%avg_sensible_lc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_sensible_lc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SENSIBLE_LC :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon averaged leaf to canopy air sensible heat flux','[W/m2]','ipoly') 
      end if
 
      if (associated(cgrid%avg_sensible_wc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_sensible_wc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SENSIBLE_WC :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon averaged wood to canopy air sensible heat flux','[W/m2]','ipoly') 
      end if
     
      if (associated(cgrid%avg_qwshed_vg)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_qwshed_vg,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_QWSHED_VG :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon averaged internal energy flux of water shed from vegetation to ground','[W/m2]','ipoly') 
      end if
       
      if (associated(cgrid%avg_qintercepted)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_qintercepted,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_QINTERCEPTED :11:hist:anal') 
         call metadata_edio(nvar,igr,&
         'Polygon averaged internal energy flux of intercepted precipitation by vegetation','[W/m2]','ipoly') 
      end if
       
      if (associated(cgrid%avg_qthroughfall)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_qthroughfall,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_QTHROUGHFALL :11:hist:anal') 
         call metadata_edio(nvar,igr,&
             'Polygon averaged internal energy flux of throughfall precipitation','[W/m2]','ipoly') 
      end if
     
      if (associated(cgrid%avg_sensible_gc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_sensible_gc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SENSIBLE_GC :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon averaged sensible heat flux ground to canopy air','[W/m2]','ipoly') 
      end if
      
      if (associated(cgrid%avg_sensible_ac)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_sensible_ac,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SENSIBLE_AC :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon averaged sensible heat flux atmosphere  to canopy','[W/m2]','ipoly') 
      end if

      if (associated(cgrid%avg_runoff_heat)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_runoff_heat,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_RUNOFF_HEAT :11:hist:anal') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      ! ---------------------------------------------
      
      if (associated(cgrid%avg_gpp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_gpp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_GPP :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon Average GPP','[umol/m2/s]','ipoly') 
      end if

!      if (associated(cgrid%avg_nppleaf)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%avg_nppleaf,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'AVG_NPPLEAF :11:hist:anal:opti') 
!         call metadata_edio(nvar,igr,'Polygon Average NPP leaf','[kgC/m2/day]','ipoly') 
!      end if

!      if (associated(cgrid%avg_nppfroot)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%avg_nppfroot,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'AVG_NPPFROOT :11:hist:anal:opti') 
!         call metadata_edio(nvar,igr,'Polygon Average NPP froot','[kgC/m2/day]','ipoly') 
!      end if
            
!      if (associated(cgrid%avg_nppsapwood)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%avg_nppsapwood,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'AVG_NPPSAPWOOD :11:hist:anal:opti') 
!         call metadata_edio(nvar,igr,'Polygon Average NPP sapwood','[kgC/m2/day]','ipoly') 
!      end if
      
!      if (associated(cgrid%avg_nppcroot)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%avg_nppcroot,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'AVG_NPPCROOT :11:hist:anal:opti') 
!         call metadata_edio(nvar,igr,'Polygon Average NPP croot','[kgC/m2/day]','ipoly') 
!      end if

!      if (associated(cgrid%avg_nppseeds)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%avg_nppseeds,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'AVG_NPPSEEDS :11:hist:anal:opti') 
!         call metadata_edio(nvar,igr,'Polygon Average NPP seeds','[kgC/m2/day]','ipoly') 
!      end if


!      if (associated(cgrid%avg_nppwood)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%avg_nppwood,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'AVG_NPPWOOD :11:hist:anal:opti') 
!         call metadata_edio(nvar,igr,'Polygon Average NPP wood','[kgC/m2/day]','ipoly') 
!      end if
      
!      if (associated(cgrid%avg_nppdaily)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%avg_nppdaily,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'AVG_NPPDAILY :11:hist:anal:opti') 
!         call metadata_edio(nvar,igr,'Polygon Average NPP daily','[kgC/m2/day]','ipoly') 
!      end if

      if (associated(cgrid%lai)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%lai,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'LAI :11:hist:anal:dail') 
         call metadata_edio(nvar,igr,'Polygon  LAI','[m2/m2]','ipoly') 
      end if

      if (associated(cgrid%wpa)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%wpa,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'WPA :11:hist:anal:dail') 
         call metadata_edio(nvar,igr,'Polygon wood projected area','[m2/m2]','ipoly') 
      end if

      if (associated(cgrid%avg_lma)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_lma,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_LMA :11:hist:dail') 
         call metadata_edio(nvar,igr,'Polygon LMA','[NA]','ipoly') 
      end if

      if (associated(cgrid%wai)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%wai,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'WAI :11:hist:anal:dail') 
         call metadata_edio(nvar,igr,'Polygon  wood area index','[m2/m2]','ipoly') 
      end if

      if (associated(cgrid%avg_leaf_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_leaf_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_LEAF_RESP :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon Average Leaf Respiration','[umol/m2/s]','ipoly') 
      end if

      if (associated(cgrid%avg_root_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_root_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_ROOT_RESP :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon Average Root Respiration','[umol/m2/s]','ipoly') 
      end if

      if (associated(cgrid%avg_growth_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_growth_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_GROWTH_RESP :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Growth Respiration','[umol/m2/s]','ipoly') 
      end if

      if (associated(cgrid%avg_storage_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_storage_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_STORAGE_RESP :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Storage Respiration','[umol/m2/s]','ipoly') 
      end if

      if (associated(cgrid%avg_vleaf_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_vleaf_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_VLEAF_RESP :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Virtual Leaf Respiration','[umol/m2/s]','ipoly') 
      end if

      if (associated(cgrid%avg_plant_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_plant_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_PLANT_RESP :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon Average Plant Respiration','[umol/m2/s]','ipoly') 
      end if

      if (associated(cgrid%avg_growth_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_growth_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_GROWTH_RESP :11:opti') 
         call metadata_edio(nvar,igr,'Polygon Average Growth Respiration','[umol/m2/s]','ipoly') 
      end if

      if (associated(cgrid%avg_storage_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_storage_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_STORAGE_RESP :11:opti') 
         call metadata_edio(nvar,igr,'Polygon Average Storage Respiration','[umol/m2/s]','ipoly') 
      end if

      if (associated(cgrid%avg_vleaf_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_vleaf_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_VLEAF_RESP :11:opti') 
         call metadata_edio(nvar,igr,'Polygon Average VLeaf Respiration','[umol/m2/s]','ipoly') 
      end if

      if (associated(cgrid%avg_htroph_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_htroph_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_HTROPH_RESP :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon Average Heterotrohic Respiration','[umol/m2/s]','ipoly') 
      end if

      if (associated(cgrid%avg_leaf_drop)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_leaf_drop,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_LEAF_DROP :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Leaf loss','[kgC/m2/yr]','ipoly') 
      end if

      if (associated(cgrid%avg_leaf_maintenance)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_leaf_maintenance,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_LEAF_MAINTENANCE :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Leaf maintenance cost','[kgC/m2/day]','ipoly') 
      end if

      if (associated(cgrid%avg_root_maintenance)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_root_maintenance,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_ROOT_MAINTENANCE :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average fine root maintenance cost','[kgC/m2/yr]','ipoly') 
      end if


         !!! added for NACP intercomparison (MCD)
      if (associated(cgrid%avg_sfcw_depth)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_sfcw_depth,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SNOWDEPTH :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Poly Avg. Snow Depth ','[m]','ipoly') 
      end if
      if (associated(cgrid%avg_sfcw_energy)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_sfcw_energy,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SNOWENERGY :11:hist:anal') 
         call metadata_edio(nvar,igr,'Poly Avg. Snow Energy ','[J/kg]','ipoly') 
      end if
      if (associated(cgrid%avg_sfcw_mass)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_sfcw_mass,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SNOWMASS :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Poly Avg. Snow Mass (SWE) ','[kg/m2]','ipoly') 
      end if
      if (associated(cgrid%avg_sfcw_tempk)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_sfcw_tempk,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SNOWTEMP :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Poly Avg. Snow Temperature','[K]','ipoly') 
      end if
      if (associated(cgrid%avg_sfcw_fracliq)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_sfcw_fracliq,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SNOWFRACLIQ :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Poly Avg. Snow liquid fraction','[proportion]','ipoly') 
      end if
      if (associated(cgrid%avg_bdead)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_bdead,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_BDEAD :11:hist:opti:dail:anal:year') 
         call metadata_edio(nvar,igr,'Poly Avg. Biomass - structural','[kgC/m2]','ipoly') 
      end if
      if (associated(cgrid%avg_bstorage)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_bstorage,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_BSTORAGE :11:hist:anal:year') 
         call metadata_edio(nvar,igr,'Poly Avg. Biomass - storage','[kgC/m2]','ipoly') 
      end if
      if (associated(cgrid%avg_bseeds)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_bseeds,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_BSEEDS :11:hist:anal:year') 
         call metadata_edio(nvar,igr,'Poly Avg. Biomass - seeds','[kgC/m2]','ipoly') 
      end if

      if (associated(cgrid%avg_balive)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_balive,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_BALIVE :11:hist:anal:year') 
         call metadata_edio(nvar,igr,'Poly Avg. Biomass -- living','[kgC/m2]','ipoly') 
      end if

      if (associated(cgrid%avg_bleaf)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_bleaf,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_BLEAF :11:hist:anal:year') 
         call metadata_edio(nvar,igr,'Poly Avg. Biomass -- leaf','[kgC/m2]','ipoly') 
      end if

      if (associated(cgrid%avg_broot)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_broot,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_BROOT :11:hist:anal:year') 
         call metadata_edio(nvar,igr,'Poly Avg. Biomass -- fine roots','[kgC/m2]','ipoly') 
      end if

      if (associated(cgrid%avg_bsapwood)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_bsapwood,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_BSAPWOOD :11:hist:anal:year') 
         call metadata_edio(nvar,igr,'Poly Avg. Biomass -- sapwood','[kgC/m2]','ipoly') 
      end if

      if (associated(cgrid%avg_fsc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_fsc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_FSC :11:hist:anal:opti:dail:year') 
         call metadata_edio(nvar,igr,'Poly Avg. Fast Soil Carbon','[kg/m2]','ipoly') 
      end if
      if (associated(cgrid%avg_ssc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_ssc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SSC :11:hist:anal:opti:dail:year') 
         call metadata_edio(nvar,igr,'Poly Avg. Slow Soil Carbon','[kg/m2]','ipoly') 
      end if
      if (associated(cgrid%avg_stsc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_stsc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_STSC :11:hist:anal:opti:year:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Structural Soil Carbon','[kg/m2]','ipoly') 
      end if


      if (associated(cgrid%avg_fsn)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_fsn,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_FSN :11:hist:anal:opti:dail:year') 
         call metadata_edio(nvar,igr,'Poly Avg. Fast Soil Nitrogen','[kg/m2]','ipoly') 
      end if
      if (associated(cgrid%avg_msn)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_msn,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_MSN :11:hist:anal:opti:dail:year') 
         call metadata_edio(nvar,igr,'Poly Avg. Mineralized Soil Carbon','[kg/m2]','ipoly') 
      end if




       !-------- TOTAL CARBON AND NITROGEN POOLS  ---------------
       ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
      if (associated(cgrid%Cleaf)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Cleaf,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'CLEAF :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Leaf Carbon','?','ipoly') 
      end if
      if (associated(cgrid%Croot)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Croot,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'CROOT :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Fine Root Carbon','?','ipoly') 
      end if
      if (associated(cgrid%Cstore)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%cstore,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'CSTORE :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Storage/TNC Carbon','?','ipoly') 
      end if
      if (associated(cgrid%Ccwd)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Ccwd,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'CCWD :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Coarse Woody Debris Carbon','?','ipoly') 
      end if
      if (associated(cgrid%Nleaf)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Nleaf,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NLEAF :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Leaf Nitrogen','?','ipoly') 
      end if
      if (associated(cgrid%Ndead)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Ndead,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NDEAD :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Wood Nitrogen','?','ipoly') 
      end if
      if (associated(cgrid%Nroot)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Nroot,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NROOT :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Fine Root Nitrogen','?','ipoly') 
      end if
      if (associated(cgrid%Nstore)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Nstore,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NSTORE :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Storage Nitrogen','?','ipoly') 
      end if
      if (associated(cgrid%Ncwd)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Ncwd,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NCWD :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Leaf Carbon','?','ipoly') 
      end if

       !-------- TOTAL CARBON AND NITROGEN FLUX  ---------------
       ! Added by MCD for NCEAS/FACE intercomparison (Apr 7 2009)
      if (associated(cgrid%Cleaf_grow)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Cleaf_grow,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'CLEAF_GROW :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Leaf Carbon Growth','?','ipoly') 
      end if
      if (associated(cgrid%Croot_grow)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Croot_grow,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'CROOT_GROW :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Fine Root Carbon Growth','?','ipoly') 
      end if
      if (associated(cgrid%Cdead_grow)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Cdead_grow,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'CDEAD_GROW :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Wood Carbon Growth','?','ipoly') 
      end if
      if (associated(cgrid%Cstore_grow)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Cstore_Grow,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'CSTORE_GROW :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Storage Carbon growth','?','ipoly') 
      end if
      if (associated(cgrid%Cleaf_litter_flux)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Cleaf_litter_flux,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'CLEAF_LITTER_FLUX :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Leaf Litter Carbon Flux','?','ipoly') 
      end if
      if (associated(cgrid%Croot_litter_flux)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Croot_litter_flux,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'CROOT_LITTER_FLUX :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Fine Root Litter Carbon Flux','?','ipoly') 
      end if
      if (associated(cgrid%Ccwd_flux)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Ccwd_flux,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'CCWD_FLUX :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Coarse Woody Debris Carbon FLUX','?','ipoly') 
      end if
      if (associated(cgrid%Nleaf_grow)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Nleaf_grow,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NLEAF_GROW :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Leaf Nitrogen growth','?','ipoly') 
      end if
      if (associated(cgrid%Nroot_grow)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Nroot_grow,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NROOT_GROW :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Fine Root Nitrogen Growth','?','ipoly') 
      end if
      if (associated(cgrid%Ndead_grow)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Ndead_grow,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NDEAD_GROW :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Wood Nitrogen growth','?','ipoly') 
      end if
      if (associated(cgrid%Nstore_grow)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Nstore_grow,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NSTORE_GROW :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Storage Nitrogen Growth','?','ipoly') 
      end if
      if (associated(cgrid%Nleaf_litter_flux)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Nleaf_litter_flux,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NLEAF_LITTER_FLUX :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Leaf litter Nitrogen flux','?','ipoly') 
      end if
      if (associated(cgrid%Nroot_litter_flux)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Nroot_litter_flux,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NROOT_LITTER_FLUX :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. fine root litter Nitrogen flux','?','ipoly') 
      end if
      if (associated(cgrid%Ncwd_flux)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Ncwd_flux,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NCWD_FLUX :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. CWD Nitrogen flux','?','ipoly') 
      end if
      if (associated(cgrid%Nbiomass_uptake)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Nbiomass_uptake,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NBIOMASS_UPTAKE :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Nitrogen Uptake','?','ipoly') 
      end if
      if (associated(cgrid%Ngross_min)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Ngross_min,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NGROSS_MIN :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Gross Nitrogen mineralization','?','ipoly') 
      end if
      if (associated(cgrid%Nnet_min)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%Nnet_min,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'NNET_MIN :11:dail') 
         call metadata_edio(nvar,igr,'Poly Avg. Net Nitrogen mineralization','?','ipoly') 
      end if

      ! ----------------------------------------------
      
      if (associated(cgrid%avg_nir_beam)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_nir_beam,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_NIR_BEAM:11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Averaged Incident Near Infrared Beam Radiation','[W/m2]','ipoly') 
      end if
      
      if (associated(cgrid%avg_nir_diffuse)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_nir_diffuse,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_NIR_DIFFUSE :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Averaged Incident Near Infrared Diffuse Radiation','[W/m2]','ipoly') 
      end if
      
      if (associated(cgrid%avg_par_beam)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_par_beam,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_PAR_BEAM :11:hist:opti:anal') 
         call metadata_edio(nvar,igr,'Polygon Averaged Incident Beam Photosynthetically Active Radiation','[W/m2]','ipoly') 
      end if
      
      if (associated(cgrid%avg_par_diffuse)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_par_diffuse,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_PAR_DIFFUSE :11:hist:opti:anal') 
         call metadata_edio(nvar,igr,'Polygon Averaged Incident Diffuse Photosynthetically Active Radiation','[W/m2]','ipoly') 
      end if
      
      if (associated(cgrid%avg_atm_tmp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_atm_tmp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_ATM_TMP :11:hist:opti:anal') 
         call metadata_edio(nvar,igr,'Polygon Averaged Atmospheric Temperature at Reference Height','[K]','ipoly') 
      end if
      
      if (associated(cgrid%avg_atm_shv)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_atm_shv,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_ATM_SHV :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Averaged Atmospheric Specific Humidity at Reference Height','[kg/kg]','ipoly') 
      end if
      
      if (associated(cgrid%avg_rshort)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_rshort,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_RSHORT :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Total Incident Shortwave Radiation','[W/m2]','ipoly') 
      end if
      
      if (associated(cgrid%avg_rshort_diffuse)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_rshort_diffuse,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_RSHORT_DIFFUSE :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Diffuse Incident Shortwave Radiation','[W/m2]','ipoly') 
      end if
      
      if (associated(cgrid%avg_rlong)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_rlong,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_RLONG :11:hist:opti:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Incident Longwave Radiation','[W/m2]','ipoly') 
      end if
      
      if (associated(cgrid%avg_pcpg)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_pcpg,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_PCPG :11:hist:opti:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Total Precipitation Rate','[kg/m2/s]','ipoly') 
      end if
      
      if (associated(cgrid%avg_qpcpg)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_qpcpg,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_QPCPG:11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Precipitation Internal Energy Deposition Rate','[W/m2/s]','ipoly') 
      end if
      
      if (associated(cgrid%avg_dpcpg)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_dpcpg,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_DPCPG :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Total Precipitation Depth Rate ','[mm/m2/s]','ipoly') 
      end if
      
      if (associated(cgrid%avg_vels)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_vels,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_VELS :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Wind Magnitude (with instability correction)','[m/s]','ipoly') 
      end if
      
      if (associated(cgrid%avg_atm_prss)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_atm_prss,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_ATM_PRSS :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Atmospheric Pressure at Ref. Height','[Pa]','ipoly') 
      end if
      
      if (associated(cgrid%avg_exner)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_exner,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_EXNER :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Exner Correction','[????]','ipoly') 
      end if
      
      if (associated(cgrid%avg_geoht)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_geoht,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_GEOHT :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Geopotential of Met. Forcing Refernce Height','[m]','ipoly') 
      end if
      
      if (associated(cgrid%avg_atm_co2)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_atm_co2,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_ATM_CO2 :11:hist:opti:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Atmospheric CO2 Concentration at Ref. Height','[ppm]','ipoly') 
      end if
  
      if (associated(cgrid%avg_albedo)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_albedo,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_ALBEDO :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon Average Surface Albedo','[W/W]','ipoly') 
      end if
      
      if (associated(cgrid%avg_albedo_beam)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_albedo_beam,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_ALBEDO_BEAM :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon Average Surface direct Albedo','[W/W]','ipoly') 
      end if
      
      if (associated(cgrid%avg_albedo_diffuse)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_albedo_diffuse,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_ALBEDO_DIFFUSE :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon Average Surface diffuse Albedo','[W/W]','ipoly') 
      end if
      
      if (associated(cgrid%avg_rlong_albedo)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_rlong_albedo,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_RLONG_ALBEDO :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Albedo for Longwave Radiation','[W/m2]','ipoly') 
      end if
          
      if (associated(cgrid%avg_rlongup)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_rlongup,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_RLONGUP :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Upwelling Longwave Radiation','[W/m2]','ipoly') 
      end if

      if (associated(cgrid%max_leaf_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%max_leaf_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MAX_LEAF_TEMP :11:anal') 
         call metadata_edio(nvar,igr,'Temp of the hottest cohort in the polygon','[K]','ipoly') 
      end if     

      if (associated(cgrid%min_leaf_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%min_leaf_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MIN_LEAF_TEMP :11:anal') 
         call metadata_edio(nvar,igr,'Temp of the coldest cohort in the polygon','[K]','ipoly') 
      end if

      if (associated(cgrid%max_wood_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%max_wood_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MAX_WOOD_TEMP :11:anal') 
         call metadata_edio(nvar,igr,'Temp of the hottest cohort in the polygon','[K]','ipoly') 
      end if     

      if (associated(cgrid%min_wood_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%min_wood_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MIN_WOOD_TEMP :11:anal') 
         call metadata_edio(nvar,igr,'Temp of the coldest cohort in the polygon','[K]','ipoly') 
      end if

      if (associated(cgrid%max_soil_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%max_soil_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MAX_SOIL_TEMP :11:anal') 
         call metadata_edio(nvar,igr,'Temp of the hottest soil layer in the polygon','[K]','ipoly') 
      end if    

      if (associated(cgrid%min_soil_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%min_soil_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MIN_SOIL_TEMP :11:anal') 
         call metadata_edio(nvar,igr,'Temp of the coldest soil layer in the polygon','[K]','ipoly') 
      end if    
      
      if (associated(cgrid%avg_leaf_energy)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_leaf_energy,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_LEAF_ENERGY :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Internal Energy of Vegetation','[J/m2]','ipoly') 
      end if

      if (associated(cgrid%avg_leaf_hcap)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_leaf_hcap,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_LEAF_HCAP :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average of Vegetation heat capacity','[J/m2/K]','ipoly') 
      end if

      if (associated(cgrid%avg_leaf_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_leaf_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_LEAF_TEMP :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Temperature of Vegetation','[K]','ipoly') 
      end if

      if (associated(cgrid%avg_leaf_fliq)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_leaf_fliq,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_LEAF_FLIQ :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Liquid fraction of Vegetation Sfc Water','[--]','ipoly') 
      end if
      
      if (associated(cgrid%avg_leaf_water)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_leaf_water,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_LEAF_WATER :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Resident Leaf Surface Water','[kg/m2]','ipoly') 
      end if
      
      if (associated(cgrid%avg_wood_energy)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_wood_energy,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_WOOD_ENERGY :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Internal Energy of wood','[J/m2]','ipoly') 
      end if

      if (associated(cgrid%avg_wood_hcap)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_wood_hcap,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_WOOD_HCAP :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average of wood heat capacity','[J/m2/K]','ipoly') 
      end if

      if (associated(cgrid%avg_wood_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_wood_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_WOOD_TEMP :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Wood of Vegetation','[K]','ipoly') 
      end if

      if (associated(cgrid%avg_wood_fliq)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_wood_fliq,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_WOOD_FLIQ :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Liquid fraction of Vegetation Sfc Water','[--]','ipoly') 
      end if
      
      if (associated(cgrid%avg_wood_water)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_wood_water,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_WOOD_WATER :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Resident Wood Surface Water','[kg/m2]','ipoly') 
      end if
      
      if (associated(cgrid%avg_can_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_can_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_CAN_TEMP :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon Average Temperature of Canopy Air Space','[K]','ipoly') 
      end if
      
      if (associated(cgrid%avg_can_shv)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_can_shv,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_CAN_SHV :11:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon Average Specific Humidity of Canopy Air','[kg/kg]','NA') 
      end if
      
      if (associated(cgrid%avg_can_co2)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_can_co2,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_CAN_CO2 :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average CO2 mixing ratio of Canopy Air','[umol/mol]','NA') 
      end if
      
      if (associated(cgrid%avg_can_rhos)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_can_rhos,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_CAN_RHOS :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Density of Canopy Air','[kg/m3]','NA') 
      end if
      
      if (associated(cgrid%avg_can_prss)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_can_prss,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_CAN_PRSS :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Canopy Air Pressure','[Pa]','NA') 
      end if
      
      if (associated(cgrid%avg_can_theta)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_can_theta,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_CAN_THETA :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Canopy Air Potential temperature','[K]','NA') 
      end if
      
      if (associated(cgrid%avg_can_theiv)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_can_theiv,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_CAN_THEIV :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Canopy Air ice-vapour equiv. pot. temp.','[K]','NA') 
      end if
      
      if (associated(cgrid%avg_can_depth)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_can_depth,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_CAN_DEPTH :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Canopy height','[m]','NA') 
      end if
      
      if (associated(cgrid%avg_soil_wetness)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_soil_wetness,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SOIL_WETNESS :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Soil Wetness RELATIVE TO WILTING POINT','[m3/m3]','ipoly') !relative to wilting point
      end if
      
      if (associated(cgrid%avg_skin_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_skin_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SKIN_TEMP :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Temperature of all surfaces','[K]','ipoly') 
      end if
      
      if (associated(cgrid%avg_available_water)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_available_water,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_AVAILABLE_WATER :11:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Available water','[K]','ipoly') 
      end if
      
      ! Daily and monthly variables. Note that all these variables need to be stored at the
      ! history file, because the averaging can be resumed...
      
      if(associated(cgrid%dmean_gpp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_gpp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_GPP :11:hist:dail') 
         call metadata_edio(nvar,igr,'Polygon Average Daily Integrated Gross Primary Productivity','[kgC/m2/yr]','ipoly') 
      end if
      
!      if (associated(cgrid%dmean_nppleaf)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%dmean_nppleaf,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'DMEAN_NPPLEAF :11:hist:dail') 
!         call metadata_edio(nvar,igr,'Polygon Average Daily NPP leaf','[kgC/m2/yr]','ipoly') 
!      end if

!      if (associated(cgrid%dmean_nppfroot)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%dmean_nppfroot,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'DMEAN_NPPFROOT :11:hist:dail') 
!         call metadata_edio(nvar,igr,'Polygon Average Daily NPP froot','[kgC/m2/yr]','ipoly') 
!      end if

            
!      if (associated(cgrid%dmean_nppsapwood)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%dmean_nppsapwood,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'DMEAN_NPPSAPWOOD :11:hist:dail') 
!         call metadata_edio(nvar,igr,'Polygon Average Daily NPP sapwood','[kgC/m2/yr]','ipoly') 
!      end if

!      if (associated(cgrid%dmean_nppcroot)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%dmean_nppcroot,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'DMEAN_NPPCROOT :11:hist:dail') 
!         call metadata_edio(nvar,igr,'Polygon Average Daily NPP croot','[kgC/m2/yr]','ipoly') 
!      end if

!      if (associated(cgrid%dmean_nppseeds)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%dmean_nppseeds,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'DMEAN_NPPSEEDS :11:hist:dail') 
!         call metadata_edio(nvar,igr,'Polygon Average Daily NPP seeds','[kgC/m2/yr]','ipoly') 
!      end if

!      if (associated(cgrid%dmean_nppwood)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%dmean_nppwood,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'DMEAN_NPPWOOD :11:hist:dail') 
!         call metadata_edio(nvar,igr,'Polygon Average Daily NPP wood','[kgC/m2/yr]','ipoly') 
!      end if
      
!      if (associated(cgrid%dmean_nppdaily)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%dmean_nppdaily,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'DMEAN_NPPDAILY :11:hist:dail') 
!         call metadata_edio(nvar,igr,'Polygon Average Daily NPP daily','[kgC/m2/yr]','ipoly') 
!      end if

      
      if(associated(cgrid%dmean_nee)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_nee,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_NEE :11:hist:dail') 
         call metadata_edio(nvar,igr,'Polygon Average Daily Integrated Net Ecosystem Exchange','[kgC/m2/yr]','ipoly') 
      end if
      
      
      if(associated(cgrid%dmean_pcpg)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_pcpg,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_PCPG :11:hist:dail') 
         call metadata_edio(nvar,igr,'total daily precipitation','[kg/m2/day]','ipoly') 
      end if

      if(associated(cgrid%dmean_runoff)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_runoff,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_RUNOFF :11:hist:dail') 
         call metadata_edio(nvar,igr,'total daily surface runoff','[kg/m2/day]','ipoly') 
      end if

      if(associated(cgrid%dmean_drainage)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_drainage,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_DRAINAGE :11:hist:dail') 
         call metadata_edio(nvar,igr,'total daily water flux through lower soil layer','[kg/m2/day]','ipoly') 
      end if

      if(associated(cgrid%dmean_vapor_ac)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_vapor_ac,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_VAPOR_AC :11:hist:dail') 
         call metadata_edio(nvar,igr,'total daily vapor flux atm->canopy','[kg/m2/day]','ipoly') 
      end if


      if(associated(cgrid%dmean_vapor_gc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_vapor_gc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_VAPOR_GC :11:hist:dail') 
         call metadata_edio(nvar,igr,'total daily vapor flux ground->canopy','[kg/m2/day]','ipoly') 
      end if


      if(associated(cgrid%dmean_vapor_lc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_vapor_lc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_VAPOR_LC :11:hist:dail') 
         call metadata_edio(nvar,igr,'total daily vapor flux leaf->canopy','[kg/m2/day]','ipoly') 
      end if

      if(associated(cgrid%dmean_vapor_wc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_vapor_wc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_VAPOR_WC :11:hist:dail') 
         call metadata_edio(nvar,igr,'total daily vapor flux wood->canopy','[kg/m2/day]','ipoly') 
      end if

      if(associated(cgrid%dmean_evap)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_evap,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_EVAP :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean (leaf+soil) evaporation','[kg/m2/s]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_transp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_transp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_TRANSP :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_sensible_lc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_sensible_lc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_SENSIBLE_LC :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_sensible_wc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_sensible_wc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_SENSIBLE_WC :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_sensible_gc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_sensible_gc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_SENSIBLE_GC :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_sensible_ac)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_sensible_ac,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_SENSIBLE_AC :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_plresp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_plresp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_PLRESP :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_rh)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_rh,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_RH :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_leaf_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_leaf_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_LEAF_RESP :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_root_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_root_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_ROOT_RESP :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_growth_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_growth_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_GROWTH_RESP :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_storage_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_storage_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_STORAGE_RESP :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_vleaf_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_vleaf_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_VLEAF_RESP :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_nep)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_nep,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_NEP :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%dmean_fs_open)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_fs_open,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_FS_OPEN :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_fsw)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_fsw,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_FSW :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_fsn)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_fsn,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_FSN :11:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_can_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_can_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_CAN_TEMP :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean canopy temperature','[K]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_can_shv)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_can_shv,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_CAN_SHV :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean canopy specific humidity','[kg/kg]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_can_co2)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_can_co2,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_CAN_CO2 :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean canopy CO2 mixing ratio','[umol/mol]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_can_rhos)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_can_rhos,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_CAN_RHOS :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean canopy air density','[kg/m3]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_can_prss)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_can_prss,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_CAN_PRSS :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean canopy air pressure','[Pa]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_can_theta)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_can_theta,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_CAN_THETA :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean canopy air potential temperature','[K]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_can_theiv)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_can_theiv,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_CAN_THEIV :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean canopy air theta_Eiv','[K]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_gnd_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_gnd_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_GND_TEMP :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean of ground temperature','[K]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_gnd_shv)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_gnd_shv,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_GND_SHV :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean of ground specific humidity','[kg/kg]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_leaf_energy)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_leaf_energy,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_LEAF_ENERGY :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean vegetation internal energy','[J/m2]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_leaf_water)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_leaf_water,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_LEAF_WATER :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean vegetation surface water','[kg/m2]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_leaf_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_leaf_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_LEAF_TEMP :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean vegetation temperature','[K]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_leaf_hcap)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_leaf_hcap,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_LEAF_HCAP :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean vegetation heat capacity','[J/m2/K]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_wood_energy)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_wood_energy,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_WOOD_ENERGY :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean wood internal energy','[J/m2]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_wood_water)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_wood_water,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_WOOD_WATER :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean wood surface water','[kg/m2]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_wood_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_wood_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_WOOD_TEMP :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean wood temperature','[K]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_wood_hcap)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_wood_hcap,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_WOOD_HCAP :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean vegetation heat capacity','[J/m2/K]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_atm_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_atm_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_ATM_TEMP :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean air temperature','[K]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_atm_shv)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_atm_shv,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_ATM_SHV :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean air specific humidity','[kg/kg]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_atm_co2)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_atm_co2,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_ATM_CO2 :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean air specific humidity','[kg/kg]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_atm_prss)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_atm_prss,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_ATM_PRSS :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean air pressure','[ Pa]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_atm_vels)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_atm_vels,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_ATM_VELS :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean wind speed','[m/s]','ipoly') 
      end if

      if(associated(cgrid%dmean_rshort)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_rshort,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_RSHORT :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean shortwave radiation','[w/m2]','ipoly') 
      end if

      if(associated(cgrid%dmean_rshort_diff)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_rshort_diff,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_RSHORT_DIFF :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean diffuse shortwave radiation','[w/m2]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_rlong)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_rlong,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_RLONG :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean longwave radiation','[w/m2]','ipoly') 
      end if

      if(associated(cgrid%dmean_rshort_gnd)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_rshort_gnd,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_RSHORT_GND :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean ground abs. shortwave radiation','[w/m2]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_rlong_gnd)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_rlong_gnd,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_RLONG_GND :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean ground abs. longwave radiation','[w/m2]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_albedo)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_albedo,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_ALBEDO :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean albedo','[---]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_albedo_beam)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_albedo_beam,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_ALBEDO_BEAM :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean direct albedo','[---]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_albedo_diffuse)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_albedo_diffuse,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_ALBEDO_DIFFUSE :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean diffuse albedo','[---]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_rlong_albedo)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_rlong_albedo,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_RLONG_ALBEDO :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean longwave albedo','[---]','ipoly') 
      end if
      
      if(associated(cgrid%dmean_rlongup)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_rlongup,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_RLONGUP :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean longwave emission from ground','[---]','ipoly') 
      end if

      if(associated(cgrid%dmean_co2_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_co2_residual,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_CO2_RESIDUAL :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean of residual canopy CO2 ','[umol/m2/s]','ipoly') 
      end if

      if(associated(cgrid%dmean_water_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_water_residual,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_WATER_RESIDUAL :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean of residual water ','[kg/m2/s]','ipoly') 
      end if

      if(associated(cgrid%dmean_energy_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_energy_residual,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_ENERGY_RESIDUAL :11:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean of residual water ','[J/m2/s]','ipoly') 
      end if

      if(associated(cgrid%mmean_gpp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_gpp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_GPP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

!      if (associated(cgrid%mmean_nppleaf)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%mmean_nppleaf,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'MMEAN_NPPLEAF :11:hist:mont:dcyc') 
!         call metadata_edio(nvar,igr,'Polygon Average Monthly NPP leaf','[kgC/m2/yr]','ipoly') 
!      end if

!      if (associated(cgrid%mmean_nppfroot)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%mmean_nppfroot,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'MMEAN_NPPFROOT :11:hist:mont:dcyc') 
!         call metadata_edio(nvar,igr,'Polygon Average Monthly NPP froot','[kgC/m2/yr]','ipoly') 
!      end if
            
!      if (associated(cgrid%mmean_nppsapwood)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%mmean_nppsapwood,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'MMEAN_NPPSAPWOOD :11:hist:mont:dcyc') 
!         call metadata_edio(nvar,igr,'Polygon Average Monthly NPP sapwood','[kgC/m2/yr]','ipoly') 
!      end if

!      if (associated(cgrid%mmean_nppcroot)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%mmean_nppcroot,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'MMEAN_NPPCROOT :11:hist:mont:dcyc') 
!         call metadata_edio(nvar,igr,'Polygon Average Monthly NPP croot','[kgC/m2/yr]','ipoly') 
!      end if

!      if (associated(cgrid%mmean_nppseeds)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%mmean_nppseeds,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'MMEAN_NPPSEEDS :11:hist:mont:dcyc') 
!         call metadata_edio(nvar,igr,'Polygon Average Monthly NPP seeds','[kgC/m2/yr]','ipoly') 
!      end if

!      if (associated(cgrid%mmean_nppwood)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%mmean_nppwood,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'MMEAN_NPPWOOD :11:hist:mont:dcyc') 
!         call metadata_edio(nvar,igr,'Polygon Average Monthly NPP wood','[kgC/m2/yr]','ipoly') 
!      end if
      
!      if (associated(cgrid%mmean_nppdaily)) then
!         nvar=nvar+1
!         call vtable_edio_r(npts,cgrid%mmean_nppdaily,nvar,igr,init,cgrid%pyglob_id, &
!              var_len,var_len_global,max_ptrs,'MMEAN_NPPDAILY :11:hist:mont:dcyc') 
!         call metadata_edio(nvar,igr,'Polygon Average Monthly NPP daily','[kgC/m2/yr]','ipoly') 
!      end if
      
      if(associated(cgrid%mmean_nee)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_nee,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_NEE :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%mmean_evap)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_evap,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_EVAP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly Mean (leaf+soil) Evaporation Rate','[kg/m2/s]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_transp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_transp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_TRANSP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmean_sensible_ac)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_sensible_ac,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_SENSIBLE_AC :11:hist:mont:dcyc')
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
      end if

      if(associated(cgrid%mmean_sensible_gc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_sensible_gc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_SENSIBLE_GC :11:hist:mont:dcyc')
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
      end if

      if(associated(cgrid%mmean_sensible_lc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_sensible_lc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_SENSIBLE_LC :11:hist:mont:dcyc')
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
      end if

      if(associated(cgrid%mmean_sensible_wc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_sensible_wc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_SENSIBLE_WC :11:hist:mont:dcyc')
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
      end if

      if(associated(cgrid%mmean_vapor_ac)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_vapor_ac,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_VAPOR_AC :11:hist:mont:dcyc')
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
      end if

      if(associated(cgrid%mmean_vapor_gc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_vapor_gc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_VAPOR_GC :11:hist:mont:dcyc')
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
      end if

      if(associated(cgrid%mmean_vapor_lc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_vapor_lc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_VAPOR_LC :11:hist:mont:dcyc')
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
      end if

      if(associated(cgrid%mmean_vapor_wc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_vapor_wc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_VAPOR_WC :11:hist:mont:dcyc')
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
      end if
      
      if(associated(cgrid%mmean_nep)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_nep,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_NEP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmean_fs_open)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_fs_open,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_FS_OPEN :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%mmean_fsw)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_fsw,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_FSW :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%mmean_fsn)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_fsn,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_FSN :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmean_plresp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_plresp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_PLRESP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%mmean_rh)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_rh,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_RH :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%mmean_leaf_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_leaf_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_LEAF_RESP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%mmean_root_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_root_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_ROOT_RESP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%mmean_growth_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_growth_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_GROWTH_RESP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%mmean_storage_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_storage_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_STORAGE_RESP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%mmean_vleaf_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_vleaf_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_VLEAF_RESP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmean_can_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_can_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_CAN_TEMP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean canopy temperature','[K]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_can_shv)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_can_shv,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_CAN_SHV :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean canopy specific humidity','[kg/kg]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_can_co2)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_can_co2,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_CAN_CO2 :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean canopy CO2 mixing ratio','[umol/mol]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_can_rhos)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_can_rhos,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_CAN_RHOS :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean canopy air density','[kg/m3]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_can_prss)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_can_prss,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_CAN_PRSS :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean canopy air pressure','[Pa]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_can_theta)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_can_theta,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_CAN_THETA :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean canopy air potential temperature','[K]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_can_theiv)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_can_theiv,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_CAN_THEIV :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean canopy air theta_Eiv','[K]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_gnd_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_gnd_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_GND_TEMP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean ground temperature','[K]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_gnd_shv)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_gnd_shv,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_GND_SHV :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean ground specific humidity','[kg/kg]','ipoly') 
      end if

      if(associated(cgrid%mmean_leaf_energy)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_leaf_energy,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_LEAF_ENERGY :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean leaf internal energy','[J/m2]','ipoly') 
      end if

      if(associated(cgrid%mmean_leaf_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_leaf_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_LEAF_TEMP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean leaf temperature','[K]','ipoly') 
      end if

      if(associated(cgrid%mmean_leaf_water)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_leaf_water,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_LEAF_WATER :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean leaf water ','[kg/m2]','ipoly') 
      end if

      if(associated(cgrid%mmean_leaf_hcap)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_leaf_hcap,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_LEAF_HCAP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean leaf heat capacity','[J/m2/K]','ipoly') 
      end if

      if(associated(cgrid%mmean_wood_energy)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_wood_energy,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_WOOD_ENERGY :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean wood internal energy','[J/m2]','ipoly') 
      end if

      if(associated(cgrid%mmean_wood_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_wood_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_WOOD_TEMP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean wood temperature','[K]','ipoly') 
      end if

      if(associated(cgrid%mmean_wood_water)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_wood_water,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_WOOD_WATER :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean wood water ','[kg/m2]','ipoly') 
      end if

      if(associated(cgrid%mmean_wood_hcap)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_wood_hcap,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_WOOD_HCAP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean wood heat capacity','[J/m2/K]','ipoly') 
      end if

      if(associated(cgrid%mmean_rshort)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_rshort,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_RSHORT :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean downwelling solar radiation','[w/m2]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_rshort_diff)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_rshort_diff,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_RSHORT_DIFF :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean incoming diffuse solar radiation','[w/m2]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_rlong)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_rlong,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_RLONG :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean downwelling longwave radiation','[w/m2]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_rshort_gnd)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_rshort_gnd,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_RSHORT_GND :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean ground abs. solar radiation','[w/m2]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_rlong_gnd)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_rlong_gnd,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_RLONG_GND :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean ground abs. longwave radiation','[w/m2]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_albedo)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_albedo,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_ALBEDO :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean albedo','[---]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_albedo_beam)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_albedo_beam,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_ALBEDO_BEAM :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean direct albedo','[---]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_albedo_diffuse)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_albedo_diffuse,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_ALBEDO_DIFFUSE :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean diffuse albedo','[---]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_rlong_albedo)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_rlong_albedo,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_RLONG_ALBEDO :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean longwave albedo','[---]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_rlongup)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_rlongup,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_RLONGUP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean longwave emission from ground','[---]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_atm_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_atm_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_ATM_TEMP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean air temperature','[K]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_atm_shv)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_atm_shv,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_ATM_SHV :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean air specific humidity','[kg/kg]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_atm_co2)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_atm_co2,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_ATM_CO2 :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean air specific humidity','[kg/kg]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_atm_prss)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_atm_prss,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_ATM_PRSS :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean air pressure','[ Pa]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_atm_vels)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_atm_vels,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_ATM_VELS :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean wind speed','[m/s]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_pcpg)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_pcpg,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_PCPG :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean precipitation rate','[kg/m2/s]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_runoff)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_runoff,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_RUNOFF :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean runoff','[kg/m2/s]','ipoly') 
      end if
      
      if(associated(cgrid%mmean_drainage)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_drainage,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_DRAINAGE :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean drainage','[kg/m2/day]','ipoly') 
      end if

      if(associated(cgrid%mmean_co2_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_co2_residual,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_CO2_RESIDUAL :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean of residual canopy CO2 ','[umol/m2/s]','ipoly') 
      end if

      if(associated(cgrid%mmean_water_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_water_residual,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_WATER_RESIDUAL :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean of residual water ','[kg/m2/s]','ipoly') 
      end if

      if(associated(cgrid%mmean_energy_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_energy_residual,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_ENERGY_RESIDUAL :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean of residual energy ','[J/m2/s]','ipoly') 
      end if



      if(associated(cgrid%mmsqu_gpp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_gpp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_GPP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmsqu_leaf_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_leaf_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_LEAF_RESP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmsqu_root_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_root_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_ROOT_RESP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmsqu_plresp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_plresp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_PLRESP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmsqu_nee)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_nee,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_NEE :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmsqu_rh)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_rh,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_RH :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmsqu_sensible_ac)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_sensible_ac,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_SENSIBLE_AC :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmsqu_sensible_lc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_sensible_lc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_SENSIBLE_LC :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmsqu_sensible_wc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_sensible_wc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_SENSIBLE_WC :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmsqu_sensible_gc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_sensible_gc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_SENSIBLE_GC :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmsqu_evap)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_evap,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_EVAP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmsqu_transp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_transp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_TRANSP :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmsqu_vapor_ac)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_vapor_ac,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_VAPOR_AC :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmsqu_vapor_lc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_vapor_lc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_VAPOR_LC :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmsqu_vapor_wc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_vapor_wc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_VAPOR_WC :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmsqu_vapor_gc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmsqu_vapor_gc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMSQU_VAPOR_GC :11:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!

      return     
   end subroutine filltab_edtype_p11
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine will fill the pointer table with the polygon-level variables         !
   ! (edtype) that have two dimensions (ndcycle,npolygons) and are real (type -11).        !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_edtype_m11(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      use ed_var_tables, only : vtable_edio_r  & ! sub-routine
                              , metadata_edio  ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype), target        :: cgrid
      integer     , intent(in)    :: init
      integer     , intent(in)    :: igr
      integer     , intent(in)    :: var_len
      integer     , intent(in)    :: max_ptrs
      integer     , intent(in)    :: var_len_global
      integer     , intent(inout) :: nvar
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: npts
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 2-D block, with diurnal cycle.  All variables must have the number !
      ! of points defined by npts.                                                         !
      !------------------------------------------------------------------------------------!
      npts = cgrid%npolygons * ndcycle
      
      if(associated(cgrid%qmean_pcpg)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_pcpg,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_PCPG :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%qmean_runoff)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_runoff,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_RUNOFF :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%qmean_drainage)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_drainage,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_DRAINAGE :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%qmean_evap)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_evap,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_EVAP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%qmean_transp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_transp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_TRANSP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%qmean_vapor_ac)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_vapor_ac,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_VAPOR_AC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%qmean_vapor_gc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_vapor_gc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_VAPOR_GC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_vapor_lc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_vapor_lc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_VAPOR_LC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_vapor_wc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_vapor_wc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_VAPOR_WC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_sensible_ac)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_sensible_ac,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_SENSIBLE_AC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_sensible_gc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_sensible_gc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_SENSIBLE_GC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_sensible_lc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_sensible_lc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_SENSIBLE_LC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_sensible_wc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_sensible_wc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_SENSIBLE_WC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_nee)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_nee,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_NEE :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_gpp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_gpp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_GPP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_nep)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_nep,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_NEP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_plresp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_plresp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_PLRESP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_rh)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_rh,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_RH :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_leaf_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_leaf_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_LEAF_RESP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_root_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_root_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_ROOT_RESP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_fs_open)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_fs_open,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_FS_OPEN :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_fsw)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_fsw,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_FSW :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_fsn)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_fsn,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_FSN :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_can_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_can_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_CAN_TEMP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_can_shv)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_can_shv,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_CAN_SHV :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_can_co2)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_can_co2,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_CAN_CO2 :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_can_rhos)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_can_rhos,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_CAN_RHOS :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_can_prss)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_can_prss,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_CAN_PRSS :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_can_theta)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_can_theta,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_CAN_THETA :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_can_theiv)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_can_theiv,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_CAN_THEIV :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_gnd_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_gnd_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_GND_TEMP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_gnd_shv)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_gnd_shv,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_GND_SHV :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_leaf_energy)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_leaf_energy,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_LEAF_ENERGY :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_leaf_water)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_leaf_water,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_LEAF_WATER :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_leaf_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_leaf_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_LEAF_TEMP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_leaf_hcap)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_leaf_hcap,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_LEAF_HCAP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_wood_energy)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_wood_energy,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_WOOD_ENERGY :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_wood_water)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_wood_water,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_WOOD_WATER :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_wood_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_wood_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_WOOD_TEMP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_wood_hcap)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_wood_hcap,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_WOOD_HCAP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_atm_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_atm_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_ATM_TEMP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_rshort)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_rshort,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_RSHORT :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_rshort_diff)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_rshort_diff,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_RSHORT_DIFF :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_rlong)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_rlong,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_RLONG :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_rshort_gnd)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_rshort_gnd,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_RSHORT_GND :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_rlong_gnd)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_rlong_gnd,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_RLONG_GND :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%qmean_albedo)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_albedo,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_ALBEDO :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%qmean_albedo_beam)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_albedo_beam,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_ALBEDO_BEAM :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%qmean_albedo_diffuse)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_albedo_diffuse,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_ALBEDO_DIFFUSE :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%qmean_rlong_albedo)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_rlong_albedo,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_RLONG_ALBEDO :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%qmean_rlongup)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_rlongup,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_RLONGUP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_atm_shv)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_atm_shv,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_ATM_SHV :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_atm_co2)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_atm_co2,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_ATM_CO2 :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_atm_prss)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_atm_prss,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_ATM_PRSS :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_atm_vels)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_atm_vels,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_ATM_VELS :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_gpp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_gpp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_GPP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_leaf_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_leaf_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_LEAF_RESP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_root_resp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_root_resp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_ROOT_RESP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_plresp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_plresp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_PLRESP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_nee)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_nee,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_NEE :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_nep)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_nep,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_NEP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_rh)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_rh,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_RH :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_sensible_ac)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_sensible_ac,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_SENSIBLE_AC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_sensible_lc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_sensible_lc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_SENSIBLE_LC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_sensible_wc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_sensible_wc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_SENSIBLE_WC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_sensible_gc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_sensible_gc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_SENSIBLE_GC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_evap)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_evap,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_EVAP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_transp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_transp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_TRANSP :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_vapor_ac)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_vapor_ac,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_VAPOR_AC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_vapor_lc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_vapor_lc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_VAPOR_LC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_vapor_wc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_vapor_wc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_VAPOR_WC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmsqu_vapor_gc)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmsqu_vapor_gc,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMSQU_VAPOR_GC :-11:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!


      return     
   end subroutine filltab_edtype_m11
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine will fill the pointer table with the polygon-level variables         !
   ! (edtype) that have two dimensions (nzg,npolygons) and are integer (type 120).         !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_edtype_p120(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      use ed_var_tables, only : vtable_edio_i  & ! sub-routine
                              , metadata_edio  ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype), target        :: cgrid
      integer     , intent(in)    :: init
      integer     , intent(in)    :: igr
      integer     , intent(in)    :: var_len
      integer     , intent(in)    :: max_ptrs
      integer     , intent(in)    :: var_len_global
      integer     , intent(inout) :: nvar
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: npts
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 2-D block, with soil layers.  All variables must have the number   !
      ! of points defined by npts.                                                         !
      !------------------------------------------------------------------------------------!
      npts = cgrid%npolygons * nzg

      if (associated(cgrid%ntext_soil)) then
         nvar=nvar+1
         call vtable_edio_i(npts,cgrid%ntext_soil,nvar,igr,init,cgrid%pyglob_id            &
                           ,var_len,var_len_global,max_ptrs                                &
                           ,'NTEXT_SOIL :120:hist:anal:dail:mont:dcyc:year')  
         call metadata_edio(nvar,igr,'Polygon mode soil class','OGE2 Class','ipoly-ngz')

      end if


      return     
   end subroutine filltab_edtype_p120
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine will fill the pointer table with the polygon-level variables         !
   ! (edtype) that have two dimensions (ndcycle,npolygons) and are real (type 12).         !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_edtype_p12(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      use ed_var_tables, only : vtable_edio_r  & ! sub-routine
                              , metadata_edio  ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype), target        :: cgrid
      integer     , intent(in)    :: init
      integer     , intent(in)    :: igr
      integer     , intent(in)    :: var_len
      integer     , intent(in)    :: max_ptrs
      integer     , intent(in)    :: var_len_global
      integer     , intent(inout) :: nvar
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: npts
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 2-D block, with soil layers.  All variables must have the number   !
      ! of points defined by npts.                                                         !
      !------------------------------------------------------------------------------------!
      npts = cgrid%npolygons * nzg

      
      if (associated(cgrid%avg_smoist_gg)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_smoist_gg,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SMOIST_GG :12:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon averaged soil moisture flux,layer nzg is flux with CAS','[kg/m2/s]','ipoly-nzg') 
      end if
      
      if (associated(cgrid%avg_transloss)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_transloss,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_TRANSLOSS :12:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon averaged soil moisture sink to transpiration','[kg/m2/s]','ipoly-nzg') 
      end if

      if (associated(cgrid%aux_s)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%aux_s,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AUX_S :12:hist:anal') 
         call metadata_edio(nvar,igr,'Soil layer discretized, auxillary variable, see rk4_derivs.f90','[user-defined]','ipoly-nzg') 
      end if

      
      if (associated(cgrid%avg_sensible_gg)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_sensible_gg,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SENSIBLE_GG :12:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon averaged soil sensible heat flux,layer nzg is flux with CAS ','[W/m2]','ipoly-nzg') 
      end if

      if (associated(cgrid%avg_soil_energy)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_soil_energy,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SOIL_ENERGY :12:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Volumetric Soil Water','[m/m]','ipoly - nzg') 
      end if
      
      if (associated(cgrid%avg_soil_water)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_soil_water,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SOIL_WATER :12:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon Average Volumetric Soil Water','[m/m]','ipoly - nzg') 
      end if
      
      if (associated(cgrid%avg_soil_mstpot)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_soil_mstpot,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SOIL_MSTPOT :12:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon Average Soil Moisture Potential','[m]','ipoly - nzg') 
      end if
      
      if (associated(cgrid%avg_soil_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_soil_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SOIL_TEMP :12:hist:anal:opti') 
         call metadata_edio(nvar,igr,'Polygon Average Soil Temperature','[K]','ipoly - nzg') 
      end if

      if (associated(cgrid%avg_soil_fracliq)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_soil_fracliq,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_SOIL_FRACLIQ :12:hist:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Soil Fraction Liquid','[proportion]','ipoly - nzg') 
      end if
      

      if (associated(cgrid%avg_soil_rootfrac)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_soil_rootfrac,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'POLY_SOIL_ROOTFRAC :12:hist:mont:year') 
         call metadata_edio(nvar,igr,'Polygon Average Root Fraction','[pdf]','ipoly - nzg') 
      end if


      if(associated(cgrid%dmean_soil_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_soil_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_SOIL_TEMP :12:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_soil_water)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_soil_water,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_SOIL_WATER :12:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_soil_mstpot)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_soil_mstpot,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_SOIL_MSTPOT :12:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%dmean_transloss)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_transloss,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_TRANSLOSS :12:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%mmean_soil_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_soil_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_SOIL_TEMP :12:hist:mont:dcyc')
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
      end if

      if(associated(cgrid%mmean_soil_water)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_soil_water,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_SOIL_WATER :12:hist:mont:dcyc')
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
      end if

      if(associated(cgrid%mmean_soil_mstpot)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_soil_mstpot,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_SOIL_MSTPOT :12:hist:mont:dcyc')
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
      end if
      
      if(associated(cgrid%mmean_transloss)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_transloss,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_TRANSLOSS :12:hist:mont:dcyc')
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA')
      end if
      
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!



      return     
   end subroutine filltab_edtype_p12
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine will fill the pointer table with the polygon-level variables         !
   ! (edtype) that have three dimensions (nzg,ndcycle,npolygons) and are real (type -12).  !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_edtype_m12(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      use ed_var_tables, only : vtable_edio_r  & ! sub-routine
                              , metadata_edio  ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype), target        :: cgrid
      integer     , intent(in)    :: init
      integer     , intent(in)    :: igr
      integer     , intent(in)    :: var_len
      integer     , intent(in)    :: max_ptrs
      integer     , intent(in)    :: var_len_global
      integer     , intent(inout) :: nvar
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: npts
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 3-D block, with soil depth class and diurnal cycle.  All vari-     !
      ! ables must have the number of points defined by npts.                              !
      !------------------------------------------------------------------------------------!
      npts = cgrid%npolygons *  nzg * ndcycle

      if(associated(cgrid%qmean_soil_temp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_soil_temp,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_SOIL_TEMP :-12:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_soil_water)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_soil_water,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_SOIL_WATER :-12:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(cgrid%qmean_soil_mstpot)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%qmean_soil_mstpot,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'QMEAN_SOIL_MSTPOT :-12:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!


      return     
   end subroutine filltab_edtype_m12
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine will fill the pointer table with the polygon-level variables         !
   ! (edtype) that have two dimensions (n_pft,npolygons) and are real (type 14).           !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_edtype_p14(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      use ed_var_tables, only : vtable_edio_r  & ! sub-routine
                              , metadata_edio  ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype), target        :: cgrid
      integer     , intent(in)    :: init
      integer     , intent(in)    :: igr
      integer     , intent(in)    :: var_len
      integer     , intent(in)    :: max_ptrs
      integer     , intent(in)    :: var_len_global
      integer     , intent(inout) :: nvar
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: npts
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 2-D block, with PFT class.  All variables must have the number of  !
      ! points defined by npts.                                                            !
      !------------------------------------------------------------------------------------!
      npts = cgrid%npolygons *  n_pft


      if(associated(cgrid%lai_pft)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%lai_pft,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'LAI_PFT :14:hist:anal:dail') 
         call metadata_edio(nvar,igr,'Leaf Area Index','[m2/m2]','NA') 
      end if
      
      if(associated(cgrid%wpa_pft)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%wpa_pft,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'WPA_PFT :14:hist:anal:dail') 
         call metadata_edio(nvar,igr,'Wood Projected Area','[m2/m2]','NA') 
      end if
      
      if(associated(cgrid%wai_pft)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%wai_pft,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'WAI_PFT :14:hist:anal:dail') 
         call metadata_edio(nvar,igr,'Wood Area Index','[m2/m2]','NA') 
      end if
      
      if(associated(cgrid%mmean_lai_pft)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_lai_pft,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_LAI_PFT :14:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%mmean_wpa_pft)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_wpa_pft,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_WPA_PFT :14:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%mmean_wai_pft)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_wai_pft,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_WAI_PFT :14:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%bseeds_pft)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%bseeds_pft,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'BSEEDS_PFT :14:hist:dail:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if(associated(cgrid%agb_pft)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%agb_pft,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AGB_PFT :14:hist:dail:mont:dcyc') 
         call metadata_edio(nvar,igr,'Above-ground biomass by PFT','[kgC/m2]','NA') 
      end if
      
      if(associated(cgrid%ba_pft)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%ba_pft,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'BA_PFT :14:hist:dail:mont:dcyc') 
         call metadata_edio(nvar,igr,'Basal area by PFT','[cm2/m2]','NA') 
      end if

      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!


      return     
   end subroutine filltab_edtype_p14
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine will fill the pointer table with the polygon-level variables         !
   ! (edtype) that have two dimensions (n_dbh,npolygons) and are real (type 16).           !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_edtype_p16(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      use ed_var_tables, only : vtable_edio_r  & ! sub-routine
                              , metadata_edio  ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype), target        :: cgrid
      integer     , intent(in)    :: init
      integer     , intent(in)    :: igr
      integer     , intent(in)    :: var_len
      integer     , intent(in)    :: max_ptrs
      integer     , intent(in)    :: var_len_global
      integer     , intent(inout) :: nvar
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: npts
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 2-D block, with size (DBH) class.  All variables must have the     !
      ! number of points defined by npts.                                                  !
      !------------------------------------------------------------------------------------!
      npts = cgrid%npolygons *  n_dbh

      if(associated(cgrid%dmean_gpp_dbh)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%dmean_gpp_dbh,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_GPP_DBH :16:hist:dail') 
         call metadata_edio(nvar,igr,'Polygon Averaged by DBH, Daily Integrated Gross Primary Production' &
              ,'[kgC/m2/yr]','ipoly - ndbh') 
      end if

      if(associated(cgrid%mmean_gpp_dbh)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%mmean_gpp_dbh,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_GPP_DBH :16:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!

      return     
   end subroutine filltab_edtype_p16
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine will fill the pointer table with the polygon-level variables         !
   ! (edtype) that have two dimensions (13 months,npolygons) and are real (type 19).       !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_edtype_p19(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      use ed_var_tables, only : vtable_edio_r  & ! sub-routine
                              , metadata_edio  ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype), target        :: cgrid
      integer     , intent(in)    :: init
      integer     , intent(in)    :: igr
      integer     , intent(in)    :: var_len
      integer     , intent(in)    :: max_ptrs
      integer     , intent(in)    :: var_len_global
      integer     , intent(inout) :: nvar
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: npts
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 2-D block, with 13 months.  All variables must have the            !
      ! number of points defined by npts.                                                  !
      !------------------------------------------------------------------------------------!
      npts = cgrid%npolygons * 13

      if(associated(cgrid%workload)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%workload,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'WORKLOAD :19:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Disturbance Rates','[NA]','NA') 
      end if
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!


      return     
   end subroutine filltab_edtype_p19
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine will fill the pointer table with the polygon-level variables         !
   ! (edtype) that have three dimensions (n_dbh,n_pft,npolygons) and are real (type 146).  !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_edtype_p146(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      use ed_var_tables, only : vtable_edio_r  & ! sub-routine
                              , metadata_edio  ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype), target        :: cgrid
      integer     , intent(in)    :: init
      integer     , intent(in)    :: igr
      integer     , intent(in)    :: var_len
      integer     , intent(in)    :: max_ptrs
      integer     , intent(in)    :: var_len_global
      integer     , intent(inout) :: nvar
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: npts
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 3-D block, with PFT and size class.  All variables must have the   !
      ! number of points defined by npts.                                                  !
      !------------------------------------------------------------------------------------!
      npts = cgrid%npolygons * n_pft * n_dbh

      if (associated(cgrid%basal_area)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%basal_area,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'BASAL_AREA :146:hist:anal:dail:mont:dcyc:year') 
         call metadata_edio(nvar,igr,'Polygon basal area profile','[cm2/m2]','ipoly - n_dbh - n_pft')
      end if
      
      if (associated(cgrid%agb)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%agb,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AGB :146:hist:anal:dail:mont:dcyc:year') 
         call metadata_edio(nvar,igr,'Polygon above ground biomass profile','[kgC/m2]','ipoly - n_dbh - n_pft')
      end if
      
      if (associated(cgrid%pldens)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%pldens,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'PLDENS :146:hist:anal:dail:dcyc:mont:year') 
         call metadata_edio(nvar,igr,'Polygon plant density profile','[#/m2]','ipoly - n_dbh - n_pft')
      end if
      
      if (associated(cgrid%bseeds)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%bseeds,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'BSEEDS :146:hist:anal:dail:dcyc:mont:year') 
         call metadata_edio(nvar,igr,'Polygon seed biomass','[kgC/m2]','ipoly - n_dbh - n_pft')
      end if
      
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!


      return     
   end subroutine filltab_edtype_p146
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine will fill the pointer table with the polygon-level variables         !
   ! (edtype) that have three dimensions (3,4,npolygons) and are real (type 199).          !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_edtype_p199(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      use ed_var_tables, only : vtable_edio_r  & ! sub-routine
                              , metadata_edio  ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype), target        :: cgrid
      integer     , intent(in)    :: init
      integer     , intent(in)    :: igr
      integer     , intent(in)    :: var_len
      integer     , intent(in)    :: max_ptrs
      integer     , intent(in)    :: var_len_global
      integer     , intent(inout) :: nvar
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: npts
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 3-D block, with LAI and 4 vars.  All variables must have the       !
      ! number of points defined by npts.                                                  !
      !------------------------------------------------------------------------------------!
      npts = cgrid%npolygons * 3 * 4

      if (associated(cgrid%avg_lai_ebalvars)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%avg_lai_ebalvars,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'AVG_LAI_EBALVARS :199:anal') 
         call metadata_edio(nvar,igr,'Polygon Average Energy Balance Variables','[variable]','ipoly - 4 - 3') 
      end if
      
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!



      return     
   end subroutine filltab_edtype_p199
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine will fill the pointer table with the polygon-level variables         !
   ! (edtype) that have three dimensions (n_dist_types,n_dist_types,npolygons) and are     !
   ! real (type 155).                                                                      !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_edtype_p155(cgrid,igr,init,var_len,var_len_global,max_ptrs,nvar)
      use ed_var_tables, only : vtable_edio_r  & ! sub-routine
                              , metadata_edio  ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype), target        :: cgrid
      integer     , intent(in)    :: init
      integer     , intent(in)    :: igr
      integer     , intent(in)    :: var_len
      integer     , intent(in)    :: max_ptrs
      integer     , intent(in)    :: var_len_global
      integer     , intent(inout) :: nvar
      !----- Local variables. -------------------------------------------------------------!
      integer                     :: npts
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 3-D block, with disturbance transitions.  All variables must have  !
      ! the number of points defined by npts.                                              !
      !------------------------------------------------------------------------------------!
      npts = cgrid%npolygons * n_dist_types * n_dist_types

      if(associated(cgrid%disturbance_rates)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cgrid%disturbance_rates,nvar,igr,init,cgrid%pyglob_id, &
              var_len,var_len_global,max_ptrs,'DISTURBANCE_RATES :155:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Disturbance Rates','[NA]','NA') 
      end if
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!


      return     
   end subroutine filltab_edtype_p155
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine filltab_polygontype(igr,ipy,init)

      use ed_var_tables, only : vtable_edio_r & ! sub-routine
                              , vtable_edio_i & ! sub-routine
                              , metadata_edio ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer, intent(in) :: init
      integer, intent(in) :: igr
      integer, intent(in) :: ipy
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer :: cpoly
      integer                    :: var_len
      integer                    :: max_ptrs
      integer                    :: var_len_global
      integer                    :: nvar
      integer                    :: npts
      !------------------------------------------------------------------------------------!



      !----- Check whether this polygon is allocated. -------------------------------------!
      if (.not.associated(edgrid_g(igr)%polygon(ipy)%sipa_id)) then
         write (unit=*,fmt='(a)')       '------------------------------------------------'
         write (unit=*,fmt='(a)')       ' Returning from filltab_polygontype'
         write (unit=*,fmt='(a,1x,i5)') ' IGR  = ',igr
         write (unit=*,fmt='(a,1x,i5)') ' IPY  = ',ipy
         write (unit=*,fmt='(a,1x,i5)') ' INIT = ',init
         write (unit=*,fmt='(a)')       '------------------------------------------------'
         return
      end if
      !------------------------------------------------------------------------------------!



      !----- Assign a pointer to the current polygon. -------------------------------------!
      cpoly => edgrid_g(igr)%polygon(ipy)
      !------------------------------------------------------------------------------------!



      !----- Define some variables to be used by all variables. ---------------------------!
      var_len        = cpoly%nsites
      var_len_global = edgrid_g(igr)%nsites_global
      max_ptrs       = edgrid_g(igr)%npolygons
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Initialise the variable number, which should be unique for each variable.      !
      ! nioglobal has the current position after all the parent structure variables have   !
      ! been positioned in the variable table.                                             !
      !------------------------------------------------------------------------------------!
      nvar=nioglobal+niogrid
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !       This part should have only 1-D vectors (nsites).  Notice that they all use   !
      ! npts = cpoly%nsites.  Add only variables of types 20 and 21 here.                  !
      !------------------------------------------------------------------------------------!
      npts = cpoly%nsites


      if (associated(cpoly%sipa_id)) then
         nvar=nvar+1
         call vtable_edio_i(npts,cpoly%sipa_id,nvar,igr,init,cpoly%siglob_id, &
              var_len,var_len_global,max_ptrs,'SIPA_ID :20:hist:dail:mont:dcyc:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if (associated(cpoly%sipa_n)) then
         nvar=nvar+1
         call vtable_edio_i(npts,cpoly%sipa_n,nvar,igr,init,cpoly%siglob_id, &
              var_len,var_len_global,max_ptrs,'SIPA_N :20:hist:dail:mont:dcyc:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if (associated(cpoly%patch_count)) then
         nvar=nvar+1
         call vtable_edio_i(npts,cpoly%patch_count,nvar,igr,init,cpoly%siglob_id, &
              var_len,var_len_global,max_ptrs,'PATCH_COUNT :20:hist:dail:mont:dcyc:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if (associated(cpoly%sitenum)) then
         nvar=nvar+1
         call vtable_edio_i(npts,cpoly%sitenum,nvar,igr,init,cpoly%siglob_id, &
              var_len,var_len_global,max_ptrs,'SITENUM :20:hist:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%lsl)) then
         nvar=nvar+1
           call vtable_edio_i(npts,cpoly%lsl,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'LSL_SI :20:hist:dail:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if   
      
      if (associated(cpoly%num_landuse_years)) then
         nvar=nvar+1
           call vtable_edio_i(npts,cpoly%num_landuse_years,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'NUM_LANDUSE_YEARS :20:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%hydro_next)) then
         nvar=nvar+1
           call vtable_edio_i(npts,cpoly%hydro_next,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'HYDRO_NEXT :20:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%hydro_prev)) then
         nvar=nvar+1
           call vtable_edio_i(npts,cpoly%hydro_prev,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'HYDRO_PREV :20:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%plantation)) then
         nvar=nvar+1
           call vtable_edio_i(npts,cpoly%plantation,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'PLANTATION_SI :20:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if 

      if (associated(cpoly%agri_stocking_pft)) then
         nvar=nvar+1
           call vtable_edio_i(npts,cpoly%agri_stocking_pft,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'AGRI_STOCKING_PFT :20:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%plantation_stocking_pft)) then
         nvar=nvar+1
           call vtable_edio_i(npts,cpoly%plantation_stocking_pft,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'PLANTATION_STOCKING_PFT :20:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%nat_dist_type)) then
         nvar=nvar+1
           call vtable_edio_i(npts,cpoly%nat_dist_type,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'NAT_DIST_TYPE :20:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%area)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%area,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'AREA_SI:21:hist:dail:mont:dcyc:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%patch_area)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%patch_area,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'PATCH_AREA:21:hist:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%elevation)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%elevation,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'ELEVATION :21:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%slope)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%slope,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'SLOPE :21:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%aspect)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%aspect,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'ASPECT :21:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%TCI)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%TCI,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'TCI :21:hist') 
         call metadata_edio(nvar,igr,'Topographic convergence index','[NA]','NA') 
      end if      

      if (associated(cpoly%pptweight)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%TCI,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'pptweight :21:hist') 
         call metadata_edio(nvar,igr,'precip lapse weighting','[NA]','NA') 
      end if      

      if (associated(cpoly%moist_W)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%moist_W,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'MOIST_W :21:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%moist_f)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%moist_f,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'MOIST_F :21:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if  

      if (associated(cpoly%moist_tau)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%moist_tau,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'MOIST_TAU :21:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%moist_zi)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%moist_zi,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'MOIST_ZI :21:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if 

      if (associated(cpoly%baseflow)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%baseflow  ,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'BASEFLOW_SI :21:hist') 
         call metadata_edio(nvar,igr,'loss of water from site to watershed discharge','[kg/m2/s]','NA') 
      end if 

      if (associated(cpoly%min_monthly_temp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%min_monthly_temp,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'MIN_MONTHLY_TEMP :21:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%agri_stocking_density)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%agri_stocking_density,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'AGRI_STOCKING_DENSITY :21:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%plantation_stocking_density)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%plantation_stocking_density,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'PLANTATION_STOCKING_DENSITY :21:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%primary_harvest_memory)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%primary_harvest_memory,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'PRIMARY_HARVEST_MEMORY :21:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%secondary_harvest_memory)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%secondary_harvest_memory,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'SECONDARY_HARVEST_MEMORY:21:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%fire_disturbance_rate)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%fire_disturbance_rate,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'FIRE_DISTURBANCE_RATE :21:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%ignition_rate)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%ignition_rate,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'IGNITION_RATE :21:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%nat_disturbance_rate)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%nat_disturbance_rate,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'NAT_DISTURBANCE_RATE :21:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if (associated(cpoly%rad_avg)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%rad_avg,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'RAD_AVG :21:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if


      if(associated(cpoly%dmean_co2_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cpoly%dmean_co2_residual,nvar,igr,init,cpoly%siglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_CO2_RESIDUAL_SI :21:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean of residual canopy CO2 ','[umol/m2/day]','ipoly') 
      end if

      if(associated(cpoly%dmean_water_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cpoly%dmean_water_residual,nvar,igr,init,cpoly%siglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_WATER_RESIDUAL_SI :21:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean of residual water ','[kg/m2/day]','ipoly') 
      end if

      if(associated(cpoly%dmean_energy_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cpoly%dmean_energy_residual,nvar,igr,init,cpoly%siglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_ENERGY_RESIDUAL_SI :21:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean of residual water ','[J/m2/day]','ipoly') 
      end if

      if(associated(cpoly%mmean_co2_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cpoly%mmean_co2_residual,nvar,igr,init,cpoly%siglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_CO2_RESIDUAL_SI :21:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean of residual canopy CO2 ','[umol/m2/day]','ipoly') 
      end if

      if(associated(cpoly%mmean_water_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cpoly%mmean_water_residual,nvar,igr,init,cpoly%siglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_WATER_RESIDUAL_SI :21:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean of residual water ','[kg/m2/day]','ipoly') 
      end if

      if(associated(cpoly%mmean_energy_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cpoly%mmean_energy_residual,nvar,igr,init,cpoly%siglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_ENERGY_RESIDUAL_SI :21:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean of residual energy ','[J/m2/day]','ipoly') 
      end if

      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 2-D block, with dimensions being nsites and nzg.  Make sure to     !
      ! include only variables of types 22 and 220 here, as they will all use the same     !
      ! npts.                                                                              !
      !------------------------------------------------------------------------------------!
      npts = cpoly%nsites * nzg

  
      if (associated(cpoly%ntext_soil)) then
         nvar=nvar+1
           call vtable_edio_i(npts,cpoly%ntext_soil,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'NTEXT_SOIL_SI :220:hist:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 2-D block, with dimensions being nsites and n_pft.  Make sure to   !
      ! include only variables of type 24 here, as they will all use the same npts.        !
      !------------------------------------------------------------------------------------!
      npts = cpoly%nsites * n_pft

      if (associated(cpoly%lai_pft)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%lai_pft,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'LAI_PFT_SI :24:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if (associated(cpoly%wpa_pft)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%wpa_pft,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'WPA_PFT_SI :24:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if (associated(cpoly%wai_pft)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%wai_pft,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'WAI_PFT_SI :24:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%green_leaf_factor)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%green_leaf_factor,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'GREEN_LEAF_FACTOR :24:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%leaf_aging_factor)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%leaf_aging_factor,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'LEAF_AGING_FACTOR :24:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 2-D block, with dimensions being nsites and n_dist_types.  Make    !
      ! sure to include only variables of type 25 here, as they will all use the same      !
      ! npts.                                                                              !
      !------------------------------------------------------------------------------------!
      npts = cpoly%nsites * n_dist_types

      if (associated(cpoly%loss_fraction)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%loss_fraction,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'LOSS_FRACTION :25:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 2-D block, with dimensions being nsites and 12 months.  Make sure  !
      ! to include only variables of type 29 here, as they will all use the same npts.     !
      !------------------------------------------------------------------------------------!
      npts = cpoly%nsites * 12

      if (associated(cpoly%lambda_fire)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%lambda_fire,nvar,igr,init,cpoly%siglob_id         &
                             ,var_len,var_len_global,max_ptrs,'LAMBDA_FIRE :29:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 3-D block, with dimensions being n_pft and n_dbh.  Make sure to    !
      ! include only variables of type 246 here, as they will all use the same npts.       !
      !------------------------------------------------------------------------------------!
      npts = cpoly%nsites * n_pft * n_dbh

      if (associated(cpoly%basal_area)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%basal_area,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'BASAL_AREA_SI :246:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%agb)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%agb,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'AGB_SI :246:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%pldens)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%pldens,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'PLDENS_SI :246:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%bseeds)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%bseeds,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'BSEEDS_SI :246:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%basal_area_growth)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%basal_area_growth,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'BASAL_AREA_GROWTH :246:hist:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%agb_growth)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%agb_growth,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'AGB_GROWTH :246:hist:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%basal_area_mort)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%basal_area_mort,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'BASAL_AREA_MORT :246:hist:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%basal_area_cut)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%basal_area_cut,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'BASAL_AREA_CUT :246:year:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%agb_mort)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%agb_mort,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'AGB_MORT :246:hist:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%agb_cut)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%agb_cut,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'AGB_CUT :246:hist:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%avg_soil_rootfrac)) then
         nvar=nvar+1
         call vtable_edio_r(npts,cpoly%avg_soil_rootfrac,nvar,igr,init,cpoly%siglob_id, &
              var_len,var_len_global,max_ptrs,'SITE_SOIL_ROOTFRAC :12:hist:mont:year') 
         call metadata_edio(nvar,igr,'Site Average Root Fraction','[kg m-3]','ipoly - nzg') 
      end if


      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     This is the 3-D block, with dimensions being nsites and n_dist_types,          !
      ! n_dist_types (from and to).  Make sure to include only variables of type 255 here, !
      ! as they will all use the same npts.                                                !
      !------------------------------------------------------------------------------------!
      npts = cpoly%nsites * n_dist_types * n_dist_types

      if (associated(cpoly%disturbance_memory)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%disturbance_memory,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'DISTURBANCE_MEMORY :255:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(cpoly%disturbance_rates)) then
         nvar=nvar+1
           call vtable_edio_r(npts,cpoly%disturbance_rates,nvar,igr,init,cpoly%siglob_id, &
           var_len,var_len_global,max_ptrs,'DISTURBANCE_RATES_SI :255:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!


      !----- Save the number of site-level (polygontype) variables that go to the output. -!
      if (init == 0) niopoly=nvar-niogrid-nioglobal
      !------------------------------------------------------------------------------------!

      return
   end subroutine filltab_polygontype
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine will fill the variable table with the sitetype variables (patch- !
   ! -level).                                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_sitetype(igr,ipy,isi,init)
     use soil_bgc, only: filltab_sitetype_sbgc
      use ed_var_tables, only : vtable_edio_r  & ! sub-routine
                              , vtable_edio_i  & ! sub-routine
                              , metadata_edio  ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer       , intent(in) :: init
      integer       , intent(in) :: igr
      integer       , intent(in) :: ipy
      integer       , intent(in) :: isi
      !----- Local variables. -------------------------------------------------------------!
      type(sitetype), pointer    :: csite
      integer                    :: var_len
      integer                    :: max_ptrs
      integer                    :: var_len_global
      integer                    :: nvar
      integer                    :: npts
      !------------------------------------------------------------------------------------!



      !----- Assign a pointer to the current polygon. -------------------------------------!
      csite => edgrid_g(igr)%polygon(ipy)%site(isi)
      !------------------------------------------------------------------------------------!



      !----- Define some variables to be used by all variables. ---------------------------!
      var_len        = csite%npatches
      var_len_global = edgrid_g(igr)%npatches_global
      max_ptrs       = get_nsites(edgrid_g(igr))
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Initialise the variable number, which should be unique for each variable.      !
      ! nioglobal has the current position after all the parent structure variables have   !
      ! been positioned in the variable table.                                             !
      !------------------------------------------------------------------------------------!
      nvar=nioglobal+niogrid+niopoly
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !       This part should have only 1-D vectors (npatches).  Notice that they all use !
      ! npts = csite%npatches.  Add only variables of types 30 and 31 here.                !
      !------------------------------------------------------------------------------------!
      npts = csite%npatches


      if (associated(csite%paco_id)) then
         nvar=nvar+1
         call vtable_edio_i(npts,csite%paco_id,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'PACO_ID :30:hist:anal:dail:mont:dcyc:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%paco_n)) then
         nvar=nvar+1
           call vtable_edio_i(npts,csite%paco_n,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'PACO_N :30:hist:anal:dail:mont:dcyc:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if (associated(csite%dist_type)) then
         nvar=nvar+1
           call vtable_edio_i(npts,csite%dist_type,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'DIST_TYPE :30:hist:anal:dail:mont:dcyc:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%plantation)) then
         nvar=nvar+1
           call vtable_edio_i(npts,csite%plantation,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'PLANTATION :30:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%nlev_sfcwater)) then
         nvar=nvar+1
           call vtable_edio_i(npts,csite%nlev_sfcwater,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'NLEV_SFCWATER :30:hist:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%fuse_flag)) then
         nvar=nvar+1
           call vtable_edio_i(npts,csite%fuse_flag,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'FUSE_FLAG :30:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%age)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%age,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'AGE :31:hist:anal:dail:mont:dcyc:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%area)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%area,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'AREA :31:hist:anal:dail:mont:dcyc:year') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if (associated(csite%sum_dgd)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%sum_dgd,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'SUM_DGD :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%sum_chd)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%sum_chd,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'SUM_CHD :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(allocated(csite%sbgc%fast_soil_C))then
         call filltab_sitetype_sbgc(nvar, npts, csite%sbgc, igr, init, csite%paglob_id,  &
              var_len, var_len_global, max_ptrs)
      endif

      if (associated(csite%can_theiv)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%can_theiv,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CAN_THEIV :31:hist') 
         call metadata_edio(nvar,igr,'Canopy air ice-vapour equivalent potential temperature','[K]','NA') 
      end if

      if (associated(csite%can_temp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%can_temp,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CAN_TEMP :31:hist') 
         call metadata_edio(nvar,igr,'Canopy air temperature','[K]','NA') 
      end if

      if (associated(csite%can_shv)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%can_shv,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CAN_SHV :31:hist') 
         call metadata_edio(nvar,igr,'Canopy air specific humidity','[kg_H2O/kg_Air]','NA') 
      end if

      if (associated(csite%can_co2)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%can_co2,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CAN_CO2 :31:hist') 
         call metadata_edio(nvar,igr,'Canopy air CO2 mixing ratio','[umol/mol]','NA') 
      end if
   
      if (associated(csite%can_rhos)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%can_rhos,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CAN_RHOS :31:hist') 
         call metadata_edio(nvar,igr,'Canopy air Density','[kg/m3]','NA') 
      end if
   
      if (associated(csite%can_prss)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%can_prss,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CAN_PRSS :31:hist') 
         call metadata_edio(nvar,igr,'Canopy air Pressure','[Pa]','NA') 
      end if
   
      if (associated(csite%can_theta)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%can_theta,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CAN_THETA :31:hist') 
         call metadata_edio(nvar,igr,'Canopy air Potential temperature','[K]','NA') 
      end if
     
      if (associated(csite%can_depth)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%can_depth,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CAN_DEPTH :31:hist') 
         call metadata_edio(nvar,igr,'Canopy depth','[m]','NA') 
      end if
     
      if (associated(csite%opencan_frac)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%opencan_frac,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'OPENCAN_FRAC :31:hist')
           call metadata_edio(nvar,igr,'Fraction of open canopy','[---]','NA') 
      end if
     
      if (associated(csite%ggbare)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ggbare,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'GGBARE :31:hist')
           call metadata_edio(nvar,igr,'Bare ground conductance','[m/s]','NA') 
      end if

      if (associated(csite%ggveg)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ggveg,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'GGVEG :31:hist')
           call metadata_edio(nvar,igr,'Vegetated ground conductance','[m/s]','NA') 
      end if

      if (associated(csite%ggnet)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ggnet,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'GGNET :31:hist')
           call metadata_edio(nvar,igr,'Net ground conductance','[m/s]','NA') 
      end if

      if (associated(csite%ggsoil)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ggsoil,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'GGSOIL :31:hist')
           call metadata_edio(nvar,igr,'Soil conductance for evaporation','[m/s]','NA') 
      end if

      if (associated(csite%lambda_light)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%lambda_light,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'LAMBDA_LIGHT :31:hist') 
         call metadata_edio(nvar,igr,'Light extinction','[m2/,2]','NA') 
      end if
     
      if (associated(csite%dmean_lambda_light)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%dmean_lambda_light,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'DMEAN_LAMBDA_LIGHT :31:hist:dail') 
         call metadata_edio(nvar,igr,'Light extinction','[m2/,2]','NA') 
      end if
     
      if (associated(csite%mmean_lambda_light)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mmean_lambda_light,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'MMEAN_LAMBDA_LIGHT :31:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Light extinction','[m2/,2]','NA') 
      end if

      if (associated(csite%lai)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%lai,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'LAI_PA :31:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%wpa)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%wpa,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'WPA_PA :31:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%wai)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%wai,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'WAI_PA :31:hist:dail') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%ground_shv)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ground_shv,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'GROUND_SHV :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%ground_ssh)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ground_ssh,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'GROUND_SSH :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%ground_temp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ground_temp,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'GROUND_TEMP :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%ground_fliq)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ground_fliq,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'GROUND_FLIQ :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%rough)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rough,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'ROUGH :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if


      if (associated(csite%mean_rh)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mean_rh,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'MEAN_RH :31:hist') 
         call metadata_edio(nvar,igr,'Heterotrophic respiration','[umol/m2/s]','ipatch') 
      end if

      if (associated(csite%dmean_rh)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%dmean_rh,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'DMEAN_RH_PA :31:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean of heterotrophic respiration','[umol/m2/s]','ipatch') 
      end if

      if (associated(csite%mmean_rh)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mmean_rh,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'MMEAN_RH_PA :31:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean of heterotrophic respiration','[umol/m2/s]','ipatch') 
      end if

      if (associated(csite%dmean_albedo)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%dmean_albedo,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'DMEAN_ALBEDO_PA :31:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean of albedo','[---]','ipatch') 
      end if

      if (associated(csite%dmean_albedo_beam)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%dmean_albedo_beam,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'DMEAN_ALBEDO_BEAM_PA :31:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean of beam albedo','[---]','ipatch') 
      end if

      if (associated(csite%dmean_albedo_diffuse)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%dmean_albedo_diffuse,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'DMEAN_ALBEDO_DIFFUSE_PA :31:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean of diffuse albedo','[---]','ipatch') 
      end if

      if (associated(csite%mmean_albedo)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mmean_albedo,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'MMEAN_ALBEDO_PA :31:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean of albedo','[---]','ipatch') 
      end if

      if (associated(csite%mmean_albedo_beam)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mmean_albedo_beam,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'MMEAN_ALBEDO_BEAM_PA :31:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean of beam albedo','[---]','ipatch') 
      end if

      if (associated(csite%mmean_albedo_diffuse)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mmean_albedo_diffuse,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'MMEAN_ALBEDO_DIFFUSE_PA :31:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'Monthly mean of diffuse albedo','[---]','ipatch') 
      end if

      if (associated(csite%mean_nep)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mean_nep,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'MEAN_NEP :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%wbudget_loss2atm)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%wbudget_loss2atm,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'WBUDGET_LOSS2ATM :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%wbudget_denseffect)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%wbudget_denseffect,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'WBUDGET_DENSEFFECT :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%wbudget_precipgain)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%wbudget_precipgain,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'WBUDGET_PRECIPGAIN :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%wbudget_loss2runoff)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%wbudget_loss2runoff,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'WBUDGET_LOSS2RUNOFF :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%wbudget_loss2drainage)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%wbudget_loss2drainage,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'WBUDGET_LOSS2DRAINAGE :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%wbudget_initialstorage)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%wbudget_initialstorage,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'WBUDGET_INITIALSTORAGE :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%wbudget_residual)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%wbudget_residual,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'WBUDGET_RESIDUAL :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%ebudget_loss2atm)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ebudget_loss2atm,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'EBUDGET_LOSS2ATM :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%ebudget_denseffect)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ebudget_denseffect,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'EBUDGET_DENSEFFECT :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%ebudget_loss2runoff)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ebudget_loss2runoff,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'EBUDGET_LOSS2RUNOFF :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%ebudget_loss2drainage)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ebudget_loss2drainage,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'EBUDGET_LOSS2DRAINAGE :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%ebudget_netrad)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ebudget_netrad,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'EBUDGET_NETRAD :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%ebudget_precipgain)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ebudget_precipgain,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'EBUDGET_PRECIPGAIN :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%ebudget_initialstorage)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ebudget_initialstorage,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'EBUDGET_INITIALSTORAGE :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%ebudget_residual)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ebudget_residual,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'EBUDGET_RESIDUAL :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%co2budget_initialstorage)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%co2budget_initialstorage,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CO2BUDGET_INITIALSTORAGE :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%co2budget_residual)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%co2budget_residual,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CO2BUDGET_RESIDUAL :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%co2budget_loss2atm)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%co2budget_loss2atm,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CO2BUDGET_LOSS2ATM :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%co2budget_denseffect)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%co2budget_denseffect,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CO2BUDGET_DENSEFFECT :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%co2budget_gpp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%co2budget_gpp,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CO2BUDGET_GPP :31:hist') 
         call metadata_edio(nvar,igr,'Patch total gross primary productivity per timestep','[umol/m2/dtlsm]','NA') 
      end if


      if (associated(csite%avg_daily_temp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%avg_daily_temp,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'AVG_DAILY_TEMP :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if  

      if (associated(csite%co2budget_plresp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%co2budget_plresp,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CO2BUDGET_PLRESP :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%co2budget_rh)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%co2budget_rh,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CO2BUDGET_RH :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%today_A_decomp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%today_A_decomp,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'TODAY_A_DECOMP :31:hist') 
         call metadata_edio(nvar,igr,'NOT A DIAGNOSTIC-WILL ZERO AT END OF DAY','[NA]','NA') 
      end if

      if (associated(csite%today_Af_decomp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%today_Af_decomp,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'TODAY_AF_DECOMP :31:hist') 
         call metadata_edio(nvar,igr,'NOT A DIAGNOSTIC-WILL ZERO AT END OF DAY','[NA]','NA') 
      end if

      if (associated(csite%dmean_A_decomp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%dmean_A_decomp,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'DMEAN_A_DECOMP :31:hist:dail') 
         call metadata_edio(nvar,igr,'A factor for decomposition','[--]','ipatch') 
      end if

      if (associated(csite%dmean_Af_decomp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%dmean_Af_decomp,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'DMEAN_AF_DECOMP :31:hist:dail') 
         call metadata_edio(nvar,igr,'A factor for decomposition, including N','[--]','ipatch') 
      end if

      if (associated(csite%mmean_A_decomp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mmean_A_decomp,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'MMEAN_A_DECOMP :31:hist:dcyc:mont') 
         call metadata_edio(nvar,igr,'Monthly mean of A factor for decomposition','[--]','ipatch') 
      end if

      if (associated(csite%mmean_Af_decomp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mmean_Af_decomp,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'MMEAN_AF_DECOMP :31:hist:dcyc:mont') 
         call metadata_edio(nvar,igr,'Monthly mean A factor for decomposition, including N','[--]','ipatch') 
      end if


      if (associated(csite%veg_rough)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%veg_rough,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'VEG_ROUGH :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%veg_height)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%veg_height ,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'VEG_HEIGHT :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%veg_displace)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%veg_displace ,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'VEG_DISPLACE :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%fsc_in)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%fsc_in,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'FSC_IN :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%ssc_in)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ssc_in,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'SSC_IN :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%ssl_in)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ssl_in,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'SSL_IN :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%fsn_in)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%fsn_in,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'FSN_IN :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if



      if (associated(csite%fsp_in)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%fsp_in,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'FSP_IN :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      if (associated(csite%stsp_in)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%stsp_in,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'STSP_IN :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if





      if (associated(csite%total_plant_P_uptake)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%total_plant_P_uptake,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'TOTAL_PLANT_P_UPTAKE :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      if (associated(csite%total_plant_nitrogen_uptake)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%total_plant_nitrogen_uptake,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'TOTAL_PLANT_NITROGEN_UPTAKE :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%mineralized_N_loss)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mineralized_N_loss,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'NMIN_LOSS :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if


      if (associated(csite%mineralized_N_input)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mineralized_N_input,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'NMIN_INPUT :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%rshort_g)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rshort_g,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RSHORT_G :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%rshort_g_beam)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rshort_g_beam,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RSHORT_G_BEAM :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%rshort_g_diffuse)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rshort_g_diffuse,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RSHORT_G_DIFFUSE :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%rlong_g)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rlong_g,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RLONG_G :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%rlong_g_surf)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rlong_g_surf,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RLONG_G_SURF :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%rlong_g_incid)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rlong_g_incid,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RLONG_G_INCID :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%rlong_s)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rlong_s,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RLONG_S :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%rlong_s_surf)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rlong_s_surf,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RLONG_S_SURF :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%rlong_s_incid)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rlong_s_incid,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RLONG_S_INCID :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%albedo)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%albedo,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'ALBEDO :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%albedo_beam)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%albedo_beam,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'ALBEDO_BEAM :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%albedo_diffuse)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%albedo_diffuse,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'ALBEDO_DIFFUSE :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%rlongup)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rlongup,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RLONGUP :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%rlong_albedo)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rlong_albedo,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RLONG_ALBEDO :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%total_sfcw_depth)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%total_sfcw_depth,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'TOTAL_SFCW_DEPTH :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%snowfac)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%snowfac,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'SNOWFAC :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%A_decomp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%A_decomp,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'A_DECOMP :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%f_decomp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%f_decomp,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'F_DECOMP :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%rh)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rh,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RH :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%cwd_rh)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%cwd_rh,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CWD_RH :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%plant_ag_biomass)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%plant_ag_biomass,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'PLANT_AG_BIOMASS :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if (associated(csite%mean_wflux)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mean_wflux,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'MEAN_WFLUX :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%mean_latflux)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mean_latflux,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'MEAN_LATFLUX :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%mean_hflux)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mean_hflux,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'MEAN_HFLUX :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%mean_runoff)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mean_runoff,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'MEAN_RUNOFF :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%mean_qrunoff)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mean_qrunoff,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'MEAN_QRUNOFF :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if (associated(csite%htry)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%htry,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'HTRY :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      
      if (associated(csite%avg_rk4step)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%avg_rk4step,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'AVG_RK4STEP :31:hist:anal') 
         call metadata_edio(nvar,igr,'Average time step used by Runge-Kutta','[s]','ipatch') 
      end if
      
      if (associated(csite%avg_available_water)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%avg_available_water,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'AVG_AVAILABLE_WATER_PA :31:hist:anal') 
         call metadata_edio(nvar,igr,'Average available water for transpiration','[kg/m2]','ipatch') 
      end if
      
      if (associated(csite%dmean_rk4step)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%dmean_rk4step,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'DMEAN_RK4STEP :31:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean time step used by Runge-Kutta','[s]','ipatch') 
      end if
      
      if (associated(csite%mmean_rk4step)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%mmean_rk4step,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'MMEAN_RK4STEP :31:hist:dcyc:mont') 
         call metadata_edio(nvar,igr,'Monthly mean time step used by Runge-Kutta','[s]','ipatch') 
      end if
      
      if (associated(csite%ustar)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ustar,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'USTAR :31:hist') 
           call metadata_edio(nvar,igr,'Patch level friction velocity','[m/s]','ipatch') 
      end if

      if (associated(csite%tstar)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%tstar,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'TSTAR :31:hist') 
         call metadata_edio(nvar,igr,'patch level heat transfer atm->canopy','[K]','ipatch') 
      end if

      if (associated(csite%qstar)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%qstar,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'QSTAR :31:hist') 
         call metadata_edio(nvar,igr,'patch level vapor transfer atm->canopy','[kg/kg]','ipatch') 
      end if

      if (associated(csite%cstar)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%cstar,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CSTAR :31:hist') 
         call metadata_edio(nvar,igr,'patch level co2 transfer atm->canopy','[ppm?]','ipatch') 
      end if

      if (associated(csite%zeta)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%zeta,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'ZETA :31:hist') 
         call metadata_edio(nvar,igr,'patch level height over Obukhov length','[---]','ipatch') 
      end if

      if (associated(csite%ribulk)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%ribulk,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RIBULK :31:hist') 
         call metadata_edio(nvar,igr,'patch level bulk Richardson number','[---]','ipatch') 
      end if
      
      if (associated(csite%upwp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%upwp,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'UPWP :31:hist') 
         call metadata_edio(nvar,igr,'Vertical flux of U-direction momentum','[kg m^-1 s^-2]','ipatch') 
      end if

      if (associated(csite%tpwp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,csite%tpwp,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'TPWP :31:hist') 
         call metadata_edio(nvar,igr,'Vertical flux of Heat','[kg K m^-2 s^-1]','ipatch')
      end if

      if (associated(csite%qpwp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,csite%qpwp,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'QPWP :31:hist') 
         call metadata_edio(nvar,igr,'Vertical flux of Moisture','[kg m^-2 s% -1]','ipatch')
      end if

      if (associated(csite%cpwp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,csite%cpwp,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'CPWP :31:hist') 
         call metadata_edio(nvar,igr,'Vertical flux of Carbon','[kg umolC mol^-1 m^-2 s-1]','ipatch')
      end if
      
      if (associated(csite%wpwp)) then
         nvar=nvar+1
         call vtable_edio_r(npts,csite%wpwp,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'WPWP :31:hist') 
         call metadata_edio(nvar,igr,'Vertical flux of W-direction momentum','[kg m^-1 s^-2]','ipatch')
      end if
      
      if (associated(csite%par_l_max)) then
         nvar=nvar+1
         call vtable_edio_r(npts,csite%par_l_max,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'PAR_L_MAX :31:hist') 
         call metadata_edio(nvar,igr,'Maximum PAR - not an output variable','[W/m2]','ipatch')
      end if
      
      if (associated(csite%par_l_beam_max)) then
         nvar=nvar+1
         call vtable_edio_r(npts,csite%par_l_beam_max,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'PAR_L_BEAM_MAX :31:hist') 
         call metadata_edio(nvar,igr,'Maximum direct PAR - not an output variable' &
                           ,'[W/m2]','ipatch')
      end if
      
      if (associated(csite%par_l_diffuse_max)) then
         nvar=nvar+1
         call vtable_edio_r(npts,csite%par_l_diffuse_max,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'PAR_L_DIFFUSE_MAX :31:hist') 
         call metadata_edio(nvar,igr,'Maximum diffuse PAR - not an output variable' &
                           ,'[W/m2]','ipatch')
      end if

      if (associated(csite%avg_rshort_gnd)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%avg_rshort_gnd,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'AVG_RSHORT_GND_PA :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%avg_rlong_gnd)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%avg_rlong_gnd,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'AVG_RLONG_GND_PA :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%avg_rlongup)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%avg_rlongup,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'AVG_RLONGUP_PA :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%avg_albedo)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%avg_albedo,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'AVG_ALBEDO_PA :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%avg_albedo_beam)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%avg_albedo_beam,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'AVG_ALBEDO_BEAM_PA :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%avg_albedo_diffuse)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%avg_albedo_diffuse,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'AVG_ALBEDO_DIFFUSE_PA :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%avg_rlong_albedo)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%avg_rlong_albedo,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'AVG_RLONG_ALBEDO_PA :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%avg_carbon_ac)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%avg_carbon_ac,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'AVG_CARBON_AC_PA :31:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if(associated(csite%dmean_co2_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,csite%dmean_co2_residual,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_CO2_RESIDUAL_PA :31:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean of residual canopy CO2 ','[umol/m2/day]','ipoly') 
      end if

      if(associated(csite%dmean_water_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,csite%dmean_water_residual,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_WATER_RESIDUAL_PA :31:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean of residual water ','[kg/m2/day]','ipoly') 
      end if

      if(associated(csite%dmean_energy_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,csite%dmean_energy_residual,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'DMEAN_ENERGY_RESIDUAL_PA :31:hist:dail') 
         call metadata_edio(nvar,igr,'Daily mean of residual water ','[J/m2/day]','ipoly') 
      end if

      if(associated(csite%mmean_co2_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,csite%mmean_co2_residual,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_CO2_RESIDUAL_PA :31:hist:dcyc:mont') 
         call metadata_edio(nvar,igr,'Monthly mean of residual canopy CO2 ','[umol/m2/day]','ipoly') 
      end if

      if(associated(csite%mmean_water_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,csite%mmean_water_residual,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_WATER_RESIDUAL_PA :31:hist:dcyc:mont') 
         call metadata_edio(nvar,igr,'Monthly mean of residual water ','[kg/m2/day]','ipoly') 
      end if

      if(associated(csite%mmean_energy_residual)) then
         nvar=nvar+1
         call vtable_edio_r(npts,csite%mmean_energy_residual,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'MMEAN_ENERGY_RESIDUAL_PA :31:hist:dcyc:mont') 
         call metadata_edio(nvar,igr,'Monthly mean of residual energy ','[J/m2/day]','ipoly') 
      end if

      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!






      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !       This part should have only 2-D vectors with dimensions npatches and ndcycle. !
      !  Notice that they all use the same npts.  Here you should only add variables of    !
      ! types -31.                                                                         !
      !------------------------------------------------------------------------------------!
      npts = csite%npatches * ndcycle

      if (associated(csite%qmean_rh)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%qmean_rh,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'QMEAN_RH_PA :-31:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%qmean_albedo)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%qmean_albedo,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'QMEAN_ALBEDO_PA :-31:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%qmean_albedo_beam)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%qmean_albedo_beam,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'QMEAN_ALBEDO_BEAM_PA :-31:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%qmean_albedo_diffuse)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%qmean_albedo_diffuse,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'QMEAN_ALBEDO_DIFFUSE_PA :-31:hist:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!






      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !       This part should have only 2-D vectors with dimensions npatches and nzg.     !
      !  Notice that they all use the same npts.  Here you should only add variables of    !
      ! types 32 and 320.                                                                  !
      !------------------------------------------------------------------------------------!
      npts = csite%npatches * nzg

      if (associated(csite%soil_energy)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%soil_energy,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'SOIL_ENERGY_PA :32:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%soil_water)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%soil_water,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'SOIL_WATER_PA :32:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%current_paw)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%current_paw,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CURRENT_PAW_PA :32:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%soil_tempk)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%soil_tempk,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'SOIL_TEMPK_PA :32:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[K]','ipatch : nzg') 
      end if

      if (associated(csite%soil_fracliq)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%soil_fracliq,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'SOIL_FRACLIQ_PA :32:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[normalized]','ipatch : nzg') 
      end if

      if (associated(csite%rootdense)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rootdense,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'PATCH_ROOT_DENSITY :32:hist') 
         call metadata_edio(nvar,igr,'Patch level root density with depth','[kg/m3]',&
              'ipatch : nzg') 
      end if

      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !       This part should have only 2-D vectors with dimensions npatches and nzs.     !
      !  Notice that they all use the same npts.  Here you should only add variables of    !
      ! type 33.                                                                           !
      !------------------------------------------------------------------------------------!
      npts = csite%npatches * nzs


      if (associated(csite%sfcwater_mass)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%sfcwater_mass,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'SFCWATER_MASS :33:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%sfcwater_energy)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%sfcwater_energy,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'SFCWATER_ENERGY :33:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%sfcwater_depth)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%sfcwater_depth,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'SFCWATER_DEPTH :33:hist:opti') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','m') 
      end if

      if (associated(csite%rshort_s)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rshort_s,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RSHORT_S :33:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%rshort_s_beam)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rshort_s_beam,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RSHORT_S_BEAM :33:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%rshort_s_diffuse)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%rshort_s_diffuse,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'RSHORT_S_DIFFUSE :33:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%sfcwater_tempk)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%sfcwater_tempk,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'SFCWATER_TEMPK :33:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%sfcwater_fracliq)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%sfcwater_fracliq,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'SFCWATER_FRACLIQ :33:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!






      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !       This part should have only 2-D vectors with dimensions npatches and n_pft.   !
      !  Notice that they all use the same npts.  Here you should only add variables of    !
      ! type 34.                                                                           !
      !------------------------------------------------------------------------------------!
      npts = csite%npatches * n_pft

      if (associated(csite%A_o_max)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%A_o_max,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'A_O_MAX :34:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%A_c_max)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%A_c_max,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'A_C_MAX :34:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if

      if (associated(csite%repro)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%repro,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'REPRO_PA :34:hist:mont:dcyc') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!






      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !       This part should have only 2-D vectors with dimensions npatches and n_dbh.   !
      !  Notice that they all use the same npts.  Here you should only add variables of    !
      ! type 36.                                                                           !
      !------------------------------------------------------------------------------------!
      npts = csite%npatches * n_dbh

      if (associated(csite%co2budget_gpp)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%co2budget_gpp,nvar,igr,init,csite%paglob_id, &
           var_len,var_len_global,max_ptrs,'CO2BUDGET_GPP_DBH :36:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!






      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !       This part should have only 3-D vectors with dimensions npatches and n_dbh.   !
      !  Notice that they all use the same npts.  Here you should only add variables of    !
      ! type 346.                                                                          !
      !------------------------------------------------------------------------------------!
      npts = csite%npatches * n_pft * n_dbh

      if (associated(csite%cumlai_profile)) then
         nvar=nvar+1
           call vtable_edio_r(npts,csite%cumlai_profile,nvar,igr,init,csite%paglob_id,     &
           var_len,var_len_global,max_ptrs,'CUMLAI_PROFILE :346:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!






      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !       This part should have only 3-D vectors with dimensions n_stoma_atts and      !
      ! n_dbh.  Notice that they all use the same npts.  Here you should only add vari-    !
      ! ables of type 316.                                                                 !
      !------------------------------------------------------------------------------------!
      npts = csite%npatches * n_stoma_atts * n_pft

      if (associated(csite%old_stoma_vector_max)) then
         nvar=nvar+1
         call vtable_edio_r(npts,csite%old_stoma_vector_max,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'OLD_STOMA_VECTOR_MAX :316:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      !------------------------------------------------------------------------------------!
      npts = csite%npatches * nzg * 10

      if (associated(csite%past_paw)) then
         nvar=nvar+1
         call vtable_edio_r(npts,csite%past_paw,nvar,igr,init,csite%paglob_id, &
              var_len,var_len_global,max_ptrs,'PAST_PAW_PA :317:hist') 
         call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
      end if
      !------------------------------------------------------------------------------------!



      !----- Save the number of patch-level (sitetype) variables that go to the output. ---!
      if (init == 0) niosite=nvar-niopoly-niogrid-nioglobal
      !------------------------------------------------------------------------------------!

      return
      
   end subroutine filltab_sitetype
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine filltab_patchtype(igr,ipy,isi,ipa,init)

      use ed_var_tables, only : vtable_edio_r & ! sub-routine
                              , vtable_edio_i & ! sub-routine
                              , metadata_edio ! ! sub-routine

      use cohort_state, only: filltab_patchtype_state
      use cohort_phen, only: filltab_patchtype_phen
      use cohort_mort, only: filltab_patchtype_mort
      use cohort_resp, only: filltab_patchtype_resp
      use cohort_photo, only: filltab_patchtype_photo
      use cohort_rad, only: filltab_patchtype_rad
      use cohort_therm, only: filltab_patchtype_therm

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer          , intent(in) :: init
      integer          , intent(in) :: igr
      integer          , intent(in) :: ipy
      integer          , intent(in) :: isi
      integer          , intent(in) :: ipa
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer    :: cpoly
      type(sitetype)   , pointer    :: csite
      type(patchtype)  , pointer    :: cpatch
      integer                       :: ip
      integer                       :: np
      integer                       :: is
      integer                       :: ns
      integer                       :: iy
      integer                       :: ny
      integer                       :: var_len
      integer                       :: max_ptrs
      integer                       :: var_len_global
      integer                       :: nvar
      integer                       :: npts
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      The condition exists where some patches may have no cohorts, particularly     !
      ! over deserts.  Therefore it is possible that the first patch in the site has no    !
      ! cohorts, but the others might have some cohorts.  Therefore, we want to make sure  !
      ! we create mapping tables if "any" of the patches have cohorts, so we will provide  !
      ! for this conditions when we make the mapping tables.  The condition has to arise   !
      ! where there are no cohorts at all in any of the patches in any of the polygons on  !
      ! the current machine to not initialize these vectors.                               !
      !------------------------------------------------------------------------------------!
      if (init==0) then
         ny = edgrid_g(igr)%npolygons

         polyloop: do iy=1,ny
            cpoly => edgrid_g(igr)%polygon(iy)
            ns = cpoly%nsites
            siteloop: do is=1,ns
               csite => cpoly%site(is)
               np = csite%npatches
               patchloop: do ip=1,np
                  cpatch => edgrid_g(igr)%polygon(iy)%site(is)%patch(ip)
                  if (cpatch%ncohorts > 0) exit polyloop
               end do patchloop
            end do siteloop
         end do polyloop

      else
         cpatch => edgrid_g(igr)%polygon(ipy)%site(isi)%patch(ipa)
      end if
      !------------------------------------------------------------------------------------!


      !----- Nothing else to do here in case this is an empty patch. ----------------------!
      if (cpatch%ncohorts == 0) return
      !------------------------------------------------------------------------------------!


      !----- Define some variables to be used by all variables. ---------------------------!
      var_len        = cpatch%ncohorts
      var_len_global = edgrid_g(igr)%ncohorts_global
      max_ptrs       = edgrid_g(igr)%npatches_global
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Initialise the variable number, which should be unique for each variable.      !
      ! nioglobal has the current position after all the parent structure variables have   !
      ! been positioned in the variable table.                                             !
      !------------------------------------------------------------------------------------!
      nvar=nioglobal+niogrid+niopoly+niosite
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !       This part should have only 1-D vectors, with dimension ncohorts.  Notice     !
      ! that they all use the same npts.  Here you should only add variables of type 40    !
      ! and 41.                                                                            !
      !------------------------------------------------------------------------------------!
      npts = cpatch%ncohorts

      call filltab_patchtype_state(nvar, npts, cpatch%costate, igr, init, cpatch%coglob_id, &
           var_len, var_len_global, max_ptrs)
      call filltab_patchtype_phen(nvar, npts, cpatch%cophen, igr, init, cpatch%coglob_id, &
           var_len, var_len_global, max_ptrs)
      call filltab_patchtype_mort(nvar, npts, cpatch%comort, igr, init, cpatch%coglob_id, &
           var_len, var_len_global, max_ptrs, n_mort)
      call filltab_patchtype_resp(nvar, npts, cpatch%coresp, igr, init, cpatch%coglob_id, &
           var_len, var_len_global, max_ptrs, ndcycle)
      call filltab_patchtype_photo(nvar, npts, cpatch%cophoto, igr, init, cpatch%coglob_id, &
           var_len, var_len_global, max_ptrs, ndcycle)
      call filltab_patchtype_rad(nvar, npts, cpatch%corad, igr, init, cpatch%coglob_id, &
           var_len, var_len_global, max_ptrs)
      call filltab_patchtype_therm(nvar, npts, cpatch%cotherm, igr, init, cpatch%coglob_id, &
           var_len, var_len_global, max_ptrs)


      return
   end subroutine filltab_patchtype
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
                   ! ========= UTILITITY FUNCTIONS =========== !
!==========================================================================================!
!==========================================================================================!
  function get_nsites(cgrid)

    implicit none
    integer :: get_nsites
    integer :: ipy,isi
    type(edtype),target           :: cgrid
    type(polygontype),pointer     :: cpoly

    get_nsites = 0

    do ipy=1,cgrid%npolygons
       cpoly=>cgrid%polygon(ipy)
       do isi=1,cpoly%nsites
          get_nsites = get_nsites + 1
       end do
    end do

    return
  end function get_nsites
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  function get_npatches(cgrid)

    implicit none
    integer :: get_npatches
    integer :: ipy,isi,ipa
    type(edtype),target           :: cgrid
    type(polygontype),pointer     :: cpoly
    type(sitetype),pointer        :: csite

    get_npatches = 0

    do ipy=1,cgrid%npolygons
       cpoly=>cgrid%polygon(ipy)
       do isi=1,cpoly%nsites
          csite=>cpoly%site(isi)
          do ipa=1,csite%npatches
             get_npatches = get_npatches + 1
          end do
       end do
    end do
    return
  end function get_npatches
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  function get_ncohorts(cgrid)

    implicit none
    integer :: get_ncohorts
    integer :: ipy,isi,ipa,ico
    type(edtype),target           :: cgrid
    type(polygontype),pointer     :: cpoly
    type(sitetype),pointer        :: csite
    type(patchtype),pointer       :: cpatch

    get_ncohorts = 0

    do ipy=1,cgrid%npolygons
       cpoly=>cgrid%polygon(ipy)
       do isi=1,cpoly%nsites
          csite=>cpoly%site(isi)
          do ipa=1,csite%npatches
             cpatch=>csite%patch(ipa)
             get_ncohorts = get_ncohorts + cpatch%ncohorts
          end do
       end do
    end do
    return
  end function get_ncohorts
!==========================================================================================!
!==========================================================================================!
end module ed_state_vars




