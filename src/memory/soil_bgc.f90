! Soil biogeochemistry module for ED2.

Module soil_bgc
  implicit none

  !====================
  ! Soil Carbon
  !====================

  ! Soil carbon, fast pool (kgC/m2)
  real, dimension(:) :: fast_soil_C 

  ! Soil carbon, slow pool (kgC/m2)
  real , dimension(:) :: slow_soil_C 
  
  ! Soil carbon, structural pool (kgC/m2)
  real , dimension(:) :: structural_soil_C
  
  !---------------
  ! This one is a legacy variable.  Probably can delete this.
  ! Soil lignin , structural pool (kgC/m2)
  real , dimension(:) :: structural_soil_L
  !---------------
  
  !====================
  ! Soil Carbon
  !====================

  ! Soil nitrogen, mineralized pool (kgN/m2)
  real , dimension(:) :: mineralized_soil_N  
  
  
  !====================
  ! Soil Carbon
  !====================

  ! Soil phosphorus , fast pool (kgP/m2)
  real , dimension(:) :: fast_soil_P  
  ! Soil phosphorus , structural pool (kgP/m2)
  real , dimension(:) :: struct_soil_P  
  ! Soil phosphorus , mineralized pool (kgP/m2)
  real , dimension(:) :: miner_soil_P 
  ! Soil phosphorus , slow pool (kgP/m2)
  real , dimension(:) :: slow_soil_P 
  







     ! Soil nitrogen , fast pool (kg/m2)
     real , dimension(:) :: fast_soil_N 



end Module soil_bgc

