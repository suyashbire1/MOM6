!> Initialization of the boundary-forced-basing configuration
module BFB_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type
use MOM_variables, only : thermo_var_ptrs
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use MOM_verticalGrid, only : verticalGrid_type
implicit none ; private

#include <MOM_memory.h>

public BFB_set_coord
public BFB_initialize_thickness, BFB_initialize_thickness_varlayth
public BFB_initialize_sponges_southonly
public BFB_initialize_topography

!> Unsafe model variable
!! \todo Remove this module variable
logical :: first_call = .true.

contains

!> Initialize topography.
subroutine BFB_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m
  real                                            :: efold, rct, ebdepth
  real                                            :: westlon,eastlon, lenlon
  real, parameter                     :: piby180 = 4.0*atan(1.0)/180.0
  character(len=40)  :: mdl = "BFB_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

!!$  call MOM_error(FATAL, &
!!$   "USER_initialization.F90, USER_initialize_topography: " // &
!!$   "Unmodified user routine called - you must edit the routine to use it")

  ! Slope at the eastern boundary, default is 2m/km
  call get_param(param_file, mdl, "TOPO_EFOLDING_SCALE", efold, &
          "E-folding scale of the topography at the eastern boundary on the continental shelf", units="degrees", default=2.5)
  call get_param(param_file, mdl, "DEPTH_EB", ebdepth, &
       "Depth at the eastern boundary", &
       units="m", default=100.0)
  call get_param(param_file, mdl, "LENLON", lenlon, &
                 "The longitudinal length of the domain.", units="degrees")
  call get_param(param_file, mdl, "WESTLON", westlon, &
                 "The western longitude of the domain.", units="degrees", default=0.0)

  eastlon = westlon + lenlon
  do j=js,je; do i=is,ie
      D(i,j) = min(max_depth,exp(-efold*G%geoLonT(i,j)))
  enddo; enddo

  if (first_call) call write_BFB_log(param_file)

end subroutine BFB_initialize_topography

!> This subroutine sets up the sponges for the southern bouundary of the domain. Maximum damping occurs
!! within 2 degrees lat of the boundary. The damping linearly decreases northward over the next 2 degrees.
subroutine BFB_initialize_sponges_southonly(G, use_temperature, tv, param_file, CSp, h)
  type(ocean_grid_type), intent(in) :: G    !< The ocean's grid structure
  logical,               intent(in) :: use_temperature !< If true, temperature and salinity are used as
                                            !! state variables.
  type(thermo_var_ptrs), intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  type(sponge_CS),       pointer    :: CSp  !< A pointer to the sponge control structure
  real, dimension(NIMEM_, NJMEM_, NKMEM_), &
                         intent(in) :: h    !< Layer thicknesses, in H (usually m or kg m-2)
  !call MOM_error(FATAL, &
  ! "BFB_initialization.F90, BFB_initialize_sponges: " // &
  ! "Unmodified user routine called - you must edit the routine to use it")

  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! A temporary array for eta.
  real :: Idamp(SZI_(G),SZJ_(G))    ! The inverse damping rate, in s-1.

  real :: H0(SZK_(G))
  real :: min_depth, D_aby
  real :: damp, e_dense, slat, wlon, lenlat, lenlon, nlat
  character(len=40)  :: mdl = "BFB_initialize_sponges_southonly" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  eta(:,:,:) = 0.0 ; Idamp(:,:) = 0.0

!  Here the inverse damping time, in s-1, is set. Set Idamp to 0     !
!  wherever there is no sponge, and the subroutines that are called  !
!  will automatically set up the sponges only where Idamp is positive!
!  and mask2dT is 1.                                                   !

!   Set up sponges for DOME configuration
  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)

  call get_param(param_file, mdl, "SOUTHLAT", slat, &
                 "The southern latitude of the domain.", units="degrees")
  call get_param(param_file, mdl, "LENLAT", lenlat, &
                 "The latitudinal length of the domain.", units="degrees")
  call get_param(param_file, mdl, "WESTLON", wlon, &
                 "The western longitude of the domain.", units="degrees", default=0.0)
  call get_param(param_file, mdl, "LENLON", lenlon, &
                 "The longitudinal length of the domain.", units="degrees")
  nlat = slat + lenlat
  call get_param(param_file, mdl, "D_ABYSS", D_aby, &
                 "Depth at which abyssal layer starts", units="m", default=1500.0)
!  do k=1,nz ; H0(k) = -G%max_depth * real(k-1) / real(nz) ; enddo
!  do k=1,nz ; H0(k) = -G%max_depth * real(k-1) / real(nz-1) ; enddo ! Use for meridional thickness profile initialization
  do k=1,nz ; H0(k) = -D_aby * real(k-1) / real(nz-1) ; enddo

  ! Use for meridional thickness profile initialization
!  do k=1,nz ; H0(k) = -G%max_depth * real(k-1) / real(nz-1) ; enddo
  do i=is,ie; do j=js,je
    if (G%geoLatT(i,j) < slat+2.0) then ; damp = 1.0
    elseif (G%geoLatT(i,j) < slat+4.0) then
       damp = 1.0*(slat+4.0-G%geoLatT(i,j))/2.0
    else ; damp = 0.0
    endif

    ! These will be streched inside of apply_sponge, so they can be in
    ! depth space for Boussinesq or non-Boussinesq models.

    ! This section is used for uniform thickness initialization
    ! except at the topography
    do k = 1,nz
       if (H0(k) < -G%bathyT(i,j)) then
          eta(i,j,k) = -G%bathyT(i,j)
       else
          eta(i,j,k) = H0(k)
       endif
    enddo

    ! The below section is used for meridional temperature profile thickness initiation
    ! do k = 1,nz; eta(i,j,k) = H0(k); enddo
    ! if (G%geoLatT(i,j) > 40.0) then
    !   do k = 1,nz
    !     eta(i,j,k) = -G%Angstrom_z*(k-1)
    !   enddo
    ! elseif (G%geoLatT(i,j) > 20.0) then
    !   do k = 1,nz
    !     eta(i,j,k) = min(H0(k) + (G%geoLatT(i,j) - 20.0)*(G%max_depth - nz*G%Angstrom_z)/20.0, -(k-1)*G%angstrom_z)
    !   enddo
    ! endif
    eta(i,j,nz+1) = -G%bathyT(i,j)

    if (G%bathyT(i,j) > min_depth) then
      Idamp(i,j) = damp/86400.0
    else ; Idamp(i,j) = 0.0 ; endif
  enddo ; enddo

!  This call sets up the damping rates and interface heights.
!  This sets the inverse damping timescale fields in the sponges.    !
  call initialize_sponge(Idamp, eta, G, param_file, CSp)

!   Now register all of the fields which are damped in the sponge.   !
! By default, momentum is advected vertically within the sponge, but !
! momentum is typically not damped within the sponge.                !

  if (first_call) call write_BFB_log(param_file)

end subroutine BFB_initialize_sponges_southonly

!> This subroutine specifies the vertical coordinate in terms of temperature at the surface and at the bottom.
!! This case is set up in such a way that the temperature of the topmost layer is equal to the SST at the
!! southern edge of the domain. The temperatures are then converted to densities of the top and bottom layers
!! and linearly interpolated for the intermediate layers.
subroutine BFB_set_coord(Rlay, g_prime, GV, param_file, eqn_of_state)
  real, dimension(NKMEM_), intent(out) :: Rlay !< Layer potential density.
  real, dimension(NKMEM_), intent(out) :: g_prime !< The reduced gravity at
                                                  !! each interface, in m2 Z-1 s-2.
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  type(param_file_type),   intent(in)  :: param_file !< A structure to parse for run-time parameters
  type(EOS_type),          pointer     :: eqn_of_state !< Integer that selects the
                                                     !! equation of state.
  ! Local variables
  real                                 :: drho_dt, SST_s, T_bot, rho_top, rho_bot
  integer                              :: k, nz
  character(len=40)  :: mdl = "BFB_set_coord" ! This subroutine's name.

  call get_param(param_file, mdl, "DRHO_DT", drho_dt, &
          "Rate of change of density with temperature.", &
           units="kg m-3 K-1", default=-0.2)
  call get_param(param_file, mdl, "SST_S", SST_s, &
          "SST at the suothern edge of the domain.", units="C", default=20.0)
  call get_param(param_file, mdl, "T_BOT", T_bot, &
                 "Bottom Temp", units="C", default=5.0)

  rho_top = GV%rho0 + drho_dt*(SST_s-T_bot)
  rho_bot = GV%rho0 + drho_dt*(T_bot-T_bot)
  nz = GV%ke

  do k = 1,nz
    Rlay(k) = (rho_bot - rho_top)/(nz-1)*real(k-1) + rho_top
    if (k >1) then
      g_prime(k) = (Rlay(k) - Rlay(k-1)) * GV%g_Earth/GV%rho0
    else
      g_prime(k) = GV%g_Earth
    endif
    !Rlay(:) = 0.0
    !g_prime(:) = 0.0
  enddo

  if (first_call) call write_BFB_log(param_file)

end subroutine BFB_set_coord

subroutine BFB_initialize_thickness(h, G, param_file)
  real, intent(out), dimension(NIMEM_,NJMEM_, NKMEM_) :: h
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! A temporary array for eta.
  real :: H0(SZK_(G))
  real :: D_aby
  character(len=40)  :: mod = "BFB_initialize_thickness" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call get_param(param_file, mod, "D_ABYSS", D_aby, &
                 "Depth at which abyssal layer starts", units="m", default=1500.0)

  eta(:,:,:) = 0.0
  do k=1,nz ; H0(k) = -D_aby * real(k-1) / real(nz-1) ; enddo
  do i=is,ie; do j=js,je
    do k = 1,nz
       if (H0(k) < -G%bathyT(i,j)) then
          eta(i,j,k) = -G%bathyT(i,j)
       else
          eta(i,j,k) = H0(k)
       endif
    enddo
    eta(i,j,nz+1) = -G%bathyT(i,j)
    do k = 1,nz; h(i,j,k) = eta(i,j,k) - eta(i,j,k+1); enddo
  enddo; enddo

 if (first_call) call write_BFB_log(param_file)

end subroutine BFB_initialize_thickness

subroutine BFB_initialize_thickness_varlayth(h, G, param_file)
  real, intent(out), dimension(NIMEM_,NJMEM_, NKMEM_) :: h
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! A temporary array for eta.
  real :: H0(SZK_(G))
  real :: hmin, hmax
  real, parameter :: pi = 4.0*atan(1.0)
  character(len=40)  :: mod = "BFB_initialize_thickness_varlayth" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call get_param(param_file, mod, "HMIN", hmin, &
                 "Initial thickness of thinnest layer", units="m")
  eta(:,:,:) = 0.0
  hmax = (pi*G%max_depth - real(nz)*hmin*(pi-2.0))/2/real(nz)
  do k=1,nz ; H0(k) = 2*nz*(hmin-hmax)/pi*(cos((k-1)*pi/2/nz)-1.0) + (k-1)*hmin; enddo
  do i=is,ie; do j=js,je
   do k = 1,nz; eta(i,j,k) = -H0(k); enddo
    eta(i,j,nz+1) = -G%max_depth
    do k = 1,nz; h(i,j,k) = eta(i,j,k) - eta(i,j,k+1); enddo
  enddo; enddo

 if (first_call) call write_BFB_log(param_file)

end subroutine BFB_initialize_thickness_varlayth


!> Write output about the parameter values being used.
subroutine write_BFB_log(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure indicating the
                                                  !! open file to parse for model
                                                  !! parameter values.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "BFB_initialization" ! This module's name.

  call log_version(param_file, mdl, version)
  first_call = .false.

end subroutine write_BFB_log

end module BFB_initialization
