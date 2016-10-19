module MOM_diagnostics

!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, February 2001                                  *
!*                                                                     *
!*    This subroutine calculates any requested diagnostic quantities   *
!*  that are not calculated in the various subroutines.  Diagnostic    *
!*  quantities are requested by allocating them memory.                *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, CoriolisBu                            *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, bathyT                                *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_coms,              only : reproducing_sum
use MOM_diag_mediator,     only : post_data, post_data_1d_k
use MOM_diag_mediator,     only : register_diag_field, register_scalar_field
use MOM_diag_mediator,     only : diag_ctrl, time_type, safe_alloc_ptr
use MOM_domains,           only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,           only : To_North, To_East
use MOM_EOS,               only : calculate_density, int_density_dz
use MOM_error_handler,     only : MOM_error, FATAL, WARNING
use MOM_file_parser,       only : get_param, log_version, param_file_type
use MOM_forcing_type,      only : forcing
use MOM_grid,              only : ocean_grid_type
use MOM_interface_heights, only : find_eta
use MOM_spatial_means,     only : global_area_mean, global_layer_mean, global_volume_mean
use MOM_variables,         only : thermo_var_ptrs, ocean_internal_state, p3d
use MOM_variables,         only : accel_diag_ptrs, cont_diag_ptrs, BT_cont_type
use MOM_verticalGrid,      only : verticalGrid_type
use MOM_wave_speed,        only : wave_speed, wave_speed_CS

implicit none ; private

#include <MOM_memory.h>

public calculate_diagnostic_fields
public register_time_deriv
public find_eta
public MOM_diagnostics_init
public MOM_diagnostics_end


type, public :: diagnostics_CS ; private
  logical :: split                 ! If true, use split time stepping scheme.
  type(diag_ctrl), pointer :: diag ! structure to regulate diagnostics timing

  ! following arrays store diagnostics calculated here and unavailable outside.

  ! following fields have nz+1 levels.
  real, pointer, dimension(:,:,:) :: &
    e => NULL(), &   ! interface height (metre)
    e_D => NULL()    ! interface height above bottom (metre)

  ! following fields have nz layers.
  real, pointer, dimension(:,:,:) :: &
    du_dt => NULL(), &    ! net i-acceleration in m/s2
    dv_dt => NULL(), &    ! net j-acceleration in m/s2
    dh_dt => NULL(), &    ! thickness rate of change in (m/s) or kg/(m2*s)

    h_Rlay => NULL(),    & ! layer thicknesses in layered potential density
                           ! coordinates, in m (Bouss) or kg/m2 (non-Bouss)
    uh_Rlay   => NULL(), & ! zonal and meridional transports in layered 
    vh_Rlay   => NULL(), & ! potential rho coordinates: m3/s(Bouss) kg/s(non-Bouss)
    uhGM_Rlay => NULL(), & ! zonal and meridional Gent-McWilliams transports in layered 
    vhGM_Rlay => NULL()    ! potential density coordinates, m3/s (Bouss) kg/s(non-Bouss)

  ! following fields are 2-D.
  real, pointer, dimension(:,:) :: &
    cg1 => NULL(),       & ! first baroclinic gravity wave speed, in m s-1
    Rd1 => NULL(),       & ! first baroclinic deformation radius, in m
    cfl_cg1 => NULL(),   & ! CFL for first baroclinic gravity wave speed, nondim
    cfl_cg1_x => NULL(), & ! i-component of CFL for first baroclinic gravity wave speed, nondim
    cfl_cg1_y => NULL()    ! j-component of CFL for first baroclinic gravity wave speed, nondim

  ! arrays to hold diagnostics in the layer-integrated energy budget.
  ! all except KE have units of m3 s-3 (when Boussinesq).
  real, pointer, dimension(:,:,:) :: &
    KE        => NULL(), &  ! KE per unit mass, in m2 s-2
    dKE_dt    => NULL(), &  ! time derivative of the layer KE
    PE_to_KE  => NULL(), &  ! potential energy to KE term
    KE_CorAdv => NULL(), &  ! KE source from the combined Coriolis and
                            ! advection terms.  The Coriolis source should be
                            ! zero, but is not due to truncation errors.  There
                            ! should be near-cancellation of the global integral
                            ! of this spurious Coriolis source.
    KE_adv     => NULL(),&  ! KE source from along-layer advection
    KE_visc    => NULL(),&  ! KE source from vertical viscosity
    KE_horvisc => NULL(),&  ! KE source from horizontal viscosity
    KE_dia     => NULL(),&  ! KE source from diapycnal diffusion
    diag_tmp3d => NULL()    ! 3D re-usable arrays for diagnostics

  real, pointer, dimension(:,:,:) :: &
    hfv         => NULL(),&  ! coriolis xtwa term
    hpfu        => NULL(),&  ! PG xtwa term
    huwb        => NULL(),&  ! diabatic xtwa term
    huuxpt      => NULL(),&  ! zonal advection xtwa term
    huvymt      => NULL(),&  ! meridional advection xtwa term
    hdudtvisc   => NULL(),&  ! vertical viscosity xtwa term
    hdiffu      => NULL(),&  ! horizontal viscosity xtwa term

    hmfu        => NULL(),&  ! coriolis ytwa term
    hpfv        => NULL(),&  ! PG ytwa term
    hvwb        => NULL(),&  ! diabatic ytwa term
    huvxpt      => NULL(),&  ! zonal advection ytwa term
    hvvymt      => NULL(),&  ! meridional advection ytwa term
    hdvdtvisc   => NULL(),&  ! vertical viscosity ytwa term
    hdiffv      => NULL()    ! horizontal viscosity ytwa term

  real, pointer, dimension(:,:,:) :: &
    h_Cu         => NULL(),&
    huu_Cu       => NULL(),&
    hv_Cu        => NULL(),&
    hw_Cu        => NULL(),&
    hwb_Cu       => NULL(),&
    pfu_masked   => NULL(),&
    esq_Cu       => NULL(),&
    e_Cu         => NULL(),&
    epfu         => NULL(),&

    h_Cv         => NULL(),&
    hvv_Cv       => NULL(),&
    hu_Cv        => NULL(),&
    hw_Cv        => NULL(),&
    hwb_Cv       => NULL(),&
    pfv_masked   => NULL(),&
    esq_Cv       => NULL(),&
    e_Cv         => NULL(),&
    epfv         => NULL()

  ! diagnostic IDs
  integer :: id_e              = -1, id_e_D            = -1
  integer :: id_du_dt          = -1, id_dv_dt          = -1
  integer :: id_col_ht         = -1, id_dh_dt          = -1
  integer :: id_KE             = -1, id_dKEdt          = -1
  integer :: id_PE_to_KE       = -1, id_KE_Coradv      = -1
  integer :: id_KE_adv         = -1, id_KE_visc        = -1
  integer :: id_KE_horvisc     = -1, id_KE_dia         = -1
  integer :: id_uh_Rlay        = -1, id_vh_Rlay        = -1
  integer :: id_uhGM_Rlay      = -1, id_vhGM_Rlay      = -1
  integer :: id_h_Rlay         = -1, id_Rd1            = -1     
  integer :: id_Rml            = -1, id_Rcv            = -1
  integer :: id_cg1            = -1, id_cfl_cg1        = -1
  integer :: id_cfl_cg1_x      = -1, id_cfl_cg1_y      = -1
  integer :: id_temp_int       = -1, id_salt_int       = -1
  integer :: id_mass_wt        = -1, id_col_mass       = -1
  integer :: id_masscello      = -1, id_masso          = -1
  integer :: id_thetaoga       = -1, id_soga           = -1
  integer :: id_sosga          = -1, id_tosga          = -1
  integer :: id_temp_layer_ave = -1, id_salt_layer_ave = -1
  integer :: id_pbo            = -1
  integer :: id_thkcello       = -1, id_rhoinsitu      = -1
  integer :: id_rhopot0        = -1, id_rhopot2        = -1

  integer :: id_hfv            = -1, id_hmfu           = -1
  integer :: id_hpfu           = -1, id_hpfv           = -1
  integer :: id_huwb           = -1, id_hvwb           = -1
  integer :: id_huuxpt         = -1, id_huvxpt         = -1
  integer :: id_huvymt         = -1, id_hvvymt         = -1
  integer :: id_hdudtvisc      = -1, id_hdvdtvisc      = -1
  integer :: id_hdiffu         = -1, id_hdiffv         = -1

  integer :: id_h_Cu           = -1, id_h_Cv           = -1
  integer :: id_huu_Cu         = -1, id_hvv_Cv         = -1
  integer :: id_hv_Cu          = -1, id_hu_Cv          = -1
  integer :: id_hw_Cu          = -1, id_hw_Cv          = -1
  integer :: id_hwb_Cu         = -1, id_hwb_Cv         = -1
  integer :: id_epfu           = -1, id_epfv           = -1
  integer :: id_e_Cu           = -1, id_e_Cv           = -1
  integer :: id_esq_Cu         = -1, id_esq_Cv         = -1
  integer :: id_pfu_masked     = -1, id_pfv_masked     = -1

  type(wave_speed_CS), pointer :: wave_speed_CSp => NULL()  

  ! pointers used in calculation of time derivatives
  type(p3d) :: var_ptr(MAX_FIELDS_)
  type(p3d) :: deriv(MAX_FIELDS_)
  type(p3d) :: prev_val(MAX_FIELDS_)
  integer   :: nlay(MAX_FIELDS_)
  integer   :: num_time_deriv = 0
  
  ! for group halo pass
  type(group_pass_type) :: pass_KE_uv 

end type diagnostics_CS

contains

subroutine calculate_diagnostic_fields(u, v, h, uh, vh, tv, ADp, CDp, fluxes, &
                                       dt, G, GV, CS, eta_bt)
  type(ocean_grid_type),                     intent(inout) :: G
  type(verticalGrid_type),                   intent(in)    :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: uh
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: vh
  type(thermo_var_ptrs),                     intent(in)    :: tv
  type(accel_diag_ptrs),                     intent(in)    :: ADp
  type(cont_diag_ptrs),                      intent(in)    :: CDp
  type(forcing),                             intent(in)    :: fluxes
  real,                                      intent(in)    :: dt
  type(diagnostics_CS),                      intent(inout) :: CS
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)    :: eta_bt

! Diagnostics not more naturally calculated elsewhere are computed here. 

! Arguments: 
!  (in)      u   - zonal velocity component (m/s)
!  (in)      v   - meridional velocity component (m/s)
!  (in)      h   - layer thickness,  meter(Bouss)  kg/m^2(non-Bouss)
!  (in)      uh  - transport through zonal faces = u*h*dy, m3/s(Bouss) kg/s(non-Bouss)
!  (in)      vh  - transport through meridional faces = v*h*dx, m3/s(Bouss) kg/s(non-Bouss)
!  (in)      tv  - structure pointing to various thermodynamic variables.
!  (in)      ADp - structure with pointers to accelerations in momentum equation
!  (in)      CDp - structure with pointers to terms in continuity equation
!  (in)      dt  - time difference in s since the last call to this subroutine
!  (in)      G   - ocean grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS  - control structure returned by a previous call to diagnostics_init
!  (in,opt)  eta_bt - An optional barotropic variable that gives the "correct"
!                     free surface height (Boussinesq) or total water column
!                     mass per unit area (non-Boussinesq).  This is used to
!                     dilate the layer thicknesses when calculating interface
!                     heights, in m or kg m-2.

  integer i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb

  ! coordinate variable potential density, in kg m-3.
  real :: Rcv(SZI_(G),SZJ_(G),SZK_(G)) 
                                   
  ! tmp array for surface properties 
  real :: surface_field(SZI_(G),SZJ_(G)) 
  real :: pressure_1d(SZI_(G)) ! Temporary array for pressure when calling EOS
  real :: wt, wt_p

  ! squared Coriolis parameter at to h-points (1/s2)
  real :: f2_h      

  ! magnitude of the gradient of f (1/(m*s))
  real :: mag_beta
 
  ! frequency squared used to avoid division by 0 (1/s2)
  ! value is roughly (pi / (the age of the universe) )^2.
  real, parameter :: absurdly_small_freq2 = 1e-34 
    
  integer :: k_list

  real, dimension(SZK_(G)) :: temp_layer_ave, salt_layer_ave
  real :: thetaoga, soga, masso, tosga, sosga

  is  = G%isc  ; ie   = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq  = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  nz  = G%ke   ; nkmb = GV%nk_rho_varies

  ! smg: is the following robust to ALE? It seems a bit opaque.   
  ! If the model is NOT in isopycnal mode then nkmb=0. But we need all the
  ! following diagnostics to treat all layers as variable density, so we set
  ! nkmb = nz, on the expectation that loops nkmb+1,nz will not iterate.
  ! This behavior is ANSI F77 but some compiler options can force at least
  ! one iteration that would break the following one-line workaround!
  if (nkmb==0) nkmb = nz

  if (loc(CS)==0) call MOM_error(FATAL, &
         "calculate_diagnostic_fields: Module must be initialized before used.")

  call calculate_derivs(dt, G, CS)

  if (ASSOCIATED(CS%e)) then
    call find_eta(h, tv, GV%g_Earth, G, GV, CS%e, eta_bt)
    if (CS%id_e > 0) call post_data(CS%id_e, CS%e, CS%diag)
  endif

  if (ASSOCIATED(CS%e_D)) then
    if (ASSOCIATED(CS%e)) then
      do k=1,nz+1 ; do j=js,je ; do i=is,ie
        CS%e_D(i,j,k) = CS%e(i,j,k) + G%bathyT(i,j)
      enddo ; enddo ; enddo
    else
      call find_eta(h, tv, GV%g_Earth, G, GV, CS%e_D, eta_bt)
      do k=1,nz+1 ; do j=js,je ; do i=is,ie
        CS%e_D(i,j,k) = CS%e_D(i,j,k) + G%bathyT(i,j)
      enddo ; enddo ; enddo
    endif

    if (CS%id_e_D > 0) call post_data(CS%id_e_D, CS%e_D, CS%diag)
  endif

  ! mass per area of grid cell (for Bouss, use Rho0)
  if (CS%id_masscello > 0) then
    do k=1,nz; do j=js,je ; do i=is,ie 
       CS%diag_tmp3d(i,j,k) = GV%H_to_kg_m2*h(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_masscello, CS%diag_tmp3d, CS%diag)
  endif

  ! mass of liquid ocean (for Bouss, use Rho0)
  if (CS%id_masso > 0) then
    do k=1,nz; do j=js,je ; do i=is,ie
       CS%diag_tmp3d(i,j,k) = GV%H_to_kg_m2*h(i,j,k)*G%areaT(i,j)
    enddo ; enddo ; enddo
    masso = (reproducing_sum(sum(CS%diag_tmp3d,3)))
    call post_data(CS%id_masso, masso, CS%diag)
  endif

  ! diagnose thickness of grid cells (meter)
  if (CS%id_thkcello > 0) then

    ! thkcello = h for Boussinesq 
    if (GV%Boussinesq) then 
      call post_data(CS%id_thkcello, GV%H_to_m*h, CS%diag)

    ! thkcello = dp/(rho*g) for non-Boussinesq
    else 
      do j=js,je

        ! pressure loading at top of surface layer (Pa) 
        if(ASSOCIATED(fluxes%p_surf)) then
          do i=is,ie
            pressure_1d(i) = fluxes%p_surf(i,j) 
          enddo
        else
          do i=is,ie
            pressure_1d(i) = 0.0
          enddo
        endif 

        do k=1,nz
          ! pressure for EOS at the layer center (Pa)
          do i=is,ie
            pressure_1d(i) = pressure_1d(i) + 0.5*(GV%g_Earth*GV%H_to_kg_m2)*h(i,j,k)
          enddo
          ! store in-situ density (kg/m3) in diag_tmp3d
          call calculate_density(tv%T(:,j,k),tv%S(:,j,k), pressure_1d, &
                                 CS%diag_tmp3d(:,j,k), is, ie-is+1, tv%eqn_of_state)
          ! cell thickness = dz = dp/(g*rho) (meter); store in diag_tmp3d
          do i=is,ie
            CS%diag_tmp3d(i,j,k) = (GV%H_to_kg_m2*h(i,j,k))/CS%diag_tmp3d(i,j,k)
          enddo
          ! pressure for EOS at the bottom interface (Pa)
          do i=is,ie
            pressure_1d(i) = pressure_1d(i) + 0.5*(GV%g_Earth*GV%H_to_kg_m2)*h(i,j,k)
          enddo
        enddo ! k

      enddo ! j
      call post_data(CS%id_thkcello, CS%diag_tmp3d, CS%diag)
    endif
  endif

  ! volume mean potential temperature 
  if (CS%id_thetaoga>0) then
    thetaoga = global_volume_mean(tv%T, h, G, GV)
    call post_data(CS%id_thetaoga, thetaoga, CS%diag)
  endif

  ! area mean SST 
  if (CS%id_tosga > 0) then
    do j=js,je ; do i=is,ie 
       surface_field(i,j) = tv%T(i,j,1)
    enddo ; enddo 
    tosga = global_area_mean(surface_field, G)
    call post_data(CS%id_tosga, tosga, CS%diag)
  endif

  ! volume mean salinity 
  if (CS%id_soga>0) then
    soga = global_volume_mean(tv%S, h, G, GV)
    call post_data(CS%id_soga, soga, CS%diag)
  endif

  ! area mean SSS 
  if (CS%id_sosga > 0) then
    do j=js,je ; do i=is,ie 
       surface_field(i,j) = tv%S(i,j,1)
    enddo ; enddo 
    sosga = global_area_mean(surface_field, G)
    call post_data(CS%id_sosga, sosga, CS%diag)
  endif

  ! layer mean potential temperature 
  if (CS%id_temp_layer_ave>0) then
    temp_layer_ave = global_layer_mean(tv%T, h, G, GV)
    call post_data_1d_k(CS%id_temp_layer_ave, temp_layer_ave, CS%diag)
  endif

  ! layer mean salinity 
  if (CS%id_salt_layer_ave>0) then
    salt_layer_ave = global_layer_mean(tv%S, h, G, GV)
    call post_data_1d_k(CS%id_salt_layer_ave, salt_layer_ave, CS%diag)
  endif

  call calculate_vertical_integrals(h, tv, fluxes, G, GV, CS)

  if ((CS%id_Rml > 0) .or. (CS%id_Rcv > 0) .or. ASSOCIATED(CS%h_Rlay) .or. &
      ASSOCIATED(CS%uh_Rlay) .or. ASSOCIATED(CS%vh_Rlay) .or. &
      ASSOCIATED(CS%uhGM_Rlay) .or. ASSOCIATED(CS%vhGM_Rlay)) then

    if (associated(tv%eqn_of_state)) then
      pressure_1d(:) = tv%P_Ref
!$OMP parallel do default(none) shared(tv,Rcv,is,ie,js,je,nz,pressure_1d)
      do k=1,nz ; do j=js-1,je+1
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pressure_1d, &
                               Rcv(:,j,k), is-1, ie-is+3, tv%eqn_of_state)
      enddo ; enddo
    else ! Rcv should not be used much in this case, so fill in sensible values.
      do k=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
        Rcv(i,j,k) = GV%Rlay(k)
      enddo ; enddo ; enddo
    endif
    if (CS%id_Rml > 0) call post_data(CS%id_Rml, Rcv, CS%diag)
    if (CS%id_Rcv > 0) call post_data(CS%id_Rcv, Rcv, CS%diag)

    if (ASSOCIATED(CS%h_Rlay)) then
      k_list = nz/2
!$OMP parallel do default(none) shared(is,ie,js,je,nz,nkmb,CS,Rcv,h,GV) &
!$OMP                          private(wt,wt_p) firstprivate(k_list)
      do j=js,je 
        do k=1,nkmb ; do i=is,ie
          CS%h_Rlay(i,j,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do i=is,ie
          CS%h_Rlay(i,j,k) = h(i,j,k)
        enddo ; enddo 
        do k=1,nkmb ; do i=is,ie
          call find_weights(GV%Rlay, Rcv(i,j,k), k_list, nz, wt, wt_p)
          CS%h_Rlay(i,j,k_list)   = CS%h_Rlay(i,j,k_list)   + h(i,j,k)*wt
          CS%h_Rlay(i,j,k_list+1) = CS%h_Rlay(i,j,k_list+1) + h(i,j,k)*wt_p
        enddo ; enddo 
      enddo

      if (CS%id_h_Rlay > 0) call post_data(CS%id_h_Rlay, CS%h_Rlay, CS%diag)
    endif

    if (ASSOCIATED(CS%uh_Rlay)) then
      k_list = nz/2
!$OMP parallel do default(none) shared(Isq,Ieq,js,je,nz,nkmb,Rcv,CS,GV,uh) &
!$OMP                          private(wt,wt_p) firstprivate(k_list)
      do j=js,je
        do k=1,nkmb ; do I=Isq,Ieq
          CS%uh_Rlay(I,j,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do I=Isq,Ieq
          CS%uh_Rlay(I,j,k) = uh(I,j,k)
        enddo ; enddo 
        k_list = nz/2
        do k=1,nkmb ; do I=Isq,Ieq
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i+1,j,k)), k_list, nz, wt, wt_p)
          CS%uh_Rlay(I,j,k_list)   = CS%uh_Rlay(I,j,k_list)   + uh(I,j,k)*wt
          CS%uh_Rlay(I,j,k_list+1) = CS%uh_Rlay(I,j,k_list+1) + uh(I,j,k)*wt_p
        enddo ; enddo 
      enddo

      if (CS%id_uh_Rlay > 0) call post_data(CS%id_uh_Rlay, CS%uh_Rlay, CS%diag)
    endif

    if (ASSOCIATED(CS%vh_Rlay)) then
      k_list = nz/2
!$OMP parallel do default(none)  shared(Jsq,Jeq,is,ie,nz,nkmb,Rcv,CS,GV,vh) &
!$OMP                          private(wt,wt_p) firstprivate(k_list)
      do J=Jsq,Jeq
        do k=1,nkmb ; do i=is,ie
          CS%vh_Rlay(i,J,k) = 0.0
        enddo ; enddo 
        do k=nkmb+1,nz ; do i=is,ie
          CS%vh_Rlay(i,J,k) = vh(i,J,k)
        enddo ; enddo
        do k=1,nkmb ; do i=is,ie
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i,j+1,k)), k_list, nz, wt, wt_p)
          CS%vh_Rlay(i,J,k_list)   = CS%vh_Rlay(i,J,k_list)   + vh(i,J,k)*wt
          CS%vh_Rlay(i,J,k_list+1) = CS%vh_Rlay(i,J,k_list+1) + vh(i,J,k)*wt_p
        enddo ; enddo
      enddo

      if (CS%id_vh_Rlay > 0) call post_data(CS%id_vh_Rlay, CS%vh_Rlay, CS%diag)
    endif

    if (ASSOCIATED(CS%uhGM_Rlay) .and. ASSOCIATED(CDp%uhGM)) then
      k_list = nz/2
!$OMP parallel do default(none) shared(Isq,Ieq,js,je,nz,nkmb,Rcv,CDP,CS,GV) &
!$OMP                          private(wt,wt_p) firstprivate(k_list)
      do j=js,je
        do k=1,nkmb ; do I=Isq,Ieq
          CS%uhGM_Rlay(I,j,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do I=Isq,Ieq
          CS%uhGM_Rlay(I,j,k) = CDp%uhGM(I,j,k)
        enddo ; enddo
        do k=1,nkmb ; do I=Isq,Ieq
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i+1,j,k)), k_list, nz, wt, wt_p)
          CS%uhGM_Rlay(I,j,k_list)   = CS%uhGM_Rlay(I,j,k_list)   + CDp%uhGM(I,j,k)*wt
          CS%uhGM_Rlay(I,j,k_list+1) = CS%uhGM_Rlay(I,j,k_list+1) + CDp%uhGM(I,j,k)*wt_p
        enddo ; enddo
      enddo

      if (CS%id_uh_Rlay > 0) call post_data(CS%id_uhGM_Rlay, CS%uhGM_Rlay, CS%diag)
    endif

    if (ASSOCIATED(CS%vhGM_Rlay) .and. ASSOCIATED(CDp%vhGM)) then
      k_list = nz/2
!$OMP parallel do default(none) shared(is,ie,Jsq,Jeq,nz,nkmb,CS,CDp,Rcv,GV) &
!$OMP                          private(wt,wt_p) firstprivate(k_list)
      do J=Jsq,Jeq
        do k=1,nkmb ; do i=is,ie
          CS%vhGM_Rlay(i,J,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do i=is,ie
          CS%vhGM_Rlay(i,J,k) = CDp%vhGM(i,J,k)
        enddo ; enddo
        do k=1,nkmb ; do i=is,ie
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i,j+1,k)), k_list, nz, wt, wt_p)
          CS%vhGM_Rlay(i,J,k_list)   = CS%vhGM_Rlay(i,J,k_list)   + CDp%vhGM(i,J,k)*wt
          CS%vhGM_Rlay(i,J,k_list+1) = CS%vhGM_Rlay(i,J,k_list+1) + CDp%vhGM(i,J,k)*wt_p
        enddo ; enddo
      enddo

      if (CS%id_vhGM_Rlay > 0) call post_data(CS%id_vhGM_Rlay, CS%vhGM_Rlay, CS%diag)
    endif
  endif

  if (associated(tv%eqn_of_state)) then
    if (CS%id_rhopot0 > 0) then
      pressure_1d(:) = 0.
!$OMP parallel do default(none) shared(tv,Rcv,is,ie,js,je,nz,pressure_1d)
      do k=1,nz ; do j=js,je
        call calculate_density(tv%T(:,j,k),tv%S(:,j,k),pressure_1d, &
                               Rcv(:,j,k),is,ie-is+1, tv%eqn_of_state)
      enddo ; enddo
      if (CS%id_rhopot0 > 0) call post_data(CS%id_rhopot0, Rcv, CS%diag)
    endif
    if (CS%id_rhopot2 > 0) then
      pressure_1d(:) = 2.E7 ! 2000 dbars
!$OMP parallel do default(none) shared(tv,Rcv,is,ie,js,je,nz,pressure_1d)
      do k=1,nz ; do j=js,je
        call calculate_density(tv%T(:,j,k),tv%S(:,j,k),pressure_1d, &
                               Rcv(:,j,k),is,ie-is+1, tv%eqn_of_state)
      enddo ; enddo
      if (CS%id_rhopot2 > 0) call post_data(CS%id_rhopot2, Rcv, CS%diag)
    endif
    if (CS%id_rhoinsitu > 0) then
!$OMP parallel do default(none) shared(tv,Rcv,is,ie,js,je,nz,pressure_1d,h,GV)
      do j=js,je
        pressure_1d(:) = 0. ! Start at p=0 Pa at surface
        do k=1,nz
          pressure_1d(:) =  pressure_1d(:) + 0.5 * h(:,j,k) * GV%H_to_Pa ! Pressure in middle of layer k
          call calculate_density(tv%T(:,j,k),tv%S(:,j,k),pressure_1d, &
                                 Rcv(:,j,k),is,ie-is+1, tv%eqn_of_state)
          pressure_1d(:) =  pressure_1d(:) + 0.5 * h(:,j,k) * GV%H_to_Pa ! Pressure at bottom of layer k
        enddo
      enddo
      if (CS%id_rhoinsitu > 0) call post_data(CS%id_rhoinsitu, Rcv, CS%diag)
    endif
  endif

  if ((CS%id_cg1>0) .or. (CS%id_Rd1>0) .or. (CS%id_cfl_cg1>0) .or. &
      (CS%id_cfl_cg1_x>0) .or. (CS%id_cfl_cg1_y>0)) then
    call wave_speed(h, tv, G, GV, CS%cg1, CS%wave_speed_CSp)
    if (CS%id_cg1>0) call post_data(CS%id_cg1, CS%cg1, CS%diag)
    if (CS%id_Rd1>0) then
!$OMP parallel do default(none) shared(is,ie,js,je,G,CS) &
!$OMP                          private(f2_h,mag_beta)
      do j=js,je ; do i=is,ie
        ! Blend the equatorial deformation radius with the standard one.
        f2_h = absurdly_small_freq2 + 0.25 * &
            ((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
             (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I,J-1)**2))
        mag_beta = sqrt(0.5 * ( &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2 + &
             ((G%CoriolisBu(I,J-1)-G%CoriolisBu(I-1,J-1)) * G%IdxCv(i,J-1))**2) + &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 + &
             ((G%CoriolisBu(I-1,J)-G%CoriolisBu(I-1,J-1)) * G%IdyCu(I-1,j))**2) ))
        CS%Rd1(i,j) = CS%cg1(i,j) / sqrt(f2_h + CS%cg1(i,j) * mag_beta)

      enddo ; enddo
      call post_data(CS%id_Rd1, CS%Rd1, CS%diag)
    endif
    if (CS%id_cfl_cg1>0) then
      do j=js,je ; do i=is,ie
        CS%cfl_cg1(i,j) = (dt*CS%cg1(i,j)) * (G%IdxT(i,j) + G%IdyT(i,j))
      enddo ; enddo
      call post_data(CS%id_cfl_cg1, CS%cfl_cg1, CS%diag)
    endif
    if (CS%id_cfl_cg1_x>0) then
      do j=js,je ; do i=is,ie
        CS%cfl_cg1_x(i,j) = (dt*CS%cg1(i,j)) * G%IdxT(i,j)
      enddo ; enddo
      call post_data(CS%id_cfl_cg1_x, CS%cfl_cg1_x, CS%diag)
    endif
    if (CS%id_cfl_cg1_y>0) then
      do j=js,je ; do i=is,ie
        CS%cfl_cg1_y(i,j) = (dt*CS%cg1(i,j)) * G%IdyT(i,j)
      enddo ; enddo
      call post_data(CS%id_cfl_cg1_y, CS%cfl_cg1_y, CS%diag)
    endif
  endif

  if (dt > 0.0) then
    if (CS%id_du_dt>0) call post_data(CS%id_du_dt, CS%du_dt, CS%diag)

    if (CS%id_dv_dt>0) call post_data(CS%id_dv_dt, CS%dv_dt, CS%diag)

    if (CS%id_dh_dt>0) call post_data(CS%id_dh_dt, CS%dh_dt, CS%diag)

    call calculate_energy_diagnostics(u, v, h, uh, vh, ADp, CDp, G, CS)
    call calculate_twa_diagnostics(u, v, h, uh, vh, ADp, CDp, G, GV, CS)
  endif

end subroutine calculate_diagnostic_fields


subroutine find_weights(Rlist, R_in, k, nz, wt, wt_p)
  real,     intent(in)    :: Rlist(:), R_in
  integer,  intent(inout) :: k
  integer,  intent(in)    :: nz
  real,     intent(out)   :: wt, wt_p

  ! This subroutine finds location of R_in in an increasing ordered
  ! list, Rlist, returning as k the element such that
  !  Rlist(k) <= R_in < Rlist(k+1), and where wt and wt_p are the linear
  ! weights that should be assigned to elements k and k+1.

  integer :: k_upper, k_lower, k_new, inc

  ! First, bracket the desired point.
  if ((k < 1) .or. (k > nz)) k = nz/2

  k_upper = k ; k_lower = k ; inc = 1
  if (R_in < Rlist(k)) then
    do
      k_lower = max(k_lower-inc, 1)
      if ((k_lower == 1) .or. (R_in >= Rlist(k_lower))) exit
      k_upper = k_lower
      inc = inc*2
    end do
  else
    do
      k_upper = min(k_upper+inc, nz)
      if ((k_upper == nz) .or. (R_in < Rlist(k_upper))) exit
      k_lower = k_upper
      inc = inc*2
    end do
  endif

  if ((k_lower == 1) .and. (R_in <= Rlist(k_lower))) then
    k = 1 ; wt = 1.0 ; wt_p = 0.0
  else if ((k_upper == nz) .and. (R_in >= Rlist(k_upper))) then
    k = nz-1 ; wt = 0.0 ; wt_p = 1.0
  else
    do
      if (k_upper <= k_lower+1) exit
      k_new = (k_upper + k_lower) / 2
      if (R_in < Rlist(k_new)) then
        k_upper = k_new
      else
        k_lower = k_new
      endif
    end do

!   Uncomment this as a code check
!    if ((R_in < Rlist(k_lower)) .or. (R_in >= Rlist(k_upper)) .or. (k_upper-k_lower /= 1)) &
!      write (*,*) "Error: ",R_in," is not between R(",k_lower,") = ", &
!        Rlist(k_lower)," and R(",k_upper,") = ",Rlist(k_upper),"."
    k = k_lower
    wt = (Rlist(k_upper) - R_in) / (Rlist(k_upper) - Rlist(k_lower))
    wt_p = 1.0 - wt

  endif

end subroutine find_weights

subroutine calculate_vertical_integrals(h, tv, fluxes, G, GV, CS)
  type(ocean_grid_type),                    intent(inout) :: G
  type(verticalGrid_type),                  intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h
  type(thermo_var_ptrs),                    intent(in)    :: tv
  type(forcing),                            intent(in)    :: fluxes
  type(diagnostics_CS),                     intent(inout) :: CS

! Subroutine calculates vertical integrals of several tracers, along
! with the mass-weight of these tracers, the total column mass, and the
! carefully calculated column height.

! Arguments: 
!  (in)      h  - layer thickness: metre (Bouss) or kg/ m2 (non-Bouss)
!  (in)      tv - structure pointing to thermodynamic variables
!  (in)      fluxes - a structure containing the surface fluxes.
!  (in)      G  - ocean grid structure
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - control structure returned by a previous call to diagnostics_init

  real, dimension(SZI_(G), SZJ_(G)) :: &
    z_top, &  ! Height of the top of a layer or the ocean, in m.
    z_bot, &  ! Height of the bottom of a layer (for id_mass) or the
              ! (positive) depth of the ocean (for id_col_ht), in m.
    mass, &   ! integrated mass of the water column, in kg m-2.  For
              ! non-Boussinesq models this is rho*dz. For Boussiensq
              ! models, this is either the integral of in-situ density
              ! (rho*dz for col_mass) or reference dens (Rho_0*dz for mass_wt).
    btm_pres,&! The pressure at the ocean bottom, or CMIP variable 'pbo'.
              ! This is the column mass multiplied by gravity plus the pressure
              ! at the ocean surface.
    dpress, &    ! Change in hydrostatic pressure across a layer, in Pa.
    tr_int    ! vertical integral of a tracer times density,
              ! (Rho_0 in a Boussinesq model) in TR kg m-2.
  real    :: IG_Earth  ! Inverse of gravitational acceleration, in s2 m-1.

  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  
  if (CS%id_mass_wt > 0) then
    do j=js,je ; do i=is,ie ; mass(i,j) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      mass(i,j) = mass(i,j) + GV%H_to_kg_m2*h(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_mass_wt, mass, CS%diag)
  endif

  if (CS%id_temp_int > 0) then
    do j=js,je ; do i=is,ie ; tr_int(i,j) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      tr_int(i,j) = tr_int(i,j) + (GV%H_to_kg_m2*h(i,j,k))*tv%T(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_temp_int, tr_int, CS%diag)
  endif

  if (CS%id_salt_int > 0) then
    do j=js,je ; do i=is,ie ; tr_int(i,j) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      tr_int(i,j) = tr_int(i,j) + (GV%H_to_kg_m2*h(i,j,k))*tv%S(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_salt_int, tr_int, CS%diag)
  endif

  if (CS%id_col_ht > 0) then
    call find_eta(h, tv, GV%g_Earth, G, GV, z_top)
    do j=js,je ; do i=is,ie
      z_bot(i,j) = z_top(i,j) + G%bathyT(i,j)
    enddo ; enddo
    call post_data(CS%id_col_ht, z_bot, CS%diag)
  endif

  if (CS%id_col_mass > 0 .or. CS%id_pbo > 0) then
    do j=js,je ; do i=is,ie ; mass(i,j) = 0.0 ; enddo ; enddo
    if (GV%Boussinesq) then
      if (associated(tv%eqn_of_state)) then
        IG_Earth = 1.0 / GV%g_Earth
!       do j=js,je ; do i=is,ie ; z_bot(i,j) = -P_SURF(i,j)/GV%H_to_Pa ; enddo ; enddo
        do j=js,je ; do i=is,ie ; z_bot(i,j) = 0.0 ; enddo ; enddo
        do k=1,nz
          do j=js,je ; do i=is,ie
            z_top(i,j) = z_bot(i,j)
            z_bot(i,j) = z_top(i,j) - GV%H_to_m*h(i,j,k)
          enddo ; enddo
          call int_density_dz(tv%T(:,:,k), tv%S(:,:,k), &
                              z_top, z_bot, 0.0, GV%H_to_kg_m2, GV%g_Earth, &
                              G%HI, G%HI, tv%eqn_of_state, dpress)
          do j=js,je ; do i=is,ie
            mass(i,j) = mass(i,j) + dpress(i,j) * IG_Earth
          enddo ; enddo
        enddo
      else
        do k=1,nz ; do j=js,je ; do i=is,ie
          mass(i,j) = mass(i,j) + (GV%H_to_m*GV%Rlay(k))*h(i,j,k)
        enddo ; enddo ; enddo
      endif
    else
      do k=1,nz ; do j=js,je ; do i=is,ie
        mass(i,j) = mass(i,j) + GV%H_to_kg_m2*h(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (CS%id_col_mass > 0) then
      call post_data(CS%id_col_mass, mass, CS%diag)
    endif
    if (CS%id_pbo > 0) then
      do j=js,je ; do i=is,ie ; btm_pres(i,j) = 0.0 ; enddo ; enddo
      ! 'pbo' is defined as the sea water pressure at the sea floor
      !     pbo = (mass * g) + pso
      ! where pso is the sea water pressure at sea water surface
      ! note that pso is equivalent to fluxes%p_surf
      do j=js,je ; do i=is,ie
        btm_pres(i,j) = mass(i,j) * GV%g_Earth
        if (ASSOCIATED(fluxes%p_surf)) then
          btm_pres(i,j) = btm_pres(i,j) + fluxes%p_surf(i,j)
        endif
      enddo ; enddo
      call post_data(CS%id_pbo, btm_pres, CS%diag)
    endif
  endif

end subroutine calculate_vertical_integrals


subroutine calculate_energy_diagnostics(u, v, h, uh, vh, ADp, CDp, G, CS)
  type(ocean_grid_type),                     intent(inout) :: G
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: uh
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: vh
  type(accel_diag_ptrs),                     intent(in)    :: ADp
  type(cont_diag_ptrs),                      intent(in)    :: CDp
  type(diagnostics_CS),                      intent(inout) :: CS

! This subroutine calculates terms in the mechanical energy budget.

! Arguments:
!  (in)      u   - zonal velocity component (m/s)
!  (in)      v   - meridional velocity componnent (m/s)
!  (in)      h   - layer thickness: metre(Bouss) of kg/m2(non-Bouss)
!  (in)      uh  - transport through zonal faces=u*h*dy: m3/s (Bouss) kg/s(non-Bouss)
!  (in)      vh  - transport through merid faces=v*h*dx: m3/s (Bouss) kg/s(non-Bouss)
!  (in)      ADp - structure pointing to accelerations in momentum equation
!  (in)      CDp - structure pointing to terms in continuity equations
!  (in)      G   - ocean grid structure
!  (in)      CS  - control structure returned by a previous call to diagnostics_init

  real :: KE_u(SZIB_(G),SZJ_(G))
  real :: KE_v(SZI_(G),SZJB_(G))
  real :: KE_h(SZI_(G),SZJ_(G))

  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  do j=js-1,je ; do i=is-1,ie
    KE_u(I,j) = 0.0 ; KE_v(i,J) = 0.0
  enddo ; enddo

  if (ASSOCIATED(CS%KE)) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      CS%KE(i,j,k) = ((u(I,j,k)*u(I,j,k) + u(I-1,j,k)*u(I-1,j,k)) + &
          (v(i,J,k)*v(i,J,k) + v(i,J-1,k)*v(i,J-1,k)))*0.25
      ! DELETE THE FOLLOWING...  Make this 0 to test the momentum balance,
      ! or a huge number to test the continuity balance.
      ! CS%KE(i,j,k) *= 1e20
    enddo ; enddo ; enddo
    if (CS%id_KE > 0) call post_data(CS%id_KE, CS%KE, CS%diag)
  endif

  if(.not.G%symmetric) then
    if(ASSOCIATED(CS%dKE_dt) .OR. ASSOCIATED(CS%PE_to_KE) .OR. ASSOCIATED(CS%KE_CorAdv) .OR. &
       ASSOCIATED(CS%KE_adv) .OR. ASSOCIATED(CS%KE_visc)  .OR. ASSOCIATED(CS%KE_horvisc).OR. &
       ASSOCIATED(CS%KE_dia) ) then
        call create_group_pass(CS%pass_KE_uv, KE_u, KE_v, G%Domain, To_North+To_East)
    endif
  endif

  if (ASSOCIATED(CS%dKE_dt)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k)*G%dxCu(I,j)*CS%du_dt(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k)*G%dyCv(i,J)*CS%dv_dt(i,J,k)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = CS%KE(i,j,k)*CS%dh_dt(i,j,k)
      enddo ; enddo
      if (.not.G%symmetric) &      
         call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%dKE_dt(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) * &
            (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_dKEdt > 0) call post_data(CS%id_dKEdt, CS%dKE_dt, CS%diag)
  endif

  if (ASSOCIATED(CS%PE_to_KE)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k)*G%dxCu(I,j)*ADp%PFu(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k)*G%dyCv(i,J)*ADp%PFv(i,J,k)
      enddo ; enddo
      if (.not.G%symmetric) &
         call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%PE_to_KE(i,j,k) = 0.5 * G%IareaT(i,j) * &
            (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_PE_to_KE > 0) call post_data(CS%id_PE_to_KE, CS%PE_to_KE, CS%diag)
  endif

  if (ASSOCIATED(CS%KE_CorAdv)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k)*G%dxCu(I,j)*ADp%CAu(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k)*G%dyCv(i,J)*ADp%CAv(i,J,k)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = -CS%KE(i,j,k) * G%IareaT(i,j) * &
            (uh(I,j,k) - uh(I-1,j,k) + vh(i,J,k) - vh(i,J-1,k))
      enddo ; enddo
      if (.not.G%symmetric) &
         call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%KE_CorAdv(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) * &
            (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_Coradv > 0) call post_data(CS%id_KE_Coradv, CS%KE_Coradv, CS%diag)
  endif

  if (ASSOCIATED(CS%KE_adv)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k)*G%dxCu(I,j)*ADp%gradKEu(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k)*G%dyCv(i,J)*ADp%gradKEv(i,J,k)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = -CS%KE(i,j,k) * G%IareaT(i,j) * &
            (uh(I,j,k) - uh(I-1,j,k) + vh(i,J,k) - vh(i,J-1,k))
      enddo ; enddo
      if (.not.G%symmetric) &
         call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%KE_adv(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) * &
            (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_adv > 0) call post_data(CS%id_KE_adv, CS%KE_adv, CS%diag)
  endif

  if (ASSOCIATED(CS%KE_visc)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k)*G%dxCu(I,j)*ADp%du_dt_visc(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k)*G%dyCv(i,J)*ADp%dv_dt_visc(i,J,k)
      enddo ; enddo
      if (.not.G%symmetric) &
         call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%KE_visc(i,j,k) = 0.5 * G%IareaT(i,j) * &
            (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_visc > 0) call post_data(CS%id_KE_visc, CS%KE_visc, CS%diag)
  endif

  if (ASSOCIATED(CS%KE_horvisc)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k)*G%dxCu(I,j)*ADp%diffu(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k)*G%dyCv(i,J)*ADp%diffv(i,J,k)
      enddo ; enddo
      if (.not.G%symmetric) &
         call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%KE_horvisc(i,j,k) = 0.5 * G%IareaT(i,j) * &
            (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_horvisc > 0) call post_data(CS%id_KE_horvisc, CS%KE_horvisc, CS%diag)
  endif

  if (ASSOCIATED(CS%KE_dia)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k)*G%dxCu(I,j)*ADp%du_dt_dia(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k)*G%dyCv(i,J)*ADp%dv_dt_dia(i,J,k)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = CS%KE(i,j,k) * &
            (CDp%diapyc_vel(i,j,k) - CDp%diapyc_vel(i,j,k+1))
      enddo ; enddo
      if (.not.G%symmetric) &
         call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%KE_dia(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) * &
            (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_dia > 0) call post_data(CS%id_KE_dia, CS%KE_dia, CS%diag)
  endif

end subroutine calculate_energy_diagnostics

subroutine calculate_twa_diagnostics(u, v, h, uh, vh, ADp, CDp, G, GV, CS)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in)    :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in)    :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: h
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in)    :: uh
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in)    :: vh
  type(accel_diag_ptrs),                  intent(in)    :: ADp
  type(cont_diag_ptrs),                   intent(in)    :: CDp
  type(ocean_grid_type),                  intent(inout) :: G
  type(verticalGrid_type),                intent(in)    :: GV
  type(diagnostics_CS),                   intent(inout) :: CS

! This subroutine calculates terms in the twa budget.

! Arguments:
!  (in)      u   - zonal velocity component (m/s)
!  (in)      v   - meridional velocity componnent (m/s)
!  (in)      h   - layer thickness: metre(Bouss) of kg/m2(non-Bouss)
!  (in)      uh  - transport through zonal faces=u*h*dy: m3/s (Bouss) kg/s(non-Bouss)
!  (in)      vh  - transport through merid faces=v*h*dx: m3/s (Bouss) kg/s(non-Bouss)
!  (in)      ADp - structure pointing to accelerations in momentum equation
!  (in)      CDp - structure pointing to terms in continuity equations
!  (in)      G   - ocean grid structure
!  (in)      CS  - control structure returned by a previous call to diagnostics_init

  real :: dwd, uhxCu,vhyCu, uhxCv, vhyCv
  real :: hmin
  real, dimension(SZK_(G)+1) :: wdatui, wdatvi
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: ishqlarge
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: wdatui1,wdatvi1

  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  do k = 1,nz
    do J=Jsq,Jeq ; do I=Isq,Ieq
      hmin = min(h(i,j,k), h(i+1,j,k), h(i,j+1,k), h(i+1,j+1,k))
      ishqlarge(I,J,k) = ceiling(abs(hmin-GV%Angstrom_z)/hmin)
    enddo ; enddo
  enddo

  if (ASSOCIATED(CS%h_Cu)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        CS%h_Cu(I,j,k) = 0.5*(h(i,j,k) + h(i+1,j,k))*ishqlarge(I,J,k)
      enddo ; enddo
    enddo
    if (CS%id_h_Cu > 0) call post_data(CS%id_h_Cu, CS%h_Cu, CS%diag)
  endif

  if (ASSOCIATED(CS%h_Cv)) then
    do k=1,nz
      do j=Jsq,Jeq ; do i=is,ie
        CS%h_Cv(i,J,k) = 0.5*(h(i,j,k) + h(i,j+1,k))*ishqlarge(I,J,k)
      enddo ; enddo
    enddo
    if (CS%id_h_Cv > 0) call post_data(CS%id_h_Cv, CS%h_Cv, CS%diag)
  endif

  if (ASSOCIATED(CS%huu_Cu)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        CS%huu_Cu(I,j,k) = CS%h_Cu(I,j,k)*u(I,j,k)*u(I,j,k)*ishqlarge(I,J,k)
      enddo ; enddo
    enddo
    if (CS%id_huu_Cu > 0) call post_data(CS%id_huu_Cu, CS%huu_Cu, CS%diag)
  endif

  if (ASSOCIATED(CS%hvv_Cv)) then
    do k=1,nz
      do J=Jsq,Jeq ; do i=is,ie
        CS%hvv_Cv(i,J,k) = CS%h_Cv(i,j,k)*v(i,J,k)*v(i,J,k)*ishqlarge(I,J,k)
      enddo ; enddo
    enddo
    if (CS%id_hvv_Cv > 0) call post_data(CS%id_hvv_Cv, CS%hvv_Cv, CS%diag)
  endif

  if (ASSOCIATED(CS%hv_Cu)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        CS%hv_Cu(I,j,k) = CS%h_Cu(I,j,k)*0.25* &
          (v(i,j,k)+v(i,j-1,k)+v(i+1,j,k)+v(i+1,j-1,k))
      enddo ; enddo
    enddo
    if (CS%id_hv_Cu > 0) call post_data(CS%id_hv_Cu, CS%hv_Cu, CS%diag)
  endif

  if (ASSOCIATED(CS%hu_Cv)) then
    do k=1,nz
      do J=Jsq,Jeq ; do i=is,ie
        CS%hu_Cv(i,J,k) = CS%h_Cv(i,J,k)*0.25* &
          (u(i,j,k)+u(i,j+1,k)+u(i-1,j,k)+u(i-1,j+1,k))
      enddo ; enddo
    enddo
    if (CS%id_hu_Cv > 0) call post_data(CS%id_hu_Cv, CS%hu_Cv, CS%diag)
  endif

!  if (ASSOCIATED(CS%hw_Cu)) then
!    do j=js,je ; do I=Isq,Ieq
!      wdatui(1) = 0.0
!      do k=2,nz
!        wdatui(k) = 0.5*(CDp%diapyc_vel(i,j,k)+CDp%diapyc_vel(i+1,j,k))*GV%g_prime(k)
!        CS%hw_Cu(I,j,k-1) = 0.5*(wdatui(k) + wdatui(k-1))*ishqlarge(I,J,k-1)
!      enddo
!      wdatui(nz+1) = 0.0
!      CS%hw_Cu(I,j,nz) = 0.5*(wdatui(nz+1) + wdatui(nz))*ishqlarge(I,J,nz)
!    enddo ; enddo
!    if (CS%id_hw_Cu > 0) call post_data(CS%id_hw_Cu, CS%hw_Cu, CS%diag)
!  endif
!
!  if (ASSOCIATED(CS%hw_Cv)) then
!    do J=Jsq,Jeq ; do i=is,ie
!      wdatvi(1) = 0.0
!      do k=2,nz
!        wdatvi(k) = 0.5*(CDp%diapyc_vel(i,j,k)+CDp%diapyc_vel(i,j+1,k))*GV%g_prime(k)
!        CS%hw_Cv(i,J,k-1) = 0.5*(wdatvi(k) + wdatvi(k-1))*ishqlarge(I,J,k-1)
!      enddo
!      wdatvi(nz+1) = 0.0
!      CS%hw_Cv(i,J,nz) = 0.5*(wdatvi(nz+1) + wdatvi(nz))*ishqlarge(I,J,nz)
!    enddo ; enddo
!    if (CS%id_hw_Cv > 0) call post_data(CS%id_hw_Cv, CS%hw_Cv, CS%diag)
!  endif
  if (ASSOCIATED(CS%hw_Cu)) then
    do j=js,je ; do I=Isq,Ieq
      do k=1,nz-1
        CS%hw_Cu(I,j,k) = 0.25*((CDp%diapyc_vel(i,j,k)+CDp%diapyc_vel(i+1,j,k))*GV%g_prime(k)&
            +(CDp%diapyc_vel(i,j,k+1)+CDp%diapyc_vel(i+1,j,k+1))*GV%g_prime(k+1))*ishqlarge(I,J,k)
      enddo
      CS%hw_Cu(I,j,nz) = 0.25*(CDp%diapyc_vel(i,j,nz)+CDp%diapyc_vel(i+1,j,nz))&
          *GV%g_prime(nz)*ishqlarge(I,J,nz)
    enddo ; enddo
    if (CS%id_hw_Cu > 0) call post_data(CS%id_hw_Cu, CS%hw_Cu, CS%diag)
  endif

  if (ASSOCIATED(CS%hw_Cv)) then
    do J=Jsq,Jeq ; do i=is,ie
      do k=1,nz-1
        CS%hw_Cv(i,J,k) = 0.25*((CDp%diapyc_vel(i,j,k)+CDp%diapyc_vel(i,j+1,k))*GV%g_prime(k)&
            +(CDp%diapyc_vel(i,j,k+1)+CDp%diapyc_vel(i,j+1,k+1))*GV%g_prime(k+1))*ishqlarge(I,J,k)
      enddo
      CS%hw_Cv(i,J,nz) = 0.25*(CDp%diapyc_vel(i,j,nz)+CDp%diapyc_vel(i+1,j,nz))&
          *GV%g_prime(nz)*ishqlarge(I,J,nz)
    enddo ; enddo
    if (CS%id_hw_Cv > 0) call post_data(CS%id_hw_Cv, CS%hw_Cv, CS%diag)
  endif

  if (ASSOCIATED(CS%hwb_Cu)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        CS%hwb_Cu(I,j,k) = 0.25*(CDp%diapyc_vel(i,j,k) + CDp%diapyc_vel(i+1,j,k) + &
          CDp%diapyc_vel(i,j,k+1) + CDp%diapyc_vel(i+1,j,k+1))*ishqlarge(I,J,k)
      enddo ; enddo
    enddo
    if (CS%id_hwb_Cu > 0) call post_data(CS%id_hwb_Cu, CS%hwb_Cu, CS%diag)
  endif

  if (ASSOCIATED(CS%hwb_Cv)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        CS%hwb_Cv(i,J,k) = 0.25*(CDp%diapyc_vel(i,j,k) + CDp%diapyc_vel(i,j+1,k) + &
          CDp%diapyc_vel(i,j,k+1) + CDp%diapyc_vel(i,j+1,k+1))*ishqlarge(I,J,k)
      enddo ; enddo
    enddo
    if (CS%id_hwb_Cv > 0) call post_data(CS%id_hwb_Cv, CS%hwb_Cv, CS%diag)
  endif

  if (ASSOCIATED(CS%esq_Cu)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        CS%esq_Cu(i,j,k) = 0.25*(CS%e(i,j,k)*CS%e(i,j,k)+CS%e(i,j,k+1)*CS%e(i,j,k+1)&
          +CS%e(i+1,j,k)*CS%e(i+1,j,k)+CS%e(i+1,j,k+1)*CS%e(i+1,j,k+1))*ishqlarge(I,J,k)
      enddo ; enddo
    enddo
    if (CS%id_esq_Cu > 0) call post_data(CS%id_esq_Cu, CS%esq_Cu, CS%diag)
  endif

  if (ASSOCIATED(CS%esq_Cv)) then
    do k=1,nz
      do J=Jsq,Jeq ; do i=is,ie
        CS%esq_Cv(i,j,k) = 0.25*(CS%e(i,j,k)*CS%e(i,j,k)+CS%e(i,j,k+1)*CS%e(i,j,k+1)&
          +CS%e(i,j+1,k)*CS%e(i,j+1,k)+CS%e(i,j+1,k+1)*CS%e(i,j+1,k+1))*ishqlarge(I,J,k)
      enddo ; enddo
    enddo
    if (CS%id_esq_Cv > 0) call post_data(CS%id_esq_Cv, CS%esq_Cv, CS%diag)
  endif

  if (ASSOCIATED(CS%e_Cu) .OR. ASSOCIATED(CS%epfu)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        CS%e_Cu(I,j,k) = 0.25*(CS%e(i,j,k)+CS%e(i+1,j,k)+&
          CS%e(i,j,k+1)+CS%e(i+1,j,k+1))*ishqlarge(I,J,k)
        CS%epfu(I,j,k) = ADp%PFu(I,j,k)*CS%e_Cu(I,j,k)
      enddo ; enddo
    enddo
    if (CS%id_e_Cu > 0) call post_data(CS%id_e_Cu, CS%e_Cu, CS%diag)
    if (CS%id_epfu > 0) call post_data(CS%id_epfu, CS%epfu, CS%diag)
  endif

  if (ASSOCIATED(CS%e_Cv) .OR. ASSOCIATED(CS%epfv)) then
    do k=1,nz
      do J=Jsq,Jeq ; do i=is,ie
        CS%e_Cv(i,J,k) = 0.25*(CS%e(i,j,k)+CS%e(i,j+1,k)+&
          CS%e(i,j,k+1)+CS%e(i,j+1,k+1))*ishqlarge(I,J,k)
        CS%epfv(i,J,k) = ADp%PFv(i,J,k)*CS%e_Cv(i,J,k)
      enddo ; enddo
    enddo
    if (CS%id_e_Cv > 0) call post_data(CS%id_e_Cv, CS%e_Cv, CS%diag)
    if (CS%id_epfv > 0) call post_data(CS%id_epfv, CS%epfv, CS%diag)
  endif

  if (ASSOCIATED(CS%pfu_masked)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        CS%pfu_masked(i,J,k) = ADp%PFu(I,j,k)*ishqlarge(I,J,k)
      enddo ; enddo
    enddo
    if (CS%id_pfu_masked > 0) call post_data(CS%id_pfu_masked, CS%pfu_masked, CS%diag)
  endif

  if (ASSOCIATED(CS%pfv_masked)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        CS%pfv_masked(i,J,k) = ADp%PFv(I,j,k)*ishqlarge(I,J,k)
      enddo ; enddo
    enddo
    if (CS%id_pfv_masked > 0) call post_data(CS%id_pfv_masked, CS%pfv_masked, CS%diag)
  endif

  if (ASSOCIATED(CS%hfv)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        CS%hfv(I,j,k) = CS%h_Cu(I,j,k)*(ADp%CAu(I,j,k) - ADp%gradKEu(I,j,k) - ADP%rv_x_v(I,j,k))
      enddo ; enddo
    enddo
    if (CS%id_hfv > 0) call post_data(CS%id_hfv, CS%hfv, CS%diag)
  endif
  
  if (ASSOCIATED(CS%hpfu)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        CS%hpfu(i,J,k) = CS%h_Cu(I,j,k)*ADp%PFu(I,j,k)
      enddo ; enddo
    enddo
    if (CS%id_hpfu > 0) call post_data(CS%id_hpfu, CS%hpfu, CS%diag)
  endif

  if (ASSOCIATED(CS%huwb)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        dwd = 0.5*(CDp%diapyc_vel(i,j,k+1) - CDp%diapyc_vel(i,j,k) &
          + CDp%diapyc_vel(i+1,j,k+1) - CDp%diapyc_vel(i+1,j,k))*ishqlarge(I,J,k)
        CS%huwb(I,j,k) = -CS%h_Cu(I,j,k)*ADp%du_dt_dia(I,j,k) - dwd*u(I,j,k)
      enddo ; enddo
    enddo
    if (CS%id_huwb > 0) call post_data(CS%id_huwb, CS%huwb, CS%diag)
  endif
  
  if (ASSOCIATED(CS%huuxpt)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        uhxCu = 0.5*((CDp%uh(I,j,k) - CDp%uh(I-1,j,k))*G%IareaT(i,j) &
                + (CDp%uh(I+1,j,k) - CDp%uh(I,j,k))*G%IareaT(i+1,j))*ishqlarge(I,J,k) 
        CS%huuxpt(I,j,k) = CS%h_Cu(I,j,k)*ADp%gradKEu(I,j,k) - uhxCu*u(I,j,k)
      enddo ; enddo
    enddo
    if (CS%id_huuxpt > 0) call post_data(CS%id_huuxpt, CS%huuxpt, CS%diag)
  endif
  
  if (ASSOCIATED(CS%huvymt)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        vhyCu = 0.5*((CDp%vh(i,j,k) - CDp%vh(i,j-1,k))*G%IareaT(i,j) &
                + (CDp%vh(i+1,j,k) - CDp%vh(i+1,j-1,k))*G%IareaT(i+1,j))*ishqlarge(I,J,k)
        CS%huvymt(I,j,k) = CS%h_Cu(I,j,k)*ADp%rv_x_v(I,j,k) - vhyCu*u(I,j,k)
      enddo ; enddo
    enddo
    if (CS%id_huvymt > 0) call post_data(CS%id_huvymt, CS%huvymt, CS%diag)
  endif

  if (ASSOCIATED(CS%hdudtvisc)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        CS%hdudtvisc(I,j,k) = CS%h_Cu(I,j,k)*ADp%du_dt_visc(I,j,k)
      enddo ; enddo
    enddo
    if (CS%id_hdudtvisc > 0) call post_data(CS%id_hdudtvisc, CS%hdudtvisc, CS%diag)
  endif

  if (ASSOCIATED(CS%hdiffu)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        CS%hdiffu(I,j,k) = CS%h_Cu(I,j,k)*ADp%diffu(I,j,k)
      enddo ; enddo
    enddo
    if (CS%id_hdiffu > 0) call post_data(CS%id_hdiffu, CS%hdiffu, CS%diag)
  endif

  if (ASSOCIATED(CS%hmfu)) then
    do k=1,nz
      do j=Jsq,Jeq ; do i=is,ie
        CS%hmfu(i,J,k) = CS%h_Cv(i,J,k)*(ADp%CAv(i,J,k) - ADp%gradKEv(i,J,k) &
          - ADP%rv_x_u(i,J,k))
      enddo ; enddo
    enddo
    if (CS%id_hmfu > 0) call post_data(CS%id_hmfu, CS%hmfu, CS%diag)
  endif
  
  if (ASSOCIATED(CS%hpfv)) then
    do k=1,nz
      do J=Jsq,Jeq ; do i=is,ie
        CS%hpfv(i,J,k) = CS%h_Cv(i,J,k)*ADp%PFv(i,J,k)
      enddo ; enddo
    enddo
    if (CS%id_hpfv > 0) call post_data(CS%id_hpfv, CS%hpfv, CS%diag)
  endif

  if (ASSOCIATED(CS%hvwb)) then
    do k=1,nz
      do J=Jsq,Jeq ; do i=is,ie
        dwd = 0.5*(CDp%diapyc_vel(i,j,k+1) - CDp%diapyc_vel(i,j,k) &
          + CDp%diapyc_vel(i,j+1,k+1) - CDp%diapyc_vel(i,j+1,k))*ishqlarge(I,J,k)  
        CS%hvwb(i,J,k) = -CS%h_Cv(i,J,k)*ADp%dv_dt_dia(i,J,k) - dwd*v(i,J,k)
      enddo ; enddo
    enddo
    if (CS%id_hvwb > 0) call post_data(CS%id_hvwb, CS%hvwb, CS%diag)
  endif
  
  if (ASSOCIATED(CS%huvxpt)) then
    do k=1,nz
      do J=Jsq,Jeq ; do i=is,ie
        uhxCv = 0.5*((CDp%uh(i,j,k) - CDp%uh(i-1,j,k))*G%IareaT(i,j) &
                + (CDp%uh(i,j+1,k) - CDp%uh(i-1,j+1,k))*G%IareaT(i,j+1))*ishqlarge(I,J,k) 
        CS%huvxpt(i,J,k) = CS%h_Cv(i,J,k)*ADp%rv_x_u(i,J,k) - uhxCv*v(i,J,k)
      enddo ; enddo
    enddo
    if (CS%id_huvxpt > 0) call post_data(CS%id_huvxpt, CS%huvxpt, CS%diag)
  endif
  
  if (ASSOCIATED(CS%hvvymt)) then
    do k=1,nz
      do J=Jsq,Jeq ; do i=is,ie
        vhyCv = 0.5*((CDp%vh(i,J,k) - CDp%vh(i,J-1,k))*G%IareaT(i,j) &
                + (CDp%vh(i,J+1,k) - CDp%vh(i,J,k))*G%IareaT(i,j+1)) *ishqlarge(I,J,k)
        CS%hvvymt(i,J,k) = CS%h_Cv(i,J,k)*ADp%gradKEv(i,J,k) - vhyCv*v(i,J,k)
      enddo ; enddo
    enddo
    if (CS%id_hvvymt > 0) call post_data(CS%id_hvvymt, CS%hvvymt, CS%diag)
  endif

  if (ASSOCIATED(CS%hdvdtvisc)) then
    do k=1,nz
      do J=Jsq,Jeq ; do i=is,ie
        CS%hdvdtvisc(I,j,k) = CS%h_Cv(i,J,k)*ADp%dv_dt_visc(i,J,k)
      enddo ; enddo
    enddo
    if (CS%id_hdvdtvisc > 0) call post_data(CS%id_hdvdtvisc, CS%hdvdtvisc, CS%diag)
  endif

  if (ASSOCIATED(CS%hdiffv)) then
    do k=1,nz
      do J=Jsq,Jeq ; do i=is,ie
        CS%hdiffv(i,J,k) = CS%h_Cv(i,J,k)*ADp%diffv(i,J,k)
      enddo ; enddo
    enddo
    if (CS%id_hdiffv > 0) call post_data(CS%id_hdiffv, CS%hdiffv, CS%diag)
  endif
end subroutine calculate_twa_diagnostics

subroutine register_time_deriv(f_ptr, deriv_ptr, CS)
  real, dimension(:,:,:), target :: f_ptr, deriv_ptr
  type(diagnostics_CS),  pointer :: CS

! This subroutine registers fields to calculate a diagnostic time derivative.
! Arguments: 
!  (target)  f_ptr     - field whose derivative is taken
!  (in)      deriv_ptr - field in which the calculated time derivatives placed
!  (in)      num_lay   - number of layers in this field
!  (in)      CS        - control structure returned by previous call to diagnostics_init

  integer :: m

  if (.not.associated(CS)) call MOM_error(FATAL, &
         "register_time_deriv: Module must be initialized before it is used.")

  if (CS%num_time_deriv >= MAX_FIELDS_) then
    call MOM_error(WARNING,"MOM_diagnostics:  Attempted to register more than " // &
                   "MAX_FIELDS_ diagnostic time derivatives via register_time_deriv.")
    return
  endif

  m = CS%num_time_deriv+1 ; CS%num_time_deriv = m

  CS%nlay(m) = size(f_ptr(:,:,:),3)
  CS%deriv(m)%p => deriv_ptr
  allocate(CS%prev_val(m)%p(size(f_ptr(:,:,:),1), size(f_ptr(:,:,:),2), CS%nlay(m)) )

  CS%var_ptr(m)%p => f_ptr
  CS%prev_val(m)%p(:,:,:) = f_ptr(:,:,:)

end subroutine register_time_deriv


subroutine calculate_derivs(dt, G, CS)
  real,                  intent(in)    :: dt
  type(ocean_grid_type), intent(inout) :: G
  type(diagnostics_CS),  intent(inout) :: CS

! This subroutine calculates all registered time derivatives.
! Arguments: 
!  (in)      dt - time interval in s over which differences occur
!  (in)      G  - ocean grid structure.
!  (in)      CS - control structure returned by previous call to diagnostics_init

  integer i, j, k, m
  real Idt

  if (dt > 0.0) then ; Idt = 1.0/dt
  else ; return ; endif

  do m=1,CS%num_time_deriv
    do k=1,CS%nlay(m) ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      CS%deriv(m)%p(i,j,k) = (CS%var_ptr(m)%p(i,j,k) - CS%prev_val(m)%p(i,j,k)) * Idt
      CS%prev_val(m)%p(i,j,k) = CS%var_ptr(m)%p(i,j,k)
    enddo ; enddo ; enddo
  enddo

end subroutine calculate_derivs

subroutine MOM_diagnostics_init(MIS, ADp, CDp, Time, G, GV, param_file, diag, CS, &
                                wave_speed_CSp)
  type(ocean_internal_state), intent(in)    :: MIS
  type(accel_diag_ptrs),      intent(inout) :: ADp
  type(cont_diag_ptrs),       intent(inout) :: CDp
  type(time_type),            intent(in)    :: Time
  type(ocean_grid_type),      intent(in)    :: G
  type(verticalGrid_type),    intent(in)    :: GV
  type(param_file_type),      intent(in)    :: param_file
  type(diag_ctrl), target,    intent(inout) :: diag
  type(diagnostics_CS),       pointer       :: CS
  type(wave_speed_CS),        pointer       :: wave_speed_CSp

! Arguments
!  (in)     MIS    - For "MOM Internal State" a set of pointers to the fields and
!                    accelerations that make up the ocean's internal physical
!                    state.
!  (inout)  ADp    - structure with pointers to momentum equation terms 
!  (inout)  CDp    - structure with pointers to continuity equation terms
!  (in)     Time   - current model time
!  (in)     G      - ocean grid structure
!  (in)     GV     - The ocean's vertical grid structure.
!  (in) param_file - structure indicating the open file to parse for
!                     model parameter values
!  (in)     diag   - structure to regulate diagnostic output
!  (in/out) CS     - pointer set to point to control structure for this module

! This include declares and sets the variable "version".
#include "version_variable.h"

  character(len=40)  :: mod = "MOM_diagnostics" ! This module's name.
  real :: omega, f2_min
  character(len=48) :: thickness_units, flux_units
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz, nkml, nkbl
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, i, j

  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "MOM_diagnostics_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "SPLIT", CS%split, &
                 "Use the split time stepping if true.", default=.true.)

  if (GV%Boussinesq) then
    thickness_units = "meter" ; flux_units = "meter3 second-1"
  else
    thickness_units = "kilogram meter-2" ; flux_units = "kilogram second-1"
  endif

  CS%id_temp_layer_ave = register_diag_field('ocean_model', 'temp_layer_ave', diag%axesZL, Time, &
      'Layer Average Ocean Temperature', 'Celsius')

  CS%id_salt_layer_ave = register_diag_field('ocean_model', 'salt_layer_ave', diag%axesZL, Time, &
      'Layer Average Ocean Salinity', 'ppt')

  CS%id_masscello = register_diag_field('ocean_model', 'masscello', diag%axesTL,&
      Time, 'Mass per unit area of liquid ocean grid cell', 'kg m-2',           &
      standard_name='sea_water_mass_per_unit_area')

  CS%id_masso = register_scalar_field('ocean_model', 'masso', Time,  &
      diag, 'Mass of liquid ocean', 'kg', standard_name='sea_water_mass')

  CS%id_thkcello = register_diag_field('ocean_model', 'thkcello', diag%axesTL, Time, &
      long_name = 'Cell Thickness', standard_name='cell_thickness', units='m')

  if (((CS%id_masscello>0) .or. (CS%id_masso>0) .or. (CS%id_thkcello>0.and..not.GV%Boussinesq)) &
      .and. .not.ASSOCIATED(CS%diag_tmp3d)) then
    call safe_alloc_ptr(CS%diag_tmp3d,isd,ied,jsd,jed,nz)
  endif

  CS%id_thetaoga = register_scalar_field('ocean_model', 'thetaoga',    &
      Time, diag, 'Global Mean Ocean Potential Temperature', 'Celsius',&
      standard_name='sea_water_potential_temperature')

  CS%id_soga = register_scalar_field('ocean_model', 'soga', &
      Time, diag, 'Global Mean Ocean Salinity', 'ppt',      &
      standard_name='sea_water_salinity')

  CS%id_tosga = register_scalar_field('ocean_model', 'sst_global', Time, diag,&
      long_name='Global Area Average Sea Surface Temperature',                &
      units='degC', standard_name='sea_surface_temperature',                  &   
      cmor_field_name='tosga', cmor_standard_name='sea_surface_temperature',  &
      cmor_units='degC', cmor_long_name='Sea Surface Temperature')

  CS%id_sosga = register_scalar_field('ocean_model', 'sss_global', Time, diag,&
      long_name='Global Area Average Sea Surface Salinity',                   &
      units='ppt', standard_name='sea_surface_salinity',                      &    
      cmor_field_name='sosga', cmor_standard_name='sea_surface_salinity',     &
      cmor_units='ppt', cmor_long_name='Sea Surface Salinity')

  CS%id_e = register_diag_field('ocean_model', 'e', diag%axesTi, Time, &
      'Interface Height Relative to Mean Sea Level', 'meter')
  if (CS%id_e>0) call safe_alloc_ptr(CS%e,isd,ied,jsd,jed,nz+1)

  CS%id_e_D = register_diag_field('ocean_model', 'e_D', diag%axesTi, Time, &
      'Interface Height above the Seafloor', 'meter')
  if (CS%id_e_D>0) call safe_alloc_ptr(CS%e_D,isd,ied,jsd,jed,nz+1)

  CS%id_Rml = register_diag_field('ocean_model', 'Rml', diag%axesTL, Time, &
      'Mixed Layer Coordinate Potential Density', 'kg meter-3')

  CS%id_Rcv = register_diag_field('ocean_model', 'Rho_cv', diag%axesTL, Time, &
      'Coordinate Potential Density', 'kg meter-3')

  CS%id_rhopot0 = register_diag_field('ocean_model', 'rhopot0', diag%axesTL, Time, &
      'Potential density referenced to surface', 'kg meter-3')
  CS%id_rhopot2 = register_diag_field('ocean_model', 'rhopot2', diag%axesTL, Time, &
      'Potential density referenced to 2000 dbar', 'kg meter-3')
  CS%id_rhoinsitu = register_diag_field('ocean_model', 'rhoinsitu', diag%axesTL, Time, &
      'In situ density', 'kg meter-3')

  CS%id_du_dt = register_diag_field('ocean_model', 'dudt', diag%axesCuL, Time, &
      'Zonal Acceleration', 'meter second-2')
  if ((CS%id_du_dt>0) .and. .not.ASSOCIATED(CS%du_dt)) then
    call safe_alloc_ptr(CS%du_dt,IsdB,IedB,jsd,jed,nz)
    call register_time_deriv(MIS%u, CS%du_dt, CS)
  endif

  CS%id_dv_dt = register_diag_field('ocean_model', 'dvdt', diag%axesCvL, Time, &
      'Meridional Acceleration', 'meter second-2')
  if ((CS%id_dv_dt>0) .and. .not.ASSOCIATED(CS%dv_dt)) then
    call safe_alloc_ptr(CS%dv_dt,isd,ied,JsdB,JedB,nz)
    call register_time_deriv(MIS%v, CS%dv_dt, CS)
  endif

  CS%id_dh_dt = register_diag_field('ocean_model', 'dhdt', diag%axesTL, Time, &
      'Thickness tendency', trim(thickness_units)//" second-1")
  if ((CS%id_dh_dt>0) .and. .not.ASSOCIATED(CS%dh_dt)) then
    call safe_alloc_ptr(CS%dh_dt,isd,ied,jsd,jed,nz)
    call register_time_deriv(MIS%h, CS%dh_dt, CS)
  endif

  ! layer thickness variables 
  !if (GV%nk_rho_varies > 0) then
    CS%id_h_Rlay = register_diag_field('ocean_model', 'h_rho', diag%axesTL, Time, &
        'Layer thicknesses in pure potential density coordinates', thickness_units)
    if (CS%id_h_Rlay>0) call safe_alloc_ptr(CS%h_Rlay,isd,ied,jsd,jed,nz)

    CS%id_uh_Rlay = register_diag_field('ocean_model', 'uh_rho', diag%axesCuL, Time, &
        'Zonal volume transport in pure potential density coordinates', flux_units)
    if (CS%id_uh_Rlay>0) call safe_alloc_ptr(CS%uh_Rlay,IsdB,IedB,jsd,jed,nz)

    CS%id_vh_Rlay = register_diag_field('ocean_model', 'vh_rho', diag%axesCvL, Time, &
        'Meridional volume transport in pure potential density coordinates', flux_units)
    if (CS%id_vh_Rlay>0) call safe_alloc_ptr(CS%vh_Rlay,isd,ied,JsdB,JedB,nz)

    CS%id_uhGM_Rlay = register_diag_field('ocean_model', 'uhGM_rho', diag%axesCuL, Time, &
        'Zonal volume transport due to interface height diffusion in pure potential &
        &density coordinates', flux_units)
    if (CS%id_uhGM_Rlay>0) call safe_alloc_ptr(CS%uhGM_Rlay,IsdB,IedB,jsd,jed,nz)

    CS%id_vhGM_Rlay = register_diag_field('ocean_model', 'vhGM_rho', diag%axesCvL, Time, &
        'Meridional volume transport due to interface height diffusion in pure &
        &potential density coordinates', flux_units)
    if (CS%id_vhGM_Rlay>0) call safe_alloc_ptr(CS%vhGM_Rlay,isd,ied,JsdB,JedB,nz)
  !endif


  ! terms in the kinetic energy budget
  CS%id_KE = register_diag_field('ocean_model', 'KE', diag%axesTL, Time, &
      'Layer kinetic energy per unit mass', 'meter2 second-2')
  if (CS%id_KE>0) call safe_alloc_ptr(CS%KE,isd,ied,jsd,jed,nz)

  CS%id_dKEdt = register_diag_field('ocean_model', 'dKE_dt', diag%axesTL, Time, &
      'Kinetic Energy Tendency of Layer', 'meter3 second-3')
  if (CS%id_dKEdt>0) call safe_alloc_ptr(CS%dKE_dt,isd,ied,jsd,jed,nz)

  CS%id_PE_to_KE = register_diag_field('ocean_model', 'PE_to_KE', diag%axesTL, Time, &
      'Potential to Kinetic Energy Conversion of Layer', 'meter3 second-3')
  if (CS%id_PE_to_KE>0) call safe_alloc_ptr(CS%PE_to_KE,isd,ied,jsd,jed,nz)

  CS%id_KE_Coradv = register_diag_field('ocean_model', 'KE_Coradv', diag%axesTL, Time, &
      'Kinetic Energy Source from Coriolis and Advection', 'meter3 second-3')
  if (CS%id_KE_Coradv>0) call safe_alloc_ptr(CS%KE_Coradv,isd,ied,jsd,jed,nz)

  CS%id_KE_adv = register_diag_field('ocean_model', 'KE_adv', diag%axesTL, Time, &
      'Kinetic Energy Source from Advection', 'meter3 second-3')
  if (CS%id_KE_adv>0) call safe_alloc_ptr(CS%KE_adv,isd,ied,jsd,jed,nz)

  CS%id_KE_visc = register_diag_field('ocean_model', 'KE_visc', diag%axesTL, Time, &
      'Kinetic Energy Source from Vertical Viscosity and Stresses', 'meter3 second-3')
  if (CS%id_KE_visc>0) call safe_alloc_ptr(CS%KE_visc,isd,ied,jsd,jed,nz)

  CS%id_KE_horvisc = register_diag_field('ocean_model', 'KE_horvisc', diag%axesTL, Time, &
      'Kinetic Energy Source from Horizontal Viscosity', 'meter3 second-3')
  if (CS%id_KE_horvisc>0) call safe_alloc_ptr(CS%KE_horvisc,isd,ied,jsd,jed,nz)

  CS%id_KE_dia = register_diag_field('ocean_model', 'KE_dia', diag%axesTL, Time, &
      'Kinetic Energy Source from Diapycnal Diffusion', 'meter3 second-3')
  if (CS%id_KE_dia>0) call safe_alloc_ptr(CS%KE_dia,isd,ied,jsd,jed,nz)


  ! gravity wave CFLs 
  CS%id_cg1 = register_diag_field('ocean_model', 'cg1', diag%axesT1, Time, &
      'First baroclinic gravity wave speed', 'meter second-1')
  CS%id_Rd1 = register_diag_field('ocean_model', 'Rd1', diag%axesT1, Time, &
      'First baroclinic deformation radius', 'meter')
  CS%id_cfl_cg1 = register_diag_field('ocean_model', 'CFL_cg1', diag%axesT1, Time, &
      'CFL of first baroclinic gravity wave = dt*cg1*(1/dx+1/dy)', 'nondim')
  CS%id_cfl_cg1_x = register_diag_field('ocean_model', 'CFL_cg1_x', diag%axesT1, Time, &
      'i-component of CFL of first baroclinic gravity wave = dt*cg1*/dx', 'nondim')
  CS%id_cfl_cg1_y = register_diag_field('ocean_model', 'CFL_cg1_y', diag%axesT1, Time, &
      'j-component of CFL of first baroclinic gravity wave = dt*cg1*/dy', 'nondim')

  CS%wave_speed_CSp => wave_speed_CSp
  if ((CS%id_cg1>0) .or. (CS%id_Rd1>0) .or. (CS%id_cfl_cg1>0) .or. &
      (CS%id_cfl_cg1_x>0) .or. (CS%id_cfl_cg1_y>0)) then
    call safe_alloc_ptr(CS%cg1,isd,ied,jsd,jed)
    if (CS%id_Rd1>0)       call safe_alloc_ptr(CS%Rd1,isd,ied,jsd,jed)
    if (CS%id_cfl_cg1>0)   call safe_alloc_ptr(CS%cfl_cg1,isd,ied,jsd,jed)
    if (CS%id_cfl_cg1_x>0) call safe_alloc_ptr(CS%cfl_cg1_x,isd,ied,jsd,jed)
    if (CS%id_cfl_cg1_y>0) call safe_alloc_ptr(CS%cfl_cg1_y,isd,ied,jsd,jed)
  endif

  CS%id_mass_wt = register_diag_field('ocean_model', 'mass_wt', diag%axesT1, Time,  &
      'The column mass for calculating mass-weighted average properties', 'kg m-2')

  CS%id_temp_int = register_diag_field('ocean_model', 'temp_int', diag%axesT1, Time,                &
      'Density weighted column integrated potential temperature', 'degC kg m-2',                    &
      cmor_field_name='opottempmint',                                                               &
      cmor_long_name='integral_wrt_depth_of_product_of_sea_water_density_and_potential_temperature',&
      cmor_units='degC kg m-2', cmor_standard_name='Depth integrated density times potential temperature')

  CS%id_salt_int = register_diag_field('ocean_model', 'salt_int', diag%axesT1, Time,   &
      'Density weighted column integrated salinity', 'ppt kg m-2',                     &
      cmor_field_name='somint',                                                        &
      cmor_long_name='integral_wrt_depth_of_product_of_sea_water_density_and_salinity',&
      cmor_units='ppt kg m-2', cmor_standard_name='Depth integrated density times salinity')

  CS%id_col_mass = register_diag_field('ocean_model', 'col_mass', diag%axesT1, Time, &
      'The column integrated in situ density', 'kg m-2')

  CS%id_col_ht = register_diag_field('ocean_model', 'col_height', diag%axesT1, Time, &
      'The height of the water column', 'm')
  CS%id_pbo = register_diag_field('ocean_model', 'pbo', diag%axesT1, Time, &
      long_name='Sea Water Pressure at Sea Floor', standard_name='sea_water_pressure_at_sea_floor', &
      units='Pa')

  ! terms in the twa budget
  CS%id_hfv = register_diag_field('ocean_model', 'twa_hfv', diag%axesCuL, Time, &
      'coriolis xtwa term', 'meter second-2')
  call safe_alloc_ptr(CS%hfv,IsdB,IedB,jsd,jed,nz)
  CS%id_hpfu = register_diag_field('ocean_model', 'twa_hpfu', diag%axesCuL, Time, &
      'PG xtwa term', 'meter second-2')
  call safe_alloc_ptr(CS%hpfu,IsdB,IedB,jsd,jed,nz)
  CS%id_huwb = register_diag_field('ocean_model', 'twa_huwb', diag%axesCuL, Time, &
      'diabatic xtwa term', 'meter second-2')
  call safe_alloc_ptr(CS%huwb,IsdB,IedB,jsd,jed,nz)
  CS%id_huuxpt = register_diag_field('ocean_model', 'twa_huuxpt', diag%axesCuL, Time, &
      'zonal advection xtwa term', 'meter second-2')
  call safe_alloc_ptr(CS%huuxpt,IsdB,IedB,jsd,jed,nz)
  CS%id_huvymt = register_diag_field('ocean_model', 'twa_huvymt', diag%axesCuL, Time, &
      'meridional advection xtwa term', 'meter second-2')
  call safe_alloc_ptr(CS%huvymt,IsdB,IedB,jsd,jed,nz)
  CS%id_hdudtvisc = register_diag_field('ocean_model', 'twa_hdudtvisc', diag%axesCuL, Time, &
      'verical viscous xtwa term', 'meter second-2')
  call safe_alloc_ptr(CS%hdudtvisc,IsdB,IedB,jsd,jed,nz)
  CS%id_hdiffu = register_diag_field('ocean_model', 'twa_hdiffu', diag%axesCuL, Time, &
      'horizontal viscous xtwa term', 'meter second-2')
  call safe_alloc_ptr(CS%hdiffu,IsdB,IedB,jsd,jed,nz)
  CS%id_hmfu = register_diag_field('ocean_model', 'twa_hmfu', diag%axesCvL, Time, &
      'coriolis ytwa term', 'meter second-2')
  call safe_alloc_ptr(CS%hmfu,isd,ied,JsdB,JedB,nz)
  CS%id_hpfv = register_diag_field('ocean_model', 'twa_hpfv', diag%axesCvL, Time, &
      'PG ytwa term', 'meter second-2')
  call safe_alloc_ptr(CS%hpfv,isd,ied,JsdB,JedB,nz)
  CS%id_hvwb = register_diag_field('ocean_model', 'twa_hvwb', diag%axesCvL, Time, &
      'diabatic ytwa term', 'meter second-2')
  call safe_alloc_ptr(CS%hvwb,isd,ied,JsdB,JedB,nz)
  CS%id_huvxpt = register_diag_field('ocean_model', 'twa_huvxpt', diag%axesCvL, Time, &
      'zonal advection ytwa term', 'meter second-2')
  call safe_alloc_ptr(CS%huvxpt,isd,ied,JsdB,JedB,nz)
  CS%id_hvvymt = register_diag_field('ocean_model', 'twa_hvvymt', diag%axesCvL, Time, &
      'meridional advection ytwa term', 'meter second-2')
  call safe_alloc_ptr(CS%hvvymt,isd,ied,JsdB,JedB,nz)
  CS%id_hdvdtvisc = register_diag_field('ocean_model', 'twa_hdvdtvisc', diag%axesCvL, Time, &
      'vertical viscous ytwa term', 'meter second-2')
  call safe_alloc_ptr(CS%hdvdtvisc,isd,ied,JsdB,JedB,nz)
  CS%id_hdiffv = register_diag_field('ocean_model', 'twa_hdiffv', diag%axesCvL, Time, &
      'horizontal viscous ytwa term', 'meter second-2')
  call safe_alloc_ptr(CS%hdiffv,isd,ied,JsdB,JedB,nz)

  CS%id_h_Cu = register_diag_field('ocean_model', 'h_Cu', diag%axesCuL, Time, &
      'h at Cu points', 'meter')
  call safe_alloc_ptr(CS%h_Cu,IsdB,IedB,jsd,jed,nz)
  CS%id_h_Cv = register_diag_field('ocean_model', 'h_Cv', diag%axesCvL, Time, &
      'h at Cv points', 'meter')
  call safe_alloc_ptr(CS%h_Cv,isd,ied,JsdB,JedB,nz)
  CS%id_huu_Cu = register_diag_field('ocean_model', 'huu_Cu', diag%axesCuL, Time, &
      'huu at Cu points', 'meter3 second-1')
  call safe_alloc_ptr(CS%huu_Cu,IsdB,IedB,jsd,jed,nz)
  CS%id_hvv_Cv = register_diag_field('ocean_model', 'hvv_Cv', diag%axesCvL, Time, &
      'hvv at Cv points', 'meter3 second-1')
  call safe_alloc_ptr(CS%hvv_Cv,isd,ied,JsdB,JedB,nz)
  CS%id_hv_Cu = register_diag_field('ocean_model', 'hv_Cu', diag%axesCuL, Time, &
      'hv at Cu points', 'meter2 second-1')
  call safe_alloc_ptr(CS%hv_Cu,IsdB,IedB,jsd,jed,nz)
  CS%id_hu_Cv = register_diag_field('ocean_model', 'hu_Cv', diag%axesCvL, Time, &
      'hu at Cv points', 'meter')
  call safe_alloc_ptr(CS%hu_Cv,isd,ied,JsdB,JedB,nz)
  CS%id_hw_Cu = register_diag_field('ocean_model', 'hw_Cu', diag%axesCuL, Time, &
      'hw at Cu points', 'meter2 second-1')
  call safe_alloc_ptr(CS%hw_Cu,IsdB,IedB,jsd,jed,nz)
  CS%id_hw_Cv = register_diag_field('ocean_model', 'hw_Cv', diag%axesCvL, Time, &
      'hw at Cv points', 'meter2 second-1')
  call safe_alloc_ptr(CS%hw_Cv,isd,ied,JsdB,JedB,nz)
  CS%id_hwb_Cu = register_diag_field('ocean_model', 'hwb_Cu', diag%axesCuL, Time, &
      'hwb at Cu points', 'meter2 second-1')
  call safe_alloc_ptr(CS%hwb_Cu,IsdB,IedB,jsd,jed,nz)
  CS%id_hwb_Cv = register_diag_field('ocean_model', 'hwb_Cv', diag%axesCvL, Time, &
      'hwb at Cv points', 'meter2 second-1')
  call safe_alloc_ptr(CS%hwb_Cv,isd,ied,JsdB,JedB,nz)
  CS%id_esq_Cu = register_diag_field('ocean_model', 'esq_Cu', diag%axesCuL, Time, &
      'e**2 at Cu points', 'meter2')
  call safe_alloc_ptr(CS%esq_Cu,IsdB,IedB,jsd,jed,nz)
  CS%id_esq_Cv = register_diag_field('ocean_model', 'esq_Cv', diag%axesCvL, Time, &
      'e**2 at Cv points', 'meter2')
  call safe_alloc_ptr(CS%esq_Cv,isd,ied,JsdB,JedB,nz)
  CS%id_e_Cu = register_diag_field('ocean_model', 'e_Cu', diag%axesCuL, Time, &
      'e at Cu points', 'meter')
  call safe_alloc_ptr(CS%e_Cu,IsdB,IedB,jsd,jed,nz)
  CS%id_e_Cv = register_diag_field('ocean_model', 'e_Cv', diag%axesCvL, Time, &
      'e at Cv points', 'meter')
  call safe_alloc_ptr(CS%e_Cv,isd,ied,JsdB,JedB,nz)
  CS%id_epfu = register_diag_field('ocean_model', 'epfu', diag%axesCuL, Time, &
      'epfu at Cu points', 'meter2 second-1')
  call safe_alloc_ptr(CS%epfu,IsdB,IedB,jsd,jed,nz)
  CS%id_epfv = register_diag_field('ocean_model', 'epfv', diag%axesCvL, Time, &
      'epfv at Cv points', 'meter2 second-1')
  call safe_alloc_ptr(CS%epfv,isd,ied,JsdB,JedB,nz)
  CS%id_pfu_masked = register_diag_field('ocean_model', 'pfu_masked', diag%axesCuL, Time, &
      'pfu_masked at Cu points', 'meter2 second-1')
  call safe_alloc_ptr(CS%pfu_masked,IsdB,IedB,jsd,jed,nz)
  CS%id_pfv_masked = register_diag_field('ocean_model', 'pfv_masked', diag%axesCvL, Time, &
      'pfv_masked at Cv points', 'meter2 second-1')
  call safe_alloc_ptr(CS%pfv_masked,isd,ied,JsdB,JedB,nz)


  call set_dependent_diagnostics(MIS, ADp, CDp, G, CS)

end subroutine MOM_diagnostics_init


subroutine set_dependent_diagnostics(MIS, ADp, CDp, G, CS)
  type(ocean_internal_state), intent(in)    :: MIS
  type(accel_diag_ptrs),      intent(inout) :: ADp
  type(cont_diag_ptrs),       intent(inout) :: CDp
  type(ocean_grid_type),      intent(in)    :: G
  type(diagnostics_CS),       pointer       :: CS

! This subroutine sets up diagnostics upon which other diagnostics depend.
! Arguments: 
!  (in)      MIS - For "MOM Internal State" a set of pointers to the fields and
!                  accelerations making up ocean internal physical state.
!  (inout)   ADp - structure pointing to accelerations in momentum equation
!  (inout)   CDp - structure pointing to terms in continuity equation
!  (in)      G   - ocean grid structure
!  (in)      CS -  pointer to the control structure for this module

  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (ASSOCIATED(CS%dKE_dt) .or. ASSOCIATED(CS%PE_to_KE) .or. &
      ASSOCIATED(CS%KE_CorAdv) .or. ASSOCIATED(CS%KE_adv) .or. &
      ASSOCIATED(CS%KE_visc) .or. ASSOCIATED(CS%KE_horvisc) .or. &
      ASSOCIATED(CS%KE_dia)) &
    call safe_alloc_ptr(CS%KE,isd,ied,jsd,jed,nz)

  if (ASSOCIATED(CS%dKE_dt)) then
    if (.not.ASSOCIATED(CS%du_dt)) then
      call safe_alloc_ptr(CS%du_dt,IsdB,IedB,jsd,jed,nz)
      call register_time_deriv(MIS%u, CS%du_dt, CS)
    endif
    if (.not.ASSOCIATED(CS%dv_dt)) then
      call safe_alloc_ptr(CS%dv_dt,isd,ied,JsdB,JedB,nz)
      call register_time_deriv(MIS%v, CS%dv_dt, CS)
    endif
    if (.not.ASSOCIATED(CS%dh_dt)) then
      call safe_alloc_ptr(CS%dh_dt,isd,ied,jsd,jed,nz)
      call register_time_deriv(MIS%h, CS%dh_dt, CS)
    endif
  endif

  if (ASSOCIATED(CS%KE_adv)) then
    call safe_alloc_ptr(ADp%gradKEu,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%gradKEv,isd,ied,JsdB,JedB,nz)
  endif

  if (ASSOCIATED(CS%KE_visc)) then
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
  endif

  if (ASSOCIATED(CS%KE_dia)) then
    call safe_alloc_ptr(ADp%du_dt_dia,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%dv_dt_dia,isd,ied,JsdB,JedB,nz)
  endif

  if (ASSOCIATED(CS%uhGM_Rlay)) call safe_alloc_ptr(CDp%uhGM,IsdB,IedB,jsd,jed,nz)
  if (ASSOCIATED(CS%vhGM_Rlay)) call safe_alloc_ptr(CDp%vhGM,isd,ied,JsdB,JedB,nz)

  if (ASSOCIATED(CS%hfv)) then
    call safe_alloc_ptr(CS%h_Cu,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%gradKEu,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%rv_x_v,IsdB,IedB,jsd,jed,nz)
  endif

  if (ASSOCIATED(CS%huwb)) then
    call safe_alloc_ptr(CS%h_Cu,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(CDp%diapyc_vel,isd,ied,jsd,jed,nz+1)
    call safe_alloc_ptr(ADp%du_dt_dia,IsdB,IedB,jsd,jed,nz)
  endif

  if (ASSOCIATED(CS%huuxpt)) then
    call safe_alloc_ptr(CS%h_Cu,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%gradKEu,IsdB,IedB,jsd,jed,nz)
  endif

  if (ASSOCIATED(CS%huvymt)) then
    call safe_alloc_ptr(CS%h_Cu,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%rv_x_v,IsdB,IedB,jsd,jed,nz)
  endif

  if (ASSOCIATED(CS%hdudtvisc)) then
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(CS%h_Cu,IsdB,IedB,jsd,jed,nz)
  endif

  if (ASSOCIATED(CS%hdiffu)) then
    call safe_alloc_ptr(ADp%diffu,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(CS%h_Cu,IsdB,IedB,jsd,jed,nz)
  endif

  if (ASSOCIATED(CS%hmfu)) then
    call safe_alloc_ptr(CS%h_Cv,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%gradKEv,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%rv_x_u,isd,ied,JsdB,JedB,nz)
  endif

  if (ASSOCIATED(CS%hvwb)) then
    call safe_alloc_ptr(CS%h_Cv,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(CDp%diapyc_vel,isd,ied,jsd,jed,nz+1)
    call safe_alloc_ptr(ADp%dv_dt_dia,isd,ied,JsdB,JedB,nz)
  endif

  if (ASSOCIATED(CS%huvxpt)) then
    call safe_alloc_ptr(CS%h_Cv,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%rv_x_u,isd,ied,JsdB,JedB,nz)
  endif

  if (ASSOCIATED(CS%hvvymt)) then
    call safe_alloc_ptr(CS%h_Cv,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%gradKEv,isd,ied,JsdB,JedB,nz)
  endif

  if (ASSOCIATED(CS%hdvdtvisc)) then
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(CS%h_Cv,isd,ied,JsdB,JedB,nz)
  endif

  if (ASSOCIATED(CS%hdiffv)) then
    call safe_alloc_ptr(ADp%diffv,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(CS%h_Cv,isd,ied,JsdB,JedB,nz)
  endif

  if (ASSOCIATED(CS%hw_Cu)) then
    call safe_alloc_ptr(CS%h_Cu,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(CDp%diapyc_vel,isd,ied,jsd,jed,nz+1)
  endif

  if (ASSOCIATED(CS%hw_Cv)) then
    call safe_alloc_ptr(CS%h_Cv,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(CDp%diapyc_vel,isd,ied,jsd,jed,nz+1)
  endif

  if (ASSOCIATED(CS%hwb_Cu)) then
    call safe_alloc_ptr(CDp%diapyc_vel,isd,ied,jsd,jed,nz+1)
  endif

  if (ASSOCIATED(CS%hwb_Cv)) then
    call safe_alloc_ptr(CDp%diapyc_vel,isd,ied,jsd,jed,nz+1)
  endif

  if (ASSOCIATED(CS%epfu)) then
    call safe_alloc_ptr(CS%e_Cu,IsdB,IedB,jsd,jed,nz)
  endif

  if (ASSOCIATED(CS%epfv)) then
    call safe_alloc_ptr(CS%e_Cv,isd,ied,JsdB,JedB,nz)
  endif

  if (ASSOCIATED(CS%e_Cv) .OR. ASSOCIATED(CS%e_Cu)&
    .OR. ASSOCIATED(CS%esq_Cu) .OR. ASSOCIATED(CS%esq_Cv)) then
    call safe_alloc_ptr(CS%e,isd,ied,jsd,jed,nz+1)
  endif
end subroutine set_dependent_diagnostics


subroutine MOM_diagnostics_end(CS, ADp)
  type(diagnostics_CS),   pointer       :: CS
  type(accel_diag_ptrs),  intent(inout) :: ADp
  integer :: m
    
  if (ASSOCIATED(CS%e))          deallocate(CS%e)
  if (ASSOCIATED(CS%e_D))        deallocate(CS%e_D)
  if (ASSOCIATED(CS%KE))         deallocate(CS%KE)
  if (ASSOCIATED(CS%dKE_dt))     deallocate(CS%dKE_dt)
  if (ASSOCIATED(CS%PE_to_KE))   deallocate(CS%PE_to_KE)
  if (ASSOCIATED(CS%KE_Coradv))  deallocate(CS%KE_Coradv)
  if (ASSOCIATED(CS%KE_adv))     deallocate(CS%KE_adv)
  if (ASSOCIATED(CS%KE_visc))    deallocate(CS%KE_visc)
  if (ASSOCIATED(CS%KE_horvisc)) deallocate(CS%KE_horvisc)
  if (ASSOCIATED(CS%KE_dia))     deallocate(CS%KE_dia)
  if (ASSOCIATED(CS%dv_dt))      deallocate(CS%dv_dt)
  if (ASSOCIATED(CS%dh_dt))      deallocate(CS%dh_dt)
  if (ASSOCIATED(CS%du_dt))      deallocate(CS%du_dt)
  if (ASSOCIATED(CS%h_Rlay))     deallocate(CS%h_Rlay)
  if (ASSOCIATED(CS%uh_Rlay))    deallocate(CS%uh_Rlay)
  if (ASSOCIATED(CS%vh_Rlay))    deallocate(CS%vh_Rlay)
  if (ASSOCIATED(CS%uhGM_Rlay))  deallocate(CS%uhGM_Rlay)
  if (ASSOCIATED(CS%vhGM_Rlay))  deallocate(CS%vhGM_Rlay)
  if (ASSOCIATED(CS%diag_tmp3d)) deallocate(CS%diag_tmp3d)

  if (ASSOCIATED(ADp%gradKEu))    deallocate(ADp%gradKEu)
  if (ASSOCIATED(ADp%gradKEu))    deallocate(ADp%gradKEu)
  if (ASSOCIATED(ADp%du_dt_visc)) deallocate(ADp%du_dt_visc)
  if (ASSOCIATED(ADp%dv_dt_visc)) deallocate(ADp%dv_dt_visc)
  if (ASSOCIATED(ADp%du_dt_dia))  deallocate(ADp%du_dt_dia)
  if (ASSOCIATED(ADp%dv_dt_dia))  deallocate(ADp%dv_dt_dia)
  if (ASSOCIATED(ADp%du_other))   deallocate(ADp%du_other)
  if (ASSOCIATED(ADp%dv_other))   deallocate(ADp%dv_other)

  if (ASSOCIATED(CS%hfv))         deallocate(CS%hfv)
  if (ASSOCIATED(CS%hpfu))        deallocate(CS%hpfu)
  if (ASSOCIATED(CS%huwb))        deallocate(CS%huwb)
  if (ASSOCIATED(CS%huuxpt))      deallocate(CS%huuxpt)
  if (ASSOCIATED(CS%huvymt))      deallocate(CS%huvymt)
  if (ASSOCIATED(CS%hdudtvisc))   deallocate(CS%hdudtvisc)
  if (ASSOCIATED(CS%hdiffu))      deallocate(CS%hdiffu)
  if (ASSOCIATED(CS%hmfu))        deallocate(CS%hmfu)
  if (ASSOCIATED(CS%hpfv))        deallocate(CS%hpfv)
  if (ASSOCIATED(CS%hvwb))        deallocate(CS%hvwb)
  if (ASSOCIATED(CS%huvxpt))      deallocate(CS%huvxpt)
  if (ASSOCIATED(CS%hvvymt))      deallocate(CS%hvvymt)
  if (ASSOCIATED(CS%hdvdtvisc))   deallocate(CS%hdvdtvisc)
  if (ASSOCIATED(CS%hdiffv))      deallocate(CS%hdiffv)

  if (ASSOCIATED(CS%h_Cu))        deallocate(CS%h_Cu)
  if (ASSOCIATED(CS%huu_Cu))       deallocate(CS%huu_Cu)
  if (ASSOCIATED(CS%hv_Cu))       deallocate(CS%hv_Cu)
  if (ASSOCIATED(CS%hw_Cu))       deallocate(CS%hw_Cu)
  if (ASSOCIATED(CS%hwb_Cu))      deallocate(CS%hwb_Cu)

  if (ASSOCIATED(CS%h_Cv))        deallocate(CS%h_Cv)
  if (ASSOCIATED(CS%hvv_Cv))       deallocate(CS%hvv_Cv)
  if (ASSOCIATED(CS%hu_Cv))       deallocate(CS%hu_Cv)
  if (ASSOCIATED(CS%hw_Cv))       deallocate(CS%hw_Cv)
  if (ASSOCIATED(CS%hwb_Cv))      deallocate(CS%hwb_Cv)

  do m=1,CS%num_time_deriv ; deallocate(CS%prev_val(m)%p) ; enddo

  deallocate(CS)

end subroutine MOM_diagnostics_end

end module MOM_diagnostics
