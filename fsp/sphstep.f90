!> Subroutines for advancing the simulation in time.
!! ODE solution logic intertwined with updates of smoothed
!! particle properties.
!!
!! Leapfrog algorithm modified from:
!! Yasmin Melean, Leonardo Di G Sigalotti, and Anwar Hasmy,
!! “On the SPH tensile instability in forming viscous liquid drops,” 
!! Computer Physics Communications 157, no. 3 (2004): 191-200.
!!
!! \todo add sound speed calculations

 !if(rho >= 2) then 
 !  if(unstable_density .eq. 0) then
 !    continue
 !  else if(unstable_density .eq. 1) then
 !    print *,rho
 !    stop 'density in unstable region'
 !  else if(unstable_density .eq. 2) then
 !    t = (u+adash*rho)/kbdash
 !    p = 0.0
 !  endif
 !else

module sphstep

use global
use kernel
use particle
use reader
use sphforce
use sphforce_ab
use density
use simulation_box
use neighbour_list
use thermostat
use boundary
use collision
use splib

implicit none

contains

!> Aggregates the call to find pair seperations,
!! compress list, and compute dv, to allow integrations to be
!! expressed more compactly
subroutine compute_pairs(nl,box,p,bnl,bp)
  type (nlist_data_type), intent(inout) :: nl
  type (box_data_type) :: box
  type (particle_type), intent(inout) :: p
  type (nlist_data_type), optional, dimension(:), intent(inout) :: bnl
  type (particle_type),optional, dimension(:), intent(inout) :: bp
  integer i,ni
  i=0

  if(present(bp))then
    do i=1,nbounds
      call find_pair_separations(bnl(i),box,p%x,bp(i)%x)
      call compress_nlist(bnl(i))
      ni = bnl(i)%nip
      call calc_dv_ab(bnl(i)%dv(1:ni,:),bnl(i)%ilist(1:ni,:),p%v,bp(i)%v &
        & ,p%n,bnl(i)%nip,dimn)
    end do
  endif
  
  call find_pair_separations(nl,box,p%x)
  call compress_nlist(nl)
  ni = nl%nip
  call calc_dv(nl%dv(1:ni,:),nl%ilist(1:ni,:),p%v,p%n,ni,dimn)

  if(validate) call validate_particles(p)
  if(validate) call validate_nlist(nl)
  if(verbose) print *,'Found pair seperations and compressed...'
end subroutine compute_pairs


!> Aggregates the function calls for computation of densities for particles.
subroutine compute_density(nl,p,bnl,bp)
  type (nlist_data_type), intent(in) :: nl
  type (particle_type), intent(inout) :: p
  type (nlist_data_type), optional, dimension(:), intent(inout) :: bnl
  type (particle_type),optional, dimension(:), intent(inout) :: bp
  integer i
  i=0
  
  call sum_density(nl,p%n,p%sml,p%m,nl%w, p%rho,p%grad_rho,nl%dwdx)
  call sum_density(nl,p%n,p%sml_long,p%m,nl%w_long,p%rho_lr,p%grad_rho &
  &               ,nl%dwdx_long)

  if (present(bp)) then
    do i=1,nbounds
      call boundary_density(bnl(i),p,bp(i))
    end do
  endif

  if(validate) call validate_particles(p)
  if(validate) call validate_nlist(nl)
  if(verbose) print *,'Densities computed...'
end subroutine compute_density


!> Executes one simulation step,
!! using a modified leapfrog method 
subroutine sph_step_leapfrog(nl,p,dt,box,bnl,bp)
  implicit none

  type (nlist_data_type), intent(inout) :: nl
  type (particle_type), intent(inout) :: p
  double precision :: dt
  type (box_data_type) :: box
  type (particle_type),optional, dimension(nbounds), intent(inout) :: bp
  type (nlist_data_type),optional, dimension(nbounds), intent(inout) :: bnl
  
  double precision, dimension(p%n) :: tdsdt
  integer :: i,j,d,n

  !> Amount moved by each particle this step.
  double precision, dimension(p%n,p%dimn) :: dr
  double precision :: dedt_tot
  double precision :: u_env_half

  !> Temporary particle structure for leapfrog method.
  type (particle_type) :: p_cent
 
  n = p%n
 
  p_cent = p

  if(verbose) print *, 'taking another step'
  
  if(validate) then
      call validate_particles(p)
      call validate_particles(p_cent)
      call validate_nlist(nl)
      if(present(bp)) then
        do i=1,nbounds
          call validate_particles(bp(i))
        enddo 
      endif
  endif

 
  ! LEAPFROG HALF STEP
  ! -----------------------
  p%xhalf = p%xhalf + p%rdot * dt
  p_cent%x = p%xhalf
  p_cent%v = p%v + 0.5 * p%a * dt
  p_cent%u = p%u + 0.5 * p%dedt * dt

  ! Apply boundary conditions
  ! ------------------------- 
  call apply_boundary(p_cent%x,p_cent%v,box)
  call apply_pbcs(box, p_cent%x)

  ! Compute pair seperatations
  ! --------------------------
  if(present(bp)) then
    call compute_pairs(nl,box,p_cent,bnl,bp)
  else
    call compute_pairs(nl,box,p_cent)
  endif
  if(validate) call validate_nlist(nl)


  ! Half step Collisions
  ! ---------- 
  if(collisions .eq. 1) then
    call collide_particles(nl,p_cent)
  endif 

  p_cent%rdot(:,:) = p_cent%v(:,:) 
  if(velocity_avg .eq. 1) call v_smooth(p_cent,nl)

  ! Half step kernels
  ! ---------------------------------
  if(verbose) print *,'Calculating kernel for all pairs'
  call calc_kernels_all(p_cent,nl)
  if(present(bnl)) then
    do i=1,nbounds
      call calc_kernels_ab(p_cent,bp(i),bnl(i))
    end do
  endif

  ! Half step densities
  ! -------------------
  if(verbose) print *,'Kernel calculation finished - calculating densities'
  call compute_density(nl,p_cent,bnl,bp)

 
  ! Half step temperatures and thermostat
  ! -------------------------------------
  if(verbose) print *,'Densities ok, updating temperatures'

  do i=1,n
    call update_temperature(p_cent%u(i),p_cent%temp(i),p_cent%rho(i), &
    & p_cent%grad_rho(i,:) )
  end do
 
  if(thermostat_particles) then 
    call apply_thermostat(n,p_cent%temp,p%thermostat,p_cent%u, &
    & p_cent%rho,p_cent%grad_rho,u_env_half,thermostat_type)
  endif

  ! Half step force
  ! NB. calc_sphforce only changes the particle acceleration and dedt
  ! ---------------
  if (eos_type .eq. 1 .or. eos_type .eq. 2 .or. eos_type .eq. 4) then
    call calc_sphforce(nl,n,dimn,nl%drij,p_cent,nl%dwdx,tdsdt)  
  else if (eos_type .eq. 3) then
    if(verbose) print *,'Calculating forces for cohesive smoothing'
    call calc_sphforce(nl,n,dimn,nl%drij,p_cent,nl%dwdx,tdsdt,nl%dwdx_long) 
  end if

  if (present(bp)) then
    do i=1,nbounds
      if(verbose) print *,'Applying sph force from boundary particles'
      call calc_sphforce_ab(bnl(i),p_cent,bp(i),tdsdt,dedt_tot)
      call calc_sphforce_boundary(bnl(i),p_cent)
      if(validate) call validate_particles(bp(i))
    end do
  endif
    
 
  ! Half step body force
  ! --------------------
  if(debug) print *,'Gravitation.'
  call apply_gravity(p_cent%a)
  if(validate) call validate_particles(p_cent)

  ! Half step smooth velocities
  ! ---------------------------
  if(debug) print *,'Velocity Smoothing'
  p_cent%rdot(:,:) = p_cent%v(:,:) 
  if(velocity_avg .eq. 1) call v_smooth(p_cent,nl)
  if(debug) print *,'Velocity Smoothing Finished'


  ! Advance with centred rates of change
  ! Advance x to xhalf and use this for collisions etc
  ! ------------------------------------
  do i=1,n
    do d=1, dimn
      dr(i,d) = p_cent%rdot(i,d)*dt
      p%x(i,d) = p%xhalf(i,d) + dr(i,d)
      p%v(i,d) = p%v(i,d) + p_cent%a(i,d)*dt        
    end do
    p%dedt(i) = p_cent%dedt(i)
    p%u(i) = p%u(i) + p_cent%dedt(i)*dt
  end do
  p%a = p_cent%a
  if(debug) print *,'Advanced...'

  ! Apply boundaries
  ! -------------------------------
  if(debug) print *,'Boundaries...'
  call apply_boundary(p%x,p%v,box)
  call apply_pbcs(box, p%x)
  if(debug) print *,'Boundaried...'

  ! Increment neighbour list and compute smoothed properties
  ! --------------------------------------------------------
  call increment_nlist_drtot(nl,dr)
  if(debug) print *,'Incremented...'

  ! Reform neighbour list if required
  ! ---------------------------------
  if(reform_nlist_now(nl)) then
    print *,'Reforming neighbour list - this may be expensive.'
    call form_nlist(nl, box, p%x)
  endif  
  
  if(present(bnl)) then
    do i=1,nbounds
      if(reform_nlist_now(bnl(i))) then
        call form_nlist(bnl(i), box, p%x,bp(i)%x)
      endif
    end do
  endif

  ! Recompute pairs
  ! -----------------------------------
  call compute_pairs(nl,box,p,bnl,bp)
  if(validate) call validate_nlist(nl)

  ! Collisions
  ! ----------
  if(collisions .eq. 1) then
    call collide_particles(nl,p)
  endif 

  ! Kernels 
  ! -------
  call calc_kernels_all(p,nl)
  if(present(bnl)) then
    do i=1,nbounds
      call calc_kernels_ab(p,bp(i),bnl(i))
    end do
  endif

  ! Density
  ! --------
  call compute_density(nl,p,bnl,bp)

  p%rdot(:,:) = p%v(:,:) 
  if(velocity_avg .eq. 1) call v_smooth(p,nl)

  ! Update temperatures with new energy
  ! -----------------------------------
  if(verbose) print *,'Updating internal energy and temperature'
  do i=1,n
    call update_temperature(p%u(i),p%temp(i),p%rho(i),p%grad_rho(i,:))
  end do
  
  if(present(bp)) then
    do i=1,nbounds
      do j=1,bp(i)%n
        bp(i)%u(j) = bp(i)%u(j) + bp(i)%dedt(j)*dt
        call update_temperature(bp(i)%u(j),bp(i)%temp(j),bp(i)%rho(j), &
        & bp(i)%grad_rho(j,:))
      enddo
    enddo
  endif

  if(thermostat_boundaries) then 
    do i=1,nbounds
      call apply_thermostat(bp(i)%n,bp(i)%temp,bp(i)%thermostat, &
      & bp(i)%u,bp(i)%rho,bp(i)%grad_rho,u_env,2)
      if(validate) call validate_particles(bp(i))
      call increment_nlist_drtot(bnl(i),dr)
    end do
  endif

  if(thermostat_particles) then  
    call apply_thermostat(n,p%temp,p%thermostat,p%u,p%rho,p%grad_rho &
    & ,u_env,thermostat_type)
  endif

  ! Compute smoothed properties
  ! ---------------------------
  call calc_smoothed_properties(nl,p,nl%w)

  if(verbose) print *, 'Leapfrog step finished - Hooray!'
           
end subroutine sph_step_leapfrog


!> This subroutine executes one simulation step using an 'improved euler' 
!! method.
!! This may also be known as a second order Runge-Kutta predictor corrector. 
!! It is quite possible no more accurate than leapfrog, and is still here for 
!! error checking and historical purposes.
subroutine sph_step_imp_euler(nl,p,dt,box,bnl,bp)
  implicit none

  type (nlist_data_type), intent(inout) :: nl
  type (particle_type), intent(inout) :: p
  integer :: n
  double precision :: dt
  type (box_data_type) :: box
  type (particle_type),optional, dimension(nbounds), intent(inout) :: bp
  type (nlist_data_type),optional, dimension(nbounds), intent(inout) :: bnl 

  double precision, dimension(p%n) :: tdsdt
  
  !> Array to hold the amount moved this step.
  double precision, dimension(p%n,p%dimn) :: dr
  integer :: i,j,d

  ! Temporary particle object for improved euler method.
  type (particle_type) :: ptmp
  double precision :: dedt_tot
  double precision :: u_env_tmp
  
  n = p%n
  if(verbose) print *, 'Taking another step'
 
  ! PREDICTOR STEP
  ! --------------

  ! Predictor force calculation
  ! ---------------------------
  ! nb. calc_sphforce only changes 
  ! the particle acceleration and dedt
  if (eos_type .eq. 1 .or. eos_type .eq. 2 .or. eos_type .eq. 4) then
    call calc_sphforce(nl,n,dimn,nl%drij,p,nl%dwdx,&
    & tdsdt)
  else if (eos_type .eq. 3) then
    call calc_sphforce(nl,n,dimn,nl%drij,p,nl%dwdx, &
    & tdsdt,nl%dwdx_long)
  end if

  if (present(bp)) then
    do i=1,nbounds
      if(verbose) print *,'Applying sph force from boundary particles'
      call calc_sphforce_ab(bnl(i),p,bp(i),tdsdt,dedt_tot)
      call calc_sphforce_boundary(bnl(i),p)
      if(validate) call validate_particles(bp(i))
    end do
  endif 

  ! Predictor body force
  ! --------------------
  call apply_gravity(p%a)

  ! Predictor smooth v
  ! ------------------
  if(debug) print *,'calling v_smooth'
  p%rdot(:,:) = p%v(:,:) 
  if(velocity_avg .eq. 1) call v_smooth(p,nl)
  if(validate) call validate_particles(p)

  ptmp = p
  if(validate) call validate_particles(ptmp)
 
  ! Predictor project state
  ! -----------------------
  do i=1,n
    do d=1, dimn
      ptmp%x(i,d) = p%x(i,d) + p%rdot(i,d)*dt
      ptmp%v(i,d) = p%v(i,d) + p%a(i,d)*dt
    end do
    ptmp%u(i) = p%u(i) + p%dedt(i)*dt
  end do
  
  if(verbose) print *,'Predictor step ended'
  if(validate) call validate_particles(ptmp)


  ! Projected Boundaries
  ! --------------------
  call apply_boundary(ptmp%x,ptmp%v,box)
  call apply_pbcs(box,ptmp%x)

  ! Projected Pairs
  ! ---------------
  call compute_pairs(nl,box,ptmp,bnl,bp)

  ! Projected Collisions
  ! --------------------
  if(collisions .eq. 1) then
    call collide_particles(nl,ptmp)
  endif

  ! Projected Kernels
  ! ------------------
  call calc_kernels_all(ptmp,nl)
  if(present(bnl)) then
    do i=1,nbounds
      call calc_kernels_ab(ptmp,bp(i),bnl(i))
    end do
  endif 
  if(validate) call validate_particles(ptmp)

  ! Projected density estimate 
  ! --------------------------
  !call sum_density(nl,n,ptmp%sml,ptmp%m,nl%w,ptmp%rho,ptmp%grad_rho, &
  !&                nl%dwdx)
  call compute_density(nl,ptmp,bnl,bp)

  ! Projected temperature
  ! ---------------------
  do i=1,n
    call update_temperature(ptmp%u(i),ptmp%temp(i),ptmp%rho(i),&
    &    ptmp%grad_rho(i,:))
  end do
   
  if(thermostat_particles) then  
    call apply_thermostat(n,ptmp%temp,ptmp%thermostat,ptmp%u, &
    & ptmp%rho,ptmp%grad_rho &
    & ,u_env_tmp,thermostat_type)
  endif
 
  if(verbose) print *,'Temperatures set'
  if(validate) call validate_particles(ptmp)

  ! Projected Forces
  ! ----------------
  if (eos_type .eq. 1 .or. eos_type .eq. 2 .or. eos_type .eq. 4) then
    call calc_sphforce(nl,n,dimn,nl%drij,ptmp,nl%dwdx,tdsdt) 
  else if (eos_type .eq. 3) then
    call calc_sphforce(nl,n,dimn,nl%drij,ptmp,nl%dwdx,tdsdt,nl%dwdx_long)
  end if

  if (present(bp)) then
    do i=1,nbounds
      if(verbose) print *,'Applying sph force from boundary particles'
      call calc_sphforce_ab(bnl(i),ptmp,bp(i),tdsdt,dedt_tot)
      call calc_sphforce_boundary(bnl(i),ptmp)
      if(validate) call validate_particles(bp(i))
    end do
  endif

  ! Projected body force
  ! --------------------
  call apply_gravity(ptmp%a)

  if(verbose) print *,'Forces calculated'
  if(validate) call validate_particles(ptmp)

  ! Projected velocity smoothing
  ! ----------------------------
  ptmp%rdot(:,:) = ptmp%v(:,:) 
  if(velocity_avg .eq. 1) call v_smooth(ptmp,nl)
  if(validate) call validate_particles(ptmp)
  if(verbose) print *,'Set rdot, stepping forward in time...'

  ! Corrector timestep
  ! ------------------
  do i=1,n
    do d=1, dimn
      dr(i,d) = (0.5) * ( p%rdot(i,d) + ptmp%rdot(i,d) )*dt
      p%x(i,d) = p%x(i,d) + dr(i,d)
      p%v(i,d) = p%v(i,d) + (0.5) * ( p%a(i,d) + ptmp%a(i,d) )*dt
    end do
    p%u(i) = p%u(i) + (0.5) * ( p%dedt(i) + ptmp%dedt(i) )*dt
  end do
  
  if(validate) call validate_particles(p)
   

  ! Finalise boundary conditions
  ! ----------------------------  
  call apply_pbcs(box,p%x)
  call apply_boundary(p%x,p%v,box)
  call increment_nlist_drtot(nl,dr)


  ! Construct neighbour list for the next step
  ! ------------------------------------------
  if(reform_nlist_now(nl)) then
    print *,'Forming nlist - costing time'
    call form_nlist(nl, box, p%x)
  endif

  if(present(bnl)) then
    do i=1,nbounds
      if(reform_nlist_now(bnl(i))) then
        call form_nlist(bnl(i), box, p%x,bp(i)%x)
      endif
    end do
  endif 

  if(verbose) print *, 'Formed neighbour list'
  if(validate) call validate_particles(p)


  ! Finalise pair separations and minimum image
  ! ---------------------------------------
  call compute_pairs(nl,box,p,bnl,bp)
 

  ! Finalise collisions
  ! --------------------
  if(collisions .eq. 1) call collide_particles(nl,p)
  

  ! Finalise kernels and density
  ! -----------------------------
  call calc_kernels_all(p,nl)
   if(present(bnl)) then
    do i=1,nbounds
      call calc_kernels_ab(p,bp(i),bnl(i))
    end do
  endif  
  call compute_density(nl,p,bnl,bp)
  if(validate) call validate_particles(p)


  ! Finalise temperature
  ! ---------------------
  if(verbose) print *,'Updating temperatures'
  do i=1,n
    call update_temperature(p%u(i),p%temp(i),p%rho(i),p%grad_rho(i,:))
  end do

  if(present(bp)) then
    do i=1,nbounds
      do j=1,bp(i)%n
        bp(i)%u(j) = bp(i)%u(j) + bp(i)%dedt(j)*dt
        call update_temperature(bp(i)%u(j),bp(i)%temp(j),bp(i)%rho(j), &
        & bp(i)%grad_rho(j,:))
      enddo
    enddo
  endif

  if(thermostat_particles) then  
    call apply_thermostat(n,p%temp,p%thermostat,p%u,p%rho,p%grad_rho &
    & ,u_env,thermostat_type)
  endif

  if(thermostat_boundaries) then 
    do i=1,nbounds
        call apply_thermostat(bp(i)%n,bp(i)%temp,bp(i)%thermostat, &
        & bp(i)%u,bp(i)%rho,bp(i)%grad_rho,u_env,2)
        if(validate) call validate_particles(bp(i))
      call increment_nlist_drtot(bnl(i),dr)
    end do
  endif

  if(validate) call validate_particles(p)
  if(validate) call validate_nlist(nl)

  if(verbose) print *,'Updated temperatures'

  call calc_smoothed_properties(nl,p,nl%w)

  if(verbose) print *, 'Improved Euler step finished. Bravo!'
end subroutine sph_step_imp_euler


!> Executes a simulation step using a fourth order Runge-Kutta method. 
!! Someone remarked to me once that this is overkill because SPH is only 
!! second order accurate, but I wasn't convinced. 
!! Bill Hoover reckons RK4 is
!! the right method to use for SPAM, and he wrote the book.
subroutine sph_step_runge_kutta(nl,p,dt,box,bnl,bp)
  implicit none

  ! this subroutine executes one simulation step.
  type (nlist_data_type), intent(inout) :: nl
  type (particle_type), intent(inout) :: p
  integer :: n
  double precision :: dt
  !double precision, dimension(:) :: w
  !double precision, intent(inout), dimension(:,:) :: dwdx
  type (box_data_type) :: box
  double precision, dimension(p%n) :: tdsdt
  
  type (particle_type),optional, dimension(nbounds), intent(inout) :: bp
  type (nlist_data_type),optional, dimension(nbounds), intent(inout) :: bnl

  !double precision, dimension(n) :: dedt
  integer :: i,j,d

  double precision, dimension(p%n,p%dimn) :: dr
  double precision :: dedt_tot
  double precision :: u_env_tmp

  !> temporary particle for runge kutta method
  type (particle_type) :: ptmp
  
  !> temporary arrays for rates of change for RK method.
  double precision, dimension(p%n,p%dimn) :: dxk1, dxk2, dxk3, dxk4
  double precision, dimension(p%n,p%dimn) :: dvk1, dvk2, dvk3, dvk4
  double precision, dimension(p%n) ::  duk1, duk2, duk3, duk4
    
  ptmp = p
  n = p%n

  if(verbose) print *, 'Taking another RK4 step...'
  

  ! K1 STEP
  ! -------

  ! K1 Forces
  ! ---------
  if (eos_type .eq. 3) then
      call calc_sphforce(nl,n,dimn,nl%drij,p,nl%dwdx, &
                       & tdsdt,nl%dwdx_long)
  else if (eos_type .eq. 1 .or. eos_type .eq. 2 .or. eos_type .eq. 4) then
      call calc_sphforce(nl,n,dimn,nl%drij,p,nl%dwdx,&
                       & tdsdt)
  end if

  if (present(bp)) then
    do i=1,nbounds
      if(verbose) print *,'Applying sph force from boundary particles'
      call calc_sphforce_ab(bnl(i),p,bp(i),tdsdt,dedt_tot)
      call calc_sphforce_boundary(bnl(i),p)
      if(validate) call validate_particles(bp(i))
    end do
  endif
  
  call apply_gravity(p%a)


  ! K1 velocity smoothing
  ! ---------------------
  p%rdot(:,:) = p%v(:,:) 
  if(velocity_avg .eq. 1) call v_smooth(p,nl)

 
  ! K1 step forward
  ! ---------------
  do i=1,n
    do d=1, dimn
      dxk1(i,d) = p%rdot(i,d)*dt
      dvk1(i,d) = p%a(i,d)*dt
    end do
    duk1(i) = p%dedt(i)*dt
  end do

  ! END K1 STEP
  ! -----------

  
  ! K2 STEP
  ! -------
  do i=1,n
    do d=1,dimn
      ptmp%x(i,d) = p%x(i,d) + dxk1(i,d)/2.0
      ptmp%v(i,d) = p%v(i,d) + dvk1(i,d)/2.0
    end do
    ptmp%u(i) = p%u(i) + duk1(i)/2.0
  end do

  call apply_pbcs(box, ptmp%x)
  call apply_boundary(ptmp%x,ptmp%v,box)


  call compute_pairs(nl,box,ptmp,bnl,bp)

  if(collisions .eq. 1) then
    call collide_particles(nl,ptmp)
  endif 

  call calc_kernels_all(ptmp,nl)
   if(present(bnl)) then
    do i=1,nbounds
      call calc_kernels_ab(ptmp,bp(i),bnl(i))
    end do
  endif
 
  ! K2 Density estimate based on predictor position
  ! -----------------------------------------------
  call compute_density(nl,ptmp,bnl,bp)
  !call sum_density(nl, n, ptmp%sml, ptmp%m, nl%w, ptmp%rho,ptmp%grad_rho,nl%dwdx)

  ! K2 update temperature
  ! ---------------------
  do i=1,n
    call update_temperature(ptmp%u(i),ptmp%temp(i),ptmp%rho(i) &
    &   ,ptmp%grad_rho(i,:))
  end do

  if(thermostat_particles) then  
    call apply_thermostat(n,ptmp%temp,ptmp%thermostat,ptmp%u, &
    & ptmp%rho,ptmp%grad_rho &
    & ,u_env_tmp,thermostat_type)
  endif

  ! K2 Forces
  ! ---------
  if (eos_type .eq. 1 .or. eos_type .eq. 2 .or. eos_type .eq. 4) then
    call calc_sphforce(nl,n,dimn,nl%drij,ptmp,nl%dwdx,& 
                & tdsdt)
  else if (eos_type .eq. 3) then
    call calc_sphforce(nl,n,dimn,nl%drij,ptmp,nl%dwdx,&
            & tdsdt,nl%dwdx_long)
  end if

  if (present(bp)) then
    do i=1,nbounds
      if(verbose) print *,'Applying sph force from boundary particles'
      call calc_sphforce_ab(bnl(i),ptmp,bp(i),tdsdt,dedt_tot)
      call calc_sphforce_boundary(bnl(i),ptmp)
      if(validate) call validate_particles(bp(i))
    end do
  endif

  call apply_gravity(ptmp%a)

  ! K2 velocity smoothing
  ! ---------------------
  ptmp%rdot(:,:) = ptmp%v(:,:)
  if(velocity_avg .eq. 1) call v_smooth(ptmp,nl)
 

  ! K2 step forward
  ! ---------------
  do i=1,n
    do d=1, dimn
      dxk2(i,d) = ptmp%rdot(i,d) * dt !/ 2.0
      dvk2(i,d) = ptmp%a(i,d) * dt !/ 2.0
    end do
    !update internal energy
    duk2(i) = ptmp%dedt(i)*dt !/2.0
  end do

  ! END K2 STEP
  ! -----------

  
  ! K3 STEP
  ! -------

  do i=1,n
    do d=1,dimn
      ptmp%x(i,d) = p%x(i,d) + dxk2(i,d) / 2.0
      ptmp%v(i,d) = p%v(i,d) + dvk2(i,d) / 2.0
    end do
    ptmp%u(i) = p%u(i) + duk2(i) / 2.0
  end do

  call apply_pbcs(box, ptmp%x)
  call apply_boundary(ptmp%x,ptmp%v,box)

 
  ! K3 Pairs
  ! -------- 
  call compute_pairs(nl,box,ptmp,bnl,bp)

  if(collisions .eq. 1) then
    call collide_particles(nl,ptmp)
  endif 

  call calc_kernels_all(ptmp,nl)
   if(present(bnl)) then
    do i=1,nbounds
      call calc_kernels_ab(ptmp,bp(i),bnl(i))
    end do
  endif
 
  ! K3 density estimate based on predictor position
  ! -----------------------------------------------
  call compute_density(nl,ptmp,bnl,bp)
  !call sum_density(nl, n, ptmp%sml, ptmp%m, nl%w, ptmp%rho,ptmp%grad_rho,nl%dwdx)

  ! K3 update temperature
  ! ---------------------
  do i=1,n
     call update_temperature(ptmp%u(i),ptmp%temp(i),ptmp%rho(i),ptmp%grad_rho(i,:))
  end do

  if(thermostat_particles) then  
    call apply_thermostat(n,ptmp%temp,ptmp%thermostat,ptmp%u, &
    & ptmp%rho,ptmp%grad_rho &
    & ,u_env_tmp,thermostat_type)
  endif

  ! K3 Forces
  ! ---------
  if (eos_type .eq. 1 .or. eos_type .eq. 2 .or. eos_type .eq. 4) then
     call calc_sphforce(nl,n,dimn,nl%drij,ptmp,nl%dwdx,& 
            & tdsdt)
  else if (eos_type .eq. 3) then
      call calc_sphforce(nl,n,dimn,nl%drij,ptmp,nl%dwdx,&
            & tdsdt,nl%dwdx_long)
  end if

  if (present(bp)) then
    do i=1,nbounds
      if(verbose) print *,'Applying sph force from boundary particles'
      call calc_sphforce_ab(bnl(i),ptmp,bp(i),tdsdt,dedt_tot)
      call calc_sphforce_boundary(bnl(i),ptmp)
      if(validate) call validate_particles(bp(i))
    end do
  endif

  call apply_gravity(ptmp%a)

  ! K3 velocity smoothing
  ! ---------------------
  ptmp%rdot(:,:) = ptmp%v(:,:) 
  if(velocity_avg .eq. 1) call v_smooth(ptmp,nl)
 
  ! K3 step forward
  ! ---------------
  do i=1,n
    do d=1, dimn
      dxk3(i,d) = ptmp%rdot(i,d)*dt !/2.0
      dvk3(i,d) = ptmp%a(i,d)*dt !/2.0
    end do
    !update internal energy
    duk3(i) = ptmp%dedt(i)*dt!/2.0
  end do

  ! END K3 STEP
  ! -----------
 
 
  ! K4 STEP 
  ! -------
  do i=1,n
    do d=1,dimn
      ptmp%x(i,d) = p%x(i,d) + dxk3(i,d)
      ptmp%v(i,d) = p%v(i,d) + dvk3(i,d)
    end do
    ptmp%u(i) = p%u(i) + duk3(i)
  end do

  call apply_pbcs(box, ptmp%x)
  call apply_boundary(ptmp%x,ptmp%v,box)

  call compute_pairs(nl,box,ptmp,bnl,bp)

  if(collisions .eq. 1) then
    call collide_particles(nl,ptmp)
  endif 

  call calc_kernels_all(ptmp,nl)
   if(present(bnl)) then
    do i=1,nbounds
      call calc_kernels_ab(ptmp,bp(i),bnl(i))
    end do
  endif
 
  ! K4 density estimate based on predictor position
  ! -----------------------------------------------
  call compute_density(nl,ptmp,bnl,bp)
  !call sum_density(nl, n, ptmp%sml, ptmp%m, nl%w, ptmp%rho,ptmp%grad_rho,nl%dwdx)

  ! K4 update temperature
  ! ---------------------
  do i=1,n
    call update_temperature(ptmp%u(i),ptmp%temp(i),ptmp%rho(i), &
    & ptmp%grad_rho(i,:))
  end do

  if(thermostat_particles) then  
    call apply_thermostat(n,ptmp%temp,ptmp%thermostat,ptmp%u, &
    & ptmp%rho,ptmp%grad_rho &
    & ,u_env_tmp,thermostat_type)
  endif

  ! K4 Forces
  ! ---------
  if (eos_type .eq. 1 .or. eos_type .eq. 2 .or. eos_type .eq. 4) then
    call calc_sphforce(nl,n,dimn,nl%drij,ptmp,nl%dwdx,& 
            & tdsdt)
  else if (eos_type .eq. 3) then
    call calc_sphforce(nl,n,dimn,nl%drij,ptmp,nl%dwdx,&
            & tdsdt,nl%dwdx_long)
  end if

  if (present(bp)) then
    do i=1,nbounds
      if(verbose) print *,'Applying sph force from boundary particles'
      call calc_sphforce_ab(bnl(i),ptmp,bp(i),tdsdt,dedt_tot)
      call calc_sphforce_boundary(bnl(i),ptmp)
      if(validate) call validate_particles(bp(i))
    end do
  endif
  
  call apply_gravity(ptmp%a)

  ptmp%rdot(:,:) = ptmp%v(:,:) 
  if(velocity_avg .eq. 1) call v_smooth(ptmp,nl)
  
  do i=1,n
    do d=1, dimn
      dxk4(i,d) = ptmp%rdot(i,d)*dt
      dvk4(i,d) = ptmp%a(i,d)*dt
    end do
    !update internal energy
    duk4(i) = ptmp%dedt(i)*dt
  end do
  ! END K4 STEP


  ! FINAL STEP
  do i=1,n
    do d=1, dimn
      dr(i,d) = (1.0/6.0)*dxk1(i,d) + (1.0/3.0)*dxk2(i,d) + (1.0/3.0)*dxk3(i,d) + (1.0/6.0)*dxk4(i,d)
      p%x(i,d) = p%x(i,d) + dr(i,d)
      p%v(i,d) = p%v(i,d) + (1.0/6.0)*dvk1(i,d) + (1.0/3.0)*dvk2(i,d) + (1.0/3.0)*dvk3(i,d) + (1.0/6.0)*dvk4(i,d) 
    end do
    
    p%u(i) = p%u(i) + (1.0/6.0)*duk1(i) + (1.0/3.0)*duk2(i) + (1.0/3.0)*duk3(i) + (1.0/6.0)*duk4(i)       
    
  end do

  ! END FINAL STEP
  ! --------------
  
  call increment_nlist_drtot(nl,dr)


  ! Finalise properties - reform nlist if required.
  ! -----------------------------------------------

  if(reform_nlist_now(nl)) then
    print *,'Reforming neighbour list'
    call form_nlist(nl, box, p%x)
  endif

  if(present(bnl)) then
    do i=1,nbounds
      if(reform_nlist_now(bnl(i))) then
        call form_nlist(bnl(i), box, p%x,bp(i)%x)
      endif
    end do
  endif 

  call apply_pbcs(box, p%x)
  call apply_boundary(p%x,p%v,box)

  ! Find pair separations and minimum image
  call compute_pairs(nl,box,p,bnl,bp)
  
  if(collisions .eq. 1) then
    call collide_particles(nl,p)
  endif 

  call calc_kernels_all(p,nl)
   if(present(bnl)) then
    do i=1,nbounds
      call calc_kernels_ab(p,bp(i),bnl(i))
    end do
  endif

  ! calculate densities
  call compute_density(nl,p,bnl,bp)
  !call sum_density(nl, n, p%sml, p%m, nl%w, p%rho,p%grad_rho,nl%dwdx)

   do i=1,n
     call update_temperature(p%u(i),p%temp(i),p%rho(i),p%grad_rho(i,:))
   end do 
  
  ! Finalise temperature
  ! ---------------------
  if(verbose) print *,'Updating temperatures'
  do i=1,n
    call update_temperature(p%u(i),p%temp(i),p%rho(i),p%grad_rho(i,:))
  end do

  if(present(bp)) then
    do i=1,nbounds
      do j=1,bp(i)%n
        bp(i)%u(j) = bp(i)%u(j) + bp(i)%dedt(j)*dt
        call update_temperature(bp(i)%u(j),bp(i)%temp(j),bp(i)%rho(j), &
        & bp(i)%grad_rho(j,:))
      enddo
    enddo
  endif

  if(thermostat_particles) then  
    call apply_thermostat(n,p%temp,p%thermostat,p%u,p%rho,p%grad_rho &
    & ,u_env,thermostat_type)
  endif

  if(thermostat_boundaries) then 
    do i=1,nbounds
        call apply_thermostat(bp(i)%n,bp(i)%temp,bp(i)%thermostat, &
        & bp(i)%u,bp(i)%rho,bp(i)%grad_rho,u_env,2)
        if(validate) call validate_particles(bp(i))
      call increment_nlist_drtot(bnl(i),dr)
    end do
  endif

  call calc_smoothed_properties(nl,p,nl%w)

  if(debug) print *, 'RK4 step finished you scabby dog.'

end subroutine sph_step_runge_kutta

end module sphstep

