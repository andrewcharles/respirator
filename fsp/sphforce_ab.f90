!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----|
!> Smooth Particle Mechanics code
!! Modified version of sphforce, designed for two particle systems
!! Copyright Andrew Charles, all rights reserved.


module sphforce_ab

use eos
use global
use neighbour_list
use simulation_box
use art_viscosity
use core_potential
use particle
use sphforce

public :: calc_grad_v_ab, calc_div_v_ab, calc_grad_v_os_ab, calc_pi_os_ab, calc_pi_one_ab

contains
!> Computes SPH rates of change
!! \param nl neighbour list data structure
!! \param box box data structure
!! \param n - number of particles
!! \param dimn - dimension of system - currently supports 1 or 2
!! \param dx - displacements between pairs
!! \param p particle data structure
    !! x,v,a - positions, velocities, accelerations
    !! u - internal energies
    !! m - masses
    !! dedt - rate of change of specific internal energy
    !! rho - densities
    !! sml - smoothing lengths
    !! c   - particle sound speed
    !! eta - dynamic shear viscosity
    !! zeta - dynamic bulk viscosity
    !! temp - temperatures
    !! pco - cohesive pressure for seperated eos
    !! pressure tensor
    !! symmetric traceless viscous (irreversible) pressure tensor
    !! pi_os
    !! pi_irr is the total irreversible pressure tensor
    !! isotropic viscous (irreversible) pressure tensor
    !! double precision, dimension (n) :: pi_one
    !! reversible, isotropic pressure tensor
    !! this accumulates core potential and avisc as well
    !! double precision, dimension(n,dimn,dimn) :: p_rev
    !! the isotropic part of the equilibrium pressure
    !! double precision, dimension(n) :: p_eq
!! \param tdsdt - entropy production
!! \param bn
!! \param dedt_total total change in internal energy of the boundary

subroutine calc_sphforce_ab(nl,pa,pb,tdsdt,dedt_total)
  implicit none

  type (nlist_data_type), intent(inout) :: nl
  !type (box_data_type), intent(inout) :: box
  type (particle_type), intent(inout) :: pa
  type (particle_type), intent(inout) :: pb
  !type (boundary_data_type), intent(in) :: bnd
  double precision :: dedt_total

  integer :: n
  !double precision, dimension(:,:) :: dx
  double precision, dimension(:), intent(inout)::  tdsdt

  ! INTERNAL VARIABLES
  ! symmetric traceless velocity gradient, used to calc pi_os
  double precision, dimension (pa%n,pa%dimn,pa%dimn) :: grad_v_os_a
  double precision, dimension (pb%n,pb%dimn,pb%dimn) :: grad_v_os_b

  double precision, dimension(pa%n) :: div_v_a
  double precision, dimension(pb%n) :: div_v_b

  !double precision, dimension(pa%n) :: gradgradrho_a
  !double precision, dimension(pb%n) :: gradgradrho_b
  !> running totals
  !! A keen observer will note the similarity in name with Liu and Liu's
  !! code and deduce correctly that i owe them a substantial debt.
  double precision :: he, hx, hy

  !> artifical viscosity for each interaction
  double precision :: avisc

  !> core repulsion
  !double precision, dimension(dimn) :: repulsion

  !> f is used in the monaghan heat flux calculation
  double precision :: f

  integer :: i,j,k
  integer :: ni

  !benchmarking
  integer :: start_time,end_time,gradv_time,stress_time,time_rate
  real :: time_seconds
  if(profile) then
      call system_clock(count=start_time,count_rate=time_rate)
  endif

  n = pa%n

  ! initialise tensor arrays to zero
  !pa%pi_os = 0
  !pa%pi_one = 0
  !pa%pi_irr = 0
  !pa%p_rev = 0
  !pa%p_eq = 0
  
  pb%pi_os = 0
  pb%pi_one = 0
  pb%pi_irr = 0
  pb%p_rev = 0
  pb%p_eq = 0

  ! intialise some other stuff to zero
  tdsdt = 0
  avisc = 0.0

  ! set rates of change and spatial gradients to zero
  !pa%a = 0.0
  !pa%q = 0
  !pa%dedt = 0
  !pa%grad_v = 0

  pb%a = 0.0
  pb%q = 0
  pb%dedt = 0
  pb%grad_v = 0

  if(verbose) print *, 'Calculating SPH forces '
  
  if(profile) call system_clock(count=gradv_time)

  call calc_grad_v_ab(nl,pa,pb)
  call calc_div_v_ab(pa%grad_v,div_v_a)
  call calc_div_v_ab(pb%grad_v,div_v_b)
  call calc_grad_v_os_ab(grad_v_os_a,pa%grad_v,div_v_a)
  call calc_grad_v_os_ab(grad_v_os_b,pb%grad_v,div_v_b)

  if(profile) then
    call system_clock(count = end_time)
    time_seconds = (end_time - gradv_time)/real(time_rate)
    print *,'**** div,grad v time ',i, time_seconds
  endif
    
  ! calculate pressure from equation of state
  
  !call calc_iso_pressure_ab(pa)
  !call calc_iso_pressure_ab(pb)

  !pa%p_eq = pa%p 
  !pb%p_eq = pb%p

  !call calc_gradgradrho_ab(nl,pa,pb,gradgradrho_a,gradgradrho_b) 
  
  if(verbose) print *, 'Calculating capillary pressure '
  
  !call calc_capillary_pressure_ab(pa,gradgradrho_a) 
  !call calc_capillary_pressure_ab(pb,gradgradrho_b) 

  !pa%p_rev(:,1,1) = pa%p_rev(:,1,1) + pa%p_eq(:)
  !pa%p_rev(:,2,2) = pa%p_rev(:,2,2) + pa%p_eq(:)

  !pb%p_rev(:,1,1) = pb%p_rev(:,1,1) + pb%p_eq(:)
  !pb%p_rev(:,2,2) = pb%p_rev(:,2,2) + pb%p_eq(:)
  
  if(profile) call system_clock(count=stress_time)
  
  ! calculate shear stress  
  if(verbose) print *, 'Calculating shear stress '
  !call calc_pi_os(pa%pi_os,grad_v_os_a,pa%eta)
  !call calc_pi_os(pb%pi_os,grad_v_os_b,pb%eta)

  !calculate bulk stress
  if(verbose) print *, 'Calculating bulk stress '
  !call calc_pi_one(pa%pi_one,div_v_a,pa%zeta)
  !call calc_pi_one(pb%pi_one,div_v_b,pb%zeta)

  ! Add bulk stress to shear stress
  !pa%pi_irr = pa%pi_os
  !pb%pi_irr = pb%pi_os

  !pa%pi_irr(:,1,1) = pa%pi_irr(:,1,1) + pa%pi_one(:)
  !pb%pi_irr(:,2,2) = pb%pi_irr(:,2,2) + pb%pi_one(:)

  ! calculate viscous entropy
  if(verbose) print *, 'calculating viscous entropy '

  !call calc_viscous_entropy_os(nl,p%pi_irr,n,tdsdt,grad_v_os,p%rho)

  if(profile) then
    call system_clock(count = end_time)
    time_seconds = (end_time - stress_time)/real(time_rate)
    print *,'**** shear,bulk,viscous entropy time ',i, time_seconds
  endif

  if(conduction_type .eq. 2) then
    call calc_heat_flux_ab(nl,pa,pb)
  end if

  ! calc SPH sum for pressure force
  if(verbose) print *, 'calculating pressure force '

  ni = nl%nip

  do k = 1, ni
    i = nl%ilist(k,1)
    j = nl%ilist(k,2)

    if(debug) print *, 'loop for pair',k

    ! hx and hy are running totals of P/rho**2
    hx=0.0
    hy=0.0
    !he is a running total of the pair work (reversible)
    he = 0.e0
    
    !calculate artificial viscosity
    
!    call calc_art_viscosity(nl%drij(k,:),nl%dv(k,:),(pa%rho(i)+pb%rho(j))/2 &
!    &,(pa%sml(i)+pb%sml(j))/2,(pa%c(i)+pb%c(j))/2,avisc)   
   
!    call core_force(nl%drij(k,:),nl%rijsq(k),bnd%bsig,repulsion,bnd%bcore)
    
!    p%a(i,1) = p%a(i,1) + repulsion(1)/p%m(i)
!    p%a(i,2) = p%a(i,2) + repulsion(2)/p%m(i)
    ! x coodrdinate of acceleration
    ! reversible contribution
!    hx = - (pa%p_rev(i,1,1)/pa%rho(i)**2 + pb%p_rev(j,1,1)/pb%rho(j)**2 )*nl%dwdx(k,1)&
!    &    - (pa%p_rev(i,1,2)/pa%rho(i)**2 + pb%p_rev(j,1,2)/pb%rho(j)**2 )*nl%dwdx(k,2)
          
    ! y coordinate of acceleration
    ! reversible contribution
!    hy = -( pa%p_rev(i,1,2)/pa%rho(i)**2 + pb%p_rev(j,1,2)/pb%rho(j)**2 )*nl%dwdx(k,1)&
!    &    -( pa%p_rev(i,2,2)/pa%rho(i)**2 + pb%p_rev(j,2,2)/pb%rho(j)**2 )*nl%dwdx(k,2)

    ! long range contribution
!    hx = hx -( pa%pco(i)/pa%rho(i)**2 + pb%pco(j)/pb%rho(j)**2 ) &
!       * nl%dwdx_long(k,1)
!    hy = hy -( pa%pco(i)/pa%rho(i)**2 + pb%pco(j)/pb%rho(j)**2 ) &
!       * nl%dwdx_long(k,2)

    ! artificial viscosity
!    hx = hx-avisc*nl%dwdx(k,1) + repulsion(1)
!    hy = hy-avisc*nl%dwdx(k,2) + repulsion(2)
    
    ! reversible work
!    he= he+ nl%dv(k,1)*hx
!    he= he+ nl%dv(k,2)*hy
          
    if(debug) print *,'force contribution of',hx,hy
   
    ! x coodrdinate of acceleration 
    ! irreversible contribution
!    hx = hx+(pa%pi_irr(i,1,1)/pa%rho(i)**2 + pb%pi_irr(j,1,1)/pb%rho(j)**2 ) &
!    &      * nl%dwdx(k,1) &
!    &      + (pa%pi_irr(i,1,2)/pa%rho(i)**2 + pb%pi_irr(j,1,2)/pb%rho(j)**2 ) &
!    &      * nl%dwdx(k,2)
          
    ! y coordinate of acceleration
    ! irreversible contribution
!    hy = hy + ( pa%pi_irr(i,1,2)/pa%rho(i)**2 + pb%pi_irr(j,1,2)/pb%rho(j)**2 ) &
!    &    *nl%dwdx(k,1)&
!    &       + ( pa%pi_irr(i,2,2)/pa%rho(i)**2 + pb%pi_irr(j,2,2)/pb%rho(j)**2 ) &
!    &      *nl%dwdx(k,2)
   
!    pa%a(i,1) = pa%a(i,1) + (pb%m(j)*hx + repulsion(1)/pa%m(i))
!    pb%a(j,1) = pb%a(j,1) - (pa%m(i)*hx + repulsion(1)/pb%m(j))
!    pa%a(i,2) = pa%a(i,2) + (pb%m(j)*hy + repulsion(2)/pa%m(i))
!    pb%a(j,2) = pb%a(j,2) - (pa%m(i)*hy + repulsion(2)/pb%m(j))
     
    !if(validate) call check_zero_distance(dx(k,:),drsq(k),rho,nl)
    !heat conduction
    if(debug) print *, 'calculating heat conduction for ', i, j
    !the heat conduction should be counted as part of tdsdt
    !if we ever want to calculate the entropy

    if(conduction_type .eq. 1) then
        !monaghan type conduction
        
      if(dimn .eq. 2) then  
        f = ( nl%dwdx(k,2)*nl%drij(k,2)+nl%dwdx(k,1)*nl%drij(k,1) ) / nl%rijsq(k)
      else if(dimn .eq. 1) then
        f= nl%dwdx(k,1)*nl%drij(k,1)/nl%rijsq(k)
      endif
 
      pa%dedt(i) = pa%dedt(i) - pb%m(j) * f * ( (thermalk + thermalk)* (pa%temp(i) - pb%temp(j)))/(pa%rho(i)*pb%rho(j))
      pb%dedt(j) = pb%dedt(j) + pb%m(j) * f * ( (thermalk + thermalk)* (pa%temp(i) - pb%temp(j)))/(pa%rho(i)*pb%rho(j))

 
    else if(conduction_type .eq. 2) then
      !print *,pb%m(j)*( pa%q(i,1)/(pa%rho(i)**2) + pb%q(j,1)/(pb%rho(j)**2) )*nl%dwdx(k,1)
      !Using full heat flux vector - ref Nugent and Posch
      !calculate contribution to dedt
      pa%dedt(i) = pa%dedt(i) - pb%m(j)*( pa%q(i,1)/(pa%rho(i)**2) + pb%q(j,1)/(pb%rho(j)**2) )*nl%dwdx(k,1)
      pa%dedt(i) = pa%dedt(i) - pb%m(j)*( pa%q(i,2)/(pa%rho(i)**2) + pb%q(j,2)/(pb%rho(j)**2) )*nl%dwdx(k,2)
      
      pb%dedt(j) = pb%dedt(j) + pb%m(j)*( pa%q(i,1)/(pa%rho(i)**2) + pb%q(j,1)/(pb%rho(j)**2) )*nl%dwdx(k,1)
      pb%dedt(j) = pb%dedt(j) + pb%m(j)*( pa%q(i,2)/(pa%rho(i)**2) + pb%q(j,2)/(pb%rho(j)**2) )*nl%dwdx(k,2)
    endif
    
    !needs to be symmetrised consistenty
    !pa%dedt(i) = pa%dedt(i) +  pb%m(j)*he/2
   ! pb%dedt(j) = pb%dedt(j) +  pa%m(i)*he/2

    if(validate) then
      call validate_double(pa%dedt(i),'dedt')
      call validate_double(pb%dedt(j),'dedt')
      call validate_double(pa%a(i,1),'a')
      call validate_double(pa%a(i,2),'a')
      call validate_double(pb%a(j,1),'a')
      call validate_double(pb%a(j,2),'a')
    endif

  enddo !end loop over pairs

  !failed experiment in simplifying dedt
  !! \todo investigate this further (not zeroing grad v was most likely 
  !! the issue. For realtime graphics this approach should be much cheaper i
  !! than the whole -tensor approach.
  !! PROTOTYPE THIS IN PYTHON FIRST :-)
  !! do i=1, n
  !! dedt_rev1(i) = -p_rev(i)*div_v(i)/rho(i)
  !! dedt(i) = dedt(i) + tdsdt(i)
  !! enddo
 
  ! change in specific internal energy de/dt = T ds/dt - p/rho
  do i=1, pb%n
!    pa%dedt(i) = tdsdt(i) + pa%dedt(i)
!    pb%dedt(i) = tdsdt(i) + pb%dedt(i)
    dedt_total = dedt_total + pb%dedt(i)
  enddo
 
  if(verbose) print *, 'finished calculating sphforces '
     
  if(profile) then   
    call system_clock(count = end_time)
    time_seconds = (end_time - start_time)/real(time_rate)
    print *,'**** force calculation time ',i, time_seconds
  endif

end subroutine calc_sphforce_ab

subroutine v_smooth_ab(p,nl)
  ! v_smooth is a function that smoothed rdot
  ! adds contributions from seperate systems
  ! \todo UNFINISHED
  type (particle_type), intent(inout) :: p
  type (nlist_data_type) :: nl
  integer i,j,k

  do k = 1, nl%nip
   i = nl%ilist(k,1)
   j = nl%ilist(k,2)

   if( (p%rho(i)+p%rho(j)) .eq. 0.0 ) stop 'divide by zero error in dx'
             
   p%rdot(i,:) = p%rdot(i,:) + v_eps * (p%v(j,:) - p%v(i,:)) * ((p%m(i)+p%m(j))) *nl%w(k)  / (p%rho(i)+p%rho(j)) 
   p%rdot(j,:) = p%rdot(j,:) + v_eps * (p%v(i,:) - p%v(j,:)) * ((p%m(i)+p%m(j))) *nl%w(k)  / (p%rho(i)+p%rho(j)) 

  end do


end subroutine v_smooth_ab

subroutine calc_grad_v_ab(nl,pa,pb)
  implicit none
  type (nlist_data_type) :: nl
  type (particle_type), intent(inout) :: pa
  type (particle_type), intent(inout) :: pb
  
  integer :: i,j,a,b,k
  
  !dv and dwdx have dimensions of (npairs,dimensions)
  
  !grad_v = 0;
  
  do k = 1, nl%nip
    i = nl%ilist(k,1)
    j = nl%ilist(k,2)
  
    !loop over dimensions
    do a=1,dimn
      do b=1,dimn     
        pa%grad_v(i,a,b) =  pa%grad_v(i,a,b)+(pb%m(j)/pb%rho(j))*nl%dv(k,a)*nl%dwdx(k,b)
        pb%grad_v(j,a,b) =  pb%grad_v(j,a,b)-(pa%m(i)/pa%rho(i))*nl%dv(k,a)*nl%dwdx(k,b)
      end do
    end do
  end do

end subroutine calc_grad_v_ab


subroutine calc_div_v_ab(grad_v, div_v)
  double precision, dimension(:,:,:) :: grad_v
  double precision, dimension(:) :: div_v
  integer :: i,n
   
  n = size(grad_v,1)
  do i = 1,n        
    div_v(i) = grad_v(i,1,1) + grad_v(i,2,2)        
  end do
  
end subroutine calc_div_v_ab


subroutine calc_grad_v_os_ab(grad_v_os,grad_v,div_v)
  double precision, dimension(:,:,:) :: grad_v
  double precision, dimension(:) :: div_v
  double precision, dimension(:,:,:) :: grad_v_os
  integer :: n,i,d
  n = size(grad_v,1)
  d = 2.0
  
  do i=1,n
    grad_v_os(i,1,1) = (1/2) * (grad_v(i,1,1)*2.0) - (1.0/d)*(div_v(i))   
    grad_v_os(i,1,2) = (1/2) * ( grad_v(i,1,2) + grad_v(i,2,1) )
    grad_v_os(i,2,1) = (1/2) * ( grad_v(i,2,1) + grad_v(i,1,2) )
    grad_v_os(i,2,2) = (1/2) * (grad_v(i,2,2)*2.0) - (1.0/d)*(div_v(i))
  end do

end subroutine calc_grad_v_os_ab


subroutine calc_iso_pressure_ab(p)
  implicit none
  type (particle_type), intent(inout) :: p
  integer ::  i
  if (verbose) print *,'calculating iso pressure'
 
  do i=1,p%n
    if(eos_type .eq. 1) then
        call igeos(p%rho(i), p%u(i), p%p(i), p%c(i))
      else if(eos_type .eq. 2) then
        call vdweos(p%rho(i), p%u(i), p%p(i))
      else if(eos_type .eq. 3) then
        call vdweos_attractive(p%rho(i),p%pco(i))
        call vdweos_repulsive(p%rho(i),p%u(i),p%p(i))
      !else if(eos_type .eq. 4) then
      !  call vdw_gradient(p%rho(i),p%u(i),p%p(i),p%grad_rho(i,:),p%cgrad(i))
      endif  
  end do
end subroutine calc_iso_pressure_ab


subroutine calc_heat_flux_ab(nl,pa,pb)
  type(nlist_data_type) :: nl
  type(particle_type) :: pa
  type(particle_type) :: pb
  integer i,j,k
  double precision :: m,rho,dtmp

!  pa%q(:,:)=0
  pb%q(:,:)=0
      
  do k = 1, nl%nip
    i = nl%ilist(k,1)
    j = nl%ilist(k,2)
    m = (pa%m(i)+pb%m(j))/2.
    rho = (pa%rho(i)+pb%rho(j))/2.
    dtmp = pb%temp(j) - pa%temp(i) 
    
    pa%q(i,1) = pa%q(i,1) - thermalk * ( m/rho) * dtmp * nl%dwdx(k,1)
    pa%q(i,2) = pa%q(i,2) - thermalk * ( m/rho) * dtmp * nl%dwdx(k,2)
   
    pb%q(j,1) = pb%q(j,1) - thermalk * ( m/rho ) * dtmp * nl%dwdx(k,1)
    pb%q(j,2) = pb%q(j,2) - thermalk * ( m/rho ) * dtmp * nl%dwdx(k,2) 
  end do

end subroutine calc_heat_flux_ab


subroutine calc_pi_os_ab(pi_os,grad_v_os,eta)
  implicit none
  double precision, dimension(:,:,:) :: pi_os, grad_v_os
  double precision, dimension(:) :: eta
  integer i,j,k
  
  do i=1,size(pi_os,1)
    do j=1,dimn
    do k=1,dimn
      pi_os(i,j,k) = -2*eta(i) * grad_v_os(i,j,k)
    end do
    end do
  end do

end subroutine calc_pi_os_ab

subroutine calc_pi_one_ab(pi_one,div_v,zeta)
  implicit none
  double precision, dimension(:) :: pi_one,div_v,zeta
  integer i
  do i=1,size(pi_one)
    pi_one(i) = - zeta(i) * div_v(i)
  end do
end subroutine calc_pi_one_ab


subroutine calc_viscous_entropy_os_ab(nl,p_ir,tdsdt,grad_v_os,rho)
  !calculate the symmetric traceless 
  !contribution to dedt
  !pass in pi_os
  
  type(nlist_data_type) :: nl
  double precision, dimension(:,:,:) :: grad_v_os
!  double precision, dimension(:,:) :: dv,dwdx
  double precision, dimension(:) :: tdsdt,rho
  double precision, dimension(:,:,:) :: p_ir
  integer :: i,j,k

if(entropy_type .eq. 1) then
    do i=1,size(rho)
      if(dimn .eq. 1) then
        tdsdt(i) = p_ir(i,1,1)* p_ir(i,1,1)
      else if( dimn .eq. 2) then
        tdsdt(i) = p_ir(i,1,1) * p_ir(i,1,1) + 2.e0 * p_ir(i,1,2)*p_ir(i,1,2) + p_ir(i,2,2)*p_ir(i,2,2)
      endif
   end do
 
  else if (entropy_type .eq. 2) then
   
  do k = 1, nl%nip
   i = nl%ilist(k,1)
   j = nl%ilist(k,2) 
      ! is there a factor of rho_i missing?
      tdsdt(i) = tdsdt(i) + ( p_ir(i,1,1)/rho(i)**2         &
      &        + p_ir(j,1,1)/rho(j)**2 )*grad_v_os(i,1,1)
      tdsdt(i) = tdsdt(i) + ( p_ir(i,1,2)/rho(i)**2 + p_ir(j,1,2)/rho(j)**2 )*grad_v_os(i,2,1)
      tdsdt(i) = tdsdt(i) + ( p_ir(i,2,1)/rho(i)**2 + p_ir(j,2,1)/rho(j)**2 )*grad_v_os(i,1,2)
      tdsdt(i) = tdsdt(i) + ( p_ir(i,2,2)/rho(i)**2 + p_ir(j,2,2)/rho(j)**2 )*grad_v_os(i,2,2)

      tdsdt(j) = tdsdt(j) + ( p_ir(j,1,1)/rho(j)**2 + p_ir(i,1,1)/rho(i)**2 )*grad_v_os(j,1,1)
      tdsdt(j) = tdsdt(j) + ( p_ir(j,1,2)/rho(j)**2 + p_ir(i,1,2)/rho(i)**2 )*grad_v_os(j,2,1)
      tdsdt(j) = tdsdt(j) + ( p_ir(j,1,2)/rho(j)**2 + p_ir(i,1,2)/rho(i)**2 )*grad_v_os(j,1,2)
      tdsdt(j) = tdsdt(j) + ( p_ir(j,2,2)/rho(j)**2 + p_ir(i,2,2)/rho(i)**2 )*grad_v_os(j,2,2)
    end do
       
  end if

end subroutine calc_viscous_entropy_os_ab

subroutine calc_gradgradrho_ab(nl,pa,pb,gradgradrho_a,gradgradrho_b) 

  type(nlist_data_type) :: nl
  type (particle_type), intent(inout) :: pa
  type (particle_type), intent(inout) :: pb
  
  integer i,j,k
  double precision, dimension(pa%n) :: gradgradrho_a
  double precision, dimension(pb%n) :: gradgradrho_b
  
  gradgradrho_a=0
  gradgradrho_b=0
  
  ! calculate grad grad rho   
  do k = 1, nl%nip
   i = nl%ilist(k,1)
   j = nl%ilist(k,2) 
   
   !check this section
   !is there a 1/rho missing?
   
      gradgradrho_a(i) = gradgradrho_a(i) - pa%rho(i)*pb%m(j)*( pa%grad_rho(i,1)/(pa%rho(i)**2) &
      & + pb%grad_rho(j,1)/(pb%rho(j)**2) )*nl%dwdx(k,1)
      gradgradrho_a(i) = gradgradrho_a(i) - pa%rho(i)*pb%m(j)*( pa%grad_rho(i,2)/(pa%rho(i)**2) &
      & + pb%grad_rho(j,2)/(pb%rho(j)**2) )*nl%dwdx(k,2)
      gradgradrho_b(j) = gradgradrho_b(j) + pb%rho(j)*pa%m(i)*( pa%grad_rho(i,1)/(pa%rho(i)**2) &
      & + pb%grad_rho(j,1)/(pb%rho(j)**2) )*nl%dwdx(k,1)
      gradgradrho_b(j) = gradgradrho_b(j) + pb%rho(j)*pb%m(j)*( pa%grad_rho(i,2)/(pa%rho(i)**2) &
      & + pb%grad_rho(j,2)/(pb%rho(j)**2) )*nl%dwdx(k,2)
  end do

end subroutine calc_gradgradrho_ab
  
subroutine calc_capillary_pressure_ab(p,gradgradrho) 

  type (particle_type), intent(inout) :: p
  double precision, dimension(p%n) :: gradgradrho
  double precision, dimension(p%n) :: maggradrhosq
  integer i
  
  do i=1,p%n
    maggradrhosq = abs(p%grad_rho(i,1)*p%grad_rho(i,1) + p%grad_rho(i,2)*p%grad_rho(i,2))
    p%p_rev(i,1,1) = p%p_rev(i,1,1) - p%cgrad(i)*p%rho(i)*gradgradrho(i) - p%cgrad(i)*(0.5)*maggradrhosq(i)
    p%p_rev(i,2,2) = p%p_rev(i,2,2) - p%cgrad(i)*p%rho(i)*gradgradrho(i) - p%cgrad(i)*(0.5)*maggradrhosq(i)
  
    p%p_rev(i,1,1) = p%p_rev(i,1,1) + p%cgrad(i)*p%grad_rho(i,1)**2
    p%p_rev(i,1,2) = p%p_rev(i,1,2) + p%cgrad(i)*p%grad_rho(i,1)*p%grad_rho(i,2)
    p%p_rev(i,2,1) = p%p_rev(i,2,1) + p%cgrad(i)*p%grad_rho(i,1)*p%grad_rho(i,2)
    p%p_rev(i,2,2) = p%p_rev(i,2,2) + p%cgrad(i)*p%grad_rho(i,2)**2
  enddo
  
end subroutine calc_capillary_pressure_ab


end module sphforce_ab
