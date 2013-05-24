! Applies boundaries to particles.
! handles pbcs
! Handles some aspects of boundary particles
! Handles gravity
! Not a good candidate for python wrapping because of the
! heavy use of data structures.

module boundary
use global
use eos
use neighbour_list
use kernel
use simulation_box
use particle
use core_potential
implicit none

private
public :: apply_boundary, initialise_boundary_particles, apply_gravity, set_upper_bounds, set_lower_bounds, calc_sphforce_boundary


contains

subroutine initialise_boundary_particles(bp)
  ! Initialise the boundary particle structure with some
  ! default values
  implicit none
  type(particle_type),intent(inout) :: bp

  bp%x=0
  bp%v=0
  bp%a=0
  bp%u=0
  bp%temp= 1.0   !temperature
  bp%m=1.1    !mass
  bp%rho = 1.2  !density
  bp%p  =0   !pressure
  bp%c = 0    !sound speed
  bp%eta  = 0 !viscosity
  bp%sml  = 3.0
  bp%pco = 0
  bp%sml_long = 6.0
  bp%dedt = 0
  bp%xhalf = 0
  bp%v_flow=0
  bp%t_smooth = 0
  bp%rdot =0
  bp%p_smooth = 0
  bp%q = 0
  bp%grad_rho=0.0
  bp%zeta=0
  bp%cgrad=0.0
  bp%pi_one=0.0
  
end subroutine initialise_boundary_particles

subroutine set_lower_bounds(x)
  ! Set up the positions of boundary particles such
  ! that they are in two staggered layers
  double precision, dimension(:,:) :: x
  double precision :: dx
  integer i,n,m
  n = size(x,1)
  m = n/2.
  dx = lx/m

  ! two layers
  do i=1,m
    x(i,2) = 0.0 
    x(i,1) = i*dx
  end do

  do i=1,m
    x(m+i,2)=dx/2.
    x(m+i,1)=dx/2.+(i-1)*dx
  end do

end subroutine set_lower_bounds

subroutine set_upper_bounds(x)
  ! Set up the positions of boundary particles such
  ! that they are in two staggered layers
  double precision, dimension(:,:) :: x
  double precision  :: dx
  integer i,n,m
  n = size(x,1)
  m = n/2.
  dx = lx/m
  ! first layer
  do i=1,m
    x(i,2) =ly 
    x(i,1) = i*dx
  end do
  do i=1,m
    x(m+i,2)=ly-dx/2.
    x(m+i,1)=dx/2.+(i-1)*dx
  end do
end subroutine set_upper_bounds

subroutine apply_gravity(a)
  ! Applies constant external acceleration to particles
  double precision, dimension(:,:) :: a
  a(:,2) = a(:,2) + grav
end subroutine apply_gravity  

subroutine apply_boundary(r,v,box)
  ! Reflective boundaries - which sides of the box are reflective is
  ! determined by the input configuration
  type (box_data_type), intent(in) :: box
  double precision, dimension(:,:) :: r
  double precision, dimension(:,:) :: v
  
  integer n,i
  
  n = size(r,1)
 
  do i=1,n
    if( bounds_reflect(1) ) then
      if( (r(i,1) .lt. 0)) then 
        r(i,1) = 0.0
        if(v(i,1) .lt. 0) then
          v(i,1) = -v(i,1)
        endif
      endif
    endif
    if( bounds_reflect(2) ) then
      if( r(i,1) .gt. box%boxvec(1,1)) then
        r(i,1) = box%boxvec(1,1)
        if (v(i,1) .gt. 0 ) then
        v(i,1) = -v(i,1)
        endif
      endif
    endif
    if( bounds_reflect(3) ) then
      if( r(i,2) .lt. 0.0) then
        r(i,2) = 0.0
     if (v(i,2) .lt. 0 ) then
        v(i,2) = -v(i,2)
        endif
      endif
    endif
    if( bounds_reflect(4) ) then
      if( r(i,2) .gt. box%boxvec(2,2) )  then
        r(i,2) = box%boxvec(2,2) 
      if (v(i,2) .gt. 0 ) then
        v(i,2) = -v(i,2)
      endif
      endif
    endif
  end do

end subroutine apply_boundary


subroutine calc_sphforce_boundary(bnl,p)
  ! Computes the sph force between normal particles
  ! and boundary particles.
  type (nlist_data_type), intent(inout) :: bnl
  type (particle_type), intent(inout) :: p
  double precision, dimension(2) :: repulsion
  integer :: k,i,j
  ! just core repulsion for now
  do k = 1,bnl%nip 
    i = bnl%ilist(k,1)
    j = bnl%ilist(k,2)
    if(debug) print *, 'boundary core loop for pair',k
    call core_force(repulsion,bnl%drij(k,:),bnl%rijsq(k),bsig,bcore,dimn)
    p%a(i,1) = p%a(i,1) + repulsion(1)/p%m(i)
    p%a(i,2) = p%a(i,2) + repulsion(2)/p%m(i)
  end do
end subroutine calc_sphforce_boundary




end module boundary
