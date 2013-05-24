!> Defines the smoothed particle data type. All properties including mass,
!! density and velocity are stored in this data structure.
!! the number of particles and dimension are set parameters
module particle

use global
use eos
use neighbour_list
use reader
use simulation_box
use kernel
use collision
use art_viscosity

implicit none

public :: particle_type
public :: read_particle_state, validate_particles
public :: calc_kernels_all, calc_smoothed_properties, calc_kernels_ab
public :: collide_particles, calc_art_viscosity
public :: update_temperature, update_energy, update_pressure

type particle_type
!  n  !<  number of particles
!  dimn !< dimension of system
!  x !< positions
!  v !< velocities
!  a !< accelerations
!  !> used for leapfrog integrator - positions 1/2 a timestep ago
!  xhalf 
!  q !<heat flux
!  u !< internal energy
!  m !< mass
!  rho !< density
!  rho_lr !< long range density
!  p !< pressure
!  c !< speed of sound
!  eta !< shear viscosity
!  temp !< temperature
!  co !< long range cohesive pressure
!  pco !> internal energy time deriv
!  double precision, dimension(:),allocatable :: dedt 
!  double precision, dimension(:),allocatable :: zeta !< bulk viscosity
!  !> gradient energy coefficient
!  double precision, dimension(:),allocatable :: cgrad 
!  double precision, dimension(:,:),allocatable :: grad_rho !< density gradient
!  double precision, dimension(:,:,:), allocatable :: grad_v !< velocity gradient
!  !> smoothed velocity at the particle location
!  double precision, dimension(:,:),allocatable :: v_flow
!  !> smoothed pressure at the particle location
!  double precision, dimension(:), allocatable :: p_smooth
!  !> smoothed temperature at the particle location
!  double precision, dimension(:), allocatable :: t_smooth
!
!!  !> 
  !! 
!  rdot --  smoothed velocity for moving particles 
!         ( Monaghan velocity smoothing)
!         not the same as smoothed field velocity
! cohesive pressure is used for vanderwaals equation of state
! with different smoothing lengths for the attractive and repulsive
! components
! sml -- smoothing length
! sml_long -- long smoothing length
! pressure tensor
! pi_os -- symmetric traceless viscous (irreversible) pressure tensor
! pi_irr -- pi_irr is the total irreversible pressure tensor
! pi_one -- pi_one is the isotropic viscous (irreversible) pressure tensor
! p_rev -- p_rev is the reversible, nonisotropic pressure tensor
! p_rev_lr -- long range reversible pressure accumulates core potential 
! p_eq -- the isotropic part of the equilibrium pressure
! thermostat  --

  integer:: n  
  integer:: dimn 
  double precision, dimension(:,:), allocatable :: x 
  double precision, dimension(:,:), allocatable :: v
  double precision, dimension(:,:), allocatable :: a
  double precision, dimension(:,:), allocatable :: xhalf 
  double precision, dimension(:,:), allocatable :: q
  double precision, dimension(:), allocatable :: u
  double precision, dimension(:), allocatable :: m
  double precision, dimension(:), allocatable :: rho
  double precision, dimension(:), allocatable :: rho_lr 
  double precision, dimension(:), allocatable :: p
  double precision, dimension(:), allocatable :: c
  double precision, dimension(:), allocatable :: eta
  double precision, dimension(:), allocatable :: temp
  double precision, dimension(:), allocatable :: pco
  double precision, dimension(:), allocatable :: dedt 
  double precision, dimension(:), allocatable :: zeta
  double precision, dimension(:), allocatable :: cgrad 
  double precision, dimension(:,:), allocatable :: grad_rho
  double precision, dimension(:,:,:), allocatable :: grad_v
  double precision, dimension(:,:), allocatable :: v_flow
  double precision, dimension(:), allocatable :: p_smooth
  double precision, dimension(:), allocatable :: t_smooth
  double precision, dimension(:,:), allocatable :: rdot
  double precision, dimension(:), allocatable :: sml
  double precision, dimension(:), allocatable :: sml_long
  double precision, dimension(:,:,:), allocatable :: pi_os
  double precision, dimension(:,:,:), allocatable :: pi_irr
  double precision, dimension(:), allocatable :: pi_one
  double precision, dimension(:,:,:), allocatable :: p_rev
  double precision, dimension(:,:,:), allocatable :: p_rev_lr
  double precision, dimension(:), allocatable :: p_eq
  double precision :: thermostat

end type particle_type

!! \todo
!! type smoothed_point_properties_type
!! sorry about the long name.
!! typically you would use this struct to hold the smoothed
!! (flow) values for properties at the particle positions
!! (as distinct from the values attached to the particles)
!! this could also be used to track a set of points
!! not associated with particles
!! end type smoothed_point_properties_type

contains

subroutine initialise_particles(p)
  implicit none

  type(particle_type),intent(inout) :: p
    p%x = 0
    p%v = 0
    p%a = 0
    p%u = 0        ! internal energy
    p%temp = 0.0   ! temperature
    p%m = 1        ! mass
    p%rho = 1      ! density
    p%rho_lr = 1   ! density
    p%p = 0        ! pressure
    p%c = 0        ! sound speed
    p%eta = 0      ! viscosity
    p%sml = 0
    p%pco = 0
    p%sml_long = 0
    p%dedt = 0
    p%xhalf = 0
    p%v_flow = 0
    p%t_smooth = 0
    p%rdot = 0
    p%p_smooth = 0
    p%q = 0
    p%grad_rho = 0
    p%zeta = 0
    p%cgrad = cgrad_ini
    p%thermostat = 10.0
    p%pi_one = 0.0
    p%p_eq = 0.0

end subroutine initialise_particles

subroutine create_particles(p,n,d)
  ! Allocates memory to particle structure
  implicit none
  type (particle_type), intent(inout) :: p
  integer :: error, n, d
  p%n = n
  p%dimn = d

  ! Allocate arrays for particles
  allocate(p%x(n,d), p%v(n,d),p%a(n,d),p%q(n,d),p%u(n), stat = error)
  allocate(p%rho_lr(n), stat = error)
  allocate(p%m(n),p%rho(n),p%p(n),p%pco(n),p%c(n),p%eta(n),stat = error)
  allocate(p%zeta(n),stat = error)
  allocate(p%sml(n),p%temp(n),p%sml_long(n),p%dedt(n),p%cgrad(n), stat=error)
  allocate(p%grad_rho(n,d),p%xhalf(n,d),p%rdot(n,d), stat=error)  
  ! Allocate arrays for smoothed properties at particle locations
  allocate(p%v_flow(n,d),p%p_smooth(n),p%t_smooth(n), stat=error)
  allocate(p%grad_v(n,d,d), stat=error)

  ! allocate arrays for the pressure tensor 
  allocate(p%pi_os(n,d,d), p%pi_irr(n,d,d),stat=error)
  allocate(p%pi_one(n),stat=error)
  allocate(p%p_rev(n,d,d),stat=error)
  allocate(p%p_rev_lr(n,d,d),stat=error)
  allocate(p%p_eq(n),stat=error)
  
  if(error .ne. 0) stop 'allocation error in particles'
   
end subroutine create_particles
   
!> Checks the entire particle data structure for inf and NaN
!! the validate_double subroutine will crash the program if
!! it finds anything amiss
subroutine validate_particles(p)
  type (particle_type) :: p
  integer :: i,d,e

  if(verbose) then
    print *,'Validating particles for',p%n,'particles'
  endif

  do i=1,p%n
    do d=1,p%dimn
      do e=1,p%dimn
        ! second rank tensors
        call validate_double(p%pi_os(i,d,e),'pi_os')
        call validate_double(p%pi_irr(i,d,e),'pi_irr')
        call validate_double(p%p_rev(i,d,e),'p_rev')
        call validate_double(p%p_rev_lr(i,d,e),'p_rev')
      end do
      ! first rank tensors (vectors)
      call validate_double(p%x(i,d),'x')
      call validate_double(p%v(i,d),'v')
      call validate_double(p%a(i,d),'a')
      call validate_double(p%rdot(i,d),'rdot')
      call validate_double(p%v_flow(i,d),'vflow')
      call validate_double(p%grad_rho(i,d),'gradrho')
      call validate_double(p%xhalf(i,d),'xhalf')
      call validate_double(p%q(i,d),'q')
    end do
    ! zeroth rank tensors (scalars)
    call validate_double(p%temp(i),'t')
    call validate_double(p%m(i),'m')
    call validate_double(p%rho(i),'rho')
    call validate_double(p%p(i),'p')
    call validate_double(p%c(i),'c')
    call validate_double(p%eta(i),'eta')
    call validate_double(p%sml(i),'sml')
    call validate_double(p%pco(i),'pco')
    call validate_double(p%sml_long(i),'smllong')
    call validate_double(p%dedt(i),'dedt')
    call validate_double(p%t_smooth(i),'tsm')
    call validate_double(p%p_smooth(i),'psm')
    call validate_double(p%zeta(i),'zeta')
    call validate_double(p%cgrad(i),'cgrad')
    call validate_double(p%pi_one(i),'pione')
    call validate_double(p%p_eq(i),'pieq')
  end do
end subroutine validate_particles


!> deallocates arrays
subroutine destroy_particles(p)
   type (particle_type) :: p
    if(allocated(p%x)) deallocate(p%x)
    if(allocated(p%v)) deallocate(p%v)
    if(allocated(p%a)) deallocate(p%a)
    if(allocated(p%u)) deallocate(p%u)
    if(allocated(p%m)) deallocate(p%m)
    if(allocated(p%rho)) deallocate(p%rho)
    if(allocated(p%p)) deallocate(p%p)
    if(allocated(p%pco)) deallocate(p%pco)
    if(allocated(p%c)) deallocate(p%c)
    if(allocated(p%zeta)) deallocate(p%zeta)
    if(allocated(p%eta)) deallocate(p%eta)
    if(allocated(p%sml)) deallocate(p%sml)
    if(allocated(p%temp)) deallocate(p%temp)
    if(allocated(p%xhalf)) deallocate(p%xhalf)
    if(allocated(p%v_flow)) deallocate(p%v_flow)
    if(allocated(p%p_smooth)) deallocate(p%p_smooth)
    if(allocated(p%t_smooth)) deallocate(p%t_smooth)
    if(allocated(p%q)) deallocate(p%q)
    if(allocated(p%grad_rho)) deallocate(p%grad_rho)
    if(allocated(p%grad_v)) deallocate(p%grad_v) 
    if(allocated(p%cgrad)) deallocate(p%cgrad)
    if(allocated(p%rdot)) deallocate(p%rdot)
    if(allocated(p%pi_one)) deallocate(p%pi_one)
    if(allocated(p%p_rev)) deallocate(p%p_rev)
    if(allocated(p%p_rev_lr)) deallocate(p%p_rev_lr)
end subroutine destroy_particles

subroutine rd_particle_input(p,ifile,iunit)
  implicit none  
  type(particle_type),intent(inout) :: p
  integer :: iunit
  double precision :: etain,zetain,cin
  character(len=32) :: ifile         ! input file name
  open(iunit,file=ifile,status='old',form='formatted', &
     & position='rewind')
  call read_dbl(iunit,'ETA',etain)
  call read_dbl(iunit,'ZETA',zetain)
  call read_dbl(iunit,'SOUND_SPEED',cin)
  p%eta(:) = etain
  p%zeta(:) = zetain
  p%c(:) = cin

 close(iunit)
end subroutine rd_particle_input

!> Reads particle states. Seriously.
!! Does not use the boundary state.
subroutine read_particle_state(p)
  !! \todo Reads the wrong pressure?
  
  type (particle_type) :: p
  integer :: funit, i
  funit = 1;
  print *,p%n

  open(unit=funit,file='sphstate.input', &
      & form='formatted',status='old',position='rewind')
  do i=1,p%n
    print *,'read particle',i
    read(funit,*) p%x(i,1),p%x(i,2),p%v(i,1),p%v(i,2),&
         &        p%a(i,1),p%a(i,2),p%temp(i),p%u(i), &
         &        p%m(i),p%rho(i),p%p(i),p%c(i),p%eta(i),p%sml(i)
  end do
  close(funit)
  
end subroutine read_particle_state

!> If any particle is moving so fast that it can cover a distance of 0.5h,
!! the simulation is stopped. This may be overly restrictive, especially for 
!! development work.
subroutine check_velocity(p,d,maxn,dt)
  type (particle_type) :: p
  integer :: d,maxn
  double precision :: dt
  integer :: i
  if(d==2) then
  if(maximum_velocity == 1) then
  
    do i=1,maxn
      if( ((p%v(i,2)*dt)**2 + (p%v(i,1)*dt)**2) > ((p%sml(i)/2)**2) ) then
        stop 'moving too fast'
      endif
    end do
    
  endif
  endif 
end subroutine check_velocity


subroutine calc_kernels_all(p,nl)
  ! Calculates kernels for all particles
  implicit none
  type (particle_type), intent(inout) :: p
  type (nlist_data_type) :: nl
  integer :: n, k, i, j
  double precision :: rhoav,selfdens,r0
  double precision, dimension(dimn) :: dx0,dw0

  n = p%n

  if(verbose) print *,"Calculating all kernels"

  if( adapt_smoothing == 2 ) then
    ! Use ADKE (Sigalotti et al 2005)
    ! Compute a kernel based on the initial smoothing length 

    do k = 1, nl%nip
    call calc_kernel(nl%rij(k), & 
                  & nl%drij(k,:), &
                  & smoothing_ini, &
                  & nl%w(k),nl%dwdx(k,:), &
                    kernel_type,dimn)

    call calc_kernel(nl%rij(k), &
                  & nl%drij(k,:), &
                  & co_smoothing_ini, &
                  & nl%w_long(k),nl%dwdx_long(k,:), &
                    kernel_type,dimn)
    end do

    ! use this kernel to compute a pilot density

    r0=0.0
    dx0=0.0
    dw0=0.0    
    do i=1,n
      call calc_kernel(r0,dx0,smoothing_ini,selfdens,dw0,kernel_type,dimn)
      p%rho(i) = selfdens*p%m(i)
    end do

    do k = 1, nl%nip
      i = nl%ilist(k,1)
      j = nl%ilist(k,2)
      p%rho(i) = p%rho(i) + p%m(j)*nl%w(k)
      p%rho(j) = p%rho(j) + p%m(i)*nl%w(k)
    end do
 
    do i=1,n
      rhoav = rhoav+p%rho(i)
    end do
    rhoav = rhoav/n

    ! rescale the smoothing lengths based on each particle's
    ! density ratio
    do i=1,p%n
      p%sml(i) = 1.0 * smoothing_ini * 1/sqrt(p%rho(i)/rhoav)  
      p%sml_long(i) = 1.0 * co_smoothing_ini  * 1/sqrt(p%rho(i)/rhoav)
    end do   

  endif
    
  ! Now just do a normal kernel calculation

  do k = 1, nl%nip
    i = nl%ilist(k,1)
    j = nl%ilist(k,2)

    if(debug) then
      print *,nl%nip
      print *,nl%nab
      print *, 'calculating kernel for (pair,1,2)',k, i, j
      print *, 'distance', nl%rij(k)
      print *, 'disp', nl%drij(k,:)
    end if

    call calc_kernel(nl%rij(k), nl%drij(k,:), &
    & (p%sml(i)+p%sml(j))/2,nl%w(k),nl%dwdx(k,:),kernel_type,dimn)
    call calc_kernel(nl%rij(k), nl%drij(k,:), &
    & (p%sml_long(i)+p%sml_long(j))/2.,nl%w_long(k),nl%dwdx_long(k,:),&
    & kernel_type,dimn)
  end do

  if(verbose) print *,"finished kernels"

end subroutine calc_kernels_all


subroutine calc_smoothed_properties(nl,p,w)
  implicit none
  !uses the f1 sph approximation
  ! <f(x)> = sum( f(i) * m(i)/rho(i) 
  type(nlist_data_type) :: nl
  type(particle_type) :: p
  !kernel for all pairs
  double precision, dimension(:) :: w
    
  double precision :: selfkern
  integer :: i,j,k,ni,d
  double precision :: r
  double precision, dimension(dimn) :: hv
  double precision, dimension(dimn) :: dw 

  hv=0.e0
  r=0
     
  !smoothed velocity
  p%v_flow = 0
  p%p_smooth = 0
  p%t_smooth = 0
   
  ! Self contribution
  do i=1,p%n
    call calc_kernel(r,hv,p%sml(i),selfkern,dw,kernel_type,dimn)
      do d=1,dimn
        p%v_flow(i,d) = selfkern*p%m(i)*p%v(i,d)/p%rho(i)
      end do
      p%p_smooth(i) = selfkern*p%m(i)*(p%p(i)+p%pco(i))/p%rho(i)
      p%t_smooth(i) = selfkern*p%m(i)*p%temp(i)/p%rho(i)
  end do
    
  ! Neighbour contributions
  ni = nl%nip
  do k = 1, ni
    i = nl%ilist(k,1)
    j = nl%ilist(k,2)
 
    do d=1,dimn
      p%v_flow(i,d) = p%v_flow(i,d) + ( p%m(j)/p%rho(j) ) * p%v(j,d) * w(k)
      p%v_flow(j,d) = p%v_flow(j,d) + ( p%m(i)/p%rho(i) ) * p%v(i,d) * w(k)
    end do
    p%p_smooth(i) = p%p_smooth(i) + (p%m(j)/p%rho(j)) * (p%p(j)+p%pco(j)) &
                    & * w(k)
    p%p_smooth(j) = p%p_smooth(j) + (p%m(i)/p%rho(i)) * (p%p(i)+p%pco(i)) &
                    & * w(k)
    p%t_smooth(i) = p%t_smooth(i) + (p%m(j)/p%rho(j)) * p%temp(j) * w(k)
    p%t_smooth(j) = p%t_smooth(j) + (p%m(i)/p%rho(i)) * p%temp(i) * w(k)
    end do
   
end subroutine calc_smoothed_properties


subroutine calc_kernels_ab(p,bp,bnl)
  implicit none
  type (particle_type), intent(inout) :: p
  type (particle_type), intent(inout) :: bp
  type (nlist_data_type) :: bnl
  integer :: k, i, j

  if(verbose) print *,"looping through neighbours"

  do k = 1, bnl%nip
    i = bnl%ilist(k,1)
    j = bnl%ilist(k,2)

    call calc_kernel(bnl%rij(k), bnl%drij(k,:),        &
    &              (p%sml(i)+bp%sml(j))/2,bnl%w(k),   &
    &               bnl%dwdx(k,:),kernel_type,dimn)

  end do
end subroutine calc_kernels_ab


subroutine collide_particles(nl,p)
  implicit none
  type (particle_type), intent(inout) :: p
  type (nlist_data_type), intent(inout) :: nl

  ! nl%rij -- distance
  ! nl%rijsq -- distance sq
  ! nl%drij -- displacement
  integer :: i,j,k
 
  do k = 1, nl%nip
      i = nl%ilist(k,1)
      j = nl%ilist(k,2)
      if (nl%rij(k) .lt. collide_dist) then
        call collide2d(p%v(i,:),p%v(j,:),p%m(i) &
        & ,p%m(j),nl%drij(k,:),nl%rijsq(k) )
      endif
  end do
end subroutine collide_particles

subroutine calc_art_viscosity(dr,dv,rij,rho,sml,c,wij,avisc)
  ! calculates artificial viscosity for one pair of particles
  ! dr -- position difference vector
  ! dv -- velocity difference vector
  ! rij -- pair separation distance
  ! rho -- average density of pair
  ! sml -- average smoothing length of pair
  ! c -- average speed of sound of pair
  ! avisc -- output
  ! The global parameter artificial_viscosity determines the type
  ! of artificial viscosity.
  ! References:
  ! Hoover, page 46, 90. 
  ! Liu, page 125.
  implicit none
  
  double precision, dimension(:), intent(in) :: dr
  double precision, dimension(:), intent(in) :: dv
  double precision, intent(in) :: rho
  double precision, intent(in) :: sml
  double precision, intent(in) :: c
  double precision, intent(in) :: rij
  double precision, intent(in) :: wij
  double precision, intent(inout) :: avisc
  
  if(artificial_viscosity .eq. 0) then

    avisc = 0.0  
  
  else if(artificial_viscosity .eq. 1) then
  
    call vnr_avisc(avisc,dr,dv,rho,sml,wij,rij,dimn)

  else if(artificial_viscosity .eq. 2) then
    
    call mon_beta_avisc(avisc,dr,dv,rho,sml,dimn)
  
  else
    ! We will use c in the future. Need to use it
    ! now to avoid compiler warnings.
    avisc = c * 0
    avisc = 0
  end if

end subroutine calc_art_viscosity

! Selects the appropriate temperature update function
! based on the configuration
subroutine update_temperature(u,temp,rho,grad_rho)
  double precision :: u
  double precision :: temp
  double precision :: rho
  double precision, dimension(:) :: grad_rho
  
  if(eos_type .eq. 1) then
    call igtemp(u,temp)
  else if(eos_type .eq. 2) then
    call vdwtemp(u,temp,rho)
  else if(eos_type .eq. 3) then
    call vdwtemp(u,temp,rho)
  endif
   
end subroutine update_temperature

subroutine update_energy(u,temp,rho,grad_rho)
  ! Selects appropriate energy update function based on configuration
  ! and applies it to the input variables
  double precision :: u
  double precision :: temp
  double precision :: rho
  double precision, dimension(:) :: grad_rho

  if(eos_type .eq. 1) then
     !! call igtemp(p%u(i),p%temp(i))
  else if (eos_type .eq. 2) then
    call vdwenergy(u,temp,rho)
  else if (eos_type .eq. 3) then
    call vdwenergy(u,temp,rho)
  endif

end subroutine update_energy


subroutine update_pressure(rho,u,p,pco,c,grad_rho)
    double precision rho, u, p, c, pco
    double precision, dimension(:) :: grad_rho

    if(eos_type .eq. 1) then
      call igeos(rho,u,p,c)
    else if(eos_type .eq. 2) then
      call vdweos(rho,u,p)
    else if(eos_type .eq. 3) then
      call vdweos_attractive(rho,pco)
      call vdweos_repulsive(rho,u,p)
    endif  

end subroutine update_pressure

subroutine calc_iso_pressure(p)
  !subroutine calc_iso_pressure(rho,u,p,c,grad_rho)
  implicit none
  type (particle_type), intent(inout) :: p
  integer :: i
  !double precision, dimension(:) :: rho
  !double precision, dimension(:) :: u, p, c
  !double precision, dimension(:,:) :: grad_rho
  if (verbose) print *,'calculating iso pressure'
 
  do i=1,p%n
    if(eos_type .eq. 1) then
        call igeos(p%rho(i), p%u(i), p%p(i), p%c(i))
      else if(eos_type .eq. 2) then
        call vdweos2d(p%rho(i), p%u(i), p%p(i))
      else if(eos_type .eq. 3) then
        call vdweos_repulsive2d(p%rho(i),p%u(i),p%p(i))
      !else if(eos_type .eq. 4) then
      !  call vdw_gradient(p%rho(i),p%u(i),p%p(i),p%grad_rho(i,:),p%cgrad(i))
      endif  
  end do
end subroutine calc_iso_pressure

subroutine calc_long_range_pressure(p)
  implicit none
  type (particle_type), intent(inout) :: p
  integer :: i
  !double prghecision, dimension(:) :: rho
  !double precision, dimension(:) :: u, p
 
  do i=1,p%n
      if(eos_type .eq. 3) then
        call vdweos_attractive(p%rho_lr(i),p%pco(i))
      else
        stop 'Invalid equation of state'
      endif  
  end do
end subroutine calc_long_range_pressure



end module particle

