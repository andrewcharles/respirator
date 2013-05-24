!--------1---------2---------3---------4---------5---------6---------7--------|
!! A Smoothed Particle continuum equation solver
!! Andrew Charles, RMIT 2006-2008.
!! 
!! Three dimensional version. No boundaries.
!!
!! This is not complete. Not all of the lower level subroutines support 3d
!!
!! With contributions from
!! Assoc Prof Peter Daivis, RMIT.
!!
!! Copyright Andrew Charles, RMIT University 
!! this code is released under the Computer Physics
!! Communications license. 
!! All other rights reserved.

program fsph

use global
use sphstep
use particle
use reader
use writer
use simulation_box
use neighbour_list
use thermostat
use boundary
use kernel
use adaptive_timestep

implicit none

!! global variable declarations
! p -- particle properties
! nl -- neighbour list contains pair properties 
! box -- simulation box
! tstat -- config data for the thermostat
! fp -- file pointer for input files
! error --
! input_file -- input filename

type (particle_type)    :: p          
type (particle_type), dimension(:), allocatable :: bp     
type (nlist_data_type)  :: nl         
type (nlist_data_type), dimension(:), allocatable :: bnl     
type (box_data_type)   :: box        

character(len=32) :: input_file      

! loop counters
integer :: i

! benchmarking
integer :: start_time,end_time,time_rate
real :: time_seconds

! Variable initialisations
! ------------------------
call system_clock(count=start_time,count_rate=time_rate)
input_file = 'spinput3d.txt'
call read_input(input_file)

print *,'***********************************************'
print *,'** A program for solving problems in         **'
print *,'** fluid dynamics using                      **'
print *,'** Smoothed Particle Applied Mechanics       **'
print *,'** Andrew Charles, RMIT. All rights reserved **'
print *,'** Contact ac1201@gmail.com                  **'
print *,'***********************************************'

print *,'System Parameters:'
print *,maxn,' particles'
print *,dimn,' dimensions'

! Initialise the particle data structures
! create_particles is in .... module
! initialise_particles is in ... module
call create_particles(p,maxn,dimn)
call initialise_particles(p)
call rd_particle_input(p,input_file,1)
if(validate) call validate_particles(p)
print *,'Particle data structure initialised ...'

! Set the thermostats up
!call rd_thermostat_input(input_file,1)
p%thermostat = temperature
print *,'Thermostat initialised'

! Initialise the simulation box
call rd_box_input(box,input_file,1)
call create_box(box,dimn)
call init_box(box, lx, ly)
print *,'Box initialised'

! Initialise the neighbour list
call rd_nlist_input(nl,input_file,1)
call init_nlist(nl, cutoff, p%n, dimn)
print *,'Neighbour list initialised ...'

! Initialise the system
step = 0
dt = dt_ini
call initialise(p,nl,box)

! Initial energies
call system_kinetic_energy(p,ek)
call system_internal_energy(p,ep)
ep_iso = ep + u_env
isoham_initial = ek+ep_iso
isoham = isoham_initial

if(profile) then
  call system_clock(count=end_time)
  time_seconds = (end_time - start_time)/real(time_rate)
  print *,'Initialisation time (seconds)', time_seconds
endif

!------------------
! Main Program Loop
!------------------

do i=1,maxsteps

  if(debug) print *,'****  Starting loop number ',i
  call system_clock(count=start_time)
  
  if(verbose) print *,'Writing data' 
  call write_properties(p,maxn)
  
  if( mod(step,snapfreq) .eq. 0 ) call write_state(p,step,bp)

  if(integrator .eq. 1) then
    call sph_step_imp_euler(nl,p,dt,       &
                          & box,bnl,bp)
  else if (integrator .eq. 2) then
    call sph_step_leapfrog(nl,p,dt,         &
                          & box,bnl,bp)
  else if (integrator .eq. 3) then
    call sph_step_runge_kutta(nl,p,dt,       &
                            & box,bnl,bp)
  end if

  ! compute total system energies (including boundaries)
  call system_kinetic_energy(p,ek,bp)
  call system_internal_energy(p,ep,bp)
  ep_iso = ep + u_env
  solver_energy_error = abs((ep_iso + ek - isoham)/ isoham) 
  isoham = ek+ep_iso

  !print *,'pot',ep_iso
  !print *,'ek',ek
  !print *,'sum',isoham

  if (solver_adapt_type .eq. 0) then
  else if (solver_adapt_type .eq. 1) then
    print *,solver_energy_error
    print *,isoham
!   if the error exceeds the error tolerance
!   reduce the timestep
    if (solver_energy_error .gt. solver_tolerance) then 
      dt = 0.5 * dt
      print *,"shrinking dt"
      print *,solver_energy_error
      print *,dt
    else if (solver_energy_error .lt. solver_tolerance/1000.0) then
      print *,"increasing dt"
      print *,solver_energy_error
      print *,dt
      dt = 1.01 * dt
    endif
  else if (solver_adapt_type .eq. 2) then

    call adapt_dt_velocity(p,dt)

  endif

  if(verbose) print *,'Checking velocity' 
  call check_velocity(p,p%dimn,p%n,dt)
  
  call system_clock(count = end_time)
  solver_sdt = (end_time - start_time)/real(time_rate)
  if(solver_sdt<0) print *,'diurnal boundary crossing -timing not accurate'
  write(*,'(8A)',advance='no') '+-> Step '
  write(*,'(i4)',advance="no") i
  write(*,'(16A)',advance='no') ' execution time: '
  write(*,'(F10.5)') solver_sdt

  if(debug) then
    print *,'*******************'
    print *,'**** Completed step', i
  endif

  step = step + 1

  if(confirm_step) then
    print *,'Hit enter to continue'
    read (*,*)
  endif

enddo

!---------------
! Main loop ends
!---------------

! Cleanup
! -------
call destroy_box(box)
call destroy_nlist(nl)
call destroy_particles(p)

do i=1,nbounds
  call destroy_particles(bp(i))
  call destroy_nlist(bnl(i))
end do

print *, 'I sincerely believe all memory has been freed.'
print *, 'Smooth particle program ended normally. You recieve an item [cookie]'

contains

subroutine initialise(p,nl,box,bp,bnl)
  !! Initialise the system:
  !!  - assign particle positions and temperatures
  !!  - computing initial neighbour list
  !!
  !! Arguments
  !! ---------
  !! p -- particles
  !! nl -- neighbours
  !! bx -- simulation box
  !! bp -- boundary particles
  !! bnl -- boundary neighbours
  !!
  !! Global configuration options
  !! ----------------------------
  !! if run_type is 1 load the state from file
  !! if run_type is 2 then load the state, but over-write the temperature

  use reader
  use writer
  implicit none

  type (particle_type), intent(inout) :: p
  type (nlist_data_type), intent(inout) :: nl
  type (box_data_type) :: box
  type (particle_type), dimension(:), optional, intent(inout) :: bp
  type (nlist_data_type), dimension(:), optional, intent(inout) :: bnl

  integer :: n
  integer i, q
  double precision j,k

  if (p%dimn .eq. 2) then
    ! arrange the particles in a central grid
    j=(box%boxvec(1,1)/2.0 - side/2.0)
    k=(box%boxvec(2,2)/2.0 - side/2.0)
     
    n = p%n 
    do i=1,p%n
      p%sml(i) = smoothing_ini
      p%sml_long(i) = co_smoothing_ini
      p%x(i,1) = j
      p%x(i,2) = k
      p%temp(i) = temperature 
      if(debug) print *, 'set position of',i, 'to', j,k
      j=j+spacing
      if (j .ge. side+(box%boxvec(1,1)/2 - side/2)) then
        k=k+spacing
        j=(box%boxvec(1,1)/2 - side/2)
      end if
    end do
    p%xhalf = p%x
  end if
      
  if(present(bp)) then
    do i=1,nbounds
      bp(i)%sml = smoothing_ini
      bp(i)%sml_long = co_smoothing_ini
    end do
  endif
    
  if(run_type .eq. 1) then
    print *,'Loading particle state from file...'
    call read_particle_state(p)
    p%xhalf = p%x
  else if(run_type .eq. 2) then
    print *,'Loading particle state from file...'
    call read_particle_state(p)
    p%xhalf = p%x
    print *,'Setting temperatures...'
    do i=1,p%n
      p%temp(i)=temperature
    end do 
  else if(run_type .eq. 3) then
    print *,'Loading particle state from file...'
    call read_particle_state(p)
    p%xhalf = p%x
    p%sml(:) = smoothing_ini
    p%sml_long(:) = co_smoothing_ini
  end if    

  call form_nlist(nl, box, p%x)
  if(verbose) print *,'formed initial neighbour list'

  if (present(bnl)) then 
    do i=1,nbounds
      call form_nlist(bnl(i), box, p%x, bp(i)%x)
      if(verbose) print *,'formed boundary neighbour list'
    end do
  endif
  
  if(verbose) print *,'finding pair separations'
  call find_pair_separations(nl,box,p%x)
  
  if(verbose) print *,'Compressing neighbour list'
  call compress_nlist(nl)

  if (present(bnl)) then 
    do i=1,nbounds
      call find_pair_separations(bnl(i),box,p%x,bp(i)%x)
      call compress_nlist(bnl(i))
    end do
  endif


  if(verbose) print *,'found pair separations - calculating kernel'
    
  call calc_kernels_all(p,nl)

  if(validate) call validate_particles(p)
  if(validate) call validate_nlist(nl)

  if(verbose) print *,'Calculated initial kernels, summing density...'
  if(verbose) then
      print *,size(p%sml),size(nl%nlist(:,1)),size(nl%ilist(:,1))
      !print *,size(nl%w),size(bp%rho),size(p%grad_rho),size(nl%dwdx)
      print *,n,nl%nip
  endif 
  call sum_density(nl, n, p%sml, p%m, nl%w, p%rho, p%grad_rho,nl%dwdx)

  if(verbose) print *,'Computing energy from equation of state...'
  do i=1,p%n
    call update_energy(p%u(i),p%temp(i),p%rho(i),p%grad_rho(i,:))
  end do
  
  if(present(bp)) then

    if(validate) call validate_particles(bp(1))
    if(validate) call validate_nlist(bnl(1))

    do i=1,nbounds
      if(verbose) print *,'Computing boundary energy from equation of state'
      do q=1,bp(i)%n
        call update_energy(bp(i)%u(q),bp(i)%temp(q),bp(i)%rho(q),bp(i)%grad_rho(q,:))
      end do
    end do
  endif

  do i=1,p%n
    call update_pressure(p%rho(i),p%u(i),p%p(i),p%pco(i),p%c(i),p%grad_rho(i,:))
  end do

  if(present(bp)) then
  if(verbose) print *,'computing boundary pressure from equation of state'
  do i=1,nbounds
    do q=1,bp(i)%n
      call update_pressure(bp(i)%rho(q),bp(i)%u(q),bp(i)%p(q),bp(i)%pco(q),bp(i)%c(q),bp(i)%grad_rho(q,:))
    end do
  end do
  endif
  
  if(validate) call validate_particles(p)
  !if(validate) call validate_particles(bp(1))
  if(verbose) print *,'Successfully initialised'
  

end subroutine initialise

end program fsph

