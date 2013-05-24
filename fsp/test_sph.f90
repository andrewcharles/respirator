! Exercises functionality of two dimensional modules 
!

program test_sph

use global
use sphstep
use kernel
use particle
use reader
use writer
use sphforce
use density
use simulation_box
use neighbour_list
use thermostat
use boundary
use adaptive_timestep

implicit none

! Initialise a particle object, nlist object and a box

type (particle_type)    :: p          
type (nlist_data_type)  :: nl         
type (box_data_type)   :: box        
character(len=32) :: input_file
integer :: i

print *,'TESTING SPAC - no data will be output'

input_file = 'spinput.txt'
call read_input(input_file)

call create_particles(p,maxn,dimn)
call initialise_particles(p)
call rd_particle_input(p,input_file,1)
call validate_particles(p)
p%thermostat = temperature

call rd_box_input(box,input_file,1)
call create_box(box,dimn)
call init_box(box, lx, ly)

call rd_nlist_input(nl,input_file,1)
call init_nlist(nl, cutoff, p%n, dimn)

step = 0
dt = dt_ini
call initialise(p,nl,box)

do i=1,maxsteps

  if(debug) print *,'****  Starting loop number ',i
  
  if(integrator .eq. 1) then
    call sph_step_imp_euler(nl,p,dt,       &
                          & box)
  else if (integrator .eq. 2) then
    call sph_step_leapfrog(nl,p,dt,         &
                          & box)
  else if (integrator .eq. 3) then
    call sph_step_runge_kutta(nl,p,dt,       &
                            & box)
  end if

  write(*,'(8A)',advance='no') '+-> Step '
  write(*,'(i4)',advance="no") i
  write(*,'(16A)',advance='no') ' execution time: '

  print *,'*******************'
  print *,'**** Completed step', i

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

print *, 'I sincerely believe all memory has been freed.'
print *, 'Smooth particle program ended normally. You recieve an item [cookie]'

contains

subroutine initialise(p,nl,box)
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
  integer :: n
  integer i, q
  double precision j,k

  if (p%dimn .eq. 2) then
    ! arrange the particles in a central grid
    j=(box%boxvec(1,1)/2.0 - side/2.0)
    k=(box%boxvec(2,2)/2.0 - side/2.0)
     
    do i=1,p%n
      p%sml(i) = smoothing_ini
      p%sml_long(i) = co_smoothing_ini
      p%x(i,1) = j
      p%x(i,2) = k
      p%eta=1.0
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

  if(verbose) print *,'finding pair separations'
  call find_pair_separations(nl,box,p%x)

  if(verbose) print *,'Compressing neighbour list'
  call compress_nlist(nl)
  call find_pair_separations(nl,box,p%x)
  call compress_nlist(nl)
  call find_pair_separations(nl,box,p%x)
  call compress_nlist(nl)
  call find_pair_separations(nl,box,p%x)
  call compress_nlist(nl)

  if(verbose) print *,'found pair separations - calculating kernel'
    
  call calc_kernels_all(p,nl)

  if(validate) call validate_particles(p)
  if(validate) call validate_nlist(nl)

  if(verbose) print *,'Calculated initial kernels, summing density...'
  if(verbose) then
      print *,size(p%sml),size(nl%nlist(:,1)),size(nl%ilist(:,1))
      print *,n,nl%nip
  endif 
  call sum_density(nl, p%n, p%sml, p%m, nl%w, p%rho, p%grad_rho,nl%dwdx)

  if(verbose) print *,'Computing energy from equation of state...'
  do i=1,p%n
    call update_energy(p%u(i),p%temp(i),p%rho(i),p%grad_rho(i,:))
  end do
  
  do i=1,p%n
    call update_pressure(p%rho(i),p%u(i),p%p(i),p%pco(i),p%c(i),p%grad_rho(i,:))
  end do

  if(validate) call validate_particles(p)
  if(verbose) print *,'Successfully initialised'
  if(verbose) print *,'************************'
  

end subroutine initialise

end program test_sph
