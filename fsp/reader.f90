module reader

use global
use eos

private
public :: read_input, read_dbl, read_int, read_bool

contains
  
subroutine read_input(ifile)
  ! Reads data from input file into the global configuration
  ! parameters.

  implicit none
  !character (len=12) :: key
  !double precision :: var 
  character(len=32) :: ifile         ! input file name
  integer iunit
  iunit = 1
  open(iunit,form='formatted',file=ifile,status='old',&
  &    position='rewind')
  call read_int(iunit,'NP',maxn)
  call read_int(iunit,'DIM',dimn)
  call read_int(iunit,'MAXIP',max_interact)
  call read_int(iunit,'TSTEPS',maxsteps)
  call read_dbl(iunit,'DT',dt_ini)
  call read_int(iunit,'SNAPFREQ',snapfreq)
  call read_dbl(iunit,'CUTOFF',cutoff)
  call read_bool(iunit,'DEBUG',debug)
  call read_bool(iunit,'VERBOSE',verbose)
  call read_bool(iunit,'VALIDATE',validate)
  call read_bool(iunit,'CONFIRM',confirm_step)
  call read_dbl(iunit,'XMAX',lx)
  call read_dbl(iunit,'YMAX',ly)
  call read_dbl(iunit,'SIDE',side)
  call read_dbl(iunit,'SPACING',spacing)
  call read_dbl(iunit,'START_TEMP',temperature)
  call read_int(iunit,'ENTROPY_TYPE',entropy_type)
  call read_int(iunit,'INTEGRATOR',integrator)
  call read_int(iunit,'ADAPT_DT',solver_adapt_type)
  call read_dbl(iunit,'TOLERANCE',solver_tolerance)
  call read_int(iunit,'VELOCITY_AVG',velocity_avg)
  call read_dbl(iunit,'V_EPS',v_eps)
  call read_int(iunit,'COLLIDE',collisions)
  call read_dbl(iunit,'COLLIDE_DIST',collide_dist)
  call read_int(iunit,'RUN_TYPE',run_type)
  call read_int(iunit,'ART_VISC',artificial_viscosity)
  call read_int(iunit,'CONDUCTION',conduction_type)
  call read_bool(iunit,'PROFILE',profile)
  call read_dbl(iunit,'CORE_SIZE',sigma)
  call read_dbl(iunit,'REPULSION',rcoef)
  call read_int(iunit,'RESTRICT_V',maximum_velocity)
  call read_int(iunit,'RHO_EQ',density_equation)
  call read_int(iunit,'UNSTABLE_RHO',unstable_density)

! Boundaries
  call read_int(iunit,'NBOUNDS',nbounds)
  call read_int(iunit,'NBP_1',nbp_1)
  call read_int(iunit,'NBP_2',nbp_2)

! Equation of state
  call read_int(iunit,'EOS',eos_type)
  call read_dbl(iunit,'ADASH',adash)
  call read_dbl(iunit,'BDASH',bdash)
  call read_dbl(iunit,'KBDASH',kbdash)
  call read_dbl(iunit,'CGRAD',cgrad_ini)

! Kernel
  call read_int(iunit,'KERNEL',kernel_type)
  call read_dbl(iunit,'H',smoothing_ini)
  call read_int(iunit,'ADAPT_H',adapt_smoothing)
  call read_dbl(iunit,'HCO',co_smoothing_ini)

! Thermostat
  call read_int(iunit,'THERMOSTAT_TYPE',thermostat_type)
  call read_int(iunit,'B_THERMOSTAT_TYPE',bound_thermostat_type)
  call read_bool(iunit,'THERMOSTAT_PARTICLES',thermostat_particles)
  call read_bool(iunit,'THERMOSTAT_BOUNDARY',thermostat_boundaries)
  call read_dbl(iunit,'THERMOSTAT_TEMP',thermostat_temperature)
  if(nbounds .gt. 0) then
    call read_dbl(iunit,'B1_TEMPERATURE',boundary_thermostat_temperature(1))
    call read_dbl(iunit,'B2_TEMPERATURE',boundary_thermostat_temperature(2))
  endif

! Boundary
    call read_bool(iunit,'REFLECT_XMIN',bounds_reflect(1))
    call read_bool(iunit,'REFLECT_XMAX',bounds_reflect(2))
    call read_bool(iunit,'REFLECT_YMIN',bounds_reflect(3))
    call read_bool(iunit,'REFLECT_YMAX',bounds_reflect(4))
    call read_dbl(iunit,'BOUNDARY_CORE',bcore)
    call read_dbl(iunit,'BOUNDARY_SIGMA',bsig) 
    call read_dbl(iunit,'GRAVITY',grav)

  close (iunit)
  max_interact = (maxn*(maxn-1))/2.0
end subroutine read_input


subroutine read_dbl(iunit,key,val)
    ! Searches the input file for the key and reads the value
    ! to the supplied variable when it finds it.
    implicit none
    integer iunit
    character (len=20) :: st
    character (len=*) :: key
    double precision :: val
    integer i
    rewind(iunit)
    do while (st(1:len(key)) .ne. key)
        read(iunit,*,iostat=i) st,val
        !    print *,key," ",st(1:len(key)),i
        if(i .eq. -1) then
            print *,key
            stop 'parameter not found' 
        endif
        if (st(1:len(key)) .eq. key) then
        !   print *,'I found ',st
           write (*,*) st,val
           return 
        endif
    enddo
end subroutine read_dbl

subroutine read_int(iunit,key,val)
    ! Searches the input file for the key and reads the value
    ! to the supplied variable when it finds it.
    implicit none
    integer iunit
    character (len=20) :: st
    character (len=*) :: key
    integer :: val
    integer i
    i=0
    rewind(iunit)
    do while (st(1:len(key)) .ne. key)
        read(iunit,*,iostat=i) st,val
        if(i .eq. -1) then
            print *,key,i
            stop 'parameter not found' 
        endif
        !print *,'key',key,len(key)
        if (st(1:len(key)) .eq. key) then
        !   print *,'I found ',st
           write (*,*) st,val
           return 
        endif
    enddo
end subroutine read_int

subroutine read_bool(iunit,key,val)
    ! Searches the input file for the key and reads the value
    ! to the supplied variable when it finds it.
    implicit none
    integer iunit
    character (len=20) :: st
    character (len=*) :: key
    logical :: val
    integer i,w
    rewind(iunit)
    w = len(key)
    do while (st(1:len(key)) .ne. key)
        !print *,'reading value'
        read(iunit,*,iostat=i) st
        !print *,i
        !print *,st,val
        !    print *,key,st(1:len(key)),i
        if(i .eq. -1) then
            write(*,*) key
            stop 'fatal error: parameter not found' 
        endif
        if (st(1:len(key)) .eq. key) then
        !   print *,'I found ',st
           backspace(iunit)
           read(iunit,*,iostat=i) st,val
           write(*,*) st,val
           return 
        endif
    enddo
end subroutine read_bool


end module reader 
