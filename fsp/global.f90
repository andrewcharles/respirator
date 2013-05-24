! Global parameters for configuration.

module global

implicit none
! snapfreq -- Snapshot frequency (simulation steps)
! step -- Current step count
! side -- Side of initial particle grid (length units)
! spacing -- Spacing of particles in initial grid
! temperature -- Initial temperature of particles
! lx,ly -- Size of box in x,y dimension
! velocity_avg --
! v_eps --
! rho_max --

integer :: snapfreq
logical :: replace
integer :: step                 
double precision side           
double precision :: spacing     
double precision :: temperature 
double precision lx, ly         
integer :: velocity_avg         
double precision :: v_eps      
double precision :: rho_max
parameter(rho_max=2)

integer :: nbounds    !< number of sets of boundary particles
integer :: nbp_1 !< Number of boundary particles in first system
integer :: nbp_2 !< Number of boundary particles in second system

!! SYSTEM ENERGIES
!! ep_iso -- isolated internal energy - used for  energy conservation
!! ek -- kinetic energy of system ep_iso + ek should be constant
!! u_env --  energy budget of environment
double precision :: ep_iso
double precision :: ek
double precision :: ep
double precision :: u_env
double precision :: isoham
double precision :: isoham_initial

integer :: maxn !< The maximum number of particles
integer :: dimn !< The dimension of the system
!parameter(dimn=2)

! Number of neighbour lists
! integer :: num_lists

!> what to do about unstable density
!! 1=stop, 2=restrict to stable region
integer :: unstable_density

!> equation of state type
! 1 = ideal gas,
! 2 = vanderwaals (combined)
! 3 = van der waals (seperated)
integer :: eos_type

!> density gradient constant
double precision :: cgrad_ini

!> range of repulsive potential
double precision :: sigma

!> coefficient for repulsive potential
double precision :: rcoef

double precision :: grav

!> What to do when particle velocities get beyond our resolution
! 0 = no restrictions
! 1 = stop simulation if v*dt < h/2
integer :: maximum_velocity

!> type of heat conduction
!> 1 = Monaghan style
!> 2 = Full heat flux
integer :: conduction_type

!> viscous entropy type
!! shear tensor type is probably a more appropriate name
!! 1 = Liu approximation
!! 2 = summantion (ref sigalotti)
integer :: entropy_type

!> hard collisions?
integer :: collisions
!> collision radius
double precision :: collide_dist

integer :: density_equation
! 0=summation density
! 1=continuity density

integer :: artificial_viscosity
! 0 = none
! 1 = von Neumann-Richtmyer (linear)only
! 2 = monaghan beta term

integer :: max_interact !< The maximum number of interacting pairs
double precision :: cutoff !< distance at which to truncate the force
logical :: debug   !< debug output on/off
logical :: verbose !< insane amount of verbous output on/off
logical :: profile !< time profile output on/off
logical :: validate !< mostly silent validation tests run in code

!> thermal conductivity
!! \todo looks like thermal k is set as a parameter
double precision :: thermalk
parameter(thermalk = 5.0)

logical :: confirm_step !< whether to prompt for input before each step

!> Type of run
! 0 = new run;
! 1 = from stored state sphstate.input.
integer :: run_type

!> compare monaghan and full flux heat terms (0,1)
integer :: compare_heat

!> physical and numerical constants
!! yeah yeah this is not resolved enough
double precision :: pi
parameter(pi = 3.1415926535897931)

! Solver Configuration
!
! solver_type -- Integration Method
!   1 = Improved Euler
!   2 = Leapfrog
!   3 = Fourth order Runge-Kutta
!
! solver_adapt_type -- adaptive timestep type
!   0 = None
!   1 = Energy Conservation
!   2 = Velocity Limiting
!
! maxsteps -- number of steps to advance the simulation
! dto -- initial timestep size
! dt -- current timestep size
! solver_sdt -- time taken to execute last step

integer :: integrator
integer :: solver_adapt_type
double precision :: solver_tolerance
double precision :: solver_energy_error
integer :: maxsteps
double precision :: dt_ini
double precision :: dt
double precision :: solver_sdt


! Kernel Configuration
! --------------------
! kernel_type -- kernel function
!   1: gaussian
!   2: Lucy quintic
! smoothing -- default smoothing length
! adapt_smoothing -- whether to adapt the smoothing length
!   0 = none
!   1 = h=ho(rho_o/rho)^(1/d)
!   2 = sig ADKE
! co_smoothing -- cohesive smoothing length
integer :: kernel_type
double precision :: smoothing_ini
integer :: adapt_smoothing
double precision :: co_smoothing_ini
logical :: normalise_kernel
parameter(normalise_kernel = .false.)

! Thermostat Configuration
! ------------------------
double precision :: thermostat_temperature
integer :: thermostat_type        
integer :: bound_thermostat_type        
double precision :: thermostat_scaling_factor
logical :: thermostat_particles
logical :: thermostat_boundaries
double precision, dimension(2) :: boundary_thermostat_temperature

! Boundary Configuration
! ----------------------
logical, dimension(4) :: bounds_reflect !1- xmin, 2-xmax, 3-ymin, 4-ymax
double precision :: bcore
double precision :: bsig
double precision :: bconfig

contains 

subroutine validate_double(x,description)
  character(len=*) :: description         
  double precision :: x
  if ( x == x ) then
    x=x
  else
    print *,description
    stop 'NaN detected' 
  end if 
end subroutine validate_double


end module global
