! Equations of state 
! Functional relationships between state variables
! Subroutines calculate the value of either pressure, density or
! temperature given the other two. For the ideal gas equation sound
! velocity is a parameter.
! Note the dimensional dependence of the coefficient on the
! temperature relation.

module eos

implicit none

! parameters for equation of state
double precision :: adash = 2.0
double precision :: bdash = 0.5
double precision :: kbdash = 1.0

public :: igeos, igtemp, vdweos, vdweos_attractive, vdweos_repulsive, &
        & vdwenergy, vdwtemp, calc_vdw_energy, calc_vdw_temp,         &
        & calc_vdw_cohesive_pressure, calc_vdw_hc_pressure,           &
        & eos_set_vdw_params, eos_print_vdw_params,                   &
        & vdweos2d, vdwenergy2d, vdwtemp2d

contains

subroutine eos_set_vdw_params(adash_in,bdash_in,kbdash_in)
  ! Sets the parameters in the imported eos module
  ! These are not accessible via the interface that f2py builds
  double precision, intent(in) :: adash_in 
  double precision, intent(in) :: bdash_in
  double precision, intent(in) :: kbdash_in
  adash = adash_in
  bdash = bdash_in
  kbdash = kbdash_in
end subroutine eos_set_vdw_params

subroutine eos_print_vdw_params()
  ! This function enables us to inspect the values
  print *,adash
  print *,bdash
  print *,kbdash
end subroutine eos_print_vdw_params

subroutine igeos(rho, u, p, c)
  ! Ideal gas equation of state
  ! density, internal energy, pressure, sound velocity
  implicit none
  double precision rho, u, p, c
  double precision gam
  ! air, gamma = 1.4
  gam = 1.4
  p = (gam-1) * rho * u
  c = sqrt(gam-1)*u
end subroutine igeos


subroutine igtemp(u,t)
  ! Incomplete
  implicit none
  double precision u, t, cv
  cv = 7.2
  t = u/cv
end subroutine igtemp


subroutine vdweos(rho,u,p)
  !! van der Waals equation of state - computes pressure
  !! Uses variables from global module
  !! adash, bdash and kbdash
  implicit none
  double precision rho,u,p,t

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
  t = (u+adash*rho)/kbdash
  !t = (2./3.)*(u+adash*rho)/kbdash
  p = (rho*kbdash*t)/(1-rho*bdash) - adash*rho*rho
  !endif
end subroutine vdweos

subroutine vdweos2d(rho,u,p)
  !! van der Waals equation of state - computes pressure
  !! Uses variables from global module
  !! adash, bdash and kbdash
  implicit none
  double precision rho,u,p,t

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
  t = (u+adash*rho)/kbdash
  p = (rho*kbdash*t)/(1-rho*bdash) - adash*rho*rho
  !endif
end subroutine vdweos2d


subroutine vdweos_attractive(rho,p)
  ! Attractive part of van der Waals pressure
  ! for seperated equation of state (ref Nugent and Posch, 2000)
  !f2py intent(in,out,overwrite) :: p
  implicit none
  double precision rho
  double precision, intent(inout) :: p
  p =  - adash*rho*rho
end subroutine vdweos_attractive


subroutine vdweos_repulsive(rho, u, p)
  ! Attractive part of van der Waals pressure
  ! for seperated equation of state (ref Nugent and Posch, 2000)
  ! Removed the test for unstable density - this should be done at
  ! a higher level in the code
  !f2py intent(in,out,overwrite) :: p
  implicit none
  double precision rho, u, t
  double precision, intent(inout) :: p
  t = (u+adash*rho)/kbdash
  !t = (2./3.)*(u+adash*rho)/kbdash
  p = (rho*kbdash*t)/(1-rho*bdash)
end subroutine vdweos_repulsive

subroutine vdweos_repulsive2d(rho, u, p)
  ! Attractive part of van der Waals pressure
  ! for seperated equation of state (ref Nugent and Posch, 2000)
  ! Removed the test for unstable density - this should be done at
  ! a higher level in the code
  !f2py intent(in,out,overwrite) :: p
  implicit none
  double precision rho, u, t
  double precision, intent(inout) :: p
  t = (u+adash*rho)/kbdash
  p = (rho*kbdash*t)/(1-rho*bdash)
end subroutine vdweos_repulsive2d

subroutine vdwtemp(u,t,rho)
  ! Caloric equation of state (computes temperature) for van der Waals
  ! model.
  !f2py intent(in,out,overwrite) :: t
  implicit none
  double precision, intent(in) :: u
  double precision, intent(inout) :: t
  double precision, intent(in) :: rho
  t = (u+adash*rho)/kbdash
  !t = (2./3.)*(u+adash*rho)/kbdash
end subroutine vdwtemp

subroutine vdwtemp2d(u,t,rho)
  ! Caloric equation of state (computes temperature) for van der Waals
  ! model.
  !f2py intent(in,out,overwrite) :: t
  implicit none
  double precision, intent(in) :: u
  double precision, intent(inout) :: t
  double precision, intent(in) :: rho
  t = (u+adash*rho)/kbdash
end subroutine vdwtemp2d

subroutine calc_vdw_temp(u,t,rho,n)
  ! Caloric equation of state for arrays
  !f2py intent(in,out,overwrite) :: t
  implicit none
  double precision, intent(in), dimension(n) :: rho 
  double precision, intent(inout), dimension(n) :: t
  double precision, intent(in), dimension(n) :: u 
  integer :: i,n
  do i=1,n
    call vdwtemp(u(i),t(i),rho(i))
  end do
end subroutine calc_vdw_temp

subroutine vdwenergy(u,t,rho)
  ! Caloric equation of state (computes internal energy) for van der Waals
  ! model.
  !f2py intent(in,out,overwrite) :: u
  implicit none
  ! Here is an f2py gotcha. Right:
  double precision, intent(inout) :: u
  double precision, intent(in) :: t, rho
  ! Wrong! The compiler will not return your value!
  !double precision :: u, t, rho
  u = t*kbdash -adash*rho
  !u = (3./2.)*t*kbdash -adash*rho
end subroutine vdwenergy

subroutine vdwenergy2d(u,t,rho)
  ! Caloric equation of state (computes internal energy) for van der Waals
  ! model.
  !f2py intent(in,out,overwrite) :: u
  implicit none
  ! Here is an f2py gotcha. Right:
  double precision, intent(inout) :: u
  double precision, intent(in) :: t, rho
  ! Wrong! The compiler will not return your value!
  !double precision :: u, t, rho
  u = t*kbdash -adash*rho
end subroutine vdwenergy2d

subroutine calc_vdw_energy(u,t,rho,n)
  ! Caloric equation of state for arrays
  !f2py intent(in,out,overwrite) :: u
  implicit none
  double precision, intent(inout), dimension(n) :: u 
  double precision, intent(in), dimension(n) :: rho 
  double precision, intent(in), dimension(n) :: t
  integer :: i,n
  do i=1,n
    call vdwenergy(u(i),t(i),rho(i))
  end do
end subroutine calc_vdw_energy


subroutine calc_vdw_hc_pressure(p,rho,u,n)
  ! Compute van der Waals equation of state for the whole system
  ! Just the hard core bit
  implicit none
  double precision, intent(inout), dimension(n) :: p 
  !f2py intent(in,out,overwrite) :: p
  double precision, intent(in), dimension(n) :: rho 
  double precision, intent(in), dimension(n) :: u 
  integer, intent(in) :: n
  integer :: i

  do i=1,n
    call vdweos_repulsive(rho(i),u(i),p(i))
  end do

end subroutine calc_vdw_hc_pressure


subroutine calc_vdw_cohesive_pressure(pco,rho,u,n)
  ! Compute van der Waals equation of state for the whole system
  ! Just the hard core bit
  implicit none
  double precision, intent(inout), dimension(n) :: pco 
  !f2py intent(in,out,overwrite) :: pco 
  double precision, intent(in), dimension(n) :: rho 
  double precision, intent(in), dimension(n) :: u 
  integer, intent(in) :: n

  integer :: i

  do i=1,n
    call vdweos_attractive(rho(i),pco(i))
  end do

end subroutine calc_vdw_cohesive_pressure


end module eos
