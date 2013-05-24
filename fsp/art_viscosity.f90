! Calculates artificial viscosity for pairs of particles.
! 1. Simple von Neumann-Richtmeyer with linear and quadratic terms.
! 2. Monaghan's beta term viscosity.
! All primitive arguments, no dependencies.

module art_viscosity

implicit none
public :: mon_beta_avisc, vnr_avisc 

contains

subroutine vnr_avisc(avisc,dr,dv,rho,sml,wij,rij,ndim)
  ! Simple von Neumann Richtmeyer quadratic artifical viscosity
  ! See Hoover for reference.
  ! To include a linear term, the following expression may be
  ! useful:
  ! avisc = (alpha * dx * c * rho * gradv) 
  ! + (alpha2 * drsq * rho * gradv*gradv)
  implicit none
  double precision, intent(inout) :: avisc
  double precision, dimension(ndim), intent(in) :: dr
  double precision, dimension(ndim), intent(in) :: dv
  double precision, intent(in) :: rho
  double precision, intent(in) :: sml
  double precision, intent(in) :: rij
  double precision, intent(in) :: wij
  integer, intent(in) :: ndim
  double precision :: alpha, alpha2
  double precision :: div_v 
  integer :: d
  alpha = 1
  alpha2 = 1
  div_v = 0
  do d=1,ndim
    div_v = div_v + dv(d)*dr(d)
  end do
  if( div_v < 0.0 ) then
    div_v = div_v * wij / rij
    avisc =  (alpha2 * sml * sml * rho * div_v  *div_v)
  else
    avisc = 0
  endif
end subroutine vnr_avisc

subroutine mon_beta_avisc(avisc,dr,dv,rho,sml,ndim)
  ! Monaghan beta term for SPH artifical viscosity
  ! The speed of sound is required to include the
  ! alpha term. It is considered (see Liu) that the alpha
  ! term is not required because we calculate shear 
  ! viscosity directly.
  ! The full monaghan artifical viscosity is
  ! avisc = (beta*phi*phi - alpha*phi*c)/rho. 
  implicit none
  double precision, dimension(ndim), intent(in) :: dr
  double precision, dimension(ndim), intent(in) :: dv
  double precision, intent(in) :: rho
  double precision, intent(in) :: sml
  double precision, intent(inout) :: avisc
  double precision :: vdotr
  double precision :: drsq, phi
  double precision :: beta, kippa
  integer :: d,ndim

  ndim = size(dr)

  beta = 1.0
  kippa = 0.1
  vdotr = 0.0
  drsq = 0.0

  avisc = 0.0
  do d=1,ndim
    vdotr = vdotr + dv(d)*dr(d)
    drsq = drsq + dr(d)*dr(d)
  end do
  if( vdotr < 0.0 ) then
    phi = sml * vdotr /(drsq + (kippa*sml)**2) 
    avisc = (beta*phi*phi)/rho
  end if
end subroutine mon_beta_avisc

end module art_viscosity
