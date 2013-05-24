! Calculates a 'core potential' force for particle pairs
! See Hoover's SPAM book, p238. The core potential is a somewhat
! ad-hoc means of preventing undesirable particle clumping and 
! interpenetration.
! A weakness of this approach is that PV work attributible to this
! interaction is not accounted for.

module core_potential

public :: core_force

! sigma -- range of repulsive potential
! rcoef -- coefficient for repulsive potential
!          NEGATIVE for repulsion

contains

! Calculates the repulsion for one pair of particles
! dr is rj - ri
! so a repulsive force if rj > ri is in the direction of
! dr for j and -dr for i
subroutine core_force(repulsion,dr,drsq,sigma,rcoef,ndim)
  double precision, dimension(ndim) :: dr
  double precision :: drsq
  !f2py intent(in,out,overwrite) :: repulsion
  double precision, dimension(ndim), intent(inout) :: repulsion
  double precision :: sigma_sq
  double precision :: sigma
  double precision :: dx
  double precision :: rcoef
  double precision :: forcemag
  integer :: d, ndim
 
  if (sigma .eq. 0.0) then
    repulsion = 0.0
    return
  endif
  
  sigma_sq = sigma * sigma
 
  if ( drsq .le. (sigma_sq) )  then
    dx=sqrt(drsq)
    forcemag = rcoef*((2*dx)/sigma_sq) * ( ((-drsq/sigma_sq) + 1)**3 )
    do d=1,ndim
      repulsion(d) = forcemag * dr(d)/dx
    end do 
  else
    repulsion = 0
  endif

end subroutine core_force

subroutine calc_coreforce3d(ilist,x,a,m,dr,dx,sigma,rcoef,n,ni)
  ! ilist -- integer array of pairs 
  ! particle variables
      ! x,a -- positions, accelerations (overwritten)
      ! m -- masses
  
  ! dr -- distances between pairs
  ! dx -- displacements between pairs
  ! n -- number of particles
  ! ni -- number of interacting pairs

  implicit none
  integer, intent(in), dimension(ni,2) :: ilist
  double precision, intent(in), dimension(n,3) :: x 
  !f2py intent(in,out,overwrite) :: a
  double precision, intent(inout), dimension(n,3) :: a
  double precision, intent(in), dimension(n) :: m
  double precision, intent(in), dimension(ni) :: dr 
  double precision, intent(in), dimension(ni,3) :: dx 
  double precision, intent(in) :: rcoef
  double precision, intent(in) :: sigma
  integer :: ni
  integer :: n
  integer :: ndim

  !> core repulsion
  double precision, dimension(3) :: repulsion

  integer :: i,j,k

  ndim = 3

  ! set rate of change to zero
  a = 0.0

  do k = 1, ni
    i = ilist(k,1)
    j = ilist(k,2)

    ! Calculate core force
    call core_force(repulsion,dx(k,:),dr(k)*dr(k),sigma,rcoef,ndim)
 
    !! Core repulsion is a straight force so we need to divide by mass.
    
    a(i,1) = a(i,1) + repulsion(1)/m(i)
    a(j,1) = a(j,1) - repulsion(1)/m(j)
    
    a(i,2) = a(i,2) + repulsion(2)/m(i)
    a(j,2) = a(j,2) - repulsion(2)/m(j)

    a(i,3) = a(i,3) + repulsion(3)/m(i)
    a(j,3) = a(j,3) - repulsion(3)/m(j)

  enddo !end loop over pairs

end subroutine calc_coreforce3d


end module core_potential
