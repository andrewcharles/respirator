! Adaptive timestep for the SPAC fortran code.
!
! Andrew Charles 2008

module adaptive_timestep
use global
use particle

implicit none
! Velocity based adaptive timestep.
! Given the particle velocities and the smoothing length
! Set the timestep such that no particle may move more than 0.25h
! For simplicity we use the particle itself's smoothing length
!
! We scale the timestep such that the particle will move only 0.25h
! This is actually pretty far, and probably needs to be restricted
! further.

contains

subroutine adapt_dt_velocity(p,dt)
  type (particle_type) :: p
  double precision :: dt
  double precision :: drsq,visq,v,vbsq
  integer :: i
  integer :: b

  vbsq = 0.0 
  ! Find the fastest moving particle 
  ! Sure some sort of sort based search would be faster.
  ! If this bugs any of you comp sci geeks then send me a patch.
  do i=1,p%n
    visq = p%v(i,2)**2 + p%v(i,1)**2
    if (visq .gt. vbsq) then 
      vbsq = visq
      b = i
    endif
  enddo

  drsq = (p%v(b,2)*dt)**2 + (p%v(b,1)*dt)**2 
  if( drsq > ((p%sml(b)/4.)**2) ) then
      v = sqrt(vbsq)
      dt = (p%sml(b)/4.)/v
  endif
    

end subroutine adapt_dt_velocity


end module adaptive_timestep
