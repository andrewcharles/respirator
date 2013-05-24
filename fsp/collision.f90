! Incomplete implementation of hard-sphere elastic collision between
! particles

module collision

implicit none
! A collision reassigns momentum among particles that are too close together
! A hard collision is just an instantanous force. 
! The subroutine changes v directly.        
contains

subroutine collide2d(va,vb,ma,mb,dx,rsq)
  ! Collides two particles, elastically reflecting
  ! the normal component of their velocities.
  ! Does not test for closeness, it is assumed this is
  ! done before calling.
  implicit none
  double precision, intent(inout), dimension(2) :: va,vb
  double precision, intent(in) :: ma,mb
  double precision, intent(in), dimension(2) :: dx
  double precision, intent(in) :: rsq
  double precision :: vidotr
  double precision :: vjdotr 
  double precision :: vtix, vtiy, vnix, vniy
  double precision :: vtjx, vtjy, vnjx, vnjy
  double precision :: tmp
  
  ! Compute the dot products to determine direction
  ! of travel
  vidotr = va(1)*dx(1) + va(2)*dx(2)
  vjdotr = vb(1)*dx(1) + vb(2)*dx(2)

  ! If ther particles are moving away from each other 
  ! our work here is done
  if( (vidotr .lt. 0) .and. (vjdotr .gt. 0) ) return

  ! Compute tangential and normal components of velocity
  vtix = (vidotr/rsq) * dx(1) 
  vtiy = (vidotr/rsq) * dx(2)
  vnix = va(1) - vtix
  vniy = va(2) - vtiy
  vtjx = (vjdotr/rsq) * dx(1) 
  vtjy = (vjdotr/rsq) * dx(2)
  vnjx = vb(1) - vtjx
  vnjy = vb(2) - vtjy     
! Transfer tangential component of momentum
  tmp = vtix
  vtix = vtjx *(mb/ma)
  vtjx = tmp * (ma/mb)
  tmp = vtiy
  vtiy = vtjy *(mb/ma)
  vtjy = tmp * (ma/mb)
! Convert back to xy frame
  va(1) = vtix + vnix 
  va(2) = vtiy + vniy
  vb(1) = vtjx + vnjx
  vb(2) = vtjy + vnjy

end subroutine collide2d


subroutine collide3d(va,vb,ma,mb,dx,rsq)
  ! Collides two particles, elastically reflecting
  ! the normal component of their velocities.
  ! Does not test for closeness, it is assumed this is
  ! done before calling.
  !
  ! vidotscl -- vidotr/rsq
  !
  implicit none
  double precision, intent(inout), dimension(3) :: va,vb
  double precision, intent(in), dimension(3) :: dx
  double precision, intent(in) :: ma,mb
  double precision, intent(in) :: rsq
  double precision :: vidotr
  double precision :: vjdotr 
  double precision :: vtix, vtiy, vtiz, vnix, vniy, vniz
  double precision :: vtjx, vtjy, vtjz, vnjx, vnjy, vnjz
  double precision :: tmp, vdotscl
  
  ! Compute the dot products to determine direction
  ! of travel
  vidotr = va(1)*dx(1) + va(2)*dx(2) + va(3)*dx(3)
  vjdotr = vb(1)*dx(1) + vb(2)*dx(2) + vb(3)*dx(3)

  ! If ther particles are moving away from each other 
  ! our work here is done
  if( (vidotr .lt. 0) .and. (vjdotr .gt. 0) ) return

  ! Compute tangential and normal components of velocity
  vdotscl = (vidotr/rsq)
  vtix = vdotscl * dx(1) 
  vtiy = vdotscl * dx(2)
  vtiz = vdotscl * dx(3)
  vnix = va(1) - vtix
  vniy = va(2) - vtiy
  vniz = va(3) - vtiz
  vdotscl = (vjdotr/rsq)
  vtjx = vdotscl * dx(1) 
  vtjy = vdotscl * dx(2)
  vtjz = vdotscl * dx(3)
  vnjx = vb(1) - vtjx
  vnjy = vb(2) - vtjy
  vnjz = vb(3) - vtjz
  ! Transfer tangential component of momentum
  tmp = vtix
  vtix = vtjx *(mb/ma)
  vtjx = tmp * (ma/mb)
  tmp = vtiy
  vtiy = vtjy *(mb/ma)
  vtjy = tmp * (ma/mb)
  tmp = vtiz
  vtiz = vtjz * (mb/ma)
  vtjz = tmp * (ma/mb)
  ! Convert back to xy frame
  va(1) = vtix + vnix 
  va(2) = vtiy + vniy
  va(3) = vtiz + vniz
  vb(1) = vtjx + vnjx
  vb(2) = vtjy + vnjy
  vb(3) = vtjz + vnjz

end subroutine collide3d








end module collision
