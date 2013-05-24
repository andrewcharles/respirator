module fpairsep
public :: compute_pairsep, pairsep_minimage, rect_minimage, gen_minimage, compute_pairsep2d

contains

subroutine compute_pairsep(dr,rij,nlist,r,n,ni,ndim)
  ! Calculates velocity difference using ilist.
  !f2py intent(in,out) :: dr
  double precision, intent(inout), dimension(ni,ndim) :: dr
  !f2py intent(in,out) :: rij
  double precision, intent(inout), dimension(ni) :: rij
  integer, intent(in), dimension(ni,2) :: nlist
  !f2py intent(in) :: nlist
  double precision, intent(in), dimension(n,ndim) :: r
  !f2py intent(in) :: r
  double precision :: drsq
  integer, intent(in) :: n,ni,ndim
  integer k,i,j,d

  do k=1,ni
    i = nlist(k,1)
    j = nlist(k,2)
    do d=1,ndim
      dr(k,d) = r(j,d) - r(i,d)
    end do
    drsq = dr(k,1)*dr(k,1) + dr(k,2)*dr(k,2) + dr(k,3)*dr(k,3) 
    rij(k) = sqrt(drsq)
  end do
end subroutine compute_pairsep

subroutine compute_pairsep2d(dr,rij,nlist,r,n,ni,ndim)
  ! Calculates velocity difference using ilist.
  !f2py intent(in,out) :: dr
  double precision, intent(inout), dimension(ni,ndim) :: dr
  !f2py intent(in,out) :: rij
  double precision, intent(inout), dimension(ni) :: rij
  integer, intent(in), dimension(ni,2) :: nlist
  !f2py intent(in) :: nlist
  double precision, intent(in), dimension(n,ndim) :: r
  !f2py intent(in) :: r
  double precision :: drsq
  integer, intent(in) :: n,ni,ndim
  integer k,i,j,d

  do k=1,ni
    i = nlist(k,1)
    j = nlist(k,2)
    do d=1,ndim
      dr(k,d) = r(j,d) - r(i,d)
    end do
    drsq = dr(k,1)*dr(k,1) + dr(k,2)*dr(k,2)
    rij(k) = sqrt(drsq)
  end do
end subroutine compute_pairsep2d

subroutine pairsep_minimage(dr,rij,nlist,lx,ly,lz,r,n,ni,ndim)
  implicit none
  ! Calculates velocity difference using ilist.
  double precision, intent(inout), dimension(ni,ndim) :: dr
  double precision, intent(inout), dimension(ni) :: rij
  double precision, intent(in) :: lx,ly,lz
  !f2py intent(in,out,overwrite) :: dr
  integer, intent(in), dimension(ni,2) :: nlist
  !f2py intent(in) :: nlist
  double precision, intent(in), dimension(n,ndim) :: r
  !f2py intent(in) :: v
  double precision :: drsq, volume
  integer, intent(in) :: n,ni,ndim
  integer k,i,j,d

  do k=1,ni
    i = nlist(k,1)
    j = nlist(k,2)
    do d=1,ndim
      dr(k,d) = r(j,d) - r(i,d)
    end do
    !rij(k) = sqrt(drsq)
  end do

  call rect_minimage(dr,lx,ly,lz,ni,ndim)

  do k=1,ni
    drsq = dr(k,1)*dr(k,1) + dr(k,2)*dr(k,2) + dr(k,3)*dr(k,3) 
    rij(k) = sqrt(drsq)
  end do

end subroutine pairsep_minimage

subroutine rect_minimage(dr,lx,ly,lz,ni,ndim)
  implicit none
  ! Calculates velocity difference using ilist.
  double precision, intent(inout), dimension(ni,ndim) :: dr
  double precision, intent(in) :: lx,ly,lz
  !f2py intent(in,out,overwrite) :: dr
  double precision :: volume
  integer, intent(in) :: ni,ndim
  integer k

  do k=1,ni
    if (dr(k,1) > lx/2.) then
        dr(k,1) = dr(k,1) - lx
    endif
    if (dr(k,2) > ly/2.) then
        dr(k,2) = dr(k,2) - ly
    endif
    if (dr(k,3) > lz/2.) then
        dr(k,3) = dr(k,3) - lz
    endif
    if (dr(k,1) < -lx/2.) then
        dr(k,1) = dr(k,1) + lx 
    endif
    if (dr(k,2) < -ly/2.) then
        dr(k,2) = dr(k,2) + ly
    endif
    if (dr(k,3) < -lz/2.) then
        dr(k,3) = dr(k,3) + lz
    endif
  end do

end subroutine rect_minimage


subroutine gen_minimage(dr,lx,ly,lz,ni,ndim)
  implicit none
  ! Calculates velocity difference using ilist.
  double precision, intent(inout), dimension(ni,ndim) :: dr
  !f2py intent(in,out,overwrite) :: dr
  double precision, intent(in) :: lx,ly,lz
  integer, intent(in) :: ni,ndim
  double precision :: volume
  integer k,i,j,d
  double precision, dimension(ni,ndim) :: shu
  double precision, dimension(ni,ndim) :: rboxv
  double precision, dimension(ndim,ndim) :: boxvec
  boxvec = 0.0
  boxvec(1,1) = lx
  boxvec(2,2) = ly
  boxvec(3,3) = lz

  volume =  boxvec(1,1)*(boxvec(2,2)*boxvec(3,3) - &
                         & boxvec(2,3)*boxvec(3,2))  &
         & -boxvec(1,2)*(boxvec(2,1)*boxvec(3,3) - &
                         & boxvec(2,3)*boxvec(3,1))  &
         & +boxvec(1,3)*(boxvec(2,1)*boxvec(3,2) - &
                         & boxvec(2,2)*boxvec(3,1))
  
    ! Evaluate position in terms of box vectors.
    rboxv(:,1) =  (dr(:,1)*boxvec(2,2)*boxvec(3,3) &
       & - dr(:,2)*boxvec(2,1)*boxvec(3,3))/volume
    rboxv(:,2) = (-dr(:,1)*boxvec(1,2)*boxvec(3,3) &
       & + dr(:,2)*boxvec(1,1)*boxvec(3,3))/volume
    rboxv(:,3) = dr(:,3)/boxvec(3,3)

  ! this integer cast may be expensive
  shu = int(2.0*rboxv - int(rboxv))

  do i = 1, 3 ! pbcs in all directions
   dr(:,i) = dr(:,i) - shu(:,1)*boxvec(1,i) - shu(:,2)*boxvec(2,i) &
                     & - shu(:,3)*boxvec(3,i)
  enddo

end subroutine gen_minimage



end module
