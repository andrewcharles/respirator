module splib
! -----------------------
! Smooth Particle Library
! -----------------------
! A library of subroutines that take only primitives as their arguments.
! These are intended to be portable implementations of algorithms that
! can be used across different SPH codes.

! Conventions:
!   - where a subroutine only works on two or three dimensional data,
! this is obvious from the title.

! F2PY has been used to expose these subroutines to python.
! Peterson, P. (2009) 'F2PY: a tool for connecting Fortran and Python
! programs', Int. J. Computational Science and Engineering.
! Vol.4, No. 4, pp.296-305.
!  http://cens.ioc.ee/~pearu/papers/IJCSE4.4_Paper_8.pdf

! Comments with a leading f2py are additional information required to
! build the python interface. These are not optional and in some cases
! the interface will not work without them.

! These functions are mostly extracted from sphforce, which will become
! a module that has only the single force loop subroutine.
! Andrew Charles, September 2009

use kernel

public :: calc_dv, calc_grad_v, calc_grad_v3d,           &
  & calc_dv_ab, calc_pi_os2d, calc_pi_os3d,              &
  & calc_pi_one, calc_viscous_entropy_liu,               &
  & calc_viscous_entropy_os_full_2d,                     &
  & calc_grad_v_os3d,calc_grad_v_os2d,                   &
  & calc_div_v2d,                                        &
  & calc_div_v3d,                                        &
  & calc_viscous_entropy_full_2d,                        &
  & calc_pv_work_full_2d,                                &
  & calc_capillary_pressure2d,                           &
  & calc_capillary_pressure3d,                           &
  & calc_heat_flux_2d,                                   &
  & calc_heat_flux_3d,                                   &
  & calc_laplacian3d,             &
  & calc_laplacian2d,             &
  & calc_laplacian3d_quot,        &
  & calc_laplacian3d_no

contains


subroutine calc_dv(dv,nlist,v,n,ni,ndim)
  ! Calculates velocity difference using ilist.
  double precision, intent(inout), dimension(ni,ndim) :: dv
  !f2py intent(in,out,overwrite) :: dv
  integer, intent(in), dimension(ni,2) :: nlist
  !f2py intent(in) :: nlist
  double precision, intent(in), dimension(n,ndim) :: v
  !f2py intent(in) :: v
  integer, intent(in) :: n,ni,ndim
  integer k,i,j,d

  do k=1,ni
    i = nlist(k,1)
    j = nlist(k,2)
    do d=1,ndim
      dv(k,d) = v(j,d) - v(i,d)
    end do
  end do
end subroutine calc_dv 


subroutine calc_dv_ab(dv,nlist,v,v2,n,ni,ndim)
    ! Calculates velocity difference using ilist for two velocity
    ! vectors (i.e. two systems).
    double precision, intent(out), dimension(ni,ndim) :: dv
    integer, intent(in), dimension(ni,2) :: nlist
    double precision, dimension(n,ndim) :: v
    double precision, dimension(n,ndim) :: v2
    integer, intent(in) :: n, ni, ndim
    integer k,i,j,d
 
    do k=1, ni
      i = nlist(k,1)
      j = nlist(k,2)
      do d=1,ndim
        dv(k,d) = v2(j,d) - v(i,d)
      enddo
    end do

end subroutine calc_dv_ab

subroutine calc_grad_v(grad_v,ilist,dwdx,dv,m,rho,n,ndim,ni)
  ! Compute the velocity gradient for 1,2 or 3 dimensions
  ! dwdx is expected to be dwdx[i] w[ij]
  ! dv is expected to be v[j] - v[i]
  implicit none
  integer, dimension(ni,2), intent(in):: ilist
  integer :: i,j,a,b,k,n,ndim,ni
  double precision, intent(inout), dimension(n,ndim,ndim) :: grad_v
  ! f2py intent(in,out) grad_v
  double precision, intent(in),  dimension(ni,ndim) :: dwdx
  double precision, intent(in),  dimension(ni,ndim) :: dv
  double precision, intent(in),  dimension(n) :: m, rho
  double precision :: mjdvk1,midvk1,mjdvk2,midvk2,mjdvk3,midvk3
  
  grad_v = 0.0
 
  do k = 1, ni
    i = ilist(k,1)
    j = ilist(k,2) 
    do a=1,ndim
      do b=1,ndim
        grad_v(i,a,b) =  grad_v(i,a,b) + m(j) * dv(k,a) * dwdx(k,b)
        grad_v(j,a,b) =  grad_v(j,a,b) + m(i) * dv(k,a) * dwdx(k,b)
       end do
    end do
  end do
  do i=1,n
    grad_v(i,:,:) =  grad_v(i,:,:)/rho(i)
  end do
end subroutine calc_grad_v

subroutine calc_grad_v3d(grad_v,ilist,dwdx,dv,m,rho,n,ndim,ni)
  ! Compute the velocity gradient for 1,2 or 3 dimensions
  ! dwdx is expected to be dwdx[i] w[ij]
  ! dv is expected to be v[j] - v[i]
  implicit none
  integer, dimension(ni,2), intent(in):: ilist
  integer :: i,j,a,b,k,n,ndim,ni
  double precision, intent(inout), dimension(n,ndim,ndim) :: grad_v
  ! f2py intent(in,out) grad_v
  double precision, intent(in),  dimension(ni,ndim) :: dwdx
  double precision, intent(in),  dimension(ni,ndim) :: dv
  double precision, intent(in),  dimension(n) :: m, rho
  double precision :: mjdvk1,midvk1,mjdvk2,midvk2,mjdvk3,midvk3
  
  grad_v = 0.0
 
  ! You could do it like the generic one above
  ! But you don't want to, because it is of the order of 10% more
  ! expensive, at least in this application with the GNU compiler
  !
  ! Instead it is done like:

  do k = 1, ni
    i = ilist(k,1)
    j = ilist(k,2)
    
    mjdvk1 = m(j) * dv(k,1)
    midvk1 = m(i) * dv(k,1)
    mjdvk2 = m(j) * dv(k,2)
    midvk2 = m(i) * dv(k,2)
    mjdvk3 = m(j) * dv(k,3)
    midvk3 = m(i) * dv(k,3)
    
    grad_v(i,1,1) =  grad_v(i,1,1) + mjdvk1  * dwdx(k,1)
    grad_v(i,1,2) =  grad_v(i,1,2) + mjdvk1  * dwdx(k,2)
    grad_v(i,1,3) =  grad_v(i,1,3) + mjdvk1  * dwdx(k,3)
    grad_v(j,1,1) =  grad_v(j,1,1) + midvk1  * dwdx(k,1)
    grad_v(j,1,2) =  grad_v(j,1,2) + midvk1  * dwdx(k,2)
    grad_v(j,1,3) =  grad_v(j,1,3) + midvk1  * dwdx(k,3)

    grad_v(i,2,1) =  grad_v(i,2,1) + mjdvk2  * dwdx(k,1)
    grad_v(i,2,2) =  grad_v(i,2,2) + mjdvk2  * dwdx(k,2)
    grad_v(i,2,3) =  grad_v(i,2,3) + mjdvk2  * dwdx(k,3)
    grad_v(j,2,1) =  grad_v(j,2,1) + midvk2  * dwdx(k,1)
    grad_v(j,2,2) =  grad_v(j,2,2) + midvk2  * dwdx(k,2)
    grad_v(j,2,3) =  grad_v(j,2,3) + midvk2  * dwdx(k,3)

    grad_v(i,3,1) =  grad_v(i,3,1) + mjdvk3  * dwdx(k,1)
    grad_v(i,3,2) =  grad_v(i,3,2) + mjdvk3  * dwdx(k,2)
    grad_v(i,3,3) =  grad_v(i,3,3) + mjdvk3  * dwdx(k,3)
    grad_v(j,3,1) =  grad_v(j,3,1) + midvk3  * dwdx(k,1)
    grad_v(j,3,2) =  grad_v(j,3,2) + midvk3  * dwdx(k,2)
    grad_v(j,3,3) =  grad_v(j,3,3) + midvk3  * dwdx(k,3)
  
  end do
  do i=1,n
    grad_v(i,:,:) =  grad_v(i,:,:)/rho(i)
  end do

end subroutine calc_grad_v3d

subroutine calc_grad_v_os2d(grad_v_os,grad_v,div_v,n,ndim)
  implicit none
  double precision, intent(inout), dimension(n,ndim,ndim) :: grad_v_os
  ! f2py intent(in,out) grad_v_os
  double precision, intent(in), dimension(n,ndim,ndim) :: grad_v
  double precision, intent(in), dimension(n) :: div_v
  integer :: i,j,k, n, ndim
 
  ! 2D case left as example
  !do i=1,n
  !  grad_v_os(i,1,1) = (0.5) * (grad_v(i,1,1)*2) - (1/ndim)*(div_v(i))   
  !  grad_v_os(i,1,2) = (0.5) * ( grad_v(i,1,2) + grad_v(i,2,1) )
  !  grad_v_os(i,2,1) = (0.5) * ( grad_v(i,2,1) + grad_v(i,1,2) )
  !  grad_v_os(i,2,2) = (0.5) * (grad_v(i,2,2)*2) - (1/ndim)*(div_v(i))
  !end do
 
  grad_v_os = 0.

  !do i=1,n
  !  do j=1,ndim
  !    do k=1,ndim
  !      grad_v_os(i,j,k) = 0.5*(grad_v(i,j,k) + grad_v(i,k,j))
  !    end do
  !    grad_v_os(i,j,j) = grad_v_os(i,j,j) - (1./ndim)*(div_v(i))
  !  end do
  !end do

  do i=1,n
    grad_v_os(i,1,1) = 0.5*(grad_v(i,1,1)+grad_v(i,1,1))-(1./2.)*(div_v(i))
    grad_v_os(i,1,2) = 0.5*(grad_v(i,1,2)+grad_v(i,2,1))
    grad_v_os(i,2,1) = 0.5*(grad_v(i,2,1)+grad_v(i,1,2))
    grad_v_os(i,2,2) = 0.5*(grad_v(i,2,2)+grad_v(i,2,2))-(1./2.)*(div_v(i))
  end do

end subroutine calc_grad_v_os2d

subroutine calc_grad_v_os3d(grad_v_os,grad_v,div_v,n,ndim)
  implicit none
  double precision, intent(inout), dimension(n,ndim,ndim) :: grad_v_os
  ! f2py intent(in,out) grad_v_os
  double precision, intent(in), dimension(n,ndim,ndim) :: grad_v
  double precision, intent(in), dimension(n) :: div_v
  integer :: i,j,k, n, ndim
 
  ! 2D case left as example
  !do i=1,n
  !  grad_v_os(i,1,1) = (1/2) * (grad_v(i,1,1)*2) - (1/ndim)*(div_v(i))   
  !  grad_v_os(i,1,2) = (1/2) * ( grad_v(i,1,2) + grad_v(i,2,1) )
  !  grad_v_os(i,2,1) = (1/2) * ( grad_v(i,2,1) + grad_v(i,1,2) )
  !  grad_v_os(i,2,2) = (1/2) * (grad_v(i,2,2)*2) - (1/ndim)*(div_v(i))
  !end do
 
  grad_v_os = 0.

  !do i=1,n
  !  do j=1,ndim
  !    do k=1,ndim
  !      grad_v_os(i,j,k) = 0.5*(grad_v(i,j,k) + grad_v(i,k,j))
  !    end do
  !    grad_v_os(i,j,j) = grad_v_os(i,j,j) - (1./ndim)*(div_v(i))
  !  end do
  !end do

  do i=1,n
    grad_v_os(i,1,1) = 0.5*(grad_v(i,1,1)+grad_v(i,1,1))-(1./3.)*(div_v(i))
    grad_v_os(i,1,2) = 0.5*(grad_v(i,1,2)+grad_v(i,2,1))
    grad_v_os(i,1,3) = 0.5*(grad_v(i,1,3)+grad_v(i,3,1))

    grad_v_os(i,2,1) = 0.5*(grad_v(i,2,1)+grad_v(i,1,2))
    grad_v_os(i,2,2) = 0.5*(grad_v(i,2,2)+grad_v(i,2,2))-(1./3.)*(div_v(i))
    grad_v_os(i,2,3) = 0.5*(grad_v(i,2,3)+grad_v(i,3,2))

    grad_v_os(i,3,1) = 0.5*(grad_v(i,3,1)+grad_v(i,1,3))
    grad_v_os(i,3,2) = 0.5*(grad_v(i,3,2)+grad_v(i,2,3))
    grad_v_os(i,3,3) = 0.5*(grad_v(i,3,3)+grad_v(i,3,3))-(1./3.)*(div_v(i))
  end do

end subroutine calc_grad_v_os3d

subroutine calc_div_v2d(grad_v,div_v,n,ndim)
  ! Computes the divergence in two dimensions
  ! This is just the trace of the velocity gradient
  double precision, intent(inout), dimension(n) :: div_v
  double precision, intent(in), dimension(n,ndim,ndim) :: grad_v
  integer :: i,j,n,ndim
 
  do i=1,n
    div_v(i) = grad_v(i,1,1) + grad_v(i,2,2)
  end do

end subroutine calc_div_v2d

subroutine calc_div_v3d(grad_v,div_v,n,ndim)
  ! Computes the divergence in one, two, three dimensions
  ! This is just the trace of the velocity gradient
  double precision, intent(inout), dimension(n) :: div_v
  double precision, intent(in), dimension(n,ndim,ndim) :: grad_v
  integer :: i,j,n,ndim
 
  ! old code
  !do i=1,n
  !  div_v(i) = 0.0
  !  do j=1,ndim
  !    div_v(i) = div_v(i) + grad_v(i,j,j)
  !  end do
  !end do

  do i=1,n
    div_v(i) = grad_v(i,1,1) + grad_v(i,2,2) + grad_v(i,3,3)
  end do

end subroutine calc_div_v3d

subroutine calc_pi_os2d(pi_os,grad_v_os,eta,n,ndim)
  ! pi_os -- returned symmetric traceless component of pressure tensor
  ! grad_v_os -- symmetric traceless velocity gradient
  ! n-- number of particles
  ! ndim -- number of dimensions
  implicit none
  double precision, intent(inout), dimension(n,ndim,ndim) :: pi_os
  ! f2py intent(in,out) pi_os
  double precision, intent(in), dimension(n,ndim,ndim) :: grad_v_os
  double precision, intent(in), dimension(n) :: eta
  integer, intent(in) :: n,ndim
  integer i,j,k

  do i=1,n
    pi_os(i,1,1) = -2.0 * eta(i) * grad_v_os(i,1,1)
    pi_os(i,1,2) = -2.0 * eta(i) * grad_v_os(i,1,2)
    pi_os(i,2,1) = -2.0 * eta(i) * grad_v_os(i,2,1)
    pi_os(i,2,2) = -2.0 * eta(i) * grad_v_os(i,2,2)
  end do

end subroutine calc_pi_os2d


subroutine calc_pi_os3d(pi_os,grad_v_os,eta,n,ndim)
  ! pi_os -- returned symmetric traceless component of pressure tensor
  ! grad_v_os -- symmetric traceless velocity gradient
  ! n-- number of particles
  ! ndim -- number of dimensions
  implicit none
  double precision, intent(inout), dimension(n,ndim,ndim) :: pi_os
  ! f2py intent(in,out) pi_os
  double precision, intent(in), dimension(n,ndim,ndim) :: grad_v_os
  double precision, intent(in), dimension(n) :: eta
  integer, intent(in) :: n,ndim
  integer i,j,k

  ! old code
  !do i=1,n
  !  do j=1,ndim
  !  do k=1,ndim
  !    pi_os(i,j,k) = -2.0 * eta(i) * grad_v_os(i,j,k)
  !  end do
  !  end do
  !end do

  do i=1,n
    pi_os(i,1,1) = -2.0 * eta(i) * grad_v_os(i,1,1)
    pi_os(i,1,2) = -2.0 * eta(i) * grad_v_os(i,1,2)
    pi_os(i,1,3) = -2.0 * eta(i) * grad_v_os(i,1,3)
    pi_os(i,2,1) = -2.0 * eta(i) * grad_v_os(i,2,1)
    pi_os(i,2,2) = -2.0 * eta(i) * grad_v_os(i,2,2)
    pi_os(i,2,3) = -2.0 * eta(i) * grad_v_os(i,2,3)
    pi_os(i,3,1) = -2.0 * eta(i) * grad_v_os(i,3,1)
    pi_os(i,3,2) = -2.0 * eta(i) * grad_v_os(i,3,2)
    pi_os(i,3,3) = -2.0 * eta(i) * grad_v_os(i,3,3)
  end do

end subroutine calc_pi_os3d


subroutine calc_pi_one(pi_one,div_v,zeta,n)
  ! Computes the contribution to the pressure tensor of the
  ! bulk viscosity. This is isotropic, and linear in the
  ! velocity gradient.
  ! n -- number of particles
  ! pi_one -- returned isotropic non-equilibrium pressre
  implicit none
  double precision, intent(inout), dimension(n) :: pi_one
  ! f2py intent(in,out) pi_one
  double precision, dimension(n) :: div_v,zeta
  integer n,i
  do i=1,n
    pi_one(i) =  - zeta(i) * div_v(i)
  end do
end subroutine calc_pi_one


subroutine calc_viscous_entropy_liu(p_ir,tdsdt,n,ndim)
  ! A simpler, cheaper way of computing tdsdt, however
  ! I have not tested this at all.
  ! TODO: Test this subroutine
  ! DO NOT USE
  integer :: n
  double precision, dimension(:) :: tdsdt
  double precision, dimension(:,:,:) :: p_ir
  integer :: ndim
  integer :: i
  do i=1,n
    if(ndim .eq. 1) then
      tdsdt(i) = p_ir(i,1,1)* p_ir(i,1,1)
    else if(ndim .eq. 2) then
      tdsdt(i) = p_ir(i,1,1) * p_ir(i,1,1) + 2.e0 * &
      &  p_ir(i,1,2)*p_ir(i,1,2) + p_ir(i,2,2)*p_ir(i,2,2)
    else if(ndim .eq. 3) then
      stop '3D not implemented in liu viscous entropy'
    endif
 end do
end subroutine calc_viscous_entropy_liu


subroutine calc_viscous_entropy_os_full_2d(ilist,p_ir,tdsdt,grad_v_os,rho,n,ni)
  ! The full computation - this may be unneccessarily complex
  integer, dimension(ni,2), intent(in):: ilist
  double precision, dimension(n,2,2) :: p_ir
  double precision, intent(inout), dimension(n) :: tdsdt
  double precision, dimension(n,2,2) :: grad_v_os
  double precision, dimension(n) :: rho
  integer :: n, ni
  integer :: i,j,k 
  do k = 1, ni
    i = ilist(k,1)
    j = ilist(k,2) 

    tdsdt(i) = tdsdt(i) + ( p_ir(i,1,1)/rho(i)**2 + p_ir(j,1,1)/rho(j)**2 ) &
      & * grad_v_os(i,1,1)
    tdsdt(i) = tdsdt(i) + ( p_ir(i,1,2)/rho(i)**2 + p_ir(j,1,2)/rho(j)**2 ) &
      & *grad_v_os(i,2,1)
    tdsdt(i) = tdsdt(i) + ( p_ir(i,2,1)/rho(i)**2 + p_ir(j,2,1)/rho(j)**2 ) &
      & *grad_v_os(i,1,2)
    tdsdt(i) = tdsdt(i) + ( p_ir(i,2,2)/rho(i)**2 + p_ir(j,2,2)/rho(j)**2 ) &
      &*grad_v_os(i,2,2)

    tdsdt(j) = tdsdt(j) + ( p_ir(j,1,1)/rho(j)**2 + p_ir(i,1,1)/rho(i)**2 ) &
      & *grad_v_os(j,1,1)
    tdsdt(j) = tdsdt(j) + ( p_ir(j,1,2)/rho(j)**2 + p_ir(i,1,2)/rho(i)**2 ) &
      & *grad_v_os(j,2,1)
    tdsdt(j) = tdsdt(j) + ( p_ir(j,1,2)/rho(j)**2 + p_ir(i,1,2)/rho(i)**2 ) &
      & *grad_v_os(j,1,2)
    tdsdt(j) = tdsdt(j) + ( p_ir(j,2,2)/rho(j)**2 + p_ir(i,2,2)/rho(i)**2 ) &
      & *grad_v_os(j,2,2)

  end do

  do i = 1,n
    tdsdt(i) = 0.5 * tdsdt(i)
  end do

end subroutine calc_viscous_entropy_os_full_2d


subroutine calc_viscous_entropy_full_2d(ilist,p_ir,tdsdt,grad_v,rho,n,ni)
  integer, dimension(ni,2), intent(in):: ilist
  double precision, dimension(n,2,2) :: p_ir
  double precision, intent(inout), dimension(n) :: tdsdt
  double precision, dimension(n,2,2) :: grad_v
  double precision, dimension(n) :: rho
  integer :: n, ni
  integer :: i,j,k 
  do k = 1, ni
    i = ilist(k,1)
    j = ilist(k,2) 

    tdsdt(i) = tdsdt(i) + ( p_ir(i,1,1)/rho(i)**2 + p_ir(j,1,1)/rho(j)**2 ) &
      & * grad_v(i,1,1)
    tdsdt(i) = tdsdt(i) + ( p_ir(i,1,2)/rho(i)**2 + p_ir(j,1,2)/rho(j)**2 ) &
      & *grad_v(i,2,1)
    tdsdt(i) = tdsdt(i) + ( p_ir(i,2,1)/rho(i)**2 + p_ir(j,2,1)/rho(j)**2 ) &
      & *grad_v(i,1,2)
    tdsdt(i) = tdsdt(i) + ( p_ir(i,2,2)/rho(i)**2 + p_ir(j,2,2)/rho(j)**2 ) &
      &*grad_v(i,2,2)

    tdsdt(j) = tdsdt(j) + ( p_ir(j,1,1)/rho(j)**2 + p_ir(i,1,1)/rho(i)**2 ) &
      & *grad_v(j,1,1)
    tdsdt(j) = tdsdt(j) + ( p_ir(j,1,2)/rho(j)**2 + p_ir(i,1,2)/rho(i)**2 ) &
      & *grad_v(j,2,1)
    tdsdt(j) = tdsdt(j) + ( p_ir(j,1,2)/rho(j)**2 + p_ir(i,1,2)/rho(i)**2 ) &
      & *grad_v(j,1,2)
    tdsdt(j) = tdsdt(j) + ( p_ir(j,2,2)/rho(j)**2 + p_ir(i,2,2)/rho(i)**2 ) &
      & *grad_v(j,2,2)

  end do

  do i = 1,n
    tdsdt(i) = 0.5 * tdsdt(i)
  end do

end subroutine calc_viscous_entropy_full_2d

subroutine calc_pv_work_full_2d(ilist,p,de,grad_v,rho,n,ni)
  ! This is not currently used
  integer, dimension(ni,2), intent(in):: ilist
  double precision, dimension(n,2,2) :: p
  double precision, intent(inout), dimension(n) :: de
  double precision, dimension(n,2,2) :: grad_v
  double precision, dimension(n) :: rho
  integer :: n, ni
  integer :: i,j,k 
  do k = 1, ni
    i = ilist(k,1)
    j = ilist(k,2) 

    de(i) = de(i) + ( p(i,1,1)/rho(i)**2 + p(j,1,1)/rho(j)**2 ) &
      & * grad_v(i,1,1)
    de(i) = de(i) + ( p(i,1,2)/rho(i)**2 + p(j,1,2)/rho(j)**2 ) &
      & * grad_v(i,2,1)
    de(i) = de(i) + ( p(i,2,1)/rho(i)**2 + p(j,2,1)/rho(j)**2 ) &
      & * grad_v(i,1,2)
    de(i) = de(i) + ( p(i,2,2)/rho(i)**2 + p(j,2,2)/rho(j)**2 ) &
      & * grad_v(i,2,2)

    de(j) = de(j) + ( p(j,1,1)/rho(j)**2 + p(i,1,1)/rho(i)**2 ) &
      & * grad_v(j,1,1)
    de(j) = de(j) + ( p(j,1,2)/rho(j)**2 + p(i,1,2)/rho(i)**2 ) &
      & * grad_v(j,2,1)
    de(j) = de(j) + ( p(j,1,2)/rho(j)**2 + p(i,1,2)/rho(i)**2 ) &
      & * grad_v(j,1,2)
    de(j) = de(j) + ( p(j,2,2)/rho(j)**2 + p(i,2,2)/rho(i)**2 ) &
      & * grad_v(j,2,2)

  end do

  do i = 1,n
    de(i) = 0.5 * de(i)
  end do

end subroutine calc_pv_work_full_2d

subroutine calc_heat_flux_2d(ilist,q,rho,m,tmp,dwdx_jij,thermalk,n,ni,ndim)
  ! Computes the heat flux for a set of SPAM particles.
  integer, dimension(ni,ndim), intent(in):: ilist
  double precision, dimension(n,ndim), intent(inout) :: q
  ! f2py intent(in,out) q
  double precision, dimension(n) :: rho,m,tmp
  double precision, dimension(ni,ndim) :: dwdx_jij
  double precision :: thermalk
  integer n,ni,ndim
  double precision :: rhoij,mij,tij
  integer i,j,k

  q=0
      
  do k = 1, ni
    i = ilist(k,1)
    j = ilist(k,2)
    rhoij = (rho(i)+rho(j))/2.
    mij = (m(i)+m(j))/2. 
    tij = tmp(j) - tmp(i) 
    
    q(i,1) = q(i,1) - thermalk * (mij/rhoij) * tij * -dwdx_jij(k,1)
    q(i,2) = q(i,2) - thermalk * (mij/rhoij) * tij * -dwdx_jij(k,2)
   
    q(j,1) = q(j,1) - thermalk * (mij/rhoij) * -tij * dwdx_jij(k,1)
    q(j,2) = q(j,2) - thermalk * (mij/rhoij) * -tij * dwdx_jij(k,2) 
  end do

end subroutine calc_heat_flux_2d

subroutine calc_heat_flux_3d(q,ilist,rho,m,tmp,dwdx_jij,thermalk,n,ni)
  ! Computes the heat flux for a set of SPAM particles.
  integer, dimension(ni,2), intent(in):: ilist
  double precision, dimension(n,3), intent(inout) :: q
  ! f2py intent(in,out) q
  double precision, dimension(n) :: rho,m,tmp
  double precision, dimension(ni,3) :: dwdx_jij
  double precision :: thermalk
  double precision :: rhoij,mij,tij,mt_on_rho
  integer i,j,k,n,ni

  q=0
      
  do k = 1, ni
    i = ilist(k,1)
    j = ilist(k,2)
    rhoij = (rho(i)+rho(j))/2.
    mij = (m(i)+m(j))/2. 
    tij = tmp(j) - tmp(i) 
    mt_on_rho = thermalk * (mij/rhoij) * tij

    ! this sign may actually be wrong for jij and assuming iij
    q(i,1) = q(i,1) + mt_on_rho * dwdx_jij(k,1)
    q(i,2) = q(i,2) + mt_on_rho * dwdx_jij(k,2)
    q(i,3) = q(i,3) + mt_on_rho * dwdx_jij(k,3)
   
    q(j,1) = q(j,1) + mt_on_rho * dwdx_jij(k,1)
    q(j,2) = q(j,2) + mt_on_rho * dwdx_jij(k,2) 
    q(j,3) = q(j,3) + mt_on_rho * dwdx_jij(k,3) 
  end do

end subroutine calc_heat_flux_3d

subroutine calc_laplacian2d(ilist,rho,m,grad_rho,dwdx,gdotgrho,n,ni)
  ! ilist -- indices of interacting pairs
  ! ni -- number of interacting pairs
  ! rho -- particle densities
  ! m -- particle masses
  ! n -- number of particles (size of rho, m, first dimension of grad_rho)
  ! grad_rho -- density gradient of particles
  ! dwdx -- kernel gradient of particles
  ! TODO - clarify whether this is grad_j_dwij or grad_i_dwij
  ! Have to assume this is grad_j_dwij (i.e the gradient at the
  ! position of particle j due to particle i)?
  ! this doesn't actually look right...

  ! Compute the capillary part of the pressure tensor.
  integer, dimension(ni,2), intent(in):: ilist
  double precision, dimension(n), intent(in) :: rho
  double precision, dimension(n), intent(in) :: m
  double precision, dimension(n,2), intent(in) :: grad_rho
  double precision, dimension(ni,2), intent(in) :: dwdx
  ! f2py intent(in,out) gdotgrho
  double precision, intent(inout), dimension(n) :: gdotgrho
  integer, intent(in) :: n, ni
 
  ! Internal variables
  integer i,j,k
  double precision :: rhosq_i,rhosq_j,ggx,ggy,iso_term,rhoimj,rhojmi,rhoij,mij
  
  gdotgrho=0
  
  ! Calculate grad grad rho   
  ! -----------------------
  do k=1, ni
    i = ilist(k,1)
    j = ilist(k,2) 
    rhosq_i = rho(i)**2
    rhosq_j = rho(j)**2
    
    ! PRODUCT RULE SYMMETRISED
    rhoij = (rho(i)+rho(j))/2.
    mij = ((m(i)+m(j))/2.)

    ggx = (mij/rhoij) * (grad_rho(j,1) - grad_rho(i,1)) * dwdx(k,1)
    ggy = (mij/rhoij) * (grad_rho(j,2) - grad_rho(i,2)) * dwdx(k,2)

    ! These are the same sign because gri - grj reverses sign
    ! when indices are reversed, but so does dwdxij
    gdotgrho(i) = gdotgrho(i) - ggx
    gdotgrho(i) = gdotgrho(i) - ggy
    gdotgrho(j) = gdotgrho(j) - ggx 
    gdotgrho(j) = gdotgrho(j) - ggy

  end do

end subroutine calc_laplacian2d

subroutine calc_capillary_pressure2d(ilist,rho,m,grad_rho,dwdx,p_rev,cgrad,n,ni,ndim)
  ! ilist -- indices of interacting pairs
  ! ni -- number of interacting pairs
  ! rho -- particle densities
  ! m -- particle masses
  ! n -- number of particles (size of rho, m, first dimension of grad_rho)
  ! grad_rho -- density gradient of particles
  ! dwdx -- kernel gradient of particles
  ! TODO - clarify whether this is grad_j_dwij or grad_i_dwij
  ! p_rev -- the capillary pressure is accumulated into this array
  ! cgrad -- gradient coefficient (constant)

  ! Compute the capillary part of the pressure tensor.
  integer, dimension(ni,ndim), intent(in):: ilist
  double precision, dimension(n), intent(in) :: rho
  double precision, dimension(n), intent(in) :: m
  double precision, dimension(n,ndim), intent(in) :: grad_rho
  double precision, dimension(ni,ndim), intent(in) :: dwdx
  ! f2py intent(in,out) p_rev
  double precision, dimension(n,ndim,ndim), intent(inout) :: p_rev
  double precision, intent(in) :: cgrad
  integer, intent(in) :: n, ni, ndim
 
  ! Internal variables
  integer i,j,k
  double precision, dimension(n) :: gradgradrho
  double precision, dimension(n) :: maggradrhosq
  double precision :: rhosq_i,rhosq_j,ggx,ggy,rhoimj,rhojmi
  
  gradgradrho=0
  
  ! Calculate grad grad rho   
  ! -----------------------
  do k=1, ni
    i = ilist(k,1)
    j = ilist(k,2) 

    rhosq_i = rho(i)**2
    rhosq_j = rho(j)**2
    rhoimj = rho(i) * m(j)
    rhojmi = rho(j) * m(i)

    ggx = ( grad_rho(i,1)/(rhosq_i) + grad_rho(j,1)/(rhosq_j) )*dwdx(k,1)
    ggy = ( grad_rho(i,2)/(rhosq_i) + grad_rho(j,2)/(rhosq_j) )*dwdx(k,2)

    gradgradrho(i) = gradgradrho(i) - rhoimj * ggx
    gradgradrho(i) = gradgradrho(i) - rhoimj * ggy
    gradgradrho(j) = gradgradrho(j) + rhojmi * ggx
    gradgradrho(j) = gradgradrho(j) + rhojmi * ggy
 
!   gradgradrho(i) = gradgradrho(i) + rho(i)*m(j)*( grad_rho(i,1)/(rho(i)**2) &
!      &            + grad_rho(j,1)/(rho(j)**2) )*dwdx(k,1)
!   gradgradrho(i) = gradgradrho(i) + rho(i)*m(j)*( grad_rho(i,2)/(rho(i)**2) &
!!      &            + grad_rho(j,2)/(rho(j)**2) )*dwdx(k,2)
!   gradgradrho(j) = gradgradrho(j) - rho(j)*m(i)*( grad_rho(i,1)/(rho(i)**2) &
!      &            + grad_rho(j,1)/(rho(j)**2) )*dwdx(k,1)
!   gradgradrho(j) = gradgradrho(j) - rho(j)*m(i)*( grad_rho(i,2)/(rho(i)**2) &
!      &            + grad_rho(j,2)/(rho(j)**2) )*dwdx(k,2)

  end do
  
  do i=1,n
    maggradrhosq = abs(grad_rho(i,1)*grad_rho(i,1) &
      &          + grad_rho(i,2)*grad_rho(i,2))
    p_rev(i,1,1) = p_rev(i,1,1) - cgrad*rho(i)*gradgradrho(i) &
      &          - cgrad*(0.5)*maggradrhosq(i)
    p_rev(i,2,2) = p_rev(i,2,2) - cgrad*rho(i)*gradgradrho(i) &
      &          - cgrad*(0.5)*maggradrhosq(i)
  
    p_rev(i,1,1) = p_rev(i,1,1) + cgrad*grad_rho(i,1)**2
    p_rev(i,1,2) = p_rev(i,1,2) + cgrad*grad_rho(i,1)*grad_rho(i,2)
    p_rev(i,2,1) = p_rev(i,2,1) + cgrad*grad_rho(i,1)*grad_rho(i,2)
    p_rev(i,2,2) = p_rev(i,2,2) + cgrad*grad_rho(i,2)**2
  enddo

end subroutine calc_capillary_pressure2d




subroutine calc_laplacian3d(ilist,rho,m,grad_rho,dwdx,gdotgrho,n,ni)
  ! ilist -- indices of interacting pairs
  ! ni -- number of interacting pairs
  ! rho -- particle densities
  ! m -- particle masses
  ! n -- number of particles (size of rho, m, first dimension of grad_rho)
  ! grad_rho -- density gradient of particles
  ! dwdx -- kernel gradient of particles
  ! TODO - clarify whether this is grad_j_dwij or grad_i_dwij
  ! Have to assume this is grad_j_dwij (i.e the gradient at the
  ! position of particle j due to particle i)?
  ! this doesn't actually look right...

  ! Compute the capillary part of the pressure tensor.
  integer, dimension(ni,2), intent(in):: ilist
  double precision, dimension(n), intent(in) :: rho
  double precision, dimension(n), intent(in) :: m
  double precision, dimension(n,3), intent(in) :: grad_rho
  double precision, dimension(ni,3), intent(in) :: dwdx
  ! f2py intent(in,out) gdotgrho
  double precision, intent(inout), dimension(n) :: gdotgrho
  integer, intent(in) :: n, ni
 
  ! Internal variables
  integer i,j,k
  double precision :: rhosq_i,rhosq_j,ggx,ggy,ggz,iso_term,rhoimj,rhojmi,rhoij,mij
  
  gdotgrho=0
  
  ! Calculate grad grad rho   
  ! -----------------------
  do k=1, ni
    i = ilist(k,1)
    j = ilist(k,2) 
    rhosq_i = rho(i)**2
    rhosq_j = rho(j)**2
    
    ! QUOTIENT RULE SYMMETRISED
    !rhoimj = rho(i) * m(j)
    !rhojmi = rho(j) * m(i)
    !ggx = ( grad_rho(i,1)/(rhosq_i) + grad_rho(j,1)/(rhosq_j) )*dwdx(k,1)
    !ggy = ( grad_rho(i,2)/(rhosq_i) + grad_rho(j,2)/(rhosq_j) )*dwdx(k,2)
    !ggz = ( grad_rho(i,3)/(rhosq_i) + grad_rho(j,3)/(rhosq_j) )*dwdx(k,3)
    !gdotgrho(i) = gdotgrho(i) - m(j) * ggx!-
    !gdotgrho(i) = gdotgrho(i) - m(j) * ggy!-
    !gdotgrho(i) = gdotgrho(i) - m(j) * ggz !-
    !gdotgrho(j) = gdotgrho(j) + m(i) * ggx 
    !gdotgrho(j) = gdotgrho(j) + m(i) * ggy
    !gdotgrho(j) = gdotgrho(j) + m(i) * ggz

    ! PRODUCT RULE SYMMETRISED
    rhoij = (rho(i)+rho(j))/2.
    mij = ((m(i)+m(j))/2.)

    ggx = (mij/rhoij) * (grad_rho(j,1) - grad_rho(i,1)) * dwdx(k,1)
    ggy = (mij/rhoij) * (grad_rho(j,2) - grad_rho(i,2)) * dwdx(k,2)
    ggz = (mij/rhoij) * (grad_rho(j,3) - grad_rho(i,3)) * dwdx(k,3)

    ! These are the same sign because gri - grj reverses sign
    ! when indices are reversed, but so does dwdxij
    gdotgrho(i) = gdotgrho(i) - ggx!-
    gdotgrho(i) = gdotgrho(i) - ggy!-
    gdotgrho(i) = gdotgrho(i) - ggz !-
    gdotgrho(j) = gdotgrho(j) - ggx 
    gdotgrho(j) = gdotgrho(j) - ggy
    gdotgrho(j) = gdotgrho(j) - ggz

    ! UNSYMMETRISED
  end do

end subroutine calc_laplacian3d

subroutine calc_laplacian3d_quot(ilist,rho,m,grad_rho,dwdx,gdotgrho,n,ni)
  ! ilist -- indices of interacting pairs
  ! ni -- number of interacting pairs
  ! rho -- particle densities
  ! m -- particle masses
  ! n -- number of particles (size of rho, m, first dimension of grad_rho)
  ! grad_rho -- density gradient of particles
  ! dwdx -- kernel gradient of particles
  ! TODO - clarify whether this is grad_j_dwij or grad_i_dwij
  ! Have to assume this is grad_j_dwij (i.e the gradient at the
  ! position of particle j due to particle i)?
  ! this doesn't actually look right...

  ! Compute the capillary part of the pressure tensor.
  integer, dimension(ni,2), intent(in):: ilist
  double precision, dimension(n), intent(in) :: rho
  double precision, dimension(n), intent(in) :: m
  double precision, dimension(n,3), intent(in) :: grad_rho
  double precision, dimension(ni,3), intent(in) :: dwdx
  ! f2py intent(in,out) gdotgrho
  double precision, intent(inout), dimension(n) :: gdotgrho
  integer, intent(in) :: n, ni
 
  ! Internal variables
  integer i,j,k
  double precision :: rhosq_i,rhosq_j,ggx,ggy,ggz,iso_term,rhoimj,rhojmi,rhoij,mij
  
  gdotgrho=0
  
  ! Calculate grad grad rho   
  ! -----------------------
  do k=1, ni
    i = ilist(k,1)
    j = ilist(k,2) 
    rhosq_i = rho(i)**2
    rhosq_j = rho(j)**2
    
    ! QUOTIENT RULE SYMMETRISED
    rhoimj = rho(i) * m(j)
    rhojmi = rho(j) * m(i)
    ggx = ( grad_rho(i,1)/(rhosq_i) + grad_rho(j,1)/(rhosq_j) )*dwdx(k,1)
    ggy = ( grad_rho(i,2)/(rhosq_i) + grad_rho(j,2)/(rhosq_j) )*dwdx(k,2)
    ggz = ( grad_rho(i,3)/(rhosq_i) + grad_rho(j,3)/(rhosq_j) )*dwdx(k,3)
    gdotgrho(i) = gdotgrho(i) - rhoimj * ggx
    gdotgrho(i) = gdotgrho(i) - rhoimj * ggy
    gdotgrho(i) = gdotgrho(i) - rhoimj * ggz
    gdotgrho(j) = gdotgrho(j) + rhojmi * ggx 
    gdotgrho(j) = gdotgrho(j) + rhojmi * ggy
    gdotgrho(j) = gdotgrho(j) + rhojmi * ggz

  end do

end subroutine calc_laplacian3d_quot

subroutine calc_laplacian3d_no(ilist,rho,m,grad_rho,dwdx,gdotgrho,n,ni)
  ! ilist -- indices of interacting pairs
  ! ni -- number of interacting pairs
  ! rho -- particle densities
  ! m -- particle masses
  ! n -- number of particles (size of rho, m, first dimension of grad_rho)
  ! grad_rho -- density gradient of particles
  ! dwdx -- kernel gradient of particles
  ! no symmetrisation

  ! Compute the capillary part of the pressure tensor.
  integer, dimension(ni,2), intent(in):: ilist
  double precision, dimension(n), intent(in) :: rho
  double precision, dimension(n), intent(in) :: m
  double precision, dimension(n,3), intent(in) :: grad_rho
  double precision, dimension(ni,3), intent(in) :: dwdx
  ! f2py intent(in,out) gdotgrho
  double precision, intent(inout), dimension(n) :: gdotgrho
  integer, intent(in) :: n, ni
 
  ! Internal variables
  integer i,j,k
  double precision :: rhosq_i,rhosq_j,ggx,ggy,ggz,iso_term,rhoimj,rhojmi,rhoij,mij,factor
  
  gdotgrho=0
  
  ! Calculate grad grad rho   
  ! -----------------------
  do k=1, ni
    i = ilist(k,1)
    j = ilist(k,2) 
    
    ! No SYMMETRISED
    factor = m(j)/rho(j)
    gdotgrho(i) = gdotgrho(i) - factor * grad_rho(j,1) * dwdx(k,1)
    gdotgrho(i) = gdotgrho(i) - factor * grad_rho(j,2) * dwdx(k,2)
    gdotgrho(i) = gdotgrho(i) - factor * grad_rho(j,3) * dwdx(k,3)
    gdotgrho(j) = gdotgrho(j) + factor * grad_rho(i,1) * dwdx(k,1)
    gdotgrho(j) = gdotgrho(j) + factor * grad_rho(i,2) * dwdx(k,2)
    gdotgrho(j) = gdotgrho(j) + factor * grad_rho(i,3) * dwdx(k,3)

  end do

end subroutine calc_laplacian3d_no


subroutine calc_capillary_pressure3d(ilist,rho,m,grad_rho,dwdx,p_rev,cgrad,n,ni)
  ! ilist -- indices of interacting pairs
  ! ni -- number of interacting pairs
  ! rho -- particle densities
  ! m -- particle masses
  ! n -- number of particles (size of rho, m, first dimension of grad_rho)
  ! grad_rho -- density gradient of particles
  ! dwdx -- kernel gradient of particles
  ! This is grad_j_dwij
  ! p_rev -- the capillary pressure is accumulated into this array
  ! cgrad -- gradient coefficient (constant)

  ! Compute the capillary part of the pressure tensor.
  integer, dimension(ni,2), intent(in):: ilist
  double precision, dimension(n), intent(in) :: rho
  double precision, dimension(n), intent(in) :: m
  double precision, dimension(n,3), intent(in) :: grad_rho
  double precision, dimension(ni,3), intent(in) :: dwdx
  double precision, dimension(n,3,3), intent(inout) :: p_rev
  ! f2py intent(in,out) p_rev
  double precision, intent(in) :: cgrad
  integer, intent(in) :: n, ni
 
  ! Internal variables
  integer i,j,k
  double precision, dimension(n) :: gradgradrho
  double precision :: rhosq_i,rhosq_j,ggx,ggy,ggz,iso_term,rhoimj,rhojmi
  
  gradgradrho=0
  
  ! Calculate grad grad rho   
  ! -----------------------
  do k=1, ni
    i = ilist(k,1)
    j = ilist(k,2) 
    rhosq_i = rho(i)**2
    rhosq_j = rho(j)**2
    rhoimj = rho(i) * m(j)
    rhojmi = rho(j) * m(i)

    ggx = ( grad_rho(i,1)/(rhosq_i) + grad_rho(j,1)/(rhosq_j) )*dwdx(k,1)
    ggy = ( grad_rho(i,2)/(rhosq_i) + grad_rho(j,2)/(rhosq_j) )*dwdx(k,2)
    ggz = ( grad_rho(i,3)/(rhosq_i) + grad_rho(j,3)/(rhosq_j) )*dwdx(k,3)

    gradgradrho(i) = gradgradrho(i) - rhoimj * ggx
    gradgradrho(i) = gradgradrho(i) - rhoimj * ggy
    gradgradrho(i) = gradgradrho(i) - rhoimj * ggz
    gradgradrho(j) = gradgradrho(j) + rhojmi * ggx
    gradgradrho(j) = gradgradrho(j) + rhojmi * ggy
    gradgradrho(j) = gradgradrho(j) + rhojmi * ggz

  end do
  
  do i=1,n
    ! this quantity is 
    ! c * rho(i) * gradgradrho(i)
    ! magnitude(gradrho squared) * c * 0.5 
    iso_term = cgrad    *                                 &
      &      ( rho(i) * gradgradrho(i) +                  &
      &        0.5 * abs( grad_rho(i,1) * grad_rho(i,1)   &
      &                 + grad_rho(i,2) * grad_rho(i,2)   &
      &                 + grad_rho(i,3) * grad_rho(i,3)))

    p_rev(i,1,1) = p_rev(i,1,1) - iso_term
    p_rev(i,2,2) = p_rev(i,2,2) - iso_term
    p_rev(i,3,3) = p_rev(i,3,3) - iso_term

    p_rev(i,1,1) = p_rev(i,1,1) + cgrad*(grad_rho(i,1)**2)
    p_rev(i,1,2) = p_rev(i,1,2) + cgrad*grad_rho(i,1)*grad_rho(i,2)
    p_rev(i,1,3) = p_rev(i,1,3) + cgrad*grad_rho(i,1)*grad_rho(i,3)

    p_rev(i,2,1) = p_rev(i,2,1) + cgrad*grad_rho(i,2)*grad_rho(i,1)
    p_rev(i,2,2) = p_rev(i,2,2) + cgrad*(grad_rho(i,2)**2)
    p_rev(i,2,3) = p_rev(i,2,3) + cgrad*grad_rho(i,2)*grad_rho(i,3)

    p_rev(i,3,1) = p_rev(i,3,1) + cgrad*grad_rho(i,3)*grad_rho(i,1)
    p_rev(i,3,2) = p_rev(i,3,2) + cgrad*grad_rho(i,3)*grad_rho(i,2)
    p_rev(i,3,3) = p_rev(i,3,3) + cgrad*(grad_rho(i,3)**2)

  enddo

end subroutine calc_capillary_pressure3d


end module splib
