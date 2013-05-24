! Subroutines for calculating the SPH kernel and kernel gradient
! for pairs of particles, and for all pairs in the system.
! 
! For compatibility with python bindings, all kernels must
! have the signature:
!   subroutine kernel(r,dx,h,w,dwdx,ndim)
!
! The kernel gradient returned here is the gradient at distance r,
! and displacement dx from the kernel center.
! i.e. if the kernel gradient contribution of particle i on particle j
! is desired, the call must be lucy_kernel(r,(rj - ri)...) where rj-ri
! is the vector from rj to ri

module kernel
!use global

implicit none
public :: calc_kernel, lucy_kernel, gauss_kernel, debrun_kernel, &
&   smoothing_kernels,adapt_smoothing_length, density_sum,       &
& adapt_kernel_width,add_boundary_density

!double precision :: pi = 4.0 * atan(1.0)
double precision :: pi = 3.14159265358979323846

contains

subroutine lucy_kernel(r,dx,h,w,dwdx,ndim)
  ! q -- normalisation factor
  ! r --
  ! dx -- distance to kernel center
  implicit none
  double precision, intent(in) :: r
  double precision, intent(in), dimension(ndim) :: dx
  double precision, intent(in) :: h
  double precision, intent(inout) :: w
  double precision, intent(inout), dimension(ndim) :: dwdx
  double precision :: q
  integer :: d, ndim

  if(ndim .eq. 1) then
    q = 5./(4. * h)
  else if(ndim .eq. 2) then
    q = 5./(pi * (h**2))
  else if(ndim .eq. 3) then
    q = 105./(16. * pi * (h**3))
  endif

  if( r < h ) then
    w = q * ( (-3.*(r**4)/(h**4)) + (8.*(r**3)/(h**3)) &
    & - (6.*(r**2)/(h**2)) + 1. )
    
    if(r .eq. 0) then
      dwdx=0.0
    else
      do d=1,ndim
        dwdx(d) = q * ( (-12./(h**4))*(r**3) &
        &       + (24./(h**3))*(r**2) - (12.*r/(h**2)) )* dx(d)/r
      end do
    end if
  else
    do d=1,ndim
      dwdx(d) = 0
    end do
      w = 0
  end if
      
end subroutine lucy_kernel


subroutine gauss_kernel(r,dx,h,w,dwdx,ndim)
  ! A Gaussian smooth particle kernel function
  ! This one should work for 1, 2 or 3 dimensions.
  implicit none
  double precision, intent(in) :: r
  double precision, intent(in), dimension(ndim) :: dx
  double precision, intent(in) :: h
  double precision, intent(inout), dimension(ndim) :: dwdx
  double precision, intent(inout) :: w
  integer d,ndim
  double precision factor, q
  q = r/h
  factor = 1.e0 / ((h*sqrt(pi))**ndim) 
  w = factor * exp(-q*q)
  do d=1,ndim
    dwdx(d) = w * ( (-2.e0) * dx(d)/(h*h))
  end do
end subroutine gauss_kernel


subroutine debrun_kernel(r,dx,h,w,dwdx,ndim)
  ! Standing on the shoulders of
  ! Muller et al
  ! Particle-Based Fluid Simulation for Interactive Applications
  ! Values for 1D and 2D were computed in Maple - not sure how
  ! accurate they are. I'd like to revisit the calculation in Sage
  ! one of these days but you know how it is with kids.
  ! w = ( (5/(pi * h**2)) *(1 + 3*r/h)*((1-r/h) )**3)
  ! dw =  q * ((h**3) - 3*(h*h)*r + 3*(r*r)*h - (r**3))
  implicit none
  double precision, intent(in) :: r
  double precision, intent(in), dimension(ndim) :: dx
  double precision, intent(in) :: h
  double precision, intent(inout), dimension(ndim) :: dwdx
  double precision, intent(inout) :: w
  integer i,ndim
  double precision q
  if(ndim .eq. 2) then
    q = (10.0/(pi*(h**5)))
  else if(ndim .eq. 3) then
    q = (15.0/(pi*(h**6)))
  else if(ndim .eq. 1) then
    q = ( ( (-7.*(h**4)/4) + 2*(h**4) )**(-1) )/ 2 
  endif

  if( r < h ) then
    w =  q * ((h - r)**3)
    if(r .eq. 0) then
      dwdx = 0.0
    else
      do i=1,ndim
        dwdx(i) = q * (-3.0*(h**2) + 6.0*h*r - 3.0*(r**2)) * (dx(i)/r)
      end do
    end if
  else
    do i=1,ndim
      dwdx(i) = 0
    end do
    w = 0
  end if
end subroutine debrun_kernel


subroutine debrun2d(r,dx,h,w,dwdx) 
  ! w = ( (5/(pi * h**2)) *(1 + 3*r/h)*((1-r/h) )**3)
  ! dw =  q * ((h**3) - 3*(h*h)*r + 3*(r*r)*h - (r**3))
  implicit none
  double precision, intent(in) :: r
  double precision, intent(in), dimension(2) :: dx
  double precision, intent(in) :: h
  double precision, intent(inout), dimension(2) :: dwdx
  double precision, intent(inout) :: w
  integer i
  double precision q
  if( r < h ) then
    q = (20.0/(pi*(h**6)))
    w =  q * ((h - r)**3)
    if(r .eq. 0) then
      dwdx = 0.0
    else
      do i=1,2
        dwdx(i) = q * (-3.0*(h**2) + 6.0*h*r - 3.0*(r**2)) * (dx(i)/r)
      end do
    end if
  else
    do i=1,2
      dwdx(i) = 0
    end do
    w = 0
  end if
end subroutine debrun2d


subroutine calc_kernel(r,dx,h,w,dwdx,kernel_type,ndim)
  ! Calculates the kernel and kernel derivative for a displacement
  ! The kernel gradient is calculated at displacment dx. If using this
  ! function to compute the kernel gradient with respect to two particles
  ! this should be the displacement from the center of the kernel to
  ! the position you want the gradient at. By convention we tend to store
  ! grad_i W_ij, so dx in this case is r_i - r_j
  !
  ! r -- distance 
  ! dx -- x,y,z displacement
  ! dwdx -- first spatial derivative of kernel
  ! w -- kernel value
  ! h -- smoothing lengths of each particle
  ! kernel_type -- selects the form of kernel
  ! ndim -- number of dimensions
  
  implicit none
  double precision, intent(in) :: r
  double precision, intent(in), dimension(ndim) :: dx
  double precision, intent(in) :: h
  double precision, intent(inout), dimension(ndim) :: dwdx
  double precision, intent(inout) :: w
  integer kernel_type, ndim

  ! Gaussian Kernel
  if(kernel_type .eq. 1) then
    call gauss_kernel(r,dx,h,w,dwdx,ndim)
  ! Lucy Kernel
  else if(kernel_type .eq. 2 .or. kernel_type .eq. 4) then
    call lucy_kernel(r,dx,h,w,dwdx,ndim)
  else if(kernel_type .eq. 3) then  
    ! Debrun spiky
    ! gamma_g := 20/(Pi*h^6);
    ! w_spiky_2d := (r,h) -> piecewise((r<(h)), gamma_g * (h-r)^3, 0);
    call debrun_kernel(r,dx,h,w,dwdx,ndim)
  endif

end subroutine calc_kernel


subroutine smoothing_kernels(w,dwdx,rij,drij,nlist,sml,kernel_type,n,ni,ndim)
  ! Calculates kernel and kernel gradients for a given set of
  ! particle displacements and smoothing lengths
  ! w -- output kernel values for each pair
  ! dwdx -- output kernel gradient for each pair
  ! rij -- pair distances
  ! drij -- pair displacements
  ! nlist - integer of pair indices
  ! h -- smoothing length
  ! ni -- number of pairs
  ! ndim -- dimensionality

  implicit none

  double precision, intent(out), dimension(ni) :: w
  double precision, intent(out), dimension(ni,ndim) :: dwdx
  double precision, intent(in), dimension(ni) :: rij
  !f2py intent(in,inplace) :: rij
  double precision, intent(in), dimension(ni,ndim) :: drij
  !f2py intent(in,inplace) :: drij
  double precision, intent(in), dimension(n) :: sml 
  !f2py intent(in,inplace) :: sml
  integer, intent(in), dimension(ni,2) :: nlist
  
  integer :: n, k, i, j, ni, ndim, kernel_type

  do k = 1, ni
    i = nlist(k,1)
    j = nlist(k,2)

    call calc_kernel(rij(k), drij(k,:), &
    & (sml(i) + sml(j)) / 2., w(k), dwdx(k,:), kernel_type,ndim)

  end do

end subroutine smoothing_kernels


subroutine adapt_kernel_width(sml,smli,rho,n)
  ! Use Sigalotti ADKE to adapt smoothing lengths
  ! Use ADKE (Sigalotti et al 2005)
  ! rho --  initial pilot densities
  ! h -- smoothing length
  ! hi -- initial smoothing length
  ! n -- number of particles
  ! ndim -- dimensionality

  implicit none
  double precision, intent(inout), dimension(n) :: sml 
  double precision, intent(in), dimension(n) :: rho
  double precision, intent(in), dimension(n) :: smli 
  double precision :: rhoav, eps
  integer :: n, i

  ! Compute a kernel based on the initial smoothing length 
  rhoav = 0.0
  eps = -0.5

  do i=1,n
    rhoav = rhoav + rho(i)
  end do
  rhoav = rhoav/n

  ! rescale the smoothing lengths based on each particle's
  ! density ratio
  ! here we use the arithmetic mean...
  do i=1,n
    sml(i) =  smli(i) * ( (rho(i)/rhoav)**(eps)  )
  end do   

end subroutine adapt_kernel_width


subroutine adapt_smoothing_length(w,dwdx,rij,drij,nlist,m,sml,smli &
  & ,kernel_type,n,ni,ndim)
  
  ! Use Sigalotti ADKE to adapt smoothing lengths
  ! Use ADKE (Sigalotti et al 2005)
  ! w -- output kernel values for each pair
  ! dwdx -- output kernel gradient for each pair
  ! rij -- pair distances
  ! drij -- pair displacements
  ! nlist - integer of pair indices
  ! h -- smoothing length
  ! hi -- initial smoothing length
  ! n -- number of particles
  ! ni -- number of pairs
  ! ndim -- dimensionality

  ! This may be redundant - no longer called?

  implicit none
  double precision, intent(inout), dimension(ni) :: w
  double precision, intent(inout), dimension(ni,ndim) :: dwdx
  double precision, intent(in), dimension(ni) :: rij
  double precision, intent(in), dimension(ni,ndim) :: drij
  integer, intent(in), dimension(ni,2) :: nlist
  double precision, intent(inout), dimension(n) :: sml 
  double precision, intent(in), dimension(n) :: m
  
  double precision, dimension(n) :: rho
  double precision :: rhoav,selfdens,r0,smli
  integer :: n, k, i, j, ni, ndim, kernel_type
  double precision, dimension(ndim) :: dx0,dw0

  ! Compute a kernel based on the initial smoothing length 

  do k = 1, ni
  call calc_kernel(rij(k), & 
                & drij(k,:), &
                & smli, &
                & w(k),dwdx(k,:), &
                  kernel_type,ndim)
  end do

  ! use this kernel to compute a pilot density
  r0 = 0.0
  dx0 = 0.0
  dw0 = 0.0
  rhoav = 0.0
  do i=1,n
    call calc_kernel(r0,dx0,smli,selfdens,dw0,kernel_type,ndim)
    rho(i) = selfdens * m(i)
  end do

  do k = 1, ni
    i = nlist(k,1)
    j = nlist(k,2)
    rho(i) = rho(i) + m(j) * w(k)
    rho(j) = rho(j) + m(i) * w(k)
  end do

  do i=1,n
    rhoav = rhoav + rho(i)
  end do
  rhoav = rhoav/n

  ! rescale the smoothing lengths based on each particle's
  ! density ratio
  do i=1,n
    sml(i) = 1.0 * smli * 1/sqrt(rho(i)/rhoav)  
  end do   

end subroutine adapt_smoothing_length

subroutine density_sum(rho, grad_rho, nlist, sml, mass, w, &
& grad_w, kernel_type, &
& n, ni, ndim)
  !! Calculates density and spatial gradient of density using
  !! the summation density method. This is unlike the one in
  !! density.f90 as the kernel normalisation has been removed.
  !! grad_w is assumed to be the neighbour list module grad_w
  !! dwdx[i,j] -- kernel gradient at position of j with respect to particle i
  !! nabla_j W_{ij}
  implicit none
  !f2py intent(in,out,overwrite) rho
  double precision, intent(inout), dimension(n) :: rho
  !f2py intent(in,out) grad_rho
  double precision, intent(inout), dimension(n,ndim) :: grad_rho
  integer, intent(in), dimension(ni,2) :: nlist
  double precision, intent(in), dimension(n) :: sml
  double precision, intent(in), dimension(n) ::  mass
  double precision, intent(in), dimension(ni) :: w
  double precision, intent(in), dimension(ni,ndim) :: grad_w
  integer :: n, ni, ndim, kernel_type
  
  double precision :: mij
  double precision,dimension(ndim) :: zdist
  integer i,j,k,d,q
  double precision selfkern, selfkg(ndim), r, wi(n)

  zdist = 0.e0
  selfkg = 0.e0
  selfkern = 0.e0
  r = 0.e0

  do i=1,n
    do d=1,ndim
      grad_rho(i,d) = 0.e0
    end do
  end do

  ! Integrate the density over space
  do i=1,n
    call calc_kernel(r,zdist,sml(i),selfkern,selfkg,kernel_type,ndim)
    rho(i) = selfkern * mass(i)
  end do

  ! SPH sum for rho
  do k = 1, ni
    i = nlist(k,1)
    j = nlist(k,2)
    rho(i) = rho(i) + mass(j)*w(k)
    rho(j) = rho(j) + mass(i)*w(k)
   
    !mij = (mass(j)+mass(i))/2.
    ! Density gradient
    ! Symmetrise?
    do d=1,ndim
      grad_rho(i,d) = grad_rho(i,d) - mass(j) * grad_w(k,d)
      grad_rho(j,d) = grad_rho(j,d) + mass(i) * grad_w(k,d)
      !grad_rho(i,d) = grad_rho(i,d) - mij * grad_w(k,d)
      !grad_rho(j,d) = grad_rho(j,d) + mij * grad_w(k,d)

      !ggx = (mij/rhoij) * (grad_rho(j,1) - grad_rho(i,1)) * dwdx(k,1)
      !ggy = (mij/rhoij) * (grad_rho(j,2) - grad_rho(i,2)) * dwdx(k,2)

    end do
  end do

end subroutine density_sum

subroutine add_boundary_density(rho, grad_rho, nlist, sml, bmass, w, grad_w, &
& kernel_type, n, ni, ndim)
  ! Add the boundary contribution to the density
  ! bmass - boundary particle masses
  ! w, gradw are the kernels and kernel gradients for the particle/boundary
  ! interaction

  implicit none
  double precision, intent(inout), dimension(n) :: rho
  !f2py intent(in,out) rho
  double precision, intent(inout), dimension(n,ndim) :: grad_rho
  !f2py intent(in,out) grad_rho
  integer, intent(in), dimension(ni,2) :: nlist
  double precision, intent(in), dimension(n) :: sml
  double precision, intent(in), dimension(n) ::  bmass
  double precision, intent(in), dimension(ni) :: w
  double precision, intent(in), dimension(ni,ndim) :: grad_w
  integer :: n, ni, ndim, kernel_type
  integer i,j,k,d
 
  do k = 1, ni
    i = nlist(k,1)
    j = nlist(k,2)
    rho(i) = rho(i) + bmass(j) * w(k)
    !Density gradient
    do d=1,ndim
      grad_rho(i,d) = grad_rho(i,d) + bmass(j) * grad_w(k,d)
    end do
  end do

end subroutine add_boundary_density

!subroutine adapt_bandwidth(sml,rho,rho_ini,n,ndim)
  ! adapt the smoothing length given the
  ! smoothing lengths for all particles,
  ! densities, and pilot densities
  !f2py intent(in,out,overwrite) :: sml
!  double precision, intent(inout), dimension(n) :: sml(:)
!  double precision, intent(in), dimension(n) :: rho(:)
!  double precision, intent(in), dimension(n) :: rho_ini(:)
!  integer :: n, ndim
!  integer :: i
!  n = size(rho)
!   do i=1,n
!     sml(i) = sml(i)*( rho_ini(i)/rho(i) )**(1.0/ndim)
!   end do
!end subroutine adapt_bandwidth


end module kernel
