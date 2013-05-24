!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----|
! Force module for Smooth Particle Mechanics code SPAC
! Subroutines for computing the pressure tensor, forces and rates of
! change for the particles.
!
! References and Influences:
! Liu and Liu
! Centre for Advanced Computing in Engineering and Science,
! "Smoothed Particle Hydrodynamics, a meshfree particle method" 
! World Scientific, 2003.
!
! William Hoover
! Smooth Particle Applied Mechanics - the state of the art
! World Scientific, 2005
!
! Copyright Andrew Charles and RMIT University, all rights reserved.

module sphforce3d

use eos
!use global
use art_viscosity
use core_potential
!use particle
use splib

double precision :: thermalk = 1.0
double precision :: cgrad = 0.0
double precision :: core_sigma = 0.0
double precision :: core_rcoef = 0.0
double precision :: inner_cutoff = 0.1

contains

subroutine set_vdw_params(adash_in,bdash_in,kbdash_in)
  ! Sets the parameters in the imported eos module
  ! These are not accessible via the interface that f2py builds
  double precision, intent(in) :: adash_in 
  double precision, intent(in) :: bdash_in
  double precision, intent(in) :: kbdash_in
  adash = adash_in
  bdash = bdash_in
  kbdash = kbdash_in
end subroutine set_vdw_params

subroutine print_vdw_params()
  ! This function enables us to inspect the values
  print *,adash
  print *,bdash
  print *,kbdash
end subroutine print_vdw_params

subroutine calc_sphforce3d(ilist,x,v,a,   &
  &   p_rev,p_rev_lr,pi_irr,              &
  &   grad_rho,grad_v,                    &
  &   u,dedt,m,rho,rho_lr,temp,q,         &
  &   c,eta,zeta,                         &
  &   dv,dr,dx,                           &
  &   sml,sml_long,w,dwdx,dwdx_long,      &
  &   avisc_on,                           &
  &   n,ni)
!  &   adash_in,bdash_in,kbdash_in,        &

  ! Computes SPH rates of change
  ! ilist -- integer array of pairs 
  ! particle variables
      ! x,v,a -- positions, velocities, accelerations
      ! p, pco -- repulsive and cohesive pressure from equation of state
      ! pressure tensor
        ! p_rev -- p_rev is the reversible, nonisotropic pressure tensor
        ! p_rev_lr -- long range reversible pressure
        ! pi_irr -- the total irreversible pressure tensor
        ! Other components are temporary internal variables
      ! u -- internal energies
      ! dedt -- rate of change of specific internal energy
      ! m -- masses
      ! rho -- densities
      ! rho_lr -- density computed with long smoothing length
      ! temp -- temperatures
      ! q -- heat flux
      ! c -- particle sound speed
      ! cgrad -- density gradient coefficient
      ! eta -- dynamic shear viscosity
      ! zeta -- dynamic bulk viscosity
  
  ! dv -- pair difference in velocity
  ! dr -- distances between pairs
  ! dx -- displacements between pairs
  
  ! sml -- smoothing lengths
  ! sml_long -- smoothing lengths
  ! dwdx -- an array of kernel derivatives and
  ! dwdx_long -- long range kernel derivatives for seperated
  ! avisc_on -- whether artificial viscosity is on or not
  ! n -- number of particles

  implicit none
  integer, intent(in), dimension(ni,2) :: ilist
  double precision, intent(in), dimension(n,3) :: x 
  double precision, intent(in), dimension(n,3) :: v
  double precision, intent(inout), dimension(n,3) :: a
  double precision, intent(inout), dimension(n,3,3) :: pi_irr
  double precision, intent(inout), dimension(n,3,3) :: p_rev
  double precision, intent(inout), dimension(n,3,3) :: p_rev_lr
  double precision, intent(inout), dimension(n,3) :: q
  double precision, intent(in), dimension(n) :: u
  double precision, intent(in), dimension(n) :: m
  double precision, intent(inout), dimension(n) :: rho
  double precision, intent(inout), dimension(n) :: rho_lr 
  double precision, intent(in), dimension(n) :: c
  double precision, intent(in), dimension(n) :: eta
  double precision, intent(inout), dimension(n) :: temp
  double precision, intent(inout), dimension(n) :: dedt 
  double precision, intent(in), dimension(n) :: zeta

  ! dr -- pair distance
  ! dx -- pair displacement
  ! dv -- velocity difference
  double precision, intent(in), dimension(ni) :: dr 
  double precision, intent(in), dimension(ni,3) :: dx 
  double precision, intent(in), dimension(ni,3) :: dv

  double precision, intent(in), dimension(n,3) :: grad_rho
  !double precision, intent(in), dimension(n,3) :: grad_rho_lr ! no longer used
  double precision, intent(inout), dimension(n,3,3) :: grad_v
  double precision, intent(in), dimension(n) :: sml
  double precision, intent(in), dimension(n) :: sml_long
  double precision, intent(in), dimension(ni) :: w
  double precision, intent(in), dimension(ni,3) :: dwdx
  double precision, dimension(ni,3) :: neg_dwdx
  double precision, intent(in), dimension(ni,3) :: dwdx_long
  logical, intent(in) :: avisc_on
  integer :: n
  integer :: ni
  integer :: ndim

  ! INTERNAL VARIABLES
  ! grad_v_os -- symmetric traceless velocity gradient, used to calc pi_os
  ! div_v -- velocity divergence
  ! pi_os -- symmetric traceless viscous (irreversible) pressure tensor
  ! pi_one -- isotropic viscous (irreversible) pressure tensor
  ! pi_eq -- reversible, isotropic pressure tensor 
  !   (accumulates core potential and avisc as well)
  double precision, dimension (n,3,3) :: grad_v_os
  double precision, dimension(n) :: div_v

  double precision, dimension(n,3,3) :: pi_os
  double precision, dimension(n) :: pi_one
  double precision, dimension(n) :: p_eq
  
  double precision, dimension(n) :: p
  double precision, dimension(n) :: pco 
  
  ! running totals
  double precision :: he, hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz, hx,hy,hz

  !> average density
  double precision :: rho_ij,c_ij,sml_ij

  double precision ::  rhosq_i
  double precision ::  rhosq_j
  double precision ::  rhosq_lr_i
  double precision ::  rhosq_lr_j

  !> artifical viscosity for each interaction
  double precision :: avisc, dux, duy, duz
!  double precision :: adash_in,bdash_in,kbdash_in

  !> core repulsion
  double precision, dimension(3) :: repulsion

  !> f is used in the monaghan heat flux calculation
  double precision :: f
  integer :: i,j,k

  !call set_vdw_params(adash_in,bdash_in,kbdash_in)
  !print *,kbdash
  !print *,cgrad

  ndim = 3

  ! initialise tensor arrays to zero
  pi_os = 0.0
  pi_one = 0.0
  pi_irr = 0.0
  p_rev = 0.0
  p_rev_lr = 0.0
  p_eq = 0.0
  p = 0.0
  pco = 0.0

  ! intialise some other stuff to zero
  avisc = 0.0

  ! set rates of change and spatial gradients to zero
  a = 0.0
  q = 0.0
  dedt = 0.0
  grad_v = 0.0

  if(ni .eq. 0) return
  neg_dwdx = -dwdx
  ! We pass negative dwdx because calc_grad_v expects grad_w_i_dij
  call calc_grad_v3d(grad_v,ilist,neg_dwdx,dv,m,rho,n,ndim,ni)
  call calc_div_v3d(grad_v,div_v,n,ndim)
  call calc_grad_v_os3d(grad_v_os,grad_v,div_v,n,ndim)

  ! Calculate pressure from equation of state
  call calc_vdw_hc_pressure(p,rho,u,n)
  call calc_vdw_cohesive_pressure(pco,rho_lr,u,n)

  p_eq = p

  ! use dwdx_long for capillary pressure
  ! todo clarfy sign of dwdx and write explanation here 
  ! 26 Feb 2012 flipped this sign.
  ! Also need to flip sign of the force, as the
  ! resulting force seems to be correctly oriented.
  ! neg_dwdx_long = -dwdx_long
  ! todo this really needs to be clarified throught the code and in the docs
  ! unless otherwise noted, dwdx is
  ! \Nabla_j W_ij
  ! This is the form it is returned from the kernel module by default
  call calc_capillary_pressure3d(ilist,rho,m,grad_rho &
    & ,dwdx,p_rev,cgrad,n,ni)

  p_rev(:,1,1) = p_rev(:,1,1) + p(:)
  p_rev(:,2,2) = p_rev(:,2,2) + p(:)
  p_rev(:,3,3) = p_rev(:,3,3) + p(:)
  !print *,pco
  !print *,p_rev_lr(:,:,:)
  p_rev_lr(:,1,1) = p_rev_lr(:,1,1) + pco(:)
  p_rev_lr(:,2,2) = p_rev_lr(:,2,2) + pco(:)
  p_rev_lr(:,3,3) = p_rev_lr(:,3,3) + pco(:)
  !print *,p_rev_lr(:,:,:)
  
  ! calculate shear stress  
  call calc_pi_os3d(pi_os,grad_v_os,eta,n,ndim)
   
  !calculate bulk stress
  call calc_pi_one(pi_one,div_v,zeta,n)

  ! Add bulk stress to shear stress
  pi_irr = pi_os
  pi_irr(:,1,1) = pi_irr(:,1,1) + pi_one(:)
  pi_irr(:,2,2) = pi_irr(:,2,2) + pi_one(:)
  pi_irr(:,3,3) = pi_irr(:,3,3) + pi_one(:)

  call calc_heat_flux_3d(q,ilist,rho,m,temp,neg_dwdx,thermalk,n,ni)

  ! calc SPH sum for pressure force

  do k = 1, ni
    i = ilist(k,1)
    j = ilist(k,2)

    ! hx and hy are running totals of the pair force
    hx = 0.e0
    hy = 0.e0
    hz = 0.e0
    hxx = 0.e0
    hxy = 0.e0
    hxz = 0.e0
    hyx = 0.e0
    hyy = 0.e0
    hyz = 0.e0
    hzx = 0.e0
    hzy = 0.e0
    hzz = 0.e0

    ! he is a running total of the pair work (reversible)
    he = 0.e0
    
    rho_ij = (rho(i) + rho(j)) / 2.0 
    c_ij = (c(i) + c(j)) / 2.0
    sml_ij = (sml(i) + sml(j)) / 2.0

    ! Calculate artificial viscosity
    if (avisc_on) then
      call mon_beta_avisc(avisc,dx(k,:),dv(k,:),rho_ij,sml_ij,ndim)
    endif

    ! Calculate core force
    call core_force(repulsion,dx(k,:),dr(k)*dr(k),core_sigma,core_rcoef,ndim)


    ! Note that dwdx is grad_j Wij, so the h's are directed so as
    ! to contribute to particle j's acceleration
   
    rhosq_i = rho(i)**2
    rhosq_j = rho(j)**2
    rhosq_lr_i = rho_lr(i)**2
    rhosq_lr_j = rho_lr(j)**2

    ! Short range reversible contribution
    ! hx == (Previ/rho^2 + Prevj/rhoj^2)* dwdx
    ! (P dot w)i = sum_j Pij Wj 
    hxx = ( p_rev(i,1,1) / rhosq_i               &
      &   + p_rev(j,1,1) / rhosq_j ) * dwdx(k,1)
    
    hxy = ( p_rev(i,1,2) / rhosq_i               &
      &   + p_rev(j,1,2) / rhosq_j ) * dwdx(k,2)    
    
    hxz = ( p_rev(i,1,3) / rhosq_i               &
      &   + p_rev(j,1,3) / rhosq_j ) * dwdx(k,3)    
          
    ! y coordinate of acceleration
    ! Short range reversible contribution
    hyx = ( p_rev(i,2,1) / rhosq_i               &
      &   + p_rev(j,2,1) / rhosq_j ) * dwdx(k,1)   
      
    hyy = ( p_rev(i,2,2) / rhosq_i               &
      &   + p_rev(j,2,2) / rhosq_j ) * dwdx(k,2)
    
    hyz = ( p_rev(i,2,3) / rhosq_i               &
      &   + p_rev(j,2,3) / rhosq_j ) * dwdx(k,3)

    ! Short range z acceleration
    hzx = ( p_rev(i,3,1) / rhosq_i               &
      &   + p_rev(j,3,1) / rhosq_j ) * dwdx(k,1)   
      
    hzy = ( p_rev(i,3,2) / rhosq_i               &
      &   + p_rev(j,3,2) / rhosq_j ) * dwdx(k,2)
    
    hzz = ( p_rev(i,3,3) / rhosq_i               &
      &   + p_rev(j,3,3) / rhosq_j ) * dwdx(k,3)

    ! Long range reversible contribution to acceleration
    ! x 
    hxx = hxx + (p_rev_lr(i,1,1)/rhosq_lr_i     &
      &       +  p_rev_lr(j,1,1)/rhosq_lr_j )   &
      &       * dwdx_long(k,1)                        
    hxy = hxy + (p_rev_lr(i,1,2)/rhosq_lr_i     &
      &       +  p_rev_lr(j,1,2)/rhosq_lr_j )   &
      &       * dwdx_long(k,2)
    hxz = hxz + (p_rev_lr(i,1,3)/rhosq_lr_i     &
      &       +  p_rev_lr(j,1,3)/rhosq_lr_j )   &
      &       * dwdx_long(k,3)                        

    ! y 
    hyx = hyx + (p_rev_lr(i,2,1)/rhosq_lr_i     &
      &       +  p_rev_lr(j,2,1)/rhosq_lr_j )   &
      &       * dwdx_long(k,1)                        
    hyy = hyy + (p_rev_lr(i,2,2)/rhosq_lr_i     & 
      &       +  p_rev_lr(j,2,2)/rhosq_lr_j )   &
      &       * dwdx_long(k,2)
    hyz = hyz + (p_rev_lr(i,2,3)/rhosq_lr_i     & 
      &       +  p_rev_lr(j,2,3)/rhosq_lr_j )   &
      &       * dwdx_long(k,3)

    ! z 
    hzx = hzx + (p_rev_lr(i,3,1)/rhosq_lr_i     &
      &       +  p_rev_lr(j,3,1)/rhosq_lr_j )   &
      &       * dwdx_long(k,1)                        
    hzy = hzy + (p_rev_lr(i,3,2)/rhosq_lr_i     & 
      &       +  p_rev_lr(j,3,2)/rhosq_lr_j )   &
      &       * dwdx_long(k,2)
    hzz = hzz + (p_rev_lr(i,3,3)/rhosq_lr_i     & 
      &       +  p_rev_lr(j,3,3)/rhosq_lr_j )   &
      &       * dwdx_long(k,3)
    
    ! Add artificial viscosity to the reversible pressure bit
    ! by doing this the artificial viscosity is
    ! conservative of energy, not dissipative.
    hxx = hxx + avisc * dwdx(k,1) 
    hyy = hyy + avisc * dwdx(k,2)
    hzz = hzz + avisc * dwdx(k,3)

    ! x coodrdinate of acceleration 
    ! irreversible contribution
    hxx = hxx + (pi_irr(i,1,1)/(rhosq_i)   &
      &       +  pi_irr(j,1,1)/(rhosq_j) ) &
      &       * dwdx(k,1)      
    hxy = hxy + (pi_irr(i,1,2)/(rhosq_i)   &
      &       +  pi_irr(j,1,2)/(rhosq_j) ) &
      &       * dwdx(k,2)
    hxz = hxz + (pi_irr(i,1,3)/(rhosq_i)   &
      &       +  pi_irr(j,1,3)/(rhosq_j) ) &
      &       * dwdx(k,3)
          
    ! y coordinate of acceleration
    ! irreversible contribution
    hyx = hyx + (pi_irr(i,2,1)/(rhosq_i)     &
      &       +  pi_irr(j,2,1)/(rhosq_j) )   &
      &       * dwdx(k,1)   
    hyy = hyy + (pi_irr(i,2,2)/(rhosq_i)     &
      &       +  pi_irr(j,2,2)/(rhosq_j) )   &
      &       * dwdx(k,2)
    hyz = hyz + (pi_irr(i,2,3)/(rhosq_i)     &
      &       +  pi_irr(j,2,3)/(rhosq_j) )   &
      &       * dwdx(k,3)

    ! z coord irreversible
    hzx = hzx + (pi_irr(i,3,1)/(rhosq_i)     &
      &       +  pi_irr(j,3,1)/(rhosq_j) )   &
      &       * dwdx(k,1)   
    hzy = hzy + (pi_irr(i,3,2)/(rhosq_i)     &
      &       +  pi_irr(j,3,2)/(rhosq_j) )   &
      &       * dwdx(k,2)
    hzz = hzz + (pi_irr(i,3,3)/(rhosq_i)     &
      &       +  pi_irr(j,3,3)/(rhosq_j) )   &
      &       * dwdx(k,3)

    !! hx is not a force, it is the
    !! term in the sph momentum equation sans the mass
   
    !! The sign here is reversed, because dwdx(k,:) is actually
    !! grad_j W_ij = -grad_i W_ij
    !! In the SPH equations of motion, grad_i W_ij appears.
    !! Reversing this convention throughout the code would be
    !! time consuming, and may not even make sense for modules
    !! like neighbour list.
    !! Core repulsion is a straight force so we need to divide by mass.
    hx = hxx + hxy + hxz
    hy = hyy + hyx + hyz
    hz = hzz + hzy + hzx

    !! If the particles are closer than inner_cutoff
    !! Only apply the core repulsion
    if(dr(i) .lt. inner_cutoff) then
      hx = 0.0
      hy = 0.0
      hz = 0.0
    endif 
    
    a(i,1) = a(i,1) + (m(j)*(hx) + repulsion(1)/m(i))
    a(j,1) = a(j,1) - (m(i)*(hx) + repulsion(1)/m(j))
    
    a(i,2) = a(i,2) + (m(j)*(hy) + repulsion(2)/m(i))
    a(j,2) = a(j,2) - (m(i)*(hy) + repulsion(2)/m(j))

    a(i,3) = a(i,3) + (m(j)*(hz) + repulsion(3)/m(i))
    a(j,3) = a(j,3) - (m(i)*(hz) + repulsion(3)/m(j))

    !looking for nans
    !print *,'f90',a(i,:)
    !print *,'f90',a(j,:)

    !! If the particles are closer than inner_cutoff
    !! Only apply the core repulsion
    if(dr(i) .lt. inner_cutoff) then
      a(i,1) = 0.0
      a(i,2) = 0.0
      a(i,3) = 0.0
      a(j,1) = 0.0
      a(j,2) = 0.0
      a(j,3) = 0.0
    endif 

    ! both reversible (PV) work and work due to the irreversible parts
    ! of the pressure tensor (shear and bulk viscosity)
    ! P/rho^2 dot vij
    ! he accumulates du/dt terms without the multiplication
    ! by particle mass
    ! dv here is v[j] - v[i]
    ! hxx uses grad[j] W[i,j]
    ! So the he term is signed such that it is added to dedt.
    
    ! Slower he
    !he = he + (hxx * dv(k,1)) &
    !  &     + (hxy * dv(k,2)) &
    !  &     + (hxz * dv(k,3)) &
    !  &     + (hyx * dv(k,1)) &
    !  &     + (hyy * dv(k,2)) &
    !  &     + (hyz * dv(k,3)) &
    !  &     + (hzx * dv(k,1)) &
    !  &     + (hzy * dv(k,2)) &
    !  &     + (hzz * dv(k,3))

    ! Faster he
    he = he + (hx * dv(k,1)) &
      &     + (hy * dv(k,2)) &
      &     + (hz * dv(k,3)) 

    !needs to be symmetrised consistenty
    dedt(i) = dedt(i) + m(j)*he/2.
    dedt(j) = dedt(j) + m(i)*he/2.

    ! Heat conduction
    ! ---------------

    ! subroutine compute_conduction(
    ! only ful heat flux conduction present for now.
    ! other varieties should be in a subroutine.
    ! probably to use monaghan style we want to
    ! use this shortcut to compute the heat flux and
    ! do this step the same anyway.
  
    !Using full heat flux vector - ref Nugent and Posch
    !calculate contribution to dedt
    ! The sign here does not need to be opposite because the heat
    ! flux code takes grad[j]W[i,j]
    dux = ( q(i,1)/rhosq_i + q(j,1)/rhosq_j ) * dwdx(k,1)
    duy = ( q(i,2)/rhosq_i + q(j,2)/rhosq_j ) * dwdx(k,2)
    duz = ( q(i,3)/rhosq_i + q(j,3)/rhosq_j ) * dwdx(k,3)
    dedt(i) = dedt(i) - m(j) * (dux + duy + duz) 
    !dedt(i) = dedt(i) - m(j) * duy 
    !dedt(i) = dedt(i) - m(j) * duz 
    dedt(j) = dedt(j) + m(i) * (dux + duy + duz)
    !dedt(j) = dedt(j) + m(i) * duy
    !dedt(j) = dedt(j) + m(i) * duz
    
  enddo !end loop over pairs

  !failed experiment in simplifying dedt
  !! \todo investigate this further (not zeroing grad v was most likely 
  !! the issue. For realtime graphics this approach should be much cheaper i
  !! than the whole -tensor approach.
  !! Prototype in python first.
  !! do i=1, n
  !! dedt_rev1(i) = -p_rev(i)*div_v(i)/rho(i)
  !! dedt(i) = dedt(i) + tdsdt(i)
  !! enddo

end subroutine calc_sphforce3d
end module sphforce3d
