!---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----|
! Force module for Smooth Particle Mechanics code SPAC
! Subroutines for computing the pressure tensor, forces and rates of
! change for the particles.
! This module is heavily borrowed from code published in:
!
! Smoothed Particle Hydrodynamics, a meshfree particle method.
! Liu and Liu
! Centre for Advanced Computing in Engineering and Science,
! a meshfree particle method, World Scientific, 2003.
! and 
! William Hoover
! Smooth Particle Applied Mechanics - the state of the art
! World Scientific, 2005
! Copyright Andrew Charles, all rights reserved.

module sphforce

use eos
use global
use neighbour_list
use simulation_box
use art_viscosity
use core_potential
use particle
use boundary
use splib

public :: v_smooth

contains

subroutine calc_sphforce(nl,n,ndim,dx,p,dwdx,tdsdt,dwdx_long)
  ! Computes SPH rates of change
  ! nl -- neighbour list data structure
  ! box -- box data structure
  ! n -- number of particles
  ! dimn -- dimension of system - currently supports 1 or 2
  ! dx -- displacements between pairs
  ! p -- particle data structure
      ! x,v,a -- positions, velocities, accelerations
      ! u -- internal energies
      ! m -- masses
      ! dedt -- rate of change of specific internal energy
      ! rho -- densities
      ! sml -- smoothing lengths
      ! c -- particle sound speed
      ! eta -- dynamic shear viscosity
      ! zeta -- dynamic bulk viscosity
      ! temp -- temperatures
      ! pco -- cohesive pressure for seperated eos
      ! pressure tensor
      ! symmetric traceless viscous (irreversible) pressure tensor
      ! pi_os
      ! pi_irr is the total irreversible pressure tensor
      ! isotropic viscous (irreversible) pressure tensor
      ! double precision, dimension (n) :: pi_one
      ! reversible, isotropic pressure tensor
      ! this accumulates core potential and avisc as well
      ! double precision, dimension(n,dimn,dimn) :: p_rev
      ! the isotropic part of the equilibrium pressure
      ! double precision, dimension(n) :: p_eq
  ! tdsdt -- entropy production
  ! dwdx -- an array of kernel derivatives and
  ! dwdx_long -- long range kernel derivatives for seperated

  implicit none

  type (nlist_data_type), intent(inout) :: nl
  integer :: n
  integer :: ndim
  double precision, dimension(:,:), intent(in) :: dx
  type (particle_type), intent(inout) :: p
  double precision, intent(in), dimension(:,:) :: dwdx
  double precision, dimension(:), intent(inout)::  tdsdt
  double precision, dimension(:,:), optional :: dwdx_long

  ! optional boundary structures
  !type (nlist_data_type), intent(inout), optional :: bnl
  !type (particle_type), intent(inout), optional :: bp
  
  ! INTERNAL VARIABLES
  ! symmetric traceless velocity gradient, used to calc pi_os
  double precision, dimension (n,ndim,ndim) :: grad_v_os

  double precision, dimension(n) :: div_v
  
  ! interaction list
  integer, dimension(nl%nip,2) :: ilist

  !> running totals
  !! A keen observer will note the similarity in name with Liu and Liu's
  !! code and deduce correctly that i owe them a substantial debt.
  double precision :: he, hxx, hxy, hyx, hyy

  !> average density
  double precision :: rho_ij,c_ij,sml_ij

  !> artifical viscosity for each interaction
  double precision :: avisc

  !> core repulsion
  double precision, dimension(ndim) :: repulsion

  !> f is used in the monaghan heat flux calculation
  double precision :: f

  integer :: i,j,k
  integer :: ni

  !benchmarking
  integer :: start_time,end_time,gradv_time,stress_time,time_rate
  real :: time_seconds
  if(profile) then
      call system_clock(count=start_time,count_rate=time_rate)
  endif


  ! initialise tensor arrays to zero
  p%pi_os = 0
  p%pi_one = 0
  p%pi_irr = 0
  p%p_rev = 0
  p%p_rev_lr = 0
  p%p_eq = 0
  p%p = 0
  p%pco = 0

  ! intialise some other stuff to zero
  tdsdt = 0.0
  avisc = 0.0

  ! set rates of change and spatial gradients to zero
  p%a = 0.0
  p%q = 0
  p%dedt = 0
  p%grad_v = 0

  ni = nl%nip
  ! Create a local copy of nl%ilist to pass to subroutines
  do k = 1,ni
    ilist(k,1) = nl%ilist(k,1) 
    ilist(k,2) = nl%ilist(k,2)
  enddo

  if(debug) print *, 'calculating sphforces '
 
  ! Correct form for an nlist loop: 
  ! ni = nl%nip
  ! do k = 1, ni
  !   i = nl%ilist(k,1)
  !   j = nl%ilist(k,2)
  ! enddo

  if(profile) call system_clock(count=gradv_time)

  !print *,nl%nip
  !print *,ni
  !print *,size(dwdx(:,1))
  !print *,size(ilist(:,1))
  if(ni .eq. 0) return

  ! We pass negative dwdx because calc_grad_v expects 
  call calc_grad_v(p%grad_v,ilist(1:ni,:),-dwdx(1:ni,:),nl%dv(1:ni,:),p%m, &
    & p%rho,p%n,ndim,nl%nip)
  call calc_div_v2d(p%grad_v,div_v,p%n,ndim)
  call calc_grad_v_os2d(grad_v_os,p%grad_v,div_v,p%n,ndim)

  if(validate) then
    call validate_double(p%grad_v(1,1,2),'gv')
  endif

  if(profile) then
    call system_clock(count = end_time)
    time_seconds = (end_time - gradv_time)/real(time_rate)
    print *,'**** div,grad v time ',time_seconds
  endif
    
  ! calculate pressure from equation of state
  call calc_iso_pressure(p) !%rho,p%u,p%p,p%c,p%grad_rho)
  call calc_long_range_pressure(p) !%n,p%rho_lr,p%u,p%pco)

  !p_eq = p%p + p%pco
  p%p_eq = p%p

  ! use dwdx_long for capillary pressure 
!subroutine calc_capillary_pressure(ilist,rho,m,grad_rho,dwdx,p_rev,cgrad,n,ni,ndim)
  call calc_capillary_pressure2d(ilist(1:ni,:),p%rho_lr,p%m,p%grad_rho &
    & ,dwdx_long(1:ni,:),p%p_rev_lr,cgrad_ini,p%n,ni,ndim)

  !do i=1,size(p_eq)
  p%p_rev(:,1,1) = p%p_rev(:,1,1) + p%p(:)
  p%p_rev(:,2,2) = p%p_rev(:,2,2) + p%p(:)
  !end do

  p%p_rev_lr(:,1,1) = p%p_rev_lr(:,1,1) + p%pco(:)
  p%p_rev_lr(:,2,2) = p%p_rev_lr(:,2,2) + p%pco(:)
  
  if(profile) call system_clock(count=stress_time)
  
  ! calculate shear stress  
  call calc_pi_os2d(p%pi_os,grad_v_os,p%eta,p%n,ndim)
   
   if(validate) then
      call validate_double(p%pi_os(1,1,2),'pios')
   endif

  !calculate bulk stress
  call calc_pi_one(p%pi_one,div_v,p%zeta,p%n)

  ! Add bulk stress to shear stress
  p%pi_irr = p%pi_os
  ! do i=1,size(pi_one)
  p%pi_irr(:,1,1) = p%pi_irr(:,1,1) + p%pi_one(:)
  p%pi_irr(:,2,2) = p%pi_irr(:,2,2) + p%pi_one(:)
  ! end do

  ! calculate viscous entropy
  ! if(verbose) print *, 'calculating viscous entropy '
  ! call calc_viscous_entropy_os(ilist(1:ni,:),p%pi_irr,tdsdt,p%grad_v,p%rho &
  !  & ,n,ni,ndim)

  if(profile) then
    call system_clock(count = end_time)
    time_seconds = (end_time - stress_time)/real(time_rate)
    print *,'**** shear,bulk,viscous entropy time ',i, time_seconds
  endif

  if(conduction_type .eq. 2) then
    !Nugent and posch full heat flux style
    !calculate heat flux vector
    !try passing it smoothed temperature?
    call calc_heat_flux_old(nl,p%q,p%rho,p%m,p%temp,-nl%dwdx)
  end if

  ! calc SPH sum for pressure force
  if(verbose) print *, 'calculating pressure force '

  ni = nl%nip

  if(validate) call check_zero_distance(nl%rijsq(1:ni),ilist(1:ni,:),nl%nip)
  !! \todo make rhoij etc pair computed properties
  !! \todo pass an averaged smoothing length consistently
  do k = 1, ni
    i = ilist(k,1)
    j = ilist(k,2)

    if(debug) print *, 'loop for pair',k

    ! hx and hy are running totals of the pair force
    hxx = 0.e0
    hxy = 0.e0
    hyx = 0.e0
    hyy = 0.e0
    ! he is a running total of the pair work (reversible)
    he = 0.e0
    
    !calculate artificial viscosity
    !subroutine calc_art_viscosity(dr,dv,rho,sml,c,avisc)
    rho_ij = (p%rho(i)+p%rho(j))/2.0 
    c_ij =(p%c(i)+p%c(j))/2.0
    sml_ij = (p%sml(i)+p%sml(j))/2.0

    call calc_art_viscosity(nl%drij(k,:),nl%dv(k,:),nl%rij(k),rho_ij &
      ,sml_ij,c_ij,nl%w(k),avisc)

    call core_force(repulsion,nl%drij(k,:),nl%rijsq(k),sigma,rcoef,ndim)
 
    ! Note that dwdx is grad_j Wij, so the h's are directed so as
    ! to contribute to particle j's acceleration
    
    ! Short range reversible contribution
    ! hx == (Previ/rho^2 + Prevj/rhoj^2)* dwdx
    ! (P dot w)i = sum_j Pij Wj 
    hxx = ( p%p_rev(i,1,1)/p%rho(i)**2               &
      &   + p%p_rev(j,1,1)/p%rho(j)**2 )* dwdx(k,1)
      
    hxy = ( p%p_rev(i,1,2)/p%rho(i)**2               &
      &   + p%p_rev(j,1,2)/p%rho(j)**2 )* dwdx(k,2)    
          
    ! y coordinate of acceleration
    ! Short range reversible contribution
    hyx = ( p%p_rev(i,2,1)/p%rho(i)**2               &
      &   + p%p_rev(j,2,1)/p%rho(j)**2 )* dwdx(k,1)   
      
    hyy = ( p%p_rev(i,2,2)/p%rho(i)**2               &
      &   + p%p_rev(j,2,2)/p%rho(j)**2 )* dwdx(k,2)

    ! Long range reversible contribution to acceleration
    if(eos_type .eq. 3) then
      if( present(dwdx_long)) then
        hxx = hxx + (p%p_rev_lr(i,1,1)/p%rho_lr(i)**2     &
          &       +  p%p_rev_lr(j,1,1)/p%rho_lr(j)**2 )   &
          &       * dwdx_long(k,1)                        
        hxy = hxy + (p%p_rev_lr(i,1,2)/p%rho_lr(i)**2     &
          &       +  p%p_rev_lr(j,1,2)/p%rho_lr(j)**2 )   &
          &       * dwdx_long(k,2)                        

        hyx = hyx + (p%p_rev_lr(i,2,1)/p%rho_lr(i)**2     &
          &       +  p%p_rev_lr(j,2,1)/p%rho_lr(j)**2 )   &
          &       * dwdx_long(k,1)                        
        hyy = hyy + (p%p_rev_lr(i,2,2)/p%rho_lr(i)**2     & 
          &       +  p%p_rev_lr(j,2,2)/p%rho_lr(j)**2 )   &
          &       * dwdx_long(k,2)

      endif
    endif
    
    ! Add artificial viscosity to the reversible pressure bit
    ! by doing this the artificial viscosity is
    ! conservative of energy, not dissipative.
    hxx = hxx + avisc * dwdx(k,1) 
    hyy = hyy + avisc * dwdx(k,2)

    !print *,'hxx',hxx
    !print *,'hxy',hxy
    !print *,'hyx',hyx
    !print *,'hzz',hyy

    if(validate) then
      call validate_double(hxx,'hxx  ')
      call validate_double(he,'he  ')
    endif
   
    ! x coodrdinate of acceleration 
    ! irreversible contribution
    hxx = hxx + (p%pi_irr(i,1,1)/(p%rho(i)**2)   &
      &       +  p%pi_irr(j,1,1)/(p%rho(j)**2) )  &
      &       * dwdx(k,1)      
    hxy = hxy + (p%pi_irr(i,1,2)/(p%rho(i)**2)   &
      &       +  p%pi_irr(j,1,2)/(p%rho(j)**2) ) &
      &       * dwdx(k,2)
          
    ! y coordinate of acceleration
    ! irreversible contribution
    hyx = hyx + ( p%pi_irr(i,2,1)/(p%rho(i)**2)    &
      &       + p%pi_irr(j,2,1)/(p%rho(j)**2) ) &
      &       * dwdx(k,1)   
    hyy = hyy + ( p%pi_irr(i,2,2)/(p%rho(i)**2)    &
      &       + p%pi_irr(j,2,2)/(p%rho(j)**2) ) &
      &       * dwdx(k,2)

    if(validate) then
      call validate_double(dwdx(k,1),'gotcha')
    endif
   
    !! hx is not a force, it is the
    !! term in the sph momentum equation sans the mass
    !! repulsion is a straight force so we need to divide by mass.
   
    !! The sign here is reversed, because dwdx(k,:) is actually
    !! grad_j W_ij = -grad_i W_ij
    !! In the SPH equations of motion, grad_i W_ij appears.
    !! Reversing this convention throughout the code would be
    !! time consuming, and may not even make sense for modules
    !! like neighbour list.

    p%a(i,1) = p%a(i,1) + (p%m(j)*(hxx+hxy) + repulsion(1)/p%m(i))
    p%a(j,1) = p%a(j,1) - (p%m(i)*(hxx+hxy) + repulsion(1)/p%m(j))
    
    p%a(i,2) = p%a(i,2) + (p%m(j)*(hyy+hyx) + repulsion(2)/p%m(i))
    p%a(j,2) = p%a(j,2) - (p%m(i)*(hyy+hyx) + repulsion(2)/p%m(j))


    ! both reversible (PV) work and work due to the irreversible parts
    ! of the pressure tensor (shear and bulk viscosity)
    ! P/rho^2 dot vij
    ! he accumulates du/dt terms without the multiplication
    ! by particle mass
    ! dv here is v[j] - v[i]
    ! hxx uses grad[j] W[i,j]
    ! So the he term is signed such that it is added to dedt.
    he = he + (hxx * nl%dv(k,1)) &
      &     + (hxy * nl%dv(k,2)) &
      &     + (hyx * nl%dv(k,1)) &
      &     + (hyy * nl%dv(k,2))

    !needs to be symmetrised consistenty
    p%dedt(i) = p%dedt(i) + p%m(j)*he/2.
    p%dedt(j) = p%dedt(j) + p%m(i)*he/2.

    ! Super sleuth October 16
    if( p%a(i,1) .gt. 10000) then
      print *,'a',p%a(i,1)
      print *,'rij',nl%rij(k)
      print *,dwdx(k,1)
      print *,k
      print *,p%p_rev_lr(i,:,:)
      print *,p%rho(i)
    end if

    ! Heat conduction
    ! ---------------
    if(debug) print *, 'calculating heat conduction for ', i, j
    !the heat conduction should be counted as part of tdsdt
    !if we ever want to calculate the entropy

    ! subroutine compute_conduction(

    if(conduction_type .eq. 1) then
        !monaghan type conduction
        
      if(dimn .eq. 2) then  
        f = ( dwdx(k,2)*dx(k,2)+dwdx(k,1)*dx(k,1) ) / nl%rijsq(k)
      else if(dimn .eq. 1) then
        f = dwdx(k,1)*dx(k,1)/nl%rijsq(k)
      endif
 
      p%dedt(i) = p%dedt(i) - p%m(j) * f * ( (thermalk + thermalk) & 
                * (p%temp(i) - p%temp(j)))/(p%rho(i)*p%rho(j))
      p%dedt(j) = p%dedt(j) + p%m(j) * f * ( (thermalk + thermalk) &
                * (p%temp(i) - p%temp(j)))/(p%rho(i)*p%rho(j))
  
    else if(conduction_type .eq. 2) then
      !Using full heat flux vector - ref Nugent and Posch
      !calculate contribution to dedt
      ! The sign here does not need to be opposite because the heat
      ! flux code takes grad[j]W[i,j]
      p%dedt(i) = p%dedt(i) - p%m(j)*( p%q(i,1)/(p%rho(i)**2) &
                + p%q(j,1)/(p%rho(j)**2) )*dwdx(k,1)
      p%dedt(i) = p%dedt(i) - p%m(j)*( p%q(i,2)/(p%rho(i)**2) &
                + p%q(j,2)/(p%rho(j)**2) )*dwdx(k,2)
      
      p%dedt(j) = p%dedt(j) + p%m(i)*( p%q(i,1)/(p%rho(i)**2) &
                + p%q(j,1)/(p%rho(j)**2) )*dwdx(k,1)
      p%dedt(j) = p%dedt(j) + p%m(i)*( p%q(i,2)/(p%rho(i)**2) &
                + p%q(j,2)/(p%rho(j)**2) )*dwdx(k,2)
    endif
    

    if(validate) then
      call validate_double(p%dedt(i),'dedt')
      call validate_double(p%dedt(j),'dedt')
      call validate_double(p%a(i,1),'a')
      call validate_double(p%a(i,2),'a')
      call validate_double(p%a(j,1),'a')
      call validate_double(p%a(j,2),'a')
    endif

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

  ! change in specific internal energy de/dt = T ds/dt - p/rho
  !do i=1, n
  !  p%dedt(i) =  p%dedt(i) - tdsdt(i)
  !enddo

  if(verbose) print *, 'Finished calculating sphforces.'
     
  if(profile) then   
    call system_clock(count = end_time)
    time_seconds = (end_time - start_time)/real(time_rate)
    print *,'**** force calculation time ', time_seconds
  endif

end subroutine calc_sphforce

subroutine v_smooth(p,nl)
  ! v_smooth is a function that smoothes rdot
  type (particle_type), intent(inout) :: p
  type (nlist_data_type) :: nl
  integer i,j,k

  do k = 1, nl%nip
   i = nl%ilist(k,1)
   j = nl%ilist(k,2)

   if( (p%rho(i)+p%rho(j)) .eq. 0.0 ) stop 'divide by zero error in dx'
            
  ! removed mass:  
   p%rdot(i,:) = p%rdot(i,:) + v_eps * (p%v(j,:) - p%v(i,:)) &
             & * ((p%m(i)+p%m(j))) *nl%w(k)  / (p%rho(i)+p%rho(j)) 
   p%rdot(j,:) = p%rdot(j,:) + v_eps * (p%v(i,:) - p%v(j,:)) &
             & * ((p%m(i)+p%m(j))) *nl%w(k)  / (p%rho(i)+p%rho(j)) 

  end do

end subroutine v_smooth


subroutine calc_heat_flux_old(nl,q,rho,m,tmp,dwdx_jij)
  ! Computes the heat flux for a set of SPAM particles.
  type(nlist_data_type) :: nl
  double precision, dimension(:,:), intent(inout) :: q
  double precision, dimension(:) :: rho,m,tmp
  double precision, dimension(:,:) :: dwdx_jij
  double precision :: rhoij,mij,tij
  integer i,j,k

  q=0
      
  do k = 1, nl%nip
    i = nl%ilist(k,1)
    j = nl%ilist(k,2)
    rhoij = (rho(i)+rho(j))/2.
    mij = (m(i)+m(j))/2. 
    tij = tmp(j) - tmp(i) 
    
    q(i,1) = q(i,1) + thermalk * (mij/rhoij) * tij * dwdx_jij(k,1)
    q(i,2) = q(i,2) + thermalk * (mij/rhoij) * tij * dwdx_jij(k,2)
    q(j,1) = q(j,1) + thermalk * (mij/rhoij) * tij * dwdx_jij(k,1)
    q(j,2) = q(j,2) + thermalk * (mij/rhoij) * tij * dwdx_jij(k,2) 
  end do

end subroutine calc_heat_flux_old


subroutine calc_viscous_entropy_os(ilist,p_ir,tdsdt,grad_v,rho,n,ni,ndim)
  ! This is now a wrapper for subroutines in sphlib
  ! calculate the symmetric traceless 
  ! contribution to dedt
  ! pass in pi_os
  ! Uses the global configuration variable viscous_entropy_type to
  ! determine which subroutine it calls.
  ! type(nlist_data_type) :: nl changed to ilist
  integer, dimension(ni,2), intent(in):: ilist
  integer :: n
  double precision, dimension(n,ndim,ndim) :: grad_v
  double precision, dimension(n) :: tdsdt,rho
  double precision, dimension(n,ndim,ndim) :: p_ir
  integer:: ni, ndim

  if(debug) print *,n
  if(entropy_type .eq. 1) then
    ! NB this fragment largely untested
    call calc_viscous_entropy_liu(p_ir,tdsdt,n,dimn)
  else if (entropy_type .eq. 2) then
    !call calc_viscous_entropy_os_full_2d(ilist,p_ir,tdsdt,grad_v_os,rho,n,ni)
    call calc_viscous_entropy_full_2d(ilist,p_ir,tdsdt,grad_v,rho,n,ni)
    ! and now the bulk
  endif
end subroutine calc_viscous_entropy_os


subroutine check_zero_distance(drsq,ilist,ni)
  ! Given an array of squared interparticle distances,
  ! indexed by ilist, check that none of the distances
  ! are zero
  ! drsq -- squared distances
  ! ilist -- indices of each interacting pair
  ! ni -- number of interacting pairs (not given by size of ilist)
  double precision, dimension(:) :: drsq
  integer, dimension(ni,2) :: ilist
  integer ni,k,i,j
    
  do k = 1, ni
      i = ilist(k,1)
      j = ilist(k,2)
      if( drsq(k).eq. 0.0 ) then
        print *,k
        stop 'zero distance between particles'
      endif
  end do
end subroutine check_zero_distance


subroutine check_zero_density(rho)
  ! Check that no particles have accidentaly been assigned zero density
  ! Zero densities can cause hard to locate divide by zero errors! 
  double precision, dimension(:) :: rho
  integer :: i
  do i = 1, size(rho)
    if(rho(i) .eq. 0) then
      print *,'density of', i, 'is zero'
      stop 'divide by zero error in density'
    endif
  end do
end subroutine check_zero_density



end module sphforce

