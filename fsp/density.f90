!> Computes particle densities by either summation or continuity density
!! methods.
!! \detail Continuity method is not implemented yet.
!! \todo Implement continuity density
!! \todo ad-hoc stopping at density=2.0 should be configurable as this only
!! applies to vdw eos

module density

use global
use neighbour_list
use kernel
use particle

implicit none

contains

!> Calculates density and spatial gradient of density using
!! the summation density method
!! Here the grad_w being passed in is gradw_ji
subroutine sum_density(nl, n, sml, mass, w, rho, grad_rho, grad_w)
  implicit none

  type (nlist_data_type), intent(in) :: nl
  integer :: n
  double precision sml(:)
  double precision mass (:), rho(:)
  double precision, dimension(:) :: w
  double precision, dimension(:,:) :: grad_w
  double precision, dimension(:,:) :: grad_rho
  double precision :: rho_ini(n)
  integer i,j,k,d,q
  double precision selfdens, hv(dimn), r, wi(n)

  if(debug) print *,'calculating densities'
  hv=0.e0
  r=0.0
  rho_ini = rho

  do i=1,n
    do d=1,dimn
      grad_rho(i,d) = 0
    end do
  end do

  !integrate the density over space

  do i=1,n
     if(debug) print *,'calculating selfdensity for',i
    call calc_kernel(r,hv,sml(i),selfdens,hv,kernel_type,dimn)
    rho(i) = selfdens*mass(i)
    if(debug) then
      print *,'selfdensity',rho(i)
    endif
  end do

  !SPH sum for rho

  do k = 1, nl%nip
    i = nl%ilist(k,1)
    j = nl%ilist(k,2)

    rho(i) = rho(i) + mass(j)*w(k)
    rho(j) = rho(j) + mass(i)*w(k)
    
    !Density gradient
    do d=1,dimn
      ! This is negative because it 
      ! is assumed this is grad_w_ji
      grad_rho(i,d) = grad_rho(i,d) - mass(j) * grad_w(k,d)
      grad_rho(j,d) = grad_rho(j,d) + mass(i) * grad_w(k,d)
    end do
    if(debug) then
      print *,'contributions of',mass(j)*w(k), 'and', mass(i)*w(k), 'to density'
    endif
  end do

  if(adapt_smoothing .eq. 1) then
    call adapt_bandwidth(sml,rho,rho_ini)
  endif  
 
  !normalise
  !used for normalising see liu  p115

  if(normalise_kernel) then
    do k=1, nl%nip
      i = nl%ilist(k,1)
      j = nl%ilist(k,2)
      wi(i)=wi(i) + mass(j)/rho(j) * w(k)
      wi(j)=wi(j) + mass(i)/rho(i) * w(k)
    end do

    wi(i) = selfdens*mass(i)/rho(i)
    do i=1,n
      rho(i) = rho(i)/wi(i)
    end do
  endif

  do i=1,n
    if (rho(i) .gt. rho_max ) then
      if (unstable_density .eq. 0) then
        continue 
      else if (unstable_density .eq. 1) then
        print *,'Unstable density detected'
        print *,'density ',i,' is ',rho(i)
        print *, 'mass of', i, 'is', mass(i)
        if(debug) then
          print *, 'contributions of'
          do k=1, nl%nip
            q = nl%ilist(k,1)
            j = nl%ilist(k,2)
            if (q .eq. i) then
              print *,q,j, 'mass', mass(j), 'times smoothing', w(k)
            endif
          enddo
        endif
        stop 'unstable density'
      else if (unstable_density .eq. 2) then
        continue
      endif

    endif
  end do

  if(debug) then
     do i=1,n
       print *, i, 'density', rho(i)
       if (rho(i) .eq. 0.0) then
         print *, 'density of',i, 'is zero'
         print *, 'mass of', i, 'is', mass(i)
         print *, 'contributions of'
         do k=1, nl%nip
           q=nl%iindex(k)
           j=nl%jindex(k)
           if (q .eq. i) then
             print *,q,j,'mass', mass(j), 'times smoothing', w(k)
           endif
         enddo
         stop 'particle has zero density'
       endif
     end do
  end if

  if(verbose) print *,'densities calculated'

end subroutine sum_density

subroutine adapt_bandwidth(sml,rho,rho_ini)
  ! adapt the smoothing length given the
  ! smoothing lengths for all particles,
  ! densities, and pilot densities
  !f2py intent(in,out,overwrite) :: sml
  ! todo: it's not clear if this is even correct
  double precision, intent(inout) :: sml(:)
  double precision, intent(in) :: rho(:)
  double precision, intent(in) :: rho_ini(:)
  integer :: n,i
  n = size(rho)
   do i=1,n
     sml(i) = sml(i)*( rho_ini(i)/rho(i) )**(1.0/dimn)
   end do
end subroutine adapt_bandwidth


!> Adds the boundary particles contribution to 
!! particle density
subroutine boundary_density(nl,p,bp)
  type (nlist_data_type), intent(inout) :: nl
  type (particle_type), intent(inout) :: p
  type (particle_type), intent(inout) :: bp
  integer i,j,k,d
 
  do k = 1, nl%nip
    i = nl%ilist(k,1)
    j = nl%ilist(k,2)
    p%rho(i) = p%rho(i) + bp%m(j)*nl%w(k)
    !Density gradient
    do d=1,p%dimn
      p%grad_rho(i,d) = p%grad_rho(i,d) + bp%m(j) * nl%dwdx(k,d)
    end do
  end do

end subroutine boundary_density

subroutine continuity_density

end subroutine continuity_density

end module density
