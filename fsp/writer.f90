!> Writes system state out to file.
!!
!! Naive ASCII format for bughunting and compatibility

module writer 

use global
use particle
use neighbour_list
use system_properties
use thermostat

private
public :: write_state, write_properties

contains

!> Writes the system state in a simple ASCII format, one line per particle.
!! \param p Particle
!! \param step Step number
!! \param 

!! Output format for a 2 dimensional system is as follows:
!! column 1 - particle x positions 
!! column 2 - y positions 
!! column 3 - x velocities
!! column 4 - y velocities
!! column 5 - x accelerations
!! column 6 - y acelerations
!! column 7 - temperatures
!! column 8 - internal energies
!! column 9 - masses
!! column 10 - densities
!! column 11 - pressures
!! column 12 - speeds of sound (not currently computed)
!! column 13 - viscosity
!! column 14 - smoothing lengths
!! column 15 - heat flux vector, x direction
!! column 16 - heat flux vector, y direction
!! column 17 - smoothed velocity, x
!! column 18 - smoothed velocity, y
!! column 19 - smoothed pressure
!! column 20 - smoothed temperature
!! column 21 - density gradient x
!! column 22 - density gradient y
!! column 23 - 0 for normal particle, 1 for boundary
!! column 24 - long smoothing lengths

subroutine write_state(p,step,bp)
  implicit none
  integer :: fp
  integer, intent(in) :: step
  type (particle_type), intent(in) :: p
  type (particle_type), intent(in), dimension(:), optional :: bp
  integer i,j
  
  character(len=8) :: kbchar
  fp=1
  if (replace) then
    write(kbchar,'(8A)') 'sphstate'
  else
    write(kbchar,'(i8.8)') step
  end if
  open(fp,file='sphstate.'//kbchar)
  
  do i=1,p%n
    if( p%dimn .eq. 1) then
      write (fp,'(9(E20.12))') p%x(i,1),p%v(i,1),p%a(i,1),p%temp(i), &
             &  p%u(i),p%m(i),p%rho(i),p%c(i),p%eta(i),p%sml(i)
      else if( p%dimn .eq. 2) then
        write (fp,'(22(E20.12),1i2,1(E20.12))')                            &
             &                p%x(i,1),p%x(i,2),p%v(i,1),p%v(i,2),         &
             &                p%a(i,1),p%a(i,2),p%temp(i),p%u(i),          &
             &                p%m(i),p%rho(i),p%p(i)+p%pco(i),p%c(i),      &
             &                p%eta(i),p%sml(i),p%q(i,1),p%q(i,2),         &
             &                p%v_flow(i,1),p%v_flow(i,2),p%p_smooth(i),   &
             &                p%t_smooth(i),p%grad_rho(i,1),               &
             &                p%grad_rho(i,2),1,p%sml_long(i)
      end if
  end do

  if(present(bp)) then
    do j=1,nbounds
      do i=1,bp(j)%n
        if( bp(j)%dimn .eq. 1) then
          write (fp,'(9(E20.12))') bp(j)%x(i,1),bp(j)%v(i,1),bp(j)%a(i,1) &
                 & ,bp(j)%temp(i)                                         &
                 & ,bp(j)%u(i),bp(j)%m(i),bp(j)%rho(i),bp(j)%c(i)         &
                 & ,bp(j)%eta(i),bp(j)%sml(i)
        else if( bp(j)%dimn .eq. 2) then
          write (fp,'(22(E20.12),1i2,1(E20.12))') bp(j)%x(i,1),bp(j)%x(i,2), &
                 &  bp(j)%v(i,1),                                            &
                 &  bp(j)%v(i,2),                                            &
                 &  bp(j)%a(i,1),bp(j)%a(i,2),bp(j)%temp(i),bp(j)%u(i),      &
                 &  bp(j)%m(i),bp(j)%rho(i),bp(j)%p(i)+bp(j)%pco(i),         &
                 &  bp(j)%c(i),                                              &
                 &  bp(j)%eta(i),bp(j)%sml(i),bp(j)%q(i,1),bp(j)%q(i,2),     &
                 &  bp(j)%v_flow(i,1),bp(j)%v_flow(i,2),bp(j)%p_smooth(i),   &
                 &  bp(j)%t_smooth(i),bp(j)%grad_rho(i,1),                   &
                 &  bp(j)%grad_rho(i,2),0,bp(j)%sml_long(i)
        end if
      end do
    end do
  endif

  close(fp); 
end subroutine write_state

subroutine write_properties(p,n)
  !< writes a line to the system property output file
  !! column 1 - kinetic energy
  !! column 2 - potential energy 
  !! column 3 - average temperature
  !! column 4 - average density
  !! column 5 - thermostat energy 
  !! column 6 - total energy with thermostat contribution subtracted 
  !! column 7 - timestep size (for this step)
  !! column 8 - system time elapsed for this step

  type(particle_type), intent(in) :: p
  integer, intent(in) :: n

  integer :: fp
  double precision :: av_temp,av_dens
  !ek = 0.0
  !ep=0.0
  !ep_iso=0.0
  av_temp=0.0
  av_dens=0.0
  fp =1
  !call system_kinetic_energy(p,ek)
  !call system_internal_energy(p,ep)
  !call isolated_internal_energy(ep,tstat%u_env,ep_iso)
  call average_scalar_property(p%temp,n,av_temp)
  call average_scalar_property(p%rho,n,av_dens)

  open(fp,file='properties.output',position='append')
  write (fp,'(8(E20.12))') ek,ep,av_temp,av_dens,u_env,isoham,dt,solver_sdt
  close(fp)

end subroutine write_properties


subroutine calc_radial_distance(x,l,w,side,radial)
  !< position, box x dimension, box y dimension, distribution side
  double precision, dimension(:) :: x
  double precision l,w,side,radial

  !should check dimensions
  radial = sqrt( (x(1)- (l-side)/2)**2 + ( x(2) - (w-side)/2)**2 )

end subroutine calc_radial_distance

subroutine write_neighbours(nl,step)
  implicit none
  type (nlist_data_type), intent(inout) :: nl
  integer, intent(in) :: step

  integer fp !file pointer
  integer i,j,ni,k

  character(len=8) :: kbchar
  write(kbchar,'(i8.8)') step
  open(fp,file='neighbour.'//kbchar)

  ni = nl%nip

  do k=1, ni
   i =nl%iindex(k)
   j =nl%jindex(k)
    write (fp,*) i,j
  enddo

  close(fp);

end subroutine write_neighbours


end module writer 
