module system_properties

use particle

contains

subroutine system_kinetic_energy(p,kinetic_energy,pb)
  ! calculates the total kinetic energy associated
  ! with the movement of the sph particles
  type (particle_type), intent(inout) :: p
  double precision, intent(inout) :: kinetic_energy
  type (particle_type), intent(inout), dimension(:), optional :: pb
  integer i,d,j
  double precision :: vsq

  kinetic_energy = 0.0
  do i=1,p%n
    vsq = 0.0
    do d=1,p%dimn
      vsq = vsq + p%v(i,d)**2
    end do 
    kinetic_energy = kinetic_energy + p%m(i)*vsq/2
  end do
   
  if(present(pb)) then
    do i=1,nbounds
      do j=1,pb(i)%n
        vsq = 0.0
        do d=1,p%dimn
          vsq = vsq + pb(i)%v(j,d)**2
        end do
        kinetic_energy = kinetic_energy + pb(i)%m(j)*vsq/2
      end do
    end do
  end if

end subroutine system_kinetic_energy

subroutine system_internal_energy(p,utot,pb)

  type (particle_type), intent(inout) :: p
  double precision, intent(inout) :: utot
  type (particle_type), intent(inout), dimension(:), optional :: pb
  integer i,j 

  utot = 0
  do i=1,p%n
    utot = utot + p%u(i)*p%m(i)
  end do

  if(present(pb))then
    do i=1,nbounds
      do j=1,pb(i)%n
        utot = utot + pb(i)%u(j)*pb(i)%m(j)
      end do
    end do
  end if

end subroutine system_internal_energy


subroutine total_internal_energy(u,n,u_tot,m)
  double precision, dimension(:) :: u,m
  integer :: n,i
  double precision :: u_tot
  ! remember u is specific internal energy 
 
  u_tot = 0
  do i=1,n
    u_tot = u_tot + u(i)*m(i)
  end do
end subroutine total_internal_energy


subroutine isolated_internal_energy(u_tot,u_env,u_iso)
  !calculated the internal energy
  !not associated with the thermostat
  double precision :: u_tot,u_env,u_iso
  u_iso = u_tot - u_env
end subroutine isolated_internal_energy


subroutine average_temperature(prop,n,av_prop)
  double precision, dimension(:) :: prop
  integer :: n
  double precision :: av_prop
  integer i
  
  av_prop = 0 
  do i=1,n
    av_prop = av_prop + prop(i)
  end do
  av_prop = av_prop/n  
end subroutine average_temperature

subroutine average_scalar_property(prop,n,av_prop)
  double precision, dimension(:) :: prop
  integer :: n
  double precision :: av_prop
  integer i
  
  av_prop = 0 
  do i=1,n
    av_prop = av_prop + prop(i)
  end do
  av_prop = av_prop/n  
end subroutine average_scalar_property

end module system_properties
