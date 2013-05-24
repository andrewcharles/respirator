! Unsurprisingly this is a thermostat. Uses scaling to keep
! the average temperature constant.
! Keeps track of the amount of energy added/removed by the thermostat.
module thermostat

use global
use particle
use system_properties
use eos

contains

subroutine apply_thermostat(n,temp,target_temp,u,rho,grad_rho,uenv,therm_type)
  ! Applies the thermostat rule to particle temperatures
  ! uenv -- amount of net system internal energy
  !         added by thermostat (not specific internal energy)
  ! grad_rho -- density gradient
  double precision, dimension(:) :: temp,u,rho
  double precision, dimension(:,:) :: grad_rho 
  integer therm_type
  integer n,i
  double precision :: uenv
  double precision :: tav, scalef, t_initial, u_initial,target_temp

  if(therm_type==1) then
    ! scaling thermostat
    ! calculate average temperature
    call average_temperature(temp,n,tav)

    ! calculate scaling factor
    scalef = target_temp/tav

    ! rescale temperature
    do i=1,n
      t_initial = temp(i)
      u_initial = u(i)
      temp(i) = scalef*temp(i)

      ! update energies to match new temperature
      call update_energy(u(i),temp(i),rho(i),grad_rho(i,:))
      uenv = uenv - (u(i) - u_initial)
    
    end do

  else if(therm_type==2) then
    ! Absolute thermostat
    do i=1,n
      t_initial = temp(i)
      u_initial = u(i)
      temp(i) = target_temp
      ! update energies to match new temperature
      call update_energy(u(i),temp(i),rho(i),grad_rho(i,:))
      uenv = uenv + (u_initial - u(i))
    end do
 
  end if

end subroutine apply_thermostat



end module thermostat
