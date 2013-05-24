!> Unused subroutine - vestige of previous implementation of the density
!! gradient contribution to the pressure tensor.
subroutine vdw_gradient(rho,u,p,grad_rho,cgrad)
  double precision :: u,t,rho
  double precision, dimension(:) :: grad_rho
  double precision, intent(out) :: p
  double precision :: cgrad
  ! calculate the pressure
  ! if(rho >= 2) then 
  !   if(unstable_density .eq. 0) then
  !     continue
  !   else if(unstable_density .eq. 1) then
  !    stop 'unstable density' 
  !   else if(unstable_density .eq. 2) then
  !     rho = 1.99
  !   endif
  ! endif
  if(rho == 0.0) stop 'density is zero'
  t = (u+adash*rho + cgrad*(grad_rho(1)**2)/rho + cgrad*(grad_rho(2)**2)/rho )/kbdash
  !t = (u+adash*rho  )/kbdash
  p = (rho*kbdash*t)/(1-rho*bdash) - adash*rho*rho - cgrad*(grad_rho(1)**2) - cgrad*(grad_rho(2)**2)
end subroutine vdw_gradient


!> Unused subroutine - vestige of previous implementation of the density
!! gradient contribution to the pressure tensor.
subroutine vdw_gradient_energy(u,t,rho,grad_rho)
  double precision :: u,t,rho
  double precision, dimension(:) :: grad_rho
  double precision :: cgrad
  cgrad = 1.0
  ! calculate the internal energy
  ! this is an experimental equation
  
  u = t*kbdash - adash*rho - cgrad*(grad_rho(1)**2)/rho - cgrad*(grad_rho(2)**2)/rho
  
end subroutine vdw_gradient_energy


!> Unused subroutine - vestige of previous implementation of the density
!! gradient contribution to the pressure tensor.
subroutine vdw_gradient_temp(u,t,rho,grad_rho)
  double precision :: u,t,rho
  double precision, dimension(:) :: grad_rho 
  double precision :: cgrad
  cgrad = 1.0
  ! calculate the temperature
  ! if(verbose) print *,size(grad_rho)
  t = (u+adash*rho + cgrad*(grad_rho(1)**2)/rho + cgrad*(grad_rho(2)**2)/rho ) /kbdash
  !omit the gradient term contribution to energy and temperature
  !t = (u+adash*rho ) /kbdash
end subroutine vdw_gradient_temp
