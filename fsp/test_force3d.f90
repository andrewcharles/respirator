! We need a fortran program to do the testing of this ambitious
! subroutine.
! Having though about it, I am wondering if declaring the data
! it the module level is not more convenient for the wrapping.
! This means we can assign the data directly from python, and
! run the subroutine without copying?

program test_force3d
use sphforce3d

implicit none
! Declare the required variables and populate with some
! known data.
!ilist = 1,2 1,3 1,4 2,3 2,4 3,4
integer, dimension(6,2) :: ilist
!type (particle_type)    :: p
integer :: n,d,ni
!character(len=32) :: input_file      

!n = 10
!d = 3
!input_file = 'sph3d.txt'

!call create_particles(p,n,d)
!call initialise_particles(p)
!call rd_particle_input(p,input_file,1)
!call rd_nlist_input(nl,input_file,1)
!call init_nlist(nl, cutoff, p%n, dimn)
!call read_input(input_file)

double precision, dimension(4,3) :: x 
double precision, dimension(4,3) :: v
double precision, dimension(4,3) :: a
double precision, dimension(4,3,3) :: pi_irr
double precision, dimension(4,3,3) :: p_rev
double precision, dimension(4,3,3) :: p_rev_lr
double precision, dimension(4,3) :: q
double precision, dimension(4) :: u
double precision, dimension(4) :: m
double precision, dimension(4) :: rho
double precision, dimension(4) :: rho_lr 
double precision, dimension(4) :: c
double precision, dimension(4) :: eta
double precision, dimension(4) :: temp
double precision, dimension(4) :: dedt 
double precision, dimension(4) :: zeta
double precision, dimension(4,3) :: grad_rho
double precision, dimension(4,3,3) :: grad_v
double precision, dimension(4) :: sml
double precision, dimension(4) :: sml_long

double precision, dimension(6) :: dr 
double precision, dimension(6,3) :: dx 
double precision, dimension(6,3) :: dv
double precision, dimension(6) :: w
double precision, dimension(6,3) :: dwdx
double precision, dimension(6,3) :: dwdx_long

integer :: start_time,end_time,time_rate,i
real :: time_seconds

adash = 2.0
bdash = 0.5
kbdash = 1.0
n = 4
ni = 6
d = 3
x(:,1) = (/ 1.0, 2.0, 3.0, 4.0 /)
x(:,2) = (/ 0.0, 0.0, 0.0, 0.0 /)
x(:,3) = (/ 0.0, 0.0, 0.0, 0.0 /)
v(:,1) = (/ 0.0, 0.0, 0.0, 0.0 /)
v(:,2) = (/ 0.0, 0.0, 0.0, 0.0 /)
v(:,3) = (/ 0.0, 0.0, 0.0, 0.0 /)
m = (/ 1.0, 1.0, 1.0, 1.0 /)
a = 0
pi_irr = 0 
p_rev = 0
p_rev_lr = 0
m = 1.0
rho = 1.0
grad_rho = 0.0
rho_lr = 0.5
c = 0.0
eta = 1.0
zeta = 1.0
temp = 1.0
u = 0.1
cgrad = 1.0
ilist(:,1) = (/ 1, 1, 1, 2, 2, 3 /)
ilist(:,2) = (/ 2, 3, 4, 3, 4, 4 /)
sml = 1.0
sml_long = 2.0
dr = (/1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)
dx(:,1) = (/1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)
dx(:,2) = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
dx(:,3) = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

dv(:,1) = (/1.0, -1.0, 0.0, 2.1, 0.0, 0.0/)
dv(:,2) = (/0.0,  0.0, 0.0, 0.0, 0.0, 0.0/)
dv(:,3) = (/0.0,  0.0, 0.0, 0.0, 0.0, 0.0/)

w(:) = (/1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)
dwdx(:,1) = (/1.0, -1.0, 0.0, 2.1, 0.0, 0.0/)
dwdx(:,2) = (/0.0,  0.0, 0.0, 0.0, 0.0, 0.0/)
dwdx(:,3) = (/0.0,  0.0, 0.0, 0.0, 0.0, 0.0/)
dwdx_long = 0.5 * dwdx


call system_clock(count=start_time,count_rate=time_rate)
! Call the force subroutine

do i=1,1000
call calc_sphforce3d(ilist,x,v,a,p_rev,p_rev_lr,pi_irr,      &
  &   grad_rho,grad_v,u,dedt,m,rho,rho_lr,temp,q,     &
  &   c,eta,zeta, &
!  &   adash,bdash,kbdash,
  &   dv,dr,dx,  &
  &   sml,sml_long,w,dwdx,dwdx_long,.True.,n,ni)
enddo

call system_clock(count=end_time)
time_seconds = (end_time - start_time)/real(time_rate)
print *,'Call time (seconds)', time_seconds

! Inspect the output. Not concerned about numbers at this stage
! so much as working.
! print *,a
! print *,dedt
! Profit. 

end program test_force3d
