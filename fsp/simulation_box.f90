!***********************************************************************
   module simulation_box
!***********************************************************************
!  Module for simulation box procedures and types.
!
!  Note: box vectors are stored with first index for vector number and
!  second index for component number. (Three row vectors.)
!
!  Note: boxvec0 is the initial box - used by elongational and bulk
!  flow simulations.
!
!  Note: allocatables in defined types requires 
!
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
use global
use reader

implicit none

private

public :: apply_pbcs, create_box, destroy_box, &
      & init_box, minimum_image, rd_box_input, rd_box_state, &
      & store_box_state, volume, wr_box 

public :: box_data_type

type box_data_type
  double precision, dimension(:,:), allocatable :: boxvec, boxvec0
  double precision :: volume, angle
  integer :: ndim, origin  ! (origin: 0=bottom left, 1=centre)
  integer :: state  ! current box state
                 ! (0=initial, 1=sheared, 2=rotated and elongated)
  integer :: pbctype  ! periodic boundary conditions
                 ! (0=all directions, 
                 !  1=all except y direction
                 !  2=none )
end type box_data_type

contains

!> Apply pbcs
subroutine apply_pbcs(box,r)
  implicit none
  type (box_data_type) :: box
  double precision, dimension(:,:) :: r
  double precision, dimension(size(r,1),box%ndim) :: rsh
  double precision, dimension(size(r,1),box%ndim) :: rb, shu
  double precision :: r0
  integer :: i

  if(box%pbctype .eq. 2) return
  select case(box%origin)
  case(0)
   r0 = 0.5
  case(1)
   r0 = 0.0
  end select

  rb = rboxv(r,box)

  do i = 1, box%ndim
   rb(:,i) = rb(:,i) - r0
  enddo

  shu = int(2.0*(rb) - int(rb))
  if(box%pbctype .eq. 1) shu(:,2) = 0.d0

  do i = 1, size(r,2)
   rsh(:,i) = shu(:,1)*box%boxvec(1,i) + shu(:,2)*box%boxvec(2,i)
   if(box%ndim .eq. 3) rsh(:,i) = rsh(:,i) + shu(:,3)*box%boxvec(3,i)
   r(:,i) = r(:,i) - rsh(:,i)
  enddo

end subroutine apply_pbcs

subroutine check_box(box, rcut)
  implicit none
  type (box_data_type) :: box
  double precision :: rcut

  if(2.0*rcut .ge. maxval(box%boxvec)) stop 'box is less than 2*rcut'

end subroutine check_box


!> Allocates all dynamic memory the box needs
subroutine create_box(box,ndim)
  implicit none
  type (box_data_type) :: box
  integer :: error,ndim
  box%ndim = ndim

  ! allocate box vector array
  allocate(box%boxvec(box%ndim,box%ndim), stat=error)
  allocate(box%boxvec0(box%ndim,box%ndim), stat=error)
  if(error .ne. 0) stop 'allocation error in create_box'

end subroutine create_box


!> deallocates all memory for box
subroutine destroy_box(box)
  implicit none
  type (box_data_type) :: box
  integer :: error

  if(allocated(box%boxvec)) deallocate(box%boxvec, box%boxvec0, stat=error)

end subroutine destroy_box

!> Initialise box variables
subroutine init_box(box, lx, ly, lz)
  implicit none
  type (box_data_type) :: box
  double precision :: lx, ly
  double precision, optional :: lz

  ! initalize box data
  box%boxvec = 0.0
  box%angle = 0.0

  box%boxvec(1,1) = lx
  box%boxvec(2,2) = ly
  if(present(lz)) box%boxvec(3,3) = lz
  box%volume = volume(box%boxvec)

end subroutine init_box

!> works for 2 or 3-dimensions
subroutine minimum_image(box, rij)

  implicit none
  type (box_data_type) :: box
  double precision, dimension(:,:), intent(inout) :: rij

  double precision, dimension(size(rij,1),box%ndim) :: shu
  double precision, dimension(size(rij,1),box%ndim) :: rb
  
  !double precision, dimension(nl%nab,box%ndim) :: shu
  !double precision, dimension(nl%nab,box%ndim) :: rb

  !! Q. Might be worth explicitly only working on up to nl%nab members?
  !! A. This will not work because sometimes we take the min image
  !! of the full npairs, other times of the compressed list
  !! for now just do the unneccesay computation, and consider passing
  !! a size paramter

  integer :: i

  if(verbose) print *,'Computing minimum image, calling box vectors'
  if (box%pbctype .eq. 2) return

  rb = rboxv(rij,box)

  if(debug) print *,'Minimum image checkpoint 2'

  ! this integer cast may be expensive
  shu = int(2.0*rb - int(rb))

  select case(box%ndim)
  case(2)
    if(box%pbctype .eq. 0) then
  !    do d = 1, size(rij,2) ! pbcs in all directions
  !     rij(:,d) = rij(:,d) - shu(:,1)*box%boxvec(1,d) - shu(:,2)*box%boxvec(2,d)
     rij(:,1) = rij(:,1) - shu(:,1)*box%boxvec(1,1) - shu(:,2)*box%boxvec(2,1)
     rij(:,2) = rij(:,2) - shu(:,1)*box%boxvec(1,2) - shu(:,2)*box%boxvec(2,2)
  !    enddo
    else if(box%pbctype .eq. 1) then ! no pbcs in y direction 
     rij(:,1) = rij(:,1) - shu(:,1)*box%boxvec(1,1)
    else if(box%pbctype .eq. 2) then ! no pbcs
      return
    end if

  case(3)
    if(box%pbctype .eq. 0) then
      do i = 1, size(rij,2) ! pbcs in all directions
       rij(:,i) = rij(:,i) - shu(:,1)*box%boxvec(1,i) - shu(:,2)*box%boxvec(2,i) &
                         & - shu(:,3)*box%boxvec(3,i)
      enddo
    else if(box%pbctype .eq. 1) then
      !! no pbcs in y direction - rectangular box assumed
      rij(:,1) = rij(:,1) - shu(:,1)*box%boxvec(1,1) - shu(:,3)*box%boxvec(3,1)
      rij(:,3) = rij(:,3) - shu(:,1)*box%boxvec(1,3) - shu(:,3)*box%boxvec(3,3)
    else if(box%pbctype .eq. 2) then ! no pbcs
      return
    end if

  case default
     stop 'invalid ndim in minimum image'
  end select

  if(verbose) print *,'completed minimum image'

  end subroutine minimum_image

  !***********************************************************************
  subroutine rd_box_input(box,ifile,iunit)
  !***********************************************************************
  ! Reads box variables from input file.
  !
  !***********************************************************************
  implicit none

  type(box_data_type) :: box
  integer :: iunit
  character(len=32) :: ifile         ! input file name
  open(iunit,file=ifile,status='old',form='formatted', &
     & position='rewind')
  call read_int(iunit,'NDIM',box%ndim)
  call read_int(iunit,'ORIGIN',box%origin)
  call read_int(iunit,'PBCTYPE',box%pbctype)
  close(iunit)

  end subroutine rd_box_input
  !***********************************************************************
  subroutine rd_box_state(box,iunit)
  !***********************************************************************
  implicit none
  type (box_data_type), intent(out) :: box
  integer, intent(in) :: iunit

  read(iunit) box%ndim, box%origin, box%boxvec, box%boxvec0, &
          & box%state, box%angle

  box%volume = volume(box%boxvec)

  end subroutine rd_box_state
  !***********************************************************************
  subroutine store_box_state(box,iunit)
  !***********************************************************************
  implicit none
  type (box_data_type), intent(in) :: box
  integer, intent(in) :: iunit

  write(iunit) box%ndim, box%origin, box%boxvec, box%boxvec0, &
           & box%state, box%angle

  end subroutine store_box_state
  !***********************************************************************
  function volume(boxvec)
  ! Evaluate triple product of box vectors to find box volume.
  !***********************************************************************
  implicit none
  double precision, dimension(:,:) :: boxvec
  double precision :: volume
  integer :: ndim

  ndim = size(boxvec,1)

  select case(ndim)
  case(2)
  volume =  boxvec(1,1)*boxvec(2,2) &
         & -boxvec(1,2)*boxvec(2,1)
  case(3)
  volume =  boxvec(1,1)*(boxvec(2,2)*boxvec(3,3) - &
                         & boxvec(2,3)*boxvec(3,2))  &
         & -boxvec(1,2)*(boxvec(2,1)*boxvec(3,3) - &
                         & boxvec(2,3)*boxvec(3,1))  &
         & +boxvec(1,3)*(boxvec(2,1)*boxvec(3,2) - &
                         & boxvec(2,2)*boxvec(3,1))
  end select

  end function volume

  function rboxv(r,box)
    ! Evaluate position in terms of box vectors.
    implicit none
    double precision, dimension(:,:) :: r
    type (box_data_type) :: box
    double precision, dimension(size(r,1),size(r,2)) :: rboxv
    integer :: ndim

    if(verbose) print *,'rboxv start'
    ndim = size(r,2)

    select case(ndim)
    case(3)
    rboxv(:,1) =  (r(:,1)*box%boxvec(2,2)*box%boxvec(3,3) &
       & - r(:,2)*box%boxvec(2,1)*box%boxvec(3,3))/box%volume
    rboxv(:,2) = (-r(:,1)*box%boxvec(1,2)*box%boxvec(3,3) &
       & + r(:,2)*box%boxvec(1,1)*box%boxvec(3,3))/box%volume
    rboxv(:,3) = r(:,3)/box%boxvec(3,3)

    case(2)
    if(verbose) print *,'rboxv: Two dimensions.'

    rboxv(:,1) =  (r(:,1)*box%boxvec(2,2) &
       & - r(:,2)*box%boxvec(2,1))/box%volume
    rboxv(:,2) = (-r(:,1)*box%boxvec(1,2) &
       & + r(:,2)*box%boxvec(1,1))/box%volume

    case(1)
        rboxv(:,1) =  (r(:,1)/box%boxvec(1,1))

    case default
    stop 'invalid dimension in function rboxv'
    end select

    if(verbose) print *,'rboxv completed'

  end function rboxv

  !***********************************************************************
  subroutine wr_box(box,ofile,ounit)
  ! writes box data to the output file
  !***********************************************************************
  implicit none
  type(box_data_type) :: box
  integer :: ounit
  character(len=32) :: ofile      !output file name

  open(ounit,file=ofile,status='old',form='formatted',position='append')
  write (ounit,*) 'Box data'
  write (ounit,*)
  write (ounit,*) 'Dimensionality            = ', box%ndim
  write (ounit,*) 'Box origin type           = ', box%origin
  write (ounit,*) 'Box vectors               = ', box%boxvec
  write (ounit,*) 'Box angle                 = ', box%angle
  write (ounit,*)
  close(ounit)

  end subroutine wr_box
  !***********************************************************************






  end module simulation_box
