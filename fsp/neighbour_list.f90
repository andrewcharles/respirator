!> Creates a neighbour list for 2-body forces.
!! based on Wmd2 nlist by Peter Daivis.
!! Further development by Andrew Charles
!! The usual procedure for using the module's verlet list is as follows:
!! 1. form iindex, jindex by brute force or cell
!! 2. form nlist by masking by tolerance shell
!! 3. form ilist by masking by cutoff radius

! COPYRIGHT NOTICE
! Copyright RMIT University, Peter Daivis and Andrew Charles, 2008.
! This code is made available under the Computer Physics Communication License.

module neighbour_list
use simulation_box
use global
use reader

implicit none

private

public :: compress_nlist, destroy_nlist, find_pair_separations, &
        & form_nlist, increment_nlist_drtot, init_nlist, &
        & rd_nlist_input, reform_nlist_now,  &
        & calc_dr,validate_nlist

public :: nlist_data_type

type :: nlist_data_type
    
  ! shell -- gives us the thickness around the cutoff radius
  !          to allow when building the nlist neighbour list
  ! search_rsq -- radius squared around each particle to include in the
  !               nlist list
  ! dsqmax -- used for deciding whether to rebuild the list
  ! drtot -- for rebuild decision
  ! drini -- for rebuild decision
  ! old_rijsq -- for rebuild decision. Set to the inital displacement 
  !             between pairs 

  double precision :: shell 
  double precision :: search_rsq 
  double precision :: dsqmax 
  double precision, dimension(:,:), allocatable :: drtot
  double precision, dimension(:,:), allocatable ::drini 
  double precision, dimension(:), allocatable :: old_rijsq 

  ! when compress_list is called, ilist contains the compressed list
  ! compression in this case is just removing pairs outside of cutoff
  ! nlist -- neighbour list with all pairs to be searched
  ! ilist -- neighbour list with non-interacting pairs removed
  ! iindex -- first indices of uncompressed list 
  ! jindex -- second indices of uncompressed list 

  integer, dimension(:,:), allocatable :: nlist
  integer, dimension(:,:), allocatable :: ilist
  integer, dimension(:), allocatable :: iindex
  integer, dimension(:), allocatable :: jindex 

  !> imask is used to prevent pairs from interacting
  logical, dimension(:,:), allocatable :: imask

  integer :: np !< number of particles
  integer :: np2 !< number of particles in second set of particles
  integer :: max_interact !< maximum interactions in this list
  integer :: npairs !< number of neighbours in the raw brute or cell list (iindex)
  integer :: nab !< number of neighbours after search radius exclusion (nlist)
  integer :: nip !< number of interacting pairs after cutoff exclusion (ilist)
  integer :: mnp !< maximum neighbours per particle   

  !> nmask is used to screen pairs that are outside of interaction range 
  !! both nmask and cmask are allocated the maximum number of pairs
  !! usage (truncation) is given by nab for nmask and nip for cmask
  logical, dimension(:), allocatable :: nmask
  !> cmask is used to mask pairs outside cutoff range
  logical, dimension(:), allocatable :: cmask
 
  ! Pair properties
  ! These are always given wrt to ilist, i.e. the cutoff
  ! truncated list, with the exception of drij and rijsq which
  ! are temporarily calculated at each level

  ! dv -- velocity difference
  ! rij -- distance 
  ! rijsq -- distance squared
  ! w -- kernel
  ! w_long -- long range kernel
  ! drij -- pair displacement r_j - r_i
  ! dwdx -- kernel gradient at position of j with respect to particle i
  ! Note that this is grad_j W_ij, so when used in the sph equations it
  ! needs to be flipped.
  ! dwdx_long -- long kernel gradient
  double precision, dimension(:,:), allocatable :: dv 
  double precision, dimension(:), allocatable :: rij 
  double precision, dimension(:), allocatable :: rijsq 
  double precision, dimension(:), allocatable :: w 
  double precision, dimension(:), allocatable :: w_long 
  double precision, dimension(:,:), allocatable :: drij 
  double precision, dimension(:,:), allocatable :: dwdx 
  double precision, dimension(:,:), allocatable :: dwdx_long 

end type nlist_data_type

contains

!***********************************************************************
 subroutine init_nlist(nl, rcut, np, ndim, np2)
!***********************************************************************
  ! np2 gives the number of particles in system 2
  type (nlist_data_type), intent(inout) :: nl
  double precision :: rcut
  integer :: np, ndim
  integer, optional :: np2 
  integer :: error
  !double precision :: factor

  if(verbose) print *,'Initialising nlist'

  ! do some initializations
  if( present(np2)) then
    nl%np2 = np2
  endif
  nl%np = np
  nl%search_rsq = (nl%shell+rcut)**2
 
  ! Set max interactions to brute force maximum 
  if( present(np2)) then
    max_interact = np*np2
    nl%max_interact = np*np2
  else
    max_interact = ((np*(np-1))/2)
    nl%max_interact = ((np*(np-1))/2)
  endif

  if(debug) print *,'Max interactions',max_interact

  nl%dsqmax = 2.0*nl%shell*(rcut + 0.5*nl%shell)

  allocate(nl%w(nl%max_interact), stat=error)
  allocate(nl%dwdx(nl%max_interact,ndim), stat=error)
  allocate(nl%w_long(nl%max_interact), stat=error)
  allocate(nl%dwdx_long(nl%max_interact,ndim), stat=error)
  allocate(nl%rij(nl%max_interact), stat=error)
  allocate(nl%drij(nl%max_interact,ndim) , stat=error)
  allocate(nl%dv(nl%max_interact,ndim), stat=error)
  allocate(nl%rijsq(nl%max_interact), stat=error)
  allocate(nl%iindex(nl%max_interact),stat=error)
  allocate(nl%jindex(nl%max_interact), stat=error)
  allocate(nl%nmask(nl%max_interact), stat=error)
  allocate(nl%cmask(nl%max_interact), stat=error)
  allocate(nl%nlist(nl%max_interact,2), stat=error)
  allocate(nl%ilist(nl%max_interact,2), stat=error)
      
  if(error .ne. 0) stop 'allocation error in main neighbour list'

  if(debug) print *,'Allocations allocated'
  if (present(np2)) then
    allocate(nl%imask(np,np2), stat = error)
  else
    if(debug) print *,'imask setup'
    allocate(nl%imask(np,np), stat = error)
    if(debug) print *,'imask allocated'
  end if
  nl%imask(:,:) = .true.
  !nl%imask(:,0) = .false.
  !nl%imask(0,:) = .false.
  if(debug) print *,'imask setup complete'
  if(verbose) print *,'Finished initialising nlist'

end subroutine init_nlist

subroutine form_nlist(nl, box, r, rb)
  ! builds the neighbour list with range searchrij
  ! nl -- neighbour list
  ! box -- simulation box
  ! r -- particle positions
  ! rb -- postions of second set of particles (for ab interactions)
  implicit none

  type (nlist_data_type), intent(inout) :: nl
  type (box_data_type), intent(inout) :: box
  double precision, dimension(:,:) :: r
  double precision, optional, dimension(:,:) :: rb 
  double precision :: drsq

  integer :: i, j, error, ndim
  integer k

  ndim = size(r,2)
 
  nl%npairs = nl%max_interact 
  nl%nab = nl%max_interact
  nl%nip = nl%max_interact

  ! form iindex and jindex arrays
  if(verbose) print *,'Starting neighbour list formation ...'

  ! all-particle neighbour list formation
  if(verbose) print *,'form_nlist: calling all-pairs method'
  if (present(rb)) then
    call all_pairs_ab(nl)
  else 
    call all_pairs_method(nl)
  endif

  if (present(rb)) then
    do k=1, nl%npairs
      nl%drij(k,:) = rb(nl%jindex(k),:) - r(nl%iindex(k),:)
    end do
  else
    do k=1, nl%npairs
      nl%drij(k,:) = r(nl%jindex(k),:) - r(nl%iindex(k),:)
    end do
  endif

  ! diagnostics
  ! at this point the neighbour lists will contain
  ! all pairs
  if(debug) then 
    print *,'number of pairs:',nl%nip 
      do i=1, nl%nip
        print *,'pair data - dr,indices,pos',i
        print *,nl%drij(i,1),nl%drij(i,2)
        print *,nl%iindex(i),nl%jindex(i)
        print *,r(nl%iindex(i),:)
        print *,r(nl%jindex(i),:)
      end do
  endif
     
  if(debug) print *,'form_nlist: calling minimum image'

  ! find minimum image separations
  call minimum_image(box, nl%drij)

  if(verbose) print *,'form_nlist: minimum image complete'

  ! eliminate pairs outside the tolerance shell and find total number of
  ! pairs in list
  
  if(verbose) then
    print *,'form_nlist: nmask is', size(nl%iindex)
    !print *,size(nl%drij)
    !print *,size(nl%nmask)
    !print *,nl%drij(12,2)
    !print *,nl%drij(12,1)
    !print *,nl%nip
    !print *,max_interact
  endif

  do k=1, nl%npairs
    drsq = nl%drij(k,1)**2 + nl%drij(k,2)**2
    if( drsq .eq. 0.0 ) then
      print *,nl%iindex(k),nl%jindex(k),k
      stop 'zero distance in nlist'
    endif 
    nl%nmask(k) = drsq .le. nl%search_rsq
  end do
  
  if(verbose) print *,'form_nlist: creating list'
 
  ! nab gives the number in the tolerance shell
  ! nip can never be greater than this
  nl%nab = count(nl%nmask)
  nl%nip = count(nl%nmask)
            
  if(nl%nab .eq. 0) then
    if(verbose) print *,'No particles in neighbour list - increase tolerance'
  endif

  ! create neighbour list
  j=1
  do i= 1, nl%npairs
    if(nl%nmask(i)) then
      nl%nlist(j,1)=nl%iindex(i)
      nl%nlist(j,2)=nl%jindex(i)
      nl%drij(j,1) = nl%drij(i,1)
      nl%drij(j,2) = nl%drij(i,2)
      nl%rij(j)=nl%rij(i)
      nl%rijsq(j)=nl%rijsq(i)
      if(debug) then
        if (nl%iindex(i) .eq. 0) stop 'zero index in formnlist - i'
        if (nl%jindex(i) .eq. 0) stop 'zero index in formnlist - j'
      endif
      j=j+1
    end if
  end do
  if (debug) print *,nl%drij(1,1)
  
  if(verbose) print *,'form_nlist: allocating dr totals'
        
    if(allocated(nl%drtot)) deallocate(nl%drtot)
    if(allocated(nl%drini)) deallocate(nl%drini)
    if(allocated(nl%old_rijsq)) deallocate(nl%old_rijsq)
    allocate(nl%drtot(nl%nab,ndim), stat=error)
    allocate(nl%drini(nl%nab,ndim), stat=error)
    allocate(nl%old_rijsq(nl%nab), stat=error)
    
    if(error .ne. 0) stop 'allocation error in form_nlist'
   
    ! set drini to the initial seperation between
    ! pairs. 
 
    do i = 1, nl%nab
      nl%drtot(i,:) = nl%drij(i,:)
      nl%drini(i,:) = nl%drij(i,:)
    enddo
    
    if(verbose) print *,'forming nlist: setting old_rijsq'
        
    do i= 1, nl%nab
      nl%old_rijsq(i) = ( nl%drtot(i,1)**2 +nl%drtot(i,2)**2 ) 
    end do
  
  if (verbose) print *,'form_nlist: complete :-)'
  if (debug) print *,nl%drij(i,1)

end subroutine form_nlist

!> Generates an all-pairs neighbour list in nl%iindex and nl%jindex
subroutine all_pairs_method(nl)
  implicit none
  type (nlist_data_type), intent(inout) :: nl
  integer :: i, j, np, ip, ni
  
  if(debug) print *, 'all pairs'
  
  np = nl%np
  ni = (np*(np-1)/2)
 
  ip = 0
  do i=1, np
    do j=i+1,np
      if(nl%imask(i,j)) then
        ip = ip + 1
        nl%iindex(ip) = i 
        nl%jindex(ip) = j 
      endif
      if(debug) then
        if (i .eq. j)  stop 'particle selfinteration'
        if (i .eq. 0) stop 'zero index'
        if (j .eq. 0) stop 'zero index'
      endif
    enddo
  enddo
  nl%npairs = ip

end subroutine all_pairs_method


subroutine all_pairs_ab(nl)
  !< all pairs for different particle systems
  implicit none
  type (nlist_data_type), intent(inout) :: nl
  integer :: i, j, ip
  
  if(debug) print *, 'all pairs ab'
 
  ip = 0
  do i=1, nl%np
    do j=1,nl%np2
      !if(nl%imask(i,j)) then
        ip = ip + 1
        nl%iindex(ip) = i 
        nl%jindex(ip) = j 
      !endif
    enddo
  enddo
  nl%npairs = ip

end subroutine all_pairs_ab

!> Computes the seperation between pairs in the neighbour list.
subroutine find_pair_separations(nl,box,r,rb)
  implicit none
  type (nlist_data_type), intent(inout) :: nl
  type (box_data_type), intent(inout) :: box
  !> positions of particles
  double precision, dimension(:,:) :: r
  double precision, optional, dimension(:,:) :: rb
  integer i,k

  ! find particle separations
  if(verbose) print *, 'subroutine find_pair_separations start'
  if(debug) print *,'rij is',size(nl%drij(:,1)),size(nl%drij,2)
  if(debug) print *,'pairs is',size(nl%nlist(:,1))
  if(debug) print *,nl%drij
  if(debug) print *,nl%drij(1,1)
 
  if (present(rb)) then
    do k=1,nl%nab
      nl%drij(k,:) = rb(nl%nlist(k,2),:) - r(nl%nlist(k,1),:)
    end do
  else
    do k=1,nl%nab
      nl%drij(k,:) = r(nl%nlist(k,2),:) - r(nl%nlist(k,1),:)
    end do
  end if

  if (verbose) print *, 'calling minimum_image '
  call minimum_image(box, nl%drij)

  if (debug) print *, 'pairsep - minimum image call successful'
  
  ! find interacting pairs
  do k=1,nl%nab
    nl%rijsq(k) = nl%drij(k,1)*nl%drij(k,1) + nl%drij(k,2)*nl%drij(k,2)

    ! Error checking
    if(nl%nlist(k,1) .eq. 0) then
      print *,'index',k,'is zero'
      print *,nl%drij(k,1),nl%drij(k,2)
      stop 'an index cannot be zero!'
    end if
    if(nl%rijsq(k) .eq. 0.0) then
      print *,'zero distance b/w pairs detected in find_pair_separation - crash'
      print *,'pair number',k
      print *,'indices',nl%nlist(k,1),nl%nlist(k,2)
      print *,nl%drij(k,1), nl%drij(k,2)
      stop 'zero distance'
    end if
  end do
  
  if (debug) print *, 'Pair-sep tests passed ok'

  nl%cmask = .false.

  do i=1, nl%nab
    nl%cmask(i) = nl%rijsq(i) .le. cutoff*cutoff
    if(nl%cmask(i)) then
        nl%rij(i) = sqrt(nl%rijsq(i))
    end if
  end do
  nl%nip = count(nl%cmask)

end subroutine find_pair_separations


!> Removes pairs outside of interaction range for nl%iilist.
!! this results in two nlists: nl%nlist which has all pairs in
!! the shell radius
!! nl%ilist, which has all pairs in interaction range 
subroutine compress_nlist(nl)
  !,dx,dr,drsq)

  type (nlist_data_type), intent(inout) :: nl
  !double precision, dimension(:,:) :: dx
  !double precision, dimension(:) :: drsq
  !double precision, dimension(:) :: dr
  integer :: i, j
  
  !dx = nl%dx
  !drsq = nl%drsq
  !dr = nl%dr


  if (debug) print *, 'compress nlist start'
  nl%nab = count(nl%nmask)
  nl%nip = count(nl%cmask)
  if (verbose) print *,'Compress nlist :cmask',count(nl%cmask) &
  & ,'nmask',count(nl%nmask)
  
  ! create compressed interaction list
  j=1
  do i=1,nl%nab
    if(nl%cmask(i)) then
      nl%ilist(j,1) = nl%nlist(i,1)
      nl%ilist(j,2) = nl%nlist(i,2)
      if (debug) print *,'set ilist',j,'from',i
      ! print *,'set ilist',j,'from',i,nl%nlist(i,1),nl%nlist(i,2)
      j=j+1
    end if
  end do

  if(validate) call validate_nlist(nl)

  if(debug) print *,'compress_nlist: interactions =',nl%nip
  if(debug) print *,'compress_nlist: mask size =',size(nl%cmask)
  if(debug) print *,'compress_nlist: mask count =',count(nl%cmask)
  if(debug) print *,'compress_nlist: ab count =',nl%nab
  
  ! repack the pair seperations
  ! and all the other neighbourly properties
   
  j=0
  do i=1,nl%nab
    if(nl%cmask(i)) then
      j=j+1
      nl%drij(j,1) = nl%drij(i,1)
      nl%drij(j,2) = nl%drij(i,2)
      nl%rij(j)=nl%rij(i)
      nl%rijsq(j)=nl%rijsq(i)
      if (debug) print *,'compressing...'
      if (debug) print *,size(nl%drij,1)
      if (debug) print *,size(nl%drij,2)
      if (debug) print *,i
      if (debug) print *,nl%drij(i,1)

      if (nl%rij(j) == 0.0)  then
        if (debug) print *,'zero compressing...'
          print *,'zero distance detected - crashing'   
          print *,nl%rij(j),i,j
          print *,nl%rij(i)
          stop 'zero distance in nlist'
      endif
      if (debug) print *,'compressing...'
      if (debug) print *,nl%drij(1,1)
      if (debug) print *,'set ilist props',j,'from',i,nl%drij(i,:)
    end if
  end do

  if (nl%nip .ne. j) then
    print *,nl%nip,j
    stop 'wtf'
    endif

  if (validate) call validate_nlist(nl)
  if (debug) print *, 'compress nlist end'

end subroutine compress_nlist
    
!> Calculates distances between neighbours using ilist.
subroutine calc_dr(nl,r,box,rb)
  type (nlist_data_type), intent(inout) :: nl
  type (box_data_type) :: box
  double precision, dimension(:,:) :: r 
  double precision, optional, dimension(:,:) :: rb 
  integer k,i,j,d

  if (present(rb)) then
    do k = 1, nl%nip
        i = nl%ilist(k,1)
        j = nl%ilist(k,2)
        nl%drij(k,:) = rb(j,:) - r(i,:)
        nl%rijsq(k) = (nl%drij(k,1)*nl%drij(k,1) + nl%drij(k,2)*nl%drij(k,2))
        nl%rij(k) = sqrt(nl%rijsq(k))
    end do
  else
    do k = 1, nl%nip
        i = nl%ilist(k,1)
        j = nl%ilist(k,2)
        do d=1,dimn
          nl%drij(k,d) = r(j,d) - r(i,d)
        enddo
        nl%rijsq(k) = (nl%drij(k,1)*nl%drij(k,1) + nl%drij(k,2)*nl%drij(k,2))
        nl%rij(k) = sqrt(nl%rijsq(k))
    end do
  end if

    !! \todo it's not entirely clear why we get the sqrt twice......
    
  call minimum_image(box,nl%drij)

end subroutine calc_dr 



subroutine destroy_nlist(nl)
  type (nlist_data_type), intent(inout) :: nl

! destroy neighbour list
  if(allocated(nl%drtot)) deallocate(nl%drtot)
  if(allocated(nl%old_rijsq)) deallocate(nl%old_rijsq)
  if(allocated(nl%imask)) deallocate(nl%imask)
  if(allocated(nl%nlist)) deallocate(nl%nlist)
  if(allocated(nl%ilist)) deallocate(nl%ilist)
  if(allocated(nl%iindex)) deallocate(nl%iindex)
  if(allocated(nl%jindex)) deallocate(nl%jindex)
  if(allocated(nl%w)) deallocate(nl%w)
  if(allocated(nl%w_long)) deallocate(nl%w_long)
  if(allocated(nl%dwdx)) deallocate(nl%dwdx)
  if(allocated(nl%dwdx_long)) deallocate(nl%dwdx_long)
  if(allocated(nl%rij)) deallocate(nl%rij)
  if(allocated(nl%rijsq)) deallocate(nl%rijsq)
  if(allocated(nl%drij)) deallocate(nl%drij)
  if(allocated(nl%dv)) deallocate(nl%dv)
  if(allocated(nl%iindex)) deallocate(nl%iindex)
  if(allocated(nl%jindex)) deallocate(nl%jindex)
  if(allocated(nl%nmask)) deallocate(nl%nmask)
  if(allocated(nl%cmask)) deallocate(nl%cmask)
  
end subroutine destroy_nlist

subroutine increment_nlist_drtot(nl,dr)
  type (nlist_data_type), intent(inout) :: nl
  !> dr is an n-particles by dim array
  !! the step subroutines fill it and pass it.
  double precision, dimension(:,:) :: dr
  integer i

  if(debug) then
    print *,'nl drtot size',size(nl%drtot)
    print *,'dr size',size(dr)
  endif

  do i = 1, nl%nab
    nl%drtot(i,:) = nl%drtot(i,:) + dr(nl%nlist(i,1),:) + dr(nl%nlist(i,2),:)
  enddo

end subroutine increment_nlist_drtot

subroutine rd_nlist_input(nl,ifile,iunit)
  ! Reads neighbour list variables from input file.
  type (nlist_data_type), intent(out) :: nl
  integer :: iunit
  character(len=32) :: ifile         ! input file name
  open(iunit,file=ifile,status='old',form='formatted', &
     & position='rewind')
  call read_dbl(iunit,'SHELL',nl%shell)
  call read_int(iunit,'MAXPN',nl%mnp)
  close(iunit)

end subroutine rd_nlist_input

function reform_nlist_now(nl)
  ! Returns a logical determining whether to reform the list now or not
  implicit none
  logical :: reform_nlist_now
  type (nlist_data_type), intent(in) :: nl
  double precision,dimension(:),allocatable :: total
  integer i,error
  allocate(total(size(nl%drtot,1)),stat=error)
  if(error .ne. 0) stop 'allocation error in reform_nlist_now'
  if(validate) call validate_nlist(nl)

  do i=1,size(nl%drtot,1)
    total(i) = nl%drtot(i,1)**2 + nl%drtot(i,2)**2
  end do
  
  if(validate) call validate_nlist(nl)
  
  reform_nlist_now = (maxval(abs(total(:) - nl%old_rijsq(:))) .ge. nl%dsqmax)
 
  deallocate(total)

end function reform_nlist_now
  
subroutine validate_nlist(nl)
  !> Checks the entire neighbour list 
  !! checks that all indices are in the required range
  type (nlist_data_type), intent(in) :: nl
  integer :: k

  if(verbose) print *,'validating nlist for',nl%npairs,nl%nab,nl%nip

  do k=1,nl%npairs
    if (nl%iindex(k) > nl%np) stop 'pairs index too large'
    if (nl%jindex(k) > nl%np) stop 'pairs index too large'
    if (nl%iindex(k) .eq. 0 ) stop 'pairs index zero'
    if (nl%jindex(k) .eq. 0) stop 'pairs index zero'
  end do
  
  do k=1,nl%nab
    if (nl%nlist(k,1) > nl%np) stop 'nlist index out of range'
    if (nl%nlist(k,2) > nl%np) stop 'nlist index out of range'
    if (nl%nlist(k,1) .eq. 0) stop 'nlist index out of range'
    if (nl%nlist(k,2) .eq. 0) stop 'nlist index out of range'
  end do
  
  do k=1,nl%nip
    if (nl%ilist(k,1) > nl%np) stop 'ilist index too large'
    if (nl%ilist(k,2) > nl%np) stop 'ilist index too large'
    if (nl%ilist(k,1) .eq. 0) then
      print *,'ilist index zero',nl%ilist(k,1),nl%ilist(k,2),k
      print *,nl%nab
      print *,nl%nip
      stop 'Neighbour list failed validation'
    endif
    if (nl%ilist(k,2) .eq. 0) stop 'ilist index zero'
  end do

  ! validate
  do k=1,nl%nip
    call validate_double(nl%drij(k,1),'rxij____')
    call validate_double(nl%drij(k,2),'ryij____')
    call validate_double(nl%rij(k),   'rij_____')
    call validate_double(nl%rijsq(k), 'rijsq___')
  end do

  ! check allocatables
  if(.not. allocated(nl%drtot)) stop 'drtot not allocated'
  if(.not. allocated(nl%drij)) stop 'drij is not allocated'

end subroutine validate_nlist

end module neighbour_list


