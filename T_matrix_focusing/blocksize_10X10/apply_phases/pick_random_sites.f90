


! this subroutine is based on Fisher–Yates shuffle methods
! the essential part : the do loop is adopted from Rosetta code
! see http://rosettacode.org/wiki/Knuth_shuffle
!pick from an array a of n elements m elements
! the first time when the subroutine is running, it use system time as seed
! the second time it use the new seed = old seed + 999
subroutine Shuffle_and_pick(a, n, b, m)
  integer :: n, m
  integer :: a(n)
  integer :: b(m)
  integer :: i, randpos, temp
  real :: r
  real, external :: ran4
  integer, dimension(8) :: system_time ! use it as a seed
  integer, save :: save_time
  !integer, save :: counter 
  integer :: seed
  character(len=10):: c1, c2
  !data counter /0/
  
  !if (counter==0) then
    !counter = 1
!    call DATE_AND_TIME(values=system_time)
!    save_time = 1000*system_time(7) + system_time(8)
!    seed(1) = save_time
    ! seed is the system time in milliseconds
    !CALL RANDOM_SEED(size = 1)
    open(11,file='fort.13')
    read(11,*) c1, c2, seed
    close(11)
    write(*,*) "seed =", seed
    !close(13)
    !call system("mv fort.10 seed.txt")
  !else
    !seed(1) = save_time + 999
  !endif
  
!  call random_seed(put=seed)
  
  do i = size(a), 2, -1
    !call random_number(r)
    r = ran4(seed)
    randpos = int(r * i) + 1
    temp = a(randpos)
    a(randpos) = a(i)
    a(i) = temp
  end do
  b(1:m) = a(1:m)
end subroutine Shuffle_and_pick



!program Knuth_Shuffle
!  implicit none 
!  integer :: i, n
!  integer, dimension(101) :: a, bins = 0, initial = (/ (n, n=1,101) /) 
!  integer, dimension(10) :: b
  
!    ! the first time
!    a = initial
! 	  call Shuffle_and_pick(a,100,b,10)
!    write(*,"(10(i4))") b
    
!    ! the second time
!    a = initial
! 	  call Shuffle_and_pick(a,100,b,10)
!    write(*,"(10(i4))") b
!end program Knuth_Shuffle
