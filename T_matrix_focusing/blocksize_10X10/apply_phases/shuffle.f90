program Knuth_Shuffle
  implicit none 
  integer, parameter :: reps = 1000000
  integer :: i, n
  integer, dimension(100) :: a, bins = 0, initial = (/ (n, n=1,100) /) 
  integer, dimension(10) :: b
  
!  do i = 1, 100 !reps
    a = initial
 	  call Shuffle_and_pick(a,100,b,10)
    write(*,"(10(i4))") b
    !where (a == initial) bins = bins + 1  ! skew tester
!  end do
!  write(*, "(10(i8))") bins
! prints  100382  100007   99783  100231  100507   99921   99941  100270  100290  100442

end program Knuth_Shuffle


! this subroutine is based on Fisher–Yates shuffle methods
! the essential part : the do loop is adopted from Rosetta code
! see http://rosettacode.org/wiki/Knuth_shuffle
!pick from an array a of n elements m elements
subroutine Shuffle_and_pick(a, n, b, m)
  integer :: n, m
  integer :: a(n)
  integer :: b(m)
  integer :: i, randpos, temp
  real :: r
  integer, dimension(8) :: system_time ! use it as a seed
  integer :: seed(1)
  
  call DATE_AND_TIME(values=system_time)
  seed(1) = 1000*system_time(7) + system_time(8)
  ! seed is the system time in milliseconds
  !CALL RANDOM_SEED(size = 1)
  write(*,*) seed
  !call random_seed()
  call random_seed(put=seed)
  
  do i = size(a), 2, -1
    call random_number(r)
    randpos = int(r * i) + 1
    temp = a(randpos)
    a(randpos) = a(i)
    a(i) = temp
  end do
  b(1:m) = a(1:m)
end subroutine Shuffle_and_pick
