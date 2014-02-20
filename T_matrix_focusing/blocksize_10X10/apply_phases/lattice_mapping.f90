! map a 1d array to a 2d array
! note the storage of a matrix in fortran is column-wise
! given the index of an element in the 1d array
! if we reshape that 1d_array to a 2d_array
! 2d_array=reshape(1d_array, /(nx, ny)/)
! we want to know the index (xi, yi) of that element in the 2d_array

subroutine map_1d_to_2d(index_in_1d, nx, ny, xi, yi)
  implicit none
  integer :: nx, ny
  integer :: xi, yi, tmp
  integer :: index_in_1d
  
  if (MOD(index_in_1d, ny)==0) then
    xi = index_in_1d/ny
    yi = ny
  else
    xi = index_in_1d/ny + 1
    yi = MOD(index_in_1d, ny)
  endif 
  tmp = xi
  xi = yi
  yi = tmp   
end subroutine map_1d_to_2d
! the result is not what you expect, it first goes through the first row
! so we swap xi, yi in the end

! given the 2d index (xi, yi), output the corresponding 1d index
subroutine map_2d_to_1d(index_in_1d, nx, ny, xi, yi)
  implicit none
  integer :: nx, ny
  integer :: xi, yi, tmp
  integer :: index_in_1d
  
  index_in_1d = xi + (yi-1)*ny
  
end subroutine map_2d_to_1d


! there is a more elegant way to do the mapping
! and the result is what you expect, it first goes through the first column
! but this is slow compared with "map_1d_to_2d"
! so it is not recommended for usage
subroutine Oned_to_2d(index_in_1d, nx, ny, xi, yi)
  implicit none
  integer :: nx, ny
  integer :: xi, yi
  integer :: i
  integer :: index_in_1d
  integer :: loc(2)
  integer :: array_1d(nx*ny)
  integer,allocatable, save :: array_2d(:,:)
  integer, save :: counter
  data counter / 0 /
  
  if (counter==0) then
    counter = 1
    do i=1, nx*ny
      array_1d(i)=i
    end do
    allocate(array_2d(nx,ny))
    array_2d = reshape(array_1d, (/ nx, ny /) )
  endif
  
  loc = minloc(ABS(array_2d - index_in_1d)) 
  xi = loc(1)
  yi = loc(2)   
end subroutine Oned_to_2d


!program test_mapping
!  implicit none
!  integer :: nx, ny, xi, yi
!  integer :: index_in_1d
!  integer :: index_in_1d2
!  integer :: xi2, yi2
!  integer :: i
  
!  nx=101
!  ny=101
!  do i=1, nx*ny
!    index_in_1d = i
!    call map_1d_to_2d(index_in_1d, nx, ny, xi, yi)
!    call Oned_to_2d(index_in_1d, nx, ny, xi2, yi2)
!    write(33,*) i, xi, yi
!    write(44,*) i, xi2, yi2
!    call map_2d_to_1d(index_in_1d2, nx, ny, xi, yi)
!    write(55,*) i, index_in_1d-index_in_1d2
!  end do
  
!end program
