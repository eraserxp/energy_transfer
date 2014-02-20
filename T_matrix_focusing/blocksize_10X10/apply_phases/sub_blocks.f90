! divide the 2D lattices into sub blocks
! the size of each block is bx*by
! nb_x is the number of blocks along x axis
! nb_y is the number of blocks along y axis
! TNX(Y): total number of lattice sites along x (y) axis 
! nb_x = TNX/bx
! nb_y = TNX/by
subroutine generate_blocks(bx, by, nb_x, nb_y, &
                           TNX, TNY, block_coordinates)
  implicit none
  integer :: bx, by
  integer :: nb_x, nb_y
  integer :: TNX, TNY 
  integer :: block_coordinates(nb_x, nb_y, 4)
  ! ( use nb_x, nb_y ) to label a block
  ! block_coordinates(nb_x, nb_y, 1:2) (x,y) coordinates for
  ! the upper left corner of the block (nb_x, nb_y).
  ! block_coordinates(nb_x, nb_y, 3:4) (x,y) coordinates for
  ! the lower right corner of the block (nb_x, nb_y). 
  integer :: i, j
  integer :: points_x(nb_x) !x coordinates for the left boundaries of all blocks
  integer :: points_y(nb_y) !y coordinates for the top boundaries of all blocks 
  
  do i=1, nb_x
    points_x(i) = 1 + (i-1)*bx
  end do
  
  do i=1, nb_y
    points_y(i) = 1 + (i-1)*by
  end do  
  
  do i=1, nb_x
    do j=1, nb_y
      block_coordinates(i,j,1) = points_x(i) 
      ! x coordinate for the upper left corner of block (i,j)
      block_coordinates(i,j,2) = points_y(j)
      ! y coordinate for the upper left corner of block (i,j)
      block_coordinates(i,j,3) = points_x(i) + bx - 1 
      ! x coordinate for the lower right corner of block (i,j)
      block_coordinates(i,j,4) = points_y(j) + by - 1
      ! y coordinate for the lower right corner of block (i,j)
      if (i==nb_x) then
        block_coordinates(i,j,3) = points_x(i) + bx
      endif
      if (j==nb_y) then
        block_coordinates(i,j,4) = points_y(j) + by
      endif
    end do
  end do
end subroutine generate_blocks



!program test_subroutine
!  implicit none
!  integer, parameter :: TNX=101, TNY=101 
!  integer, parameter :: bx=10, by=10
!  integer, parameter :: nb_x=TNX/bx, nb_y=TNY/by
!  integer :: i, j
!  integer :: block_coordinates(nb_x, nb_y, 4)
  
!  call  generate_blocks(bx, by, nb_x, nb_y, &
!                        TNX, TNY, block_coordinates) 
  
!  do i=1, nb_x
!    do j=1, nb_y
!      write(46,*) i, j, block_coordinates(i,j,1), block_coordinates(i,j,2)
!      write(46,*) i, j, block_coordinates(i,j,3), block_coordinates(i,j,4)      
!    end do
!  end do
!end program test_subroutine
