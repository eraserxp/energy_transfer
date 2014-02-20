! the module generate random vacancy sites
module random_vacancy
  implicit none
  integer, save :: n_vacancy ! the number of vacancy sites 
  integer, allocatable, save :: vacancy_site(:) ! store the vacancy sites 
  
  contains 
  
  subroutine generate_vacancy(nx, ny, vacancy_percentage)
    implicit none
    integer :: nx, ny, neq
    integer :: i
    integer :: xi, yi
    real*8 :: vacancy_percentage
    integer, allocatable :: all_sites(:)
    integer :: index_in_1d 
    
    
    !write(*,*) "random vacancy 1"
    n_vacancy = INT(vacancy_percentage*nx*ny) + 1
    vacancy_percentage = 100*float(n_vacancy)/(nx*ny)
    !open(23, file="vacancy_sites.txt")
    write(23,*) "The percentage of vacancy is ", vacancy_percentage, "%"
    write(23,*) "There is ", n_vacancy, "vacancy sites in the crystal"
    
    neq = nx*ny
    allocate(all_sites(neq))    
    ! generate the random vacancy sites
    do i=1, neq
      all_sites(i) = i
    end do
    allocate(vacancy_site(n_vacancy))
    call Shuffle_and_pick(all_sites, neq, vacancy_site, n_vacancy)
    
    deallocate(all_sites)
    
    write(23,*) "These vacancy sites are:"
    do i=1, n_vacancy
      index_in_1d = vacancy_site(i)
      call map_1d_to_2d(index_in_1d, nx, ny, xi, yi)
      !call Oned_to_2d(index_in_1d, nx, ny, xi, yi)
      write(23,*) i, "  (", xi, yi, ")"       
    end do
    close(23) 
  end subroutine generate_vacancy
end module random_vacancy
