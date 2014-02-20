



!=======================================================================================     
      program wavepacket
        use random_vacancy
        implicit none 
        DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer, parameter :: nx = 101 ! number of molecules in x axis
        integer, parameter :: ny = 101 ! number of molecules in y axis
        integer, parameter :: xmin = -nx/2
        integer, parameter :: xmax =  nx/2     
        integer, parameter :: ymin = -ny/2
        integer, parameter :: ymax =  ny/2   
        integer, parameter :: xfocus = 0
        integer, parameter :: yfocus = 0          
        double precision, parameter :: time_interval = 1.D-3 !2.5D-3/(2*pi)/2 !1.D-5
        double precision, parameter :: total_time = 500*time_interval !0.5D-3 !0.5/(2*pi)/2 !.D-3
        integer :: distance
        integer :: n_iter
        integer :: neq 
        integer :: i,j,k,L
        integer :: xi, yi
        integer :: xloc, yloc
        integer :: loc(2)
        integer :: bx, by ! bx*by represent the size of small blocks in the 2D lattices
        double complex, allocatable :: Y_tmp(:, :) ! Y(-500,500)
        !double complex, allocatable :: Y_tmp2(:,:)
        double complex :: Y(nx*ny) ! a vector of coefficients
        !double complex :: c_mat(xmin:xmax, ymin:ymax)        
        double precision :: p(nx*ny) ! a vector of probabilities
        double precision :: p_mat(xmin:xmax, ymin:ymax)
        double precision :: t, tout

        double precision :: summation
        integer :: flag1, flag2
        character(len=15) :: file_name, file_name2
        character(len=6) :: number_string
        double precision :: Ja, alpha, beta
        double precision :: vacancy_percentage
        integer :: width_x, width_y
        common /flags/ flag1, flag2 !determine how theta changes with time
        common /vacancy_percentage/ vacancy_percentage

        ! specify the size of small blocks
        bx = 10
        by = 10
        
        ! specify the number of percentage here
        vacancy_percentage = 0.0 !0.1 ! 5%
        
        if (ABS(vacancy_percentage) > 1.D-10 ) then
          ! generate the random vacancy sites
          call generate_vacancy(nx, ny, vacancy_percentage)
        endif
        
                
        flag1 = 1  ! theta = 90 degree
        flag2 = 2  ! phi = 0 degree
        neq = nx*ny
           
        call system("rm *.dat *.txt")
        call initialize_molecule_information
        allocate(Y_tmp(xmin:xmax, ymin:ymax))

        ! specify the width of the initial wavepacket
        !width_x = 100
        !width_y = 100
        
        ! the initial coefficient is a Gaussian distribution
        ! C(nx, ny) = EXP( -nx^2/(2*width_x^2) - ny^2/(2*width_y^2) )
        ! it is not normalized
!        do i=xmin, xmax
!          do j=ymin,ymax
!            Y_tmp(i,j) = EXP( -i**2/float(2*width_x**2) - j**2/float(2*width_y**2) )
!            write(88,*) i, j, Y_tmp(i,j)
!          end do
!        end do
        
        !************ STARTING FROM A LOCAL EXCITATION AT THE CENTER ****************
        Y_tmp = (0.D0, 0.D0)
        Y_tmp(0,0) = (1.D0, 0.D0)
        

        
        Y = reshape( Y_tmp, (/neq/) )
        
        
        if (ABS(vacancy_percentage) > 1.D-10 ) then
          ! adding the vacancy sites
          do i=1, n_vacancy
            Y(vacancy_site(i)) = (0.D0, 0.D0) 
          end do 
        endif
                    
        Y_tmp = reshape(Y, (/nx, ny/))
        
        !Y_tmp(xmin:xmax, ymin:ymax) = Y_tmp2(1:nx, 1:ny)
        !deallocate(Y_tmp2)
        
        
        ! normalize the initial wavefunction
        p = ABS(Y)**2
        summation = SUM(p)

!======== record the initial probability distribution===========
        open(unit=10,file="p_0.dat")
        write(10,*) "# time = 0" 
!*******************************************************************************************!
        do xi= xmin, xmax                                                                   !
          do yi = ymin, ymax                                                                !
            write(10,*) xi, yi, abs(Y_tmp(xi,yi))**2/summation ! Y_tmp has already been normalized    !
          end do                                                                            !
        end do                                                                              !
!*******************************************************************************************!        
        close(10)
        
!======== record the initial coefficients =====================
!        open(unit=10,file="c_0.dat")
!        write(10,*) "# time = 0" 
!*******************************************************************************************!
!        do xi= xmin, xmax                                                                   !
!          do yi = ymin, ymax                                                                !
!            write(10,*) xi, yi, Y_tmp(xi,yi) ! Y_tmp has already been normalized            !
!          end do                                                                            !
!        end do                                                                              !
!*******************************************************************************************!        
!        close(10)  
              
! deallocate Y_tmp to save memory
        deallocate(Y_tmp)

        

        T = 0.D0
        TOUT = time_interval
        n_iter = 2*INT(total_time/time_interval)
        
        
        open(29, file="center_locations.dat")

! ************** TIME EVOLUTION CALCULATIONS **************************
        do i=1, n_iter
          !write(*,*) "main 5"
          CALL solve_ode(neq,Y,T,TOUT) 
          !write(*,*) "main 6"
          ! c represents the probabilities at each lattice site
          p = ABS(Y)**2 ! calculate the modulus of a complex number
! ABS can accept an array and return an array, the operation is done element-wise
          summation = SUM(p)
          p = p/summation ! element-wise operation
          

          !c_mat = reshape(Y, (/ nx, ny / )) ! matrix version of Y
          p_mat = reshape(p, (/ nx, ny / )) ! matrix version of c


! find the range of indices of A for the next iteration 
          loc = MAXLOC(p_mat) ! maxloc give the locations with respect to xmin and ymin
          write(*,*) "time = ", t, "s"
          write(*,*) i, MAXVAL(p_mat)
          write(*,*) "===================================="
          
          xloc = loc(1) 
          yloc = loc(2) 

          write(29,*) tout, xloc, yloc


           write(number_string, "(i5)") i 
! ADJUSTL(STRING) will left adjust a string by removing leading spaces. 
! Spaces are inserted at the end of the string as needed. 
           number_string = adjustl(number_string)
           file_name = "p_" // trim(number_string) // ".dat"
           file_name = trim(file_name)
           !file_name2 = "c_" // trim(number_string) // ".dat"
           !file_name2 = trim(file_name2)           
! record probability distribution           
           open(unit=10,file=file_name)
           write(10,*) "# time = ", tout            
!*******************************************************************************************!
           do xi= xmin, xmax                                                                !
             do yi = ymin, ymax                                                             !
               write(10,*) xi, yi, p_mat(xi,yi) ! Y_tmp has already been normalized !
             end do                                                                         !
           end do                                                                           !
!*******************************************************************************************!
           close(10)
           
! record coefficients           
!           open(unit=20,file=file_name2)
!           write(20,*) "# time = ", tout            
!*******************************************************************************************!
!           do xi= xmin, xmax                                                                !
!             do yi = ymin, ymax                                                             !
!               write(20,*) xi, yi, c_mat(xi,yi) ! Y_tmp has already been normalized !
!             end do                                                                         !
!           end do                                                                           !
!*******************************************************************************************!
!           close(20) 
                     
  		   TOUT = TOUT + time_interval

        end do
        
        close(29)
        
        ! record the final wavefunctions
        ! normalize the wavefunction
        Y = Y/sqrt(summation)
        open(15, file='wavefunction_T.bin', form='unformatted')
        write(15) Y
        close(15)
        
      end program










