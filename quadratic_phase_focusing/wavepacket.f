



!=======================================================================================     
      program wavepacket
        implicit none 
        DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer, parameter :: nx = 101 ! number of molecules in x axis
        integer, parameter :: ny = 101 ! number of molecules in y axis
        integer, parameter :: xmin = -nx/2
        integer, parameter :: xmax =  nx/2     
        integer, parameter :: ymin = -ny/2
        integer, parameter :: ymax =  ny/2   
        integer, parameter :: xfocus = -20
        integer, parameter :: yfocus = -20          
        double precision, parameter :: time_interval = 2.5D-3/(2*pi)/2 !1.D-5
        double precision, parameter :: total_time = 0.5/(2*pi)/2 !.D-3
        integer :: distance
        integer :: n_iter
        integer :: neq 
        integer :: i,j,k,L
        integer :: xi, yi
        integer :: xloc, yloc
        integer :: loc(2)
        double complex, allocatable :: Y_tmp(:, :) ! Y(-500,500)
        double complex :: Y(nx*ny) ! a vector of coefficients
        double complex :: c_mat(xmin:xmax, ymin:ymax)        
        double precision :: p(nx*ny) ! a vector of probabilities
        double precision :: p_mat(xmin:xmax, ymin:ymax)
        double precision :: t, tout

        double precision :: summation
        integer :: flag1, flag2
        character(len=10) :: file_name, file_name2
        character(len=5) :: number_string
        double precision :: Ja, alpha, beta

        common /flags/ flag1, flag2 !determine how theta changes with time


        flag1 = 1
        flag2 = 2 
        neq = nx*ny   
        call system("rm *.dat *.txt")
        call initialize_molecule_information
        allocate(Y_tmp(xmin:xmax, ymin:ymax))
        call initial_condition_2d(Y_tmp, nx, ny) ! one wave packet  

!       adding the phases needed for each molecules
        Ja = 22.8257862398621*2*pi !khz
        beta = -25*Ja/12
        alpha = 1/ (4*beta*total_time)
        do i = xmin, xmax
          do j= ymin, ymax
            distance = (i-xfocus)**2 + (j-yfocus)**2
            Y_tmp(i,j) = Y_tmp(i,j)*EXP(-(0.D0, 1.D0)*alpha*distance)
          end do
        end do


            
        Y = reshape( Y_tmp, (/neq/) )
        
!======== record the initial probability distribution===========
        open(unit=10,file="p_0.dat")
        write(10,*) "# time = 0" 
!*******************************************************************************************!
        do xi= xmin, xmax                                                                   !
          do yi = ymin, ymax                                                                !
            write(10,*) xi, yi, abs(Y_tmp(xi,yi))**2 ! Y_tmp has already been normalized    !
          end do                                                                            !
        end do                                                                              !
!*******************************************************************************************!        
        close(10)
        
!======== record the initial coefficients =====================
        open(unit=10,file="c_0.dat")
        write(10,*) "# time = 0" 
!*******************************************************************************************!
        do xi= xmin, xmax                                                                   !
          do yi = ymin, ymax                                                                !
            write(10,*) xi, yi, Y_tmp(xi,yi) ! Y_tmp has already been normalized            !
          end do                                                                            !
        end do                                                                              !
!*******************************************************************************************!        
        close(10)  
              
! deallocate Y_tmp to save memory
        deallocate(Y_tmp)

        

        T = 0.D0
        TOUT = time_interval
        n_iter = 2*INT(total_time/time_interval)
        
        
        open(29, file="center_locations.dat")

        do i=1, n_iter
          CALL solve_ode(neq,Y,T,TOUT) 
          ! c represents the probabilities at each lattice site
          p = ABS(Y)**2 ! calculate the modulus of a complex number
! ABS can accept an array and return an array, the operation is done element-wise
          summation = SUM(p)
          p = p/summation ! element-wise operation
          

          c_mat = reshape(Y, (/ nx, ny / )) ! matrix version of Y
          p_mat = reshape(p, (/ nx, ny / )) ! matrix version of c


! find the range of indices of A for the next iteration 
          loc = MAXLOC(p_mat) ! maxloc give the locations with respect to xmin and ymin
          write(*,*) "time = ", t, "s"
          write(*,*) i, MAXVAL(p_mat)
          write(*,*) " "
          
          xloc = loc(1) 
          yloc = loc(2) 

          write(29,*) tout, xloc, yloc


           write(number_string, "(i5)") i 
! ADJUSTL(STRING) will left adjust a string by removing leading spaces. 
! Spaces are inserted at the end of the string as needed. 
           number_string = adjustl(number_string)
           file_name = "p_" // trim(number_string) // ".dat"
           file_name = trim(file_name)
           file_name2 = "c_" // trim(number_string) // ".dat"
           file_name2 = trim(file_name2)           
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
           open(unit=20,file=file_name2)
           write(20,*) "# time = ", tout            
!*******************************************************************************************!
           do xi= xmin, xmax                                                                !
             do yi = ymin, ymax                                                             !
               write(20,*) xi, yi, c_mat(xi,yi) ! Y_tmp has already been normalized !
             end do                                                                         !
           end do                                                                           !
!*******************************************************************************************!
           close(20) 
                     
  		   TOUT = TOUT + time_interval

        end do

        close(29)
      end program










