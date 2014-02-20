



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
        integer, parameter :: bx=10, by=10 
        ! bx*by represent the size of small blocks in the 2D lattices 
        integer, parameter :: nb_x=nx/bx, nb_y=ny/by
        integer, parameter :: neq = nx*ny       
        double precision, parameter :: time_interval = 1.D-3 !2.5D-3/(2*pi)/2 !1.D-5
        double precision, parameter :: total_time = 500*time_interval !0.5/(2*pi)/2 !.D-3
        
        !integer :: distance
        integer :: n_iter

        integer :: i,j,k,L
        integer :: xi, yi
        integer :: xloc, yloc
        integer :: loc(2)


        double complex, allocatable :: Y_tmp(:, :) ! Y(-500,500)
        !double complex, allocatable :: Y_tmp2(:,:)
        double complex :: Y(nx*ny) ! a vector of coefficients
        double complex :: wavefunction_T(nx*ny) 
        !the wavefunction calculated from evolution starting from a local excitation
        !double complex :: c_mat(xmin:xmax, ymin:ymax)   
        double complex :: phases(nb_x,nb_y)  !the phases added to different blocks   
        double precision :: p(nx*ny) ! a vector of probabilities
        double precision :: p_mat(xmin:xmax, ymin:ymax)
        double precision :: t, tout

        double precision :: summation
        integer :: flag1, flag2
        integer :: block_coordinates(nb_x, nb_y, 4)
        integer :: index_in_1d
        integer :: x_index_min, x_index_max, y_index_min, y_index_max
        character(len=13) :: file_name, file_name2
        character(len=6) :: number_string
        double precision :: Ja, alpha, beta
        double precision :: vacancy_percentage
        integer :: width_x, width_y
        common /flags/ flag1, flag2 !determine how theta changes with time
        common /vacancy_percentage/ vacancy_percentage

 
        
        ! specify the number of percentage here
        vacancy_percentage = 0.0 !0.1 ! 5%
        
        if (ABS(vacancy_percentage) > 1.D-10 ) then
          ! generate the random vacancy sites
          call generate_vacancy(nx, ny, vacancy_percentage)
        endif
        
                
        flag1 = 1  ! theta = 90 degree
        flag2 = 2  ! phi = 0 degree
           
        call system("rm *.dat *.txt")
        call initialize_molecule_information
        allocate(Y_tmp(xmin:xmax, ymin:ymax))
        !write(*,*) "1"

        !************* START WITH A WAVEPACKET ************************************
!        ! specify the width of the initial wavepacket
!        width_x = 100
!        width_y = 100
        
!        ! the initial coefficient is a Gaussian distribution
!        ! C(nx, ny) = EXP( -nx^2/(2*width_x^2) - ny^2/(2*width_y^2) )
!        ! it is not normalized
!        do i=xmin, xmax
!          do j=ymin,ymax
!            Y_tmp(i,j) = EXP( -i**2/float(2*width_x**2) - j**2/float(2*width_y**2) )
!            write(88,*) i, j, Y_tmp(i,j)
!          end do
!        end do
        !****************************************************************************



        !************* START WITH A PLANEWAVE ************************************

        ! the initial coefficient is a plane wave
        do i= 1, nx*ny
            Y(i) = 1.D0/sqrt(dfloat(nx*ny))
        end do
        !****************************************************************************
        !write(*,*) "2"
        
        ! obtain the phases from "wavefunction_T.bin"
        open(160,file='wavefunction_T.bin',form='unformatted')
        read(160) wavefunction_T  ! or read(160) potentialC
        close(160)
        !write(*,*) "3"
        
        call  generate_blocks(bx, by, nb_x, nb_y, NX, NY, block_coordinates)
        !write(*,*) "4"
        
        do i=1, nb_x
          do j=1, nb_y
            
            x_index_min = block_coordinates(i,j,1)
            y_index_min = block_coordinates(i,j,2)
            x_index_max = block_coordinates(i,j,3)
            y_index_max = block_coordinates(i,j,4)
            !write(*,*) "5"
            
            !*************obtaining the phase needed for block (i,j)**********
            phases(i,j) = (0.D0, 0.D0)
            do xi = x_index_min, x_index_max
              do yi = y_index_min, y_index_max
                !write(*,*) "6"
                call map_2d_to_1d(index_in_1d, nx, ny, xi, yi)
                phases(i,j) = phases(i,j) + wavefunction_T(index_in_1d)
              end do
            end do
            
            !normalize the phases
            phases(i,j)=phases(i,j)/ABS(phases(i,j))
            ! complex conjugate
            phases(i,j)=CONJG(phases(i,j))
            !**************************************************************
            
            !*************adding the phase needed for block (i,j)**********
            do xi = x_index_min, x_index_max
              do yi = y_index_min, y_index_max
                call map_2d_to_1d(index_in_1d, nx, ny, xi, yi)
                ! don't apply phase to the center of the lattices
                if ( (xi==nx/2).and.(yi==ny/2) ) then
                  Y(index_in_1d) = Y(index_in_1d)
                else
                  Y(index_in_1d) = phases(i,j)*Y(index_in_1d)  
                endif
              end do
            end do 
            !**************************************************************           
          end do
        end do
        
                
        
        
        if (ABS(vacancy_percentage) > 1.D-10 ) then
          ! adding the vacancy sites
          do i=1, n_vacancy
            Y(vacancy_site(i)) = (0.D0, 0.D0) 
          end do 
        endif
                    
        Y_tmp = reshape(Y, (/nx, ny/))
        
        
        
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
      end program










