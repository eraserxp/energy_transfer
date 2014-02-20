



!=======================================================================================     
      program momentum_kick
        use some_parameters
        implicit none 
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer, parameter :: neq = 201  !1001  ! neq must be odd
        integer, parameter :: ntime = 101
        double precision :: time_interval2
        integer :: i,j,k,L
        double complex :: Y_tmp(-neq/2:neq/2) ! Y(-500,500)
        double complex ::    Y(neq)
        double complex ::    Y_k(neq)
        double complex ::    x_complex_image(ntime,neq)
        double precision ::  x_real_image(ntime,neq)
        double precision ::  c_k(-neq/2:neq/2)
        double complex ::    k_complex_image(ntime,neq)
        double precision ::  k_real_image(ntime,neq)
        double precision :: t, tout, t_tmp
        double precision :: c(neq), summation
		    integer :: flag
        character(len=10) :: file_name, file_name2
        character(len=5) :: number_string
        double precision :: ka
        double precision, external :: gaussian_func
        common /flag/ flag !determine how theta changes with time

        flag = 1 
        !displacement =  460      
        call system("rm -rf *.dat")
        write(*,*) "Obtaining parameters from files"
        call initialize_molecule_information
        
        call initial_condition(Y_tmp,neq) ! one wave packet
        Y(1:neq) = Y_tmp(-neq/2:neq/2) ! Y_tmp is already normalized
        x_complex_image(1,1:neq) = Y(1:neq)
        
        !Y = cshift(Y,displacement)
        c = ABS(Y)**2
        summation = SUM(c)
        c = c/summation
        x_real_image(1,1:neq) = c(1:neq)
        
        open(11,file='x_0.dat')
        do j=1, neq
          write(11,*) j, c(j)
        end do
        close(11)

        call from_x_to_k(Y, neq, Y_k, c_k)
        k_complex_image(1,1:neq) = Y_k(-neq/2:neq/2)
        k_real_image(1,1:neq) = c_k(-neq/2:neq/2) 
               
        open(12,file='k_0.dat')
        do k=-neq/2, neq/2
          write(12,*) k*pi/(neq/2), C_k(k)
        end do 
        close(12)

        

        T = 0.D0
        time_interval2 = time_interval
        write(*,*) t_var, time_interval
        call time_in_atomic_unit(time_interval2)        
        TOUT = time_interval2
        
        do i=1,  n_iterate 
          write(*,*) i
          !CALL solve_ode(NEQ,Y,T,TOUT) 
          CALL solve_schrodinger(NEQ,Y,T,TOUT)
          
            
          if ( MOD(i,10)==0 ) then
            c = ABS(Y)**2 ! calculate the modulus of a complex number
! ABS can accept an array and return an array, the operation is done element-wise
            summation = SUM(c)
            c = c/summation ! element-wise operation
            
            call from_x_to_k(Y, neq, Y_k, c_k) ! after this call, Y_k,c_k is already normalized           
            !----------- produce complex image array -------------------
              x_complex_image(i/10+1,1:neq) = Y(1:neq)/sqrt(summation)
              k_complex_image(i/10+1,1:neq) = Y_k(1:neq) 
            !-----------------------------------------------------------

            !----------- produce real image array -------------------
              x_real_image(i/10+1,1:neq) = c(1:neq)
              k_real_image(i/10+1,1:neq) = c_k(1:neq)
            !-----------------------------------------------------------            

          
            write(number_string, "(i5)") i/10
! ADJUSTL(STRING) will left adjust a string by removing leading spaces. 
! Spaces are inserted at the end of the string as needed. 
            number_string = adjustl(number_string)
            file_name = "x_" // trim(number_string) // ".dat"
            file_name = trim(file_name)
            open(unit=10,file=file_name)
            write(10,*) "# time = ", tout/4.134137D16 
            do j=1, neq
              write(10,*) j, c(j)
            end do
            close(10)
           
            file_name2 = "k_" // trim(number_string) // ".dat"
            file_name2 = trim(file_name2)
            open(unit=19,file=file_name2)
            write(19,*) "# time = ", tout/4.134137D16 
            do k=-neq/2, neq/2
              write(19,*) k*pi/(neq/2), C_k(k)
            end do 
            close(19)
          endif 
                  
  		  TOUT = TOUT + time_interval2
        end do

!*********************************************************************************
!*********************************************************************************
!        do i=101, 1000
!          CALL solve_ode(NEQ,Y,T,TOUT) 
!          c = ABS(Y)**2 ! calculate the modulus of a complex number
!! ABS can accept an array and return an array, the operation is done element-wise
!          summation = SUM(c)
!          c = c/summation ! element-wise operation

!           write(number_string, "(i5)") i 
!! ADJUSTL(STRING) will left adjust a string by removing leading spaces. 
!! Spaces are inserted at the end of the string as needed. 
!           number_string = adjustl(number_string)
!           file_name = "x_" // trim(number_string) // ".dat"
!           file_name = trim(file_name)
!           open(unit=10,file=file_name)
!           write(10,*) "# time = ", tout 
!           do j=1, neq
!             write(10,*) j-1-neq/2, c(j)
!           end do
!           close(10)
!           TOUT = TOUT + time_interval2
!        end do
!*********************************************************************************
!*********************************************************************************

        
        ! the shape of the wavepacket in k space
!        do L=1,500
!          ka=L*3.14/500
!          write(20,*) ka, gaussian_func(ka)
!        end do


        open(160,file='x_complex_image.bin',form='unformatted')
        write(160) x_complex_image  ! or read(160) potentialC
        close(160)

        open(160,file='k_complex_image.bin',form='unformatted')
        write(160) k_complex_image  ! or read(160) potentialC
        close(160)  
        
        open(160,file='x_real_image.bin',form='unformatted')
        write(160) x_real_image  ! or read(160) potentialC
        close(160)

        open(160,file='k_real_image.bin',form='unformatted')
        write(160) k_real_image  ! or read(160) potentialC
        close(160)
                              
      end program










