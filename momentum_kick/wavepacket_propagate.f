



!=======================================================================================     
      program wavepacket_1d
        implicit none 
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer, parameter :: neq = 1001  !1001  ! neq must be odd
        integer, parameter :: timepoint = 2601
        integer, parameter :: mdim = neq !/2 + 100
        double precision, parameter :: time_interval = 1.D-6 !1.D-5
        double precision, parameter :: time_interval2 = 1.D-6
        double precision :: image_array(timepoint,mdim)
        integer :: i,j,k,L
        double complex :: Y_tmp(-neq/2:neq/2) ! Y(-500,500)
        double complex :: Y(neq)
!        double precision :: initial_phase(-neq/2:neq/2)
!        double precision :: phase_shift(-neq/2:neq/2)
        double precision :: c_k(-neq/2:neq/2)
        double precision :: t, tout
        double precision :: c(neq), summation
		integer :: flag
        character(len=10) :: file_name, file_name2
        character(len=5) :: number_string
        character :: flag_string
        character(len=20) :: fname
        double precision :: ka
        double precision, external :: gaussian_func
        double precision, external :: theta_vs_time
        common /flag/ flag !determine how theta changes with time

        write(*,*) "Specifying the way that field angle changes:"
        write(*,*) "1. remains 90 degree"
        write(*,*) "2. remains 0 degree"
        write(*,*) "3. 0 degree (0-0.3ms); increase to theta_c (0.3-1.4ms); theta_c (>1.4ms)"
        write(*,*) "4. DACOS(SQRT(2.D0/3.D0))"
        write(*,*) "5. 0(0-0.3ms); increase to theta_c(0.3-1.3ms); theta_c(1.3-1.6ms); decrease to 0(1.6-2.6ms)"
        write(*,*) "6. 0(0-0.4ms); increase to theta_c(0.4-1.4ms); theta_c(1.4-2ms); increase to 90(2-3ms)"
        write(*,*) "Input your choice:"
        read(*,*) flag
        !flag = 5        
        call system("rm -rf *.dat *.txt")
        call initialize_molecule_information
        
!        call initial_condition(Y_tmp,neq) ! one wave packet
!        Y(1:neq) = Y_tmp(-neq/2:neq/2)
!        ! adjust the peak position
!        Y = cshift(Y,shift=400)

        open(160,file='Y_final.bin',form='unformatted')
        write(160) Y
        close(160)

        summation = SUM(abs(Y)**2)
        image_array(1,1:mdim) = ABS(Y(1:mdim))**2/summation
        !write(99,*) "# time = 0 s"
        !write(99,*) "#molecule index    population"
        !do i=1, neq
          !write(99,*) i-neq/2-1, ABS(Y(i))**2/summation
        !end do
        !call system("mv fort.99 x_0.dat")



        T = 0.D0
        TOUT = time_interval
        write(flag_string, "(I1)") flag
        do i=1, timepoint-1
          write(*,*) i

          CALL solve_ode(NEQ,Y,T,TOUT) 
        
          c = ABS(Y)**2 ! calculate the modulus of a complex number
! ABS can accept an array and return an array, the operation is done element-wise
          summation = SUM(c)
          c = c/summation ! element-wise operation
          image_array(i+1,1:mdim) = c(1:mdim)
!          call from_x_to_k(Y, neq, c_k)
          
!           write(number_string, "(i5)") i 
!! ADJUSTL(STRING) will left adjust a string by removing leading spaces. 
!! Spaces are inserted at the end of the string as needed. 
!           number_string = adjustl(number_string)
!           file_name = "x_" // trim(number_string) // ".dat"
!           file_name = trim(file_name)
!           open(unit=10,file=file_name)
!           write(10,*) "# time = ", tout 
!           write(10,*) "#molecule index    population"
!           do j=1, neq
!             write(10,*) j-1-neq/2, c(j)
!           end do
!           close(10)
           
!           file_name2 = "k_" // trim(number_string) // ".dat"
!           file_name2 = trim(file_name2)
!           open(unit=19,file=file_name2)
!           write(19,*) "# time = ", tout 
!           do k=-neq/2, neq/2
!             write(19,*) k*pi/(neq/2), C_k(k)
!           end do 
!           close(19)
                     

           open(76,file="flag"//flag_string//"_theta_vs_time.txt", position='append')
           write(76,*) TOUT, theta_vs_time(TOUT, flag)*180/Pi
           close(76)
           TOUT = TOUT + time_interval
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

        fname="flag"//flag_string//".bin"
        !call realimage(timepoint,neq,image_array,fname)
        open(160,file=fname,form='unformatted')
        write(160) image_array
        close(160)
        
      end program










