



!=======================================================================================     
      program wavepacket
        implicit none 

        integer, parameter :: nx = 101 ! number of molecules in x axis
        integer, parameter :: ny = 101 ! number of molecules in y axis
        double precision, parameter :: time_interval = 1.D-5
        integer :: neq 
        integer :: i,j,k,L
        integer :: xi, yi
        integer :: xindex, yindex
        integer :: xmin, ymin, xmax, ymax
        integer :: xmin_new, ymin_new, xmax_new, ymax_new
        integer :: xloc, yloc
        integer :: loc(2)
        double complex, allocatable :: Y_tmp(:, :) ! Y(-500,500)
        double precision, allocatable :: A(:,:)
        double complex :: Y(nx*ny)
        double complex, allocatable :: Y_new(:,:)
        double complex, allocatable :: Y_old(:,:)
        double precision :: t, tout
        double precision, allocatable :: c(:) 
        double precision :: summation
		integer :: flag1, flag2
        character(len=10) :: file_name
        character(len=5) :: number_string
        double precision :: ka

        common /flags/ flag1, flag2 !determine how theta changes with time


        flag1 = 7
        flag2 = 2 
        neq = nx*ny       

        call initialize_molecule_information
        allocate(Y_tmp(-nx/2:ny/2, -ny/2:ny/2))
        call initial_condition_2d(Y_tmp, nx, ny) ! one wave packet      
        Y = reshape( Y_tmp, (/neq/) )
! record the initial probability distribution
        open(unit=10,file="t_0.dat")
        write(10,*) "# time = 0" 
        do xi = -500, 500
          do yi = -500, 500
            if ( ((xi >= -50).and.(xi <=50)).and.((yi >= -50).and.(yi<=50)) ) then
              write(10,*) xi, yi, abs(Y_tmp(xi,yi))**2 ! Y_tmp has already been normalized
            else
              write(10,*) xi, yi, "0"
            endif
          end do
          write(10,*) " "
        end do
        close(10)
! deallocate Y_tmp to save memory
        deallocate(Y_tmp)

        

        T = 0.D0
        TOUT = time_interval


        xmin = -50 
        xmax = 50
        ymin = -50
        ymax = 50
        allocate(A(xmin:xmax, ymin:ymax))
        allocate(Y_old(xmin:xmax, ymin:ymax))

        open(29, file="center_locations.dat")

        do i=1, 300
          CALL solve_ode(neq,Y,T,TOUT) 
          allocate(c(nx*ny))
          c = ABS(Y)**2 ! calculate the modulus of a complex number
! ABS can accept an array and return an array, the operation is done element-wise
          summation = SUM(c)
          c = c/summation ! element-wise operation
          
          if (i/=1) then
! xmin .. acquire their values from last iteraction
            allocate(A(xmin:xmax, ymin:ymax))
            allocate(Y_old(xmin:xmax, ymin:ymax))
          endif

          Y_old = reshape(Y, (/ 101, 101/ ))
          A = reshape(c, (/ 101, 101/ )) 
          deallocate(c) ! to save memory

! find the range of indices of A for the next iteration 
          loc = MAXLOC(A) ! maxloc give the locations with respect to xmin and ymin
          write(*,*) MAXVAL(A)
          xloc = loc(1) + xmin -1
          yloc = loc(2) + ymin -1

          write(29,*) tout, xloc, yloc
     
          xmin_new = xloc - 50
          xmax_new = xloc + 50
          ymin_new = yloc - 50
          ymax_new = yloc + 50
          

           write(number_string, "(i5)") i 
! ADJUSTL(STRING) will left adjust a string by removing leading spaces. 
! Spaces are inserted at the end of the string as needed. 
           number_string = adjustl(number_string)
           file_name = "t_" // trim(number_string) // ".dat"
           file_name = trim(file_name)
           open(unit=10,file=file_name)
           write(10,*) "# time = ", tout 
           do j=-500, 500
              do k=-500, 500
                if ( ((j >= xmin).and.(j<=xmax)).and.((k >= ymin).and.(k<=ymax)) ) then
                  write(10,*) j, k, A(j, k) 
                else
                  write(10,*) j, k, "0"
                endif
              end do
              write(10,*) " "
           end do
           close(10)
  		   TOUT = TOUT + time_interval

          deallocate(A)

!-------------------- shift the calculation zone---------------------------------------
          allocate( Y_new(xmin_new:xmax_new, ymin_new:ymax_new) )
! set the value of Y to be (0,0)
          Y_new = (0.D0, 0.D0)
          do xindex = xmin_new, xmax_new
            do yindex = ymin_new, ymax_new
              if ( ((xindex >= xmin).and.(xindex<=xmax)).and.((yindex >= ymin).and.(yindex<=ymax)) ) then
                Y_new(xindex,yindex) = Y_old(xindex,yindex)           
              end if
            end do
          end do
          Y = reshape(Y_new, (/ nx*ny /))
          deallocate(Y_new)
!--------------------------------------------------------------------------------------
          deallocate(Y_old)

! assign the values for the next iteration
           xmin = xmin_new
           xmax = xmax_new
           ymin = ymin_new
           ymax = ymax_new


        end do

        close(29)
      end program










