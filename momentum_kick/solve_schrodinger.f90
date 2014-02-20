! in order to use it, dvode_f90_m.f90 is needed
! blas95 is needed
! solving the Schrodiner-like equation using the real ode solver
! it is not as fast as solve_complex_ode.f, but is more stable

      MODULE provide_hamiltonian
        use blas95
        implicit none
        !double complex, allocatable :: Ham(:,:) ! complex ham is more time-consuming
        double precision, allocatable :: Ham(:,:) !usually hamiltonian is real
        CONTAINS

        subroutine Fex(neq, T, Y_tmp, YDOT)
! you have to provide the Hamiltonian in this subroutine
! if the hamiltonian (or part of it) is independent of time, 
! you may want to calculate it once and save it for later use
          use some_parameters
          implicit none
          integer :: neq ! neq must be even
          integer :: n   ! n is neq/2
          double precision :: T
          DOUBLE precision :: Y_tmp(NEQ), YDOT(NEQ)
! don't define any new variable inside this subroutine, may cause segmentation fault  
          
          external :: wavepacket_ham_AC  !form_ham_matrix ! user provided
        
          n = neq/2 
           
          allocate(Ham(n,n))  
          call wavepacket_ham_AC(n, Ham, T)
                  
!          If (T <= 50*3*10**(-9) *4.1341373374D16 ) then
!            call wavepacket_ham_AC(n, Ham, T) 
!          else
!            call form_matrix_1d(n,Ham,T)
!          endif
          
          ! for complex ham, more time-consuming        
          !YDOT(1:n) = MATMUL(REAL(Ham,kind=8),Y_tmp(n+1:neq)) + MATMUL(AIMAG(Ham),Y_tmp(1:n))
          !YDOT(n+1:neq) = MATMUL(AIMAG(Ham),Y_tmp(n+1:neq)) - MATMUL(REAL(Ham,kind=8),Y_tmp(1:n))
           
          ! for real ham
          call gemv( Ham, Y_tmp(n+1:2*n), YDOT(1:n) )
          call gemv( Ham, -Y_tmp(1:n), YDOT(n+1:2*n) ) 
                       
          deallocate(Ham)
          return 
        END SUBROUTINE Fex

        subroutine JAC(neq,T,Y,ML,MU,PD,NROWPD)
          implicit none
          double precision :: T, Y(neq), PD(NROWPD,NEQ)
          integer :: neq, ML, MU, NROWPD
          integer :: n, i
          double precision :: half_H(neq/2,neq/2)
          external :: wavepacket_ham_AC 
          
          n = neq/2
          PD = 0.D0
          call wavepacket_ham_AC(n, half_H, T) 
          

!          If (T <= 50*3*10**(-9) *4.1341373374D16 ) then
!            call wavepacket_ham_AC(n, Half_H, T) 
!          else
!            call form_matrix_1d(n,Half_H,T)
!          endif
         
          do i = n+1, neq
            PD(1:n,i) = half_H(1:n,i-n)
          end do 
          
          do i = 1, n
            PD(n:neq,i) = -half_H(1:n,i)
          end do 
                
        end subroutine JAC
        
      END MODULE provide_hamiltonian

!******************************************************************



! solving the schrodinger equation using dvode
!=======================================================================================     
      subroutine solve_schrodinger(n,Y,t,tout)
        USE DVODE_F90_M
        USE provide_hamiltonian
        implicit none 
        integer :: n
        integer :: neq 
        double precision :: t, tout
        double complex :: Y(n)
        double precision :: Y_tmp(2*n)
        DOUBLE PRECISION :: ATOL, RTOL
        INTEGER :: ITASK, ISTATE, ISTATS, IOUT, IERROR

        TYPE (VODE_OPTS) :: OPTIONS  

        RTOL = 1.D-6
        ATOL = 1.D-6
        ITASK = 1
        ISTATE = 1
        NEQ = 2*n

        Y_tmp(1:n) = REAL(Y)
        Y_tmp(n+1:2*n) = AIMAG(Y)
        
        OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE.,&
                   ABSERR=ATOL,RELERR=RTOL, &
                   USER_SUPPLIED_JACOBIAN=.TRUE.)

        CALL DVODE_F90(FEX,NEQ,Y_tmp,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JAC)
        
        Y(1:n) = Y_tmp(1:n) + DCMPLX(0.D0,1.D0)*Y_tmp(n+1:2*n)

      end subroutine solve_schrodinger

    
!      program test
!        implicit none
!        integer, parameter :: n = 2
!        double precision, parameter :: dt = 1.D-2
!        double precision :: t, tout
!        double complex :: Y(n)
!        integer :: i
        
!        t = 0
!        tout = t + dt
!        Y(1) = (1.D0,0.D0)
!        Y(2) = (0.D0,0.D0)
!        write(10,*) tout, ABS( Y(1) )**2, ABS( Y(2) )**2  
              
!        do i = 1, 1000
!          call solve_schrodinger(n,Y,t,tout)
!          write(10,*) tout, ABS( Y(1) )**2, ABS( Y(2) )**2
!          tout = tout + dt
!        end do
!        call system('xmgrace -nxy fort.10')
!      end program test


!      subroutine form_ham_matrix(n, Ham,T)
!        implicit none
!        integer :: n
!        double precision :: T
!        double complex :: Ham(n,n)
      
!        Ham(1,1) = DCMPLX(4.D0,0.D0)
!        Ham(1,2) = DCMPLX(1.D0,0.D0)
!        Ham(2,1) = DCMPLX(1.D0,0.D0)
!        Ham(2,2) = DCMPLX(4.D0,0.D0)       
!      end subroutine form_ham_matrix




























