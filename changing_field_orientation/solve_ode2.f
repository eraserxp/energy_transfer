

      module provide_hamiltonian
        implicit none
        double precision,allocatable :: ham(:,:)
        contains
!==========================================================================================
        subroutine fex(neq, T, Y, YDOT, RPAR, IPAR)
 
          implicit none
          integer :: neq
          double precision :: T
          DOUBLE complex :: Y(NEQ), YDOT(NEQ), RPAR !real parameters
          integer :: IPAR !integer parameters
          external :: form_hamiltonian

          allocate(ham(neq,neq))
          call form_matrix_2d(neq, Ham, T)
!          call form_matrix_1d(neq,Ham,T)
          YDOT = RPAR*MATMUL(Ham,Y)
          deallocate(ham)
        END SUBROUTINE

!===========================================================================================
        SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
          implicit none
          integer :: NEQ
          integer :: NRPD
          integer :: IPAR
          integer :: ML, MU
          DOUBLE COMPLEX Y(NEQ), PD(NRPD,NEQ), RPAR
          DOUBLE PRECISION T
          RETURN
        END SUBROUTINE

      end module provide_hamiltonian





      subroutine solve_ode(NEQ,Y,T,TOUT)
        use provide_hamiltonian
        implicit none
        integer :: neq
        integer :: IPAR
        integer :: ITOL,ITASK,ISTATE,IOPT,LZW,LRW,LIW,MF
        DOUBLE COMPLEX :: Y(NEQ), ZWORK(15*NEQ + 10), RPAR
        DOUBLE PRECISION ::  ATOL, RTOL, RWORK(20 + NEQ + 10), T, TOUT
        integer :: IWORK(30 + 10)

        
        ITOL = 1
        RTOL = 1.D-10
        ATOL = 1.D-10
        ITASK = 1
        ISTATE = 1
        IOPT = 0
        LZW = 15*NEQ + 10 
        LRW = 20 + NEQ + 10 
        LIW = 30 + 10
        MF = 10
        RPAR = DCMPLX(0.0D0,1.0D0)

        CALL ZVODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,
     &             ZWORK,LZW,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR) ! segementation fault occurs here
!        write(*,*) "============================="
      end subroutine











     















!=======================================================================================     
!      program test_solve_ode
!        implicit none 
!        integer,parameter :: neq = 7
!        integer :: i,j
!        double complex :: Y(neq)
!        double precision :: t, tout
!        double precision :: c1, c2, c3, c4, c5, c6, c7, summation
        
!        Y(1) = (1.D0,0.D0)
!        Y(2:neq) = (0.D0, 0.D0)
        

!        do i=1, 100001
!          T = (i-1)*1.D-4
!          TOUT = i*1.D-4
!          CALL solve_ode(NEQ,Y,T,TOUT)
!          summation = 0.D0
!          c1 = real(y(1))**2 + aimag(y(1))**2
!          c2 = real(y(2))**2 + aimag(y(2))**2
!          c3 = real(y(3))**2 + aimag(y(3))**2
!          c4 = real(y(4))**2 + aimag(y(4))**2
!          c5 = real(y(5))**2 + aimag(y(5))**2
!          c6 = real(y(6))**2 + aimag(y(6))**2
!          c7 = real(y(7))**2 + aimag(y(7))**2
!          do j=1, neq
!             summation = summation + real(y(j))**2 + aimag(y(j))**2
!          end do
!          summation = c1 + c2 + c3 + c4 + c5 + c6 + c7
!          c1 = c1/summation
!          c2 = c2/summation
!          c3 = c3/summation
!          c4 = c4/summation
!          c5 = c5/summation
!          c6 = c6/summation
!          c7 = c7/summation
          
!          write(19,'(8(2XE18.10))') TOUT, c1, c2, c3, c4, c5, c6, c7    		
!        end do

!        call system("xmgrace -nxy fort.19")
!      end program






