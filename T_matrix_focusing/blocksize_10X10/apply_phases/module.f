



!==========================================================================================
      module molecule_information
! to save the values for dipole, electric_field, rot_const, lattice_constant
        implicit none
        double precision, save :: dipole
        double precision, save :: electric_field
        double precision, save :: rot_const
        double precision, save :: lattice_constant
      end module molecule_information

!=======================================================================================
      subroutine initialize_molecule_information
        use molecule_information

      	rot_const=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	Dipole=5.529D0 
        Lattice_Constant=4.D-7         
        electric_field = 1.D5 ! V/m

        CALL Dipole_in_atomic_unit(Dipole)
        CALL Length_in_Bohr(Lattice_Constant)
        CALL Hz_to_atomic_unit(rot_const)
        CALL Field_in_atomic_unit(electric_field)
      end subroutine initialize_molecule_information 














































C**********************************************************************
C
C
C
      FUNCTION THRJ(F1,F2,F3,G1,G2,G3)
C  the input variables should be of double precision types ----added by Ping
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SMALL CHANGES 31 JUL 95 (SG)
      SAVE MUNG,X,Y
      PARAMETER (MXIX=302)
      DIMENSION X(MXIX),Y(MXIX)
      DATA MUNG/0/
      IF (MUNG.EQ.21) GO TO 69
      MUNG = 21
      X(1) = 0.D0
      DO 100 I = 1, MXIX-1
      A = I
      X(I+1) = LOG(A) +X(I)
      Y(I+1) = LOG(A)
  100 CONTINUE
   69 IF(F1-ABS(G1)) 1,13,13
   13 IF(F2-ABS(G2))1,14,14
   14 IF(F3-ABS(G3))1,15,15
   15 SUM=F1+F2+F3
      NSUM=SUM+.001D0
      IF(SUM-NSUM)2,2,1
    1 THRJ=0.D0
      RETURN
    2 IF(ABS(G1+G2+G3)-1.D-08)3,3,1
    3 IF(F1+F2-F3)1,4,4
    4 IF(F1+F3-F2)1,5,5
    5 IF(F2+F3-F1)1,6,6
    6 J1=2.D0*F3+2.001D0
      J2=F1+F2-F3+1.001D0
      J3=F1-F2+F3+1.001D0
      J4=-F1+F2+F3+1.001D0
      J5=F1+F2+F3+2.001D0
      J6=F1+G1+1.001D0
      J7=F1-G1+1.001D0
      J8=F2+G2+1.001D0
      J9=F2-G2+1.001D0
      J10=F3+G3+1.001D0
      J11=F3-G3+1.001D0
      IF(J5.GT.MXIX) THEN
        WRITE(6,601) J5,MXIX
  601   FORMAT(' *** DIMENSION ERROR IN THRJ - INDEX.GT.MXIX',2I5)
        STOP
      ENDIF
      R=0.5D0*(Y(J1)+X(J2)+X(J3)+X(J4)-X(J5)
     1+X(J6)+X(J7)+X(J8)+X(J9)+X(J10)+X(J11))
      SUM=0.D0
      F=-1
      KZ=-1
    7 KZ=KZ+1
      F=-F
      J1=KZ+1
      J2=F1+F2-F3-KZ+1.001D0
      IF(J2)20,20,8
    8 J3=F1-G1-KZ+1.001D0
      IF(J3)20,20,9
    9 J4=F2+G2-KZ+1.001D0
      IF(J4)20,20,10
   10 J5=F3-F2+G1+KZ+1.001D0
      IF(J5)7,7,11
   11 J6=F3-F1-G2+KZ+1.001D0
      IF(J6)7,7,12
   12 JMAX=MAX(J1,J2,J3,J4,J5,J6)
      IF(JMAX.GT.MXIX) THEN
        WRITE(6,601) JMAX,MXIX
        STOP
      ENDIF
      S=-(X(J1)+X(J2)+X(J3)+X(J4)+X(J5)+X(J6))
      SUM=SUM+F*EXP(R+S)
      GO TO 7
   20 INT=ABS(F1-F2-G3)+0.0001D0
      VAL=((-1.D0)**INT)*SUM/SQRT(2.D0*F3+1.D0)
      IF(ABS(VAL).LE.1.D-6) VAL=0.D0
      THRJ=VAL
      RETURN
      END
C

C
C***********************************************************************



C******* DIAGONALIZE THE REAL SYMMETRIC MATRIX *******************
C
C            ..... good for small matrices .....
C
      SUBROUTINE JACOBI(A,N,NP,D,V,NROT)
C  elements of A above the diagonal are destroyed

Cf2py intent(in) N
Cf2py intent(in) NP
Cf2py intent(in,out) A 
Cf2py intent(out) D
Cf2py depend(NP) D
Cf2py intent(out) V
Cf2py depend(NP) V
Cf2py intent(out) NROT

      IMPLICIT REAL*8 (A-H,O-Z) 
      PARAMETER (NMAX=250)
      DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)
      DO 12 IP=1,N
        DO 11 IQ=1,N
          V(IP,IQ)=0.
11      CONTINUE
        V(IP,IP)=1.
12    CONTINUE
      DO 13 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.
13    CONTINUE
      NROT=0
      DO 24 I=1,50
        SM=0.
        DO 15 IP=1,N-1
          DO 14 IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
14        CONTINUE
15      CONTINUE
        IF(SM.EQ.0.)RETURN
        IF(I.LT.4)THEN
          TRESH=0.2*SM/N**2
        ELSE
          TRESH=0.
        ENDIF
        DO 22 IP=1,N-1
          DO 21 IQ=IP+1,N
            G=100.*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP)))
     *         .AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ)=0.
            ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(ABS(H)+G.EQ.ABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5*H/A(IP,IQ)
                T=1./(ABS(THETA)+SQRT(1.+THETA**2))
                IF(THETA.LT.0.)T=-T
              ENDIF
              C=1./SQRT(1+T**2)
              S=T*C
              TAU=S/(1.+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=0.
              DO 16 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
16            CONTINUE
              DO 17 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
17            CONTINUE
              DO 18 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
18            CONTINUE
              DO 19 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
19            CONTINUE
              NROT=NROT+1
            ENDIF
21        CONTINUE
22      CONTINUE
        DO 23 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.
23      CONTINUE
24    CONTINUE
      PAUSE '50 iterations should never happen'
      RETURN
      END


C**** SORT THE EIGENVALUES AND VECTORS OF THE DIAGONALIZED MATRIX *********
C
      SUBROUTINE EIGSRT(D,V,N,NP)
C
C  Modified by R. Krems to sort the eigenvalues and 
C  eigenvectors in the ascending rather then descending 
C  order
C               April 20, 2000, Goteborg
C
C

Cf2py intent(in,out) D
Cf2py intent(in,out) V
Cf2py intent(in) N
Cf2py intent(in) NP

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION D(NP),V(NP,NP)
      DO 13 I=1,N-1
        K=I
        P=D(I)
        DO 11 J=I+1,N
          IF(D(J).LE.P)THEN
            K=J
            P=D(J)
          ENDIF
11      CONTINUE
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO 12 J=1,N
            P=V(J,I)
            V(J,I)=V(J,K)
            V(J,K)=P
12        CONTINUE
        ENDIF
13    CONTINUE
      RETURN
      END
C********************************************************************      



!***********************************************************************
! input A = matrix to be diagonalized
! ouput A = orthonormal eigenvector matrix
! N ----- dimension of A
! W ----- eigenvalue array
! if the subroutine fails, it will generate a error message: lapack_eig.err
      SUBROUTINE lapack_eig(A, N, W)

Cf2py intent(in,out) A
Cf2py intent(in) N
Cf2py intent(out) W

      INTEGER          N
      INTEGER          LDA
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 100000000  ) 
! LWMAX >= 1 + 6*N + 2*N**2 
! don't delete () in the above line
!
!     .. Local Scalars ..
      INTEGER          INFO, LWORK, LIWORK
!
!     .. Local Arrays ..
      INTEGER          IWORK( LWMAX )
      DOUBLE PRECISION :: A( N, N )  
      double precision :: W( N ), WORK( LWMAX ) 

!
!     .. External Subroutines ..
      EXTERNAL         DSYEVD

!
!     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
!
      LDA=N
!
!     Query the optimal workspace.
!
      LWORK = -1
      LIWORK = -1
      CALL DSYEVD( 'Vectors', 'Upper', N, A, LDA, W, WORK, LWORK,
     $             IWORK, LIWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      LIWORK = MIN( LWMAX, IWORK( 1 ) )
!
!     Solve eigenproblem.
!
      CALL DSYEVD( 'Vectors', 'Upper', N, A, LDA, W, WORK, LWORK,
     $             IWORK, LIWORK, INFO )
!
!     Check for convergence.
!
!      open(unit=13,file="lapack_eig.err")
      IF( INFO.GT.0 ) THEN
         WRITE(13,*)'The algorithm failed to compute eigenvalues.'
!         close(13)
         STOP
      END IF
      end subroutine 
!**********************************************************************
! slove A*X=B
! dimension of A : N by N
! dimension of B : N by 1
! in output, B is the solution X 
      SUBROUTINE lapack_solve(A, N, B)

Cf2py intent(in,out) A
Cf2py intent(in) N
Cf2py intent(in,out) B

      INTEGER          N, NRHS
      PARAMETER        ( NRHS = 1 ) 
! don't delete the ()
      INTEGER          LDA, LDB
!
!     .. Local Scalars ..
      INTEGER          INFO
!
!     .. Local Arrays ..
      INTEGER          IPIV( N )
      DOUBLE PRECISION A( N, N ), B( N, NRHS )

!     .. External Subroutines ..
      EXTERNAL         DGESV

!
      LDA = N
      LDB = N
!
!     Solve the equations A*X = B.
!
      CALL DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!     Check for the exact singularity.
!
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The diagonal element of the triangular factor of A,'
         WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
         WRITE(*,*)'A is singular; the solution could not be computed.'
         STOP
      END IF
      END subroutine







!******************Convert Debye to atomic units*************************
      SUBROUTINE Dipole_in_atomic_unit(Dipole)
Cf2py intent(in,out) Dipole
      	IMPLICIT NONE
      	DOUBLE PRECISION :: Dipole
      	Dipole=Dipole*0.3934302014076827
      END SUBROUTINE

!*****************Convert meter to atomic unit****************************
      SUBROUTINE Length_in_Bohr(Lattice_Constant)
Cf2py intent(in,out) Lattice_Constant
      	IMPLICIT NONE
      	DOUBLE PRECISION :: lattice_constant
      	Lattice_Constant=Lattice_Constant*1.889726133921252D10
      END SUBROUTINE	

!*************Convert from energy atomic unit  to kHz**********************
      SUBROUTINE Energy_in_kHz(energy_in_atomic_unit)
Cf2py intent(in,out) energy_in_atomic_unit
      	IMPLICIT NONE
      	DOUBLE PRECISION :: energy_in_atomic_unit
        energy_in_atomic_unit=energy_in_atomic_unit*6.57968392072144D12
      END SUBROUTINE     

!*************Convert from energy atomic unit  to kHz**********************
      SUBROUTINE Energy_in_MHz(energy_in_atomic_unit)
Cf2py intent(in,out) energy_in_atomic_unit
      	IMPLICIT NONE
      	DOUBLE PRECISION :: energy_in_atomic_unit
        energy_in_atomic_unit=energy_in_atomic_unit*6.57968392072144D9
      END SUBROUTINE 
 
!*****************Convert meter to atomic unit****************************
      SUBROUTINE Field_in_atomic_unit(electric_field)
Cf2py intent(in,out) electric_field
      	IMPLICIT NONE
      	DOUBLE PRECISION :: electric_field
      	electric_field=electric_field*1.944690567144141D-12
      END SUBROUTINE	
 

!*****************Convert Hz to atomic unit****************************
      SUBROUTINE Hz_to_atomic_unit(rot_const)
Cf2py intent(in,out) rot_const
      	IMPLICIT NONE
      	DOUBLE PRECISION :: rot_const
      	rot_const=rot_const*2.418884324306202D-17
      END SUBROUTINE	







!*****************find new rotational eigenstates (in terms of old ones)***********
      SUBROUTINE New_EigStates(NewEigValue, coefficient_array,
     &                      Dipole, electric_field,
     &                      rot_const,Number_State,
     &                      N, M) 

Cf2py intent(out) NewEigValue 
Cf2py intent(out) coefficient_array
Cf2py intent(in) Dipole 
Cf2py intent(in) electric_field  
Cf2py intent(in) rot_const 
Cf2py intent(out)number_state
Cf2py intent(in) N
Cf2py intent(in) M
		          
      	IMPLICIT NONE
      	DOUBLE PRECISION :: Dipole, electric_field, rot_const      	
      	INTEGER :: N, M
        DOUBLE PRECISION, ALLOCATABLE :: TotalHamiltonian(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: EigValue(:)
		double precision :: NewEigValue
        DOUBLE PRECISION :: coefficient_array(Number_State)
      	INTEGER :: number_state
        external :: lapack_eig, Matrix_for_new_EigStates 

        coefficient_array = 0.d0 
!assign initial values for the matrix, otherwise it may cause errors in some system
 
        ALLOCATE( TotalHamiltonian(Number_State,Number_State) ) 
       	ALLOCATE( EigValue(Number_State) )
!     .. Form the matrix ..		
       	call Matrix_for_new_EigStates(TotalHamiltonian,Number_State,N,M,
     &                  Dipole,electric_field,rot_const)

!      .... solve the eigenvalue and eigenvector ...	
        call lapack_eig(TotalHamiltonian,number_state, EigValue)

!      assign values for output
        NewEigValue=EigValue(N-ABS(M)+1)
      	coefficient_array(1:Number_State) =
     &            TotalHamiltonian(1:Number_State, N-ABS(M)+1)
       		
       	DEALLOCATE(TotalHamiltonian,EigValue)
      END SUBROUTINE

!*****************alternative for the above subroutine (don't use lapack)***********
      SUBROUTINE New_EigStates2(NewEigValue, coefficient_array,
     &                      Dipole, electric_field,
     &                      rot_const,Number_State,
     &                      N, M)   

Cf2py intent(out) NewEigValue 
Cf2py intent(out) coefficient_array
Cf2py intent(in) Dipole 
Cf2py intent(in) electric_field  
Cf2py intent(in) rot_const 
Cf2py intent(out)number_state
Cf2py intent(in) N
Cf2py intent(in) M
   
      	IMPLICIT NONE
      	DOUBLE PRECISION :: Dipole, electric_field, rot_const
      	INTEGER :: number_state
      	INTEGER :: N, M
      	INTEGER ::  nrot        
!required by Jacobi subroutine 
        DOUBLE PRECISION, ALLOCATABLE :: TotalHamiltonian(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: EigVector(:,:) 
!for recording the eigenvector
        DOUBLE PRECISION, ALLOCATABLE :: EigValue(:)
		double precision :: NewEigValue
        DOUBLE PRECISION :: coefficient_array(100)
        external :: Matrix_for_new_EigStates, jacobi, eigsrt 

       number_state= 50 !100!30
        coefficient_array = 0.d0  !assign initial values for the matrix, otherwise it may cause errors in some system
        

       	ALLOCATE( TotalHamiltonian(Number_State,Number_State) )
       	ALLOCATE( EigValue(Number_State) )
       	ALLOCATE( EigVector(Number_State,Number_State) ) !for recording the eigenvector
		
       	call Matrix_for_new_EigStates(TotalHamiltonian,Number_State,N,M,
     &                  Dipole,electric_field,rot_const)

       	CALL Jacobi(TotalHamiltonian,Number_State,Number_State,
     &               EigValue,EigVector,nrot)

       	CALL EIGSRT(EigValue,EigVector,Number_State,number_state)


        NewEigValue=EigValue(N-ABS(M)+1)
!        write(*,*) NewEigValue
      	coefficient_array(1:Number_State)=EigVector(1:Number_State, N-ABS(M)+1)
       		
       	DEALLOCATE(TotalHamiltonian,EigValue,EigVector)
      END subroutine
  
 	
!***********form matrix for calculating the new rotational eigenstates***************     
      SUBROUTINE Matrix_for_new_EigStates(MatrixName,MatrixDIM,N,M,
     &                      Dipole, electric_field,
     &                      rot_const) 
!NOTE: MatrixDIM > N

Cf2py intent(out) MatrixName
Cf2py intent(in) MatrixDIM
Cf2py intent(in) N
Cf2py intent(in) M
Cf2py intent(in) Dipole 
Cf2py intent(in) Electric_Field  
Cf2py intent(in) rot_const
	 
	    IMPLICIT NONE
        DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        INTEGER :: I, J, K, L, MatrixDIM, N, M,  NStart, NEndWith
        DOUBLE PRECISION,allocatable :: InteractionHamiltonian(:,:) 
!the dimension shouldn't exceed 100
        DOUBLE PRECISION,allocatable :: MoleculeHamiltonian(:,:)  
! default: index starts from 1 to 100
        DOUBLE PRECISION :: MatrixName(1:MatrixDIM,1:MatrixDIM)
        DOUBLE PRECISION :: Dipole, electric_field 
        DOUBLE PRECISION :: rot_const
        double precision, external :: THRJ

!        write(20,*) "Dipole = ", Dipole
!        write(20,*) "Rotational constant =", rot_const  
!        write(20,*) "Electric Field = ", electric_field        
        
		NStart=ABS(M)
	    NEndWith=ABS(M) + MatrixDIM -1
	    
		allocate( InteractionHamiltonian
     &                 (NStart : NEndWith, NStart : NEndWith) )
		allocate( MoleculeHamiltonian
     &                 (NStart : NEndWith, NStart : NEndWith) )
		
	    IF ( N < ABS(M) ) THEN 
	      WRITE(*,*) "Wrong! N should be larger than or equal to M."
	    ELSE IF (MatrixDIM <= N) THEN
	    WRITE(*,*) "Wrong! Matrix Dimension should be larger than N"
	    ELSE       
          DO J=NStart,NEndWith  
!           form the matrix for interaction Hamiltonian
            DO I=NStart,J
              InteractionHamiltonian(I,J)=(-1)*Dipole*(-1)**(M)
     &        *electric_field*DSQRT((2.d0*I+1.d0)*(2.d0*J+1.d0))               
     &        *THRJ(DFLOAT(I),1.D0,DFLOAT(J),-DFLOAT(M),0.D0,DFLOAT(M))
     &        *THRJ(DFLOAT(I),1.D0,DFLOAT(J),0.D0,0.D0,0.D0)
	          InteractionHamiltonian(J,I)=InteractionHamiltonian(I,J) 
! since the matrix is symmetric
               END DO
            END DO
!            write(*,*) InteractionHamiltonian 
        DO L=ABS(M),NEndWith 
          DO K=ABS(M),NEndWith
            IF (L==K) THEN 
              MoleculeHamiltonian(L,K)=2*Pi*rot_const*(K+1)*K
            ELSE 
			  MoleculeHamiltonian(L,K)=0.D0   
            END IF
          END DO
        END DO

	    MatrixName(1:MatrixDIM, 1:MatrixDIM)
     &	    = InteractionHamiltonian(NStart:NEndWith, NStart:NEndWith) 
     &        + MoleculeHamiltonian(NStart:NEndWith, NStart:NEndWith)
	   END IF
      
	  deallocate(InteractionHamiltonian, MoleculeHamiltonian)
	  
	  
	  RETURN
      END SUBROUTINE




!***************************************************************************** 
      SUBROUTINE RotSplit(electric_field,NewEigValue00,NewEigValue10,
     & NewEigValue20,NewEigValue11,NewEigValue21,NewEigValue22)

Cf2py intent(in) electric_field
Cf2py intent(out) NewEigValue00 
Cf2py intent(out) NewEigValue10
Cf2py intent(out) NewEigValue20
Cf2py intent(out) NewEigValue11
Cf2py intent(out) NewEigValue21
Cf2py intent(out) NewEigValue22

      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
      DOUBLE PRECISION :: Dipole, rot_const
      DOUBLE PRECISION :: electric_field
      DOUBLE PRECISION :: NewEigValue00,NewEigValue10, NewEigValue20,
     & NewEigValue11,NewEigValue21,NewEigValue22
      DOUBLE PRECISION :: coefficient_array00(100) 
! the dimension is smaller than 100
      DOUBLE PRECISION :: coefficient_array10(100)
      DOUBLE PRECISION :: coefficient_array20(100)
      DOUBLE PRECISION :: coefficient_array11(100)
      DOUBLE PRECISION :: coefficient_array21(100)
      DOUBLE PRECISION :: coefficient_array22(100)
      INTEGER :: final_state_number00, final_state_number10
      INTEGER :: final_state_number20, final_state_number11
      INTEGER :: final_state_number21, final_state_number22
      external :: New_EigStates 

      Dipole=5.529D0*0.393430307D0
      rot_const=(11.7998D9/2.D0)*2.418884324306202D-17
      write(20,*) "+++++++++++++++++++++++++++++++++++++++++++++++"
      write(20,*) "Dipole = ", Dipole
      write(20,*) "Rotational constant =", rot_const
      write(20,*) "Electric Field = ", electric_field
      write(20,*) "+++++++++++++++++++++++++++++++++++++++++++++++"
      coefficient_array00=0.d0
      coefficient_array10=0.d0
      coefficient_array20=0.d0
      coefficient_array11=0.d0
      coefficient_array21=0.d0
      coefficient_array22=0.d0
      call Field_in_atomic_unit(electric_field)
      CALL New_EigStates(NewEigValue00,coefficient_array00,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number00,
     &                      0, 0)
      CALL New_EigStates(NewEigValue10,coefficient_array10,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number10,
     &                      1, 0)
      CALL New_EigStates(NewEigValue20,coefficient_array20,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number20,
     &                      2, 0)     
      CALL New_EigStates(NewEigValue11,coefficient_array11,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number11,
     &                      1, 1)
      CALL New_EigStates(NewEigValue21,coefficient_array21,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number21,
     &                      2, 1) 
      CALL New_EigStates(NewEigValue22,coefficient_array22,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number22,
     &                      2, 2)     
     
      

       NewEigValue00 = NewEigValue00*219474.6314 
       NewEigValue10 = NewEigValue10*219474.6314
       NewEigValue20 = NewEigValue20*219474.6314 
       NewEigValue11 = NewEigValue11*219474.6314 
       NewEigValue21 = NewEigValue21*219474.6314 
       NewEigValue22 = NewEigValue22*219474.6314 
  
       
      END subroutine
!*************************************************************************  








!*******************calculate the half_sum1 of Cos[mka]/m^3*****************
! k --- wave vector, a --- lattice constant
! num_mole = Number of molecules
      FUNCTION half_sum1(num_mole, ka)
Cf2py intent(in) num_mole
Cf2py intent(in) ka
      	IMPLICIT NONE
      	DOUBLE PRECISION :: half_sum1,ka 
      	INTEGER :: num_mole, I
      	half_sum1=0.D0
      	DO I=1, num_mole/2 
           !do not use DFLOAT(I**3), may exceed the limit for integer 
      		half_sum1=half_sum1 + ( DCOS(I*ka) )/(DFLOAT(I)**3) 
      	END DO
      	RETURN
      END FUNCTION	
!*******************calculate the half_sum1 of 1/m^3*****************
      FUNCTION half_sum2(num_mole)
Cf2py intent(in) num_mole
      	IMPLICIT NONE
      	DOUBLE PRECISION :: half_sum2
      	INTEGER :: num_mole, I
      	half_sum2=0.D0
      	DO I=1, num_mole/2
      		half_sum2=half_sum2 + 1.d0/(DFLOAT(I)**3) 
!do not use DFLOAT(I**3), may exceed the limit for integer
      	END DO
      	RETURN
      END FUNCTION
	
!***********this is for calculating the 2d dispersion curve***************************
      FUNCTION half_sum1_2d(num_mole, kax,kay)
Cf2py intent(in) num_mole
Cf2py intent(in) kax
Cf2py intent(in) kay
      	IMPLICIT NONE
      	DOUBLE PRECISION :: half_sum1_2d,kax, kay 
      	INTEGER :: num_mole, I, J, K, L
      	half_sum1_2d=0.D0
      	DO I=1, num_mole/2 
           !do not use DFLOAT(I**3), may exceed the limit for integer 
      		half_sum1_2d=half_sum1_2d +  DCOS(I*kax)/(DFLOAT(I)**3) 
      	END DO

      	DO L=1, num_mole/2 
           !do not use DFLOAT(I**3), may exceed the limit for integer 
      		half_sum1_2d = half_sum1_2d +  DCOS(L*kay)/(DFLOAT(L)**3) 
      	END DO

! the following is the cross term, related to both kax,kay
        do j=1, num_mole/2
          do k=1, num_mole/2
            half_sum1_2d = half_sum1_2d 
     &          + 2*DCOS(J*kax)*DCOS(K*kay)/( SQRT(DFLOAT(J)**2 + DFLOAT(K)**2) )**3
          end do
        end do
      	RETURN
      END FUNCTION	










	
!***************************************************************************      
      FUNCTION AminusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)

Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_

	    IMPLICIT NONE
      	DOUBLE PRECISION :: AminusBminus
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
        double precision, external :: thrj
      	AminusBminus=THRJ(1.D0,1.D0,2.D0,-1.D0,-1.D0,2.D0)*
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),-1.D0,DFLOAT(MA_))*
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),-1.D0,DFLOAT(MB_))
        RETURN
      END  FUNCTION
 

!****************************************************************************      	      
      FUNCTION AplusBplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)

Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_

	    IMPLICIT NONE
      	DOUBLE PRECISION :: AplusBplus
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
        DOUBLE PRECISION, external :: thrj
      	AplusBplus=THRJ(1.D0,1.D0,2.D0,1.D0,1.D0,-2.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),1.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),1.D0,DFLOAT(MB_))
     
        RETURN
      END  FUNCTION      

!****************************************************************************      
      FUNCTION AplusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)

Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_

	    IMPLICIT NONE
      	DOUBLE PRECISION :: AplusBminus
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
        DOUBLE PRECISION, external :: thrj
      	AplusBminus=THRJ(1.D0,1.D0,2.D0,1.D0,-1.D0,0.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),1.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),-1.D0,DFLOAT(MB_))
     
        RETURN
      END  FUNCTION      

!****************************************************************************      
      FUNCTION AminusBplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)

Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_

	    IMPLICIT NONE
      	DOUBLE PRECISION :: AminusBplus
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
        DOUBLE PRECISION, external :: thrj
      	AminusBplus=THRJ(1.D0,1.D0,2.D0,-1.D0,1.D0,0.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),-1.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),1.D0,DFLOAT(MB_))
     
        RETURN
      END   FUNCTION     


!****************************************************************************      
      FUNCTION A0B0(NA,MA,NB,MB,NA_,MA_,NB_,MB_)

Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_

	    IMPLICIT NONE
      	DOUBLE PRECISION :: A0B0
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
        DOUBLE PRECISION, external :: thrj
      	A0B0=THRJ(1.D0,1.D0,2.D0,0.D0,0.D0,0.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),0.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),0.D0,DFLOAT(MB_))
        RETURN
      END   FUNCTION   
  
!****************************************************************************      
      FUNCTION A0Bminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)

Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_

	    IMPLICIT NONE
      	DOUBLE PRECISION :: A0Bminus
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
        DOUBLE PRECISION, external :: thrj
      	A0Bminus=THRJ(1.D0,1.D0,2.D0,0.D0,-1.D0,1.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),0.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),-1.D0,DFLOAT(MB_))
     
        RETURN
      END  FUNCTION      

!****************************************************************************      
      FUNCTION A0Bplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)

Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_

	    IMPLICIT NONE
      	DOUBLE PRECISION :: A0Bplus
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
        DOUBLE PRECISION, external :: thrj
      	A0Bplus=THRJ(1.D0,1.D0,2.D0,0.D0,1.D0,-1.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),0.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),1.D0,DFLOAT(MB_))
     
        RETURN
      END  FUNCTION      
      
!****************************************************************************      
      FUNCTION AminusB0(NA,MA,NB,MB,NA_,MA_,NB_,MB_)

Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_

	    IMPLICIT NONE
      	DOUBLE PRECISION :: AminusB0
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
        DOUBLE PRECISION, external :: thrj
      	AminusB0=THRJ(1.D0,1.D0,2.D0,-1.D0,0.D0,1.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),-1.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),0.D0,DFLOAT(MB_))
     
        RETURN
      END  FUNCTION 

!****************************************************************************      
      FUNCTION AplusB0(NA,MA,NB,MB,NA_,MA_,NB_,MB_)

Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_

	    IMPLICIT NONE
      	DOUBLE PRECISION :: AplusB0
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
        DOUBLE PRECISION, external :: thrj
      	AplusB0=THRJ(1.D0,1.D0,2.D0,1.D0,0.D0,-1.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),1.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),0.D0,DFLOAT(MB_))
     
        RETURN
      END  FUNCTION 

!*******************************dipole-dipole interaction************************     
      FUNCTION DDInteraction_angle(Dipole,lattice_constant,
     &                    NA,MA,NB,MB,NA_,MA_,NB_,MB_,theta,phi)

Cf2py intent(in) Dipole
Cf2py intent(in) lattice_constant
Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_
Cf2py intent(in) theta
Cf2py intent(in) phi

	    IMPLICIT NONE
        DOUBLE complex :: DDInteraction_angle
      	DOUBLE PRECISION, EXTERNAL :: AminusBminus, AplusBplus,thrj, 
     & 	                              AplusBminus, AminusBplus, A0B0,
     &                                AplusB0,AminusB0,A0Bplus,A0Bminus
        INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_

        DOUBLE PRECISION :: Dipole, Lattice_Constant,prefactor  
        double precision :: theta, phi 
        double complex :: imag_unit
        
        imag_unit = (0.D0,1.D0)                          
       
!        write(*,*) theta, phi
!        write(*,*) (DSIN(theta)**2) *exp(-2.D0*phi*imag_unit)*AminusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)
!        write(*,*) AminusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_) 
        prefactor=-DSQRT(5.D0)*( (-1)**(MA+MB) )        
     & *( (Dipole*Dipole) / (2.D0*Lattice_Constant**3) )      
     & *DSQRT( (2.D0*NA+1.D0)*(2.D0*NA_+1.D0)*(2.D0*NB+1.D0)*(2.D0*NB_+1.D0) )     
     & *THRJ(DFLOAT(NA), 1.D0,DFLOAT(NA_),0.D0,0.D0,0.D0)     
     & *THRJ(DFLOAT(NB), 1.D0,DFLOAT(NB_),0.D0,0.D0,0.D0) 

        DDInteraction_angle = prefactor    
     & *(   3.D0* (DSIN(theta)**2) *exp(-2.D0*phi*imag_unit)*AminusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)      
     &    -6.D0*dsin(theta)*dcos(theta)*exp(-phi*imag_unit)*
     &         ( A0Bminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_) + AminusB0(NA,MA,NB,MB,NA_,MA_,NB_,MB_) )    
     &    + DSQRT(6.D0)* (3*dcos(theta)**2 -1) *(  AplusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)     
     &         + AminusBplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_) +  A0B0(NA,MA,NB,MB,NA_,MA_,NB_,MB_)  )  
     &    + 6.D0*dsin(theta)*dcos(theta)*exp(phi*imag_unit)*
     &         ( A0Bplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_) + AplusB0(NA,MA,NB,MB,NA_,MA_,NB_,MB_) )
     &    + 3.D0*dsin(theta)**2*exp(2.D0*phi*imag_unit)*AplusBplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_) 
     &     )
!       write(*,*) prefactor*DSQRT(6.D0)* (3*dcos(theta)**2 -1)*A0B0(NA,MA,NB,MB,NA_,MA_,NB_,MB_)   
      
        RETURN
      END  FUNCTION

!*******************************dipole-dipole interaction************************     
      FUNCTION DDInteraction_angle2(Dipole,lattice_constant,
     &                    NA,MA,NB,MB,NA_,MA_,NB_,MB_,theta)

Cf2py intent(in) Dipole
Cf2py intent(in) lattice_constant
Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_
Cf2py intent(in) theta

	    IMPLICIT NONE
        DOUBLE precision :: DDInteraction_angle2
      	DOUBLE PRECISION, EXTERNAL :: AminusBminus, AplusBplus,thrj, 
     & 	                              AplusBminus, AminusBplus, A0B0,
     &                                AplusB0,AminusB0,A0Bplus,A0Bminus
        INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_

        DOUBLE PRECISION :: Dipole, Lattice_Constant,prefactor  
        double precision :: theta 
                         
       
!        write(*,*) theta, phi
!        write(*,*) (DSIN(theta)**2) *exp(-2.D0*phi*imag_unit)*AminusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)
!        write(*,*) AminusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_) 
        prefactor=-DSQRT(5.D0)*( (-1)**(MA+MB) )        
     & *( (Dipole*Dipole) / (2.D0*Lattice_Constant**3) )      
     & *DSQRT( (2.D0*NA+1.D0)*(2.D0*NA_+1.D0)*(2.D0*NB+1.D0)*(2.D0*NB_+1.D0) )     
     & *THRJ(DFLOAT(NA), 1.D0,DFLOAT(NA_),0.D0,0.D0,0.D0)     
     & *THRJ(DFLOAT(NB), 1.D0,DFLOAT(NB_),0.D0,0.D0,0.D0) 

        DDInteraction_angle2 = prefactor *       
     &     DSQRT(6.D0)* (3*dcos(theta)**2 -1) *A0B0(NA,MA,NB,MB,NA_,MA_,NB_,MB_) 
!       write(*,*) prefactor*DSQRT(6.D0)* (3*dcos(theta)**2 -1)*A0B0(NA,MA,NB,MB,NA_,MA_,NB_,MB_)   
      
        RETURN
      END  FUNCTION






!***************************calculate matrix elements in electric field*********************
      FUNCTION DDInteract_in_Field_angle(Dipole,electric_field,rot_const,Lattice_Constant,
     &                             NA1,MA1,NB1,MB1,NA2,MA2,NB2,MB2,theta,phi)  
!<B-ground-A-excited|V|A-ground-B-excited>

Cf2py intent(in) Dipole
Cf2py intent(in) electric_field
Cf2py intent(in) rot_const
Cf2py intent(in) Lattice_Constant
Cf2py intent(in) NA1
Cf2py intent(in) MA1
Cf2py intent(in) NB1
Cf2py intent(in) MB1
Cf2py intent(in) NA2
Cf2py intent(in) MA2
Cf2py intent(in) NB2
Cf2py intent(in) MB2
Cf2py intent(in) theta
Cf2py intent(in) phi

	    IMPLICIT NONE
		INTEGER, PARAMETER :: number_state= 50 !100 !30
      	INTEGER :: NA1, MA1, NB1, MB1, NA2,MA2,NB2,MB2
      	INTEGER :: I, J, K, L
      	INTEGER :: final_state_number1, final_state_number2
      	INTEGER :: final_state_number3, final_state_number4      	      	
      	DOUBLE PRECISION :: Dipole, electric_field, rot_const 
        double precision :: Lattice_Constant
      	DOUBLE complex :: DDInteract_in_Field_angle
        DOUBLE PRECISION :: coefficient_array1(number_state)
        DOUBLE PRECISION :: coefficient_array2(number_state)
        DOUBLE PRECISION :: coefficient_array3(number_state)
        DOUBLE PRECISION :: coefficient_array4(number_state)
        DOUBLE PRECISION :: NewEigValue1,NewEigValue2,NewEigValue3,
     &                       NewEigValue4
        double precision :: theta, phi
      	DOUBLE complex, EXTERNAL :: DDInteraction_angle
        external :: New_EigStates 
		
		final_state_number1=number_state
		final_state_number2=number_state
		final_state_number3=number_state
		final_state_number4=number_state
		
		
        CALL New_EigStates(NewEigValue1,coefficient_array1,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number1,
     &                      NA1, MA1)
     
!        Write(40,*) coefficient_array1(1:final_state_number1)
        
        CALL New_EigStates(NewEigValue2,coefficient_array2,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number2,
     &                      NB1, MB1)
     
        CALL New_EigStates(NewEigValue3,coefficient_array3,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number3,
     &                      NA2, MA2)
    
        
        CALL New_EigStates(NewEigValue4,coefficient_array4,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number4,
     &                      NB2, MB2)     
        DDInteract_in_Field_angle=(0.D0,0.D0)
        DO I=1, final_state_number1
          DO J=1, final_state_number2
            DO K=1, final_state_number3
              DO L=1, final_state_number4
        		DDInteract_in_Field_angle=DDInteract_in_Field_angle + 
     &               coefficient_array1(I)*coefficient_array2(J)
     &   	        *coefficient_array3(K)*coefficient_array4(L)	
     &              *DDInteraction_angle( Dipole, Lattice_Constant,
     &                         ABS(MA1)+I-1,MA1,ABS(MB1)+J-1,MB1,    
     &                         ABS(MA2)+K-1,MA2,ABS(MB2)+L-1,MB2,theta,phi )
                END DO
              END DO   
            END DO
          END DO
        
        RETURN
      END  FUNCTION
!*******************************************************************************************


!***************************calculate matrix elements in electric field*********************
      FUNCTION DDInteract_in_Field_angle2(Dipole,electric_field,rot_const,Lattice_Constant,
     &                             NA1,MA1,NB1,MB1,NA2,MA2,NB2,MB2,theta)  
!<B-ground-A-excited|V|A-ground-B-excited>

Cf2py intent(in) Dipole
Cf2py intent(in) electric_field
Cf2py intent(in) rot_const
Cf2py intent(in) Lattice_Constant
Cf2py intent(in) NA1
Cf2py intent(in) MA1
Cf2py intent(in) NB1
Cf2py intent(in) MB1
Cf2py intent(in) NA2
Cf2py intent(in) MA2
Cf2py intent(in) NB2
Cf2py intent(in) MB2
Cf2py intent(in) theta
Cf2py intent(in) phi

	    IMPLICIT NONE
		INTEGER, PARAMETER :: number_state= 50 !100 !30
      	INTEGER :: NA1, MA1, NB1, MB1, NA2,MA2,NB2,MB2
      	INTEGER :: I, J, K, L
      	INTEGER :: final_state_number1, final_state_number2
      	INTEGER :: final_state_number3, final_state_number4      	      	
      	DOUBLE PRECISION :: Dipole, electric_field, rot_const 
        double precision :: Lattice_Constant
      	DOUBLE precision :: DDInteract_in_Field_angle2
        DOUBLE PRECISION :: coefficient_array1(number_state)
        DOUBLE PRECISION :: coefficient_array2(number_state)
        DOUBLE PRECISION :: coefficient_array3(number_state)
        DOUBLE PRECISION :: coefficient_array4(number_state)
        DOUBLE PRECISION :: NewEigValue1,NewEigValue2,NewEigValue3,
     &                       NewEigValue4
        double precision :: theta
      	DOUBLE precision, EXTERNAL :: DDInteraction_angle2
        external :: New_EigStates 
		
		final_state_number1=number_state
		final_state_number2=number_state
		final_state_number3=number_state
		final_state_number4=number_state
		
		
        CALL New_EigStates(NewEigValue1,coefficient_array1,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number1,
     &                      NA1, MA1)
     
!        Write(40,*) coefficient_array1(1:final_state_number1)
        
        CALL New_EigStates(NewEigValue2,coefficient_array2,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number2,
     &                      NB1, MB1)
     
        CALL New_EigStates(NewEigValue3,coefficient_array3,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number3,
     &                      NA2, MA2)
    
        
        CALL New_EigStates(NewEigValue4,coefficient_array4,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number4,
     &                      NB2, MB2)     
        DDInteract_in_Field_angle2=0.D0
        DO I=1, final_state_number1
          DO J=1, final_state_number2
            DO K=1, final_state_number3
              DO L=1, final_state_number4
        		DDInteract_in_Field_angle2=DDInteract_in_Field_angle2 + 
     &               coefficient_array1(I)*coefficient_array2(J)
     &   	        *coefficient_array3(K)*coefficient_array4(L)	
     &              *DDInteraction_angle2( Dipole, Lattice_Constant,
     &                         ABS(MA1)+I-1,MA1,ABS(MB1)+J-1,MB1,    
     &                         ABS(MA2)+K-1,MA2,ABS(MB2)+L-1,MB2,theta)
                END DO
              END DO   
            END DO
          END DO
        
        RETURN
      END  FUNCTION
!*******************************************************************************************


!*******************************************************************************************
      function dipole_dipole_1001(dipole, electric_field, rot_const, lattice_constant)

Cf2py intent(in) Dipole
Cf2py intent(in) electric_field
Cf2py intent(in) rot_const
Cf2py intent(in) Lattice_Constant
Cf2py intent(in) NA1
Cf2py intent(in) MA1
Cf2py intent(in) NB1
Cf2py intent(in) MB1
Cf2py intent(in) NA2
Cf2py intent(in) MA2
Cf2py intent(in) NB2
Cf2py intent(in) MB2

	    IMPLICIT NONE
		INTEGER, PARAMETER :: number_state= 50 !100 !30
      	INTEGER :: NA1, MA1, NB1, MB1, NA2,MA2,NB2,MB2
      	INTEGER :: I, J, K, L
      	INTEGER :: final_state_number1, final_state_number2
      	INTEGER :: final_state_number3, final_state_number4      	      	
      	DOUBLE PRECISION :: Dipole, electric_field, rot_const 
        double precision :: Lattice_Constant
      	DOUBLE precision :: dipole_dipole_1001
        DOUBLE PRECISION :: coefficient_array1(number_state)
        DOUBLE PRECISION :: coefficient_array2(number_state)
        DOUBLE PRECISION :: coefficient_array3(number_state)
        DOUBLE PRECISION :: coefficient_array4(number_state)
        DOUBLE PRECISION :: NewEigValue1,NewEigValue2,NewEigValue3,NewEigValue4
        double precision :: prefactor
      	DOUBLE precision, EXTERNAL :: thrj
        external :: New_EigStates 
		
! because the matrix element here is <10|<00| V |00>|10>
        NA1 = 1
        MA1 = 0
        NA2 = 0
        MA2 = 0

        NB1 = 0
        MB1 = 0
        NB2 = 1
        MB2 = 0        

		final_state_number1=number_state
		final_state_number2=number_state
		final_state_number3=number_state
		final_state_number4=number_state
		
		
        CALL New_EigStates(NewEigValue1,coefficient_array1,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number1,
     &                      NA1, MA1)                          ! I
     

        
        CALL New_EigStates(NewEigValue2,coefficient_array2,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number2,
     &                      NB1, MB1)                          ! J
     
        CALL New_EigStates(NewEigValue3,coefficient_array3,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number3,
     &                      NA2, MA2)                          ! K
    
        
        CALL New_EigStates(NewEigValue4,coefficient_array4,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number4,
     &                      NB2, MB2)                          ! L

        dipole_dipole_1001 = 0.D0
        prefactor = dipole**2/(lattice_constant**3)

        DO I=1, final_state_number1
          DO J=1, final_state_number2
            DO K=1, final_state_number3
              DO L=1, final_state_number4
                dipole_dipole_1001 = dipole_dipole_1001  
     &                 + dsqrt( (2.D0*(I-1) + 1)*(2.D0*(J-1) + 1)*(2.D0*(K-1) + 1)*(2.D0*(L-1) + 1) ) 
!                         sqrt can only accept real variable (don't use integer)
     &                 *coefficient_array1(I)*coefficient_array2(J)
     &   	           *coefficient_array3(K)*coefficient_array4(L)
     &                 *( thrj( 1.D0*(I-1) ,1.D0,1.D0*(K -1),0.D0,0.D0,0.D0)**2 ) 
     &                 *( thrj( 1.D0*(J-1) ,1.D0,1.D0*(L -1),0.D0,0.D0,0.D0)**2 ) 

                END DO
              END DO   
            END DO
          END DO
        dipole_dipole_1001 = prefactor*dipole_dipole_1001
        RETURN
      END  FUNCTION
!*******************************************************************************************








!*******************************************************************************************
      function dipole_dipole_1001_angle(dipole, electric_field, rot_const, lattice_constant,theta)
Cf2py intent(in) Dipole
Cf2py intent(in) electric_field
Cf2py intent(in) rot_const
Cf2py intent(in) Lattice_Constant
Cf2py intent(in) theta

        implicit none
        double precision :: dipole_dipole_1001_angle
        double precision :: dipole, electric_field, rot_const, lattice_constant, theta
        double precision,save :: dd
        integer :: counter
        double precision, external :: dipole_dipole_1001
        data counter / 0 /
        
        if (counter .NE. 1) then
          dd = dipole_dipole_1001(dipole, electric_field, rot_const, lattice_constant)
          counter = 1
        endif

        dipole_dipole_1001_angle = (1 - 3*dcos(theta)**2)*dd
        return
      end function dipole_dipole_1001_angle


      






      
!*******************************dipole-dipole interaction************************     
      FUNCTION DDInteraction(Dipole,lattice_constant,
     &                    NA,MA,NB,MB,NA_,MA_,NB_,MB_)

Cf2py intent(in) Dipole
Cf2py intent(in) lattice_constant
Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_

	    IMPLICIT NONE
        DOUBLE PRECISION :: DDInteraction
      	DOUBLE PRECISION, EXTERNAL :: AminusBminus, AplusBplus,thrj, 
     & 	                              AplusBminus, AminusBplus, A0B0
        INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_

        DOUBLE PRECISION :: Dipole, Lattice_Constant,prefactor                             
       
        
        prefactor=-3.D0*DSQRT(5.D0)*( (-1)**((-1)*(MA+MB)) )        
     & *( (Dipole*Dipole) / (Lattice_Constant)**(3) )      
     & *DSQRT( (2.D0*NA+1.D0)*(2.D0*NA_+1.D0)*(2.D0*NB+1.D0)*(2.D0*NB_+1.D0) )     
     & *THRJ(DFLOAT(NA), 1.D0,DFLOAT(NA_),0.D0,0.D0,0.D0)     
     & *THRJ(DFLOAT(NB), 1.D0,DFLOAT(NB_),0.D0,0.D0,0.D0) 

        DDInteraction = prefactor    
     & *( 0.5*AplusBplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)      
     &    + 0.5*AminusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)     
     &    - DSQRT(1.D0/6.D0)*AplusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)     
     &    - DSQRT(1.D0/6.D0)*AminusBplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)     
     &    - DSQRT(1.D0/6.D0)*A0B0(NA,MA,NB,MB,NA_,MA_,NB_,MB_)  )
      
        RETURN
      END  FUNCTION

!***************************calculate matrix elements in electric field*********************
      FUNCTION DDInteract_in_Field(Dipole,electric_field,rot_const,Lattice_Constant,
     &                             NA1,MA1,NB1,MB1,NA2,MA2,NB2,MB2)  
!<B-ground-A-excited|V|A-ground-B-excited>

Cf2py intent(in) Dipole
Cf2py intent(in) electric_field
Cf2py intent(in) rot_const
Cf2py intent(in) Lattice_Constant
Cf2py intent(in) NA1
Cf2py intent(in) MA1
Cf2py intent(in) NB1
Cf2py intent(in) MB1
Cf2py intent(in) NA2
Cf2py intent(in) MA2
Cf2py intent(in) NB2
Cf2py intent(in) MB2

	    IMPLICIT NONE
		INTEGER, PARAMETER :: number_state=50 !100 !30
      	INTEGER :: NA1, MA1, NB1, MB1, NA2,MA2,NB2,MB2
      	INTEGER :: I, J, K, L
      	INTEGER :: final_state_number1, final_state_number2
      	INTEGER :: final_state_number3, final_state_number4      	      	
      	DOUBLE PRECISION :: Dipole, electric_field, rot_const 
        double precision :: Lattice_Constant
      	DOUBLE PRECISION :: DDInteract_in_Field
        DOUBLE PRECISION :: coefficient_array1(number_state)
        DOUBLE PRECISION :: coefficient_array2(number_state)
        DOUBLE PRECISION :: coefficient_array3(number_state)
        DOUBLE PRECISION :: coefficient_array4(number_state)
        DOUBLE PRECISION :: NewEigValue1,NewEigValue2,NewEigValue3,
     &                       NewEigValue4
      	DOUBLE PRECISION, EXTERNAL :: DDInteraction
        external :: New_EigStates 
		
		final_state_number1=number_state
		final_state_number2=number_state
		final_state_number3=number_state
		final_state_number4=number_state
		
		
        CALL New_EigStates(NewEigValue1,coefficient_array1,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number1,
     &                      NA1, MA1)
     
!        Write(40,*) coefficient_array1(1:final_state_number1)
        
        CALL New_EigStates(NewEigValue2,coefficient_array2,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number2,
     &                      NB1, MB1)
     
        CALL New_EigStates(NewEigValue3,coefficient_array3,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number3,
     &                      NA2, MA2)
    
        
        CALL New_EigStates(NewEigValue4,coefficient_array4,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number4,
     &                      NB2, MB2)     
        DDInteract_in_Field=0.D0
        DO I=1, final_state_number1
          DO J=1, final_state_number2
            DO K=1, final_state_number3
              DO L=1, final_state_number4
        		DDInteract_in_Field=DDInteract_in_Field + 
     &               coefficient_array1(I)*coefficient_array2(J)
     &   	        *coefficient_array3(K)*coefficient_array4(L)	
     &              *DDInteraction( Dipole, Lattice_Constant,
     &                         ABS(MA1)+I-1,MA1,ABS(MB1)+J-1,MB1,    
     &                         ABS(MA2)+K-1,MA2,ABS(MB2)+L-1,MB2 )
                END DO
              END DO   
            END DO
          END DO
        
        RETURN
      END  FUNCTION
!*******************************************************************************************






!*************************************************************************
! the unit of exciton10_energy is kHz
! this is the dispersion energy
! exciton10_energy is an array of excitonic energies for different ka
      Subroutine Exciton10_dispersion(electric_field, num_mole,exciton10_energy)

Cf2py intent(in) electric_field
Cf2py intent(in) num_mole
Cf2py intent(out) exciton10_energy

	    IMPLICIT NONE
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        double precision :: interaction_energy
        double precision :: exciton10_energy(num_mole/2 + 1)  
      	INTEGER ::  I 
      	INTEGER :: num_mole
      	DOUBLE PRECISION :: rot_const,electric_field
        double precision :: summation(num_mole/2 + 1)
      	DOUBLE PRECISION :: Dipole, Lattice_Constant
        double precision :: ka(num_mole/2 + 1)
        external :: Dipole_in_atomic_unit, Length_in_Bohr,
     &              Hz_to_atomic_unit, Energy_in_kHz
        double precision, external :: half_sum1
        double precision, external :: DDInteract_in_Field

      	rot_const=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	Dipole=5.529D0 
        Lattice_Constant=4.D-7         
        
        CALL Dipole_in_atomic_unit(Dipole)
        CALL Length_in_Bohr(Lattice_Constant)
        CALL Hz_to_atomic_unit(rot_const)

        interaction_energy=DDInteract_in_Field(Dipole,electric_field,rot_const,Lattice_Constant,1,0,0,0,0,0,1,0)                
        
        call Energy_in_kHz(interaction_energy) !convert atomic unit to kHz 
        do i=1, num_mole/2 + 1
           ka(i) = (i-1)*2*pi/num_mole
          summation(i)=2*half_sum1(num_mole, ka(i)) 
        end do

        exciton10_energy=interaction_energy*summation 
              
      END SUBROUTINE


!*************************************************************************
! the unit of exciton10_energy is kHz
! this is the dispersion energy
! exciton10_energy is an array of excitonic energies for different ka
      Subroutine Exciton10_dispersion_angle(electric_field, num_mole,exciton10_energy,theta)

Cf2py intent(in) electric_field
Cf2py intent(in) num_mole
Cf2py intent(out) exciton10_energy
Cf2py intent(in)  theta

	    IMPLICIT NONE
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        double precision :: interaction_energy
        double precision :: exciton10_energy(num_mole/2 + 1)
        double precision :: theta  
      	INTEGER ::  I 
      	INTEGER :: num_mole
      	DOUBLE PRECISION :: rot_const,electric_field
        double precision :: summation(num_mole/2 + 1)
      	DOUBLE PRECISION :: Dipole, Lattice_Constant
        double precision :: ka(num_mole/2 + 1)
        external :: Dipole_in_atomic_unit, Length_in_Bohr,
     &              Hz_to_atomic_unit, Energy_in_kHz
        double precision, external :: half_sum1
        double precision, external :: dipole_dipole_1001_angle


      	rot_const=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	Dipole=5.529D0 
        Lattice_Constant=4.D-7         
        
        CALL Dipole_in_atomic_unit(Dipole)
        CALL Length_in_Bohr(Lattice_Constant)
        CALL Hz_to_atomic_unit(rot_const)

        interaction_energy = dipole_dipole_1001_angle(dipole, electric_field, rot_const, lattice_constant,theta)                  
        call Energy_in_kHz(interaction_energy) !convert atomic unit to kHz 
 

        do i=1, num_mole/2 + 1
          ka(i) = (i-1)*2*pi/num_mole
          summation(i)=2*half_sum1(num_mole, ka(i)) 
        end do
 

        exciton10_energy= interaction_energy*summation 
              
      END SUBROUTINE









!*************************************************************************
! the unit of exciton10_energy is kHz
! this is the dispersion energy
! exciton10_energy is an array of excitonic energies for different ka
      subroutine Exciton10_dispersion_2d(electric_field, num_mole)

Cf2py intent(in) electric_field
Cf2py intent(in) num_mole


	    IMPLICIT NONE
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer, parameter :: fileid = 10
        double precision :: interaction_energy
        double precision :: kax, kay,delta_ka  
      	INTEGER ::  I, J 
      	INTEGER :: num_mole
      	DOUBLE PRECISION :: rot_const,electric_field
      	DOUBLE PRECISION :: Dipole, Lattice_Constant
        character(len=20) :: filename
        logical :: alive
        external :: Dipole_in_atomic_unit, Length_in_Bohr,
     &              Hz_to_atomic_unit, Energy_in_kHz
        double precision, external :: half_sum1_2d
        double precision, external :: DDInteract_in_Field
        
      	rot_const=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	Dipole=5.529D0 
        Lattice_Constant=4.D-7 
        kax = 0.D0
        kay = 0.D0        
        delta_ka = 2*pi/num_mole        

        CALL Dipole_in_atomic_unit(Dipole)
        CALL Length_in_Bohr(Lattice_Constant)
        CALL Hz_to_atomic_unit(rot_const)

        interaction_energy=DDInteract_in_Field(Dipole,electric_field,rot_const,Lattice_Constant,1,0,0,0,0,0,1,0)                
        
        call Energy_in_kHz(interaction_energy) !convert atomic unit to kHz 

        write(*,*) "filename: "
        read(*,"(A20)") filename
        inquire(file=filename, exist=alive)
        if (alive) then
          open(unit=fileid, file=filename)
          close(fileid, status="DELETE")
        endif
        open(unit=fileid, file=filename) 
        do i = 1, num_mole/2 + 1
          kax = (i -1)*delta_ka 
          do j=1, num_mole/2 + 1
             kay = (j-1)*delta_ka
             write(unit=fileid, FMT=*) kax, kay, interaction_energy*2*half_sum1_2d(num_mole,kax,kay) 
          end do
          write(unit=fileid,FMT=*) "      " 
        end do 
        close(unit=fileid)
      END SUBROUTINE



!*********************************************************************************************************
      subroutine Exciton10_dispersion_2d_angle(electric_field, num_mole, theta, phi)

Cf2py intent(in) electric_field
Cf2py intent(in) num_mole
Cf2py intent(in) theta
Cf2py intent(in) phi

	    IMPLICIT NONE
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer, parameter :: fileid = 10
        double precision :: interaction_energy
        double precision :: kax, kay,delta_ka  
      	INTEGER ::  I, J 
      	INTEGER :: num_mole ! in the current case, it is 101. It is the number of molecules along one side of the crystal
        integer :: xmax,xmin,ymax,ymin
        integer :: xindex,yindex
        double precision :: theta, phi
      	DOUBLE PRECISION :: rot_const,electric_field
      	DOUBLE PRECISION :: Dipole, Lattice_Constant
        double precision :: distance
        double complex :: summation
        double precision :: cos_alpha, sin_alpha
        character(len=30) :: filename
        logical :: alive
        external :: Dipole_in_atomic_unit, Length_in_Bohr,
     &              Hz_to_atomic_unit, Energy_in_kHz
!        double precision, external :: half_sum1_2d
        double precision, external :: DDInteract_in_Field
        
      	rot_const=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	Dipole=5.529D0 
        Lattice_Constant=4.D-7 
        kax = 0.D0
        kay = 0.D0        
        delta_ka = 2*pi/num_mole        
        xmax = num_mole/2
        xmin = -xmax
        ymax = num_mole/2
        ymin = -ymax

        CALL Dipole_in_atomic_unit(Dipole)
        CALL Length_in_Bohr(Lattice_Constant)
        CALL Hz_to_atomic_unit(rot_const)

        interaction_energy=DDInteract_in_Field(Dipole,electric_field,rot_const,Lattice_Constant,1,0,0,0,0,0,1,0)                
        
        call Energy_in_kHz(interaction_energy) !convert atomic unit to kHz 

        write(*,*) "filename: "
        read(*,"(A30)") filename
        inquire(file=filename, exist=alive)
        if (alive) then
          open(unit=fileid, file=filename)
          close(fileid, status="DELETE")
        endif
        open(unit=fileid, file=filename) 
        do i = -num_mole/2, num_mole/2
          kax = i*delta_ka 
          do j=-num_mole/2, num_mole/2
             kay = j*delta_ka
             summation = (0.D0, 0.D0)
!---------------------------------------------------------------------------------------------------
             do xindex = xmin, xmax
               do yindex = ymin, ymax
                 if ((xindex==0).and.(yindex==0)) then
                    summation = summation + 0.D0 ! molecule cannot interact with itself
                 else
                   distance = sqrt(DFLOAT(xindex)**2 + DFLOAT(yindex)**2)
                   cos_alpha = ABS(DFLOAT(xindex))/distance
                   sin_alpha = ABS(DFLOAT(yindex))/distance 
                   if ( (xindex>=0 .and. yindex >=0) .or. (xindex<0.and. yindex<0) ) then
                     summation = summation + (  1 - 3*(dcos(theta)**2)*
     &                                       (dcos(phi)*cos_alpha + dsin(phi)*sin_alpha)**2  )*
     &                           exp( (0.D0,1.D0)*(kax*xindex + kay*yindex) )
     &                                     /distance**3
                   else
                     summation = summation + (  1 - 3*(dcos(theta)**2)*
     &                                       (dcos(phi)*cos_alpha - dsin(phi)*sin_alpha)**2  )*
     &                           exp( (0.D0,1.D0)*(kax*xindex + kay*yindex) )
     &                                     /distance**3
                   endif
                 endif
               end do
             end do
!---------------------------------------------------------------------------------------------------
             write(unit=fileid, FMT=*) kax, kay, interaction_energy*REAL(summation) 
          end do
          write(unit=fileid,FMT=*) "      " 
        end do 
        close(unit=fileid)
      END SUBROUTINE




!*************************************************************************
! the output exciton20_energy is in atomic unit
! exciton20_energy is the dispersion energy (purely excitonic energy)
! exciton20_energy is an array of excitonic energies for different ka
      Subroutine Exciton20_dispersion(electric_field,num_mole,exciton20_energy)

Cf2py intent(in) electric_field
Cf2py intent(in) num_mole
Cf2py intent(out) exciton20_energy



	    IMPLICIT NONE
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        double precision :: interaction_energy
        double precision :: exciton20_energy(num_mole/2 + 1)  
      	INTEGER ::  I 
      	INTEGER :: num_mole
      	DOUBLE PRECISION :: rot_const,electric_field
        double precision :: summation(num_mole/2 + 1)
      	DOUBLE PRECISION :: Dipole, Lattice_Constant
        double precision :: ka(num_mole/2 + 1)
        external :: Dipole_in_atomic_unit, Length_in_Bohr,
     &              Hz_to_atomic_unit, Energy_in_kHz 
        double precision, external :: half_sum1
        double precision, external :: DDInteract_in_Field

      	rot_const=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	Dipole=5.529D0 
        Lattice_Constant=4.D-7         
        
        CALL Dipole_in_atomic_unit(Dipole)
        CALL Length_in_Bohr(Lattice_Constant)
        CALL Hz_to_atomic_unit(rot_const)

        
        interaction_energy=DDInteract_in_Field(Dipole,electric_field,rot_const,Lattice_Constant,2,0,0,0,0,0,2,0)                
        call Energy_in_kHz(interaction_energy) !convert atomic unit to kHz
        
        do i=1, num_mole/2 + 1
           ka(i) = (i-1)*2*pi/num_mole
           summation(i)=2*half_sum1(num_mole, ka(i)) 
        end do

       exciton20_energy=interaction_energy*summation 
              
      END SUBROUTINE







!********************calculate the matrix element <f|V|i> as a function of k*******
      Subroutine M_vs_k(electric_field, num_mole,m_k_array) 
!unit of M_vs_k is: kHz  Electric_field :: atomic unit

Cf2py intent(in) electric_field
Cf2py intent(in) num_mole
Cf2py intent(out) m_k_array

        IMPLICIT NONE
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
      	DOUBLE PRECISION :: m_k_array(num_mole/2 + 1)
      	DOUBLE PRECISION :: ka(num_mole/2 + 1)
        double precision :: summation(num_mole/2 + 1)
      	double precision :: matrix_element
      	INTEGER :: num_mole,i
      	DOUBLE PRECISION :: rot_const,electric_field
      	DOUBLE PRECISION :: Dipole, Lattice_Constant
        external :: Dipole_in_atomic_unit, Length_in_Bohr,
     &              Hz_to_atomic_unit, Energy_in_kHz 
        double precision, external :: DDInteract_in_Field
        double precision, external :: half_sum1


      	rot_const=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	Dipole=5.529D0 
        Lattice_Constant=4.D-7 
                   
        CALL Dipole_in_atomic_unit(Dipole)
        CALL Length_in_Bohr(Lattice_Constant)
        CALL Hz_to_atomic_unit(rot_const)
     
      	matrix_element = DDInteract_in_Field(Dipole,electric_field,rot_const,Lattice_Constant,
     &                             2,0,0,0,1,0,1,0)
        CALL Energy_in_kHz(matrix_element)

        do i=1, num_mole/2 + 1
           ka(i) = (i-1)*2*pi/num_mole
           summation(i)=2*half_sum1(num_mole, ka(i)) 
        end do

        m_k_array=matrix_element*summation
        return            
      END subroutine
    
!*******calculate E_{k1} + E_{k2} - E_{k3} excluding the excitonic energy**********
! the return energy is in the unit of kHz
      Function delta_energy_exclude(electric_field, num_mole) 
!unit of electric field: atomic unit 

Cf2py intent(in) electric_field
Cf2py intent(in) num_mole
        
        implicit none
        double precision :: electric_field, delta_energy_exclude
        integer :: num_mole
        double precision,external :: E1G,E2G,GPSA
        external :: Energy_in_kHz

        delta_energy_exclude= 2*( E1G(electric_field) + GPSA(electric_field,num_mole, 1) )  
     &                       - E2G(electric_field) - GPSA(electric_field,num_mole, 2)  
!       here delta_energy is in atomic unit
        CALL Energy_in_kHz(delta_energy_exclude)
        return 
! return energy in kHz
      end function



!*******energy difference between |10> and |00> at certain electric field************
      FUNCTION E1G(electric_field) 
!unit of electric field: atomic units

Cf2py intent(in) electric_field
        
        IMPLICIT NONE
        DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer, parameter :: Number_state = 100 !30
        DOUBLE PRECISION :: Dipole, rot_const
        DOUBLE PRECISION :: electric_field, E1G
        DOUBLE PRECISION :: NewEigValue00,NewEigValue10
        DOUBLE PRECISION :: coefficient_array00(Number_state)
        DOUBLE PRECISION :: coefficient_array10(Number_state)
        INTEGER :: final_state_number00, final_state_number10
        external :: New_EigStates 
      
        final_state_number00 = Number_state 
        final_state_number10 = Number_state 
        Dipole=5.529D0*0.393430307D0 !convert to atomic unit
        rot_const=(11.7998D9/2.D0)*2.418884324306202D-17 !convert to atomic unit
        CALL New_EigStates(NewEigValue00,coefficient_array00,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number00,
     &                      0, 0)
        CALL New_EigStates(NewEigValue10,coefficient_array10,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number10,
     &                      1, 0)
        E1G=NewEigValue10 - NewEigValue00
        return
      END function

!*******energy difference between |20> and |00> at certain electric field************
      FUNCTION E2G(electric_field) 
!unit of electric field: atomic unit

Cf2py intent(in) electric_field
       
        IMPLICIT NONE
        DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer, parameter :: Number_state = 100 !30
        DOUBLE PRECISION :: Dipole, rot_const
        DOUBLE PRECISION :: electric_field, E2G
        DOUBLE PRECISION :: NewEigValue00,NewEigValue20
        DOUBLE PRECISION :: coefficient_array00(Number_state)
        DOUBLE PRECISION :: coefficient_array20(Number_state)
        INTEGER :: final_state_number00, final_state_number20
        external :: New_EigStates

        final_state_number00= Number_state 
        final_state_number20 = Number_state 
        Dipole=5.529D0*0.393430307D0 !convert to atomic unit
        rot_const=(11.7998D9/2.D0)*2.418884324306202D-17 !convert to atomic unit
        CALL New_EigStates(NewEigValue00,coefficient_array00,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number00,
     &                      0, 0)
        CALL New_EigStates(NewEigValue20,coefficient_array20,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number20,
     &                      2, 0)
        E2G=NewEigValue20 - NewEigValue00
        return
      END FUNCTION
!************************************************************************************** 
  
!*********gas phase shift for N=A exciton, A can be 1 or 2***************************
      Function GPSA(electric_field,num_mole, A)

Cf2py intent(in) electric_field
Cf2py intent(in) num_mole
Cf2py intent(in) A

        IMPLICIT NONE
        integer :: num_mole 
        integer :: A
        double precision :: dipole, electric_field, rot_const,Lattice_Constant
        double precision :: GPSA       
        external :: Dipole_in_atomic_unit, Length_in_Bohr,
     &              Hz_to_atomic_unit, Energy_in_kHz 
        double precision, external :: half_sum2
        double precision, external :: DDInteract_in_Field
        

      	rot_const=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	Dipole=5.529D0 
        Lattice_Constant=4.D-7 
           
!***************unit conversion*******************************************        
        CALL Dipole_in_atomic_unit(Dipole)
        CALL Length_in_Bohr(Lattice_Constant)
        CALL Hz_to_atomic_unit(rot_const)
        
        GPSA=(  DDInteract_in_Field( Dipole,electric_field,rot_const,Lattice_Constant,
     &                             A,0,0,0,A,0,0,0 )
     & - DDInteract_in_Field( Dipole,electric_field,rot_const,Lattice_Constant,
     &                             0,0,0,0,0,0,0,0)  )
     &            *2*half_sum2(num_mole)

        return ! unit of GPSA: atomic unit
      end FUNCTION

!*************************Form the evolution Hamiltonian matrix*****************************
      subroutine form_evolution_matrix(electric_field, 
     &            num_mole, evolution_matrix) 

Cf2py intent(in) electric_field
Cf2py intent(in) num_mole
Cf2py intent(out) evolution_matrix  
     
        implicit none
        double precision, parameter :: Pi=3.141592653589793115997963468d0  
        double precision :: electric_field, record
        integer :: num_mole
        double precision :: evolution_matrix(num_mole/2 + 2, num_mole/2 + 2)
        integer :: i, j, K
        double precision ::delta_energy(num_mole/2 + 1)
        double precision :: exciton10_energy(num_mole/2 + 1) 
        double precision :: exciton20_energy(num_mole/2 + 1)
        double precision :: m_k_array(num_mole/2 + 1)
        external :: Exciton10_dispersion, Exciton20_dispersion, M_vs_k 
        double precision, external :: delta_energy_exclude 
        evolution_matrix = 0.D0 

        CALL Exciton10_dispersion(electric_field,num_mole,exciton10_energy) 
        CALL Exciton20_dispersion(electric_field,num_mole,exciton20_energy)
        CALL M_vs_k(electric_field, num_mole,m_k_array)

        record =delta_energy_exclude(electric_field,num_mole) 

        do i=1, num_mole/2 + 1  
          delta_energy(i) = 2*exciton10_energy(i) - exciton20_energy(i) + record
        end do

        do j=2, num_mole/2 + 2 
          evolution_matrix(j,j) = delta_energy(j-1) 
        end do
        
        do K=2,num_mole/2 + 2
           evolution_matrix(1, K) = 2*m_k_array(K-1)/num_mole 
           evolution_matrix(K, 1) = evolution_matrix(1, K) 
        end do

        evolution_matrix(1,2) = evolution_matrix(1,2)/2  ! k=0 case (kronecker delta)
        evolution_matrix(2,1) = evolution_matrix(2,1)/2  ! k=0 case 
      

        evolution_matrix = 2000*Pi*evolution_matrix
      end subroutine



!*****************find the root of func(x) in range (x1,x2)*************************
      FUNCTION zriddr(func,x1,x2,xacc)

Cf2py intent(in) func
Cf2py intent(in) x1
Cf2py intent(in) x2
Cf2py intent(in) xacc

      INTEGER MAXIT
      double precision :: zriddr,x1,x2,xacc,func,UNUSED
      PARAMETER (MAXIT=6000,UNUSED=-1.11E30)
      EXTERNAL func
CU    USES func
      INTEGER j
      double precision :: fh,fl,fm,fnew,s,xh,xl,xm,xnew
      fl=func(x1)
      fh=func(x2)
      if((fl.gt.0..and.fh.lt.0.).or.(fl.lt.0..and.fh.gt.0.))then
        xl=x1
        xh=x2
        zriddr=UNUSED
        do 11 j=1,MAXIT
          xm=0.5*(xl+xh)
          fm=func(xm)
          s=sqrt(fm**2-fl*fh)
          if(s.eq.0.)return
          xnew=xm+(xm-xl)*(sign(1.d0,fl-fh)*fm/s)
          if (abs(xnew-zriddr).le.xacc) return
          zriddr=xnew
          fnew=func(zriddr)
          if (fnew.eq.0.) return
          if(sign(fm,fnew).ne.fm) then
            xl=xm
            fl=fm
            xh=zriddr
            fh=fnew
          else if(sign(fl,fnew).ne.fl) then
            xh=zriddr
            fh=fnew
          else if(sign(fh,fnew).ne.fh) then
            xl=zriddr
            fl=fnew
          else
            pause 'never get here in zriddr'
          endif
          if(abs(xh-xl).le.xacc) return
11      continue
        pause 'zriddr exceed maximum iterations'
      else if (fl.eq.0.) then
        zriddr=x1
      else if (fh.eq.0.) then
        zriddr=x2
      else
        pause 'root must be bracketed in zriddr'
      endif
      return
      END function


!***********************************************************************************
      function func(electric_field)
! accept electric_field in atomic unit
Cf2py intent(in) electric_field
        
        implicit none
        double precision :: func
        double precision :: electric_field
        integer, parameter :: num_mole = 5000
        double precision :: exciton10_energy(num_mole/2 + 1)
        double precision :: exciton20_energy(num_mole/2 + 1)
        external:: Exciton10_dispersion, Exciton20_dispersion
        double precision, external :: delta_energy_exclude       

!        call Field_in_atomic_unit(electric_field) ! don't change variables in function
        call Exciton10_dispersion(electric_field,num_mole, exciton10_energy)
        call Exciton20_dispersion(electric_field,num_mole, exciton20_energy) 
        func = delta_energy_exclude(electric_field, num_mole) 
     &                      + 2*exciton10_energy(1) ! can change the index to let different k in resonance
     &                      - exciton20_energy(1)
        return   
      end function


!*************************************************************************************
      subroutine find_electric_field(start_value, end_value, electric_field)

Cf2py intent(in) start_value
Cf2py intent(in) end_value
Cf2py intent(out) electric_field

        implicit none
        double precision :: electric_field
        double precision :: start_value, end_value ! in atomic unit
        double precision :: xacc 
        double precision, external :: func ! func = 0 is the equation to solve
        double precision, external :: zriddr ! root-finding subroutine

        xacc = 1.D-22
        
        electric_field = zriddr(func,start_value,end_value,xacc)
        
      end subroutine 



!*******************************************************************************************
      subroutine probability_dis(time, eig_vector_matrix, eig_value,initial_coefficients, coefficients_squared)

Cf2py intent(in) time
Cf2py intent(in) eig_vector_matrix
Cf2py intent(in) eig_value
Cf2py intent(in) initial_coefficients
Cf2py intent(out) coefficients_squared

        implicit none
        double precision, parameter :: Pi=3.141592653589793115997963468d0
        integer, parameter :: num_mole = 5000
        integer :: i, j, K
        double precision ::  time
        double precision :: coefficients_squared(num_mole/2+2)
        double precision :: coefficients_real(num_mole/2+2)
        double precision :: coefficients_imaginary(num_mole/2+2)
        double precision :: eig_vector_matrix(num_mole/2+2, num_mole/2 + 2)
        double precision :: eig_vector_matrix_t(num_mole/2+2, num_mole/2 + 2)     
        double precision :: eig_value(num_mole/2+2)
        double precision :: initial_coefficients(num_mole/2+2)
        double precision :: exponential_Dt_real(num_mole/2+2, num_mole/2+2)
        double precision :: exponential_Dt_imaginary(num_mole/2+2, num_mole/2+2)
        double precision :: sum_coefficient
        

        exponential_Dt_real = 0.D0
        exponential_Dt_imaginary = 0.D0 
        sum_coefficient = 0.D0

        eig_vector_matrix_t = TRANSPOSE(eig_vector_matrix) 
        
        do i=1, num_mole/2 + 2
          exponential_Dt_real(i,i) = COS(eig_value(i)*time)
          exponential_Dt_imaginary(i,i) = SIN(eig_value(i)*time)
        end do

        coefficients_real = MATMUL(  eig_vector_matrix, 
     &                         MATMUL( exponential_Dt_real, 
     &         MATMUL(eig_vector_matrix_t, initial_coefficients) )  )

        coefficients_imaginary = MATMUL(  eig_vector_matrix, 
     &                              MATMUL( exponential_Dt_imaginary, 
     &         MATMUL(eig_vector_matrix_t, initial_coefficients) )  )

  
        do j=1,num_mole/2 + 2
          coefficients_squared(j) = (coefficients_real(j))**2 + (coefficients_imaginary(j))**2
          sum_coefficient = sum_coefficient + coefficients_squared(j) 
        end do
! normalization
       do K=1, num_mole/2 + 2 
        coefficients_squared(K)=coefficients_squared(K)/sum_coefficient
       end do
   
      end subroutine probability_dis





!*******************************************************************************************
      subroutine form_matrix_1d(neq,H,T)
Cf2py intent(in) neq
Cf2py intent(in,out) H
Cf2py intent(in) T

        use molecule_information
        implicit none
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer :: neq
        double precision :: H(neq,neq)
        double precision :: T
        double precision, allocatable,save :: H_save(:,:)
        double precision :: prefactor
        double precision :: theta
        integer :: I,J,K
        integer :: counter,error,flag
        double precision, external :: theta_vs_time
        double precision, external :: dipole_dipole_1001
        common /flag/ flag
        data counter / 0 /


        if (counter .NE. 1) then
           allocate(H_save(neq,neq),stat=error)
           if (error/=0) then
             write(*,*) "Allocation of H1 failed"
           endif
           counter = 1
        
           H_save = 0.D0

           do i=1,neq-1
             H_save(1,i+1) = 1.D0/(DFLOAT(i)**3)
           end do
 
           do j=2, neq
             H_save(j,j+1:neq) = H_save(1,2:neq+1-j)
           end do

           do K=1,neq
             H_save(K:neq,K) = H_save(K,K:neq)
           end do 
      
           prefactor = dipole_dipole_1001(dipole,electric_field, rot_const,lattice_constant) 
           call energy_in_kHz(prefactor)
           prefactor = -prefactor*1000*2*pi
           H_save = prefactor*H_save
        endif
        theta = theta_vs_time(T,flag)

        H = (1 - 3*dcos(theta)**2)*H_save 
        
      end subroutine form_matrix_1d 

!********************************one wavepacket**********************************************
      subroutine initial_condition(Y_tmp,neq)
Cf2py intent(in,out) Y_tmp
Cf2py intent(in) neq

        implicit none
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer :: neq
        double complex :: Y_tmp(-neq/2:neq/2)
        double precision :: ka
        double precision :: normalization_constant
        integer :: i,j
        double precision, external :: gaussian_func
        Y_tmp = (0.D0,0.D0)
        do i = -neq/2, neq/2 ! i ---- molecule index
          do j = -neq/2, neq/2
            ka = j*2*pi/neq
            Y_tmp(i) = Y_tmp(i) + gaussian_func(ka)*exp(ka*i*DCMPLX(0.D0,1.D0))
          end do
        end do
!     normalization of the coefficients
        normalization_constant = SUM(ABS(Y_tmp)**2)
        Y_tmp = Y_tmp/normalization_constant
      end subroutine initial_condition


!********************************one wavepacket in 2d***************************************
      subroutine initial_condition_2d(Y_tmp,nx,ny)
Cf2py intent(in,out) Y_tmp
Cf2py intent(in) nx
Cf2py intent(in) ny

        implicit none
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer :: nx, ny
        double complex :: Y_tmp(-nx/2:nx/2, -ny/2:ny/2)
        double precision :: kax, kay
        double precision :: normalization_constant
        integer :: i,j, k, l
        double precision, external :: gaussian_func_2d
        Y_tmp = (0.D0,0.D0)
        do i = -nx/2, nx/2 ! i ---- molecule index
          do j = -ny/2, ny/2
             do k = -nx/2, nx/2
               do L = -ny/2, ny/2
                 kax = K*2*pi/nx
                 kay = L*2*pi/ny
                 Y_tmp(i,j) = Y_tmp(i,j) + gaussian_func_2d(kax,kay)
     &                        *exp( (kax*i + kay*j)*DCMPLX(0.D0,1.D0) )
               end do
             end do
          end do
        end do
!  no need to normalize it 
!     normalization of the coefficients
        normalization_constant = SQRT( SUM(ABS(Y_tmp)**2) )
        Y_tmp = Y_tmp/normalization_constant
      end subroutine initial_condition_2d






!***************************two wavepackets************************************************
      subroutine initial_condition2(Y_tmp,neq)
Cf2py intent(in,out) Y_tmp
Cf2py intent(in) neq

        implicit none
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer :: neq
        double complex :: Y_tmp(-neq/2:neq/2)
        double precision :: ka
        double precision :: normalization_constant
        integer :: i,j
        double complex, external :: twogaussian_func
        Y_tmp = (0.D0,0.D0)
        do i = -neq/2, neq/2 ! i ---- molecule index
          do j = -neq/2, neq/2
            ka = j*2*pi/neq
            Y_tmp(i) = Y_tmp(i) + twogaussian_func(ka)*exp(ka*i*DCMPLX(0.D0,1.D0))
          end do
        end do
!     normalization of the coefficients
        normalization_constant = SUM(ABS(Y_tmp)**2)
        Y_tmp = Y_tmp/normalization_constant
      end subroutine initial_condition2





!*************************************************************************************
      function gaussian_func(ka)
Cf2py intent(in) ka

        implicit none
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        double precision, parameter :: variance = 0.005D0 !0.005D0
        double precision, parameter :: central_ka = pi/2 
        double precision :: ka
        double precision :: gaussian_func

        gaussian_func = (1./sqrt(2*pi*variance))*exp( -(ka- central_ka)**2/(2*variance) )
        return
      end function gaussian_func



!*************************************************************************************
      function twogaussian_func(ka)
Cf2py intent(in) ka

        implicit none
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        double precision, parameter :: variance = 1.D-3 !0.005D0
        double precision, parameter :: central_ka1 = -pi/2 !pi/2 
        double precision, parameter :: central_ka2 = pi/2
        double precision :: ka
        double complex :: twogaussian_func

        twogaussian_func = (1./sqrt(2*pi*variance))*exp( -(ka- central_ka1)**2/(2*variance) )
     &                + (1./sqrt(2*pi*variance))*exp( -(ka- central_ka2)**2/(2*variance) )
!     &                  *exp(DCMPLX(0.D0,1.D0)*pi) 
        return
      end function twogaussian_func



!*************************************************************************************
      function gaussian_func_2d(kax,kay)
Cf2py intent(in) kax
Cf2py intent(in) kay

        implicit none
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        double precision, parameter :: variance = 0.001D0 !0.005D0
        double precision, parameter :: central_kax = 0.D0 !pi/2
        double precision, parameter :: central_kay = 0.D0 !pi/2 
        double precision :: kax, kay
        double precision :: gaussian_func_2d

        gaussian_func_2d =(1./(2*pi*variance))*exp(-(kax- central_kax)**2/(2*variance))
     &                     *exp( -(kay- central_kay)**2/(2*variance) ) 
        return
      end function gaussian_func_2d



!*******************************************************************************************
      subroutine form_matrix_2d(neq,Ham,T)
        use molecule_information
        use random_vacancy
        implicit none
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer :: neq, nx, ny
        double precision :: Ham(neq,neq)
        double precision :: T
        double precision :: H1(neq,neq), H_sin(neq,neq), H_cos(neq,neq)
        double precision, allocatable, save :: H_tmp(:,:)
        integer :: vacancy
        integer :: I,J
        integer :: xi,yi,xj,yj,neq_sqrt
        double precision :: distance
        integer :: error1, error2, error3
        double precision :: theta
        double precision :: phi
        double precision :: alpha
        double precision :: cos_alpha
        double precision :: sin_alpha
        double precision :: prefactor
        double precision :: vacancy_percentage
        integer, save :: counter
        double precision, external :: theta_vs_time
        double precision, external :: phi_vs_time
        double precision, external :: dipole_dipole_1001
        external :: energy_in_kHz
!       flag1 determine how theta would change
!       flag2 determine how phi would change
        integer :: flag1, flag2
        common /flags/ flag1, flag2
        common /vacancy_percentage/ vacancy_percentage
        data counter / 0 /

        !write(*,*) "============================="
        if (counter .NE. 1) then
           counter = 1
!           allocate(H1(neq,neq),stat=error1)
!             if (error1/=0) then
!               write(10,*) "Allocation of H1 failed"
!             endif
!           allocate(H_sin(neq,neq),stat=error2)
!             if (error2/=0) then
!               write(20,*) "Allocation of H_sin failed"
!             endif
!           allocate(H_cos(neq,neq),stat=error3)
!             if (error3/=0) then
!               write(30,*) "Allocation of H_cos failed"
!             endif

          allocate(H_tmp(neq,neq),stat=error1)
          if (error1/=0) then
            write(10,*) "Allocation of H_tmp failed"
            call system("mv fort.10 error.txt")
          endif
          
          
          H1 = 0.D0
          H_sin = 0.D0
          H_cos = 0.D0
          neq_sqrt = int(sqrt(real(neq)))
          nx = neq_sqrt
          ny = neq_sqrt
          
          do i=1,neq-1
            do j=i+1, neq

!              yi = i/neq_sqrt + 1  
!              xi = mod(i,neq_sqrt)
!              if (xi==0) then
!                xi = neq_sqrt
!                yi = i/neq_sqrt
!              endif 
              call map_1d_to_2d(i, nx, ny, xi, yi)
              !call Oned_to_2d(i, nx, ny, xi, yi)
              
!              yj = j/neq_sqrt + 1
!              xj = mod(j,neq_sqrt)
!              if (xj==0) then
!                xj = neq_sqrt
!                yj = j/neq_sqrt
!              endif 
              call map_1d_to_2d(j, nx, ny, xj, yj)
              !call Oned_to_2d(j, nx, ny, xi, yi)
              
              distance = sqrt(real((xi-xj)**2 + real(yi-yj)**2))           
              H1(j,i) = 1.D0/distance**3  
              H1(i,j) = H1(j,i) ! make H1 symmetric
             
              alpha = ATAN2(DFLOAT(yj-yi), DFLOAT(xj-xi))
              H_sin(j,i) = SIN(alpha)
              H_cos(j,i) = COS(alpha)
 
!              sin_alpha = ABS(yi-yj)/distance
!              cos_alpha = ABS(xi-xj)/distance
!              H_sin(j,i) = sin_alpha
!              H_cos(j,i) = cos_alpha

!             make them symmetric 
              H_sin(i,j) = H_sin(j,i)
              H_cos(i,j) = H_cos(j,i)

            
             end do
           end do 

           prefactor = dipole_dipole_1001(dipole, electric_field, rot_const, lattice_constant) 
           call energy_in_kHz(prefactor)
           !prefactor = 22.52 ! kHz
           !prefactor = -prefactor*1000*2*pi
           write(66,*) "J(a) = ", prefactor, "kHz"
           call system("mv fort.66 Ja.txt")
           prefactor = prefactor*2*pi
           H1 = prefactor*H1
 
!        endif
        

          ! it is assumed that theta and phi are independent of time in this calculation
          theta = theta_vs_time(T, flag1) 
          phi = phi_vs_time(T, flag2)
          H_tmp = (1 - 3*(dcos(theta)**2)*(dcos(phi)*H_cos + dsin(phi)*H_sin)**2 )*H1
        
        
          ! consider the vacancy sites by setting all the relevant elements to be zero
          ! if the vacancy_percentage = 0.0, do nothing
          if (ABS(vacancy_percentage) > 1.D-10 ) then
            do i=1, n_vacancy
              vacancy = vacancy_site(i)
              H_tmp(vacancy,:) = 0.D0
              H_tmp(:,vacancy) = 0.D0
            end do
          endif
          
        endif
        
        Ham = H_tmp 
        !write(*,*) "============================="                
        return
      end subroutine form_matrix_2d 





 



!*******************************************************************************************
      function theta_vs_time(time,flag1)
! specify how the angle between the direction of electric field and the intermolecular axis changes as 
! a function of time
Cf2py intent(in) time
Cf2py intent(in) flag1

        implicit none 
        double precision :: time
        integer :: flag1
        double precision :: theta_vs_time
        double precision, parameter :: pi = 3.141592653589793115997963468d0 
        double precision,parameter :: theta_critical = DACOS(1.D0/SQRT(3.D0))
        double precision, parameter :: time_start = 1.D-4
        double precision, parameter :: time_interval = 1.D-4

        select case(flag1)
!--------------------------------------------------------------------------------------------
          case(1)
            theta_vs_time = pi/2
!--------------------------------------------------------------------------------------------
          case(2)
            theta_vs_time = 0.D0
!--------------------------------------------------------------------------------------------
          case(3)
            IF ((time .GE. 0.D0) .AND. (time .LT. time_start) ) then
              theta_vs_time = 0.D0
            else if ((time.GE.time_start).AND.(time .LE.(time_start + time_interval))) then
              theta_vs_time = theta_critical*( DSIN(pi*(time - time_start)/(2*time_interval)) )**2
            else if (time .GT.(time_start + time_interval) ) then
              theta_vs_time = theta_critical
            else
              write(*,*) "the input time is not valid"
              stop 
            end if
!--------------------------------------------------------------------------------------------
          case(4)
            theta_vs_time = DACOS(SQRT(2.D0/3.D0))
!--------------------------------------------------------------------------------------------
          case(5)
            if ((time >= 0.D0) .AND. (time < 5.D-4) ) then
              theta_vs_time = 0.D0
            else if ((time >= 5.D-4).and.(time <= 1.5D-3)) then
              theta_vs_time =theta_critical*( DSIN(pi*(time - 5.D-4)/(2*1.D-3)) )**2
            else if ((time > 1.5D-3).and.(time < 2.D-3)) then
              theta_vs_time = theta_critical
            else if ((time >= 2.D-3).and.(time <= 3.D-3)) then
              theta_vs_time = theta_critical*( DCOS(pi*(time - 2.D-3)/(2*1.D-3)) )**2
            else
              theta_vs_time = 0.D0
            endif 
!-------------------------------------------------------------------------------------------- 
          case(6)
            if ((time >= 0.D0) .AND. (time < 5.D-4) ) then
              theta_vs_time = 0.D0
            else if ((time >= 5.D-4).and.(time <= 1.5D-3)) then
              theta_vs_time =theta_critical*( DSIN(pi*(time - 5.D-4)/(2*1.D-3)) )**2
            else if ((time > 1.5D-3).and.(time < 2.D-3)) then
              theta_vs_time = theta_critical
            else if ((time >= 2.D-3).and.(time <= 3.D-3)) then
              theta_vs_time = (pi/2-theta_critical)*( DSIN(pi*(time - 2.D-3)/(2*1.D-3)) )**2 + theta_critical
            else
              theta_vs_time = pi/2
            endif 
!--------------------------------------------------------------------------------------------
          case(7)
              theta_vs_time = pi/20

          case default
             write(*,*) "Invalid time"
        end select
        return
      end function theta_vs_time   









!*******************************************************************************************
      function phi_vs_time(time,flag2)
! specify how the angle between the direction of electric field and the intermolecular axis changes as 
! a function of time
Cf2py intent(in) time
Cf2py intent(in) flag2

        implicit none 
        double precision :: time
        integer :: flag2
        double precision :: phi_vs_time
        double precision, parameter :: pi = 3.141592653589793115997963468d0 
        double precision,parameter :: critical_angle = DACOS(1.D0/SQRT(3.D0))
        double precision, parameter :: time_start = 1.D-4
        double precision, parameter :: time_interval = 1.D-4
        select case(flag2)
          case(1)
            phi_vs_time = pi/2

          case(2)
            phi_vs_time = 0.D0

          case(3)
            phi_vs_time = critical_angle

          case(4)
            phi_vs_time = pi/4
          case(5)
            phi_vs_time = pi/10
          case default
            write(*,*) "Invalid time"
        end select
        return
      end function phi_vs_time 










