



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












!*****************find new rotational eigenstates (in terms of old ones)***********
! accept the number of rotational levels for the calculation from outside
      SUBROUTINE New_EigStates(NewEigValue, coefficient_array,
     &                      Dipole, electric_field,
     &                      rot_const,Number_State,
     &                      N, M) 

		          
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
      END SUBROUTINE New_EigStates





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
      END SUBROUTINE Matrix_for_new_EigStates






!***********form matrix for calculating the new rotational eigenstates***************     
      SUBROUTINE Matrix_for_new_EigStates_AC(MatrixName,MatrixDIM,N,M,
     &                      Dipole, DC_field, intensity,
     &                      alpha_parallel, alpha_perpendicular,
     &                      rot_const) 
!NOTE: MatrixDIM > N

	 
	    IMPLICIT NONE
        DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        INTEGER :: I, J, K, L, MatrixDIM, N, M,  NStart, NEndWith
        DOUBLE PRECISION,allocatable :: InteractionHamiltonian(:,:) 
!the dimension shouldn't exceed 100
        DOUBLE PRECISION,allocatable :: MoleculeHamiltonian(:,:)  
! default: index starts from 1 to 100
        DOUBLE PRECISION :: MatrixName(1:MatrixDIM,1:MatrixDIM)
        DOUBLE PRECISION :: Dipole, DC_field 
        double precision :: intensity
        double precision :: alpha_parallel, alpha_perpendicular
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
     &        *DC_field*DSQRT((2.d0*I+1.d0)*(2.d0*J+1.d0))               
     &        *THRJ(DFLOAT(I),1.D0,DFLOAT(J),-DFLOAT(M),0.D0,DFLOAT(M))
     &        *THRJ(DFLOAT(I),1.D0,DFLOAT(J),0.D0,0.D0,0.D0)
     & !  the above is the contribution from DC field
     & !  the below is the contribution from AC field
     &      - (1.D0/6)*(-1)**(M) * DSQRT((2.d0*I+1.d0)*(2.d0*J+1.d0))
     &        *(alpha_parallel - alpha_perpendicular)*intensity
     &        *THRJ(DFLOAT(I),2.D0,DFLOAT(J),-DFLOAT(M),0.D0,DFLOAT(M))
     &        *THRJ(DFLOAT(I),2.D0,DFLOAT(J),0.D0,0.D0,0.D0)
     
	          InteractionHamiltonian(J,I)=InteractionHamiltonian(I,J) 
! since the matrix is symmetric
               END DO
            END DO
!            write(*,*) InteractionHamiltonian 
        DO L=ABS(M),NEndWith 
          DO K=ABS(M),NEndWith
            IF (L==K) THEN 
              MoleculeHamiltonian(L,K)=2*Pi*rot_const*(K+1)*K
     & ! below is the contribution from AC field        
     &           -(1.D0/4)*alpha_perpendicular*intensity
     &           -(1.D0/12)*(alpha_parallel - alpha_perpendicular)
     &             *intensity       
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
      END SUBROUTINE Matrix_for_new_EigStates_AC





!***************************************************************************** 
      SUBROUTINE RotSplit(electric_field,NewEigValue00,NewEigValue10,
     & NewEigValue20,NewEigValue11,NewEigValue21,NewEigValue22)
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
  
       
      END subroutine RotSplit
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
     &                 + dsqrt( (2.D0*NA1 + 1)*(2.D0*NA2 + 1)*(2.D0*NB1 + 1)*(2.D0*NB2 + 1) ) 
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
      END  FUNCTION DDInteract_in_Field
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
        
        !call Energy_in_kHz(interaction_energy) !convert atomic unit to kHz 
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
        !call Energy_in_kHz(interaction_energy) !convert atomic unit to kHz 
 

        do i=1, num_mole/2 + 1
          ka(i) = (i-1)*2*pi/num_mole
          summation(i)=2*half_sum1(num_mole, ka(i)) 
        end do
 

        exciton10_energy= interaction_energy*summation 
              
      END SUBROUTINE














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
           !call energy_in_kHz(prefactor)
           prefactor = -prefactor*1000*2*pi
           H_save = prefactor*H_save
        endif
        theta = theta_vs_time(T,flag)

        H = (1 - 3*dcos(theta)**2)*H_save 
        
      end subroutine form_matrix_1d 



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
          case default
             write(*,*) "Invalid time"
        end select
        return
      end function theta_vs_time   







!********************************one wavepacket**********************************************
! convert the gaussian wavepacket in k-space to wavepacket in real space
      subroutine initial_condition(Y_tmp,neq)
        implicit none
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer :: neq ! it is also the number of molecules in the crystal
        double complex :: Y_tmp(-neq/2:neq/2)
        double precision :: ka !wavevector
        double precision :: normalization_constant
        integer :: i,j
        double precision, external :: gaussian_func
        
        Y_tmp = (0.D0,0.D0)
        
        do i = -neq/2, neq/2 ! i ---- molecule index
        ! calculate C(n) for each molecule        
          do j = -neq/2, neq/2
            ! Fourier transform: C(n) = sum_{k} f(k)*exp(i*k*n), ignoring the prefactor in fourier transform
            ka = j*pi/(neq/2)
            Y_tmp(i) = Y_tmp(i) + gaussian_func(ka)*exp(ka*i*DCMPLX(0.D0,1.D0))
          end do  
          !write(17,*) i, Y_tmp(i)       
        end do
        
!     normalization of the coefficients
        normalization_constant = SUM(ABS(Y_tmp)**2)
        Y_tmp = Y_tmp/sqrt(normalization_constant)
      end subroutine initial_condition




! this subroutine gives the initial shape of the wavepacket in k-space
!*************************************************************************************
      function gaussian_func(ka)
        implicit none
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        double precision, parameter :: variance = 0.005D0 !0.005D0
        double precision, parameter :: central_ka = 0.D0 !pi/2 
        double precision :: ka
        double precision :: gaussian_func

        gaussian_func = (1./sqrt(2*pi*variance))*exp( -(ka- central_ka)**2/(2*variance) )
        return
      end function gaussian_func




!*************************************************************************************
      subroutine from_x_to_k(Y, neq, Y_k, C_k)   !initial_phase, phase_shift)
        implicit none
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
      	integer :: neq
        double complex :: Y(neq)
        double complex :: Y_k(-neq/2:neq/2) 
!        double precision :: initial_phase(-neq/2:neq/2)
!        double precision :: phase_shift(-neq/2:neq/2)
        double precision :: C_k(-neq/2:neq/2)
        integer :: ka_n
        integer :: n
        double precision :: ka
        double precision :: summation
        
!        do n=-neq/2, neq/2
!          phase_shift(n) = ATAN2( AIMAG(Y(n+neq/2+1)), REAL( Y(n+neq/2+1) )  ) - initial_phase(n)
!        end do
        
        Y_k = (0.D0, 0.D0)
        do ka_n = -neq/2, neq/2
          ka = ka_n*pi/(neq/2)
          do n = -neq/2, neq/2
            Y_k(ka_n) = Y_k(ka_n) + Y(n+ neq/2 + 1)*exp( -(0.D0,1.D0)*ka*n )
          end do
!          phase_shift(ka_n) = ATAN2( AIMAG(Y_k(ka_n)), REAL( Y_k(ka_n) )  ) - initial_phase(ka_n)         
        end do
        C_k = abs(Y_k)**2
        summation = SUM(C_k)
        C_k = C_k/summation
        Y_k = Y_k/sqrt(summation)
      end subroutine from_x_to_k


























!*****************find new rotational eigenstates (in terms of old ones)***********
! for |n=0,m=0> and |n=1,m=0> states
! accept the number of rotational levels for the calculation from outside
      SUBROUTINE New_EigStates_for_n01(En10, En00, coefficient_array_n0,
     &                      coefficient_array_n1,
     &                      Dipole, electric_field,
     &                      rot_const,Number_State) 

      	IMPLICIT NONE
      	DOUBLE PRECISION :: Dipole, electric_field, rot_const      	
      	INTEGER :: N, M, i
        DOUBLE PRECISION, ALLOCATABLE :: TotalHamiltonian(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: EigValue(:)
		double precision :: En10 !in kHz
		double precision :: En00
        DOUBLE PRECISION :: coefficient_array_n0(Number_State)
        DOUBLE PRECISION :: coefficient_array_n1(Number_State)        
      	INTEGER :: number_state
        external :: lapack_eig, Matrix_for_new_EigStates 

        N=0 
        M=0
!        write(11,*) electric_field
        ALLOCATE( TotalHamiltonian(Number_State,Number_State) ) 
       	ALLOCATE( EigValue(Number_State) )
!     .. Form the matrix ..		
       	call Matrix_for_new_EigStates(TotalHamiltonian,Number_State,N,M,
     &                  Dipole,electric_field,rot_const)

!      .... solve the eigenvalue and eigenvector ...	
        call lapack_eig(TotalHamiltonian,number_state, EigValue)

        En10 = EigValue(N-ABS(M)+1+1)
        En00 = EigValue(N-ABS(M)+1)
        !call energy_in_kHz(En10)  
        !call energy_in_kHz(En00)              
        ! calculate the coefficients of |n=0,m=0> in terms of the old rotational states        
      	coefficient_array_n0(1:Number_State) =
     &            TotalHamiltonian(1:Number_State, N-ABS(M)+1)
     
        ! calculate the coefficients of |n=1,m=0> in terms of the old rotational states     
        coefficient_array_n1(1:Number_State) =
     &            TotalHamiltonian(1:Number_State, N-ABS(M)+1+1)

        !for test purposes
!        do i=1, Number_State
!          write(13,*) i, coefficient_array_n0(i), coefficient_array_n1(i)  
!        end do
!       	write(13,*) "=============================="
!        write(13,*) " "
!        write(13,*) " "
!        write(13,*) " "
!        write(13,*) " "
	
       	DEALLOCATE(TotalHamiltonian,EigValue)
      END SUBROUTINE New_EigStates_for_n01
!**********************************************************************************************




!*****************find new rotational eigenstates (in terms of old ones)***********
! DC field with AC field
! for |n=0,m=0> and |n=1,m=0> states
! accept the number of rotational levels for the calculation from outside
!      SUBROUTINE New_EigStates_for_n01_AC(En10, En00, coefficient_array_n0,
!     &                      coefficient_array_n1,
!     &                      Dipole, DC_field, intensity,
!     &                      alpha_parallel, alpha_perpendicular, ! polarizability
!     &                      rot_const,Number_State) 

!      	IMPLICIT NONE
!      	DOUBLE PRECISION :: Dipole, DC_field, rot_const
!      	double precision :: intensity      	! laser intensity, related to the AC field
!      	double precision :: alpha_parallel, alpha_perpendicular
!      	INTEGER :: N, M, i
!        DOUBLE PRECISION, ALLOCATABLE :: TotalHamiltonian(:,:)
!        DOUBLE PRECISION, ALLOCATABLE :: EigValue(:)
!		double precision :: En10 !in kHz
!		double precision :: En00
!        DOUBLE PRECISION :: coefficient_array_n0(Number_State)
!        DOUBLE PRECISION :: coefficient_array_n1(Number_State)        
!      	INTEGER :: number_state
!        external :: lapack_eig 
!        external :: Matrix_for_new_EigStates_AC 

!        N=0 
!        M=0
!!        write(11,*) electric_field
!        ALLOCATE( TotalHamiltonian(Number_State,Number_State) ) 
!       	ALLOCATE( EigValue(Number_State) )
!!     .. Form the matrix ..		
!       	call Matrix_for_new_EigStates_AC(TotalHamiltonian,Number_State,N,M,
!     &                  Dipole, DC_field, intensity,
!     &                      alpha_parallel, alpha_perpendicular,rot_const)

!!      .... solve the eigenvalue and eigenvector ...	
!        call lapack_eig(TotalHamiltonian,number_state, EigValue)

!        En10 = EigValue(N-ABS(M)+1+1)
!        En00 = EigValue(N-ABS(M)+1)
!        !call energy_in_kHz(En10)  
!        !call energy_in_kHz(En00)              
!        ! calculate the coefficients of |n=0,m=0> in terms of the old rotational states        
!      	coefficient_array_n0(1:Number_State) =
!     &            TotalHamiltonian(1:Number_State, N-ABS(M)+1)
     
!        ! calculate the coefficients of |n=1,m=0> in terms of the old rotational states     
!        coefficient_array_n1(1:Number_State) =
!     &            TotalHamiltonian(1:Number_State, N-ABS(M)+1+1)

!       	DEALLOCATE(TotalHamiltonian,EigValue)
!      END SUBROUTINE New_EigStates_for_n01_AC
!**********************************************************************************************

      SUBROUTINE New_EigStates_for_n01_AC(En10, En00, coefficient_array_n0,
     &                      coefficient_array_n1,
     &                      Dipole, DC_field, intensity,
     &                      alpha_parallel, alpha_perpendicular, ! polarizability
     &                      rot_const,Number_State) 
        use molecule_information, only : electric_field
      	IMPLICIT NONE
      	double precision, parameter :: relative_error = 1.D-30
      	DOUBLE PRECISION :: Dipole, DC_field, rot_const  
      	double precision :: alpha_parallel, alpha_perpendicular
      	double precision :: intensity    	
      	INTEGER :: N, M, i
        DOUBLE PRECISION, ALLOCATABLE :: TotalHamiltonian(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: EigValue(:)
		double precision :: En10 !in au
		double precision, save :: En10_s ! s is for the fixed field
		double precision :: En00
		double precision, save :: En00_s
        DOUBLE PRECISION :: coefficient_array_n0(Number_State)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: coefficient_array_n0_s(:)
        DOUBLE PRECISION :: coefficient_array_n1(Number_State) 
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: coefficient_array_n1_s(:)       
      	INTEGER :: number_state
      	integer, save :: counter
        external :: lapack_eig, Matrix_for_new_EigStates 
        data counter / 0 /

        N=0 
        M=0
        if (counter == 0) then
          counter = 1        
          ALLOCATE( TotalHamiltonian(Number_State,Number_State) ) 
       	  ALLOCATE( EigValue(Number_State) )
       	  ALLOCATE( coefficient_array_n0_s(Number_State) )
       	  ALLOCATE( coefficient_array_n1_s(Number_State) )
!     .. Form the matrix ..		
       	  call Matrix_for_new_EigStates(TotalHamiltonian,Number_State,N,M,
     &                  Dipole,electric_field,rot_const)

!      .... solve the eigenvalue and eigenvector ...	
          call lapack_eig(TotalHamiltonian,number_state, EigValue)

          En10_s = EigValue(N-ABS(M)+1+1)
          En00_s = EigValue(N-ABS(M)+1)
           
        ! calculate the coefficients of |n=0,m=0> in terms of the old rotational states        
      	  coefficient_array_n0_s(1:Number_State) =
     &            TotalHamiltonian(1:Number_State, N-ABS(M)+1)
     
        ! calculate the coefficients of |n=1,m=0> in terms of the old rotational states     
          coefficient_array_n1_s(1:Number_State) =
     &            TotalHamiltonian(1:Number_State, N-ABS(M)+1+1)

          DEALLOCATE(TotalHamiltonian,EigValue)
        endif
        
        if ( ABS(Intensity)< ABS(relative_error) )  then      
          En10 = En10_s
          En00 = En00_s
          coefficient_array_n0 = coefficient_array_n0_s
          coefficient_array_n1 = coefficient_array_n1_s
        else
          ALLOCATE( TotalHamiltonian(Number_State,Number_State) ) 
       	  ALLOCATE( EigValue(Number_State) ) 
!     .. Form the matrix ..		
        	call Matrix_for_new_EigStates_AC(TotalHamiltonian,Number_State,N,M,
     &                  Dipole, DC_field, intensity,
     &                      alpha_parallel, alpha_perpendicular,rot_const)

!      .... solve the eigenvalue and eigenvector ...	
          call lapack_eig(TotalHamiltonian,number_state, EigValue)

          En10 = EigValue(N-ABS(M)+1+1)
          En00 = EigValue(N-ABS(M)+1)
                    
        ! calculate the coefficients of |n=0,m=0> in terms of the old rotational states        
      	  coefficient_array_n0(1:Number_State) =
     &            TotalHamiltonian(1:Number_State, N-ABS(M)+1)
     
        ! calculate the coefficients of |n=1,m=0> in terms of the old rotational states     
          coefficient_array_n1(1:Number_State) =
     &            TotalHamiltonian(1:Number_State, N-ABS(M)+1+1)
        	
       	  DEALLOCATE(TotalHamiltonian,EigValue)
       	endif
      END SUBROUTINE New_EigStates_for_n01_AC







!**********************************************************************************************
      subroutine coefficient_matrix(fixed_field, varying_field, n0_coeff, n1_coeff, 
     &                              Energy_array, num_mol, num_rot)
         use molecule_information       
         implicit none
         integer :: num_mol
         integer :: num_rot
         double precision :: fixed_field
         double precision :: varying_field(num_mol) 
         double precision :: En10, En00 
         double precision :: En10_array(num_mol)
         double precision :: En00_array(num_mol)
         double precision :: En00sum
         double precision :: Energy_array(num_mol)       
         double precision :: n0_coeff(num_rot, num_mol)
         double precision :: n1_coeff(num_rot, num_mol) 
         double precision :: coefficient_array_n0(num_rot)
         double precision :: coefficient_array_n1(num_rot) 
         integer :: n           
        

         
         do n=1, num_mol         
           call New_EigStates_for_n01(En10, En00, coefficient_array_n0,
     &                      coefficient_array_n1,
     &                      Dipole, fixed_field+varying_field(n),
     &                      rot_const,num_rot) 		
           n0_coeff(1:num_rot, n) =  coefficient_array_n0(1:num_rot) 
           n1_coeff(1:num_rot, n) =  coefficient_array_n1(1:num_rot) 
           En10_array(n) = En10
           En00_array(n) = En00 
!           write(8,*) n, Energy_array(n)
         end do
         
         En00sum = SUM(En00_array)
         do n=1, num_mol
           Energy_array(n) = En00sum - En00_array(n) + En10_array(n)
         end do
         
      end subroutine  coefficient_matrix
!*******************************************************************************************




!**********************************************************************************************
      subroutine coefficient_matrix_AC(fixed_field, intensity_array,n0_coeff, n1_coeff, 
     &                              Energy_array, num_mol, num_rot)
         use molecule_information       
         implicit none
         integer :: num_mol
         integer :: num_rot
         double precision :: fixed_field
         double precision :: intensity_array(num_mol)
         double precision :: En10, En00 
         double precision :: En10_array(num_mol)
         double precision :: En00_array(num_mol)
         double precision :: En00sum
         double precision :: Energy_array(num_mol)       
         double precision :: n0_coeff(num_rot, num_mol)
         double precision :: n1_coeff(num_rot, num_mol) 
         double precision :: coefficient_array_n0(num_rot)
         double precision :: coefficient_array_n1(num_rot) 
         integer :: n  
         double precision, external :: hartree_to_kHz         
        

         !write(*,*) "b1"
         do n=1, num_mol 
           !write(*,*) "b2"        
           call New_EigStates_for_n01_AC(En10, En00, coefficient_array_n0,
     &                      coefficient_array_n1,
     &                      Dipole, fixed_field, intensity_array(n),
     &                      alpha_parallel, alpha_perpendicular,
     &                      rot_const,num_rot) 		
           n0_coeff(1:num_rot, n) =  coefficient_array_n0(1:num_rot) 
           n1_coeff(1:num_rot, n) =  coefficient_array_n1(1:num_rot) 
           En10_array(n) = En10
           En00_array(n) = En00 
           write(8,*) "E10-E00=", hartree_to_kHz(En10-En00)/1.D6, "GHz"
         end do
         
         !write(*,*) "b3"
         En00sum = SUM(En00_array)
         do n=1, num_mol
           Energy_array(n) = En00sum - En00_array(n) + En10_array(n)
         end do
         
      end subroutine  coefficient_matrix_AC
!*******************************************************************************************
















!*******************************************************************************************
      function dipole_dipole_1001_AB(dipole, rot_const, lattice_constant, num_rot,
     &                               coeff_A_n0, coeff_A_n1, coeff_B_n0, coeff_B_n1)

	    IMPLICIT NONE
      	INTEGER :: I, J, K, L 
      	integer :: num_rot  
      	double precision, external :: dd
      	double precision, allocatable, save :: dd_array(:,:,:,:)	      	
      	DOUBLE PRECISION :: Dipole, rot_const , Lattice_Constant
      	DOUBLE precision :: dipole_dipole_1001_AB
        DOUBLE PRECISION :: coeff_A_n0(num_rot)
        DOUBLE PRECISION :: coeff_A_n1(num_rot)
        DOUBLE PRECISION :: coeff_B_n0(num_rot)
        DOUBLE PRECISION :: coeff_B_n1(num_rot)
        double precision :: prefactor
      	DOUBLE precision, EXTERNAL :: thrj
      	integer, save :: counter
      	data counter /0/
      	
      	if (counter==0) then
      	  allocate( dd_array(num_rot, num_rot, num_rot, num_rot) )
      	  counter = 1
          DO I=1, num_rot
            DO J=1, num_rot
              DO K=1, num_rot
                DO L=1, num_rot
                  dd_array(i,j,k,l) = dd(dipole, lattice_constant, num_rot,i,j,k,l)
                END DO
              END DO   
            END DO
          END DO         
		endif

        dipole_dipole_1001_AB = 0.D0

        DO I=1, num_rot
          DO J=1, num_rot
            DO K=1, num_rot
              DO L=1, num_rot
                dipole_dipole_1001_AB = dipole_dipole_1001_AB  
     &                 + coeff_A_n1(I)*coeff_A_n0(J)
     &   	           *coeff_B_n0(K)*coeff_B_n1(L)
     &                 *dd_array(i,j,k,l) 
              END DO
            END DO   
          END DO
        END DO

        RETURN
      END  FUNCTION dipole_dipole_1001_AB
!*******************************************************************************************


!*******************************************************************************************
      function dd(dipole, lattice_constant, num_rot,i,j,k,l)
	    IMPLICIT NONE
      	INTEGER :: NA1, MA1, NB1, MB1, NA0, MA0, NB0, MB0
      	integer :: i, j, k, l
      	integer :: num_rot  	      	
      	DOUBLE PRECISION :: Dipole, Lattice_Constant
      	DOUBLE precision :: dd
        double precision :: prefactor
      	DOUBLE precision, EXTERNAL :: thrj

		
! because the matrix element here is <10|<00| V |00>|10>
! <NA1,MA1|<NB0,MB0| V |NB1,MB1> |NA0,MA0> 
        NA1 = 1 ! A molecule in |n=1,m=0> state
        MA1 = 0
        NA0 = 0 ! A molecule in |n=0,m=0> state
        MA0 = 0

        NB1 = 1 ! B molecule in |n=1,m=0> state
        MB1 = 0
        NB0 = 0 ! B molecule in |n=0,m=0> state
        MB0 = 0        

		
        dd = 0.D0
        prefactor = dipole**2/(lattice_constant**3)

        dd = dsqrt( (2.D0*NA1 + 1)*(2.D0*NA0 + 1)*(2.D0*NB1 + 1)*(2.D0*NB0 + 1) ) 
!                         sqrt can only accept real variable (don't use integer)
     &                 *( thrj( 1.D0*(I-1) ,1.D0,1.D0*(J -1),0.D0,0.D0,0.D0)**2 ) 
     &                 *( thrj( 1.D0*(K-1) ,1.D0,1.D0*(L -1),0.D0,0.D0,0.D0)**2 ) 


        dd = prefactor*dd
        RETURN
      END  FUNCTION dd
!*******************************************************************************************






!*****************find the root of func(x) in range (x1,x2)*************************
      FUNCTION zriddr(func,x1,x2,xacc)

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
      END function zriddr


      


      
      
      
!*******************************************************************************************
      subroutine wavepacket_ham(neq,H,T)
        use some_parameters
        use molecule_information
        use blas95
        implicit none
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
      	double precision, save :: A_k  !electric field
      	integer, parameter :: num_rot = 10 ! the number of rotational levels included in the calculation
      	double precision :: phase_shift 
        integer :: neq ! the number of molecules in the crystal
        integer :: num_mol ! the number of molecules
        double precision :: H(neq,neq) 
        double precision :: H_tmp(neq, neq) 
        double precision :: basis_trans(neq,neq)  
        double precision :: basis_trans_inv(neq,neq)   
        double precision :: n0_coeff(num_rot,neq) ! each column represent the coefficients for 1-20 rotational states
        double precision :: n1_coeff(num_rot,neq) 
        double precision :: coeff_A_n0(num_rot) !for molecule A in |00> state
        double precision :: coeff_A_n1(num_rot) !for molecule A in |10> state
        double precision :: coeff_B_n0(num_rot) !for molecule B in |00> state 
        double precision :: coeff_B_n1(num_rot) !for molecule B in |10> state
        double precision, allocatable, save :: Energy_array_t0(:) !energies with respect the defined zero energy at t=0
        double precision :: Energy_array(neq)                             
        double precision :: T
        double precision, allocatable, save :: n0_coeff_t0(:,:) 
        double precision, allocatable, save :: n1_coeff_t0(:,:) 
        double precision :: prefactor
        double precision :: fixed_field ! set it equal to electric_field in molecule_information
        double precision :: varying_field(neq) !different molecules feel different electric fields at time t
        double precision :: varying_field_t0(neq)
        integer :: n, m ! molecule index
        integer :: i, j
        double precision, external :: dipole_dipole_1001_AB !dipole-dipole interaction between molecules A and B
        integer, save :: icounter
        data icounter /0/
        
        double precision, external :: mu_e_dressed_state

        
         !write(*,*) "==========="  
         
      
        num_mol = neq
        fixed_field = electric_field
        
        if (icounter == 0) then
          icounter = 1
          write(*,*) "==========="
          phase_shift = pi/2
          !time_scale = 1.D-6 !second
          !call find_Ak(A_k, t_var, phase_shift)
!          A_k = 6.626D-34/(2*5.529*3.33564D-30 *t_var)  !V/m
          ! scale A_k by a factor of 0.3731
          A_k = A_k*0.3731
          write(*,*) "A suitable A_k has been found"
          write(99,*) "A_k = ", (A_k/1.944690567144141D-12)/1.D5, "kV/cm"
          !call Field_in_atomic_unit(A_k)
          varying_field_t0 = 0.D0
          allocate( n0_coeff_t0(num_rot,neq) )
          allocate( n1_coeff_t0(num_rot,neq) )
          allocate(Energy_array_t0(neq))
          call coefficient_matrix_AC(fixed_field, 
     &                              n0_coeff_t0, n1_coeff_t0, 
     &                              Energy_array_t0, num_mol, num_rot)  
!          zero_energy = Energy_array_t0(1)       
        endif
        

        do i=-neq/2, neq/2
          varying_field(i+neq/2+1) = A_k*(-i)*( DSIN(pi*t/(t_var*4.134137D16) )**2 ) !*100
!          write(7,*) i, A_k, varying_field(i+neq/2+1)
        end do

        !write(23,*) t, varying_field(neq)
        !write(7,*) t, mu_e_dressed_state( varying_field(neq) )
        
 

        call coefficient_matrix(fixed_field, varying_field, n0_coeff, n1_coeff, 
     &                           Energy_array, num_mol, num_rot)

        ! diagonal part
!        H(1,1) = 0.D0 ! set as the zero energy 

        do n=1, neq
          H(n,n) = Energy_array(n) - Energy_array_t0(n)
          !write(21,*) n, H(n,n)
        end do
        
        ! off-diagonal part
        do n = 1, neq ! molecule A
          do m = n+1, neq ! molecule B
          coeff_A_n0(1:num_rot) = n0_coeff(1:num_rot, n)
          coeff_A_n1(1:num_rot) = n1_coeff(1:num_rot, n)
          coeff_B_n0(1:num_rot) = n0_coeff(1:num_rot, m)
          coeff_B_n1(1:num_rot) = n1_coeff(1:num_rot, m)

          H(n,m) = dipole_dipole_1001_AB(dipole, rot_const, lattice_constant, num_rot,
     &                               coeff_A_n0, coeff_A_n1, coeff_B_n0, coeff_B_n1) / abs(n-m)**3 

          !write(22,*) n, m, H(n,m)
          H(m,n) = H(n,m)                  
          end do
        end do
        

!        call form_transform_matrix(n0_coeff, n1_coeff, n0_coeff_t0, n1_coeff_t0, neq, num_rot, 
!     &                                 basis_trans, basis_trans_inv)
!        call gemm(H, basis_trans, H_tmp)
!        call gemm(basis_trans_inv, H_tmp, H)


!        do i=1, neq
!          do j=1, neq
!            write(14,*) i, j, H(i,j)
!          end do
!        end do
!        write(14,*) "====================================="
!        write(14,*) " "
        
      end subroutine wavepacket_ham 
!*******************************************************************************************






!*******************************************************************************************
      subroutine wavepacket_ham_AC1(neq,H,T)
        use some_parameters
        use molecule_information
        use blas95
        implicit none
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
      	double precision, parameter :: d0=0.000364653 ! meter, the distance between molecule 1 and beam center
      	integer, parameter :: light_region = 100
      	double precision, allocatable, save :: intensity(:)
      	integer :: displacement
      	double precision :: intensity_array(neq)
      	double precision :: I0(neq) ! initial intensity is zero
      	integer, parameter :: num_rot = 15 ! the number of rotational levels included in the calculation
      	double precision :: phase_shift 
        integer :: neq ! the number of molecules in the crystal
        integer :: num_mol ! the number of molecules
        double precision :: H(neq,neq) 
        double precision :: H_tmp(neq, neq) 
        double precision :: basis_trans(neq,neq)  
        double precision :: basis_trans_inv(neq,neq)   
        double precision :: n0_coeff(num_rot,neq) ! each column represent the coefficients for 1-20 rotational states
        double precision :: n1_coeff(num_rot,neq) 
        double precision :: coeff_A_n0(num_rot) !for molecule A in |00> state
        double precision :: coeff_A_n1(num_rot) !for molecule A in |10> state
        double precision :: coeff_B_n0(num_rot) !for molecule B in |00> state 
        double precision :: coeff_B_n1(num_rot) !for molecule B in |10> state
        double precision, allocatable, save :: Energy_array_t0(:) !energies with respect the defined zero energy at t=0
        double precision :: Energy_array(neq)                             
        double precision :: T
        double precision, allocatable, save :: n0_coeff_t0(:,:) 
        double precision, allocatable, save :: n1_coeff_t0(:,:) 
        double precision :: prefactor
        double precision :: fixed_field ! set it equal to electric_field in molecule_information
        integer :: n, m ! molecule index
        integer :: i, j
        double precision, external :: dipole_dipole_1001_AB !dipole-dipole interaction between molecules A and B
        integer, save :: icounter
        data icounter /0/
        
        double precision, external :: mu_e_dressed_state
        double precision, external :: second_to_au
        double precision, external :: intensity_func
        common /displace/ displacement
         !write(*,*) "==========="  
         
      
        num_mol = neq
        fixed_field = electric_field
        
        if (icounter == 0) then
          icounter = 1
          write(*,*) "==========="
          !phase_shift = pi/2
          !time_scale = 1.D-6 !second
          !call find_Ak(A_k, t_var, phase_shift)
!          A_k = 6.626D-34/(2*5.529*3.33564D-30 *t_var)  !V/m
          ! scale A_k by a factor of 0.3731
          !A_k = A_k*0.3731
          !write(*,*) "A suitable A_k has been found"
          !write(99,*) "A_k = ", (A_k/1.944690567144141D-12)/1.D5, "kV/cm"
          !call Field_in_atomic_unit(A_k)
          !varying_field_t0 = 0.D0
          !intensity = 15*pi/( (alpha_parallel-alpha_perpendicular)*second_to_au(t_var) )        !2*5.15D4 ! W/cm^2
          allocate(intensity(neq))
          do i=1,neq
            intensity(i) = intensity_func(d0 + i*400*1.D-9 )
            !write(18,*) i, intensity(i)
            call intensity_in_au( intensity(i) )
          end do
          !call intensity_in_au(intensity)
          I0 = 0.D0
          allocate( n0_coeff_t0(num_rot,neq) )
          allocate( n1_coeff_t0(num_rot,neq) )
          allocate(Energy_array_t0(neq))
          call coefficient_matrix_AC(fixed_field, I0,
     &                               n0_coeff_t0, n1_coeff_t0, 
     &                               Energy_array_t0, num_mol, num_rot) 
          !write(*,*) "a1" 
!          zero_energy = Energy_array_t0(1)       
        endif
        

        do i=1, 100 !neq !-neq/2, neq/2
          intensity_array(i) = Intensity(i)*( DSIN(pi*t/(t_var*4.134137D16) )**2 ) !*100
!          write(7,*) i, A_k, varying_field(i+neq/2+1)
        end do

        !write(23,*) t, varying_field(neq)
        !write(7,*) t, mu_e_dressed_state( varying_field(neq) )
        
 
          !write(*,*) "a2" 
        call coefficient_matrix_AC(fixed_field, intensity_array, 
     &                             n0_coeff, n1_coeff, 
     &                             Energy_array, num_mol, num_rot)
          !write(*,*) "a3" 
        ! diagonal part
        H(1,1) = 0.D0 ! set as the zero energy 

        do n=2, neq
          H(n,n) = Energy_array(n) - Energy_array(1)
          !write(21,*) n, H(n,n)
        end do
        
        ! off-diagonal part
        do n = 1, neq ! molecule A
          do m = n+1, neq ! molecule B
          coeff_A_n0(1:num_rot) = n0_coeff(1:num_rot, n)
          coeff_A_n1(1:num_rot) = n1_coeff(1:num_rot, n)
          coeff_B_n0(1:num_rot) = n0_coeff(1:num_rot, m)
          coeff_B_n1(1:num_rot) = n1_coeff(1:num_rot, m)

          H(n,m) = dipole_dipole_1001_AB(dipole, rot_const, lattice_constant, num_rot,
     &                               coeff_A_n0, coeff_A_n1, coeff_B_n0, coeff_B_n1) / abs(n-m)**3 

          !write(22,*) n, m, H(n,m)
          H(m,n) = H(n,m)                  
          end do
        end do
        
      end subroutine wavepacket_ham_AC1 
!*******************************************************************************************


!*******************************************************************************************
      subroutine wavepacket_ham_AC(neq,H,T)
        use some_parameters
        use molecule_information
        use blas95
        implicit none
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
      	double precision, parameter :: d0= 5.D-6 ! meter, the distance between molecule 1 and beam center
      	!integer, parameter :: light_region = 150
      	!integer, parameter :: imin = 425
      	!integer, parameter :: imax = imin + light_region-1
      	double precision, allocatable, save :: intensity(:)
      	integer :: displacement
      	double precision :: intensity_array(neq)
      	double precision :: I0(neq) ! initial intensity is zero
      	integer, parameter :: num_rot = 15 ! the number of rotational levels included in the calculation
        integer :: neq ! the number of molecules in the crystal
        integer :: num_mol ! the number of molecules
        double precision :: H(neq,neq)  
        double precision :: n0_coeff(num_rot,neq) ! each column represent the coefficients for 1-20 rotational states
        double precision :: n1_coeff(num_rot,neq) 
        double precision :: coeff_A_n0(num_rot) !for molecule A in |00> state
        double precision :: coeff_A_n1(num_rot) !for molecule A in |10> state
        double precision :: coeff_B_n0(num_rot) !for molecule B in |00> state 
        double precision :: coeff_B_n1(num_rot) !for molecule B in |10> state
        double precision, allocatable, save :: Energy_array_t0(:) !energies with respect the defined zero energy at t=0
        double precision :: Energy_array(neq)                             
        double precision :: T
        double precision, allocatable, save :: n0_coeff_t0(:,:) 
        double precision, allocatable, save :: n1_coeff_t0(:,:) 
        double precision :: fixed_field
        double precision, save :: dd_fixed
        integer :: n, m ! molecule index
        integer :: i, j
        double precision, external :: dipole_dipole_1001_AB !dipole-dipole interaction between molecules A and B
        integer, save :: icounter
        double precision, external :: intensity_func
        double precision, external :: hartree_to_kHz
        data icounter /0/

        
 
         
      
        num_mol = neq
        fixed_field = electric_field
        
        if (icounter == 0) then
          icounter = 1
          write(*,*) "==========="
!          allocate(intensity(light_region))
!          do i=1,light_region
!            intensity(i) = intensity_func(d0 + i*400*1.D-9 )
!            call intensity_in_au( intensity(i) )
!          end do
          allocate(intensity(neq))
          do i=1,neq
            intensity(i) = intensity_func(d0 + i*400*1.D-9 )
            call intensity_in_au( intensity(i) )
          end do
          I0 = 0.D0
          allocate( n0_coeff_t0(num_rot,neq) )
          allocate( n1_coeff_t0(num_rot,neq) )
          allocate(Energy_array_t0(neq))
          call coefficient_matrix_AC(fixed_field, I0,
     &                               n0_coeff_t0, n1_coeff_t0, 
     &                               Energy_array_t0, num_mol, num_rot) 
          n = 1
          m = 2
          coeff_A_n0(1:num_rot) = n0_coeff_t0(1:num_rot, n)
          coeff_A_n1(1:num_rot) = n1_coeff_t0(1:num_rot, n)
          coeff_B_n0(1:num_rot) = n0_coeff_t0(1:num_rot, m)
          coeff_B_n1(1:num_rot) = n1_coeff_t0(1:num_rot, m)
          dd_fixed =  dipole_dipole_1001_AB(dipole, rot_const, lattice_constant, num_rot,
     &                               coeff_A_n0, coeff_A_n1, coeff_B_n0, coeff_B_n1) !for two molecules in the fixed field 
     
          write(9,*) "J(a)=", hartree_to_kHz(dd_fixed), "kHz"  
        endif
        
        !intensity_array = 0.D0
        !intensity_array(imin:imax) = Intensity(1:light_region) *( DSIN(pi*t/(t_var*4.134137D16) )**2 )
        intensity_array = Intensity *( DSIN(pi*t/(t_var*4.134137D16) )**2 )
        
 

        call coefficient_matrix_AC(fixed_field, intensity_array, 
     &                             n0_coeff, n1_coeff, 
     &                             Energy_array, num_mol, num_rot)

        H = 0.D0
        
        do n=1, neq
          H(n,n) = Energy_array(n) - Energy_array_t0(n)
          !write(71,*) n, hartree_to_kHz( H(n,n) )
        end do
        
        ! off-diagonal part
        !do n = 1, neq ! molecule A
          !do m = n+1, neq ! molecule B
        
        ! only consider the nearest three neighbors
        do n = 1, neq ! molecule A
          do m = n+1, neq ! molecule B   
            if (abs(m-n)<=3) then 
            !-------------------------------------------------------------------------------------       
              !if ( ((n>=imin).and.(n<=imax)) .and. ((m>=imin).and.(m<=imax))  ) then
             
                coeff_A_n0(1:num_rot) = n0_coeff(1:num_rot, n)
                coeff_A_n1(1:num_rot) = n1_coeff(1:num_rot, n)
                coeff_B_n0(1:num_rot) = n0_coeff(1:num_rot, m)
                coeff_B_n1(1:num_rot) = n1_coeff(1:num_rot, m)

                H(n,m) = dipole_dipole_1001_AB(dipole, rot_const, lattice_constant, num_rot,
     &                               coeff_A_n0, coeff_A_n1, coeff_B_n0, coeff_B_n1) / abs(n-m)**3 

                H(m,n) = H(n,m)
             
!              else if ( ((n>=imin).and.(n<=imax)) .and. ((m<imin).or.(m>imax))  ) then
!               else if ( ((n>=imin).and.(n<=imax)) .and. (m>imax)  ) then
!                coeff_A_n0(1:num_rot) = n0_coeff(1:num_rot, n)
!                coeff_A_n1(1:num_rot) = n1_coeff(1:num_rot, n)
!                coeff_B_n0(1:num_rot) = n0_coeff_t0(1:num_rot, m)
!                coeff_B_n1(1:num_rot) = n1_coeff_t0(1:num_rot, m)

!                H(n,m) = dipole_dipole_1001_AB(dipole, rot_const, lattice_constant, num_rot,
!     &                               coeff_A_n0, coeff_A_n1, coeff_B_n0, coeff_B_n1) / abs(n-m)**3 

!                H(m,n) = H(n,m)
            
!              else if ( ((n>imax).and.(n<imax+50)).and.((m>imax).and.(m<imax+50)) ) then
!                H(n,m) = dd_fixed/ABS(n-m)**3
!                H(m,n) = H(n,m)
              else
                H(n,m) = 0.D0
              endif
            !-------------------------------------------------------------------------------------  
            !endif                  
          end do
        end do
                
      end subroutine wavepacket_ham_AC 
!*******************************************************************************************









      subroutine matrix_inverse(A, A_inverse, n)
        use lapack95
        implicit none
        integer :: n
        double precision :: A(n,n)
        double precision :: A_inverse(n,n)
        integer :: ipiv(n)
        
        A_inverse = A
        call getrf( A_inverse, ipiv )
        call getri( A_inverse, ipiv ) ! A_inverse overwritten by the inverse
        
      end subroutine matrix_inverse





      Function intensity_func(z)
        implicit none
        DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        double precision, parameter :: lambda = 1000*1.D-9 ! wavelength is 800 nm
        double precision, parameter :: w0 = 5*1.D-6 ! beam waist is 5 micron meter
        double precision :: I0 = 1.D7 ! W/cm^2, intensity at the center
        double precision :: intensity_func
        double precision :: z
        double precision :: wz
        double precision :: zr
        
        zr = Pi*w0**2/lambda
        wz = w0*sqrt(1 + (z/zr)**2)
        intensity_func = I0*(w0/wz)**2
        return
      end function intensity_func
