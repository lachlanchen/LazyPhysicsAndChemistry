      PROGRAM OptHeHF

C Program to calculate the l=0 ground state of the helium
C atom (l is angular momentum quantum number). 
C Program described in "Computational Physics", J. M. Thijssen
C Care has been taken to optimise the program using the symmetry
C of the various matrix elements.
C Section 4.3.2
C Program written by J. M. Thijssen, June 1998

      include "globHe"
C The file 'globHe' contains global variable declarations 
C and the IMPLICIT NONE statement. See this file for details

      DOUBLE PRECISION OldEner, Energy, Precision
      PARAMETER (Precision = 1.D-11)

C Initialise: calculate all matrices which do not change
C during the self-consistency iterations.....
      CALL Initialise 
      CALL CalcOverlap
      CALL CalcVMatrix(SMatrix, VMatrix, N, NSize)
      CALL CalcHMatrix
      CALL BuildSuper

      OldEner = -1.D0
      Energy = 0.D0
          
C Self-consistency loop.....
      DO WHILE(ABS(OldEner-Energy).GT. Precision)
        OldEner = Energy
        CALL BuildG
        CALL CalcFockMatrix
        CALL DiagFock(Energy)
      ENDDO

      END


      SUBROUTINE DiagFock(Energy)

      include "globHe"

      DOUBLE PRECISION Energy, 
     .       Diag(NSize), Eigen(NSize, NSize)

      INTEGER I, J

      CALL LowdinDiag (VMatrix, FMatrix, Eigen, Diag, N, NSize)
      DO I=1, N
        C(I) = Eigen(I, 1)
      END DO
      CALL BuildDensMat
      Energy = Diag(1)
      DO I=1,N
        DO J=1, I-1
          Energy = Energy + DensMat(I,J)*Hamilton(I,J)
        ENDDO
        Energy = Energy + 0.5D0*DensMat(I,I)*Hamilton(I,I)
      ENDDO
      WRITE (6,*) 'The total energy is ....... ', Energy
      END         





      SUBROUTINE Initialise

      INCLUDE "globHe"

      INTEGER I

      DATA Alpha/0.298073,1.242567,5.782948,38.47497/

      PI = 4.D0*atan(1.D0)

C Set the initial eigenvector C equal to
C (1, 1, 1, 1). Other choices are possible.
      DO I=1, N
        C(I) = 1.D0
      ENDDO
      
      END



      SUBROUTINE CalcOverlap
C Calculate the overlap matrix

      INCLUDE "globHe"

      INTEGER R, S

      DOUBLE PRECISION Factor

      DO R=1,N
        DO S=1,R-1
          Factor = Pi/(Alpha(R)+Alpha(S))
          SMatrix(R,S) = Factor*DSQRT(Factor)
        ENDDO
        Factor    = Pi/(Alpha(R)+Alpha(R))
        SMatrix(R,R)    =  Factor*DSQRT(Factor)
      ENDDO
      END



      SUBROUTINE CalcHMatrix
C Calculate the Hamiltonian matrix

      INCLUDE "globHe"

      INTEGER R, S

      DOUBLE PRECISION A, B, Kinet11, Coul11

      DO R=1,N
        DO S=1,R-1
          A = Alpha(R)
          B = Alpha(S)
          Hamilton(R,S) = Kinet11(A,B)+Coul11(A,B)
c          Hamilton(S,R) = Hamilton(R,S)
        ENDDO
        A = Alpha(R)
        Hamilton(R,R) = Kinet11(A,A)+Coul11(A,A)
      ENDDO
      END





      DOUBLE PRECISION FUNCTION Kinet11 (A, B)
C Kinetic matrix element. The subscript 11 indicates that the
C basis functions are centred on the same nucleus.

      INCLUDE "globHe"      
      
      DOUBLE PRECISION A, B, Alph, Factor

      Alph     = A+B
      Factor   = Pi/Alph
      Kinet11   = 3.D0*Factor*DSQRT(Factor)*A*B/Alph
      END



      DOUBLE PRECISION FUNCTION Coul11 (A, B)
C Coulomb matrix element. The subscript 11 indicates that the
C basis functions are centred on the same nucleus.
C The charge Z=2 of the helium nucleus is taken into account

      INCLUDE "globHe"

      DOUBLE PRECISION A, B

      Coul11 =  -4.D0*Pi/(A+B)

      END




      SUBROUTINE WriteMatrix (Matrix, FileName, N, NSize)
C This routine is used to write matrices in a file for checking
C purposes. Calls to this routine may be removed from the program.

      IMPLICIT NONE

      INTEGER N, NSize, I, J

      DOUBLE PRECISION Matrix(NSize, NSize)

      CHARACTER*16 FileName

      OPEN (11, File=FileName)
      REWIND(11)

      WRITE (11, '(4F12.6)') ((Matrix(I, J), I=1, N), J=1, N)

      CLOSE(11)
      END



      SUBROUTINE BuildDensMat
C Calculates the normalies density matrix D=2 C C^T/(C^T C),
C where C is the eigenvector corresponding to the ground state.

      INCLUDE "globHe"

      INTEGER R, S
      DOUBLE PRECISION Norm

      Norm = 0.D0
      DO R=1,N
        DO S=1,R-1
          DensMat(R, S) = 2.D0*C(R)*C(S)
          Norm = Norm + DensMat(R,S)*SMatrix(R,S)
        END DO
        DensMat(R,R) = 2.D0*C(R)*C(R)
        Norm = Norm + 0.5D0*DensMat(R,R)*SMatrix(R,R)
      END DO
      Norm = 1.D0/Norm
      DO R=1, N
        DO S=1, R
          DensMat(R,S) = DensMat(R,S)*Norm
        END DO
      END DO
      END





      SUBROUTINE BuildG

C Calculation of the matrix 
C G(T, U) = 0.5 Sum_{R,S} D(R, S) <RT| g |SU>
C The symmetries of the two-electron matrix elements has been exploited
C as much as possible in order to speed up this step, resulting however
C in some book-keeping. 
C This is done only to illustrate the method: for the small number 
C of basis functions used here, this procedure does not result in 
C a noticeable speed-up, and might even result in a somewhat 
C slower program. 
C of the matrix elements: only the symmetry r<->s has 
C been used.

      include "globHe"

      REAL*8 SuperElem


      INTEGER R, S, T, U, MaxU

      DO R=1,N
         DO S=1,R
            GMatrix(R,S) = 0.D0
         ENDDO
      ENDDO

      DO R=1,N
        DO S=1,R
          DO T=1,R
            IF (T .LT. R) THEN
              MaxU = T
            ELSE
              MaxU = S
            ENDIF
            DO U=1,MaxU
              SuperElem = QMatrix(R,S,T,U)
              IF (TypeMatrix(R,S,T,U).EQ.1) THEN
                GMatrix(R,R) = GMatrix(R,R)+0.5*SuperElem*
     .                                          DensMat(R,R)
              ELSE IF (TypeMatrix(R,S,T,U).EQ.2) THEN
                GMatrix(R,R) = GMatrix(R,R)+0.5D0*SuperElem*
     .                                          DensMat(T,T)
                GMatrix(T,T) = GMatrix(T,T)+0.5D0*SuperElem*
     .                                            DensMat(R,R)
              ELSE IF (TypeMatrix(R,S,T,U).EQ.3) THEN
                GMatrix(R,R) = GMatrix(R,R)+SuperElem*
     .                                      DensMat(T,U)
                GMatrix(T,U) = GMatrix(T,U)+0.5D0*SuperElem*
     .                                      DensMat(R,R)
              ELSE IF (TypeMatrix(R,S,T,U).EQ.4) THEN
                GMatrix(T,T) = GMatrix(T,T)+SuperElem*
     .                                      DensMat(R,S)
                GMatrix(R,S) = GMatrix(R,S)+0.5D0*SuperElem*
     .                                            DensMat(T,T)
              ELSE IF (TypeMatrix(R,S,T,U).EQ.5) THEN
                GMatrix(R,S) = GMatrix(R,S)+SuperElem*
     .                                         DensMat(R,S)
              ELSE IF (TypeMatrix(R,S,T,U).EQ.6) THEN
                GMatrix(R,S) = GMatrix(R,S)+SuperElem*
     .                                          DensMat(T,U)
                GMatrix(T,U) = GMatrix(T,U)+SuperElem*
     .                                          DensMat(R,S)
              END IF
            ENDDO
          ENDDO
        ENDDO
      ENDDO   
      END


      SUBROUTINE CalcFockMatrix
C Calculate Fock matrix: F = H + G

      INCLUDE "globHe"

      INTEGER  R, S

      DO R = 1, N
        DO S = 1, R
          FMatrix(R,S) = Hamilton(R,S) + GMatrix(R,S)
        ENDDO
      ENDDO

      END





      SUBROUTINE BuildSuper
C Calculation of the two-electron matrix elements 
C <rt|g|su> Only elements with r>=s are calculated.

      INCLUDE "globHe"

      DOUBLE PRECISION HF1111

      INTEGER R, S, T, U, MaxU, DetermineType

      DO R=1,N
        DO S=1,R
          DO T=1,R
            IF (T .LT. R) THEN
              MaxU = T
            ELSE
              MaxU = S
            ENDIF
            DO U=1,MaxU
              TypeMatrix (R,S,T,U) = DetermineType (R,S,T,U)
              QMatrix(R,S,T,U) = HF1111(R,S,T,U)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      END


      INTEGER FUNCTION DetermineType(R,S,T,U)
C This function is called when the two-electron matrix elements are calculated.
C The aim is to determine what the relation between the indices r,s,t,u
C is, so that GMatrix can be filled in an efficient way (see BuildG)

      IMPLICIT NONE 

      INTEGER R, S, T, U, Type

      TYPE = 0
      IF ((R .EQ. S) .AND. (T .EQ. U) .AND. (R.EQ.T)) THEN
         Type = 1
      ELSEIF ((R .EQ. S) .AND. (T .EQ. U)) THEN
         Type = 2
      ELSEIF ((R .EQ. S) .AND. (T .GT. U)) THEN
         Type = 3
      ELSEIF ((R .GT. S) .AND. (T .EQ. U)) THEN
         Type = 4
      ELSEIF ((R .GT. S) .AND. (T .GT. U) .AND. 
     .        (R.EQ.T) .AND. (S.EQ.U))THEN
         Type = 5
      ELSE
         Type = 6
      ENDIF

      DetermineType = Type

      END




      DOUBLE PRECISION FUNCTION HF1111 (R, S, T, U)
C Four-electron matrix element. The subscript 1111 indicates
C that all orbitals are centred on the same nucleus.

      INCLUDE "globHe"

      DOUBLE PRECISION A, B, PiFac

      INTEGER R, S, T, U

      PiFac = 2*PI*PI*DSQRT(PI)

      A = Alpha(R) + Alpha(S)
      B = Alpha(T) + Alpha(U)

      HF1111 = PiFac/DSQRT(A+B)/A/B

      END


