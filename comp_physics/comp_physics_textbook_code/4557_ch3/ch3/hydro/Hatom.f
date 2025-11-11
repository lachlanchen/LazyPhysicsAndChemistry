      PROGRAM Hatom

C Program to calculate the l=0 spectrum of the hydrogen
C atom (l is angular momentum quantum number). 
C Program described in "Computational Physics", J. M. Thijssen
C Section 3.2.2
C Program written by J. M. Thijssen, June 1998
C
C PROGRAM USES ROUTINE GenEig IN FILE 'geneig.f'
C
C      INCLUDE 'globHatom'
C The file 'globHatom' contains global variable declarations 
C and the IMPLICIT NONE statement. See this file for details

      CALL Initialise 
      CALL Overlap
      CALL CalcHMatrix
      CALL CalcSpectrum
      END



      SUBROUTINE Initialise
C Initialise Gaussian exponents Alpha and 
C scale them according to the nuclear charge

      INCLUDE 'globHatom'
      INTEGER I

      Alpha(1) = 13.00773D0
      Alpha(2) = 1.962079D0
      Alpha(3) = 0.444529D0
      Alpha(4) = 0.1219492D0
      PI = 4.D0*atan(1.D0)
      Z = -1.D0
      DO WHILE (Z.LE.0.D0) 
        WRITE (6, *) 'Give nuclear charge Z'
        READ (5, *) Z
        IF (Z.LE.0.D0) THEN
          WRITE (6,*) 'Give a positive number'
        END IF
      END DO

      DO I=1, 4
        Alpha(I) = Alpha(I)*(Z*Z)
      END DO

      END





      SUBROUTINE Overlap
C Calculate overlap matrix

      INCLUDE 'globHatom'

      INTEGER R, S

      DOUBLE PRECISION Factor

      DO R=1,N
        DO S=1,R-1
          Factor = Pi/(Alpha(R)+Alpha(S))
          SMatrix(R,S) = Factor*DSQRT(Factor)
          SMatrix(S,R) = SMatrix(R,S)
        ENDDO
        Factor    = Pi/(Alpha(R)+Alpha(R))
        SMatrix(R,R)    =  Factor*DSQRT(Factor)
      ENDDO

      CALL WriteMatrix (SMatrix, "overlap", NSize, NSize)
      END



      SUBROUTINE CalcHMatrix
C Calculate Hamiltonian 'HMatrix'
C This routine calls the routines Kinet11 and Coul11 which 
C calculate the kinetic and Coulomb potential contribution 
C respectively.

      INCLUDE 'globHatom'

      INTEGER R, S

      DOUBLE PRECISION A, B, Kinet11, Coul11

      DO R=1,N
        DO S=1,R-1
          A = Alpha(R)
          B = Alpha(S)
          HMatrix(R,S) = Kinet11(A,B)+Coul11(A,B)
          HMatrix(S,R) = HMatrix(R,S)
        ENDDO
        A = Alpha(R)
        HMatrix(R,R) = Kinet11(A,A)+Coul11(A,A)
      ENDDO

      CALL WriteMatrix (HMatrix, "hamilton", NSize, NSize)
      END





      DOUBLE PRECISION FUNCTION Kinet11 (A, B)
C Kinetic matrix element. The subscript 11 indicates that the
C basis functions are centred on the same nucleus.

      INCLUDE 'globHatom'      
      
      DOUBLE PRECISION A, B, Alph, Factor

      Alph     = A+B
      Factor   = Pi/Alph
      Kinet11   = 3.D0*Factor*DSQRT(Factor)*A*B/Alph
      END





      DOUBLE PRECISION FUNCTION Coul11 (A, B)
C Coulomb matrix element. The subscript 11 indicates that the
C basis functions are centred on the same nucleus.

      INCLUDE 'globHatom'

      DOUBLE PRECISION A, B

      Coul11 = -2.D0*Z*Pi/(A+B)

      END



      SUBROUTINE WriteMatrix (Matrix, FileName, N, NSize)
C This routine is used to write matrices in a file for checking
C purposes. Calls to this routine may be removed from the program.

      IMPLICIT NONE

      INTEGER N, NSize, I, J

      DOUBLE PRECISION Matrix(NSize, NSize)

      CHARACTER*7 FileName

      OPEN (11, File=FileName)
      REWIND(11)

      WRITE (11, '(8F12.6)') ((Matrix(I, J), I=1, N), J=1, N)

      CLOSE(11)
      END



      SUBROUTINE CalcSpectrum
C The actual calculation. The call to GenEig solves the
C generalised eigenvalue problem. The eigenvalue is then output
C to screen and the variational groundstate wave function is
C written to the file 'WaveFunc'

      INCLUDE 'globHatom'

      DOUBLE PRECISION Diag(NSize), R, R2, Val
      INTEGER I, J

      CALL GenEig (HMatrix, SMatrix, N, NSize, Diag)
      WRITE (6, '("Lowest eigenvalue is", F12.8)') Diag(1)
      WRITE (6, '("Exact value         ", F12.8)') -Z*Z*0.5D0
      
C Write ground state wave function to file 'WaveFunc':   
      OPEN (10, File='WaveFunc')
      DO I=0, 100
        R = I/100.D0*4.D0/Z
        Val = 0.D0
        R2 = R*R
        DO J=1, N
          Val = Val + HMatrix(J, 1)*EXP(-Alpha(J)*R2)
        END DO
        WRITE (10, '(2F10.6)') R, Val
      END DO
      CLOSE (10)
      END
