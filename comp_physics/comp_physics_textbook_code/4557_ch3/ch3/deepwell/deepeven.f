      PROGRAM DeepEven

C Program to calculate the spectrum of the even states in the
C infinitely deep potential well. 
C Program described in "Computational Physics", J. M. Thijssen
C Section 3.2.1
C Program written by J. M. Thijssen, June 1998

      INCLUDE 'globvarc'

      N = MaxN*2
      DO WHILE (N.GT.MaxN)
        WRITE (6,*) 'Give nr. of basis states'
        READ (5, *) N
        IF (N.GT.MaxN) WRITE(6,*) 'N must be smaller than ', MaxN
      END DO

      CALL CalcSpectrum
      END 


      SUBROUTINE CalcSpectrum
C The hamilton and overlap matrices are filled, and the generalised
C eigenvalue problem is solved. The eigenvalues are output to
C the screen and the eigenvectors to a file. 
C In this program, only the even basis functions are treated. 

      INCLUDE 'globvarc'

      INTEGER I, J, K, IPLusJ, PtNum
      PARAMETER (PtNum=60)

      DOUBLE PRECISION SMatrix(0:MaxN, 0:MaxN), HMatrix(0:MaxN, 0:MaxN),
     .                 Diag(0:MaxN), X, Res, PI

C Fill upper triangles of HMatrix and SMatrix
      DO I=0,N-1
        DO J=I,N-1
          IPlusJ = 2*(I+J)
          SMatrix(I,J) = 2.D0/(IPlusJ+5.D0) - 4.D0/(IPlusJ+3.D0) +
     .                   2.D0/(IPlusJ+1.D0)

          HMatrix(I,J) = 8*(-1.D0+2*I+2*J+8.D0*I*J)/(IPlusJ+3.D0)/
     .                            (IPlusJ+1.D0)/(IPlusJ-1.D0)
        ENDDO
      ENDDO

C The next call is to a LAPACK routine which solves the generalised eigenvalue 
C problem

      CALL GenEig(HMatrix, SMatrix, N, MaxN+1, Diag)

C Output the variational eigenvalues to the screen  together with the exact ones
      PI = 4.D0*DATAN(1.D0)
      WRITE (6,*) 'Variational     Exact'
      DO I=0, N-1
        WRITE (6,'(5F12.4)') Diag(I), PI*PI*(2*I+1)*(2*I+1)*0.25D0
      ENDDO

C Write variational eigenfunctions
      OPEN (10, File='EvenVecs')
      DO I=0, N-1
        DO J=-PtNum, PtNum
          X = 1.D0*J/PtNum
          Res = 0.D0
          DO K=0, N-1
            Res = Res+X**(2*K)*(X**2-1)*HMatrix(K, I)
          END DO
          WRITE (10, '(2F15.8)') X, Res
        END DO
        WRITE (10,*)
      END DO
      CLOSE(10)

      END




