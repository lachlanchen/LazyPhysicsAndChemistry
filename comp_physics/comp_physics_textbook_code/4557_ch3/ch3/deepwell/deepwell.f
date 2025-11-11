      PROGRAM DeepWell

C Program to calculate the spectrum of the infinitely deep
C potential well. 
C Program described in "Computational Physics", J. M. Thijssen
C Section 3.2.1
C Program written by J. M. Thijssen, June 1998
C
C PROGRAM USES ROUTINE GenEig IN FILE 'geneig.f'
C

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
C There is no distinction between odd and even basis functions. 
C It is more efficient to make this distinction as it leads to 
C two N/2xN/2 matrices rather than a single NxN matrix. 
C Remember the (generalised) eigenvalue problem is an order N^3 
C problem!

      INCLUDE 'globvarc'

      INTEGER I, J, K, IPLusJ, PtNum
      PARAMETER (PtNum=60)

      DOUBLE PRECISION SMatrix(0:MaxN, 0:MaxN), HMatrix(0:MaxN, 0:MaxN),
     .                 Diag(0:MaxN), X, Res, PI

C Fill upper triangles of HMatrix and SMatrix
      DO I=0,N-1
        DO J=0,N-1
          IPlusJ = I+J
          IF (MOD(IPlusJ,2) .EQ. 0) THEN
            SMatrix(I,J) = 2.D0/(IPlusJ+5.D0) - 4.D0/(IPlusJ+3.D0) +
     .                     2.D0/(IPlusJ+1.D0)

            HMatrix(I,J) = 8*(I+J+2.D0*I*J-1)/(IPlusJ+3.D0)/
     .                        (IPlusJ+1.D0)/(IPlusJ-1.D0)
          ELSE
            SMatrix(I,J) = 0.D0
            HMatrix(I,J) = 0.D0
          ENDIF
        ENDDO
      ENDDO

C The next call is to a LAPACK routine which solves the generalised eigenvalue 
C problem

      CALL GenEig(HMatrix, SMatrix, N, MaxN+1, Diag)

C Output the variational eigenvalues to the screen  together with the exact ones
      PI = 4.D0*DATAN(1.D0)
      WRITE (6,*) 'Variational     Exact'
      DO I=0, N-1
        WRITE (6,'(5F12.4)') Diag(I), PI*PI/4*(I+1)*(I+1)
      ENDDO

C Write variational eigenfunctions
      OPEN (10, File='EigVecs')
      DO I=0, N-1
        DO J=-PtNum, PtNum
          X = 1.D0*J/PtNum
          Res = 0.D0
          DO K=0, N-1
            Res = Res+X**K*(X**2-1)*HMatrix(K, I)
          END DO
          WRITE (10, '(2F15.8)') X, Res
        END DO
      END DO
      CLOSE(10)

      END


