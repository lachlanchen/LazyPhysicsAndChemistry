      PROGRAM H2HF

C Program to calculate the l=0 ground state of the h2 molecule
C atom (l is angular momentum quantum number). 
C Program described in "Computational Physics", J. M. Thijssen
C Section 4.3.2
C Program written by J. M. Thijssen, June 1998-2016
C The program calculates the total energy for 100 distances 
C between 1 and 2 a.u. The results are stored in the file `e_vs_dist.dat'
C Thanks to Dario Mitnik for pointing out an error in the previous version

      include "globH2"
C The file 'globH2' contains global variable declarations 
C and the IMPLICIT NONE statement. See this file for details

      DOUBLE PRECISION OldEner, Energy, Precision
      INTEGER I
      PARAMETER (Precision = 1.D-11)

C Initialise: calculate all matrices which do not change
C during the self-consistency iterations.....
      CALL Initialise 
      DO I=1, 100
        Dist = 1.0D0 + 1.D-2*(I-1)
        CALL CalcOverlap
        CALL CalcVMatrix(SMatrix, VMatrix, 2*N, NSize)
        CALL CalcHMatrix
        CALL BuildSuper
        OldEner = -1.D0
        Energy = 0.D0
C Self-consistency loop.....
        DO WHILE(ABS(OldEner-Energy).GT. Precision)
          CALL BuildDensMat
          OldEner = Energy
          CALL BuildG
          CALL CalcFockMatrix
          CALL DiagFock(Energy)
        END DO
        WRITE (6,*) 'The total energy is ....... ', Energy + 1/dist
        WRITE (8,*) Dist, Energy+1.D0/Dist
      ENDDO
      CLOSE(8)
      END


      SUBROUTINE DiagFock(Energy)

      include "globH2"

      DOUBLE PRECISION Energy, 
     .       Diag(NSize), Eigen(NSize, NSize)

      INTEGER I, J

      CALL LowdinDiag (VMatrix, FMatrix, Eigen, Diag, 2*N, NSize)
      DO I=1, 2*N
        C(I) = Eigen(I, 1)
      END DO
      CALL BuildDensMat
      Energy = 0.D0  
      DO I=1,2*N
         DO J=1, I-1
            Energy = Energy + DensMat(I,J)*(HMatrix(I,J)+FMatrix(I,J))
         ENDDO
         Energy = Energy+0.5D0*DensMat(I,I)*(HMatrix(I,I)+FMatrix(I,I))
      ENDDO
      END         





      SUBROUTINE Initialise

      INCLUDE "globH2"

      INTEGER I

      DATA Alpha/13.00773D0,1.962079D0,0.444529D0,0.1219492D0,
     .           13.00773D0,1.962079D0,0.444529D0,0.1219492D0/


      PI = 4.D0*atan(1.D0)

C Set the initial eigenvector C equal to
C (1, 1, 1, 1). Other choices are possible.
      DO I=1, 2*N
        C(I) = 1.D0
      ENDDO
      OPEN (8, file='e_vs_dist.dat')
     
      END



      SUBROUTINE CalcOverlap
C Calculate the overlap matrix, S(r,s); r>s

      INCLUDE "globH2"

      INTEGER R, S

      DOUBLE PRECISION Factor, A, B

      DO R=1,N
        DO S=1,R-1
          A = Alpha(R)
          B = Alpha(S)
          Factor = Pi/(A+B)
          SMatrix(R,S) = Factor*DSQRT(Factor)
          SMatrix(R+N,S+N) = SMatrix(R,S)
          SMatrix(R+N,S) = SMatrix(R,S)*EXP(-A*B/(A+B)*Dist*Dist)
          SMatrix(S+N,R) = SMatrix(R+N,S)
        ENDDO
        A = Alpha(R)
        Factor    = 0.5D0*Pi/A
        SMatrix(R,R)    =  Factor*DSQRT(Factor)
        SMatrix(R+N,R+N) = SMatrix(R,R)
        SMatrix(R+N,R) = SMatrix(R,R)*EXP(-0.5D0*A*Dist*Dist)
      ENDDO
      END



      SUBROUTINE CalcHMatrix
C Calculate the HMatrixian matrix H(r,s); r>s

      INCLUDE "globH2"

      INTEGER R, S

      DOUBLE PRECISION A, B, Kinet11, Coul11, Kinet12, Coul12

      DO R=1,N
        DO S=1,R-1
          A = Alpha(R)
          B = Alpha(S)
          CALL Kinet(A, B, SMatrix(R+N, S), Kinet11, Kinet12)
          CALL Coul(A, B, SMatrix(R+N, S), Coul11, Coul12)
          HMatrix(R,S) = Kinet11+Coul11
          HMatrix(R+N,S+N) = HMatrix(R,S)
          HMatrix(R+N,S) = Kinet12+Coul12
          HMatrix(S+N,R) = HMatrix(R+N,S)
        ENDDO
        A = Alpha(R)
        CALL Kinet(A, A, SMatrix(R+N, R), Kinet11, Kinet12)
        CALL Coul(A, A, SMatrix(R+N, R), Coul11, Coul12)
        HMatrix(R,R) = Kinet11+Coul11
        HMatrix(R+N,R+N) = HMatrix(R,R)
        HMatrix(R+N,R) = Kinet12+Coul12
      ENDDO
      END





      SUBROUTINE Kinet (A, B, SElem, Kinet11, Kinet12)
C Kinetic matrix element. The subscript 11 indicates that the
C basis functions are centred on the same nucleus, 12 means that
C the kinetic matrix element is calculated for two basis functions
C centred at different nuclei.

      INCLUDE "globH2"      
      
      DOUBLE PRECISION A, B, K0, Kinet11, Kinet12,
     .                 SElem

      K0 = A*B/(A+B)
      Kinet11 = 3*K0*(PI/(A+B))**1.5D0
      Kinet12 = K0*(3-2*K0*Dist*Dist)*SElem
      END


      DOUBLE PRECISION FUNCTION F (X)

      include "globH2"

      DOUBLE PRECISION X, FFac
C FFac = sqrt(PI)/2
      DATA FFac/0.88622693D0/

      IF (X.EQ.0.D0) THEN
        F = 1.d0
      ELSE
        F = FFac/SQRT(X)*erf(SQRT(x))
      ENDIF
      END





      SUBROUTINE Coul (A, B, SElem, Coul11, Coul12)
C Coulomb matrix element. The subscript 11 indicates that the
C basis functions are centred on the same nucleus, 12 means that
C the kinetic matrix element is calculated for two basis functions

      INCLUDE "globH2"

      DOUBLE PRECISION A, B, Coul11, Coul12, K0, K1, t1, t2, t, F,
     .                 SElem
      
      K0 = PI/(A+B)
      t = (A+B)*Dist*Dist
      Coul11 = -2*K0*(1+F(t))

      K1 = SElem/SQRT(K0)
      t1 = A*A*Dist*Dist/(A+B)
      t2 = B*B*Dist*Dist/(A+B)
      Coul12 = -2*K1*(F(t1)+F(t2))

      END







      SUBROUTINE BuildDensMat
C Calculates the normalies density matrix D=2 C C^T/(C^T C),
C where C is the eigenvector corresponding to the ground state.

      INCLUDE "globH2"

      INTEGER R, S
      DOUBLE PRECISION Norm

      Norm = 0.D0
      DO R=1,2*N
        DO S=1,R-1
          DensMat(R, S) = 2.D0*C(R)*C(S)
          Norm = Norm + DensMat(R,S)*SMatrix(R,S)
        END DO
        DensMat(R,R) = 2.D0*C(R)*C(R)
        Norm = Norm + 0.5D0*DensMat(R,R)*SMatrix(R,R)
      END DO
      Norm = 1.D0/Norm
      DO R=1, 2*N
        DO S=1, R
          DensMat(R,S) = DensMat(R,S)*Norm
        END DO
      END DO
      END





      SUBROUTINE BuildG

C Calculation of the matrix 
C G(T, U) = 0.5 Sum_{R,S} D(R, S) <RT| g |SU>
C Only very limited use has been made of the symmetry 
C of the matrix elements: only the symmetries r<->s and
C t<->u have been used.

      INCLUDE "globH2"

      INTEGER R, S, T, U

      DO R=1,2*N
        DO S=1,R
          GMatrix(R,S) = 0.D0
          DO T=1, 2*N
            DO U=1, T-1
              GMatrix(R, S) = GMatrix(R,S) + DensMat(T,U)
     .                       *QMatrix(R,S,T,U)
            END DO
            GMatrix(R, S) = GMatrix(R,S) + 0.5D0*DensMat(T,T)
     .                     *QMatrix(R,S,T,T)
          END DO
        END DO
      END DO
!      print *, 'g'
!      print '(8F10.5)', GMatrix
!      STOP
      END



      SUBROUTINE CalcFockMatrix
C Calculate Fock matrix: F = H + G,  F(r,s); r>s

      INCLUDE "globH2"

      INTEGER  R, S

      DO R = 1, 2*N
        DO S = 1, R
          FMatrix(R,S) = HMatrix(R,S) + GMatrix(R,S)
        ENDDO
      ENDDO
      END





      SUBROUTINE BuildSuper
C Calculation of the two-electron matrix elements 
C <rt|g|su> Only elements with r>=s are calculated.

      INCLUDE "globH2"

      DOUBLE PRECISION MatElem, A, B, CC, D, P1, P2, P3,
     .                 P, Q, tt, Lambda, f, g

      INTEGER R, S, T, U, MaxU

      DO R=1,2*N
        DO S=1,R
          DO T=1,R
            IF (T .LT. R) THEN
              MaxU = T
            ELSE
              MaxU = S
            ENDIF
            DO U=1,MaxU
              A = Alpha(R) 
              B = Alpha(S) 
              CC = Alpha(T) 
              D = Alpha(U) 
              P1 = A+B
              P2 = CC+D
              P3 = P1+P2
              P = ((2*((R-1)/4)-1)*A+(2*((S-1)/4)-1)*B)/(A+B)*Dist/2
              Q = ((2*((T-1)/4)-1)*CC+(2*((U-1)/4)-1)*D)/(CC+D)*Dist/2
              tt = P1*P2/P3*(P-Q)**2
              Lambda = 2*SQRT(P1*P2/PI/P3)
              G = F(tt)
              MatElem = SMatrix(R,S)*SMatrix(T,U)*Lambda*G
              QMatrix(R,S,T,U) = MatElem
              QMatrix(R,S,U,T) = MatElem
              QMatrix(T,U,R,S) = MatElem
              QMatrix(T,U,S,R) = MatElem
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      END


 


