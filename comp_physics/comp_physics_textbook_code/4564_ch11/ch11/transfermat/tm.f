      PROGRAM tm
C*****************************************************************
C* This program calculates the largest eigenvalue of the         *
C* transfer matrix of the 2D Ising model for a strip width       *
C* and coupling constant to be specified by the user.            *
C* Program written by Jos Thijssen	                         *
C* Summer 2000                                                   *
C* This program is described in section 11.3   of the book       *
C* "Computational Physics" by Jos Thijssen,                      *
C* Cambridge University Press 1999                               *
C*****************************************************************

C Main program

      include "globals.f"

      CALL Initialise
      CALL DiagTransMat

      END


      SUBROUTINE Initialise
C Initialise problem. Ask strip width and coupling constant J/(K_b T)
C Initialise array ExpR with nearest neighbour Boltzman weights.
      include "globals.f"

      INTEGER I

      print *, 'give strip width'
      read *, Width
      MatSize = 2**Width
      print *, 'Geef J'
      read *, J
      DO I = -4, 4
        ExpR(I) = EXP(J*I*0.5D0)
        print *, expr(i)
      ENDDO
      END

      SUBROUTINE MultVec(OldPhi, NewPhi)
      include "globals.f"
C This is the heart of the program. In this routine, an arbitrary 
C input vector OldPhi is multiplied by the Ising transfer matrix and      
C the resulting vector is stored in NewPhi.

      DOUBLE PRECISION OldPhi(0:MaxSize), NewPhi(0:MaxSize), 
     .                 InterMedPhi(0:MaxSize)

      INTEGER Factor, SecFac, N, M, N1, N2, N3, Count
C Initialise NewPhi to zero
C The `old' vector OldPhi is locally copied to InterMedPhi
      DO M = 0, MatSize-1
        NewPhi(M) = 0.D0
        InterMedPhi(M) = OldPhi(M)
      ENDDO

C First stage: the first (lowest) new spin is added.
      DO N1 = 0, MatSize/8-1
        SecFac = 2**(Width-1)
        N = N1*4
        NewPhi(N)   = InterMedPhi(N)*ExpR(4)+InterMedPhi(N+1)*ExpR(-4)
        NewPhi(N+SecFac) =    InterMedPhi(N+SecFac)*ExpR(2)+
     .                           InterMedPhi(N+SecFac+1)*ExpR(-2)
        NewPhi(N+1) = InterMedPhi(N+1)+InterMedPhi(N)
        NewPhi(N+1+SecFac) =  InterMedPhi(N+1+SecFac)*ExpR(2)+
     .                           InterMedPhi(N+SecFac)*ExpR(-2)
        NewPhi(N+2) = InterMedPhi(N+2)*ExpR(2)+InterMedPhi(N+3)*ExpR(-2)
        NewPhi(N+2+SecFac) =  InterMedPhi(N+2+SecFac)+
     .                           InterMedPhi(N+3+SecFac)
        NewPhi(N+3) = InterMedPhi(N+3)*ExpR(2)+InterMedPhi(N+2)*ExpR(-2)
        NewPhi(N+3+SecFac) =  InterMedPhi(N+3+SecFac)*ExpR(4)+
     .                           InterMedPhi(N+2+SecFac)*ExpR(-4)
      ENDDO

C Second stage: spin 2 to M-1 are added
      Factor = 4
      DO Count = 1, Width-2
        DO M=0, MatSize-1
          InterMedPhi(M) = NewPhi(M)
        ENDDO
        Factor = Factor*2
        SecFac = Factor/8
        DO N1 = 0, MatSize/Factor-1
          DO N2 = 0, SecFac-1
            N3 = N1*Factor+N2
            NewPhi (N3)          = InterMedPhi(N3)*ExpR(4)+
     .                                InterMedPhi(N3+2*SecFac)*ExpR(-2)
            NewPhi (N3+1*SecFac) = InterMedPhi(N3+1*SecFac)*ExpR(2)+
     .                                InterMedPhi(N3+3*SecFac)*ExpR(-4)
            NewPhi (N3+2*SecFac) = InterMedPhi(N3+2*SecFac)+
     .                                InterMedPhi(N3)*ExpR(-2)
            NewPhi (N3+3*SecFac) = InterMedPhi(N3+3*SecFac)*ExpR(2)+
     .                                InterMedPhi(N3+1*SecFac)
            NewPhi (N3+4*SecFac) = InterMedPhi(N3+4*SecFac)*ExpR(2)+ 
     .                                InterMedPhi(N3+6*SecFac)
            NewPhi (N3+5*SecFac) = InterMedPhi(N3+5*SecFac)+
     .                                InterMedPhi(N3+7*SecFac)*ExpR(-2)
            NewPhi (N3+6*SecFac) = InterMedPhi(N3+6*SecFac)*ExpR(2)+
     .                                InterMedPhi(N3+4*SecFac)*ExpR(-4)
            NewPhi (N3+7*SecFac) = InterMedPhi(N3+7*SecFac)*ExpR(4)+ 
     .                                InterMedPhi(N3+5*SecFac)*ExpR(-2)
          ENDDO
        ENDDO
      ENDDO

C Stage three: last spin is added
      DO M=0, MatSize-1
        InterMedPhi(M) = NewPhi(M)
      ENDDO
      Factor = 2**(Width-3)
      SecFac = Factor*2
      DO N2=0, Factor-1
        N3 = N2*2
        NewPhi(N3) = InterMedPhi(N3)*ExpR(4)+InterMedPhi(N3+2*SecFac)
        NewPhi(N3+1) = InterMedPhi(N3+1)*ExpR(2) + 
     .                     InterMedPhi(N3+1+2*SecFac)*ExpR(-2)
        NewPhi(N3+2*SecFac) = InterMedPhi(N3+2*SecFac)+InterMedPhi(N3)*
     .                                                     ExpR(-4)
        NewPhi(N3+2*SecFac+1) = InterMedPhi(N3+2*SecFac+1)*ExpR(2)+
     .                             InterMedPhi(N3+1)*ExpR(-2)
        NewPhi(N3+SecFac) = InterMedPhi(N3+SecFac)*ExpR(2)+
     .                          InterMedPhi(N3+3*SecFac)*ExpR(-2)
        NewPhi(N3+SecFac+1) = InterMedPhi(N3+SecFac+1)+
     .                          InterMedPhi(N3+3*SecFac+1)*ExpR(-4)
        NewPhi(N3+3*SecFac) = InterMedPhi(N3+3*SecFac)*ExpR(2) +
     .                           InterMedPhi(N3+SecFac)*ExpR(-2)
        NewPhi(N3+3*SecFac+1) = InterMedPhi(N3+3*SecFac+1)*ExpR(4) +
     .                             InterMedPhi(N3+SecFac+1)
      ENDDO
      END




      SUBROUTINE DiagTransMat
C Lanczos diagonalisation of the transfer matrix.
C Uses BLAS routines and functions DNRM2, DDOT, DSCAL
C and LAPACK routine DSTEV

      INCLUDE "globals.f"

      INTEGER Index, Level, LevelCount, MaxLevel, IterNum, INFO

      PARAMETER (MaxLevel = 100)

      DOUBLE PRECISION  Phi(MaxSize), PhiPrev(MaxSize), 
     .        Matrix (MaxLevel, MaxLevel), Diag(MaxLevel), 
     .        OffDiag(MaxLevel), Norm, RealRand, B, MaxEig,
     .        Precision, TmpDiag(MaxLevel), 
     .        TmpOffDiag(MaxLevel), PrevEig,
     .        DNRM2, Work(2*MaxLevel-2), Dummy

      PARAMETER (Precision = 1E-7)

      CALL InitRand(44357)

      print *, 'L-value is equal to', Width

      Norm = 0.D0
C Initialise Phi to a random vector. Better initialisations 
C are possible
      DO Index=1, MatSize
         Phi(Index) = (RealRand())-0.5D0
         PhiPrev(Index) = 0.D0
      ENDDO
      Norm = DNRM2(MatSize, Phi, 1)
      Norm = 1.D0/Norm
      CALL DSCAL(MatSize, Norm, Phi, 1)
      B = 0.D0
      OffDiag(1) = 0.D0

      Level = 100

C Build the tridiagonal Lanczos matrix. Tridiagonalise at each step 
C and stop when the largest eigenvalue doesn't change anymore...
      PrevEig = 1.0D0
      LevelCount = 1
      DO WHILE (LevelCount .LE. Level)
        CALL CalcNew(Phi,PhiPrev,Diag,OffDiag,LevelCount,MaxLevel,B)
        DO Index = 1, LevelCount
          TmpDiag(Index) = Diag(Index)
          TmpOffDiag(Index) = OffDiag(Index)
        ENDDO
        CALL DSTEV('n', LevelCount, TmpDiag, TmpOffDiag, Dummy,
     .             1, Work, INFO)
        IF (INFO.LT.0) THEN
          print *, 'wrong arg nr. ', -INFO, 'to routine DSTEV'
          STOP
        ELSE IF (INFO.GT.0) THEN
          print *, 'no convergence reached by routine DSTEV'
        ENDIF
        DO Index = 1, LevelCount
           IF (TmpDiag(Index).GT. MaxEig) THEN
              MaxEig = TmpDiag(Index)
           ENDIF
        ENDDO
        print *, MaxEig
        IF (ABS(MaxEig/PrevEig-1.0D0).LT.Precision) THEN
           IterNum = LevelCount
           LevelCount = Level+100
        ENDIF
        PrevEig = MaxEig
        LevelCount = LevelCount + 1
      ENDDO

      print *, 'Number of iterations:', IterNum
      print *, 'Max. eigenvalue: ', MaxEig
      print *, '1/Width^2 vs Log (MaxEig/Width)'
      WRITE (6, '(2F12.6)') 1.D0/Width/Width, 
     .                        LOG(MaxEig)/Width

      END




      SUBROUTINE CalcNew(Phi, PhiPrev, 
     .                   Diag, OffDiag, LevelCount, MaxLevel, B)
C Calculate the next Lanczos vector
C Uses BLAS routines and functions DDOT, DAXPY, DSCAL and DCOPY
      include "globals.f"

      INTEGER Index, LevelCount, MaxLevel, M

      DOUBLE PRECISION Phi (MaxSize), PhiPrev(MaxSize),
     .       PhiNext(MaxSize), A, B, TmpB, 
     .       NormOf, Diag(MaxLevel), OffDiag(MaxLevel), 
     .       DDOT, DNRM2

      CALL MultVec (Phi, PhiNext)
      A = DDOT(MatSize, PhiNext, 1, Phi, 1)
      Diag(LevelCount) = A
      A = -A
      CALL DAXPY(MatSize, A, Phi, 1, PhiNext, 1)
      B = -B
      CALL DAXPY(MatSize, B, PhiPrev, 1, PhiNext, 1)

      B = DNRM2 (MatSize, PhiNext, 1)
      OffDiag(LevelCount) = B
      TmpB = 1.D0/B
      CALL DSCAL(MatSize, TmpB, PhiNext, 1)
      CALL DCOPY(MatSize, Phi, 1, PhiPrev, 1)
      CALL DCOPY(MatSize, PhiNext, 1, Phi, 1)
      END
      
