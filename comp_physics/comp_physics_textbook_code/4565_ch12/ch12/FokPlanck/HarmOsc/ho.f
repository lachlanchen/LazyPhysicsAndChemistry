      PROGRAM HarmOscFP
C*****************************************************************
C* This program estimates the ground state energy of the         *
C* harmonic oscillator using a variational quantum Monta Carlo   *
C* procedure based on the Fokker-Planck equation.                *
C* Note that in the algorithm described on page 316, the equation*
C* p = \psi_T(R')/\psi_T(R) must be replaced by                  *
C* p = [\psi_T(R')/\psi_T(R)]^2                                  *
C*                                                               *
C* Program written by Jos Thijssen	                         *
C* Spring 1999                                                   *
C* This program is described in section 12.2.5 of the book       *
C* "Computational Physics" by Jos Thijssen,                      *
C* Cambridge University Press 1999                               *
C*****************************************************************

      DOUBLE PRECISION Alpha, Energy

      INTEGER InitStep, WalkNum, StepNum

      CALL InitWalkers(Alpha, InitStep, StepNum, WalkNum)

      CALL Metropolis(Energy, Alpha, InitStep, StepNum, WalkNum)
      print *, Energy
      CLOSE(8)
      END



      SUBROUTINE InitWalkers(Alpha, InitStep, StepNum, WalkNum)
C Initialise the parameters of the procedure and initialse the 
C random walkers to random initial positions. 

      include "globHOFP"

      INTEGER K, InitStep, StepNum, WalkNum

      DOUBLE PRECISION RealRand, Alpha

      print *, 'Give nr of steps and nr of equilibration steps'
      read *, StepNum, InitStep
      print *, 'Give nr of Walkers'
      read *, WalkNum
      IF (WalkNum.GT.MaxWalkNum) THEN
        print *, 'This number should be smaller than', MaxWalkNum
        STOP
      END IF 
      print *, 'Give exponential parameter alpha'
      read *, Alpha
      CALL InitRand(4537)
      DO K=1, WalkNum
        Walker(K) = RealRand()-0.5D0
      ENDDO
      OPEN (8, File='ground_state.dat')
      END





      SUBROUTINE CalcLocal (ELocal, X, Alpha)
C Calculate the local energy for a walker at position X

      include "globHOFP"

      DOUBLE PRECISION ELocal, Alpha, X

      ELocal = Alpha + X**2*(0.5D0-2.D0*Alpha**2)
      END
      


      SUBROUTINE Metropolis (Energy, Alpha, InitStep, StepNum, WalkNum)
C The metropolis procedure. StepNum steps are carried out. In one step,
C all WalkNum walkers are moved to new positions. The first InitStep
C steps are used for equilibration. The energy is evaluated for the 
C variational parameter Alpha.  

      include "globHOFP"

      DOUBLE PRECISION Energy, Alpha, X, NewX, ELocal,
     .                 GaussWidth, Rand1, Rand2, RandArr(MaxWalkNum), 
     .                 RealRand, DeltaT, OldF, NewF, Acc, Exponent 

      INTEGER K, Step, NAcc, I, Hist(-100:100), NetStep,
     .        InitStep, StepNum, WalkNum
      
      Energy = 0.D0
      NAcc = 0
      DeltaT = 0.1D0
      GaussWidth = SQRT(DeltaT)
      DO I=-100,100
        Hist(I) = 0
      END DO
      DO Step = 1, StepNum
        DO K=1, WalkNum/2
C WalkNum should be EVEN!!!!
          CALL ExpRand(Rand1, Rand2)
          RandArr(2*K-1) = Rand1
          RandArr(2*K)   = Rand2
        END DO
        DO K=1, WalkNum
          X    = Walker(K)
          OldF = -4*Alpha*X
          NewX = Walker(K) + RandArr(K)*GaussWidth + 0.5D0*OldF*DeltaT
          NewF = -4*Alpha*NewX
          Exponent =  0.5D0*(OldF+NewF)*(0.25D0*DeltaT*(OldF-NewF)-
     .                NewX+X)
          Acc = EXP(-2*Alpha*(NewX**2-X**2))*EXP(Exponent)
          IF (RealRand().LT.Acc) THEN
            NAcc = NAcc+1
            Walker(K) = NewX
            X = NewX
          END IF
          IF (Step.GT.InitStep) THEN
            CALL CalcLocal(ELocal, X, Alpha)
            Energy = Energy + ELocal
            IF (ABS(RandArr(K)).LT.2.D0) THEN
              Hist(NINT(RandArr(K)*10)) = Hist(NINT(RandArr(K)*10))+1
            END IF
          END IF
        ENDDO
C Dynamic adaptation of the width of the Gaussian distribution
        IF (MOD(Step,100) .EQ. 0) THEN
          DeltaT = DeltaT*NAcc/(90.D0*WalkNum)
          GaussWidth = SQRT(DeltaT)
          NAcc = 0
        END IF
      ENDDO
      NetStep = StepNum-InitStep
      Energy = Energy/(NetStep*WalkNum)
      DO I=-20,20
        write (8, '(F10.5 F10.5)') I*0.1D0, 
     .                             Hist(I)/DBLE(WalkNum*NetStep)
      ENDDO
      END
      

