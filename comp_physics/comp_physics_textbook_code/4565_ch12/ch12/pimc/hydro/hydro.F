C Thir program is an implementation of the path-integral Monte Carlo
C procedure for the hydrogen atom. Because of the singularity of the 
C Coulomb potential at r=0, the cumulant version of the potential is
C used. This potential is produced by another program, gencumulant, to be 
C found in the subdir "cumulant", and is written to a file CoulPot.
C This program will give an error when CoulPot is not present.
C A sample input file is given: "input"
      PROGRAM Hydro

      include "globvar.f"

      CALL Initialise
      CALL Simulate
      CALL YieldRes
#ifdef Plot
      CALL EndPlot
#endif

      END



      SUBROUTINE Initialise

      include "globvar.f"

      INTEGER I, K
      DOUBLE PRECISION RealRand, Pos, Theta, Phi, PI
      LOGICAL Accept

      print *,"give size (=M), Delta tau, MaxStep (MCS), 
     . InitStep (for relaxation), MaxDisp (max displacement in MC)
     . , HistSize (75)"

      read *, size, epsilon, MaxStep, InitStep, MaxDisp, HistSize
      InvEps = 1.d0/Epsilon/Epsilon

      CALL InitRand(453711)
      PI = 4.d0*DATAN(1.d0)

      DO I=0, Size-1
        Accept = .FALSE.
        DO WHILE (.NOT.(Accept))
          Pos = 3.0D0*RealRand()
          Accept = (Pos*Pos*EXP(-2*Pos)*EXP(2.d0).GT.RealRand())
        ENDDO
        Theta = PI*RealRand()
        Phi = 2*PI*RealRand()
        Chain(I,1) = Pos*COS(Phi)*SIN(Theta)
        Chain(I,2) = Pos*SIN(Phi)*SIN(Theta)
        Chain(I,3) = Pos*COS(Theta)
      ENDDO
      TotE = 0.D0
      OPEN (7, File='gs_density.dat')
#ifdef Plot
      CALL InitPlot('lightblue', 700,700,'out.ps', 1)
      CALL PutStopButton
      CALL ReadPotArray
#endif
      END


      SUBROUTINE ReadPotArray
C Reads the cumulant potential from the file "CoulPot". This potential
C should be constructed for the same epsilon (=Delta tau) as used 
C in this program!
      include "globvar.f"
      INTEGER I, J, K, MaxK, KK, cnt
      OPEN(11,FILE='CoulPot')
      cnt = 0
      DO I=0, 49
        DO J=I, 49
          IF ((I.EQ.0) .OR. (J.EQ.0)) THEN
            MaxK = 0
          ELSE
            MaxK = 20
          ENDIF
          DO K=0, MaxK
            cnt = cnt + 1
            READ(11,*)PotArray(I, J, K)
            PotArray(J,I,K) = PotArray(I,J,K)
            IF (MaxK.EQ.0) THEN
              DO KK=1, 20
                PotArray(I, J, KK) = PotArray(I, J, 0)
                PotArray(J, I, KK) = PotArray(I, J, 0)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      print *, 'Have read file CoulPot ...'
      print *, 'Approaching equilibrium....'
      END



      DOUBLE PRECISION FUNCTION CumPot(Pos1, Pos2)
C Generates the cumulant potential from the file "CoulPot" which has been 
C read in in routine ReadPotArray
      include "globvar.f"
      DOUBLE PRECISION Pos1(3), Pos2(3), Norm1, Norm2, IP, CosTheta,
     .       RestI, RestIC, RestJ, RestJC, RestK, RestKC, Result
      INTEGER I, J, K
      REAL X, erf

      Norm1 = SQRT(Pos1(1)**2+Pos1(2)**2+Pos1(3)**2)
      Norm2 = SQRT(Pos2(1)**2+Pos2(2)**2+Pos2(3)**2)
      IP = Pos1(1)*Pos2(1)+Pos1(2)*Pos2(2)+Pos1(3)*Pos2(3)
      CosTheta = IP/Norm1/Norm2
      I = INT(Norm1*12.5D0)
      RestI = 12.5D0*Norm1-I
      RestIC = 1.d0-RestI
      J = INT(Norm2*12.5D0)
      RestJ = 12.5D0*Norm2-J
      RestJC = 1.d0-RestJ
      K = INT(10.D0*(1.d0-CosTheta))
      RestK = 10.D0*(1.d0-CosTheta)-K
      RestKC = 1.d0-RestK
      IF ((I.LE.49) .AND.(J.LE.49) .AND.(K.LE.20)) THEN
        Result = RestIC*RestJC*RestKC*PotArray(I, J, K) +
     .           RestI *RestJC*RestKC*PotArray(I+1, J, K) +
     .           RestIC *RestJ*RestKC*PotArray(I, J+1, K) +
     .           RestIC *RestJC*RestK*PotArray(I, J, K+1) +
     .           RestIC*RestJ*RestK*PotArray(I, J+1, K+1) +
     .           RestI*RestJC*RestK*PotArray(I+1, J, K+1) +
     .           RestI*RestJ*RestKC*PotArray(I+1, J+1, K) +
     .           RestI*RestJ*RestK*PotArray(I+1, J+1, K+1)
      ELSE
        X = Norm1
        Result = -erf(X)/Norm1
        X = Norm2
        Result = Result-erf(X)/Norm2
        Result = 0.5D0*Result*Epsilon
      ENDIF

      CumPot = Result
      END



      DOUBLE PRECISION FUNCTION Potential(Pos, LPos, RPos)
C Calls the routine CumPot which returns the cumulant potential
      include "globvar.f"
      DOUBLE PRECISION Pos(3), LPos(3), RPos(3),  Tau, R1(3), R2(3),
     .       CoulDens, CumPot, Norm1, Norm2, Pot, Pot1
      INTEGER K

      Tau= Epsilon
      DO K=1, 3
        R1(K) = Pos(K)
        R2(K) = LPos(K)
      ENDDO
      Pot = CumPot(R1, R2)
      DO K=1, 3
        R2(K) = RPos(K)
      ENDDO
      Potential = Pot+CumPot(R1, R2)
      END



      SUBROUTINE CalcSprings (OldPos, NewPos, LPos, RPos, DeltE, I)
C 'Kinetic' energy (energy of the springs) is calculated and returned in DeltE
      include "globvar.f"

      INTEGER K, I

      DOUBLE PRECISION DeltE, OldPos(3), NewPos(3), OldLeft, OldRight,
     .       NewLeft, NewRight, RPos(3), LPos(3)

      OldLeft  = 0.d0
      OldRight = 0.d0
      NewLeft  = 0.d0
      NewRight = 0.d0
      DO K=1, 3 
        OldLeft  = OldLeft  + (OldPos(K)-LPos(K))**2
        OldRight = OldRight + (OldPos(K)-RPos(K))**2
        NewLeft  = NewLeft  + (NewPos(K)-LPos(K))**2
        NewRight = NewRight + (NewPos(K)-RPos(K))**2
      ENDDO
      DeltE = NewLeft+NewRight-OldLeft-OldRight
      END




      SUBROUTINE Simulate
C Main simulation: Metropolis Monte Carlo of the chain
      include "globvar.f"

      INTEGER Step, I, IntRand, PartNum, K, acc

      DOUBLE PRECISION RealRand, DeltE, OldPos(3), NewPos(3), Potential, 
     .       DeltPot, OldR, NewR, TransFac, GuidFunc, LPos(3),
     .       RPos(3), DeltE2

      CALL UpdateHist
      CALL YieldRes2
      DO Step = 1, MaxStep
        DO I=1, Size
          PartNum = IntRand(Size)
          DO K=1, 3
            OldPos(K) = Chain(PartNum,K)
            NewPos(K) = OldPos(K)+MaxDisp*(2*RealRand()-1.d0)
            RPos(K)  = Chain(MOD(PartNum+1,Size),K)
            LPos(K)  = Chain(MOD(PartNum-1+Size,Size),K)
          ENDDO
          CALL CalcSprings (OldPos, NewPos, LPos, RPos, DeltE, PartNum) 
          DeltE = DeltE*0.5D0*InvEps
          OldR = SQRT(OldPos(1)**2+ OldPos(2)**2 + OldPos(3)**2)
          NewR = SQRT(NewPos(1)**2+ NewPos(2)**2 + NewPos(3)**2)
          DeltPot = Potential(NewPos, LPos, RPos)-
     .              Potential(OldPos, LPos, RPos)
          DeltE2 = 0.5D0/Epsilon-(DeltE-DeltPot/Epsilon)
C DeltE2 is used for calculating the total energy...
          DeltE = Epsilon*DeltE+DeltPot
C Note that epsilon is already included in DeltPot 
          TransFac = EXP(-DeltE)
          IF (TransFac.GT.RealRand()) THEN
            DO K=1, 3
              Chain(PartNum,K) = NewPos(K)
            ENDDO
            acc = acc+1
          ENDIF
        ENDDO
        IF (Step .GT. InitStep) THEN
         CALL CalcEner
         TotE = TotE+Energy
          IF (MOD(Step, 100) .EQ. 0) THEN
            print *, '<E>=', TotE/(Step-InitStep)
            IF (Step.GT.InitStep) CALL YieldRes2
          END IF
        END IF
        IF ((Step .GT. InitStep-100).AND.(MOD(Step,20).eq.0)) THEN
          CALL UpdateHist
        ENDIF
      ENDDO
      END



      SUBROUTINE CalcEner
      include "globvar.f"
C Only the direct estimator is used, as the virial estimator requires
C taking the derivative of the cumulant potential. You are invited to
C implement the virial estimator yourself. 
      DOUBLE PRECISION Pos(3), Pot
      INTEGER :: I, IP1

      Energy = 0.D0
      Pot = 0.D0
      DO I=0, Size-1
        IP1 = MOD(I+1, Size)
        DO K=1, 3
          Energy = Energy + (Chain(I, K) - Chain(IP1, K))**2/
     .                      (Epsilon*Epsilon*Size)*0.5D0
          Pos(K) = Chain(I,K)
        END DO
        Pot = Pot + CumPot(Pos, Pos)/Size/Epsilon
      END DO
      Energy = Energy - Pot
      Energy = 1.5D0/Epsilon-Energy
      END 



      SUBROUTINE UpdateHist
C Histogram of the ground state density distrib
C Maximum distance is 4 au. 
      include "globvar.f"

      INTEGER I, Count

      DOUBLE PRECISION R

      DO I=0, Size-1
        R = SQRT(Chain(I,1)**2+Chain(I,2)**2+Chain(I,3)**2)
        Count = NINT(R/4.D0*HistSize*0.5D0)
        IF (Count.LE.HistSize) THEN
          Histogram(Count) = Histogram(Count)+1
        ENDIF
      ENDDO
      END



      SUBROUTINE YieldRes

      include "globvar.f"
 
      INTEGER Count

      DO Count = 0, HistSize
        WRITE (7, '(2F12.6)') Count*4.D0/0.5D0/HistSize, 
     .               dble(Histogram(Count))/MaxStep
      ENDDO
      END



      SUBROUTINE YieldRes2
C Plots the numerical density distribution at 1/T = Epsilon*Size
C against the analytic ground state (shown in red)

      include "globvar.f"
 
      INTEGER Count, Max
      DOUBLE PRECISION Norm, X

      Max = 0
      DO Count = 0, HistSize
        IF (Histogram(Count).GT. Max) THEN
          Max = Histogram(Count)
        ENDIF
      ENDDO
#ifdef Plot
      CALL SetNamedBackground('lightblue')
      CALL Framing (0.d0, 0.d0, 
     .              8.D0, DBLE(Max))
      CALL SetPoint (0.d0, 0.d0)
      DO Count = 0, HistSize
        CALL DrawTo(Count*8.D0/HistSize, 
     .              DBLE(Histogram(Count)))
      ENDDO
      CALL SetFastColor(3)
      Norm = Max*exp(2.d0)
      CALL SetPoint (0.d0, 0.d0)
      DO Count = 1, HistSize
        x = Count*8.D0/HistSize
        CALL DrawTo (x, x*x*exp(-2*x)*Norm)
      ENDDO
      CALL SetNamedColor("black")
#endif
      END



      
      

        
      
