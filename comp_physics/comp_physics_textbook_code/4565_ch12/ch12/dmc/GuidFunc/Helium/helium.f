      PROGRAM HeliumDMC
C*****************************************************************
C* This program estimates the ground state energy of the         *
C* helium atom using a diffusion Monte Carlo                     *
C* procedure with a guide function.                              *
C*                                                               *
C* Program written by Jos Thijssen	                         *
C* Spring 1999                                                   *
C* This program is described in section 12.3.3 of the book       *
C* "Computational Physics" by Jos Thijssen,                      *
C* Cambridge University Press 1999                               *
C*****************************************************************

      DOUBLE PRECISION Alpha, Energy, DeltaT

      INTEGER InitStep, StepNum, TargetNum

      CALL InitWalkers(Alpha, DeltaT, InitStep, StepNum, TargetNum)

      CALL Metropolis(Energy, Alpha, DeltaT, InitStep, StepNum, 
     .                TargetNum)
      print *, Energy
      CLOSE(8)
      END



      SUBROUTINE InitWalkers(Alpha, DeltaT, InitStep, StepNum, 
     .                       TargetNum)
C Initialise the parameters of the procedure and initialse the 
C 3D random walkers to random initial positions. 

      include "globHe"

      INTEGER K, L, I, InitStep, StepNum, TargetNum

      DOUBLE PRECISION RealRand, Alpha, DeltaT

      print *, 'Give nr of steps and nr of equilibration steps'
      read *, StepNum, InitStep
      print *, 'Give target nr of Walkers'
      read *, TargetNum
      IF (TargetNum.GT.MaxWalkNum) THEN
        print *, 'This number should be smaller than', MaxWalkNum
        STOP
      END IF 
      print *, 'Give exponential parameter alpha'
      read *, Alpha
      print *, 'Give time step'
      read *, DeltaT
      CALL InitRand(4537)
      DO L=1, 3
        DO I=1, PartNum
          DO K=1, TargetNum
            Walker(K,I,L) = RealRand()-0.5D0
          END DO
        END DO
      ENDDO
      OPEN (8, FILE='energy.dat')
      END


      SUBROUTINE CalcDist (Pos, R1, R2, R12, R1DR12, R2DR12)
C Calculates distances R1, R2 etc and inner products R1DR12 etc 
C (See other helium programs)

      include "globHe"

      DOUBLE PRECISION Pos(3,PartNum), R1, R2, R12, R1DR12, R2DR12

      INTEGER L

      R1 = SQRT(Pos(1,1)**2+Pos(2,1)**2+Pos(3,1)**2)
      R2 = SQRT(Pos(1,2)**2+Pos(2,2)**2+Pos(3,2)**2)
      R12 = SQRT((Pos(1,1)-Pos(1,2))**2+
     .           (Pos(2,1)-Pos(2,2))**2+
     .           (Pos(3,1)-Pos(3,2))**2)
      R1DR12 = 0.d0
      R2DR12 = 0.d0
      DO L=1, 3
        R1DR12 = R1DR12+Pos(L,1)*(Pos(L,1)-Pos(L,2))
        R2DR12 = R2DR12+Pos(L,2)*(Pos(L,1)-Pos(L,2))
      ENDDO
      END



      SUBROUTINE CalcPsi(Psi, Alpha, R1, R2, R12)
C Calculates the trial wave function for walker K. R1 is
C the distance of electron 1 to the nucleus, and similar
C for R2. R12 is the relative separation between the electrons.
C The R's are returned by this routine. Alpha is the 
C variational parameter. 

      include "globHe"

      DOUBLE PRECISION Alpha, Psi, R1, R2, R12, 
     .                 Denom
 
      Denom = 1.d0/(1.d0+Alpha*R12)
      Psi = EXP(-2*(R1+R2))*EXP(R12*0.5D0*Denom)
      END


      SUBROUTINE CalcF(Pos, F, Alpha, R1, R2, R12)
C Calculates the trial wave function for walker K. R1 is
C the distance of electron 1 to the nucleus, and similar
C for R2. R12 is the relative separation between the electrons.
C The R's are returned by this routine. Alpha is the 
C variational parameter. 

      include "globHe"

      DOUBLE PRECISION Alpha, R1, R2, R12, 
     .                 Pos(3, PartNum), F(3,2), Denom, DU
      INTEGER L
 
      Denom = 1.d0/(1.d0+Alpha*R12)
      DU =  0.5D0*Denom*Denom
      DO L=1, 3
        F(L,1) = 2.d0*(-2*Pos(L,1)/R1-DU*(Pos(L,1)-Pos(L,2))/R12)
        F(L,2) = 2.d0*(-2*Pos(L,2)/R2-DU*(Pos(L,2)-Pos(L,1))/R12)
      END DO
      END



      SUBROUTINE CalcUs(R, Alpha, U12, DU12, DDU12)
C Value and derivative of exponent of correlation term in 
C the Pade-Jastrow wave function
      DOUBLE PRECISION R, Alpha, U12, DU12, DDU12, Denom
      
      Denom = 1.0d0/(1+Alpha*R)
      U12 = 0.5D0*Denom*R
      DU12 = -0.5D0*R*Denom*Denom*Alpha+0.5D0*Denom
      DDU12 = R*Denom*Denom*Denom*(Alpha**2)-Denom*Denom*Alpha
      END



      SUBROUTINE CalcVs(R1, R2, Alpha, V1, V2, DV1, DV2, DDV1, DDV2)
      DOUBLE PRECISION R1, R2, Alpha, V1, V2, DV1, DV2, DDV1, DDV2
C Value and derivative of exponent of hartree 1-electron wave functions
C in the Pade-Jastrow wave function
      
      V1 = -2*R1
      V2 = -2*R2
      DV1 = -2.D0
      DV2 = -2.D0
      DDV1 = 0.D0
      DDV2 = 0.D0
      END



      SUBROUTINE CalcLocal (ELocal, R1, R2, R12, R1DR12, R2DR12, Alpha)

      include "globHe"

      DOUBLE PRECISION ELocal, R1, R2, R12, ETemp, 
     .       U12, DU12, DDU12, V1, V2, DV1, DV2, DDV1, DDV2,
     .       Alpha, R1DR12, R2DR12


C Preliminary calculations....
      CALL CalcUs(R12, Alpha, U12, DU12, DDU12)
      CALL CalcVs(R1, R2, Alpha, V1, V2, DV1, DV2, DDV1, DDV2)

      ETemp = -0.5*(2.d0*(DV1/R1+DV2/R2)+DV1**2+DV2**2+
     .             DDV1+DDV2+
     .             R1DR12 *2.D0* DU12/R12/R1*DV1-
     .             R2DR12 *2.D0* DU12/R12/R2*DV2+
     .             4.d0*DU12/R12+2.D0*DU12**2+2.d0*DDU12)
      ETemp = ETemp-2.d0*(1/R1+1/R2) +1.d0/R12
      ELocal = ETemp
      END
      



      SUBROUTINE Metropolis (Energy, Alpha, DeltaT, InitStep, StepNum, 
     .                       TargetNum)
C The metropolis procedure. StepNum steps are carried out. In one step,
C all WalkNum walkers are moved to new positions. The first InitStep
C steps are used for equilibration. The guide function is parameterised
C by the real number Alpha.
C We move both particles at each step. It is also possible 
C to move a single particle at each step, as was done in the
C Fokker-Planck program for Helium

      include "globHe"

      DOUBLE PRECISION Energy, Alpha, Pos(3, PartNum), 
     .                 NewPos(3, PartNum), DeltaT, GaussWidth, 
     .                 RealRand, RandArr(3,2), 
     .                 OldPsi, NewPsi, R1, R2, R12,
     .                 R1DR12, R2DR12, 
     .                 OldF(3, PartNum), NewF(3, PartNum),  
     .                 Acc, ELocOld, ELocNew, EZero, Q, Weight

      INTEGER K, L, Step, NAcc, NetStep, NewWN, TargetNum,
     .        InitStep, StepNum, WalkNum, NCopies, KK, NTrials

      LOGICAL Alive(MaxWalkNum)
      
      Energy = 0.D0
      NAcc = 0
      NTrials = 0
      WalkNum = TargetNum
      DO K=1, WalkNum
        Alive(K) = .TRUE.
      END DO
      EZero = -2.85D0
      GaussWidth = SQRT(DeltaT)
      DO Step = 1, StepNum
        DO K=1, WalkNum
          DO L=1, 3
            CALL ExpRand(R1, R2)
            RandArr(L,1) = R1*GaussWidth
            RandArr(L,2) = R2*GaussWidth
          END DO
          DO L=1, 3
            Pos(L, 1)    = Walker(K, 1, L)
            Pos(L, 2)    = Walker(K, 2, L)
          END DO
          CALL CalcDist(Pos, R1, R2, R12, R1DR12, R2DR12)
          CALL CalcF(Pos, OldF, Alpha, R1, R2, R12)
          CALL CalcPsi(OldPsi, Alpha, R1, R2, R12)
          CALL CalcLocal (ELocOld, R1, R2, R12, R1DR12, R2DR12, Alpha)
          DO L=1, 3
            NewPos(L,1) = Pos(L,1)+RandArr(L,1)+0.5D0*OldF(L,1)*DeltaT
            NewPos(L,2) = Pos(L,2)+RandArr(L,2)+0.5D0*OldF(L,2)*DeltaT
          END DO
          CALL CalcDist(NewPos, R1, R2, R12, R1DR12, R2DR12)
          CALL CalcF(NewPos, NewF, Alpha, R1, R2, R12)
          CALL CalcPsi(NewPsi, Alpha, R1, R2, R12)
          CALL CalcLocal (ELocNew, R1, R2, R12, R1DR12, R2DR12, Alpha)
          Q = 0.D0
          DO L=1, 3
            Q = Q + 0.5D0*((OldF(L,1)+NewF(L,1))*(0.25D0*DeltaT*
     .                     (OldF(L,1)-NewF(L,1))-NewPos(L,1)+Pos(L,1))+
     .                     (OldF(L,2)+NewF(L,2))*(0.25D0*DeltaT*
     .                     (OldF(L,2)-NewF(L,2))-NewPos(L,2)+Pos(L,2)))
          ENDDO
          Q = EXP(Q)*(NewPsi/OldPsi)**2
          Acc = min(1.d0, Q)
          NTrials = NTrials + 1
c          print *, acc, (NewPsi/OldPsi)**2
          IF (Q.GT.RealRand()) THEN
            Nacc = Nacc+1
            Weight = EXP(-0.5D0*DeltaT*(ELocOld+ELocNew-2*EZero))
            NCopies = INT(Weight+RealRand())
            DO L=1, 3
              Walker(K, 1, L) = NewPos(L, 1)
              Walker(K, 2, L) = NewPos(L, 2)
            END DO
          ELSE
            NCopies = 1
          END IF         
C Append NCopies-1 copies of Walker(K,:) to the end of the array Walker
          DO KK=1, NCopies-1
c            print *, WalkNum
            IF (WalkNum.LT.MaxWalkNum) THEN
              WalkNum = WalkNum+1
              DO L=1, 3
                Walker(WalkNum, 1, L) = NewPos(L, 1)
                Walker(WalkNum, 2, L) = NewPos(L, 2)
              END DO
              Alive(WalkNum) = .TRUE.
            ELSE
              print *, 'Too many walkers'
              STOP
            END IF
          END DO
C Kill the walker if NCopies equals 0
          IF (NCopies.EQ.0) THEN
            Alive(K) = .FALSE.
          END IF
        END DO
C Compactify the array Walker by removing all dead walkers
          NewWN = 0
          DO K=1, WalkNum
            IF (Alive(K)) THEN
              NewWN = NewWN+1
            DO L=1, 3
              Walker(NewWN, 1, L) = Walker(K, 1, L)
              Walker(NewWN, 2, L) = Walker(K, 2, L)
            END DO
            Alive(NewWN)=.TRUE.
          ENDIF
        END DO
        WalkNum = NewWN
        IF (MOD(Step,10).EQ.0) THEN
          print *, Step, walknum, ezero
        END IF
        IF (MOD(Step,100).EQ.0) THEN
c          DeltaT = DeltaT*Nacc/(0.9D0*NTrials)
          print *, Nacc, DeltaT
          Nacc = 0
          NTrials = 0
        END IF
c Estimate new value of EZero
        Ezero= EZero+0.2D0*LOG(DBLE(TargetNum)/WalkNum)
        IF (Step.GT.InitStep) THEN
          Energy = Energy+EZero
          WRITE (8,*) EZero
        END IF        
      END DO
      NetStep = StepNum-InitStep
      Energy = Energy/(NetStep)
      END
