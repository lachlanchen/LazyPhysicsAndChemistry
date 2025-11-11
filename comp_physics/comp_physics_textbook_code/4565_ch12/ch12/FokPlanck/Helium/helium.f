      PROGRAM HeliumVMC
C*****************************************************************
C* This program estimates the ground state energy of the         *
C* helium atom using a variational quantum Monte Carlo           *
C* procedure based on the Fokker-Planc equation.                 *
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
      
      END



      SUBROUTINE InitWalkers(Alpha, InitStep, StepNum, WalkNum)
C Initialise the parameters of the procedure and initialse the 
C 3D random walkers to random initial positions. 

      include "globHe"

      INTEGER K, L, I, InitStep, StepNum, WalkNum

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
      DO L=1, 3
        DO I=1, PartNum
          DO K=1, WalkNum
            Walker(K,I,L) = RealRand()-0.5D0
          END DO
        END DO
      ENDDO
      END





      SUBROUTINE CalcPsiAndF(Pos, I, Psi, F, Alpha, R1, R2, R12)
C Calculates the trial wave function for walker K. R1 is
C the distance of electron 1 to the nucleus, and similar
C for R2. R12 is the relative separation between the electrons.
C The R's are returned by this routine. Alpha is the 
C variational parameter. 

      include "globHe"

      DOUBLE PRECISION Alpha, Psi, R1, R2, R12, 
     .                 Pos(3, PartNum), F(3), Denom, DU
      INTEGER I, L
 
      R1 = SQRT(Pos(1,1)**2+Pos(2,1)**2+Pos(3,1)**2)
      R2 = SQRT(Pos(1,2)**2+Pos(2,2)**2+Pos(3,2)**2)
      R12 = SQRT((Pos(1,1)-Pos(1,2))**2+
     .           (Pos(2,1)-Pos(2,2))**2+
     .           (Pos(3,1)-Pos(3,2))**2)
      Denom = 1.d0/(1.d0+Alpha*R12)
      Psi = EXP(-2*(R1+R2))*EXP(R12*0.5D0*Denom)
      DU =  -0.5D0*R12*Denom*Denom*Alpha+0.5D0*Denom
      IF (I.EQ.1) THEN
        DO L=1, 3
          F(L) = 2.d0*(-2*Pos(L,1)/R1+DU*(Pos(L,1)-Pos(L,2))/R12)
        END DO
      ELSE
        DO L=1, 3
          F(L) = 2.d0*(-2*Pos(L,2)/R2+DU*(Pos(L,2)-Pos(L,1))/R12)
        END DO
      END IF
      END



      SUBROUTINE CalcUs(R, Alpha, U12, DU12, DDU12)
      DOUBLE PRECISION R, Alpha, U12, DU12, DDU12, Denom
      
      Denom = 1.0d0/(1+Alpha*R)
      U12 = 0.5D0*Denom*R
      DU12 = -0.5D0*R*Denom*Denom*Alpha+0.5D0*Denom
      DDU12 = R*Denom*Denom*Denom*(Alpha**2)-Denom*Denom*Alpha
      END



      SUBROUTINE CalcVs(R1, R2, Alpha, V1, V2, DV1, DV2, DDV1, DDV2)
      DOUBLE PRECISION R1, R2, Alpha, V1, V2, DV1, DV2, DDV1, DDV2
      
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
      



      SUBROUTINE Metropolis (Energy, Alpha, InitStep, StepNum, WalkNum)
C The metropolis procedure. StepNum steps are carried out. In one step,
C all WalkNum walkers are moved to new positions. The first InitStep
C steps are used for equilibration. The energy is evaluated for the 
C variational parameter Alpha.  

      include "globHe"

      DOUBLE PRECISION Energy, Alpha, Pos(3, PartNum), 
     .                 NewPos(3, PartNum), DeltaT, GaussWidth, 
     .                 ELocal, MaxStep, RealRand, RandArr(6), 
     .                 OldPsi, NewPsi, R1, R2, R12, S1, S2, S12,
     .                 R1DR12, R2DR12, OldF(3), NewF(3), Exponent, Acc 

      INTEGER K, L, Step, NAcc, I, NetStep, 
     .        InitStep, StepNum, WalkNum
      
      Energy = 0.D0
      NAcc = 0
      DeltaT = 0.5D0
      GaussWidth = SQRT(DeltaT)
      DO Step = 1, StepNum
        DO K=1, WalkNum
          DO L=1, 3
            CALL ExpRand(R1, R2)
            RandArr(2*L-1) = R1*GaussWidth
            RandArr(2*L) = R2*GaussWidth
          END DO
          DO I=1, PartNum
            DO L=1, 3
              Pos(L, I)     = Walker(K, I, L)
              Pos(L, 3-I)   = Walker(K, 3-I, L)
              NewPos(L, 3-I)= Walker(K, 3-I, L)
            END DO
            CALL CalcPsiAndF(Pos,   I, OldPsi, OldF, Alpha, R1, R2, R12)
            DO L=1, 3
              NewPos(L, I) = Pos(L, I) + RandArr(3*(I-1)+L) +
     .                       0.5D0*OldF(L)*DeltaT
            END DO
            CALL CalcPsiAndF(NewPos,I, NewPsi, NewF, Alpha, S1, S2, S12)
            Exponent = 0.0D0
            DO L=1, 3
              Exponent = Exponent + 0.5D0*(OldF(L)+NewF(L))*
     .                   (0.25D0*DeltaT*(OldF(L)-NewF(L))
     .                    -NewPos(L,I)+Pos(L,I))
            END DO
            Acc = (NewPsi/OldPsi)**2*EXP(Exponent)
            IF (RealRand().LT.Acc) THEN
              NAcc = NAcc+1
              DO L=1, 3
                Walker(K, I, L) = NewPos(L, I)
                Pos(L, I)       = NewPos(L, I)
              END DO
              R1  = S1
              R2  = S2
              R12 = S12
            END IF
            IF (Step.GT.InitStep) THEN
              R1DR12 = Pos(1, 1)*(Pos(1, 1)-Pos(1, 2))+
     .                 Pos(2, 1)*(Pos(2, 1)-Pos(2, 2))+
     .                 Pos(3, 1)*(Pos(3, 1)-Pos(3, 2))
              R2DR12 = Pos(1, 2)*(Pos(1, 1)-Pos(1, 2))+
     .                 Pos(2, 2)*(Pos(2, 1)-Pos(2, 2))+
     .                 Pos(3, 2)*(Pos(3, 1)-Pos(3, 2))
              CALL CalcLocal(ELocal, R1, R2, R12, R1DR12, R2DR12, Alpha)
              Energy = Energy + ELocal
            END IF
          ENDDO
        END DO
C Dynamic adaptation of the time step -- 90 % of the moves should 
C be accepted
        IF (MOD(Step,100) .EQ. 0) THEN
          DeltaT = DeltaT*NAcc/(90.D0*WalkNum)
c          print *, DeltaT, NAcc, 90*WalkNum
          GaussWidth = SQRT(DeltaT)
          NAcc = 0
        END IF
      END DO
      NetStep = StepNum-InitStep
      Energy = Energy/(2*NetStep*WalkNum)
      END
      

