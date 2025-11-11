      PROGRAM Ho3
C Path integral Monte Carlo. See section 12.4 of "Computational 
C Physics", second edition, by Jos Thijssen, Cambridge University
C Press 2007. 
C Harmonic oscillator in 3 D
C The program plots the acceptance rate
C No effort is made to dynamically adjust the maximum displacement 
C to reach a predefined acceptance rate

      include "globvar.f"

      CALL Initialise
      CALL Simulate
      CALL YieldRes

      END



      SUBROUTINE Initialise

      include "globvar.f"

      INTEGER I, K

      print *, "give size, epsilon, MaxStep, MaxDisp, HistSize"

      read *, size, epsilon, MaxStep, MaxDisp, HistSize
C Size: Chain Length
C Epsilon: Delta beta
C MaxStep: Nr of MCS
C InitStep: Nr of Equilibration steps
C MaxDisp: Maximum displacement in MC step
      InvEps = 1.d0/Epsilon/Epsilon

      CALL InitRand(453711)

      DO I=0, Size-1
        DO K=1, 3
          Chain(I,K) = 1.0d0
        ENDDO
      ENDDO
      Energy = 0.D0
      TotE = 0.D0
      CALL CalcChainEner
      END



      SUBROUTINE CalcChainEner

      include "globvar.f"

      INTEGER I

      DOUBLE PRECISION Potential, R
      
      Energy = 0.D0
      DO I=0, Size-1
        R = Chain(I,1)*Chain(I,1)+Chain(I,2)*Chain(I,2)+
     .      Chain(I,3)*Chain(I,3)
        R = SQRT(R)
C This is the energy according to Eq. (12.84)
        Energy = Energy +
     .        2*Potential(R)
      ENDDO
      END


      DOUBLE PRECISION FUNCTION Potential(R)
      DOUBLE PRECISION R
     
      Potential = 0.5*R*R
      END




      SUBROUTINE CalcSprings (OldPos, NewPos, DeltE, I)

      include "globvar.f"

      INTEGER K, I

      DOUBLE PRECISION DeltE, OldPos(3), NewPos(3), OldLeft, OldRight,
     .       NewLeft, NewRight, RPos(3), LPos(3), Fac1, Fac2

      OldLeft  = 0.d0
      OldRight = 0.d0
      NewLeft  = 0.d0
      NewRight = 0.d0
      Fac1 = 1.d0 - 0.5D0*Epsilon
      Fac2 = 1.d0 + 0.5D0*Epsilon
      DO K=1, 3 
        RPos(K)  = Chain(MOD(I+1,Size),K)
        LPos(K)  = Chain(MOD(I-1+Size,Size),K)
        OldLeft  = OldLeft  + (Fac1*OldPos(K)-Fac2*LPos(K))**2
        OldRight = OldRight + (Fac2*OldPos(K)-Fac1*RPos(K))**2
        NewLeft  = NewLeft  + (Fac1*NewPos(K)-Fac2*LPos(K))**2
        NewRight = NewRight + (Fac2*NewPos(K)-Fac1*RPos(K))**2
      ENDDO
      DeltE = NewLeft+NewRight-OldLeft-OldRight
      END
     



      SUBROUTINE Simulate
C Metropolis MC steps
      include "globvar.f"

      INTEGER Step, I, IntRand, PartNum, K, acc

      DOUBLE PRECISION RealRand, DeltE, OldPos(3), NewPos(3), Potential, 
     .       DeltPot, OldR, NewR, TransFac, GuidFunc

      DO Step = 1, MaxStep
        DO I=1, Size
          PartNum = IntRand(Size)
          DO K=1, 3
            OldPos(K) = Chain(PartNum,K)
            NewPos(K) = OldPos(K)+MaxDisp*(2*RealRand()-1.d0)
          ENDDO
          CALL CalcSprings (OldPos, NewPos, DeltE, PartNum) 
          DeltE = DeltE*0.5D0*InvEps
          TransFac = EXP(-Epsilon*DeltE)
          IF (TransFac.GT.1.d0) THEN
            DO K=1, 3
              Chain(PartNum,K) = NewPos(K)
            ENDDO
            acc = acc+1
            Energy = Energy + 2*DeltPot
          ELSE
            IF (TransFac.GT.RealRand()) THEN
              DO K=1, 3
                Chain(PartNum,K) = NewPos(K)
              ENDDO
              acc = acc+1
              Energy = Energy + 2*DeltPot
            ENDIF
          ENDIF
        ENDDO
        TotE = TotE+Energy
        CALL UpdateHist
        if (Mod(step,100).eq.0) then
           print *, dble(acc)/dble(Step)/Size
        endif
      ENDDO
      END


      SUBROUTINE UpdateHist

      include "globvar.f"

      INTEGER I, Count

      DOUBLE PRECISION R

      DO I=0, Size-1
        R = SQRT(Chain(I,1)**2+Chain(I,2)**2+Chain(I,3)**2)
        Count = NINT(R/MaxDisp*HistSize*0.05D0)
        Histogram(Count) = Histogram(Count)+1
      ENDDO
      END

      SUBROUTINE YieldRes

      include "globvar.f"
 
      INTEGER Count

      DO Count = 0, HistSize
        WRITE (7, '(2F12.6)') Count*MaxDisp/0.05D0/HistSize, 
     .               dble(Histogram(Count))/MaxStep
      ENDDO
      print *, TotE/MaxStep/Size
      END



      
      

        
      
