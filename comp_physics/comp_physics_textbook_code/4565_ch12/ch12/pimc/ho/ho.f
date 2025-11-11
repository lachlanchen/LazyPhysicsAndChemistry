      PROGRAM Ho
C Path integral Monte Carlo. See section 12.4 of "Computational 
C Physics", second edition, by Jos Thijssen, Cambridge University
C Press 2007. 
C Harmonic oscillator in 1 D
C The program prints the acceptance rate
C No effort is made to dynamically adjust the maximum displacement 
C to reach a predefined acceptance rate

      include "globvar1.f"

      CALL Initialise
      CALL Simulate
      CALL YieldRes
      CLOSE(7)
      CLOSE(8)
      END



      SUBROUTINE Initialise

      include "globvar1.f"

      INTEGER I, K, Count
      DOUBLE PRECISION RealRand
C Size: Chain Length
C Epsilon: Delta beta
C MaxStep: Nr of MCS
C InitStep: Nr of Equilibration steps
C MaxDisp: Maximum displacement in MC step
      print *,"give size, epsilon, MaxStep, InitStep, MaxDisp, 
     .         HistSize"

      read *, size, epsilon, MaxStep, InitStep, MaxDisp, 
     .        HistSize
      InvEps = 1.d0/Epsilon**2

      CALL InitRand(453711)

      DO I=0, Size-1
        Chain(I) = RealRand()
      ENDDO
      Energy = 0.D0
      TotE = 0.D0
C Calculate energy before starting simulation
      CALL CalcChainEner
C Histogram will contain the ground state prob distr
      DO Count =-HistSize, HistSize
        Histogram(Count)=0
      ENDDO
      OPEN (7, File='gs_distr.dat')
      OPEN (8, File='Energy.dat')
      END



      SUBROUTINE CalcChainEner

      include "globvar1.f"

      INTEGER I, IP1

      DOUBLE PRECISION Potential, R
C This routine calculates the total energy, according to
C the virial estimator, see Eq. (12.84).
      Energy = 0.D0
      DO I=0, Size-1
        R = Chain(I)
        IP1 = MOD(I+1, Size)
C Direct Estimator:
C        Energy = Energy +0.5D0*InvEps*
C     .               (Chain(IP1)-Chain(I))**2/Size-
C     .               0.5*R**2/Size
C Virial estimator:
        Energy = Energy + 
     .        2*Potential(R)
      ENDDO
C Virial estimator:
      Energy = Energy/Size
C Direct Estimator:
C      Energy = 0.5D0/Epsilon - Energy
      END


      DOUBLE PRECISION FUNCTION Potential(R)
      DOUBLE PRECISION R
     
      Potential = 0.5d0*R*R
      END




      SUBROUTINE CalcSprings (OldPos, NewPos, DeltE, I)
C Calculates energy difference for a step from OldPos to NewPos

      include "globvar1.f"

      INTEGER I

      DOUBLE PRECISION DeltE, OldPos, NewPos, OldLeft, OldRight,
     .       NewLeft, NewRight, RPos, LPos, Fac1, Fac2,
     .       Potential

      RPos  = Chain(MOD(I+1,Size))
      LPos  = Chain(MOD(I-1+Size,Size))
      OldLeft  = (OldPos-LPos)**2
      OldRight = (OldPos-RPos)**2
      NewLeft  = (NewPos-LPos)**2
      NewRight = (NewPos-RPos)**2
      DeltE = NewLeft+NewRight-OldLeft-OldRight
      END



      SUBROUTINE Simulate
C Standard Metropolis procedure for the chain
      include "globvar1.f"

      INTEGER Step, I, IntRand, PartNum, acc

      DOUBLE PRECISION RealRand, DeltE, OldPos, NewPos, Potential, 
     .       DeltPot, TransFac

      acc = 0
      DO Step = 1, MaxStep
        DO I=0, Size-1
          PartNum = IntRand(Size)
          OldPos = Chain(PartNum)
          NewPos = OldPos+MaxDisp*(2*RealRand()-1.d0)
          CALL CalcSprings (OldPos, NewPos, DeltE, PartNum) 
          DeltE = DeltE*0.5D0*InvEps
          DeltE = DeltE + Potential(NewPos)-Potential(OldPos)
          TransFac = EXP(-Epsilon*DeltE)
          IF (TransFac.GT.RealRand()) THEN
            Chain(PartNum) = NewPos
            acc = acc+1
          ENDIF
        ENDDO
        IF (Step.GT.InitStep) THEN
          CALL UpdateHist
          CALL CalcChainEner
          TotE = TotE+Energy
          WRITE (8, '(F12.7)') Energy
          if (Mod(step,100).eq.0) then
             print *, 'acc', dble(acc)/dble(Step)/Size
          endif
        ENDIF
      ENDDO
      END


      SUBROUTINE UpdateHist

      include "globvar1.f"

      INTEGER I, Count

      DOUBLE PRECISION R

      DO I=0, Size-1
        R = Chain(I)
        Count = NINT(R/MaxDisp*HistSize*0.10D0)
        IF (ABS(Count)<HistSize) THEN
          Histogram(Count) = Histogram(Count)+1
        END IF
      ENDDO
      END

      SUBROUTINE YieldRes
! Print histogram for the ground state dens distr
      include "globvar1.f"
 
      INTEGER Count

      DO Count = -HistSize, HistSize
        WRITE (7, '(2F12.6)') Count*MaxDisp/0.10D0/HistSize, 
     .               dble(Histogram(Count))/MaxStep
      ENDDO
      print *, TotE/(MaxStep-InitStep)
      END



      
      

        
      
