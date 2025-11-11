      PROGRAM GenCumul

C*****************************************************************
C* This program generates the effective Coulomb potential between*
C* succesive time steps in the cumulant approximation.           *
C* The program contains the routine CoulDens in which the        *
C* integral of Eq. (12.83a) is evaluated, using a Simpson        *
C* procedure.                                                    *
C* Therefore, the external file Simpson.f is needed.             *
C* Output of this program is a table of Coulomb potentials.      *
C* This table is 3-D: the potential is output for different r_1, *
C* r_2 and angles theta between the two positions.               *
C* The procedure is briefly described in section 12.4.2 of the   *
C* book "Computational Physics" by Jos Thijssen,                 *
C* Cambridge University Press 1999                               *
C* See also the README file in this directory.                   *
C*****************************************************************

      include "globcumul.f"

      DOUBLE PRECISION Result, CoulDens, Y
      INTEGER I, J, K, Select, KMax, cnt
    
      WRITE (6,*), 'Give  Delta tau'
      WRITE (6,*) 'Delta tau is the "time step" used in the'
      WRITE (6,*) 'path integral discretisation'
      READ (5,*) tau
      WRITE (6,*) 'The potential is ouput to a file called'
      WRITE (6,*) '"CoulPot".'

      OPEN (UNIT=7,  FILE='CoulPot')
C Initialise weights for simpson's alternative extended rule
      WW(1) = 17.D0/48.D0
      WW(2) = 59.D0/48.D0
      WW(3) = 43.D0/48.D0
      WW(4) = 49.D0/48.D0
C Initialise PI
      PI = 4.D0*DATAN (1.D0)

C Parameters I and J run over different lengths r_1 and r_2
      cnt = 0
      DO I=0, 49
        print *, i
        DO J=I, 49
          IF ((I.EQ.0) .OR. (J.EQ.0)) THEN
            KMax = 0
          ELSE
            KMax = 20
          ENDIF
C The loop over K generates KMax (=20) values of Cos Theta between 
C -1 and 1
          R1 = I/12.5D0
          R2 = J/12.5d0
          DO K =  0, KMax
            CosTheta = 1.D0-K/10.d0
            Y = CoulDens()
C Y is the cumulant approximation to the potential, which is 
C output to file
            cnt = cnt + 1 
            WRITE (7, '(F12.8)') Y
          ENDDO
        ENDDO
      ENDDO
      print *, 'cnt', cnt
      END





       DOUBLE PRECISION FUNCTION CoulDens()

C Calculation of the Coulomb density matrix using the cumulant expansion.
C In this formulation the density matrix is given as a one-dimensional
C integral which is calculated using Romberg's integration method.

       INCLUDE "globcumul.f"
       EXTERNAL Vs
       DOUBLE PRECISION Result, Vs, X, Y, Result2, 
     .                  Eps, tmed, r
       INTEGER Select, I
       DATA Eps/1.D-12/

C Some preparation
       R12 = SQRT(R1*R1+R2*R2-2*R1*R2*CosTheta)
       IF ((R12.GT.Eps).AND.(R1.GT.Eps)) THEN
         CosAlpha = 0.5D0*(R12*R12+R1*R1-R2*R2)/(R12*R1)
       ELSE
         CosAlpha = 1.D0
       END IF

C Now, the integral is calculated...
C Gauss Legendre is used with 50 points, for which abscissas and 
C weights are given in the file 'gaussleg'
c       CALL Romberg(Vs, 0.D0, Tau, Result)
        CALL Simpson(WW, 0.D0, Tau, 400, Vs, Result)
c       CALL GaussLeg(N, XX, WW, 0.D0, Tau, Vs, Result2)
c       IF (ABS(Result-Result2).GT.1.D-6) THEN
c         print *, 'gevaar', ABS(Result-Result2)
c       END IF

C Result is corrected for added 1/sqrt(t)-like terms
       IF (DABS(R1).GT.1.D-12) THEN
         CoulDens = Result
       ELSE IF(DABS(R2).GT.1.D-12) THEN
          CoulDens = Result-2*sqrt(2.d0*Tau/PI)
       ELSE 
          CoulDens = Result-4*sqrt(2.d0*Tau/PI)
       END IF
    
       END



       DOUBLE PRECISION FUNCTION Vs(TMed)
C Returns the integrand. Normally, this is erf(r/sqrt(2*sigma)), but
C for sigma ->0 or r-> 0, limiting cases have been put

       INCLUDE "globcumul.f"
       DOUBLE PRECISION TMed, Sigma, 
     .        RMed, TVs, D
       REAL X, erf

       D = TMed/Tau*R12
       RMed = R1*R1+D*D-2*R1*D*CosAlpha
       IF (DABS(RMed).LT.1.D-12) THEN
         RMed = 0.D0
       ELSE
         RMed = SQRT(RMed)
       END IF
       Sigma = (Tau-TMed)*TMed/Tau

C Distinguish between (i) both r1 and r2 nonzero, (ii) only r1 zero,
C both r1 and r2 zero...

       IF (DABS(R1).GT.1.D-12) THEN
         IF ((ABS(TMed).LT.1.D-12).OR.(ABS(TMed-Tau).LT.1.D-12)) THEN
           TVs = -1.D0/RMed
         ELSE IF (DABS(RMed).LT.1.D-12) THEN
           TVs = -SQRT(2.D0/PI/Sigma)
         ELSE
           X = RMed/SQRT(2.D0*Sigma)
           TVs  = -DBLE(erf(X))/RMed
         END IF
       ELSE IF (DABS(R2).GT.1.D-12) THEN
         IF (DABS(TMed).LT.1.D-12) THEN
           TVs = 0.0D0
         ELSE IF (DABS(Tau-TMed).LT. 1.D-12) THEN
           TVs = -1.D0/RMed + SQRT(2.D0/(PI*Tau))
         ELSE
           X = RMed/SQRT(2.D0*Sigma)
           TVs  = -DBLE(erf(X))/RMed + SQRT(2.D0/(PI*TMed))
         ENDIF
       ELSE
         IF ((ABS(TMed).LT.1.D-12).OR.(ABS(TMed-Tau).LT.1.D-12)) THEN
           TVs = SQRT(2.D0/(PI*Tau))
         ELSE
           TVs = -SQRT(2.D0/(PI*Sigma)) + SQRT(2.D0/(PI*TMed))+
     .            SQRT(2.D0/(PI*(Tau-TMed)))
         END IF
       END IF
       Vs = TVs
       END 

