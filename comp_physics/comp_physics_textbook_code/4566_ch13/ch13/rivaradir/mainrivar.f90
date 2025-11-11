PROGRAM Main
! Main program for a finite element elastic deformation calculation
! involving local refinement
! Some parameters of the calculation are set in the routine InitFEM, which
! is part of the file 'triangfem.f90'. Others should are set at the top
! of the file 'rivara.f90'
USE RIVARA
USE triangfem
USE rivardat
IMPLICIT NONE


CALL InitGrid
CALL InitFEM
DO ! Main relaxation loop
  CALL InitFEMStep
  CALL StressCalc(Displacements, 2*VertPtNum)
  CALL CG(Displacements(1:2*VertPtNum), 2*VertPtNum, BodyForces(1:2*VertPtNum),StiffnessMul)
  CALL StressCalc(Displacements, 2*VertPtNum)
  CALL CG(StressVec(1:3*VertPtNum), 3*VertPtNum, R_Sigma, StressMul)
  CALL CheckStress()
  IF (Accurate) Exit
  CALL DeAllocTriang()
END DO
  CALL DrawFem(3)
print *, 'algorithm converged!'
CALL EndPlot()

CONTAINS

SUBROUTINE CG(X, N, RHS, Multiply)
! Conjugate gradients
  INTEGER :: N, Cnt, I
  DOUBLE PRECISION, INTENT(INOUT) :: X(N), RHS(N)
  DOUBLE PRECISION :: R(N), Z(N), P(N), Q(N), &
                      Rho, NewRho, Alpha, Beta, Error, Tmp, RealRand
  DOUBLE PRECISION, PARAMETER :: MaxRes = 1.D-10
  LOGICAL :: First, MaskArr(N)
INTERFACE
  SUBROUTINE Multiply (X, Y, N)
    INTEGER :: N
    DOUBLE PRECISION :: X(N), Y(N) 
  END SUBROUTINE Multiply
END INTERFACE

! This routine calculates the solution of AX=b and 
! stores the result in X
! SparseMult is a routine which multiplies the matrix A by some arbitrary vector X

  X = 0.D0
  R = RHS
!  print *, 'res 0', N!SQRT(SUM(R*R)), N
  Cnt = 0
  First = .TRUE.
  Error = 10*MaxRes
  DO I=1, N/2
    MaskArr(2*I-1:2*I) = .NOT.FixedPoint(I)
  END DO

  DO WHILE (Error>MaxRes.AND.Cnt<3000)
    Cnt = Cnt + 1
    Z = R
    NewRho = DOT_PRODUCT(R, Z)
    IF (First) THEN
      P = Z
      First = .FALSE.
    ELSE
      Beta = NewRho/Rho
      P = Z + Beta*P
    END IF
    CALL Multiply (P, Q, N)
    Tmp = DOT_PRODUCT(P, Q)
    alpha = NewRho/Tmp
    X = X+Alpha*P
    R = R - alpha*Q
    Rho = NewRho
    Error = SQRT(SUM(R*R, Mask=MaskArr))
!    print '(2F10.5)', X
!    print *, Cnt, error, beta, alpha, newrho
  END DO
  CALL Multiply (X, Q, N)
  R = Q - RHS
!  print *, 'res', SQRT(SUM(R*R)), Cnt
END SUBROUTINE CG

END PROGRAM Main
