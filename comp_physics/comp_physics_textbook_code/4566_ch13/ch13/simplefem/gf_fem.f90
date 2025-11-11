PROGRAM FEM
IMPLICIT NONE
! Standard finite element (FEM) calculation for an elastic beam in 2D
! attached to a stiff wall.
! The triangular elements are constructed from a structured rectangular 
! grid. 
! The first HighFree vertices of the grid can move; those with 
! indices larger than HighFree are fixed at the wall; therefore 
! they are not allowed to move
! Program written by Jos Thijssen 2004-2007


INTEGER :: VertPtNum, TriangNum, HighFree
DOUBLE PRECISION :: Force(2), TestMat(4,4), GlobForce = 0.03D0, TimeStep = 0.02D0
DOUBLE PRECISION :: MassRho
INTEGER, PARAMETER :: Nx = 50, Ny=6


TYPE TriangleType
  INTEGER :: VertIndex(3)
  DOUBLE PRECISION :: Stiffness(6,6), AMatrix(6,6)
  INTEGER :: FreeNum                    ! Number of non-fixed vertices
END TYPE TriangleType

DOUBLE PRECISION, ALLOCATABLE :: Displacements(:), Momenta(:), RHS(:), MassVec(:)
DOUBLE PRECISION, ALLOCATABLE :: Locations(:,:)
TYPE (TriangleType), ALLOCATABLE :: Triangles(:)
DOUBLE PRECISION, ALLOCATABLE :: BodyForces(:)
INTEGER, ALLOCATABLE :: TractBound(:)   ! Allocatable list of traction boundary pts
INTEGER, ALLOCATABLE :: DispBound(:)    ! Allocatable list of displacement boundary pts
DOUBLE PRECISION :: E, nu

DOUBLE PRECISION :: CMatrix(3,3)

CALL Initialise
CALL DrawFem()
print *, HighFree, NX*NY
CALL CG(Displacements(1:2*HighFree), 2*HighFree,BodyForces(1:2*HighFree))
CALL DrawFem()
CALL Endplot()

CONTAINS


SUBROUTINE Initialise
! Set up the FEM grid and elementary initialisation
! The elastic constants are also initialised
INTEGER :: IX, IY, TriangCount, LowCount, HighCount, Count
INTEGER :: Structured (Nx, Ny) ! Structured grid, where the positions 
                               ! of the grid points are directly related to Nx, Ny
DOUBLE PRECISION :: h = 1.0D0, Xll, Yll, Xur, Yur, ScaleConst

ScaleConst = 2**(1.D0/6)
MassRho = 1/(h*h*0.5*SQRT(3.D0))
h = h*ScaleConst
Force(1) = 0.0D0
Force(2) = -GlobForce
TriangNum = 2*(NX-1)*(NY-1)
VertPtNum = NX*NY
ALLOCATE (Triangles(TriangNum))
ALLOCATE (Displacements(2*VertPtNum))
ALLOCATE (Momenta(2*VertPtNum))
ALLOCATE (RHS(2*VertPtNum))
ALLOCATE (BodyForces(2*VertPtNum))
ALLOCATE (Locations(VertPtNum, 2))
ALLOCATE (MassVec(VertPtNum))
BodyForces = 0
HighFree = NX*NY-NY
LowCount = 0
HighCount = HighFree
DO IX = 1, NX
  DO IY = 1, NY
    IF (.NOT.OnBoundary(IX,IY)) THEN
      LowCount = LowCount + 1
      Count = LowCount 
    ELSE
      HighCount = HighCount + 1
      Count = HighCount
    END IF
    Structured(IX,IY) = Count
    Locations(Count,:) = (/ (IX-1)*h, (IY-1)*h /)
  END DO
END DO
TriangCount = 0
DO IX = 1, NX-1
  DO IY = 1, NY-1
    TriangCount = TriangCount + 1
    Triangles(TriangCount)%VertIndex(:) = (/ Structured(IX,IY), &
                                             Structured(IX+1, IY), &
                                             Structured(IX, IY+1) /)
    TriangCount = TriangCount + 1
    Triangles(TriangCount)%VertIndex(:) = (/ Structured(IX+1,IY), &
                                             Structured(IX+1, IY+1), &
                                             Structured(IX, IY+1) /)
  END DO
END DO
E=428.8D0
Nu = 0.33333333D0
CMatrix(1,:) = (/1.D0,   Nu, 0.D0/)
CMatrix(2,:) = (/  Nu, 1.D0, 0.D0/)
CMatrix(3,:) = (/0.D0, 0.D0, 0.5D0*(1.D0-Nu)/)
CMatrix = CMatrix*E/(1-Nu*Nu)
CMatrix(1,:) = (/1.D0-Nu,   Nu, 0.D0/)
CMatrix(2,:) = (/  Nu, 1.D0-Nu, 0.D0/)
CMatrix(3,:) = (/0.D0, 0.D0, 0.5D0*(1.D0-2*Nu)/)
CMatrix = CMatrix*E/((1+Nu)*(1-2*Nu))
MassVec = 0.D0
DO TriangCount = 1, TriangNum
  CALL Find_Stiffness(Triangles(TriangCount))
END DO
Displacements = 0.D0
CALL InitPlot('lightblue', NX*15,NY*120, 'out.ps', 1)
Xll = -h; Yll = -55*h
Xur = 60*h; Yur = NY*h+5*h
CALL Framing (Xll, Yll, Xur, Yur);
CALL PutStopButton()
END SUBROUTINE Initialise



LOGICAL FUNCTION OnBoundary(IX, IY)
! Check whether point is on the lattice edge

INTEGER :: IX, IY

OnBoundary = (Ix == 1)
END FUNCTION OnBoundary


SUBROUTINE Find_Stiffness(Triangle)
! The stiffness matrix and some other matrices are constructed 
! for the triangle under consideration

TYPE (TriangleType) :: Triangle
DOUBLE PRECISION :: Locs(3,2), Area2
DOUBLE PRECISION :: b(3), c(3),  BMatrix(3,6), TempMatrix(3,6)
INTEGER :: I, J, K, L

DO I=1, 3
  Locs(I,:) = Locations(Triangle%VertIndex(I),:)
END DO

Area2 = Locs(2,1)*Locs(3,2)-Locs(3,1)*Locs(2,2)+&
        Locs(3,1)*Locs(1,2)-Locs(1,1)*Locs(3,2)+&
        Locs(1,1)*Locs(2,2)-Locs(2,1)*Locs(1,2)
b(1) = Locs(2,2)-Locs(3,2)
c(1) = Locs(3,1)-Locs(2,1)
b(2) = Locs(3,2)-Locs(1,2)
c(2) = Locs(1,1)-Locs(3,1)
b(3) = Locs(1,2)-Locs(2,2)
c(3) = Locs(2,1)-Locs(1,1)

DO I=1, 3
  BMatrix(:,2*I-1) = (/ b(I), 0.D0, c(I) /)
  BMatrix(:,2*I) = (/ 0.D0, c(I), b(I) /)
END DO
BMatrix = BMatrix/Area2
Triangle%FreeNum = 3
TempMatrix(1:3,1:6) = MATMUL(CMatrix, BMatrix(1:3,1:6)) 
Triangle%Stiffness = MATMUL(TRANSPOSE(BMatrix(1:3,1:6)), TempMatrix(1:3,1:6))
Triangle%Stiffness = Triangle%Stiffness*0.5D0*Area2 ! Surface integral
DO I=1, 3
  BodyForces(2*Triangle%VertIndex(I)-1) = Area2*Force(1)/3.0D0
  BodyForces(2*Triangle%VertIndex(I))   = Area2*Force(2)/3.0D0
END DO
J=0
DO I=1, 3
  J = J + 1
  L=0
  DO K=1, 3
    L=L+1
    IF (J==L) THEN
      Triangle%AMatrix(2*J, 2*L) = Area2/12.D0
      Triangle%AMatrix(2*J-1, 2*L-1) = Area2/12.D0
    ELSE
      Triangle%AMatrix(2*J, 2*L) = Area2/24.D0
      Triangle%AMatrix(2*J-1, 2*L-1) = Area2/24.D0
    END IF
  END DO
END DO
END SUBROUTINE Find_Stiffness


SUBROUTINE DrawFem
! Draw the current configuration
INTEGER :: TC, VertIndex
DOUBLE PRECISION :: X1, Y1, X2, Y2, X3, Y3

CALL SetNamedBackground('lightblue')
DO TC = 1, TriangNum
  VertIndex = Triangles(TC)%VertIndex(1)
  X1 = Locations(VertIndex,1)+Displacements(2*VertIndex-1)
  Y1 = Locations(VertIndex,2)+Displacements(2*VertIndex)
  VertIndex = Triangles(TC)%VertIndex(2)
  X2 = Locations(VertIndex,1)+Displacements(2*VertIndex-1)
  Y2 = Locations(VertIndex,2)+Displacements(2*VertIndex)
  VertIndex = Triangles(TC)%VertIndex(3)
  X3 = Locations(VertIndex,1)+Displacements(2*VertIndex-1)
  Y3 = Locations(VertIndex,2)+Displacements(2*VertIndex)
CALL Draw(X1, Y1, X2, Y2)
CALL DrawTo (X3, Y3)
CALL DrawTo (X1, Y1)
END DO
CALL Draw(X1, Y1, X2, Y2)
CALL DrawTo (X3, Y3)
CALL DrawTo (X1, Y1)
END SUBROUTINE DrawFem
  



SUBROUTINE StiffnessMul(OldDisp, NewDisp, N)
!Multiplication of the old displacement by the stiffness matrix 

INTEGER :: TC, VC, VertIndex, N, J

DOUBLE PRECISION :: OldNodeVals(6), NewNodeVals(6)
DOUBLE PRECISION :: OldDisp(N), NewDisp(N)

NewDisp = 0.D0
DO TC=1, TriangNum
  OldNodeVals =0.D0
  J=0
  DO VC=1, 3
    VertIndex = Triangles(TC)%VertIndex(VC)
    IF (VertIndex<=HighFree) THEN
      J = J + 1
      OldNodeVals(2*J-1:2*J) = OldDisp(2*VertIndex-1:2*VertIndex)
    END IF
  END DO
  NewNodeVals = MATMUL(Triangles(TC)%Stiffness(1:2*J,1:2*J), OldNodeVals(1:2*J))
  J = 0
  DO VC=1, 3
    VertIndex = Triangles(TC)%VertIndex(VC)
    IF (VertIndex<=HighFree) THEN
      J = J + 1
      NewDisp(2*VertIndex-1:2*VertIndex) = NewDisp(2*VertIndex-1:2*VertIndex)+&
                                           NewNodeVals(2*J-1:2*J)
    END IF
  END DO
END DO
END SUBROUTINE StiffnessMul





      
SUBROUTINE CG(X, N, RHS)
! This routine calculates the solution of AX=b and 
! stores the result in X
! Multiply is a routine which multiplies the matrix A by some arbitrary vector X
INTEGER :: N, Cnt, I
DOUBLE PRECISION, INTENT(INOUT) :: X(N), RHS(N)
DOUBLE PRECISION :: R(N), Z(N), P(N), Q(N), &
                    Rho, NewRho, Alpha, Beta, Error, Tmp
DOUBLE PRECISION, PARAMETER :: MaxRes = 1.D-14
LOGICAL :: First


X = 0.D0
R = RHS
 Cnt = 0
First = .TRUE.
Error = 10*MaxRes
DO WHILE ((Error>MaxRes).AND.(Cnt<1000))
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
  CALL StiffnessMul (P, Q, N)
  Tmp = DOT_PRODUCT(P, Q)
  alpha = NewRho/Tmp
  X = X+Alpha*P
  R = R - alpha*Q
  Rho = NewRho
  Error = SQRT(SUM(R*R))
END DO
CALL StiffnessMul (X, Q, N)
R = Q - RHS
END SUBROUTINE CG





END PROGRAM FEM





