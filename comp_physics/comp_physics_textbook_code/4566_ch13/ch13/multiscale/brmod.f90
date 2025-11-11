MODULE BrMod
IMPLICIT NONE
! Elastic deformation of a plane structure.
! In this program, a finite element calculation is carried out for a triangular mesh.
! We construct the finite element metrices B from the three corner points.
! Per node we have its offset position and we calculate its displacement
! All displacements are stored in a system vector, r
! The stiffness matrix is a 3x3 matrix. It originates from the linear 
! relations between the six degrees of freedom (the u's and the v's at the corners) 
! which are not independent.

INTEGER :: VertPtNum, TriangNum, HighFree, NativeHigh
DOUBLE PRECISION :: Force(2), TestMat(4,4), Fac2, ScaleConst
DOUBLE PRECISION :: MassRho
INTEGER, PARAMETER :: Nx = 24, Ny=23


TYPE TriangleType
  INTEGER :: VertIndex(3)
  DOUBLE PRECISION :: Stiffness(6,6), AMatrix(6,6)
  INTEGER :: FreeNum                    ! Number of non-fixed vertices
  LOGICAL :: Border
END TYPE TriangleType

DOUBLE PRECISION, ALLOCATABLE :: Displacements(:), Momenta(:), RHS(:), MassVec(:)
DOUBLE PRECISION, ALLOCATABLE :: Locations(:,:)
TYPE (TriangleType), ALLOCATABLE :: Triangles(:)
DOUBLE PRECISION, ALLOCATABLE :: BodyForces(:)
INTEGER, ALLOCATABLE :: TractBound(:)   ! Allocatable list of traction boundary pts
INTEGER, ALLOCATABLE :: DispBound(:)    ! Allocatable list of displacement boundary pts
DOUBLE PRECISION :: E, nu, Lx, Ly

DOUBLE PRECISION :: CMatrix(3,3)


CONTAINS

INTEGER FUNCTION OutputNy()
  
  OutputNY = Ny
END FUNCTION OutputNy

SUBROUTINE Initialise
INTEGER :: IX, IY, TriangCount, LowCount, HighCount, Count
INTEGER :: Structured (Nx+1, Ny)
DOUBLE PRECISION :: h = 0.996891D0, Xll, Yll, Xur, Yur, h2

ScaleConst = 2**(1.D0/6)
h = h*ScaleConst
h2 = 0.5*h*sqrt(3.d0)
MassRho = Nx*Ny/((Nx-1)*(Ny-1)*h*h2)
Force(1) = 0.0D0
Force(2) = -0.0D0
TriangNum = 2*NX*(NY-1)
VertPtNum = (NX+1)*NY
Ly = h2*Ny
ALLOCATE (Triangles(TriangNum))
ALLOCATE (Displacements(2*VertPtNum))
ALLOCATE (Momenta(2*VertPtNum))
ALLOCATE (RHS(2*VertPtNum))
ALLOCATE (BodyForces(2*VertPtNum))
ALLOCATE (Locations(VertPtNum, 2))
ALLOCATE (MassVec(VertPtNum))
BodyForces = 0
NativeHigh = (NX)*NY-(NY+1)/2
HighFree = (NX+1)*NY-(NY+1)/2
LowCount = 0
HighCount = HighFree
print *, 'In fem, Ny is ', Ny
DO IX = 1, NX+1
  DO IY = 1, NY
    IF (.NOT.OnBoundary(IX,IY)) THEN
      LowCount = LowCount + 1
      Count = LowCount 
    ELSE
      HighCount = HighCount + 1
      Count = HighCount
    END IF
    Structured(IX,IY) = Count
    IF (MOD(IY,2)==1) THEN
      Locations(Count,:) = (/ (IX-1)*h, (IY-1)*h2 /)
    ELSE
      Locations(Count,:) = (/ (IX-0.5D0)*h, (IY-1)*h2 /)
    END IF
  END DO
END DO
! CALL CopyMDToFem (MDPos, IY)
TriangCount = 0
DO IX = 1, NX
  DO IY = 1, NY-1
    TriangCount = TriangCount + 1
    IF (MOD(IY,2).EQ.1) THEN
      print *, ix, nx
      Triangles(TriangCount)%VertIndex(:) = (/ Structured(IX,IY), &
                                               Structured(IX+1, IY), &
                                               Structured(IX, MOD(IY,Ny)+1) /)
      IF (IX==NX) THEN 
        Triangles(TriangCount)%Border = .TRUE.
      ELSE
        Triangles(TriangCount)%Border = .FALSE.
      END IF
      TriangCount = TriangCount + 1
      Triangles(TriangCount)%VertIndex(:) = (/ Structured(IX+1,IY), &
                                               Structured(IX+1, MOD(IY,Ny)+1), &
                                               Structured(IX, MOD(IY,Ny)+1) /)
      IF (IX==NX) THEN 
        Triangles(TriangCount)%Border = .TRUE.
      ELSE
        Triangles(TriangCount)%Border = .FALSE.
      END IF
    ELSE
      Triangles(TriangCount)%VertIndex(:) = (/ Structured(IX,IY), &
                                               Structured(IX+1, MOD(IY,Ny)+1), &
                                               Structured(IX, MOD(IY,Ny)+1) /)
      IF (IX==NX) THEN 
        Triangles(TriangCount)%Border = .TRUE.
      ELSE
        Triangles(TriangCount)%Border = .FALSE.
      END IF
      TriangCount = TriangCount + 1
      Triangles(TriangCount)%VertIndex(:) = (/ Structured(IX,IY), &
                                               Structured(IX+1, IY), &
                                               Structured(IX+1, MOD(IY,Ny)+1) /)
      IF (IX==NX) THEN 
        Triangles(TriangCount)%Border = .TRUE.
      ELSE
        Triangles(TriangCount)%Border = .FALSE.
      END IF
    END IF
  END DO
END DO
CMatrix = 0.D0
E=80.0D0
Nu = 0.3333333333D0
CMatrix(1,:) = (/1.D0,   Nu, 0.D0/)
CMatrix(2,:) = (/  Nu, 1.D0, 0.D0/)
CMatrix(3,:) = (/0.D0, 0.D0, 0.5D0*(1.D0-Nu)/)
CMatrix = CMatrix*E/(1-Nu*Nu)
MassVec = 0.D0
DO TriangCount = 1, TriangNum
  CALL Find_Stiffness(Triangles(TriangCount))
END DO
Displacements = 0.D0
CALL InitPlot('lightblue', NX*80,NY*40, 'out.ps', 1)
Xll = -h; Yll = -h
Xur = 100*h; Yur = NY*h
CALL Framing (Xll, Yll, Xur, Yur);
CALL PutStopButton()
END SUBROUTINE Initialise


LOGICAL FUNCTION OnBoundary(IX, IY)

INTEGER :: IX, IY

OnBoundary = ((Ix == 1).AND.(MOD(IY,2)==1))
END FUNCTION OnBoundary


SUBROUTINE Find_Stiffness(Triangle)

TYPE (TriangleType) :: Triangle
DOUBLE PRECISION :: Locs(3,2), Area2
DOUBLE PRECISION :: b(3), c(3),  BMatrix(3,6), TempMatrix(3,6)
INTEGER :: I, J, K, L
LOGICAL :: BorI, BorJ

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
DO K=1, 3
  IF (ABS(B(k)-Ly)<ABS(B(k))) THEN
    B(k) = B(k) - Ly
  ELSE IF (ABS(B(k)+Ly)<ABS(B(k))) THEN
    B(K) = B(K) + Ly
  END IF
END DO
Area2 = c(3)*b(2)-c(2)*b(3)

J = 0
DO I=1, 3
    BMatrix(:,2*I-1) = (/ b(I), 0.D0, c(I) /)
    BMatrix(:,2*I) = (/ 0.D0, c(I), b(I) /)
END DO
BMatrix = BMatrix/Area2
! Triangle%FreeNum = J
TempMatrix(1:3,1:2*3) = MATMUL(CMatrix, BMatrix(1:3,1:2*3)) 
Triangle%Stiffness = MATMUL(TRANSPOSE(BMatrix(1:3,1:2*3)), TempMatrix(1:3,1:2*3))
Triangle%Stiffness = Triangle%Stiffness*0.5D0*Area2 ! Surface integral
IF (Triangle%Border) THEN
  Triangle%Stiffness = Triangle%Stiffness*0.5D0
  print *, 'piep'
END IF
DO I=1, 3
  BodyForces(2*Triangle%VertIndex(I)-1) = Area2*Force(1)/6.0D0
  BodyForces(2*Triangle%VertIndex(I))   = Area2*Force(2)/6.0D0
  MassVec(Triangle%VertIndex(I)) = MassVec(Triangle%VertIndex(I))+Area2*MassRho/6.D0
  ! ie, one third of the triangle's mass is moved to each vertex
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
INTEGER :: TC, VertIndex
DOUBLE PRECISION :: X1, Y1, X2, Y2, X3, Y3
LOGICAL :: Border

!CALL SetNamedBackground('lightblue')
CALL SetNamedColor('black')
DO TC = 1, TriangNum
  Border = .FALSE.
  VertIndex = Triangles(TC)%VertIndex(1)
  IF (VertIndex>NativeHigh .AND. VertIndex<=HighFree) Border = .TRUE.
  X1 = Locations(VertIndex,1)+Displacements(2*VertIndex-1)
  Y1 = Locations(VertIndex,2)+Displacements(2*VertIndex)
  VertIndex = Triangles(TC)%VertIndex(2)
  IF (VertIndex>NativeHigh .AND. VertIndex<=HighFree) Border = .TRUE.
  X2 = Locations(VertIndex,1)+Displacements(2*VertIndex-1)
  Y2 = Locations(VertIndex,2)+Displacements(2*VertIndex)
  VertIndex = Triangles(TC)%VertIndex(3)
  IF (VertIndex>NativeHigh .AND. VertIndex<=HighFree) Border = .TRUE.
  X3 = Locations(VertIndex,1)+Displacements(2*VertIndex-1)
  Y3 = Locations(VertIndex,2)+Displacements(2*VertIndex)
  IF (Border) CALL SetNamedColor('red')
  CALL Draw(X1, Y1, X2, Y2)
  CALL DrawTo (X3, Y3)
  CALL DrawTo (X1, Y1)
END DO
END SUBROUTINE DrawFem


SUBROUTINE DrawFem2
INTEGER :: TC, VertIndex, IX, IY, MaxPos
DOUBLE PRECISION :: Aver1, Aver2, Mass1, Mass2, MaxVal

! CALL SetNamedBackground('lightblue')
CALL SetNamedColor('black')
VertIndex = 0
CALL SetPoint(0.D0, 0.D0)
MaxVal = 0.D0
MaxPos = 0
DO IX=1, NX
  Aver1 = 0.D0
  Aver2 = 0.D0
  Mass1=0.D0; Mass2=0.D0
  DO IY = 1, NY
    IF (.NOT.OnBoundary(IX,IY)) THEN
      VertIndex = VertIndex + 1
      IF (MOD(IY,2).EQ.0) THEN
        Mass2 = Mass2 + MassVec(VertIndex)
        Aver2 = Aver2 + Displacements(2*VertIndex)
      ELSE
        Mass1 = Mass1 + MassVec(VertIndex)
        Aver1 = Aver1 + Displacements(2*VertIndex)
      END IF
    END IF
  END DO
  IF ((Aver2+Aver1)/(Mass2+Mass1)>MaxVal) THEN
    MaxVal = (Aver2+Aver1)/(Mass2+Mass1)
    MaxPos = IX
  END IF
  IF (Mass2>0) THEN
    CALL DrawTo(DBLE(IX)*1.12246D0-0.5D0, 10000*(Aver2+Aver1)/(Mass2+Mass1))
  END IF
END DO
END SUBROUTINE DrawFem2
  

  


SUBROUTINE StiffnessMul(OldDisp, NewDisp, N)

INTEGER :: TC, VC, VertIndex, N, J

DOUBLE PRECISION :: OldNodeVals(6), NewNodeVals(6)
DOUBLE PRECISION :: OldDisp(N), NewDisp(N)

NewDisp = 0.D0
!print *, triangnum, nativehigh, highfree
DO TC=1, TriangNum
  OldNodeVals =0.D0
  J=0
  DO VC=1, 3
    VertIndex = Triangles(TC)%VertIndex(VC)
!    IF (VertIndex<=HighFree) THEN
      J = J + 1
      OldNodeVals(2*J-1:2*J) = OldDisp(2*VertIndex-1:2*VertIndex)
!    END IF
  END DO
  NewNodeVals = MATMUL(Triangles(TC)%Stiffness(1:2*J,1:2*J), OldNodeVals(1:2*J))
!  IF (TC>2*(Nx-1)*(Ny-1)) THEN
!    print *, Triangles(TC)%AMatrix(1:2*J,1:2*J)
!  END IF
  J = 0
  DO VC=1, 3
    VertIndex = Triangles(TC)%VertIndex(VC)
!    IF (VertIndex<=HighFree) THEN
      J = J + 1
      NewDisp(2*VertIndex-1:2*VertIndex) = NewDisp(2*VertIndex-1:2*VertIndex)+&
                                           NewNodeVals(2*J-1:2*J)
!    END IF
  END DO
END DO
END SUBROUTINE StiffnessMul



SUBROUTINE MassMul(OldDisp, NewDisp, N)

INTEGER :: TC, VC, VertIndex, N, J

DOUBLE PRECISION :: OldNodeVals(6), NewNodeVals(6)
DOUBLE PRECISION, INTENT(IN) :: OldDisp(N)
DOUBLE PRECISION, INTENT(OUT) :: NewDisp(N)

NewDisp = 0.D0
print *, triangnum, nativehigh, highfree
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
  NewNodeVals(1:2*J) = MATMUL(Triangles(TC)%AMatrix(1:2*J,1:2*J), OldNodeVals(1:2*J))
  IF (TC>2*NativeHigh) THEN
    print *, Triangles(TC)%AMatrix(1:2*J,1:2*J)
  END IF
!  NewNodeVals = 0.002*OldNodeVals
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
END SUBROUTINE MassMul
      


SUBROUTINE InitDyn(DeltaT)

DOUBLE PRECISION :: X, h, DeltaT
INTEGER :: StepCount, IY, I

Displacements(1:2*HighFree) = 0.D0
Fac2=0.5D0*DeltaT**2
Displacements = 0.D0
Momenta = -0.D0
CALL StiffnessMul(Displacements(1:2*HighFree), RHS(1:2*HighFree), 2*HighFree)
RHS(1:2*HighFree) = -RHS(1:2*HighFree)+BodyForces(1:2*HighFree)
DO I=1, HighFree
  RHS(2*I-1:2*I) = RHS(2*I-1:2*I)/MassVec(I)
END DO
END SUBROUTINE InitDyn

SUBROUTINE DoFEMStep1(DeltaT)
  DOUBLE PRECISION :: DeltaT
  INTEGER :: I

  Displacements(1:2*NativeHigh) = Displacements(1:2*NativeHigh)+&
                                DeltaT*Momenta(1:2*NativeHigh)+&
                                Fac2*RHS(1:2*NativeHigh)
  Momenta(1:2*NativeHigh) = Momenta(1:2*NativeHigh) + 0.5D0*DeltaT*RHS(1:2*NativeHigh)
!  DO I=1, NativeHigh
!     print *, I, NativeHigh
!     print '(4F14.10)', Displacements(2*I), Momenta(2*I), RHS(2*I), Locations(I, 1)
!  END DO
END SUBROUTINE DoFEMStep1

SUBROUTINE CalcFEMForce
  RHS = 0.D0
  CALL StiffnessMul(Displacements(1:2*HighFree), RHS(1:2*HighFree), 2*HighFree)
  RHS(1:2*HighFree) = -RHS(1:2*HighFree)+BodyForces(1:2*HighFree)
!  print *, RHS(2*NativeHigh+1:2*HighFree)
END SUBROUTINE CalcFEMForce

SUBROUTINE DoFEMStep2(DeltaT)
  DOUBLE PRECISION :: DeltaT
  INTEGER :: I
  DO I=1, NativeHigh! HighFree
    RHS(2*I-1:2*I) = RHS(2*I-1:2*I)/MassVec(I)
!    print *, massvec
  END DO
  Momenta(1:2*NativeHigh) = Momenta(1:2*NativeHigh) + 0.5D0*DeltaT*RHS(1:2*NativeHigh)
END SUBROUTINE DoFEMStep2

      
SUBROUTINE CG(X, N, RHS)
INTEGER :: N, Cnt, I
DOUBLE PRECISION, INTENT(INOUT) :: X(N), RHS(N)
DOUBLE PRECISION :: R(N), Z(N), P(N), Q(N), &
                    Rho, NewRho, Alpha, Beta, Error, Tmp, RealRand
DOUBLE PRECISION, PARAMETER :: MaxRes = 1.D-14
LOGICAL :: First

! This routine calculates the solution of AX=b and 
! stores the result in X
! SparseMult is a routine which multiplies the matrix A by some arbitrary vector X

X = 0.D0
R = RHS
!print *, 'res 0', SQRT(SUM(R*R)), N
Cnt = 0
First = .TRUE.
Error = 10*MaxRes
DO WHILE (Error>MaxRes)
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
  CALL MassMul (P, Q, N)
  Tmp = DOT_PRODUCT(P, Q)
  alpha = NewRho/Tmp
  X = X+Alpha*P
  R = R - alpha*Q
  Rho = NewRho
  Error = SQRT(SUM(R*R))
!  print *, error
END DO
CALL MassMul (X, Q, N)
R = Q - RHS
!print *, 'res', SQRT(SUM(R*R)), Cnt
END SUBROUTINE CG





END MODULE BrMod





