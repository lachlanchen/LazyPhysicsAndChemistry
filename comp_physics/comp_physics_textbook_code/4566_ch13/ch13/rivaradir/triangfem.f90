MODULE triangfem
! This module contains the handling of the stiffness and other matrices, and the
! displacement vectors
USE rivardat
USE Rivara
IMPLICIT NONE

TYPE(TriangleType), POINTER :: TriangHead
DOUBLE PRECISION, ALLOCATABLE :: BodyForces(:), StressVec(:), R_sigma(:)
DOUBLE PRECISION :: h, Force(2), CMatrix(3,3), CInv(3,3), Nu, E, MaxErr

CONTAINS

SUBROUTINE InitFEM
! General initialisation of parameters
  Force(1) = 0.08D0
  Force(2) = -0.0010D0
  MaxErr = 1.D-3
  ALLOCATE (BodyForces(2*MaxVertNum))
  ALLOCATE (Displacements(2*MaxVertNum))
  ALLOCATE (R_Sigma(3*MaxVertNum))
  ALLOCATE (StressVec(3*MaxVertNum))
  Displacements = 0.D0
  E=1.1D0
  Nu = 0.3D0
  CMatrix(1,:) = (/1.D0,   Nu, 0.D0/)
  CMatrix(2,:) = (/  Nu, 1.D0, 0.D0/)
  CMatrix(3,:) = (/0.D0, 0.D0, 0.5D0*(1.D0-Nu)/)
  CMatrix = CMatrix*E/(1-Nu*Nu)
  CInv(1,:) = (/1.D0,  -Nu, 0.D0/)
  CInv(2,:) = (/ -Nu, 1.D0, 0.D0/)
  CInv(3,:) = (/0.D0, 0.D0, 2.D0*(1.D0+Nu)/)
  CInv = CInv/E
END SUBROUTINE InitFEM


SUBROUTINE InitFEMStep()
  CALL FindTriangles()! Triangles are found from the structures list of grid points
  CALL FindAllStiffnesses() ! Construct stiffness matrices
END SUBROUTINE InitFEMStep

SUBROUTINE DrawFem(Color)
! Draw the configuration
  TYPE(TriangleType), Pointer :: TriangPtr
  INTEGER :: VertIndex, I, Color
  DOUBLE PRECISION :: X1, Y1, X2, Y2, X3, Y3, LL
  TriangPtr => TriangHead
  CALL SetNamedBackground ('lightblue')
  CALL SetFastColor(0)

  CALL FillRectangle(-0.5*hx, -0.5*hy, 0.D0, Ly)
  LL = (INT(0.5D0*Nx)-1)*hx
  CALL FillRectangle(0.0D0, -0.5*hy, LL, 0.D0)
  CALL SetFastColor(Color)
  DO WHILE (ASSOCIATED(TriangPtr)) 
    VertIndex = TriangPtr%VertIndex(1)
    X1 = Vertices(VertIndex)%Location(1)+Displacements(2*VertIndex-1)
    Y1 = Vertices(VertIndex)%Location(2)+Displacements(2*VertIndex)
    VertIndex = TriangPtr%VertIndex(2)
    X2 = Vertices(VertIndex)%Location(1)+Displacements(2*VertIndex-1)
    Y2 = Vertices(VertIndex)%Location(2)+Displacements(2*VertIndex)
    VertIndex = TriangPtr%VertIndex(3)
    X3 = Vertices(VertIndex)%Location(1)+Displacements(2*VertIndex-1)
    Y3 = Vertices(VertIndex)%Location(2)+Displacements(2*VertIndex)
    CALL Draw(X1, Y1, X2, Y2)
    CALL DrawTo (X3, Y3)
    CALL DrawTo (X1, Y1)
    TriangPtr => TriangPtr%Next
  END DO
END SUBROUTINE DrawFem


SUBROUTINE Find_Stiffness(Triangle)
! Find stiffness matrix for the triangle. Also, body forces etcetera are set.
  TYPE (TriangleType) :: Triangle
  TYPE (VertexType) :: CurVertex
  DOUBLE PRECISION :: Locs(3,2), Area2
  DOUBLE PRECISION :: b(3), c(3),  BMatrix(3,6), TempMatrix(3,6), FullB(3,6)
  INTEGER :: I, J
  DO I=1, 3
    CurVertex = Vertices(Triangle%VertIndex(I))
    Locs(I,:) = CurVertex%Location(:)
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
    Triangle%AMatrix(I,I) = Area2/12.D0
    DO J=1, 3
      IF (J/=I) THEN
        Triangle%AMatrix(I,J) = Area2/24.D0
      END IF
    END DO
  END DO
  Triangle%CBMatrix = 0.D0
  J = 0
  DO I=1, 3
    IF (.NOT.FixedPoint(Triangle%VertIndex(I))) THEN
      J = J + 1
      BMatrix(:,2*J-1) = (/ b(I), 0.D0, c(I) /)
      BMatrix(:,2*J) = (/ 0.D0, c(I), b(I) /)
    END IF
    FullB(:,2*I-1) = (/ b(I), 0.D0, c(I) /)
    FullB(:,2*I) = (/ 0.D0, c(I), b(I) /)
  END DO
  BMatrix = BMatrix/Area2
  Triangle%Area = 0.5D0*Area2
  Triangle%FreeNum = J
  Triangle%CBMatrix(1:3,1:2*J) = MATMUL(CMatrix, BMatrix(1:3,1:2*J)) 
  TempMatrix(1:3,1:2*J) = MATMUL(CMatrix, BMatrix(1:3,1:2*J)) 
  Triangle%Stiffness = MATMUL(TRANSPOSE(BMatrix(1:3,1:2*J)), TempMatrix(1:3,1:2*J))
  Triangle%Stiffness = Triangle%Stiffness*0.5D0*Area2 ! Surface integral
  DO I=1, 3
    IF (FixedPoint(Triangle%VertIndex(I))) THEN
      BodyForces(2*Triangle%VertIndex(I)-1:2*Triangle%VertIndex(I)) = 0.D0
    ELSE
      BodyForces(2*Triangle%VertIndex(I)-1) = Area2*Force(1)/3.0D0
      BodyForces(2*Triangle%VertIndex(I))   = Area2*Force(2)/3.0D0
    END IF
  END DO
END SUBROUTINE Find_Stiffness


SUBROUTINE StiffnessMul(OldDisp, NewDisp, N)
! Multiply the old displacement vector by the new one
  INTEGER :: TC, VC, VertIndex, N, J
  DOUBLE PRECISION :: OldNodeVals(6), NewNodeVals(6)
  DOUBLE PRECISION :: OldDisp(N), NewDisp(N)
  TYPE(TriangleType), Pointer :: TriangPtr
  INTEGER :: TriangCount = 0

  NewDisp = 0.D0  
  TriangPtr => TriangHead
  DO WHILE (ASSOCIATED(TriangPtr))
    OldNodeVals =  0.D0
    J=0
    DO VC=1, 3
      VertIndex = TriangPtr%VertIndex(VC)
      IF (.NOT.FixedPoint(VertIndex)) THEN
        J = J + 1
        OldNodeVals(2*J-1:2*J) = OldDisp(2*VertIndex-1:2*VertIndex)
      END IF
    END DO
    NewNodeVals = MATMUL(TriangPtr%Stiffness(1:2*J,1:2*J), OldNodeVals(1:2*J))
    J = 0
    DO VC=1, 3
      VertIndex = TriangPtr%VertIndex(VC)
      IF (.NOT.FixedPoint(VertIndex)) THEN
        J = J + 1
        NewDisp(2*VertIndex-1:2*VertIndex) = NewDisp(2*VertIndex-1:2*VertIndex)+&
                                             NewNodeVals(2*J-1:2*J)
      END IF
    END DO
    TriangPtr => TriangPtr%Next
    TriangCount = TriangCount + 1
  END DO
END SUBROUTINE StiffnessMul


SUBROUTINE StressMul(OldStress, NewStress, N)
! Multiply the stress vector by the AMatrix
  INTEGER :: TC, VC, VertIndex, N, J, VC2
  DOUBLE PRECISION :: OldStressVals(3,3), NewStressVals(3,3)
  DOUBLE PRECISION :: OldStress(N), NewStress(N)
  TYPE(TriangleType), Pointer :: TriangPtr

  TriangPtr => TriangHead
  NewStress = 0.D0
  DO WHILE (ASSOCIATED(TriangPtr))
    OldStressVals =  0.D0
    DO VC=1, 3
      VertIndex = TriangPtr%VertIndex(VC)
      OldStressVals(VC,:) = OldStress(3*VertIndex-2:3*VertIndex)
    END DO
    NewStressVals = 0.D0
    DO J=1, 3
      NewStressVals(:,J) = MATMUL(TriangPtr%AMatrix,OldStressVals(:,J))
    END DO
    DO VC=1, 3
      VertIndex = TriangPtr%VertIndex(VC)
      NewStress(3*VertIndex-2:3*VertIndex) = NewStress(3*VertIndex-2:3*VertIndex)+&
                                             NewStressVals(VC,:)
    END DO
    TriangPtr => TriangPtr%Next
  END DO
END SUBROUTINE StressMul



SUBROUTINE FindTriangles()
! From the list of vertices, that contains information about the neighbours
! of each vertex, construct the list of triangles
  INTEGER :: VertCount, TriangCount, VIC, VI1, VI2, EdgeCount
  TYPE(NeighbourListType), Pointer :: NbPtr1, NbPtr2, FirstPtr
  TYPE(TriangleType), Pointer :: TriangPtr, NextTriang

  TriangCount = 0
  EdgeCount = 0
  print *, 'vertptnum is', vertptnum
  DO VertCount = 1, VertPtNum
    VIC = Vertices(VertCount)%MyNum
    NbPtr1 => Vertices(VertCount)%NeighbourList
    FirstPtr => NbPtr1
    NbPtr2 => NbPtr1%Next
    DO WHILE (ASSOCIATED(NbPtr2) .AND. NbPtr2%Num /= FirstPtr%Num)
      EdgeCount = EdgeCount + 1
      VI1 = NbPtr1%Num; VI2 = NbPtr2%Num
      IF ((VI1>VIC).AND.(VI2>VIC)) THEN
        IF (TriangCount == 0) THEN
          ALLOCATE(TriangHead)
          NULLIFY (TriangHead%Next)
          NULLIFY (TriangHead%Previous)
          TriangPtr => TriangHead
        ELSE
          ALLOCATE (TriangPtr%Next)
          NextTriang => TriangPtr%Next
          NextTriang%Previous => TriangPtr
          TriangPtr => TriangPtr%Next
        END IF
        TriangPtr%VertIndex(1) = VIC
        TriangPtr%VertIndex(2) = VI1
        TriangPtr%VertIndex(3) = VI2
        NULLIFY (TriangPtr%Next)
        TriangCount = TriangCount + 1
      END IF
      NbPtr1 => NbPtr2
      NbPtr2 => NbPtr2%Next
    END DO
    VI1 = NbPtr1%Num; VI2 = FirstPtr%Num
    IF (.NOT.(BoundaryPoint(VertCount))) THEN
      IF ((VI1>VIC).AND.(VI2>VIC)) THEN
        IF (TriangCount == 0) THEN
          ALLOCATE(TriangHead)
          TriangPtr => TriangHead
        ELSE
          ALLOCATE (TriangPtr%Next)
          NextTriang => TriangPtr%Next
          NextTriang%Previous => TriangPtr
          TriangPtr => TriangPtr%Next
        END IF
        TriangPtr%VertIndex(1) = VIC
        TriangPtr%VertIndex(2) = VI1
        TriangPtr%VertIndex(3) = VI2
        NULLIFY (TriangPtr%Next)
        TriangCount = TriangCount + 1
      END IF
    END IF
  END DO
  print *, 'I found ', TriangCount, 'Triangles'
END SUBROUTINE FindTriangles


SUBROUTINE DeAllocTriang
! Triangle no longer needed? Then deallocate
  TYPE(TriangleType), Pointer :: TriangPtr, NextTriang
  INTEGER :: TriangCount = 0

  TriangPtr => TriangHead
  DO WHILE (ASSOCIATED(TriangPtr%Next)) 
    NextTriang => TriangPtr%Next
    DEALLOCATE(TriangPtr)
    TriangPtr => NextTriang
    TriangCount = TriangCount + 1
  END DO
END SUBROUTINE DeAllocTriang

SUBROUTINE ScanTriangles
! Scan through the list of triangles 
  TYPE(TriangleType), Pointer :: TriangPtr
  INTEGER :: TriangCount = 0

  TriangPtr => TriangHead
  DO WHILE (ASSOCIATED(TriangPtr%Next)) 
    TriangPtr => TriangPtr%Next
    TriangCount = TriangCount + 1
  END DO
! NOTE one triangle less scanned; for ALL triangles: DO WHILE (ASSOCIATED(TriangPtr)) !!
  TriangPtr => TriangPtr%Previous
  DO WHILE (ASSOCIATED(TriangPtr)) 
    TriangPtr => TriangPtr%Previous
    TriangCount = TriangCount + 1
  END DO
  print *, 'I scanned through ', TriangCount, ' triangles'
END SUBROUTINE ScanTriangles

SUBROUTINE CheckStress
! Calculate the error based on the stress distribution
  TYPE(TriangleType), Pointer :: TriangPtr
  DOUBLE PRECISION :: StressAver(3), SA2(3), TotStressVec(3), TotStressTRiang(3),&
                      ErrorEst, TempVec(3), TempVec2(3,3), TempVec3(3,3)
  INTEGER :: VC, I, VertIndex, J
  LOGICAL :: Neglect

  TotStressVec = 0.D0
  TotStressTriang = 0.D0
  DO I=1, VertPtNum
    TotStressVec = TotStressVec+StressVec(3*I-2:3*I)
  END DO
  Accurate = .TRUE.
  TriangPtr => TriangHead
  DO WHILE (ASSOCIATED(TriangPtr))
    Neglect = .FALSE.
    TempVec = TriangPtr%Area*MATMUL(CInv, TriangPtr%Stress)
    ErrorEst = DOT_PRODUCT(TempVec, TriangPtr%Stress)
    DO VC = 1, 3
      VertIndex = TriangPtr%VertIndex(VC)
      ErrorEst = ErrorEst-2*DOT_PRODUCT(StressVec(3*VertIndex-2:3*VertIndex), &
                                        TempVec)/3.D0
      TempVec2(VC,:) = MatMul(CInv, StressVec(3*VertIndex-2:3*VertIndex))
    END DO
    DO J=1, 3
      TempVec3(:,J) = MATMUL(TriangPtr%AMatrix, TempVec2(:, J))
    END DO
    DO VC=1, 3
      VertIndex = TriangPtr%VertIndex(VC)
      DO J=1, 3
        ErrorEst = ErrorEst+TempVec3(VC,J)*StressVec(3*VertIndex-3+J)
      END DO
    END DO
!    print *, 'errorest', errorest
    IF ((ErrorEst>MaxErr).AND. (.NOT.Neglect)) THEN
      Accurate = .FALSE.
      CALL RefineTriangle(TriangPtr)
    END IF
    TriangPtr => TriangPtr%Next
  END DO
END SUBROUTINE CheckStress


SUBROUTINE StressCalc(Displacements, N)
! Calculate the stress on the grid
  TYPE(TriangleType), Pointer :: TriangPtr
  INTEGER :: VC, I, N, VertIndex, J
  DOUBLE PRECISION :: NodeVals(6)
  DOUBLE PRECISION :: Displacements(N)

  R_Sigma = 0.D0
  TriangPtr => TriangHead
  DO WHILE (ASSOCIATED(TriangPtr))
    DO VC=1, 3
      VertIndex = TriangPtr%VertIndex(VC)
      NodeVals(2*VC-1:2*VC) = Displacements(2*VertIndex-1:2*VertIndex)
    END DO
    J = TriangPtr%FreeNum
    TriangPtr%Stress = MATMUL(TriangPtr%CBMatrix,NodeVals)
    DO VC=1, 3
      VertIndex = TriangPtr%VertIndex(VC)
      R_Sigma(3*VertIndex-2:3*VertIndex) = R_Sigma(3*VertIndex-2:3*VertIndex)+&
                                TriangPtr%Stress(:)*TriangPtr%Area/3.D0
    END DO
    TriangPtr => TriangPtr%Next
  END DO
END SUBROUTINE StressCalc



SUBROUTINE FindAllStiffnesses
! Routine which loops over all triangles and calls Find_Stiffness 
! which calculates the stiffness matrix for a single triangle
  TYPE(TriangleType), Pointer :: TriangPtr
  TriangPtr => TriangHead
  DO WHILE (ASSOCIATED(TriangPtr)) 
    CALL Find_Stiffness(TriangPtr)
    TriangPtr => TriangPtr%Next
  END DO
END SUBROUTINE FindAllStiffnesses


END MODULE TriangFEM