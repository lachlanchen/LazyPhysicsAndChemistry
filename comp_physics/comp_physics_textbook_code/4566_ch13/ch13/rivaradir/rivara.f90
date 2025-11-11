MODULE Rivara
USE rivardat
IMPLICIT NONE
! In this unit, a triangular grid is set up according to the structure
! and method suggested by Rivara (1984)

INTEGER :: MaxVertNum, VertPtNum
INTEGER :: Nx = 60
INTEGER :: Ny ! is determined below
DOUBLE PRECISION :: Lx = 15.0, Ly=2.0, hx, hy


TYPE(VertexType), ALLOCATABLE :: Vertices(:)
LOGICAL, ALLOCATABLE :: BoundaryPoint(:), FixedPoint(:)


CONTAINS

SUBROUTINE InitGrid
! Set up the grid and initialise
! Also, the neighbour lists are set up in this routine. These lists contain
! the information on all neighbours of each vertex.
! Lots of book-keeping.
  INTEGER :: IX, IY, TriangCount, LowCount, HighCount, Count
  DOUBLE PRECISION :: Xll, Yll, Xur, Yur
  TYPE (NeighbourListType), Pointer :: NbPtr, First

  MaxVertNum = 2000 ! Might be read from input!!!!!!!!! Indeed, not dynamically
                    ! allocated, I'm afraid ;-)
  hx = Lx/(Nx-1);
  Ny = NINT(Nx*Ly/Lx)
  hy = Ly/(Ny-1)
  IF (NY<3) THEN
    NY = 3
    NX = NINT(Lx/Ly*Ny)
    hx = Lx/(Nx-1);
    hy = Ly/(Ny-1)
  END IF  
  VertPtNum = NX*NY
  ALLOCATE (Vertices(MaxVertNum))
  ALLOCATE (BoundaryPoint(MaxVertNum))
  ALLOCATE (FixedPoint(MaxVertNum))
  FixedPoint = .FALSE.
  Count = 0
  DO IX = 2, NX-1
    DO IY = 2, NY-1
      Count = VertexIndex(Ix, IY)
      BoundaryPoint(Count) = .FALSE.
      Vertices(Count)%Location(:) = (/ (IX-1)*hx, (IY-1)*hy /)
      Vertices(Count)%MyNum = Count
      CALL MakeFirstNeighbour(Vertices(Count), IX+1,IY)
      NbPtr => Vertices(Count)%NeighbourList
      CALL MakeNextNeighbour (NbPtr, IX, IY+1)
      CALL MakeNextNeighbour (NbPtr, IX-1, IY+1)
      CALL MakeNextNeighbour (NbPtr, IX-1, IY)
      CALL MakeNextNeighbour (NbPtr, IX, IY-1)
      CALL MakeNextNeighbour (NbPtr, IX+1, IY-1)
      First => Vertices(Count)%NeighbourList
      NbPtr%Next =>  First
      First%Previous => NbPtr
    END DO
  END DO
  DO IX = 2, NX-1
    Count = VertexIndex(IX, 1)
    BoundaryPoint(Count) = .TRUE.
    Vertices(Count)%Location(:) = (/ (IX-1)*hx, 0.D0 /)
    Vertices(Count)%MyNum = Count
    IF (IX<=NX/2) FixedPoint(Count) = .TRUE.
    CALL MakeFirstNeighbour(Vertices(Count), IX+1,1)
    NbPtr => Vertices(Count)%NeighbourList
    CALL MakeNextNeighbour(NbPtr, IX, 2)
    CALL MakeNextNeighbour(NbPtr, IX-1, 2)
    CALL MakeNextNeighbour(NbPtr, IX-1, 1)
    Count = VertexIndex(IX, NY)
    BoundaryPoint(Count) = .TRUE.
    Vertices(Count)%Location(:) = (/ (IX-1)*hx, (NY-1)*hy /)
    Vertices(Count)%MyNum = Count
    CALL MakeFirstNeighbour(Vertices(Count), IX-1,NY)
    NbPtr => Vertices(Count)%NeighbourList
    CALL MakeNextNeighbour(NbPtr, IX, NY-1)
    CALL MakeNextNeighbour(NbPtr, IX+1, NY-1)
    CALL MakeNextNeighbour(NbPtr, IX+1, NY)
  END DO
  DO IY = 2, NY-1
    Count = VertexIndex(1, IY)
    BoundaryPoint(Count) = .TRUE.
    FixedPoint(Count) = .TRUE.
    Vertices(Count)%Location(:) = (/ 0.D0, (IY-1)*hy /)
    Vertices(Count)%MyNum = Count
    CALL MakeFirstNeighbour(Vertices(Count), 1,IY-1)
    NbPtr => Vertices(Count)%NeighbourList
    CALL MakeNextNeighbour(NbPtr, 2, IY-1)
    CALL MakeNextNeighbour(NbPtr, 2, IY)
    CALL MakeNextNeighbour(NbPtr, 1, IY+1)
    Count = VertexIndex(NX, IY)
    BoundaryPoint(Count) = .TRUE.
    Vertices(Count)%Location(:) = (/ (NX-1)*hx, (IY-1)*hy /)
    Vertices(Count)%MyNum = Count
    CALL MakeFirstNeighbour(Vertices(Count), NX,IY+1)
    NbPtr => Vertices(Count)%NeighbourList
    CALL MakeNextNeighbour(NbPtr, NX-1, IY+1)
    CALL MakeNextNeighbour(NbPtr, NX-1, IY)
    CALL MakeNextNeighbour(NbPtr, NX, IY-1)
  END DO
  Count = VertexIndex(1, 1)
  BoundaryPoint(Count) = .TRUE.
  Vertices(Count)%Location(:) = (/ 0.D0, 0.D0 /)
  Vertices(Count)%MyNum = Count
  FixedPoint(Count) = .TRUE.
  CALL MakeFirstNeighbour(Vertices(Count), 2, 1)
  NbPtr => Vertices(Count)%NeighbourList
  CALL MakeNextNeighbour(NbPtr, 1, 2)
  Count = VertexIndex(NX, NY)
  BoundaryPoint(Count) = .TRUE.
  Vertices(Count)%Location(:) = (/ (NX-1)*hx, (NY-1)*hy /)
  Vertices(Count)%MyNum = Count
  CALL MakeFirstNeighbour(Vertices(Count), NX-1, NY)
  NbPtr => Vertices(Count)%NeighbourList
  CALL MakeNextNeighbour(NbPtr, NX, NY-1)
  Count = VertexIndex(NX, 1)
  BoundaryPoint(Count) = .TRUE.
  Vertices(Count)%Location(:) = (/ (NX-1)*hx, 0.D0 /)
  Vertices(Count)%MyNum = Count
  CALL MakeFirstNeighbour(Vertices(Count), NX, 2)
  NbPtr => Vertices(Count)%NeighbourList
  CALL MakeNextNeighbour(NbPtr, NX-1, 2)
  CALL MakeNextNeighbour(NbPtr, NX-1, 1)
  Count = VertexIndex(1, NY)
  BoundaryPoint(Count) = .TRUE.
  Vertices(Count)%Location(:) = (/ 0.D0, (NY-1)*hy /)
  Vertices(Count)%MyNum = Count
  FixedPoint(Count) = .TRUE.
  CALL MakeFirstNeighbour(Vertices(Count), 1, NY-1)
  NbPtr => Vertices(Count)%NeighbourList
  CALL MakeNextNeighbour(NbPtr, 2, NY-1)
  CALL MakeNextNeighbour(NbPtr, 2, NY)
  CALL InitPlot('lightblue', 1400,400, 'out.ps', 1)
  Xll = -hx
  Yll = -0.8D0*Lx
  Xur = 1.6D0*NX*hx 
  Yur = 0.8D0*Lx
  CALL Framing (Xll, Yll, Xur, Yur);
  CALL PutStopButton()
!  CALL DrawFirst()
END SUBROUTINE InitGrid



SUBROUTINE MakeFirstNeighbour (Vertex, IX, IY)
! Creates the first element of the neighbourList of Vertex
  TYPE (VertexType) :: Vertex
  INTEGER, INTENT(IN) :: IX, IY

  ALLOCATE(Vertex%NeighbourList)
  NULLIFY(Vertex%NeighbourList%Previous)
  NULLIFY(Vertex%NeighbourList%Next)
  Vertex%NeighbourList%Num = VertexIndex(IX, IY)
END SUBROUTINE MakeFirstNeighbour



SUBROUTINE MakeNextNeighbour (NbPtr, IX, IY)
! Creates the following vertex
  TYPE (NeighbourListType), Pointer :: NbPtr, CurPtr
  INTEGER, INTENT(IN) :: IX, IY

  CurPtr => NbPtr
  ALLOCATE(NbPtr%Next)
  NbPtr => NbPtr%Next
  NbPtr%Previous => CurPtr
  NbPtr%Num = VertexIndex(IX, IY)
  NULLIFY (NbPtr%Next)
END SUBROUTINE MakeNextNeighbour


SUBROUTINE DrawFirst
! Draw initial configuration
  INTEGER :: VertIndex, r, g, b, VC, J, I
  TYPE (NeighbourListType), POINTER :: NbPtr, First
  DOUBLE PRECISION :: X1, Y1, X2, Y2, X3, Y3, Polygon(6), StressAver

  CALL SetNamedBackground('lightblue')
  DO I=1, VertPtNum
    X1 = Vertices(I)%Location(1); Y1 = Vertices(I)%Location(2)
    VertIndex = Vertices(I)%NeighbourList%Num;
    X2 = Vertices(VertIndex)%Location(1); Y2 = Vertices(VertIndex)%Location(2)
    CALL Draw(X1, Y1, X2, Y2)
    NbPtr => Vertices(I)%NeighbourList%Next
    First => NbPtr
    DO
      VertIndex = NbPtr%Num;
      X3 = Vertices(VertIndex)%Location(1); Y3 = Vertices(VertIndex)%Location(2)
      CALL Draw(X1, Y1, X3, Y3)
      NbPtr => NbPtr%Next
      IF (.NOT.ASSOCIATED(NbPtr) .OR. NbPtr%Num == First%Num) Exit
    END DO
  END DO
END SUBROUTINE DrawFirst


INTEGER FUNCTION VertexIndex(IX, IY)
! This function returns the index which the vertex at position 
! IX, IY in the long list called Vertices. Simply based on 
! the fact that this list was set up in a structured way.

  INTEGER :: IX, IY

  VertexIndex = ((IX-1)*NY+IY)
END FUNCTION VertexIndex

SUBROUTINE RefineTriangle(TriangPtr)
! This is an important routine which carries out the 
! refinement of a triangle where the error estimator is too big
  TYPE(TriangleType), Pointer :: TriangPtr
  INTEGER :: VI1, VI2

  IF (.NOT.Refined(TriangPtr)) THEN
    CALL FindLongestEdge (TriangPtr, VI1, VI2)
    CALL AddNeighbour(VI1, VI2)
  END IF
END SUBROUTINE RefineTriangle

LOGICAL FUNCTION Refined(TriangPtr)
! Tells you whether the triangle has been refined or not.
  TYPE (TriangleType), Pointer :: TriangPtr
  INTEGER :: VC1, VC2, VI1, VI2

  Refined = .FALSE.
  DO VC1 = 1, 3
    VI1 = TriangPtr%VertIndex(VC1)
    DO VC2 = VC1+1, 3
      VI2 = TriangPtr%VertIndex(VC2)
      IF (NoNeighbours(VI1, VI2)) THEN
        Refined = .TRUE.
      END IF
    END DO
  END DO
END FUNCTION Refined

LOGICAL FUNCTION NoNeighbours(VI1, VI2)
! Are VI1 and VI2 neighbours?
  INTEGER :: VI1, VI2
  TYPE(NeighbourListType), Pointer :: NbPtr, First

  NoNeighbours = .TRUE.
  NbPtr => Vertices(VI1)%NeighbourList
  First => NbPtr
  DO    
    IF (.NOT.ASSOCIATED(NbPtr) .OR. NbPtr%Num == VI2) EXIT
    NbPtr => NbPtr%Next
    IF (ASSOCIATED(NbPtr) .AND. NbPtr%Num == First%Num) EXIT
  END DO
  IF (ASSOCIATED(NbPtr) .AND. NbPtr%Num == VI2) THEN
    NoNeighbours = .FALSE.
  END IF
  NbPtr => Vertices(VI2)%NeighbourList
  First => NbPtr
  DO
    IF (.NOT.ASSOCIATED(NbPtr) .OR. NbPtr%Num == VI1) EXIT
    NbPtr => NbPtr%Next
    IF (ASSOCIATED(NbPtr) .AND. NbPtr%Num == First%Num) EXIT
  END DO
  IF (ASSOCIATED(NbPtr) .AND. NbPtr%Num == VI1) THEN
    IF (NoNeighbours) print *, 'something is wrong in NoNeighbours'
    NoNeighbours = .FALSE.
  END IF
END FUNCTION NoNeighbours
  


SUBROUTINE FindLongestEdge (TriangPtr, VI1, VI2)
! Finds the longest edge of the triangle TriangPtr and returns the
! end points of that edge as VI1 and VI2
  TYPE(TriangleType), Pointer :: TriangPtr
  INTEGER :: VI1, VI2, VC, VC1, VC2, TVI1, TVI2
  DOUBLE PRECISION :: LongDist

  LongDist = 0.D0
  DO VC=1, 3
    VC1 = MOD(VC,3)+1
    VC2 = MOD(VC-2+3,3)+1
    TVI1 = TriangPtr%VertIndex(VC1)
    TVI2 = TriangPtr%VertIndex(VC2)
    IF (Distance(TVI1, TVI2)>LongDist) THEN
      VI1 = TVI1
      VI2 = TVI2
      LongDist = Distance(TVI1, TVI2)
    END IF
  END DO
END SUBROUTINE FindLongestEdge

DOUBLE PRECISION FUNCTION Distance(VI1, VI2)
  INTEGER :: VI1, VI2
  DOUBLE PRECISION :: D(2)
  
  D = Vertices(VI1)%Location-Vertices(VI2)%Location
  Distance = D(1)*D(1)+D(2)*D(2)
END FUNCTION Distance


SUBROUTINE CheckVertex(VI)
! Consistency check, for debugging purposes
  INTEGER :: VI
  TYPE(NeighbourListType), Pointer :: NeighPtr, First
  NeighPtr => Vertices(VI)%NeighbourList
  First => NeighPtr
  print *, 'I am ', VI
  DO 
    IF (.NOT.Associated(NeighPtr)) Exit
    print *, NeighPtr%Num
    NeighPtr => NeighPtr%Next
    IF (ASSOCIATED(NeighPtr) .AND. First%Num == NeighPtr%Num) Exit
  END DO
  IF (ASSOCIATED(NeighPtr)) THEN
    print *, NeighPtr%Num
  END IF
  IF (ASSOCIATED(First%Previous)) THEN
    NeighPtr => First%Previous
    print *, 'last', NeighPtr%Num
  END IF
  PAUSE
END SUBROUTINE CheckVertex


SUBROUTINE FindNeighbour(VIC, VIN, NeighPtr)
! From a vertex with index VIC, find the element of its neighbourlist corresponding 
! to its neighbour VIN
  INTEGER :: VIC, VIN
  TYPE(NeighbourListType), Pointer :: NeighPtr, First

  NeighPtr => Vertices(VIC)%NeighbourList
  First => NeighPtr
  DO 
    IF (.NOT.Associated(NeighPtr) .OR. NeighPtr%Num == VIN) Exit
    NeighPtr => NeighPtr%Next
    IF (ASSOCIATED(NeighPtr) .AND. First%Num == NeighPtr%Num) Exit
  END DO  
END SUBROUTINE FindNeighbour


SUBROUTINE FindVI3_4(VI1, VI2, VI3, VI4)
! For two vertices VI1 and VI2, that define an edge where a new point is
! to be created, find the other two points to which this new point should be
! connected. So VI3 and VI4 are both neighbours of VI1 and VI2.
  INTEGER :: VI1, VI2, VI3, VI4
  TYPE(NeighbourListType), Pointer :: NeighPtr, PrevNbPtr

  CALL FindNeighbour (VI1, VI2, NeighPtr)
  IF (ASSOCIATED(NeighPtr%Previous)) THEN
    PrevNbPtr => NeighPtr%Previous
    VI4 = PrevNbPtr%Num
  ELSE
    VI4 = 0
  END IF
  IF (ASSOCIATED(NeighPtr%Next)) THEN
    PrevNbPtr => NeighPtr%Next
    VI3 = PrevNbPtr%Num
  ELSE
    VI3 = 0
  END IF
END SUBROUTINE FindVI3_4


SUBROUTINE CheckRefine(VI1, VI2, VI3, VI4, RefRes, RV1, RV2)
! Check whether the edge on which a new point is to be created,
! is indeed the longest of its two adjacent triangles

  INTEGER :: VI1, VI2, VI3, VI4, RefRes, RV1, RV2
  DOUBLE PRECISION :: D12, D13, D14, D23, D24

  RefRes = 0

  D12 = Distance(VI1, VI2)
  IF (VI3 /=0) THEN
    D13 = Distance(VI1, VI3)
    D23 = Distance(VI2, VI3)
  END IF
  IF (VI4 /=0) THEN
    D14 = Distance(VI1, VI4)
    D24 = Distance(VI2, VI4)
  END IF
  IF (VI3*VI4 /= 0) THEN
    IF (D12<D13 .OR. D12<D23) THEN
      RefRes = 3 ! The 'wrong' triangle has vertices vi1, vi2 and vi3
      RV1 = VI3
      IF (D13>D23) THEN
        RV2 = VI1
      ELSE
        RV2 = VI2
      END IF
    ELSE IF (D12<D14 .OR. D12<D24) THEN
           RefRes = 4
           RV1 = VI4
           IF (D14>D24) THEN
             RV2 = VI1
           ELSE
             RV2 = VI2
           END IF
         END IF
    END IF
END SUBROUTINE CheckRefine


RECURSIVE SUBROUTINE AddNeighbour(VI1, VI2)
  INTEGER, INTENT(IN) :: VI1, VI2
  INTEGER :: VI3, VI4, NewNeighbour, T1, T2, T3, T4, VNum, RefRes, RV1, RV2
  TYPE(NeighbourListType), Pointer :: NeighPtr, PrevNbPtr, First, Temp

  CALL FindVI3_4(VI1, VI2, VI3, VI4)
  CALL CheckRefine(VI1, VI2, VI3, VI4, RefRes, RV1, RV2)
  IF (RefRes /= 0) THEN
    CALL AddNeighbour (RV1, RV2)
  END IF
  CALL FindVI3_4(VI1, VI2, VI3, VI4)
  VertPtNum = VertPtNum + 1
  Vertices(VertPtNum)%Location = 0.5D0*(Vertices(VI1)%Location &
                                       +Vertices(VI2)%Location)
  Displacements(2*VertPtNum-1:2*VertPtNum) = 0.5D0*(Displacements(2*VI1-1:2*VI1)&
                                       +Displacements(2*VI2-1:2*VI2))
  NewNeighbour = VertPtNum
  Vertices(NewNeighbour)%MyNum = NewNeighbour
! Redirect neighbour from VI2 to newneighbour
  CALL FindNeighbour (VI2, VI1, NeighPtr)
  IF (NeighPtr%Num /=VI1) print *, 'alarm 2'
  NeighPtr%Num = NewNeighbour
!Redirect neighbour from VI1 to newneighbour  
  CALL FindNeighbour (VI1, VI2, NeighPtr)
  IF (NeighPtr%Num /= VI2) print *, 'alarm 1'
  NeighPtr%Num = NewNeighbour
!Create neighbours of new vertex
! Check whether new point is on the boundary; fill its value for BoundaryPoint
  BoundaryPoint(NewNeighbour) = (BoundaryPoint(VI1).AND.BoundaryPoint(VI2).AND.VI3*VI4==0)
! Check whether new point is fixed; fill its value for FixedPoint
  FixedPoint(NewNeighbour) = (BoundaryPoint(NewNeighbour).AND.FixedPoint(VI1).AND.FixedPoint(VI2))
  IF (BoundaryPoint(NewNeighbour)) THEN
    VNum = 3
    IF (VI4 == 0) THEN
      T1 = VI2; T2 = VI3; T3 = VI1 
    ELSE 
      IF (VI3 == 0) THEN
        T1 = VI1; T2 = VI4; T3 = VI2
      ELSE
        print *, 'problem, only 2 neighbours', VI3, VI4, BoundaryPoint(NewNeighbour)
      END IF
    ENDIF
  ELSE
    VNum = 4
    T1 = VI1; T2 = VI4; T3 = VI2; T4 = VI3 !T's number the four points anti-clockwise
  END IF
! Now we construct the new links containing the new-born point
  ALLOCATE(Vertices(NewNeighbour)%NeighbourList)
  NeighPtr => Vertices(NewNeighbour)%NeighbourList
  NeighPtr%Num = T1
  NULLIFY (NeighPtr%Previous)
  ALLOCATE (NeighPtr%Next)
  PrevNbPtr => NeighPtr
  NeighPtr => NeighPtr%Next
  NeighPtr%Previous => PrevNbPtr
  PrevNbPtr%Next => NeighPtr
  NeighPtr%Num = T2
  ALLOCATE (NeighPtr%Next)
  PrevNbPtr => NeighPtr
  NeighPtr => NeighPtr%Next
  NeighPtr%Previous => PrevNbPtr
  PrevNbPtr%Next => NeighPtr
  NeighPtr%Num = T3
  NULLIFY(NeighPtr%Next)
  IF (VNum == 4) THEN ! generic case, four points: end points of longest edge
    ALLOCATE (NeighPtr%Next) ! and the two neighbours of both endpoints
    PrevNbPtr => NeighPtr
    NeighPtr => NeighPtr%Next
    NeighPtr%Previous => PrevNbPtr
    PrevNbPtr%Next => NeighPtr
    NeighPtr%Num = T4
    NULLIFY(NeighPtr%Next)
  END IF
  IF (.NOT.BoundaryPoint(NewNeighbour)) THEN
    First => Vertices(NewNeighbour)%NeighbourList
    First%Previous => NeighPtr
    NeighPtr%Next => First
  END IF
! Redirect neighbour from VI3 to newneighbour
  IF (VI3/=0) THEN
    CALL FindNeighbour(VI3, VI1, NeighPtr)
    IF (NeighPtr%Num /=VI1) then ! consistency check
      print *, 'alarm 3', vi1, vi2, vi3, vi4, newneighbour; stop
    end if
    ALLOCATE(Temp)
    Temp%Previous => NeighPtr
    Temp%Next => NeighPtr%Next
    Temp%Num = NewNeighbour
    PrevNbPtr => NeighPtr
    NeighPtr  => NeighPtr%Next
    PrevNbPtr%Next => Temp
    NeighPtr%Previous => Temp
  END IF
  IF (VI4/=0) THEN
    CALL FindNeighbour(VI4, VI2, NeighPtr)
    IF (NeighPtr%Num /=VI2) then ! consistency check
     print *, 'alarm 4', NeighPtr%Num, VI2; stop
    end if
    ALLOCATE(Temp)
    Temp%Previous => NeighPtr
    Temp%Next => NeighPtr%Next
    Temp%Num = NewNeighbour
    PrevNbPtr => NeighPtr
    NeighPtr  => NeighPtr%Next
    PrevNbPtr%Next => Temp
    NeighPtr%Previous => Temp
  END IF
END SUBROUTINE AddNeighbour



END MODULE RIVARA