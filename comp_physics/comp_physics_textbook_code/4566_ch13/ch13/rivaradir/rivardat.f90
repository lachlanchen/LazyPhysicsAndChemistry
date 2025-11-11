MODULE rivardat

TYPE NeighbourListType
  TYPE(NeighbourListType), Pointer :: Previous, Next
  INTEGER :: Num
END TYPE NeighbourListType

TYPE VertexType 
  INTEGER :: MyNum
  TYPE (NeighbourListType), Pointer :: NeighbourList
  DOUBLE PRECISION :: Location(2)
END TYPE VertexType

TYPE TriangleType
  TYPE (TriangleType), POINTER :: Previous, Next
  INTEGER :: VertIndex(3)
  DOUBLE PRECISION :: Stiffness(6,6)
  DOUBLE PRECISION :: CBMatrix(3,6)
  DOUBLE PRECISION :: AMatrix(3,3)
  DOUBLE PRECISION :: Stress(3)
  DOUBLE PRECISION :: Area
  INTEGER :: FreeNum   ! Number of non-fixed vertices
END TYPE TriangleType

LOGICAL :: Accurate
DOUBLE PRECISION, ALLOCATABLE :: Displacements(:)

END MODULE rivardat

