 PROGRAM FEM
 IMPLICIT NONE
 ! Dynamics of an (in)elastic beam suspended on a wall, bending under the 
 ! action of gravity.
 ! A finite element calculation is carried out for a triangular mesh.
 ! We construct the finite element metrices B from the three corner points.
 ! Per node we have its offset position and we calculate its displacement
 ! All displacements are stored in a system vector, r
 ! The stiffness matrix is a 3x3 matrix. It originates from the linear 
 ! relations between the six degrees of freedom (the u's and the v's at the corners) 
 ! which are not independent.
 ! The dynamics is solved using the Verlet algorithm
 ! Program written by Jos Thijssen, 2004-2006
 
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
 CALL Dynamics()
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
 Force(2) = -0.0D0
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
 !Multiplication of the old displacement by the stiffness matrix to
 !obtain the new displacements
 
 INTEGER :: TC, VC, VertIndex, N, J
 
 DOUBLE PRECISION :: OldNodeVals(6), NewNodeVals(6)
 DOUBLE PRECISION :: OldDisp(N), NewDisp(N)
 
 NewDisp = 0.D0
 DO TC=1, TriangNum
   OldNodeVals =0.D0
   J=0
   DO VC=1, 3
     VertIndex = Triangles(TC)%VertIndex(VC)
     J = J + 1
     OldNodeVals(2*J-1:2*J) = OldDisp(2*VertIndex-1:2*VertIndex)
   END DO
   NewNodeVals = MATMUL(Triangles(TC)%Stiffness(1:2*J,1:2*J), OldNodeVals(1:2*J))
   J = 0
   DO VC=1, 3
     VertIndex = Triangles(TC)%VertIndex(VC)
     J = J + 1
     NewDisp(2*VertIndex-1:2*VertIndex) = NewDisp(2*VertIndex-1:2*VertIndex)+&
                                            NewNodeVals(2*J-1:2*J)
   END DO
 END DO
 END SUBROUTINE StiffnessMul
 
 
 
 SUBROUTINE MassMul(OldDisp, NewDisp, N)
 ! This routine multiplies the old displacement by the 
 ! mass matrix. It is called by CG (Conjugate gradients) in order to
 ! solve the implicit equation occurring when solven for the new displacements
 ! in Verlet's algorithm
 
 INTEGER :: TC, VC, VertIndex, N, J
 
 DOUBLE PRECISION :: OldNodeVals(6), NewNodeVals(6)
 DOUBLE PRECISION, INTENT(IN) :: OldDisp(N)
 DOUBLE PRECISION, INTENT(OUT) :: NewDisp(N)
 
 NewDisp = 0.D0
 DO TC=1, TriangNum
   OldNodeVals =0.D0
   J=0
   DO VC=1, 3
     VertIndex = Triangles(TC)%VertIndex(VC)
     J = J + 1
     OldNodeVals(2*J-1:2*J) = OldDisp(2*VertIndex-1:2*VertIndex)
   END DO
   NewNodeVals(1:2*J) = MATMUL(Triangles(TC)%AMatrix(1:2*J,1:2*J), OldNodeVals(1:2*J))
   J = 0
   DO VC=1, 3
     VertIndex = Triangles(TC)%VertIndex(VC)
     J = J + 1
     NewDisp(2*VertIndex-1:2*VertIndex) = NewDisp(2*VertIndex-1:2*VertIndex)+&
                                            NewNodeVals(2*J-1:2*J)
   END DO
 END DO
 END SUBROUTINE MassMul
       
 
 
 SUBROUTINE Dynamics
 ! This routine carries out the Verlet algorithm for the dynamical
 ! equations for the grid points.
 
 DOUBLE PRECISION :: DeltaT, X, h, TT, Fac2, AvPos, AvLoc, AvY, &
                     eps11, eps22, stress, force, TotWidth, nnu
 INTEGER :: StepCount, IY, I
 
 DeltaT = TimeStep
 h = 2.D0**(1.D0/6.D0)
 TotWidth = (Nx-1)*h
 DO I=1, 2*HighFree
   Displacements(2*I-1) = Locations(I,1)/TotWidth*4.323D0*1.D-4
   Displacements(2*I) = 0.D0
 END DO
 Momenta(1:2*HighFree) = 0.D0
 Fac2 = Fac2*DeltaT**2*0.5D0
 Force = GlobForce
 CALL CalcForces(Force)
 DO StepCount = 1, 120000
   TT = StepCount*DeltaT
   Momenta(1:2*HighFree) = Momenta(1:2*HighFree)+0.5D0*DeltaT*RHS(1:2*HighFree)
   Displacements(1:2*HighFree) =  Displacements(1:2*HighFree)+DeltaT*Momenta(1:2*HighFree)
   CALL CalcForces(Force)
   Momenta(1:2*HighFree) = Momenta(1:2*HighFree)+0.5D0*DeltaT*RHS(1:2*HighFree)
   Momenta = Momenta*(1-DeltaT*0.003D0)! Artificial damping
   IF ((MOD(StepCount, 600).EQ.0).AND.(StepCount>90)) THEN
     CALL DrawFem()
     AvPos = 0.D0
     AvLoc = 0.D0
     AvY = 0.D0
     DO IY=1, NY
       AvY = AvY+Displacements(2*(HighFree-(8+IY)*NY))
       AvLoc = AvLoc + Locations(HighFree-NY+IY,1)
       AvPos = AvPos + Displacements(2*(HighFree-NY+IY)-1)
     END DO
     AvPos = AvPos/Ny; AvY = AvY/Ny; AvLoc = AvLoc/Ny
     stress = Ny*force/((NY-1)*h)
     eps11 = avpos/((NX-1)*h)
     eps22 = 2*avy/((NY-1)*h)
     nnu = -eps22/(eps11-eps22)
   END IF
 END DO
 END SUBROUTINE Dynamics
 
 
 SUBROUTINE CalcForces(force)
 ! Calculate the force vectors occurring in Verlet's algorithm
   INTEGER :: I, IY
   DOUBLE PRECISION :: Force
 
   CALL StiffnessMul(Displacements(1:2*HighFree), RHS(1:2*HighFree), 2*HighFree)
   DO IY = 1, NY
     RHS(2*(HighFree-NY+IY)-1) = RHS(2*(HighFree-NY+IY)-1)
   END DO
   DO I=1, HighFree
     RHS(2*I-1) = -RHS(2*I-1)/MassVec(I)
     RHS(2*I) = -RHS(2*I)/MassVec(I)-Force
   END DO
 END SUBROUTINE CalcForces
 
       
 SUBROUTINE CG(X, N, RHS)
 INTEGER :: N, Cnt, I
 DOUBLE PRECISION, INTENT(INOUT) :: X(N), RHS(N)
 DOUBLE PRECISION :: R(N), Z(N), P(N), Q(N), &
                     Rho, NewRho, Alpha, Beta, Error, Tmp
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
 
 
 
 
 
 END PROGRAM FEM
 
 
 
 
 
