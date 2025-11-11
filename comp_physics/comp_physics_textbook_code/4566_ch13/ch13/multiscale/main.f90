PROGRAM Main
USE BrMod
USE MD2D
IMPLICIT NONE

INTEGER :: I, MainNY, IY
DOUBLE PRECISION :: DeltaT=0.002D0, TT, aver
DOUBLE PRECISION, ALLOCATABLE :: MDPos(:,:), FEMPos(:,:), MDBForce(:,:)

CALL Initialise
Call InitMD
CALL DrawFem()
MainNy = OutputNy()
print *, 'mainny', mainny
DO IY=1, 2*MainNY
  FemPart(IY,:) = Locations(NativeHigh-2*MainNy+IY,:) +&
     (/Displacements(2*(NativeHigh-2*MainNy+IY)-1),&
       Displacements(2*(NativeHigh-2*MainNy+IY))/)
END DO
DO IY=1, MainNY
  Displacements(2*(NativeHigh+IY)-1) = Pos(BorderPart(IY),1)-&
                                    Locations(NativeHigh+IY,1)
  Displacements(2*(NativeHigh+IY)) = Pos(BorderPart(IY),2)-&
                                    Locations(NativeHigh+IY,2)
END DO
CALL InitDyn(DeltaT)
CALL InitSim(DeltaT)
TT = 0
DO I=1, 35000
  TT = TT + DeltaT
  DO IY=HighFree+1, VertPtNum
    Displacements(2*IY) = 0.001D0*EXP(-(TT-3.0)**2*2.D0)
  END DO
  IF (MOD(I,10)==0) THEN
    CALL SetNamedBackground('lightblue')
    CALL DrawFem2()
  END IF
  CALL DoFEMStep1(DeltaT)
  CALL DoMDStep1(I, DeltaT)
  DO IY=1, 2*MainNY
    FemPart(IY,:) = Locations(NativeHigh-2*MainNy+IY,:) +&
       (/Displacements(2*(NativeHigh-2*MainNy+IY)-1),&
         Displacements(2*(NativeHigh-2*MainNy+IY))/)
  END DO
  DO IY = 1, MainNY
    Displacements(2*(NativeHigh+IY)-1) = Pos(BorderPart(IY),1)-&
                                      Locations(NativeHigh+IY,1)
    Displacements(2*(NativeHigh+IY)) = Pos(BorderPart(IY),2)-&
                                      Locations(NativeHigh+IY,2)
  END DO
  CALL CalcFEMForce
  CALL CalcMDForce(MDForce)
  DO IY = 1, 2*MainNY
    RHS(2*(NativeHigh-2*MainNy+IY)-1) = RHS(2*(NativeHigh-2*MainNy+IY)-1)+FEMForce(IY,1)
    RHS(2*(NativeHigh-2*MainNy+IY)) = RHS(2*(NativeHigh-2*MainNy+IY))+FEMForce(IY,2)
  END DO
  DO IY = 1, MainNY
    MDForce(BorderPart(IY),1) = MDForce(BorderPart(IY),1)+RHS(2*(NativeHigh+IY)-1)
    MDForce(BorderPart(IY),2) = MDForce(BorderPart(IY),2)+RHS(2*(NativeHigh+IY))
  END DO
  CALL DoFEMStep2(DeltaT)
  CALL DoMDStep2(I, DeltaT)
  IF (TT<5.D0) THEN
    Momenta = Momenta*0.999D0
    Momentum = Momentum*0.999D0
  END IF
END DO
CALL DrawFem()
CALL Endplot()

CONTAINS
  SUBROUTINE SendFEMPosToMD
  END SUBROUTINE SendFEMPosToMD

  SUBROUTINE SendMDPosToFEM
  END SUBROUTINE SendMDPosToFEM

  SUBROUTINE SendFEMFToMD
  END SUBROUTINE SendFEMFToMD

  SUBROUTINE SendMDFToFEM
  END SUBROUTINE SendMDFToFEM

END PROGRAM Main