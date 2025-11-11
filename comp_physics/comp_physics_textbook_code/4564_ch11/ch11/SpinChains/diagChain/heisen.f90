! Programming exercise, section 11.5.1 of the textbook `computational physics',
! second edition, Cambridge 2007 
! Jos Thijssen
! 2006-2007
! Straightforward lanczos diagonalisation for the spin-half Heisenberg chain.
! Note: the present version works for ifort, not for gfortran. This is however
! easily fixed by explitly including the lanczos code into the CONTAINS section
! and removing the subroutine argument MultPsi, replacing it in the lanczos
! code by SuperPsi
PROGRAM HEISEN ! With Periodic BC!!!!
  IMPLICIT NONE
  INTEGER, PARAMETER :: MS=2, L=12
  INTEGER :: Size, MaxNCV, ldv, Nev, MaxNv, I ! Spin-1/2
  DOUBLE PRECISION :: HMat(0:MS-1,0:MS-1,0:MS-1,0:MS-1), &
                      EigVals(MS**L), EigVecs(MS**L,2), Psi0(MS**L), HPsi(MS**L)
  DOUBLE PRECISION, DIMENSION(0:MS-1,0:MS-1) :: Sx, Sy, Sz, SPl, SMi
  DOUBLE PRECISION :: JX, JZ, Norm ! Couplings of the Heisenberg model

  Size = MS**L
  Jz = 1.d0
  Jx = 1.D0
  CALL FillHMat
  maxncv=25
  ldv=Size
  Nev = 2
  CALL Random_Number(Psi0)
  CALL Lanczos(Psi0, Size, 50, SuperPsi)
 
CONTAINS

SUBROUTINE SuperPsi(Psi, HPsi, Size)
! Multiplication of H with Psi
  INTEGER :: Size
  DOUBLE PRECISION :: Psi(Size), HPsi(Size)
  INTEGER :: s1, s2, s1p, s2p, I, S, SP, N1, N2,MSI, MSLI, L
  L = NINT(LOG(DBLE(Size))/Log(DBLE(MS)))
  HPsi = 0.D0
  MSI=1
  MSLI = MS**(L-2)
  DO i=0, L-2
    DO N1=0, MSI-1
      DO N2=0, MSLI-1
        DO s1 = 0, MS-1
          DO s2 = 0, MS-1
            DO s1p = 0, MS-1
              DO s2p = 0, MS-1
                S  = N1+N2*MSI*MS**2+s1*MSI + s2*MSI*MS+1
                SP = N1+N2*MSI*MS**2+s1p*MSI+ s2p*MSI*MS+1
                HPsi(SP) = HPsi(SP)+HMat(s1p, s2p, s1, s2)*Psi(S)
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    MSI = MSI*MS ! MSI = MS**I
    MSLI = MSLI/MS ! MSLI = MS**(L-I-2)
  END DO
!  return ! Remove for PBC!!!!
  DO N1=0, MS**(L-2)-1
    DO s1 = 0, MS-1
      DO s2 = 0, MS-1
        DO s1p = 0, MS-1
          DO s2p = 0, MS-1
            S  = N1*MS+s1*MS**(L-1)+s2+1
            SP = N1*MS+s1p*MS**(L-1)+s2p+1
            HPsi(SP)=HPsi(SP)+ HMat(s1p, s2p, s1, s2)*Psi(S)
          END DO
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE SuperPsi



SUBROUTINE FillHMat
!  Fills the coupling hamiltonian for two neighbouring sites.
  INTEGER :: s1, s2, s1p, s2p
  DOUBLE PRECISION :: FirstPl, FirstMi, FirstZ, PreFac, Mz

  Sz = 0.D0; SMi=0.D0; SPl=0.D0
  DO I=0, Ms-1
    Sz(I,I) = DBLE(Ms-1-2*I)*0.5D0
  END DO
  DO I=0, Ms-2
    Mz = DBLE(-Ms+1+2*I)*0.5D0
    PreFac = SQRT((MS-1)*(MS+1)*0.25D0-Mz*(Mz+1))
    SPl(I, I+1) = PreFac
    SMi(I+1,I)  = PreFac
  END DO
  DO s1=0, MS-1
    DO s1p=0, MS-1
      FirstPl = SPl(s1p,s1)
      FirstMi = SMi(s1p,s1)
      FirstZ = Sz(s1p,s1)
      DO s2=0, MS-1
        DO s2p=0, MS-1
          HMat(s1p,s2p,s1,s2) = 0.5*Jx*(FirstPl*SMi(s2p,s2)+FirstMi*SPl(s2p,s2))+&
                               Jz*(FirstZ*Sz(s2p,s2))
        END DO
      END DO
    END DO
  END DO
! Check the symmetry
  DO s1=0, MS-1
    DO s1p=0, MS-1
      DO s2=0, MS-1
        DO s2p=0, MS-1
          IF (ABS(HMat(s1,s1p,s2,s2p)-hmat(s2,s2p,s1,s1p))>1.D-8) THEN
            print *, 'Disaster: hamiltonian non-hermitian!'
            stop
          end if
        end do
      end do
    end do
  end do
  print '(6F10.4)', HMat
END SUBROUTINE FillHMat



END PROGRAM Heisen
    
    