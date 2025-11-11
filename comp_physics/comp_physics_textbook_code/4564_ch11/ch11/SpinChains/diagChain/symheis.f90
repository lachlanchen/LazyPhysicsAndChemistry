! Programming exercise, section 11.5.2 of the textbook `computational physics',
! second edition, Cambridge 2007 
! Jos Thijssen
! 2006-2007
! Diagonalisation for the spin-half Heisenberg chain, taking symmetry into account
! gfortran does not yet accept allocatables within derived types Hence this
! programs does not work under gfortran
PROGRAM SymHeis
  IMPLICIT NONE
  INTEGER, PARAMETER :: MS=2, L=12
  INTEGER :: Size, MaxNCV, ldv, Nev, MaxNv, I
  DOUBLE PRECISION :: HMat(0:MS-1,0:MS-1,0:MS-1,0:MS-1), &
                      EigVals(MS**L), EigVecs(MS**L,2), Psi0(MS**L), HPsi(MS**L)
  DOUBLE PRECISION, DIMENSION(0:MS-1,0:MS-1) :: SPl, SMi, Sz
  DOUBLE PRECISION :: JX, JZ, Norm, PI ! Couplings of the Heisenberg model
  TYPE SECTORMATS
    DOUBLE PRECISION, ALLOCATABLE :: Mat(:,:)
  END TYPE
  TYPE(SECTORMATS) :: SectorHam(0:MS*L-1, L), SectorBas(0:MS*L-1, L)
  LOGICAL :: Unused(MS**L)

  Size = MS**L
  PI = 4.D0*ATAN(1.D0)
  Jz = 1.d0
  Jx = 1.D0
  CALL FillHMat
  CALL ConstructSymHam(.TRUE.)
CONTAINS

SUBROUTINE ConstructSymHam(Parity)
! Carefully read the text of section 11.5.2 in order to follow this routine.
! The problem is mostly book-keeping.
  LOGICAL, INTENT(IN) :: Parity
  INTEGER :: I, HelpI, LastSpin, Sz, LC, K, KK, BVC, J, UniqNum, &
             count, WorkSize, INFO, MaxKK, HelpJ, SCount, Period, LocalK
  INTEGER, PARAMETER :: MaxUniqNum=MS**L, MaxWorkSize=500
  DOUBLE PRECISION, DIMENSION(MS**L) :: Psi, HelpPsi, Psi2
  INTEGER :: BasVecCnt(0:MS*L-1, L), BasVecNum(0:MS*L-1, L), &
             UniqHead(MaxUniqNum), UniqKK(MaxUniqNum), hELPi2
  DOUBLE PRECISION :: DWork(MaxWorkSize), Diag(MaxWorkSize), Factor

  Unused = .TRUE.
  UniqNum = 0
  count = 0
  DO I = 1, MS**L
    IF (Unused(I)) THEN
      Unused(I) = .FALSE.
      UniqNum = UniqNum+1
      Sz = 0
      HelpI = I
      DO LC=1, L
        Sz = Sz+MOD(HelpI-1, MS)
        HelpI = (HelpI-1)/MS+1
      END DO
      count = count + 1
      IF (UniqNum>MaxUniqNum) THEN
        print *, 'too many uniques; allocate more memory', UniqNum
        STOP
      END IF
      UniqHead(UniqNum) = I
      HelpI = I
      DO KK=1, L-1 ! Calculate period of basis state I (=UniqHead(UniqNum))
        LastSpin=(HelpI-1)/MS**(L-1)
        HelpI = MOD(((HelpI-1)*MS + LastSpin),MS**L)+1 ! Translation
        IF (HelpI == I) THEN
          Period = L/KK ! Is period of state in k-space
          BasVecCnt(Sz, 1:(L/2+1):Period) = BasVecCnt(Sz, 1:(L/2+1):Period) + 2
          BasVecCnt(Sz, 1) = BasVecCnt(Sz, 1)-1
          UniqKK(UniqNum) = KK
          EXIT
        END IF
        IF (Unused(HelpI)) THEN 
          Unused(HelpI) = .FALSE.
        END IF
        IF ((KK==L-1))THEN      ! i.e. if state not periodic
          BasVecCnt(Sz, 1:L/2+1) = BasVecCnt(Sz, 1:L/2+1) + 2
          BasVecCnt(Sz, 1) = BasVecCnt(Sz, 1)-1
          UniqKK(UniqNum) = L
        END IF
      END DO
    END IF
  END DO
  BasVecNum = BasVecCnt
! Now allocate memory for basis states and Hamiltonians in sectors
  DO Sz = 0, L*(MS-1)
    DO K=1, L/2+1
      ALLOCATE (SectorHam(Sz,K).Mat(BasVecCnt(Sz,K), BasVecCnt(Sz,K)))
      ALLOCATE (SectorBas(Sz,K).Mat(MS**L, BasVecCnt(Sz,K)))
      print *, Sz, K, BasVecCnt(Sz,K)
    END DO
  END DO
  BasVecCnt = 0
  DO J=1, UniqNum
    Psi = 0.D0
    I = UniqHead(J)
    Psi(I) = 1.D0
! `Measure' spin of this sector
    Sz = 0
    DO LC=1, L
      Sz = Sz+MOD(I-1, MS)
      I = (I-1)/MS+1
    END DO
    I = UniqHead(J)
    HelpI = I
    Period = L/UniqKK(J)
    BasVecCnt(Sz,1:(L/2+1):Period) = BasVecCnt(Sz, 1:(L/2+1):Period)+1
    print *, 'b'
    DO K=1, L/2+1, Period
      IF (BasVecNum(Sz,K)>0) THEN
        BVC = BasVecCnt(Sz,K)
        IF (Parity) THEN
          SectorBas(Sz,K).Mat(:,BVC) = Psi
         IF (K/=1) THEN
           BasVecCnt(Sz,K) = BasVecCnt(Sz,K)+1
         END IF
        ELSE
          SectorBas(Sz,K).Mat(:,BVC) = 0.D0 ! Psi
        END IF
      END IF
    END DO
    DO KK=1, UniqKK(J)-1! We started already with HelpPsi = Psi
      LastSpin = (HelpI-1)/MS**(L-1)
      HelpI = MOD(((HelpI-1)*MS + LastSpin),MS**L)+1
      DO K=1, L/2+1, Period ! VERT
        IF (K==1) THEN
          BVC = BasVecCnt(Sz, K) 
          Factor = COS(2*KK*(K-1)*PI/DBLE(L))
          SectorBas(Sz,K).Mat(HelpI,BVC) = SectorBas(Sz,K).Mat(HelpI,BVC) +&
                            Factor
        ELSE
          BVC = BasVecCnt(Sz,K)-1
          Factor = COS(2*KK*(K-1)*PI/DBLE(L))
          SectorBas(Sz,K).Mat(HelpI,BVC) = SectorBas(Sz,K).Mat(HelpI,BVC) +&
                            Factor
          Factor = SIN(2*KK*(K-1)*PI/DBLE(L)) 
          BVC = BasVecCnt(Sz,K)
           SectorBas(Sz,K).Mat(HelpI,BVC) = SectorBas(Sz,K).Mat(HelpI,BVC) +&
                            Factor
         END IF
      END DO
    END DO
  END DO
  IF (.FALSE.) THEN
  DO Sz = 0, L*(MS-1)
    DO K=1, L/2+1
      DO BVC = 1, BasVecCnt(Sz,k)
        IF (Parity) THEN
          Unused = .TRUE.
          Psi = 0.D0
          DO I = 1, MS**L
              HelpI = I
              CALL Reflect(HelpI)
              IF ((Sz==1).AND.(K==2).AND.(ABS(SectorBas(Sz,K).Mat(I, BVC))>1.D-12)) THEN
                print *, 'I'
                CALL ShowSpins(I)
                print *, 'HelpI'
                CALL ShowSpins(HelpI)
              END IF
              Psi(HelpI) = SectorBas(Sz,K).Mat(I, BVC)
          END DO
          SectorBas(Sz,K).Mat(:, BVC) = (Psi+SectorBas(Sz,K).Mat(:, BVC))
        END IF
      END DO
      CALL Gram_Schmidt(SectorBas(Sz,K).Mat, BVC, L*(MS-1))
      BasVecCnt(Sz, K) = BVC
    END DO
  END DO
  END IF
  Count= 0
  BasVecNum =BasVecCnt
  DO Sz = 0, MS*L-1
    DO K=1, L/2+1
      BVC = BasVecNum(Sz, K)
      IF (BVC>0) THEN
        HelpI =0
        DO I=1, BVC
          Psi = SectorBas(Sz,K).Mat(:,I)
          IF (SQRT(SUM(Psi*Psi))>1.d-14) THEN
            HelpI = HelpI + 1
            PSI = Psi/SQRT(SUM(Psi*Psi))
            CALL SuperPsi(Size, Psi, HelpPsi)
            SectorHam(Sz, K).Mat(HelpI,HelpI) = DOT_PRODUCT(Psi, HelpPsi)
            HelpJ = 0
            Psi2 = Psi
            DO J=1, I-1
              Psi = SectorBas(Sz,K).Mat(:,J)
              IF (SQRT(SUM(Psi*Psi))>1.d-14) THEN
                HelpJ = HelpJ + 1
                IF (HelpJ>HelpI) THEN
                  print *, 'problem in routine ConstructSymHam'
                  stop
                END IF
                IF (DOT_PRODUCT(Psi, Psi2)>1.D-10) THEN
                  print *, 'dependent basis', Sz, K, HelpI, HelpJ
                  STOP
                END IF
                PSI = Psi/SQRT(SUM(Psi*Psi))
                SectorHam(Sz, K).Mat(HelpI,HelpJ) = DOT_PRODUCT(Psi, HelpPsi)
                SectorHam(Sz, K).Mat(HelpJ,HelpI) = &
                                            SectorHam(Sz, K).Mat(HelpI,HelpJ)
              END IF
            END DO
          END IF
        END DO
        print *, 'Sz, K, Size', Sz, K, HelpI
        IF (HelpI>0) THEN 
        CALL DSYEV('n', 'u', HelpI, SectorHam(Sz, K).Mat, BVC, Diag, DWork, MaxWorkSize, INFO) 
          Count = Count+HelpI
          print '(6F10.5)', Diag(1:HelpI)
        END IF
      END IF
   END DO
  END DO
  DO Sz = 0, L*(MS-1)
    DO K=1, L/2+1
      DEALLOCATE (SectorHam(Sz,K).Mat)
      DEALLOCATE (SectorBas(Sz,K).Mat)
    END DO
  END DO
  print *, 'nr of states:', Count
END SUBROUTINE ConstructSymHam


SUBROUTINE ShowSpins(I)
INTEGER :: HelpI, LC, I
  HELPI = I
  DO LC=1, L
    print *, MOD(HelpI-1+MS, MS)
    HelpI = (HelpI-1)/MS+1
  END DO
END SUBROUTINE ShowSpins


SUBROUTINE Reflect(I)
INTEGER :: HelpI, LC, I
  HELPI = 1
  DO LC=1, L
    HelpI = HelpI+MOD(I-1+MS, MS)*MS**(L-LC)
    I = (I-1)/MS+1
  END DO
  I = HelpI
END SUBROUTINE Reflect

  SUBROUTINE Gram_Schmidt(Vectors, Number, Dimen)
  IMPLICIT NONE
  INTEGER :: Number, Dimen
  DOUBLE PRECISION, INTENT(INOUT) :: Vectors(Dimen, Number)
  INTEGER Iorb1, Iorb2, NewCount
  DOUBLE PRECISION :: IP
  LOGICAL :: OnOf(Number)

  NewCount = 1
  OnOf = .TRUE.
  DO Iorb1 = 1, Number
    DO Iorb2 = 1, Iorb1-1
      IF (OnOf(Iorb2)) THEN
        IP = DOT_PRODUCT(Vectors(:, Iorb2),Vectors(:, Iorb1))
        Vectors(:, Iorb1) = Vectors(:, IOrb1)-IP*Vectors(:, Iorb2)
      END IF
    END DO
    IP = DOT_PRODUCT(Vectors(:, Iorb1),Vectors(:, Iorb1))
    IF (IP>1.D-10) THEN
      IP = 1/SQRT(IP)
      Vectors(:,NewCount) = Vectors(:,NewCount)*IP
      NewCount = NewCount + 1
    ELSE
      OnOf(Iorb1) = .FALSE.
    END IF
    Number = NewCount
  END DO

  END SUBROUTINE Gram_Schmidt

SUBROUTINE SuperPsi(Size, Psi, HPsi)
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
                S  = N1+N2*MSI*MS**2+s1*MSI + s2 *MSI*MS+1
                SP = N1+N2*MSI*MS**2+s1p*MSI+ s2p*MSI*MS+1
!		print *, 'm', S, SP
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
  DO N1=0, MS**(L-2)-1
    DO s1 = 0, MS-1
      DO s2 = 0, MS-1
        DO s1p = 0, MS-1
          DO s2p = 0, MS-1
            S  = N1*MS+s1*MS**(L-1)+s2+1
            SP = N1*MS+s1p*MS**(L-1)+s2p+1
!	    print *, s, sp
            HPsi(SP)=HPsi(SP)+ HMat(s1p, s2p, s1, s2)*Psi(S)
          END DO
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE SuperPsi
  
SUBROUTINE FillHMat
  INTEGER :: s1, s2, s1p, s2p
  DOUBLE PRECISION :: FirstPl, FirstMi, FirstZ, PreFac, Mz
  Sz = 0.D0; SMi=0.D0; SPl=0.D0
  DO I=0, Ms-1
    Sz(I,I) = (Ms-1-2*I)*0.5D0
  END DO
  DO I=0, Ms-2
    Mz = DBLE(-Ms+1+2*I)*0.5D0
    PreFac = SQRT((MS-1)*(MS+1)*0.25D0-Mz*(Mz+1))
    SPl(I+1, I) = PreFac
    SMi(I,I+1)  = PreFac
  END DO
  print *, 'Sz'
  print '(6F10.4)', Sz
  print *, 'S+'
  print '(6F10.4)', SPl
  print *, 'S-'
  print '(6F10.4)', SMi
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
  DO s1=0, MS-1
    DO s1p=0, MS-1
      DO s2=0, MS-1
        DO s2p=0, MS-1
          IF (ABS(HMat(s1,s1p,s2,s2p)-hmat(s2,s2p,s1,s1p))>1.D-8) THEN
            print *, 'wat een ramp'
            stop
          end if
        end do
      end do
    end do
  end do
  print '(6F10.4)', HMat
END SUBROUTINE FillHMat



END PROGRAM SymHeis
    
    