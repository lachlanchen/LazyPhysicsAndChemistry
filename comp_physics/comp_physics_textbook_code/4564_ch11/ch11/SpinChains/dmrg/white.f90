! Programming exercise, section 11.7 of the textbook `computational physics',
! second edition, Cambridge 2007 
! Jos Thijssen
! 2006-2007
! Straightforward DMRG program for the spin-half Heisenberg chain.

PROGRAM WHITE
  IMPLICIT NONE
  INTEGER, PARAMETER :: M=40, MS=2, L=12 ! S is the `working size' of 
  ! the algorithm. MS =2 corresponds to spin S=1/2, L is the initial chain length.
  INTEGER, PARAMETER :: MaxDensMatSize= MS*M
  INTEGER :: Size, I, LeftHalfSize 
  DOUBLE PRECISION :: EV
  DOUBLE PRECISION :: HMat(0:MS-1,0:MS-1,0:MS-1,0:MS-1), &
                      EigVals(MS**L), EigVecs(MS**L,2), Psi0(MS**L), HPsi(M,M)
  DOUBLE PRECISION, DIMENSION(MS,MS) :: ESMinus, ESPlus, ESZ
  DOUBLE PRECISION, DIMENSION(M*MS,M*MS) :: HL, SMinusM, SPlusM, SZM
  DOUBLE PRECISION :: JX, JZ, Norm ! Couplings of the Heisenberg model
  DOUBLE PRECISION, DIMENSION (M, MS, MS, M) :: Psi
  DOUBLE PRECISION, DIMENSION(MaxDensMatSize,MaxDensMatSize) :: DensMat
  
  Size = MS**L
  Jz = 0.25d0
  Jx = 0.25D0
  CALL FillHMat
  CALL FillSpinMats
  CALL Random_Number(Psi0)
  CALL Random_Number(Psi0)
  Norm = SQRT(SUM(Psi0*Psi0))
  Psi0 = Psi0/Norm
  CALL Lanczos(Psi0, Size, 75)
  LeftHalfSize = MS**(L/2)
  print *, 'LeftHalfSize =', LeftHalfSize
  CALL PreCalcDensMat(Psi0,LeftHalfSize) !Preliminary calculations
  CALL PreCalcMats(LeftHalfSize)
  CALL ReduceDensMat(LeftHalfSize) ! Use the eigenvectormatrix V to reduce the !
                                   ! size of the density matrix
  CALL Random_Number(Psi)	! First calculation for full chain ground state
  CALL BigLanczos(Psi, 100, EV)
  CALL DMRG(Psi)
CONTAINS


SUBROUTINE DMRG(Psi)
! Main machinery. 
  DOUBLE PRECISION :: Psi(M, MS, MS, M), EV, PrevEV
  INTEGER :: CurrentSize, Iter, MaxIter

  MaxIter = 10000
  PrevEV = 1.D0
  CurrentSize = M*MS
  DO Iter = 1, MaxIter
    CALL CalcDensMat(Psi) 	! Trace out E
    CALL CalcMats		! Calculate relevant matrices: H on S
    CALL ReduceDensMat(CurrentSize) ! Reduce size of rho_S
    CALL WaveFunctionTransform(Psi, CurrentSize) ! Transform psi according to V
    CALL BigLanczos(Psi, 100, EV) ! Find current ground state
    print *, 'iter', Iter, 0.5*(EV-PrevEV), 0.5*(EV-PrevEV)+log(2.d0)-0.25D0
    WRITE (8,*) Iter, 0.5*(EV-PrevEV)
    PrevEV = EV
  END DO
END SUBROUTINE DMRG


SUBROUTINE CalcDensMat(Psi)
! Straightforward calculation of the density matrix, rho = |psi><psi|
  INTEGER :: sl1, sl2, sr, ml1,ml2,mr, I1,I2
  DOUBLE PRECISION :: Psi(M, MS, MS, M)
  Psi = Psi/SQRT(SUM(Psi*Psi))
  DensMat = 0.D0
  DO sl1=1, MS
    DO sl2=1, MS
      DO sr = 1, MS
        DO ml1=1, M
          DO ml2=1, M
            DO mr=1, M
              I1 = (sl1-1)*M+ml1
              I2 = (sl2-1)*M+ml2
              DensMat(I1, I2) = DensMat(I1,I2) + &
                 Psi(ml1, sl1, sr, mr)*Psi(ml2, sl2, sr, mr)
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE CalcDensMat
 

SUBROUTINE CalcMats
! Calculate matrix representations of various operators within current basis
  INTEGER :: m1, m2, s1, s2, I, J,worksize, INFO
  DOUBLE PRECISION :: Psi(M, MS), NewPsi(M, MS), NewHL(M*MS,M*MS)
!, diag(m*Ms), dwork(3*M*MS)

  worksize= 3*M*MS
  NewHL = 0.D0
  DO s1 = 1, MS
    DO m1 = 1, M
      Psi = 0.D0
      I = (s1-1)*M+m1
      Psi(m1,s1) = 1.D0
      CALL SystemHPsi(Psi, NewPsi) 
      DO s2=1, MS
        DO m2=1, M
          J=(s2-1)*M+m2
          NewHL(J,I) = NewPsi(m2,s2)
        END DO
      END DO
    END DO
  END DO
  HL=NewHL
  SPlusM = 0.D0
  SMinusM =0.D0
  SZM = 0.D0
  DO m1 = 1, M
    I = m1
    SZM(I,I) = 1.D0
    SZM(I+M,I+M) = -1.D0
    SPlusM(I,I+M) = 2.D0
    SMinusM(I+M, I) = 2.D0
  END DO   
END SUBROUTINE CalcMats    


SUBROUTINE SysBorderInt(SPsi, NewSPsi, SpinMat1, SpinMatM2, JCoupl)
! Calculates the part of H Psi which involves the border of S.
  INTEGER :: ml, mr, mlp, mrp, sl, sr, mp, sp
  DOUBLE PRECISION, DIMENSION(M,MS) :: SPsi, NewSPsi, IntSPsi1
  DOUBLE PRECISION, DIMENSION(MS,MS) :: SpinMat1
  DOUBLE PRECISION, DIMENSION(M*MS,M*MS) :: SpinMatM2
  DOUBLE PRECISION :: JCoupl
  IntSPsi1 = 0.D0
  DO sl=1, MS ! First multiplication of 
    DO ml = 1, M
      DO sp = 1, MS
        IntSPsi1(ml,sl) = IntSPsi1(ml,sl) + &
                             SpinMat1(sl,sp)*SPsi(ml,sp)
      END DO
    END DO
  END DO
  DO sl=1, MS
    DO ml = 1, M
      DO mp = 1, M
        NewSPsi(ml,sl) = NewSPsi(ml,sl) + &
                        JCoupl*SpinMatM2(ml,mp)*IntSPsi1(mp,sl)
      END DO
    END DO
  END DO
END SUBROUTINE SysBorderInt



SUBROUTINE SystemHPsi(SPsi, NewSPsi)
! H Psi calculated on S
  INTEGER :: ml, mlp, sl, sp, mp
  DOUBLE PRECISION, DIMENSION(M,MS) :: SPsi, NewSPsi, IntSPsi
  
  NewSPsi = 0.D0
  DO sl = 1, MS ! Multiply with system and environment Hamiltonians
    DO ml = 1, M
      DO mp = 1, M
        NewSPsi(ml,sl) = NewSPsi(ml,sl) + &
                                   HL(ml,mp)*SPsi(mp,sl)
      END DO  
    END DO
  END DO
  !return
  CALL SysBorderInt(SPsi, NewSPsi, ESMinus, SPlusM, 0.5*Jx)
  CALL SysBorderInt(SPsi, NewSPsi, ESPlus, SMinusM, 0.5*Jx)
  CALL SysBorderInt(SPsi, NewSPsi, ESZ, SZM, Jz)
END SUBROUTINE SystemHPsi








SUBROUTINE BorderInt(SPsi, NewSPsi, SpinMat1, SpinMat2, SpinMatM1, SpinMatM2, JCoupl)
! Part of Homailtonian-wavefunction multiplication, which takes the 
! interactions at the boundary between S and the new central spins into account.
  INTEGER :: ml, mr, mlp, mrp, sl, sr, mp, sp
  DOUBLE PRECISION, DIMENSION(M,MS,MS,M) :: SPsi, NewSPsi, IntSPsi1, IntSPsi2
  DOUBLE PRECISION, DIMENSION(MS,MS) :: SpinMat1, SpinMat2
  DOUBLE PRECISION, DIMENSION(M*MS,M*MS) :: SpinMatM1, SpinMatM2
  DOUBLE PRECISION :: JCoupl
  IntSPsi1 = 0.D0
  IntSPsi2 = 0.D0
  DO sl=1, MS 
    DO sr=1, MS
      DO mr = 1, M
        DO ml = 1, M
          DO sp = 1, MS
            IntSPsi1(ml,sl,sr,mr) = IntSPsi1(ml,sl,sr,mr) + &
                             SpinMat1(sl,sp)*SPsi(ml,sp,sr,mr)
            IntSPsi2(ml,sl,sr,mr) = IntSPsi2(ml,sl,sr,mr) + &
                             SPinMat2(sr,sp)*SPsi(ml,sl,sp,mr)
          END DO
        END DO
      END DO
    END DO
  END DO
  DO sl=1, MS
    DO sr=1, MS
      DO mr = 1, M
        DO ml = 1, M
          DO mp = 1, M
            NewSPsi(ml,sl,sr,mr) = NewSPsi(ml,sl,sr,mr) + &
                        JCoupl*SpinMatM2(ml,mp)*IntSPsi1(mp,sl,sr,mr)
            NewSPsi(ml,sl,sr,mr) = NewSPsi(ml,sl,sr,mr) + &
                        JCoupl*SpinMatM1(mr,mp)*IntSPsi2(ml,sl,sr,mp)
          END DO
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE BorderInt



SUBROUTINE CentralInt(SPsi, NewSPsi, SpinMat1, SpinMat2, JCoupl)
! Part of Homailtonian-wavefunction multiplication, which takes the 
! interactions for the central spins into account.
  INTEGER :: ml, mr, mlp, mrp, sl, sr, mp, sp
  DOUBLE PRECISION, DIMENSION(M,MS,MS,M) :: SPsi, NewSPsi, IntSPsi
  DOUBLE PRECISION, DIMENSION(MS,MS) :: SpinMat1, SpinMat2
  DOUBLE PRECISION :: JCoupl
  IntSPsi = 0.D0
  DO sl=1, MS ! First multiplication of 
    DO sr=1, MS
      DO mr = 1, M
        DO ml = 1, M
          DO sp = 1, MS
            IntSPsi(ml,sl,sr,mr) = IntSPsi(ml,sl,sr,mr) + &
                             SpinMat1(sl,sp)*SPsi(ml,sp,sr,mr)
          END DO
        END DO
      END DO
    END DO
  END DO
  DO sl=1, MS
    DO sr=1, MS
      DO mr = 1, M
        DO ml = 1, M
          DO sp = 1, MS
            NewSPsi(ml,sl,sr,mr) = NewSPsi(ml,sl,sr,mr) + &
                        JCoupl*SpinMat2(sr,sp)*IntSPsi(ml,sl,sp,mr)
          END DO
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE CentralInt



SUBROUTINE SuperPsi2(SPsi, NewSPsi)
! Hamiltonian-wavefunction multiplication for DMRG. Basis is split into 
! the states labelled ml and mr of the left and right block,
! and the central spins sl and sr.
  INTEGER :: ml, mr, mlp, mrp, sl, sr, sp, mp
  DOUBLE PRECISION, DIMENSION(M,MS,MS,M) :: SPsi, NewSPsi, IntSPsi
  
  NewSPsi = 0.D0
  DO sl = 1, MS ! Multiply with system and environment Hamiltonians
    DO sr = 1, MS
      DO mr = 1, M
        DO ml = 1, M
          DO mp = 1, M
            NewSPsi(ml,sl,sr,mr) = NewSPsi(ml,sl,sr,mr) + &
                                   HL(ml,mp)*SPsi(mp,sl,sr,mr)
            NewSPsi(ml,sl,sr,mr) = NewSPsi(ml,sl,sr,mr) + &
                HL(mr,mp)*SPsi(ml,sl,sr,mp) !Note: HL=HR
          END DO  
        END DO
      END DO
    END DO
  END DO
  CALL BorderInt(SPsi, NewSPsi, ESMinus, ESPlus, SMinusM, SPlusM, 0.5*Jx)
  CALL BorderInt(SPsi, NewSPsi, ESPlus, ESMinus, SPlusM, SMinusM, 0.5*Jx)
  CALL BorderInt(SPsi, NewSPsi, ESZ, ESZ, SZM, SZM, Jz)
  CALL CentralInt(SPsi, NewSPsi, ESMinus, ESPlus, 0.5*Jx)
  CALL CentralInt(SPsi, NewSPsi, ESPlus, ESMinus, 0.5*Jx)
  CALL CentralInt(SPsi, NewSPsi, ESZ, ESZ, Jz)
END SUBROUTINE SuperPsi2

  
  

SUBROUTINE SuperPsi(Size, Psi, HPsi)
! Multiplication of the wavefunction by the Hamiltonian, without reference to
! the DMRG blocks. Basis is the spin-up and down basis on each site.
  INTEGER :: Size, L
  DOUBLE PRECISION :: Psi(Size), HPsi(Size)
  INTEGER :: s1, s2, s1p, s2p, I, S, SP, N1, N2, MSI, MSLI
  
  L = NINT(log(DBLE(Size))/log(DBLE(MS)))
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
! The following part contains the couplings between the last and 
! first spins. Only present for periodic chains
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
! Fills the interaction matrix between neighbouring spins
  INTEGER :: s1, s2, s1p, s2p
  DOUBLE PRECISION :: FirstMi, FirstPl, FirstZ
  DOUBLE PREcISION, DIMENSION (0:MS-1, 0:MS-1) :: SPl, Smi, Sz
  Spl(0,:) = (/0.d0, 2.d0 /)
  Spl(1,:) = (/0.d0, 0.d0 /)
  Smi(0,:) = (/0.d0, 0.d0 /)
  Smi(1,:) = (/2.d0, 0.d0 /)
  Sz(0,:) = (/1.d0, 0.d0 /)
  Sz(1,:) = (/0.d0, -1.d0 /)
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
!  print '(4F10.4)', HMat
END SUBROUTINE FillHMat





SUBROUTINE FillSpinMats
! Matrices for S+, S- and Sz
  INTEGER :: s, sp

  ESPlus(1,:)  = (/0.d0, 2.d0/)
  ESPlus(2,:)  = (/0.D0, 0.D0/)
  ESMinus(1,:) = (/0.d0, 0.d0/)
  ESMinus(2,:) = (/2.D0, 0.D0/)
  ESZ(1,:)     = (/1.d0, 0.d0/)
  ESZ(2,:)     = (/0.D0, -1.D0/)
  
END SUBROUTINE FillSpinMats



DOUBLE PRECISION FUNCTION Big_Dot_Product(Psi1, Psi2)
  INTEGER :: ml, mr, mlp, mrp, sl, sr, slp, srp
  DOUBLE PRECISION, DIMENSION(M,MS,MS,M), INTENT(IN)::Psi1, Psi2

  Big_Dot_Product = SUM(Psi1*Psi2)
END FUNCTION Big_Dot_Product


SUBROUTINE BigLanczos(phi0,mm, EV)
  INTEGER, INTENT(IN)::mm
  DOUBLE PRECISION, DIMENSION(M,MS,MS,M), INTENT(INOUT)::phi0
  INTEGER::i, INFO, tm
  DOUBLE PRECISION, DIMENSION(M,MS,MS,M) ::bphi,Aphi
  DOUBLE PRECISION, DIMENSION(0:mm-1,M,MS,MS,M) ::phi
  DOUBLE PRECISION, DIMENSION(0:Mm-1) :: a,b, d, e
  DOUBLE PRECISION, DIMENSION(1:2*Mm-2) :: WrkSpace
  DOUBLE PRECISION, DIMENSION(Mm,Mm) :: Z
  DOUBLE PRECISION :: EV, preva0
  tm =mm
  preva0 = 1.d8
  phi0 = phi0/SQRT(SUM(Phi0*Phi0))
  phi(0,:,:,:,:)=phi0
  CALL SuperPsi2(phi(0,:,:,:,:), APhi)
  a(0)=big_dot_product(phi(0,:,:,:,:),Aphi)
  bphi=Aphi-a(0)*phi(0,:,:,:,:)
  b(0)=SQRT(SUM(bphi*bphi))
  phi(1,:,:,:,:)=bphi/b(0)
  b(0) = big_dot_product(phi(1,:,:,:,:),APhi)
  DO i=1,mm-1
    CALL SuperPsi2(phi(i,:,:,:,:),APhi)
    a(i)=big_dot_product(phi(i,:,:,:,:),Aphi)
    bphi=Aphi-(a(i)*phi(i,:,:,:,:))-(b(i-1)*phi(i-1,:,:,:,:))
    b(i)=SQRT(SUM(bphi*bphi))
    IF (i<mm-1) phi(i+1,:,:,:,:)=bphi/b(i)
    b(i) = big_dot_product(phi(i+1,:,:,:,:), APhi)
    d = a
    e = b
    CALL DSTEV('V',i, d , e(0:i) , Z, mm,  wrkspace , info )
    IF (abs((preva0-d(0))/d(0))<1.d-11) THEN
      tm = i
      a(0)=d(0)
      EXIT
    ELSE
      preva0 = d(0)
    END IF
  ENDDO
!  print *, tm, mm
!  print *, a(0)
  Phi0 = 0.D0
  DO I=0, tm-1
    phi0 = Phi0 + Z(I+1,1)*Phi(I,:,:,:,:)
  END DO
  EV = a(0)
END SUBROUTINE BigLanczos


SUBROUTINE PreCalcMats(Size)
!Prelim calculations for HL (Hamiltonian of left half) and
!local intercation matrices
!Size should be the basis size of the LEFT HALF!!!!!
  INTEGER :: I, L, MSI, J, JP, Size
  DOUBLE PRECISION :: Psi(Size), HPsi(Size)

  L = NINT(LOG(DBLE(Size))/LOG(DBLE(MS)))
  
  MSI = MS**L
  SPlusM = 0.D0
  SMinusM = 0.D0
  SZM = 0.D0
  DO I=1, MSI
    Psi = 0.D0
    Psi(I) = 1.d0
    CALL SuperPsi(MSI,Psi, HPsi)
    HL(1:MSI,I) = HPsi(1:MSI)
  END DO
  MSI = MS**(L-1)
  DO I=1, MSI
    J = I
    JP = J+MSI
    SPlusM(J,JP) = 2.D0
    SPlusM(J,J)  = 0.D0
    SPlusM(JP,JP)= 0.D0
    SZM(J,J)  =  1.D0
    SZM(JP,JP)= -1.D0
    SZM(J,JP) =  0.D0
    SZM(JP,J) =  0.D0
    SMinusM(JP,J) =2.D0
    SMinusM(J,J)  =0.D0
    SMinusM(JP,JP)=0.D0
  END DO      
END SUBROUTINE PreCalcMats    
    



SUBROUTINE PreCalcDensMat(Psi0, MSI)
! Preliminary calculation of density matrix
  INTEGER :: MSI, I1, I2, J
  DOUBLE PRECISION :: Psi0(MSI*MSI)
  Psi0 = Psi0/SQRT(SUM(Psi0*Psi0))
  DensMat = 0.D0
  DO I1=0, MSI-1
    DO I2=0, MSI-1
      DO J=0, MSI-1
        DensMat(I1+1,I2+1) = DensMat(I1+1, I2+1) + Psi0(I1+1+J*MSI)*Psi0(I2+1+J*MSI)
      END DO
    END DO
  END DO
END SUBROUTINE PreCalcDensMat




SUBROUTINE ReduceDensMat(CurSize)
! This is a cxrucial routine, in which the density matrix is 
! used to reduce the Hamiltonian on S and the central spin 
! matrices to a smaller size.
  INTEGER :: CurSize, WorkSize, INFO
  DOUBLE PRECISION :: DWork(3*CurSize-1), Diag(CurSize)
  DOUBLE PRECISION :: TempMat(M,M)
  WorkSize=(3*CurSize-1)
  IF (CurSize<M) THEN
    print *,'problem in reducedensmat: size of densmat too small'
    STOP
  END IF
  CALL DSYEV('v', 'u', CurSize, DensMat, MaxDensMatSize, Diag, DWork, WorkSize, INFO)
  print *, 'denstrace', CurSize, SUM(diag), SUM(diag(CurSize-M+1:CurSize))
  IF (INFO/=0) THEN
    PRINT *, 'error diagonalising density matrix'
    STOP
  END IF
  TempMat = MatMul(TRANSPOSE(DensMat(1:CurSize,CurSize-M+1:CurSize)), &
       MATMUL(HL(1:CurSize,1:CurSize),DensMat(1:CurSize,CurSize-M+1:CurSize)))
  HL(1:M,1:M) = TempMat
  TempMat = MatMul(TRANSPOSE(DensMat(1:CurSize,CurSize-M+1:CurSize)), &
                   MATMUL(SPlusM(1:CurSize,1:CurSize),DensMat(1:CurSize,CurSize-M+1:CurSize)))
  SPlusM(1:M,1:M) = TempMat
  TempMat = MatMul(TRANSPOSE(DensMat(1:CurSize,CurSize-M+1:CurSize)), &
                   MATMUL(SMinusM(1:CurSize,1:CurSize),DensMat(1:CurSize,CurSize-M+1:CurSize)))
  SMinusM(1:M,1:M) = TempMat
  TempMat = MatMul(TRANSPOSE(DensMat(1:CurSize,CurSize-M+1:CurSize)), &
                   MATMUL(SZM(1:CurSize,1:CurSize),DensMat(1:CurSize,CurSize-M+1:CurSize)))
  SZM(1:M,1:M) = TempMat
END SUBROUTINE ReduceDensMat 




SUBROUTINE WaveFunctionTransform(Psi, CurSize)
! The ground state wave function is transformed to the new basis
! defined by the eigenvectors with the largest eigenvalues of rho
  DOUBLE PRECISION, DIMENSION(M, MS, MS, M) :: Psi, NewPsi
  DOUBLE PRECISION, DIMENSION(M, MS, M) :: IntPsi
  INTEGER :: sl, sr, ml, mlp, mr, mrp, I, CurSize
  IntPsi =0.D0 ! Intermediate
  NewPsi = 0.D0 ! New wavefunction
  DO sl = 1, MS
    DO sr = 1, MS
      DO ml = 1, M
        DO mlp =1, M
          DO mr =1, M
            I = (sl-1)*M + ml
            IntPsi(mlp, sr, mr) = IntPsi(mlp, sr, mr) + DensMat(I,CurSize-M+mlp)*Psi(ml, sl, sr, mr)
          END DO
        END DO
      END DO
    END DO
  END DO
  DO sl = 1, MS
    DO sr = 1, MS
      DO ml = 1, M
        DO mr = 1, M
          DO mrp = 1, M
            I = (sr-1)*M + mr
            NewPsi(ml, sl, sr, mr) = NewPsi(ml, sl, sr, mr) + DensMat(I,CurSize-M+mrp)*IntPsi(ml,sl,mrp)
          END DO
        END DO
      END DO
    END DO
  END DO
  Psi = NewPsi
END SUBROUTINE WaveFunctionTransform


SUBROUTINE Lanczos(phi0,n,m)
! Straightforward Lanczos for the chain without reference to 
! S and E
  INTEGER, INTENT(IN)::n,m
  DOUBLE PRECISION, DIMENSION(n), INTENT(INOUT)::phi0
  INTEGER::i, INFO
  DOUBLE PRECISION, DIMENSION(1:n) ::bphi,Aphi
  DOUBLE PRECISION, DIMENSION(0:m-1,1:n) ::phi
  DOUBLE PRECISION, DIMENSION(0:M-1) :: a,b
  DOUBLE PRECISION, DIMENSION(1:2*M-2) :: WrkSpace
  DOUBLE PRECISION, DIMENSION(M,M) :: Z
  phi(0,:)=phi0
  CALL SuperPsi(n, phi(0,:), APhi)
  a(0)=dot_product(phi(0,:),Aphi)
  bphi=Aphi-a(0)*phi(0,:)
  b(0)=SQRT(SUM(bphi*bphi))
  phi(1,:)=bphi/b(0)
  b(0) = dot_product(phi(1,:),APhi)
  DO i=1,m-1
    CALL SuperPsi(n,phi(i,:),APhi)
    a(i)=dot_product(phi(i,:),Aphi)
    bphi=Aphi-(a(i)*phi(i,:))-(b(i-1)*phi(i-1,:))
    b(i)=SQRT(SUM(bphi*bphi))
    IF (abs(b(i))<1.D-19) EXIT
    IF (i<m-1) phi(i+1,:)=bphi/b(i)
    b(i) = dot_product(phi(i+1,:), APhi)
  ENDDO
  print *, 'b'
  CALL DSTEV('V',m, a , b(0:m-2) , Z, m,  wrkspace , info )
  print *, 'alanc',info
  print *, a(0:10)
  DO I=1, N
    phi0(I) = DOT_PRODUCT(Z(1:M,1),Phi(0:M-1,I))
  END DO
END SUBROUTINE Lanczos

END PROGRAM White
    
    