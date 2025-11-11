! Programming exercise, section 11.6 of the textbook `computational physics',
! second edition, Cambridge 2007 
! Jos Thijssen
! 2006-2007
PROGRAM QMRG
  IMPLICIT NONE
  INTEGER, PARAMETER :: M=8 ! This is the parameter m of the paper; should be a multiple of 4.
  DOUBLE PRECISION, ALLOCATABLE :: H_00(:,:), H_10(:,:), H_01(:,:), H_11(:,:), T(:,:), &
                                   BigH_00(:,:), BigH_10(:,:), BigH_01(:,:), BigH_11(:,:),BigT(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: EigVecs(:,:), EigVals(:)
  
  CALL Initialise
  CALL RunQMRG
  
CONTAINS
  
SUBROUTINE Initialise

  INTEGER :: MM

  WRITE (*,*) 'Give number of states M'
  READ *, MM
  ALLOCATE (H_00(2*M,2*M))
  ALLOCATE (H_10(2*M,2*M))
  ALLOCATE (H_01(2*M,2*M))
  ALLOCATE (H_11(2*M,2*M))
  ALLOCATE (T(2*M,2*M))
  ALLOCATE (EigVecs(M, 2*M))
  ALLOCATE (EigVals(2*M))
  ALLOCATE (BigH_00(2*M, 2*M))
  ALLOCATE (BigH_10(2*M, 2*M))
  ALLOCATE (BigH_01(2*M, 2*M))
  ALLOCATE (BigH_11(2*M, 2*M))
  ALLOCATE (BigT(2*M, 2*M))
  H_00(1,:) = (/ 2.d0, -1.d0/)
  H_00(2,:) = (/-1.d0,  2.d0/)
  H_10(1,:) = (/ 1.d0, -1.d0/)
  H_10(2,:) = (/-1.d0,  2.d0/)
  H_01(1,:) = (/ 2.d0, -1.d0/)
  H_01(2,:) = (/-1.d0,  1.d0/)
  H_11(1,:) = (/ 1.d0, -1.d0/)
  H_11(2,:) = (/-1.d0,  1.d0/)
  T(1,:) = (/  0.d0,  0.d0/)
  T(2,:) = (/ -1.d0,  0.d0/)
END SUBROUTINE Initialise


SUBROUTINE RunQMRG
  INTEGER :: Level, CurDim, NewCurDim, MaxLevel=11, I
  CurDim = 2 ! Only first iteration, later twice as big
  DO Level = 1,MaxLevel
    CALL Diagonalise (H_00(1:CurDim,1:CurDim), EigVecs(1:M,1:CurDim), 1, CurDim)
    CALL Diagonalise (H_01(1:CurDim,1:CurDim), EigVecs(1:M,1:CurDim), 2, CurDim)
    CALL Diagonalise (H_10(1:CurDim,1:CurDim), EigVecs(1:M,1:CurDim), 3, CurDim)
    CALL Diagonalise (H_11(1:CurDim,1:CurDim), EigVecs(1:M,1:CurDim), 4, CurDim)
! The eigenvalues and vectors of H_00 and H_11 contain representative eigenvalues of the entire system
! whose quality increase every iteration
    CALL GramSchmidt(EigVecs(1:M,1:CurDim), CurDim, NewCurDim, M) ! Basis is truncated to M/4+orth
                                                    ! NewCurdim is returned with new value!!
    CALL Transform(H_00(1:CurDim,1:CurDim), EigVecs(1:NewCurdim,1:CurDim), NewCurDim, CurDim) ! Results in 
                                                 ! NewCurdim by NewCurDim matrices
    CALL Transform(H_01(1:CurDim,1:CurDim), EigVecs(1:NewCurdim,1:CurDim), NewCurDim, CurDim)
    CALL Transform(H_10(1:CurDim,1:CurDim), EigVecs(1:NewCurdim,1:CurDim), NewCurDim, CurDim)
    CALL Transform(H_11(1:CurDim,1:CurDim), EigVecs(1:NewCurdim,1:CurDim), NewCurDim, CurDim)
    CALL Transform(   T(1:CurDim,1:CurDim), EigVecs(1:NewCurdim,1:CurDim), NewCurDim, CurDim)
    print *, 'curdim, newcurdim', curdim, newcurdim
    CurDim = NewCurDim
    CALL MakeBigMat(BigH_00(1:2*CurDim, 1:2*CurDim), H_00(1:CurDim, 1:CurDim), &
                    H_00(1:CurDim, 1:CurDim), T(1:CurDim, 1:CurDim), CurDim)
    CALL MakeBigMat(BigH_01(1:2*CurDim, 1:2*CurDim), H_00(1:CurDim, 1:CurDim), &
                    H_01(1:CurDim, 1:CurDim), T(1:CurDim, 1:CurDim), CurDim)
    CALL MakeBigMat(BigH_10(1:2*CurDim, 1:2*CurDim), H_10(1:CurDim, 1:CurDim), &
                    H_00(1:CurDim, 1:CurDim), T(1:CurDim, 1:CurDim), CurDim)
    CALL MakeBigMat(BigH_11(1:2*CurDim, 1:2*CurDim), H_10(1:CurDim, 1:CurDim), &
                    H_01(1:CurDim, 1:CurDim), T(1:CurDim, 1:CurDim), CurDim)
    CALL MakeBigT(BigT(1:2*CurDim, 1:2*CurDim), T(1:CurDim, 1:CurDim), CurDim)
    H_00 = BigH_00
    H_01 = BigH_01
    H_10 = BigH_10
    H_11 = BigH_11
    T = BigT
    CurDim = 2*CurDim
  END DO
END SUBROUTINE RunQMRG

SUBROUTINE TransForm(H, EigVecs, NewCurDim, CurDim)
  INTEGER :: CurDim, NewCurDim
  DOUBLE PRECISION :: H(CurDim,CurDim), EigVecs(NewCurDim, CurDim)
  H(1:NewCurDim, 1:CurDim) = MatMul(EigVecs, H)
  H(1:NewCurDim,1:NewCurDim) = MatMul(H(1:NewCurDim,1:CurDim), TRANSPOSE(EigVecs))
END SUBROUTINE TransForm
  

SUBROUTINE Diagonalise(Ham, EigVecs, Num, Size)
  INTEGER, INTENT(IN) :: Size, Num
  DOUBLE PRECISION :: Ham(Size,Size), EigVecs(M,Size), TempHam(Size, Size), &
                      EigVals(Size)
  INTEGER :: LDWork, INFO, I
  DOUBLE PRECISION :: Work(3*Size-1)
  LDWork=3*Size-1
  TempHam = Ham
  CALL DSYEV('v', 'u', Size, TempHam, Size, EigVals, Work, LDWork, INFO)
  IF (INFO/=0) THEN
    print *, INFO
    STOP
  END IF
  print *,'eigvals', EigVals(1:6)
  EigVecs(M/4*(Num-1)+1:M/4*Num, 1:Size) = TRANSPOSE(TempHam(1:Size, 1:M/4)) 
! Fill EigVec with correct eigenvectors
END SUBROUTINE Diagonalise



SUBROUTINE GramSchmidt(EigVecs, CurDim, NewCurDim, M)
  INTEGER :: CurDim, M, NewCurDim, I, J, Count, ListI(M)
  DOUBLE PRECISION :: EigVecs(M, CurDim), Norm

  Count = 0
  DO I=1, M
    DO J=1, I-1
      EigVecs(I,:) = EigVecs(I,:) - Dot_Product(EigVecs(I,:), EigVecs(J,:))*EigVecs(J,:)
    END DO
    Norm = SQRT(Dot_Product(EigVecs(I,:), EigVecs(I,:)))
    IF (Norm>1.D-6) THEN
      Count = Count+1
      ListI(Count) = I
      EigVecs(I,:) = EigVecs(I,:)/Norm
    ELSE
    END IF
  END DO
  NewCurDim = Count
  DO Count = 1, NewCurDim
    EigVecs(Count,:) = EigVecs(ListI(Count), :)
  END DO
END SUBROUTINE GramSchmidt

SUBROUTINE MakeBigMat(BigMat, ULMat, LRMat, T, CurDim)
  INTEGER :: CurDim
  DOUBLE PRECISION :: BigMat(2*CurDim, 2*CurDim), ULMat(CurDim, CurDim), &
                      LRMat(CurDim,CurDim), T(CurDim, CurDim)
  BigMat(1:CurDim, 1:CurDim) = ULMat
  BigMat(CurDim+1:2*CurDim, CurDim+1:2*CurDim) = LRMat
  BigMat(1:CurDim, CurDim+1:2*CurDim) = T
  BigMat(CurDim+1:2*CurDim, 1:CurDim) = TRANSPOSE(T)
END SUBROUTINE MakeBigMat

SUBROUTINE MakeBigT(BigMat, T, CurDim)
  INTEGER :: CurDim
  DOUBLE PRECISION :: BigMat(2*CurDim, 2*CurDim), T(CurDim, CurDim)
  BigMat(1:CurDim, 1:CurDim) = 0.D0
  BigMat(CurDim+1:2*CurDim, CurDim+1:2*CurDim) = 0.D0
  BigMat(1:CurDim, CurDim+1:2*CurDim) = 0.D0
  BigMat(CurDim+1:2*CurDim, 1:CurDim) = T
END SUBROUTINE MakeBigT

END PROGRAM QMRG  
  
  