

      SUBROUTINE GenEig(HMatrix, SMatrix, N, MaxN, Diag)

C Driver subroutine for solving the generalised eigenvalue problem using the
C LAPACK subroutine DSYGV. The generalised eigenvalue problem is of the form
C H C = lambda S C, where H is HMatrix, S is SMatrix. These should both be 
C symmetric. 
C Parameters: HMatrix: matrix with leading dimension MaxN
C             SMatrix: matrix with leading dimension MaxN
C             N:       only the NxN blocks of HMatrix and SMatrix are used
C             MaxN:    leading dimension of HMatrix and SMatrix
C             Diag:    array containing the eigenvalues. Its size should be 
C                      at least N.

      IMPLICIT NONE
      INTEGER MaxWork, MaxN, N, WorkDim, ErrInfo
      PARAMETER (MaxWork = 100000)
      DOUBLE PRECISION SMatrix(MaxN, MaxN), HMatrix(MaxN, MaxN),
     .                 Work(MaxWork), Diag(MaxN)

      WorkDim = MaxWork
      CALL DSYGV(1, 'v', 'u', N, HMatrix, MaxN, SMatrix, MaxN, Diag, 
     .           Work, WorkDim, ErrInfo)
      IF (ErrInfo.LT.0) THEN
        WRITE (6,*) 'Wrong value of argument ', -ErrInfo, 
     .              ' in subroutine call to DSYGV'
        STOP
      ELSE IF ((ErrInfo.GT.0).AND.(ErrInfo.LE.N)) THEN
        WRITE (6,*) 'Convergence failure routine DSYGV'
        STOP
      ELSE IF (ErrInfo.GT.N) THEN
        WRITE (6,*) 'Leading minor of order', ErrInfo-N,
     .              ' is not positive definite'
        STOP
      ENDIF
      END
