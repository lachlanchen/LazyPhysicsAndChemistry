      SUBROUTINE Simpson(Weights, Left, Right, Num, Func, Result)
      IMPLICIT NONE
      EXTERNAL Func
      DOUBLE PRECISION Func, Left, Right, Result, Weights(4), X,Step
      INTEGER Num, I
      Result = 0.D0
      Step = (Right-Left)/(Num-1)
      DO I=1, 4
        X = Left + Step*DBLE(I-1)
        Result = Result + Func(X)*Weights(I)
      END DO
      DO I=5, Num-4
        X = Left + Step*DBLE(I-1)
        Result = Result + Func(X)
      END DO
      DO I=Num-3, Num
        X = Left + Step*DBLE(I-1)
        Result = Result + Weights(Num-I+1)*Func(X)
      END DO
      Result = Result*Step
      END

