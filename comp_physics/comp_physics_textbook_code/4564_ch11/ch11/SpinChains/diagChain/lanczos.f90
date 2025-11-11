SUBROUTINE Lanczos(phi0,n,m, MultPsi)
 
  INTEGER, INTENT(IN)::n,m
  DOUBLE PRECISION, DIMENSION(n), INTENT(INOUT)::phi0
  INTEGER::i, INFO
  DOUBLE PRECISION, DIMENSION(1:n) ::bphi,Aphi
  DOUBLE PRECISION, DIMENSION(0:m-1,1:n) ::phi
  DOUBLE PRECISION, DIMENSION(0:M-1) :: a,b
  DOUBLE PRECISION, DIMENSION(1:2*M-2) :: WrkSpace
  DOUBLE PRECISION, DIMENSION(M,M) :: Z
INTERFACE
  SUBROUTINE MultPsi (X, Y, N)
    INTEGER :: N
    DOUBLE PRECISION :: X(N), Y(N) 
  END SUBROUTINE MultPsi
END INTERFACE

  phi(0,:)=phi0/SQRT(dot_product(phi0,phi0))
  CALL MultPsi(phi(0,:), APhi, N)
  a(0)=dot_product(phi(0,:),Aphi)
  bphi=Aphi-a(0)*phi(0,:)
  b(0)=SQRT(SUM(bphi*bphi))
  phi(1,:)=bphi/b(0)
  b(0) = dot_product(phi(1,:),APhi)
  DO i=1,m-1
    CALL MultPsi(phi(i,:),APhi, N)
    a(i)=dot_product(phi(i,:),Aphi)
    bphi=Aphi-(a(i)*phi(i,:))-(b(i-1)*phi(i-1,:))
    b(i)=SQRT(SUM(bphi*bphi))
    IF (abs(b(i))<1.D-9) EXIT
    IF (i<m-1) phi(i+1,:)=bphi/b(i)
    b(i) = dot_product(phi(i+1,:), APhi)
  ENDDO
  CALL DSTEV('V',m, a , b(0:m-2) , Z, m,  wrkspace , info )
  print *, 'a',info
  print *, a(0:10)
  return
! The rest is for checking.....
  DO I=1, N
    phi0(I) = DOT_PRODUCT(Z(1:M,1),Phi(0:M-1,I))
  END DO
  CALL MultPsi(phi0, APhi, N)
  DO I=1, 4
  print *, Phi0(I), APhi(I), APhi(I)/Phi0(I), M
  END DO
  print *, APhi(1:4)/Phi0(1:4)
ENDSUBROUTINE Lanczos
