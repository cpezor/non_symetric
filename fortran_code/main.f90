PROGRAM generador
  !Program that computes QR itaration for a symmetric matrix.
  USE tools
  USE placa, ONLY: construir, write_diag
  USE transform, ONLY: householder, hessenberg
  USE givens, ONLY: chase, wilkinson

  IMPLICIT NONE
  ! Variable
  ! INTEGER :: dis, m, n, iter, var, i, j
  INTEGER :: dis, m, n, iter, i
  REAL(dp), ALLOCATABLE, DIMENSION(:)   :: X, D1, D2, EIG_VAL
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: A, B, P, H, Iden, EIG_VEC
  REAL(dp) :: mu, norma, tol

  ! dis=6  ! Numero de discretizacion.

  ! m=dis-2

  ! ALLOCATE(A(m**2,m**2))  ! Matriz de discretizacion
  ! A(:,:) = 0.D0

  ! CALL construir(A,m) ! Matriz de discretizacion

  CALL read_matrix(A,'datos/matriz_prueba.dat')
  n=SIZE(A,1)

  ALLOCATE(EIG_VAL(n),X(n), D1(n-1), D2(n-1))  ! D1, D2 sub_diag elements  
  EIG_VAL(:) = 0.D0
  X(:) = 0.D0
  D1(:) = 0.D0
  D2(:) = 0.D0

  ALLOCATE(EIG_VEC(n,n), B(n,n), P(n,n), H(n,n), Iden(n,n))  
  EIG_VEC(:,:) = 0.D0
  B(:,:) = 0.D0
  P(:,:) = 0.D0
  H(:,:) = 0.D0
  Iden(:,:) = 0.D0

  ! CALL write_matrix(A,'discre.dat')
  B = A ! Save
  ! Householder
  CALL eye(EIG_VEC)
  CALL eye(P)
  CALL eye(Iden)

  CALL hessenberg(B,P) ! work
  H = B

  tol = 1.D0*1e-5
  iter=0
  norma=1000.D0 ! inicio
  DO WHILE ( norma>=tol .AND. iter<=2000)

    CALL sub_diag(H,D1)
    CALL wilkinson(H(n-1:n,n-1:n),mu)
    H= H- mu*Iden
    CALL chase(H,P)
    H = H + mu*Iden
    CALL sub_diag(H,D2)
    norma=NORM2(D2-D1)

    iter= iter+1
  END DO

  ! DO WHILE ( norma>=tol .AND. iter<=2000)

  !   CALL sub_diag(H,D1)
  !   mu=H(n,n)
  !   H= H- mu*Iden
  !   CALL chase(H,P)
  !   H = H + mu*Iden
  !   CALL sub_diag(H,D2)
  !   norma=NORM2(D2-D1)

  !   iter= iter+1
  ! END DO

  CALL diag(H,EIG_VAL)
  EIG_VEC=P
  CALL sorting(EIG_VAL,EIG_VEC)
  ! Resultados

  CALL write_array(EIG_VAL,'autovalores.dat')
  CALL write_matrix(EIG_VEC,'autovector.dat')
  CALL print_matrix(A)

  DEALLOCATE(A,B,X,P,H,Iden,D1,D2)

END PROGRAM generador
