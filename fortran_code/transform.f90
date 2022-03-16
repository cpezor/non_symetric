MODULE transform
  USE tools, ONLY: dp, eye
  IMPLICIT NONE
  ! Module to compute a bidiagonal form of a matrix.
  PRIVATE

  ! Global variables
  INTEGER ::  i, k

  PUBLIC  ::  householder, hessenberg
CONTAINS
  SUBROUTINE householder(B,U,beta)
    ! Subroutine that compute householder beta and u

    ! Parameters
    REAL(dp) :: B(:,:)
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: U
    REAL(dp) :: beta

    ! Varibles
    INTEGER ::  n
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: X
    REAL(dp)  ::  norma_x, norma_u, valor_max

    n = SIZE(B,1)

    ALLOCATE(U(n,1))  
    U(:,:)=0.D0

    ALLOCATE(X(n))  
    X(:)=0.D0

    DO i=1, n
      X(i)=B(i,1)
    END DO

    valor_max=MAXVAL(DABS(X))
    X = (1.D0/valor_max)*X
    norma_x=NORM2(X)
    X(1) = X(1) + SIGN(norma_x,X(1))
    norma_u=NORM2(X)
    beta = 2.D0/(norma_u**2)

    DO i=1, n
      U(i,1)=X(i) ! First column
    END DO

    DEALLOCATE(X) ! WORK
  END SUBROUTINE householder

  SUBROUTINE hessenberg(B,P)
    IMPLICIT NONE
    ! Compute the Hessenberg form

    ! Parameters
    REAL(dp) :: B(:,:)
    REAL(dp) :: P(:,:)

    ! Varibles
    INTEGER n
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: Z
    REAL(dp) :: beta

    n = SIZE(B,1)
    CALL eye(P)

    DO k=1, n-2
      CALL householder(B(k+1:n,k:k), Z, beta) ! Allocate Z

      B(k+1:n,k:n) = B(k+1:n,k:n) - (beta)*MATMUL(Z,MATMUL(TRANSPOSE(Z),B(k+1:n,k:n)))

      B(1:n,k+1:n) = B(1:n,k+1:n) - (beta)*MATMUL(MATMUL(B(1:n,k+1:n),Z),TRANSPOSE(Z))


      P(1:n,k+1:n) = P(1:n,k+1:n) - (beta)*MATMUL(MATMUL(P(1:n,k+1:n),Z),TRANSPOSE(Z))

      DEALLOCATE(Z) ! WORK
    END DO
  END SUBROUTINE hessenberg

END MODULE transform
