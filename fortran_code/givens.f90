MODULE givens
  USE tools, ONLY: dp, eye, print_matrix
  IMPLICIT NONE
  ! Module to compute a bidiagonal form of a matrix.
  PRIVATE

  ! Global variables
  ! INTEGER ::  i

  PUBLIC  ::  givensparams, chase, wilkinson
CONTAINS

  SUBROUTINE GIVENSPARAMS (f, g, c1, s1)
    ! Subroutine wich return the value of cos and sen for xi and xj.
    IMPLICIT NONE

    REAL(dp), INTENT(in) :: f, g  
    REAL(dp) :: c1, s1   

    REAL(dp) :: u,t,tol                     ! Dummy value

    tol=1.D0*1e-6

    IF ( DABS(g-0.D0) < tol ) THEN                           ! If xj=0, create an ident matrix
      c1=1.D0
      s1=0.D0
    ELSE IF ( DABS(f-0.D0) < tol ) THEN        ! Create sen and cos.
      c1=0.D0
      s1=DSIGN(1.D0,g)
    ELSE IF (DABS(f) > DABS(g)) THEN
      t=g/f
      u=DSIGN(1.D0,f)*DSQRT(1 + (t**2))
      c1=1/u
      s1=c1*t
    ELSE
      t=f/g
      u=DSIGN(1.D0,g)*DSQRT(1 + (t**2))
      s1=1/u
      c1=s1*t

    END IF
  END SUBROUTINE 

  SUBROUTINE chase(H,P)
    IMPLICIT NONE
    REAL(dp) :: H(:,:)
    REAL(dp) :: P(:,:)
    INTEGER :: n, k, i

    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: param, X1, X2, Y1, Y2, Z1, Z2

    n=SIZE(H,1) !Dimension of H
    
    ! ALLOCATE(param(n-1,2),X1(1,n),X2(1,n),Y1(n,1),Y2(n,1),Z1(1,n),Z2(1,n))  
    ALLOCATE(param(n-1,2),X1(1,n),X2(1,n),Y1(n,1),Y2(n,1),Z1(n,1),Z2(n,1))  

    ! P = TRANSPOSE(P)
    param(:,:)=0.D0
    X1(:,:) = 0.D0
    X2(:,:) = 0.D0
    Y1(:,:) = 0.D0
    Y2(:,:) = 0.D0
    Z1(:,:) = 0.D0
    Z2(:,:) = 0.D0

    ! Left mult
    DO k=1, n-1
      CALL givensparams(H(k,k),H(k+1,k),param(k,1),param(k,2))

      DO i=1, n ! Rows of B
        X1(1,i)=H(k,i)
        X2(1,i)=H(k+1,i)
      END DO

      DO i=1, n ! Left
        H(k,i) = param(k,1)*X1(1,i) + param(k,2)*X2(1,i)
        H(k+1,i) = param(k,1)*X2(1,i) - param(k,2)*X1(1,i)
      END DO

    END DO


    ! Right mul
    DO k=1, n-1
      
      Z1=P(1:n,k:k)
      Z2=P(1:n,k+1:k+1)

      DO i=1, n ! Columns of B
        Y1(i,1)=H(i,k)
        Y2(i,1)=H(i,k+1)
      END DO

      DO i=1, n ! Right
        H(i,k)= param(k,1)*Y1(i,1) + param(k,2)*Y2(i,1)
        H(i,k+1)= param(k,1)*Y2(i,1) - param(k,2)*Y1(i,1)

      END DO

      P(1:n,k:k)= param(k,1)*Z1 + param(k,2)*Z2
      P(1:n,k+1:k+1)= param(k,1)*Z2 - param(k,2)*Z1

    END DO

    ! CALL print_matrix(param)
    ! WRITE(*,*) ''

    DEALLOCATE(param,X1,X2,Y1,Y2,Z1,Z2)
  END SUBROUTINE chase

  SUBROUTINE wilkinson(M,w)
    ! Compute the Wilkinson shift


    ! Constants
    REAL(dp), DIMENSION(2,2) :: M
 
    ! Variables
    REAL(dp)  ::  a, b, c, d, E1, E2, dis, w
 
    a=M(1,1)
    b=M(1,2)
    c=M(2,1)
    d=M(2,2)
    ! WRITE(*,*) a, b, c ,d
 
    dis = (a+d)**2 - 4.D0*( (a*d) - (b*c) )
 
    ! If-else of dis
    IF (dis > 0.D0) THEN
      E1= (a+d+ DSQRT(dis))/2.D0
      E2= (a+d- DSQRT(dis))/2.D0
    ELSE IF (DABS(dis-0.D0) < 1.D0*1e-4) THEN
      E1 = (a+d)/2
      E2 = E1
    ELSE
      E1=(a+d)/2
      E2=DSQRT(-1.D0*dis)/2.D0
   END IF
    ! Wilkinson
    IF (DABS(E1-M(2,2)) < DABS(E2 - M(2,2))) THEN
      w=E1
    ELSE
      w=E2
    END IF

  END SUBROUTINE wilkinson

END MODULE givens
