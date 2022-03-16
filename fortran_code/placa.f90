MODULE placa
  USE tools, only: DP
  IMPLICIT NONE
  ! Module that generate a matrix for a plate.

  PRIVATE
  ! Global constants
  INTEGER :: i, j


  PUBLIC  ::  construir, interior, zona_fila, zona_columna, zona_especial, write_diag
CONTAINS
  SUBROUTINE construir(A,m)
    IMPLICIT NONE
    ! Construccion total de la matriz A

    ! Argument
    INTEGER, INTENT(in) :: m
    REAL(dp), DIMENSION(m**2,m**2) :: A

    CALL interior(A,m)
    CALL zona_fila(A,m)
    CALL zona_columna(A,m)
    CALL zona_especial(A,m)

  END SUBROUTINE construir

  SUBROUTINE interior(A,m)
    IMPLICIT NONE
    ! Construccion del interior de la matriz A.

    ! Argument
    INTEGER, INTENT(in) :: m
    REAL(dp), DIMENSION(m**2,m**2) :: A

    ! Variable
    REAL(dp), DIMENSION(m**2) :: V

    ! Variable
    
    ! Molecula interna completa. A(i,j)
    DO j=3, m-2       ! Constante en la fila
      DO i=3, m-2     ! Varia en la columna

        V(:) = 0.D0   ! Inicializamos

        V((i-2)+((j)- 1)*m)=1.D0     !X(1,3)
        V((i+2)+((j)- 1)*m)=1.D0     !X(5,3)
        V((i)+((j-2)- 1)*m)=1.D0     !X(3,1)
        V((i)+((j+2)- 1)*m)=1.D0     !X(3,5)

        V((i-1)+((j-1)- 1)*m)=2.D0   !X(2,2)
        V((i-1)+((j+1)- 1)*m)=2.D0   !X(2,4)
        V((i+1)+((j-1)- 1)*m)=2.D0   !X(4,2)
        V((i+1)+((j+1)- 1)*m)=2.D0   !X(4,4)

        V((i+1)+((j)- 1)*m)= -8.D0   !X(4,3)
        V((i-1)+((j)- 1)*m)= -8.D0   !X(2,3)
        V((i)+((j-1)- 1)*m)= -8.D0   !X(3,2)
        V((i)+((j+1)- 1)*m)= -8.D0   !X(3,4)

        V((i)+((j)- 1)*m)=20.D0      !X(3,3)

        A(:,i + (j-1)*m)=V(:)       ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...
      END DO
    END DO
  END SUBROUTINE interior

  SUBROUTINE zona_fila(A,m)
    IMPLICIT NONE
    ! Construccion de la zona por filas

    ! Argument
    INTEGER, INTENT(in) :: m
    REAL(dp),DIMENSION(m**2,m**2) :: A

    ! Variable
    REAL(dp), DIMENSION(m**2) :: V

    ! Zona 1
    j=1               ! Primera fila
      DO i=3, m-2     ! Varia en la columna

        V(:) = 0.D0   ! Inicializamos

        V((i-2)+((j)- 1)*m)=1.D0    !X(i-2,j)
        V((i+2)+((j)- 1)*m)=1.D0    !X(i+2,j)
        ! V((i)+((j-2)- 1)*m)=1.D0    !X(i,j-2)
        V((i)+((j+2)- 1)*m)=1.0     !X(i,j+2)

        ! V((i-1)+((j-1)- 1)*m)=2.D0  !X(i-1,j-1)
        V((i-1)+((j+1)- 1)*m)=2.0   !X(i-1,j+1)
        ! V((i+1)+((j-1)- 1)*m)=2.D0  !X(i+1,j-1)
        V((i+1)+((j+1)- 1)*m)=2.0   !X(i+1,j+1)

        V((i+1)+((j)- 1)*m)= -8.D0  !X(i+1,j)
        V((i-1)+((j)- 1)*m)= -8.D0  !X(i-1,j)
        ! V((i)+((j-1)- 1)*m)= -8.D0  !X(i,j-1)
        V((i)+((j+1)- 1)*m)= -8.0   !X(i,j+1)

        V((i)+((j)- 1)*m)=20.0      !X(3,3)

        A(:,i + (j-1)*m)=V(:)       ! Fila (j-1)*m+i
      END DO

    ! Zona 5
    j=2               ! 2da fila
      DO i=3, m-2     ! Varia en la columna

        V(:) = 0.D0   ! Inicializamos

        V((i-2)+((j)- 1)*m)=1.D0    !X(i-2,j)
        V((i+2)+((j)- 1)*m)=1.D0    !X(i+2,j)
        ! V((i)+((j-2)- 1)*m)=1.D0    !X(i,j-2)
        V((i)+((j+2)- 1)*m)=1.0     !X(i,j+2)

        V((i-1)+((j-1)- 1)*m)=2.D0  !X(i-1,j-1)
        V((i-1)+((j+1)- 1)*m)=2.0   !X(i-1,j+1)
        V((i+1)+((j-1)- 1)*m)=2.D0  !X(i+1,j-1)
        V((i+1)+((j+1)- 1)*m)=2.0   !X(i+1,j+1)

        V((i+1)+((j)- 1)*m)= -8.D0  !X(i+1,j)
        V((i-1)+((j)- 1)*m)= -8.D0  !X(i-1,j)
        V((i)+((j-1)- 1)*m)= -8.D0  !X(i,j-1)
        V((i)+((j+1)- 1)*m)= -8.0   !X(i,j+1)

        V((i)+((j)- 1)*m)=20.0      !X(3,3)

        A(:,i + (j-1)*m)=V(:)       ! Fila (j-1)*m+i
      END DO

    ! Zona 7
    j=m-1             ! m-1 fila
      DO i=3, m-2     ! Varia en la columna

        V(:) = 0.D0   ! Inicializamos

        V((i-2)+((j)- 1)*m)=1.D0    !X(i-2,j)
        V((i+2)+((j)- 1)*m)=1.D0    !X(i+2,j)
        V((i)+((j-2)- 1)*m)=1.D0    !X(i,j-2)
        ! V((i)+((j+2)- 1)*m)=1.0     !X(i,j+2)

        V((i-1)+((j-1)- 1)*m)=2.D0  !X(i-1,j-1)
        V((i-1)+((j+1)- 1)*m)=2.0   !X(i-1,j+1)
        V((i+1)+((j-1)- 1)*m)=2.D0  !X(i+1,j-1)
        V((i+1)+((j+1)- 1)*m)=2.0   !X(i+1,j+1)

        V((i+1)+((j)- 1)*m)= -8.D0  !X(i+1,j)
        V((i-1)+((j)- 1)*m)= -8.D0  !X(i-1,j)
        V((i)+((j-1)- 1)*m)= -8.D0  !X(i,j-1)
        V((i)+((j+1)- 1)*m)= -8.0   !X(i,j+1)

        V((i)+((j)- 1)*m)=20.0      !X(3,3)

        A(:,i + (j-1)*m)=V(:)       ! Fila (j-1)*m+i
      END DO

    ! Zona 3
    j=m               ! Ultima fila
      DO i=3, m-2     ! Varia en la columna

        V(:) = 0.D0   ! Inicializamos

        V((i-2)+((j)- 1)*m)=1.D0    !X(i-2,j)
        V((i+2)+((j)- 1)*m)=1.D0    !X(i+2,j)
        V((i)+((j-2)- 1)*m)=1.D0    !X(i,j-2)
        ! V((i)+((j+2)- 1)*m)=1.0     !X(i,j+2)

        V((i-1)+((j-1)- 1)*m)=2.D0  !X(i-1,j-1)
        ! V((i-1)+((j+1)- 1)*m)=2.0   !X(i-1,j+1)
        V((i+1)+((j-1)- 1)*m)=2.D0  !X(i+1,j-1)
        ! V((i+1)+((j+1)- 1)*m)=2.0   !X(i+1,j+1)

        V((i+1)+((j)- 1)*m)= -8.D0  !X(i+1,j)
        V((i-1)+((j)- 1)*m)= -8.D0  !X(i-1,j)
        V((i)+((j-1)- 1)*m)= -8.D0  !X(i,j-1)
        ! V((i)+((j+1)- 1)*m)= -8.0   !X(i,j+1)

        V((i)+((j)- 1)*m)=20.0      !X(i,j)

        A(:,i + (j-1)*m)=V(:)       ! Fila (j-1)*m+i
      END DO
  END SUBROUTINE zona_fila

  SUBROUTINE zona_columna(A,m)
    IMPLICIT NONE
    ! Construccion de las zonas por columnas

    ! Argument
    INTEGER, INTENT(in) :: m
    REAL(dp), DIMENSION(m**2,m**2) :: A

    ! Variable
    REAL(dp), DIMENSION(m**2) :: V

    ! Zona 2
    i= m            ! Columna m
    DO j=3, m-2     ! Constante en la fila

      V(:) = 0.D0   ! Inicializamos

      V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      ! V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      ! V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      ! V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)       ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...
    END DO

    ! Zona 6
    i= m-1            ! Columna m-1
    DO j=3, m-2     ! Constante en la fila

      V(:) = 0.D0   ! Inicializamos

      V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...
    END DO

    ! Zona 8
    i= 2            ! Columna 2
    DO j=3, m-2     ! Variable en la fila

      V(:) = 0.D0   ! Inicializamos

      ! V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...
    END DO

    ! Zona 4
    i= 1            ! Columna 1
    DO j=3, m-2     ! Variable en la fila

      V(:) = 0.D0   ! Inicializamos

      ! V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      ! V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      ! V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      ! V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...
    END DO
  END SUBROUTINE zona_columna

  SUBROUTINE zona_especial(A,m)
    IMPLICIT NONE
    ! Construccion de las zonas alpha, beta, theta y gamma

    ! Argument
    INTEGER, INTENT(in) :: m
    REAL(dp), DIMENSION(m**2,m**2) :: A

    ! Variable
    REAL(dp), DIMENSION(m**2) :: V

    ! Zona alpha 1
    V(:) = 0.D0   ! Inicializamos

    i= 1            ! Columna 1
    j= 1            ! Fila 1

      ! V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      ! V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      ! V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      ! V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      ! V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      ! V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      ! V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    ! Zona alpha 2
    V(:) = 0.D0   ! Inicializamos

    i= 2            ! Columna 2
    j= 1            ! Fila 1

      ! V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      ! V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      ! V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      ! V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      ! V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    ! Zona alpha 3
    V(:) = 0.D0   ! Inicializamos

    i= 1            ! Columna 1
    j= 2            ! Fila 2

      ! V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      ! V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      ! V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      ! V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      ! V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    ! Zona alpha 4
    V(:) = 0.D0   ! Inicializamos

    i= 2            ! Columna 2
    j= 2            ! Fila 2

      ! V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      ! V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    !------------------------------------------------------------------

    ! Zona beta 1
    V(:) = 0.D0   ! Inicializamos

    i=m-1         !Columna
    j=1           !Fila

      V((i-2)+((j)- 1)*m)=1.D0      !X(i-2,j)
      ! V((i+2)+((j)- 1)*m)=1.D0    !X(i+2,j)
      ! V((i)+((j-2)- 1)*m)=1.D0    !X(i,j-2)
      V((i)+((j+2)- 1)*m)=1.D0      !X(i,j+2)

      ! V((i-1)+((j-1)- 1)*m)=2.D0  !X(i-1,j-1)
      V((i-1)+((j+1)- 1)*m)=2.D0    !X(i-1,j+1)
      ! V((i+1)+((j-1)- 1)*m)=2.D0  !X(i+1,j-1)
      V((i+1)+((j+1)- 1)*m)=2.D0    !X(i+1,j+1)

      V((i+1)+((j)- 1)*m)= -8.D0    !X(i+1,j)
      V((i-1)+((j)- 1)*m)= -8.D0    !X(i-1,j)
      ! V((i)+((j-1)- 1)*m)= -8.D0  !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0    !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0       !X(i,j)

      A(:, i + (j-1)*m)=V(:)         ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    ! Zona beta 2
    V(:) = 0.D0   ! Inicializamos

    i=m           ! Columna
    j=1           ! Fila

      V((i-2)+((j)- 1)*m)=1.D0      !X(i-2,j)
      ! V((i+2)+((j)- 1)*m)=1.D0    !X(i+2,j)
      ! V((i)+((j-2)- 1)*m)=1.D0    !X(i,j-2)
      V((i)+((j+2)- 1)*m)=1.D0      !X(i,j+2)

      ! V((i-1)+((j-1)- 1)*m)=2.D0  !X(i-1,j-1)
      V((i-1)+((j+1)- 1)*m)=2.D0    !X(i-1,j+1)
      ! V((i+1)+((j-1)- 1)*m)=2.D0  !X(i+1,j-1)
      ! V((i+1)+((j+1)- 1)*m)=2.D0  !X(i+1,j+1)

      ! V((i+1)+((j)- 1)*m)= -8.D0  !X(i+1,j)
      V((i-1)+((j)- 1)*m)= -8.D0    !X(i-1,j)
      ! V((i)+((j-1)- 1)*m)= -8.D0  !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0    !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0       !X(i,j)

      A(:,i + (j-1)*m)=V(:)         ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    ! Zona beta 3
    V(:) = 0.D0   ! Inicializamos
    i=m-1         ! Columna
    j=2           ! Fila

      V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      ! V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    !! Zona beta 4
    V(:) = 0.D0   ! Inicializamos

    i=m           ! Columna
    j=2           ! fila

      V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      ! V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      ! V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      ! V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      ! V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    !!------------------------------------------------------------------

    ! Zona theta 1
    V(:) = 0.D0   ! Inicializamos

    i=1           ! Columna
    j=m-1         ! Fila

      ! V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      ! V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      ! V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      ! V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    ! Zona theta 2
    V(:) = 0.D0   ! Inicializamos

    i=2           ! Columna
    j=m-1         ! Fila

      ! V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    ! Zona theta 3
    V(:) = 0.D0   ! Inicializamos

    i=1           ! Columna
    j=m           ! Fila

      ! V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      ! V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      ! V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      ! V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      ! V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      ! V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    ! Zona theta 4
    V(:) = 0.D0   ! Inicializamos

    i=2
    j=m

      ! V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      ! V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      ! V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      ! V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    !!------------------------------------------------------------------

    ! Zona gamma 1
    V(:) = 0.D0   ! Inicializamos

    i=m-1         ! Columna
    j=m-1         ! Fila

      V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    ! Zona gamma 2
    V(:) = 0.D0   ! Inicializamos

    i=m           ! Columna
    j=m-1         ! Fila

      V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      ! V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      ! V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      ! V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    ! Zona gamma 3
    V(:) = 0.D0   ! Inicializamos

    i=m-1         ! Columna
    j=m           ! Fila

      V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      ! V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      ! V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      ! V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    ! Zona gamma 4
    V(:) = 0.D0   ! Inicializamos

    i=m           ! Columna
    j=m           ! Fila

      V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
      ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
      V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
      ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

      V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
      ! V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
      ! V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
      ! V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

      ! V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
      V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
      V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
      ! V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

      V((i)+((j)- 1)*m)=20.D0      !X(i,j)

      A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    !------------------------------------------------------------------

  END SUBROUTINE zona_especial

  SUBROUTINE write_diag(A,filename)
    ! Make the write of a matrix

    ! Constants
    CHARACTER (len=*) :: filename
    REAL(dp) :: A(:,:)

    ! Variables
    INTEGER ::  m, n, unidad

    m=SIZE(A,1) ! Rows of A
    n=SIZE(A,2) ! Columns of A

    OPEN(newunit=unidad, file=filename, status='replace', action='write')

    ! Write this form (matrix)
    DO i=1, m
      WRITE(unidad,'(ES14.4)')A(i,i)
    END DO

    CLOSE(unidad)
  END SUBROUTINE write_diag

END MODULE placa

    !! Zona beta 1
    !i= m-1            ! Columna m-1
    !j= 1            ! Fila 1

    !  V(:) = 0.D0   ! Inicializamos

    !  V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
    !  ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
    !  ! V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
    !  V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

    !  ! V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
    !  V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
    !  ! V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
    !  V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

    !  V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
    !  V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
    !  ! V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
    !  V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

    !  V((i)+((j)- 1)*m)=20.D0      !X(i,j)

    !  A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    !! Zona beta 2
    !i= m            ! Columna m
    !j= 1            ! Fila 1

    !  V(:) = 0.D0   ! Inicializamos

    !  V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
    !  ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
    !  ! V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
    !  V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

    !  ! V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
    !  V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
    !  ! V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
    !  ! V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

    !  ! V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
    !  V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
    !  ! V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
    !  V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

    !  V((i)+((j)- 1)*m)=20.D0      !X(i,j)

    !  A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    !! Zona beta 3
    !i= m-1            ! Columna m-1
    !j= 2              ! Fila 2

    !  V(:) = 0.D0   ! Inicializamos

    !  V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
    !  ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
    !  ! V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
    !  V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

    !  V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
    !  V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
    !  V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
    !  V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

    !  V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
    !  V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
    !  V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
    !  V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

    !  V((i)+((j)- 1)*m)=20.D0      !X(i,j)

    !  A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    !! Zona beta 4
    !i= m              ! Columna m
    !j= 2              ! Fila 2

    !  V(:) = 0.D0   ! Inicializamos

    !  V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
    !  ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
    !  ! V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
    !  V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

    !  V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
    !  V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
    !  ! V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
    !  ! V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

    !  ! V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
    !  V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
    !  V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
    !  V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

    !  V((i)+((j)- 1)*m)=20.D0      !X(i,j)

    !  A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    !!------------------------------------------------------------------

    !! Zona theta 1
    !i= 1              ! Columna 1
    !j= m-1            ! Fila m-1

    !  V(:) = 0.D0   ! Inicializamos

    !  ! V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
    !  V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
    !  V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
    !  ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

    !  ! V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
    !  ! V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
    !  V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
    !  V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

    !  V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
    !  ! V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
    !  V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
    !  V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

    !  V((i)+((j)- 1)*m)=20.D0      !X(i,j)

    !  A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    !! Zona theta 2
    !i= 2              ! Columna 2
    !j= m-1            ! Fila m-1

    !  V(:) = 0.D0   ! Inicializamos

    !  ! V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
    !  V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
    !  V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
    !  ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

    !  V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
    !  V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
    !  V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
    !  V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

    !  V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
    !  V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
    !  V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
    !  V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

    !  V((i)+((j)- 1)*m)=20.D0      !X(i,j)

    !  A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    !! Zona theta 3
    !i= 1              ! Columna 1
    !j= m              ! Fila m

    !  V(:) = 0.D0   ! Inicializamos

    !  ! V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
    !  V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
    !  V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
    !  ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

    !  ! V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
    !  ! V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
    !  V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
    !  ! V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

    !  V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
    !  ! V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
    !  V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
    !  ! V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

    !  V((i)+((j)- 1)*m)=20.D0      !X(i,j)

    !  A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    !! Zona theta 4
    !i= 2              ! Columna 2
    !j= m              ! Fila m

    !  V(:) = 0.D0   ! Inicializamos

    !  ! V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
    !  V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
    !  V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
    !  ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

    !  V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
    !  ! V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
    !  V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
    !  ! V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

    !  V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
    !  V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
    !  V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
    !  ! V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

    !  V((i)+((j)- 1)*m)=20.D0      !X(i,j)

    !  A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    !!------------------------------------------------------------------


    ! ! Zona gamma 1
    ! i= m-1              ! Columna m-1
    ! j= m-1              ! Fila m-1

    !   V(:) = 0.D0   ! Inicializamos

    !   V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
    !   ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
    !   V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
    !   ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

    !   V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
    !   V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
    !   V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
    !   V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

    !   V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
    !   V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
    !   V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
    !   V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

    !   V((i)+((j)- 1)*m)=20.D0      !X(i,j)

    !   A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    ! ! Zona gamma 2
    ! i= m                ! Columna m-1
    ! j= m-1              ! Fila m-1

    !   V(:) = 0.D0   ! Inicializamos

    !   V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
    !   ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
    !   V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
    !   ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

    !   V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
    !   V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
    !   ! V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
    !   ! V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

    !   ! V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
    !   V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
    !   V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
    !   V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

    !   V((i)+((j)- 1)*m)=20.D0      !X(i,j)

    !   A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    ! ! Zona gamma 3
    ! i= m-1                ! Columna m-1
    ! j= m                  ! Fila m-1

    !   V(:) = 0.D0   ! Inicializamos

    !   V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
    !   ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
    !   V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
    !   ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

    !   V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
    !   ! V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
    !   V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
    !   ! V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

    !   V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
    !   V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
    !   V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
    !   ! V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

    !   V((i)+((j)- 1)*m)=20.D0      !X(i,j)

    !   A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...

    ! ! Zona gamma 4
    ! i= m                ! Columna m-1
    ! j= m                  ! Fila m-1

    !   V(:) = 0.D0   ! Inicializamos

    !   V((i-2)+((j)- 1)*m)=1.D0     !X(i-2,j)
    !   ! V((i+2)+((j)- 1)*m)=1.D0     !X(i+2,j)
    !   V((i)+((j-2)- 1)*m)=1.D0     !X(i,j-2)
    !   ! V((i)+((j+2)- 1)*m)=1.D0     !X(i,j+2)

    !   V((i-1)+((j-1)- 1)*m)=2.D0   !X(i-1,j-1)
    !   ! V((i-1)+((j+1)- 1)*m)=2.D0   !X(i-1,j+1)
    !   ! V((i+1)+((j-1)- 1)*m)=2.D0   !X(i+1,j-1)
    !   ! V((i+1)+((j+1)- 1)*m)=2.D0   !X(i+1,j+1)

    !   ! V((i+1)+((j)- 1)*m)= -8.D0   !X(i+1,j)
    !   V((i-1)+((j)- 1)*m)= -8.D0   !X(i-1,j)
    !   V((i)+((j-1)- 1)*m)= -8.D0   !X(i,j-1)
    !   ! V((i)+((j+1)- 1)*m)= -8.D0   !X(i,j+1)

    !   V((i)+((j)- 1)*m)=20.D0      !X(i,j)

    !   A(:,i + (j-1)*m)=V(:)        ! Fila 2m+3 A(1,2m+3), A(2,2m+3)...
