MODULE tools
  IMPLICIT NONE
  ! Module to define precision of real variables.

  PRIVATE

  ! Module constants
  INTEGER            ::  i, j
  INTEGER, PARAMETER :: dp  = SELECTED_REAL_KIND(P=10, R=30)
  INTEGER, PARAMETER :: sp  = SELECTED_REAL_KIND(P=5, R=15)
  INTEGER, PARAMETER :: dp_alt  = KIND(0.D0)

  REAL(dp), PARAMETER :: pi  = 4.D0*DATAN(1.D0)

  PUBLIC  ::  DP, SP, pi, read_matrix, write_matrix, print_matrix, write_array, read_bin, write_bin, eye, sub_diag, diag, sorting
CONTAINS
  SUBROUTINE sub_diag(A,D)
    ! Diagonal terms of a matrix

    ! Constants
    REAL(dp):: A(:,:)
    REAL(dp):: D(:)

    ! Varibles
    INTEGER ::  m

    m= SIZE(A,1)

    D(:)=0.D0

    DO i=1, m-1
      D(i)=A(i+1,i)
    END DO
  END SUBROUTINE sub_diag

  SUBROUTINE diag(A,D)
    ! Diagonal terms of a matrix

    ! Constants
    REAL(dp):: A(:,:)
    REAL(dp):: D(:)

    ! Varibles
    INTEGER ::  m

    m= SIZE(A,1)

    D(:)=0.D0

    DO i=1, m
      D(i)=A(i,i)
    END DO
  END SUBROUTINE diag

  SUBROUTINE eye(A)
    ! Make an identidy matrix

    ! Constants
    REAL(dp):: A(:,:)

    ! Varibles
    INTEGER ::  m, n

    m= SIZE(A,1)
    n= SIZE(A,2)

    A(:,:)=0.D0

    DO i=1, m
      A(i,i)=1.D0
    END DO
  END SUBROUTINE eye

  SUBROUTINE print_matrix(A)
    ! Print a matrix on a terminal

    ! Constants
    REAL(dp):: A(:,:)

    ! Varibles
    INTEGER ::  m, n

    m= SIZE(A,1)
    n= SIZE(A,2)

    DO i=1, m
     WRITE(*,'(20F7.2)') (A(i,j), j=1,n)
    END DO
  END SUBROUTINE print_matrix

  SUBROUTINE read_matrix(A,filename)
    ! Make the read of a matrix

    ! Constants
    CHARACTER (len=*) :: filename
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: A

    ! Variables
    INTEGER ::  m, n, unidad

    OPEN(newunit=unidad, file=filename, status='old', action='read')
    READ(unidad,*) m, n

    ALLOCATE(A(m,n)) ! WORK
    A(:,:) = 0.D0

    ! Write this form (matrix for COLUMNS)
    DO j=1,n
      DO i=1, m
        READ(unidad,*)A(i,j)
      END DO
    END DO
   CLOSE(unidad)

  END SUBROUTINE read_matrix

  SUBROUTINE write_matrix(A,filename)
    ! Make the write of a matrix

    ! Constants
    CHARACTER (len=*) :: filename
    REAL(dp) :: A(:,:)

    ! Variables
    INTEGER ::  m, n, unidad

    m=SIZE(A,1) ! Rows of A
    n=SIZE(A,2) ! Columns of A

    OPEN(newunit=unidad, file=filename, status='replace', action='write')
    WRITE(unidad,'(2I5)') m, n

    ! Write this form (matrix)
    DO j=1, n
      DO i=1, m
        WRITE(unidad,'(ES14.4)')A(i,j)
      END DO
    END DO
    CLOSE(unidad)
  END SUBROUTINE write_matrix

  SUBROUTINE write_array(A,filename)
    ! Make the write of a matrix

    ! Constants
    CHARACTER (len=*) :: filename
    REAL(dp) :: A(:)

    ! Variables
    INTEGER ::  m, unidad

    m=SIZE(A) ! Rows of A

    OPEN(newunit=unidad, file=filename, status='replace', action='write')
    WRITE(unidad,'(I5)') m

    ! Write this form (matrix)
    DO i=1, m
      WRITE(unidad,'(ES14.4)')A(i)
    END DO
    CLOSE(unidad)
  END SUBROUTINE write_array

  SUBROUTINE read_bin(A,filename)
    ! Make the read of a matrix (element an element for COLUMNS)

    ! Constants
    CHARACTER (len=*) :: filename
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: A

    ! Variables
    INTEGER ::  m, n, unidad

    OPEN(newunit=unidad, file=filename, status='old', form='unformatted', action='read')
    READ(unidad) m 
    READ(unidad) n

    ALLOCATE(A(m,n)) ! WORK
    A(:,:) = 0.D0

    ! Write this form (matrix for COLUMNS)
    DO j=1,n
      DO i=1, m
        READ(unidad)A(i,j)
      END DO
    END DO
   CLOSE(unidad)
  END SUBROUTINE read_bin

  SUBROUTINE write_bin(A,filename)
    ! Make the write of a matrix (binaries and Columns form)

    ! Constants
    CHARACTER (len=*) :: filename
    REAL(dp) :: A(:,:)

    ! Variables
    INTEGER ::  m, n, unidad

    m=SIZE(A,1) ! Rows of A
    n=SIZE(A,2) ! Columns of A

    OPEN(newunit=unidad, file=filename, status='replace', form='unformatted', action='write')
    WRITE(unidad) m 
    WRITE(unidad) n

    ! Write this form (matrix for COLUMNS)
    DO j=1,n
      DO i=1, m
        WRITE(unidad)A(i,j)
      END DO
    END DO
    CLOSE(unidad)
  END SUBROUTINE write_bin

  SUBROUTINE sorting(D,P)
    IMPLICIT NONE
    ! Make a sorting algorithm

    ! Parameters
    REAL(dp) :: P(:,:)
    REAL(dp) :: D(:)

    ! Varaibles
    INTEGER :: n, pos, i
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: column_current
    REAL(dp) :: current

    n=SIZE(P,1)
    ALLOCATE(column_current(n)) ! WORK
    column_current(:)=0.D0

    DO i= 2, n
      column_current= P(:,i)
      current = D(i)
      pos = i

      DO WHILE( pos>1 .AND. D(pos-1) < current)
        D(pos) = D(pos-1)
        P(:,pos) = P(:,pos-1)
        pos=pos-1
      END DO
      D(pos) = current
      P(:,pos) = column_current
    END DO
    DEALLOCATE(column_current) ! WORK

  END SUBROUTINE sorting
END MODULE tools
