!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                VERSION                        !
!Martin Sebastian Zöllner: Ver. 05.07.2016.14:23!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM Fock_Over_Calc_TM
  IMPLICIT NONE
  INTEGER          :: dime,ierr,i,j
  CHARACTER(LEN=31) :: r
  ierr = Calc(number_aos())
  CONTAINS

  FUNCTION number_aos() RESULT(n_mo)
    IMPLICIT NONE
    INTEGER :: n_mo
    CHARACTER(LEN=31) :: r
    OPEN(UNIT=100, FILE='mos', ACTION='READ')
    DO i = 1,3
      READ(100,*) r
    ENDDO
    READ(100,*) r,r,r,r
    READ(r(7:31), '(1i10)') n_mo
    CLOSE(100)
    RETURN
  END FUNCTION

  FUNCTION Calc(dime) RESULT(worked)
    IMPLICIT NONE
    INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(14)
    INTEGER           :: i, j, worked,dime
    REAL(KIND=db)     :: fock(dime,dime), &
                         over(dime,dime), &
                         eigc(dime,dime), &
                         eige(dime,dime)
    REAL(KIND=db)     :: temp1(dime,dime),temp2(dime,dime),&
                         InvC(dime,dime)
    WRITE(*,*) 'Number of atomic orbitals:',dime
    eige    = 0
    over    = 0
    fock    = 0
    temp1   = 0
    temp2   = 0

    OPEN(UNIT=100, FILE='mos',ACTION='READ')
    OPEN(UNIT=111, FILE='coef',ACTION='WRITE')
    OPEN(UNIT=222, FILE='eigenvalues.dat',ACTION='WRITE')
      DO i = 1,3
        READ(100,'(1A)') r
      ENDDO
      DO i = 1,dime
        READ(100,*) r,r,r
        READ(r(12:31),'(1d20.14)') eige(i,i)
        READ(100,'(4d20.14)') eigc(:,i)
      END DO

      DO i = 1,dime
	write(222,*) (eige(i,j), j=1,dime)
      END DO

      DO i = 1,dime
	write(111,*) (eigc(i,j), j=1,dime)
      END DO

    CLOSE(100)

    Invc = inv(eigc) 
    CALL DGEMM('T','N',dime,dime,dime,1.0D0,Invc,dime,Invc,dime,0.0D0,over,dime)
    CALL DGEMM('N','N',dime,dime,dime,1.0D0,over,dime,eigc,dime,0.0D0,temp1,dime)
    CALL DGEMM('N','N',dime,dime,dime,1.0D0,eige,dime,Invc,dime,0.0D0,temp2,dime)
    CALL DGEMM('N','N',dime,dime,dime,1.0D0,temp1,dime,temp2,dime,0.0D0,fock,dime)

    OPEN(UNIT=1,FILE='hamiltonian.1',ACTION='WRITE')
    OPEN(UNIT=2,FILE='overlap',ACTION='WRITE')
    DO i = 1, dime
      WRITE(1,'(4d20.10)') fock(:,i)
      WRITE(2,'(4d20.10)') over(:,i)
    END DO
    CLOSE(1)
    CLOSE(2)
    worked = -1
    RETURN
  END FUNCTION Calc
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !KOMMENTAR: MARTIN SEBASTIAN ZÖLLNER, 02.06.2016       !
    ! http://fortranwiki.org/fortran/show/Matrix+inversion !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
  FUNCTION inv(A) RESULT(Ainv)
    INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(14)
    REAL(KIND=db), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=db), DIMENSION(SIZE(A,1),SIZE(A,2)) :: Ainv
    REAL(KIND=db), DIMENSION(SIZE(A,1))  :: work  ! work array for LAPACK
    INTEGER,DIMENSION(SIZE(A,1)) :: ipiv   ! pivot indices
    INTEGER :: n, info
    ! External procedures defined in LAPACK
    EXTERNAL DGETRF
    EXTERNAL DGETRI
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = SIZE(A,1)
    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    CALL DGETRF(n, n, Ainv, n, ipiv, info)
    IF (info /= 0) THEN
       STOP 'Matrix is numerically singular!'
    END IF
    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    CALL DGETRI(n, Ainv, n, ipiv, work, n, info)
    IF (info /= 0) THEN
       STOP 'Matrix inversion failed!'
    END IF
  END FUNCTION inv
END PROGRAM
