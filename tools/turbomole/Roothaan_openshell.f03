!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                VERSION                        !
!Martin Sebastian Zöllner: Ver. 05.07.2016.14:23!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM Fock_Over_Calc_TM
  IMPLICIT NONE
  INTEGER          :: dime,ierr,i
  CHARACTER(LEN=31) :: r
  ierr = Calc(number_aos())
  CONTAINS

  FUNCTION number_aos() RESULT(n_mo)
    IMPLICIT NONE
    INTEGER :: n_mo_a,n_mo_b,n_mo
    CHARACTER(LEN=31) :: r
    OPEN(UNIT=100, FILE='alpha', ACTION='READ')
    DO i = 1,3
      READ(100,*) r
    ENDDO
    READ(100,*) r,r,r,r
    READ(r(7:31), '(1i5)') n_mo_a
    CLOSE(100)
    OPEN(UNIT=100, FILE='beta', ACTION='READ')
    DO i = 1,3
      READ(100,*) r
    ENDDO
    READ(100,*) r,r,r,r
    READ(r(7:31), '(1i5)') n_mo_b
    CLOSE(100)
    IF(n_mo_a == n_mo_b) THEN
      n_mo = n_mo_a
      RETURN
    ELSE
      STOP("Number of alpha orbitals does not equals the number of beta orbitals")
    ENDIF
  END FUNCTION


  FUNCTION Calc(dime) RESULT(worked)
    IMPLICIT NONE
    INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(14)
    INTEGER           :: i, j, worked,dime
    REAL(KIND=db)     :: fock_a(dime,dime), &
                         fock_b(dime,dime), &
                         over_a(dime,dime), &
                         over_b(dime,dime), &
                         eigc_a(dime,dime), &
                         eigc_b(dime,dime), &
                         eige_a(dime,dime), &
                         eige_b(dime,dime)
    REAL(KIND=db)     :: temp1(dime,dime),temp2(dime,dime),&
                         InvC(dime,dime)
    WRITE(*,*) 'Number of alpha / beta orbitals:', dime
    eige_a = 0
    eige_b = 0
    eigc_a = 0
    eigc_b = 0
    fock_a = 0
    fock_b = 0
    temp1  = 0
    temp2  = 0
    over_a = 0
    over_b = 0
OPEN(UNIT=111, FILE='coef.1',ACTION='WRITE')
OPEN(UNIT=112, FILE='coef.2',ACTION='WRITE')
OPEN(UNIT=222, FILE='eigenvalues.1.dat',ACTION='WRITE')
OPEN(UNIT=333, FILE='eigenvalues.2.dat',ACTION='WRITE')    
    OPEN(UNIT=100, FILE='alpha',ACTION='READ')
    OPEN(UNIT=200, FILE='beta',ACTION='READ')
      DO i = 1,3
        READ(100,'(1A)') r
        READ(200,'(1A)') r
      ENDDO
      DO i = 1,dime
        READ(100,*) r,r,r
        READ(r(12:31),'(1d20.14)') eige_a(i,i)
        READ(200,*) r,r,r
        READ(r(12:31),'(1d20.14)') eige_b(i,i)
        READ(100,'(4d20.14)') eigc_a(:,i)
        READ(200,'(4d20.14)') eigc_b(:,i)
!pb write coefs
      ENDDO


      DO i = 1,dime
        write(111,*) (eigc_a(i,j), j=1,dime)
      END DO
      DO i = 1,dime
        write(112,*) (eigc_b(i,j), j=1,dime)
      END DO
DO i = 1, dime
 write(222,*) (eige_a(i,j), j=1,dime)
 write(333,*) (eige_b(i,j), j=1,dime)
ENDDO

    CLOSE(222)
    CLOSE(333)
    CLOSE(100)
    CLOSE(200)

    Invc = inv(eigc_a) 
    CALL DGEMM('T','N',dime,dime,dime,1.0D0,Invc,dime,Invc,dime,0.0D0,over_a,dime)
    CALL DGEMM('N','N',dime,dime,dime,1.0D0,over_a,dime,eigc_a,dime,0.0D0,temp1,dime)
    CALL DGEMM('N','N',dime,dime,dime,1.0D0,eige_a,dime,Invc,dime,0.0D0,temp2,dime)
    CALL DGEMM('N','N',dime,dime,dime,1.0D0,temp1,dime,temp2,dime,0.0D0,fock_a,dime)
    Invc = inv(eigc_b) 
    CALL DGEMM('T','N',dime,dime,dime,1.0D0,Invc,dime,Invc,dime,0.0D0,over_b,dime)
    CALL DGEMM('N','N',dime,dime,dime,1.0D0,over_b,dime,eigc_b,dime,0.0D0,temp1,dime)
    CALL DGEMM('N','N',dime,dime,dime,1.0D0,eige_b,dime,Invc,dime,0.0D0,temp2,dime)
    CALL DGEMM('N','N',dime,dime,dime,1.0D0,temp1,dime,temp2,dime,0.0D0,fock_b,dime)
    OPEN(UNIT=100,FILE='hamiltonian.1',ACTION='WRITE')
    OPEN(UNIT=200,FILE='overlap.1',ACTION='WRITE')
    OPEN(UNIT=300,FILE='hamiltonian.2',ACTION='WRITE')
    OPEN(UNIT=400,FILE='overlap.2',ACTION='WRITE')
    DO i = 1, dime
      WRITE(100,'(4d20.10)') fock_a(:,i)
      WRITE(200,'(4d20.10)') over_a(:,i)
      WRITE(300,'(4d20.10)') fock_b(:,i)
      WRITE(400,'(4d20.10)') over_b(:,i)
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
