MODULE utilities_ewma
CONTAINS

FUNCTION std(data, work_arr) RESULT(sigma)
USE globalVars_par
USE REAL_PRECISION
IMPLICIT NONE

!
! std CALCULATES THE STANDARD DEVIATION FOR THE DATA IN VARIALBLE data
!
   REAL(KIND=R8), INTENT(IN) :: data(:)
   INTEGER  :: m, i
   REAL(KIND=R8) :: mean, var, sigma
   type(ewmaWA), intent(inout) :: work_arr

   m = size(data, 1)
   if (m==0) then
      print *, 'm=0'
   endif
   if (m==1) then
      sigma = 0
      return
   endif
   
   mean = (SUM(data(1:m)))/dble(m)
   work_arr%dev(1:m)  =  data(1:m) - mean
   work_arr%dev(1:m) = (/ (work_arr%dev(i)*work_arr%dev(i),  i=1,m) /)
   var = SUM(work_arr%dev(1:m))/dble(m-1)  ! because R uses N-1. Otherwise, it would be divided by N only.
   sigma = SQRT(var)

RETURN
END FUNCTION std

SUBROUTINE lsfit(t, u, K, alphaStar, IERR, work_arr)
!
! lsfit CALCULATES THE LEAST SQUARES FIT FOR THE DATA GIVEN BY VECTORS T AND U
!
! Input parameters:
!      t:            x-coordinates
!      u:            y-coordinates (observations)
!      K:            number of harmonics
!
! Output parameters:
!      alphaStar:    2K+1 harmonic regression coefficients after refining
!      IERR:         
!---------------------------------------------------------------
USE globalVars_par
USE REAL_PRECISION
IMPLICIT NONE

   INTEGER(kind=2), INTENT(IN) :: K
   REAL(KIND=R8), INTENT(IN) :: t(:), u(:)
   REAL(KIND=R8), intent(inout) :: alphaStar(:)
   type(ewmaWA), intent(inout) :: work_arr
   INTEGER             :: i, j, lwork, DeAllocateStatus
   INTEGER(kind=2)     :: ccol, ncols
   INTEGER             :: info, IERR, sz, M, mn
   REAL(KIND=R8)       :: sigma, tau1, detR1, detXtX
   REAL(KIND=R8), DIMENSION(size(t,1),2*K+1) :: X
   REAL(KIND=R8), DIMENSION(:),  ALLOCATABLE :: work
   ! the variables below may be useful if true conditioning of the matrix
   ! is to be tested. Determinant is not quite the correct indicatior
   ! of matrix conditioning.
   !INTEGER             :: svd_lwork, svd_liwork, nlvl, SMLSIZ, rank
   !REAL(KIND=R8)       :: rcond
   !REAL(kind=r8), DIMENSION(:,:), ALLOCATABLE :: svd_tmp_mat
   !REAL(kind=r8), DIMENSION(:), ALLOCATABLE   :: svd_work, svd_iwork 
   !INTEGER :: ILAENV
   !EXTERNAL ILAENV
   INTERFACE
      SUBROUTINE DGELS(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
      USE REAL_PRECISION
      CHARACTER:: TRANS
      INTEGER::   INFO, LDA, LDB, LWORK, M, N, NRHS
      REAL(KIND=R8):: A( LDA, * ), B( LDB, * ), WORK( * )
      END SUBROUTINE DGELS
   END INTERFACE
   INTERFACE
      SUBROUTINE dgeqrf(M, N, A, LDA, TAU, WORK, LWORK, INFO)
      USE REAL_PRECISION
      INTEGER :: M, N, LDA, lwork, info
      REAL(KIND=R8) :: A(LDA, *), tau(*), work(*)
      END SUBROUTINE dgeqrf 
   END INTERFACE
   
   M = size(t, 1)

   ncols = 2*K + 1

   X(:,1) = (/ (1, i = 1,M) /)
   DO j = 1,K
     ccol = K+j+1
     X(:, j+1) = (/ (SIN(dble(j) * t(i)), i=1,M) /)
     X(:, ccol) = (/ (COS(dble(j) * t(i)), i=1,M) /)
   ENDDO

   !(1) checking the determinant (conditioning!) of the matrix
   !was originally designed for dgels
   lwork = min(M, int(ncols, kind=4)) +  max(max(M, INT(ncols, KIND=4)), 1) * nb  
   ALLOCATE (work(lwork))  !TODO: cross check if dgeqrf needs a different value
   work_arr%tmp_mat(1:M, 1:ncols) = X
   mn = min(M, int(ncols,kind=4))
   CALL dgeqrf(M, INT(ncols, KIND=4), work_arr%tmp_mat(1:M, 1:ncols), M, &
               work_arr%rtmp_vec1(1:mn),  work(1:lwork), lwork, info)
   ! determinant of R.
   detR1 = 1
   DO i = 1,mn   ! len(u) > ncols, so mn = min(M,ncols) = ncols, for us.
      detR1 = detR1 * work_arr%tmp_mat(i,i)
   END DO
   ! X = QR = [Q1, Q2]*[R1; 0] = Q1R1, Q, Q1, Q2 orthogonal. So, det(XtX) = det(R1tR1).  
   ! R1 ix ncols x ncols. So det(R1tR1) = det(R1)*det(R1).
   detXtX = detR1 * detR1
   IF (abs(detXtX) < 0.001) THEN 
        alphaStar = (/ (-2222, i=1,ncols) /)
        IERR = 1
        DEALLOCATE (work, stat = DeAllocateStatus)
        RETURN 
   ENDIF

   !(3) uses QR factorization to solve the system, way more stable
   !    drawback of DGELS: It is assumed that A has full rank. 
   work_arr%tmp_mat(1:M, 1:ncols) = X
   work_arr%rtmp_vec1(1:M) = u
   CALL DGELS('N', M, int(ncols, kind=4), 1, work_arr%tmp_mat(1:M, 1:ncols), &
            &   M, work_arr%rtmp_vec1(1:M), M, work(1:lwork), lwork, info)
   alphaStar = work_arr%rtmp_vec1(1:ncols)   ! alpha

   work_arr%rtmp_vec1(1:M) = MATMUL(X, alphaStar)   !utilde
   work_arr%rtmp_vec2(1:M) = u - work_arr%rtmp_vec1(1:M)     !Ealpha
   IF (MAXVAL(work_arr%rtmp_vec2(1:M)) > 10.0E6  .OR. MINVAL(work_arr%rtmp_vec2(1:M)) < -10.0E6) THEN 
        alphaStar = (/ (-2222, i=1, ncols) /)
        IERR = 1
        DEALLOCATE (work, stat = DeAllocateStatus)
        RETURN
   ENDIF
   sigma = std(work_arr%rtmp_vec2(1:M), work_arr)

   tau1 = xbarlim1 * sigma
   sz = 0
   do i=1,m
      if (ABS(work_arr%rtmp_vec2(i)) < tau1) then
         sz = sz + 1
         work_arr%itmp_vec1(sz) = i                 ! I
      endif
   enddo
   IF (sz <= ncols) THEN
        alphaStar = (/ (-2222, i=1, ncols) /)
        IERR = 1
        DEALLOCATE (work, stat = DeAllocateStatus)
        RETURN
   ENDIF

   work_arr%tmp_mat(1:sz, :) = X(work_arr%itmp_vec1(1:sz), :)    ! Xsub
   work_arr%tmp_mat(1:sz,1:ncols) = X(work_arr%itmp_vec1(1:sz), :)    ! Xsub
   work_arr%rtmp_vec1(1:sz) = u(work_arr%itmp_vec1(1:sz))
   CALL DGELS('N', sz, int(ncols, kind=4), 1, work_arr%tmp_mat(1:sz,1:ncols), sz ,   &
             &   work_arr%rtmp_vec1(1:sz), sz, work(1:lwork), lwork, info)
   alphaStar = work_arr%rtmp_vec1(1:ncols)
  
   DEALLOCATE (work, stat = DeAllocateStatus)

END SUBROUTINE LSFIT


SUBROUTINE getResiduals(alphaStar, t, D, K, M, ierr_lsfit, & 
                      sigmaIhat, EstarAlphastar, Ibar, sz, nullFlag, work_arr)
!
!Residuals calculation for Step 1 of the pseudocode.
!    Inputs:
!        alphaStar:         refined coefficients determined with repeated regression
!        t:                 the vector of present timepoints
!        D:                 the vector of present observations
!        K:                 number of harmonics used
!        M:                 length of training period
!        ierr_lsfit:        error flag (takes value 0 or 1) from lsfit.
!                           1 if matrix determinant was too small, or, 
!                             if outlier free training data is too small less, or,
!                             if error vector Ealpha had any value of too large magnitude.
!                           0 otherwise.
!        xbarlimit1:        parameter to determine 'outliers' in the initial 
!                           phase of determining 
!                           refined coefficients.
!        xbarlimit2:        Parameter to calculate the threshold to decide outliers
!                           when refined coefficients are being used.
!        
!    Outputs:
!        IBar:              array
!                           of indices of outlier-free time series
!        sigmaIhat:         scalar.
!                           standard deviation of outlier-free training data errors
!        nullFlag:          Booloean.
!                           True if further processing cannot be carried out. 
!                           (E.g., if outlier-free-training data is of length
!                           less than 2K+1, or, if ierr_lsfit was 1.).
!                           False otherwise.
!        EstarAlphastar:    array of residuals (errors) for the entire time series
!        sz:                Only the first sz elements of Ibar are to be used.
!                           (Ibar was allocated in the main code to avoid
!                           local allocations etc.)
    
USE globalVars_par
USE REAL_PRECISION
IMPLICIT NONE
   INTEGER, INTENT(IN) :: M, ierr_lsfit
   INTEGER(kind=2), INTENT(IN) :: K
   REAL(KIND=R8), INTENT(IN) :: t(:), D(:), alphaStar(:)
   INTEGER, INTENT(INOUT) ::  IBar(:)
   REAL(KIND=R8), INTENT(INOUT) ::  EstarAlphastar(:)
   type(ewmaWA), intent(inout) :: work_arr
   LOGICAL, INTENT(INOUT)     ::  nullFlag
   REAL(KIND=R8), INTENT(OUT)    ::  sigmaIhat
   INTEGER, INTENT(OUT) ::  sz
   REAL(KIND=R8)     :: sigma2, tau1, tau2 
   INTEGER     :: ncols, i, j, ccol, szhat, S

   IF (ierr_lsfit == 1) then
       nullFlag = .TRUE.
       sigmaIhat = 0
       RETURN
   ENDIF

   S = size(t, 1)
   ncols = 2*K + 1
   work_arr%tmp_mat(:,1) = (/ (1, i=1,S) /)
   DO j = 1,K
     ccol = K+j+1
     work_arr%tmp_mat(1:S, j+1) = (/ (SIN(dble(j) * t(i)), i=1,S) /)
     work_arr%tmp_mat(1:S, ccol) = (/ (COS(dble(j) * t(i)), i=1,S) /)
   ENDDO

   ! EstarAlphastar  ! to be returned, will be needed later 
   work_arr%rtmp_vec1(1:S) = D - MATMUL(work_arr%tmp_mat(1:S,:), alphaStar)  
   if (MAXVAL(work_arr%rtmp_vec1(1:M)) > 10.0E6 .OR. &
           MINVAL(work_arr%rtmp_vec1(1:M)) < -10.0E6) THEN
       nullFlag = .TRUE.
       sigmaIhat = 0
       RETURN
   endif

   EstarAlphastar(1:S) = work_arr%rtmp_vec1(1:S)

   if (M == 0) then
      print *, 'M = ', 0 
   endif

   sigma2 = std(work_arr%rtmp_vec1(1:M), work_arr)
   tau1 = xbarlim1 * sigma2
   tau2 = xbarlim2 * sigma2

   sz = 0
   DO i = 1,M
      IF ((abs(work_arr%rtmp_vec1(i)) < tau1) .and. (D(i) > ewmacd_lowthresh)) then
         sz = sz + 1
         work_arr%itmp_vec1(sz) = i   ! Ibar       
         work_arr%itmp_vec2(sz) = i   ! Ihat
      ENDIF
   ENDDO
   szhat = sz
   IF (szhat <= ncols) THEN
      nullFlag = .TRUE.
      sigmaIhat = 0
      RETURN
   ENDIF
   DO i = M+1, S
      IF ((abs(work_arr%rtmp_vec1(i)) < tau2) .and. (D(i) > ewmacd_lowthresh)) then
          sz = sz + 1
          work_arr%itmp_vec1(sz) = i  ! Ibar continued    ! to be returned, will be needed later 
      ENDIF
   ENDDO

   Ibar(1:sz) = work_arr%itmp_vec1(1:sz)

   work_arr%rtmp_vec2(1:szhat) = (/ (work_arr%rtmp_vec1(work_arr%itmp_vec2(i)), i=1,szhat) /)   !Estar_alphastar_Ihat 
   sigmaIhat = std(work_arr%rtmp_vec2(1:szhat), work_arr)    ! to be returned, will be needed later 

END SUBROUTINE getResiduals


SUBROUTINE getControlLimits(sigma, len_Ibar, tau)
!
!    Step 2 of the pseudocode.
!    
!    A control chart is built.
!    Inputs:
!        sigma, len_Ibar:  Statistical information about the 
!                          'outlier-processed' timeseries.
!
!    Output:
!        tau:    control chart limits vector.

USE globalVars_par
USE REAL_PRECISION
IMPLICIT NONE
   REAL(KIND=R8), INTENT(IN)  :: sigma
   INTEGER, INTENT(IN)  :: len_Ibar
   REAL(KIND=R8), INTENT(INOUT) :: tau(:)
   INTEGER   ::   i
   REAL(KIND=R8)   :: sl, f, a, p, b
  
   sl = sigma*L
   f = lam/(2.0 - lam)
   a = 1 - lam
   DO i=1,len_Ibar
      p = a**(2*i)
      b = 1.0 - p
      tau(i) = mu + sl * SQRT(f * b)
   ENDDO

END SUBROUTINE getControlLimits


SUBROUTINE getEWMA(Ibar, EstarAlphastar, z, len_Ibar)
!
!    Step 3 of the pseudocode.
!    
!    EWMA of residuals resulting from harmonic regression is calculated.
!
!    Inputs:
!        Ibar:               indices of the outlier free timeseries
!        Estar_alphastar:    numpy array of the residuals of the time series
!        len_Ibar:           length of Ibar
!
!    Output:
!        z:                  EWMA of residuals resulting from harmonic regression

USE globalVars_par
USE REAL_PRECISION
IMPLICIT NONE

  INTEGER, INTENT(IN) :: Ibar(:), len_Ibar
  REAL(KIND=R8), INTENT(IN) :: EstarAlphastar(:)
  REAL(KIND=R8), INTENT(INOUT) :: z(:)
  INTEGER :: i
  REAL(KIND=R8) ::  a

  !if(size(Ibar) == 0) then    !This case won't arise because size(Ihat) > 2*K+1
  if(len_Ibar == 0) then
    print *, "Invalid Ibar"
    STOP
  endif

  !print *, "In SUBROUTINE getEWMA"
  
  z = (/ (0, i=1,len_Ibar)  /)
  z(1) = EstarAlphastar(Ibar(1))
  
  a = 1.0 - lam
  DO i=2, len_Ibar
     z(i) = a * z(i-1) + lam * EstarAlphastar(Ibar(i))
  ENDDO

END SUBROUTINE getEWMA

SUBROUTINE flagHistory(z, tau, f, len_Ibar)
!
!    Step 4, flag history calculation of the pseudocode.
!    
!    Flag history is claculated here, based on the 'persistence' of any nonzero signals.
!    
!    Inputs:
!        z:         EWMA of residuals resulting from harmonic regression
!        tau:       Control chart vector
!        len_Ibar:  lenght of Ibar (i.e., number of indices of outlier free timeseries)
!
!    Output:
!        f:     flag history

USE REAL_PRECISION
IMPLICIT NONE

  INTEGER   ::  i, tmp
  REAL(KIND=R8), DIMENSION(:), INTENT(IN) :: z, tau
  INTEGER, INTENT(IN) :: len_Ibar
  REAL(KIND=R8), INTENT(INOUT) :: f(:)

  f = (/ (0, i=1,len_Ibar) /)
  DO i=1,len_Ibar
     IF (z(i) /= 0) THEN
        tmp = SIGN(1, z(i))
     ELSE
        tmp = 0
     ENDIF
     !f(i) = INT(tmp * FLOOR(ABS(z(i)/tau(i))))
     f(i) = tmp * (ABS(z(i)/tau(i)))
  END DO
     !Can we not vectorize this?:
     !f = (/ (SIGN(1,z(i)) * (ABS(z(i)/tau(i))),  i=1,len_Ibar) /)
     !The problem with this is SIGN(1,0) is return as 1, not 0. We want 0.
 
END SUBROUTINE flagHistory

SUBROUTINE persistenceCounting(f, persistenceVec, work_arr, len_Ibar)
!
!    Step 4, persistence calculation of the pseudocode
!    Inputs:
!        f:              the flag history
!        len_Ibar:       no. of points in outlier free time series
!        
!    Output:
!        persistenceVec: list of integers
!                        each element indicates how long the current sign (+ or -) lasted.

USE REAL_PRECISION
USE globalVars_par
IMPLICIT NONE

  REAL(KIND=R8), INTENT(IN) :: f(:)
  INTEGER, INTENT(IN) :: len_Ibar
  REAL(KIND=R8), INTENT(INOUT)  ::  persistenceVec(:)
  type(ewmaWA), intent(inout) :: work_arr
  INTEGER :: imax, i, tmp3hi, tmp3lo  !,mx
  REAL(KIND=R8) :: mx 
  
  ! sign of f
  DO i=1,len_Ibar
     IF (f(i) > 0) THEN
        work_arr%itmp_vec1(i) = 1      ! f_sgn
     ELSEIF (f(i) < 0) THEN
        work_arr%itmp_vec1(i) = -1
     ELSE
        work_arr%itmp_vec1(i) = 0 
     ENDIF
  ENDDO

  ! Count consecutive dates in which directions are sustained
  work_arr%itmp_vec2(1:len_Ibar) = (/ (0, i=1,len_Ibar)  /)  
  DO i = 1,len_Ibar
    tmp3lo = 0
    tmp3hi = 0
    DO WHILE (i - tmp3lo >  0)
       IF (work_arr%itmp_vec1(i) - work_arr%itmp_vec1(i - tmp3lo) == 0) THEN
          tmp3lo = tmp3lo + 1
       ELSE
          EXIT 
       ENDIF
    ENDDO
    DO WHILE (i + tmp3hi <= len_Ibar)
        IF (work_arr%itmp_vec1(i + tmp3hi) - work_arr%itmp_vec1(i) == 0) THEN
           tmp3hi = tmp3hi + 1
        ELSE
           EXIT
        ENDIF
    ENDDO
    work_arr%itmp_vec2(i) = tmp3lo + tmp3hi - 1
  ENDDO

  persistenceVec = (/ (0, i=1,len_Ibar)  /)  ! tmp4
  mx = 0
  imax = -1
  DO i=1,len_Ibar
     IF (work_arr%itmp_vec2(i) >= persistence) THEN
        persistenceVec(i) = f(i)
        mx = f(i)
        imax = i
     ELSE
        IF (imax == -1) THEN
           persistenceVec(i) = 0
        ELSE
           persistenceVec(i) = mx
        ENDIF
     ENDIF
  ENDDO

END SUBROUTINE persistenceCounting

SUBROUTINE summarize(pixel_x, pixel_y, jumpValsOrgSten, presInd, method, Sfinal, S)
!
!    Summarizing.
!    
!    Missing data points are addressed here.
!    Breakpoints are determined here.
!    
!    Inputs:
!        pixel_x:        global index
!        pixel_y:
!        jumpValsOrgSten:
!        presInd:
!        method:         string ('annual_mean' or 'on-the-fly')
!                        For getting breakpoints, choose 'on-the-fly'
!    
!    Output:
!        ewma_summary:   global 3D array
!                        positions (pixel_x, pixel_y, :) are filled up
        
USE REAL_PRECISION
USE globalVars_par
IMPLICIT NONE
   
  INTEGER(KIND=8), INTENT(IN) :: pixel_x, pixel_y
  INTEGER(KIND=2), INTENT(IN) :: jumpValsOrgSten(:), presInd(:)
  INTEGER, INTENT(IN) ::  method, S
  INTEGER :: i, tt, Sfinal

  tt = 1
  DO i = 1, Sfinal
     ewma_summary(presInd(i), pixel_x, pixel_y) = REAL(jumpValsOrgSten(tt), KIND= 4)
     tt = tt + 1
  END DO                         ! missing and 'outlier' timestamps still have -2222

  IF (ewma_summary(1, pixel_x, pixel_y) == -2222) THEN  ! if the very first obs was missing/outlier, set it to zero
     ewma_summary(1, pixel_x, pixel_y) = 0 
  ENDIF

  IF (ewma_summary(1, pixel_x, pixel_y).ne.(-2222)) THEN   ! replace the -2222s with last available ewma value
     DO i = 2,S
       IF (ewma_summary(i, pixel_x, pixel_y) == -2222) THEN
          ewma_summary(i, pixel_x, pixel_y) = ewma_summary(i-1, pixel_x, pixel_y)
       ENDIF
     END DO
  ENDIF

  i = method

END SUBROUTINE summarize

SUBROUTINE summarize_residuals(pixel_x, pixel_y, residualsPresentSten, presInd, method, &
                  &             Sfinal)
USE REAL_PRECISION
USE globalVars_par
IMPLICIT NONE
  
  INTEGER(KIND=8) :: pixel_x, pixel_y
  REAL(KIND = 4), INTENT(IN) :: residualsPresentSten(:)
  INTEGER(KIND=2), INTENT(IN) :: presInd(:)
  INTEGER, INTENT(IN) ::  method, Sfinal
  INTEGER :: i, tt

  tt = 1
  DO i = 1, Sfinal
     ewma_residuals(presInd(i), pixel_x, pixel_y) = residualsPresentSten(tt)
     tt = tt + 1
  END DO                         ! missing and 'outlier' timestamps will have -2222

  i = method

END SUBROUTINE summarize_residuals

END MODULE utilities_ewma
