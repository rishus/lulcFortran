MODULE utilities_bfast
USE REAL_PRECISION

CONTAINS

SUBROUTINE cumsum(arr)
USE REAL_PRECISION
IMPLICIT NONE

   REAL(KIND=R8), INTENT(INout) :: arr(:)
   INTEGER :: n, i
   !REAL(KIND=R8), optional, INTENT(INOUT) :: cumsum

   n = SIZE(arr)
   IF (n .EQ. 0) RETURN
   DO i = 2, n
      arr(i) = arr(i-1) + arr(i)
   END DO

END SUBROUTINE cumsum

SUBROUTINE regress(t, y, model, fit, coeffs, K, work_arr)
USE globalVars_par   !only for the definition of bfastWA and numHarmonics?
IMPLICIT NONE
   !If model is harmonic, K MUST be supplied
   REAL(KIND=R8), INTENT(IN) :: t(:), y(:)
   CHARACTER(LEN=6), INTENT(IN) :: model
   INTEGER :: ncols, m, i, lwork, info, DeAllocateStatus, j
   REAL(KIND=R8), DIMENSION(:,:),  ALLOCATABLE :: A
   REAL(KIND=R8), DIMENSION(:),  ALLOCATABLE :: work
   INTEGER(kind=2), OPTIONAL, INTENT(IN) :: K
   REAL(KIND=R8), OPTIONAL, INTENT(INOUT) :: fit(:), coeffs(:)
   type(bfastWA), intent(inout) :: work_arr
   INTERFACE
      SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
      USE REAL_PRECISION
      CHARACTER:: TRANS
      INTEGER::   INFO, LDA, LDB, LWORK, M, N, NRHS
      REAL(KIND=R8):: A( LDA, * ), B( LDB, * ), WORK( * )
      END SUBROUTINE DGELS
   END INTERFACE

   m = SIZE(t)
   IF (model .EQ. "linear") THEN
      ncols = 2
   ELSE IF (model .EQ. "harmon") THEN
      ncols = 2*K + 1
   ELSE
      PRINT *, "MODEL = ", model
      PRINT *, "model not supplied"
      STOP
   ENDIF

   ALLOCATE (A(m, ncols))
   A(:,:) = 0
   A(:, 1) = (/ (1.0_R8, i = 1, m)  /)
   IF (model .EQ. "linear") THEN
      A(:, 2) = (/ (t(i), i = 1, m) /) 
   ELSEIF (model .EQ. "harmon") THEN
      DO j = 1, K
         A(:, 2*j) = (/ (COS(dble(j) * t(i)), i=1,m) /)
         A(:, 2*j+1) = (/ (SIN(dble(j) * t(i)), i=1,m) /)
      END DO
   ELSE
      PRINT *, "Model not supported"
      STOP
   ENDIF
  
   lwork = min(m, ncols) +  max(max(m, ncols), 1) * nbDGELS
   ALLOCATE (work(lwork))

   work_arr%tmp_mat(1:m, 1:ncols) = A
   work_arr%rtmp_vec1(1:m) = y
   CALL DGELS('N', m, ncols, 1, work_arr%tmp_mat, size(work_arr%tmp_mat, 1),  &
           &      work_arr%rtmp_vec1, size(work_arr%rtmp_vec1, 1), work, lwork, info)
   coeffs = work_arr%rtmp_vec1(1:ncols)
   DO i = 1, m
      fit(i) = DOT_PRODUCT(A(i, :), coeffs) 
   END DO   
   DEALLOCATE (A, work, stat = DeAllocateStatus)

END SUBROUTINE regress

SUBROUTINE OLSMosum(t, obs, process, model, h, K, work_arr)
USE globalVars_par
IMPLICIT NONE

   ! K is the degree of the fitting polynomial.
   ! For a linear fit, K = 1.
   ! For a trigonometric fit, K = number of harmonics we want to use.
   REAL(KIND=R8), INTENT(IN) :: t(:), obs(:), h
   CHARACTER(LEN=6), OPTIONAL, INTENT(IN) :: model
   INTEGER(kind=2), OPTIONAL, INTENT(IN) :: K
   INTEGER :: Sfinal, dofResid, i, nh
   REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: fit, residuals, coeffs
   REAL(KIND=R8) :: sigma
   REAL(KIND=R8), INTENT(INOUT) :: process(:)
   type(bfastWA), intent(inout) :: work_arr

   Sfinal = SIZE(t)

   ALLOCATE (fit(Sfinal), residuals(Sfinal))

   IF ( model .EQ. "linear") THEN
      ALLOCATE (coeffs(2))
      CALL regress(t, obs, "linear", fit, coeffs, K, work_arr)
   ELSEIF (model .EQ. "harmon") THEN
      IF (.NOT. PRESENT(K)) THEN
         PRINT *, "In SUBROUTINE OLSMosum. K MUST be supplied for harmonic model."
         PRINT *, "K is the number of harmonics we want to USE (and not the number of columns)."
         STOP
      ENDIF
      ALLOCATE (coeffs(2*K + 1))
      CALL regress(t, obs, "harmon", fit, coeffs, K, work_arr)
   ELSE
      PRINT *, "Model not supported"
   ENDIF
   residuals = obs - fit
   dofResid = Sfinal - 2   ! that's what it is for least squares fit. Ref., wikipedia
   !we're taking the sample standard deviation
   sigma = SQRT(SUM( (/ (residuals(i)*residuals(i), i=1,Sfinal) /) ) / dofResid)
   nh = FLOOR(Sfinal*h)
   CALL cumsum(residuals)  !the R code has a zero in front of residuals
   process = residuals(nh:Sfinal) - (/ 0.0_R8, residuals(1: Sfinal-nh) /)
   process = process / (sigma*SQRT(REAL(Sfinal, KIND=R8)))

END SUBROUTINE OLSMosum

SUBROUTINE criticalValueTables(method, critValTable, startRow, ENDRow)
USE globalVars_par 
IMPLICIT NONE

   CHARACTER(LEN=26) :: method
   INTEGER, INTENT(IN) :: startRow, ENDRow
   REAL(KIND=R8), INTENT(INOUT), allocatable :: critValTable(:,:)

   IF (method .EQ. "brownian bridge increments") THEN
      !Brownian bridge increments with "max" 
      critValTable = BBI(startRow:ENDRow,:)
   END IF
   !Brownian bridge with maxL2

   !Brownian bridge with meanL2 

END SUBROUTINE criticalValueTables

SUBROUTINE fillLine2(x, xEndpts, yEndpts, filledPoints, slopeOut)
IMPLICIT NONE
  ! Differs from fillLine only in the fact that this one 
  ! (1) evaluates only one segment at a time  
  ! (2) input x is the x-coordinates at which we want the evaluation done
  ! (3) output y is the computed y-coordinates cors to x
  REAL(KIND=R8), INTENT(IN):: x(:), xEndpts(:), yEndpts(:) !xEndpts, yEndpts are just 2 element arrays consisting of the starting and ENDing points of the desired line
  INTEGER :: i
  REAL(KIND=R8) :: slope
  REAL(KIND=R8), OPTIONAL, INTENT(OUT) :: slopeOut
  REAL(kind=R8), INTENT(INOUT) :: filledPoints(:)

  slope = (yEndpts(2) - yEndpts(1))/(xEndpts(2) - xEndpts(1))
  filledPoints = (/ (yEndpts(1) + slope * (x(i) - xEndpts(1)), i = 1, size(x))/)
  IF (PRESENT(slopeOut)) THEN
     slopeOut = slope
  ENDIF

END SUBROUTINE fillLine2

SUBROUTINE calcpValue(x, method, pval, k, h, fnalin)
IMPLICIT NONE
   ! x is the value at which the pValue has to be evaluated.
   ! method ... right now only ""brownian bridge increments" is supported
   ! pval is the output to be returned, again, a scalar
   ! k is the number of columns in the 'process'
   ! h is ... i don't know ... looks like sth to do with grid spacing. it's a
   ! parameter.
   ! fnalin ... right now only "max" is supported.
   REAL(KIND=R8), INTENT(IN) :: x, h
   INTEGER(kind=2) :: k
   CHARACTER(LEN=3), OPTIONAL, INTENT(IN) :: fnalin
   CHARACTER(LEN=3) :: fnal
   CHARACTER(LEN=26) :: method
   INTEGER :: start, ENDR, i, tableNcols
   REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: tableH, tableP, tableipl
   REAL(KIND=R8), DIMENSION(:,:), allocatable :: critValTable
   REAL(KIND=R8), INTENT(OUT) :: pval

   IF (.NOT. PRESENT(fnalin)) THEN
       fnal = "max"
   ELSE
       fnal = fnalin
   ENDIF
   ALLOCATE(critValTable(10,4), tableH(10), tableP(4+1), tableipl(4+1))  !why 10?
   IF (method .EQ. "brownian bridge increments") THEN
      IF (k > 6) THEN
         k = 6
      ENDIF
      start = (k-1)*10+1
      ENDR = k*10
      IF (fnal .EQ. "max") THEN
         CALL criticalValueTables("brownian bridge increments",critValTable, start, ENDR)
         tableNcols = SIZE(critValTable, 2)  !redundant as of now becuz we've allocated fixen no. of cols = 4 in/for the table.
         tableH = (/ (0.05*i, i=1,10) /)
         tableP = (/ 1.0_R8, 0.1_R8, 0.05_R8, 0.02_R8, 0.01_R8 /)
         tableipl = (/ (0, i=1,tableNcols+1) /)
         start = 1
         DO WHILE ((((tableH(start) <= h).and.(tableH(start+1) >= h)) .EQV.  .FALSE.) .and. (start < 10))
             start = start +1
         END DO
         ENDR = start + 1
         DO i = 1,tableNcols
            CALL fillLine2((/h /), tableH(start: ENDR), critValTable(start:endR,i), tableipl(i+1:i+1))
         END DO

         !Determine the interval in which x lies.
         IF ( x > tableipl(5)) THEN
            pval = tableP(5)
            RETURN
         ENDIF
         IF ( x < tableipl(1)) THEN
            pval = tableP(1)
            RETURN
         ENDIF
         start = 1
         DO WHILE ((((tableipl(start) <= x).and.(tableipl(start+1) >= x)) .EQV. .FALSE.) .and. (start < 5))
             start = start +1
         END DO
         ENDR = start + 1
         CALL fillLine2((/x/), tableipl(start:ENDR), tableP(start:endR), tableH(1:1))
         pval = tableH(1)
      ENDIF
   ENDIF
END SUBROUTINE calcpValue

SUBROUTINE recres(vec_timestamps, vec_obs, begin_idx, model, deg, & 
                &   work_arr, recResidual, Sfinal, ncols)
USE globalVars_par
IMPLICIT NONE

   REAL(KIND=R8), INTENT(IN) :: vec_timestamps(:), vec_obs(:)
   INTEGER, INTENT(IN) :: begin_idx, Sfinal, ncols
   INTEGER(kind=2), INTENT(IN) :: deg
   CHARACTER(LEN=6), INTENT(IN) :: model

   !INTEGER :: Sfinal, ncols, curr_idx, j, i, info !, rank
   INTEGER :: curr_idx, j, i, info !, rank
   REAL(KIND=R8), DIMENSION(Sfinal, ncols) :: mat_design
   REAL(KIND=R8) :: vec_fitobs(Sfinal), vec_fitcoefs(ncols)
   INTEGER, DIMENSION(ncols) :: ipiv
   REAL(KIND=R8) :: fr
   LOGICAL :: check

   REAL(KIND=R8), INTENT(INOUT) :: recResidual(:)
   type(bfastWA), intent(inout) :: work_arr 
   INTERFACE 
       SUBROUTINE DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       USE REAL_PRECISION
       INTEGER :: N, NRHS, LDA, IPIV(*), LDB, INFO
       REAL(KIND=R8) :: A(LDA, *), B(LDB, *) 
       END SUBROUTINE DGESV
   END INTERFACE

!   if (model == 'harmon')then
!     print *, 'in recres, Sfinal = ', Sfinal, 'model = ', model
!   endif
   IF (model == 'linear') THEN
      mat_design(:, 1) = 1
      mat_design(:, 2) = vec_timestamps
   ELSEIF (model == 'harmon') THEN
      mat_design(:, 1) = 1
      DO j = 1, deg
         mat_design(:, 2*j) = (/ (COS(dble(j) * vec_timestamps(i)), i=1,Sfinal) /)
         mat_design(:, 2*j+1) = (/ (SIN(dble(j) * vec_timestamps(i)), i=1,Sfinal) /)
      END DO
   ELSE 
      PRINT *, "model not supported"
      STOP
   ENDIF

   check = .TRUE.
   ! begin_idx is the firstmost index where the residual will be calculated.
   ! Begin_idx must be >=  ncols+1 
   ! The lastmost residual will be calculate for the index Sfinal.
   ! So residuals will be calculated from indices begin_idx to Sfinal.
   ! curr_idx is the point where the estimation window ends.
   ! So curr_idx will run from begin_idx-1 to Sfinal-1.
   DO curr_idx = begin_idx-1, Sfinal-1
      IF (check .eqv. .TRUE.) THEN
          ! the matrix XtX:  (X.t)_{ncols x curr_idx } * X_{curr_idx x ncols}
          work_arr%tmp_mat(1:ncols, 1:ncols) =   &
                    &  MATMUL(TRANSPOSE(mat_design(1:curr_idx,:)), mat_design(1:curr_idx,:))

          ! Solve (XtX) vec_fitcoefs = Xt vec_obs
          work_arr%rtmp_vec1(1:ncols) =  &
                  MATMUL(TRANSPOSE(mat_design(1:curr_idx, :)), vec_obs(1:curr_idx))
          CALL DGESV(ncols, 1, work_arr%tmp_mat(1:ncols, 1:ncols), ncols, ipiv, work_arr%rtmp_vec1,  &
                   &   size(work_arr%rtmp_vec1,1), info)
          vec_fitcoefs = work_arr%rtmp_vec1(1:ncols)

          ! Solve (XtX) y = xt
          work_arr%tmp_mat(1:ncols, 1:ncols) =   &
                    &  MATMUL(TRANSPOSE(mat_design(1:curr_idx,:)), mat_design(1:curr_idx,:))
          work_arr%rtmp_vec2(1:ncols) = mat_design(curr_idx+1, :)
          CALL DGESV(ncols, 1, work_arr%tmp_mat(1:ncols, 1:ncols), ncols, ipiv,  &
                  work_arr%rtmp_vec2(1:ncols), ncols, info)
          
          !rank = 2
!          IF (rank > 2) then
!             print *, 'rank > 2. Cannot proceed. Cross check recresid.'
!             STOP
!          ENDIF
!          IF (rank == 0) then
!             print *, 'rank = 0 in recres.'
!             STOP
!          ENDIF
      ENDIF

      fr = 1 + DOT_PRODUCT(mat_design(curr_idx+1, :), work_arr%rtmp_vec2(1:ncols))
      !if (model == 'harmon') then
!          print *, 'reached here in recres 5', curr_idx, 'fr=', fr
      !endif
      if (fr <0) then
         print *, 'fr < 0'
         print *, model, 'model'
         print *, 'info = ', info
         do j = 1,ncols
            print *, work_arr%tmp_mat(j, 1:ncols)
         end do
         print *, mat_design(curr_idx+1, :)
         print *, 'COEFFS =', work_arr%rtmp_vec2(1:ncols)
         stop
      endif
      recResidual(curr_idx+1) = (vec_obs(curr_idx+1) -  &
              DOT_PRODUCT(mat_design(curr_idx+1,:), vec_fitcoefs))/SQRT(fr)
   END DO

END SUBROUTINE recres 

SUBROUTINE buildDynPrTable(RSStri_full, Sfinal, numBrks, matPos, vecBrkPts)
USE globalVars_par
IMPLICIT NONE

   !numBrks is the number of internal breakpoints.
   !vecBrkPts is of dimension numBrks, i.e., t1 ans tS are NOT included in it.

   REAL(KIND=R8), intent(in) :: RSStri_full(:,:)
   INTEGER(kind=2), INTENT(IN) ::  numBrks, Sfinal
   INTEGER :: beginidx, ENDidx, potidxbegin, potidxend, idx,j, nbs, tmp(1)
   INTEGER :: brkpt_spacing
   INTEGER, INTENT(INOUT) :: matPos(:,:), vecBrkPts(:)
   REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE ::  matCost
   REAL(KIND=R8), DIMENSION(:), ALLOCATABLE ::  vecCost
   INTEGER:: last_brkpt_pos, curr_brkpt_pos, i

   brkpt_spacing = INT(FLOOR(Sfinal*bfast_brkptSpacingProp))
   ALLOCATE (matCost(numBrks, Sfinal))
   ALLOCATE(vecCost(Sfinal))
   matCost(:,:) = 0.0_R8
   matPos = -1 !no splits
   matCost(1,:) = RSStri_full(1,1:Sfinal) !return minimum cost index from here

   !make sure n > numBrks+1 * h
   DO nbs = 2,numBrks
      beginIdx = nbs*brkpt_spacing  !for nbs=2, h=10, beginIdx=20
      ENDIdx = Sfinal-brkpt_spacing      !                 endIdx = n-10 = 100-90 = 90 (say)
      potIdxBegin = (nbs-1)*brkpt_spacing   !for nbs=2,h=10, potIdxBegin=10
      DO idx = beginIdx,ENDIdx     !matCost gets filled at entries cors to idx colums.
         potIdxEnd = idx - brkpt_spacing       !potIdxEnd= idx-10. So vecCost gets filled from 10 to idx-10 for brk2, 20 to idx-10 for brk3, and so on.
         vecCost = 9999
         DO j = potIdxBegin, potIdxEnd           ! j=10 to idx-10
            vecCost(j) = matCost(nbs-1,j) + RSStri_full(j+1,idx)
         END DO
         matCost(nbs,idx) = MINVAL(vecCost)
         tmp = MINLOC(vecCost)
         matPos(nbs,idx) = tmp(1)
      END DO
   END DO
   !Note 1:  add the cost of segment(idx+1, Sfinal)
   do idx = numBrks*brkpt_spacing,Sfinal-brkpt_spacing
       matCost(numBrks, idx) = matCost(numBrks, idx) + RSStri_full(idx+1, Sfinal) !This is the last row-- what happens for numBrks == 1
   enddo
   tmp = MINLOC(matCost(numBrks,numBrks*brkpt_spacing:Sfinal-brkpt_spacing)) !the first break point
   last_brkpt_pos = tmp(1) + numBrks *brkpt_spacing -1
   curr_brkpt_pos = last_brkpt_pos
   vecBrkPts(numBrks) = last_brkpt_pos
   i = numBrks
   DO WHILE (i > 1)
      curr_brkpt_pos = matPos(i,curr_brkpt_pos)
      i = i-1
      vecBrkPts(i) = curr_brkpt_pos
   enddo

END SUBROUTINE buildDynPrTable

SUBROUTINE getBrkPts(vec_timestamps, vec_obs, model, deg, numBrks, vec_breakpoints, &
                work_arr, RSStri_full, Sfinal)
! work_arr is needed in recres.
! RSStri_full is needed in buildDynPrTable.
USE REAL_PRECISION
USE globalVars_par
IMPLICIT NONE

   !numBrks is the number of internal breakpoints. 
   !vec_breakpoints is of dimension numbrks,i.e., t1 & tS are NOT listed in it.
   REAL(KIND=R8), INTENT(IN) :: vec_timestamps(:), vec_obs(:)
   INTEGER(kind=2), INTENT(IN) :: deg, numBrks
   INTEGER, INTENT(IN) :: Sfinal
   CHARACTER(LEN=6) :: model 
 
   INTEGER :: i, j, brkpt_spacing 
   INTEGER :: ncols
   REAL(KIND=R8), DIMENSION(Sfinal) :: recResidual
   INTEGER, DIMENSION(numBrks+1,Sfinal) :: matPos
   INTEGER, INTENT(INOUT) :: vec_breakpoints(:)
   type (bfastWA), intent(inout) :: work_arr
   REAL(KIND=R8), INTENT(INOUT) :: RSStri_full(:,:)
  
   IF (model .eq. "linear") THEN
      ncols = 2
   ELSEIF (model .eq. "harmon") THEN
      ncols = INT(2*deg) + 1
   ELSE
      PRINT *, "model not supported"
   ENDIF

!   print *, 'in getBrkPts, num_valid_obs =', Sfinal, 'i.e.', SIZE(vec_obs,1), 'i.e.', SIZE(vec_timestamps,1)
!   print *, 'in getBrkPts, size of RSStri_full=', SIZE(RSStri_full, 1), SIZE(RSStri_full, 2)
   RSStri_full(1:Sfinal,1:Sfinal) = 0
   brkpt_spacing = INT(FLOOR(Sfinal*bfast_brkptSpacingProp))
!   print *, 'brkpt_spacing =', brkpt_spacing
   DO i = 1,Sfinal-brkpt_spacing+1 
      recResidual = 0 
      !print *, model, 'fit in (', i, ',', Sfinal,')'
      CALL recres(vec_timestamps(i:Sfinal), vec_obs(i:Sfinal), &
               &      ncols+1, model, deg, work_arr, recResidual(i:Sfinal),  &
               &      Sfinal-i+1, ncols )
      recResidual = (/ (recResidual(j)**2, j=1,SIZE(recResidual) )/)

      CALL cumsum(recResidual)
      RSStri_full(i, i:Sfinal) = recResidual(i:Sfinal)
   END DO

   CALL buildDynPrTable(RSStri_full, INT(Sfinal, kind=2), numBrks, matPos, vec_breakpoints)

END SUBROUTINE getBrkPts

SUBROUTINE hammingDistance(vec1, vec2, dist)
IMPLICIT NONE

   INTEGER, INTENT(IN) :: vec1(:), vec2(:)
   INTEGER :: i, len_vec1, len_vec2
   INTEGER, INTENT(OUT) :: dist  

   len_vec1 = SIZE(vec1, 1)
   len_vec2 = SIZE(vec2, 1)
   if ((len_vec1 == 0) .and. (len_vec2 == 0)) then
      dist = 0
      return 
   elseif (len_vec1 == 0) then
      dist = len_vec2
      return
   elseif (len_vec2 == 0) then
      dist = len_vec1
      return
   endif

   dist = 0
   if (len_vec1 > len_vec2) then
       do i = 1,len_vec2 
          if (vec1(i) /= vec2(i)) then
             dist = dist + 1
          endif
        end do
   elseif (len_vec2 > len_vec1) then
        do i = 1, len_vec1
           if (vec1(i) /= vec2(i)) then
              dist = dist + 1
           endif
        end do
   else
        do i = 1, len_vec1
           if (vec1(i) /= vec2(i)) then
              dist = dist + 1
           endif
        end do
   endif
   
   dist = dist + abs(len_vec1 - len_vec2)

   return

END SUBROUTINE hammingDistance

END MODULE utilities_sept16_par
