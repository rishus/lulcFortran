MODULE bfast_fn_par
CONTAINS

SUBROUTINE bfast(imgIndex)
use globalVars_par
use utilities_bfast
use real_precision
implicit none

! ************* input vars **********************
INTEGER(kind = 8), intent(in) :: imgIndex(:)
! ***********************************************

! ************* utility variables ********************
!general arrays and other variables
INTEGER(KIND=8) :: num_pixels, pixel, pixel_x, pixel_y  !just like ewmacd
INTEGER :: i, ncols, len_training_vec, num_valid_obs, pctr
REAL(KIND=R8), DIMENSION(num_obs) ::  vec_timestamps    ! wasn't needed in ewmacd 
REAL(KIND=R8), DIMENSION(num_obs) ::  vec_timestamps_pres, vec_obs_pres  !t and D of ewmacd
REAL(KIND=R8), DIMENSION(num_obs) :: vec_u, vec_util
REAL(KIND=R8), DIMENSION(num_obs) :: vec_fit_pres
REAL(KIND=R8), DIMENSION(num_obs, num_obs) :: RSStri_full
integer (kind=2), dimension(num_obs) :: presInd
type(bfastWA) :: work_arr

!arrays for linear part
INTEGER, DIMENSION(bfast_numBrks+2) :: vec_trend_brkpts, vec_trend_brkpts_old
REAL(KIND=R8), DIMENSION(bfast_numBrks+1, 2) ::  linear_fit_coeffs
 
!arrays for seasonal part
REAL(KIND=R8), DIMENSION(num_obs) ::  vec_timestamps_seasonal
REAL(KIND=R8), DIMENSION(num_obs) ::  vec_timestamps_seasonal_pres
INTEGER, DIMENSION(bfast_numBrks+2) :: vec_seasonal_brkpts, vec_seasonal_brkpts_old 
REAL(KIND=R8), DIMENSION(2*bfast_numHarmonics + 1) ::  seasonal_fit_coeffs

!arrays for 'work'
REAL(KIND=R8), DIMENSION(num_obs) ::  process
INTEGER :: len_process

!other variables (parameters)
INTEGER(kind=2)  ::  numTrendBrks, numSeasonalBrks, numTrendSegs, numSeasonalSegs
INTEGER(kind=2)  ::  numTrendBrksFinal, numSeasonalBrksFinal 
REAL(KIND = R8) :: tauTrend, tauSeason, pvalTrend, pvalSeason
INTEGER  :: iter, brkpt_spacing
INTEGER :: startPoint, endPoint, j, hamTrend, hamSeason
! ***********************************************

! ********* PURPOSE ****************
!     Code for BFAST algorithm
!     Supporting routines are all in utilities_bfast.f90
! **********************************

numTrendBrks = bfast_numBrks ! this is internal number of breaks.If we do count the boundary pts
numSeasonalBrks = bfast_numBrks ! as breaks, then there are a total of numbrks + 2 breaks /
tauTrend = bfast_pvalThresh
tauSeason = bfast_pvalThresh

ALLOCATE (work_arr%tmp_mat(num_obs, 2*bfast_numHarmonics+1),  &
          work_arr%rtmp_vec1(num_obs), work_arr%rtmp_vec2(num_obs),  &
          work_arr%itmp_vec1(num_obs), work_arr%itmp_vec2(num_obs))

num_pixels = SIZE(imgIndex, 1)
ncols = 2*bfast_numHarmonics + 1

DO pixel = 1, num_pixels

     pixel_x = mod(imgIndex(pixel), NumCols)
     IF (pixel_x == 0) THEN
         pixel_x = NumCols
     ENDIF
     pixel_y = INT(CEILING(REAL(imgIndex(pixel))/NumCols))

     ! ********* develop the presInd vector *******************
     pctr = 0
     presInd = 0
     len_training_vec = 0
     DO i = 1, num_obs
        IF (input_mat(i, pixel_x, pixel_y) > ewmacd_lowthresh) THEN
           pctr = pctr + 1
           presInd(pctr) = int(i, kind=2)
           IF ((tyeardoy(i,1) >= trainingStart) .and. (tyeardoy(i,1) < trainingEnd)) THEN
              len_training_vec = len_training_vec + 1
           ENDIF
        ENDIF
     END DO
     num_valid_obs = pctr

     IF (len_training_vec <= 2.0 * ewmacd_numHarmonics+1 ) THEN
        CYCLE     !note the use of ewmacd harmonics here
     ENDIF
    ! **************** prepare data *************************
    
    brkpt_spacing = INT (FLOOR (num_valid_obs * bfast_brkptSpacingProp))

    vec_timestamps = (/ (tyeardoy(i,1) + (REAL(tyeardoy(i,2),kind=r8)/365.0_R8), i=1,num_obs)  /)  
    vec_timestamps_pres(1:num_valid_obs) = (/ (vec_timestamps(presInd(i)), i=1,num_valid_obs) /)
    vec_timestamps_seasonal = (/ (2.0_R8 * Pi * REAL(tyeardoy(i,2),KIND=R8)/365.0_R8, i=1,num_obs) /)
    vec_timestamps_seasonal_pres(1:num_valid_obs) = (/(vec_timestamps_seasonal(presInd(i)), i=1,num_valid_obs) /)
    vec_obs_pres(1:num_valid_obs) = (/ (input_mat(presInd(i), pixel_x, pixel_y), i = 1,num_valid_obs) /)

    ! ********************** actual processing starts here ***********
    !initialize
    vec_fit_pres = 0

    len_process = num_valid_obs - brkpt_spacing + 1

    vec_fit_pres = 0
    vec_trend_brkpts = 0
    vec_seasonal_brkpts = 0
    vec_trend_brkpts_old = 1
    vec_seasonal_brkpts_old = 1
    CALL hammingDistance(vec_trend_brkpts, vec_trend_brkpts_old, hamTrend) 
    iter = 0

    DO WHILE (((hamTrend /= 0).or.(hamSeason /= 0)) .and.  (iter < bfast_maxIter))
         vec_trend_brkpts_old = vec_trend_brkpts
         vec_seasonal_brkpts_old = vec_seasonal_brkpts
      
         !'adjust' the data by deseasoning
         vec_u(1:num_valid_obs) = vec_obs_pres(1:num_valid_obs)- vec_fit_pres(num_valid_obs)
      
         !get OLS-Mosum statistics
         CALL OLSMosum(vec_timestamps_pres(1:num_valid_obs), vec_u(1:num_valid_obs),   &
                       process(1:len_process), "linear", bfast_brkptSpacingProp, int(1, kind=2),   &
                       work_arr) 
         CALL calcpValue(MAXVAL(ABS(process(1:len_process))), "brownian bridge increments",   &
                &        pvalTrend, bfast_numColsProcess, bfast_brkptSpacingProp, "max")
       
         !get breakpoints in trend
         if (pvalTrend <= tauTrend) then
             numTrendBrksFinal = numTrendBrks + 2
             vec_trend_brkpts(1) = 1
             vec_trend_brkpts(numTrendBrks+2) = num_valid_obs
             ! numBrks argument in getBrkPts is the internal number of breakpoints
             
             CALL getBrkPts(vec_timestamps_pres(1:num_valid_obs), vec_u(1:num_valid_obs),   &
                     "linear", INT(1, KIND=2), numTrendBrks,    &
                     vec_trend_brkpts(2:numTrendBrks+1), work_arr,  &
                     RSStri_full(1:num_valid_obs, 1:num_valid_obs), num_valid_obs)
             numTrendSegs = numTrendBrks + 1 !which is = numTrendBrksFinal-1; numTrendBrks is the internal brk pts.
         else
             numTrendBrksFinal = 2
             vec_trend_brkpts(1) = 1
             vec_trend_brkpts(2) = num_valid_obs
             numTrendSegs = 1
         endif
         !do piecewise linear approximation
         DO i = 1, numTrendSegs
            startPoint = vec_trend_brkpts(i)
            endPoint = vec_trend_brkpts(i+1)
            IF ( i == numTrendSegs) THEN
               CALL regress(vec_timestamps_pres(startPoint:endPoint), &
                      &     vec_u(startPoint:endPoint),"linear",   &
                      &     vec_fit_pres(startPoint:endPoint), linear_fit_coeffs(i, 1:2), & 
                      &     int(0,kind=2), work_arr)
            ELSE
               CALL regress(vec_timestamps_pres(startPoint:endPoint-1),   &
                      &     vec_u(startPoint:endPoint-1),"linear",   &
                      &     vec_fit_pres(startPoint:endPoint-1), linear_fit_coeffs(i, 1:2), &
                      &     int(0, kind=2), work_arr)
           ENDIF 
         END DO 
      
         !adjust the data by detrending (use the most recently calculated, i.e., updated, trend model)
         vec_util(1:num_valid_obs) = vec_obs_pres(1:num_valid_obs) - vec_fit_pres(1:num_valid_obs)
      
         !get OLS-Mosum statistics for seasonal 
         CALL OLSMosum(vec_timestamps_seasonal_pres(1:num_valid_obs),   &
                       vec_util(1:num_valid_obs), process(1:len_process), "harmon",  & 
                 &   bfast_brkptSpacingProp, bfast_numHarmonics, work_arr) 
         CALL calcpValue(MAXVAL(ABS(process(1:len_process))), "brownian bridge increments", &
                         &  pvalSeason, bfast_numColsProcess, bfast_brkptSpacingProp, "max")
         
         !get breakpoints in season
         if (pvalSeason <= tauSeason) then
             numSeasonalBrksFinal = numSeasonalBrks + 2
             vec_seasonal_brkpts(1) = 1
             vec_seasonal_brkpts(numSeasonalBrks+2) = num_valid_obs
             CALL getBrkPts(vec_timestamps_seasonal_pres(1:num_valid_obs),   &
                            vec_util(1:num_valid_obs), "harmon", &
                            bfast_numHarmonics, numSeasonalBrks,   &
                            vec_seasonal_brkpts(2:numSeasonalBrks+1),  &
                            work_arr, RSStri_full(1:num_valid_obs, 1:num_valid_obs), &
                            num_valid_obs) 
             numSeasonalSegs = numSeasonalBrks + 1  ! internal brks + 1
         else
             numSeasonalBrksFinal = 2
             vec_seasonal_brkpts(1) = 1
             vec_seasonal_brkpts(2) = num_valid_obs
             numSeasonalSegs = 1
         endif
         !do piecewise harmonic approximation
         DO i = 1, numSeasonalSegs
              startPoint = vec_seasonal_brkpts(i)
              endPoint = vec_seasonal_brkpts(i+1)
              IF (i== numSeasonalSegs) THEN
                 CALL regress(vec_timestamps_seasonal_pres(startPoint:endPoint),   &
                      &     vec_util(startPoint:endPoint), "harmon",        &
                      & vec_fit_pres(startPoint:endPoint),seasonal_fit_coeffs, &
                      & bfast_numHarmonics, work_arr)
              ELSE
                 CALL regress(vec_timestamps_seasonal_pres(startPoint:endPoint-1),   &
                      &     vec_util(startPoint:endPoint-1),"harmon",   &
                      & vec_fit_pres(startPoint:endPoint-1), seasonal_fit_coeffs, &
                      & bfast_numHarmonics, work_arr)
              ENDIF
         END DO
         
         ! get the Hamming distance between the previous brkpts and current brkpts.
         CALL hammingDistance(vec_trend_brkpts, vec_trend_brkpts_old, hamTrend)
         CALL hammingDistance(vec_seasonal_brkpts, vec_seasonal_brkpts_old, hamSeason)
         
         iter = iter+1
    ENDDO

    ! fill in the missing points. We just want the trend. No interest in season.
    !bfast_summary(:, pixel_x, pixel_y) = 0 !we already initialized. no need here.
    work_arr%tmp_mat(1:num_obs,1) = 1
    work_arr%tmp_mat(1:num_obs,2) = vec_timestamps(1:num_obs)
    DO j = 1, numTrendSegs
          if (j== numTrendSegs) then
             startPoint = vec_timestamps_pres(vec_trend_brkpts(j))
             endPoint = vec_timestamps(num_obs)
             FORALL (i=1:num_obs, (vec_timestamps(i) >= startPoint)) &
              &   bfast_summary(i, pixel_x, pixel_y) = DOT_PRODUCT(work_arr%tmp_mat(i, 1:2), linear_fit_coeffs(j,1:2))
          else
             startPoint = vec_timestamps_pres(vec_trend_brkpts(j))
             endPoint = vec_timestamps(vec_trend_brkpts(j+1))
             FORALL (i=1:num_obs, (vec_timestamps(i) >= startPoint .and. & 
                     &   vec_timestamps(i) < endPoint))  &
                   bfast_summary(i, pixel_x, pixel_y) = DOT_PRODUCT(work_arr%tmp_mat(i, 1:2), &
                                  &  linear_fit_coeffs(j,1:2))
          endif 
    END DO
END DO

DEALLOCATE (work_arr%tmp_mat, work_arr%rtmp_vec1, work_arr%rtmp_vec2,  &
          work_arr%itmp_vec1, work_arr%itmp_vec2)
END SUBROUTINE bfast

END MODULE bfast_fn_par
