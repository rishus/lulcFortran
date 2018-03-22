MODULE ltr_fn_par
CONTAINS

SUBROUTINE ltr(imgIndex)
use globalVars_par
use real_precision
implicit none

! ************* input vars **********************
INTEGER(kind = 8), intent(in) :: imgIndex(:)
! ***********************************************

! ************* utility vars ********************
!general arrays and other variables
INTEGER(KIND=8) :: num_pixels, pixel, pixel_x, pixel_y !just like ewmacd
REAL(KIND=R8),DIMENSION(num_obs) :: vec_timestamps !wasn't needed in ewmacd 
REAL(KIND=R8),DIMENSION(num_obs) :: vec_timestamps_pres, vec_obs_pres !t and D of ewmacd
INTEGER :: i, len_training_vec, num_valid_obs, pctr, mpnp1, num_days_gone
INTEGER :: nSegs, nParams, prev_model_ind
integer(kind=2), dimension(num_obs) :: presInd, vec_timestamps_edited
integer :: nVerts_firstModel, numInitVerts1, numInitVerts2
integer :: nVerts_thisModel, nVerts_prevModel, this_model_ind
integer, dimension(num_obs) :: newInitBrkpts, initBrkpts, updatedBrkpts
LOGICAL :: accepted, accept
type(ltrWA) :: work_arr
type(bestTrace) :: thisModelTrace, thisModelTrace_alternate

!arrays for 'work'
REAL(KIND=R8), DIMENSION(num_obs) :: vec_obs_despiked
REAL(KIND=R8), DIMENSION(ltr_mu+ltr_nu+1) :: all_fstats
TYPE(ltr_model), DIMENSION(ltr_mu + ltr_nu + 1) :: my_models
TYPE(fittingStats_ltr) :: thisModelStats
real :: slope, intercept
integer, dimension(ltr_mu + ltr_nu+1) :: vecTrendBrks, brkptsGI
integer, dimension(ltr_mu + ltr_nu+1, 2) :: brkPtYrDoy

!other variables (parameters)
INTEGER :: bestModelInd, num_models_evaluated, bestModelIndTmp(1:1)
INTEGER :: order, numTotalVertices, numInternalVertices, splines_l, ktimesl
INTEGER :: ctntyCondsPerVertex, sum_ctntyCondsPerInternalVertex, n_dim
INTEGER :: num_models, left, right
REAL(kind = r8) :: x1, y1, x2, y2
! ***********************************************

! ********* PURPOSE **************
!   Main code for algorithm LandTrendR
!   Supporting files are all in utilities_ltr.f90, bsplines.f90, and BRATIO.f90
!   Output is written to ltr_summary
! ***********************************************


!numTrendBrks = ltr_max_numBrks ! this is internal number of breaks.

ALLOCATE (work_arr%tmp_mat(num_obs, 5),  &
          work_arr%rtmp_vec1(num_obs), work_arr%rtmp_vec2(num_obs),  &
          work_arr%itmp_vec1(num_obs), work_arr%itmp_vec2(num_obs))
num_pixels = SIZE(imgIndex, 1)

vec_obs_despiked = mdv
DO pixel = 1, num_pixels

     pixel_x = mod(imgIndex(pixel), NumCols)
     IF (pixel_x == 0) THEN
         pixel_x = NumCols
     ENDIF
     pixel_y = INT(CEILING(REAL(imgIndex(pixel), KIND=R8)/NumCols))  

     ! ********* develop the presInd vector *******************
     pctr = 0
     presInd = 0
     len_training_vec = 0
     num_days_gone = 0
     vec_timestamps_edited(1) = tyeardoy(1, 2)
     DO i = 1, num_obs
        IF (i > 1) THEN    
           IF (tyeardoy(i,1) /= tyeardoy(i-1, 1)) THEN
              IF (mod(INT(tyeardoy(i,1), KIND=4),4) == 0 ) then
                 num_days_gone = num_days_gone + 366
              ELSE
                 num_days_gone = num_days_gone + 365
              ENDIF
              vec_timestamps_edited(i)=int(num_days_gone+INT(tyeardoy(i,2),KIND=4), kind=2)
           ELSE
              vec_timestamps_edited(i)=int(num_days_gone+INT(tyeardoy(i,2),KIND=4), kind=2)
           ENDIF
        ENDIF
        ! the present index part
        IF (input_mat(i, pixel_x, pixel_y) > ewmacd_lowthresh) THEN
           pctr = pctr + 1
           presInd(pctr) = int(i, kind=2)
           IF ((tyeardoy(i,1) >= trainingStart) .and. (tyeardoy(i,1)  &
                   < trainingEnd)) THEN
              len_training_vec = len_training_vec + 1
           ENDIF
        ENDIF
     END DO
     num_valid_obs = pctr

     IF (len_training_vec <= 2.0 * ewmacd_numHarmonics+1 ) THEN
        CYCLE     !note the use of ewmacd harmonics here
     ENDIF
    ! **************** prepare data *************************
     
    vec_timestamps  = vec_timestamps_edited
    vec_timestamps_pres(1:num_valid_obs) = &
            (/ (vec_timestamps(presInd(i)), i=1,num_valid_obs) /)
    vec_obs_pres(1:num_valid_obs) = &
            (/ (input_mat(presInd(i), pixel_x, pixel_y), i = 1,num_valid_obs)/)

    ! ********************** actual processing starts here ******

    mpnp1 = ltr_mu + ltr_nu + 1

    DO i = 1, mpnp1
        ALLOCATE (my_models(i)%yfit(num_valid_obs), &
                  my_models(i)%vertices(mpnp1), &
                  my_models(i)%vertYvals(mpnp1), &
                  my_models(i)%slopes(mpnp1))
    END DO

    !despike
    CALL despike(vec_obs_pres(1:num_valid_obs), num_valid_obs, &
                 vec_obs_despiked(1:num_valid_obs))

    !Find potential breakpoints: ltr_mu, ltr_nu, ltr_distwtfactor are global vars
    !Vertices correspond to indices of the presInd vector
    CALL getInitVerts(vec_timestamps_pres(1:num_valid_obs), &
                      vec_obs_despiked(1:num_valid_obs), num_valid_obs, &
                      mpnp1, initBrkpts, numInitVerts1, &
                      work_arr)

    !Cull by angle change. Output: newInitVerts
    CALL cullByAngleChange(vec_timestamps_pres(1:num_valid_obs), &
                       vec_obs_despiked(1:num_valid_obs), &
                       num_valid_obs, initBrkpts(1:numInitVerts1), numInitVerts1, &
                       newInitBrkpts(1:ltr_mu+2), numInitVerts2)

    !Fit trajectory in this model. That is, find first model. Also get it's goodness.
    ALLOCATE (thisModelTrace%vertices(numInitVerts2), &
              thisModelTrace%vertYvals(numInitVerts2), &
              thisModelTrace%yFitVals(num_valid_obs), &
              thisModelTrace%slopes(numInitVerts2-1))

    ! Output: thisModelTrace
    CALL findBestTrace(vec_timestamps_pres(1:num_valid_obs), &
                       vec_obs_despiked(1:num_valid_obs), num_valid_obs, &
                       newInitBrkpts(1:numInitVerts2), numInitVerts2,  &
                       numInitVerts2-1, thisModelTrace, work_arr)

    !Note: IDL code uses an unnecessarily complicated expression for 
    !nSegs :o :( :/. But we do agree with them on the no. of parameters.
    nSegs = numInitVerts2 - 1  !TODO: Why this?
    nParams = nSegs
    !Output: thisModelStats. 
    !        Attributes -- 17 attributes in all. Except yfit, all others are scalars.
    !                      Further, yfit is just a copy of yFitVals from thisModelTrace.
    CALL calcFittingStats(vec_obs_despiked(1:num_valid_obs), &
                          thisModelTrace%yFitVals(1:num_valid_obs), & 
                          num_valid_obs, nParams, thisModelStats)

    !put together first model
    nVerts_firstModel = numInitVerts2
    ! number of vertices for this model
    my_models(1)%numVertices = nVerts_firstModel
    ! vertices for this first model
    my_models(1)%vertices(1:nVerts_firstModel) = thisModelTrace%vertices(1:nVerts_firstModel)
    ! fit obtained using anchored regression. Variable 'y' in the ps.
    my_models(1)%yfit(1:num_valid_obs) = thisModelTrace%yFitVals(1:num_valid_obs)
    ! yfit values AT the vertices
    my_models(1)%vertYvals(1:nVerts_firstModel) = thisModelTrace%vertYvals(1:nVerts_firstModel)
    ! slopes for this first model
    my_models(1)%slopes(1:nVerts_firstModel-1) = thisModelTrace%slopes(1:nVerts_firstModel - 1)
    ! statistics for this first model
    my_models(1)%fstat = thisModelStats%fstat
    my_models(1)%p_of_f = thisModelStats%p_of_f
    my_models(1)%aicc = thisModelStats%aicc

    ! Simplify (reduce) model, one vertex at a time, to generate \mu different models.
    prev_model_ind = 1
    this_model_ind = 2
    nVerts_prevModel = nVerts_firstModel
    DO WHILE (nSegs > 1)
       nVerts_thisModel = nVerts_prevModel -1
       ! Step 1: update vertices
       ! Output: updatedVertices, stored in updatedBrkpts(1:nVerts_thisModel)
       CALL takeOutWeakest(my_models(prev_model_ind),  &
                           ltr_recoveryThreshold, &
                           vec_timestamps_pres(1:num_valid_obs), &
                           vec_obs_despiked(1:num_valid_obs), &
                           my_models(prev_model_ind)%vertices(1:nVerts_prevModel), &
                           my_models(prev_model_ind)%vertYvals(1:nVerts_prevModel), &
                           num_valid_obs, nVerts_prevModel, nVerts_thisModel,&
                           updatedBrkpts(1:nVerts_thisModel))
      ! Step 2: update trace (i.e., vertYals, yFit, slopes, and a copy of vertices)
      ! Output: thisModelTrace. 
      !         It has 4 attributes-- 
      !         1. vertices:  no. of vertices keeps reducing
      !         2. vertYvals: no. of vertYvals also keeps reducing
      !         3. yFitVals:  no. of yFitVals stays same
      !         4. slopes:    no. of slopes also keeps reducing
      ! Be careful about the indices to be used for each attributes.
      CALL findBestTrace(vec_timestamps_pres(1:num_valid_obs), &
                         vec_obs_despiked(1:num_valid_obs), num_valid_obs, &
                         updatedBrkpts(1:nVerts_thisModel),  &
                         nVerts_thisModel, nVerts_thisModel-1, thisModelTrace, work_arr)

      ! Step 3: update statistics (everything except yfit is scalar. yfit is only a copy
      !         of the yfit from trace)
      nSegs = nVerts_thisModel - 1
      nParams = nSegs
      CALL calcFittingStats(vec_obs_despiked(1:num_valid_obs), &
                           thisModelTrace%yFitVals, num_valid_obs, nParams, &
                           thisModelStats)
      ! Step 4: Update my_models, prev_model_ind
      my_models(this_model_ind)%numVertices = nVerts_thisModel
      my_models(this_model_ind)%vertices(1:nVerts_thisModel) = &
          thisModelTrace%vertices(1:nVerts_thisModel)  ! same as updatedBrkpts(1:nVerts_thisModel)
      my_models(this_model_ind)%vertYvals(1:nVerts_thisModel) = &
              thisModelTrace%vertYvals(1:nVerts_thisModel)
      my_models(this_model_ind)%slopes(1:nVerts_thisModel-1) = &
              thisModelTrace%slopes(1:nVerts_thisModel-1)
      my_models(this_model_ind)%yfit(1:num_valid_obs) = thisModelTrace%yFitVals
      my_models(this_model_ind)%fstat = thisModelStats%fstat
      my_models(this_model_ind)%p_of_f = thisModelStats%p_of_f
      my_models(this_model_ind)%aicc = thisModelStats%aicc

      prev_model_ind = prev_model_ind + 1
      this_model_ind = this_model_ind + 1
      nVerts_prevModel = nVerts_thisModel

     END DO

     num_models = this_model_ind - 1
     
     ! Pick best model
     accepted = .FALSE.
     all_fstats(1:num_models) = (/ (my_models(i)%fstat, i=1, num_models) /)
     num_models_evaluated = 0
     DO WHILE ((accepted .EQV. .FALSE.) .and. (num_models_evaluated < num_models))
        num_models_evaluated = num_models_evaluated + 1
        CALL pickBestModel(my_models, num_models, bestModelInd)
        IF (bestModelInd /= -1) THEN
           CALL checkSlopes(my_models(bestModelInd),   &
                   my_models(bestModelInd)%numVertices, accept)
           IF (accept .EQV. .FALSE.) THEN
              my_models(bestModelInd)%p_of_f = 10000
              my_models(bestModelInd)%fstat  = 0
              my_models(bestModelInd)%aicc    = 0
              accepted = .FALSE.
           ELSE
              accepted = .TRUE.
           ENDIF
        ELSE
           bestModelIndTmp = MINLOC(all_fstats(1:num_models))
           bestModelInd = bestModelIndTmp(1)
           accepted = .TRUE.
        ENDIF
     ENDDO
     !**********************************************************
     !**********************************************************
     !If no good fit found, try the MPFITFN approach
     IF (my_models(bestModelInd)%p_of_f > ltr_pval) THEN
        !redo the whole model generation part, this time with 
        !Levenberg-Marquardt algorithm based fitting.
        !restart from vertices determined by vetVerts3
        !i.e., the vertices POST- cullByAngleChange.
        !my_models is not useless. Use it to store the alternate
        !models. Dimension of the array attributes are still the same.

        !Step 4: Fit trajectory using the marquardt approach. That's a
        !global method in contrast to the above local method.
        ALLOCATE (thisModelTrace_alternate%vertices(numInitVerts2), &
                  thisModelTrace_alternate%vertYvals(numInitVerts2), &
                  thisModelTrace_alternate%yFitVals(num_valid_obs), &
                  thisModelTrace_alternate%slopes(numInitVerts2-1))
        order = 2             ! order = degree + 1. k stands for order
        ! Assuming that the vertices are labelled as 
        ! \xi_1, \xi_2, ...., \xi_{l+1}. Then,
        numTotalVertices = numInitVerts2
        numInternalVertices = numTotalVertices - 2
        splines_l = numTotalVertices - 1 ! or, numInternalVertices + 1
        ktimesl = order * splines_l
        ctntyCondsPerVertex = 1  !becuz landtrend only wants continuous.
                             !no imposition on derivatives.
        sum_ctntyCondsPerInternalVertex = numInternalVertices
        n_dim = ktimesl - sum_ctntyCondsPerInternalVertex  !dim of spline space

        CALL findBestTrace_alternate(vec_timestamps_pres(1:num_valid_obs), &
                       vec_obs_despiked(1:num_valid_obs), num_valid_obs, & 
                       newInitBrkpts(1:ltr_mu+2), numInitVerts2, &
                       order, n_dim, thisModelTrace_alternate)
        nSegs = numInitVerts2 - 1
        nParams = nSegs 
        CALL calcFittingStats(vec_obs_despiked(1:num_valid_obs), &
                              thisModelTrace_alternate%yFitVals(1:num_valid_obs), &
                              num_valid_obs, nParams, thisModelStats)
        !put together first alternate model
        nVerts_firstModel = numInitVerts2
        ! number of vertices for this model
        my_models(1)%numVertices = nVerts_firstModel
        ! vertices for this first model
        my_models(1)%vertices(1:nVerts_firstModel) =  &
                        thisModelTrace_alternate%vertices(1:nVerts_firstModel)
        ! fit obtained using anchored regression. Variable 'y' in the ps.
        my_models(1)%yfit(1:num_valid_obs) = &
                            thisModelTrace_alternate%yFitVals(1:num_valid_obs)
        ! yfit values AT the vertices
        my_models(1)%vertYvals(1:nVerts_firstModel) = &
                       thisModelTrace_alternate%vertYvals(1:nVerts_firstModel)
        ! slopes for this first model
        my_models(1)%slopes(1:nVerts_firstModel-1) = &
                      thisModelTrace_alternate%slopes(1:nVerts_firstModel - 1)
        ! statistics for this first model
        my_models(1)%fstat = thisModelStats%fstat
        my_models(1)%p_of_f = thisModelStats%p_of_f
        my_models(1)%aicc = thisModelStats%aicc

        prev_model_ind = 1
        this_model_ind = 2
        nVerts_prevModel = nVerts_firstModel
        DO WHILE (nSegs > 1)
            nVerts_thisModel = nVerts_prevModel -1
            ! Step 1: update vertices
            ! Output: updatedVertices, stored in updatedBrkpts(1:nVerts_thisModel)
            CALL takeOutWeakest_alternate(vec_timestamps_pres(1:num_valid_obs), &
                       vec_obs_despiked(1:num_valid_obs), &
                       num_valid_obs,  &
                       my_models(prev_model_ind)%vertices(1:nVerts_prevModel), &
                       my_models(prev_model_ind)%vertYvals(1:nVerts_prevModel), &
                       nVerts_prevModel, updatedBrkpts(1:nVerts_thisModel))

            numTotalVertices = nVerts_thisModel
            numInternalVertices = numTotalVertices - 2
            splines_l = numTotalVertices - 1 ! or, numInternalVertices + 1
            ktimesl = order * splines_l
            ctntyCondsPerVertex = 1  !becuz landtrend only wants continuous.
            !no imposition on derivatives.
            sum_ctntyCondsPerInternalVertex = numInternalVertices
            n_dim = ktimesl - sum_ctntyCondsPerInternalVertex  !dim of spline space
           
             CALL findBestTrace_alternate(vec_timestamps_pres(1:num_valid_obs), &
                       vec_obs_despiked(1:num_valid_obs), num_valid_obs, & 
                       updatedBrkpts(1:nVerts_thisModel), nVerts_thisModel, &
                       order, n_dim, thisModelTrace_alternate)
             nSegs = nVerts_thisModel - 1
             nParams = nSegs 
             CALL calcFittingStats(vec_obs_despiked(1:num_valid_obs), &
                          thisModelTrace_alternate%yFitVals(1:num_valid_obs), &
                              num_valid_obs, nParams, thisModelStats)

             ! number of vertices for this model
             my_models(this_model_ind)%numVertices = nVerts_thisModel
             ! vertices for this first model
             my_models(this_model_ind)%vertices(1:nVerts_thisModel) =  &
                     thisModelTrace_alternate%vertices(1:nVerts_thisModel)
             ! fit obtained using anchored regression. Variable 'y' in the ps.
             my_models(this_model_ind)%yfit(1:num_valid_obs) = &
                     thisModelTrace_alternate%yFitVals(1:num_valid_obs)
             ! yfit values AT the vertices
             my_models(this_model_ind)%vertYvals(1:nVerts_thisModel) = &
                     thisModelTrace_alternate%vertYvals(1:nVerts_thisModel)
             ! slopes for this first model
             my_models(this_model_ind)%slopes(1:nVerts_thisModel-1) = &
                     thisModelTrace_alternate%slopes(1:nVerts_thisModel - 1)
             ! statistics for this first model
             my_models(this_model_ind)%fstat = thisModelStats%fstat
             my_models(this_model_ind)%p_of_f = thisModelStats%p_of_f
             my_models(this_model_ind)%aicc = thisModelStats%aicc

             prev_model_ind = prev_model_ind + 1
             this_model_ind = this_model_ind + 1
             nVerts_prevModel = nVerts_thisModel
        END DO

        num_models = this_model_ind - 1

        ! Pick best model
        accepted = .FALSE.
        all_fstats(1:num_models) = (/ (my_models(i)%fstat, i=1, num_models) /)
        num_models_evaluated = 0
        DO WHILE ((accepted .EQV. .FALSE.) .and. (num_models_evaluated < num_models))
           num_models_evaluated = num_models_evaluated + 1
           CALL pickBestModel(my_models, num_models, bestModelInd)
           IF (bestModelInd /= -1) THEN
              CALL checkSlopes(my_models(bestModelInd),   &
                      my_models(bestModelInd)%numVertices, accept)
              IF (accept .EQV. .FALSE.) THEN
                 my_models(bestModelInd)%p_of_f = 10000
                 my_models(bestModelInd)%fstat  = 0
                 my_models(bestModelInd)%aicc    = 0
                 accepted = .FALSE.
              ELSE
                 accepted = .TRUE.
              ENDIF
           ELSE
              bestModelIndTmp = MINLOC(all_fstats(1:num_models))
              bestModelInd = bestModelIndTmp(1)
              accepted = .TRUE.
           ENDIF
        ENDDO
        DEALLOCATE (thisModelTrace_alternate%vertices, &
        thisModelTrace_alternate%vertYvals, thisModelTrace_alternate%yFitVals,&
        thisModelTrace_alternate%slopes)
        
     ENDIF

     ! Summarize
     ! (1) Trend summary:
     ltr_summary(presInd(1:num_valid_obs), pixel_x, pixel_y) = real(my_models(bestModelInd)% &
                                        yfit(1:num_valid_obs), kind=4)
     ! using the array initVerts itself here for storing the final vertices.
     ! these are wrt the pres indices, though
     numTotalVertices = my_models(bestModelInd)%numVertices
     vecTrendBrks(1:numTotalVertices) =  my_models(bestModelInd)%vertices(1:numTotalVertices)
     brkptsGI(1:numTotalVertices) = presInd(vecTrendBrks(1:numTotalVertices))
     brkptsGI(1) = 1
     brkptsGI(numTotalVertices) = num_obs
     DO i = 1, numTotalVertices 
        brkPtYrDoy(i, 1:2) = tyeardoy(brkptsGI(i), 1:2)
     end do
     left = 1
     right = 2
     DO i = presInd(1), presInd(num_valid_obs)
       ! Check if this location was present. If it was, we already have the value; no need
       ! for further calculation.
       IF ( ANY(presInd(1:num_valid_obs) == i ) .EQV. .TRUE.) THEN
          CYCLE
       ELSE
       ! Locate left, right st. vts(i) \in [vts_pres(left), vts_pres(right))
       ! (locate the interval in which this x lies.)
       DO WHILE (vec_timestamps(i) >= vec_timestamps_pres(vecTrendBrks(right)))
              left = left + 1
              right = right + 1
              if (right >= numTotalVertices) then
                 EXIT
              endif
          END DO
         ! Fetch the (x,y) coordinates of the segment representing this interval.
          x1 = vec_timestamps_pres(vecTrendBrks(left))
          x2 = vec_timestamps_pres(vecTrendBrks(right))
          y1 = my_models(bestModelInd)%vertYvals(left)
          y2 = my_models(bestModelInd)%vertYvals(right)
          slope = real(y2-y1)/real(x2-x1)
          intercept = real(y1) - slope * real(x1)
          ltr_summary(i, pixel_x, pixel_y) = slope * real(vec_timestamps(i)) + intercept
       ENDIF
     END DO

     ! the right end:
     x1 = vec_timestamps_pres(vecTrendBrks(left))
     x2 = vec_timestamps_pres(vecTrendBrks(right))
     y1 = my_models(bestModelInd)%vertYvals(left)
     y2 = my_models(bestModelInd)%vertYvals(right)
     slope = real(y2-y1)/real(x2-x1)
     intercept = real(y1) - slope * real(x1)
     DO i = presInd(num_valid_obs)+1 ,num_obs
          ltr_summary(i, pixel_x, pixel_y) = slope * real(vec_timestamps(i)) + intercept
     END DO

     ! the left end:
     x1 = vec_timestamps_pres(vecTrendBrks(1))
     x2 = vec_timestamps_pres(vecTrendBrks(2))
     y1 = my_models(bestModelInd)%vertYvals(1)
     y2 = my_models(bestModelInd)%vertYvals(2)
     slope = real(y2-y1)/real(x2-x1)
     intercept = real(y1) - slope * real(x1)
     DO i = 1, presInd(1)-1
          ltr_summary(i, pixel_x, pixel_y) = slope * real(vec_timestamps(i)) + intercept
     END DO

     DO i = 1, mpnp1
         DEALLOCATE (my_models(i)%yfit, my_models(i)%vertices, &
                   my_models(i)%vertYvals,  my_models(i)%slopes)
     END DO
     DEALLOCATE (thisModelTrace%vertices, thisModelTrace%vertYvals, &
               thisModelTrace%yFitVals, thisModelTrace%slopes)

END DO

DEALLOCATE (work_arr%tmp_mat, work_arr%rtmp_vec1, work_arr%rtmp_vec2,  &
          work_arr%itmp_vec1, work_arr%itmp_vec2)

RETURN
END SUBROUTINE ltr

END MODULE ltr_fn_par

