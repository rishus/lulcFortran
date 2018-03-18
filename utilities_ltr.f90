MODULE utilities_ltr
USE REAL_PRECISION

CONTAINS

SUBROUTINE regress(x, y, work_arr, fit, slope, const)
USE globalVars_par
IMPLICIT NONE
   REAL(KIND=R8), INTENT(IN) :: x(:), y(:)
   INTEGER :: ncols, m, i, lwork, info, DeAllocateStatus
   REAL(KIND=R8), DIMENSION(:,:),  ALLOCATABLE :: A(:,:)
   REAL(KIND=R8), DIMENSION(:),  ALLOCATABLE :: work(:)
   REAL(KIND=R8), INTENT(OUT) :: const
   REAL(KIND=R8), INTENT(INOUT) :: fit(:), slope
   type(ltrWA), intent(inout) :: work_arr
   INTERFACE
      SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
      USE REAL_PRECISION
      CHARACTER:: TRANS
      INTEGER::   INFO, LDA, LDB, LWORK, M, N, NRHS
      REAL(KIND=R8):: A( LDA, * ), B( LDB, * ), WORK( * )
      END SUBROUTINE DGELS
   END INTERFACE

   ncols = 2
   m = SIZE(x)
   ALLOCATE (A(m, ncols))
   A(:, 1) = (/ (1, i = 1, m)  /)
   A(:, 2) = (/ (x(i), i = 1, m) /) 
  
   lwork = min(m, ncols) +  max(max(m, ncols), 1) * nb
   ALLOCATE (work(lwork))

   work_arr%tmp_mat(1:m, 1:ncols) = A
   work_arr%rtmp_vec1(1:m) = y
   CALL DGELS('N', m, ncols, 1, work_arr%tmp_mat, size(work_arr%tmp_mat, 1), &
           work_arr%rtmp_vec1, size(work_arr%rtmp_vec1, 1), work, lwork, info)
   const = work_arr%rtmp_vec1(1)
   slope = work_arr%rtmp_vec1(2)
   
   DO i = 1, m
      fit(i) = const + slope * x(i)
   END DO   
   DEALLOCATE (A, work, stat = DeAllocateStatus)

END SUBROUTINE regress

SUBROUTINE despike(vals, Sfinal, vec_obs_updated)
USE globalVars_par
IMPLICIT NONE

   ! parameters needed in this routine
   ! ltr_despikeTol is a global variable.
   REAL(KIND=R8), INTENT(IN) :: vals(:)
   INTEGER, INTENT(IN) :: Sfinal
   INTEGER :: ctr, i, prop_ind_tmp(1:1), ctr2, prop_ind(Sfinal)
   REAL(KIND=R8) :: prop, md
   REAL(KIND=R8), DIMENSION(Sfinal) :: fwd_diffs, bkwd_diffs, central_diffs
   REAL(KIND=R8), DIMENSION(Sfinal) :: correction, prop_correction
   REAL(KIND=R8), DIMENSION(:), INTENT(INOUT) :: vec_obs_updated

   vec_obs_updated = vals(1:Sfinal)   !(/ (vals(ctr), ctr=1,n) /)
   IF (ltr_despikeTol > 1) THEN
      ltr_despikeTol = 0.9
   ENDIF
   prop = 1.0  ! ltr_despikeTol + 1.0 : TODO check this
   ctr = 1
   !fwd_diffs, bkwd_diffs, central_diffs (1) and (Sfinal) never really get used.
   DO WHILE ((prop.GT. ltr_despikeTol) .and. (ctr <= Sfinal))
      fwd_diffs(1) = 0
      fwd_diffs(2:Sfinal) = (/(vec_obs_updated(i+1)-vec_obs_updated(i), &
                                                         i=1, Sfinal-1) /)
      bkwd_diffs(Sfinal) = 0
      bkwd_diffs(1:Sfinal-1) = (/ (vec_obs_updated(i) - vec_obs_updated(i-1), &
                                                         i=2,Sfinal) /)
      central_diffs(1) = 0
      central_diffs(2:Sfinal-1) = (/(vec_obs_updated(i+1)- vec_obs_updated(i-1),&
                                                         i = 2,Sfinal-1)/)
      central_diffs(Sfinal) = 0

      prop_correction = 0
      correction = 0
      DO i = 2, Sfinal-1  !no correction at the ends
         md = MAX(ABS(fwd_diffs(i)), ABS(bkwd_diffs(i)))
         ! what is both are zero?
         IF (md <= EPSILON(FLOAT(1))) THEN
            md = central_diffs(i)  !md = *** ? or correction(i) = ***?
         ENDIF
         !what if central diffs is also 0? Then we have a 0/0 form.
         !then correction = 0 becuz rhs of correction formula is 0.
         IF (md == 0) THEN
            correction(i) = 0
         ELSE
            prop_correction(i) = 1.0 - (ABS(central_diffs(i))/md)
            correction(i) = prop_correction(i) * 0.5 * (vec_obs_updated(i+1) - &
                    &            2.0*vec_obs_updated(i) + vec_obs_updated(i-1))
         ENDIF
      END DO

      prop =  MAXVAL(prop_correction(2:Sfinal-1))
      prop_ind_tmp = MAXLOC(prop_correction(2:Sfinal-1))
      ctr2 = 1
      prop_ind(1) = prop_ind_tmp(1) + 1
      DO i = prop_ind(1)+1, Sfinal-1
         IF (prop_correction(i) == prop) THEN
             ctr2 = ctr2 + 1
             prop_ind(ctr2) = i
         ENDIF
      END DO
      DO i = 1, ctr2
         vec_obs_updated(prop_ind(i)) = vec_obs_updated(prop_ind(i)) + &
                                       correction(prop_ind(i))
      END DO
      ctr = ctr + 1
  END DO

END SUBROUTINE despike

SUBROUTINE scoreSegments(vec_timestamps, vec_obs, vertices, &
                                   vertexCount, leny, segmentScores, work_arr)
USE globalVars_par
IMPLICIT NONE

  ! Global variables needed here:
  ! ltrWA deifinition
  ! the score of a segment is basically the error in that fit.
  INTEGER, INTENT(IN) :: vertexCount, vertices(:), leny
  REAL(KIND=R8), INTENT(IN) :: vec_timestamps(:), vec_obs(:)
  TYPE(ltrWA) :: work_arr
  INTEGER :: i, span, startInd, endInd, j
  REAL(KIND=R8), DIMENSION(leny) :: yfit
  REAL(KIND = R8) :: dummy1, dummy2
  REAL(KIND=R8), INTENT(INOUT) :: segmentScores(:)
 
  DO i=1,vertexCount-1
     startInd = vertices(i)
     endInd = vertices(i+1)
     span = endInd -  startInd +1  !total number of pts in [t_i, t_{i+1}]
     IF (span > 2) THEN 
        !if we've done desawtooth, it's possible that all of the values in a
        !segment have the same value, in which case regress would choke. 
        !So deal with that.
        IF (MAXVAL(vec_obs(startInd:endInd)) - &
                MINVAL(vec_obs(startInd:endInd)) > 0) THEN
           CALL regress(vec_timestamps(startInd:endInd),   &
                   vec_obs(startInd:endInd), work_arr, yfit(1:span), dummy1, &
                   dummy2)
        ELSE
           yfit(1:span) = vec_obs(startInd : endInd)
        ENDIF
        work_arr%rtmp_vec2(1:span) = vec_obs(startInd:endInd) - yfit(1:span)
        work_arr%rtmp_vec2(1:span) = &
                (/ (work_arr%rtmp_vec2(j)* work_arr%rtmp_vec2(j), j=1,span)/)
        segmentScores(i) = SUM(work_arr%rtmp_vec2(1:span)) / span
     ELSE
        ! if vertex_i and vertex_{i+1} are really really right next to each other
        ! then, we can't split that interval further anyways. So set the MSE for 
        ! such an interval equal to 0, then the other interval will get chosen
        ! for the next steps (of splitting etc).
        ! This situation will arise, for example, when there is a 'sharp' drop
        ! in values. Check out the figures in CG paper.
        segmentScores(i) = 0
     ENDIF
  END DO

END SUBROUTINE scoreSegments

SUBROUTINE splitSeries(vec_timestamps, vec_obs, numObs, endSegment, &
                           distTest, ok, maxdiffInd, work_arr)
USE globalVars_par
IMPLICIT NONE
  ! Global variables needed:
  ! ltrWA definition 
! given an x and y split the series into two smaller segments. 
! However, invoke a rule where there can be no series where the value decreases
! (implication of recovery) for only 1 or 2 years --
! This will help with the overfitting, since this is not really a prominent type
! of phenomenon, and if it's minor anyway a coarser fit would make more sense. 
! Endsegment is a flag set to .true. if this segment is at the end of the time period,
! when we don't invoke the recovery restriction rule -- it will get thrown out later
! if it's really extreme.

  REAL(KIND=R8), INTENT(IN) :: vec_timestamps(:), vec_obs(:)
  INTEGER, INTENT(IN) :: numObs
  LOGICAL, INTENT(IN) :: endSegment, distTest
  TYPE(ltrWA) :: work_arr
  REAL(KIND=R8), DIMENSION(numObs) :: yfit, diff
  INTEGER :: maxTmpInd(1:1)
  REAL(KIND=R8) :: dummy1, dummy2
  INTEGER, INTENT(OUT) :: maxdiffInd   ! maxdiff will be a vertex (the 
                                       ! vertex this subroutine picks up)
  LOGICAL, intent(OUT) :: ok

  ok = .FALSE.
  maxDiffInd = 0
  CALL regress(vec_timestamps, vec_obs, work_arr, yfit, dummy1, dummy2)
  diff = ABS(yfit - vec_obs)
  diff(1) = 0
  diff(numObs) = 0  ! end points are already vertices. So take them out of consideration.

  IF ((distTest ).AND.(endSegment)) THEN
      IF (vec_obs(numObs).LE. vec_obs(numObs-1)) THEN
         diff(numObs-1) = 0    ! We know there are 3 segments at least, 
                               ! because assured in the calling program.
      ENDIF
  ENDIF
  
  maxTmpInd = MAXLOC(diff)
  maxdiffInd = maxTmpInd(1)

  ! if maxdiff is 1, then we know there were no segments that met the rule of no 1-yr recovery
  IF (maxdiffInd > 1) THEN
     ! if maxDiff occurs at an interior index, then we will accept it.
     ! Otherwise we'll look for sth else.
     ! we will use 0 for 'false' and 1 for 'true'
     ok = .TRUE.  ! TODO: confirm whether it should be 0 or 1
  ENDIF

END SUBROUTINE splitSeries

SUBROUTINE getInitVerts(vec_timestamps, vec_obs, numObs, mpnp1, &
                    initVertices, initVertCount, work_arr)
USE globalVars_par
IMPLICIT NONE
   ! Global variables needed:
   ! ltrWA definition
   ! ltr_distwtfactor
   !TODO: 2. I THINK the pseudocode wants us to calculate MSE only in the newest
   !         interval found. But this subroutine finds MSE of all existing
   !         intervals. The latter sounds like a much better idea. 
   !vec_timestamps, vec_obs: arrays of length numObs, allocated smwhere above in 
   !                         the calling program
   !initVertices:            array of length maxCount, allocated in the calling program.
   !mpnp1:                   maximum number of (initial) vertices allowed
   !initVertCount:           actual number of (initial) vertices found. 
   !                         initVertCount <= mpnp1

   REAL(KIND=R8), INTENT(IN) :: vec_timestamps(:), vec_obs(:)
   INTEGER, INTENT(IN) :: mpnp1, numObs
   REAL(KIND=R8), DIMENSION(mpnp1-1) :: mses !num segs = num verts - 1
   INTEGER :: currNumVerts, i, ctr, newVertInd, max_mse_ind(1:1), mmi
   INTEGER :: numMaxVerts, vertexFlags(numObs)
   REAL(KIND=R8) :: maxMSE
   LOGICAL :: isFirstSeg, isLastSeg, ok
   INTEGER, INTENT(INOUT) :: initVertices(:)
   INTEGER, INTENT(OUT) :: initVertCount
   TYPE(ltrWA) :: work_arr

   numMaxVerts = MIN(mpnp1, numObs-2)  !internal
   ! set initVertices to the end points
   vertexFlags = 0
   vertexFlags(1) = 1
   vertexFlags(numObs) = 1
   initVertices(1) = 1
   initVertices(2) = numObs
   currNumVerts = 2
   mses = 0
 
   !there is a step to 'grab a few with big changes in convolve value'
   !but it looks redundant to me and am, therefore, ignoring it.
   DO WHILE (currNumVerts .LT. numMaxVerts)
      !The score of a segment is basically the mse in that fit.
      !Send in the entire domain and score EACH segment.
      CALL scoreSegments(vec_timestamps, vec_obs, &
                         initVertices(1:currNumVerts), &
                         currNumVerts, numObs, mses(1:currNumVerts-1), &
                         work_arr)
      ok = .FALSE.
      DO WHILE (ok .EQV. .FALSE.)
          ! Amongst ALL segments, find the segment that has highest MSE
          maxMSE = MAXVAL(mses(1:currNumVerts-1))
          IF (maxMSE .EQ. 0) THEN
             !bail if we can't get high enough without breaking recovery rule
             !currNumVerts and initVertices is already set to whatever in the 
             !previous iteration
             initVertCount = currNumVerts
             RETURN
          ENDIF
          max_mse_ind = MAXLOC(mses(1:currNumVerts-1))
          mmi = max_mse_ind(1)
          isFirstSeg = (mmi .EQ. 1)
          isLastSeg = (mmi .EQ.  currNumVerts - 1)
          !just use distweightfactor to determine if disturbance should be
          !considered in initial segments.
          !Take the interval with HIGHEST MSE and split it into two at the
          !point of highest deviation. So this 'highest mse' interval is
          !sent into splitSeries. Both endpoints are included.
          CALL splitSeries(vec_timestamps(initVertices(mmi): initVertices(mmi+1)), &
                           vec_obs(initVertices(mmi): initVertices(mmi+1)), &
                           SIZE(vec_obs(initVertices(mmi): initVertices(mmi+1))), &
                           isLastSeg, (ltr_distwtfactor /= 0), ok, newVertInd, &
                           work_arr)
          IF (ok .eqv. .FALSE.) THEN   
             !means, this interval yielded a preexisting vertex. So this segment is useless. 
             mses(mmi) = 0    !look at the next best option
          ENDIF
      END DO

      !the vertex picked by split series, but adjusted so that it corresponds
      !to the global indices. Fortran does not have SORT. So, another strategy
      !has to be used to maintain a sorted list of vertices.
      vertexFlags(newVertInd + initVertices(mmi)-1) = 1
      initVertices(1) = 1
      ctr = 1
      DO i = 2, numObs
         IF (vertexFlags(i) == 1) THEN
            ctr = ctr + 1
            initVertices(ctr) = i
         ENDIF
      END DO
      currNumVerts = ctr  !which shd be same as currNumVerts = currNumVerts+1 

      IF (currNumVerts .GT. 20) THEN             !TODO: Why 20?
         initVertCount = currNumVerts
         RETURN
      ENDIF
   END DO

   initVertCount = currNumVerts

END SUBROUTINE getInitVerts

SUBROUTINE angleDiff(xcoords, ycoords, yrange, distweightfactor,  angDiff)
USE globalVars_par
IMPLICIT NONE
! Global variables needed here:
! Pi ... just one parameter
!
! Note that this is NOT the same as compareAngles. compareAngles gets used 
! later, separately, while collapsing segments in landtrendr_point_run.
!
!   1. Need three points -- middle point is the one that gets the score
!       the others are the ones preceding and following
!   2. Note that ycoords needs to be scaled to the range of the whole
!       trajectory for this to really be meaningful
!   3. Distweightfactor helps determine how much weight is given to angles
!       that precede a disturbance.  If set to 0, the angle difference is
!       passed straighton.
!   4. If disturbance (positive), give more weight

  REAL(KIND=R8), INTENT(IN) :: xcoords(:), ycoords(:), yrange 
  REAL(KIND=R8) :: angle1, angle2, mx, mn, scaler, distweightfactor
  REAL(KIND = R8), INTENT(OUT) :: angDiff

  IF ((SIZE(xcoords) .LT. 3) .OR. (SIZE(ycoords) .LT. 3)) THEN
     angDiff = -1
  ENDIF
  angDiff = 0.0_r8  
  IF (distweightfactor .EQ. 0 ) THEN
  ENDIF
  angle1 = ATAN2((ycoords(2) - ycoords(1)),(xcoords(2) - xcoords(1))) * 180.0 / Pi
  angle2 = ATAN2((ycoords(3) - ycoords(2)),(xcoords(3) - xcoords(2))) * 180.0 / Pi
  !TODO: Check correctness of the lines from this point on
  mx = MAXVAL((/angle1, angle2/))
  mn = MINVAL((/angle1, angle2/))
  scaler = MAXVAL((/REAL(0, KIND=R8), &
          ((ycoords(3) - ycoords(2))*distweightfactor)/yrange /)) + 1.0
  angDiff = (mx-mn) * scaler
!  angDiff = MAXVAL((/ ABS(angle1), ABS(angle2) /) )*scaler
  !have no idea why this is there in the IDL code.
  if (abs(angDiff) <= 0.00001 ) then
     angDiff = 0.0
  endif

END SUBROUTINE angleDiff

SUBROUTINE get_next(CA_currVertInd, nextArray, nextVert)
IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: CA_currVertInd, nextArray(:)
    INTEGER, INTENT(OUT) :: nextVert
    
    nextVert = nextArray(CA_currVertInd)

END SUBROUTINE get_next

SUBROUTINE get_prev(CA_currVertInd, nextArray, CA_prevVertInd)
IMPLICIT NONE
    
    ! size of nextArray is one less than the size of currVerts array.
    ! CA_currVertInd minimum value can be 2, max value can be nVerts.
    ! nextArray(1) will never be -1.
    INTEGER, INTENT(IN) :: CA_currVertInd, nextArray(:)
    INTEGER, intent(out) :: CA_prevVertInd
    
    if (CA_currVertInd == 2) then
       CA_prevVertInd = 1
       return
    endif

    CA_prevVertInd = CA_currVertInd - 1
    do while ( nextArray(CA_prevVertInd) /= -1 )
       CA_prevVertInd = CA_prevVertInd - 1
    end do

    RETURN
END SUBROUTINE get_prev

SUBROUTINE delPos(currVertInd, nextArray)
IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: nextArray(:)
    INTEGER, INTENT(IN) :: currVertInd
    integer :: prev

    prev = currVertInd - 1
    do while (nextArray(prev) /= -1)   ! eventually, prev will become =1. then 
                                       ! nextArray(1) /= -1 will not be satisfied.
       prev = prev - 1
    end do
    nextArray(prev) = nextArray(currVertInd)
    nextArray(currVertInd) = -1

    RETURN
END SUBROUTINE delPos

SUBROUTINE cullByAngleChange(vec_timestamps, vec_obs, numObs, origVerts, nOrigVerts, &
                &  finalVerts, nFinalVerts)
USE globalVars_par
IMPLICIT NONE
   ! Global parameters being used:
   ! ltr_mu  
   ! ltr_distwtfactor 
   !
   ! Culls by angle change
   !orgiVerts:  origVerts(1:nOrigVerts) is the set of starting vertices,
   !            including 1st and last vertices.
   !mu:         maximum number of segments allowed   
   !nu:         number of extra vertices allowed on top of mu+1 while 
   !            forming the initial stencil
   !finalVerts: finalVerts(1:nFinalVerts) is the final set of vertices we will use

   REAL(KIND=R8), INTENT(IN) :: vec_timestamps(:), vec_obs(:)
   INTEGER, INTENT(IN) :: origVerts(:), nOrigVerts, numObs
   INTEGER, DIMENSION(nOrigVerts) :: currVerts 
   INTEGER :: nAngles, nToRemove, endVertMarker, endAngMarker
   INTEGER :: minVertInd, minAngInd, i, itmp(1:1), nVerts
   REAL(KIND=R8) :: yr, yrScaled
   REAL(KIND=R8), DIMENSION(numObs) :: yScaled
   REAL(KIND=R8), DIMENSION(nOrigVerts-2) :: angles !1st & last vertices not included
   INTEGER, INTENT(INOUT) :: finalVerts(:)  !dimension shd be mu+1
   INTEGER, INTENT(OUT) :: nFinalVerts      !should be mu+1, mu=num segs
   
   nVerts = nOrigVerts
   nToRemove = nVerts - (ltr_mu+1)  ! shd turn out to be (mu+nu+1)-(mu+1) = nu
                                
   !if we don't have enough, just return what was given get initial slope ratios
   !across all points
   IF ((nToRemove .LT. 0) .OR. (nVerts .LE. 3)) THEN
        nFinalVerts = nVerts
        finalVerts(1:MINVAL((/nVerts, ltr_mu+1/))) = origVerts(1:MINVAL((/nVerts, ltr_mu+1/)))
        RETURN
   ENDIF

   !for angles, need to make it scaled in y like we see in display
   yr = MAXVAL(vec_obs) - MINVAL(vec_obs)
   yScaled = ((vec_obs - MINVAL(vec_obs))/yr) * &
                   (MAXVAL(vec_timestamps) - MINVAL(vec_timestamps))
   yrScaled = MAXVAL(yscaled) - MINVAL(yscaled)

   ! calculate all angles (interior vertices only)
   currVerts = origVerts
   angles = 0
!   print *, currVerts(1:nVerts)
   DO i = 2, nVerts - 1
      ! the first element in angles cors to the angle
      ! between the first segment and the second segment
      CALL angleDiff(vec_timestamps(currVerts(i-1:i+1)),  &
                     yscaled(currVerts(i-1:i+1)), &
                     yrScaled, ltr_distwtfactor, angles(i-1))
                                                                 
   END DO
!   print *, angles(1:nVerts-2)

   !Now go through and iteratively remove them
   endVertMarker = nVerts              !keeps track of how many good ones are left 
   nAngles = nVerts - 2  !1st and last vertex won't have any cors angleDiffs
   endAngMarker = nAngles  !currently excludes the 1st and last vertex

   DO i=1, nToRemove                 ! nToRemove = nu
      !pick from the slope diffs in play (not the ones at the end, 
      !which are shifted there from prior iterations)
      itmp = MINLOC(angles(1:endAngMarker))
      minAngInd = itmp(1)
      minVertInd = minAngInd + 1          
!      print *, '---- i = ', i, '--------'
!      print *, 'currVerts', currVerts(1:endVertMarker)
!      print *, 'angles', angles(1:endAngMarker)
!      print *, '                     '

      IF (minVertInd .EQ. (endVertMarker-1)) THEN
         ! drop this vertex. update the endVertMarker
         currVerts(endVertMarker-1) = currVerts(endVertMarker) 
         endVertMarker = endVertMarker - 1

         ! drop this angle. update the endAngMarker
         endAngMarker = endAngMarker - 1

         !recalculate the angle for the vertex to the left of the dropped one:
         CALL angleDiff(vec_timestamps(currVerts(endVertMarker-2:endVertMarker)), &
                      & yscaled(currVerts(endVertMarker-2:endVertMarker)), &
                        & yrScaled, ltr_distwtfactor, angles(endAngMarker))

      ELSEIF (minVertInd .EQ. 2) THEN
         ! drop this vertex. update the endVertMarker
         currVerts(2:endVertMarker-1) = currVerts(3:endVertMarker) 
         endVertMarker = endVertMarker - 1

         ! drop this angle. update the endAngMarker
         angles(minAngInd : endAngMarker-1) = angles(minAngInd + 1 : endAngMarker)
         endAngMarker = endAngMarker - 1

         !recalculate the angle to the right of the one taken out:
         CALL angleDiff(vec_timestamps(currVerts(1:3)),  &
                        yscaled(currVerts(1:3)), yrScaled, &
                        ltr_distwtfactor, angles(1))

      ELSE
         ! drop this vertex. update the endVertMarker
!         if (minVertInd >= endVertMarker) then
!                 print *,"bad condition"
!                 stop
!         end if 
!         if (minVertInd >= nVerts) then
!                 print *,"bad condition 1"
!                 stop
!                 
!         end if 
!         print *, 'minVertInd =', minVertInd
!         print *, 'endVertMarker = ', endVertMarker
!         print *, 'minAngInd', minAngInd
!         print *, 'endAngMarker =', endAngMarker
!         print *, '            '
         currVerts(minVertInd:endVertMarker-1) = currVerts(minVertInd+1:endVertMarker) 
         endVertMarker = endVertMarker - 1
         ! drop this angle. update the endAngMarker
         angles(minAngInd : endAngMarker-1) = angles(minAngInd + 1 : endAngMarker)
         endAngMarker = endAngMarker - 1

         !recalculate the angle to the left of the one taken out:
         CALL angleDiff(vec_timestamps(currVerts(minVertInd-2 : ((minVertInd+1)-1))), &
                      & yscaled(currVerts(minVertInd-2 : ((minVertInd+1)-1))), &
                        & yrScaled, ltr_distwtfactor, angles(minAngInd-1))

         !recalculate the angle to the right of the one taken out:
!         print *, vec_timestamps( currVerts(minVertInd-1:((minVertInd+2)-1)) )
         CALL angleDiff(vec_timestamps( currVerts(minVertInd-1:((minVertInd+2)-1)) ), &
                      & yscaled( currVerts(minVertInd-1:((minVertInd+2)-1)) ), &
                      & yrScaled, ltr_distwtfactor, angles(minAngInd))

      ENDIF

   END DO
   nFinalVerts = endVertMarker
   finalVerts(1:nFinalVerts) = currVerts(1:endVertMarker)

    RETURN
END SUBROUTINE cullByAngleChange

SUBROUTINE anchoredRegression(xvals, yvals, numPts, yanchorval, fit, slope)
IMPLICIT NONE
  !the idl codes have this routine in helper directory.
  !do a simple least-squares regression, but anchor it so that the 
  !zeroth element has the value "yanchorval"

  REAL(KIND=R8), INTENT(IN) :: xvals(:), yvals(:), yanchorval
  INTEGER, INTENT(IN) :: numPts
  REAL(KIND=R8), DIMENSION(numPts) :: x, y
  REAL(KIND=R8) :: xy, xx
  REAL(kind=R8), INTENT(OUT) :: slope
  REAL(kind=R8), INTENT(INOUT) :: fit(:)

  x = xvals - xvals(1)
  y = yvals - yanchorval

  xy = DOT_PRODUCT(x, y)
  xx = DOT_PRODUCT(x, x)
  slope = xy/xx
  fit = slope * x + yanchorval

END SUBROUTINE anchoredRegression

SUBROUTINE pickBetterFit(vec_obs, fit1, fit2, choice)
IMPLICIT NONE
  !This assumes same number of parameters to have created the yfits
  REAL(KIND=R8), INTENT(IN) :: vec_obs(:), fit1(:), fit2(:)
  REAL(KIND=R8) :: mse1, mse2
  INTEGER :: choice

  mse1 = DOT_PRODUCT((vec_obs - fit1), (vec_obs - fit1))
  mse2 = DOT_PRODUCT((vec_obs - fit2), (vec_obs - fit2))
  IF (mse1 .LT. mse2) THEN
     choice = 1
  ELSE
     choice = 2
  ENDIF

END SUBROUTINE pickBetterFit

SUBROUTINE fillLine(x, xEndpts, yEndpts, filledLine, slope)
IMPLICIT NONE
  !xEndpts, yEndpts are just 2 element arrays consisting of the starting
  !and ending points of the desired line
  REAL(KIND=R8), INTENT(IN) :: x(:), xEndpts(:), yEndpts(:)   
  REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: xvals
  INTEGER :: i, maxIdx, minIdx
  REAL(KIND=R8), INTENT(OUT) :: slope
  REAL(kind=R8), INTENT(INOUT) :: filledLine(:)

  slope = (yEndpts(2) - yEndpts(1))/(xEndpts(2) - xEndpts(1))

  i = 1
  ! in case where there are multiple images per year
  DO WHILE (x(i) /= xEndpts(1))
        i = i+1
  END DO
  minIdx = i
  DO WHILE (x(i) /= xEndPts(2))
        i = i + 1
  END DO
  maxIdx = i

  ALLOCATE (xvals(maxidx - minidx + 1))
  xvals = x(minidx : maxidx)
  filledLine = (/ (yEndpts(1) + slope * (xvals(i) - xvals(1)), i = 1, size(xvals))/)
  DEALLOCATE (xvals)

END SUBROUTINE fillLine

SUBROUTINE findBestTrace(vec_timestamps, vec_obs, numObs, &
                         currVerts, nVerts, nSegs, bt, work_arr)
USE globalVars_par
IMPLICIT NONE
! Global parameters needed here:
! ltr_WA ... it is an inout array.
! The definition for bestTrace.
!
!Given set of vertices (x-vals), find the the combo of vertex y-vals
!that results in best fit for each segment x and y are the original values.
!vertices: is the list of vertices (in terms of array position, not the x-value.
!nVerts:   is the number of vertices -- passed in to avoid allocation.
!nSegs:    is the number of segments -- passed in just to save calc time.
!This is used only on the first run, with all of the segments.  From
!here on out, we just eliminate each one and calc the vals

  REAL(KIND=R8), INTENT(IN) :: vec_timestamps(:), vec_obs(:)
  INTEGER, INTENT(IN) :: currVerts(:)
  INTEGER, INTENT(IN) :: nSegs, nVerts, numObs
  TYPE(ltrWA) :: work_arr
  REAL(KIND=R8), DIMENSION(numObs) :: fillway, anchway
  INTEGER :: seg, choice, v1, v2, i
  REAL(KIND=R8) :: fillSlope, regressSlope, anchSlope, const
  TYPE(bestTrace), INTENT(INOUT) :: bt 

  !TODO:  (becuz there was 'TODO' written in the idl code)

  ! initialize bt%vertYvals to be same as vec_obs
  bt%vertYvals(1:nVerts) = (/ (vec_obs(currVerts(i)), i=1,nVerts) /)
  ! linear fit in the interval [t_1, t_{v_1}]
  seg = 1
  v1 = currVerts(1)  !which should be 1, anyways.
  v2 = currVerts(2)
  CALL fillLine( vec_timestamps(v1:v2), &
                 (/ vec_timestamps(v1), vec_timestamps(v2) /),  &
                 (/ bt%vertYvals(1), bt%vertYvals(2) /), fillway(v1:v2), fillSlope)

  !just using anchway to store output here in order to avoid 
  !allocating one more array.
  CALL regress(vec_timestamps(v1:v2), vec_obs(v1:v2), work_arr, &
                            anchway(v1:v2), regressSlope, const)
  CALL pickBetterFit(vec_obs(v1:v2), fillway(v1:v2), anchway(v1:v2), choice)

  IF (choice .EQ. 1) THEN
     ! at the vertices, yvals resulting from dotway will be the same as original yvals
     bt%vertYVals(seg) = fillway(v1)
     bt%vertYVals(seg+1) = fillway(v2)
     bt%slopes(seg) = fillSlope
     bt%yFitVals(v1:v2) = fillway(v1:v2)
  ELSE
     ! yvals resulting from regression may be different from the original yvals at the vertices
     bt%vertYVals(seg) = anchway(v1)
     bt%vertYVals(seg+1) = anchway(v2)
     bt%slopes(seg) = regressSlope
     bt%yFitVals(v1:v2) = anchway(v1:v2)
  END IF

  DO seg = 2, nSegs 
     v1 = currVerts(seg)
     v2 = currVerts(seg+1)
     CALL fillLine(vec_timestamps(v1:v2),  &
                   (/ vec_timestamps(v1), vec_timestamps(v2) /), &
                   (/ bt%vertYVals(seg), bt%vertYVals(seg+1) /),  &
                   fillway(v1:v2), fillSlope)

     CALL anchoredRegression(vec_timestamps(v1:v2), vec_obs(v1:v2), &
                              v2-v1+1, bt%vertYvals(seg), anchway(v1:v2), anchSlope)
     CALL pickBetterFit(vec_obs(v1:v2), fillway(v1:v2), anchway(v1:v2), choice)
     IF (choice .EQ. 1) THEN
        ! at the vertices, yvals resulting from dotway will be the same
        ! as original yvals
        bt%yFitVals(v1+1:v2) = fillway(v1+1:v2)
        ! doesn't quite need to be updated in this case; will stay same only.
        !bt%vertYVals(seg+1) = fillway(v2) 
        bt%slopes(seg) = fillSlope
     ELSE
        ! yvals resulting from regression may be different from the 
        ! original yvals at the vertices
        bt%yFitVals(v1+1:v2) = anchway(v1+1:v2)
        bt%vertYVals(seg+1) = anchway(v2)
        bt%slopes(seg) = anchSlope
     END IF
  END DO
  bt%vertices(1:nVerts) = currVerts(1:nVerts)

END SUBROUTINE findBestTrace

SUBROUTINE calcFittingStats(vec_obs, vec_fitted, nObs, nParams, modelStats)
USE globalVars_par
USE stat_utils
IMPLICIT NONE
  ! Global parameters needed here:
  ! The definition of the derived type fittingStats_ltr
  REAL(KIND=R8), INTENT(IN) :: vec_obs(:), vec_fitted(:)
  INTEGER, INTENT(IN) :: nParams, nObs 
  REAL(KIND=R8), DIMENSION(nObs) :: observed, fitted
  INTEGER  ::  nPred, ok, ierr
  REAL(KIND=R8) :: ubar, ss_mean, abs_diff, ss_resid, X1_squared, X2_squared
  REAL(KIND=R8) :: infinity, AICC, AIC, ms_resid, ms_regr, f !, fRegr 
  REAL(KIND=R8) :: total_variance, residual_variance, adjusted_rsquared
  REAL :: Ix, pOfF
  INTEGER :: dof1, dof2
  TYPE(fittingStats_ltr), INTENT(OUT) :: modelStats
!  REAL :: HUGE
!  EXTERNAL :: HUGE
!   INTERFACE
!      SUBROUTINE bratio( a, b, x, y, w, w1, ierr )
!      USE REAL_PRECISION
!      real :: a, b, x, y, w, w1
!      integer :: ierr
!      END SUBROUTINE bratio
!   END INTERFACE

  infinity = 1000000000000.0000005   !HUGE(1000000000000.0000005)
  nPred = SIZE(vec_fitted)
  IF (nObs /= nPred) THEN
     ok = 0
     !TODO: modelStats components?
  ENDIF
  
  observed = vec_obs
  fitted = vec_fitted
 
!  !first, take out the points in the record where the observed and predicted are
!  !identical -- these are nodes in the fitting that should not be allowed to 
!                inflate the fit
!  ctr = 0
!  DO i = 1, nObs
!     IF ((observed(i) - fitted(i)) .LT. 0.00001) THEN
!        ctr = ctr + 1
!        exacts(ctr) = 1
!     ENDIF
!  END DO
!  nExacts = ctr

  ubar = SUM(observed)/nObs               ! mean in 'observed' group
  ss_mean = DOT_PRODUCT((observed - ubar), (observed - ubar))
  abs_diff = SUM(ABS(observed - fitted))
  ss_resid = DOT_PRODUCT((observed - fitted), (observed - fitted))

  X1_squared = ss_mean - ss_resid
  X2_squared = ss_resid

  dof1 = nParams  
  dof2 = nObs - (nParams + 1)   
  
  IF ((dof2 .LE. 0) .or. (X2_squared < 0.000001) .or. (dof1 ==0 )) THEN 
        modelStats%ok = 0
        modelStats%mean_u = ubar
        modelStats%sum_of_squares = ss_mean
        modelStats%sum_of_squares_resid = X2_squared
        modelStats%sum_of_squares_regressor = X1_squared
        modelStats%dof1 = dof1
        modelStats%dof2 = dof2
        modelStats%fStat = 0
        modelStats%p_of_f = 1.0
        modelStats%aicc = 0
        modelStats%residual_variance = 0
        modelStats%total_variance = 0
        modelStats%adjusted_rsquare = 0
        modelStats%ms_resid = 0
        modelStats%ms_regr = 0            ! ms denotes mean squared value
        modelStats%yfit = vec_fitted
        modelStats%abs_diff = abs_diff
  ELSE
        residual_variance = ss_resid/dof2
        total_variance = ss_mean/(nObs - 1)
        adjusted_rsquared = 1 - (residual_variance/total_variance)

        ms_regr = X1_squared/dof1
        IF ((ms_regr .LT. 0.00001 ) .or. (X2_squared < 0.000001) & 
            .or. (dof1 == 0))  THEN
           f = 0.0001
        ELSE
           f = X1_squared * dof2 / ( X2_squared * dof1 )  !(X1_squared/dof1)/(X2_squared/dof2)
        ENDIF

        ms_resid = X2_squared/dof2
        ! Get the probability that F > f, i.e., Q(f| d1, d2)
        ! Ix is the ratio of incomplete beta and complete beta.
        ! It can also be calculated using DBETAI.f etc
        ! but BRATIO uses a newer algorithm which is better.
        CALL BRATIO(REAL(0.5*dof2), REAL(0.5*dof1), &
                    REAL(dof2/(dof2 + dof1*f)),  &
                    REAL(1.0 - (dof2/(dof2 + dof1*f))), &
                    Ix, pOfF, ierr)
!        CALL BRATIO(REAL(1.0), REAL(2.0), &
!                    REAL(0.25),  &
!                    REAL(0.75), &
!                    Ix, pOfF, ierr)

        AIC = (2.0 * nParams) + (nObs * ALOG(REAL(X2_squared/nObs)))
        AICC = AIC + ((2.0 * nParams *(nParams+1))/(nObs - nParams - 1))
        IF (AICC > infinity) THEN
           AICC = -1
        ENDIF
        modelStats%ok = 1
        modelStats%mean_u = ubar
        modelStats%sum_of_squares = ss_mean
        modelStats%sum_of_squares_resid = X2_squared
        modelStats%sum_of_squares_regressor = X1_squared
        modelStats%dof1 = dof1
        modelStats%dof2 = dof2
        modelStats%fstat = f
        modelStats%p_of_f = pOfF
        modelStats%AICC = AICC
        modelStats%residual_variance = residual_variance
        modelStats%total_variance = total_variance
        modelStats%adjusted_rsquare = adjusted_rsquared !terms from Jongman et al. pg 37
        modelStats%ms_regr = ms_regr    
        modelStats%ms_resid = ms_resid   
        modelStats%yfit = vec_fitted
        modelStats%abs_diff = abs_diff
  ENDIF  

END SUBROUTINE calcFittingStats


SUBROUTINE fillFromVertices(x, nVerts, vertices, vertVals, yfit, slopes)
IMPLICIT NONE
    
    REAL(KIND=R8), INTENT(IN) :: x(:), vertVals(:)
    INTEGER, INTENT(IN) :: nVerts, vertices(:)
    INTEGER :: i, nSegments
    REAL(KIND=R8), INTENT(INOUT) :: yfit(:), slopes(:)
!fillLine(x, xEndpts, yEndpts, filledLine, slope)
    nSegments = nVerts - 1
    DO i = 1, nSegments
       CALL fillLine(x, (/x(vertices(i)), x(vertices(i+1))/), & 
                    (/vertVals(i), vertVals(i+1)/), &
                    yfit(vertices(i):vertices(i+1)), slopes(i))
    END DO

END SUBROUTINE fillFromVertices

SUBROUTINE takeOutWeakest(currModel, threshold, vec_timestamps, vec_obs, v, &
                vertYVals, nObs, nCurrVerts, nUpdatedVerts, updatedVerts)
USE globalVars_par
IMPLICIT NONE
   ! Global parameters needed here:
   ! The definition for the derived type ltr_model
   !
   ! currModel is the model being evaluated. It has:
   !      -- yfit:       fitted values over all timepoints
   !      -- vertices:   global indices (in present data?) of vertex positions
   !      -- fstat:      F-statistic of this model
   !      -- p_of_f:     
   !      -- aicc:
   !      -- vertYVals:  fitted values at vertices
   !      -- slopes:     slope of each segment in this model
   !
   ! nCurrVerts is the number of vertices in the current (incoming) model.
   ! vec_timestamps is the x-coordinate vector
   ! vec_obs is the y-coordinate vector
   ! nObs is the number of timepoints
   ! v is the set of vertices in the incoming model, same as currModel%vertices
   ! vertYVals is same as currModel%vertYVals
   ! nCurrVerts is len(v)
   ! Exactly 1 vertex will get dropped. 
   ! updatedVerts has nCurrVerts-1 elements.
   !
   ! We will only USE the currModel here. We won't update it. It gets updated
   ! in the calling program. Here, only the updated vertices are gathered
   ! and returned as the vector updatedVerts.

   TYPE(ltr_model), INTENT(IN) :: currModel !will get allocated and calculated in
                                        !the calling program when it calls 
                                        !calcFittingStats. Note that, currModel 
                                        !will just an intersecting set of
                                        !fittingStats and bestTrace. 
   REAL(KIND=R8), INTENT(IN) :: threshold, vec_timestamps(:), vec_obs(:)
   INTEGER, INTENT(IN) ::  nCurrVerts, nObs, v(:)
   INTEGER :: nSlopes, nNegatives, i , k
   INTEGER, DIMENSION(nCurrVerts-1):: negatives
   REAL(KIND=R8),DIMENSION(nCurrVerts-1) :: scaledSlopes !no. segs is one less than no. verts
   REAL(KIND=R8) :: range_of_values
   LOGICAL :: runMse
   REAL(KIND=R8), DIMENSION(nCurrVerts):: MSE !One element per vertex-drop 
                                        !MSE is actually calculated only for interior vertices.
                                        !The first two indices are set to a huge number.
   INTEGER :: j(1:1), weakest_segment, weakest_vert
   REAL(KIND=R8), DIMENSION(nObs) :: yfit, uPP
   INTEGER :: vleft, vright, nUpdatedVerts
   REAL(KIND=R8) :: leftx, lefty, rightx, righty, fillSlope
   INTEGER, INTENT(INOUT) :: updatedVerts(:)
   REAL(KIND=R8), INTENT(INOUT) :: vertYVals(:)

   nSlopes = currModel%numVertices - 1
   yfit = currModel%yfit
   ! yfit is used for 
   ! (i)  calculating the range_of_vals
   ! (ii) IF no seg is found, then for calculating the MSEs.
   !we operate under the knowledge that disturbance is always considered to
   !have a negative (positive) slope, and recovery a positive (negative) slope
   !based on NDVI (band5) type indicators).
   nNegatives = 0
   scaledSlopes = 0.0
   !print *, 'nSplopes =', nSlopes
   DO i = 1, nSlopes
!      IF ((currModel%slopes(i) < 0).AND.(currModel%slopes(i) /= -1)) THEN
      IF ((currModel%slopes(i) > 0).AND.(-1.0*currModel%slopes(i) /= -1)) THEN
         nNegatives = nNegatives + 1
         negatives(nNegatives) = i
      ENDIF
   END DO
   ! only negatives(1:nNegatives) will be relevant. Rest is garbage.
   range_of_values = MAXVAL(yfit) - MINVAL(yfit)
   
   ! again, only scaledSlopes(1:nNegatives) will be relevant. Rest is garbage.
   IF (nNegatives .GT. 0) THEN
      scaledSlopes(1:nNegatives) = &
              ABS(currModel%slopes(negatives(1:nNegatives)))/range_of_values
   ELSE
      scaledSlopes(1:nNegatives) = threshold -1 !set so it won't be GT threshold
   ENDIF
   runMSE = .TRUE.

   !If none of the segments is neg, then scSlp is all less that tau, i.e.,
   !maxval(scSlp) < tau. So, basically, this if-block will execute only if
   !there is atleast one segment with negative (and /= -1) slope, is it?
   !Summarizing, if a neg slope is found, we choose that vertex as the one
   !to be dropped.
   IF (MAXVAL(scaledSlopes(1:nNegatives)).GT.threshold) THEN
        !Note that scaled slopes has the absolute values of negative slopes.
        j = MAXLOC(scaledSlopes(1:nNegatives)) 
        weakest_segment = negatives(j(1))
        weakest_vert = weakest_segment + 1 
        !The violator is a segment -- which vertex to remove? Since we are
        !tracking through time we assume that it is the latter vertex that
        !is causing the problem and take it out. This will be violated only
        !if there are spikes in brightness that are not fixed by desawtooth,
        !but that situation would be no better removing the first vertex
        !anyway, so we stick with this approach since it will take out more
        !shadow problems. the only now interpolate to get rid of this point,
        !so it doesn't mess up the fits later.
        IF (weakest_vert .EQ. nCurrVerts) THEN
           !since the violating point was at end, need to run mse instead
           !after fixing
           yfit(v(weakest_vert)) = yfit(v(weakest_vert)-1)
           vertYVals(weakest_vert) = yfit(v(weakest_vert)) !These two dont make sense either
           runMSE = .TRUE.
        ELSE 
           ! No need for MSE in this case. Just drop the vertex.

!           thisx = vec_timestamps(v(weakest_vert))
!           leftx = vec_timestamps(v(weakest_vert-1))
!           rightx = vec_timestamps(v(weakest_vert+1))
!           lefty = currModel%vertYVals(weakest_vert-1)
!           righty = currModel%vertYVals(weakest_vert+1)
!           tmp = (righty - lefty)/(rightx - leftx)  !slope 
!           currModel%yfit(v(weakest_vert)) = (thisx-leftx)*tmp + lefty
           !the above calculation make no sense!
           !why are we changing the fit values in the current model????!!!!!
           !the idea is to generate an updated set of vertices
           !and then a new (simpler) model using those vertices.
           !that happens in findBestTrace, after exiting this routine.
           !There's no sense in altering the current model.

           updatedVerts(1:weakest_vert-1) = v(1:weakest_vert-1)
           updatedVerts(weakest_vert: nCurrVerts-1) = v(weakest_vert+1: nCurrVerts)
           ! So, MSE will not be calculated anymore.
           runMSE = .FALSE.
        ENDIF
   ENDIF

   !If no neg. segment is found, or, the neg. segment cors to the last vertex, then
   !choose the vertex to be dropped based on MSE.
   IF (runMSE .EQV. .TRUE.) THEN
  
     MSE(1) = 1000000000000.0000005   !HUGE(1000000000000.0000005)
     ! actually, can't we just have MSE of length nCurrVerts-2?
     MSE(nCurrVerts) = 1000000000000.0000005   !HUGE(1000000000000.0000005)
     DO i = 2, nCurrVerts - 1
         vleft = v(i-1)
         vright = v(i+1)
         leftx =  vec_timestamps(vleft)
         rightx = vec_timestamps(vright)
         lefty = yfit(vleft)   !per pseudocode, these shd be the fitted values
         righty = yfit(vright)  !at these vertices that were obtained in modelFit.
         ! calculate the segment that will result from dropping vertex i.
         ! the new values of the potential fit are stored in yfit itself
         ! at the corresponding indices.
         CALL fillLine(vec_timestamps(vleft:vright), &
                      (/ leftx, rightx /), (/ lefty, righty /), &
                      uPP(vleft:vright), fillSlope)
         MSE(i) = SUM((/ ((uPP(k) - vec_obs(k))* (uPP(k) - vec_obs(k))/(rightx-leftx), k = vleft,vright) /))
     END DO

     j = MINLOC(MSE(2:nCurrVerts-1)) !min location wrt indices of original stencil will be j + 1 
     weakest_vert = j(1) + 1
     updatedVerts(1:weakest_vert-1) = v(1:weakest_vert-1)
     updatedVerts(weakest_vert: nCurrVerts-1) = v(weakest_vert+1: nCurrVerts)
   ENDIF
   nUpdatedVerts = nCurrVerts - 1 ! cud we possibly need any check for this?
   !we will do newYFit later after getting out of this subroutine.
   !Separation of responsibilities: This routine will only generate
   !the new vertices. A separate subroutine (findBestTrace) will 
   !generate the new model using the updated vertices.
   !I also think that the IDL code is very buggy.
   !CALL fillFromVertices(vec_timestamps, nRemainVerts, remainVerts, &
   !                        remVertVals, newYFit, newSlopes)

END SUBROUTINE takeOutWeakest

SUBROUTINE pickBestModel(my_models, numModels, bestModelInd)
USE globalVars_par
IMPLICIT NONE
   ! Global variables needed:
   ! ltr_model definition
   ! ltr_useFstat 
   ! ltr_bestModelProp
   ! ltr_pval
   ! currModel is the model being evaluated. It has:
   !      -- yfit:       fitted values over all timepoints
   !      -- vertices:   global indices (in present data?) of vertex positions
   !      -- fstat:      F-statistic of this model
   !      -- p_of_f:     
   !      -- aicc:
   !      -- vertYVals:  fitted values at vertices
   !      -- slopes:     slope of each segment in this model
   !
   TYPE(ltr_model), INTENT(IN) :: my_models(:)
   INTEGER, INTENT(IN) :: numModels
   INTEGER :: ctr, passTheTest, i
   REAL(KIND=R8) :: tau, mn, mx
   INTEGER, DIMENSION(numModels) :: goodModelsInd
   INTEGER, INTENT(OUT) :: bestModelInd
   
   !Corner case 1: check on useFstat. Not sure why it's needed, though.
   !I don't see myself NOT supplying the value of this parameter.
!   IF (useFstat .EQ. null) THEN
!      useFstat = 0
!   ENDIF
   
   !Corner case 2: ensure that there is at least one model with 1 or more segment.
   !This also looks redundant to me.
   ctr = 0
   DO i = 1, numModels
      IF (SIZE(my_models(i)%slopes) /= 0) THEN
         ctr = ctr + 1
      ENDIF
   END DO
   IF (ctr .EQ. 0) THEN
      STOP
   ENDIF

   !Now pick the best one
   IF (ltr_useFstat == 0) THEN
      !This is how we ideally want to choose our model. F-statistic in
      !comibination with it's p-value is more indicative.
      !Want to settle on a model that has low p of F-statistic
      mn = MINVAL((/ (my_models(i)%p_of_f, i=1, numModels) /))
      tau = (2.0_R8 - ltr_bestModelProp) * mn
      ctr = 0
      DO i = 1, numModels
         IF (my_models(i)%p_of_f .LE. tau) THEN
            ctr = ctr + 1
            goodModelsInd(ctr) = i
         ENDIF
      END DO
      bestModelInd = goodModelsInd(1)
      !Question: If we do have to take the first instance of a good model as
      !          our bestModel, then why are we collecting indices of all the
      !          good models? As soon as the first good model is encountered,
      !          we can just accept that and return!
   ELSE
      ! Compromise: We'll make a choice based on F-statistic only. But, at least, we make
      ! sure that we have some models that have p of F-statistic lower that a
      ! certain threshold (pval).
      passTheTest = 0
      DO i = 1, numModels
         IF (my_models(i)%p_of_f .LE. ltr_pval) THEN
             passTheTest = passTheTest + 1
         ENDIF
      END DO
      !Corner case
      IF (passTheTest .EQ. 0) THEN
         bestModelInd = -1
         RETURN
      ENDIF
      !Now make a decision based on F-statistic only. What to do? :(
      mx = MAXVAL( (/  (my_models(i)%fstat, i=1,numModels ) /))
      tau = ltr_bestModelProp * mx
      DO i=1,numModels
         IF (my_models(i)%fstat .GT. tau) THEN
            ctr = ctr + 1
            goodModelsInd(ctr) = i
         ENDIF
      END DO
      bestModelInd = goodModelsInd(1)
   ENDIF

END SUBROUTINE pickBestModel

SUBROUTINE checkSlopes(model, nVerts, accept)
USE globalVars_par
IMPLICIT NONE
   ! Global variables needed:
   ! ltr_model
   ! ltr_recoveryThreshold
   !
   ! model is the model being evaluated. It has:
   !      -- numVerts:   number of vertices in the model
   !      -- yfit:       fitted values over all timepoints
   !      -- vertices:   global indices (in present data?) of vertex positions
   !      -- fstat:      F-statistic of this model
   !      -- p_of_f:     
   !      -- aicc:
   !      -- vertYVals:  fitted values at vertices
   !      -- slopes:     slope of each segment in this model
   !
   !Given one model, look at its slopes.
   !Filter out if recovery happens quicker than quickest disturbance --
   !a value-free way to get unreasonable things out.
   !but of course don't do if all we have is recovery. No way to tell for sure then.
   TYPE(ltr_model), INTENT(IN) :: model
   INTEGER, INTENT(IN) :: nVerts
   INTEGER ::  i, nPositives, nSlopes
   INTEGER, DIMENSION(nVerts-1) :: positives
   REAL(KIND=R8) :: rangeOfVals
   REAL(KIND=R8), DIMENSION(nVerts-1) :: scaledSlopes
   LOGICAL, INTENT(OUT) :: accept

   accept = .TRUE.  ! 1 is for yes, 0 is for no
   nSlopes = nVerts -1
   !all of these operate under the knowledge that disturbance is
   !always considered to have a positive (should be negative for NDVI)
   !slope, and recovery a negative slope (based on band5 type indicators).
   !Always make sure that the adjustment factor that happens
   !upstream of find_segments6 ensures that recovery is negative
   rangeOfVals = MAXVAL(model%yfit) - MINVAL(model%yfit)
   nPositives = 0
   DO i=1,nSlopes
      IF (model%slopes(i) > 0) THEN
         nPositives = nPositives + 1
         positives(nPositives) = i
      ENDIF
   END DO
   IF (nPositives .GT. 0) THEN 
      scaledSlopes(1:nPositives) = ABS(model%slopes(positives(1:nPositives)))/rangeOfVals
      IF (MAXVAL(scaledSlopes(1:nPositives)) .GT. ltr_recoveryThreshold) THEN
          accept = .FALSE.
          RETURN
      ENDIF
   ENDIF

END SUBROUTINE checkSlopes

SUBROUTINE findBestTrace_alternate(vec_timestamps, vec_obs, Sfinal, &
                                   currVerts, numCurrVerts, &
                                   order, n_dim, bt)
USE globalVars_par
USE bsplines
IMPLICIT NONE
     ! Global variables needed:
     ! bestTrace definition
     ! 
!    this routine is for the marquardt approach. Err ... actually, 
!    we are using Bsplines
!    "x values"
!    "y values"
!    "currVerts"
!    "x coords of vertices"
!    "y coords of vertices"
!    "y-fit values"
!    "slopes is needed for later analysis in the main algorithm. So 
!     we store and return those as well"

    REAL(KIND=R8), INTENT(IN) :: vec_timestamps(:), vec_obs(:)
    INTEGER, INTENT(IN) :: currVerts(:)
    INTEGER, INTENT(IN) :: order, n_dim, Sfinal, numCurrVerts
    integer :: i, thisVert, nextVert, left
    !INTEGER  :: order, n_dim, Sfinal, numCurrVerts
    real, dimension(Sfinal) :: vec_ts, vec_obsr
    REAL, DIMENSION(order, n_dim) :: Q
    REAL, DIMENSION(n_dim) :: DIAG, bcoeff
    REAL, DIMENSION(Sfinal) :: weights
    REAL(kind=r8), DIMENSION(Sfinal) :: fit
    REAL, DIMENSION(numCurrVerts+2) :: knots
!    REAL, DIMENSION(100) :: weights
!    REAL, DIMENSION(7) :: knots
!    REAL, DIMENSION(:, :), allocatable :: Q
!    REAL, DIMENSION(:), allocatable :: DIAG, bcoeff
    type(bestTrace), intent(inout) :: bt


    ! TODO: figure out the thing between out knot sequence here vs n_dim+K
    vec_ts = (/ (real(vec_timestamps(i), kind=4), i=1, Sfinal ) /)
    vec_obsr = (/ (real(vec_obs(i), kind=4), i=1, Sfinal) /)
    knots(1) = vec_ts(1)
    knots(2:numCurrVerts+1) = vec_ts(currVerts(1:numCurrVerts)) 
    knots(numCurrVerts+2) = vec_ts(Sfinal)
    !order = 2
    !n_dim = 5
    !Sfinal = 100
    !vec_ts = (/ (real(i), i=1, Sfinal ) /)
    !vec_obsr = (/ (3.0*sin(2.0*3.148*real(i)/10.0) + &
    !               2.0*cos(2.0*3.148*real(i)/10.0) + &
    !               0.5*real(i/10.0), i=1,100)  /)
!
!    knots(1) = 1.0
!    knots(2) = 1.0
!    knots(3) = 25.0
!    knots(4) = 50.0
!    knots(5) = 75.0
!    knots(6) = 100.0
!    knots(7) = 100.0
    
!    allocate (Q(order, n_dim), diag(n_dim), bcoeff(n_dim))
    weights = 1.0
    Q = 0.0
    DIAG = 0.0
    !print *, knots
    CALL l2appr(knots, n_dim, order, Sfinal, vec_ts, vec_obsr, &
                     weights, Q, DIAG, bcoeff)
    !print *, 'bcoeffs:' , bcoeff
    fit = 0.0
    left = order
    tsloop: DO i = 1, Sfinal
       do while ((vec_ts(i) >= knots(left+1)) .or. (knots(left-1) &
                                                  >= knots(left+1)))
           if (left >= n_dim) then
              fit(i) = 0.0
              cycle tsloop
           endif
           left = left + 1
       end do
       ! bvalue will return a single precision number
       fit(i) =  real(bvalue(knots, bcoeff, n_dim, order, &
                    vec_ts(i), 0, left, 0), kind=r8)
    END DO tsloop
    ! figure out n_dim vs Sfinal
    fit(Sfinal) = fit(Sfinal-1) !TODO: this is cheating. Figure out if that
                                ! last breakpoint processing can be safely
                                ! updated to use left = t_{n_dim-1} 
    bt%vertices(1:numCurrVerts) = currVerts(1:numCurrVerts)
    DO i=1,numCurrVerts-1
       thisVert = currVerts(i)
       nextVert = currVerts(i+1)
       bt%vertYVals(i) = fit(thisVert)
       bt%slopes(i) = (fit(nextVert) - fit(thisVert))/  & 
                     (vec_timestamps(nextVert) - vec_timestamps(thisVert))
    END DO
    bt%vertYVals(numCurrVerts) = fit(currVerts(numCurrVerts))
    bt%yFitVals(1:Sfinal) = fit(1:Sfinal)

END SUBROUTINE findBestTrace_alternate

SUBROUTINE takeOutWeakest_alternate(vec_timestamps, vec_obs, nObs, v, &
                                    vertVals, nCurrVerts, updatedVerts)

    ! vertVals are the fitted y-values at the vertices
    REAL(KIND=R8), INTENT(IN) :: vec_timestamps(:), vec_obs(:), vertVals(:)
    INTEGER, INTENT(IN) :: v(:), nCurrVerts, nObs
    REAL(KIND=R8), DIMENSION(nCurrVerts) :: MSE !there will be one value per vertex (drop)
                                            !The very first and the very last MSE values
                                            !are set to HUGE cuz we dont want to drop
                                            !these two vertices.
    REAL(KIND=R8), DIMENSION(nObs) :: uPP
    INTEGER :: i, j(1:1), weakest_vert, vleft, vright
    REAL(KIND=R8) :: fillSlope, leftx, rightx, lefty, righty
    INTEGER, INTENT(INOUT) :: updatedVerts(:)

    MSE(1) = 1000000000000.0000005   !HUGE(1000000000000.0000005)
    MSE(nCurrVerts) = 1000000000000.0000005   !HUGE(1000000000000.0000005)
    DO i = 2, nCurrVerts - 1
        vleft = v(i-1)
        vright = v(i+1)
        leftx =  vec_timestamps(vleft)
        rightx = vec_timestamps(vright)
        lefty = vertVals(i-1)   !per pseudocode, these shd be the fitted values
        righty = vertVals(i+1)  !at these vertices that were obtained in modelFit.
        ! calculate the segment that will result from dropping vertex i.
        ! the new values of the potential fit are stored in uPP itself
        ! at the corresponding indices.
        CALL fillLine(vec_timestamps(vleft:vright), &
                       (/ vec_timestamps(vleft), vec_timestamps(vright) /), &
                       (/ lefty, righty /), &
                     uPP(vleft:vright), fillSlope)


        MSE(i) = DOT_PRODUCT((uPP(vleft:vright) - vec_obs(vleft:vright)), &
                             (uPP(vleft:vright) - vec_obs(vleft:vright)))
    END DO
    j = MINLOC(MSE(2:nCurrVerts-1)) !min location wrt indices of original stencil will be j + 1 
    weakest_vert = j(1) + 1
    updatedVerts(1:weakest_vert-1) = v(1:weakest_vert-1)
    updatedVerts(weakest_vert: nCurrVerts-1) = v(weakest_vert+1: nCurrVerts)
    
END SUBROUTINE takeOutWeakest_alternate

END MODULE utilities_ltr

