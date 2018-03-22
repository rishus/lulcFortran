MODULE ewmacd_fn_par
CONTAINS

SUBROUTINE ewmacd(imgIndex)  !pixel_x is a k-element vector. So is pixel_y.
USE globalVars_par
USE utilities_ewma
USE REAL_PRECISION
IMPLICIT NONE

! ************* input vars **********************
INTEGER(kind = 8), intent(in) :: imgIndex(:)
! ***********************************************

! ************* utility vars ********************
integer (kind=8) :: num_pixels, pixel, pixel_x, pixel_y
INTEGER :: i, ncols, len_training_vec, Sfinal, pctr
INTEGER :: ierr_lsfit, aerr, derr, sIbar
REAL(KIND=R8), DIMENSION(num_obs) ::  t, D
REAL(KIND=R8)  :: sigmaIhat  
REAL(KIND=R8), DIMENSION(2*ewmacd_numHarmonics+1) :: alphaStar
REAL(KIND=R8), DIMENSION(num_obs) :: EstarAlphastar, tau, z, f, persistenceVec
INTEGER, DIMENSION(num_obs) :: Ibar
INTEGER (KIND = 2), DIMENSION(num_obs) :: jump_vals_presSten, presInd
!REAL (KIND = 4), DIMENSION(num_obs) :: residuals_presSten  ! uncomment only residuals are needed
LOGICAL :: nullFlag
type(ewmaWA) :: work_arr
character(256) :: my_errmsg
! ***********************************************

! ********* PURPOSE **************
!   Main code for algorithm EWMACD
!   Supporting files are all in utilities_ewma.f90
!   Output is written to the array ewma_summary
! ***********************************************

num_pixels = SIZE(imgIndex, 1)
ncols = 2*ewmacd_numHarmonics+1

ALLOCATE (work_arr%tmp_mat(num_obs, 2*ewmacd_numHarmonics+1), &
               &  work_arr%rtmp_vec1(num_obs), work_arr%rtmp_vec2(num_obs),   &
               &  work_arr%itmp_vec1(num_obs), work_arr%itmp_vec2(num_obs), &
               &  work_arr%dev(num_obs), STAT=aerr)

DO pixel = 1, num_pixels

     pixel_x = mod(imgIndex(pixel), NumCols)
     IF (pixel_x == 0) THEN
         pixel_x = NumCols
     ENDIF
     pixel_y = INT(CEILING(REAL(imgIndex(pixel))/NumCols))
    
    ! ************* develop the presInd vector ***********************
    pctr = 0
    presInd = 0
    len_training_vec = 0
    !IF (reverse .eqv. .FALSE.) THEN
         DO i = 1, num_obs
            IF (input_mat(i, pixel_x, pixel_y) > ewmacd_lowthresh) THEN
               pctr = pctr + 1
               presInd(pctr) = int(i, kind=2)
               IF ((tyeardoy(i,1) >= trainingStart) .and. (tyeardoy(i,1) < trainingEnd)) THEN
                  len_training_vec = len_training_vec + 1
               ENDIF
            ENDIF
         END DO
    !ENDIF
    Sfinal = pctr

    IF (len_training_vec <= 2.0 * ewmacd_numHarmonics+1 ) THEN
       CYCLE
    ENDIF

    ! *************** prepare data ***********************************
    
    t(1:Sfinal) = (/ (tyeardoy(presInd(i),2)*2.0_R8*Pi/365.0_R8, i=1,Sfinal) /)
    
    D(1:Sfinal) = (/ (input_mat(presInd(i), pixel_x, pixel_y), i=1,Sfinal  ) /)
    
    nullFlag = .FALSE.
    ierr_lsfit = 0
    
    ! ****************** actual processing starts here ********************

    ! alphaStar gets allocated
    CALL LSFIT(t(1:len_training_vec), D(1:len_training_vec),   &
             &     ewmacd_numHarmonics, alphaStar, ierr_lsfit, work_arr)

    ! EstarAlphastar, Ibar get allocated; Ibar is size sz
    CALL getResiduals(alphaStar, t(1:Sfinal), D(1:Sfinal),  &
                      ewmacd_numHarmonics, len_training_vec, ierr_lsfit, sigmaIhat,  &
                      EstarAlphastar, Ibar, sIbar, nullFlag, work_arr)
    !only the first Sfinal elements of EstarAlphastar and first sIbar elements of Ibar are relevant

    IF (nullFlag .eqv. .FALSE.) THEN
    
       ! tau gets allocated; tau is size sz
       CALL getControlLimits(sigmaIhat, sIbar, tau(1:sIbar))
    
       ! z gets allocated; z is size sz
       CALL getEWMA(Ibar(1:sIbar), EstarAlphastar(1:Sfinal), z(1:sIbar), sIbar)
    
       ! f gets allocated; f is size sz
       CALL flagHistory(z(1:sIbar), tau(1:sIbar), f(1:sIbar), sIbar)
    
       ! persistenceVec gets allocated; again size sz
       CALL persistenceCounting(f(1:sIbar), persistenceVec(1:sIbar), work_arr, sIbar)
    
       jump_vals_presSten(1:Sfinal) = INT(mdv, KIND=2)  ! present data = 'good' data + 'outlier' data
       jump_vals_presSten(Ibar) = INT(persistenceVec(1:sIbar), KIND=2)  ! remember that Ibar is relative to 'Sfinal'
    
       CALL summarize(pixel_x, pixel_y, jump_vals_presSten(1:Sfinal),  &
                     &          presInd(1:Sfinal), 3, Sfinal, num_obs)
    
!       ! Uncomment IF residuals are desired
!       residuals_presSten(1:Sfinal) = mdv
!       DO i = 1, sIbar
!          residuals_presSten(Ibar(i)) = REAL(EstarAlphastar(Ibar(i)), KIND=4) ! EstarAlphastar is of size Sfinal
!       END DO                                                        ! However, we are retaining only the values 
!                                                                   ! cors to 'good' data. Outlier timestamps
!                                                                   ! will have -2222.
!       CALL summarize_residuals(pixel_x, pixel_y, residuals_presSten(1:Sfinal), presInd(1:Sfinal), 3, Sfinal)
!    
!       ! Uncomment IF coefficients are needed
!       ewma_coeffs(1:ncols, pixel_x, pixel_y) = (/ (REAL(alphaStar(i), KIND = 4), i= 1,ncols) /)
    
    ENDIF

END DO

DEALLOCATE (work_arr%tmp_mat, work_arr%rtmp_vec1, work_arr%rtmp_vec2,  &
            work_arr%itmp_vec1, work_arr%itmp_vec2, work_arr%dev, &
            STAT=derr, errmsg=my_errmsg)

END SUBROUTINE ewmacd

END MODULE ewmacd_fn_par
