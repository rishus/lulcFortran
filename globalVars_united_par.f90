MODULE globalVars_par
USE REAL_PRECISION
IMPLICIT NONE

!http://www.cs.uwm.edu/~cs151/Bacon/Lecture/HTML/ch06s09.html
! INTEGER(KIND=2): 16 bit (2byte) signed integer, -32768 to 32768
! INTEGER(KIND=4): 32 bit (4byte) signed integer, 2,147,483,648 to +2,147,483,647
! INTEGER(KIND=8): 64 bit (8byte) signed integer, +/-9.22 x 1018
! REAL(KIND=4)   : 32 bit (4byte ieee float) floating point,+/-(1.1754 x 10-38 to 3.4028 x 1038)
! REAL(KIND=8)   : 64 bit (8byte iee float) floating point,+/-(2.2250x10^-308 to 1.7976 x 10^308)
! REAL(KIND=R8)  : If R8=8, 64 bit (8byte ieee float), ensures what the computer supports. 
!                   (Make similar arrangement for integers also.)

! shared (or, read only) variables
! declared here. Assigned in the main program. Any subprogram
! that uses this module will have access to these variables.
!********* common ********************
INTEGER (KIND = 2), ALLOCATABLE :: input_mat(:,:,:), tyeardoy(:,:)  
LOGICAL :: reverse
INTEGER ::  num_obs
INTEGER (KIND = 8) :: NumRows, NumCols
REAL(KIND=R8), PARAMETER :: Pi = 3.141592653589793_R8

!*********for EWMACD **********************
REAL(KIND=4), ALLOCATABLE :: ewma_summary(:,:,:), ewma_residuals(:,:,:), ewma_coeffs(:,:,:)
INTEGER(KIND=2) :: ewmacd_numHarmonics
REAL(KIND=R8)  :: xbarlim1, xbarlim2, ewmacd_lowthresh, lam, L, mu
INTEGER :: trainingStart, trainingEnd, persistence
INTEGER :: nb

! work arrays. Must be threadprivate (or private?).
! as of now, I am declaring and allocating these in the main program only. 
! By this module, the subprogram will only have the definition of this 
! type bfastWA. The actual work_arr will have to be passed as an argument to it.
type ewmaWA
INTEGER, DIMENSION(:), ALLOCATABLE :: itmp_vec1, itmp_vec2
REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: rtmp_vec1, rtmp_vec2, dev
REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE :: tmp_mat
end type ewmaWA

!*********** for BFAST *******************
REAL(KIND=4), ALLOCATABLE :: bfast_summary(:,:,:)  
INTEGER(KIND=2) :: bfast_numHarmonics, bfast_numBrks, bfast_numColsProcess, bfast_maxIter
REAL(KIND=R8)  :: bfast_brkptSpacingProp, bfast_pvalThresh
INTEGER ::  nbDGELS
REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE :: BBI

type bfastWA
INTEGER, DIMENSION(:), ALLOCATABLE :: itmp_vec1, itmp_vec2
REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: rtmp_vec1, rtmp_vec2
REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE :: tmp_mat !may be do RWork_mat also here itself?!
end type bfastWA

!*********** for LANDTRENDR *******************
REAL(KIND=4), ALLOCATABLE :: ltr_summary(:,:,:), ltr_brkpt_summary(:,:,:)
INTEGER :: ltr_mu, ltr_nu, ltr_useFstat
REAL(KIND=R8) :: ltr_despikeTol, ltr_distwtfactor, ltr_recoveryThreshold
REAL(KIND=R8) :: ltr_bestModelProp, ltr_pval

type ltrWA
INTEGER, DIMENSION(:), ALLOCATABLE :: itmp_vec1, itmp_vec2
REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: rtmp_vec1, rtmp_vec2
REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE :: tmp_mat !may be do RWork_mat also here itself?!
end type ltrWA

type ltr_model
    REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: yfit
    INTEGER, DIMENSION(:), ALLOCATABLE :: vertices
    REAL(KIND=R8) :: fstat, p_of_f, aicc
    REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: vertYVals, slopes
    INTEGER :: numVertices
end type ltr_model

type fittingStats_ltr
    INTEGER ::  ok
    REAL(KIND=R8) :: mean_u
    REAL(KIND=R8) :: sum_of_squares 
    REAL(KIND=R8) :: sum_of_squares_resid
    REAL(KIND=R8) :: sum_of_squares_regressor
    INTEGER :: dof1
    INTEGER :: dof2
    REAL(KIND=R8) :: residual_variance 
    REAL(KIND=R8) :: total_variance 
    REAL(KIND=R8) ::  adjusted_rsquare 
    REAL(KIND=R8) :: fstat
    REAL(KIND=R8) :: p_of_f
    REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: yfit 
    REAL(KIND=R8) :: ms_regr 
    REAL(KIND=R8) :: ms_resid 
    REAL(KIND=R8) :: aicc
    REAL(KIND=R8) :: abs_diff 
end type fittingStats_ltr

TYPE bestTrace
    INTEGER, DIMENSION(:), ALLOCATABLE :: vertices
    ! This is an array. Use a pointer
    REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: vertYVals
    ! This is an array. Use a pointer 
    REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: yfitVals
    ! This is an array. Use a pointer(:) 
    REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: slopes
end type bestTrace



END MODULE globalVars_par
