program EWMACD_main
USE REAL_PRECISION
USE globalVars_par
!USE ewmacd_fn_par
!USE bfast_fn_par
USE ltr_fn_par
!$ USE OMP_LIB, ONLY: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM, OMP_GET_WTIME

        implicit none

        INTEGER, parameter :: INTKIND = 2, REALKIND = 4, S_LEN = 100, STRINGBUFFER = 500, NEG1 = -1
        INTEGER, parameter :: BAND_INFO_FIELDS = 3    ! (since num cols in info csv file is 3)
        INTEGER, parameter :: UTM_LINE_NUMBER = 16
        INTEGER, parameter :: CELL_SIZE = 30
        INTEGER, parameter :: HALF_CELL_SIZE = 15

        INTEGER, parameter :: INPUT_FILE = 11
        INTEGER, parameter :: INPUT_HEADER_FILE = 12
        INTEGER, parameter :: INPUT_DATES_FILE = 13
        INTEGER, parameter :: OUTPUT_FILE = 14
        INTEGER, parameter :: OUTPUT_FILE2 = 15
        INTEGER, parameter :: OUTPUT_FILE3 = 16
        
        INTEGER (KIND = 2) :: b    ! NumBands
        INTEGER (KIND = 2) :: CurrHeaderLineNum = 8
        INTEGER (KIND = 2) :: MemoryDivisor = 0
        INTEGER (KIND = 4) :: ULX = 0
        INTEGER (KIND = 4) :: ULY = 0
        INTEGER (KIND = 4) :: ULX_Center = 0
        INTEGER (KIND = 4) :: ULY_Center = 0
        INTEGER (KIND = 4) :: CurrEasting
        INTEGER (KIND = 4) :: CurrNorthing
       ! INTEGER (KIND = 8) :: NumRows, NumCols 
        INTEGER (KIND = 8) :: NumExtraRows
        INTEGER (KIND = 8) :: NumRowsChunk = 0
        INTEGER (KIND = 4) :: ix
        INTEGER (KIND = 8) :: i, pixel_x, pixel_y, tmp, imgIndex
        INTEGER (KIND = 8) :: NumPixels, numValidPixels, chunk_size

        INTEGER (KIND = 2), ALLOCATABLE :: input_mat_old(:,:,:), z (:,:)  ! for generating txt output file
        !INTEGER (KIND = 4) :: currInd, time
        INTEGER  :: numArgs, derr
        
        CHARACTER (LEN = 100), DIMENSION(:), ALLOCATABLE :: args
        CHARACTER (LEN = S_LEN) :: InputFileName
        CHARACTER (LEN = S_LEN) :: InputHeaderFileName
        CHARACTER (LEN = S_LEN) :: InputDatesFileName
        CHARACTER (LEN = S_LEN) :: OutputFileName, OutputFileName2
        CHARACTER (LEN = S_LEN) :: OutputFilePrefix, OutputFileName3
        CHARACTER (LEN = STRINGBUFFER) TempString, TempString3
        
        LOGICAL :: OS = .FALSE.
        LOGICAL :: outputFileExists
        REAL(KIND=R8) :: START, FINISH
        !INTEGER :: len_training_vec, Sfinal, pctr
        INTEGER :: pctr
!        INTEGER (KIND=2), DIMENSION(:), ALLOCATABLE :: presInd

!        type (myWA) :: work_arr
        integer :: ilaenv
        EXTERNAL ILAENV 

!        PRINT *, 'On this computer, R8=', R8
        ! *********** Read hdr file filename, memory divisor ************* 
        ! *********** (by default = 1), band csv filename.   ************* 
        numArgs =  command_argument_count()
        ALLOCATE (args(numArgs))
        DO ix = 1, numArgs
           CALL GET_COMMAND_ARGUMENT(ix, args(ix))
        END DO
        InputFileName = args(1)
        InputDatesFileName = args(2)
        OutputFilePrefix = args(3)
        OutputFileName = TRIM(OutputFilePrefix) // "_summary"
        OutputFileName2 = TRIM(OutputFilePrefix) // "_resids"
        OutputFileName3 = TRIM(OutputFilePrefix) // "_coeffs"
        InputHeaderFileName = TRIM(InputFileName) // ".hdr"  ! '//' concatenates the two strings
        
        !print *, "Memory divisor (enter a one if none needed): "
        !read (unit =*, FMT = *) MemoryDivisor
        MemoryDivisor = 1      
        !print *, "Input band info CSV file name: "
        !read (unit = *, FMT = *) InputDatesFileName

        ! *********************************************************    


        ! ************ from the hdr file, read image details *************

        open (unit = INPUT_HEADER_FILE, file = InputHeaderFileName, status = "old", action = "read")
        INQUIRE (unit = INPUT_HEADER_FILE, Opened=OS)

        OS = .FALSE.

        read (unit = INPUT_HEADER_FILE, fmt = "(A)") TempString
        
        read (unit = INPUT_HEADER_FILE, fmt = "(A)") TempString

        read (unit = INPUT_HEADER_FILE, fmt = "(A)") TempString
        
        read (unit = INPUT_HEADER_FILE, fmt = "(A)") TempString

        read (unit = INPUT_HEADER_FILE, fmt = "(A10,I4)") TempString, NumCols

        read (unit = INPUT_HEADER_FILE, fmt = "(A10,I4)") TempString, NumRows

        read (unit = INPUT_HEADER_FILE, fmt = "(A10,I3)") TempString, b

        do while (CurrHeaderLineNum .LT. UTM_LINE_NUMBER)
                read (unit = INPUT_HEADER_FILE, fmt = "(A)") TempString
                CurrHeaderLineNum = CurrHeaderLineNum + 1_2
        end do
    
        read (unit = INPUT_HEADER_FILE, fmt = "(A31, I6, A6, I7)") TempString, ULX, TempString3, ULY
        TempString = TRIM(TempString)
        TempString3 = TRIM(TempString3)

        ULY_Center = ULY
        ULX_Center = ULX

        CurrNorthing = ULY_Center
        CurrEasting = ULX_Center
        NumExtraRows = MOD(NumRows, INT(MemoryDivisor, KIND= 8))

        CLOSE(INPUT_HEADER_FILE)

        ! *********************************************************    
                
        NumPixels = NumRows*NumCols
        NumRowsChunk = NumRows
        print *, 'NumRows = ', NumRows, ', NumCols =', NumCols
        ALLOCATE (input_mat_old(NumCols, NumRows, b))
        ALLOCATE (input_mat(b, INT(NumCols/2), NumRows))
        ALLOCATE (z(NumPixels, 4*b))
        !ALLOCATE (ewma_summary (b, NumCols, NumRows))
        !ALLOCATE (ewma_residuals (b, NumCols, NumRows))
        !ALLOCATE (ewma_coeffs (5, NumCols, NumRows))
        !ALLOCATE (bfast_summary (b, NumCols, NumRows))
        !ALLOCATE (ltr_summary (b, NumCols, NumRows))
        !ALLOCATE (ltr_brkpt_summary (b, NumCols, NumRows))
  
       ! ******* Read the time stamps from band info csv file **********
       ALLOCATE (tyeardoy (b,BAND_INFO_FIELDS))

       open (unit = INPUT_DATES_FILE, file = InputDatesFileName,  status = "old", action = "read")
       INQUIRE (unit = INPUT_DATES_FILE, Opened=OS)

       OS = .FALSE.

       tyeardoy = 0

       read (unit = INPUT_DATES_FILE, fmt = "(A)") TempString
       do i = 1, b
             read (unit = INPUT_DATES_FILE, fmt = *) tyeardoy (i,1), tyeardoy (i,2), tyeardoy (i,3)
             !print *, tyeardoy (i,1), tyeardoy (i,2), tyeardoy (i,3)
       end do     
        
       !******** parameters for EWMACD **************
       ewmacd_numHarmonics = 2
       reverse = .FALSE.
       xbarlim1 = 1.5
       xbarlim2 = 20
       ewmacd_lowthresh = 1
       trainingStart = 2009
       trainingEnd = 2012
       L = 3
       mu = 0
       lam = 0.5
       persistence = 10
       !ewma_summary = -2222
       !ewma_residuals = -2222
       !ewma_coeffs = -2222

       !******** parameters for BFAST **************
       bfast_numBrks = 2
       bfast_numHarmonics = 1
       bfast_brkptSpacingProp = 0.15
       bfast_numColsProcess = 1
       bfast_pvalThresh = 0.1
       bfast_maxIter = 2
       ! the critical value table for use in bfast
       ALLOCATE (BBI(60, 4))
       BBI = reshape(  (/ 0.7552, 0.9809, 1.1211, 1.217, 1.2811, 1.3258, 1.3514, &
                 & 1.3628, 1.361, 1.3751, 0.7997, 1.0448, 1.203, 1.3112, 1.387,   &
                 & 1.4422, 1.4707, 1.4892, 1.4902, 1.5067, 0.825, 1.0802, 1.2491, &
                 & 1.3647, 1.4449, 1.5045, 1.5353, 1.5588, 1.563, 1.5785, 0.8414, &
                 & 1.1066, 1.2792, 1.3973, 1.4852, 1.5429, 1.5852, 1.6057, 1.6089, &
                 & 1.6275, 0.8541, 1.1247, 1.304, 1.425, 1.5154, 1.5738, 1.6182,   &
                 & 1.646, 1.6462, 1.6644, 0.8653, 1.1415, 1.3223, 1.4483, 1.5392, &
                 & 1.6025, 1.6462, 1.6697, 1.6802, 1.6939, 0.8017, 1.0483, 1.2059, &
                 & 1.3158, 1.392, 1.4448, 1.4789, 1.4956, 1.4976, 1.5115, 0.8431, &
                 & 1.1067, 1.2805, 1.4042, 1.4865, 1.5538, 1.59, 1.6105, 1.6156,   &
                 & 1.6319, 0.8668, 1.1419, 1.3259, 1.4516, 1.5421, 1.6089, 1.656, &
                 & 1.6751, 1.6828, 1.6981, 0.8828, 1.1663, 1.3533, 1.4506, 1.5791, &
                 & 1.6465, 1.6927, 1.7195, 1.7245, 1.7435, 0.8948, 1.1846, 1.3765, &
                 & 1.5069, 1.6077, 1.677, 1.7217, 1.754, 1.7574, 1.7777, 0.9048, &
                 & 1.1997, 1.3938, 1.5305, 1.6317, 1.7018, 1.7499, 1.7769, 1.7889, &
                 & 1.8052, 0.8444, 1.1119, 1.2845, 1.4053, 1.4917, 1.5548, 1.5946, &
                 & 1.6152, 1.621, 1.6341, 0.8838, 1.1654, 1.3509, 1.4881, 1.5779, &
                 & 1.653, 1.6953, 1.7206, 1.7297, 1.7455, 0.904, 1.1986, 1.3951, &
                 & 1.5326, 1.6322, 1.7008, 1.751, 1.7809, 1.7901, 1.8071, 0.9205, &
                 & 1.2217, 1.4212, 1.5593, 1.669, 1.742, 1.7941, 1.8212, 1.8269, &
                 & 1.8495, 0.9321, 1.2395, 1.444, 1.5855, 1.6921, 1.7687, 1.8176, &
                 & 1.8553, 1.8615, 1.8816, 0.9414, 1.253, 1.4596, 1.61, 1.7139, &
                 & 1.793, 1.8439, 1.8763, 1.8932, 1.9074, 0.8977, 1.1888, 1.3767, &
                 & 1.5131, 1.6118, 1.6863, 1.7339, 1.7572, 1.7676, 1.7808, 0.9351, &
                 & 1.2388, 1.4362, 1.5876, 1.693, 1.7724, 1.8223, 1.8559, 1.8668, &
                 & 1.8827, 0.9519, 1.27, 1.482, 1.6302, 1.747, 1.8143, 1.8756, 1.9105, &
                 & 1.919, 1.9395, 0.9681, 1.2918, 1.5013, 1.6536, 1.7741, 1.8573, &
                 & 1.914, 1.945, 1.9592, 1.9787, 0.9799, 1.3088, 1.5252, 1.6791, &
                 & 1.7967, 1.8837, 1.9377, 1.9788, 1.9897, 2.0085, 0.988, 1.622, &
                 & 1.5392, 1.7014, 1.8154, 1.9061, 1.9605, 1.9986, 2.0163, 2.0326 /), &
                 & (/ 60, 4/))
       !bfast_summary = 0
       !******** parameters for LandTrendR **************
       ltr_despikeTol = 0.5
       ltr_mu = 6
       ltr_nu = 3
       ltr_distwtfactor = 2.0
       ltr_recoveryThreshold = 1.0
       ltr_bestModelProp = 0.5
       ltr_pval = 0.05
       ltr_useFstat = 0  !0 means 'use p_of_f', '1' is for 'dont use p_of_f'
                         !So if use_fstat = 0, the code will use p_of_f.
       !ltr_summary = -2222
       !ltr_brkpt_summary = 0.0

       ! ***********************************************************************

       ! ************* Read the band values from the input file *********
       open (INPUT_FILE, file = InputFileName, form = "unformatted" , & 
                &   access = 'stream', status = "old", action = "read") 

       inquire(file = outputFileName, exist = OutputFileExists) 
  
       if(OutputFileExists) THEN
             open (OUTPUT_FILE, file = outputFileName, form = "unformatted" , &
                                access = 'stream', status = "old", action = "write")
       else if(.NOT.OutputFileExists) THEN
                open (OUTPUT_FILE, file = outputFileName, form = "unformatted" , &
                                access = 'stream', status = "new", action = "write")
       end if

       inquire(file = outputFileName2, exist = OutputFileExists) 

       if(OutputFileExists) THEN
             open (OUTPUT_FILE2, file = outputFileName2, form = "unformatted" , &
                                access = 'stream', status = "old", action = "write")
       else if(.NOT.OutputFileExists) THEN
                open (OUTPUT_FILE2, file = outputFileName2, form = "unformatted" , &
                                access = 'stream', status = "new", action = "write")
       end if

       inquire(file = outputFileName3, exist = OutputFileExists) 

       if(OutputFileExists) THEN
             open (OUTPUT_FILE3, file = outputFileName3, form = "unformatted" , &
                                access = 'stream', status = "old", action = "write")
       else if(.NOT.OutputFileExists) THEN
                open (OUTPUT_FILE3, file = outputFileName3, form = "unformatted" , &
                                access = 'stream', status = "new", action = "write")
       end if

       read (unit=INPUT_FILE) input_mat_old
       print *, 'going to read input'    
      !put input matrix in the order congenial to our need: for a fixed pixel,
      !we want all the timestamps in one place, i.e., time contiguous because we
      !are going to be processing one time series at a time.
      !Also, sort out valid pixels
       START = OMP_GET_WTIME()
       !$OMP PARALLEL DO  &
       !$OMP& SHARED(input_mat, input_mat_old, NumRows, NumCols) &
       !$OMP& PRIVATE(pixel_x, pixel_y, imgIndex)  &
       !$OMP& SCHEDULE(DYNAMIC)  &
       !$OMP& DEFAULT (NONE)
       DO imgIndex = 1,  NumPixels
            !pixel_x is the horizontal shift. So, the column number.
            !pixel_y is the vertical shift (downwards). So the row number.
            !So, for a 50x50 image:
            !pixel_x = mod(60, 50) = 10 th column
            !pixel_y = int(ceiling(60/50)) = 2  nd row
            pixel_x = mod(imgIndex, NumCols)
            IF (pixel_x == 0) THEN
                pixel_x = NumCols
            ENDIF
            pixel_y = INT(CEILING(REAL(imgIndex, KIND=R8)/NumCols))
            if (pixel_x <= FLOOR(NumCols/2.0)) then
               input_mat(:, pixel_x, pixel_y) = input_mat_old(pixel_x, pixel_y, :)
            endif
!            if ((pixel_x == 156) .and. (pixel_y==112)) then
!               print *, input_mat(:, pixel_x, pixel_y)
!            endif
       END DO
       !$OMP END PARALLEL DO
       FINISH = OMP_GET_WTIME()
       print *, FINISH - START

       deallocate (input_mat_old, z)

       num_obs = size(tyeardoy(:,1), 1)
       nb = ILAENV(1, 'DGELS', 'N', num_obs, 2*ewmacd_numHarmonics+1, -1, -1)
       nbDGELS = nb

       chunk_size = 100
       !CALL OMP_SET_DYNAMIC(.FALSE.)
       print *, 'going to parallel call ltr now '
       START = OMP_GET_WTIME()
       !$OMP PARALLEL DO  &
       !
       !The SHARED list contains shared variables across all threads.
       !$OMP& SHARED(input_mat, tyeardoy, &
       !$OMP&            reverse, num_obs, &
       !$OMP&            NumRows, NumCols, NumPixels, chunk_size, &
       !
       !$OMP&            ewmacd_numHarmonics, &
       !$OMP&            ewma_summary, ewma_residuals, ewma_coeffs, & 
       !$OMP&            xbarlim1, xbarlim2, ewmacd_lowthresh, lam, L, mu, &
       !$OMP&            trainingStart, trainingEnd, persistence, nb, &
       !
       !$OMP&            bfast_summary, &
       !$OMP&            bfast_numBrks, bfast_numHarmonics, &
       !$OMP&            bfast_brkptSpacingProp,  bfast_numColsProcess, &
       !$OMP&            bfast_pvalThresh, bfast_maxIter, nbDGELS, BBI,  &
       !
       !$OMP&            ltr_summary, ltr_brkpt_summary, &
       !$OMP&            ltr_despikeTol, ltr_mu, ltr_nu, ltr_distwtfactor, &
       !$OMP&            ltr_recoveryThreshold, ltr_bestModelProp, ltr_pval, &
       !$OMP&            ltr_useFstat)  &
       !The PRIVATE list contains uninitialized variables every thread should have a 
       !private copy of.
       !$OMP& PRIVATE(imgIndex, tmp)  &
       !$OMP& SCHEDULE(DYNAMIC, 1)  &
       !$OMP& DEFAULT (NONE)
         imgIndexloop: do imgIndex = 1, NumPixels, chunk_size
!             CALL ewmacd((/(tmp, tmp = imgIndex, MIN(NumPixels, imgIndex+chunk_size-1))/))

!             CALL bfast((/(tmp, tmp = imgIndex, MIN(NumPixels, imgIndex+chunk_size-1))/))

             CALL ltr((/ (tmp, tmp = imgIndex, MIN(NumPixels, imgIndex+chunk_size-1))/))
          end do imgIndexloop
       !$OMP END PARALLEL DO
       FINISH=OMP_GET_WTIME()
       print *, 'omp time=', FINISH-START
!       DEALLOCATE (work_arr%tmp_mat, work_arr%rtmp_vec1, work_arr%rtmp_vec2,   &
!                 &  work_arr%itmp_vec1, work_arr%itmp_vec2, work_arr%dev, STAT = derr)
   
       !leaving the output in time continuous order: will need THREE new variables otherwise.

       ! write the binary file of the output
!       write(OUTPUT_FILE)   ewma_summary
!       write(OUTPUT_FILE2) ewma_residuals
!       write(OUTPUT_FILE3) ewma_coeffs 
!       write(OUTPUT_FILE)   ltr_summary
      
       close (OUTPUT_FILE)
       close (OUTPUT_FILE2)
       close (OUTPUT_FILE3)
       close (INPUT_FILE)
    
       deallocate (input_mat, tyeardoy)
!       deallocate (ewma_summary, ewma_coeffs, ewma_residuals, STAT=derr)
!       deallocate (bfast_summary, STAT=derr)
!       deallocate (ltr_summary, STAT=derr)
!       deallocate (ltr_brkpt_summary, STAT=derr)
       deallocate (args)
       deallocate (BBI)

end program EWMACD_main 
