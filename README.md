# lulcFortran
Fortran codes for lulc

Input: A binary file written in BSQ format.

Each image has M rows and N columns.
There are T timestamps (i.e., T images) in all.

BSQ format: The intuitive, simplest format where each line of the data is followed immediately by the next line in the same time point (or, spectral band).

            [111, 121, 131, 141, ...., 1NT, 
             211, 221, 231, 241, ...., 2N1,
             .
             M11, M21, M31, M41, ...., MN1;
             112, 122, 132, 142, ...., 1N2,
             212, 222, 232, 242, ...., 2N2,
             .
             M12, M22, M32, M42, ...., MN2;
             .
             .
             M1T, M2T, M3T, M4T, ...., MNT]

 
Output: Binary files written in BIL (time sequential) format. 

BIL format: All timepoints of first (11) pixel, then all timepoints of second (12) pixel,...., all timepoints of last (MN) pixel.

[111, 112, 113, 114, ...., 11T,
 
 121, 122, 123, 124, ...., 12T,
             .
             .
             .
 
 MN1, MN2, MN3, MN4, ...., MNT]
