# lulcFortran
Fortran codes for lulc

Inputs: (i) A binary file written in BSQ format (integer spectral values), (ii) A header file for the binary file, and (iii) a csv file with dates information.

Each image is an M x N matrix (M rows and N columns).
There are S images (i.e., S timestamps) in all.


BSQ format: The intuitive, simplest format where each line of the data is followed immediately by the next line in the same time point (or, spectral band).

            [111, 121, 131, 141, ...., 1N1,
             211, 221, 231, 241, ...., 2N1,
             .
             M11, M21, M31, M41, ...., MN1;
             112, 122, 132, 142, ...., 1N2,
             212, 222, 232, 242, ...., 2N2,
             .
             M12, M22, M32, M42, ...., MN2;
             .
             .
             11S, 12S, 13S, 14S, ...., 1NS,
             .
             .
             M1T, M2T, M3T, M4T, ...., MNS]
 
Output: Binary files written in BIL (time sequential) format. 

BIL format: All timepoints of first (11) pixel, then all timepoints of second (12) pixel,...., all timepoints of last (MN) pixel.
           
           [111, 112, 113, 114, ...., 11T,
            121, 122, 123, 124, ...., 12T,
             .
             .
             .
            MN1, MN2, MN3, MN4, ...., MNT]

Notes:
1) BRATIO.f90 and bsplines.f are part of standard libraries. 
I copied these codes from the splines codes I received from Dr. L.T. Watson.
2) Topmost code is poly\_united\_par.f90.
