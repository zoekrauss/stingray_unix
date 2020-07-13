# stingray_unix
Fortran source code for Stingray ray-tracing algorithm written in ascii I/O format by the Wilcock Lab Group

STINGRAY UNIX HOW-TO

The Stingray main code, which is Fortran, can be found in /stingray_src.F

Due to the complications of Fortran compiling and deprecation within Matlab, Zoe & William rewrote this Fortran source code to a simple unix command file, title stingray_src_unix.F

Instead of using the Matlab Mex file method, I have instead written in stingray.m to transfer all the variables used by the Fortran code to ascii files. Then, stingray.m calls the unix executable, and within that unix program (stingray_src_unix), these ascii files are manually read in and dimensionalized, and then the Fortran loop is called (located at the end of stingray_src_unix). 
NOTE: everytime the stingray_src_unix.F file is edited, it needs to be recompiled. So, if you have edited the Fortran code, once you are ready to run it, you need to delete the previous unix executable stingray_src_unix and recompile the stingray_src_unix.F file into a new executable. 


NOTE:
stingray_src_unix.F begins by specifying the dimensions of all the variables used. Since many of the variables are dimensioned USING input variables, this creates an issue- the program needs to use the input variables before we define them by reading them out of the ascii files we have created. 
To work around this, I wrote an “include” file that defines “parameters”, not variables (so just single numbers), that are way larger than the variables will ever be, in order to ensure that the dimensions of all of the variables are large enough. This does not cause any problems later by making variables too big, as it only stores temporary memory on the computer, but it allows us to use the ascii file method. 
The file can be found: Stingray/source/stingray_basedims.F
There is a “catch” function that will create an error if by chance these parameters in the include file are not large enough. In this case, the included file, Stingray/source/stingray_basedims.F, will need to be edited to make the parameters larger than the dimensions of the incoming variables (nx, ny, nz, nodes, etc.)

The logic of how this works:

1) A matlab script written to organize all input files and call the Stingray program, such as P_Crustal_Stingray_run.m, calls to the function stingray.m
     -->   This needs to be in the same directory as an srInput and srOutput folder
2) Stingray.m writes all necessary variables to ascii files
3) Stingray.m executes the unix file, stingray_src_unix
4) The unix file reads all the necessary variables out of the ascii files
5) The unix file executes the Fortran code
     -->   The loop is the second part of the unix file.
6) The unix file writes the output variables to ascii files, and returns to stingray.m
7) Stingray.m reads the output variables and uses them to create an srRays structure


To compile the stingray_src_unix.F code into a unix executable, Zoe uses the Intel iFort compiler, specifically with this command in the Matlab command window:
!/opt/intel/compilers_and_libraries_2020.0.166/mac/bin/intel64/ifort -o stingray_src_unix stingray_src_unix.F
The filepath part of this command will change depending on where in your computer the compiler is saved.

