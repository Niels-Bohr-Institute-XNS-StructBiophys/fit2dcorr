# fit2dcorr
SAXS 2D > 1D data reduction software (wrapper for Fit2D software):

	Should run on most Linux distros and Windows.
	Macintosh not tested yet, most probably only tiny code changes necessary.

	Intended as an easy terminal interface for data reduction with Fit2D.
	All images readable by Fit2D can be processed.
	Especially useful for batch processing and when using tiff-files from Pilatus detector images from SAXSLAB cameras.

	Integration is based on Fit2D's standard SAXS / GISAXS -> INTEGRATE (av == 1) module or the more flexible SAXS / GISAXS -> CAKE -> INTEGRATE (av == 0) module.
	For batch processing of large file sequences parallelization can be applied via OpenMP.

	Features:
		Use of masks created in Fit2D
		Absolute units calibration based on instrumental calibration factor, sample thickness, sample transmission and exposure time
		Errorbars dI based on Poisson approach
		Batch processing of image files
		Subtraction of a background / buffer file from all sample files
		Support for SAXSLAB's Pilatus tiffs and their xml entries to automatize processing through reading many parameters automatically, if provided in the tiff file
		Azimuthally averaging within (partial) ring sectors (I vs Q) or radially averaging along (partial) ring profiles (I vs phi) via mode av == 0 (CAKE -> INTEGRATE)
		Different scales for Q axis


Included files:

	Files				Comment

	fit2dcorr			Linux executable (compiled with gcc on Ubuntu 16.04 LTS 64bit using libtiff, 64 bit executable)
	fit2dcorr.exe			Windows executable (compiled with gcc using MinGW and libtiff on Win7 64bit, 32bit-executable)
	fit2dcorr.cpp			Source code

	Makefile			Makefile (for Linux)
	compile_fit2dcorr.sh		Shell script calling Makefile (for Linux)

	README.md			This README.md
	LICENSE				GNU General Public License v3.0
	HOWTO_RUN.md			Information on how-to run the program and its arguments, with examples

	fit2d/fit2d_saxs_av_0.mac	Fit2D macro file
	fit2d/fit2d_saxs_av_1.mac	Fit2D macro file
	fit2d/<Fit2D executables>	Fit2D executables to be downloaded by user (see below)

	test/winmingw_vs_linux/*	Test files to compare results between Linux and Windows version


Dependencies:

	libtiff
	(openmp)


Installation

	Linux (64bit)

		Installation, using compiled executables:

			Install libtiff* package (e.g. Debian-based):

				sudo apt-get install libtiff5

			Clone git package (or download & extract zip-archive):

				git clone https://github.com/Niels-Bohr-Institute-XNS-StructBiophys/fit2dcorr.git
				cd fit2dcorr/

			Download Fit2D executables from http://ftp.esrf.eu/pub/expg/FIT2D/ to ../fit2d/ and make them executable:
			(note that both versions must be downloaded, the newer version 18 has a bug what affects the av == 0, for which the older version 12 must be used)

				fit2d_12_081_i686_linux2.4.20
				fit2d_beta_18_002_Debian7_intel64

				cd fit2d/
				wget "http://ftp.esrf.eu/pub/expg/FIT2D/fit2d_12_081_i686_linux2.4.20" "http://ftp.esrf.eu/pub/expg/FIT2D/fit2d_beta_18_002_Debian7_intel64"
				chmod +x fit2d_12_081_i686_linux2.4.20 fit2d_beta_18_002_Debian7_intel64
				cd ..

			Run fit2dcorr executable file in terminal:
	
				./fit2dcorr
				./fit2dcorr -h(elp)
				./fit2dcorr [<spec> <arg>]


		Installation, compilation from source code:

			Install libtiff*(-dev) packages (e.g. Debian-based):

				sudo apt-get install libtiff5 libtiff5-dev

			Clone git package (see above)

			Download Fit2D executables (see above)

			Compiling:

				to compile it with gcc/g++ with OpenMP support and optimization level 3, use the provided Makefile and compile_fit2dcorr.sh script

				chmod +x compile_fit2dcorr.sh
				./compile_fit2dcorr.sh g++ 3 NONE

			Run executable file in terminal (see above)

			Test run in terminal:

				../../../fit2dcorr -av 1 -err 1 +mask ../masks_fit2d/mask.msk +bc 245.33 271.6 +f im_0049365_caz.tiff +sdd 266.4 +lambda auto +pix_size auto -abs_units 0.04764 0.1 0.67512 auto -nonnegative -qscale Q_A-1


	Windows 7 (MinGW + gcc + libtiff, 32 bit (works on 64bit, too) )

		Installation, compilation from source code:

			Download and install MinGW packages:
				https://sourceforge.net/projects/mingw/files/Installer/
				->mingw-get-setup.exe

				mingw32-base
				mingw32-gcc-g++
				mingw32-pthreads-w32
				(NB: https://stackoverflow.com/questions/39185248/eclipse-mingw-c-cannot-find-lpthread )

				Optionally apply some checks if gcc works in Windows cmd:
				gcc --version
				gcc -dM -E - <NUL:
				( -> __WIN32__ macro is known in MinGW )

			Add path for the gcc compiler (and other MinGW binaries) to the Windows $PATH ( e.g. C:\Program Files (x86)\MinGW\bin )

			Download and install libtiff ( e.g. under C:\Program Files (x86)\GnuWin32 ):
				http://gnuwin32.sourceforge.net/packages/tiff.htm
				-> tiff-3.8.2-1.exe

			Add path for libtiff to the Windows $PATH ( e.g. C:\Program Files (x86)\GnuWin32\bin )
			( $PATH is then ending on sth like ...;C:\Program Files (x86)\MinGW\bin;C:\Program Files (x86)\GnuWin32\bin )

			Clone git package (or download & extract zip-archive):

				git clone https://github.com/Niels-Bohr-Institute-XNS-StructBiophys/fit2dcorr.git
				cd fit2dcorr/

			Download Fit2D executables from http://ftp.esrf.eu/pub/expg/FIT2D/ to ../fit2d/ and make them executable:
			(note that both versions must be downloaded, the newer version 18 has a bug what affects the av == 0, for which the older version 12 must be used)

				fit2d_12_077_i686_WXP.exe
				fit2d_beta_18_002_Windows7_intel32.exe

			Compiling:

				e.g. ( -I(nclude) and -L(ibrary) paths for libtiff might differ ! ) :
				g++ -O3 -Wall -g3 -o fit2dcorr fit2dcorr.cpp -I"C:\Program Files (x86)\GnuWin32\include" -L"C:\Program Files (x86)\GnuWin32\lib" -ltiff -lm -fopenmp

			Run executable file in terminal (see above)

			Test run in Windows cmd:

				..\..\..\fit2dcorr.exe -av 1 -err 1 +mask ..\masks_fit2d\mask.msk +bc 245.33 271.6 +f im_0049365_caz.tiff +sdd 266.4 +lambda auto +pix_size auto -abs_units 0.04764 0.1 0.67512 auto -nonnegative -qscale Q_A-1


			Remarks:
			-> Path to the directory to fit2dcorr and Fit2D etc should not contain any spaces, since they cannot be read by Fit2D !
			-> currently same source file as under Linux
			-> works with direct.h since included in MinGW (getcwd)
			-> works with OPENMP_STDOUT macros (_iob)
			-> as under Linux OPENMP_STDOUT must be declared when using default(none) declaration !
			-> compilation and running fit2dcorr is done in Windows shell (cmd.exe)
			-> works with __WIN32__ macro
			-> only one instance of Fit2D can be run (not parallel) in macro-mode with popen.
			-> MinGW generated 32 bit executables cannot be run in Cygwin(64)


	Windows 7 (Cygwin + gcc + libtiff, 32/64bit)

		!!! NOT WORKING YET !!!

		Installation, compilation from source code:

			Remarks:

				-> allows compiling and calling / working in a nice unix-like shell 
				-> req. __CYGWIN__ macro
				-> does not work with OPENMP_STDOUT macros in parallel sessions
				-> stdout (w/o OMP declaration) could be used and OPENMP_STDOUT removed
				-> incompatible with MinGW mode and maybe Linux, too, since there must be declared ( at least when using default(none) declaration in pragma omp parallel )

				-> in future try to use stdout w/o declaration and omitting default(none) and skip MinGW tree for Windows (keep only __CYGWIN__ and __linux__)
				-> in future same source file as with linux but maybe not as with MinGW

				-> might be working with parallel sessions of popen (Fit2D)

				-> problem with process_chi_files() must be solved

			Download and install Cygwin 64 bit packages (gcc-g++ and libtiff) from https://www.cygwin.com/ :
				-> setup-x86_64.exe

				"gcc-package"
				libtiff-devel
				libtiff6

			Clone git package (or download & extract zip-archive):

				git clone https://github.com/Niels-Bohr-Institute-XNS-StructBiophys/fit2dcorr.git
				cd fit2dcorr/

			Download Fit2D executables from http://ftp.esrf.eu/pub/expg/FIT2D/ to ../fit2d/ and make them executable:
			(note that both versions must be downloaded, the newer version 18 has a bug what affects the av == 0, for which the older version 12 must be used)

				fit2d_12_077_i686_WXP.exe
				fit2d_beta_18_002_Windows7_intel32.exe

			Compiling:

				g++ -O3 -Wall -g3 -o fit2dcorr fit2dcorr.cpp -ltiff -lm -fopenmp


	Macintosh

		!!! NEVER TRIED !!!

		Installation, compilation from source code:

			Install libtiff*(-dev) packages:

				???

			Clone git package (or download & extract zip-archive):

				git clone https://github.com/Niels-Bohr-Institute-XNS-StructBiophys/fit2dcorr.git
				cd fit2dcorr/

			Download Fit2D executables from http://ftp.esrf.eu/pub/expg/FIT2D/ to ../fit2d/ and make them executable:
			(note that both versions must be downloaded, the newer version 18 has a bug what affects the av == 0, for which the older version 12 must be used)

				fit2d_12_080_G3_MacOSX10.3.5
				fit2d_beta_18_002_MacOSX_7_5_intel64

				cd fit2d/
				wget "http://ftp.esrf.eu/pub/expg/FIT2D/fit2d_12_080_G3_MacOSX10.3.5" "http://ftp.esrf.eu/pub/expg/FIT2D/fit2d_beta_18_002_MacOSX_7_5_intel64"
				chmod +x fit2d_12_080_G3_MacOSX10.3.5 fit2d_beta_18_002_MacOSX_7_5_intel64
				cd ..

			Compiling:

				???

			Run executable file in terminal (see above)


