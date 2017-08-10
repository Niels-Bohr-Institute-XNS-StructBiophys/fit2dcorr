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

	fit2dcorr			Linux executable
	fit2dcorr.exe			Windows executable
	fit2dcorr.cpp			source code

	Makefile			Makefile
	compile_fit2dcorr.sh		shell script calling Makefile
	README				this README

	fit2d/fit2d_saxs_av_0.mac	Fit2D macro file
	fit2d/fit2d_saxs_av_1.mac	Fit2D macro file


Dependencies:

	libtiff
	(openmp)


Installation, using compiled executables:

	Install libtiff* package:

		Linux:
			sudo apt-get install libtiff5

	Clone git package (or download & extract zip-archive):

		git clone https://github.com/Niels-Bohr-Institute-XNS-StructBiophys/fit2dcorr.git
		cd fit2dcorr/

	Download Fit2D executables from http://ftp.esrf.eu/pub/expg/FIT2D/ to ../fit2d/:
	(note that both versions must be downloaded, the newer version 18 has a bug what affects the av == 0, for which the older version 12 must be used)

		Linux:
			fit2d_12_081_i686_linux2.4.20
			fit2d_beta_18_002_Debian7_intel64

			cd fit2d/
			wget "http://ftp.esrf.eu/pub/expg/FIT2D/fit2d_12_081_i686_linux2.4.20" "http://ftp.esrf.eu/pub/expg/FIT2D/fit2d_beta_18_002_Debian7_intel64"
			cd ..

		Windows:
			fit2d_12_077_i686_WXP.exe
			fit2d_beta_18_002_Windows7_intel32.exe

		Macintosh:
			fit2d_12_080_G3_MacOSX10.3.5
			fit2d_beta_18_002_MacOSX_7_5_intel64

			cd fit2d/
			wget "http://ftp.esrf.eu/pub/expg/FIT2D/fit2d_12_080_G3_MacOSX10.3.5" "http://ftp.esrf.eu/pub/expg/FIT2D/fit2d_beta_18_002_MacOSX_7_5_intel64"
			cd ..

	Run executable file in terminal 

		./fit2dcorr
		./fit2dcorr -h(elp)
		./fit2dcorr [<spec> <arg>]



Installation, compilation from source code:

	Install libtiff*(-dev) packages:

		Linux:
			sudo apt-get install libtiff5 libtiff5-dev

		Windows:


		MacOS:
			???


	Clone git package (see above)

	Download Fit2D executables (see above)

	Compiling:
		Linux:
			to compile it with gcc/g++ with OpenMP support and optimization level 3, use the provided Makefile and compile_fit2dcorr.sh script

			./compile_fit2dcorr.sh g++ 3 NONE

		Windows:
			to compile on win7 with libtiff and gcc/MSVC compiler

		Macintosh:
			???

	Run executable file in terminal (see above)









