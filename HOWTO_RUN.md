
fit2dcorr -- usage


specifier with "-" are optionally, those with "+" are mandatory:


	-av <av>					average mode:
							0 -> SAXS/GISAXS in Fit2D, without or Poisson-like errorbars
							     uses CAKE -> INTEGRATE module
							     allows more flexibility with regard to azimuthal and radial range
							1 -> SAXS/GISAXS in Fit2D, without or Poisson-like errorbars (default)
							     uses standard INTEGRATE module
							2 -> azimuthal averaging without Fit2D, without / with errorbars (not implemented yet)

	-err <err>					errorbar mode:
							0 -> no error bars (produces -1 numbers in 3rd column)
							1 -> Poisson-style (default)


	+mask <maskfile>				apply mask from Fit2D maskfile to all images in filelist
							e.g. -mask 300k_20Hz.msk

	+seq <first_file last_file min_pos max_pos>	sequence of input files with a fixed run number scheme between min_pos and max_pos
							e.g. +seq latest_0002367_craw.tiff latest_0003689_craw.tiff 8 14
							     +seq 002367.tif 003689.tif 1 6
							files in the sequence that don't exist will not be considered, a warning message will be posted
	OR
	+f <file_1 ... file_n>				list of several tif-files, not necessarily with a run number scheme
							e.g. 002367.tiff 002345.tiff 1-a6.tif test.tif

	+bc <bc[0] bc[1]				beam center x-y-coordinates applied for all images in filelist [pix]
							e.g. +bc 270.5 89.6
							can be read automatically for SAXSLAB tifs if xml entry is included
							e.g. +bc auto

	+sdd <sdd>					sample-detector-distance applied for all images in filelist [mm]
							e.g. +sdd 290
							can be read automatically for SAXSLAB tifs if xml entry is included
							e.g. +sdd auto

	+lambda <lambda>				wavelength applied for all images in filelist [Angstroem]
							e.g. +lambda 1.5418
							can be read automatically for SAXSLAB tifs if xml entry is included
							e.g. +lambda auto

	+pix_size <pixsz_x pixsz_y>			pixel sizes of detector applied for all images in filelist [microns]
							e.g. +pix_size 172 172
							can be read automatically for SAXSLAB tifs if xml entry / tif-header is included
							e.g. +pix_size auto

	-subtract <file_b vol_fract_b>			option to subtract a background file file_b scaled by volume fraction vol_fract_b from all files in filelist
							the subtraction is done with the 1D azimuthally averaged data, i.e. NOT on the 2D images (and then averaged)

	-abs_units <cf d T t (d_b T_b t_b) >		absolute units calibration applied for all images in filelist
							where: calibration factor cf [cps], sample thickness d [cm], Transmission T [0...1] and exposure time t [s]
							d, T, and t can be read automatically for SAXSLAB tifs if xml entries are included
							e.g. -abs_units 0.7 0.1 0.2456 3600
							     -abs_units 0.7 auto auto auto
							     -abs_units 0.7 auto 0.1 auto
							     -abs_units 0.7 0.1 0.2456 auto
							if -subtract option is used also d_b, T_b and t_b for the background file must be appended !
							e.g. -abs_units 0.7 auto auto auto auto auto auto
							     -abs_units 0.7 0.1 0.234 auto 0.1 0.221 auto

	-qscale <Qscale>				Q_nm-1 for "Q [1/nm]" (default), Q_A-1 for "Q [1/A]", s_nm-1 for "s [1/nm]", s_A-1 for "s [1/A]"

	-l <ranges>					list(s) of lines to skip in all averaged files, multiple instances are possible !
							CSV lists and / or range specifications possible
							e.g. -l 1:1:3
							     -l 3,4,5,6,235,236
							     -l 2:3:14
							     -l 3,4,5,6,235,236 -l 100:1:121
							note that other options like -nonnegative (considered first when writing output) might affect line numbers too!

	-rad_st <rad_st>				start inner radius for radial distance [pix] applied for all images in filelist
							default is 0.0 [pix], for modes av == 0 & 2 only

	-rad_end <rad_end>				end outer radius for radial distance [pix] applied for all images in filelist
							will be automatically calculated by default from the most distant corner of the image with respect to the bc
	OR
	-max_2theta <max_2theta>			maximum 2theta angle applied for all images in filelist
							will be automatically calculated by default from the most distant corner of the image with respect to the bc

	-rad_bins <rad_bins>				number of radial bins [integer] applied for all images in filelist
							will be automatically calculated by default from the most distant corner of the image to the bc and the pix_size
							for -av 0 and -rad_bins 1 the azimuthal intensity is exported as chi file for a (partial) ring / sector that,
							is defined by rad_st, rad_end (or max_2theta), azi_st, azi_end and azi_bins

	-azi_st <azi_st>				start azimuthal angle [deg] applied for all images in filelist
							default is 0.0 [deg], for modes av == 0 & 2 only

	-azi_end <azi_end>				end azimuthal angle [deg] applied for all images in filelist
							default is 360.0 [deg], for modes av == 0 & 2 only

	-azi_bins <azi_bins>				number of azimuthal bins [integer] applied for all images in filelist
							by default will be computed from azi_st and azi_end in 1 [deg] steps, for modes av == 0 & 2 only

	-pol_fac <pol_fac>				polarisation factor
							default is 0.99 for lab sources (for synchrotrons 0.95 might be better)

	-nonnegative					export only lines in *.chi-files where intensities are strictly positive (>0)
							in case of -subtract option this affects all lines where sample, background and difference are >0

	-openmp						OpenMP parallelization support, by default off

	-v						verbose mode, if used, temporary files like *.{chi,tif}_{avg,tot} will not be deleted to facilitate error analyis / debugging

	-mac <macrofile>				user-defined Fit2D macro-file, for modes av == 0 & 1 only


examples:

	(1) standard data reduction for a tiff-file with absolute intensity calibration using default mode (av == 1, err == 1) requesting a Q-scale in units of [1/nm]

	    fit2dcorr +mask mask.msk +f test.tiff +sdd 1023.3 +bc 270.4 245.9 +lambda 1.5418 +pix_size 172 172 -abs_units 23.4 0.1 0.65 3600 -qscale Q_nm-1

	(2) same as (1) but for for a list of SAXSLAB tiff-files where the different transmissions T and exposure times t can be read automatically

	    fit2dcorr +mask mask.msk +f file_1.tiff file_2.tiff test.tiff +sdd 1023.3 +bc 270.4 245.9 +lambda 1.5418 +pix_size 172 172 -abs_units 23.4 0.1 auto auto

	(3) data reduction of a sequence of SAXSLAB tiff-files with absolute intensity calibration and background subtraction,
	    images are SAXSLAB data format, containing the data d T t lambda pix_size as xml entries, that are automatically read,
	    integration shall be done in an angular sector (70-280 [deg], 60-140 [pix]) with 25 q-bins, what requires mode av == 0,
	    errorbars in subtracted file are computed by error propagation from the sample and background files,
	    all non-positive data points will be skipped as well as the first 5 points,
	    parallelization is used via OpenMP (depending on compilation and OS)

	    fit2dcorr -av 0 +mask mask.msk +seq im_0049301_caz.tiff im_0049632_caz.tiff 6 10 +sdd 1023.3 +bc 270.4 245.9 +lambda auto +pix_size auto -subtract buffer.tiff 0.99 -abs_units 23.4 auto auto auto auto auto auto -rad_st 60.0 -rad_end 140.0 -azi_st 70.0 -azi_end 280.0 -rad_bins 25 -l 1:1:5 -nonnegative -openmp

	(4) data reduction of a list of SAXSLAB tiff-files to extract the scattered intensity along a Debye-Scherrer ring sector,
	    normalization of intensity only by d(=0.1 [cm] fix) T t, that are automatically read
	    most parameters are read automatically, since they were written to the tiff-files during beamtime,
	    integration shall be done in a ring sector (0-360 [deg], 100-140 [pix]) with 120 angular-bins (3 [deg] steps), what requires mode av == 0,
	    rad_bins == 1 enforces intensity vs azimuthal angle,
	    no errobars are calculated (err == 0)

	    fit2dcorr -av 0 -err 0 +mask mask.msk +f file_1.tiff file_2.tiff file_3.tiff file_4.tiff +sdd auto +bc auto +lambda auto +pix_size auto -abs_units 1.0 0.1 auto auto -rad_st 120.0 -rad_end 140.0 -rad_bins 1


compilation:

	on Linux OS install package libtiff and call

	./compile_fit2dcorr.sh g++ 3 NONE

	to compile it with gcc/g++ with OpenMP support and optimization level 3, using the provided Makefile and compile_fit2dcorr.sh script

