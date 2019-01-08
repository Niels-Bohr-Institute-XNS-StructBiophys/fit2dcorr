test
/***************************************************************************
 *   Copyright (C) 2010-2017                                               *
 *                                                                         *
 *   Author: Martin Schmiele                                               *
 *   Email:  martin.schmiele@nbi.ku.dk                                     *
 *                                                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

using namespace std;

#ifdef WIN64
#define __WINDOWS__
#endif

#ifdef WIN32
#define __WINDOWS__
#endif

#ifdef __WINDOWS__
#include <Windows.h>
#include <direct.h>
#define OPENMP_STDOUT _iob
#define exp_str " "
#define linux_isdef false
#endif

#ifdef __linux__
#include <unistd.h>
#define OPENMP_STDOUT stdout
#define exp_str ""
#define linux_isdef true
#endif

#include <iostream>
#include <fstream>
#include <sstream>
// #include <string>
#include <algorithm>
#include <string.h>
#include <math.h>
#include <stdlib.h>
// #include <stdio.h>
#include <vector>
#include "tiffio.h"
#include <typeinfo>

// #ifndef M_PI
// #define M_PI 3.14159265358979323846
// #endif

#include <omp.h> /* OpenMP support */

class fit2dcorr
{
  
  public:

	vector<unsigned int> linelist ; /* -l <list of line numbers to skip> */

	string macrofile ; /* macro filename (full path) */
	bool mac_isdef ; /* -mac */

	string maskfile ; /* mask filename (full or relative path) */
	bool mask_isdef ; /* +mask */

	string fit2dfile ; /* Fit2D executable file (full path) */

	bool file_isdef ; /* +f */
	vector<string> file_list ; /* filename list to be processed */
	vector<string> fit2dfile_list ; /* list of latest Fit2D output files, av < 2 only */
 	vector<string> fit2dcall_list ; /* list of latest Fit2D calls, av < 2 only */

	bool seq_isdef ; /* +seq */
	vector<string> seqfile_list ; /* list containing the pairs of filenames defining the sequences */
	vector<unsigned int> seqpos_list ; /* list containing pairs of starting and terminating numbers in the seqfile_list defining the sequence */

	/* using a calibration factor for a specific instrument setting; sample thickness, transmission and exposure time */
	bool abs_units_isdef ;  /* -abs_units */ 
	bool cf_isauto, thickness_isauto, transmission_isauto, exptime_isauto ; /* flags if parameters have have been set */
	bool cf_is_i0, cf_b_is_i0 ;
	double cf_def, thickness_def, transmission_def, exptime_def ; /* values used for all files */
	vector<double> cf_list ; /* list of calibration factors [1/cps] */
	vector<double> thickness_list ; /* list of thicknesses [cm] */
	vector<double> transmission_list ; /* list of transmissions */
	vector<double> exptime_list ; /* list of exptimes [t] */
	vector<int> phi0_list ; /* list of flux [cps] */

	bool subtract_isdef ; /* -subtract */
	string file_b ; /* background file */
	double vol_fract_b ;
	bool cf_b_isauto, thickness_b_isauto, transmission_b_isauto, exptime_b_isauto ;
	double cf_b_def, thickness_b_def, transmission_b_def, exptime_b_def ;
	double cf_b, thickness_b, transmission_b, exptime_b ;
	int phi0_b ;

	double sdd ; /* sample-to-detector distance */
	bool sdd_isdef, sdd_isauto ; /* +sdd */

	double bc[2] ; /* x,y pixel coordinates for the beam center */
	bool bc_isdef, bc_isauto ; /* +bc */

	double lambda ; /* wavelength in Angstroem */
	bool lambda_isdef, lambda_isauto ; /* +lambda */

	double pix_size[2] ; /* pixel geometry */
	bool pix_size_isdef, pix_size_isauto ; /* +pix_size */

	string x_scale, x_scale_def ;
	double x_scale_fac, x_scale_fac_def ;

	string y_scale ;

	string executable_directory ;
	string working_directory ;
	string cmdbind ;
	string folderbind ;

	unsigned int av, av_def, max_av ; /* -av */
	unsigned int err, err_def, max_err[3] ; /* -err */

	bool is_intensity_conserved ;


	unsigned int array_size[2] ;
	unsigned int image_size[2] ;

	double pix_x_r_out, pix_y_r_out ;

	double max_2theta, min_2theta ; /* maximum and minimum theta angle for Q-range */
	bool max_2theta_isdef ; /* -max_2theta */
	double max_Q, min_Q ;

	double azi_st, azi_st_def ;
	double azi_end, azi_end_def ;
	unsigned int azi_bins ;
	bool azi_bins_isdef, azi_st_isdef, azi_end_isdef ; /* -azi_bins, -azi_st, -azi_end */

	bool ring_profile_isdef ;
	bool nonnegative_isdef ; /* -nonnegative */

	double rad_st, rad_st_def, rad_end, rad_max ;
	unsigned int rad_bins ;
	bool rad_st_isdef, rad_end_isdef, rad_bins_isdef ;
	double pol_fac, pol_fac_def ; /* -pol_fac polarisation factor */

	bool openmp_isdef ; /* -openmp */
	bool verbose_isdef ; /* -v */

	string presp ;

	static const unsigned int defsigns = 4096 ;

	/* for export of function call parameters */
	char **evarg ;
	unsigned int ecarg ;


	/* constructor */
	fit2dcorr() 
	{
		/* default qscale */
		this->x_scale_def = "Q_nm-1" ;
		this->x_scale_fac_def = 1.0 ;

		/* modes for averaging */
		this->av_def = 1 ;
		this->max_av = 2 ;

		/* errorbar handling for the modes */
		this->err_def = 1 ;
		this->max_err[0] = 1 ; /* av == 0 */
		this->max_err[1] = 1 ; /* av == 1 */
		this->max_err[2] = -1 ; /* av == 2 not implemented yet */

		/* defaults for the integration */
		this->azi_st_def = 0.0 ;
		this->azi_end_def = 360.0 ;

		this->rad_st_def = 0.0 ;

		this->pol_fac_def = 0.99 ;
 	}


	/* get sample (and background) thicknesses d, transmissions T and exptimes t from PILATUS/SAXSLAB xml tif-files or from user input */
	void get_cf_d_T_t()
	{
		string str ;
		bool interupt_isdef = false ;

		for ( unsigned int i=0; i<file_list.size(); ++i)
		{
			if ( cf_isauto )
			{
				str = read_xml_entry_from_tif( file_list[i], "saxsconf_Izero") ;
				if ( str.length() == 0 )
				{
					fprintf( OPENMP_STDOUT, "Warning: No I0 flux found from xml entries in SAXSLAB tif file %s.\n", file_list[i].c_str() ) ;
					cf_list.push_back( 0.0 ) ;

					interupt_isdef = true ;
				}
				else
				{
					/* scale flux with solid angle to obtain cf = sdd^2 / I0 / pix_area [1/cps]  */
					phi0_list.push_back( (int) strtol( str.c_str(), NULL, 10 ) ) ;
					cf_list.push_back( 1.0e6 * ( sdd * sdd ) / (double) phi0_list[i] / ( pix_size[0] * pix_size[1] ) ) ;
				}
			}
			else
			{
				if ( cf_is_i0 ) { cf_list.push_back( 1.0e6 * ( sdd * sdd ) / cf_def / ( pix_size[0] * pix_size[1] ) ) ; }
				else { cf_list.push_back( cf_def ) ; }
			}

			if ( thickness_isauto )
			{
				str = read_xml_entry_from_tif( file_list[i], "sample_thickness") ;
				if ( str.length() == 0 )
				{
					fprintf( OPENMP_STDOUT, "Warning: No thickness found from xml entries in SAXSLAB tif file %s.\n", file_list[i].c_str() ) ;
					thickness_list.push_back( 0.0 ) ;
					interupt_isdef = true ;
				}
				else
				{
					thickness_list.push_back( strtod( str.c_str(), NULL ) ) ;
				}
			}
			else { thickness_list.push_back( thickness_def ) ; }

			if ( transmission_isauto )
			{
				str = read_xml_entry_from_tif( file_list[i], "sample_transfact") ;
				if ( str.length() == 0 )
				{
					fprintf( OPENMP_STDOUT, "Warning: No transmission found from xml entries SAXSLAB tif in file %s.\n", file_list[i].c_str() ) ;
					transmission_list.push_back( 0.0 ) ;
					interupt_isdef = true ;
				}
				else
				{
					transmission_list.push_back( strtod( str.c_str(), NULL ) ) ;
				}
			}
			else { transmission_list.push_back( transmission_def ) ; }

			if ( exptime_isauto )
			{
				str = read_xml_entry_from_tif( file_list[i], "det_exposure_time") ;
				if ( str.length() != 0 ) { exptime_list.push_back( strtod( str.c_str(), NULL) ) ; }
				else // fallback to PILATUS header
				{
					fprintf( OPENMP_STDOUT, "Could not read exptime from SAXSLAB header in file %s, try to read it from PILATUS header instead.\n", file_list[i].c_str() ) ;
					exptime_list.push_back( get_exposure_time_from_PILATUS_tif( file_list[i]) ) ;
				}
			}
			else { exptime_list.push_back( exptime_def ) ; }
		}

		if ( subtract_isdef )
		{
			if ( cf_b_isauto )
			{
				str = read_xml_entry_from_tif( file_b, "saxsconf_Izero") ;
				if ( str.length() == 0 )
				{
					fprintf( OPENMP_STDOUT, "Warning: No I0 flux found from xml entries in SAXSLAB tif background file %s.\n", file_b.c_str() ) ;
					phi0_b = -1 ;
					cf_b = 0.0 ;
					interupt_isdef = true ;
				}
				else
				{
					/* scale flux with solid angle to obtain cf = sdd^2 / I0 / pix_area [1/cps]  */
					phi0_b = (int) strtol( str.c_str(), NULL, 10 ) ;
					cf_b = 1.0e6 * ( sdd * sdd ) / (double) phi0_b / ( pix_size[0] * pix_size[1] ) ;
				}
			}
			else
			{
				if ( cf_b_is_i0 ) { cf_b = 1.0e6 * ( sdd * sdd ) / cf_b_def / ( pix_size[0] * pix_size[1] )  ; }
				else { cf_b = cf_b_def ; }
			}

			if ( thickness_b_isauto )
			{
				str = read_xml_entry_from_tif( file_b, "sample_thickness") ;
				if ( str.length() == 0 )
				{
					fprintf( OPENMP_STDOUT, "Warning: No thickness found from xml entries in SAXSLAB tif background file %s.\n", file_b.c_str() ) ;
					thickness_b = 0.0 ;
					interupt_isdef = true ;
				}
				else
				{
					thickness_b = strtod( str.c_str(), NULL) ;
				}
			}
			else { thickness_b = thickness_b_def ; }

			if ( transmission_b_isauto )
			{
				str = read_xml_entry_from_tif( file_b, "sample_transfact") ;
				if ( str.length() == 0 )
				{
					fprintf( OPENMP_STDOUT, "Warning: No transmission found from xml entries in SAXSLAB tif background file %s.\n", file_b.c_str() ) ;
					transmission_b = 0.0 ;
					interupt_isdef = true ;
				}
				else
				{
					transmission_b = strtod( str.c_str(), NULL) ;
				}
			}
			else { transmission_b = transmission_b_def ; }

			if ( exptime_b_isauto )
			{
				str = read_xml_entry_from_tif( file_b, "det_exposure_time") ;
				if ( str.length() != 0 ) { exptime_b = strtod( str.c_str(), NULL) ; }
				else // fallback to PILATUS header
				{
					fprintf( OPENMP_STDOUT, "Read exptime for background file %s from PILATUS header.\n", file_b.c_str() ) ;
					exptime_b = get_exposure_time_from_PILATUS_tif( file_b) ;
				}
			}
			else { exptime_b = exptime_b_def ; }
		}
		
		string decision ;
		/* ask whether to proceed with data reduction if important warnings did arise */
		if ( interupt_isdef )
		{
			fprintf( stdout, "Some thicknesses and/or transmissions could not be read automatically and are set to zero. \n") ;
			fprintf( stdout, "Do you want to proceed [y] or exit [n] (recommended)? :\t");
			while (true)
			{
				cin >> decision ;
				if ( !decision.compare("y") ) { break ; }
				else if ( !decision.compare("n") ) { fprintf( stdout, "Terminating fit2dcorr.\n") ; exit(0) ; }
				else { fprintf( stdout, "Type 'y' or 'n':\t") ; }
			}
		}
	}


	/* from +seq option(s) the file_list is deduced */
	/* if files do not exist in the sequence they won't be included to file_list */
	void deduce_file_list_from_sequence()
	{
		FILE* file ;
		char data[defsigns] ;
		char sdummy[defsigns] ;
		unsigned int i_seq[2] ; /* first and last sequence integer */

		string strdummy, filename ;

		/* loop over all defined sequences */
		for ( unsigned int i=0; i < (seqpos_list.size() / 2) ; ++i)
		{
			/* fill i_seq from seqfile_list */
			for ( unsigned int j=0; j<2; ++j)
			{
				strdummy = seqfile_list[2*i+j].substr( seqpos_list[2*i] - 1, seqpos_list[2*i+1] - seqpos_list[2*i] + 1) ;
				stringstream ss(strdummy) ;
				if ( (ss >> i_seq[j]).fail() )
				{
					fprintf( stdout, "Error: Could not deduce start or terminating number from seqfile_list[%u] %s. Exit.\n", 2*i+j, seqfile_list[2*i+j].c_str()) ; 
					exit(1) ;
				}
			}

			if ( i_seq[0] > i_seq[1] ) { fprintf( stdout, "Error: Start number (%u) must be less equal than the end number (%u) in the sequence. Exit.\n", i_seq[0], i_seq[1]) ; exit(1) ; }

			strdummy = seqfile_list[2*i] ;
			sprintf( sdummy, "%%0%uu", seqpos_list[2*i+1] - seqpos_list[2*i] + 1) ;

			for ( unsigned int j=i_seq[0]; j<=i_seq[1]; ++j)
			{
				sprintf( data, sdummy, j) ;
				filename = strdummy.replace( seqpos_list[2*i] - 1, seqpos_list[2*i+1] - seqpos_list[2*i] + 1, string(data) ) ;

				/* check if file exists */
				if ( ( file = fopen( filename.c_str(), "r")) == NULL)
				{
					fprintf( stdout, "Warning: File %s in a defined file sequence is missing. Continue.\n", filename.c_str()) ;
				}
				fclose(file) ;

				file_list.push_back( filename ) ;
			}
		}

		/* if file_list is completely empty, terminate with an error */
		if ( file_list.empty() )
		{
			fprintf( stdout, "Error: None of the files in the defined file sequences does exist. Exit.\n") ;
			exit(1) ;
		}
	}


	void usage()
	{
		/* 
		   compilation:
		   GCC Linux
		   ./compile_fit2dcorr.sh g++ 3 NONE

		   GCC Windows
		   -> limited OpenMP support, popen cannot run in parallel in run_fit2d() !!!
		   using libtiff.dll.a
		   g++ -O3 -Wall -g3 -o fit2dcorr fit2dcorr.cpp -Iinclude -Llib -ltiff.dll -lm -fopenmp
		   equivalently, using libtiff.lib
		   g++ -O3 -Wall -g3 -o fit2dcorr fit2dcorr.cpp -Iinclude -Llib -ltiff -lm -fopenmp

		   MSVC in MSVS 2010 Express
		   -> no OpenMP support with Express Edition !!!
		   -> use _popen instead of popen
		   -> changes must be applied for GetModuleFileName
		   -> ...
		   -> better don't use MSVC ...
		*/

		/*  copy of usage from terminal */
		/*
		fit2dcorr -- usage


		specifiers with "-" are optionally, those with "+" are mandatory:


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


			+mask <maskfile>				apply mask from Fit2D maskfile to all images in file_list
									e.g. -mask 300k_20Hz.msk

			+seq <first_file last_file min_pos max_pos>	sequence of input files with a fixed run number scheme between min_pos and max_pos
									e.g. +seq latest_0002367_craw.tiff latest_0003689_craw.tiff 8 14
									     +seq 002367.tif 003689.tif 1 6
									files in the sequence that don't exist will not be considered, a warning message will be posted
									it is possible to use multiple instances of +seq option and to use +seq in combination with +f

			+f <file_1 ... file_n>				list of several tif-files, not necessarily with a run number scheme
									e.g. 002367.tiff 002345.tiff 1-a6.tif test.tif
									it is possible to use multiple instances of +f option and to use +f in combination with +seq

			+bc <bc[0] bc[1]				beam center y and z coordinates applied for all images in file_list [pix]
									e.g. +bc 270.5 89.6
									with auto argument, for SAXSLAB tifs it will be automatically read from each file if the xml entry is set: beamcenter_actual
									e.g. +bc auto

			+sdd <sdd>					sample-detector-distance applied for all images in file_list [mm]
									e.g. +sdd 290
									with auto argument, for SAXSLAB tifs it will be automatically read from each file if the xml entry is set: detector_dist
									e.g. +sdd auto

			+lambda <lambda>				wavelength applied for all images in file_list [Angstroem]
									e.g. +lambda 1.5418
									with auto argument, for SAXSLAB tifs it will be automatically read from each file if the xml entry is set: wavelength
									e.g. +lambda auto

			+pix_size <pixsz_y pixsz_z>			pixel sizes of detector applied for all images in file_list [microns]
									e.g. +pix_size 172 172
									with auto argument, for SAXSLAB tifs it will be automatically read from each file if the xml entry is set: pix_size
									with auto argument, for PILATUS tifs it will be automatically read from each file (via # Pixel_size)
									e.g. +pix_size auto

			-subtract <file_b vol_fract_b>			option to subtract a background file file_b scaled by volume fraction vol_fract_b from all files in file_list
									the subtraction is done with the 1D azimuthally averaged data, i.e. NOT on the 2D images (and then averaged)

			-abs_units <cf d T t (cf_b d_b T_b t_b)>	absolute units calibration applied for all images in file_list
									where:
									     cf (cf_b) - sample (background) calibration factor [1/cps]
									     d (d_b) - sample (background) thickness [cm]
									     T (T_b) - sample (background) transmission [0...1]
									     t (t_b) - sample (background) exposure time [s]
									with auto argument, for SAXSLAB tifs cf will be automatically calculated (with pix_size and sdd) from each file if the xml entry is set: saxsconf_Izero
									with auto argument, for SAXSLAB tifs d will be automatically read from each file if the xml entry is set: sample_thickness
									with auto argument, for SAXSLAB tifs d will be automatically read from each file if the xml entry is set: sample_transfact)
									with auto argument, for SAXSLAB tifs d will be automatically read from each file if the xml entry is set: det_exposure_time
									with auto argument, for PILATUS tifs t will be automatically read from each file (via # Exposure_time)
									e.g. -abs_units 0.7 0.1 0.2456 3600
									     -abs_units 0.7 auto auto auto
									     -abs_units 0.7 auto 0.1 auto
									     -abs_units 0.7 0.1 0.2456 auto
									     -abs_units auto auto auto auto
									if -subtract option is used also cf_b, d_b, T_b and t_b for the background file must be appended !
									e.g. -abs_units 0.72 auto auto auto 0.69 auto auto auto
									     -abs_units 0.72 0.1 0.234 auto 0.69 0.1 0.221 auto
									     -abs_units auto 0.15 auto auto auto 0.15 auto auto
									     -abs_units auto auto auto auto auto auto auto auto

			-cf_is_i0					CF value for sample files in -abs_units option is I0 and not CF, thus scaling must be applied from I0 to CF
			-cf_b_is_i0					CF value for backgr files in -abs_units option is I0 and not CF, thus scaling must be applied from I0 to CF

			-qscale <Qscale>				Q_nm-1 for "Q [1/nm]" (default), Q_A-1 for "Q [1/A]", s_nm-1 for "s [1/nm]", s_A-1 for "s [1/A]"

			-l <ranges>					list(s) of lines to skip in all averaged files, multiple instances are possible !
									comma separated lists and / or range specifications possible
									e.g. -l 1:1:3
									     -l 3,4,5,6,235,236
									     -l 2:3:14
									     -l 3,4,5,6,235,236 -l 100:1:121
									note that other options like -nonnegative (considered first when writing output) might affect line numbers too!

			-rad_st <rad_st>				start inner radius for radial distance [pix] applied for all images in file_list
									default is 0.0 [pix], for modes av == 0 & 2 only

			-rad_end <rad_end>				end outer radius for radial distance [pix] applied for all images in file_list
									will be automatically calculated by default from the most distant corner of the image with respect to the bc
			OR
			-max_2theta <max_2theta>			maximum 2theta angle applied for all images in file_list
									will be automatically calculated by default from the most distant corner of the image with respect to the bc

			-rad_bins <rad_bins>				number of radial bins [integer] applied for all images in file_list
									will be automatically calculated by default from the most distant corner of the image to the bc and the pix_size
									for -av 0 and -rad_bins 1 the azimuthal intensity is exported as chi file for a (partial) ring / sector that,
									is defined by rad_st, rad_end (or max_2theta), azi_st, azi_end and azi_bins

			-azi_st <azi_st>				start azimuthal angle [deg] applied for all images in file_list
									default is 0.0 [deg], for modes av == 0 & 2 only

			-azi_end <azi_end>				end azimuthal angle [deg] applied for all images in file_list
									default is 360.0 [deg], for modes av == 0 & 2 only

			-azi_bins <azi_bins>				number of azimuthal bins [integer] applied for all images in file_list
									by default will be computed from azi_st and azi_end in 1 [deg] steps, for modes av == 0 & 2 only

			-pol_fac <pol_fac>				polarisation factor
									default is 0.99 for lab sources (for synchrotrons 0.95 might be better)

			-nonnegative					export only lines in *.chi-files where intensities are strictly positive (>0)
									in case of -subtract option this affects all lines where sample, background and difference are >0

			-openmp						OpenMP parallelization support, by default off

			-v						verbose mode, if used, temporary files like *.{chi,tif}_{avg,tot} will not be deleted to facilitate error analysis / debugging

			-mac <macrofile>				user-defined Fit2D macro-file, for modes av == 0 & 1 only


		examples:

			(1) standard data reduction for a tiff-file with absolute intensity calibration using default mode (av == 1, err == 1) requesting a Q-scale in units of [1/nm]

			    fit2dcorr +mask mask.msk +f test.tiff +sdd 1023.3 +bc 270.4 245.9 +lambda 1.5418 +pix_size 172 172 -abs_units 23.4 0.1 0.65 3600 -qscale Q_nm-1

			(2) same as (1) but for for a list of SAXSLAB tiff-files where the different transmissions T and exposure times t can be read automatically

			    fit2dcorr +mask mask.msk +f file_1.tiff file_2.tiff test.tiff +sdd 1023.3 +bc 270.4 245.9 +lambda 1.5418 +pix_size 172 172 -abs_units 23.4 0.1 auto auto

			(3) data reduction of a sequence of SAXSLAB tiff-files with absolute intensity calibration and background subtraction,
			    images are SAXSLAB data format, containing the data I0 d T t lambda pix_size as xml entries, that are automatically read,
			    integration shall be done in an angular sector (70-280 [deg], 60-140 [pix]) with 25 q-bins, what requires mode av == 0,
			    errorbars in subtracted file are computed by error propagation from the sample and background files,
			    all non-positive data points will be skipped as well as the first 5 points,
			    parallelization is used via OpenMP (depending on compilation and OS)

			    fit2dcorr -av 0 +mask mask.msk +seq im_0049301_caz.tiff im_0049632_caz.tiff 6 10 +sdd 1023.3 +bc 270.4 245.9 +lambda auto +pix_size auto -subtract buffer.tiff 0.99 -abs_units auto auto auto auto auto auto auto -rad_st 60.0 -rad_end 140.0 -azi_st 70.0 -azi_end 280.0 -rad_bins 25 -l 1:1:5 -nonnegative -openmp

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


		*/

		fprintf( stdout, "\n") ;
		fprintf( stdout, "fit2dcorr -- usage\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "specifiers with \"-\" are optionally, those with \"+\" are mandatory:\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-av <av>\t\t\t\t\taverage mode:\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t0 -> SAXS/GISAXS in Fit2D, without or Poisson-like errorbars\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     uses CAKE -> INTEGRATE module\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     allows more flexibility with regard to azimuthal and radial range\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t1 -> SAXS/GISAXS in Fit2D, without or Poisson-like errorbars (default)\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     uses standard INTEGRATE module\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t2 -> azimuthal averaging without Fit2D, without / with errorbars (not implemented yet)\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-err <err>\t\t\t\t\terrorbar mode:\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t0 -> no error bars (produces -1 numbers in 3rd column)\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t1 -> Poisson-style (default)\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t+mask <maskfile>\t\t\t\tapply mask from Fit2D maskfile to all images in file_list\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\te.g. -mask 300k_20Hz.msk\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t+seq <first_file last_file min_pos max_pos>\tsequence of input files with a fixed run number scheme between min_pos and max_pos\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\te.g. +seq latest_0002367_craw.tiff latest_0003689_craw.tiff 8 14\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     +seq 002367.tif 003689.tif 1 6\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\tfiles in the sequence that don't exist will not be considered, a warning message will be posted\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\tit is possible to use multiple instances of +seq option and to use +seq in combination with +f\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t+f <file_1 ... file_n>\t\t\t\tlist of several tif-files, not necessarily with a run number scheme\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\te.g. 002367.tiff 002345.tiff 1-a6.tif test.tif\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\tit is possible to use multiple instances of +f option and to use +f in combination with +seq\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t+bc <bc[0] bc[1]\t\t\t\tbeam center y and z coordinates applied for all images in file_list [pix]\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\te.g. +bc 270.5 89.6\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\twith auto argument, for SAXSLAB tifs it will be automatically read from each file if the xml entry is set: beamcenter_actual\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\te.g. +bc auto\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t+sdd <sdd>\t\t\t\t\tsample-detector-distance applied for all images in file_list [mm]\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\te.g. +sdd 290\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\twith auto argument, for SAXSLAB tifs it will be automatically read from each file if the xml entry is set: detector_dist\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\te.g. +sdd auto\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t+lambda <lambda>\t\t\t\twavelength applied for all images in file_list [Angstroem]\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\te.g. +lambda 1.5418\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\twith auto argument, for SAXSLAB tifs it will be automatically read from each file if the xml entry is set: wavelength\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\te.g. +lambda auto\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t+pix_size <pixsz_y pixsz_z>\t\t\tpixel sizes of detector applied for all images in file_list [microns]\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\te.g. +pix_size 172 172\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\twith auto argument, for SAXSLAB tifs it will be automatically read from each file if the xml entry is set: pix_size\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\twith auto argument, for PILATUS tifs it will be automatically read from each file (via # Pixel_size)\n") ; // both are equivalent
		fprintf( stdout, "\t\t\t\t\t\t\te.g. +pix_size auto\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-subtract <file_b vol_fract_b>\t\t\toption to subtract a background file file_b scaled by volume fraction vol_fract_b from all files in file_list\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\tthe subtraction is done with the 1D azimuthally averaged data, i.e. NOT on the 2D images (and then averaged)\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-abs_units <cf d T t (cf_b d_b T_b t_b)>\tabsolute units calibration applied for all images in file_list\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\twhere:\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     cf (cf_b) - sample (background) calibration factor [1/cps]\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     d (d_b) - sample (background) thickness [cm]\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     T (T_b) - sample (background) transmission [0...1]\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     t (t_b) - sample (background) exposure time [s]\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\twith auto argument, for SAXSLAB tifs cf will be automatically calculated (with pix_size and sdd) from each file if the xml entry is set: saxsconf_Izero\n") ; // Izero = Imon * Ieff (Ieff for PIN-Diode typ >= 1, for Detector == 1)
		fprintf( stdout, "\t\t\t\t\t\t\twith auto argument, for SAXSLAB tifs d will be automatically read from each file if the xml entry is set: sample_thickness\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\twith auto argument, for SAXSLAB tifs d will be automatically read from each file if the xml entry is set: sample_transfact)\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\twith auto argument, for SAXSLAB tifs d will be automatically read from each file if the xml entry is set: det_exposure_time\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\twith auto argument, for PILATUS tifs t will be automatically read from each file (via # Exposure_time)\n") ; // both are equivalent
		fprintf( stdout, "\t\t\t\t\t\t\te.g. -abs_units 0.7 0.1 0.2456 3600\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     -abs_units 0.7 auto auto auto\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     -abs_units 0.7 auto 0.1 auto\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     -abs_units 0.7 0.1 0.2456 auto\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     -abs_units auto auto auto auto\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\tif -subtract option is used also cf_b, d_b, T_b and t_b for the background file must be appended !\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\te.g. -abs_units 0.72 auto auto auto 0.69 auto auto auto\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     -abs_units 0.72 0.1 0.234 auto 0.69 0.1 0.221 auto\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     -abs_units auto 0.15 auto auto auto 0.15 auto auto\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     -abs_units auto auto auto auto auto auto auto auto\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-cf_is_i0\t\t\t\t\tCF value for sample files in -abs_units option is I0 and not CF, thus scaling must be applied from I0 to CF\n") ;
		fprintf( stdout, "\t-cf_b_is_i0\t\t\t\t\tCF value for backgr files in -abs_units option is I0 and not CF, thus scaling must be applied from I0 to CF\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-qscale <Qscale>\t\t\t\tQ_nm-1 for \"Q [1/nm]\" (default), Q_A-1 for \"Q [1/A]\", s_nm-1 for \"s [1/nm]\", s_A-1 for \"s [1/A]\"\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-l <ranges>\t\t\t\t\tlist(s) of lines to skip in all averaged files, multiple instances are possible !\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\tcomma separated lists and / or range specifications possible\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\te.g. -l 1:1:3\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     -l 3,4,5,6,235,236\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     -l 2:3:14\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\t     -l 3,4,5,6,235,236 -l 100:1:121\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\tnote that other options like -nonnegative (considered first when writing output) might affect line numbers too!\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-rad_st <rad_st>\t\t\t\tstart inner radius for radial distance [pix] applied for all images in file_list\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\tdefault is 0.0 [pix], for modes av == 0 & 2 only\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-rad_end <rad_end>\t\t\t\tend outer radius for radial distance [pix] applied for all images in file_list\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\twill be automatically calculated by default from the most distant corner of the image with respect to the bc\n") ;
		fprintf( stdout, "\tOR\n") ;
		fprintf( stdout, "\t-max_2theta <max_2theta>\t\t\tmaximum 2theta angle applied for all images in file_list\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\twill be automatically calculated by default from the most distant corner of the image with respect to the bc\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-rad_bins <rad_bins>\t\t\t\tnumber of radial bins [integer] applied for all images in file_list\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\twill be automatically calculated by default from the most distant corner of the image to the bc and the pix_size\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\tfor -av 0 and -rad_bins 1 the azimuthal intensity is exported as chi file for a (partial) ring / sector that,\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\tis defined by rad_st, rad_end (or max_2theta), azi_st, azi_end and azi_bins\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-azi_st <azi_st>\t\t\t\tstart azimuthal angle [deg] applied for all images in file_list\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\tdefault is 0.0 [deg], for modes av == 0 & 2 only\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-azi_end <azi_end>\t\t\t\tend azimuthal angle [deg] applied for all images in file_list\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\tdefault is 360.0 [deg], for modes av == 0 & 2 only\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-azi_bins <azi_bins>\t\t\t\tnumber of azimuthal bins [integer] applied for all images in file_list\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\tby default will be computed from azi_st and azi_end in 1 [deg] steps, for modes av == 0 & 2 only\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-pol_fac <pol_fac>\t\t\t\tpolarisation factor\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\tdefault is 0.99 for lab sources (for synchrotrons 0.95 might be better)\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-nonnegative\t\t\t\t\texport only lines in *.chi-files where intensities are strictly positive (>0)\n") ;
		fprintf( stdout, "\t\t\t\t\t\t\tin case of -subtract option this affects all lines where sample, background and difference are >0\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-openmp\t\t\t\t\t\tOpenMP parallelization support, by default off\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-v\t\t\t\t\t\tverbose mode, if used, temporary files like *.{chi,tif}_{avg,tot} will not be deleted to facilitate error analysis / debugging\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t-mac <macrofile>\t\t\t\tuser-defined Fit2D macro-file, for modes av == 0 & 1 only\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "examples:\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t(1) standard data reduction for a tiff-file with absolute intensity calibration using default mode (av == 1, err == 1) requesting a Q-scale in units of [1/nm]\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t    fit2dcorr +mask mask.msk +f test.tiff +sdd 1023.3 +bc 270.4 245.9 +lambda 1.5418 +pix_size 172 172 -abs_units 23.4 0.1 0.65 3600 -qscale Q_nm-1\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t(2) same as (1) but for for a list of SAXSLAB tiff-files where the different transmissions T and exposure times t can be read automatically\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t    fit2dcorr +mask mask.msk +f file_1.tiff file_2.tiff test.tiff +sdd 1023.3 +bc 270.4 245.9 +lambda 1.5418 +pix_size 172 172 -abs_units 23.4 0.1 auto auto\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t(3) data reduction of a sequence of SAXSLAB tiff-files with absolute intensity calibration and background subtraction,\n") ;
		fprintf( stdout, "\t    images are SAXSLAB data format, containing the data I0 d T t lambda pix_size as xml entries, that are automatically read,\n") ;
		fprintf( stdout, "\t    integration shall be done in an angular sector (70-280 [deg], 60-140 [pix]) with 25 q-bins, what requires mode av == 0,\n") ;
		fprintf( stdout, "\t    errorbars in subtracted file are computed by error propagation from the sample and background files,\n") ;
		fprintf( stdout, "\t    all non-positive data points will be skipped as well as the first 5 points,\n") ;
		fprintf( stdout, "\t    parallelization is used via OpenMP (depending on compilation and OS)\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t    fit2dcorr -av 0 +mask mask.msk +seq im_0049301_caz.tiff im_0049632_caz.tiff 6 10 +sdd 1023.3 +bc 270.4 245.9 +lambda auto +pix_size auto -subtract buffer.tiff 0.99 -abs_units auto auto auto auto auto auto auto -rad_st 60.0 -rad_end 140.0 -azi_st 70.0 -azi_end 280.0 -rad_bins 25 -l 1:1:5 -nonnegative -openmp\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t(4) data reduction of a list of SAXSLAB tiff-files to extract the scattered intensity along a Debye-Scherrer ring sector,\n") ;
		fprintf( stdout, "\t    normalization of intensity only by d(=0.1 [cm] fix) T t, that are automatically read\n") ;
		fprintf( stdout, "\t    most parameters are read automatically, since they were written to the tiff-files during beamtime,\n") ;
		fprintf( stdout, "\t    integration shall be done in a ring sector (0-360 [deg], 100-140 [pix]) with 120 angular-bins (3 [deg] steps), what requires mode av == 0,\n") ;
		fprintf( stdout, "\t    rad_bins == 1 enforces intensity vs azimuthal angle,\n") ;
		fprintf( stdout, "\t    no errobars are calculated (err == 0)\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t    fit2dcorr -av 0 -err 0 +mask mask.msk +f file_1.tiff file_2.tiff file_3.tiff file_4.tiff +sdd auto +bc auto +lambda auto +pix_size auto -abs_units 1.0 0.1 auto auto -rad_st 120.0 -rad_end 140.0 -rad_bins 1\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "compilation:\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\ton Linux OS install package libtiff and call\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\t./compile_fit2dcorr.sh g++ 3 NONE\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\tto compile it with gcc/g++ with OpenMP support and optimization level 3, using the provided Makefile and compile_fit2dcorr.sh script\n") ;
		fprintf( stdout, "\n") ;

		exit(0) ;
	}


	/* returns a particular SAXSLAB-xml-entry as string from the SAXSLAB-modified PILATUS-tif-images */
	string read_xml_entry_from_tif(string filename, string entry)
	{
		/* typical structure for SAXSLAB 300k, the #-commented lines are from Dectris:

			<binary header line> # Pixel_size 172e-6 m x 172e-6 m
			# Silicon sensor, thickness 0.000320 m
			# Exposure_time 0.100000 s
			# Exposure_period 0.100000 s
			# Tau = 0 s
			# Count_cutoff 819995 counts
			# Threshold_setting: 4024 eV
			# Gain_setting: autog (vrf = 1.000)
			# N_excluded_pixels = 6
			# Excluded_pixels: badpixel_mask.tif
			# Flat_field: FF_p3-0166-500Hz_E8048_T4024_vrf_m0p100.tif
			# Trim_file: p3-0166-500Hz_E8048_T4024.bin
			# Image_path: /data/datatemp/
			# Ratecorr_lut_directory: Continuous
			# Retrigger_mode: 1
			<SAXSLAB xml entries>
			<binary data for tif-image>

		   a similar struture without the xml can be found for the 100k PILATUS images
		*/

		/* 
		  <param name="pixelsize">0.172 0.172</param> -> "0.172 0.172"
		  <param name="vp1">-0.500000</param> -> "-0.500000"
		  <param name="source_kV">50.00</param> -> "50.00"
		  <param name="detx">1900.000000</param> -> "1900.000000"
		  <param name="dety">0.572969</param> -> "0.572969"
		*/

		string line ;
		size_t start, end ;
		ifstream FILE ( filename.c_str() ) ;
		
		entry = "<param name=\"" + entry + "\">" ;
		
		if ( FILE.is_open() )
		{
			while( FILE.good() )
			{
				getline( FILE, line) ;
				start = line.find( entry ) ;
				if ( start != string::npos)
				{
					start += entry.length() ;
					end = line.find( "</param>", start ) ;
					FILE.close() ;
					return line.substr( start, end - start) ;
				}
			}
			FILE.close() ;
			cout << "Error: Could not find xml-entry " << entry << " in file " << filename << ". Maybe you can provide the parameter as an option on the command line? Exit.\n" ;
			exit(1) ;
		}
		else 
		{
			cout << "Error: Cannot open file " << filename << ". Exit.\n" ;
			exit(1) ;
		}
	}


	/* read and check optional(-) and mandatory(+) command line arguments */
	void eval_cmd (unsigned int carg, char **varg)
	{

		int j ;
		char *strptr1, *strptr2 ; 
		unsigned int first, last, increment, n;
		char data[defsigns] ;
		FILE *file ;
		double ddummy, ddummy2 ;
 
		/* hardcopy command line call parameters to evarg */
		ecarg = carg ;
		evarg = (char**) calloc( ecarg, sizeof(char*)) ;

		for ( unsigned int i=0; i<ecarg; ++i)
		{
			evarg[i] = (char*) calloc( strlen(varg[i]) + 1, sizeof(char)) ;
			strcpy( evarg[i], varg[i]) ;
		}

		/* print help (usage) information if no arguments appear or -h(elp) option is called and exit */
		if ( carg < 2 ) { usage() ; }

		/* assign default values */
		/* 1D chi-files */
		presp = "#" ;

		av = av_def ;
		err = err_def ;

		/* qscaling */
		x_scale = x_scale_def ;
		x_scale_fac = x_scale_fac_def ;

		pol_fac = pol_fac_def  ;

		azi_st = azi_st_def ; /* av == 0 & 2 */
		azi_end = azi_end_def ; /* av == 0 & 2 */

		rad_st = rad_st_def ; /* av == 0 & 2 */

		/* values to be defined by the user */
		/* mandatory */
		seq_isdef = false ;
		file_isdef = false ;
		mask_isdef = false ;

		bc_isdef = false ;
		bc_isauto = false ;
		sdd_isdef = false ;
		sdd_isauto = false ;
		lambda_isdef = false ;
		lambda_isauto = false ;
		pix_size_isdef = false ;
		pix_size_isauto = false ;

		/* optional */
		mac_isdef = false ;

		abs_units_isdef = false ;
		cf_isauto = false ;
		cf_is_i0 = false ;
		thickness_isauto = false ;
		transmission_isauto = false ;
		exptime_isauto = false ;

		subtract_isdef = false ;
		cf_b_isauto = false ;
		cf_b_is_i0 = false ;
		thickness_b_isauto = false ;
		transmission_b_isauto = false ;
		exptime_b_isauto = false ;

		/* prior to reading of input check for -subtract option, since it affects -abs_units input arguments */
		for ( unsigned int i=1; i<carg; ++i) { if ( !strcmp(varg[i], "-subtract") ) { subtract_isdef = true ; } }

		rad_st_isdef = false ;
		rad_end_isdef = false ;
		max_2theta_isdef = false;
		rad_bins_isdef = false ;

		azi_st_isdef = false ;
		azi_end_isdef = false ;
		azi_bins_isdef = false ;

		ring_profile_isdef = false ;
		nonnegative_isdef = false ;

		openmp_isdef = false ;
		verbose_isdef = false ;



		/* read option flags from cmd */
		fprintf( stdout, "Evaluate command line arguments\n") ;
		for ( unsigned int i=1; i<carg; ++i)
		{
			/* allow only flags of the kind +/-, numbers and text without +/- in front will cause the program to terminate */
			fprintf( stdout, "%s\n", varg[i]) ; 
			fflush(stdout) ;

			if ( ( (*varg[i] == '-') || (*varg[i] == '+') ) && ( is_numeric(varg[i],false) == false ) )
			{
				switch (varg[i][1])
				{
					/* -av -azi_st -azi_end -azi_bins -abs_units */
					case 'a':
						if ( !strcmp(varg[i], "-abs_units") )
						{
							/* 
								input must be either auto or numeric
								if auto the corresponding flag will remain false and program will try to read it automatically from SAXSLAB xml records in the tif file
								if non-auto it will be tried to read input as a number (if fails error and exit) and flag will set in case of success to true
							*/

							/* calibration factor cf */
							if ( ++i<carg )
							{
								if ( strcmp( varg[i], "auto") )
								{
									if ( is_numeric(varg[i]) ) { cf_def = strtod( varg[i], NULL) ; }
									else { fit2dcorr_error(2) ; }
								}
								else { cf_isauto = true ; }
							}
							else { fit2dcorr_error(1) ; }

							/* d, T and t for sample files */
							if ( ++i<carg )
							{
								if ( strcmp( varg[i], "auto") )
								{
									if ( is_numeric(varg[i]) ) { thickness_def = strtod( varg[i], NULL) ; }
									else { fit2dcorr_error(2) ; }
								}
								else { thickness_isauto = true ; }
							}
							else { fit2dcorr_error(1) ; }

							if ( ++i<carg )
							{
								if ( strcmp( varg[i], "auto") )
								{
									if ( is_numeric(varg[i]) ) { transmission_def = strtod( varg[i], NULL) ; }
									else { fit2dcorr_error(2) ; }
								}
								else { transmission_isauto = true ;} 
							}
							else { fit2dcorr_error(1) ; }

							if ( ++i<carg )
							{
								if ( strcmp( varg[i], "auto") )
								{
									if ( is_numeric(varg[i]) ) { exptime_def = strtod( varg[i], NULL) ; }
									else { fit2dcorr_error(2) ; }
								}
								else { exptime_isauto = true ; }
							}
							else { fit2dcorr_error(1) ; }

							/* cf_b, d_b, T_b and t_b for background file */
							if ( subtract_isdef )
							{
								if ( ++i<carg )
								{
									if ( strcmp( varg[i], "auto") )
									{
										if ( is_numeric(varg[i]) ) { cf_b_def = strtod( varg[i], NULL) ; }
										else { fit2dcorr_error(2) ; }
									}
									else { cf_b_isauto = true ; }
								}
								else { fit2dcorr_error(1) ; }

								if ( ++i<carg )
								{
									if ( strcmp( varg[i], "auto") )
									{
										if ( is_numeric(varg[i]) ) { thickness_b_def = strtod( varg[i], NULL) ; }
										else { fit2dcorr_error(2) ; }
									}
									else { thickness_b_isauto = true ; }
								}
								else { fit2dcorr_error(1) ; }

								if ( ++i<carg )
								{
									if ( strcmp( varg[i], "auto") )
									{
										if ( is_numeric(varg[i]) ) { transmission_b_def = strtod( varg[i], NULL) ; }
										else { fit2dcorr_error(2) ; }
									}
									else { transmission_b_isauto = true ;} 
								}
								else { fit2dcorr_error(1) ; }

								if ( ++i<carg )
								{
									if ( strcmp( varg[i], "auto") )
									{
										if ( is_numeric(varg[i]) ) { exptime_b_def = strtod( varg[i], NULL) ; }
										else { fit2dcorr_error(2) ; }
									}
									else { exptime_b_isauto = true ; }
								}
								else { fit2dcorr_error(1) ; }
							}

							abs_units_isdef = true ;
						}
						else if ( !strcmp(varg[i], "-av") )
						{
							if ( ++i<carg )
							{
								if ( is_numeric(varg[i]) ) { av = (unsigned int) strtol( varg[i], NULL, 10) ; }
								else { fit2dcorr_error(2) ; }
							}
							else { fit2dcorr_error(1) ; }

							if ( av > max_av ) { fit2dcorr_error(4) ; }
						}
						/*
						else if ( !strcmp(varg[i], "-array_size") )
						{
							for ( j=0; j<2; ++j)
							{
								if ( ++i<carg )
								{
									if ( is_numeric(varg[i]) ) { array_size[j] = strtol( varg[i], NULL, 10) ; }
									else { fit2dcorr_error(2) ; }
								}
								else { fit2dcorr_error(1) ; }
							}
						}
						*/
						else if ( !strcmp(varg[i], "-azi_st") )
						{
							if ( ++i<carg )
							{
								if ( is_numeric(varg[i]) ) { azi_st = strtod( varg[i], NULL) ; }
								else { fit2dcorr_error(2) ; }
								azi_st_isdef = true ;
							}
							else { fit2dcorr_error(1) ; }
							break;
						}
						else if ( !strcmp(varg[i], "-azi_end") )
						{
							if ( ++i<carg )
							{
								if ( is_numeric(varg[i]) ) { azi_end = strtod( varg[i], NULL) ; }
								else { fit2dcorr_error(2) ; }
								azi_end_isdef = true ;
							}
							else { fit2dcorr_error(1) ; }
							break;
						}
						else if ( !strcmp(varg[i], "-azi_bins") )
						{
							if ( ++i<carg )
							{
								if ( is_numeric(varg[i]) ) { azi_bins = strtol( varg[i], NULL, 10) ; }
								else { fit2dcorr_error(2) ; }
								azi_bins_isdef = true ;
							}
							else { fit2dcorr_error(1) ; }
							break;
						}
						else { fit2dcorr_error(3) ; }
						break ;
					/* +bc */
					case 'b':
						if ( !strcmp( varg[i], "+bc") )
						{
							/* 
								input must be either auto or two numerics
								if auto the corresponding flag will remain false and program will try to read it automatically from SAXSLAB xml entry in tif
								if non-auto it will be tried to read input as numbers (if fails error and exit) and flag will set in case of success to true
							*/
							for ( j=0; j<2; ++j)
							{
								if ( ++i<carg )
								{
									if ( !strcmp( varg[i], "auto") ) { bc_isauto = true ; break ; }

									if ( is_numeric(varg[i]) ) { bc[j] = strtod( varg[i], NULL) ; }
									else { fit2dcorr_error(2) ; }
									bc_isdef = true ;
								}
								else { fit2dcorr_error(1) ; }
							}
							break;
						}
						else { fit2dcorr_error(3) ; }
						break ;
					case 'c':
						if ( !strcmp( varg[i], "-cf_is_i0") )
						{
							cf_is_i0 = true ;
						}
						else if ( !strcmp( varg[i], "-cf_b_is_i0") )
						{
							cf_b_is_i0 = true ;
						}
						else { fit2dcorr_error(3) ; }
						break ;
					/* -err */
					case 'e':
						if ( !strcmp(varg[i], "-err") )
						{
							if ( ++i<carg )
							{
								if ( is_numeric(varg[i]) ) { err = (unsigned int) strtol( varg[i], NULL, 10) ; }
								else { fit2dcorr_error(2) ; }
							}
							else { fit2dcorr_error(1) ; }
						}
						else { fit2dcorr_error(3) ; }
						break ;
					/* +f */
					case 'f':
						if ( !strcmp(varg[i], "+f") )
						{
							while ( ++i<carg )
							{
								if ( *varg[i] != '-' && *varg[i] != '+' )
								{
									/* check if file exists annd add to file_list */
									if ( ( file = fopen( varg[i], "r")) == NULL)
									{
										fprintf( stdout, "Error: File %s in a +f option does not exist. Exit.\n", varg[i]) ;
										exit(1) ;
									}
									fclose(file) ;

									file_list.push_back( (string)varg[i] ) ;

									file_isdef = true ;
								}
								else
								{
									--i ;
									break ;
								}
							}
							if ( file_isdef == false )
							{
								fprintf( stdout, "Error: At least one file must be provided for +f. Exit.\n") ;
								exit(1) ;
							}
						}
						else { fit2dcorr_error(3) ; }
						break ;
					/* -h(elp) */
					case 'h':
						if ( !strcmp(varg[i], "-h") ) { usage() ; }
						else if ( !strcmp(varg[i], "-help") ) { usage() ; }
						else { fit2dcorr_error(3) ; }
						break ;
					/* -l +lambda */
					case 'l':
						/* it is possible to use -l option multiple times */
						if ( !strcmp( varg[i], "-l") )
						{
							++i;
							do
							{
								strptr1 = varg[i] ;
								if ( (strptr2 = strchr( varg[i], ':')) != NULL)
								{
									/* range specification of the form 2:17:53 or 1:2:13 */
									j = (int)(strptr2-strptr1) ;
									memmove( data, strptr1, j) ;
									data[j] = 0 ;
									first = (unsigned int)strtol( data, NULL, 10) ;

									strptr1 = strptr2 + 1 ;
									if ( (strptr2 = strchr( strptr1, ':')) != NULL)
									{
										j = (int)(strptr2-strptr1) ;
										memmove( data, strptr1, j) ;
										data[j] = 0 ;
										increment = (unsigned int)strtol( data, NULL, 10) ;
									}
									else
									{
										fprintf( stdout, "Error: Missing ':' in input for line numbers.\n") ;
										exit(1) ;
									}

									strptr1 = strptr2 + 1 ;
									strptr2 = strchr( strptr1, 0);
									{
										j = (int)(strptr2-strptr1) ;
										memmove( data, strptr1, j) ;
										data[j] = 0 ;
										last = (unsigned int)strtol( data, NULL, 10) ;
									}

									if ( (last-first) % increment != 0 )
									{
										fprintf( stdout, "Error: Wrong index specification in input for line numbers.\n") ;
										exit(1) ;
									}
									n = (last-first)/increment + 1 ;

									for ( unsigned int k=0; k<n; ++k) { linelist.push_back( first + k * increment) ; }
								}
								else
								{
									/* comma separated list e.g. 2,5,6,7,19 */
									strptr1 = varg[i] ;
									while ( (strptr2 = strchr( strptr1, ',')) != NULL)
									{
										j = (int)(strptr2-strptr1) ;
										memmove( data, strptr1, j) ;
										data[j] = 0 ;

										linelist.push_back( (unsigned int)strtol( data, NULL, 10)) ;

										strptr1 = strptr2 + 1 ;
									}
									strptr2 = strchr( varg[i], 0) ;
									j = (int)(strptr2-strptr1) ;
									memmove( data, strptr1, j) ;
									data[j] = 0 ;

									linelist.push_back( (unsigned int)strtol( data, NULL, 10)) ;
								}
								++i ;
								if ( i==carg ) { break ; }
							} while ( *varg[i] != '-' && *varg[i] != '+' ) ;
							--i;
							break;
						}
						else if ( !strcmp( varg[i], "+lambda") )
						{
							if ( ++i<carg )
							{
								if ( strcmp( varg[i], "auto") )
								{
									if ( is_numeric(varg[i]) ) { lambda = strtod( varg[i], NULL) ; }
									else { fit2dcorr_error(2) ; }
									lambda_isdef = true ;
								}
								else { lambda_isauto = true ; }
							}
							else { fit2dcorr_error(1) ; }
							break;
						}
						else { fit2dcorr_error(3) ; }
						break ;
					/* +mask -mac -max_2theta */
					case 'm':
						if ( !strcmp( varg[i], "-mac") )
						{
							if ( ++i<carg )
							{
								macrofile = (string)varg[i] ;
								mac_isdef = true ;
							}
							else { fit2dcorr_error(1) ; }
							break;
						}
						else if ( !strcmp( varg[i], "+mask") )
						{
							if ( ++i<carg )
							{
								/* check if maskfile exists */
								if ( ( file = fopen( varg[i], "r")) == NULL)
								{
									fprintf( stdout, "Error: Cannot read maskfile %s. File does not exist. Exit.\n", varg[i]) ;
									exit(1) ;
								}
								fclose(file) ;

								maskfile = (string)varg[i] ;

								mask_isdef = true ;
							}
							else { fit2dcorr_error(1) ; }
							break;
						}
						else if ( !strcmp( varg[i], "-max_2theta") )
						{
							if ( ++i<carg )
							{
								if ( is_numeric(varg[i]) ) { max_2theta = strtod( varg[i], NULL) ; }
								else { fit2dcorr_error(2) ; }
								max_2theta_isdef = true ;
							}
							else { fit2dcorr_error(1) ; }
							break;
						}
						else { fit2dcorr_error(3) ; }
						break ;
					/* -nonnegative */
					case 'n':
						if ( !strcmp(varg[i], "-nonnegative") )
						{
							nonnegative_isdef = true ;
						}
						else { fit2dcorr_error(3) ; }
						break ;
					/* -openmp */
					case 'o':
						if ( !strcmp(varg[i], "-openmp") )
						{
							/* -openmp, use OpenMP parallelization */
							openmp_isdef = true ;
						}
						else { fit2dcorr_error(3) ; }
						break ;
					/* -pol_fac +pix_size */
					case 'p':
						if ( !strcmp(varg[i], "+pix_size") )
						{
							/* 
								input must be either auto or two numerics
								if auto the corresponding flag will remain false and program will try to read it automatically from SAXSLAB xml entry in tif
								if non-auto it will be tried to read input as numbers (if fails error and exit) and flag will set in case of success to true
							*/
							for ( j=0; j<2; ++j)
							{
								if ( ++i<carg )
								{
									if ( !strcmp( varg[i], "auto") ) { pix_size_isauto = true ; break ; }

									if ( is_numeric(varg[i]) ) { pix_size[j] = strtod( varg[i], NULL) ; }
									else { fit2dcorr_error(2) ; }
									pix_size_isdef = true ;
								}
								else { fit2dcorr_error(1) ; }
							}
						}
						else if ( !strcmp(varg[i], "-pol_fac") )
						{
							if ( ++i<carg )
							{
								if ( is_numeric(varg[i]) ) { pol_fac = strtod( varg[i], NULL) ; }
								else { fit2dcorr_error(2) ; }
							}
							else { fit2dcorr_error(1) ; }
							break;
						}
						else { fit2dcorr_error(3) ; }
						break ;
					/* -qscale */
					case 'q':
						if ( !strcmp( varg[i], "-qscale") )
						{
							if ( ++i<carg )
							{
								x_scale = string(varg[i]) ;
							}
							else { fit2dcorr_error(1) ; }
							break;
						}
						else { fit2dcorr_error(3) ; }
						break ;
					/* -rad_st -rad_end -rad_bins */
					case 'r':
						if ( !strcmp( varg[i], "-rad_st") )
						{
							if ( ++i<carg )
							{
								if ( is_numeric(varg[i]) ) { rad_st = strtod( varg[i], NULL) ; }
								else { fit2dcorr_error(2) ; }
								rad_st_isdef = true ;
							}
							else { fit2dcorr_error(1) ; }
							break;
						}
						else if ( !strcmp( varg[i], "-rad_end") )
						{
							if ( ++i<carg )
							{
								if ( is_numeric(varg[i]) ) { rad_end = strtod( varg[i], NULL) ; }
								else { fit2dcorr_error(2) ; }
								rad_end_isdef = true ;
							}
							else { fit2dcorr_error(1) ; }
							break;
						}
						else if ( !strcmp( varg[i], "-rad_bins") )
						{
							if ( ++i<carg )
							{
								if ( is_numeric(varg[i]) ) { rad_bins = strtol( varg[i], NULL, 10) ; }
								else { fit2dcorr_error(2) ; }
								rad_bins_isdef = true ;
							}
							else { fit2dcorr_error(1) ; }
							break;
						}
						else { fit2dcorr_error(3) ; }
						break ;
					/* +sdd +seq -subtract */
					case 's':
						if ( !strcmp( varg[i], "+seq") )
						{
							for ( j=0; j<2; ++j)
							{
								if ( ++i<carg )
								{
									if ( *varg[i] != '-' && *varg[i] != '+' )
									{
										seqfile_list.push_back( (string)varg[i] ) ;
									}
									else
									{
										/* two files must be provided */
										fprintf( stdout, "Error: For +seq two input files must be provided. Exit.\n") ;
										exit(1) ;
									}
								}
								else { fit2dcorr_error(1) ; }
							}
							for ( j=0; j<2; ++j)
							{
								if ( ++i<carg )
								{
									if ( is_numeric(varg[i]) ) { seqpos_list.push_back( (unsigned int)strtol( varg[i], NULL, 10) ) ; }
									else { fit2dcorr_error(2) ; }
								}
								else { fit2dcorr_error(1) ; }
							}
							seq_isdef = true ;
							break;
						}
						else if ( !strcmp( varg[i], "+sdd") )
						{
							if ( ++i<carg )
							{
								if ( strcmp( varg[i], "auto") )
								{
									if ( is_numeric(varg[i]) ) { sdd = strtod( varg[i], NULL) ; }
									else { fit2dcorr_error(2) ; }
									sdd_isdef = true ;
								}
								else { sdd_isauto = true ; }
							}
							else { fit2dcorr_error(1) ; }
							break;
						}
						else if ( !strcmp(varg[i], "-subtract") )
						{
							if ( ++i<carg ) { file_b = (string)varg[i] ; }
							else { fit2dcorr_error(1) ; }

							if ( ++i<carg )
							{
								if ( is_numeric(varg[i]) ) { vol_fract_b = strtod( varg[i], NULL) ; }
								else { fit2dcorr_error(2) ; }
							}
							else { fit2dcorr_error(1) ; }
							/* optionally, set subtract_isdef (again) to true */
							subtract_isdef = true ;
						}
						else { fit2dcorr_error(3) ; }
						break ;
					/* -v */
					case 'v':
						if ( !strcmp( varg[i], "-v") )
						{
							verbose_isdef = true ;
						}
						else { fit2dcorr_error(3) ; }
						break ;
					default:
						/* unknown options */
						fprintf( stdout, "Error: Unknown option %s. Exit.\n", varg[i]) ;
						exit(1) ;
						break ;
				}
			}
			else
			{
				/* unknown input */
				fprintf( stdout, "Error: Unknown input %s. Exit.\n", varg[i]) ; 
				exit(1) ; 
				break ;
			}
		}
		fprintf( stdout, "done\n\n") ;



		/* av == 0 & 1 where Fit2D is involved */
		if ( av < 2 )
		{
			/*
				define cmdbind for command line options for Fit2D for all OS :
				"-" for Windows, MacOSX and Linux e.g. -<var=...>
	
				define folderbind for all OS

				get + check path to the Fit2D executable for all OS

				get + check path to the macrofile
			*/
			char buf[defsigns] ;
			size_t bufsize = sizeof(buf) ;
			#ifdef __WINDOWS__
				/* for av == 0 use Fit2D v12 since higher versions produce sometimes bugs when using azimuthal maps (CAKE -> INTEGRATE ) */
				if ( av == 0 ) { fit2dfile = "fit2d_12_077_i686_WXP.exe" ; }
				else if ( av == 1 ) { fit2dfile = "fit2d_beta_18_002_Windows7_intel32.exe" ; }

				/* 
				If the function succeeds, the return value is the length of the string that is copied to the buffer, 
				in characters, not including the terminating null character. If the buffer is too small to hold the
				module name, the string is truncated to nSize characters including the terminating null character, 
				the function returns nSize.
				*/
				/* MSVC only, use LPTSTR (wchar_t)
				wchar_t buf_wchar[defsigns] = {0} ;
				size_t len = GetModuleFileName( NULL, buf_wchar, bufsize) ;
				http://www.gamedev.net/topic/392816-converting-from-wchar_t-to-char/
				wcstombs( buf, buf_wchar, wcslen(buf_wchar) ) ;
				strcpy( buf, buf_wchar) ;
				*/

				size_t len = GetModuleFileName( NULL, buf, bufsize) ;

				if ( len != bufsize )
				{
					buf[len] = '\0' ;
				}
				else
				{
					fprintf( stdout, "buffer too small for executable_directory. Exit.\n") ;
					exit(1) ;
				}

				/* however Fit2D has problems with path including white spaces in Windows, independent of compiler used */

				/* furthermore the cmdbind for Fit2D is not "/" as written in the manual, it is as in Linux "-" */
				/* cmdbind = "/" ; */
				cmdbind = "-" ;

				folderbind = "\\" ;
			#elif __APPLE__
				/* not tested yet */

				/* for av == 0 use Fit2D v12 since higher versions produce sometimes bugs when using azimuthal maps (CAKE -> INTEGRATE ) */
				if ( av == 0 ) { fit2dfile = "fit2d_12_080_G3_MacOSX10.3.5" ; }
				else if ( av == 1 ) { fit2dfile = "fit2d_beta_18_002_MacOSX_7_5_intel64" ; }

				if ( _NSGetExecutablePath( buf, &bufsize) != 0 )
				{
					fprintf( stdout, "buffer too small for executable_directory, need size %u. Exit,\n", bufsize) ;
					exit(1) ;
				}
				cmdbind = "-" ;
				folderbind = "/" ;
			#elif __linux__
				/* for av == 0 use Fit2D v12 since higher versions produce sometimes bugs when using azimuthal maps (CAKE -> INTEGRATE ) */
				if ( av == 0 ) { fit2dfile = "fit2d_12_081_i686_linux2.4.20" ; }
				else if ( av == 1 ) { fit2dfile = "fit2d_beta_18_002_Debian7_intel64" ; }

				ssize_t len = readlink("/proc/self/exe", buf, bufsize-1) ;
				if ( len != -1 )
				{
					buf[len] = '\0' ;
				}
				else
				{
					fprintf( stdout, "buffer too small for executable_directory, need size %ld. Exit.\n", len) ;
					exit(1) ;
				}
				cmdbind = "-" ;
				folderbind = "/" ;
			#endif
			/* http://www.cplusplus.com/reference/string/string/find_last_of/ */
			executable_directory = string(buf) ;
			executable_directory = executable_directory.substr( 0, executable_directory.find_last_of("/\\") + 1 ) ;
			fprintf( stdout, "executable_directory is %s\n", executable_directory.c_str()) ;

			/* get current working directory */
			bufsize = sizeof(buf) ;
			if ( getcwd( buf, bufsize) != NULL)
			{
				working_directory = string(buf) + folderbind ;
				fprintf( stdout, "working_directory is %s\n", working_directory.c_str()) ;
			}
			else
			{
				fprintf( stdout, "buffer too small for working_directory. Exit.\n") ;
				exit(1) ;
			}

			/* full path to Fit2D executable */
			fit2dfile = executable_directory + "fit2d" + folderbind + fit2dfile ;

			/* check if Fit2D executable exists */
			if ( ( file = fopen( fit2dfile.c_str(), "r")) == NULL) 
			{ 
				fprintf( stdout, "Error: Cannot read Fit2D executable %s. File does not exist. Exit.\n", fit2dfile.c_str()) ;
				exit(1) ;
			}
			fclose(file) ;

			/* automatically define macrofile, if not set by user */
			if ( mac_isdef == false )
			{
				if ( av == 0 ) { macrofile = "fit2d_saxs_av_0.mac" ; }
				else if ( av == 1 ) { macrofile = "fit2d_saxs_av_1.mac" ;	}

				macrofile = executable_directory + "fit2d" + folderbind + macrofile ;
			}

			/* check if macrofile exists */
			if ( ( file = fopen( macrofile.c_str(), "r")) == NULL)
			{
				fprintf( stdout, "Error: Cannot read macrofile %s. File does not exist. Exit.\n", macrofile.c_str()) ;
				exit(1) ;
			}
			fclose(file) ;
		}
		else
		{
			/* av == 2 without Fit2D */
			fprintf( stdout, "Error: av == 2 mode has not been implemented yet. Exit.\n") ;
			exit(1) ;
		}



		/*
			apply some initial input checks and look for incompatibilities in the input (options), derive some arguments automatically
		*/




		/* check if there is only one rad_bin for av == 0 -> ring profile mode */
		if ( av == 0 && rad_bins == 1 ) { ring_profile_isdef = true ; }

		/* check if chosen err mode is permitted for chosen av mode */ 
		if ( err > max_err[av] ) { fit2dcorr_error(5) ; }

		/* 1st check if all essential parameters (+...) have been defined / flagged as auto */
		if ( ! ( bc_isdef || bc_isauto ) || ! ( sdd_isdef || sdd_isauto ) || ! ( lambda_isdef || lambda_isauto ) || ! ( pix_size_isdef || pix_size_isauto ) || !mask_isdef || ! ( seq_isdef || file_isdef ) ) { fit2dcorr_error(6) ; }


		/* in case of a +seq option derive filenames from the given sequence(s) and add filenames to file_list */
		/* file existence is checked, missing files are omitted (warning message is shown) */
		if ( seq_isdef ) { deduce_file_list_from_sequence() ; }


		/* check for duplicate files in file_list (incl file_b) since this leads to erros */
		if ( subtract_isdef ) { file_list.push_back( file_b ) ; }
		for ( unsigned int i = 0; i < file_list.size(); i++)
		{
			for ( unsigned int k = i+1; k < file_list.size(); k++)
			{
				if ( file_list.at(i).compare( file_list.at(k) ) == 0 )
				{
					fprintf( stdout, "Duplicate filenames in file_list (incl file_b) found. Exit.\n") ;
					exit(1) ;
				}
			}
		}
		if ( subtract_isdef ) { file_list.pop_back() ; }




		/*
		   now file_list is defined and other parameters might be derived from 1st file in it automatically

		   read some parameters that were flagged by auto from SAXSLAB xml entries
		*/


		/* 
		   image size of PILATUS-tif (e.g. 487x195 or 487x619) will be read via libtiff,
		   it is assumed that all files in file_list have the same image dimensions !!!
		*/
		info_tif( file_list[0], image_size[0], image_size[1]) ;


		/* 
		   derive pixel geometry from 1st image in file_list in case of auto via SAXSLAB xml entry
		   it is assumed that all images in the list have the same pixel sizes as the first one
		   transform to [microns] since they are given in [m]
		*/
		if ( pix_size_isauto )
		{
			string str = read_xml_entry_from_tif( file_list[0], "det_pixel_size") ;
			if ( str.length() != 0 )
			{
				stringstream ss(str) ;
				ss >> pix_size[0] >> pix_size[1] ;
				pix_size[0] *= 1.0e+6 ; pix_size[1] *= 1.0e+6 ;
			}
			else // fallback to PILATUS header, does NOT work yet
			{
				fprintf( stdout, "Could not read pix_size from SAXSLAB header in file %s, try to read it from PILATUS header instead.\n", file_list[0].c_str()) ;
				get_pixel_size_from_PILATUS_tif( pix_size[0], pix_size[1], file_list[0]) ;
			}
			pix_size_isdef = true ;
		}

		/* 
		   derive bc from 1st image in file_list in case of auto via SAXSLAB xml entry
		   it is assumed that all images in the list have the same bc as the first one
		   note that SAXSLABs bc coordinates are not the same as for Fit2D, apply transformation
		   y = z and z = image_size[1] - y
		*/
		if ( bc_isauto )
		{
			string str = read_xml_entry_from_tif( file_list[0], "beamcenter_actual") ;
			stringstream ss(str) ;
			ss >> bc[1] >> bc[0] ;
			// could be also image_size[1] - 1 instead ... same for bc[1] 1 pix offset unclear ...
			bc[1] = image_size[1] - bc[1] ; 

			bc_isdef = true ;
		}

		/* 
		   derive lambda from 1st image in file_list in case of auto via SAXSLAB xml entry
		   it is assumed that all images in the list have the same lambda as the first one
		 */
		if ( lambda_isauto )
		{
			lambda = strtod( read_xml_entry_from_tif( file_list[0], "wavelength").c_str(), NULL) ;
			lambda_isdef = true ;
		}

		/* 
		   derive sdd from 1st image in file_list in case of auto via SAXSLAB xml entry
		   it is assumed that all images in the list have the same sdd as the first one
		*/
		if ( sdd_isauto )
		{
			sdd = strtod( read_xml_entry_from_tif( file_list[0], "detector_dist").c_str(), NULL) ;
			sdd_isdef = true ;
		}


		/* 2nd check if all essential parameters (+...) are now defined */
		if ( !bc_isdef || !sdd_isdef || !lambda_isdef || !pix_size_isdef || !mask_isdef || ! ( seq_isdef || file_isdef ) ) { fit2dcorr_error(6) ; }



		/* print file_list and some common parameters that apply for all files */
		fprintf( stdout, "\n") ;
		fprintf( stdout, "summary of some input arguments:\n") ;
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\tfiles:\n") ;
		for ( unsigned int i=0; i<file_list.size(); ++i)
		{
			fprintf( stdout, "\t%s\n", file_list[i].c_str()) ;
		}
		fprintf( stdout, "\n") ;
		fprintf( stdout, "\tmask=%s\n", maskfile.c_str()) ;
		fprintf( stdout, "\tbc=(%-.1lf, %-.1lf) [pix]\n", bc[0], bc[1]) ;
		fprintf( stdout, "\tsdd=%-.1lf [mm]\n", sdd) ;
		fprintf( stdout, "\tlambda=%-.4lf [Angstroem]\n", lambda) ;
		fprintf( stdout, "\tpix_size=%-.1lf x %-.1lf [microns x microns]\n", pix_size[0], pix_size[1]) ;
		fprintf( stdout, "\n") ;





		/* x_scale labeling and define qscale scaling factor */
		if ( x_scale.compare("Q_nm-1") == 0 ) { x_scale = "Q [1/nm]" ; x_scale_fac = 1.0 ; }
		else if ( x_scale.compare("Q_A-1") == 0 ) { x_scale = " Q [1/A]" ; x_scale_fac = 0.1 ; }
		else if ( x_scale.compare("s_nm-1") == 0 ) { x_scale = "s [1/nm]" ; x_scale_fac = 1.0 / ( 2.0 * M_PI ) ; }
		else if ( x_scale.compare("s_A-1") == 0 ) { x_scale = " s [1/A]" ; x_scale_fac = 0.1 / ( 2.0 * M_PI ) ; }
		else
		{
			fprintf( stdout, "Error: valid arguments for -qscale option are: Q_nm-1 (default), Q_A-1, s_nm-1, s_A-1. Exit.\n") ;
			exit(1) ;
		}



		/* 
			get sample (and background) thicknesses d, transmissions T and exptimes t from PILATUS/SAXSLAB xml tif-files or from user input
			-> cf_list, thickness_list, transmission_list, exptime_list
			-> additionally also cf_b, thickness_b, transmission_b, exptime_b (for -subtract)
		 */
		if ( abs_units_isdef ) { get_cf_d_T_t() ; }

		/* y_scale labeling */
		if ( abs_units_isdef ) { y_scale = "I [1/cm]" ; }
		else { y_scale = "I [a.u.]" ; }



		/* 
		   either max_2theta[deg] or rad_end[pix] can be provided by user, not both.
		   actually, av == 0 strictly accepts only rad_st and rad_end but no max_2theta
		             av == 1 strictly accepts only max_2theta
		   however, these issues will be circumvented by the program, for both av modes rad_end/max_2theta will be accepted
		*/
		if ( max_2theta_isdef && rad_end_isdef )
		{ 
			fprintf( stdout, "Error: options -max_2theta and -rad_end cannot be used concurrently. Exit.\n") ;
			exit(1) ;
		}


		/* check for use of forbidden options for av == 1 */
		if ( av == 1 )
		{
			/*  -rad_st option to change rad_st_def = 0.0 is only permitted for av == 0 & 2 */
			if ( rad_st_isdef )
			{
				fprintf( stdout, "Error: option -rad_st is supported only for -av 0. Exit.\n") ;
				exit(1) ;
			}

			/* 
			-azi_st option to change azi_st_def = 0.0 is only permitted for av == 0,2
			-azi_end option to change azi_st_def = 360.0 is only permitted for av == 0,2
			-azi_bins option to change azi_bins is only permitted for av == 0,2
			*/
			if ( azi_st_isdef )
			{
				fprintf( stdout, "Error: option -azi_st is supported only for -av 0. Exit.\n") ;
				exit(1) ;
			}
			if ( azi_end_isdef )
			{
				fprintf( stdout, "Error: option -azi_end is supported only for -av 0. Exit.\n") ;
				exit(1) ;
			}
			if ( azi_bins_isdef )
			{
				fprintf( stdout, "Error: option -azi_bins is supported only for -av 0. Exit.\n") ;
				exit(1) ;
			}
		}


		/* radial stuff: pix_x_r_out, pix_y_r_out, rad_st, rad_end, rad_bins, max_2theta, min_Q, max_Q */

		/* determine automatically rad_max [pix] by checking check all four corners with respect to the beamcenter */
		/* start with southwest corner */
		ddummy = sqrt( pow( ( 1.0 - bc[0]) * pix_size[0], 2.0) + pow( ( 1.0 - bc[1]) * pix_size[1], 2.0) ) ;
		pix_x_r_out = 1.0 ;
		pix_y_r_out = 1.0 ;
		/* check southeast corner */
		ddummy2 = sqrt( pow( ((double)image_size[0] - bc[0]) * pix_size[0], 2.0) + pow( ( 1.0 - bc[1]) * pix_size[1], 2.0) ) ;
		if ( ddummy2 > ddummy )
		{
			pix_x_r_out = (double)image_size[0] ;
			pix_y_r_out = 1.0 ;
			ddummy = ddummy2 ;
		}
		/* check northwest corner */
		ddummy2 = sqrt( pow( ( 1.0 - bc[0]) * pix_size[0], 2.0) + pow( ( (double)image_size[1] - bc[1]) * pix_size[1], 2.0) ) ;
		if ( ddummy2 > ddummy )
		{
			pix_x_r_out = 1.0 ;
			pix_y_r_out = (double)image_size[1] ;
			ddummy = ddummy2 ;
		}
		/* check northeast corner */
		ddummy2 = sqrt( pow( ( (double)image_size[0] - bc[0]) * pix_size[0], 2.0) + pow( ( (double)image_size[1] - bc[1]) * pix_size[1], 2.0) ) ;
		if ( ddummy2 > ddummy )
		{
			pix_x_r_out = (double)image_size[0] ;
			pix_y_r_out = (double)image_size[1] ;
		}
		/* determine now rad_max[pix] with the most distant pix_x_r_out & pix_y_r_out */
		rad_max = sqrt( pow( ( pix_x_r_out - bc[0] ) * pix_size[0], 2.0) + pow( ( pix_y_r_out - bc[1] ) * pix_size[1], 2.0) ) ;
		rad_max /= (double)min( pix_size[0], pix_size[1] ) ;


		/* check now if 0 <= rad_st < rad_max, otherwise print error message */
		if ( rad_st_isdef ) /* applies only to av == 0 & 2 */
		{
			/* if user-defined, check for consistency */
			if ( ( rad_st >= rad_max ) || ( rad_st < 0.0 ) )
			{
				fprintf( stdout, "Error: rad_st must be in the range 0.0 <= rad_st < rad_max = %lf [pix]. Exit.\n", rad_max) ; 
				exit(1) ;
			}
		}


		/* set rad_end[pix] */

		/* nothing to do if user-defined */
		// if ( rad_end_isdef ) { rad_end = rad_end ; /* [pix] */ }

		/* compute from max_2theta, if user-defined */
		if ( max_2theta_isdef )
		{
			rad_end = ( sdd * 1000.0 / (double)min( pix_size[0], pix_size[1] ) ) * tan( M_PI / 180.0 * max_2theta ) ; /* [pix] */
		}
		/* if no limit was provided by user (default), set it automatically to the largest possible rad */
		if ( rad_end_isdef == false && max_2theta_isdef == false ) { rad_end = rad_max ; /* [pix] */ }

		/* check now if rad_end <= rad_max, otherwise print warning message */
		if ( rad_end > rad_max )
		{
			fprintf( stdout, "Warning: user-defined rad_end (or max_2theta) value produces rad_end > rad_max. Better try rad_end <= rad_max = %.2lf [pix].\n", rad_max) ;
		}

		/* 
		   if not already provided by user, derive max_2theta (2theta and Q-scans) and max_Q from rad_end, sdd, pix_size
		   max_2theta[deg] = 180 / pi * arctan( rad_end[pix] * pix_size[microns] / ( sdd[mm] * 1000 ) )
		   max_Q[1/nm] = ( 2.0 * M_PI ) * 10.0 * ( 2.0 / lambda[Angstroem] * sin( 0.5 * max_2theta[rad] ) )
		*/
		if ( max_2theta_isdef == false )
		{
			max_2theta = 180.0 / M_PI * atan( rad_end * (double)min( pix_size[0], pix_size[1] ) / ( sdd * 1000.0 ) ) ; /* [deg] */
		}
		max_Q = ( 2.0 * M_PI ) * 10.0 * ( 2.0 / lambda * sin( 0.5 * max_2theta * M_PI / 180.0 ) ) ; /* [1/nm] */

		/* compute min_2theta [deg] and min_Q [1/nm] */
		min_2theta = 180.0 / M_PI * atan( rad_st * (double)min( pix_size[0], pix_size[1] ) / ( sdd * 1000.0 ) ) ; /* [deg] */
		min_Q = ( 2.0 * M_PI ) * 10.0 * ( 2.0 / lambda * sin( 0.5 * min_2theta * M_PI / 180.0 ) ) ; /* [1/nm] */

		/* 
		   compute now number of rad_bins automatically if not user defined,
		   Fit2D computes it as the distance to the most remote corner of the image as seen from the bc
		*/
		if ( rad_bins_isdef == false )
		{
			/*
				Tests with Fit2D and automatically computed rad_max

				quadratic pixels pix_size = (172,172), 512x512
				(x,y)=(99,49.18) sqrt((487-99)^2+(195-49.18)^2)=414.497 -> rad_bins=414 (furthermost corner Right-Up)
				(x,y)=(99,49.18) sqrt((487-99)^2+(195-49.17)^2)=414.500 -> rad_bins=415 (furthermost corner Right-Up)

				(x,y)=(99,131.93) sqrt((487-99)^2+(131.93-1)^2)=409.496 -> rad_bins=409 (furthermost corner Right-Bottom)
				(x,y)=(99,131.95) sqrt((487-99)^2+(131.95-1)^2)=409.502 -> rad_bins=410 (furthermost corner Right-Bottom)

				rectangular pixels pix_size = (200,80), 1024x1024
				(x,y)=(100,50) sqrt(((487-100)*200)^2+((195-50)*80)^2)/80 = 978.3 -> rad_bins=978 (furthermost corner Right-Up)
			*/

			/* 
			   the approach works also for rectangular pixels, cause rad_end and rad_st are given by the user / derived from the program in pixels,
			   if automatic detection takes place for rad_st and rad_end -> rad_end - rad_st = rad_max 
			 */
			rad_bins = (unsigned int)( rad_end - rad_st + 0.5 ) ;
		}


		/* azimuthal stuff: azi_bins */

		/* setup automatically azi_bins if not user-defined with a stepsize of 1 [deg] */
		if ( azi_bins_isdef == false ) /* applies only to av == 0 & 2 */
		{
			/* 
			  difference values should be within [0,360] degrees.
			  the approach works for arbitrary azi_st and azi_end respecting the absolute position on the circle
			  89 -250 -> -339 -> 21 
			  -5 -10 -> -5 -> 355
			  89 87 -> -2 -> 358
			  87 89 -> 2 -> 2
			  -90 90 -> 180 -> 180
			  270 90 -> -180 -> 180
			  181 179 -> -2 -> 358
			  179 181 -> 2 -> 2
			  -540 90 -> 630 -> 270 
			  450 540 -> 90 -> 90
			  -540 540 -> 1080 -> 360 !! okay since rad_st < rad_end, see 0 360
			  0 360 -> 360 -> 360
			  540 -540 -> -1080 -> 0 !! okay since rad_st > rad_end, see 360 0
			  360 0 -> -360 -> 0 
			*/
			ddummy = azi_end - azi_st ;
			while ( ddummy < 0.0 ) { ddummy += 360.0 ; }
			while ( ddummy > 360.0 ) { ddummy -= 360.0 ; }
			/* ddummy is positive and within [0,360] */
			azi_bins = (unsigned int)( ddummy + 0.5 ) ;
			/* by contrast Fit2D would round values <=0.5 to lower binsizes and >0.5 to higher binsizes !!! */
		}



		/* 
		  define array_size for Fit2D  automatically

		  av = 0 :		array_size[0] = max( image_size[0], rad_bins )
		  (Fit2D SAXS_azi)	array_size[1] = max( image_size[1], azi_bins )

		  av = 1 :		array_size[0] = max( image_size[0], rad_bins )
		  (Fit2D SAXS)		array_size[1] = image_size[1]

		  av = 2 :		array_size[0] = max( image_size[0], rad_bins )
		  (no Fit2D)		array_size[1] = max( image_size[1], azi_bins )

		  check that array_size[i] for Fit2D is not larger than 100000 !
		*/

		array_size[0] = max( image_size[0], rad_bins) ;

		if ( av != 1 ) { array_size[1] = max( image_size[1], azi_bins) ; }
		else { array_size[1] = image_size[1] ; }

		if ( av < 2 )
		{
			for ( unsigned int i=0; i<2; ++i)
			{
				if ( array_size[i] > 100000 )
				{
					fprintf( stdout, "Error: array_size[%d] > 100000. Exit.\n", i) ;
					exit(1) ;
				}
			}
		}




		/*
			now all of the following class-wide variables are set either with default values, by input or are derived automatically

			sdd, lambda, bc[2], pix_size[2], maskfile, file_list
			macrofile, array_size[2]
			
			azi_st, azi_end, azi_bins
			pix_x_r_out, pix_y_r_out, rad_st, rad_end, rad_bins, max_2theta, min_Q, max_Q
			pol_fac
			x_scale, x_scale_fac, y_scale
			
			optionally also:
			cf_list, thickness_list, transmission_list, exptime_list (-abs_units)
			cf_b, file_b, thickness_b, transmission_b, exptime_b, vol_fract_b (-subtract and/or -abs_units)
		*/
	}


	/* error module for frequent error messages */
	void fit2dcorr_error(unsigned int ErrorCode)
	{
		const char *MP_ERR[] =
		{
			"No error", /* 00 */
			"Missing input after previous command line option.", /* 01 */
			"Expected a numeric input argument for previous command line option.", /* 02 */
			"Previous command line option is unknown.", /* 03 */
			"av argument for -av option must be less than max_av.", /* 04 */
			"err argument for -err option must be less than max_err.", /* 05 */
			"Missing input argument. Check if minimal input arguments +bc +sdd +lambda +pix_size +mask +seq (or +f) are defined." /* 06 */
		} ;

		/* before exit flush stdout */
		fflush(stdout) ;

		if ( ErrorCode == 0 )
		{
			return ;
		}

		fprintf( stderr, "Error: fit2dcorr_error(%d): %s\n", ErrorCode, MP_ERR[ErrorCode]) ;
		exit(1) ;
	}


	/* returns true if data contains a number representation (0,1,...,9,+,-,.,e,E)
	   e.g. -1.056, 2, 2.56E+07
	   but also for e.g. E-0.4+4.098E98....-+e
	*/
	bool is_numeric(char* data, bool with_exp=true)
	{
		bool bdummy = true ;
		unsigned int i = 0 ;
		while ( data[i] != 0 )
		{
			if ( isdigit(data[i]) ) { ++i ; continue ; }
			if ( data[i] == '.' ) { ++i ; continue ; }
			if ( ( data[i] == '-' ) || ( data[i] == '+' ) ) { ++i ; continue ; }
			if ( with_exp ) {  if ( ( data[i] == 'e' ) || ( data[i] == 'E' ) ) { ++i ; continue ; } }
			bdummy = false ;
			break ;
		}
		return bdummy ;
	}



	/*
	   perform absolute intensity calibration
	   returns dSigma/dOmega[1/cm] = cf[1/cps] * (I/t)[cps] / thickness[cm] / transmission
	*/
	void get_abs_intensity( double &I, double cf, double thickness, double transmission, double exptime)
	{
		I *= cf / exptime  / thickness / transmission ;
	}



	/* 
	   returns the exposure time in seconds given in the comment section of the PILATUS tif-files
	   checks for entry with the tag "# Exposure_time <time> <unit>"
	   in the following it is assumed that the unit is always seconds (s)
	   by this, if the user uses e.g. min, h as units in the camserver this is likely to cause problems ... 
	*/
	double get_exposure_time_from_PILATUS_tif(string filename)
	{
		FILE *file ;
		char sdummy[defsigns] ;
		char data[defsigns] ;
		char *strptr1, *strptr2 ;
		bool found = false ;
		int i;
		double exptime = -1.0 ;

		/* for Windows binary reading is essential, for Unix it doesn't matter, it works also with text modus */
		if ( ( file = fopen( filename.c_str(), "rb") ) == NULL )
		{ 
			fprintf( stdout, "Error: Cannot read file %s. File does not exist. Exit.\n", filename.c_str() ) ; 
			exit(1) ;
		}

		while ( fgets( sdummy, defsigns, file) != NULL )
		{
			strptr1 = sdummy ;

			if ( *strptr1 != '#' ) { continue ; }

			++strptr1 ;

			while ( *strptr1 == ' ' ) { ++strptr1 ; }

			if ( ( strptr2 = strchr( strptr1, ' ') ) != NULL )
			{
				/* check for correct Exposure_time pattern */

				i = (int)(strptr2-strptr1) ;
				memmove( data, strptr1, i) ;
				data[i] = 0 ;

				if ( !strcmp( data, "Exposure_time") )
				{
					found = true ;
					strptr1 = strptr2 ;
					while ( *strptr1 == ' ' ) { ++strptr1 ; }

					if ( ( strptr2 = strchr( strptr1, ' ') ) != NULL )
					{
						/* read exposure time, assume value is given in seconds */

						i = (int)(strptr2-strptr1) ;
						memmove( data, strptr1, i) ;
						data[i] = 0 ;

						exptime = strtod( data, NULL) ;
					}
					else
					{
						fprintf( stdout, "Error: Missing numeric value for exposure time in PILATUS tif-file %s. Exit.\n", filename.c_str() ) ;
						exit(1) ;
					}
					/* if successfully found, escape loop */
					break ;
				}
			}
		}

		fclose(file) ;

		if ( found == false )
		{
			fprintf( stdout, "Error: Cannot find exposure time in PILATUS tif-file %s. Exit.\n", filename.c_str() ) ; 
			exit(1) ;
		}
		return exptime ;
	}



	/* 
		does not work yet cause of binary file issues
	*/
	void get_pixel_size_from_PILATUS_tif( double& pix_size_y, double& pix_size_z, string filename)
	{
		FILE *file ;
		char sdummy[defsigns] ;
		char data[defsigns] ;
		char *strptr1, *strptr2 ;
		bool found = false ;
		int i;

		/* for Windows binary reading is essential, for Unix it doesn't matter, it works also with text modus */
		if ( ( file = fopen( filename.c_str(), "rb") ) == NULL )
		{ 
			fprintf( stdout, "Error: Cannot read file %s. File does not exist. Exit.\n", filename.c_str() ) ; 
			exit(1) ;
		}

		while ( fgets( sdummy, defsigns, file) != NULL )
		{
			strptr1 = sdummy ;

			if ( *strptr1 != '#' ) { continue ; }

			++strptr1 ;

			while ( *strptr1 == ' ' ) { ++strptr1 ; }

			if ( ( strptr2 = strchr( strptr1, ' ') ) != NULL )
			{
				/* check for correct Pixel_size pattern */
				// e.g. # Pixel_size 172e-6 m x 172e-6 m
				i = (int)(strptr2-strptr1) ;
				memmove( data, strptr1, i) ;
				data[i] = 0 ;

				if ( !strcmp( data, "Pixel_size") )
				{
					found = true ;
					strptr1 = strptr2 ;
					while ( *strptr1 == ' ' ) { ++strptr1 ; }

					if ( ( strptr2 = strchr( strptr1, ' ') ) != NULL )
					{
						/* read , assume value is given in m, transform to microns */
						i = (int)(strptr2-strptr1) ;
						memmove( data, strptr1, i) ;
						data[i] = 0 ;
						pix_size_y = 1000.0 * strtod( data, NULL) ;
					}
					else
					{
						fprintf( stdout, "Error: Missing numeric value for 1st pixel size in PILATUS tif-file %s. Exit.\n", filename.c_str() ) ;
						exit(1) ;
					}

					strptr1 = strptr2 ;
					if ( ( strptr2 = strchr( strptr1, 'x') + 1 ) == NULL )
					{
						fprintf( stdout, "Error: Missing numeric value for 2nd pixel size in PILATUS tif-file %s. Exit.\n", filename.c_str() ) ;
						exit(1) ;
					}
					strptr1 = strptr2 ;
					while ( *strptr1 == ' ' ) { ++strptr1 ; }

					if ( ( strptr2 = strchr( strptr1, ' ') ) != NULL )
					{
						/* read , assume value is given in m, transform to microns */
						i = (int)(strptr2-strptr1) ;
						memmove( data, strptr1, i) ;
						data[i] = 0 ;
						pix_size_z = 1000.0 * strtod( data, NULL) ;
					}
					else
					{
						fprintf( stdout, "Error: Missing numeric value for 2nd pixel size in PILATUS tif-file %s. Exit.\n", filename.c_str() ) ;
						exit(1) ;
					}

					/* if successfully found, escape loop */
					break ;
				}
			}
		}

		fclose(file) ;

		if ( found == false )
		{
			fprintf( stdout, "Error: Cannot find pixel size in PILATUS tif-file %s. Exit.\n", filename.c_str() ) ; 
			exit(1) ;
		}
	}



	/* function delallspc deletes all spaces in the input string str */
	/* http://www.cs.bu.edu/teaching/cpp/string/array-vs-ptr/ */
	void delallspc (char *str)
	{
		char dummy[strlen(str)+1] ;

		for (; *str; ++str)
		{
			if (*str == ' ') 
			{
				strcpy( dummy, str+1) ;
				strcpy( str, dummy) ;
				--str ;
			}
		}
	}



	/* azimuthal averaging */
	void run()
	{
		/* MSVC only, http://msdn.microsoft.com/de-de/library/0fatw238%28v=vs.80%29.aspx */
		/* _set_output_format(_TWO_DIGIT_EXPONENT) ; */


		// in case of background subtraction, append background file to file_list, since processing in Fit2D must be exactly the same
		if ( subtract_isdef ) { file_list.push_back( file_b ) ; }


		if ( av < 2 ) /* using Fit2D */
		{
			/* use always azimuthal averaging and not integration in 1st run */
			is_intensity_conserved = false ;

			// err == 0 -> one run, err == 1 two runs w/ and w/o intensity conservation
			for ( unsigned int i=0 ; i <= err ; ++i)
			{
				/* create Fit2D batch command and run on cmd */
				run_fit2d() ;

				if ( av == 0 )
				{
					/* read azimuthal map from files generated from Fit2D and derive azimuthal average or ring profile */
					process_azi_tif_files() ;
				}
				else /* av == 1 */
				{
					/* apply corrections for Fit2D generated chi-files, create file headers with parameters */
					process_chi_files() ;
				}

				/* use integration in 2nd run in case of err == 1 */
				is_intensity_conserved = true ;
			}
			/* reset to default */
			is_intensity_conserved = false ;

			// in case of background subtraction, remove background file now from file_list, to avoid its processing as a sample file later !
			if ( subtract_isdef ) { file_list.pop_back() ; fit2dfile_list.pop_back() ; fit2dcall_list.pop_back() ; }

			/*
				absolute intensity scaling (if abs_units_isdef)
				background subtraction (if subtract_isdef)
				error bars (if err == 1) on basis of Poisson approach sqrt( <I> / N_bin ) == sqrt( <I>^2 / I_tot ) == <I> / sqrt ( I_tot )
				line skipping ( nonnegative_isdef and linelist)
			*/
			correct_chi_files() ;
		}
		else /* av == 2 without Fit2D */
		{
			process() ;


			// in case of background subtraction, remove background file from file_list, to avoid its processing as a sample file 
//			if ( subtract_isdef ) { file_list.pop_back() ; fit2dfile_list.pop_back() ; fit2dcall_list.pop_back() ; }


		}
	}


	/* use Fit2D for av == 0 and av == 1 */
	void run_fit2d()
	{
		FILE *file ;
		char sdummy[defsigns], sdummy2[defsigns] ; /* Fit2D call C-string */
		char buff[1048576] ; /* 1 MByte buffer for Fit2D should be fine */
		string strconservint ;
		string strsuffix ;
		size_t pos ;

		unsigned int i ;

		/* 
		prepare template string, leaving i.a. in/output files free to be altered later

		av == 1 Example 
		fit2d -key -dim512x512 -fvar#SDD=292.0 -fvar#LAMBDA=1.5418 -fvar#BCX=278.0 -fvar#BCY=91.0 -fvar#PIXSZ_X=172.0 -fvar#PIXSZ_Y=172.0 -fvar#POL_FAC=0.99 -svar#MASK=mask20110414.msk -svar#CONSERVINT=NO -svar#FILE_IN=000870.tif -svar#FILE_OUT=000870.chi -ivar#RAD_BINS=230 -fvar#MAX_2THETA=7.85 -macfit2d_saxs.mac

		av == 0 Example 
		fit2d -key -dim512x512 -fvar#SDD=292.0 -fvar#LAMBDA=1.5418 -fvar#BCX=278.0 -fvar#BCY=91.0 -fvar#PIXSZ_X=172.0 -fvar#PIXSZ_Y=172.0 -fvar#POL_FAC=0.99 -svar#MASK=mask20110414.msk -svar#CONSERVINT=NO -svar#FILE_IN=000870.tif -svar#FILE_OUT=000870_pol.tif -fvar#AZI_START=-90.0 -fvar#AZI_END=90.0 -ivar#AZI_BINS=180 -fvar#PIX_X_ROUT=487.0 -fvar#PIX_Y_ROUT=195.0 -fvar#RAD_IN=0.0 -fvar#RAD_OUT=230.0 -ivar#RAD_BINS=230 -macfit2d_saxs_azi.mac
		*/


		//
		// For av == 0 and azi_bins == 1 (ring profiles, sectors) it would be also an option to export to chi directly instead of azi_tif and then using process_chi_files() as fot av == 1 !!!
		//



		/* file suffix handling */
		if ( av == 1 ) { strsuffix = ".chi" ; } 
		else /* av == 0 */ { strsuffix = "_azi.tif" ; }

		if ( is_intensity_conserved ) { strconservint = "YES" ; strsuffix += "_tot" ; }
		else { strconservint = "NO" ; strsuffix += "_avg" ; }


		/* store common part in sdummy */
 		sprintf( sdummy,"%s %skey %sdim%ux%u %sfvar#SDD=%-.1lf %sfvar#LAMBDA=%-.4lf %sfvar#BCX=%-.1lf %sfvar#BCY=%-.1lf %sfvar#PIXSZ_X=%-.1lf %sfvar#PIXSZ_Y=%-.1lf %sfvar#POL_FAC=%-.3lf %ssvar#MASK=%s %ssvar#CONSERVINT=%s %ssvar#FILE_IN=%%s %ssvar#FILE_OUT=%%s", fit2dfile.c_str(), cmdbind.c_str(), cmdbind.c_str(), array_size[0], array_size[1], cmdbind.c_str(), sdd, cmdbind.c_str(), lambda, cmdbind.c_str(), bc[0], cmdbind.c_str(), bc[1], cmdbind.c_str(), pix_size[0], cmdbind.c_str(), pix_size[1], cmdbind.c_str(), pol_fac, cmdbind.c_str(), maskfile.c_str(), cmdbind.c_str(), strconservint.c_str(), cmdbind.c_str(), cmdbind.c_str()) ;

		/* define individual part in sdummy2 for the different av modes */
		if ( av == 1 )
		{
			sprintf( sdummy2," %sivar#RAD_BINS=%u %sfvar#MAX_2THETA=%-.4lf %smac%s", cmdbind.c_str(), rad_bins, cmdbind.c_str(), max_2theta, cmdbind.c_str(), macrofile.c_str() ) ;
		}
		else /* av == 0 */
		{
			sprintf( sdummy2," %sfvar#AZI_START=%-.1lf %sfvar#AZI_END=%-.1lf %sivar#AZI_BINS=%u %sfvar#PIX_X_ROUT=%-.2lf %sfvar#PIX_Y_ROUT=%-.2lf %sfvar#RAD_IN=%-.2lf %sfvar#RAD_OUT=%-.2lf %sivar#RAD_BINS=%u %smac%s", cmdbind.c_str(), azi_st, cmdbind.c_str(), azi_end, cmdbind.c_str(), azi_bins, cmdbind.c_str(), pix_x_r_out, cmdbind.c_str(), pix_y_r_out, cmdbind.c_str(), rad_st, cmdbind.c_str(), rad_end, cmdbind.c_str(), rad_bins, cmdbind.c_str(), macrofile.c_str() ) ;
		}

		/* join them to sdummy */
		strcat( sdummy, sdummy2) ;

		fit2dfile_list.resize( file_list.size()) ;
		fit2dcall_list.resize( file_list.size()) ;

		#pragma omp parallel if ( openmp_isdef )
		{
			if ( omp_in_parallel() )
			{
				#pragma omp single
				fprintf( stdout ,"Parallelized Run: #-CPU=%d, #-Threads=%d\n", omp_get_num_procs(), omp_get_num_threads() ) ;
			}
			else
			{
				fprintf( stdout ,"Serialized Run\n") ;
			}
		}

		/* 
		   file_list, fit2dfile_list, fit2dcall_list are class members and always shared
		   however append operations on vectors should be not be applied to avoid race conditions !!!

		   tests with serialized runs and parallelized runs (-openmp flag) give exactly the same results

		   for Windows parallel runs with popen to run Fit2D do not work !!! However otherwise it works.
		*/
		#pragma omp parallel private( sdummy2, file, buff, pos) shared( strsuffix, sdummy, OPENMP_STDOUT ) default(none) if ( openmp_isdef && linux_isdef )
		{
			/* run Fit2D in batch mode */
			#pragma omp for schedule(static)
			for ( i=0; i<file_list.size(); ++i)
			{
				/* replace *.tif(f) with strsuffix for current filename and update accordingly fit2dfile_list */
				pos = file_list[i].find(".") ;
				fit2dfile_list[i] = file_list[i].substr( 0, pos) + strsuffix ;


				/* use pattern from sdummy with placeholders for strings with INPUT + OUTPUT files and write finally to sdummy2 */
				sprintf( sdummy2, sdummy, file_list[i].c_str(), fit2dfile_list[i].c_str()) ;

				/* save the program call in the string list fit2dcall_list */
				fit2dcall_list[i] = string(sdummy2) ;


				/* run Fit2D call */

				/* NB: for MSVC use _popen and _pclose instead of popen and pclose */
				fprintf( OPENMP_STDOUT, "Running %s\n", sdummy2) ;
				if ( ( file = popen( sdummy2, "r")) == NULL ) { fprintf( OPENMP_STDOUT, "Error: An error occured when executing Fit2D. Exit.\n") ; exit(1) ; }

				/* 
					WITHOUT SCANNING THE BUFFER, FIT2D WILL NOT RUN WHEN CALLED IN 
							fit2dcorr IN A TERMINAL

					WHEN DEBUGGING WITHIN KDevelop BUFFER READING IS NOT NECESSARY

					The origin of this strange behavior remains unclear, yet !
				*/
				while ( fgets( buff, sizeof(buff), file) != NULL )
				{
					fprintf( OPENMP_STDOUT, "%s\n", buff) ;
				}

				/* close the pipe */
				pclose(file) ;
			}
		}
	}



	/* edit (open-read-close, overwrite) the chi-files generated from Fit2D for av == 1 */
	void process_chi_files()
	{
		FILE *file ;
		char sdummy[defsigns] ;
		char sdummy2[defsigns] ;
		char data[defsigns] ;

		vector<string> content ;
		char *strptr1, *strptr2 ;
		unsigned int n ;

		double x, y ;

		unsigned int i, k ;
		unsigned int lineindex ;
		unsigned int offset ; /* number of first lines to skip from Fit2D generated chi-files */
		bool offset_isdef ;

		#pragma omp parallel if ( openmp_isdef )
		{
			if ( omp_in_parallel() )
			{
				#pragma omp single
				fprintf( stdout ,"Parallelized Run: #-CPU=%d, #-Threads=%d\n", omp_get_num_procs(), omp_get_num_threads() ) ;
			}
			else
			{
				fprintf( stdout ,"Serialized Run\n") ;
			}
		}

		/* 
		   fit2dfile_list, and other class members are always shared
		   however append operations on vectors should be not be applied to avoid race conditions !!!

		   tests with serialized runs and parallelized runs (-openmp flag) give exactly the same results
		*/
		#pragma omp parallel private( k, file, x, y, lineindex, offset, offset_isdef, content, sdummy, sdummy2, strptr1, strptr2, data, n ) shared( OPENMP_STDOUT ) default(none) if ( openmp_isdef )
		{
			/* batch mode */
			#pragma omp for schedule(static)
			for ( i=0; i<fit2dfile_list.size(); ++i)
			{
				/* read i-th file into memory and correct defect line */
				if ( ( file = fopen( fit2dfile_list[i].c_str(), "r")) == NULL )
				{ 
					fprintf( stdout, "Error: Cannot read file %s. File does not exist. Exit.\n", fit2dfile_list[i].c_str() ) ; 
					exit(1) ;
				}


				lineindex = 0 ;
				offset = 3 ;
				offset_isdef = false ;
				sprintf( sdummy2, "%d\n", rad_bins) ;
				while ( fgets( sdummy, defsigns, file) != NULL )
				{
					++lineindex ;
					/* skip the first offset lines, that must be detected automatically via the rad_bins entry */
					if ( lineindex <= offset ) { continue ; }
					if ( ( lineindex == offset + 1 ) && !offset_isdef )
					{
						delallspc( sdummy ) ;
						++offset ;
						if ( !strcmp( sdummy, sdummy2) ) { offset_isdef = true ; }
						continue ;
					}


					/* all other lines store in content */
					content.push_back( string(sdummy) ) ;
				}

				fclose(file) ;

				/* rewrite i-th file */
				if ( ( file = fopen( fit2dfile_list[i].c_str(), "w")) == NULL ) { fprintf( stdout, "Error: Cannot write to file %s. Exit.\n", fit2dfile_list[i].c_str()) ; exit(1) ; }

				/* create header with relevant parameters of the following type:

				LINUX example

				# working_directory=/home/mschmiele/projects/fit2dcorr/test/
				# fit2dfile=/home/mschmiele/projects/fit2dcorr/fit2d/fit2d_12_081_i686_linux2.4.20
				# macrofile=/home/mschmiele/projects/fit2dcorr/fit2d/fit2d_saxs.mac
				#
				# in=000871.tif
				# out=000871.chi
				# mask=mask20110414.msk
				# bcx=277.0 [pix]
				# bcy=91.5 [pix]
				# sdd=293.0 [mm]
				# lambda=1.5418 [Angstroem]
				# pix_size=172.0 x 172.0 [microns x microns]
				# max_2theta=7.8255 [deg]
				# rad_bins=234
				# array_size=487 x 195 [pix x pix]
				# pol_fac=0.990
				#
				# cf=0.700 [1/cps]
				# thickness=0.100 [cm]
				# transmission=0.3200
				# exposure time=10800.000 [s]
				#
				# Program called with: /home/mschmiele/projects/fit2dcorr/fit2dcorr -av 1 +mask mask20110414.msk +bc 277 91.5 +seq 000869.tif 000889.tif +sdd 293.0 +lambda 1.5418 +pix_size 172.0 172.0 -l 1:1:3 -qscale Q_nm-1 -abs_units 0.7 0.1 0.32 -openmp 
				#
				# Fit2D called with: /home/mschmiele/projects/fit2dcorr/fit2d/fit2d_12_081_i686_linux2.4.20 -key -dim487x195 -fvar#SDD=293.0 -fvar#LAMBDA=1.5418 -fvar#BCX=277.0 -fvar#BCY=91.5 -fvar#PIXSZ_X=172.0 -fvar#PIXSZ_Y=172.0 -fvar#POL_FAC=0.990 -svar#MASK=mask20110414.msk -svar#FILE_IN=000871.tif -svar#FILE_OUT=000871.chi -ivar#RAD_BINS=234 -fvar#MAX_2THETA=7.8255 -mac/home/mschmiele/projects/fit2dcorr/fit2d/fit2d_saxs.mac
				#
				#     Q [1/nm]       I [1/cm]      dI [1/cm]
				8.31870590E-02 1.75991132E+01 0.00000000E+00
				...


				WINDOWS example

				# working_directory=C:\Users\usr\Documents\cpp\fit2dcorr_gcc\test\
				# fit2dfile=C:\Users\usr\Documents\cpp\fit2dcorr_gcc\fit2d\fit2d_12_077_i686_WXP.exe
				# macrofile=C:\Users\usr\Documents\cpp\fit2dcorr_gcc\fit2d\fit2d_saxs.mac
				#
				# in=000871.tif
				# out=000871.chi
				# mask=mask20110414.msk
				# bcx=277.0 [pix]
				# bcy=91.5 [pix]
				# sdd=293.0 [mm]
				# lambda=1.5418 [Angstroem]
				# pix_size=172.0 x 172.0 [microns x microns]
				# max_2theta=7.8255 [deg]
				# rad_bins=234
				# array_size=487 x 195 [pix x pix]
				# pol_fac=0.990
				#
				# cf=0.700 [1/cps]
				# thickness=0.100 [cm]
				# transmission=0.3200
				# exposure time=10800.000 [s]
				#
				# Program called with: ..\fit2dcorr.exe -av 1 +mask mask20110414.msk +bc 277 91.5 +seq 000869.tif 000889.tif +sdd 293.0 +lambda 1.5418 +pix_size 172.0 172.0 -l 1:1:3 -qscale Q_nm-1 -abs_units 0.7 0.1 0.32 -openmp 
				#
				# Fit2D called with: C:\Users\usr\Documents\cpp\fit2dcorr_gcc\fit2d\fit2d_12_077_i686_WXP.exe -key -dim487x195 -fvar#SDD=293.0 -fvar#LAMBDA=1.5418 -fvar#BCX=277.0 -fvar#BCY=91.5 -fvar#PIXSZ_X=172.0 -fvar#PIXSZ_Y=172.0 -fvar#POL_FAC=0.990 -svar#MASK=mask20110414.msk -svar#FILE_IN=000871.tif -svar#FILE_OUT=000871.chi -ivar#RAD_BINS=234 -fvar#MAX_2THETA=7.8255 -macC:\Users\usr\Documents\cpp\fit2dcorr_gcc\fit2d\fit2d_saxs.mac
				#
				#      Q [1/nm]        I [1/cm]       dI [1/cm]
				8.31870660E-002 1.75991132E+001 0.00000000E+000
				...

				*/

				fprintf( file, "%s working_directory=%s\n", presp.c_str(), working_directory.c_str()) ;
				fprintf( file, "%s fit2dfile=%s\n", presp.c_str(), fit2dfile.c_str()) ;
				fprintf( file, "%s macrofile=%s\n", presp.c_str(), macrofile.c_str()) ;
				fprintf( file, "%s\n", presp.c_str()) ;
				fprintf( file, "%s mask=%s\n", presp.c_str(), maskfile.c_str()) ;
				fprintf( file, "%s bcx=%-.1lf [pix]\n", presp.c_str(), bc[0]) ;
				fprintf( file, "%s bcy=%-.1lf [pix]\n", presp.c_str(), bc[1]) ;
				fprintf( file, "%s sdd=%-.1lf [mm]\n", presp.c_str(), sdd) ;
				fprintf( file, "%s lambda=%-.4lf [Angstroem]\n", presp.c_str(), lambda) ;
				fprintf( file, "%s pix_size=%-.1lf x %-.1lf [microns x microns]\n", presp.c_str(), pix_size[0], pix_size[1]) ;
				fprintf( file, "%s max_2theta=%-.4lf [deg]\n", presp.c_str(), max_2theta) ;
				fprintf( file, "%s rad_bins=%d\n", presp.c_str(), rad_bins) ;
				fprintf( file, "%s array_size=%d x %d [pix x pix]\n", presp.c_str(), array_size[0], array_size[1]) ;
				fprintf( file, "%s pol_fac=%-.3lf\n", presp.c_str(), pol_fac) ;
				fprintf( file, "%s\n", presp.c_str()) ;
				fprintf( file, "%s in=%s\n", presp.c_str(), file_list[i].c_str()) ;
				fprintf( file, "%s out=%s\n", presp.c_str(), fit2dfile_list[i].c_str()) ;
				if ( abs_units_isdef )
				{
					fprintf( file, "%s cf=%-.3lf [1/cps]", presp.c_str(), cf_list[i]) ;
					if ( cf_isauto ) { fprintf( file, " ( phi0=%d [cps] )", phi0_list[i]) ; }
					fprintf( file, "\n") ;
					fprintf( file, "%s thickness=%-.3lf [cm]\n", presp.c_str(), thickness_list[i]) ;
					fprintf( file, "%s transmission=%-.5lf\n", presp.c_str(), transmission_list[i]) ;
					fprintf( file, "%s exptime=%-.3lf [s]\n", presp.c_str(), exptime_list[i]) ;
					fprintf( file, "%s\n", presp.c_str()) ;
				}
				if ( subtract_isdef )
				{
					fprintf( file, "%s background=%s\n", presp.c_str(), file_b.c_str()) ;
					fprintf( file, "%s vol_fract_b=%-.3lf\n", presp.c_str(), vol_fract_b) ;
					if ( abs_units_isdef )
					{
						fprintf( file, "%s cf_b=%-.3lf [1/cps]", presp.c_str(), cf_b) ;
						if ( cf_b_isauto ) { fprintf( file, " ( phi0_b=%d [cps] )", phi0_b) ; }
						fprintf( file, "\n") ;
						fprintf( file, "%s thickness_b=%-.3lf [cm]\n", presp.c_str(), thickness_b) ;
						fprintf( file, "%s transmission_b=%-.5lf\n", presp.c_str(), transmission_b) ;
						fprintf( file, "%s exptime_b=%-.3lf [s]\n", presp.c_str(), exptime_b) ;
						fprintf( file, "%s\n", presp.c_str()) ;
					}
				}
				fprintf ( file, "%s offset=%d\n", presp.c_str(), offset) ;
				fprintf ( file, "%s skipped lines=", presp.c_str()) ;
				if ( linelist.empty() ) { fprintf( file, "%s", "none") ; }
				for ( k=0; k<linelist.size(); ++k) 
				{ 
					fprintf ( file, "%d", linelist[k]) ;
					if ( k < linelist.size()-1 ) { fprintf ( file, "%s", ", ") ; }
				}
				fprintf( file, "\n") ;
				fprintf( file, "%s\n", presp.c_str()) ;
				fprintf ( file, "%s Program called with: ", presp.c_str()) ;
				/* Export function call parameters ecarg */
				for ( k=0; k<ecarg; ++k) { fprintf ( file, "%s ", evarg[k]) ; }
				fprintf ( file, "\n%s\n", presp.c_str()) ;

				fprintf ( file, "%s Fit2D called with: %s\n", presp.c_str(), fit2dcall_list[i].c_str()) ;
				fprintf ( file, "%s\n", presp.c_str()) ;

				fprintf( file, "%s     %s%s       %s%s      %sd%s\n", presp.c_str(), exp_str, x_scale.c_str(), exp_str, y_scale.c_str(), exp_str, y_scale.c_str()) ;

				/* append content */
				for ( k=0; k<content.size(); ++k)
				{
					/* fprintf( file, "%s", content[k].c_str() ) ; */
					strcpy( sdummy, content[k].c_str()) ;
					strptr1 = sdummy ;

					/* first column */
					while (*strptr1 == ' ') { ++strptr1 ; }

					if ( ( strptr2 = strchr( strptr1, ' ')) == NULL ) { fprintf( stdout, "Error: Expected two columns in file %s, only one is present. Exit.\n", fit2dfile_list[i].c_str()) ; exit(1) ; }

					n = (unsigned int)(strptr2-strptr1) ;
					memmove( data, strptr1, n) ;
					data[n] = 0 ;
					if ( is_numeric(data) ) { x = strtod( data, NULL) ; }
					else { fprintf( stdout, "Error: Found entry in the first column in file %s with non-numeric type. Exit.\n", fit2dfile_list[i].c_str()) ; exit(1) ; }

					/* q-scale default is Q [1/nm] */
					x *= x_scale_fac ;

					/* second column */
					strptr1 = strptr2 ;
					while (*strptr1 == ' ') { ++strptr1 ; }

					strptr2 = strchr( strptr1, '\n') ;

					n = (unsigned int)(strptr2-strptr1) ;
					memmove( data, strptr1, n) ;
					data[n] = 0 ;
					if ( is_numeric(data) ) { y = strtod( data, NULL) ; }
					else { fprintf( stdout, "Error: Found entry in the second column in file %s with non-numeric type. Exit.\n", fit2dfile_list[i].c_str()) ; exit(1) ; }

					fprintf( file, "%.8E %.8E %.8E\n", x, y, -1.0) ;
				}

				fclose(file) ;

				/* !!! IMPORTANT : clear content !!! */
				content.clear() ;
			}
		}
	}


	/* post-processing of azi-tif files from av == 0 */
	void process_azi_tif_files()
	{
		float *raw_image ; /* read 32 bit *_azi.tif into float(32 bit) column raw_image */
		double **image ; /* transform raw_image to double(64 bit) image */
		FILE *file ;

		vector<double> x ;
		vector<double> y ;

		double ddummy ;
		unsigned int uidummy ;

		string strsuffix, filename ;
		size_t pos ;

		double dd ;

		unsigned int i, j, k ;

		unsigned int xi, eta ;



		/* file suffix handling for output files */
		strsuffix =  ".chi" ;
		if ( is_intensity_conserved ) { strsuffix += "_tot" ; }
		else { strsuffix += "_avg" ; }


		/* define either phi or s/Q scale and store it in vector x */
		if ( ring_profile_isdef )
		{
			/*
			compute phi in [deg] as in Fit2D,
			phi_k = azi_st + ( k - 0.5 ) * ( azi_end - azi_st ) / azi_bins= azi_st + ( k - 0.5 ) * dd
			(for k = 1 ... azi_bins)

			cause indexing is 0-based here, use k + 0.5 instead !!!
			by default azi_st is 0.0 and azi_end is 360.0, if it is not overridden by the user
			*/
			dd = ( azi_end - azi_st) / ( (double)azi_bins ) ;
			for ( k=0; k<azi_bins; ++k) { x.push_back( azi_st + ( (double)k + 0.5 ) * dd ) ; }
		}
		else
		{
			/*
			compute Q or s in [1/nm] or [1/A], default is Q [1/nm]

			Q_k = Q_min + ( k - 0.5 ) * ( Q_max - Q_min ) / rad_bins = Q_min + ( k - 0.5 ) * dd
			(for k = 1 ... rad_bins)

			cause indexing is 0-based here, use k + 0.5 instead !!!
			by default min_Q is 0.0 and max_Q given via rad_max, if it is not overridden by the user
			*/
			dd = ( max_Q - min_Q ) / ( (double)rad_bins ) ;
			for ( k=0; k<rad_bins; ++k)
			{
				x.push_back( x_scale_fac * ( min_Q + ( (double)k + 0.5 ) * dd ) ) ;
			}
		}


		#pragma omp parallel if ( openmp_isdef )
		{
			if ( omp_in_parallel() )
			{
				#pragma omp single
				fprintf( stdout ,"Parallelized Run: #-CPU=%d, #-Threads=%d\n", omp_get_num_procs(), omp_get_num_threads() ) ;
			}
			else
			{
				fprintf( stdout ,"Serialized Run\n") ;
			}
		}

		/* 
		   fit2dfile_list, and other class members are always shared
		   however append operations on vectors should be not be applied to avoid race conditions !!!

		   tests with serialized runs and parallelized runs (-openmp flag) give exactly the same results
		*/
 		#pragma omp parallel private( j, k, xi, eta, raw_image, image, file, filename, y, ddummy, uidummy, pos) shared( x, strsuffix, OPENMP_STDOUT ) default(none) if ( openmp_isdef )
		{
			/* allocate dynamic arrays within parallel region */
			raw_image = (float*)calloc( azi_bins * rad_bins, sizeof(float) ) ;

			/* make image an rad_bins x azi_bins array allowing easy access to columns for a common radial bin */
			image = (double**)calloc( rad_bins, sizeof(double*) ) ;
			for ( j=0; j <rad_bins; ++j)
			{
				image[j] = (double*)calloc( azi_bins, sizeof(double) ) ;
			}

			/* run Fit2D in batch mode */
 			#pragma omp for schedule(static)
			for ( i=0; i<fit2dfile_list.size(); ++i)
			{
				/* resize y to avoid potential race conditions with append */
				if ( ring_profile_isdef ) { y.resize( azi_bins ) ; }
				else { y.resize( rad_bins ) ; }

				read_tif( fit2dfile_list[i], raw_image, azi_bins, rad_bins) ;

				/*
				transform raw image such that is has the same shape and indexing scheme as displayed in Fit2D,
				intensity values are okay and must not be transformed,
				Fit2D (indexing 1-based, fit2dcorr 0-based) 
				Fit2D coordinates x+1,y+1 -> fit2dcorr x = x' and y = azi_bins - y',
				where x',y' are pixel coordinates from raw_image

				convert from float(32bit) to double(64bit)
				*/
				for ( eta = 0; eta < azi_bins; ++eta)
				{
					for ( xi = 0; xi < rad_bins; ++xi)
					{
						image[xi][eta] = (double)raw_image[ ( azi_bins - 1 - eta ) * rad_bins + xi ];
					}
				}


				if ( ring_profile_isdef ) 
				{
					/* just copy all values (incl. zero / negative) since only one radial bin (uidummy always 1) */
					for ( eta = 0; eta < azi_bins; ++eta) { y[eta] = image[0][eta] ; }
				}
				else
				{
					/* average / integrate azimuthal map */

					/*
					Notes:

					compared with av == 1 or by using azi_bins == 1 in Fit2D Cake/Integrate one 
					observes differences in the averaged intensities, typically for the first 3
					(around beam center) and last 3 points (far most corner) the intensities
					deviate up to 20 %, otherwise there is a fairly good agreement

					presumably, if the beamcenter and the corner are avoided, i.e. there enough
					pixels for the average the problem should not be observed (e.g. it must 
					sufficiently hold rad_st > 0 and rad_end < rad_max )

					the exact reason for the problem is not revealed yet. It must be related
					to the azimuthal mapping / interpolation and number of pixels
					*/
					for ( xi = 0; xi < rad_bins; ++xi)
					{
						ddummy = 0.0 ;
						uidummy = 0 ;

						for ( eta = 0; eta < azi_bins; ++eta)
						{
							if ( image[xi][eta] > 0.0 )
							{
								ddummy += image[xi][eta] ;
								++uidummy ;
							}
						}
						if ( uidummy > 0 )
						{
							/* 
								if intensity conservation YES -> keep sum (sum of all azimuthal bins is independent of azi_bins
								if intensity conservation NO -> take average (all azimuthal bins have nearly the same value, independent on azi_bins)
							*/
							if ( is_intensity_conserved ) { y[xi] = ddummy ; }
							else { y[xi] = ddummy / ( (double)uidummy ) ; } 
						}
						else { y[xi] = 0.0 ; }
					}
				}



				/* replace *.tif(f) with strsuffix for current filename */
				pos = file_list[i].find(".") ;
				filename = file_list[i].substr( 0, pos) + strsuffix ; 


				if ( ( file = fopen( filename.c_str(), "w")) == NULL) { fprintf( stdout, "Error: Cannot write to file %s. Exit.\n", filename.c_str()) ; exit(1) ; }

				/* create header with relevant parameters of the following type:

				LINUX example

				# working_directory=/home/mschmiele/projects/fit2dcorr/test/
				# fit2dfile=/home/mschmiele/projects/fit2dcorr/fit2d/fit2d_12_081_i686_linux2.4.20
				# macrofile=/home/mschmiele/projects/fit2dcorr/fit2d/fit2d_saxs_azi.mac
				#
				# in=000884.tif
				# out=000884.chi
				# mask=mask20110414.msk
				# bcx=277.0 [pix]
				# bcy=91.5 [pix]
				# sdd=293.0 [mm]
				# lambda=1.5418 [Angstroem]
				# pix_size=172.0 x 172.0 [microns x microns]
				# rad_st=0.0 [pix]
				# rad_end=234.1 [pix]
				# rad_bins=234
				# max_2theta=7.8255 [deg]
				# azi_st=0.0 [deg]
				# azi_end=360.0 [deg]
				# azi_bins=360
				# array_size=487 x 360 [pix x pix]
				# pol_fac=0.980
				#
				# cf=0.700 [1/cps]
				# thickness=0.100 [cm]
				# transmission=0.3200
				# exposure time=10800.000 [s]
				#
				# Program called with: /home/mschmiele/projects/fit2dcorr/fit2dcorr -av 0 -err 1 +mask mask20110414.msk +bc 277 91.5 +seq 000869.tif 000889.tif +sdd 293.0 +lambda 1.5418 +pix_size 172.0 172.0 -l 10:1:12 -pol_fac 0.98 -qscale s_nm-1 -abs_units 0.7 0.1 0.32 -azi_st 0.0 -azi_end 360.0 -openmp 
				#
				# Fit2D called with: /home/mschmiele/projects/fit2dcorr/fit2d/fit2d_12_081_i686_linux2.4.20 -key -dim487x360 -fvar#SDD=293.0 -fvar#LAMBDA=1.5418 -fvar#BCX=277.0 -fvar#BCY=91.5 -fvar#PIXSZ_X=172.0 -fvar#PIXSZ_Y=172.0 -fvar#POL_FAC=0.980 -svar#MASK=mask20110414.msk -svar#FILE_IN=000884.tif -svar#FILE_OUT=000884_azi.tif -fvar#AZI_START=0.0 -fvar#AZI_END=360.0 -ivar#AZI_BINS=360 -fvar#PIX_X_ROUT=487.00 -fvar#PIX_Y_ROUT=195.00 -fvar#RAD_IN=0.00 -fvar#RAD_OUT=234.12 -ivar#RAD_BINS=234 -mac/home/mschmiele/projects/fit2dcorr/fit2d/fit2d_saxs_azi.mac
				#
				#     s [1/nm]       I [1/cm]      dI [1/cm]
				1.32395706E-02 1.62458552E+01 9.25921749E-02
				...


				WINDOWS EXAMPLE

				# working_directory=C:\Users\usr\Documents\cpp\fit2dcorr_gcc\test\
				# fit2dfile=C:\Users\usr\Documents\cpp\fit2dcorr_gcc\fit2d\fit2d_12_077_i686_WXP.exe
				# macrofile=C:\Users\usr\Documents\cpp\fit2dcorr_gcc\fit2d\fit2d_saxs_azi.mac
				#
				# in=000884.tif
				# out=000884.chi
				# mask=mask20110414.msk
				# bcx=277.0 [pix]
				# bcy=91.5 [pix]
				# sdd=293.0 [mm]
				# lambda=1.5418 [Angstroem]
				# pix_size=172.0 x 172.0 [microns x microns]
				# rad_st=0.0 [pix]
				# rad_end=234.1 [pix]
				# rad_bins=234
				# max_2theta=7.8255 [deg]
				# azi_st=0.0 [deg]
				# azi_end=360.0 [deg]
				# azi_bins=360
				# array_size=487 x 360 [pix x pix]
				# pol_fac=0.980
				#
				# cf=0.700 [1/cps]
				# thickness=0.100 [cm]
				# transmission=0.3200
				# exposure time=10800.000 [s]
				#
				# Program called with: ..\fit2dcorr.exe -av 0 -err 1 +mask mask20110414.msk +bc 277 91.5 +seq 000869.tif 000889.tif +sdd 293.0 +lambda 1.5418 +pix_size 172.0 172.0 -l 10:1:12 -pol_fac 0.98 -qscale Q_nm-1 -abs_units 0.7 0.1 0.32 -azi_st 0.0 -azi_end 360.0 
				#
				# Fit2D called with: C:\Users\usr\Documents\cpp\fit2dcorr_gcc\fit2d\fit2d_12_077_i686_WXP.exe -key -dim487x360 -fvar#SDD=293.0 -fvar#LAMBDA=1.5418 -fvar#BCX=277.0 -fvar#BCY=91.5 -fvar#PIXSZ_X=172.0 -fvar#PIXSZ_Y=172.0 -fvar#POL_FAC=0.980 -svar#MASK=mask20110414.msk -svar#FILE_IN=000884.tif -svar#FILE_OUT=000884_azi.tif -fvar#AZI_START=0.0 -fvar#AZI_END=360.0 -ivar#AZI_BINS=360 -fvar#PIX_X_ROUT=487.00 -fvar#PIX_Y_ROUT=195.00 -fvar#RAD_IN=0.00 -fvar#RAD_OUT=234.12 -ivar#RAD_BINS=234 -macC:\Users\usr\Documents\cpp\fit2dcorr_gcc\fit2d\fit2d_saxs_azi.mac
				#
				#      Q [1/nm]        I [1/cm]       dI [1/cm]
				8.31866756E-002 1.62458571E+001 9.25922185E-002
				...

				*/

				fprintf( file, "%s working_directory=%s\n", presp.c_str(), working_directory.c_str()) ;
				fprintf( file, "%s fit2dfile=%s\n", presp.c_str(), fit2dfile.c_str()) ;
				fprintf( file, "%s macrofile=%s\n", presp.c_str(), macrofile.c_str()) ;
				fprintf( file, "%s\n", presp.c_str()) ;
				fprintf( file, "%s mask=%s\n", presp.c_str(), maskfile.c_str()) ;
				fprintf( file, "%s bcx=%-.1lf [pix]\n", presp.c_str(), bc[0]) ;
				fprintf( file, "%s bcy=%-.1lf [pix]\n", presp.c_str(), bc[1]) ;
				fprintf( file, "%s sdd=%-.1lf [mm]\n", presp.c_str(), sdd) ;
				fprintf( file, "%s lambda=%-.4lf [Angstroem]\n", presp.c_str(), lambda) ;
				fprintf( file, "%s pix_size=%-.1lf x %-.1lf [microns x microns]\n", presp.c_str(), pix_size[0], pix_size[1]) ;
				fprintf( file, "%s rad_st=%-.1lf [pix]\n", presp.c_str(), rad_st) ;
				fprintf( file, "%s rad_end=%-.1lf [pix]\n", presp.c_str(), rad_end) ;
				fprintf( file, "%s rad_bins=%d\n", presp.c_str(), rad_bins) ;
				fprintf( file, "%s max_2theta=%-.4lf [deg]\n", presp.c_str(), max_2theta) ;
				fprintf( file, "%s azi_st=%-.1lf [deg]\n", presp.c_str(), azi_st) ;
				fprintf( file, "%s azi_end=%-.1lf [deg]\n", presp.c_str(), azi_end) ;
				fprintf( file, "%s azi_bins=%d\n", presp.c_str(), azi_bins) ;
				fprintf( file, "%s array_size=%d x %d [pix x pix]\n", presp.c_str(), array_size[0], array_size[1]) ;
				fprintf( file, "%s pol_fac=%-.3lf\n", presp.c_str(), pol_fac) ;

				fprintf( file, "%s\n", presp.c_str()) ;
				fprintf( file, "%s in=%s\n", presp.c_str(), file_list[i].c_str()) ;
				fprintf( file, "%s out=%s\n", presp.c_str(), fit2dfile_list[i].c_str()) ;
				if ( abs_units_isdef )
				{
					fprintf( file, "%s cf=%-.3lf [1/cps]", presp.c_str(), cf_list[i]) ;
					if ( cf_isauto ) { fprintf( file, " ( phi0=%d [cps] )", phi0_list[i]) ; }
					fprintf( file, "\n") ;
					fprintf( file, "%s thickness=%-.3lf [cm]\n", presp.c_str(), thickness_list[i]) ;
					fprintf( file, "%s transmission=%-.5lf\n", presp.c_str(), transmission_list[i]) ;
					fprintf( file, "%s exptime=%-.3lf [s]\n", presp.c_str(), exptime_list[i]) ;
					fprintf( file, "%s\n", presp.c_str()) ;
				}
				if ( subtract_isdef )
				{
					fprintf( file, "%s background=%s\n", presp.c_str(), file_b.c_str()) ;
					fprintf( file, "%s vol_fract_b=%-.3lf\n", presp.c_str(), vol_fract_b) ;
					if ( abs_units_isdef )
					{
						fprintf( file, "%s cf_b=%-.3lf [1/cps]", presp.c_str(), cf_b) ;
						if ( cf_b_isauto ) { fprintf( file, " ( phi0_b=%d [cps] )", phi0_b) ; }
						fprintf( file, "\n") ;
						fprintf( file, "%s thickness_b=%-.3lf [cm]\n", presp.c_str(), thickness_b) ;
						fprintf( file, "%s transmission_b=%-.5lf\n", presp.c_str(), transmission_b) ;
						fprintf( file, "%s exptime_b=%-.3lf [s]\n", presp.c_str(), exptime_b) ;
						fprintf( file, "%s\n", presp.c_str()) ;
					}
				}
				fprintf ( file, "%s skipped lines=", presp.c_str()) ;
				if ( linelist.empty() ) { fprintf( file, "%s", "none") ; }
				for ( k=0; k<linelist.size(); ++k) 
				{ 
					fprintf ( file, "%d", linelist[k]) ;
					if ( k < linelist.size()-1 ) { fprintf ( file, "%s", ", ") ; }
				}
				fprintf( file, "\n") ;
				fprintf( file, "%s\n", presp.c_str()) ;
				fprintf( file, "%s Program called with: ", presp.c_str()) ;
				/* Export function call parameters ecarg */
				for ( k=0; k<ecarg; ++k) { fprintf( file, "%s ", evarg[k]) ; }
				fprintf( file, "\n%s\n", presp.c_str()) ;

				fprintf( file, "%s Fit2D called with: %s\n", presp.c_str(), fit2dcall_list[i].c_str()) ;
				fprintf( file, "%s\n", presp.c_str()) ;


				/* write data */
				if ( ring_profile_isdef )
				{
					fprintf( file, "%s    %s%s [deg]       %s%s      %sd%s\n", presp.c_str(), exp_str, "phi", exp_str, y_scale.c_str(), exp_str, y_scale.c_str() ) ;

					for ( eta=0; eta<azi_bins; ++eta) { fprintf( file, "%.8E %.8E %.8E\n", x[eta], y[eta], -1.0) ; }
				}
				else
				{
					fprintf( file, "%s     %s%s       %s%s      %sd%s\n", presp.c_str(), exp_str, x_scale.c_str(), exp_str, y_scale.c_str(), exp_str, y_scale.c_str()) ;

					for ( xi=0; xi<rad_bins; ++xi) { fprintf( file, "%.8E %.8E %.8E\n", x[xi], y[xi], -1.0) ; }
				}


				/* close file */
				fclose(file) ;

				/* clear y */
				y.clear() ;

				/* remove temporary *.tif-files */
				if ( !verbose_isdef ) 
				{
					if ( remove( fit2dfile_list[i].c_str() ) != 0 ) {fprintf( stdout, "Error: Could not delete file %s. Exit.\n", fit2dfile_list[i].c_str()) ; exit(1) ; }
				}
			}

			/* free (raw_)image */
			free(raw_image) ;

			for ( j=0; j<rad_bins; ++j)
			{
				free(image[j]) ;
			}
			free(image) ;
		}
	}


	/* edit the chi-files generated from Fit2D for av == 0 & 1 in case of err == 1 to generate the errorbars */
	void correct_chi_files()
	{
		FILE *file[5] ;
		string filename[5] ;
		char *strptr1, *strptr2 ;

		char data[defsigns] ; /* dummy char array to store data */
		char sdummy[4][defsigns] ; /* dummy variable to store all read lines fro the files */
		char x_string[4][defsigns] ;

		double y[4] ;

		unsigned int n, N, N_del ; /* N - number of temporary chi-files to read, N_del - number of temporary chi-files to delete */
		unsigned int i, j ;
		unsigned int lineindex, writtenlineindex ;
		bool line_isnegative ;

		string strdummy ;

		#pragma omp parallel if ( openmp_isdef )
		{
			if ( omp_in_parallel() )
			{
				#pragma omp single
				fprintf( stdout ,"Parallelized Run: #-CPU=%d, #-Threads=%d\n", omp_get_num_procs(), omp_get_num_threads() ) ;
			}
			else
			{
				fprintf( stdout ,"Serialized Run\n") ;
			}
		}


		/* write chi file, copy all #-lines and write 1st and 2nd from *.chi_avg and calculate 3rd column from *.chi_avg and *.chi_tot */


		/* 
		   fileist, thickness_list, transmission_list, exptime_list and other class members are always shared
		   however append operations on vectors should be not be applied to avoid race conditions !!!

		   tests with serialized runs and parallelized runs (-openmp flag) give exactly the same results
		*/
		#pragma omp parallel private( filename, file, N, N_del, x_string, y, lineindex, writtenlineindex, line_isnegative, sdummy, strptr1, strptr2, data, n, j, strdummy ) shared( OPENMP_STDOUT ) default(none) if ( openmp_isdef )
		{
			/* batch mode */
			#pragma omp for schedule(static)
			for ( i=0; i<file_list.size(); ++i)
			{
				/* filename handling */
				strdummy = file_list[i].substr( 0, file_list[i].find(".") ) ;

				/* define N and N_del depending on err and subtract_isdef */
				N = 1 ;
				/* write i-th *.chi file */
				filename[0] = strdummy + ".chi" ;

				/* read i-th *.chi_avg/tot-files */
				filename[1] = strdummy + ".chi_avg" ;
				if ( err == 0 )
				{
					if ( subtract_isdef ) // N == 2
					{
						N = N + 1 ;
						strdummy = file_b.substr( 0, file_b.find(".") ) ;
						filename[2] = strdummy + ".chi_avg" ;
					}
					// else N == 1
				}
				else
				{
					N = N + 1 ; // N == 2
					filename[2] = strdummy + ".chi_tot" ;

					if ( subtract_isdef ) // N == 4
					{
						N = N + 2 ;
						strdummy = file_b.substr( 0, file_b.find(".") ) ;
						filename[3] = strdummy + ".chi_avg" ;
						filename[4] = strdummy + ".chi_tot" ;
					}
				}

				/* w/o subtraction N == N_del, otherwise it is half of it */
				N_del = N ;
				if ( subtract_isdef ) { N_del /= 2 ; }



				for ( j=1; j<=N; ++j)
				{
					if ( ( file[j] = fopen( filename[j].c_str(), "r")) == NULL )
					{ 
						fprintf( stdout, "Error: Cannot read file %s. File does not exist. Exit.\n", filename[j].c_str() ) ;
						exit(1) ;
					}
				}

				if ( ( file[0] = fopen( filename[0].c_str(), "w")) == NULL )
				{
					fprintf( stdout, "Error: Cannot write to file %s. Exit.\n", filename[0].c_str() ) ; 
					exit(1) ;
				}


				lineindex = 0 ;
				writtenlineindex = 0 ;
				while ( fgets( sdummy[0], defsigns, file[1]) != NULL )
				{
					++lineindex ;

					for ( j=1; j<N; ++j)
					{
						if ( fgets( sdummy[j], defsigns, file[j+1]) == NULL )
						{
							fprintf( stdout, "Error: Line %d in file %s does not exist as in corresponding file %s. Exit.\n", lineindex, filename[j].c_str(), filename[1].c_str() ) ;
							exit(1) ;
						}
					}

					/* copy the #-comment lines (incl EOL) from sample's '*.chi_avg file to output file */
					if ( sdummy[0][0] == '#' )
					{
						fprintf( file[0], "%s", sdummy[0] ) ;
						continue ;
					}

					/* process line by line the data (read q and I / I_tot ) in all N files */
					for ( j=0; j<N; ++j)
					{
						strptr1 = sdummy[j] ;

						/* first column */
						while (*strptr1 == ' ') { ++strptr1 ; }

						if ( ( strptr2 = strchr( strptr1, ' ')) == NULL ) { fprintf( stdout, "Error: Expected three columns in file %s, only one is present. Exit.\n", filename[j+1].c_str()) ; exit(1) ; }

						n = (unsigned int)(strptr2-strptr1) ;
						memmove( x_string[j], strptr1, n) ;
						x_string[j][n] = 0 ;

						/* second column */
						strptr1 = strptr2 ;
						while (*strptr1 == ' ') { ++strptr1 ; }

						if ( ( strptr2 = strchr( strptr1, ' ')) == NULL ) { fprintf( stdout, "Error: Expected three columns in file %s, only two are present. Exit.\n", filename[j+1].c_str()) ; exit(1) ; }

						n = (unsigned int)(strptr2-strptr1) ;
						memmove( data, strptr1, n) ;
						data[n] = 0 ;
						if ( is_numeric(data) ) { y[j] = strtod( data, NULL) ; }
						else { fprintf( stdout, "Error: Found entry in the second column in line %d in file %s with non-numeric type. Exit.\n", lineindex, filename[j+1].c_str()) ; exit(1) ; }
					}

					/* check if q-values are same in all temporary chi-files */
					for ( j=1; j<N; ++j)
					{
						if ( strcmp( x_string[0], x_string[j] ) )
						{
							fprintf( stdout, "Error: x-values at line %d in file %s and the corresponding file %s do not match. Exit.\n", lineindex, filename[j+1].c_str(), filename[1].c_str() ) ;
						 	exit(1) ;
						}
					}

					/*
						POST PROCESSING

						absolute units
						errorbars
						(weighted) background subtraction

						skipping negative (<=0) lines if requested (in case of background subtraction skip lines where sample, background and difference are <=0 )
						skipping user-defined lines if requested
					 */

					/* omit intensities values <= 0 if requested */
					if ( nonnegative_isdef )
					{
						line_isnegative = false ;
						for ( j=0; j<N; ++j) { if ( y[j]<=0.0 ) { line_isnegative = true ; break ; } }
						if ( line_isnegative ) { continue ; }
					}


					/* calculate errobars for sample / background and save in y[1] / y[b+1] */
					if ( err == 1 )
					{
						/* avoid divsion by zero if y[1] == 0 (implies usually y[0] == 0, too) */
						if ( y[1]==0.0 ) { y[1] = 1.0 ; }

						/* write dI to y[1] */
						y[1] = y[0] / sqrt(y[1]) ;

						if ( subtract_isdef )
						{
							if ( y[3]==0.0 ) { y[3] = 1.0 ; }
							y[3] = y[2] / sqrt(y[3]) ;
						}
					}


					/* 
						update y[0] and y[1] in case of abs_units_isdef and/or subtract_isdef
					*/
					if ( abs_units_isdef )
					{
						get_abs_intensity( y[0], cf_list[i], thickness_list[i], transmission_list[i], exptime_list[i]) ;
						if ( err == 1 )
						{
							get_abs_intensity( y[1], cf_list[i], thickness_list[i], transmission_list[i], exptime_list[i]) ;
						}

						if ( subtract_isdef )
						{
							if ( err == 0 )
							{
								get_abs_intensity( y[1], cf_b, thickness_b, transmission_b, exptime_b) ;
								y[0] = y[0] - vol_fract_b * y[1] ;
							}
							else
							{
								get_abs_intensity( y[2], cf_b, thickness_b, transmission_b, exptime_b) ;
								y[0] = y[0] - vol_fract_b * y[2] ;

								get_abs_intensity( y[3], cf_b, thickness_b, transmission_b, exptime_b) ;
								y[1] = sqrt( y[1] * y[1] + vol_fract_b * vol_fract_b * y[3] * y[3] ) ;
							}
						}
					}
					else
					{
						if ( subtract_isdef )
						{
							if ( err == 0 )
							{
								y[0] = y[0] - vol_fract_b * y[1] ;
							}
							else
							{
								y[0] = y[0] - vol_fract_b * y[2] ;
								y[1] = sqrt( y[1] * y[1] + vol_fract_b * vol_fract_b * y[3] * y[3] ) ;
							}
						}
					}

					/* in case of err == 0, y[1] needs to stay untouched at -1 (it is not reasd from the file(s)) */ 
					if ( err == 0 ) { y[1] = -1.0 ; }


					/* in case of background subtraction check again difference */
					if ( nonnegative_isdef && subtract_isdef && ( y[0] <= 0.0 ) ) { continue ; }

					++writtenlineindex ;

					/* check if the line should be skipped */
					if ( find( linelist.begin(), linelist.end(), writtenlineindex ) != linelist.end() ) { continue ; }

					fprintf( file[0], "%s %.8E %.8E\n", x_string[0], y[0], y[1] ) ;
				}


				/* close all N+1 files and - if not in verbose mode - remove temporary sample *.chi-files (remove background files (>N/2) in the last run) */
				for ( j=0; j<=N; ++j) { fclose(file[j]) ; }

				if ( !verbose_isdef ) 
				{
					if ( ( i == file_list.size() - 1 ) && subtract_isdef ) { N_del *= 2 ; }

					for ( j=1; j<=N_del; ++j)
					{
						if ( remove( filename[j].c_str() ) != 0 ) {fprintf( stdout, "Error: Could not delete file %s. Exit.\n", filename[j].c_str()) ; exit(1) ; }
					}
				}
			}
		}
	}


	/* returns image_size from tif-file */
	void info_tif( string filename, unsigned int &width, unsigned int &height)
	{
		TIFF* tif = TIFFOpen( filename.c_str(), "r") ;
		if (tif)
		{
			/* image dimensions */
			TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height) ;
			TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width) ;
		}
		else 
		{ 
			fprintf( stdout, "Error: Cannot open file %s. Exit.\n", filename.c_str()) ;
			exit(1) ;
		}
		TIFFClose(tif) ;
	}


// 	template <class T>
// 	void read_tif( string filename, T* buf, unsigned int &height, unsigned int &width)
// 	{
// 		/* NOTE: 
// 		   in Fit2D float-tif export dialogue
// 		   Support for floating point image data is not part of "baseline TIFF", 
// 		   many programs may input this data incorrectly !
// 		*/
// 		TIFF* tif = TIFFOpen( filename.c_str(), "r") ;
// 		if (tif)
// 		{
// 			uint32 row ;
// 			uint32 col ;
// 
// 			/* image dimensions */
// 			TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height) ;
// 			TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width) ;
// 
// 			uint16 bps, spp, bpp, photo, int_float ;
// 
// 			/* TIFFTAG_BITSPERSAMPLE
// 			   32 bit tif for PILATUS and for Fit2D exported tif's
// 			   http://www.awaresystems.be/imaging/tiff/tifftags/bitspersample.html
// 			*/
// 			TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps) ;
// 
// 			/* TIFFTAG_SAMPLESPERPIXEL, 
// 			   not provided in the tif-files -> default value 0 is returned
// 			   generally:
// 			   1 for grayscale, b/w, Hue
// 			   3 for RGB
// 			   http://www.awaresystems.be/imaging/tiff/tifftags/samplesperpixel.html
// 			*/
// 			TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp) ;
// 
// 			bpp = bps * spp;
// 
// 			/* 
// 			   TIFFTAG_PHOTOMETRIC
// 			   0 = WhiteIsZero. For bilevel and grayscale images: 0 is imaged as white. (for tif's exported from Fit2D, ???)
// 			   1 = BlackIsZero. For bilevel and grayscale images: 0 is imaged as black. (for PILATUS tif)
// 			   ...
// 			   http://www.awaresystems.be/imaging/tiff/tifftags/photometricinterpretation.html
// 			*/
// 			TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &photo) ;
// 
// 			/* 
// 			   TIFFTAG_SAMPLEFORMAT
// 			   1 = unsigned integer data
// 			   2 = two's complement signed integer data  (for PILATUS tif)
// 			   3 = IEEE floating point data (for float tif's exported from Fit2D)
// 			   4 = undefined data format 
// 			   ...
// 			   http://www.awaresystems.be/imaging/tiff/tifftags/sampleformat.html
// 			*/
// 			TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &int_float) ;
// 
// 			/* get type of array buf and compare with type of the tif-file
// 			cout << "type " << typeid(buf).name() << endl ;
// 			
// 			cout << "type " << typeid(int*).name() << endl ;
// 			cout << "type " << typeid(float*).name() << endl ;
// 			cout << "type " << typeid(int).name() << endl ;
// 			cout << "type " << typeid(float).name() << endl ;
// 			*/
// 			if ( !strcmp( typeid(buf).name(), typeid(int*).name() ) && ( int_float != 2 ) ) { fprintf( stdout, "Error: Expect an int-array for tif-file %s. Exit.", filename.c_str()) ; exit(1) ; } 
// 			if ( !strcmp( typeid(buf).name(), typeid(float*).name() ) && ( int_float != 3 ) ) { fprintf( stdout, "Error: Expect a float-array for tif-file %s. Exit.", filename.c_str()) ; exit(1) ; } 
// 
// 			/* linesize = width * sizeof(pixval), for 32Bit(=4Byte) tif -> 4 * width */
// 			/* int linesize = TIFFScanlineSize(tif) ; */
// 
// // 			int32* buf ;
// // 			buf = (int32*)malloc( height * linesize ) ;
// // 			int32 maxval = 0 ;
// // 			int32 minval = 0 ;
// // 			float* buf ;
// // 			buf = (float*)malloc( height * linesize ) ;
// // 			float maxval = 0.0 ;
// // 			float minval = 0.0 ;
// 
// // 			buf = (T*)malloc( height * linesize ) ;
// 			T maxval = 0 ;
// 			T minval = 0 ;
// 
// 			for ( row = 0; row < height; row++)
// 			{
// 				TIFFReadScanline( tif, &buf[row * width], row, 0) ;
// 				printf( "%d\n", row) ;
// 				for ( col = 0; col < width; col++)
// 				{
// 					if ( buf[row * width+col] > maxval ) { maxval = buf[row * width+col] ; }
// 					if ( buf[row * width+col] < minval ) { minval = buf[row * width+col] ; }
// 					printf( "%f ", (double)buf[row * width+col]) ;
// // 					printf( "%d ", buf[row * width+col]) ;
// 				}
// 				printf("\n") ;
// 			}
// 
// 			/* for char buf -> transform later to int32 or float */
// // 			char* buf ;
// // 			buf = (char*)malloc( height * linesize ) ;
// // 			for (row = 0; row < height; row++)
// // 			{
// // 				printf("%d\n", row);
// // 				TIFFReadScanline(tif, &buf[row * linesize], row, 0);
// // 				for ( col = 0; col < linesize; col++)
// // 				{
// // 					if ( buf[row * linesize+col] > maxval ) { maxval = buf[row * linesize+col] ; }
// // 					if ( buf[row * linesize+col] < minval ) { minval = buf[row * linesize+col] ; }
// // 					printf("%d ", buf[row * linesize+col]);
// // 				}
// // 				printf("\n") ;
// // 			}
// 
// // 			/* int images TIFFTAG_SAMPLEFORMAT == 2 */
// // 			int* data = (int*) buf ;
// // 			maxval = 0 ;
// // 			minval = 0 ;
// // 
// // 			for (row = 0; row < height; row++)
// // 			{
// // 				printf("%d\n", row);
// // 				for (col = 0; col < width; col++)
// // 				{
// // 					if ( data[row * width+col] > maxval ) { maxval = data[row * width+col] ; }
// // 					if ( data[row * width+col] < minval ) { minval = data[row * width+col] ; }
// // 					printf("%d ", data[row * width+col]);
// // 				}
// // 				printf("\n");
// // 			}
// 
// 			/* <=32bit float images TIFFTAG_SAMPLEFORMAT == 3 */
// // 			float* data = (float*) buf ;
// // 			float maxval_f = 0.0 ;
// // 			float minval_f = 0.0 ;
// // 
// // 			for (row = 0; row < height; row++)
// // 			{
// // 				printf("%d\n", row);
// // 				for (col = 0; col < width; col++)
// // 				{
// // 					if ( data[row * width+col] > maxval_f ) { maxval_f = data[row * width+col] ; }
// // 					if ( data[row * width+col] < minval_f ) { minval_f = data[row * width+col] ; }
// // 					printf("%f ", data[row * width+col]);
// // 				}
// // 				printf("\n");
// // 			}
// 			TIFFClose(tif) ;
//  		}
//  		else { fprintf( stdout, "Error: Cannot open tif-file %s. Exit.", filename.c_str()) ; exit(1) ; }
// 	}


	/*
	   read tif-file and return via buf
	   apply checks:
	   -first, if image fits exactly into buf (compare with given values height, width)
	   -check for correct data format matching between buf and tif format (unsigned int, int(4 byte, 32 bit), float(4 byte, 32 bit))
	   -double i.e. 64 bit tif-files are not supported
	*/
	template <class T>
	void read_tif( string filename, T* buf, unsigned int height, unsigned int width)
	{
		/* NOTE:
		   in Fit2D float-tif export dialogue
		   Support for floating point image data is not part of "baseline TIFF",
		   many programs may input this data incorrectly !
		*/
		TIFF* tif = TIFFOpen( filename.c_str(), "r") ;
		if (tif)
		{
			unsigned int uidummy ;

			/* image dimensions */
			TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &uidummy) ;
			if ( height != uidummy )
			{
				fprintf( stdout, "Error: tif-file %s has height = %u and not %u. Exit.\n", filename.c_str(), uidummy, height) ;
				exit(1) ;
			}
			TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &uidummy) ;
			if ( width != uidummy )
			{
				fprintf( stdout, "Error: tif-file %s has width = %u and not %u. Exit.\n", filename.c_str(), uidummy, width) ;
				exit(1) ;
			}

			uint16 bps, spp, /* bpp,*/ photo, int_float ;

			/* TIFFTAG_BITSPERSAMPLE
			   32 bit tif for PILATUS and for Fit2D exported tif's
			   http://www.awaresystems.be/imaging/tiff/tifftags/bitspersample.html
			*/
			TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps) ;

			/* TIFFTAG_SAMPLESPERPIXEL, 
			   not provided in the tif-files from PILATUS and Fit2D -> default value 0 is returned
			   generally:
			   1 for grayscale, b/w, Hue
			   3 for RGB
			   http://www.awaresystems.be/imaging/tiff/tifftags/samplesperpixel.html
			*/
			TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp) ;

			/* bpp = bps * spp ; */

			/* 
			   TIFFTAG_PHOTOMETRIC
			   0 = WhiteIsZero. For bilevel and grayscale images: 0 is imaged as white. (for tif's exported from Fit2D, ???)
			   1 = BlackIsZero. For bilevel and grayscale images: 0 is imaged as black. (for PILATUS tif)
			   ...
			   http://www.awaresystems.be/imaging/tiff/tifftags/photometricinterpretation.html
			*/
			TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &photo) ;

			/* 
			   TIFFTAG_SAMPLEFORMAT
			   1 = unsigned integer data
			   2 = two's complement signed integer data  (for PILATUS tif)
			   3 = IEEE floating point data (for float tif's exported from Fit2D)
			   4 = undefined data format 
			   ...
			   http://www.awaresystems.be/imaging/tiff/tifftags/sampleformat.html
			*/
			TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &int_float) ;

			/* get type of array buf and compare with type of the tif-file */
			if ( !strcmp( typeid(buf).name(), typeid(unsigned int*).name() ) && ( int_float != 1 ) )
			{
				fprintf( stdout, "Error: Expected an unisgned int array for tif-file %s. Exit.\n", filename.c_str()) ;
				exit(1) ;
			}
			if ( !strcmp( typeid(buf).name(), typeid(int*).name() ) && ( int_float != 2 ) )
			{
				fprintf( stdout, "Error: Expected an int array for tif-file %s. Exit.\n", filename.c_str()) ;
				exit(1) ;
			}
			if ( !strcmp( typeid(buf).name(), typeid(float*).name() ) && ( int_float != 3 ) )
			{
				fprintf( stdout, "Error: Expected a float array for tif-file %s. Exit.\n", filename.c_str()) ;
				exit(1) ;
			} 

			/* read tif-file into buf */
			for ( uint32 row = 0; row < height; row++)
			{
				TIFFReadScanline( tif, &buf[row * width], row, 0) ;
			}
			TIFFClose(tif) ;
 		}
 		else { fprintf( stdout, "Error: Cannot open file %s. Exit.\n", filename.c_str()) ; exit(1) ; }
	}



	/* av == 2 */
	void process()
	{
		/*
		   without using Fit2D, not yet implemented

		   necessary steps:
		   read (PILATUS) tiff-files -> read_PILATUS_tiff()
		   read mask -> read_Fit2D_mask()
		   perform azimuthal average and generate errorbars -> azimuthal_averaging()
		   absolute intensity calibration,
		   export as chi-file
		*/
	}


	/*
	   read PILATUS tiff-file
	*/
	template <class T>
	void read_PILATUS_tiff( string filename, T* buf, unsigned int height, unsigned int width)
	{

	}

	/*
	   read Fit2D mask file
	*/
	void read_Fit2D_mask( string filename, bool* mask, unsigned int height, unsigned int width)
	{

	}

	/* azimuthal averaging without Fit2D (av == 2) */
	void azimuthal_averaging( )
	{

	}

} ;


int main (int carg, char **varg)
{
	fit2dcorr job = fit2dcorr() ;

	/* read and derive parameters from the cmd */
	job.eval_cmd( (unsigned int)carg, varg) ;

	/* run azimuthal averaging */
	job.run() ;

	return (0) ;
}
