NAME

       fnl_mapmaking

       Generating non-Gaussian simulations of the Cosmic Microwave
	   Background (CMB)


DESCRIPTION

       The purpose of this algorithm is to simulate temperature and
	   polarization realizations of the CMB that contain non-Gaussianity
	   of local type in the Healpix pixelization scheme. For a full
	   description, see https://arxiv.org/abs/0909.0009. The code is
	   written in Fortran and hybrid MPI/OpenMP parallelized for an
	   efficient execution on computer clusters.


INSTALLATION

       Prerequisites: - Healpix
                      - Intel Math Kernel Library
                      - CFITSIO library

       Prior to compilation, please edit the 'Makefile' to reflect
	   your MPI compiler wrapper. The path to all relevant components
	   of the Intel MKL, including lapack95 and blas95, should be
	   adapted. Additionally, the Healpix and CFITSIO libraries must
	   be available for linking.


USAGE

       The algorithm reads keyword - keyvalue pairs from a namelist in
       the configuration file located at './config/input.cfg'. Here,
	   the user is asked to specify

	   NR_SIMS:
	      The number of simulations to generate.

	   NR_SHELLS:
	      The number of radial nodes used in the line of sight
		  integration of curvature perturbations.

	   NSIDE:
	      Resolution parameter of Healpix maps.

	   LMAX:
	      Maximum multipole moment of the simulated maps.

	   INCLUDE_POL:
	      Boolean flag indicating if simulations should also include
		  polarization.

	   DIR_IN_W8RING:
          Directory of the Healpix pixel weights, typically
          $HEALPIX/data. Only the directory should be specified, the
          file name itself will be appended according to the Healpix
          naming convention.

	   FILE_PREFIX_COSMO:
	      Abbreviation describing the cosmology assumed in the
		  simulation.

       FILE_PREFIX_NG_AUX:
	      Directory containing all auxiliary input files specified 
		  below.

	   FILE_IN_WEIGHTS:
          File containing the numerical quadrature weights used in
		  the line of sight integration of curvature perturbations.
		  The string given by 'NR_SHELLS.fits' will be automatically
		  appended.

	   FILE_IN_WEIGHTS_MU:
          File prefix containing the numerical quadrature weights that
		  can be used to compute mu perturbations. The string given
		  by 'NR_SHELLS.fits' will be automatically appended.

       FILE_IN_COV:
          File prefix for the covariance matrices of primordial 
		  curvature perturbations. The string given by
		  'NR_SHELLS.fits' will be automatically appended.

	   FILE_IN_COV_MU:
          File prefix for the covariance matrices of mu perturbations. 
		  The string given by 'NR_SHELLS.fits' will be automatically 
		  appended.

	   FILE_OUT_SIGNAL_L:
          File naming convention for output files containing the linear
		  (Gaussian) component of the CMB realization. The
		  placeholder '%MCNUMBER%' will be replaced by the simulation
		  number.

	   FILE_OUT_SIGNAL_NL:
          File naming convention for output files containing the
		  non-Gaussian CMB map of local type. The placeholder
		  '%MCNUMBER%' will be replaced by the current simulation 
		  number.

	   FILE_OUT_SIGNAL_NNL:
          File naming convention for output files containing the
		  CMB maps for the cubic term. The placeholder '%MCNUMBER%'
		  will be replaced by the current simulation number.

	   FILE_OUT_SIGNAL_MU:
          File naming convention for output files containing maps of
		  mu perturbations. The placeholder '%MCNUMBER%' will be
		  replaced by the current simulation number.


RESTRICTIONS

       Example files to compute non-Gaussian maps for Planck cosmology
	   are provided in the 'data' directory.

	   We do not provide inputs to generate mu perturbations!
	   Corresponding code sections have to be commented out in the
	   source code prior to compilation.

	   Please cite https://arxiv.org/abs/0909.0009 if you make use of
       this algorithm in scientific publications.
