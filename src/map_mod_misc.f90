!
!     Copyright (C) 2009 - 2015 Franz Elsner
!
!     This file is part of the f_nl map making package.
!     See http://arxiv.org/abs/0909.0009 for details.
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <https://www.gnu.org/licenses/>.
!


MODULE MAP_MOD_MISC

  USE MPI_F08
  USE HEALPIX_TYPES
  USE MAP_MOD_CLASSES

  IMPLICIT NONE


CONTAINS

!-----------------------------------------------------------------------
! Parse input parameters
!-----------------------------------------------------------------------

  SUBROUTINE PARSE_INPUT_PARAMETERS(parameter_file, files, par)

    USE PIX_TOOLS,     ONLY: NSIDE2NPIX

    CHARACTER(LEN=*),                       INTENT(IN)    :: parameter_file
    TYPE(FILENAMES),                        INTENT(INOUT) :: files
    TYPE(PARAMETERS),                       INTENT(INOUT) :: par

    LOGICAL(LGT)                                          :: include_pol
    CHARACTER(LEN=FNL)                                    :: filename
    CHARACTER(LEN=FNL)                                    :: dir_in_w8ring
    CHARACTER(LEN=FNL)                                    :: file_prefix_cosmo
    CHARACTER(LEN=FNL)                                    :: file_prefix_ng_aux
    CHARACTER(LEN=FNL)                                    :: file_in_weights
    CHARACTER(LEN=FNL)                                    :: file_in_weights_mu
    CHARACTER(LEN=FNL)                                    :: file_in_cov
    CHARACTER(LEN=FNL)                                    :: file_in_cov_mu
    CHARACTER(LEN=FNL)                                    :: file_out_signal_l
    CHARACTER(LEN=FNL)                                    :: file_out_signal_nl
    CHARACTER(LEN=FNL)                                    :: file_out_signal_nnl
    CHARACTER(LEN=FNL)                                    :: file_out_signal_mu
    INTEGER(I4B)                                          :: nr_sims
    INTEGER(I4B)                                          :: nr_shells
    INTEGER(I4B)                                          :: nside
    INTEGER(I4B)                                          :: lmax

    NAMELIST / input_parameters / nr_sims, nr_shells, nside, lmax,&
         & include_pol, dir_in_w8ring, file_prefix_cosmo,&
         & file_prefix_ng_aux, file_in_weights, file_in_weights_mu,&
         & file_in_cov, file_in_cov_mu, file_out_signal_l,&
         & file_out_signal_nl, file_out_signal_nnl, file_out_signal_mu

    IF (COMMAND_ARGUMENT_COUNT() .EQ. 1) THEN
       CALL GET_COMMAND_ARGUMENT(1, filename)
       WRITE(*,'(\3X, "Loading parameters from:", 1X, A)')&
            & TRIM(filename)
    ELSE
       filename = parameter_file
    ENDIF

    OPEN(UNIT=10, FILE=filename, STATUS='OLD', ACTION='READ')
    READ(10, NML=input_parameters)
    CLOSE(10)

    ALLOCATE(par%radii__r(1:nr_shells))
    CALL par%INIT()

    par%pol       = include_pol
    par%nr_sims   = nr_sims
    par%nside     = nside
    par%npix      = NSIDE2NPIX(nside)
    par%lmax      = lmax
    par%mmax      = par%lmax
    par%nr_shells = nr_shells
    par%radii__r  = 0.0_DP

    files%in_hpx_data      = dir_in_w8ring
    files%prefix_cosmo     = file_prefix_cosmo
    files%prefix_ng_aux    = file_prefix_ng_aux
    files%in_weights       = file_in_weights
    files%in_weights_mu    = file_in_weights_mu
    files%in_covariance    = file_in_cov
    files%in_covariance_mu = file_in_cov_mu
    files%out_signal_l     = file_out_signal_l
    files%out_signal_nl    = file_out_signal_nl
    files%out_signal_nnl   = file_out_signal_nnl
    files%out_signal_mu    = file_out_signal_mu

    CALL par%sht%INIT(par)

    CALL READ_W8RINGS(par, files)


  END SUBROUTINE PARSE_INPUT_PARAMETERS

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Read quadrature weights from file
!-----------------------------------------------------------------------

  SUBROUTINE READ_W8RINGS(par, files)

    USE FITSTOOLS,     ONLY: READ_DBINTAB

    TYPE(FILENAMES),                        INTENT(IN)    :: files
    TYPE(PARAMETERS),                       INTENT(INOUT) :: par

    CHARACTER(LEN=FNL)                                    :: filename
    CHARACTER(LEN=5)                                      :: nside_string
    LOGICAL(LGT)                                          :: anynull
    REAL(DP)                                              :: nullval

    IF (par%thread_id .EQ. par%root) THEN

       WRITE(nside_string,'(I5.5)') par%nside

       filename = TRIM(files%in_hpx_data) // '/weight_ring_n'&
            & // nside_string // '.fits'

       CALL READ_DBINTAB(filename, par%sht%w8ring__ring_IQU,&
            & 2*par%nside, 1, nullval, anynull)

       par%sht%w8ring__ring_IQU = 1.0_DP + par%sht%w8ring__ring_IQU

    ENDIF

    CALL MPI_BCAST(par%sht%w8ring__ring_IQU,&
         & SIZE(par%sht%w8ring__ring_IQU), MPI_DOUBLE_PRECISION,&
         & par%root, MPI_COMM_WORLD, par%mpi_error)
    CALL CHECK_ERROR_FLAG("MPI_BCAST", par%mpi_error)


  END SUBROUTINE READ_W8RINGS

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------

  SUBROUTINE ALLOCATE_ARRAYS(weights_p, weights_mu, cov_p, cov_mu,&
       & rng, signal_l, signal_nl, signal_nnl, signal_mu, par)

    TYPE(PARAMETERS),                       INTENT(IN)    :: par
    TYPE(RND_SEEDS),                        INTENT(INOUT) :: rng
    TYPE(WEIGHTS),                          INTENT(INOUT) :: weights_p
    TYPE(WEIGHTS),                          INTENT(INOUT) :: weights_mu
    TYPE(COVARIANCE),                       INTENT(INOUT) :: cov_p
    TYPE(COVARIANCE),                       INTENT(INOUT) :: cov_mu
    TYPE(SKYMAP),                           INTENT(INOUT) :: signal_l
    TYPE(SKYMAP),                           INTENT(INOUT) :: signal_nl
    TYPE(SKYMAP),                           INTENT(INOUT) :: signal_nnl
    TYPE(SKYMAP),                           INTENT(INOUT) :: signal_mu

    CALL rng%INIT(par)
    CALL weights_p%INIT(par)
    CALL weights_mu%INIT(par)
    CALL cov_p%INIT(par)
    CALL cov_mu%INIT(par)
    CALL signal_l%INIT(par)
    CALL signal_nl%INIT(par)
    CALL signal_nnl%INIT(par)
    CALL signal_mu%INIT(par)


  END SUBROUTINE ALLOCATE_ARRAYS

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Read auxiliary files
!-----------------------------------------------------------------------

  SUBROUTINE READ_AUX_FILES(weights_p, weights_mu, cov_p, cov_mu,&
       & rng, files, par)

    USE RNGMOD,        ONLY: PLANCK_RNG, RAND_INIT, RAND_UNI

    TYPE(FILENAMES),                        INTENT(IN)    :: files
    TYPE(PARAMETERS),                       INTENT(INOUT) :: par
    TYPE(RND_SEEDS),                        INTENT(INOUT) :: rng
    TYPE(WEIGHTS),                          INTENT(INOUT) :: weights_p
    TYPE(WEIGHTS),                          INTENT(INOUT) :: weights_mu
    TYPE(COVARIANCE),                       INTENT(INOUT) :: cov_p
    TYPE(COVARIANCE),                       INTENT(INOUT) :: cov_mu

    TYPE(PLANCK_RNG)                                      :: rng_handle
    INTEGER(I4B)                                          :: i
    INTEGER(I4B)                                          :: rnd_seed

    IF (par%thread_id .EQ. par%root) THEN

       CALL SYSTEM_CLOCK(COUNT = rnd_seed)

       CALL RAND_INIT(rng_handle, rnd_seed)
       DO i = 1,rng%nr_sims
          rng%seeds__map(i)&
               & = INT(MAX_I4B*RAND_UNI(rng_handle), KIND=I4B)
       ENDDO

       CALL READ_WGT_FITS(weights_p, files, files%in_weights, par)

       CALL READ_WGT_FITS(weights_mu, files, files%in_weights_mu, par)

       CALL READ_COV_FITS(cov_p, files, files%in_covariance)

       CALL READ_COV_FITS(cov_mu, files, files%in_covariance_mu)

    ENDIF

    CALL MPI_BCAST(rng%seeds__map, rng%nr_sims, MPI_INTEGER,&
         & 0, MPI_COMM_WORLD, par%mpi_error)
    CALL CHECK_ERROR_FLAG("MPI_BCAST", par%mpi_error)

    CALL MPI_BCAST(weights_p%m__r_TE_l, SIZE(weights_p%m__r_TE_l),&
         & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, par%mpi_error)
    CALL CHECK_ERROR_FLAG("MPI_BCAST", par%mpi_error)

    CALL MPI_BCAST(weights_mu%m__r_TE_l, SIZE(weights_mu%m__r_TE_l),&
         & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, par%mpi_error)
    CALL CHECK_ERROR_FLAG("MPI_BCAST", par%mpi_error)

    CALL MPI_BCAST(cov_p%sqrt__r1_r2_l, SIZE(cov_p%sqrt__r1_r2_l),&
         & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, par%mpi_error)
    CALL CHECK_ERROR_FLAG("MPI_BCAST", par%mpi_error)
    CALL MPI_BCAST(cov_p%var__r, cov_p%nr_shells,&
         & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, par%mpi_error)
    CALL CHECK_ERROR_FLAG("MPI_BCAST", par%mpi_error)

    CALL MPI_BCAST(cov_mu%sqrt__r1_r2_l, SIZE(cov_mu%sqrt__r1_r2_l),&
         & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, par%mpi_error)
    CALL CHECK_ERROR_FLAG("MPI_BCAST", par%mpi_error)
    CALL MPI_BCAST(cov_mu%var__r, cov_mu%nr_shells,&
         & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, par%mpi_error)
    CALL CHECK_ERROR_FLAG("MPI_BCAST", par%mpi_error)


  END SUBROUTINE READ_AUX_FILES

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Read weights, FITS format
!-----------------------------------------------------------------------

  SUBROUTINE READ_WGT_FITS(wgt, files, basename, par)

    CHARACTER(LEN=*),                       INTENT(IN)    :: basename
    TYPE(FILENAMES),                        INTENT(IN)    :: files
    TYPE(PARAMETERS),                       INTENT(INOUT) :: par
    TYPE(WEIGHTS),                          INTENT(INOUT) :: wgt

    LOGICAL(LGT)                                          :: is_phi
    LOGICAL(LGT)                                          :: anynull
    CHARACTER(LEN=2)                                      :: nr_shells_string
    CHARACTER(LEN=FNL)                                    :: filename
    CHARACTER(LEN=80)                                     :: comment
    CHARACTER(LEN=80)                                     :: tmp_hdr
    CHARACTER(LEN=80)                                     :: cosmo_file
    INTEGER(I4B)                                          :: status
    INTEGER(I4B)                                          :: unit
    INTEGER(I4B)                                          :: blocksize
    INTEGER(I4B)                                          :: nr_elements
    INTEGER(I4B)                                          :: lmax_file
    INTEGER(I4B)                                          :: nr_shells_file
    REAL(DP),          DIMENSION(:,:),      ALLOCATABLE   :: tmp

    WRITE(nr_shells_string,'(I2.2)') wgt%nr_shells

    status   = 0
    unit     = 0
    filename = TRIM(files%prefix_ng_aux) // '/'&
         & // TRIM(files%prefix_cosmo) // '/' // TRIM(basename)&
         & // nr_shells_string // '.fits'

    CALL FTGIOU(unit, status)
    CALL FTOPEN(unit, filename, 0, blocksize, status)
    CALL FTGKYS(unit, 'COSMO', tmp_hdr, cosmo_file, status)
    CALL FTGKYL(unit, 'PHI', is_phi, comment, status)
    CALL FTGKYJ(unit, 'MAX-LPOL', lmax_file, comment, status)
    CALL FTGKYJ(unit, 'NR_R', nr_shells_file, comment, status)

    IF (TRIM(cosmo_file) .NE. TRIM(files%prefix_cosmo)) THEN
       WRITE(*,'(/3X, "WARNING: Possible mismatch in cosmology")')
       WRITE(*,'(3X, "Expected input: ", A)') TRIM(files%prefix_cosmo)
       WRITE(*,'(3X, "Found in file:  ", A)') TRIM(cosmo_file)
    ENDIF

    IF (.NOT. is_phi) THEN
       CALL CHECK_ERROR_FLAG('read weights: Input for&
            & potential required', 1)
    ENDIF

    IF (lmax_file .LT. wgt%lmax) THEN
       CALL CHECK_ERROR_FLAG('read weights: Input lmax incompatible', 1)
    ENDIF

    IF (nr_shells_file .NE. wgt%nr_shells) THEN
       CALL CHECK_ERROR_FLAG('read weights: Input shell&
            & number incompatible', 1)
    ENDIF

    CALL FTMNHD(unit, -1, 'R', 0, status)
    CALL FTGPVD(unit, 1, 1, nr_shells_file, 0, par%radii__r, anynull,&
         & status)

    ALLOCATE(tmp(0:lmax_file,1:nr_shells_file))
    tmp         = 0.0_DP
    nr_elements = (lmax_file+1)*nr_shells_file

    CALL FTMNHD(unit, -1, 'WEIGHTS_T', 0, status)
    CALL FTGPVD(unit, 1, 1, nr_elements, 0, tmp, anynull, status)
    wgt%m__r_TE_l(:,1,2:wgt%lmax) = TRANSPOSE(tmp(2:wgt%lmax,:))

    IF (wgt%pol) THEN

       CALL FTMNHD(unit, -1, 'WEIGHTS_E', 0, status)
       CALL FTGPVD(unit, 1, 1, nr_elements, 0, tmp, anynull, status)
       wgt%m__r_TE_l(:,2,2:wgt%lmax) = TRANSPOSE(tmp(2:wgt%lmax,:))

    ENDIF

    CALL FTCLOS(unit, status)
    CALL FTFIOU(unit, status)

    DEALLOCATE(tmp)


  END SUBROUTINE READ_WGT_FITS

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Read covariance, FITS format
!-----------------------------------------------------------------------

  SUBROUTINE READ_COV_FITS(cov, files, basename)

    CHARACTER(LEN=*),                       INTENT(IN)    :: basename
    TYPE(FILENAMES),                        INTENT(IN)    :: files
    TYPE(COVARIANCE),                       INTENT(INOUT) :: cov

    LOGICAL(LGT)                                          :: is_phi
    LOGICAL(LGT)                                          :: anynull
    CHARACTER(LEN=2)                                      :: nr_shells_string
    CHARACTER(LEN=FNL)                                    :: filename
    CHARACTER(LEN=80)                                     :: comment
    CHARACTER(LEN=80)                                     :: tmp_hdr
    CHARACTER(LEN=80)                                     :: cosmo_file
    INTEGER(I4B)                                          :: i
    INTEGER(I4B)                                          :: l
    INTEGER(I4B)                                          :: status
    INTEGER(I4B)                                          :: unit
    INTEGER(I4B)                                          :: blocksize
    INTEGER(I4B)                                          :: nr_elements
    INTEGER(I4B)                                          :: lmax_file
    INTEGER(I4B)                                          :: nr_shells_file
    REAL(DP),          DIMENSION(:,:,:),    ALLOCATABLE   :: tmp

    WRITE(nr_shells_string,'(I2.2)') cov%nr_shells

    status   = 0
    unit     = 0
    filename = TRIM(files%prefix_ng_aux) // '/'&
         & // TRIM(files%prefix_cosmo) // '/'&
         & // TRIM(basename) // nr_shells_string // '.fits'

    CALL FTGIOU(unit, status)
    CALL FTOPEN(unit, filename, 0, blocksize, status)
    CALL FTGKYS(unit, 'COSMO', tmp_hdr, cosmo_file, status)
    CALL FTGKYL(unit, 'PHI', is_phi, comment, status)
    CALL FTGKYJ(unit, 'MAX-LPOL', lmax_file, comment, status)
    CALL FTGKYJ(unit, 'NR_R', nr_shells_file, comment, status)

    IF (TRIM(cosmo_file) .NE. TRIM(files%prefix_cosmo)) THEN
       WRITE(*,'(/3X, "WARNING: Possible mismatch in cosmology")')
       WRITE(*,'(3X, "Expected input: ", A)') TRIM(files%prefix_cosmo)
       WRITE(*,'(3X, "Found in file:  ", A)') TRIM(cosmo_file)
    ENDIF

    IF (.NOT. is_phi) THEN
       CALL CHECK_ERROR_FLAG('read covariance:&
            & Input for potential required', 1)
    ENDIF

    IF (lmax_file .LT. cov%lmax) THEN
       CALL CHECK_ERROR_FLAG('read covariance:&
            & Input lmax incompatible', 1)
    ENDIF

    IF (nr_shells_file .NE. cov%nr_shells) THEN
       CALL CHECK_ERROR_FLAG('read covariance:&
            & Input shell number incompatible', 1)
    ENDIF

    ALLOCATE(tmp(0:lmax_file,1:nr_shells_file,1:nr_shells_file))
    tmp         = 0.0_DP
    nr_elements = (lmax_file+1) * nr_shells_file**2
    CALL FTMNHD(unit, -1, 'COV_P', 0, status)
    CALL FTGPVD(unit, 1, 1, nr_elements, 0, tmp, anynull, status)
    DO l = 2,cov%lmax
       cov%sqrt__r1_r2_l(:,:,l) = tmp(l,:,:)
    ENDDO

    CALL FTCLOS(unit, status)
    CALL FTFIOU(unit, status)

    DO i = 1,cov%nr_shells
       DO l = 2,cov%lmax
          cov%var__r(i) = cov%var__r(i) + REAL(2*l+1, KIND=DP)&
               & / FOURPI * cov%sqrt__r1_r2_l(i,i,l)
       ENDDO
    ENDDO

    DEALLOCATE(tmp)


  END SUBROUTINE READ_COV_FITS

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Create sampling matrices
!-----------------------------------------------------------------------

  SUBROUTINE COMPUTE_MATRICES(cov, par)

    USE LAPACK95,      ONLY: SYEVR

    TYPE(PARAMETERS),                       INTENT(INOUT) :: par
    TYPE(COVARIANCE),                       INTENT(INOUT) :: cov

    INTEGER(I4B)                                          :: i
    INTEGER(I4B)                                          :: j
    INTEGER(I4B)                                          :: error
    REAL(DP),                               EXTERNAL      :: DLAMCH
    REAL(DP),          DIMENSION(:),        ALLOCATABLE   :: eigenvalues__k
    REAL(DP),          DIMENSION(:,:),      ALLOCATABLE   :: eigenvectors__r_k
    REAL(DP),          DIMENSION(:,:),      ALLOCATABLE   :: diagonal_matrix

    ALLOCATE(eigenvectors__r_k(1:cov%nr_shells,1:cov%nr_shells))
    ALLOCATE(diagonal_matrix(1:cov%nr_shells,1:cov%nr_shells))
    ALLOCATE(eigenvalues__k(1:cov%nr_shells))
    eigenvectors__r_k = 0.0_DP
    diagonal_matrix   = 0.0_DP
    eigenvalues__k    = 0.0_DP

    CALL GET_LOOP_BOUNDS(par, 2, cov%lmax, par%thread_id)

    DO i = par%loop_min,par%loop_max

       diagonal_matrix = cov%sqrt__r1_r2_l(:,:,i)

       CALL SYEVR(diagonal_matrix, eigenvalues__k, uplo='U',&
            & z=eigenvectors__r_k, il=1, iu=cov%nr_shells,&
            & abstol=DLAMCH('S'), info=error)

       CALL CHECK_ERROR_FLAG("Lapack", error)

       WHERE(eigenvalues__k .LT. 0.0_DP) eigenvalues__k = 0.0_DP

       diagonal_matrix = 0.0_DP
       DO j = 1,cov%nr_shells
          diagonal_matrix(j,j) = DSQRT(eigenvalues__k(j))
       ENDDO

       cov%sqrt__r1_r2_l(:,:,i) = MATMUL(eigenvectors__r_k,&
            & MATMUL(diagonal_matrix, TRANSPOSE(eigenvectors__r_k)))

    ENDDO

    CALL MPI_BARRIER(MPI_COMM_WORLD, par%mpi_error)

    DO i = 0,par%nr_nodes-1

       CALL GET_LOOP_BOUNDS(par, 2, par%lmax, i)

       CALL MPI_BCAST(cov%sqrt__r1_r2_l(:,:,par%loop_min:par%loop_max),&
            & SIZE(cov%sqrt__r1_r2_l(:,:,par%loop_min:par%loop_max)),&
            & MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, par%mpi_error)

    ENDDO

    DEALLOCATE(eigenvalues__k)
    DEALLOCATE(diagonal_matrix)
    DEALLOCATE(eigenvectors__r_k)


  END SUBROUTINE COMPUTE_MATRICES

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Compute loop boundaries for each thread
!-----------------------------------------------------------------------

  SUBROUTINE GET_LOOP_BOUNDS(par, loop_start, loop_stop, thread_id)

    INTEGER(I4B),                           INTENT(IN)    :: loop_start
    INTEGER(I4B),                           INTENT(IN)    :: loop_stop
    INTEGER(I4B),                           INTENT(IN)    :: thread_id
    TYPE(PARAMETERS),                       INTENT(INOUT) :: par

    INTEGER(I4B)                                          :: nr_iterations
    REAL(DP)                                              :: index_fraction
    REAL(DP)                                              :: index_min
    REAL(DP)                                              :: index_max

    nr_iterations = loop_stop - loop_start + 1

    IF ((thread_id .LT. 0) .OR. (par%nr_nodes .LT. thread_id - 1)&
         & .OR. (nr_iterations .LT. par%nr_nodes)) THEN
       WRITE(*,'(3X, "Error: Loop indices")')
       STOP
    ENDIF

    index_fraction = REAL(nr_iterations, KIND=DP)&
         & / REAL(par%nr_nodes, KIND=DP)

    index_min = index_fraction * REAL(thread_id, KIND=DP)
    index_max = index_fraction * REAL(thread_id + 1, KIND=DP)

    par%loop_min = loop_start + INT(index_min, KIND=I4B)
    par%loop_max = loop_start + INT(index_max, KIND=I4B) - 1

    IF (thread_id .EQ. par%nr_nodes - 1) THEN
       par%loop_max = loop_stop
    ENDIF

    IF (par%loop_max .LT. par%loop_min) THEN
       WRITE(*,'(3X, "Error: Loop indices")')
       STOP
    ENDIF


  END SUBROUTINE GET_LOOP_BOUNDS

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Compute seeds for random number generator
!-----------------------------------------------------------------------

  SUBROUTINE GENERATE_RND_SEEDS(rng, sim_nr, par)

    USE RNGMOD,        ONLY: RAND_UNI, RAND_INIT, PLANCK_RNG

    TYPE(PARAMETERS),                       INTENT(IN)    :: par
    INTEGER(I4B),                           INTENT(IN)    :: sim_nr
    TYPE(RND_SEEDS),                        INTENT(INOUT) :: rng

    INTEGER(I4B)                                          :: i
    INTEGER(I4B)                                          :: l
    INTEGER(I4B)                                          :: m
    INTEGER(I4B)                                          :: start_seed
    TYPE(PLANCK_RNG)                                      :: rng_handle

    start_seed = rng%seeds__map(sim_nr - par%sim_nr_start + 1)

    CALL RAND_INIT(rng_handle, start_seed)

    DO m = 0,rng%mmax
       DO l = MAX(2,m),rng%lmax
          DO i = 1,4
             rng%seeds__nr_l_m(i,l,m)&
                  & = INT(MAX_I4B * RAND_UNI(rng_handle), KIND=I4B)
          ENDDO
       ENDDO
    ENDDO


  END SUBROUTINE GENERATE_RND_SEEDS

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Generate univariate alm vector
!-----------------------------------------------------------------------

  SUBROUTINE GENERATE_RND_VECTOR(gauss_rand__r, rng, l, m)

    USE RNGMOD,        ONLY: RAND_GAUSS, RAND_INIT, PLANCK_RNG

    TYPE(RND_SEEDS),                        INTENT(IN)    :: rng
    INTEGER(I4B),                           INTENT(IN)    :: l
    INTEGER(I4B),                           INTENT(IN)    :: m
    COMPLEX(DPC),      DIMENSION(1:),       INTENT(INOUT) :: gauss_rand__r

    INTEGER(I4B)                                          :: i
    TYPE(PLANCK_RNG)                                      :: rng_handle

    CALL RAND_INIT(rng_handle, rng%seeds__nr_l_m(1,l,m),&
         & rng%seeds__nr_l_m(2,l,m), rng%seeds__nr_l_m(3,l,m),&
         & rng%seeds__nr_l_m(4,l,m))

    IF (m .EQ. 0) THEN
       DO i = 1,rng%nr_shells
          gauss_rand__r(i) = DCMPLX(rand_gauss(rng_handle), 0.0_DP)
       ENDDO
    ELSE
       DO i = 1,rng%nr_shells
          gauss_rand__r(i) = DCMPLX(rand_gauss(rng_handle)/SQRT2,&
               & rand_gauss(rng_handle)/SQRT2)
       ENDDO
    ENDIF


  END SUBROUTINE GENERATE_RND_VECTOR

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Convert map to alms
!-----------------------------------------------------------------------

  SUBROUTINE MAP_2_ALM(map)

    USE ALM_TOOLS,     ONLY: MAP2ALM

    CLASS(SKYMAP),                          INTENT(INOUT) :: map

    IF (map%pol) THEN
       CALL MAP2ALM(map%nside, map%lmax, map%mmax,&
            & map%map__pix_IQU, map%a__TEB_l_m,&
            & map%sht%zbounds__lat, map%sht%w8ring__ring_IQU)
    ELSE
       CALL MAP2ALM(map%nside, map%lmax, map%mmax,&
            & map%map__pix_IQU(:,1), map%a__TEB_l_m,&
            & map%sht%zbounds__lat, map%sht%w8ring__ring_IQU)
    ENDIF


  END SUBROUTINE MAP_2_ALM

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Convert alms to map
!-----------------------------------------------------------------------

  SUBROUTINE ALM_2_MAP(map)

    USE ALM_TOOLS,     ONLY: ALM2MAP

    CLASS(SKYMAP),                          INTENT(INOUT) :: map

    IF (map%pol) THEN
       CALL ALM2MAP(map%nside, map%lmax, map%mmax, map%a__TEB_l_m,&
            & map%map__pix_IQU)
    ELSE
       CALL ALM2MAP(map%nside, map%lmax, map%mmax, map%a__TEB_l_m,&
            & map%map__pix_IQU(:,1))
    ENDIF


  END SUBROUTINE ALM_2_MAP

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Get next free file number
!-----------------------------------------------------------------------

  SUBROUTINE GET_FILE_NUMBER(par, files)

    TYPE(FILENAMES),                        INTENT(IN)    :: files
    TYPE(PARAMETERS),                       INTENT(INOUT) :: par

    CHARACTER(LEN=FNL)                                    :: filename
    LOGICAL(LGT)                                          :: is_present
    INTEGER(I4B)                                          :: i

    IF (par%thread_id .EQ. par%root) THEN

       i          = 0
       is_present = .TRUE.

       DO WHILE(is_present)

          i = i + 1

          filename = GENERATE_FILENAME(files%out_signal_l, i)

          INQUIRE(FILE=filename, EXIST=is_present)

          IF (i .GE. 10000) THEN
             WRITE(*,'(3X, "Error in filename")')
             STOP
          ENDIF

       ENDDO

       par%sim_nr_start = i

    ENDIF

    CALL MPI_BCAST(par%sim_nr_start, 1, MPI_INTEGER, 0,&
         & MPI_COMM_WORLD, par%mpi_error)
    CALL CHECK_ERROR_FLAG("MPI_BCAST", par%mpi_error)


  END SUBROUTINE GET_FILE_NUMBER

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Write alm to file
!-----------------------------------------------------------------------

  SUBROUTINE OUTPUT_ALM(map, rng, sim_nr, output_file)

    USE MISC_UTILS
    USE HEAD_FITS,     ONLY: WRITE_MINIMAL_HEADER, ADD_CARD
    USE FITSTOOLS,     ONLY: ALMS2FITS

    CHARACTER(LEN=*),                       INTENT(IN)    :: output_file
    TYPE(RND_SEEDS),                        INTENT(IN)    :: rng
    INTEGER(I4B),                           INTENT(IN)    :: sim_nr
    TYPE(SKYMAP),                           INTENT(INOUT) :: map

    CHARACTER(LEN=FNL)                                    :: filename
    INTEGER(I4B)                                          :: nalms
    INTEGER(I4B)                                          :: i
    INTEGER(I4B)                                          :: j
    INTEGER(I4B)                                          :: index
    REAL(DP),          DIMENSION(:,:,:),    ALLOCATABLE   :: dummy_alm

    nalms = INT((map%lmax + 1)*(map%lmax + 2)/2, KIND=I4B)

    IF (map%pol) THEN
       ALLOCATE(dummy_alm(1:nalms,1:4,1:3))
    ELSE
       ALLOCATE(dummy_alm(1:nalms,1:4,1:1))
    ENDIF
    dummy_alm = 0.0_DP

    DO i = 0,map%lmax
       DO j = 0,i
          index = INT(i*(i+1)/2)+j+1
          dummy_alm(index,1,:) = REAL(i, KIND=DP)
          dummy_alm(index,2,:) = REAL(j, KIND=DP)
          IF ((i .GE. 2) .AND. (j .LE. map%mmax)) THEN
             dummy_alm(index,3,1) = REAL(map%a__TEB_l_m(1,i,j), KIND=DP)
             dummy_alm(index,4,1) = DIMAG(map%a__TEB_l_m(1,i,j))
             IF (map%pol) THEN
                dummy_alm(index,3,2) = REAL(map%a__TEB_l_m(2,i,j), KIND=DP)
                dummy_alm(index,4,2) = DIMAG(map%a__TEB_l_m(2,i,j))
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    CALL WRITE_MINIMAL_HEADER(map%header(:,1), 'ALM', append=.FALSE.,&
         & nside=map%nside, coordsys='G',&
         & randseed=rng%seeds__map(sim_nr), units='dimensionless',&
         & nlmax=map%lmax, polar=map%pol, nmmax=map%mmax)
    CALL ADD_CARD(map%header(:,1), "COMMENT",&
         & "---------------------------------------------")
    CALL ADD_CARD(map%header(:,1), "COMMENT",&
         & "           NON-GAUSSIAN SIMULATION           ")
    CALL ADD_CARD(map%header(:,1), "COMMENT",&
         & "---------------------------------------------")

    filename = GENERATE_FILENAME(output_file, sim_nr)

    CALL ASSERT_NOT_PRESENT(filename)
    IF (map%pol) THEN
       CALL ALMS2FITS(filename, nalms, dummy_alm, 3, map%header, 80, 3)
    ELSE
       CALL ALMS2FITS(filename, nalms, dummy_alm, 3, map%header, 80, 1)
    ENDIF

    DEALLOCATE(dummy_alm)


  END SUBROUTINE OUTPUT_ALM

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Write map to file
!-----------------------------------------------------------------------

  SUBROUTINE OUTPUT_MAP(map, rng, sim_nr, output_file)

    USE MISC_UTILS,    ONLY: ASSERT_NOT_PRESENT
    USE HEAD_FITS,     ONLY: WRITE_MINIMAL_HEADER, ADD_CARD
    USE FITSTOOLS,     ONLY: WRITE_BINTAB

    CHARACTER(LEN=*),                       INTENT(IN)    :: output_file
    TYPE(RND_SEEDS),                        INTENT(IN)    :: rng
    INTEGER(I4B),                           INTENT(IN)    :: sim_nr
    TYPE(SKYMAP),                           INTENT(INOUT) :: map

    CHARACTER(LEN=FNL)                                    :: filename

    CALL WRITE_MINIMAL_HEADER(map%header(:,1), 'MAP', append=.FALSE.,&
         & nside=map%nside, ordering='RING', coordsys='G',&
         & randseed=rng%seeds__map(sim_nr), units='dimensionless',&
         & nlmax=map%lmax, polar=map%pol, nmmax=map%mmax)
    CALL ADD_CARD(map%header(:,1), "COMMENT",&
         & "---------------------------------------------")
    CALL ADD_CARD(map%header(:,1), "COMMENT",&
         & "           NON-GAUSSIAN SIMULATION           ")
    CALL ADD_CARD(map%header(:,1), "COMMENT",&
         & "---------------------------------------------")

    filename = GENERATE_FILENAME(output_file, sim_nr)

    CALL ASSERT_NOT_PRESENT(filename)
    IF (map%pol) THEN
       CALL WRITE_BINTAB(map%map__pix_IQU, map%npix, 3,&
            & map%header(:,1), 80, filename)
    ELSE
       CALL WRITE_BINTAB(map%map__pix_IQU, map%npix, 1,&
            & map%header(:,1), 80, filename)
    ENDIF


  END SUBROUTINE OUTPUT_MAP

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Convert alm to map, write to file
!-----------------------------------------------------------------------

  SUBROUTINE OUTPUT_AS_MAP(map, rng, sim_nr, output_file)

    CHARACTER(LEN=*),                       INTENT(IN)    :: output_file
    TYPE(RND_SEEDS),                        INTENT(IN)    :: rng
    INTEGER(I4B),                           INTENT(IN)    :: sim_nr
    TYPE(SKYMAP),                           INTENT(INOUT) :: map

    CALL ALM_2_MAP(map)

    CALL OUTPUT_MAP(map, rng, sim_nr, output_file)


  END SUBROUTINE OUTPUT_AS_MAP

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Write matrix to file
!-----------------------------------------------------------------------

  SUBROUTINE OUTPUT_MATRIX(filename, matrix)

    CHARACTER(LEN=*),                       INTENT(IN)    :: filename
    REAL(DP),          DIMENSION(1:,1:),    INTENT(IN)    :: matrix

    INTEGER(I4B)                                          :: i

    OPEN(UNIT=17, FILE=filename, STATUS='UNKNOWN',&
         & ACTION='WRITE', FORM='FORMATTED')
    REWIND 17

    DO i = 1,SIZE(matrix, 2)
       WRITE(17,'(10000(1X, ES23.15E3, :))') matrix(:,i)
    ENDDO

    CLOSE(17)


  END SUBROUTINE OUTPUT_MATRIX

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Construct filenames
!-----------------------------------------------------------------------

  FUNCTION GENERATE_FILENAME(file_inout, map_nr)

    USE MISC_UTILS,    ONLY: STRING

    CHARACTER(LEN=*),                       INTENT(IN)    :: file_inout
    INTEGER(I4B),                           INTENT(IN)    :: map_nr
    CHARACTER(LEN=FNL)                                    :: GENERATE_FILENAME

    CHARACTER(LEN=4)                                      :: map_nr_string

    map_nr_string = STRING(map_nr, '(I4.4)')

    GENERATE_FILENAME = REPLACE_SUBSTRING(file_inout, '%MCNUMBER%',&
         & map_nr_string)


  END FUNCTION GENERATE_FILENAME

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Replace string placeholder with variable value
!-----------------------------------------------------------------------

  FUNCTION REPLACE_SUBSTRING(string_in, string_replace, string_insert)

    CHARACTER(LEN=*),                       INTENT(IN)    :: string_in
    CHARACTER(LEN=*),                       INTENT(IN)    :: string_replace
    CHARACTER(LEN=*),                       INTENT(IN)    :: string_insert
    CHARACTER(LEN=FNL)                                    :: REPLACE_SUBSTRING

    CHARACTER(LEN=FNL)                                    :: tmp_in
    CHARACTER(LEN=FNL)                                    :: tmp_replace
    CHARACTER(LEN=FNL)                                    :: tmp_insert
    INTEGER(I4B)                                          :: i
    INTEGER(I4B)                                          :: len_in
    INTEGER(I4B)                                          :: len_replace
    INTEGER(I4B)                                          :: len_insert

    tmp_in      = ADJUSTL(string_in)
    tmp_replace = ADJUSTL(string_replace)
    tmp_insert  = ADJUSTL(string_insert)

    len_in      = LEN(TRIM(tmp_in))
    len_replace = LEN(TRIM(tmp_replace))
    len_insert  = LEN(TRIM(tmp_insert))

    DO i = 1,len_in-len_replace+1
        IF (tmp_in(i:i+len_replace-1) .EQ.&
             & tmp_replace(1:len_replace)) THEN
          REPLACE_SUBSTRING = tmp_in(1:i-1)&
               & // TRIM(tmp_insert(1:len_insert))&
               & // tmp_in(i+len_replace:len_in)
          RETURN
       ENDIF
    ENDDO

    REPLACE_SUBSTRING = string_in


  END FUNCTION REPLACE_SUBSTRING

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Clock program
!-----------------------------------------------------------------------

  SUBROUTINE CLOCK()

    USE MISC_UTILS,    ONLY: WALL_CLOCK_TIME

    LOGICAL(LGT),                           SAVE          :: firstcall = .TRUE.
    REAL(SP),                               SAVE          :: wct_start
    REAL(SP),                               SAVE          :: wct_stop

    IF (firstcall) THEN
       CALL WALL_CLOCK_TIME(wct_start)
       firstcall= .FALSE.
    ELSE
       CALL WALL_CLOCK_TIME(wct_stop)
       WRITE(*,'(3X, "Time: ", F9.2, " min")')&
            & (wct_stop - wct_start) / 60.0_DP
    ENDIF


  END SUBROUTINE CLOCK

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Check error flag
!-----------------------------------------------------------------------

  SUBROUTINE CHECK_ERROR_FLAG(code, error)

    CHARACTER(LEN=*),                       INTENT(IN)    :: code
    INTEGER(I4B),                           INTENT(IN)    :: error

    IF (error .NE. 0) THEN
       WRITE(*,'(/3X, "Error detected by ", A"; errorcode ", I3)')&
            & code, error
       STOP
    ENDIF


  END SUBROUTINE CHECK_ERROR_FLAG

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Cleanup
!-----------------------------------------------------------------------

  SUBROUTINE DEALLOCATE_ARRAYS(weights_p, weights_mu, cov_p, cov_mu,&
       & rng, signal_l, signal_nl, signal_nnl, signal_mu, par)

    TYPE(PARAMETERS),                       INTENT(INOUT) :: par
    TYPE(RND_SEEDS),                        INTENT(INOUT) :: rng
    TYPE(WEIGHTS),                          INTENT(INOUT) :: weights_p
    TYPE(WEIGHTS),                          INTENT(INOUT) :: weights_mu
    TYPE(COVARIANCE),                       INTENT(INOUT) :: cov_p
    TYPE(COVARIANCE),                       INTENT(INOUT) :: cov_mu
    TYPE(SKYMAP),                           INTENT(INOUT) :: signal_l
    TYPE(SKYMAP),                           INTENT(INOUT) :: signal_nl
    TYPE(SKYMAP),                           INTENT(INOUT) :: signal_nnl
    TYPE(SKYMAP),                           INTENT(INOUT) :: signal_mu

    CALL par%FREE()
    CALL rng%FREE()
    CALL weights_p%FREE()
    CALL weights_mu%FREE()
    CALL cov_p%FREE()
    CALL cov_mu%FREE()
    CALL signal_l%FREE()
    CALL signal_nl%FREE()
    CALL signal_nnl%FREE()
    CALL signal_mu%FREE()


  END SUBROUTINE DEALLOCATE_ARRAYS

!-----------------------------------------------------------------------

END MODULE MAP_MOD_MISC
