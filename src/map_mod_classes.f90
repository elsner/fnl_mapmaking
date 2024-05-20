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


MODULE MAP_MOD_CLASSES

  USE HEALPIX_TYPES

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Global variables
!-----------------------------------------------------------------------

  INTEGER(I4B),                             PARAMETER     :: FNL = 200

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Compilation of filenames
!-----------------------------------------------------------------------

  TYPE FILENAMES

     CHARACTER(LEN=FNL)                                   :: in_hpx_data
     CHARACTER(LEN=FNL)                                   :: prefix_cosmo
     CHARACTER(LEN=FNL)                                   :: prefix_ng_aux
     CHARACTER(LEN=FNL)                                   :: in_weights
     CHARACTER(LEN=FNL)                                   :: in_weights_mu
     CHARACTER(LEN=FNL)                                   :: in_covariance
     CHARACTER(LEN=FNL)                                   :: in_covariance_mu
     CHARACTER(LEN=FNL)                                   :: out_signal_l
     CHARACTER(LEN=FNL)                                   :: out_signal_nl
     CHARACTER(LEN=FNL)                                   :: out_signal_nnl
     CHARACTER(LEN=FNL)                                   :: out_signal_mu


  END TYPE FILENAMES

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Seeds of random number generator
!-----------------------------------------------------------------------

  TYPE RND_SEEDS

     INTEGER(I4B)                                         :: nr_sims
     INTEGER(I4B)                                         :: nr_shells
     INTEGER(I4B)                                         :: lmax
     INTEGER(I4B)                                         :: mmax
     INTEGER(I4B),     DIMENSION(:),        ALLOCATABLE   :: seeds__map
     INTEGER(I4B),     DIMENSION(:,:,:),    ALLOCATABLE   :: seeds__nr_l_m


   CONTAINS

     PROCEDURE                                            :: INIT => INIT_RND
     PROCEDURE                                            :: FREE => FREE_RND


  END TYPE RND_SEEDS

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Integration weights
!-----------------------------------------------------------------------

  TYPE WEIGHTS

     LOGICAL(LGT)                                         :: pol
     INTEGER(I4B)                                         :: lmax
     INTEGER(I4B)                                         :: nr_shells
     REAL(DP),         DIMENSION(:,:,:),    ALLOCATABLE   :: m__r_TE_l


   CONTAINS

     PROCEDURE                                            :: INIT => INIT_WGT
     PROCEDURE                                            :: FREE => FREE_WGT


  END TYPE WEIGHTS

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Square root of covariance matrix
!-----------------------------------------------------------------------

  TYPE COVARIANCE

     INTEGER(I4B)                                         :: lmax
     INTEGER(I4B)                                         :: nr_shells
     REAL(DP),         DIMENSION(:),        ALLOCATABLE   :: var__r
     REAL(DP),         DIMENSION(:,:,:),    ALLOCATABLE   :: sqrt__r1_r2_l


   CONTAINS

     PROCEDURE                                            :: INIT => INIT_COV
     PROCEDURE                                            :: FREE => FREE_COV


  END TYPE COVARIANCE

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Auxiliary data for spherical harmonic transform
!-----------------------------------------------------------------------

  TYPE S_H_TRANSFORM

     REAL(DP),         DIMENSION(:),        ALLOCATABLE   :: zbounds__lat
     REAL(DP),         DIMENSION(:,:),      ALLOCATABLE   :: w8ring__ring_IQU


   CONTAINS

     PROCEDURE                                            :: INIT => INIT_SHT
     PROCEDURE                                            :: FREE => FREE_SHT


  END TYPE S_H_TRANSFORM

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Map in pixel and harmonic space T or T & E
!-----------------------------------------------------------------------

  TYPE SKYMAP

     TYPE(S_H_TRANSFORM)                                  :: sht
     LOGICAL(LGT)                                         :: pol
     CHARACTER(LEN=80), DIMENSION(:,:),     ALLOCATABLE   :: header
     INTEGER(I4B)                                         :: lmax
     INTEGER(I4B)                                         :: mmax
     INTEGER(I4B)                                         :: nside
     INTEGER(I4B)                                         :: npix
     REAL(DP),         DIMENSION(:,:),      ALLOCATABLE   :: map__pix_IQU
     COMPLEX(DPC),     DIMENSION(:,:,:),    ALLOCATABLE   :: a__TEB_l_m


   CONTAINS

     PROCEDURE                                            :: INIT => INIT_MAP
     PROCEDURE                                            :: FREE => FREE_MAP


  END TYPE SKYMAP

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Collection of basic simulation parameters
!-----------------------------------------------------------------------

  TYPE PARAMETERS

     TYPE(S_H_TRANSFORM)                                  :: sht
     LOGICAL(LGT)                                         :: pol
     INTEGER(I4B)                                         :: nr_sims
     INTEGER(I4B)                                         :: sim_nr_start
     INTEGER(I4B)                                         :: nside
     INTEGER(I4B)                                         :: npix
     INTEGER(I4B)                                         :: lmax
     INTEGER(I4B)                                         :: mmax
     INTEGER(I4B)                                         :: nr_shells
     INTEGER(I4B)                                         :: mpi_error
     INTEGER(I4B)                                         :: root
     INTEGER(I4B)                                         :: thread_id
     INTEGER(I4B)                                         :: nr_nodes
     INTEGER(I4B)                                         :: loop_min
     INTEGER(I4B)                                         :: loop_max
     REAL(DP),         DIMENSION(:),        ALLOCATABLE   :: radii__r


   CONTAINS

     PROCEDURE                                            :: INIT => INIT_PAR
     PROCEDURE                                            :: FREE => FREE_PAR


  END TYPE PARAMETERS

!-----------------------------------------------------------------------


CONTAINS

!-----------------------------------------------------------------------
! Constructor routine: random number seeds
!-----------------------------------------------------------------------

  SUBROUTINE INIT_RND(rnd, par)

    TYPE(PARAMETERS),                       INTENT(IN)    :: par
    CLASS(RND_SEEDS),                       INTENT(INOUT) :: rnd

    rnd%nr_sims   = par%nr_sims
    rnd%nr_shells = par%nr_shells
    rnd%lmax      = par%lmax
    rnd%mmax      = par%mmax

    ALLOCATE(rnd%seeds__map(1:rnd%nr_sims))
    ALLOCATE(rnd%seeds__nr_l_m(1:4,0:rnd%lmax,0:rnd%mmax))
    rnd%seeds__map    = 0
    rnd%seeds__nr_l_m = 0


  END SUBROUTINE INIT_RND

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Destructor routine: random number seeds
!-----------------------------------------------------------------------

  SUBROUTINE FREE_RND(rnd)

    CLASS(RND_SEEDS),                       INTENT(INOUT) :: rnd

    DEALLOCATE(rnd%seeds__map)
    DEALLOCATE(rnd%seeds__nr_l_m)


  END SUBROUTINE FREE_RND

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Constructor routine: weights
!-----------------------------------------------------------------------

  SUBROUTINE INIT_WGT(wgt, par)

    TYPE(PARAMETERS),                       INTENT(IN)    :: par
    CLASS(WEIGHTS),                         INTENT(INOUT) :: wgt

    wgt%pol       = par%pol
    wgt%lmax      = par%lmax
    wgt%nr_shells = par%nr_shells

    IF (wgt%pol) THEN
       ALLOCATE(wgt%m__r_TE_l(1:wgt%nr_shells,1:2,0:wgt%lmax))
    ELSE
       ALLOCATE(wgt%m__r_TE_l(1:wgt%nr_shells,1:1,0:wgt%lmax))
    ENDIF
    wgt%m__r_TE_l = 0.0_DP


  END SUBROUTINE INIT_WGT

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Destructor routine: random number seeds
!-----------------------------------------------------------------------

  SUBROUTINE FREE_WGT(wgt)

    CLASS(WEIGHTS),                         INTENT(INOUT) :: wgt

    DEALLOCATE(wgt%m__r_TE_l)


  END SUBROUTINE FREE_WGT

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Constructor routine: covariance
!-----------------------------------------------------------------------

  SUBROUTINE INIT_COV(cov, par)

    TYPE(PARAMETERS),                       INTENT(IN)    :: par
    CLASS(COVARIANCE),                      INTENT(INOUT) :: cov

    cov%lmax      = par%lmax
    cov%nr_shells = par%nr_shells

    ALLOCATE(cov%var__r(1:cov%nr_shells))
    ALLOCATE(cov%sqrt__r1_r2_l(1:cov%nr_shells,1:cov%nr_shells,0:cov%lmax))
    cov%var__r        = 0.0_DP
    cov%sqrt__r1_r2_l = 0.0_DP


  END SUBROUTINE INIT_COV

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Destructor routine: covariance
!-----------------------------------------------------------------------

  SUBROUTINE FREE_COV(cov)

    CLASS(COVARIANCE),                       INTENT(INOUT) :: cov

    DEALLOCATE(cov%sqrt__r1_r2_l)


  END SUBROUTINE FREE_COV

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Constructor routine: SHT
!-----------------------------------------------------------------------

  SUBROUTINE INIT_SHT(sht, par)

    TYPE(PARAMETERS),                       INTENT(IN)    :: par
    CLASS(S_H_TRANSFORM),                   INTENT(INOUT) :: sht

    IF (par%nside .LE. 0 ) THEN
       WRITE(*,'(3X, "Error: Initializing SHT")')
       STOP
    ENDIF

    ALLOCATE(sht%zbounds__lat(1:2))
    IF (par%pol) THEN
       ALLOCATE(sht%w8ring__ring_IQU(1:2*par%nside,1:1))
    ELSE
       ALLOCATE(sht%w8ring__ring_IQU(1:2*par%nside,1:3))
    ENDIF

    sht%zbounds__lat     = (/ -1.0_DP, 1.0_DP /)
    sht%w8ring__ring_IQU = 0.0_DP


  END SUBROUTINE INIT_SHT

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Destructor routine: SHT
!-----------------------------------------------------------------------

  SUBROUTINE FREE_SHT(sht)

    CLASS(S_H_TRANSFORM),                   INTENT(INOUT) :: sht

    DEALLOCATE(sht%zbounds__lat)
    DEALLOCATE(sht%w8ring__ring_IQU)


  END SUBROUTINE FREE_SHT

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Constructor routine: Map
!-----------------------------------------------------------------------

  SUBROUTINE INIT_MAP(map, par, polarization)

    TYPE(PARAMETERS),                       INTENT(IN)    :: par
    LOGICAL(LGT),      OPTIONAL,            INTENT(IN)    :: polarization
    CLASS(SKYMAP),                          INTENT(INOUT) :: map

    IF (PRESENT(polarization)) THEN
       map%pol = polarization
    ELSE
       map%pol = par%pol
    ENDIF

    CALL map%sht%INIT(par)

    map%lmax                 = par%lmax
    map%mmax                 = par%mmax
    map%npix                 = par%npix
    map%nside                = par%nside
    map%sht%w8ring__ring_IQU = par%sht%w8ring__ring_IQU

    IF (.NOT. ((map%npix .GT. 0) .AND. (map%lmax .GT. 0) .AND.&
         & (map%mmax .GT. 0)) ) THEN
       WRITE(*,'(3X, "Error: Initializing map")')
       STOP
    ENDIF

    IF (map%pol) THEN
       ALLOCATE(map%header(1:80,1:3))
       ALLOCATE(map%map__pix_IQU(0:map%npix-1,1:3))
       ALLOCATE(map%a__TEB_l_m(1:3,0:map%lmax,0:map%mmax))
    ELSE
       ALLOCATE(map%header(1:80,1:1))
       ALLOCATE(map%map__pix_IQU(0:map%npix-1,1:1))
       ALLOCATE(map%a__TEB_l_m(1:1,0:map%lmax,0:map%mmax))
    ENDIF

    map%header       = ''
    map%map__pix_IQU = 0.0_DP
    map%a__TEB_l_m   = 0.0_DPC


  END SUBROUTINE INIT_MAP

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Destructor routine: Map
!-----------------------------------------------------------------------

  SUBROUTINE FREE_MAP(map)

    CLASS(SKYMAP),                          INTENT(INOUT) :: map

    DEALLOCATE(map%header)
    DEALLOCATE(map%map__pix_IQU)
    DEALLOCATE(map%a__TEB_l_m)

    CALL map%sht%FREE()


  END SUBROUTINE FREE_MAP

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Constructor routine: Parameter
!-----------------------------------------------------------------------

  SUBROUTINE INIT_PAR(par)

    CLASS(PARAMETERS),                      INTENT(INOUT) :: par

    par%pol          = .FALSE.
    par%nr_sims      = 0
    par%sim_nr_start = 0
    par%nside        = 0
    par%npix         = 0
    par%lmax         = 0
    par%mmax         = 0
    par%nr_shells    = 0
    par%mpi_error    = 0
    par%root         = 0
    par%loop_min     = 0
    par%loop_max     = 0


  END SUBROUTINE INIT_PAR

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Destructor routine: Parameter
!-----------------------------------------------------------------------

  SUBROUTINE FREE_PAR(par)

    CLASS(PARAMETERS),                      INTENT(INOUT) :: par

    IF (ALLOCATED(par%sht%zbounds__lat))&
         & DEALLOCATE(par%sht%zbounds__lat)
    IF (ALLOCATED(par%sht%w8ring__ring_IQU))&
         & DEALLOCATE(par%sht%w8ring__ring_IQU)
    IF (ALLOCATED(par%radii__r)) DEALLOCATE(par%radii__r)


  END SUBROUTINE FREE_PAR

!-----------------------------------------------------------------------

END MODULE MAP_MOD_CLASSES
