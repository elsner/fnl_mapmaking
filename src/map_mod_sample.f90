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


MODULE MAP_MOD_SAMPLE

  USE HEALPIX_TYPES
  USE MAP_MOD_MISC

  IMPLICIT NONE


CONTAINS

!-----------------------------------------------------------------------
! Generate maps
!-----------------------------------------------------------------------

  SUBROUTINE SAMPLE_MAP(signal_l, signal_nl, signal_nnl, signal_mu,&
          & cov_p, cov_mu, weights_p, weights_mu, rng, par)

    USE UDGRADE_NR
    USE PIX_TOOLS,     ONLY: NSIDE2NPIX
#ifdef SILO
    USE VISIT_INTERFACE
#endif

    TYPE(PARAMETERS),                       INTENT(IN)    :: par
    TYPE(RND_SEEDS),                        INTENT(IN)    :: rng
    TYPE(WEIGHTS),                          INTENT(IN)    :: weights_p
    TYPE(WEIGHTS),                          INTENT(IN)    :: weights_mu
    TYPE(COVARIANCE),                       INTENT(IN)    :: cov_p
    TYPE(COVARIANCE),                       INTENT(IN)    :: cov_mu
    TYPE(SKYMAP),                           INTENT(INOUT) :: signal_l
    TYPE(SKYMAP),                           INTENT(INOUT) :: signal_nl
    TYPE(SKYMAP),                           INTENT(INOUT) :: signal_nnl
    TYPE(SKYMAP),                           INTENT(INOUT) :: signal_mu

    INTEGER(I4B)                                          :: shell
    TYPE(SKYMAP)                                          :: phi_l
    TYPE(SKYMAP)                                          :: phi_mu
    TYPE(SKYMAP)                                          :: phi_transform
#ifdef SILO
    REAL(DP),          DIMENSION(:,:),      ALLOCATABLE   :: phi_l__pix_r
    REAL(DP),          DIMENSION(:,:),      ALLOCATABLE   :: phi_nl__pix_r

    ALLOCATE(phi_l__pix_r(0:par%npix-1,1:par%nr_shells))
    ALLOCATE(phi_nl__pix_r(0:par%npix-1,1:par%nr_shells))
    phi_l__pix_r  = 0.0_DP
    phi_nl__pix_r = 0.0_DP
#endif

    CALL phi_l%INIT(par, polarization=.FALSE.)
    CALL phi_mu%INIT(par, polarization=.FALSE.)
    CALL phi_transform%INIT(par, polarization=.FALSE.)
    signal_l%a__TEB_l_m   = 0.0_DPC
    signal_nl%a__TEB_l_m  = 0.0_DPC
    signal_nnl%a__TEB_l_m = 0.0_DPC
    signal_mu%a__TEB_l_m  = 0.0_DPC

    DO shell = 1,par%nr_shells

       phi_l%map__pix_IQU         = 0.0_DP
       phi_mu%map__pix_IQU        = 0.0_DP
       phi_transform%map__pix_IQU = 0.0_DP
       phi_l%a__TEB_l_m           = 0.0_DPC
       phi_mu%a__TEB_l_m          = 0.0_DPC
       phi_transform%a__TEB_l_m   = 0.0_DPC

       CALL SAMPLE_PHI(phi_l, phi_mu, cov_p, cov_mu, rng, shell)

       CALL ADD_2_SIGNAL(signal_l, phi_l, shell, weights_p)

       CALL CALCULATE_PHI_SQR(phi_transform, phi_l, phi_l, cov_p,&
            & shell)

       CALL ADD_2_SIGNAL(signal_nl, phi_transform, shell, weights_p)

#ifdef SILO
       phi_l__pix_r(:,shell)  = phi_l%map__pix_IQU(:,1)
       phi_nl__pix_r(:,shell) = phi_transform%map__pix_IQU(:,1)
#endif

       CALL CALCULATE_PHI_SQR(phi_transform, phi_mu, phi_l, cov_mu,&
            & shell)

       CALL ADD_2_SIGNAL(signal_mu, phi_transform, shell, weights_mu)

       CALL CALCULATE_PHI_CUBE(phi_transform, phi_l, phi_l, phi_l)

       CALL ADD_2_SIGNAL(signal_nnl, phi_transform, shell, weights_p)

    ENDDO

#ifdef SILO
    CALL HEALPIX_MESH_MANY(phi_l__pix_r, par%nside, par%npix,&
         & par%nr_shells, 'db_phi_l.dat', par%radii__r)
    CALL HEALPIX_MESH_MANY(phi_nl__pix_r, par%nside, par%npix,&
         & par%nr_shells, 'db_phi_nl.dat', par%radii__r)

    DEALLOCATE(phi_l__pix_r)
    DEALLOCATE(phi_nl__pix_r)
#endif

    CALL phi_l%FREE()
    CALL phi_mu%FREE()
    CALL phi_transform%FREE()


  END SUBROUTINE SAMPLE_MAP

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Sample phi
!-----------------------------------------------------------------------

  SUBROUTINE SAMPLE_PHI(phi_l, phi_mu, cov_p, cov_mu, rng, shell)

    USE BLAS95,        ONLY: DOTC

    TYPE(RND_SEEDS),                        INTENT(IN)    :: rng
    TYPE(COVARIANCE),                       INTENT(IN)    :: cov_p
    TYPE(COVARIANCE),                       INTENT(IN)    :: cov_mu
    INTEGER(I4B),                           INTENT(IN)    :: shell
    TYPE(SKYMAP),                           INTENT(INOUT) :: phi_l
    TYPE(SKYMAP),                           INTENT(INOUT) :: phi_mu

    INTEGER(I4B)                                          :: l
    INTEGER(I4B)                                          :: m
    COMPLEX(DPC),      DIMENSION(:),        ALLOCATABLE   :: gauss_rand__r

!$OMP PARALLEL PRIVATE(l, m, gauss_rand__r)

    ALLOCATE(gauss_rand__r(1:rng%nr_shells))
    gauss_rand__r = 0.0_DPC

!$OMP DO SCHEDULE(GUIDED)
    DO m = 0,phi_l%mmax
       DO l = MAX(2,m),phi_l%lmax

          CALL GENERATE_RND_VECTOR(gauss_rand__r, rng, l, m)

          phi_l%a__TEB_l_m(1,l,m)&
               & = DOTC(DCMPLX(cov_p%sqrt__r1_r2_l(:,shell,l), 0.0_DP),&
               & gauss_rand__r)

          phi_mu%a__TEB_l_m(1,l,m)&
               & = DOTC(DCMPLX(cov_mu%sqrt__r1_r2_l(:,shell,l), 0.0_DP),&
               & gauss_rand__r)

       ENDDO
    ENDDO
!$OMP ENDDO

    DEALLOCATE(gauss_rand__r)

!$OMP END PARALLEL

    CALL ALM_2_MAP(phi_l)
    CALL ALM_2_MAP(phi_mu)


  END SUBROUTINE SAMPLE_PHI

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Generate phi_l**2 map from phi_l
!-----------------------------------------------------------------------

  SUBROUTINE CALCULATE_PHI_SQR(phi_out, phi_1, phi_2, cov_p, shell)

    TYPE(COVARIANCE),                       INTENT(IN)    :: cov_p
    TYPE(SKYMAP),                           INTENT(IN)    :: phi_1
    TYPE(SKYMAP),                           INTENT(IN)    :: phi_2
    INTEGER(I4B),                           INTENT(IN)    :: shell
    TYPE(SKYMAP),                           INTENT(INOUT) :: phi_out

    phi_out%map__pix_IQU = phi_1%map__pix_IQU * phi_2%map__pix_IQU&
         & - cov_p%var__r(shell)

    CALL MAP_2_ALM(phi_out)


  END SUBROUTINE CALCULATE_PHI_SQR

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Generate phi_l**3 map from phi_l
!-----------------------------------------------------------------------

  SUBROUTINE CALCULATE_PHI_CUBE(phi_out, phi_1, phi_2, phi_3)

    TYPE(SKYMAP),                           INTENT(IN)    :: phi_1
    TYPE(SKYMAP),                           INTENT(IN)    :: phi_2
    TYPE(SKYMAP),                           INTENT(IN)    :: phi_3
    TYPE(SKYMAP),                           INTENT(INOUT) :: phi_out

    phi_out%map__pix_IQU = phi_1%map__pix_IQU&
         & * phi_2%map__pix_IQU * phi_3%map__pix_IQU

    CALL MAP_2_ALM(phi_out)


  END SUBROUTINE CALCULATE_PHI_CUBE

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Project phi_l/nl/nnl/mu maps to signal_l/nl/nnl/mu map
!-----------------------------------------------------------------------

  SUBROUTINE ADD_2_SIGNAL(signal, phi, shell, wgt)

    TYPE(WEIGHTS),                          INTENT(IN)    :: wgt
    TYPE(SKYMAP),                           INTENT(IN)    :: phi
    INTEGER(I4B),                           INTENT(IN)    :: shell
    TYPE(SKYMAP),                           INTENT(INOUT) :: signal

    INTEGER(I4B)                                          :: l
    INTEGER(I4B)                                          :: m

!$OMP PARALLEL PRIVATE(l, m)

!$OMP DO SCHEDULE(GUIDED)
    DO m = 0,signal%mmax
       DO l = MAX(2,m),signal%lmax
          signal%a__TEB_l_m(1,l,m) = signal%a__TEB_l_m(1,l,m)&
               & + DCMPLX(wgt%m__r_TE_l(shell,1,l), 0.0_DP)&
               & * phi%a__TEB_l_m(1,l,m)
       ENDDO
    ENDDO
!$OMP ENDDO NOWAIT

    IF (signal%pol) THEN

!$OMP DO SCHEDULE(GUIDED)
       DO m = 0,signal%mmax
          DO l = MAX(2,m),signal%lmax
             signal%a__TEB_l_m(2,l,m) = signal%a__TEB_l_m(2,l,m)&
                  & + DCMPLX(wgt%m__r_TE_l(shell,2,l), 0.0_DP)&
                  & * phi%a__TEB_l_m(1,l,m)
          ENDDO
       ENDDO
!$OMP ENDDO

    ENDIF

!$OMP END PARALLEL


  END SUBROUTINE ADD_2_SIGNAL

!-----------------------------------------------------------------------

END MODULE MAP_MOD_SAMPLE
