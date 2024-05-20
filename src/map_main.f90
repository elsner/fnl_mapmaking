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


PROGRAM MAP_MAIN

  USE MPI_F08
  USE HEALPIX_TYPES
  USE MAP_MOD_CLASSES
  USE MAP_MOD_MISC
  USE MAP_MOD_SAMPLE

  IMPLICIT NONE

  TYPE(FILENAMES)                                         :: files
  TYPE(PARAMETERS)                                        :: par
  TYPE(RND_SEEDS)                                         :: rng
  TYPE(WEIGHTS)                                           :: weights_p
  TYPE(WEIGHTS)                                           :: weights_mu
  TYPE(COVARIANCE)                                        :: cov_p
  TYPE(COVARIANCE)                                        :: cov_mu
  TYPE(SKYMAP)                                            :: signal_l
  TYPE(SKYMAP)                                            :: signal_nl
  TYPE(SKYMAP)                                            :: signal_nnl
  TYPE(SKYMAP)                                            :: signal_mu
  INTEGER(I4B)                                            :: sim_nr

!-----------------------------------------------------------------------
! Initialize MPI
!-----------------------------------------------------------------------

  CALL MPI_INIT(par%mpi_error)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, par%nr_nodes, par%mpi_error)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, par%thread_id, par%mpi_error)

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Parse input parameters
!-----------------------------------------------------------------------

  CALL PARSE_INPUT_PARAMETERS("./config/input.cfg", files, par)

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------

  CALL ALLOCATE_ARRAYS(weights_p, weights_mu, cov_p, cov_mu, rng,&
       & signal_l, signal_nl, signal_nnl, signal_mu, par)

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Read auxiliary files
!-----------------------------------------------------------------------

  CALL READ_AUX_FILES(weights_p, weights_mu, cov_p, cov_mu, rng,&
       & files, par)

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Compute square root of covariance matrices
!-----------------------------------------------------------------------

  CALL COMPUTE_MATRICES(cov_p, par)

  CALL COMPUTE_MATRICES(cov_mu, par)

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Compute loop boundaries for each thread
!-----------------------------------------------------------------------

  CALL GET_FILE_NUMBER(par, files)

  CALL GET_LOOP_BOUNDS(par, par%sim_nr_start, par%nr_sims&
       & + par%sim_nr_start - 1, par%thread_id)

!-----------------------------------------------------------------------


  CALL CLOCK()

!-----------------------------------------------------------------------
! Sample maps
!-----------------------------------------------------------------------

  DO sim_nr = par%loop_min,par%loop_max

     WRITE(*,'(/20X, "--- Map ", 1X, I4, " ---"/)') sim_nr


   !--------------------------------------------------------------------
   ! Initialize random seeds
   !--------------------------------------------------------------------

     CALL GENERATE_RND_SEEDS(rng, sim_nr, par)

   !--------------------------------------------------------------------



   !--------------------------------------------------------------------
   ! Generate maps
   !--------------------------------------------------------------------

     CALL SAMPLE_MAP(signal_l, signal_nl, signal_nnl, signal_mu,&
          & cov_p, cov_mu, weights_p, weights_mu, rng, par)

   !--------------------------------------------------------------------



   !--------------------------------------------------------------------
   ! Write Gaussian signal alms to file
   !--------------------------------------------------------------------

     CALL OUTPUT_AS_MAP(signal_l, rng, sim_nr, files%out_signal_l)

!     CALL OUTPUT_ALM(signal_l, rng, sim_nr, files%out_signal_l)

   !--------------------------------------------------------------------



   !--------------------------------------------------------------------
   ! Write non-Gaussian signal alms to file
   !--------------------------------------------------------------------

     CALL OUTPUT_AS_MAP(signal_nl, rng, sim_nr, files%out_signal_nl)

!     CALL OUTPUT_ALM(signal_nl, rng, sim_nr, files%out_signal_nl)

   !--------------------------------------------------------------------



   !--------------------------------------------------------------------
   ! Write non-Gaussian signal alms to file
   !--------------------------------------------------------------------

     CALL OUTPUT_AS_MAP(signal_nnl, rng, sim_nr, files%out_signal_nnl)

!     CALL OUTPUT_ALM(signal_nnl, rng, sim_nr, files%out_signal_nnl)

   !--------------------------------------------------------------------



   !--------------------------------------------------------------------
   ! Write mu signal alms to file
   !--------------------------------------------------------------------

     CALL OUTPUT_AS_MAP(signal_mu, rng, sim_nr, files%out_signal_mu)

!     CALL OUTPUT_ALM(signal_mu, rng, sim_nr, files%out_signal_mu)

   !--------------------------------------------------------------------


  ENDDO

!-----------------------------------------------------------------------

  CALL CLOCK()


!-----------------------------------------------------------------------
! Cleanup
!-----------------------------------------------------------------------

  CALL DEALLOCATE_ARRAYS(weights_p, weights_mu, cov_p, cov_mu, rng,&
       & signal_l, signal_nl, signal_nnl, signal_mu, par)

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! MPI Cleanup
!-----------------------------------------------------------------------

  CALL MPI_FINALIZE(par%mpi_error)

!-----------------------------------------------------------------------

END PROGRAM MAP_MAIN
