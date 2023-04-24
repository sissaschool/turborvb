! Copyright (C) 2022 TurboRVB group
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.

subroutine gather_avgs(wbra, ener, etot, wtotf, wbra_t, &
                       ener_t, etot_t, wtotf_t, np3, pdim)

    use allio, only: commcolrep_mpi
    use kpoints_mod, only: kaverage, decoupled_run

    integer, intent(in) :: np3, pdim
    real(8), intent(in) :: wbra, ener, wtotf(2), etot(np3)
    real(8), intent(inout) :: ener_t(pdim), wbra_t(pdim), wtotf_t(2, pdim), etot_t(np3, pdim)
    integer :: ierr

# if defined PARALLEL
    include "mpif.h"
# endif

    ener_t = 0.d0
    wbra_t = 0.d0
    wtotf_t = 0.d0
    etot_t = 0.d0

# if defined PARALLEL
    ! correlation functions
    call mpi_gather(etot, np3, MPI_DOUBLE_PRECISION, etot_t, &
                    np3, MPI_DOUBLE_PRECISION, 0, commcolrep_mpi, ierr)
    ! weights
    call mpi_gather(wtotf, 2, MPI_DOUBLE_PRECISION, wtotf_t, &
                    2, MPI_DOUBLE_PRECISION, 0, commcolrep_mpi, ierr)
    ! energy (real part)
    call mpi_gather(ener, 1, MPI_DOUBLE_PRECISION, ener_t, &
                    1, MPI_DOUBLE_PRECISION, 0, commcolrep_mpi, ierr)
    ! wbra (total weight)
    call mpi_gather(wbra, 1, MPI_DOUBLE_PRECISION, wbra_t, &
                    1, MPI_DOUBLE_PRECISION, 0, commcolrep_mpi, ierr)
# endif

    return

end subroutine gather_avgs

subroutine gather_wconfn(ist, ien, w_tot, w_loc_kp, nw, nk)

    use allio, only: rankrep, nprocrep, commrep_mpi, in1, psip
    implicit none

    integer, intent(in) :: nw, nk, ien, ist
    real(8), intent(in) :: w_tot(nw)
    real(8), intent(inout) :: w_loc_kp(nw/nk)
    integer ist_loc, ien_loc, nw_rep, ierr
# if defined PARALLEL
    include "mpif.h"
# endif

    nw_rep = nw/nk
    ist_loc = rankrep*(nw_rep/nprocrep) + 1
    ien_loc = (rankrep + 1)*(nw_rep/nprocrep)

    ! collect local weights on each pool
    ! in1 = # of walkers per process must always be equal to 1
    w_loc_kp(ist_loc:ien_loc) = w_tot(ist:ien)
    psip(1:in1) = w_loc_kp(ist_loc:ien_loc)
# if defined PARALLEL
    call mpi_gather(psip, in1, MPI_DOUBLE_PRECISION, w_loc_kp, &
                    in1, MPI_DOUBLE_PRECISION, 0, commrep_mpi, ierr)
# endif

    return

end subroutine gather_wconfn

subroutine scatter_wconfn(ist, ien, w_tot, w_loc_kp, nw, nk)

    use allio, only: rankrep, nprocrep, commrep_mpi, psip, in1
    implicit none

    integer, intent(in) :: nw, nk, ien, ist
    real(8), intent(inout) :: w_tot(nw), w_loc_kp(nw/nk)
    integer ist_loc, ien_loc, nw_rep, ierr
# if defined PARALLEL
    include "mpif.h"
# endif

    nw_rep = nw/nk
    ist_loc = rankrep*(nw_rep/nprocrep) + 1
    ien_loc = (rankrep + 1)*(nw_rep/nprocrep)

# if defined PARALLEL
    call mpi_scatter(w_loc_kp(ist_loc), in1, MPI_DOUBLE_PRECISION, &
                     psip, in1, MPI_DOUBLE_PRECISION, 0, commrep_mpi, ierr)
# endif
    w_loc_kp(ist_loc:ien_loc) = psip(1:in1)
    w_tot(ist:ien) = w_loc_kp(ist_loc:ien_loc)

    return

end subroutine scatter_wconfn
