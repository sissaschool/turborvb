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

subroutine GRAHAM_SCALAPACK(PSI, over, rn, lda, NDIME, MH, info, desc, desch, nlax, rank)

    use descriptors
    use allio, only: commrep_mpi, rankrep, nprocrep
    implicit none
    integer mh, ndime, i, ii, jj, lda, nlax, ierr, rank, dimbuf, dime, &
        nslice, nlast, proclast, ind, info
    real(8) :: over(nlax, nlax)
    real(8), dimension(:, :), allocatable :: over_part, over_new
    integer :: desch(20), ir, ic, j
    integer :: desc(descla_siz_)
    real(8) PSI(lda, *), RN(*)
    real(8), dimension(:), allocatable :: sc1, sc2, sp
    real(8), dimension(:, :), allocatable :: mat

#ifdef PARALLEL
    include 'mpif.h'
#endif

#if defined (_OPENMP) && defined (__NOOMPDFT)
    integer, external :: omp_get_max_threads
    integer old_threads
    old_threads = omp_get_max_threads()
    call omp_set_num_threads(1) ! scalar code
#endif
    if (mh .le. 0) return ! do nothing

    !          This subroutine uses gramh-schmidt orthogonalization
    !          in a metric given by the overlap matrix "over".
    !          The output PSI will be orthogonal vectors in this metric such that:

    !          sum_kl [ over(k,l)   psi(k,i) x psi(l,j)] = 1  if i=j, 0 otherwise.
    !          where  1 <= k,l <= ndime
    !          On output info=0 successful orthogonalization
    !          info > 0  MH-info linearly independent vectors found.
    !          rn(i) is the normalization of the vectors.
    !          the components of psi such that k> ndime are not touched.

    allocate (sc1(ndime), sp(mh), mat(ndime, mh), sc2(ndime))
    dime = ndime*mh
    mat = 0.d0
    sc1 = 0.d0
    sc2 = 0.d0
    sp = -1.d0
    info = 0
#ifdef     PARALLEL
    nslice = ndime/nprocrep
    proclast = nprocrep - 1
    if (nslice*nprocrep .ne. ndime) then
        nslice = ndime/proclast
        nlast = ndime - nslice*proclast
    else
        nlast = nslice
    end if
#endif

#ifdef __SCALAPACK
    ir = 0
    ic = 0
    if (desc(lambda_node_) > 0) then
        allocate (over_part(nlax, nlax))
        allocate (over_new(nlax, nlax))
        ir = desc(ilar_)
        ic = desc(ilac_)
        over_part = 0.d0
        over_new = 0.d0

        do j = 1, descla(nlac_)
            do i = 1, descla(nlar_)
                if ((j + ic - 1) <= mh) then
                    over_part(i, j) = psi((i + ir - 1), (j + ic - 1))
                end if
            end do
        end do

        call PDGEMM('N', 'N', ndime, ndime, ndime, 1.0d0, over, 1, 1,&
             & desch, over_part, 1, 1, desch, 0.0d0, over_new, 1, 1, desch)

        do j = 1, descla(nlac_)
            do i = 1, descla(nlar_)
                if ((j + ic - 1) <= mh) then
                    mat((i + ir - 1), (j + ic - 1)) = over_new(i, j)
                end if
            end do
        end do

        deallocate (over_new, over_part)

    end if
  !!  now distribute the matrix mat
#ifdef PARALLEL
    call reduce_base_real(dime, mat, commrep_mpi, -1)

    !       call  reduce_base(dime,mat,buff,dimbuf)
#endif

#endif
    !          if(rank.eq.0) then
    !
    RN(1) = dsqrt(sum(psi(1:ndime, 1)*mat(:, 1)))
    if (rn(1) .gt. 0.d0) then
        psi(1:ndime, 1) = psi(1:ndime, 1)/rn(1)
    else
        info = info + 1
        psi(1:ndime, 1) = 0.d0
    end if

    !          endif

    do I = 2, MH
        sc2 = 0.d0
        !           if(rank.eq.0) then
#ifdef PARALLEL
        ind = rankrep*nslice + 1
        SP(1:I - 1) = 0.d0
        if (rankrep .ne. proclast .and. nslice .gt. 0) then
            call DGEMV('T', nslice, I - 1, 1.d0, PSI(ind, 1), lda, mat(ind, I), 1, 0.d0, SP, 1)
        elseif (nlast .gt. 0) then
            call DGEMV('T', nlast, I - 1, 1.d0, PSI(ind, 1), lda, mat(ind, I), 1, 0.d0, SP, 1)
        end if
        call reduce_base_real(I - 1, sp, commrep_mpi, -1)

        sc1 = 0.d0

        if (rankrep .ne. proclast .and. nslice .gt. 0) then
            call DGEMV('N', nslice, I, 1.d0, PSI(ind, 1), lda, SP, 1, 0.d0, SC1(ind), 1)
        elseif (nlast .gt. 0) then
            call DGEMV('N', nlast, I, 1.d0, PSI(ind, 1), lda, SP, 1, 0.d0, SC1(ind), 1)
        end if
        call reduce_base_real(ndime, sc1, commrep_mpi, -1)
#else

        call DGEMV('T', ndime, I - 1, 1.d0, PSI, lda, mat(1, I), 1, 0.d0, SP, 1)
        call DGEMV('N', ndime, I, 1.d0, PSI, lda, SP, 1, 0.d0, SC1, 1)
#endif
        !           endif
#ifdef __SCALAPACK
        !#ifdef   PARALLEL
        !        call mpi_bcast(sc1,ndime,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        !#endif

        if (desc(lambda_node_) > 0) then

            ir = desc(ilar_)
            ic = desc(ilac_)

            do jj = 1, descla(nlac_)
                do ii = 1, descla(nlar_)
                    sc2(ii + ir - 1) = sc2(ii + ir - 1) + over(ii, jj)*sc1(jj + ic - 1)
                end do
            end do

        end if

#ifdef PARALLEL
        !        sc3=sc2
        !        call mpi_allreduce(sc3,sc2,ndime,MPI_DOUBLE_PRECISION&
        !    &,MPI_SUM,MPI_COMM_WORLD,ierr)
        call reduce_base_real(ndime, sc2, commrep_mpi, -1)
#endif

#else
        call dgemv('N', ndime, ndime, 1.d0, over, lda, sc1, 1, 0.d0, sc2, 1)
#endif
        !         if(rank.eq.0) then
        RN(i) = dsqrt(sum(sc1(:)*sc2(:)))
        if (rn(i) .gt. 0.d0) then
            psi(1:ndime, i) = -sc1(1:ndime)/rn(i)
        else
            psi(1:ndime, i) = 0.d0
            info = info + 1
        end if
        !         endif
    end do
    deallocate (sc1, sc2, sp, mat)
    !#ifdef PARALLEL
    !        call mpi_bcast(psi,lda*mh,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    !        call mpi_bcast(rn,mh,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    !#endif

#if defined (_OPENMP) && defined (__NOOMPDFT)
    call omp_set_num_threads(old_threads) ! restore the previous threads
#endif
    return
end subroutine GRAHAM_SCALAPACK

subroutine GRAHAM_SCALAPACK_complex(PSI, over, rn, lda, NDIME, MH, info, desc, desch, nlax, rank)
    use descriptors
    use allio, only: commrep_mpi, rankrep, nprocrep
    use constants, only: zone, zzero
    implicit none
    integer, intent(in) :: lda, nlax, rank, desch(20), desc(descla_siz_), ndime, mh
    integer i, ii, jj, ierr, dimbuf, dime&
            &, nslice, nlast, proclast, ind, info, ir, ic, j
    complex(8) over(nlax, nlax)
    complex(8), dimension(:, :), allocatable :: over_part, over_new
    real(8) :: RN(*)
    complex(8) PSI(lda, *)
    complex(8), dimension(:), allocatable :: sc1, sc2, sp
    complex(8), dimension(:, :), allocatable :: mat
#ifdef PARALLEL
    include 'mpif.h'
#endif

#if defined (_OPENMP) && defined (__NOOMPDFT)
    integer, external :: omp_get_max_threads
    integer old_threads
    old_threads = omp_get_max_threads()
    call omp_set_num_threads(1) ! scalar code
#endif

    if (mh .le. 0) return ! do nothing

    !          This subroutine uses gramh-schmidt orthogonalization
    !          in a metric given by the overlap matrix over
    !          The output psi will be orthogonal vectors in this metric such that:

    !          sum_kl [ over(k,l)   psi(k,i) x psi(l,j)] = 1  if i=j, 0 otherwise.
    !          where  1 <= k,l <= ndime
    !          On output info=0 successful orthogonalization
    !          info > 0  MH-info linearly independent vectors found.
    !          rn(i) is the normalization of the vectors.
    !          the components of psi such that k> ndime are not touched.

    allocate (sc1(ndime), sp(mh), mat(ndime, mh), sc2(ndime))
    dime = ndime*mh
    mat = zzero
    sc1 = zzero
    sc2 = zzero
    sp = -zone
    info = 0

#ifdef PARALLEL
    nslice = ndime/nprocrep
    proclast = nprocrep - 1
    if (nslice*nprocrep .ne. ndime) then
        nslice = ndime/proclast
        nlast = ndime - nslice*proclast
    else
        nlast = nslice
    end if
#endif

#ifdef __SCALAPACK

    ir = 0
    ic = 0
    if (desc(lambda_node_) > 0) then

        allocate (over_part(nlax, nlax))
        allocate (over_new(nlax, nlax))
        ir = desc(ilar_)
        ic = desc(ilac_)
        over_part = zzero
        over_new = zzero

        do j = 1, descla(nlac_)
            do i = 1, descla(nlar_)
                if ((j + ic - 1) <= mh) then
                    over_part(i, j) = psi((i + ir - 1), (j + ic - 1))
                end if
            end do
        end do

        call PZGEMM('N', 'N', ndime, ndime, ndime, zone, over, 1, 1,&
             & desch, over_part, 1, 1, desch, zzero, over_new, 1, 1, desch)

        do j = 1, descla(nlac_)
            do i = 1, descla(nlar_)
                if ((j + ic - 1) <= mh) then
                    mat((i + ir - 1), (j + ic - 1)) = over_new(i, j)
                end if
            end do
        end do

        deallocate (over_new, over_part)

    end if
    !  now distribute the matrix mat
#ifdef PARALLEL
    call reduce_base_complex(dime, mat, commrep_mpi, -1)
#endif

#endif

    RN(1) = dsqrt(real(sum(dconjg(psi(1:ndime, 1))*mat(:, 1))))
    if (rn(1) .gt. 0.d0) then
        psi(1:ndime, 1) = psi(1:ndime, 1)/rn(1)
    else
        info = info + 1
        psi(1:ndime, 1) = zzero
    end if

    ind = rankrep*nslice + 1

    do I = 2, MH

        sc2 = zzero
#ifdef PARALLEL
        ind = rankrep*nslice + 1
        SP(1:I - 1) = zzero
        if (rankrep .ne. proclast .and. nslice .gt. 0) then
            call ZGEMV('C', nslice, I - 1, zone, PSI(ind, 1), lda, mat(ind, I), 1, zzero, SP, 1)
        elseif (nlast .gt. 0) then
            call ZGEMV('C', nlast, I - 1, zone, PSI(ind, 1), lda, mat(ind, I), 1, zzero, SP, 1)
        end if
        call reduce_base_complex(I - 1, sp, commrep_mpi, -1)

        sc1 = zzero

        if (rankrep .ne. proclast .and. nslice .gt. 0) then
            call ZGEMV('N', nslice, I, zone, PSI(ind, 1), lda, SP, 1, zzero, SC1(ind), 1)
        elseif (nlast .gt. 0) then
            call ZGEMV('N', nlast, I, zone, PSI(ind, 1), lda, SP, 1, zzero, SC1(ind), 1)
        end if
        call reduce_base_complex(ndime, sc1, commrep_mpi, -1)
#else
        call ZGEMV('C', ndime, I - 1, zone, PSI, lda, mat(1, I), 1, zzero, SP, 1)
        call ZGEMV('N', ndime, I, zone, PSI, lda, SP, 1, zzero, SC1, 1)
#endif

#ifdef __SCALAPACK

        if (desc(lambda_node_) > 0) then

            ir = desc(ilar_)
            ic = desc(ilac_)

            do jj = 1, descla(nlac_)
                do ii = 1, descla(nlar_)
                    sc2(ii + ir - 1) = sc2(ii + ir - 1) + over(ii, jj)*sc1(jj + ic - 1)
                end do
            end do

        end if

#ifdef PARALLEL
        call reduce_base_complex(ndime, sc2, commrep_mpi, -1)
#endif

#else
        call zgemv('N', ndime, ndime, zone, over, lda, sc1, 1, zzero, sc2, 1)
#endif

        RN(i) = dsqrt(real(sum(dconjg(sc1(:))*sc2(:))))
        if (rn(i) .gt. 0.d0) then
            psi(1:ndime, i) = -sc1(1:ndime)/rn(i)
        else
            psi(1:ndime, i) = zzero
            info = info + 1
        end if
    end do
    deallocate (sc1, sc2, sp, mat)
#if defined (_OPENMP) && defined (__NOOMPDFT)
    call omp_set_num_threads(old_threads) ! restore the previous threads
#endif
    return
end subroutine GRAHAM_SCALAPACK_complex
