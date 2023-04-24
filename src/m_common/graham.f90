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

subroutine GRAHAM(PSI, over, ldo, rn, lda, NDIME, MH, info)

    implicit none
    integer mh, ndime, i, lda, ldo, info
    real(8) over(ldo, *)
    real(8) PSI(lda, *), RN(*)
    real(8), dimension(:), allocatable :: sc1, sc2, sp
    real(8), dimension(:, :), allocatable :: mat

    ! This subroutine uses gramh-schmidt orthogonalization
    ! in a metric given by the overlap matrix over
    ! The output psi will be orthogonal vectors in this metric such that:

    ! sum_kl [ over(k,l)   psi(k,i) x psi(l,j)] = 1  if i=j, 0 otherwise.
    ! where  1 <= k,l <= ndime
    ! On output info=0 successful orthogonalization
    ! info > 0  MH-info linearly independent vectors found.
    ! rn(i) is the normalization of the vectors.
    ! the components of psi such that k>ndime are not touched.

    if (mh .le. 0) return ! do nothing

    allocate (sc1(ndime), sp(mh), mat(ndime, mh), sc2(ndime))
    mat = 0.d0
    sc1 = 0.d0
    sc2 = 0.d0
    sp = -1.d0
    info = 0

    call dgemm('N', 'N', ndime, mh, ndime, 1.d0, over, ldo, psi, lda, 0.d0&
            &, mat, ndime)
    !
    RN(1) = dsqrt(sum(psi(1:ndime, 1)*mat(:, 1)))
    if (rn(1) .gt. 0.d0) then
        psi(1:ndime, 1) = psi(1:ndime, 1)/rn(1)
    else
        info = 1
        psi(1:ndime, 1) = 0.d0
    end if
    do I = 2, MH
        call DGEMV('T', ndime, I - 1, 1.d0, PSI, lda, mat(1, I), 1, 0.d0, SP, 1)
        call DGEMV('N', ndime, I, 1.d0, PSI, lda, SP, 1, 0.d0, SC1, 1)
        call dgemv('N', ndime, ndime, 1.d0, over, ldo, sc1, 1, 0.d0, sc2, 1)
        RN(i) = dsqrt(sum(sc1(:)*sc2(:)))
        if (rn(i) .gt. 0.d0) then
            psi(1:ndime, i) = -sc1(1:ndime)/rn(i)
        else
            psi(1:ndime, i) = 0.d0
            info = info + 1
        end if
    end do
    deallocate (sc1, sc2, sp, mat)
    return
end subroutine GRAHAM

subroutine GRAHAM_complex(PSI, over, ldo, rn, lda, NDIME, MH, info)

    use constants, only: zzero, zone

    implicit none
    integer mh, ndime, i, lda, ldo, info
    complex(8) over(ldo, *), PSI(lda, *)
    real(8) :: RN(*)
    complex(8), dimension(:), allocatable :: sc1, sc2, sp
    complex(8), dimension(:, :), allocatable :: mat

    !          This subroutine uses gramh-schmidt orthogonalization
    !          in a metric given by the overlap matrix over
    !          The output psi will be orthogonal vectors in this metric such that:

    !          sum_kl [ over(k,l)   psi(k,i) x psi(l,j)] = 1  if i=j, 0 otherwise.
    !          where  1 <= k,l <= ndime
    !          On output info=0 successful orthogonalization
    !          info > 0  MH-info linearly independent vectors found.
    !          rn(i) is the normalization of the vectors.
    !          the components of psi such that k> ndime are not touched.

    if (mh .le. 0) return ! do nothing

    allocate (sc1(ndime), sp(mh), mat(ndime, mh), sc2(ndime))
    mat = 0.d0
    sc1 = 0.d0
    sc2 = 0.d0
    sp = -1.d0
    info = 0

    call zgemm('N', 'N', ndime, mh, ndime, zone, over, ldo, psi, lda, zzero &
            &, mat, ndime)
    !
    RN(1) = sqrt(sum(conjg(psi(1:ndime, 1))*mat(:, 1)))
    if (rn(1) .gt. 0.d0) then
        psi(1:ndime, 1) = psi(1:ndime, 1)/rn(1)
    else
        info = 1
        psi(1:ndime, 1) = zzero
    end if
    !
    ! NB: in the case of complex vectors the scalar product has a different definition
    ! in order to keep valid its main properties:
    ! a \dot b = \sum_i a_i b_i^*
    !
    do I = 2, MH
        call ZGEMV('C', ndime, I - 1, zone, PSI, lda, mat(1, I), 1, zzero, SP, 1)
        call ZGEMV('N', ndime, I, zone, PSI, lda, SP, 1, zzero, sc1, 1)
        call ZGEMV('N', ndime, ndime, zone, over, ldo, sc1, 1, zzero, sc2, 1)
        RN(i) = sqrt(sum(conjg(sc1(:))*sc2(:)))
        if (rn(i) .gt. 0.d0) then
            psi(1:ndime, i) = -sc1(1:ndime)/rn(i)
        else
            psi(1:ndime, i) = zzero
            info = info + 1
        end if
    end do
    deallocate (sc1, sc2, sp, mat)
    return
end subroutine GRAHAM_complex
