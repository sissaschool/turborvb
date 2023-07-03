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

subroutine shells(ndim, a, cut, nshlls, rkcomp, rknorm, kmult           &
        &, nvects, mnkv, mnsh, mdim)
    !      implicit real*8 (a-h,o-z)
    implicit none
    integer ndim, nshlls, mnsh, kmult(0:*), nvects, mnkv, npts, l            &
            &, nkspan(3), i, icount(3), nzero, j, kj, jp, mdim
    real*8 a(*), rkcomp(mdim, *), rknorm(0:*)                            &
            &, x(3), cut, c2, rsq, rks, rnow

    ! computes the vectors x(ndim)=(a(1)*n(1),..,a(ndim)*n(ndim))
    ! where n(i) are integers and x(1)**2+..+x(ndim)**2.le.cut**2
    ! the vectors x(i) are stored in rkcomp in
    ! the order given by the values of their norms.
    !  also nshlls gives the number of
    ! different values of the norms ( the relative square norms
    ! differ by less than 1.e-5) and knorm(i) gives
    ! these nshlls norms and kmult(i) is the last vector
    ! whose norm is given by knorm(i). hence the total
    ! number of vectors of magnitude less than cut is kmult(nshlls).
    !

    kmult(0) = 0
    rknorm(0) = 0.d0
    c2 = cut**2
    npts = 1
    do l = 1, ndim
        nkspan(l) = (0.00001 + abs(cut/a(l)))
        !     range of search is (-nkspan(l),+nkspan(l))
        icount(l) = -nkspan(l)
        npts = (2*nkspan(l) + 1)*npts
    end do

    nvects = 0
    do i = 1, npts
        rsq = 0.d0
        ! nzero will be the first nonzero entry in icount
        nzero = 0
        do l = 1, ndim
            if (nzero .eq. 0) nzero = icount(l)
            x(l) = icount(l)*a(l)
            rsq = rsq + x(l)**2
            if (rsq .gt. c2) go to 30
        end do
        !     we only take half of the vectors and exclude the one at the origin
        if (nzero .le. 0) go to 30
        if (nvects .gt. mnkv) then
            write (*, *) ' mnkv too small ', mnkv, nvects
            stop
        end if
        ! we have found a vector
        nvects = nvects + 1

        ! go thru previous vectors. if they have a greater norm move
        ! them up one slot
        do j = 1, nvects - 1
            kj = j
            rks = 0.d0
            do l = 1, ndim
                rks = rks + rkcomp(l, j)**2
            end do
            if (rks .ge. rsq) go to 7
        end do

        kj = nvects

7       do jp = nvects, kj + 1, -1
            do l = 1, ndim
                rkcomp(l, jp) = rkcomp(l, jp - 1)
            end do
        end do

        do l = 1, ndim
            rkcomp(l, kj) = x(l)
        end do

        !
        !
        ! increase counters with carries for the next vector
30      do l = 1, ndim
            icount(l) = icount(l) + 1
            if (icount(l) .le. nkspan(l)) go to 3
            icount(l) = -nkspan(l)
        end do

3   end do

    nshlls = 0
    ! count number of different norms and find pointers
    rnow = 0.d0
    do i = 1, nvects
        rsq = 0.d0
        do l = 1, ndim
            rsq = rsq + rkcomp(l, i)**2
        end do

        if (rsq - rnow .gt. .001*rnow) nshlls = nshlls + 1
        if (nshlls .gt. mnsh) then
            write (*, *) ' mnsh too small ', mnsh, nshlls
            stop
        end if
        rnow = rsq
        rknorm(nshlls) = sqrt(rnow)
        kmult(nshlls) = i
    end do
    !      write(*,*) 'nshlls=',nshlls

end
