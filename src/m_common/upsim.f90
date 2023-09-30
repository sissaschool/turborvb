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

subroutine upsim(amat, ndim, ind, value, symmagp, ipf)
    use allio, only: nelorb_at, pfaffup, kiontot
    implicit none

    ! argument parameters
    integer, intent(in) :: ndim, ind, ipf
    real*8, intent(in) :: value
    real*8, intent(inout) :: amat(ndim, *)
    logical, intent(in) :: symmagp

    ! local variables
    integer :: i, j, ndimh

    if (ipf .eq. 1) then
        j = (ind - 1)/ndim + 1
        i = ind - ndim*(j - 1)
        amat(i, j) = value
        if (symmagp .and. j .le. ndim) amat(j, i) = value
    else
        ! if symmagp
        !   Pfaff(up,down)=A with A= A^T  and Pfaff(up,up)=Pfaff(do,do)
        !  In the complex case two cases are possible
        !  if yeshermite      --> A=A^+  and Pfaff(up,up)=Pfaff(do,do)^*
        !  if .not.yeshermite --> A=A^T  and Pfaff(up,up)=Pfaff(do,do)
        !  On top of that
        !  if pfaffup=.true. the down-down part of the pfaffian is absent

        ndimh = nelorb_at/2
        j = (ind - 1)/ndim + 1
        i = ind - ndim*(j - 1)
        amat(i, j) = value
        if (j .le. ndim) amat(j, i) = -value ! With the unpaired j>ndim
        if (symmagp .and. kiontot(i) .ne. 0 .and. kiontot(j) .ne. 0) then
            if (j .le. ndimh .and. i .le. ndimh .and. .not. pfaffup) then
                amat(i + ndimh, j + ndimh) = value
                amat(j + ndimh, i + ndimh) = -value
            elseif (j .gt. ndimh .and. j .le. nelorb_at .and. i .le. ndimh) then
                amat(j - ndimh, i + ndimh) = value
                amat(i + ndimh, j - ndimh) = -value
                !       elseif(j.le.ndimh.and.i.gt.ndimh) then
                !          amat(j+ndimh, i-ndimh)=value
                !          amat(i-ndimh, j+ndimh)=-value
            end if
        end if
    end if
    return
end subroutine upsim

!Beccato, è corretto usare nelorb_at
subroutine upsim_complex(amat, ndim, ind, value, symmagp, ipf)
    use allio, only: yes_hermite, nelorb_at, pfaffup, kiontot
    implicit none

    ! argument parameters
    integer, intent(in) :: ndim, ind, ipf
    complex*16, intent(in) :: value(*)
    complex*16, intent(inout) :: amat(ndim, *)
    logical, intent(in) :: symmagp

    ! local variables
    integer i, j, ndimh

    !      write(6,*) value
    if (ipf .eq. 1) then
        j = (ind - 1)/ndim + 1
        i = ind - ndim*(j - 1)
        amat(i, j) = value(1)
        if (symmagp .and. j .le. nelorb_at .and. i .le. nelorb_at) then
            if (yes_hermite) then
                amat(j, i) = conjg(value(1))
            else
                amat(j, i) = value(1)
            end if
        end if
    else
        ndimh = nelorb_at/2
        j = (ind - 1)/ndim + 1
        i = ind - ndim*(j - 1)
        amat(i, j) = value(1)
        if (j .le. ndim) amat(j, i) = -value(1) !   if it is unpaired element j,i is not present

        if (symmagp .and. kiontot(i) .ne. 0 .and. kiontot(j) .ne. 0) then
            if (yes_hermite) then
                if (j .le. ndimh .and. i .le. ndimh .and. .not. pfaffup) then
                    amat(i + ndimh, j + ndimh) = conjg(value(1))
                    amat(j + ndimh, i + ndimh) = -conjg(value(1))
                elseif (j .gt. ndimh .and. i .le. ndimh) then
                    amat(j - ndimh, i + ndimh) = conjg(value(1))
                    amat(i + ndimh, j - ndimh) = -conjg(value(1))
                    !          elseif(j.le.ndimh.and.i.gt.ndimh) then
                    !             amat(j+ndimh, i-ndimh)=conjg(value(1))
                    !             amat(i-ndimh, j+ndimh)=-conjg(value(1))
                end if
            else
                if (j .le. ndimh .and. i .le. ndimh .and. .not. pfaffup) then
                    amat(i + ndimh, j + ndimh) = value(1)
                    amat(j + ndimh, i + ndimh) = -value(1)
                elseif (j .gt. ndimh .and. i .le. ndimh) then
                    amat(j - ndimh, i + ndimh) = value(1)
                    amat(i + ndimh, j - ndimh) = -value(1)
                    !          elseif(j.le.ndimh.and.i.gt.ndimh) then
                    !             amat(j+ndimh, i-ndimh)=value(1)
                    !             amat(i-ndimh, j+ndimh)=-value(1)
                end if
            end if
        end if
    end if
    return
end subroutine upsim_complex
