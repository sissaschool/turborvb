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

subroutine dsktrs(uplo, n, nhrs, a, b, ldb, x, info)
    implicit none
    integer n, nhrs, ldb, info, i, j
    real*8 a(*), b(ldb, nhrs), x(nhrs, *)
    character uplo
    if (uplo .eq. 'l' .or. uplo .eq. 'L') then
        a(1:n - 1) = -a(1:n - 1)
    end if
    !       Solve   A X = B , where B and X  are  rectangular matrices n x nhrs
    !       A is a non singular skew symmetric matrix with dimension n even
    !       stored in A(i+1,i)=a(i)  i=1,..,n-1
    !       assumed  A(i,i+1)=-a(i)
    !         even x
    !         x(j,2)=b(j,1)/a(1)
    x(:, 2) = b(1, :)/a(1)
    do i = 4, n, 2
        x(:, i) = (b(i - 1, :) + a(i - 2)*x(:, i - 2))/a(i - 1)
    end do
    !         odd x
    !         x(j,n-1)=-b(j,n-2)/a(n-1)
    x(:, n - 1) = -b(n, :)/a(n - 1)
    do i = n - 3, 1, -2
        x(:, i) = -(b(i + 1, :) - a(i + 1)*x(:, i + 2))/a(i)
    end do
    do i = 1, n
        b(i, :) = x(:, i)
    end do
    !       restore value
    if (uplo .eq. 'l' .or. uplo .eq. 'L') then
        a(1:n - 1) = -a(1:n - 1)
    end if
    return
end
