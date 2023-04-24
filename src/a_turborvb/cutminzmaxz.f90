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

subroutine cutminzmaxz(ipc, contractionj, iesm, npar3bodyr_c, iesuptransbj&
        &, jbraiesm, minz, maxz, indc, alphavar, alphab, scalpar, fixpar, tpar, nmat, rank)
    implicit none
    integer iesuptransbj(*), jbraiesm(*), i, j, iesm, npar3bodyr_c   &
            &, rank, contractionj, indc, nmat, ipc
    real*8 alphavar(*), minz, maxz, cost&
            &, tpar, scalealphab, scaletry, alphab(*), scalpar(*)
    logical fixpar
#ifdef PARALLEL
    include 'mpif.h'
#endif

    !This subroutine rescale alphab if the Z are outside the interval minz< Z < maxz
    !                  Does not modify any other input
    !      if the previous Z (alphavar) are outside the interval --> not changed

    if (maxz .le. minz .or. iesm .eq. 0) return ! input does not make sense
    scalealphab = 1.d0
    if (contractionj .eq. 0) then
        do i = 1, iesm, ipc
            cost = alphavar(indc + i) + tpar*alphab(indc + i)
            if (cost .lt. minz .and. alphavar(indc + i) .ge. minz) then
                scaletry = (minz - alphavar(indc + i))/tpar/alphab(indc + i)/2.d0
                alphab(indc + i) = scaletry*alphab(indc + i)
                if (rank .eq. 0) write (6, *) ' Warning Z too small  ', i, cost, minz
                if (fixpar) then
                    scalpar(indc + i) = 0.d0
                    if (rank .eq. 0) write (6, *) ' Fixing forever this Z ', i, cost
                end if
            elseif (alphavar(indc + i) .lt. minz) then
                if (rank .eq. 0) write (6, *) ' Warning Z smaller than limit  '
                alphavar(indc + i) = minz
                alphab(indc + i) = 0.d0
                if (fixpar) then
                    scalpar(indc + i) = 0.d0
                    if (rank .eq. 0) write (6, *) ' Fixing forever this Z ', i, cost
                end if
            end if
            if (cost .gt. maxz .and. alphavar(indc + i) .le. maxz) then
                if (rank .eq. 0) write (6, *) 'Warning Z too large ', i, cost, maxz
                scaletry = (maxz - alphavar(indc + i))/tpar/alphab(indc + i)/2.d0

                alphab(indc + i) = scaletry*alphab(indc + i)
                if (fixpar) then
                    scalpar(indc + i) = 0.d0
                    if (rank .eq. 0) write (6, *) ' Fixing forever this Z ', i, cost
                end if

                !           if(scaletry.lt.scalealphab) scalealphab=scaletry
            elseif (alphavar(indc + i) .gt. maxz) then
                if (rank .eq. 0) write (6, *) ' Warning Z larger than limit  '
                alphavar(indc + i) = maxz
                alphab(indc + i) = 0.d0
                if (fixpar) then
                    scalpar(indc + i) = 0.d0
                    if (rank .eq. 0) write (6, *) ' Fixing forever this Z ', i, cost
                end if
            end if
        end do
    else
        do i = 1, npar3bodyr_c
            if (iesuptransbj(i) .ne. 0 .and. jbraiesm(i) .ne. 0) then
                j = ipc*(abs(jbraiesm(i)) - 1) + 1
                cost = alphavar(indc + j) + tpar*alphab(indc + j)
                if (cost .lt. minz .and. alphavar(indc + j) .ge. minz) then
                    scaletry = (minz - alphavar(indc + j))/tpar/alphab(indc + j)/2.d0
                    alphab(indc + j) = scaletry*alphab(indc + j)
                    !           if(scaletry.lt.scalealphab) scalealphab=scaletry
                    if (rank .eq. 0) write (6, *) ' Warning Z too small ', j, cost, minz
                    if (fixpar) then
                        scalpar(indc + j) = 0.d0
                        if (rank .eq. 0) write (6, *) ' Fixing forever this Z ', j, cost
                    end if
                elseif (alphavar(indc + j) .lt. minz) then
                    if (rank .eq. 0) write (6, *) ' Warning Z smaller than limit  '
                    alphavar(indc + j) = minz
                    alphab(indc + j) = 0.d0
                    if (fixpar) then
                        scalpar(indc + j) = 0.d0
                        if (rank .eq. 0) write (6, *) ' Fixing forever this Z ', j, cost
                    end if
                end if
                if (cost .gt. maxz .and. alphavar(indc + j) .le. maxz) then
                    scaletry = (maxz - alphavar(indc + j))/tpar/alphab(indc + j)/2.d0
                    !           if(scaletry.lt.scalealphab) scalealphab=scaletry
                    alphab(indc + j) = scaletry*alphab(indc + j)
                    if (rank .eq. 0) write (6, *) ' Warning Z too large  ', j, cost, maxz
                    if (fixpar) then
                        scalpar(indc + j) = 0.d0
                        if (rank .eq. 0) write (6, *) ' Fixing forever this Z ', j, cost
                    end if
                elseif (alphavar(indc + j) .gt. maxz) then
                    if (rank .eq. 0) write (6, *) ' Warning Z larger than limit  '
                    alphavar(indc + j) = maxz
                    alphab(indc + j) = 0.d0
                    if (fixpar) then
                        scalpar(indc + j) = 0.d0
                        if (rank .eq. 0) write (6, *) ' Fixing forever this Z ', j, cost
                    end if
                end if
                !            write(6,*) ' Inside cut i,j =',i,jbraiesm(i)
            end if
        end do
    end if
    !         if(scalealphab.ne.1.d0) call dscal(nmat,0.5d0*scalealphab,alphab,1)
#ifdef PARALLEL
!  To avoid roundoff difference between processors we bcast scalpar,alphab,alphavar
!  safely.
    call bcast_real(alphavar, nmat, 0, MPI_COMM_WORLD)
    call bcast_real(alphab, nmat, 0, MPI_COMM_WORLD)
    call bcast_real(scalpar, nmat, 0, MPI_COMM_WORLD)
#endif
    return
end
