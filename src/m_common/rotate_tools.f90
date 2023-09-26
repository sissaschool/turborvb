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

subroutine ruota_xyz(alpha, xrot, yrot, zrot, nion_1, rion_1, rion_2)
    implicit none
    real(8) u(3, 3)
    real(8), dimension(:, :), allocatable :: emme
    integer nion_1, i, j, ii

    real(8) alpha, xrot, yrot, zrot, rion_1(3, nion_1), rion_2(3, nion_1)

    call dscal(9, 0.d0, u(1, 1), 1)
    call make_u(alpha, xrot, yrot, zrot, u)

    allocate (emme(3*nion_1, 3*nion_1))
    emme = 0.d0

    do i = 1, 3*nion_1, 3
        call caricap(emme, u, i, i, 3*nion_1)
    end do
    call dgemm('N', 'N', 3*nion_1, 1, 3*nion_1, 1.d0, emme, 3*nion_1, &
            &    rion_1, 3*nion_1, 0.d0, rion_2, 3*nion_1)

    do i = 1, nion_1
        write (20, *) (rion_2(j, i), j=1, 3)
    end do

    deallocate (emme)

end subroutine ruota_xyz

subroutine ruota_molec(ipc, alpha, xrot, yrot, zrot, iesupr_1, &
                    &  dupr_1, ioptorb_1, nparam_1, nshell_1, &
                    &  ioptorb, nshell, nelorb, dupr_2)
    use constants, only: zzero, zone
    implicit none
    integer nshell_1, indpar, nelorb, nshell, shift, i, j, ii, iesupr_1, ipc
    integer ioptorb(nshell), ioptorb_1(nshell_1), nparam_1(nshell_1)

    real(8) alpha, xrot, yrot, zrot
    real(8) dupr_1(ipc*iesupr_1)
    real(8) dupr_2(ipc*iesupr_1)

    real(8) u(3, 3), emmed(5, 5), emmef(7, 7), emmeg(9, 9)
    real(8), dimension(:, :), allocatable :: emme

    allocate (emme(nelorb, nelorb))

    call dscal(9, 0.d0, u(1, 1), 1)
    call make_u(alpha, xrot, yrot, zrot, u)
    call prepemmed(u, emmed)
    call prepemmef(u, emmef)
    !          write(*,*)'emmeg'
    call prepemmeg(u, emmeg)

    !          Building matrix lambda_ij

    call build_emme(emme, nelorb, u, emmed, emmef, emmeg, ioptorb, nshell)

    dupr_2 = dupr_1
    indpar = 0
    shift = nelorb + 1

    do i = 1, nshell_1
        if (ioptorb_1(i) .ge. 1000000 .and. nparam_1(i) .eq. 2*nelorb) then
            if (ipc .eq. 1) then
                call dgemv('N', nelorb, nelorb, 1.d0, emme, nelorb, dupr_1(indpar + shift), 1&
                        &, 0.d0, dupr_2(indpar + shift), 1)
            else
                call zgemv('N', nelorb, nelorb, zone, emme, nelorb, dupr_1(2*(indpar + shift) - 1), 2&
                        &, zzero, dupr_2(2*(indpar + shift) - 1), 2)
                call zgemv('N', nelorb, nelorb, zone, emme, nelorb, dupr_1(2*(indpar + shift)), 2&
                        &, zzero, dupr_2(2*(indpar + shift)), 2)
            end if
        end if
        indpar = indpar + nparam_1(i)
    end do
    deallocate (emme)
    return
end subroutine ruota_molec

subroutine ruota_lambda(ipc, ipf, alpha, xrot, yrot, zrot, &
                    &   ix_1, iy_1, detmat_1, nnozero_1, occ_1, nelcol, &
                    &   ioptorb_1, nshell_1, nnozero_2, ix_2, iy_2, detmat_2, &
                    &   symmagp)
    use allio, only: yes_hermite
    use constants, only: zone, zzero
    !         INPUT
    !         alpha,theta,xrot,yrot,zrot
    !         nion_1
    !         occ_1
    !         ixc_1,iyc_1
    !         detmat_1
    !        OUTPUT
    !         ixc_2,iyc_2
    !         detmat_2
    implicit none
    integer ipc, ipf
    integer occ_1, nnozero_1, nnozero_2, iesswr_2
    integer nshell_1, nelcol, neldiff, &
            &            ioptorb_1(nshell_1)

    real(8) alpha, xrot, yrot, zrot
    integer ix_1(nnozero_1), iy_1(nnozero_1)
    real(8) detmat_1(*)
    integer ix_2(*), iy_2(*)
    real(8) detmat_2(*)

    integer i, j, ii
    real(8) u(3, 3), emmed(5, 5), emmef(7, 7), emmeg(9, 9)
    real(8), dimension(:, :), allocatable :: emme, lwork, lambda
    logical symmagp

    allocate (lambda(ipc*occ_1, nelcol), &
            &   emme(ipc*occ_1, occ_1), lwork(ipc*occ_1, nelcol))

    call dscal(9, 0.d0, u(1, 1), 1)
    call make_u(alpha, xrot, yrot, zrot, u)
    call prepemmed(u, emmed)
    call prepemmef(u, emmef)
    !          write(*,*)'emmeg'
    call prepemmeg(u, emmeg)

    !          Building matrix lambda_ij
    !           write(*,*)' Building matrix lambda_ij'
    call dscal(ipc*occ_1*nelcol, 0.d0, lambda, 1)
    do i = 1, nnozero_1
        !          lambda(ipc*(ix_1(i)-1)+1:ipc*ix_1(i),iy_1(i))=&
        !    &detmat_1(ipc*occ_1*(iy_1(i)-1)+ipc*(ix_1(i)-1)+1:ipc*occ_1*(iy_1(i)-1)+ipc*ix_1(i))
        if (ipc .eq. 1) then
            call upsim(lambda, occ_1, (iy_1(i) - 1)*occ_1 + ix_1(i) &
                       , detmat_1(occ_1*(iy_1(i) - 1) + ix_1(i)), symmagp, ipf)
        else
            call upsim_complex(lambda, occ_1, (iy_1(i) - 1)*occ_1 + ix_1(i) &
                               , detmat_1(2*occ_1*(iy_1(i) - 1) + 2*ix_1(i) - 1), symmagp, ipf)
        end if
    end do

    !          New  lambda = emme lambda emme^{T}

    !          write(*,*)' Building emme'
    call build_emme(emme, occ_1, u, emmed, emmef, emmeg, ioptorb_1, nshell_1)

    !           write(*,*)' After Building emme'
    if (ipc .eq. 2) then
        lwork = 0.d0
        ! turn emme to complex
        call dcopy(occ_1*occ_1, emme, 1, lwork, 2)
        call dcopy(2*occ_1*occ_1, lwork, 1, emme, 1)
    end if

    if (ipc .eq. 1) then

        call dgemm('N', 'T', occ_1, occ_1, occ_1, 1.d0, lambda, occ_1, &
                &    emme, occ_1, 0.d0, lwork, occ_1)
        !          write(*,*)' After dgemm'
        call dgemm('N', 'N', occ_1, occ_1, occ_1, 1.d0, emme, occ_1, &
                &    lwork, occ_1, 0.d0, lambda, occ_1)
    else
        call zgemm('N', 'T', occ_1, occ_1, occ_1, zone, lambda, occ_1, &
                &    emme, occ_1, zzero, lwork, occ_1)
        !          write(*,*)' After dgemm'
        call zgemm('N', 'N', occ_1, occ_1, occ_1, zone, emme, occ_1, &
                &    lwork, occ_1, zzero, lambda, occ_1)
    end if

    !          write(*,*)' After dgemm II'

    if (nelcol .gt. occ_1) then
        !         Rotating unpaired orbitals
        neldiff = nelcol - occ_1
        call dcopy(ipc*neldiff*occ_1, lambda(1, occ_1 + 1), 1, lwork, 1)
        if (ipc .eq. 1) then
            call dgemm('N', 'N', occ_1, neldiff, occ_1, 1.d0, emme, occ_1, &
                    &    lwork, occ_1, 0.d0, lambda(1, occ_1 + 1), occ_1)
        else
            call zgemm('N', 'N', occ_1, neldiff, occ_1, zone, emme, occ_1, &
                    &    lwork, occ_1, zzero, lambda(1, occ_1 + 1), occ_1)
        end if
    end if

    do i = 1, nnozero_1
        detmat_2(ipc*occ_1*(iy_1(i) - 1) + ipc*(ix_1(i) - 1) + 1:ipc*occ_1*(iy_1(i) - 1) &
                 + ipc*ix_1(i)) = lambda(ipc*(ix_1(i) - 1) + 1:ipc*ix_1(i), iy_1(i))
    end do
    !           write(*,*)'deallocate'
    deallocate (lambda, emme, lwork)

    !          write(*,*)'esco ruota lambda'
    return
end subroutine ruota_lambda

!#############################################################
!
!##############################################################

subroutine build_emme(emme, occ_1, u, emmed, emmef, emmeg, ioptorb_1, &
        &                        nshell_1)
    implicit none
    integer occ_1, nshell_1, icek, ish, i, j
    real*8 emme(occ_1, occ_1), u(3, 3), emmed(5, 5), emmef(7, 7), &
            &   emmeg(9, 9)
    integer ioptorb_1(nshell_1)

    emme = 0.d0

    ish = 1
    do icek = 1, nshell_1
        !           write(6,*)ish,ioptorb_1(icek)
        selectcase (ioptorb_1(icek))
            ! orbitali s
        case (34)
            emme(ish, ish) = 1.d0
            ish = ish + 1
        case (10)
            emme(ish, ish) = 1.d0
            ish = ish + 1
        case (16)
            emme(ish, ish) = 1.d0
            ish = ish + 1
        case (161, 131)
            emme(ish, ish) = 1.d0
            ish = ish + 1
        case (17)
            emme(ish, ish) = 1.d0
            ish = ish + 1
        case (300:399)
            emme(ish, ish) = 1.d0
            ish = ish + 1
        case (100)
            emme(ish, ish) = 1.d0
            ish = ish + 1
        case (200)
            emme(ish, ish) = 1.d0
            ish = ish + 1
        case (3000:3999)
            emme(ish, ish) = 1.d0
            ish = ish + 1
        case (1000:1099)
            emme(ish, ish) = 1.d0
            ish = ish + 1
        case (900000:1000000)
            emme(ish, ish) = 1.d0
            ish = ish + 1
! orbitali p
        case (20)
            call caricap(emme, u, ish, ish, occ_1)
            ish = ish + 3
        case (22)
            call caricap(emme, u, ish, ish, occ_1)
            ish = ish + 3
        case (36)
            call caricap(emme, u, ish, ish, occ_1)
            ish = ish + 3
        case (400:499)
            call caricap(emme, u, ish, ish, occ_1)
            ish = ish + 3
        case (103)
            call caricap(emme, u, ish, ish, occ_1)
            ish = ish + 3
        case (150)
            call caricap(emme, u, ish, ish, occ_1)
            ish = ish + 3
        case (4000:4999)
            call caricap(emme, u, ish, ish, occ_1)
            ish = ish + 3
        case (1100:1199)
            call caricap(emme, u, ish, ish, occ_1)
            ish = ish + 3
!          orbital d
        case (30:33, 37, 47, 66, 68, 84, 85, 127, 133, 147)
            do i = 0, 4
                do j = 0, 4
                    emme(ish + i, ish + j) = emmed(i + 1, j + 1)
                end do
            end do
            ish = ish + 5
        case (500:599)
            do i = 0, 4
                do j = 0, 4
                    emme(ish + i, ish + j) = emmed(i + 1, j + 1)
                end do
            end do
            ish = ish + 5
        case (5000:5999, 1200:1299, 2200:2299)
            do i = 0, 4
                do j = 0, 4
                    emme(ish + i, ish + j) = emmed(i + 1, j + 1)
                end do
            end do
            ish = ish + 5
!          case f
        case (48, 58, 70, 86, 154)
            do i = 0, 6
                do j = 0, 6
                    emme(ish + i, ish + j) = emmef(i + 1, j + 1)
                end do
            end do
            ish = ish + 7
        case (600)

            do i = 0, 6
                do j = 0, 6
                    emme(ish + i, ish + j) = emmef(i + 1, j + 1)
                end do
            end do
            ish = ish + 7
!           case g

        case (700:701, 51:53)
            do i = 0, 8
                do j = 0, 8
                    emme(ish + i, ish + j) = emmeg(i + 1, j + 1)
                end do
            end do
            ish = ish + 9

        case default
            write (6, *) 'Sorry, rotation not yet defined for orbital '     &
                    &, ioptorb_1(icek)
            stop
        end select

    end do

    if (occ_1 .eq. 2*(ish - 1)) then
    do i = 1, ish - 1
        do j = 1, ish - 1
            emme(ish - 1 + i, ish - 1 + j) = emme(i, j)
        end do
    end do
    end if

    return
end subroutine build_emme

subroutine caricap(a, u, posx, posy, occ)
    implicit none
    integer posx, posy, occ, i, j
    real*8 u(3, 3), a(occ, occ)

    do i = 1, 3
    do j = 1, 3
        a(posx + i - 1, posy + j - 1) = u(i, j)
    end do
    end do

    return
end subroutine caricap

subroutine make_u(alpha, x, y, z, u)
    implicit none
    real*8 alpha, x, y, z, uc
    real*8 u(3, 3)

    uc = 1.d0 - dcos(alpha)

    u(1, 1) = dcos(alpha) + uc*x**2
    u(1, 2) = uc*x*y - dsin(alpha)*z
    u(1, 3) = uc*x*z + dsin(alpha)*y

    u(2, 1) = uc*y*x + dsin(alpha)*z
    u(2, 2) = dcos(alpha) + uc*y**2
    u(2, 3) = uc*y*z - dsin(alpha)*x

    u(3, 1) = uc*z*x - dsin(alpha)*y
    u(3, 2) = uc*z*y + dsin(alpha)*x
    u(3, 3) = dcos(alpha) + uc*z**2

    return
end subroutine make_u

subroutine prepemmed(u, emmed)
    implicit none
    integer i, j
    real*8 mat(3, 3), vec(5), emmed(5, 5), u(3, 3), ut(3, 3), ddot
    ! on input U is the matrix defining the change of coordinates as:
    !   x'_mu = U_{mu,nu} x_nu, x_nu= U_{mu,nu} x'_mu (U is unitary)
    do i = 1, 5
        vec = 0.d0
        vec(i) = 1.d0
        call mat2d(mat, vec, -1)
        call dgemm('N', 'N', 3, 3, 3, 1.d0, u, 3, mat, 3, 0.d0, ut, 3)
        call dgemm('N', 'T', 3, 3, 3, 1.d0, ut, 3, u, 3, 0.d0, mat, 3)
        call mat2d(mat, emmed(1, i), 1)
    end do

    !     write(6,*) ' The matrix '
    !     do i=1,5
    !     write(6,*) (emmed(j,i),j=1,5)
    !     enddo
    !      write(6,*) ' check orthogonality d'
    !      do i=1,5
    !        do j=i,5
    !        write(6,*) i,j,ddot(5,emmed(1,i),1,emmed(1,j),1)
    !        enddo
    !      enddo
    !     stop
    return
end

subroutine prepemmef(u, emmef)
    implicit none
    integer i, j, jj, k
    real*8 mat(3, 3, 3), vec(7), emmef(7, 7), u(3, 3), ut(3, 3), ddot &
        , matblok(3, 3), matscra(3, 3, 3), vecscra(3), vecscra2(3)

    do i = 1, 7
        vec = 0.d0
        matblok = 0.d0
        matscra = 0
        vecscra = 0
        vecscra2 = 0

        vec(i) = 1.d0

        call mat3f(mat, vec, -1)

        do j = 1, 3
            matblok(:, :) = mat(j, :, :)
            call dgemm('N', 'N', 3, 3, 3, 1.d0, u, 3, matblok, 3, 0.d0, ut, 3)
            call dgemm('N', 'T', 3, 3, 3, 1.d0, ut, 3, u, 3, 0.d0, matblok, 3)
            matscra(j, :, :) = matblok(:, :)
        end do

        do j = 1, 3
        do jj = 1, 3
            vecscra(:) = matscra(:, j, jj)
            call dgemv('N', 3, 3, 1.d0, u, 3, vecscra(:), 1, 0.d0, vecscra2, 1)
            mat(:, j, jj) = vecscra2(:)
        end do
        end do

        call mat3f(mat, emmef(1, i), 1)

    end do

    !     write(6,*) ' The matrix '
    !      do i=1,7
    !       write(6,12) (emmef(j,i),j=1,7)
    !      enddo

    !      write(6,*) ' check orthogonality f'
    !      do i=1,7
    !        do j=i,7
!        write(6,*) i,j,ddot(7,emmef(1,i),1,emmef(1,j),1)
    !        enddo
    !      enddo
12  format(7(f10.6, 1x))

    return
end

subroutine prepemmeg(u, emmeg)
    implicit none
    integer i, j, k, jj, i1, i2, i3, i4
    real*8 mat(3, 3, 3, 3), vec(9), emmeg(9, 9), u(3, 3), ut(3, 3), ddot &
        , matblok(3, 3), matscra(3, 3, 3, 3)
    !     real(8) cost,cost1g,cost2g,cost3g,cost4g,cost5g

    !     write(6,*) ' Input matrix '
    !      do i=1,3
    !       do j=1,3
    !       write(6,*) i,j,sum(u(:,i)*u(:,j))
    !       enddo
    !      enddo

    !     cost1g=0.125d0
    !     cost2g=0.79056941504209d0
!     cost3g=0.55901699437494d0
    !     cost4g=2.09165006633518d0
    !     cost5g=0.73950997288745d0

    do i = 1, 9

        matblok = 0.d0
        matscra = 0

        vec = 0.d0
        vec(i) = 1.d0

        !       do j=1,9
        !       vec(j)=dsin(dble(j**2-j))
        !       enddo

        !       write(6,*) ' g element =',i
        !       if(i.le.1) then
        !       cost=1.d0/cost1g
        !       elseif(i.le.3) then
        !       cost=1.d0/cost2g
        !       elseif(i.le.5) then
        !       cost=1.d0/cost3g
        !       elseif(i.le.7) then
        !       cost=1.d0/cost4g
        !       else
        !       cost=1.d0/cost5g
        !       endif

        call mat2g(mat, vec, -1)

        do j = 1, 3
        do k = 1, 3
            matblok(:, :) = mat(j, k, :, :)
            call dgemm('N', 'N', 3, 3, 3, 1.d0, u, 3, matblok, 3, 0.d0, ut, 3)
            call dgemm('N', 'T', 3, 3, 3, 1.d0, ut, 3, u, 3, 0.d0, matblok, 3)
            matscra(j, k, :, :) = matblok(:, :)
        end do
        end do

        do j = 1, 3
        do k = 1, 3
            matblok(:, :) = matscra(:, :, j, k)
            call dgemm('N', 'N', 3, 3, 3, 1.d0, u, 3, matblok, 3, 0.d0, ut, 3)
            call dgemm('N', 'T', 3, 3, 3, 1.d0, ut, 3, u, 3, 0.d0, matblok, 3)
            matscra(:, :, j, k) = matblok(:, :)
        end do
        end do
        !       write(6,*) '  check symmetry ',i
        !       write(6,*) 'I',mat(1,1,1,2),mat(2,1,1,1),mat(1,2,1,1),mat(1,1,2,1)
        !       write(6,*) 'II',mat(1,1,1,3),mat(3,1,1,1),mat(1,3,1,1),mat(1,1,3,1)
        !       write(6,*) 'III',mat(1,1,2,2),mat(1,2,1,2),mat(1,2,2,1),mat(2,1,1,2),mat(2,1,2,1),mat(2,2,1,1)
        !       write(6,*) 'IV',mat(1,1,2,3),mat(1,1,3,2),mat(1,2,1,3),mat(1,3,1,2),mat(1,2,3,1),mat(1,3,2,1)&
        !                      ,mat(2,1,3,1),mat(3,1,2,1),mat(2,3,1,1),mat(3,2,1,1),mat(1,2,3,1),mat(1,3,2,1)

        !      do i1=1,3
        !       do i2=i1,3
        !        do i3=i2,3
        !          do i4=i3,3
        !          if(mat(i1,i2,i3,i4).ne.0.d0) then
        !     write(6,*) ' Non zero element =',i1,i2,i3,i4,mat(i1,i2,i3,i4)*cost
        !          endif
        !          enddo
        !        enddo
        !       enddo
        !      enddo

        call mat2g(matscra, emmeg(1, i), 1)

        !       call mat2g(mat,emmeg(1,i),-1)

        !       write(6,*) ' check matrix ',i,sum(abs(mat(:,:,:,:)-matscra(:,:,:,:)))

    end do

    !    write(6,*) ' The matrix '
    !     do i=1,9
    !      write(6,12) (emmeg(j,i),j=1,9)
    !     enddo
    !     write(6,*) ' check orthogonality g'
    !     do i=1,9
    !       do j=i,9
    !       write(6,*) i,j,ddot(9,emmeg(1,i),1,emmeg(1,j),1)
    !       enddo
    !     enddo
12  format(9(f10.6, 1x))

    return
end

subroutine mat2d(mat, vec, iopt)
    implicit none
    integer iopt
    real*8 mat(3, 3), vec(5), cost1d, cost2d, cost3d
    !     This subroutine maps a symmetric traceless matrix 3x3 to
    !     a 5 dimensional vector of the l=2 representation of the angular
    !     momentum.
    !     iopt=1   mat--> vec, iopt=-1  vec--> mat
    cost1d = 0.5d0
    cost2d = dsqrt(3.d0)/2.d0
    cost3d = dsqrt(3.d0)
    !      The mapping:
    !      Orbitald(#k) = sum_{nu,mu} A^k_{nu,mu} x_nu x_mu
    !      x_nu = x,y,z  for nu=1,2,3 respectively.

    if (iopt .eq. -1) then
        mat(1, 2) = cost3d/2.d0*vec(3)
        mat(2, 1) = mat(1, 2)
        mat(2, 3) = cost3d/2.d0*vec(4)
        mat(3, 2) = mat(2, 3)
        mat(1, 3) = cost3d/2.d0*vec(5)
        mat(3, 1) = mat(1, 3)
        mat(3, 3) = cost1d*2.d0*vec(1)
        mat(2, 2) = (-cost1d*vec(1) + cost2d*vec(2))
        mat(1, 1) = (-cost1d*vec(1) - cost2d*vec(2))
    else
        vec(3) = 2.d0/cost3d*mat(1, 2)
        vec(4) = 2.d0/cost3d*mat(2, 3)
        vec(5) = 2.d0/cost3d*mat(1, 3)
        vec(1) = -(mat(1, 1) + mat(2, 2))/(cost1d*2.d0)
        vec(2) = (mat(2, 2) - mat(1, 1))/(cost2d*2.d0)
    end if

    return
end

subroutine mat3f(mat, vec, iopt)
    implicit none
    real*8 mat(3, 3, 3), vec(7)
    real*8 cost1f, cost2f, cost3f, cost4f
    integer iopt

    cost1f = 0.5d0
    cost2f = 0.612372435695794d0
    cost3f = 1.93649167310371d0
    cost4f = 0.790569415042095d0

    if (iopt .eq. -1) then
        !  symmetric tensor rank 3
        !  the symmetry property remains under rotation

        mat(1, 1, 1) = -cost2f*vec(2) + cost4f*vec(6)
        mat(1, 1, 2) = (-cost2f*vec(3) + 3.d0*cost4f*vec(7))/3.d0
        mat(1, 2, 1) = mat(1, 1, 2)
        mat(2, 1, 1) = mat(1, 1, 2)

        mat(1, 1, 3) = (-3.d0*cost1f*vec(1) + cost3f*vec(4))/3.d0
        mat(1, 3, 1) = mat(1, 1, 3)
        mat(3, 1, 1) = mat(1, 1, 3)

        mat(1, 2, 2) = (-cost2f*vec(2) - 3.d0*cost4f*vec(6))/3.d0
        mat(2, 1, 2) = mat(1, 2, 2)
        mat(2, 2, 1) = mat(1, 2, 2)

        mat(1, 2, 3) = (2.d0*cost3f*vec(5))/6.d0
        mat(2, 3, 1) = mat(1, 2, 3)
        mat(3, 1, 2) = mat(1, 2, 3)
        mat(3, 2, 1) = mat(1, 2, 3)
        mat(2, 1, 3) = mat(1, 2, 3)
        mat(1, 3, 2) = mat(1, 2, 3)

        mat(1, 3, 3) = (4.d0*cost2f*vec(2))/3.d0
        mat(3, 1, 3) = mat(1, 3, 3)
        mat(3, 3, 1) = mat(1, 3, 3)

        mat(2, 2, 2) = (-cost2f*vec(3) - cost4f*vec(7))

        mat(2, 2, 3) = (-3.d0*cost1f*vec(1) - cost3f*vec(4))/3.d0
        mat(3, 2, 2) = mat(2, 2, 3)
        mat(2, 3, 2) = mat(2, 2, 3)

        mat(2, 3, 3) = (4.d0*cost2f*vec(3))/3.d0
        mat(3, 2, 3) = mat(2, 3, 3)
        mat(3, 3, 2) = mat(2, 3, 3)

        mat(3, 3, 3) = 2.d0*cost1f*vec(1)

    else

        vec(1) = mat(3, 3, 3)/2.d0/cost1f

        vec(2) = 3.d0*mat(1, 3, 3)/4.d0/cost2f

        vec(3) = 3.d0*mat(2, 3, 3)/4.d0/cost2f

        vec(4) = 3.d0*(mat(1, 1, 3) - mat(2, 2, 3))/2.d0/cost3f

        vec(5) = mat(1, 2, 3)*3.d0/cost3f

        vec(6) = (mat(1, 1, 1) + cost2f*vec(2))/cost4f

        vec(7) = (-mat(2, 2, 2) - cost2f*vec(3))/cost4f

    end if

    return
end

subroutine mat2g(mat, v, iopt)
    implicit none
    real(8) mat(3, 3, 3, 3), v(9)
    integer iopt, i1, i2, i3, i4, i5
    real(8) cost1g, cost2g, cost3g, cost4g, cost5g

    cost1g = 0.125d0
    cost2g = 0.79056941504209d0
    cost3g = 0.55901699437494d0
    cost4g = 2.09165006633518d0
    cost5g = 0.73950997288745d0

    !     rank 4 symmetric tensor
    !     The symmetry remain conserved under rotations.
    ! Each function can be written as:
    !  sum_{i,j,k,l}  A_i,j,k,l  x_i x_j x_k x_l
    !  Under rotation of x = U x , A transform as a symmetric tensor does.
    !  All indices are the same so the permutation symmetry of indices
    !  is preserved under rotations.

    if (iopt .eq. -1) then

        mat(1, 1, 1, 1) = 3.d0*cost1g*v(1) - cost3g*v(4) + cost5g*v(8)
        !  1

        mat(1, 1, 1, 2) = (-2.d0*cost3g*v(5) + 4.d0*cost5g*v(9))/4.d0
        mat(1, 1, 2, 1) = mat(1, 1, 1, 2)
        mat(1, 2, 1, 1) = mat(1, 1, 1, 2)
        mat(2, 1, 1, 1) = mat(1, 1, 1, 2)

        !  5

        mat(1, 1, 1, 3) = (-3.d0*cost2g*v(2) + cost4g*v(6))/4.d0
        mat(1, 1, 3, 1) = mat(1, 1, 1, 3)
        mat(1, 3, 1, 1) = mat(1, 1, 1, 3)
        mat(3, 1, 1, 1) = mat(1, 1, 1, 3)

        !  9

        mat(1, 1, 2, 2) = (6.d0*cost1g*v(1) - 6.d0*cost5g*v(8))/6.d0
        mat(1, 2, 1, 2) = mat(1, 1, 2, 2)
        mat(1, 2, 2, 1) = mat(1, 1, 2, 2)
        mat(2, 1, 1, 2) = mat(1, 1, 2, 2)
        mat(2, 1, 2, 1) = mat(1, 1, 2, 2)
        mat(2, 2, 1, 1) = mat(1, 1, 2, 2)

        ! 15

        mat(1, 1, 2, 3) = (-3.d0*cost2g*v(3) + 3.d0*cost4g*v(7))/12.d0
        mat(1, 1, 3, 2) = mat(1, 1, 2, 3)
        mat(1, 2, 1, 3) = mat(1, 1, 2, 3)
        mat(1, 3, 1, 2) = mat(1, 1, 2, 3)
        mat(1, 2, 3, 1) = mat(1, 1, 2, 3)
        mat(1, 3, 2, 1) = mat(1, 1, 2, 3)
        mat(2, 1, 3, 1) = mat(1, 1, 2, 3)
        mat(3, 1, 2, 1) = mat(1, 1, 2, 3)
        mat(2, 3, 1, 1) = mat(1, 1, 2, 3)
        mat(3, 2, 1, 1) = mat(1, 1, 2, 3)
        mat(2, 1, 1, 3) = mat(1, 1, 2, 3)
        mat(3, 1, 1, 2) = mat(1, 1, 2, 3)

        ! 27

        mat(1, 1, 3, 3) = (-24.d0*cost1g*v(1) + 6.d0*cost3g*v(4))/6.d0
        mat(1, 3, 1, 3) = mat(1, 1, 3, 3)
        mat(1, 3, 3, 1) = mat(1, 1, 3, 3)
        mat(3, 1, 1, 3) = mat(1, 1, 3, 3)
        mat(3, 1, 3, 1) = mat(1, 1, 3, 3)
        mat(3, 3, 1, 1) = mat(1, 1, 3, 3)

        ! 33

        mat(1, 2, 2, 2) = (-2.d0*cost3g*v(5) - 4*cost5g*v(9))/4.d0
        mat(2, 1, 2, 2) = mat(1, 2, 2, 2)
        mat(2, 2, 1, 2) = mat(1, 2, 2, 2)
        mat(2, 2, 2, 1) = mat(1, 2, 2, 2)

        ! 37

        mat(1, 2, 2, 3) = (-3.d0*cost2g*v(2) - 3.d0*cost4g*v(6))/12.d0
        mat(3, 2, 2, 1) = mat(1, 2, 2, 3)
        mat(2, 1, 2, 3) = mat(1, 2, 2, 3)
        mat(2, 3, 2, 1) = mat(1, 2, 2, 3)
        mat(2, 1, 3, 2) = mat(1, 2, 2, 3)
        mat(2, 3, 1, 2) = mat(1, 2, 2, 3)
        mat(1, 2, 3, 2) = mat(1, 2, 2, 3)
        mat(3, 2, 1, 2) = mat(1, 2, 2, 3)
        mat(1, 3, 2, 2) = mat(1, 2, 2, 3)
        mat(3, 1, 2, 2) = mat(1, 2, 2, 3)
        mat(2, 2, 3, 1) = mat(1, 2, 2, 3)
        mat(2, 2, 1, 3) = mat(1, 2, 2, 3)

        ! 49

        mat(1, 2, 3, 3) = (12.d0*cost3g*v(5))/12.d0
        mat(2, 1, 3, 3) = mat(1, 2, 3, 3)
        mat(3, 3, 2, 1) = mat(1, 2, 3, 3)
        mat(3, 3, 1, 2) = mat(1, 2, 3, 3)
        mat(3, 1, 3, 2) = mat(1, 2, 3, 3)
        mat(3, 2, 3, 1) = mat(1, 2, 3, 3)
        mat(3, 1, 2, 3) = mat(1, 2, 3, 3)
        mat(3, 2, 1, 3) = mat(1, 2, 3, 3)
        mat(1, 3, 2, 3) = mat(1, 2, 3, 3)
        mat(2, 3, 1, 3) = mat(1, 2, 3, 3)
        mat(1, 3, 3, 2) = mat(1, 2, 3, 3)
        mat(2, 3, 3, 1) = mat(1, 2, 3, 3)

        !  61

        mat(1, 3, 3, 3) = (4.d0*cost2g*v(2))/4.d0
        mat(3, 1, 3, 3) = mat(1, 3, 3, 3)
        mat(3, 3, 1, 3) = mat(1, 3, 3, 3)
        mat(3, 3, 3, 1) = mat(1, 3, 3, 3)

        !  65
        mat(2, 2, 2, 2) = 3.d0*cost1g*v(1) + cost3g*v(4) + cost5g*v(8)

        !  66

        mat(2, 2, 2, 3) = (-3.d0*cost2g*v(3) - cost4g*v(7))/4.d0
        mat(3, 2, 2, 2) = mat(2, 2, 2, 3)
        mat(2, 3, 2, 2) = mat(2, 2, 2, 3)
        mat(2, 2, 3, 2) = mat(2, 2, 2, 3)

        !  70

        mat(2, 2, 3, 3) = (-24.d0*cost1g*v(1) - 6.d0*cost3g*v(4))/6.d0
        mat(2, 3, 2, 3) = mat(2, 2, 3, 3)
        mat(2, 3, 3, 2) = mat(2, 2, 3, 3)
        mat(3, 2, 2, 3) = mat(2, 2, 3, 3)
        mat(3, 2, 3, 2) = mat(2, 2, 3, 3)
        mat(3, 3, 2, 2) = mat(2, 2, 3, 3)

        !  76

        mat(2, 3, 3, 3) = (4.d0*cost2g*v(3))/4.d0
        mat(3, 2, 3, 3) = mat(2, 3, 3, 3)
        mat(3, 3, 2, 3) = mat(2, 3, 3, 3)
        mat(3, 3, 3, 2) = mat(2, 3, 3, 3)

        !  80

        mat(3, 3, 3, 3) = 8.d0*cost1g*v(1)

        !  81 --> ALL DEFINED

    else
        !     First the simplest
        v(1) = mat(3, 3, 3, 3)/cost1g/8.d0
        v(2) = mat(1, 3, 3, 3)/cost2g
        v(3) = mat(2, 3, 3, 3)/cost2g
        v(5) = mat(1, 2, 3, 3)/cost3g
        !     get v(4) from :
        !     mat(1,1,1,1)=3.d0*cost1g*v(1)-cost3g*v(4)+cost5g*v(8)
        !     mat(2,2,2,2)=3.d0*cost1g*v(1)+cost3g*v(4)+cost5g*v(8)
        v(4) = (mat(2, 2, 2, 2) - mat(1, 1, 1, 1))/cost3g/2.d0
        !     get v(6) from:
        !     mat(1,1,1,3)=(-3.d0*cost2g*v(2)+cost4g*v(6))/4.d0
        v(6) = (4.d0*mat(1, 1, 1, 3) + 3.d0*cost2g*v(2))/cost4g
        !     get v(7) from:
        !     mat(2,2,2,3)=(-3.d0*cost2g*v(3)-cost4g*v(7))/4.d0
        v(7) = -(4.d0*mat(2, 2, 2, 3) + 3.d0*cost2g*v(3))/cost4g
        !     get v(8) from:
        !     mat(1,1,2,2)=cost1g*v(1)-cost5g*v(8)
        v(8) = -(mat(1, 1, 2, 2) - cost1g*v(1))/cost5g
        !     get v(9) from:
        !     mat(1,1,1,2)=(-2.d0*cost3g*v(5)+4.d0*cost5g*v(9))/4.d0
        v(9) = (4.d0*mat(1, 1, 1, 2) + 2.d0*cost3g*v(5))/4.d0/cost5g
    end if

    return
end

subroutine prep_rotate(ipc, nshell_c, nelorb, nion, ioptorb_c, kion_c, mult_c&
&, nparam_c, rion, zeta, emmep, emme, cellscale, iespbc, apbc)
    implicit none
    integer nshell_c, nelorb, nelorbt, ind, nion, nshell, i, j, k, i3, jj, apbc_sign&
    &, ind_sco, ipc
    integer ioptorb_c(nshell_c), kion_c(nshell_c), mult_c(nshell_c)&
    &, nparam_c(nshell_c)
    real*8 rion(3, nion), zeta(nion), emmep(3, 3), emmed(5, 5), emmef(7, 7)&
&, emmeg(9, 9), emme(ipc*nelorb, nelorb), riont(3), cellscale(3), dist(3)
    real*8 eps, ddot
    logical iespbc, apbc(3), found

    !      integer, dimension(:), allocatable:: ioptorb,ioptorba,cellmap,indsh
    integer, dimension(:), allocatable :: ioptorb, cellmap, indsh
    logical, dimension(:), allocatable :: notfound, notfoundo

    !         emmep is the input rotation assumed from the origin

    eps = 1d-4

    call prepemmed(emmep, emmed)
    call prepemmef(emmep, emmef)
    !          write(*,*)'emmeg'
    call prepemmeg(emmep, emmeg)

    nshell = 0
    do i = 1, nshell_c
    if (nparam_c(i) .le. 1) then
        nshell = nshell + 1
    elseif (ioptorb_c(i) .lt. 900000) then
        nshell = nshell + nparam_c(i)/2
    end if
    end do

    !         allocate(ioptorb(nshell),indsh(nshell),ioptorba(nelorb))
    allocate (ioptorb(nshell), indsh(nshell))

    ind = 0
    nelorbt = 0
    do i = 1, nshell_c
    if (nparam_c(i) .le. 1) then
        ind = ind + 1
        ioptorb(ind) = ioptorb_c(i)
        !         ioptorba(nelorbt+1:nelorbt+mult_c(i))=ioptorb_c(i)
        nelorbt = nelorbt + mult_c(i)
    elseif (ioptorb_c(i) .lt. 900000) then
        do j = 1, nparam_c(i)/2
            ind = ind + 1
            !         ioptorba(nelorbt+1:nelorbt+mult_c(i))=ioptorb_c(i)
            nelorbt = nelorbt + mult_c(i)
            ioptorb(ind) = ioptorb_c(i)
        end do
    end if
    end do

    !         write(6,*) ' nelorbt nelorb found =',nelorb,nelorbt

    allocate (cellmap(nion))

    do i = 1, nion

        call dgemv('N', 3, 3, 1.d0, emmep, 3, rion(1, i), 1, 0.d0, riont, 1)

        found = .false.
        do j = 1, nion
            dist(:) = riont(:) - rion(:, j)
            apbc_sign = 1

            if (iespbc) then

                ! Sign  due to APBC condition on each side
                do i3 = 1, 3
                    if (mod(nint(dist(i3)/cellscale(i3) + eps), 2) .ne. 0 .and. apbc(i3)) then
                        apbc_sign = -apbc_sign
                    end if
                end do

                dist(:) = dist(:) - cellscale(:)*anint(dist(:)/cellscale(:) + eps)
            end if

            if (sum(abs(dist(:))) .lt. eps .and. zeta(i) .eq. zeta(j)) then
                cellmap(i) = j*apbc_sign
                found = .true.
            end if
        end do
        if (.not. found) then
            write (6, *) ' Error symmetry not found !!!', sum(abs(dist(:)))
            stop
        end if

    end do

    !         now define emme
    allocate (notfound(nshell), notfoundo(nshell))

    indsh = 0

    notfound = .true.
    notfoundo = .true.
    ind = 0

    do i = 1, nshell_c
    if (nparam_c(i) .le. 1) then
        ind = ind + 1
        !              write(6,*) ' index updated =',ind
        nelorbt = 0
        ind_sco = 0
        do j = 1, nshell_c
            if (nparam_c(j) .le. 1) then
                ind_sco = ind_sco + 1
                if (kion_c(j) .eq. abs(cellmap(kion_c(i))) .and. ioptorb_c(i) .eq. ioptorb_c(j)&
                &.and. notfound(ind_sco) .and. notfoundo(ind)) then
                    indsh(ind) = nelorbt + 1
                    if (cellmap(kion_c(i)) .lt. 0) indsh(ind) = -indsh(ind)
                    notfound(ind_sco) = .false.
                    notfoundo(ind) = .false.
                end if
                nelorbt = nelorbt + mult_c(j)
            elseif (ioptorb_c(j) .lt. 900000) then
                do k = 1, nparam_c(j)/2
                    ind_sco = ind_sco + 1
                    nelorbt = nelorbt + mult_c(j)
                end do
            end if
        end do
    elseif (ioptorb_c(i) .lt. 900000) then
        do jj = 1, nparam_c(i)/2
            ind = ind + 1
            nelorbt = 0
            ind_sco = 0
            do j = 1, nshell_c
            if (nparam_c(j) .le. 1) then
                ind_sco = ind_sco + 1
                nelorbt = nelorbt + mult_c(j)
            elseif (ioptorb_c(j) .lt. 900000) then
                do k = 1, nparam_c(j)/2
                    ind_sco = ind_sco + 1
                    if (k .eq. jj .and. kion_c(j) .eq. abs(cellmap(kion_c(i))) .and. ioptorb_c(i) .eq.&
                    &ioptorb_c(j) .and. notfound(ind_sco) .and. notfoundo(ind)) then
                        indsh(ind) = nelorbt + 1
                        if (cellmap(kion_c(i)) .lt. 0) indsh(ind) = -indsh(ind)
                        notfound(ind_sco) = .false.
                        notfoundo(ind) = .false.
                    end if
                    nelorbt = nelorbt + mult_c(j)
                end do
            end if
            end do
        end do
    end if
    end do

    write (6, *) ' Index orbitals found ='
    do i = 1, nshell
        write (6, *) i, indsh(i)
    end do

    write (6, *) ' cellmap found '
    do i = 1, nion
        write (6, *) i, cellmap(i)
    end do

    !     stop

    call build_emmel(ipc, emme, nelorbt, emmep, emmed, emmef, emmeg, ioptorb, &
    &                        nshell, indsh)

    !         write(6,*) ' check ortho '

    !         do i=1,nelorbt,23
    !          write(6,*) i,j,ddot(nelorbt,emme(1,i),1,emme(1,i),1)
    !          do j=i+1,nelorbt,23
    !          write(6,*) i,j,ddot(nelorbt,emme(1,i),1,emme(1,j),1)
    !          enddo
    !         enddo

    !         deallocate(cellmap,ioptorb,ioptorba,indsh,notfound,notfoundo)
    deallocate (cellmap, ioptorb, indsh, notfound, notfoundo)

    return
end

subroutine build_emmel(ipc, emme, occ_1, emmep, emmed, emmef, emmeg, ioptorb_1, &
&                        nshell_1, indsh)
    implicit none
    integer occ_1, nshell_1, icek, i, j, ish, ishr, ipc
    real*8 emme(ipc*occ_1, occ_1), emmep(3, 3), emmed(5, 5), emmef(7, 7), &
&   emmeg(9, 9)
    integer ioptorb_1(nshell_1), indsh(nshell_1)
    real*8 sign
    ish = 1
    !  vanishing of emme
    emme = 0.d0

    do icek = 1, nshell_1
        !           write(6,*)ish,ioptorb_1(icek)
        ishr = abs(indsh(icek))
        sign = 1.d0
        if (indsh(icek) .lt. 0) sign = -1.d0

        selectcase (ioptorb_1(icek))
            ! orbitali s
        case (34, 80, 81)
            emme(ishr, ish) = sign
            ish = ish + 1
        case (10)
            emme(ishr, ish) = sign
            ish = ish + 1
        case (16)
            emme(ishr, ish) = sign
            ish = ish + 1
        case (161, 131)
            emme(ishr, ish) = sign
            ish = ish + 1
        case (17)
            emme(ishr, ish) = sign
            ish = ish + 1
        case (300:399)
            emme(ishr, ish) = sign
            ish = ish + 1
        case (100)
            emme(ishr, ish) = sign
            ish = ish + 1
        case (200)
            emme(ishr, ish) = sign
            ish = ish + 1
        case (3000:3999)
            emme(ishr, ish) = sign
            ish = ish + 1
        case (1000:1099)
            emme(ishr, ish) = sign
            ish = ish + 1
        case (900000:1000000)
            emme(ishr, ish) = 1.d0
            ish = ish + 1
            ! orbitali p
        case (20, 22, 50, 36, 82, 83, 400:499, 103, 150, 4000:4999, 1100:1199)
            do i = 0, 2
            do j = 0, 2
                emme(ishr + i, ish + j) = sign*emmep(i + 1, j + 1)
            end do
            end do
            ish = ish + 3

            !          orbital d
        case (30:33, 37, 47, 66, 68, 84, 85, 127, 133, 147)
            do i = 0, 4
            do j = 0, 4
                emme(ishr + i, ish + j) = sign*emmed(i + 1, j + 1)
            end do
            end do
            ish = ish + 5
        case (500:599)
            do i = 0, 4
            do j = 0, 4
                emme(ishr + i, ish + j) = sign*emmed(i + 1, j + 1)
            end do
            end do
            ish = ish + 5
        case (5000:5999, 1200:1299, 2200:2299)
            do i = 0, 4
            do j = 0, 4
                emme(ishr + i, ish + j) = sign*emmed(i + 1, j + 1)
            end do
            end do
            ish = ish + 5
            !          case f
        case (48, 58, 70, 86, 154)
            do i = 0, 6
            do j = 0, 6
                emme(ishr + i, ish + j) = sign*emmef(i + 1, j + 1)
            end do
            end do
            ish = ish + 7
        case (600:699)

            do i = 0, 6
            do j = 0, 6
                emme(ishr + i, ish + j) = sign*emmef(i + 1, j + 1)
            end do
            end do
            ish = ish + 7
            !           case g
        case (700:799, 51:55, 88)
            do i = 0, 8
            do j = 0, 8
                emme(ishr + i, ish + j) = sign*emmeg(i + 1, j + 1)
            end do
            end do
            ish = ish + 9
        case default
            write (6, *) 'Sorry, rotation not yet defined for orbital '     &
    &, ioptorb_1(icek)
            stop
        end select

    end do
    if (occ_1 .eq. 2*(ish - 1)) then
    do i = 1, ish - 1
        do j = 1, ish - 1
            emme(ish - 1 + i, ish - 1 + j) = emme(i, j)
        end do
    end do
    end if
    if (ipc .ne. 1) then
        !   Turn to complex
        do j = 1, occ_1
        do i = occ_1, 1, -1
            emme(2*i - 1, j) = emme(i, j)
            emme(i, j) = 0.d0
        end do
        end do
    end if
    return
end subroutine build_emmel
