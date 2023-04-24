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

subroutine reshuffhub(Lz, Lzr, Ltab, Ltabb, nelnion, nw, np, jbra, kel       &
        &, dist, econf, table, tabler, tabpip, winv, nel2wt, winvj, nel2wtj         &
        &, winvup, nel2upt, winvdo, nel2dot, ainvup, nel2up, wsto                 &
        &, diagfn, diag, psiln, psisn, enert, vpot, vpotreg, enertrue, diffuse      &
        &, nelkel, tmu, naccm, winvbar, nel2bar, winvjbar, winvjbarsz, nel2jbar    &
        &, nel2jbarsz, ivic, pseudolocal, nel, indt, ncore, gradtot, gradtotbar    &
        &, angle, gradpsi, gradpsibar, n_gvec, sum_q_cos_gr, sum_q_sin_gr        &
        &, rank, nproc, ierr, status, psip, skip, iessz, Lbox, psidetln              &
        &, jastrowall_ee, dimee, jastrowall_ei, dimei, indtm, yesivic, vcut, diffkin&
        &, winvfn, nel2wtfn, winvbarfn, nel2barfn, vpotsav_ee, nelsquare)
    use allio, only: ipc, nbra_cyrus, first_cyrus, dim_cyrus, queue_cyrus&
            &, yes_fastbranch, istart
    implicit none

    integer nw, np, nel2up, nel2do, nel2upt, nel2dot, jbra(*), j, kmax2       &
            &, ind, in1, Lz, nel, Ltab, Ltabb, nelkel, nel2wt, nel2wtj, nelnion, naccm(*) &
        &, nel2bar, nel2jbar, nel2jbarsz, indt, ncore, nel3, nel9, nelsquare&
            &, dimee, dimei, n_gvec, tag, tag0, ntags, nel2wtfn, nel2barfn, Lzr
    real*8 econf(*), tabpip(Ltab, *), tmu(Ltabb, *), table(Lz, *), psisn(*)        &
            &, tabler(Lzr, *), diagfn(*)                                          &
            &, wsto(*), diag(*), winv(nel2wt, *), winvj(max(nel2wtj, 1), *)           &
            &, ainvup(nel2up, *), enert(ipc, *), winvup(nel2upt, *), winvdo(nel2dot, *)    &
            &, dist(nelnion, *), vpot(*), diffuse(*), winvjbar(max(nel2jbar, 1), *)   &
            &, winvjbarsz(max(nel2jbarsz, 1), *), enertrue(*)                      &
            &, winvbar(nel2bar, *), ivic(3, max(indt, 1), nel, *), pseudolocal(nel, *)  &
            &, gradtot(*), gradtotbar(*), gradpsi(3*nel, *), gradpsibar(3*nel, *)    &
            &, angle(18*nel, *), psip(*), psiln(*), Lbox, vcut(*), diffkin(3, *)        &
            &, psidetln(*), jastrowall_ee(dimee, *), jastrowall_ei(dimei, *)        &
            &, sum_q_cos_gr(n_gvec + 1, *), sum_q_sin_gr(n_gvec, *), vpotreg(2*nel, *)   &
            &, winvbarfn(nel2barfn, *), winvfn(nel2wtfn, *), vpotsav_ee(nelsquare, *)
    real*8 kel(nelkel, *)
    integer rank, nproc, ierr, skip, nien, nrank, srank, nel2, nelt&
            &, j_proc, ind_proc, indtm(nel, *)
    logical iessz, yesivic

#ifdef PARALLEL
    include 'mpif.h'
    integer status(MPI_STATUS_SIZE), shift
#else
    integer status(*)
#endif

    !#ifdef PARALLEL
    ! total extension of the real*8 elements to be sent for one walker
    !      skip=10+2*Ltab+Lz+nel2upt+nel2dot
    !     &+6*nel+nelkel
    !      if(ncore.gt.0) skip=skip+nel+3*indt*nel+4*nel
    !#endif
    nel2 = 2*nel
    nel3 = 3*nel
    nel9 = 18*nel
    ntags = 20
    nelt = nel3 ! The number of component useful to restore the conf + angle

    in1 = nw/nproc

    !        diagonal part
    do j = 1, nw
        ind = jbra(j)
        if (j .ne. ind) then

#ifdef PARALLEL
            ! find the process (nrank) which j belongs to
            ! and find the process (srank) which ind belongs to
            nrank = (j - 1)/in1
            srank = (ind - 1)/in1

            j_proc = j - nrank*in1
            ind_proc = ind - srank*in1

            ! belonging to the same process
            if (nrank .eq. srank .and. rank .eq. nrank) then
#else
                j_proc = j
                ind_proc = ind
#endif
                if (yes_fastbranch) then
                    call dcopy(nelt, kel(1, ind_proc), 1, kel(1, j_proc), 1)
                else
                    call dcopy(nelkel, kel(1, ind_proc), 1, kel(1, j_proc), 1)
                    diag(j_proc) = diag(ind_proc)
                    diagfn(j_proc) = diagfn(ind_proc)
                    !          enert(j)=enert(ind)
                    enert(:, j_proc) = enert(:, ind_proc)
                    wsto(j_proc) = wsto(ind_proc)
                    psiln(j_proc) = psiln(ind_proc)
                    psidetln(j_proc) = psidetln(ind_proc)
                    psisn(j_proc) = psisn(ind_proc)
                    vpot(j_proc) = vpot(ind_proc)
                    !     vpotreg(j_proc)=vpotreg(ind_proc)
                    call dcopy(nel2, vpotreg(1, ind_proc), 1, vpotreg(1, j_proc), 1)
                    enertrue(j_proc) = enertrue(ind_proc)
                    naccm(j_proc) = naccm(ind_proc)
                    diffuse(j_proc) = diffuse(ind_proc)
                    vcut(j_proc) = vcut(ind_proc)
                    diffkin(:, j_proc) = diffkin(:, ind_proc)
                    gradtot(j_proc) = gradtot(ind_proc)
                    gradtotbar(j_proc) = gradtotbar(ind_proc)
                    indtm(:, j_proc) = indtm(:, ind_proc)
                    call dcopy(Ltab, tabpip(1, ind_proc), 1, tabpip(1, j_proc), 1)
                    call dcopy(Ltabb, tmu(1, ind_proc), 1, tmu(1, j_proc), 1)
                    call dcopy(Lz, table(1, ind_proc), 1, table(1, j_proc), 1)
                    call dcopy(Lzr, tabler(1, ind_proc), 1, tabler(1, j_proc), 1)
                    call dcopy(nel2wt, winv(1, ind_proc), 1, winv(1, j_proc), 1)
                    call dcopy(nel2wtfn, winvfn(1, ind_proc), 1, winvfn(1, j_proc), 1)
                    call dcopy(nelsquare, vpotsav_ee(1, ind_proc), 1, vpotsav_ee(1, j_proc), 1)
                    call dcopy(nel2upt, winvup(1, ind_proc), 1, winvup(1, j_proc), 1)
                    call dcopy(nel2dot, winvdo(1, ind_proc), 1, winvdo(1, j_proc), 1)
                    call dcopy(nel2up, ainvup(1, ind_proc), 1, ainvup(1, j_proc), 1)
                    call dcopy(nelnion, dist(1, ind_proc), 1, dist(1, j_proc), 1)
                    call dcopy(np, econf(ind_proc), in1, econf(j_proc), in1)
                    call dcopy(nel2bar, winvbar(1, ind_proc), 1, winvbar(1, j_proc), 1)
                    call dcopy(nel2barfn, winvbarfn(1, ind_proc), 1, winvbarfn(1, j_proc), 1)
                    call dcopy(nel3, gradpsi(1, ind_proc), 1, gradpsi(1, j_proc), 1)
                    call dcopy(nel3, gradpsibar(1, ind_proc), 1, gradpsibar(1, j_proc), 1)
                    call dcopy(dimee, jastrowall_ee(1, ind_proc)                        &
                            &, 1, jastrowall_ee(1, j_proc), 1)
                    call dcopy(dimei, jastrowall_ei(1, ind_proc)                        &
                            &, 1, jastrowall_ei(1, j_proc), 1)
                    if (LBox .gt. 0.d0) then
                        call dcopy(n_gvec + 1, sum_q_cos_gr(1, ind_proc), 1                      &
                                &, sum_q_cos_gr(1, j_proc), 1)
                        call dcopy(n_gvec, sum_q_sin_gr(1, ind_proc), 1                      &
                                &, sum_q_sin_gr(1, j_proc), 1)
                    end if

                    if (nel2wtj .ne. 0)                                                  &
                            &call dcopy(nel2wtj, winvj(1, ind_proc), 1, winvj(1, j_proc), 1)
                    if (nel2jbar .ne. 0) then
                        call dcopy(nel2jbar, winvjbar(1, ind_proc), 1, winvjbar(1, j_proc), 1)
                        if (iessz) call dcopy(nel2jbar, winvjbarsz(1, ind_proc)              &
                                &, 1, winvjbarsz(1, j_proc), 1)
                    end if
                end if
                if (yesivic) call dcopy(3*indt*nel, ivic(1, 1, 1, ind_proc), 1, ivic(1, 1, 1, j_proc), 1)
                if (ncore .gt. 0) then
                    call dcopy(nel, pseudolocal(1, ind_proc), 1, pseudolocal(1, j_proc), 1)
                end if
                call dcopy(nel9, angle(1, ind_proc), 1, angle(1, j_proc), 1)

                if (nbra_cyrus .gt. 0) then
                    first_cyrus(j_proc) = first_cyrus(ind_proc)
                    call dcopy(dim_cyrus, queue_cyrus(1, 1, 1, ind_proc), 1, queue_cyrus(1, 1, 1, j_proc), 1)
                end if

#ifdef PARALLEL
                ! belonging to another process: sender
            elseif (nrank .ne. srank .and. rank .eq. srank) then

                !      write(6,*) ' before send ',nrank,srank,rank

                shift = 1
                if (yes_fastbranch) then
                    call dcopy(nelt, kel(1, ind_proc), 1&
                         &, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + nelt
                else
                    call dcopy(nelkel, kel(1, ind_proc), 1&
                         &, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + nelkel
                    call dcopy(1, diag(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + 1
                    call dcopy(1, diagfn(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + 1
!          call dcopy(1,enert(ind),1,psip(shift+(ind_proc-1)*skip),1)
                    call dcopy(ipc, enert(1, ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + ipc
                    call dcopy(1, wsto(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + 1
                    call dcopy(1, psiln(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + 1
                    call dcopy(1, psidetln(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + 1
                    call dcopy(1, vpot(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + 1
                    !     call dcopy(1,vpotreg(ind_proc),1,psip(shift+(j_proc-1)*skip),1)
                    !     shift=shift+1
                    call dcopy(nel2, vpotreg(1, ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + nel2
                    call dcopy(1, enertrue(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + 1
                    call dcopy(1, diffuse(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + 1
                    call dcopy(1, vcut(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + 1
                    call dcopy(3, diffkin(1, ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + 3
                    call dcopy(1, gradtot(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + 1
                    call dcopy(1, gradtotbar(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + 1
                    call dcopy(Ltab, tabpip(1, ind_proc), 1                              &
                         &, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + Ltab
                    call dcopy(Ltabb, tmu(1, ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + Ltabb
                    call dcopy(Lz, table(1, ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + Lz
                    call dcopy(Lzr, tabler(1, ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + Lzr
                    call dcopy(nel2upt, winvup(1, ind_proc), 1                           &
                         &, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + nel2upt
                    call dcopy(nel2dot, winvdo(1, ind_proc), 1                           &
                         &, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + nel2dot
                    call dcopy(nel3, gradpsi(1, ind_proc), 1                             &
                         &, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + nel3
                    call dcopy(nel3, gradpsibar(1, ind_proc), 1                          &
                         &, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + nel3

                    if (yesivic) then
                        call dcopy(3*indt*nel, ivic(1, 1, 1, ind_proc), 1          &
                             &, psip(shift + (ind_proc - 1)*skip), 1)
                        shift = shift + 3*indt*nel
                    end if

                    if (ncore .gt. 0) then
                        call dcopy(nel, pseudolocal(1, ind_proc)                            &
                             &, 1, psip(shift + (ind_proc - 1)*skip), 1)
                        shift = shift + nel
                    end if
                end if
                call dcopy(nel9, angle(1, ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + nel9

                if (.not. yes_fastbranch) then
                    call dcopy(np, econf(ind_proc), in1, psip(shift + (ind_proc - 1)*skip), 1)
                    if (shift + np - 1 .ne. skip) then
                        write (6, *) ' Error in dimensioning skipreshuff  in main !!! ',    &
                             &rank, shift + np - 1, skip
                        call mpi_finalize(ierr)
                        stop
                    end if

                    !      write(6,*) ' before send ',nrank,srank
                    tag0 = ntags*(j_proc - 1)
                    tag = tag0 + 1
                    call mpi_send(psip(1 + (ind_proc - 1)*skip), skip                        &
                         &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                    tag = tag0 + 2
                    call mpi_send(naccm(ind_proc), 1, MPI_INTEGER                       &
                         &, nrank, tag, MPI_COMM_WORLD, ierr)
                    tag = tag0 + 3
                    call mpi_send(psisn(ind_proc), 1, MPI_DOUBLE_PRECISION             &
                         &, nrank, tag, MPI_COMM_WORLD, ierr)
                    tag = tag0 + 4
                    call mpi_send(indtm(1, ind_proc), nel, MPI_INTEGER                   &
                         &, nrank, tag, MPI_COMM_WORLD, ierr)
                    !     now sending each big vectors independently
                    tag = tag0 + 5
                    call mpi_send(winv(1, ind_proc), nel2wt                             &
                         &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                    tag = tag0 + 6
                    call mpi_send(winvfn(1, ind_proc), nel2wtfn                         &
                         &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                    tag = tag0 + 7
                    call mpi_send(ainvup(1, ind_proc), nel2up                           &
                         &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                    tag = tag0 + 8
                    call mpi_send(dist(1, ind_proc), nelnion                            &
                         &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                    tag = tag0 + 9
                    call mpi_send(winvbar(1, ind_proc), nel2bar                         &
                         &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                    tag = tag0 + 10
                    call mpi_send(winvbarfn(1, ind_proc), nel2barfn                         &
                         &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                    tag = tag0 + 11
                    call mpi_send(jastrowall_ee(1, ind_proc), dimee                     &
                         &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                    tag = tag0 + 12
                    call mpi_send(jastrowall_ei(1, ind_proc), dimei                     &
                         &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                    tag = tag0 + 13
                    call mpi_send(vpotsav_ee(1, ind_proc), nelsquare&
                         &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                    if (nel2wtj .ne. 0) then
                        tag = tag0 + 14
                        call mpi_send(winvj(1, ind_proc), nel2wtj&
                             &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                    end if
                    if (nel2jbar .ne. 0) then
                        tag = tag0 + 15
                        call mpi_send(winvjbar(1, ind_proc), nel2jbar                       &
                             &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                        if (iessz) then
                            tag = tag0 + 16
                            call mpi_send(winvjbarsz(1, ind_proc), nel2jbar                     &
                                 &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                        end if
                    end if
                    if (Lbox .gt. 0) then
                        tag = tag0 + 17
                        call mpi_send(sum_q_cos_gr(1, ind_proc), n_gvec + 1                     &
                             &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                        tag = tag0 + 18
                        call mpi_send(sum_q_sin_gr(1, ind_proc), n_gvec                     &
                             &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                    end if
                else
                    if (shift - 1 .ne. skip) then
                        write (6, *) ' Error in dimensioning skipreshuff  in main !!! ',    &
                             &rank, shift - 1, skip
                        call mpi_finalize(ierr)
                        stop
                    end if
                    tag0 = ntags*(j_proc - 1)
                    tag = tag0 + 1
                    call mpi_send(psip(1 + (ind_proc - 1)*skip), skip                        &
                         &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                end if
                if (nbra_cyrus .gt. 0) then
                    tag = tag0 + 19
                    call mpi_send(first_cyrus(ind_proc), 1, MPI_INTEGER                       &
                         &, nrank, tag, MPI_COMM_WORLD, ierr)
                    tag = tag0 + 20
                    call mpi_send(queue_cyrus(1, 1, 1, ind_proc), dim_cyrus                     &
                    &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, ierr)
                end if
                ! belonging to another process: receiver
            elseif (nrank .ne. srank .and. rank .eq. nrank) then

                !      write(6,*) ' before receive  ',nrank,srank,rank
                tag0 = ntags*(j_proc - 1)
                tag = tag0 + 1

                call mpi_recv(psip(1 + skip*(j_proc - 1)), skip                        &
                     &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)

                if (yes_fastbranch) then
                    shift = 1
                    call dcopy(nelt, psip(shift + (j_proc - 1)*skip), 1, kel(1, j_proc), 1)
                    shift = shift + nelt
                else
                    shift = 1
                    call dcopy(nelkel, psip(shift + (j_proc - 1)*skip), 1, kel(1, j_proc), 1)
                    shift = shift + nelkel
                    call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, diag(j_proc), 1)
                    shift = shift + 1
                    call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, diagfn(j_proc), 1)
                    shift = shift + 1
!          call dcopy(1,psip(shift+skip*(j_proc-1)),1,enert(j),1)
                    call dcopy(ipc, psip(shift + skip*(j_proc - 1)), 1, enert(1, j_proc), 1)
                    shift = shift + ipc
                    call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, wsto(j_proc), 1)
                    shift = shift + 1
                    call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, psiln(j_proc), 1)
                    shift = shift + 1
                    call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, psidetln(j_proc), 1)
                    shift = shift + 1
                    call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, vpot(j_proc), 1)
                    shift = shift + 1
                    !     call dcopy(1,psip(shift+skip*(j_proc-1)),1,vpotreg(j_proc),1)
                    !     shift=shift+1
                    call dcopy(nel2, psip(shift + skip*(j_proc - 1)), 1, vpotreg(1, j_proc), 1)
                    shift = shift + nel2
                    call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, enertrue(j_proc), 1)
                    shift = shift + 1
                    call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, diffuse(j_proc), 1)
                    shift = shift + 1
                    call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, vcut(j_proc), 1)
                    shift = shift + 1
                    call dcopy(3, psip(shift + skip*(j_proc - 1)), 1, diffkin(1, j_proc), 1)
                    shift = shift + 3
                    call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, gradtot(j_proc), 1)
                    shift = shift + 1
                    call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, gradtotbar(j_proc), 1)
                    shift = shift + 1
                    call dcopy(Ltab, psip(shift + skip*(j_proc - 1)), 1, tabpip(1, j_proc), 1)
                    shift = shift + Ltab
                    call dcopy(Ltabb, psip(shift + skip*(j_proc - 1)), 1, tmu(1, j_proc), 1)
                    shift = shift + Ltabb
                    call dcopy(Lz, psip(shift + skip*(j_proc - 1)), 1, table(1, j_proc), 1)
                    shift = shift + Lz
                    call dcopy(Lzr, psip(shift + skip*(j_proc - 1)), 1, tabler(1, j_proc), 1)
                    shift = shift + Lzr
                    call dcopy(nel2upt, psip(shift + skip*(j_proc - 1))                    &
                         &, 1, winvup(1, j_proc), 1)
                    shift = shift + nel2upt
                    call dcopy(nel2dot, psip(shift + skip*(j_proc - 1))                    &
                         &, 1, winvdo(1, j_proc), 1)
                    shift = shift + nel2dot
                    call dcopy(nel3, psip(shift + (j_proc - 1)*skip), 1, gradpsi(1, j_proc), 1)
                    shift = shift + nel3
                    call dcopy(nel3, psip(shift + (j_proc - 1)*skip)                       &
                         &, 1, gradpsibar(1, j_proc), 1)
                    shift = shift + nel3

                    if (yesivic) then
                        call dcopy(3*indt*nel, psip(shift + skip*(j_proc - 1))                 &
                             &, 1, ivic(1, 1, 1, j_proc), 1)
                        shift = shift + 3*indt*nel
                    end if

                    if (ncore .gt. 0) then
                        call dcopy(nel, psip(shift + skip*(j_proc - 1))                        &
                             &, 1, pseudolocal(1, j_proc), 1)
                        shift = shift + nel
                    end if
                end if
                call dcopy(nel9, psip(shift + (j_proc - 1)*skip), 1, angle(1, j_proc), 1)
                shift = shift + nel9

                if (.not. yes_fastbranch) then

                    call dcopy(np, psip(shift + (j_proc - 1)*skip), 1, econf(j_proc), in1)

                    !      write(6,*) ' before receive  ',nrank,srank

                    tag = tag0 + 2
                    call mpi_recv(naccm(j_proc), 1, MPI_INTEGER                         &
                         &, srank, tag, MPI_COMM_WORLD, status, ierr)
                    tag = tag0 + 3
                    call mpi_recv(psisn(j_proc), 1, MPI_DOUBLE_PRECISION               &
                         &, srank, tag, MPI_COMM_WORLD, status, ierr)
                    tag = tag0 + 4
                    call mpi_recv(indtm(1, j_proc), nel, MPI_INTEGER                     &
                         &, srank, tag, MPI_COMM_WORLD, status, ierr)
                    !     now receiving each big vectors independently
                    tag = tag0 + 5
                    call mpi_recv(winv(1, j_proc), nel2wt                               &
                         &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)
                    tag = tag0 + 6
                    call mpi_recv(winvfn(1, j_proc), nel2wtfn                           &
                         &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)
                    tag = tag0 + 7
                    call mpi_recv(ainvup(1, j_proc), nel2up                             &
                         &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)
                    tag = tag0 + 8
                    call mpi_recv(dist(1, j_proc), nelnion                              &
                         &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)
                    tag = tag0 + 9
                    call mpi_recv(winvbar(1, j_proc), nel2bar                           &
                         &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)
                    tag = tag0 + 10
                    call mpi_recv(winvbarfn(1, j_proc), nel2barfn                       &
                         &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)
                    tag = tag0 + 11
                    call mpi_recv(jastrowall_ee(1, j_proc), dimee                       &
                         &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)
                    tag = tag0 + 12
                    call mpi_recv(jastrowall_ei(1, j_proc), dimei                       &
                         &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)
                    tag = tag0 + 13
                    call mpi_recv(vpotsav_ee(1, j_proc), nelsquare&
                         &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)
                    if (nel2wtj .ne. 0) then
                        tag = tag0 + 14
                        call mpi_recv(winvj(1, j_proc), nel2wtj                             &
                             &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)
                    end if
                    if (nel2jbar .ne. 0) then
                        tag = tag0 + 15
                        call mpi_recv(winvjbar(1, j_proc), nel2jbar                         &
                             &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)
                        if (iessz) then
                            tag = tag0 + 16
                            call mpi_recv(winvjbarsz(1, j_proc), nel2jbar                       &
                                 &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)
                        end if
                    end if
                    if (Lbox .gt. 0) then
                        tag = tag0 + 17
                        call mpi_recv(sum_q_cos_gr(1, j_proc), n_gvec + 1                       &
                             &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)
                        tag = tag0 + 18
                        call mpi_recv(sum_q_sin_gr(1, j_proc), n_gvec                       &
                             &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)
                    end if
                end if
                if (nbra_cyrus .gt. 0) then
                    tag = tag0 + 19
                    call mpi_recv(first_cyrus(j_proc), 1, MPI_INTEGER                         &
                         &, srank, tag, MPI_COMM_WORLD, status, ierr)
                    tag = tag0 + 20
                    call mpi_recv(queue_cyrus(1, 1, 1, j_proc), dim_cyrus                         &
                            &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, status, ierr)
                end if
                !      write(6,*) ' after  receive  ',nrank,srank

            end if

#endif

        end if

    end do

    return
end subroutine reshuffhub

subroutine reshuffhub_noblock(Lz, Lzr, Ltab, Ltabb, nelnion, nw, np, jbra, kel       &
        &, dist, econf, table, tabler, tabpip, winv, nel2wt, winvj, nel2wtj         &
        &, winvup, nel2upt, winvdo, nel2dot, ainvup, nel2up, wsto                 &
        &, diagfn, diag, psiln, psisn, enert, vpot, vpotreg, enertrue, diffuse      &
        &, nelkel, tmu, naccm, winvbar, nel2bar, winvjbar, winvjbarsz, nel2jbar    &
        &, nel2jbarsz, ivic, pseudolocal, nel, indt, ncore, gradtot, gradtotbar    &
        &, angle, gradpsi, gradpsibar, n_gvec, sum_q_cos_gr, sum_q_sin_gr        &
        &, rank, nproc, ierr, status, psip, skip, iessz, Lbox, psidetln              &
        &, jastrowall_ee, dimee, jastrowall_ei, dimei, indtm, yesivic, vcut, diffkin&
        &, winvfn, nel2wtfn, winvbarfn, nel2barfn, vpotsav_ee, nelsquare)
    use allio, only: ipc, nbra_cyrus, first_cyrus, dim_cyrus, queue_cyrus, istart
    implicit none

    integer nw, np, nel2up, nel2do, nel2upt, nel2dot, jbra(*), j, k, kmax2       &
            &, ind, in1, Lz, nel, Ltab, Ltabb, nelkel, nel2wt, nel2wtj, nelnion, naccm(*) &
            &, nel2bar, nel2jbar, nel2jbarsz, indt, ncore, nel3, nel9, nrec        &
            &, dimee, dimei, n_gvec, tag, tag0, ntags, nel2wtfn, nel2barfn&
            &, Lzr, icount, nelsquare
    real*8 econf(*), tabpip(Ltab, *), tmu(Ltabb, *), table(Lz, *), psisn(*)        &
            &, tabler(Lzr, *), diagfn(*)                                          &
            &, wsto(*), diag(*), winv(nel2wt, *), winvj(max(nel2wtj, 1), *)           &
            &, ainvup(nel2up, *), enert(ipc, *), winvup(nel2upt, *), winvdo(nel2dot, *)    &
            &, dist(nelnion, *), vpot(*), diffuse(*), winvjbar(max(nel2jbar, 1), *)   &
            &, winvjbarsz(max(nel2jbarsz, 1), *), enertrue(*)                      &
            &, winvbar(nel2bar, *), ivic(3, max(indt, 1), nel, *), pseudolocal(nel, *)  &
            &, gradtot(*), gradtotbar(*), gradpsi(3*nel, *), gradpsibar(3*nel, *)    &
            &, angle(18*nel, *), psip(*), psiln(*), Lbox, vcut(*), diffkin(3, *)        &
            &, psidetln(*), jastrowall_ee(dimee, *), jastrowall_ei(dimei, *)        &
            &, sum_q_cos_gr(n_gvec, *), sum_q_sin_gr(n_gvec, *), vpotreg(2*nel, *)   &
            &, winvbarfn(nel2barfn, *), winvfn(nel2wtfn, *), vpotsav_ee(nelsquare, *)
    real*8 kel(nelkel, *)
    integer rank, nproc, ierr, skip, nien, nrank, srank, nel2&
            &, j_proc, ind_proc, indtm(nel, *)
    logical iessz, yesivic

#ifdef PARALLEL
    include 'mpif.h'
    integer, allocatable :: req(:), jrec(:), statusi(:, :)
    integer shift, status(MPI_STATUS_SIZE)
    allocate (req(20*nw), jrec(nw), statusi(MPI_STATUS_SIZE, 20*nw))
    req = 0
    jrec = 0
#else
    integer status(*)
#endif

    !#ifdef PARALLEL
    ! total extension of the real*8 elements to be sent for one walker
    !      skip=10+2*Ltab+Lz+nel2upt+nel2dot
    !     &+6*nel+nelkel
    !      if(ncore.gt.0) skip=skip+nel+3*indt*nel+4*nel
    !#endif

    nrec = 0
    nel2 = 2*nel
    nel3 = 3*nel
    nel9 = 18*nel
    ntags = 20
    in1 = nw/nproc

    icount = 0
    do j = 1, nw
        ind = jbra(j)
        if (j .ne. ind) then

#ifdef PARALLEL
            ! find the process (nrank) which j belongs to
            ! and find the process (srank) which ind belongs to
            nrank = (j - 1)/in1
            srank = (ind - 1)/in1

            j_proc = j - nrank*in1
            ind_proc = ind - srank*in1

            ! belonging to the same process
            if (nrank .eq. srank .and. rank .eq. nrank) then
#else
                j_proc = j
                ind_proc = ind
#endif

                diag(j_proc) = diag(ind_proc)
                diagfn(j_proc) = diagfn(ind_proc)
                !          enert(j)=enert(ind)
                enert(:, j_proc) = enert(:, ind_proc)
                wsto(j_proc) = wsto(ind_proc)
                psiln(j_proc) = psiln(ind_proc)
                psidetln(j_proc) = psidetln(ind_proc)
                psisn(j_proc) = psisn(ind_proc)
                vpot(j_proc) = vpot(ind_proc)
                !     vpotreg(j_proc)=vpotreg(ind_proc)
                call dcopy(nel2, vpotreg(1, ind_proc), 1, vpotreg(1, j_proc), 1)
                enertrue(j_proc) = enertrue(ind_proc)
                naccm(j_proc) = naccm(ind_proc)
                diffuse(j_proc) = diffuse(ind_proc)
                vcut(j_proc) = vcut(ind_proc)
                diffkin(:, j_proc) = diffkin(:, ind_proc)
                gradtot(j_proc) = gradtot(ind_proc)
                gradtotbar(j_proc) = gradtotbar(ind_proc)
                indtm(:, j_proc) = indtm(:, ind_proc)
                call dcopy(Ltab, tabpip(1, ind_proc), 1, tabpip(1, j_proc), 1)
                call dcopy(Ltabb, tmu(1, ind_proc), 1, tmu(1, j_proc), 1)
                call dcopy(Lz, table(1, ind_proc), 1, table(1, j_proc), 1)
                call dcopy(Lzr, tabler(1, ind_proc), 1, tabler(1, j_proc), 1)
                call dcopy(nel2wt, winv(1, ind_proc), 1, winv(1, j_proc), 1)
                call dcopy(nel2wtfn, winvfn(1, ind_proc), 1, winvfn(1, j_proc), 1)
                call dcopy(nelsquare, vpotsav_ee(1, ind_proc), 1, vpotsav_ee(1, j_proc), 1)
                call dcopy(nel2upt, winvup(1, ind_proc), 1, winvup(1, j_proc), 1)
                call dcopy(nel2dot, winvdo(1, ind_proc), 1, winvdo(1, j_proc), 1)
                call dcopy(nel2up, ainvup(1, ind_proc), 1, ainvup(1, j_proc), 1)
                call dcopy(nelnion, dist(1, ind_proc), 1, dist(1, j_proc), 1)
                call dcopy(np, econf(ind_proc), in1, econf(j_proc), in1)
                call dcopy(nel2bar, winvbar(1, ind_proc), 1, winvbar(1, j_proc), 1)
                call dcopy(nel2barfn, winvbarfn(1, ind_proc), 1, winvbarfn(1, j_proc), 1)
                call dcopy(nel3, gradpsi(1, ind_proc), 1, gradpsi(1, j_proc), 1)
                call dcopy(nel3, gradpsibar(1, ind_proc), 1, gradpsibar(1, j_proc), 1)
                call dcopy(nelkel, kel(1, ind_proc), 1, kel(1, j_proc), 1)
                call dcopy(dimee, jastrowall_ee(1, ind_proc)                        &
                        &, 1, jastrowall_ee(1, j_proc), 1)
                call dcopy(dimei, jastrowall_ei(1, ind_proc)                        &
                        &, 1, jastrowall_ei(1, j_proc), 1)
                if (LBox .gt. 0.d0) then
                    call dcopy(n_gvec + 1, sum_q_cos_gr(1, ind_proc), 1                      &
                            &, sum_q_cos_gr(1, j_proc), 1)
                    call dcopy(n_gvec, sum_q_sin_gr(1, ind_proc), 1                      &
                            &, sum_q_sin_gr(1, j_proc), 1)
                end if

                if (nel2wtj .ne. 0)                                                  &
                        &call dcopy(nel2wtj, winvj(1, ind_proc), 1, winvj(1, j_proc), 1)
                if (nel2jbar .ne. 0) then
                    call dcopy(nel2jbar, winvjbar(1, ind_proc), 1, winvjbar(1, j_proc), 1)
                    if (iessz) call dcopy(nel2jbar, winvjbarsz(1, ind_proc)              &
                            &, 1, winvjbarsz(1, j_proc), 1)
                end if
                if (yesivic) call dcopy(3*indt*nel, ivic(1, 1, 1, ind_proc), 1, ivic(1, 1, 1, j_proc), 1)
                if (ncore .gt. 0) then
                    call dcopy(nel, pseudolocal(1, ind_proc), 1, pseudolocal(1, j_proc), 1)
                end if
                call dcopy(nel9, angle(1, ind_proc), 1, angle(1, j_proc), 1)

                if (nbra_cyrus .gt. 0) then
                    first_cyrus(j_proc) = first_cyrus(ind_proc)
                    call dcopy(dim_cyrus, queue_cyrus(1, 1, 1, ind_proc), 1, queue_cyrus(1, 1, 1, j_proc), 1)
                end if

#ifdef PARALLEL
                ! belonging to another process: sender
            elseif (nrank .ne. srank .and. rank .eq. srank) then

                !      write(6,*) ' before send ',nrank,srank,rank

                shift = 1
                call dcopy(1, diag(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + 1
                call dcopy(1, diagfn(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + 1
!          call dcopy(1,enert(ind),1,psip(shift+(ind_proc-1)*skip),1)
                call dcopy(ipc, enert(1, ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + ipc
                call dcopy(1, wsto(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + 1
                call dcopy(1, psiln(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + 1
                call dcopy(1, psidetln(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + 1
                call dcopy(1, vpot(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + 1
                !     call dcopy(1,vpotreg(ind_proc),1,psip(shift+(j_proc-1)*skip),1)
                !     shift=shift+1
                call dcopy(nel2, vpotreg(1, ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + nel2
                call dcopy(1, enertrue(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + 1
                call dcopy(1, diffuse(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + 1
                call dcopy(1, vcut(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + 1
                call dcopy(3, diffkin(1, ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + 3
                call dcopy(1, gradtot(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + 1
                call dcopy(1, gradtotbar(ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + 1
                call dcopy(Ltab, tabpip(1, ind_proc), 1                              &
                     &, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + Ltab
                call dcopy(Ltabb, tmu(1, ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + Ltabb
                call dcopy(Lz, table(1, ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + Lz
                call dcopy(Lzr, tabler(1, ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + Lzr
                call dcopy(nel2upt, winvup(1, ind_proc), 1                           &
                     &, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + nel2upt
                call dcopy(nel2dot, winvdo(1, ind_proc), 1                           &
                     &, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + nel2dot
                call dcopy(nel3, gradpsi(1, ind_proc), 1                             &
                     &, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + nel3
                call dcopy(nel3, gradpsibar(1, ind_proc), 1                          &
                     &, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + nel3
                call dcopy(nelkel, kel(1, ind_proc), 1                               &
                     &, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + nelkel

                if (yesivic) then
                    call dcopy(3*indt*nel, ivic(1, 1, 1, ind_proc), 1          &
                         &, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + 3*indt*nel
                end if

                if (ncore .gt. 0) then
                    call dcopy(nel, pseudolocal(1, ind_proc)                            &
                         &, 1, psip(shift + (ind_proc - 1)*skip), 1)
                    shift = shift + nel
                end if
                call dcopy(nel9, angle(1, ind_proc), 1, psip(shift + (ind_proc - 1)*skip), 1)
                shift = shift + nel9
                call dcopy(np, econf(ind_proc), in1, psip(shift + (ind_proc - 1)*skip), 1)

                if (shift + np - 1 .ne. skip) then
                    write (6, *) ' Error in dimensioning skipreshuff  in main !!! ',    &
                         &rank, shift + np - 1, skip
                    call mpi_finalize(ierr)
                    stop
                end if

                !      write(6,*) ' before send ',nrank,srank
                icount = icount + 1
                tag0 = ntags*(j_proc - 1)
                tag = tag0 + 1
                call mpi_isend(psip(1 + (ind_proc - 1)*skip), skip                    &
              &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                icount = icount + 1
                tag = tag0 + 2
                call mpi_isend(naccm(ind_proc), 1, MPI_INTEGER                       &
                      &, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                icount = icount + 1
                tag = tag0 + 3
                call mpi_isend(psisn(ind_proc), 1, MPI_DOUBLE_PRECISION             &
                     &, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                icount = icount + 1
                tag = tag0 + 4
                call mpi_isend(indtm(1, ind_proc), nel, MPI_INTEGER                   &
                     &, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                !     now sending each big vectors independently
                icount = icount + 1
                tag = tag0 + 5
                call mpi_isend(winv(1, ind_proc), nel2wt                             &
                     &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                icount = icount + 1
                tag = tag0 + 6
                call mpi_isend(winvfn(1, ind_proc), nel2wtfn                         &
                &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                icount = icount + 1
                tag = tag0 + 7
                call mpi_isend(ainvup(1, ind_proc), nel2up                           &
                     &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                icount = icount + 1
                tag = tag0 + 8
                call mpi_isend(dist(1, ind_proc), nelnion                            &
                     &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                icount = icount + 1
                tag = tag0 + 9
                call mpi_isend(winvbar(1, ind_proc), nel2bar                         &
                     &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                icount = icount + 1
                tag = tag0 + 10
                call mpi_isend(winvbarfn(1, ind_proc), nel2barfn                         &
                     &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                icount = icount + 1
                tag = tag0 + 11
                call mpi_isend(jastrowall_ee(1, ind_proc), dimee                     &
                 &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                icount = icount + 1
                tag = tag0 + 12
                call mpi_isend(jastrowall_ei(1, ind_proc), dimei                     &
              &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                icount = icount + 1
                tag = tag0 + 13
                call mpi_isend(vpotsav_ee(1, ind_proc), nelsquare&
                &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                if (nel2wtj .ne. 0) then
                    icount = icount + 1
                    tag = tag0 + 14
                    call mpi_isend(winvj(1, ind_proc), nel2wtj                           &
                         &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                end if
                if (nel2jbar .ne. 0) then
                    tag = tag0 + 15
                    icount = icount + 1
                    call mpi_isend(winvjbar(1, ind_proc), nel2jbar                       &
                         &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                    if (iessz) then
                        tag = tag0 + 16
                        icount = icount + 1
                        call mpi_isend(winvjbarsz(1, ind_proc), nel2jbar                     &
                             &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                    end if
                end if
                if (Lbox .gt. 0) then
                    tag = tag0 + 17
                    icount = icount + 1
                    call mpi_isend(sum_q_cos_gr(1, ind_proc), n_gvec + 1                     &
                         &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                    tag = tag0 + 18
                    icount = icount + 1
                    call mpi_isend(sum_q_sin_gr(1, ind_proc), n_gvec                     &
                         &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                end if

                if (nbra_cyrus .gt. 0) then
                    tag = tag0 + 19
                    icount = icount + 1
                    call mpi_isend(first_cyrus(ind_proc), 1, MPI_INTEGER                       &
                         &, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                    tag = tag0 + 20
                    icount = icount + 1
                    call mpi_send(queue_cyrus(1, 1, 1, ind_proc), dim_cyrus                     &
                    &, MPI_DOUBLE_PRECISION, nrank, tag, MPI_COMM_WORLD, req(icount), ierr)
                end if

                !      write(6,*) ' after  send ',nrank,srank

                ! belonging to another process: receiver
            elseif (nrank .ne. srank .and. rank .eq. nrank) then

                nrec = nrec + 1
                jrec(nrec) = j
                !      write(6,*) ' before receive  ',nrank,srank,rank
                tag0 = ntags*(j_proc - 1)
                tag = tag0 + 1
                icount = icount + 1

                call mpi_irecv(psip(1 + skip*(j_proc - 1)), skip                        &
                 &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)

                icount = icount + 1
                tag = tag0 + 2
                call mpi_irecv(naccm(j_proc), 1, MPI_INTEGER                         &
                     &, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                icount = icount + 1
                tag = tag0 + 3
                call mpi_irecv(psisn(j_proc), 1, MPI_DOUBLE_PRECISION               &
                     &, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                tag = tag0 + 4
                icount = icount + 1
                call mpi_irecv(indtm(1, j_proc), nel, MPI_INTEGER                     &
                     &, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                !     now receiving each big vectors independently
                tag = tag0 + 5
                icount = icount + 1
                call mpi_irecv(winv(1, j_proc), nel2wt                               &
                     &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                tag = tag0 + 6
                icount = icount + 1
                call mpi_irecv(winvfn(1, j_proc), nel2wtfn                           &
                     &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                tag = tag0 + 7
                icount = icount + 1
                call mpi_irecv(ainvup(1, j_proc), nel2up                             &
                     &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                tag = tag0 + 8
                icount = icount + 1
                call mpi_irecv(dist(1, j_proc), nelnion                              &
                     &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                tag = tag0 + 9
                icount = icount + 1
                call mpi_irecv(winvbar(1, j_proc), nel2bar                           &
                     &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                tag = tag0 + 10
                icount = icount + 1
                call mpi_irecv(winvbarfn(1, j_proc), nel2barfn                       &
                     &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                tag = tag0 + 11
                icount = icount + 1
                call mpi_irecv(jastrowall_ee(1, j_proc), dimee                       &
                     &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                tag = tag0 + 12
                icount = icount + 1
                call mpi_irecv(jastrowall_ei(1, j_proc), dimei                       &
                     &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                tag = tag0 + 13
                icount = icount + 1
                call mpi_irecv(vpotsav_ee(1, j_proc), nelsquare&
               &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)

                if (nel2wtj .ne. 0) then
                    tag = tag0 + 14
                    icount = icount + 1
                    call mpi_irecv(winvj(1, j_proc), nel2wtj                             &
                         &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                end if
                if (nel2jbar .ne. 0) then
                    tag = tag0 + 15
                    icount = icount + 1
                    call mpi_irecv(winvjbar(1, j_proc), nel2jbar                         &
                         &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                    if (iessz) then
                        tag = tag0 + 16
                        icount = icount + 1
                        call mpi_irecv(winvjbarsz(1, j_proc), nel2jbar                       &
                             &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                    end if
                end if
                if (Lbox .gt. 0) then
                    tag = tag0 + 17
                    icount = icount + 1
                    call mpi_irecv(sum_q_cos_gr(1, j_proc), n_gvec + 1                       &
                         &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                    tag = tag0 + 18
                    icount = icount + 1
                    call mpi_irecv(sum_q_sin_gr(1, j_proc), n_gvec                       &
                         &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                end if
                !      write(6,*) ' after  receive  ',nrank,srank

                if (nbra_cyrus .gt. 0) then
                    tag = tag0 + 19
                    icount = icount + 1
                    call mpi_irecv(first_cyrus(j_proc), 1, MPI_INTEGER                         &
                         &, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                    tag = tag0 + 20
                    icount = icount + 1
                    call mpi_irecv(queue_cyrus(1, 1, 1, j_proc), dim_cyrus                         &
                            &, MPI_DOUBLE_PRECISION, srank, tag, MPI_COMM_WORLD, req(icount), ierr)
                end if

            end if

#endif

        end if

    end do

#ifdef PARALLEL
    call MPI_WAITALL(icount, req, statusi, ierr)
    deallocate (req, statusi)
!   call mpi_barrier(mpi_comm_world,ierr)

    do k = 1, nrec
        j = jrec(k)
        ind = jbra(j)

        ! find the process (nrank) which j belongs to
        ! and find the process (srank) which ind belongs to
        nrank = (j - 1)/in1
        srank = (ind - 1)/in1

        j_proc = j - nrank*in1
        ind_proc = ind - srank*in1

        ! belonging to the same process

        !      write(6,*) ' before receive  ',nrank,srank,rank
        shift = 1
        call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, diag(j_proc), 1)
        shift = shift + 1
        call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, diagfn(j_proc), 1)
        shift = shift + 1
!          call dcopy(1,psip(shift+skip*(j_proc-1)),1,enert(j),1)
        call dcopy(ipc, psip(shift + skip*(j_proc - 1)), 1, enert(1, j_proc), 1)
        shift = shift + ipc
        call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, wsto(j_proc), 1)
        shift = shift + 1
        call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, psiln(j_proc), 1)
        shift = shift + 1
        call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, psidetln(j_proc), 1)
        shift = shift + 1
        call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, vpot(j_proc), 1)
        shift = shift + 1
        !     call dcopy(1,psip(shift+skip*(j_proc-1)),1,vpotreg(j_proc),1)
        !     shift=shift+1
        call dcopy(nel2, psip(shift + skip*(j_proc - 1)), 1, vpotreg(1, j_proc), 1)
        shift = shift + nel2
        call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, enertrue(j_proc), 1)
        shift = shift + 1
        call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, diffuse(j_proc), 1)
        shift = shift + 1
        call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, vcut(j_proc), 1)
        shift = shift + 1
        call dcopy(3, psip(shift + skip*(j_proc - 1)), 1, diffkin(1, j_proc), 1)
        shift = shift + 3
        call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, gradtot(j_proc), 1)
        shift = shift + 1
        call dcopy(1, psip(shift + skip*(j_proc - 1)), 1, gradtotbar(j_proc), 1)
        shift = shift + 1
        call dcopy(Ltab, psip(shift + skip*(j_proc - 1)), 1, tabpip(1, j_proc), 1)
        shift = shift + Ltab
        call dcopy(Ltabb, psip(shift + skip*(j_proc - 1)), 1, tmu(1, j_proc), 1)
        shift = shift + Ltabb
        call dcopy(Lz, psip(shift + skip*(j_proc - 1)), 1, table(1, j_proc), 1)
        shift = shift + Lz
        call dcopy(Lzr, psip(shift + skip*(j_proc - 1)), 1, tabler(1, j_proc), 1)
        shift = shift + Lzr
        call dcopy(nel2upt, psip(shift + skip*(j_proc - 1))                    &
             &, 1, winvup(1, j_proc), 1)
        shift = shift + nel2upt
        call dcopy(nel2dot, psip(shift + skip*(j_proc - 1))                    &
             &, 1, winvdo(1, j_proc), 1)
        shift = shift + nel2dot
        call dcopy(nel3, psip(shift + (j_proc - 1)*skip), 1, gradpsi(1, j_proc), 1)
        shift = shift + nel3
        call dcopy(nel3, psip(shift + (j_proc - 1)*skip)                       &
             &, 1, gradpsibar(1, j_proc), 1)
        shift = shift + nel3
        call dcopy(nelkel, psip(shift + (j_proc - 1)*skip), 1, kel(1, j_proc), 1)
        shift = shift + nelkel

        if (yesivic) then
            call dcopy(3*indt*nel, psip(shift + skip*(j_proc - 1))                 &
                 &, 1, ivic(1, 1, 1, j_proc), 1)
            shift = shift + 3*indt*nel
        end if

        if (ncore .gt. 0) then
            call dcopy(nel, psip(shift + skip*(j_proc - 1))                        &
                 &, 1, pseudolocal(1, j_proc), 1)
            shift = shift + nel
            call dcopy(nel9, psip(shift + (j_proc - 1)*skip), 1, angle(1, j_proc), 1)
            shift = shift + nel9
        end if
        call dcopy(nel9, psip(shift + (j_proc - 1)*skip), 1, angle(1, j_proc), 1)
        shift = shift + nel9

        call dcopy(np, psip(shift + (j_proc - 1)*skip), 1, econf(j_proc), in1)

        !      write(6,*) ' before receive  ',nrank,srank

    end do
    deallocate (jrec)
#endif

    return
end subroutine reshuffhub_noblock
