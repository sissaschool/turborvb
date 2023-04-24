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

subroutine write_fort11_begin
    use allio
    implicit none
    rewind (11)
    write (6, *) ' Length record unit 12  =', nmatb
    nmatb = nnozero_c
    write (11) ngg, nbra, npow, etryr, nrest, nw, tbra, np, npbra, nmatb      &
            &, nnozeroj, yesmin, membig&
            &, iread, nmp, nfat, np3, npsov, wcort, wcorwt, nweight, ngn, iend          &
            &, indcor, parbest, iflagerr, enermin, varmin                           &
            &, iesconv, fmax, jmax, parcut, inext, pippoc, pippo, nproc, nbin, stepcg     &
            &, ncg, tmes, tstep, iskipdyn, epscuttype, change_parr, change_tpar, use_stable_tpar, nbra_cyrus
    write (11) rmax, rmaxj
    if (change_parr) write (11) parr, iesdelay, delay_changeparr
    if (change_tpar) write (11) tpar, min_running_ave, min_running_ave_energy, min_running_std, &
        min_running_std_energy, tpar_buffer_filled, error_energy_list, energy_list, &
        stop_increasing_tpar, tpar_stable_list

end subroutine write_fort11_begin

subroutine write_fort11_end
    use allio
    implicit none
    integer indtel, l, i, j, k
    indtel = nel*(indt + 1)
    if (itestr .eq. -5) then
        if (contraction .gt. 0) then
            write (11) (wconfn(j), j=1, nw), (wcorw(j), j=1, nfat)&
                    &, (dek(k), k=1, ieskin)&
                    &, (dekg(k), k=1, iesking)             &
                    &, ((velion(j, i), j=1, 3), i=1, ieskindim)                              &
                    &, ((reduce(i, j), j=1, np), i=1, ncgdim)                                &
                    &, (allowcontr(i, 1), allowcontr(i, 2), i=1, nelorb_c)
        else
            write (11) (wconfn(j), j=1, nw), (wcorw(j), j=1, nfat)&
                    &, (dek(k), k=1, ieskin)             &
                    &, (dekg(k), k=1, iesking)                                            &
                    &, ((velion(j, i), j=1, 3), i=1, ieskindim)                              &
                    &, ((reduce(i, j), j=1, np), i=1, ncgdim)
        end if
    else
        write (11) (wconfn(j), j=1, nw), (wcorw(j), j=1, nfat)&
                &, (dek(k), k=1, ieskin)             &
                &, (dekg(k), k=1, iesking)
    end if
    if (allocated(detmat_proj) .and. .not. kaverage) write (11) detmat_proj, projmat_c
    if (allocated(cov_old)) write (11) cov_old
end subroutine write_fort11_end

subroutine read_fort11_begin
    use allio
    implicit none
    real*8 tstepr, parr_read, tpar_read, min_running_ave_read, &
        min_running_ave_energy_read, min_running_std_read, min_running_std_energy_read
    logical, dimension(tpar_buffer_len) :: tpar_buffer_filled_read
    real(8), dimension(tpar_buffer_len) :: energy_list_read, error_energy_list_read
    real(8), dimension(size(tpar_stable_list)) :: tpar_stable_list_read
    integer iesconv_read, iesdelay_read, delay_changeparr_read
    logical change_parr_old, change_tpar_old, stop_increasing_tpar_read, use_stable_tpar_old
    rewind (11, err=100)
    read (11, err=100, end=100) ngenc, nbrar, npow, etryr, nrest, nwr, tbrar, npr&
            &, npbrar, nmatbr, nnozerojr, yesminr, membigr&
            &, ireadr, nmpr, nfatr, np3r, npsovr&
            &, wcort, wcorwt, nweightr, ngs, iendr, indcor, parbest                   &
            &, iflagerr, enermin, varmin, iesconv_read, fmax, jmax                        &
            &, parcutr, inext, pippoc, pippo, nprocr, nbinread, stepcg, ncgread         &
            &, tmes, tstepr, iskipdynr, epscuttyper, change_parr_old, change_tpar_old, use_stable_tpar_old, nbra_cyrus_read
    if (scalermax) then
        read (11, err=100, end=100) rmax, rmaxj
    else
        read (11, err=100, end=100)
    end if

    if (change_parr_old) then
        read (11) parr_read, iesdelay_read, delay_changeparr_read
        if (iopt .eq. 0 .and. change_parr) then
            parr = parr_read
            iesdelay = iesdelay_read
            delay_changeparr = delay_changeparr_read
        end if
    end if
    if (change_tpar_old) then
        if ((iopt .eq. 2 .or. tpar .gt. 0) .and. change_tpar) then
            read (11) tpar_read, min_running_ave_read, min_running_ave_energy_read, min_running_std_read, &
                min_running_std_energy_read, tpar_buffer_filled_read, error_energy_list_read, &
                energy_list_read, stop_increasing_tpar_read, tpar_stable_list_read
        else
            read (11)
        end if
    end if
    if (iopt .eq. 0 .and. change_tpar) then
        if (tpar .gt. 0 .and. change_tpar_old) then
            use_stable_tpar = use_stable_tpar_old
            stop_increasing_tpar = stop_increasing_tpar_read
            tpar = tpar_read
            tpar_buffer_filled = tpar_buffer_filled_read
            min_running_ave = min_running_ave_read
            min_running_ave_energy = min_running_ave_energy_read
            min_running_std = min_running_std_read
            min_running_std_energy = min_running_std_energy_read
        else
            write (6, *) ' Warning using tpar read in input |tpar| and restarting buffer '
            stop_increasing_tpar = .false.
            tpar = abs(tpar)
            tpar_buffer_filled(:) = .false.
            min_running_ave = 1.d20
            min_running_std = 0.0d0
            min_running_ave_energy = 1d20
            min_running_std_energy = 0.0d0
            error_energy_list_read = 0.d0
            energy_list_read = 0.d0
        end if
        if (change_tpar_old) then
            error_energy_list = error_energy_list_read
            energy_list = energy_list_read
            tpar_stable_list = tpar_stable_list_read
        else
            error_energy_list = 0.d0
            energy_list = 0.d0
        end if
    end if

    if (iopt .eq. 0) iesconv = iesconv_read
    if (change_tstep) tstep = tstepr
    return
100 write (6, *) "From allio ERROR fort.11 begin"
    iflagerr = 1
end subroutine read_fort11_begin

subroutine read_fort11_end
    use allio
    implicit none
    integer nelnwr, i, j, k, l, indtel

    if (nwr .ne. nw .or. nprocr .ne. nproc) then
        if (nwr .ne. nw) write (6, *) ' ERROR in read_fort11: &
                & you cannot continue with different # walkers', nw, nwr
        if (nprocr .ne. nproc) write (6, *) ' ERROR in read_fort11: &
                & you cannot continue with different # processors ', nproc, nprocr
        iflagerr = 1
    end if
    nelnwr = nel*nwr*(indt + 1)
    indtel = nel*(indt + 1)
    if (iopt .eq. 2) then
        read (11, err=100, end=100) (wconfn(j), j=1, nwr)&
                &, (wcorw(j), j=1, nfatr)
    else
        if (itestr .eq. -5) then
            if (nwr .ne. nw) then
                write (6, *) ' ERROR in read_fort11: &
                        & you cannot continue with different # walkers', nw, nwr
                iflagerr = 1
            end if
            reduce = 0.d0
            if (contraction .gt. 0) then
                if (npr .eq. np) then
                    read (11, err=100, end=100) (wconfn(j), j=1, nwr)&
                            &, (wcorw(j), j=1, nfatr)&
                            &, (dek(k), k=1, ieskin)             &
                            &, (dekg(k), k=1, iesking)                                            &
                            &, ((velion(j, i), j=1, 3), i=1, ieskindim)                              &
                            &, ((reduce(mod(i - 1, ncgdim) + 1, j), j=1, np), i=1, ncgread)               &
                            &, (allowcontr(i, 1), allowcontr(i, 2), i=1, nelorb_c)
                else
                    read (11, err=100, end=100) (wconfn(j), j=1, nwr)&
                            &, (wcorw(j), j=1, nfatr)&
                            &, (dek(k), k=1, ieskin)             &
                            &, (dekg(k), k=1, iesking)                                            &
                            &, ((velion(j, i), j=1, 3), i=1, ieskindim)
                    allowcontr = .true.
                    write (6, *) ' Warning np read different from np found !!! ', np, npr
                    if (ncgdim .gt. 1) then
                        write (6, *) ' Warning previous CG steps forgotten in read '
                        stepcg = 0
                    end if
                    reduce = 0.d0
                end if

            else
                if (npr .eq. np) then
                    read (11, err=100, end=100) (wconfn(j), j=1, nwr)&
                            &, (wcorw(j), j=1, nfatr)&
                            &, (dek(k), k=1, ieskin)             &
                            &, (dekg(k), k=1, iesking)                                            &
                            &, ((velion(j, i), j=1, 3), i=1, ieskindim)                              &
                            &, ((reduce(mod(i - 1, ncgdim) + 1, j), j=1, np), i=1, ncgread)
                else
                    read (11, err=100, end=100) (wconfn(j), j=1, nwr)&
                            &, (wcorw(j), j=1, nfatr)&
                            &, (dek(k), k=1, ieskin)             &
                            &, (dekg(k), k=1, iesking)                                            &
                            &, ((velion(j, i), j=1, 3), i=1, ieskindim)
                    write (6, *) ' Warning np read different from np found !!! ', np, npr
                    if (ncgdim .gt. 1) then
                        write (6, *) ' Warning previous CG steps forgotten in read '
                        stepcg = 0
                    end if
                    reduce = 0.d0
                end if
            end if
        else
            read (11, err=100, end=100) (wconfn(j), j=1, nwr)&
                    &, (wcorw(j), j=1, nfatr)&
                    &, (dek(k), k=1, ieskin)             &
                    &, (dekg(k), k=1, iesking)
        end if
    end if
    if (allocated(detmat_proj) .and. .not. kaverage) then
        if (iopt .ne. 3) then
            read (11, err=100, end=100) detmat_proj, projmat_c
        else
            read (11, err=100, end=100)
        end if
    end if
    if (allocated(cov_old)) read (11, err=100, end=100) cov_old
    write (6, *) ' File fort.11 is read correctly ...'
    return

100 write (6, *) " From fort11_io.f90, Error fort.11 end"
    iflagerr = 1
end subroutine read_fort11_end
