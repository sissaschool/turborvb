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

! sobroutine to automatically change tpar (the learning rate) during the

subroutine adjust_tpar(i_main, nweight, energy, error_energy, tpar, ngentry, itestr4)

    use allio, only: energy_list, error_energy_list, tpar_buffer_filled, &
                     tpar_increased, inc_tpar_frequency, stop_increasing_tpar, &
                     tpar_unstble_stop, counter_unstable_tpar, &
                     counter_unstable_energy, min_running_ave, min_running_std, &
                     min_running_ave_energy, min_running_std_energy, tpar_max, counter_unstable_err, divide_tpar, &
                     multiply_tpar, len_shorter_buffer, cut_sigma, tpar_buffer_len, times_tpar_decreased, &
                     use_stable_tpar, tpar_stable_list, len_tpar_stable_list, n_sigmas_tpar

    implicit none
    integer nweight, i_main, ngentry, itestr4, iter_step, start_index, stop_tpar
    real*8 energy, error_energy, tpar, running_ave, running_std, running_ave_energy, running_std_energy, ener_try, ener_sav

    integer index_ene, count_list, index_stable_list, half_stable_tpars

    half_stable_tpars = int(len_tpar_stable_list/2)

    if (use_stable_tpar) then
        !check if we have filled the list of "stable" values of tpar
        do index_stable_list = 1, size(tpar_stable_list)
            if (tpar_stable_list(index_stable_list) .eq. -1) exit
        end do

        write (*, *) tpar_stable_list

        ! if we have all the stable values of tpar we just use the mean of the last ten values
        if (index_stable_list .gt. size(tpar_stable_list)) then
            write (*, *) 'adjust_tpar: warning, using average of stable value for tpar', int(i_main/nweight)
            tpar = sum(tpar_stable_list(size(tpar_stable_list) - (half_stable_tpars - 1):size(tpar_stable_list)))/half_stable_tpars
            return
        end if
    end if

    count_list = 0
    do index_ene = 1, size(energy_list)
        if (tpar_buffer_filled(index_ene)) count_list = count_list + 1
    end do

    do index_ene = 1, size(energy_list)
        if (.not. tpar_buffer_filled(index_ene)) then
            energy_list(index_ene) = energy
            error_energy_list(index_ene) = error_energy
            tpar_buffer_filled(index_ene) = .true.
            exit
        end if
    end do

    ! If return_tpar just stores them, this subroutine skips the following part.
    if (count_list .lt. tpar_buffer_len - 1) then
        !    The  buffer  is filled even when count_list=size(energy_list)-1
        !    because the  previous loop will certainly  add one element in the list
        return
    end if

    running_ave_energy = sum(energy_list(1:count_list))/count_list
    running_std_energy = std(energy_list, running_ave_energy, count_list)

    if ((running_std_energy/error_energy)**2 .ge. 3) then !.and. count_list .ge. len_shorter_buffer) then
        write (*, *) 'adjust_tpar: warning, using only last elements of the buffer', int(i_main/nweight)
        running_ave_energy = sum(energy_list(tpar_buffer_len - len_shorter_buffer + 1:tpar_buffer_len))/len_shorter_buffer
        running_std_energy = std(energy_list(tpar_buffer_len - len_shorter_buffer + 1), running_ave_energy, len_shorter_buffer)
        running_ave = sum(error_energy_list(tpar_buffer_len - len_shorter_buffer + 1:tpar_buffer_len))/len_shorter_buffer
        running_std = std(error_energy_list(tpar_buffer_len - len_shorter_buffer + 1), running_ave, len_shorter_buffer)
    else
        running_ave = sum(error_energy_list)/size(error_energy_list)
        running_std = std(error_energy_list, running_ave, count_list)

    end if

    if (running_ave .lt. min_running_ave) then
        min_running_ave = running_ave
        min_running_std = running_std
    end if

    if (running_ave_energy .lt. min_running_ave_energy) then
        min_running_ave_energy = running_ave_energy
        min_running_std_energy = running_std_energy
    end if

    write (*, *) "adjust_tpar: minimum running average and std var/ene", &
        min_running_ave, min_running_std, min_running_ave_energy, min_running_std_energy

    iter_step = mod(i_main/nweight, inc_tpar_frequency) + 1

    !  When the counter is just filled  begin  to increase tpar  too
    if (count_list .eq. size(energy_list) - 1 .or. mod(int(i_main/nweight), inc_tpar_frequency) .eq. 0) then
        stop_increasing_tpar = .false.
        times_tpar_decreased = 0
    end if

    if (energy .gt. min_running_ave_energy + 3*min_running_std_energy) then
        counter_unstable_energy = counter_unstable_energy + 1
        !       if(energy.gt.min_running_ave_energy+6*min_running_std_energy) &
        !        counter_unstable_energy=counter_unstable_energy+3
    else
        counter_unstable_energy = 0
    end if

    if (error_energy .gt. min_running_ave + n_sigmas_tpar*min_running_std) then
        counter_unstable_err = counter_unstable_err + 1
        !        if (error_energy .gt. min_running_ave+6*min_running_std) &
        !       counter_unstable_err=counter_unstable_err+3
    else
        counter_unstable_err = 0
    end if

    if (times_tpar_decreased .eq. 0) then
        stop_tpar = tpar_unstble_stop
    else
        stop_tpar = tpar_unstble_stop*2
    end if

    if (counter_unstable_energy .ge. stop_tpar .or. counter_unstable_err .ge. stop_tpar) then

        if (tpar_increased) then

            tpar = tpar/divide_tpar

            if (use_stable_tpar) then
                ! we fill the list of stable values of tpar until it is filled
                if (tpar_stable_list(index_stable_list) .eq. -1) then
                    tpar_stable_list(index_stable_list) = tpar
                end if

            end if

            if (counter_unstable_energy .ge. stop_tpar .and. counter_unstable_err .lt. stop_tpar) then
                write (*, *) " adjust_tpar: energy too big, tpar decreased", tpar
            elseif (counter_unstable_err .ge. stop_tpar .and. counter_unstable_energy .lt. stop_tpar) then
                write (*, *) " adjust_tpar: variance too big, tpar decreased", tpar
            else
                write (*, *) " adjust_tpar: variance & energy too big, tpar decreased", tpar
            end if

            counter_unstable_energy = 0
            counter_unstable_err = 0
            stop_increasing_tpar = .true.
            times_tpar_decreased = times_tpar_decreased + 1
            if (times_tpar_decreased .gt. 1) then
                write (*, *) 'adjust_tpar: warning, reinitializing references at iteration ', int(i_main/nweight)
                tpar_increased = .false.
                min_running_ave_energy = (running_ave_energy + min_running_ave_energy)/2
                min_running_std_energy = running_std_energy
                min_running_ave = (running_ave + min_running_ave)/2
                min_running_std = running_std
            end if
        end if

    elseif (.not. stop_increasing_tpar) then
        times_tpar_decreased = 0
        tpar = tpar*multiply_tpar
        if (itestr4 .eq. -4) tpar = min(tpar, tpar_max)
        stop_increasing_tpar = .false.
        tpar_increased = .true.
    end if

    if (count_list .eq. size(energy_list)) then
        do index_ene = 1, size(energy_list) - 1
            energy_list(index_ene) = energy_list(index_ene + 1)
            error_energy_list(index_ene) = error_energy_list(index_ene + 1)
        end do
        energy_list(index_ene) = energy
        error_energy_list(index_ene) = error_energy
    end if
contains

    function std(vector, mean, dim)
        integer dim
        real*8 vector(dim)
        real*8 mean, sum, std
        integer i

        sum = 0
        do i = 1, dim
            sum = sum + (vector(i) - mean)**2
        end do
        std = sqrt(sum/(dim - 1))
    end function std

end subroutine adjust_tpar
