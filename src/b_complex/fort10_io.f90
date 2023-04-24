!TL off
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

subroutine read_fort10(ufort10)
    use Cell, only : cellscale, map, metric, CartesianToCrystal,case_map,chosen_map
    use allio
    use convertmod
    implicit none
    integer, intent(in) :: ufort10
    integer ioptorbcontr, i1, nel2, ilaenv, iflagpip, indshell&
            &, indorb, indpar, i, ii, j, jj, k, kk, ind, indref, indrefg, nshelljmax, &
            indi, num200, id, numvjpar, numpaired, nnozero_all,dim_jasmat
    character(2) chars
    character(3) checkpbc
    character(5) checkpbc_c
    real*8 dnrm2, rc(3), r0, sumdet, costr3
    real*8, dimension(:), allocatable :: atom_s
    integer, dimension(:), allocatable :: iatom_s
    logical notequal, done,donen, molyesp
    integer, external :: check_multioptorb, iesdr1iesd, iesdr2iesd, pointvjf, num_vjpar
    real*8, external :: jastrow_ei
    real*8 rion_ref0(3), volume_unit, alphap, betap, gammap
    real(8), external :: cond_find
    integer iesd_onebody, iesd_twobody,natoms_hyb,natoms_empty
    real*8 timeg
    real*8, external:: cclock
    integer, dimension(:), allocatable:: nozeron
    logical, external:: slaterorb
#ifdef PARALLEL
  include "mpif.h"
#endif

    ! $$$$$$$$$$$$$$$$ READING fort.10 $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    if(rank.eq.0) then
        write(6, *)
        write(6, *) ' START reading the wave function fort.10 '
        write(6, *)
    endif
    case_map = 0  ! default value of case_map
    chosen_map=.false. ! No abs(iesdrr)>=100 is assumed
    if(rankrep.eq.0) then
        ! read input data
        rewind(ufort10)
        read(ufort10, *, err = 100, end = 100) chars, checkpbc_c
        yes_tilted = .false.
        yes_complex = .false. ! initialize yes_complex to .false. in ANY case
        yes_crystal = .false. ! flag for Crystal periodic basis set
        if(checkpbc_c.eq.'PBC_T') then
            yes_tilted = .true.
            iespbc = .true.
            yes_crystal = .true.
            gamma_point = .true.
            if(rank.eq.0) write(6, *) ' Warning: using Crystal periodic basis set &
                    &definition.  Complex wave function with phase is possible now! '
        elseif(checkpbc_c.eq.'PBC_C') then
            iespbc = .true.
            yes_crystal = .true.
            gamma_point = .true.
            if(rank.eq.0) write(6, *) ' Warning: using Crystal periodic basis set &
                    &definition.  Complex wave function with phase is possible now! '
        else
            rewind(ufort10)
            read(ufort10, *, err = 100, end = 100) chars, checkpbc
            if(checkpbc.eq.'PBC') then
                iespbc = .true.
                ! defining gamma_point in any case to avoid illegal intr. in mpi_bcast.
                gamma_point = .true.
                if(rank.eq.0) write(6, *) ' Warning we are assuming PBC '
            else
                iespbc = .false.
                gamma_point = .true.
            endif
        endif
        rewind(ufort10)
    endif
#ifdef PARALLEL
  call mpi_bcast(gamma_point,1,MPI_LOGICAL,0,commrep_mpi,ierr)
  call mpi_bcast(iespbc,1,MPI_LOGICAL,0,commrep_mpi,ierr)
  call mpi_bcast(yes_crystal,1,MPI_LOGICAL,0,commrep_mpi,ierr)
  call mpi_bcast(yes_complex,1,MPI_LOGICAL,0,commrep_mpi,ierr)
  call mpi_bcast(yes_tilted,1,MPI_LOGICAL,0,commrep_mpi,ierr)
#ifdef UNREL
  !   For unreliable  networks.
  call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif

#endif
    yes_hermite = .false.

    if(iespbc) then
        if(rankrep.eq.0) then
            if(rank.eq.0) write(*, *) ' Reading celldm . . . . '
            read(ufort10, *, err = 102, end = 102)
            if(.not.yes_crystal) then
                read(ufort10, *, err = 102, end = 102) rs, celldm(2:3), phase(:)
            else
                if(yes_tilted) then
                    read(ufort10, *, err = 102, end = 102) (s2r(:, k), k = 1, 3), phase(:), phase_down(:)
                    case_map = 8
                    if(rank.eq.0) write(6, *) ' Warning new mapping used case_map=8 !!! '
                else
                    ! complex w.f. : read also a phase for down-spin electrons!
                    read(ufort10, *, err = 102, end = 102) rs, celldm(2:3), phase(:), phase_down(:)
                endif
            endif
            if(sum(abs(phase(:))).eq.0.and..not.yes_crystal) then
                gamma_point = .true.
            else
                gamma_point = .false.
            endif
        endif

#ifdef PARALLEL
   call mpi_bcast(case_map,1,MPI_INTEGER,0,commrep_mpi,ierr)
     if(manyfort10.and.(.not.molyes.or.yesdft)) then
       id=rankcolrep+1
       phase(:)=xkp(:,id)
       phase_down(:)=xkp_down(:,id)
     else
    call bcast_real(phase,3,0,commrep_mpi)
    call bcast_real(phase_down,3,0,commrep_mpi)
!   call mpi_bcast(phase,3,MPI_DOUBLE_PRECISION,0,commrep_mpi,ierr)
!   call mpi_bcast(phase_down,3,MPI_DOUBLE_PRECISION,0,commrep_mpi,ierr)

     if(nbead.gt.1.and.yesdft) then
!   xkp(:,id)=phase(:)
!   xkp_down(:,id)=phase_down(:)
!    Communication phases read to xkp and xkp_down


     call mpi_allgather(phase,3,MPI_DOUBLE_PRECISION,xkp&
           &,3,MPI_DOUBLE_PRECISION,commcolrep_mpi,ierr)
     call mpi_allgather(phase_down,3,MPI_DOUBLE_PRECISION,xkp_down&
           &,3,MPI_DOUBLE_PRECISION,commcolrep_mpi,ierr)
     if(rank.eq.0) then
        ! write file "kp_weights.dat" needed to perform averages
        open(unit=37,file='kp_info.dat',form='formatted',status='unknown',position='rewind')
        write(6,*) ' Writing k-points information on file '
        write(37,'(I6)') nk
        write(37,*) '# up spin electrons '
        do i = 1,nk
           write(37,301) i,xkp(1,i),xkp(2,i),xkp(3,i),wkp(i)
        enddo
        write(37,*) '# down spin electrons '
        do i = 1,nk
           write(37,301) i,xkp_down(1,i),xkp_down(2,i),xkp_down(3,i),wkp_down(i)
        enddo
        close(37)
301     format(I6,3F13.8,1F19.14)
     endif

     endif

     endif

     if(yes_tilted) then
     call bcast_real(s2r,9,0,commrep_mpi)
     else
     call bcast_real(celldm(2),5,0,commrep_mpi)
     endif
!    call mpi_bcast(celldm(2),2,MPI_DOUBLE_PRECISION,0,commrep_mpi,ierr)
     call mpi_bcast(gamma_point,1,MPI_LOGICAL,0,commrep_mpi,ierr)
#ifdef UNREL
     !   For unreliable  networks.
     call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
#endif


        if(yes_crystal) then
            yes_hermite = .true.
        else
            yes_hermite = .false.
        endif
        opposite_phase = .true.
        same_phase = .true.

        do i = 1, 3
            if(abs(phase(i) - nint(phase(i))).ne.0.5d0.and.abs(phase(i) - nint(phase(i))).ne.0.d0) yes_complex = .true.
            ! treat gamma and other real boundaries, like (0.5,0.5,0.5), without complex w.f.
            ! unless complex algorithm is forced on input (nshell<0)
            if(yes_crystal) then
                if(abs(phase_down(i) - nint(phase_down(i))).ne.0.5d0.and.abs(phase_down(i) - nint(phase_down(i))).ne.0.d0) yes_complex = .true.
                if((abs(phase(i) - nint(phase(i))).ne.0.5.and.&
                        &(phase_down(i) - nint(phase_down(i)).ne.nint(phase(i)) - phase(i))).or.&
                        &(abs(phase(i) - nint(phase(i))).eq.0.5.and.&
                                &abs(phase_down(i) - nint(phase_down(i))).ne.0.5)) opposite_phase = .false.
                if((abs(phase(i) - nint(phase(i))).ne.0.5.and.&
                        &(phase_down(i) - nint(phase_down(i)).ne.phase(i) - nint(phase(i)))).or.&
                        &(abs(phase(i) - nint(phase(i))).eq.0.5.and.&
                                &abs(phase_down(i) - nint(phase_down(i))).ne.0.5)) same_phase = .false.
            endif
        enddo
        if(opposite_phase.and.same_phase) same_phase = .false.
    else
        rs = 0.d0
        !    define in any case at the unit vector of the mesh and the unit volume.
        at = 0.d0
        do i = 1, 3
            at(i, i) = 1.d0
        enddo
        unit_volume = 1.d0
        same_phase = .true.
        opposite_phase = .false.
    endif

    if(rankrep.eq.0) then

        if(rank.eq.0) write(6, *) ' Reading the begin . . . . '
        read(ufort10, *, err = 104, end = 104)
        read(ufort10, *, err = 104, end = 104) nelup, nel, nion
        if(nel.lt.0) then
            nel = -nel
            ipf = 2
        else
            ipf = 1
        endif
        if(nelup.le.0) then
            pfaffup = .true.
            nelup = -nelup
        else
            pfaffup = .false.
        endif

        if(nion.lt.0) then
            yesbump = .true.
            nion = -nion
        else
            yesbump = .false.
        endif

        if(rs.lt.0.d0) then
            yeslbox = .true.
            if(yes_tilted) then
                alphap = celldm(4)
                betap = celldm(5)
                gammap = celldm(6)
                volume_unit = celldm(2) * celldm(3) * dsqrt(&
                        & 1.d0 - dcos(alphap)**2.d0 - dcos(betap)**2.d0 - dcos(gammap)**2.d0 + &
                                & 2.d0 * dcos(alphap) * dcos(betap) * dcos(gammap))
                rs = (-3.d0 / 4.d0 / Pi * rs**3 * volume_unit / nel)**(1.d0 / 3.d0)
            else
                rs = (-3.d0 / 4.d0 / Pi * rs**3 * celldm(2) * celldm(3) / nel)**(1.d0 / 3.d0)
            endif
        else
            yeslbox = .false.
        endif

        read(ufort10, *, err = 106, end = 106)
        read(ufort10, *, err = 106, end = 106) nshell, nshellj
        forcecomplex = .false.
        forcesymm = .false.
        if(nshell.lt.0) then
            if(rank.eq.0.and..not.yes_complex) write(6, *) ' Warning complex wave function  !!!'
            yes_complex = .true.
            forcecomplex = .true.
            nshell = -nshell
        endif
        if(.not.yes_complex.and.yes_crystal) then
            opposite_phase = .false.
            same_phase = .true.
            do i = 1, 3
                if(abs(phase(i) - nint(phase(i))).eq.0.5.and.abs(phase_down(i) - nint(phase_down(i))).ne.0.5) same_phase = .false.
                if(abs(phase(i) - nint(phase(i))).eq.0..and.abs(phase_down(i) - nint(phase_down(i))).ne.0) same_phase = .false.
            enddo
            done = .false.
            do i = 1, 3
                !  setting to zero anyway if multiple of 2pi
                if(phase_down(i) - nint(phase_down(i)).eq.0.d0) phase_down(i) = 0.d0
                if(phase(i) - nint(phase(i)).eq.0.d0) phase(i) = 0.d0
                if(abs(phase_down(i) - nint(phase_down(i))).eq.0.5d0) then
                    phase_down(i) = -0.5d0
                    opposite_phase = .true.
                    same_phase = .false.
                    if(rank.eq.0.and..not.done) write(6, *) ' Warning set phases opposite in the real case '
                    done = .true.
                endif
                if(abs(phase(i) - nint(phase(i))).eq.0.5d0) then
                    phase(i) = 0.5d0
                    opposite_phase = .true.
                    same_phase = .false.
                    if(rank.eq.0.and..not.done) write(6, *) ' Warning set phases opposite in the real case '
                    done = .true.
                endif
            enddo
        endif

        yes_crystalj = .false.
        !comment from K.Nakano on 8th Dec. 2022.
        !The bug is trivial and it happened only with a non-orthorhombic cell without Jas.
        !We should exclude nshellj==0 case: otherwise, yes_crystalj is set .true.,
        !even though nshellj=0, which induced the double allocation in allio.f90 (kgrid). 
        !if(nshellj.lt.0.or.yes_tilted) then
        if(nshellj.lt.0.or.(yes_tilted.and.nshellj.ne.0.)) then
            if(rank.eq.0.and.iespbc) write(6, *) ' Warning Crystal basis for Jastrow  !!!'
            if(iespbc) yes_crystalj = .true.
            nshellj = abs(nshellj)
        endif

        read(ufort10, *, err = 106, end = 106)
        if(nelup.gt.nel) then
            read(ufort10, *, err = 106, end = 106) iesdrr, iesupr, npar3body, iesswr_eagp, nnozero_eagp
        else
            iesswr_eagp = 0
            nnozero_eagp = 0
            read(ufort10, *, err = 106, end = 106) iesdrr, iesupr, npar3body
        endif
        if(abs(iesdrr).ge.100) then
        chosen_map=.true.
        case_map=abs(iesdrr)/100
         if(iesdrr.gt.0) then
         iesdrr=iesdrr-100*case_map
         else
         iesdrr=iesdrr+100*case_map
         endif
        endif
        read(ufort10, *, err = 106, end = 106)
        read(ufort10, *, err = 106, end = 106) nnozero, nnozeroj
        if(nnozero.lt.0) then
            if(yes_complex) forcesymm = .true.
            nnozero = -nnozero
            if(rank.eq.0.and.yes_complex) write(6, *) ' Warning complex symmetric geminal wave function  !!!'
        endif
        read(ufort10, *, err = 106, end = 106)
        read(ufort10, *, err = 106, end = 106) iesupind, iesmind
        read(ufort10, *, err = 106, end = 106)
        read(ufort10, *, err = 106, end = 106) iesfreer, iesswr, ieskinr, ireadmin

    endif

#ifdef PARALLEL
  call mpi_bcast(nelup,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(nel,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(ipf,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(nion,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(nshell,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(nshellj,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(iesdrr,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(iesupr,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(npar3body,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(nnozero,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(nnozero_eagp,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(nnozeroj,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(iesupind,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(iesmind,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(iesfreer,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(iesswr,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(iesswr_eagp,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(ieskinr,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(ireadmin,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(case_map,1,MPI_INTEGER,0,commrep_mpi,ierr)
  call mpi_bcast(yesbump,1,MPI_LOGICAL,0,commrep_mpi,ierr)
  call mpi_bcast(yes_complex,1,MPI_LOGICAL,0,commrep_mpi,ierr)
  call mpi_bcast(yes_hermite,1,MPI_LOGICAL,0,commrep_mpi,ierr)
  call mpi_bcast(rs,1,MPI_DOUBLE_PRECISION,0,commrep_mpi,ierr)
  call mpi_bcast(yeslbox,1,MPI_LOGICAL,0,commrep_mpi,ierr)
  call mpi_bcast(forcesymm,1,MPI_LOGICAL,0,commrep_mpi,ierr)
  call mpi_bcast(pfaffup,1,MPI_LOGICAL,0,commrep_mpi,ierr)
  call mpi_bcast(yes_crystalj,1,MPI_LOGICAL,0,commrep_mpi,ierr)
  call mpi_bcast(chosen_map,1,MPI_LOGICAL,0,commrep_mpi,ierr)
#ifdef UNREL
  !   For unreliable  networks.
  call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
#endif
    if(mod(nelup, nel).ne.0) then
        numpaired = 2 * (nelup / nel)
        nelup = mod(nelup, nel)
    else
        numpaired = 2 * (nelup / nel - 1)
        nelup = nel
    endif

    if(iesdrr.eq.-12.or.iesdrr.eq.-22.or.iesdrr.eq.-26.or.iesdrr.eq.-27.or.iesdrr.eq.-30.or.iesdrr.eq.31) then
        ipj = 2
    else
        ipj = 1
    endif


    ! first define LBox for open systems
    ! if complex w.f. only option LBox.eq.3 works for open
    if(yes_complex) then
        LBox = 3.d0
        if(.not.iespbc) LBox = -3.d0
    endif

    if(yesbump.and.rank.eq.0) then
        write(6, *) ' Warning working with bumped Gaussian orbitals !!! '
    endif

    if(yes_complex) then
        ipc = 2
        if(ipf.eq.2) then
            nbdgetri = nelup_mat * ILAENV(1, 'ZSYTRF', ' ', nelup_mat, -1, -1, -1)
            !    nbdgetri=max(nbdgetri,12*nelup_mat+nelup_mat**2)
        else
            nbdgetri = nelup * ILAENV(1, 'ZGETRI', ' ', nelup, -1, -1, -1)
        endif
    else
        ipc = 1
        if(ipf.eq.2) then
            nbdgetri = nelup_mat * ILAENV(1, 'DSYTRF', ' ', nelup_mat, -1, -1, -1)
            !    nbdgetri=max(nbdgetri,12*nelup_mat+nelup_mat**2)
        else
            nbdgetri = nelup * ILAENV(1, 'DGETRI', ' ', nelup, -1, -1, -1)
        endif
    endif

    nel3 = 3 * nel

    if(iespbc) then
        if(yes_tilted) then
            givens2r = .true.
        else
            celldm(1) = 1.d0
            if(.not.yes_tilted) celldm(4:6) = 90d0 * PI / 180.d0

            alphap = celldm(4)
            betap = celldm(5)
            gammap = celldm(6)

            omega = celldm(2) * celldm(3) * dsqrt(abs(&
                    & 1.d0 - dcos(alphap)**2.d0 - dcos(betap)**2.d0 - dcos(gammap)**2.d0 + &
                            & 2.d0 * dcos(alphap) * dcos(betap) * dcos(gammap)))

            !    omega=celldm(2)*celldm(3)
            celldm(1) = (PI * nel * 4.d0 / 3.d0 / omega)**(1.d0 / 3.d0) * rs
            ! setting up direct and reciprocal lattices
            givens2r = .false.
        endif
        call InitCell(nion, nel, yes_complex)
        if(yes_tilted) then
            rs = (3.d0 * omega / (4.d0 * nel * PI))**(1.d0 / 3.d0)
        endif
        ! compute Bravais lattice symmetry operations (general routine from QE)

        call set_sym_bl(s2r)
        if(.not.manyfort10) then
            if(rank.eq.0) then
                if(.not.yes_crystal) then
                    write(*, '(a,3f12.6)') ' Phase : ', phase(:)
                else
                    write(*, '(a,3f12.6)') ' Phase up spin  : ', phase(:)
                    write(*, '(a,3f12.6)') ' Phase down spin: ', phase_down(:)
                endif
            endif
            if(sum(abs(phase(:))).eq.0.or.yes_complex) then
                gamma_point = .true.
                if(rank.eq.0) then
                    if(.not.yes_complex) write(6, *) ' Gamma Point Caculation '
                    ! complex wf can be used with any phase
                    if(yes_complex.and.sum(abs(phase(:))).ne.0.d0) write(*, *) ' Complex phase calculation '
                    if(yes_complex.and.sum(abs(phase(:))).eq.0.d0) write(*, *) ' Gamma phase calculation with complex wf'
                endif
                ! gamma_point is true not only for true gamma point, but also for complex phases (out of gamma)
            else
                gamma_point = .false.
                if(rank.eq.0) write(6, *) ' Real boundary calculation '
                ! gamma_point is false only for phases with real wave functions (real boundaries)
            endif
        else
            gamma_point = .true.
        endif

        300  format(3X, I3, X, 3F10.6)
        ! definition LBox in case of PBC
        if(yesbump) then
            LBox = 2.d0
        elseif(yes_crystal) then ! choice for non trivial phases
            LBox = 3.d0
        else
            LBox = 1.d0
        endif

    endif

    ! complex wave function is possible only with LBox.eq.3.d0
    if(abs(LBox).ne.3.d0.and.yes_complex) &
            call error('fort10_io', ' Use keyword PBC_C in the first line to use a complex wave function! ', 1, rank)

    if(ieskint.gt.ieskinr_pos + 4) &
            &  call error('fort10_io', ' Warning ieskin > ieskinr ', -1, rank)

    if(ieskint.ne.ieskinr_pos)                                 &
            &  call error('fort10_io', ' Warning ieskinr != ieskin ', -1, rank)

    ! To be sure in any event do not compute g
    cellderiv = .false.
    if(itestr.eq.-5.and.ieskint.gt.ieskinr_pos) cellderiv = .true.

    if(iespbc) then
        if(ieskint.eq.ieskinr_pos + 2) then
            if(rank.eq.0)                                              &
                    write(6, *) 'Computing Cell derivatives b,c , constant volume '
        elseif(ieskint.eq.ieskinr_pos + 3) then
            if(rank.eq.0)                                              &
                    write(6, *) 'Calculation at constant  pressure variable cell ', pressfixed
        elseif(ieskint.eq.ieskinr_pos + 1) then
            if(rank.eq.0)                                              &
                    write(6, *) 'Calculation at constant  pressure b/a c/a fixed ', pressfixed
        endif
    endif

    if(iesdrr.eq.-8) then
        iesdr = -7
        iessz = .true.
        iesgros = .true.
    elseif(iesdrr.eq.-12) then
        iesdr = -15
        iessz = .false.
        iesgros = .false.
    elseif(iesdrr.eq.-22) then
        iesdr = -20
        iessz = .false.
        iesgros = .false.
    elseif(iesdrr.eq.-27) then
        iesdr = -21
        iessz = .false.
        iesgros = .false.
    elseif(iesdrr.eq.-30) then
        iesdr = 9
        iessz = .false.
        iesgros = .false.
    elseif(iesdrr.eq.-31) then
        iesdr = 10
        iessz = .false.
        iesgros = .false.
    elseif(iesdrr.eq.-26) then
        iesdr = -7
        iessz = .false.
        iesgros = .false.
    elseif(iesdrr.eq.-10) then
        iesdr = 0
        iessz = .true.
        iesgros = .true.
    elseif(iesdrr.eq.-11) then
        iesdr = 0
        iessz = .true.
        iesgros = .false.
    elseif(iesdrr.eq.-18) then
        iesdr = -17
        iessz = .true.
        iesgros = .true.
    elseif(iesdrr.eq.-28) then
        iesdr = -20
        iessz = .true.
        iesgros = .true.
    elseif(iesdrr.eq.-9) then
        iesdr = -7
        iessz = .true.
        iesgros = .false.
    elseif(iesdrr.eq.-16) then
        iesdr = -5
        iessz = .true.
        iesgros = .false.
    elseif(iesdrr.eq.-19) then
        iesdr = -17
        iessz = .true.
        iesgros = .false.
    elseif(iesdrr.eq.-29) then
        iesdr = -20
        iessz = .true.
        iesgros = .false.
    else
        iessz = .false.
        iesgros = .true.
        iesdr = iesdrr
    endif

    !         if(iesdr.eq.-5.or.iesdr.eq.-6.or.iesdr.eq.-7.or.iesdr.eq.-17) &
    !    &     then
    !         iessingle=.true.
    !         else
    !         iessingle=.false.
    !         endif

    if(iesinv.ne.0.and..not.iessz) then
        call error('main', ' iesinv should be zero, no spin Jastrow in fort.10 ', -1, rank)
    endif

    ALLOCATE(indtm(nel, nws))
    indtm = 0
    !    IREADMIN=1, replace all the parameters (prob. averaged by readall)
    !    with the ones in the last records
    ! definition of indt (the neighbouring positions)
    !      find the maximum number of atoms below the core radius

    istart = 1
    !      if(npsa.ge.0.and.npsa.le.3) then
    !         indt=npsa*nintpseudo
    !      else

    indt = npsamax * nintpseudo
    !      endif

    if(itestr.eq.-2.or.itestr.eq.-3) indt = indt + 1
    ! add the identity
    ! in the non local DMC the term  < x | exp{-\tau V_nl} | x' > / \tau \approx 1/\tau - < x | V_nl | x'>
    ! is sampled via a heat bath algorithm
    ! the identity is splitted into single electron contributions
    if(itestr.eq.-6.or.itestr.eq.-7.or.itestr.eq.-1.or.itestr4.lt.-20) then
        if(alat2.ne.0) then
            indt = 12 + indt
            istart = 13
        else
            indt = 6 + indt
            istart = 7
        endif
    endif

    indteff = indt
    indtupt = indt
    if(itestr.eq.6) indteff = istart - 1 !  Obsolete option
    !      Use npow=1 for locality approximation.
    if(itestr.eq.-2.or.itestr.eq.-3) indtupt = indt - 1

    !#################### ALLOCATION ##########################
    ALLOCATE(naccm(nws), jbra(nwm), jbraw(nwm), tcost(max(1, indt)))
    naccm = 0
    jbra = 1
    jbraw = 1
    tcost = 0.d0
    ALLOCATE(pseudolocal(nel * nws), ivic(3, max(indt, 1), nel * nws))
    pseudolocal = 0.d0
    ivic = 0
    if(npsa.gt.0) then
        ALLOCATE(prefactor(indt - istart + 1, nel), angle(18, nel * nws))
        yesivic = .true.
    else
        ALLOCATE(prefactor(1, 1), angle(18, nel * nws))
        yesivic = .false.
    endif
    !     initialization diagonal matrix
    angle = 0.d0
    angle(1, :) = 1.d0
    angle(5, :) = 1.d0
    angle(9, :) = 1.d0
    angle(10, :) = 1.d0
    angle(14, :) = 1.d0
    angle(18, :) = 1.d0
    !         definition ivic
    versoralat = 0.d0
    if(iesrandoma) then
        yesivic = .true.
        alat = abs(alat)
    endif

    do i = 1, 3
        do kk = 1, max(indt, 1)
            do jj = 1, nel * nws
                ivic(i, kk, jj) = 0.d0
            enddo
        enddo
    enddo
    if((itestr.eq.-6.or.itestr.eq.-7.or.itestr.eq.-1           &
            &.or.itestr4.lt.-20).and.istart.ge.7) then
        do jj = 1, nel * nws
            do i = 1, 3
                ivic(i, i, jj) = alat
            enddo
            do i = 4, 6
                ivic(i - 3, i, jj) = -alat
            enddo
            if(istart.gt.7) then
                do i = 7, 9
                    ivic(i - 6, i, jj) = alat2 * alat
                enddo
                do i = 10, 12
                    ivic(i - 9, i, jj) = -alat2 * alat
                enddo
            endif
        enddo
        do i = 1, 3
            versoralat(i, i) = alat
        enddo
        do i = 4, 6
            versoralat(i - 3, i) = -alat
        enddo
        do i = 7, 9
            versoralat(i - 6, i) = alat2 * alat
        enddo
        do i = 10, 12
            versoralat(i - 9, i) = -alat2 * alat
        enddo
    endif


    ! definition of costa
    if(alat.ne.0) costa = plat(2) / alat**2
    nnozero_c = nnozero
    nnozeroj_c = nnozeroj

    iesup_c = iesupr
    iesupr_c = iesup_c
    npar3body_c = npar3body
    npar3bodyr_c = npar3body_c

    nshellr = nshell
    nshelljr = nshellj
    nshell_c = nshell
    nshellj_c = nshellj
    ! det <--> jastrow
    ! nnozero <--> nnozeroj (root)
    ! nnozero_c <--> nnozeroj_c (contracted)
    ! iesupr <--> npar3body (root)
    ! iesup_c <--> npar3body_c (contracted)
    ! iesupr_c <--> npar3bodyr_c (contracted zeta optimization)
    ! nshellr <--> nshelljr (root)
    ! nshell <--> nshellj (root zeta opt)
    ! nshell_c <--> nshellj_c (contracted)

    nelm = nel
    neldo = nel - nelup
    neldomax = max(ipc * neldo, 1)
    nelup_mat = nelup
    nel_mat = 2 * nelup
    npar_eagp = 0
    ndiff = nelup - neldo
    if(ipf.eq.2) then
        if(mod(nelup + neldo, 2).ne.0) then
            ndiff = 1 + numpaired
        else
            ndiff = numpaired
        endif
        nelup_mat = nelup + neldo + ndiff
        nel_mat = nelup_mat
        npar_eagp = (ndiff * (ndiff - 1)) / 2
    endif
    if(yes_complex) then
        if(ipf.eq.2) then
            nbdgetri = nelup_mat * ILAENV(1, 'ZSYTRF', ' ', nelup_mat, -1, -1, -1)
            !     nbdgetri=max(nbdgetri,12*nelup_mat+nelup_mat**2)
        else
            nbdgetri = nelup * ILAENV(1, 'ZGETRI', ' ', nelup, -1, -1, -1)
        endif
    else
        if(ipf.eq.2) then
            nbdgetri = nelup_mat * ILAENV(1, 'DSYTRF', ' ', nelup_mat, -1, -1, -1)
            !     nbdgetri=max(nbdgetri,12*nelup_mat+nelup_mat**2)
        else
            nbdgetri = nelup * ILAENV(1, 'DGETRI', ' ', nelup, -1, -1, -1)
        endif
    endif


    ! ################## REAL(8) ARRAY ################################

    ALLOCATE(multranspip(iesup_c), iesuptransb(iesupr_c))

    if(yes_complex) then
        allocate(scale_c(2 * (nnozero_c + nnozero_eagp)), dup_c(2 * iesupr_c))
    else
        allocate(scale_c(nnozero_c + nnozero_eagp), dup_c(iesupr_c))
    endif
    ALLOCATE(vju_c(max(npar3bodyr_c, 1)), vju(max(npar3body_c, 1))&
            &, scalej_c(max(nnozeroj_c, 1)), multranspipj(max(npar3body_c, 1))&
            &, iesuptransbj(max(npar3bodyr_c, 1)))

    !     Initialize
    dup_c = 0.d0
    scale_c = 0.d0
    vju_c = 0.d0
    vju = 0.d0
    scalej_c = 0.d0
    multranspip = 0
    multranspipj = 0
    iesuptransb = 0
    iesuptransbj = 0

    ALLOCATE(scalejsz_c(max(nnozeroj_c, 1)))
    ALLOCATE(ddw(max(iesfreer, 1))&
            &, ddwsz(max(iesfreer, 1)), dek(ieskindim), dekg(ieskingdim), scalpar(npm))

    scalejsz_c = 0.d0
    ddw = 0.d0
    ddwsz = 0.d0
    dek = 0.d0
    dekg = 0.d0
    scalpar = 0.d0

    !###################### INTEGER ARRAY #############################

    ALLOCATE(nozero_c(nnozero_c + nnozero_eagp), jbradet(nnozero_c + nnozero_eagp), sjbradet(nnozero_c)&
            &, jbradetn(3 * (nnozero_c + nnozero_eagp)), mult_c(nshell_c), kion_c(nshell_c)&
            &, ioptorb_c(nshell_c), nparam_c(nshell_c))

    nozero_c = 0.d0
    jbradet = 0
    jbradetn = 0
    sjbradet = .false.
    mult_c = 0
    kion_c = 0
    ioptorb_c = 0
    nparam_c = 0

    ALLOCATE(jbraj(max(nnozeroj_c, 1))    &
            &, jbrajn(3 * max(nnozeroj_c, 1)), multj_c(max(nshellj_c, 1))&
            &, kionj_c(max(nshellj_c, 1)), ioptorbj_c(max(nshellj_c, 1))&
            &, nparamj_c(max(nshellj_c, 1)))

    jbraj = 0
    jbrajn = 0
    multj_c = 0
    kionj_c = 0
    ioptorbj_c = 0
    nparamj_c = 0

    memtot = 3.5d0 * nnozero_c

    ALLOCATE(jbrajsz(max(nnozeroj_c, 1)))
    ALLOCATE(itouch(npm))

    jbrajsz = 0
    itouch = 0

    ! ################## READ ALL NEC. for input ######################

    ALLOCATE (wconfw(ipc * nwm), econfw(nwm), psilnw(nws)   &
            &, rcarto(3), rcart(3), dx_old(3), dx_new(3)                           &
            &, kel(3, (indt + 1) * nel * nws), keln(3, 0:indt)   &
            &, zetar(nion), atom_number(nion), type_atom(nion), rion(3, nion)       &
            &, rionsav(3, nion), pointvj(2, nion))

    wconfw = 0.0
    econfw = 0.0
    rcarto = 0.d0
    rcart = 0.d0
    dx_old = 0.d0
    dx_new = 0.d0
    kel = 0.d0
    keln = 0.d0
    zetar = 0.d0
    atom_number = 0.d0
    rion = 0.d0
    rionsav = 0.d0


    !  if(cellderiv.and.epsder.lt.0) then
    !     ip4=7
    !     if(rank.eq.0) write(6,*) ' Warning cellderiv on '
    !  else
    ip4 = 4   ! the number of calculations beyond indt in upnewwf
    !  endif

    if(.not.membig) then
        indt4 = 0
        indt4j = 0
        if(itest.ne.2) indt4j = indt + 4   ! the jastrow is always necessary
    else
        !
        !      For optimization using Hessian no way we have to know winv
        !      The same for DMC.
        !
        indt4 = indt + ip4
        indt4j = indt + ip4
    endif

    if(yesdft) then
        ! optimize memory: arrays not required when doing DFT or tools
        allocate(gradpsibar(1, 1), gradpsi(1, 1), gradtot(1) &
                &, gradtotbar(1), gradpsiold(1, 1), gradpsibarold(1, 1)           &
                &, jastrowall_ee(1, 1, 0:0, 1), jasnew_ee(1)               &
                &, jastrowall_ei(1, 1, 1), jasnew_ei(1), rcne(1, 1))

        gradpsibar = 0.d0
        gradpsi = 0.d0
        gradtot = 0.d0
        gradtotbar = 0.d0
        gradpsiold = 0.d0
        gradpsibarold = 0.d0
        jastrowall_ee = 0.d0
        jasnew_ee = 0.d0
        jastrowall_ei = 0.d0
        jasnew_ei = 0.d0
        rcne = 0.d0

    else

        allocate(gradpsibar(3, nel * nws), gradpsi(3, nel * nws), gradtot(nws) &
                &, gradtotbar(nws), gradpsiold(3, nel), gradpsibarold(3, nel)           &
                &, jastrowall_ee(nel, nel, 0:indt4j, nws), jasnew_ee(nel)               &
                &, jastrowall_ei(nion, nel, nws), jasnew_ei(nion)                      &
                &, rcne(3, nel))

        gradtot = 0.d0
        gradtotbar = 0.d0
        gradpsiold = 0.d0
        gradpsibarold = 0.d0
        jastrowall_ee = 0.d0
        jasnew_ee = 0.d0
        jastrowall_ei = 0.d0
        jasnew_ei = 0.d0
        rcne = 0.d0

        do j = 1, 3
            do jj = 1, nel * nws
                gradpsi(j, jj) = 1.d0
                gradpsibar(j, jj) = 1.d0
            enddo
        enddo

    endif


    !##################### END ALLOCATION #############################

    !          Definition scalpar
    do i = 1, nmat
        scalpar(i) = -1.d0
    enddo

    do i = 1, ncore
        scalpar(icore(i)) = tcore(i)
    enddo

    do i = 1, npp
        itouch(i) = 1
    enddo


    ! read atomic numbers and nuclear positions

    if(rankrep.eq.0) then

        read(ufort10, *, err = 108, end = 108)
        do j = 1, nion
            read(ufort10, *, err = 108, end = 108) zetar(j), atom_number(j), (rion(i, j), i = 1, 3)
            if(atom_number(j).lt.0.d0) zetar(j) = 0.d0
        enddo

        allocate(atom_s(nion), iatom_s(nion))
        atom_s(1:nion) = atom_number(1:nion)
        call dsortx(atom_s, 1, nion, iatom_s)
        cost = 0.d0
        nmax_ion = 0
        do i = 1, nion
            indi = iatom_s(i)
            if(atom_number(indi).ne.cost) then
                nmax_ion = nmax_ion + 1
                type_atom(indi) = nmax_ion
                cost = atom_number(indi)
            else
                type_atom(indi) = nmax_ion
            endif
        enddo
        deallocate(atom_s, iatom_s)
        if(rank.eq.0) write(6, *) ' Number of different atomic species =', nmax_ion

    endif

#ifdef PARALLEL

  call mpi_bcast(nmax_ion,1,MPI_INTEGER,0,commrep_mpi,ierr)
  !   For unreliable  networks.
  call bcast_real(zetar,nion,0,commrep_mpi)
  call bcast_real(atom_number,nion,0,commrep_mpi)
  call bcast_integer(type_atom,nion,0,commrep_mpi)
  call bcast_real(rion,3*nion,0,commrep_mpi)

#ifdef UNREL
  call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif

#endif
    if(iespbc) &
            call purge_isymm(.true., nrot, isymm, sname, t_rev, rion, nion, atom_number, cellscale)

    if(iespbc) then

        if(iesdr.eq.-2) &
                call error('read_fort10', ' Jastrow = -2 not implemented for periodic systems!', 1, rank)

        if(rank.eq.0) write(6, *) " Periodic System "

        derEVp = 0.d0
        errEVp = 0.d0

        allocate(rmucos(3, nion * (max(indt, 1) + 1)), rmusin(3, nion * (max(indt, 1) + 1)))
      allocate(cosphase(0:max(indt, 1),nion), sinphase(3, 0:max(indt, 1),nion))
        ! initialization for gamma point
        rmucos = 0.d0
        rmusin = 0.d0
        cosphase = 1.d0
        sinphase = 0.d0
        Linv = 1.d0
        ! variables for AAD
        if(.not.yesdft) then
    allocate(cosphaseb(0:max(indt, 1),nion), sinphaseb(3, 0:max(indt, 1),nion))
            cosphaseb = 0.d0
            sinphaseb = 0.d0
        endif

        if(cellscale(3).ne.0.d0.and.cellscale(2).ne.0.d0) then
            kappa = kappar / min(cellscale(1), cellscale(2), cellscale(3))
        elseif(cellscale(2).ne.0.d0) then
            kappa = kappar / min(cellscale(1), cellscale(2))
        else
            kappa = kappar / cellscale(1)
        endif
     if(yes_tilted.and.ksq.ne.0.5d0) kappa = kappa/cond_find(metric(1,2),metric(1,3),metric(2,3))

        if(rank.eq.0) then

            write(6, *) ' Rs = ', rs
            write(6, *) ' Notice celldm(1) set to ', 1
            write(6, *) ' Unit of length a.u. '
            write(6, *) ' LMin = ', lmin
            write(6, *) ' Real Volume of the Cell = ', omega
            write(6, *) ' Rescaled tstep ', tstep
            write(6, *) ' celscale before ', cellscale(1:3)
            if(yes_tilted.and.ksq.ne.0.5d0) write(6,*) ' Warning increased kappa Ewald by = ',1.d0/cond_find(metric(1,2),metric(1,3),metric(2,3))

            if(cellscale(3).ne.0.d0.and.cellscale(2).ne.0.d0) then
                lmin = min(cellscale(1), cellscale(2), cellscale(3))
            elseif(cellscale(2).ne.0.d0) then
                lmin = min(cellscale(1), cellscale(2))
            else
                lmin = cellscale(1)
            endif

            if(tstep.gt.lmin / 2.0) call error('read_fort10', " WARNING ! tstep > LBox/2 ", -1, rank)
            if(rank.eq.0) write(6, *) ' Rescaling Distances . . . '

        endif

#ifdef PARALLEL
     call mpi_bcast(lmin,1,MPI_DOUBLE_PRECISION,0,commrep_mpi,ierr)
#endif
        allocate(derEV(nws), p_pulay(nws))

        derEV = 0.d0
        p_pulay = 0.d0

        !      *************** PERIODIC WORLD *****************
    else

        if(rank.eq.0) write(6, *) ' Open Boundary Conditions '
        kMax = 1
        kMax2 = 2
        ALLOCATE(rmucos(3, nion * (indt + 1)), rmusin(3, nion * (indt + 1)))
        rmucos = 0.d0
        rmusin = 0.d0
        !       write(6,*) ' After allocation '
        if(yesbump) then
            LBox = -2.d0
        elseif(yes_complex) then
            LBox = -3.d0
        else
            LBox = -1.d0
        endif
        kappa = -1.d0
        kSq = -1.d0
        !        Unit Hartree
        ! do ii=1,4
        !   ris(ii)=1.d0
        ! enddo
    endif
    ! Unit Hartree
    ris(1) = 1.d0
    ris(2) = 0.5d0
    ris(3) = 0.5d0
    ris(4) = ris(2)**2
    if(oldscaling) then
        ris(5) = 0.5d0
    else
        ris(5) = 1.d0 / sqrt(2.d0) ! This is the scaling of energy  1/s   E(Ryd)=s^2 E(Hartree)
    endif

    LBoxj = LBox
    if(LBox.eq.3.and..not.yes_crystalj) LBoxj = 1.d0   ! Jastrow with standard PBC

    ! rescaling unit energy
    tbra = tbra * ris(2)
    etry = etry / ris(2)
    shift = shift / ris(2)
    temp = temp / ris(2)
    pressfixed = pressfixed / ris(2)
    if(ldynsecond) then
        friction = friction / ris(5)
        if(.not.oldscaling) delta0 = delta0 / ris(5)
    else
        friction = friction / ris(2)
        if(.not.oldscaling) delta0 = delta0 / ris(2)
    endif

    ! two body jastrow coefficients (if iesdr=-5,-6,-7)
    allocate(costz(nion), costz3(nion), zetaq(nion))

    if(npsar.gt.0) then

        do i = 1, nion
            ind = 1
            do while(ind.le.npsar.and.kindion(ind).ne.i)
                ind = ind + 1
            enddo

            cost = 2.d0 * zetar(i)
            if(ind.le.npsar) then

                if(LBox.gt.0) then

                    if(rcutoff(ind).gt.Lmin) then
                        if(rank.eq.0) write(6, *) ' Warning core radius too large for pseudo ion: ', ind, i
                    endif
                endif

                do k = 1, nparpshell(pshell(ind), ind)
                    if(parshell(2, jpseudo(pshell(ind), ind) + k - 1).eq.1) then
                        cost = cost - 2.d0 * parshell(1, jpseudo(pshell(ind), ind) + k - 1)
                    endif
                enddo
                if(cost.eq.0.and.softcusp) then
                    costr3 = 0.d0
                    do k = 1, nparpshell(pshell(ind), ind)
                        if(parshell(2, jpseudo(pshell(ind), ind) + k - 1).eq.1) then
                            costr3 = costr3 + parshell(1, jpseudo(pshell(ind), ind) + k - 1) * &
                                    &parshell(3, jpseudo(pshell(ind), ind) + k - 1)
                        endif
                    enddo
                    costr3 = costr3 / 6.d0
                    !          write(6,*) ' r^3 singularity found =',costr3
                endif
            endif






            !              Rescaling pseudo  ADD PSEUDO PBC +4 lines
            !  In the code for simplicity of Ewald there is a global rescaling
            !   of variables  (x',y',z')= (x,y,z)/LBox
            !  The one-body potential following Pierleoni et al. PRE 68, p. 046707 (
            !  is defined as u_1(r)=  (2 z)^(3/4) u ( (2z)^{1/4} r)
            !   where u is the two-body potential satisfying the cusp condition
            !     u'(r=0)=1/2
            !  In order to satisfy the cusp elctron-ion condition it is enough that
            !     u'_1(r=0)=Z

            if(zetar(i).ne.0.d0) then
                if(cost.ne.0.d0.or..not.softcusp) then
                    costz(i) = cost**0.25d0
                    costz3(i) = cost**0.75d0
                else
                    !         soft pseudo for atom i
                    costz(i) = 1.d0
                    costz3(i) = costr3
                endif
            else
                costz(i) = 0.d0
                costz3(i) = 0.d0
            endif

        enddo
    else
        do i = 1, nion
            if(zetar(i).ne.0.d0) then
                costz(i) = (2.d0 * zetar(i))**0.25d0
                costz3(i) = (2.d0 * zetar(i)) / costz(i)
            else
                costz(i) = 0.d0
                costz3(i) = 0.d0
            endif
        enddo
    endif

    ! define zetaq (effective charge for Q caffarel)
    do i = 1, nion
        if(costz(i).ne.1.d0) then
            zetaq(i) = costz(i) * costz3(i) / 2.d0
        else
            zetaq(i) = 0.d0
        endif
    enddo


    !  definition n_body_on
    if(dnrm2(nion, costz, 1).ne.0.d0.and.iesdr.ne.0) then
        n_body_on = 1   ! one_body_on
    else
        n_body_on = 0   ! no one body
    endif
    if(n_body_on.eq.0) add_onebody2det = .false. ! we cannot add in any case
    !  read the set of constraints  for the forces
    if(rankrep.eq.0) then
        read(ufort10, *, err = 110, end = 110)
    endif

    indn = 0
    indnn = 0
    indp = 0

    iscraipsip = max(2 * nwm, 2 * (nnozero + nnozero_eagp) + 1, 2 * nnozeroj + 1, 9 * nion + 1, iesupr + 1)
    iscraipsip = max(iscraipsip, 6 * npm + 2 + 5 * npmn)

    if(symmagp.and.ipc.eq.2) then
        iscraipsip = max(iscraipsip, 6 * nnozero_c)
    else
        iscraipsip = max(iscraipsip, 2 * ipc * nnozero_c)
    endif
    iscraipsip = max(iscraipsip, 2 * nnozeroj_c + 1)

    iscraipsip = max((ipf + 1) * nelup_mat, iscraipsip)  ! for upscratch

    iscraipsip = max(iscraipsip, 18 * nion) ! for forces when used

    if(yesdft) then
        iscraipsip = max(nelup, iesupr_c + iesupind, npar3bodyr_c + iesmind, 9 * nion + 1, 2 * nnozero + 1, 2 * nnozeroj + 1)
    else
        iscraipsip = max(iscraipsip, iesupr_c + iesupind                  &
                &, npar3bodyr_c + iesmind)
    endif

    iscraipsip = max(iscraipsip, nel) ! for forces when used

    ALLOCATE(ipsip(iscraipsip))

    ipsip = 0

    if(rank.eq.0) write(6, *) ' Reading ieskin . . . . '
    !            if(rank.eq.0.and.ieskinr.eq.3*nion.and.iesking.eq.0) then
    !            write(6,*) 'Warning assuming all force components included '
    !            endif

    allocate(ion_table(ieskinr))
    indn = 0
    do i = 1, ieskinr
        ! ion index, coordinate index
        if(rankrep.eq.0) read(ufort10, *, err = 112, end = 112) ii, (ipsip(j), j = 1, 2 * abs(ii))
        !              if(ieskinr.eq.3*nion.and.iesking.eq.0) then
        !                ii=1
        !                ipsip(1)=(i-1)/3+1
        !                ipsip(2)=mod(i-1,3)+1
        !              endif
#ifdef PARALLEL
     call mpi_bcast(ii,1,MPI_INTEGER,0,commrep_mpi,ierr)
#ifdef UNREL
     !   For unreliable  networks.
     call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
     call bcast_integer(ipsip,2*abs(ii),0,commrep_mpi)
#ifdef UNREL
     !   For unreliable  networks.
     call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
#endif
     ! indp is the number of irreducible moved ionic coordinates,
        ! later as movedion
        ion_table(i)%mult = ii
        indp = abs(ii)
        allocate(ion_table(i)%ion(indp))
        allocate(ion_table(i)%comp(indp))
        if(ii.gt.0) indn = indn + indp
        do j = 1, indp
            ! indn is the number of moved ionic coordinates
            ion_table(i)%comp(j) = abs(ipsip(2 * j))
            ion_table(i)%ion(j) = ipsip(2 * j - 1)
            if(ipsip(2 * j).lt.0) ion_table(i)%ion(j) = -abs(ipsip(2 * j - 1))
            !                  indn=indn+1
            !                  ! nozeroforces a mapping table indn --> 1~3*nion
            !                  nozeroforces(indn)=abs(ipsip(2*j))+3*(abs(ipsip(2*j-1))-1)
            !                  if(ipsip(2*j-1).gt.0) then
            !                    !jbraforces is a mapping table indn --> indp with \pm sign
            !                    jbraforces(indn)=indp
            !                  elseif(ipsip(2*j-1).lt.0) then
            !                    jbraforces(indn)=-indp
            !                  endif
        enddo
        !              else
        !                indnn=indnn+1
        !                jbraforcesn(indnn)=ii
        !                do j=1,-2*ii
        !                  indnn=indnn+1
        !                  jbraforcesn(indnn)=ipsip(j)
        !                enddo
        !              endif
    enddo

    if(rank.eq.0) write(6, *) 'Number of moved ionic coordinates', indn
    movedion = indn
    stopdyn = .false.

    ! read 2 body jastrow
    if(rankrep.eq.0) then
        read(ufort10, *, err = 114, end = 114)
    endif

    if(iesdr.ne.0) then
        if(rankrep.eq.0) then
            if(rank.eq.0) write(6, *) ' Reading 2-body jastrow . . . . '
            allocate(psip(10000))
            read(ufort10, *, err = 118, end = 118) niesd, (psip(i), i = 1, abs(niesd))
            if(niesd.lt.0) then
                symmagp = .false.
                niesd = -niesd
            else
                symmagp = .true.
            endif
        endif
#ifdef PARALLEL
     call mpi_bcast(symmagp,1,MPI_LOGICAL,0,commrep_mpi,ierr)
     call mpi_bcast(niesd,1,MPI_INTEGER,0,commrep_mpi,ierr)
#endif
        if(niesd.gt.10000) call error('read_fort10', 'Too many variables in the two body Jastrow ! ', 1, rank)
        allocate(vj(niesd))
        if(rankrep.eq.0) then
            vj(1:niesd) = psip(1:niesd)
            deallocate(psip)
        endif
#ifdef PARALLEL
     call mpi_bcast(vj,niesd,MPI_DOUBLE_PRECISION,0,commrep_mpi,ierr)
#endif

    else

        if(rankrep.eq.0) then
            read(ufort10, *, err = 120, end = 120) niesd
        endif
#ifdef PARALLEL
     call mpi_bcast(niesd,1,MPI_INTEGER,0,commrep_mpi,ierr)
#endif
        if(niesd.lt.0) then
            symmagp = .false.
            niesd = 0
        else
            symmagp = .true.
        endif
        allocate(vj(1))
        vj = 0.d0
        !  deallocate(psip)
    endif
    !  setting to true before checking the symmagp case
    !   if(niesd.gt.2.and.niesd.ne.nmax_ion+1) then
    !    call error( ' read_fort10 ', ' Number of parameters 1-2body J should be equal to &
    !  & 1,2 or nmax_ion+1 ',1,rank)
    !   endif

    if(iespbc.and.yes_crystal) then
        allowed_averagek = .true.
        !  At gamma it is a useless operation
        if(sum(abs(phase(:) - nint(phase(:)))).eq.0.and.&
                &sum(abs(phase_down(:) - nint(phase_down(:)))).eq.0.and..not.manyfort10) allowed_averagek = .false.
    else
        allowed_averagek = .false.
    endif
    ! up to this line the definition of yes_hermite is the same as in makefort10
    ! just check whether the boundary conditions for opposite spins are opposite.

    if(rank.eq.0.and.allowed_averagek) &
            write(6, *) ' Warning algorithm with effective lambda for k-average '

    !   Calculation yes_hermite
    if(.not.opposite_phase) yes_hermite = .false.
    if(yes_hermite.and..not.symmagp) yes_hermite = .false.

    if(forcesymm) then
        yes_hermite = .false.
        opposite_phase = .false.
        same_phase = .true.
    endif

#ifdef PARALLEL
   if(manyfort10) then
    if(rank.eq.0) write(6,*) ' Warning forcing all to have the same &
    &  Hermite/non Hermitian case=',yes_hermite
    call mpi_bcast(yes_hermite,1,MPI_LOGICAL,0,mpi_comm_world,ierr)
    call mpi_bcast(opposite_phase,1,MPI_LOGICAL,0,mpi_comm_world,ierr)
    call mpi_bcast(same_phase,1,MPI_LOGICAL,0,mpi_comm_world,ierr)
   endif
#endif

    if(rank.eq.0) write(6, *) 'opposite phase =', opposite_phase
    if(rank.eq.0) write(6, *) 'same  phase =', same_phase
    if(rank.eq.0) write(6, *) 'yes_hermite before =', yes_hermite

    if(yes_crystal) then
        if(symmagp) then
            if(same_phase) then
                do ii = 1, 3
                    if((abs(phase(ii) - nint(phase(ii))).ne.0.5d0.and.abs(phase_down(ii) - nint(phase_down(ii))).ne.0.5d0&
                            &.and.phase(ii) - nint(phase(ii)).ne.phase_down(ii) - nint(phase_down(ii))).or.&
                            &(abs(phase(ii) - nint(phase(ii))).eq.0.5.and.abs(phase_down(ii) - nint(phase_down(ii))).ne.0.5d0).or.&
                            &(abs(phase(ii) - nint(phase(ii))).ne.0.5.and.abs(phase_down(ii) - nint(phase_down(ii))).eq.0.5d0)) then
                        call error(' read_fort10 ', 'Phases for up and down spin &
                                & electrons are not compatible with symmagp=.true.: symmetric or hermitian &
                                &  AGP!! Please check your fort.10 (e.g. use forcesymm=.true.) ', 1, rankrep)
                    endif
                enddo
            elseif(opposite_phase) then
                do ii = 1, 3
                    if((abs(phase(ii) - nint(phase(ii))).ne.0.5d0.and.abs(phase_down(ii) - nint(phase_down(ii))).ne.0.5d0&
                            &.and.phase(ii) - nint(phase(ii)).ne.nint(phase_down(ii)) - phase_down(ii)).or.&
                            &(abs(phase(ii) - nint(phase(ii))).eq.0.5.and.abs(phase_down(ii) - nint(phase_down(ii))).ne.0.5d0).or.&
                            &(abs(phase(ii) - nint(phase(ii))).ne.0.5.and.abs(phase_down(ii) - nint(phase_down(ii))).eq.0.5d0)) then
                        call error(' read_fort10 ', 'Phases for up and down spin &
                                & electrons are not compatible with symmagp=.true.: symmetric or hermitian &
                                &  AGP!! Please check your fort.10 (e.g. use forcesymm=.true.)', 1, rankrep)
                    endif
                enddo
            endif
        endif
    endif

    if(symmagp.and.ipc.eq.2.and.yes_correct) then
        allocate(dsw(4 * (iesswr + iesswr_eagp)))
    else
        allocate(dsw(ipc * (iesswr + iesswr_eagp)))
    endif
    dsw = 0.d0

    if(iesd.ne.0.and.iesd.ne.niesd) then
        if(rank.eq.0) write(6, *) ' Warning iesd changed to ', niesd
        iesd = niesd
    endif

    ! read orbitals determinantal part
    if(rankrep.eq.0) then
        read(ufort10, *, err = 122, end = 122)
    endif
    indpar = 0
    occ = 0
    contraction = 0
    nshell_max = 0
    molyes = .false.
    molyesp = .false.
    hybyes = .false.
    molecular = 0
    nelorbmax = 0
    maxshell = 1
    allocate(whereiesup(iesup_c))
    whereiesup = 0
    if(rank.eq.0) write(6, *) ' Reading det shells . . . . '
    do i = 1, nshell_c
        if(rankrep.eq.0) then
            read(ufort10, *, err = 124, end = 124) mult_c(i), nparam_c(i), ioptorb_c(i)
            if(check_multioptorb(mult_c(i), ioptorb_c(i), 'E').eq.1) then
                call error(' read_fort10 ', 'Error checking determinant orbitals multiplicity.', 1, rank)
            endif
        endif
#ifdef PARALLEL
     call mpi_bcast(mult_c(i),1,MPI_INTEGER,0,commrep_mpi,ierr)
     call mpi_bcast(nparam_c(i),1,MPI_INTEGER,0,commrep_mpi,ierr)
     call mpi_bcast(ioptorb_c(i),1,MPI_INTEGER,0,commrep_mpi,ierr)
#ifdef UNREL
     !   For unreliable  networks.
     call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
#endif

        if(mult_c(i).gt.maxshell) maxshell = mult_c(i)

        if(molyesp.and.ioptorb_c(i).lt.900000) &
                call error('read_fort10', ' You cannot use atomic orbital after molecular ! ', 1, rank)

        if(ioptorb_c(i).ge.900000) molyesp = .true.
        if(ioptorb_c(i).eq.900000) hybyes = .true.
        if(ioptorb_c(i).eq.1000000)  molecular = molecular + 1

        if(ioptorb_c(i).ne.ioptorbcontr(ioptorb_c(i), LBox, 1)) then
            contraction = contraction + 1

            if(.not.molyesp) then
                nshell_max = nshell_max + nparam_c(i) / 2
                nelorbmax = nelorbmax + mult_c(i) * nparam_c(i) / 2
            endif

        elseif(.not.molyesp) then
            nshell_max = nshell_max + 1
            nelorbmax = nelorbmax + mult_c(i)
        endif

        if(rankrep.eq.0) then
            if(yes_complex.and.nparam_c(i).gt.1) then ! complex w.f. with contracted orb
                read(ufort10, *, err = 126, end = 126) kion_c(i), (dup_c(2 * (indpar + j) - 1), j = 1, nparam_c(i) / 2), & ! reading exponents
                        (dup_c(2 * indpar + nparam_c(i) + j), j = 1, nparam_c(i))        ! reading complex coefficients
                do j = 2 * indpar + 2, 2 * indpar + nparam_c(i), 2
                    dup_c(j) = 0.d0
                enddo
            elseif(yes_complex.and.nparam_c(i).eq.1) then ! complex w.f. with uncontracted orb
                read(ufort10, *, err = 126, end = 126) kion_c(i), dup_c(2 * indpar + 1)
                dup_c(2 * indpar + 2) = 0.d0 ! exponents always real
            else ! real w.f.
                read(ufort10, *, err = 126, end = 126) kion_c(i), (dup_c(indpar + j), j = 1, nparam_c(i))
            endif
        endif

        !         write(6,*) ' Reading orbital =',i
        !         do j=1,2*nparam_c(i)
        !         write(6,*) j,dup_c(j)
        !         enddo

#ifdef PARALLEL
     call mpi_bcast(kion_c(i),1,MPI_INTEGER,0,commrep_mpi,ierr)
     if(yes_complex) then
        call mpi_bcast(dup_c(2*indpar+1),2*nparam_c(i),MPI_DOUBLE_PRECISION,0,commrep_mpi,ierr)
     else
        call mpi_bcast(dup_c(indpar+1),nparam_c(i),MPI_DOUBLE_PRECISION ,0,commrep_mpi,ierr)
     endif
#ifdef UNREL
     !   For unreliable  networks.
     call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
#endif

        do ii = 1, nparam_c(i)
            whereiesup(indpar + ii) = kion_c(i)
        enddo

        indpar = indpar + nparam_c(i)
        occ = occ + mult_c(i)
    enddo  ! end do i
#ifdef PARALLEL
        call mpi_bcast(hybyes,1,MPI_LOGICAL,0,commrep_mpi,ierr)
#endif

    !         stop

    !         calculation upper bound for maxparam
       maxparam = 1
       natoms_hyb=0
       do i=1,nion
          done=.false.
          do j=1,nshell_c
          if(kion_c(j).eq.i.and.ioptorb_c(j).eq.900000) done=.true. 
          enddo
          if(done) natoms_hyb=natoms_hyb+1
       enddo
       natoms_empty=0
       do i=1,nion
          donen=.true.
          do j=1,nshell_c
          if(kion_c(j).eq.i.and.ioptorb_c(j).lt.900000) donen=.false. 
          enddo
          if(donen) natoms_empty=natoms_empty+1
       enddo

     if(natoms_hyb.gt.0.and.natoms_hyb+natoms_empty.ne.nion) then
      call error(' read_fort10 ', 'ALL atoms should have hybrid orbitals!', 1, rank)
     endif
    do i = 1, nshell_c
        if(ioptorb_c(i).le.900000) then
            ind = 0
            do j = 1, nshell_c
                if(kion_c(j).eq.kion_c(i).and.ioptorb_c(j).le.900000) ind = ind + 1
            enddo
            if(ind.gt.maxparam) maxparam = ind
        endif
    enddo
    maxparam = maxparam + molecular

    occ_c = occ

    do i = 2, nshell_c
        if(kion_c(i).lt.kion_c(i - 1)) then
            if(rank.eq.0) then
                write(6, *) 'Warn. Ions are not in ascending order in Det !'
                write(6, *) 'Check shell det =', i
            endif
        endif
    enddo

    if(indpar.ne.iesup_c) &
            call error('read_fort10', ' Some errors in input Slater part ', 1, rank)

    if(contraction.ne.0) then

        if(rank.eq.0)  then
            write(6, *) 'USING', contraction, 'CONTRACTED DET ORBITALS'
            write(6, *) 'USING', molecular, 'MOLECULAR  DET ORBITALS'
        endif

        ! definition iesup_atom
        iesup_atom = 0
        do i = 1, nshell_c
            if(ioptorb_c(i).ne.1000000) iesup_atom = iesup_atom + nparam_c(i)
        enddo
        if(rank.eq.0) write(6, *)  ' iesup dimension  atomic basis =', iesup_atom, iesup_c

    else
        iesup_atom = iesup_c

        ALLOCATE(dupr(iesupr), duprold(iesupr))
        ALLOCATE(mult(nshell), kion(nshell), ioptorb(nshell)                &
                &, nparam(nshell))
        duprold = 0
        dupr = 0.d0
        mult = 0
        kion = 0
        ioptorb = 0
        nparam = 0

        ALLOCATE(nozero(nnozero + nnozero_eagp))
        ! storing vector
        if(yes_complex) then
            allocate(scale(2 * nnozero))
        else
            allocate(scale(nnozero))
        endif
        scale = 0.d0
        nozero = 0

        allocate(nozerodet(nnozero + nnozero_eagp))
        nozerodet = 0
        if(yes_complex) then
            do i = 1, iesup_c
                dupr(i) = dup_c(2 * i - 1)
            enddo
        else
            do i = 1, iesup_c
                dupr(i) = dup_c(i)
            enddo
        endif
        do i = 1, nshell_c
            kion(i) = kion_c(i)
            mult(i) = mult_c(i)
            ioptorb(i) = ioptorb_c(i)
            nparam(i) = nparam_c(i)
        enddo
        occtot = occ
        nshelldo = nshell
    endif

    ! read jastrow orbitals
    if(rankrep.eq.0) then
        read(ufort10, *, err = 128, end = 128)
    endif
    allocate(whereiesm(max(npar3body_c, 1)))
    whereiesm = 0
    indpar = 0
    occ = 0
    contractionj = 0
    num200 = 0
    moljyes = .false.
    nshellj_max = 0
    nelorbmaxj = 0
    molecularj = 0
    maxshellj = 1
    maxparamj = 0
    numcost = 0
    if(rank.eq.0) write(6, *) ' Reading jas shells . . . . '
    do i = 1, nshellj_c
        if(rankrep.eq.0) then
            read(ufort10, *, err = 130, end = 130) multj_c(i), nparamj_c(i), ioptorbj_c(i)
            if(check_multioptorb(multj_c(i), ioptorbj_c(i), 'E').eq.1) then
                call error(' read_fort10 ', 'Error checking Jastrow orbitals multiplicity.', 1, rank)
            endif
        endif
#ifdef PARALLEL
     call mpi_bcast(multj_c(i),1,MPI_INTEGER,0,commrep_mpi,ierr)
     call mpi_bcast(nparamj_c(i),1,MPI_INTEGER,0,commrep_mpi,ierr)
     call mpi_bcast(ioptorbj_c(i),1,MPI_INTEGER,0,commrep_mpi,ierr)
#ifdef UNREL
     !   For unreliable  networks.
     call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
#endif

        if(moljyes.and.ioptorbj_c(i).lt.900000) &
                call error('read_fort10', ' You cannot use atomic orbital after molecular ! ', 1, rank)

        if(ioptorbj_c(i).eq.200) numcost = numcost + 1
        if(maxshellj.lt.multj_c(i)) maxshellj = multj_c(i)

        if(ioptorbj_c(i).ge.900000) moljyes = .true.
        if(ioptorbj_c(i).eq.1000000)  molecularj = molecularj + 1

        if(ioptorbj_c(i).ne.ioptorbcontr(ioptorbj_c(i), LBoxj, 1)&
                &.or.(num200.gt.0.and.ioptorbj_c(i).eq.200)) then
            if(num200.eq.1.and.ioptorbj_c(i).eq.200) then
                if(rank.eq.0) write(6, *) ' Warning >1  constant orbitals--> considered contracted '
            endif
            contractionj = contractionj + 1
            if(.not.moljyes) then
                nshellj_max = nshellj_max + max(nparamj_c(i) / 2, 1)
                nelorbmaxj = nelorbmaxj + max(multj_c(i) * nparamj_c(i) / 2, 1)
            endif
        elseif(.not.moljyes) then
            nshellj_max = nshellj_max + 1
            nelorbmaxj = nelorbmaxj + multj_c(i)
        endif

        if(ioptorbj_c(i).eq.200) num200 = num200 + 1

        if(rankrep.eq.0) then
            read(ufort10, *, err = 132, end = 132) kionj_c(i), (vju_c(indpar + j), j = 1, nparamj_c(i))
        endif
#ifdef PARALLEL
     call mpi_bcast(kionj_c(i),1,MPI_INTEGER,0,commrep_mpi,ierr)
     if(nparamj_c(i).gt.0) &
          call mpi_bcast(vju_c(indpar+1),nparamj_c(i),MPI_DOUBLE_PRECISION,0,commrep_mpi,ierr)
#ifdef UNREL
     !   For unreliable  networks.
     call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
#endif
        ! ion index: kion_c (det) --> kionj_c (jastrow)
        ! orbital coeff contracted basis: dup_c (det) --> vju_c (jastrow)

        do jj = 1, nparamj_c(i)
            whereiesm(indpar + jj) = kionj_c(i)
        enddo

        indpar = indpar + nparamj_c(i)
        occ = occ + multj_c(i)
    enddo

     if(yes_sparse_choose) then
       if(int(nnozeroj,8)*max_sparse_choice.lt.int(nelorbmaxj*ipj,8)**2) then
       yes_sparse=.true.
       else
       yes_sparse=.false.
       endif
     endif


    if((contractionj.ne.0.or.iessz).and.yes_sparse.and..not.yes_sparse_choose) then
     if(rank.eq.0.and.contractionj.ne.0) write(6,*) ' Warning  yes_sparse=.true. implemented only with  uncontracted Jastrow '
     if(rank.eq.0.and.iessz) write(6,*) ' Warning  yes_sparse=.true. not implemented with  old spin  Jastrow/ use -22 or similar '
      yes_sparse=.false.
    endif

    nnozeromax = max(nnozero + nnozero_eagp, nnozeroj)
    if(yes_sparse) nnozeromax=max(nnozeromax,2*nnozeroj)

    ALLOCATE(nozeron(nnozeromax))
    nozeron = 0
    if(yes_sparse) then
    ALLOCATE(nozeroj_c(max(2*nnozeroj_c, 1)))
    else
    ALLOCATE(nozeroj_c(max(nnozeroj_c, 1)))
    endif
    nozeroj_c = 0

    !         calculation upper bound for maxparamj
    maxparamj = 1
    do i = 1, nshellj_c
        if(ioptorbj_c(i).le.900000) then
            ind = 0
            do j = 1, nshellj_c
                if(kionj_c(j).eq.kionj_c(i)                                 &
                        &      .and.ioptorbj_c(j).le.900000) ind = ind + 1
            enddo
            if(ind.gt.maxparamj) maxparamj = ind
        endif
    enddo
    maxparamj = maxparamj + molecularj

    if(numcost.gt.maxparamj) maxparamj = numcost

    occj_c = occ

    do i = 2, nshellj_c
        if(kionj_c(i).lt.kionj_c(i - 1)) then
            if(rank.eq.0) then
                write(6, *) 'Warn. Ions are not in ascending order in Jas !'
                write(6, *) 'Check shell Jas =', i
            endif
        endif
    enddo

    if(indpar.ne.npar3body_c) &
            call error('read_fort10', ' Some error in input Jastrow 3-body ', 1, rank)

    if(rank.eq.0) then
        if(contractionj.ne.0) then
            write(6, *) ' USING CONTRACTED JASTROW ORBITALS ', contractionj
        else
            write(6, *) ' USING UNCONTRACTED JASTROW ORBITALS '
        endif
    endif

    if(contractionj.ne.0) then

        if(rank.eq.0) then
            write(6, *) 'USING', contractionj, 'CONTRACTED + CONST JASTROW ORBITALS'
            write(6, *) 'USING', molecularj, 'MOLECULAR  JASTROW ORBITALS'
            write(6, *) 'USING', maxshellj, 'max shells jastrow'
            write(6, *) 'USING', maxshell, 'max shells det '
        endif

    else

        ! orbital coeff root basis (zeta): dupr (det) --> vjur (jastrow)

        nshelljmax = max(nshellj, 1)

        ALLOCATE(vjur(max(npar3body, 1)), vjurold(0:max(npar3body, 1)))
        ALLOCATE(multj(nshelljmax), kionj(nshelljmax), ioptorbj(nshelljmax)&
                &, nparamj(nshelljmax))
        vjurold = 0.d0
        vjur = 0.d0
        multj = 0
        kionj = 0
        ioptorbj = 0
        nparamj = 0
        if(yes_sparse) then 
        allocate(nozeroj(max(2*nnozeroj, 1)),nozerojder(max(nnozeroj, 1)))
        else
        allocate(nozeroj(max(nnozeroj, 1)),nozerojder(max(nnozeroj, 1)))
        endif

        ALLOCATE(scalej(max(nnozeroj, 1)))
        ALLOCATE(scalejsz(max(nnozeroj, 1)))

        scalej = 0.d0
        nozeroj = 0
        nozerojder = 0
        scalejsz = 0.d0
        ! in this case root and contracted must coincide
        do i = 1, npar3body_c
            vjur(i) = vju_c(i)
        enddo
        do i = 1, nshellj_c
            kionj(i) = kionj_c(i)
            multj(i) = multj_c(i)
            ioptorbj(i) = ioptorbj_c(i)
            nparamj(i) = nparamj_c(i)
        enddo

        occtotj = occ

    endif


    ! read occupation of determinantal orbitals
    if(rankrep.eq.0) then
        read(ufort10, *, err = 134, end = 134)
    endif

    ALLOCATE(ioccup_c(occ_c))
    ioccup_c = 0

    if(rankrep.eq.0) then
        if(rank.eq.0) write(*, *) ' Reading det occupation . . . . '
        do i = 1, occ_c
            read(ufort10, *, err = 136, end = 136) ioccup_c(i)
        enddo
    endif
#ifdef PARALLEL
  call bcast_integer(ioccup_c,occ_c,0,commrep_mpi)
#ifdef UNREL
  !   For unreliable  networks.
  call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
#endif

    nelorb_c = 0

    ind = 0
    do i = 1, nshell_c
        do j = 1, mult_c(i)
            ind = ind + 1
            if(ioccup_c(ind).eq.1.and.ioptorb_c(i).ne.1000000) nelorb_c = nelorb_c + ipf
            if(ioccup_c(ind).eq.1.and.ioptorb_c(i).eq.1000000) nelorb_c = nelorb_c + 1
        enddo
    enddo

    ! do i=1,occ_c
    !    if(ioccup_c(i).eq.1.and.ioptorb_c(i).ne.1000000) nelorb_c=nelorb_c+ipf
    !    if(ioptorb_c(i).eq.1000000) nelorb_c=nelorb_c+1
    ! enddo

    if(rank.eq.0) then
        write(6, *) 'Number of total det orbitals', occ_c
        write(6, *) 'Number of occupied det orbitals', nelorb_c
    endif

    if(contraction.ne.0) then

        !            nshell=nshell_c*maxparam ! max value for nshell (root)
        ! max value for nshell (root)
        nshell = nshell_max


        !            write(6,*) ' nshell max =',nshell
        !            nelorbmax=maxshell*nshell ! max value for ioccup (root)

        dimtranspip = iesup_c
        do ii = 2, maxshell
            dimtranspip = dimtranspip + iesup_atom
        enddo

        if(rank.eq.0) then
            write(6, *) ' Preliminary allocation ', 24d-9 * iesup_c
            write(6, *) ' transpip =', 4d-9 * dimtranspip
        endif
        memtot = memtot + 3 * iesup_c + 0.5 * dimtranspip

        if(rankrep.eq.0) then
            ALLOCATE(occshell(nshell), mu_touch(nelorbmax, occ_c))
            if(yes_complex) then
                ALLOCATE(mu_tmp(2 * nelorbmax, occ_c))
            else
                ALLOCATE(mu_tmp(nelorbmax, occ_c))
            endif
            mu_tmp = 0.d0
            occshell = 0
            mu_touch = 0.d0
        endif
        allocate(transpip(maxshell))
        allocate(transpip(1)%col(iesup_c))
        transpip(1)%col = 0
        do i1 = 2, maxshell
            allocate(transpip(i1)%col(iesup_atom))
            transpip(i1)%col = 0
        enddo

        ALLOCATE(mult(nshell), ioptorb(nshell), nparam(nshell), kion(nshell)  &
                &, ioccup(nelorbmax), ioccdo(nelorbmax)                       &
                &, iesuptrans(nshell))
        ALLOCATE(duprold(nshell), dupr(nshell))

        !             Inizialize
        duprold = 0.d0
        dupr = 0.d0
        mult = 0
        ioptorb = 0
        nparam = 0
        kion = 0
        ioccup = 0
        ioccdo = 0
        iesuptrans = 0

        if(rankrep.eq.0) then

            mu_touch = 0
            mu_tmp = 0.d0

            iesuptrans = 0
            iesuptransb = 0
            duprold = 0
            kion = 0
            mult = 0
            ioptorb = 0
            do i1 = 1, maxshell
                transpip(i1)%col = 0
            enddo
            multranspip = 0
            nshell = 0
            indpar = 0
            occ = 0
            indocc = 0

            do i = 1, nshell_c
                if(ioptorb_c(i).ne.ioptorbcontr(ioptorb_c(i), LBox, 1)&
                        &.and.ioptorb_c(i).lt.900000) then

                    do j = 1, nparam_c(i) / 2

                        !look for inequivalent shells
                        jj = 1
                        if(yes_complex) then
                            do while((dup_c(2 * (indpar + j) - 1).ne.duprold(jj).or.      &
                                    &                kion_c(i).ne.kion(jj).or.                         &
                                    &ioptorbcontr(ioptorb_c(i), LBox, j).ne.ioptorb(jj))        &
                                    &                .and.jj.le.nshell)
                                jj = jj + 1
                            enddo
                        else
                            do while((dup_c(indpar + j).ne.duprold(jj).or.      &
                                    &                kion_c(i).ne.kion(jj).or.                         &
                                    &ioptorbcontr(ioptorb_c(i), LBox, j).ne.ioptorb(jj))        &
                                    &                .and.jj.le.nshell)
                                jj = jj + 1
                            enddo
                        endif
                        ind = indpar + nparam_c(i) / 2 + j

                        if(jj.gt.nshell) then

                            nshell = nshell + 1
                            if(yes_complex) then
                                duprold(nshell) = dup_c(2 * (indpar + j) - 1)
                            else
                                duprold(nshell) = dup_c(indpar + j)
                            endif
                            ioptorb(nshell) = ioptorbcontr(ioptorb_c(i), LBox, j)
                            mult(nshell) = mult_c(i)
                            kion(nshell) = kion_c(i)
                            iesuptrans(nshell) = indpar + j
                            occshell(nshell) = occ
                            iesuptransb(indpar + j) = nshell
                            multranspip(ind) = 0

                            do kk = 1, mult(nshell)
                                if(ioccup_c(indocc + kk).ne.0) then
                                    multranspip(ind) = multranspip(ind) + 1
                                    transpip(multranspip(ind))%col(ind) = &
                                            &                      (indocc + kk - 1) * nelorbmax + occ + kk
                                endif
                                ioccup(occ + kk) = 1
                                if(yes_complex) then
                                    mu_tmp(2 * (occ + kk) - 1, indocc + kk) = dup_c(2 * ind - 1)         &
                                            & + mu_tmp(2 * (occ + kk) - 1, indocc + kk)
                                    mu_tmp(2 * (occ + kk), indocc + kk) = dup_c(2 * ind)         &
                                            & + mu_tmp(2 * (occ + kk), indocc + kk)
                                else
                                    mu_tmp(occ + kk, indocc + kk) = dup_c(ind)         &
                                            & + mu_tmp(occ + kk, indocc + kk)
                                endif
                                mu_touch(occ + kk, indocc + kk) = 1
                            enddo
                            occ = occ + mult(nshell)
                        else
                            iesuptransb(indpar + j) = jj
                            occ_tmp = occshell(jj)
                            multranspip(ind) = 0
                            do kk = 1, mult(jj)
                                if(ioccup_c(indocc + kk).ne.0) then
                                    multranspip(ind) = multranspip(ind) + 1
                                    transpip(multranspip(ind))%col(ind) = &
                                            &                      (indocc + kk - 1) * nelorbmax + occ_tmp + kk
                                endif
                                ioccup(occ_tmp + kk) = 1
                                if(yes_complex) then
                                    mu_tmp(2 * (occ_tmp + kk) - 1, indocc + kk) = dup_c(2 * ind - 1)     &
                                            & + mu_tmp(2 * (occ_tmp + kk) - 1, indocc + kk)
                                    mu_tmp(2 * (occ_tmp + kk), indocc + kk) = dup_c(2 * ind)     &
                                            & + mu_tmp(2 * (occ_tmp + kk), indocc + kk)
                                else
                                    mu_tmp(occ_tmp + kk, indocc + kk) = dup_c(ind)     &
                                            & + mu_tmp(occ_tmp + kk, indocc + kk)
                                endif
                                mu_touch(occ_tmp + kk, indocc + kk) = 1
                            enddo
                        endif

                    enddo

                    indpar = indpar + nparam_c(i)
                    indocc = indocc + mult_c(i)

                elseif(ioptorb_c(i).lt.900000) then

                    if(nparam_c(i).gt.1) then
                        if(rankrep.eq.0) write(6, *) 'Uncontracted orbital with', nparam_c(i), 'zeta!', rankcolrep
                        call error('read_fort10', ' It is supposed to be a single zeta orbital! ', 1, rank)
                    endif

                    jj = 1
                    if(yes_complex) then
                        do while((dup_c(2 * indpar + 1).ne.duprold(jj).or.              &
                                &             kion_c(i).ne.kion(jj).or.                &
                                &    ioptorbcontr(ioptorb_c(i), LBox, j).ne.ioptorb(jj)) &
                                &             .and.jj.le.nshell)
                            jj = jj + 1
                        enddo
                    else
                        do while((dup_c(indpar + 1).ne.duprold(jj).or.                &
                                &             kion_c(i).ne.kion(jj).or.                &
                                &    ioptorbcontr(ioptorb_c(i), LBox, j).ne.ioptorb(jj)) &
                                &             .and.jj.le.nshell)
                            jj = jj + 1
                        enddo
                    endif

                    if(jj.gt.nshell) then
                        nshell = nshell + 1
                        if(yes_complex) then
                            duprold(nshell) = dup_c(2 * indpar + 1)
                        else
                            duprold(nshell) = dup_c(indpar + 1)
                        endif
                        ioptorb(nshell) = ioptorbcontr(ioptorb_c(i), LBox, j)
                        mult(nshell) = mult_c(i)
                        kion(nshell) = kion_c(i)
                        iesuptrans(nshell) = indpar + 1
                        iesuptransb(indpar + 1) = nshell
                        occshell(nshell) = occ
                        if(yes_complex) then
                            do kk = 1, mult(nshell)
                                ioccup(occ + kk) = 1
                                mu_tmp(2 * (occ + kk) - 1, indocc + kk) = 1.d0
                                mu_tmp(2 * (occ + kk), indocc + kk) = 0.d0
                                mu_touch(occ + kk, indocc + kk) = 1
                            enddo
                        else
                            do kk = 1, mult(nshell)
                                ioccup(occ + kk) = 1
                                mu_tmp(occ + kk, indocc + kk) = 1.d0
                                mu_touch(occ + kk, indocc + kk) = 1
                            enddo
                        endif
                        occ = occ + mult(nshell)
                    else
                        iesuptransb(indpar + 1) = jj
                        occ_tmp = occshell(jj)
                        if(yes_complex) then
                            do kk = 1, mult(jj)
                                ioccup(occ_tmp + kk) = 1
                                mu_tmp(2 * (occ_tmp + kk) - 1, indocc + kk) = 1.d0
                                mu_tmp(2 * (occ_tmp + kk), indocc + kk) = 0.d0
                                mu_touch(occ_tmp + kk, indocc + kk) = 1
                            enddo
                        else
                            do kk = 1, mult(jj)
                                ioccup(occ_tmp + kk) = 1
                                mu_tmp(occ_tmp + kk, indocc + kk) = 1.d0
                                mu_touch(occ_tmp + kk, indocc + kk) = 1
                            enddo
                        endif
                    endif

                    indpar = indpar + nparam_c(i)
                    indocc = indocc + mult_c(i)

                endif

            enddo
        endif  ! endif rank=0


#ifdef PARALLEL
     call bcast_integer(nshell,1,0,commrep_mpi)
     call bcast_integer(iesuptrans,nshell,0,commrep_mpi)
     call bcast_integer(iesuptransb,iesupr_c,0,commrep_mpi)
     call bcast_real(duprold,size(duprold),0,commrep_mpi)
     call bcast_real(dupr,size(dupr),0,commrep_mpi)
     call bcast_integer(transpip(1)%col,iesup_c,0,commrep_mpi)
     call bcast_integer(multranspip,iesup_c,0,commrep_mpi)
     call bcast_integer(kion,nshell,0,commrep_mpi)
     call bcast_integer(mult,nshell,0,commrep_mpi)
!    call mpi_bcast(kion,nshell,MPI_INTEGER,0,commrep_mpi,ierr)
!    call mpi_bcast(mult,nshell,MPI_INTEGER,0,commrep_mpi,ierr)
     call bcast_integer(nparam,nshell,0,commrep_mpi)
     call bcast_integer(ioptorb,nshell,0,commrep_mpi)
     do i1=2,maxshell
     call bcast_integer(transpip(i1)%col,iesup_atom,0,commrep_mpi)
     enddo
     call bcast_integer(ioccup,nelorbmax,0,commrep_mpi)
     call bcast_integer(ioccdo,nelorbmax,0,commrep_mpi)
     call bcast_integer(indpar,1,0,commrep_mpi)
     call bcast_integer(occ,1,0,commrep_mpi)
#endif

        nelorb = 0
        do i = 1, occ
            if(ioccup(i).ne.0) nelorb = nelorb + 1
        enddo

        if(rank.eq.0) then
            write(6, *) 'Number of total det orbitals (root)', occ
            write(6, *) 'Number of occupied det orbitals (root)', nelorb
        endif

        nelorbh = nelorb

        nelcolh = ipf * nelorbh + ndiff
        iesupr = nshell
        ! true value of nshell
        nshellr = nshell
        occtot = occ
        do i = 1, nshell
            nparam(i) = 1
            dupr(i) = duprold(i)
        enddo
        !    It should be no longer used
        deallocate(duprold)
        indpar = nshell


        !if(rank.eq.0) then
        !   do i=1,nshell
        !      write(*,*) 'shell,iopt,mult,dupr',i,ioptorb(i),mult(i),dupr(i)
        !   enddo
        !endif

        iesupr_2 = iesupr

        do i = 1, occtot
            ioccdo(i) = ioccup(i)
        enddo
        nshelldo = nshell

    else ! contraction>0
        nelorbh = nelorb_c / ipf
        nelcolh = ipf * nelorbh + ndiff
        nelorb = nelorbh

        ALLOCATE(ioccup(occtot), ioccdo(occtot))
        ioccup = 0
        ioccdo = 0
        do i = 1, occ_c
            ioccup(i) = ioccup_c(i)
            if(ioccup(i).eq.0) then
                call error(' read_fort10 ', 'Holes  in  the occupation no longer supported with uncontracted', 1, rank)
            endif
        enddo
        do i = 1, occtot
            ioccdo(i) = ioccup(i)
        enddo
    endif ! endif contraction

    ndiffdim = ipc * ndiff * nelorbh * ipf
    indndiff = ipc * ipf * ipf * nelorbh * nelorbh + 1


    ! Initialization of direct lattice vectors for building the
    ! periodic basis set. The definition is the same employed in the
    ! Crystal DFT code. This basis can be used for both real and complex
    ! wave functions. In the case of complex wave function it can be used
    ! for an open system too. Use the keyword "PBC_C" in the first line
    ! of the wave function to enforce this option.
    ! NB: Changing epsbas only once after read_datas
    if(epsbas.gt.0.d0) lepsbas = -log(epsbas)

    if(.not.iespbc) then 
    yes_scemama=yes_scemama_open
    endif

    ! read occupation of jastrow orbitals
    if(rankrep.eq.0) then
        read(ufort10, *, err = 138, end = 138)
    endif

    ALLOCATE(ioccj_c(max(occj_c, 1)))
    ioccj_c = 0

    if(rankrep.eq.0) then
        if(rank.eq.0) write(6, *) ' Reading jas occupation . . . . '
        do i = 1, occj_c
            read(ufort10, *, err = 140, end = 140) ioccj_c(i)
        enddo
    endif
#ifdef PARALLEL
  call bcast_integer(ioccj_c,occj_c,0,commrep_mpi)
#ifdef UNREL
  !   For unreliable  networks.
  call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
#endif

    nelorbj_c = 0
    do i = 1, occj_c
        if(ioccj_c(i).eq.1) nelorbj_c = nelorbj_c + 1
    enddo

    if(rank.eq.0) then
        write(6, *) 'Number of total Jas orbitals', occj_c
        write(6, *) 'Number of occupied Jas orbitals', nelorbj_c
    endif

    ! write(6,*) ' Lbox inside fort10',LBox
    if(LBox.eq.3.and..not.yes_crystalj) then
        LBoxj = 1.d0
    else
        LBoxj = LBox
    endif

    if(contractionj.ne.0) then

        nshellj = nshellj_max
        !            nshellj=nshellj_c*maxparamj ! max value for nshellj (root)
        !            nelorbmaxj=maxshellj*nshellj ! max value for ioccj (root)
        ALLOCATE(muj_tmp(nelorbmaxj, occj_c), occshellj(nshellj)     &
                &, muj_touch(nelorbmaxj, occj_c))

        allocate(transpipj(maxshellj))
        do i1 = 1, maxshellj
            allocate(transpipj(i1)%col(npar3body_c))
            transpipj(i1)%col = 0
        enddo

        muj_touch = 0
        muj_tmp = 0.d0
        occshellj = 0

        nshelljmax = max(nshellj, 1)
        ALLOCATE(vjurold(0:nshellj), vjutouch(nshellj), vjur(nshellj)       &
                &, multj(nshelljmax), ioptorbj(nshelljmax), nparamj(nshelljmax)&
                &, kionj(nshelljmax) &
                &, ioccj(max(nelorbmaxj, 1)), iesuptransj(nshelljmax))

        iesuptransj = 0
        iesuptransbj = 0
        vjurold = 0.d0
        vjutouch = 0
        kionj = 0
        multj = 0
        nparamj = 0
        ioptorbj = 0
        do i1 = 1, maxshellj
            transpipj(i1)%col = 0
        enddo
        multranspipj = 0

        nshellj = 0
        npar3body_fill = 0
        indpar = 0
        occj = 0
        indocc = 0
        do i = 1, nshellj_c
            if(ioptorbj_c(i).ne.ioptorbcontr(ioptorbj_c(i), LBoxj, 1)&
                    &.and.ioptorbj_c(i).lt.900000) then

                do j = 1, nparamj_c(i) / 2

                    !look for inequivalent shells
                    jj = 1
                    do while((kionj_c(i).ne.kionj(jj).or.             &
                            &                vju_c(indpar + j).ne.vjurold(vjutouch(jj)).or.      &
                            &     ioptorbcontr(ioptorbj_c(i), LBoxj, j).ne.ioptorbj(jj)) &
                            &.and.jj.le.nshellj)
                        jj = jj + 1
                    enddo

                    ind = indpar + nparamj_c(i) / 2 + j

                    if(jj.gt.nshellj) then
                        nshellj = nshellj + 1
                        ioptorbj(nshellj) = ioptorbcontr(ioptorbj_c(i), LBoxj, j)
                        multj(nshellj) = multj_c(i)
                        kionj(nshellj) = kionj_c(i)
                        nparamj(nshellj) = 1
                        occshellj(nshellj) = occj

                        npar3body_fill = npar3body_fill + 1
                        vjurold(npar3body_fill) = vju_c(indpar + j)
                        vjutouch(nshellj) = npar3body_fill
                        ! nshellj can be .ne. npar3body_fill
                        iesuptransj(npar3body_fill) = indpar + j
                        iesuptransbj(indpar + j) = npar3body_fill

                        multranspipj(ind) = 0

                        do kk = 1, multj(nshellj)
                            if(ioccj_c(indocc + kk).ne.0) then
                                multranspipj(ind) = multranspipj(ind) + 1
                                transpipj(multranspipj(ind))%col(ind) = &
                                        &                      (indocc + kk - 1) * nelorbmaxj + occj + kk
                            endif
                            ioccj(occj + kk) = 1
                            muj_tmp(occj + kk, indocc + kk) = vju_c(ind)       &
                                    & + muj_tmp(occj + kk, indocc + kk)
                            muj_touch(occj + kk, indocc + kk) = 1
                        enddo
                        occj = occj + multj(nshellj)

                    else
                        iesuptransbj(indpar + j) = vjutouch(jj)
                        occ_tmp = occshellj(jj)
                        multranspipj(ind) = 0
                        do kk = 1, multj(jj)
                            if(ioccj_c(indocc + kk).ne.0) then
                                multranspipj(ind) = multranspipj(ind) + 1
                                transpipj(multranspipj(ind))%col(ind) = &
                                        &                      (indocc + kk - 1) * nelorbmaxj + occ_tmp + kk
                            endif
                            ioccj(occ_tmp + kk) = 1
                            muj_tmp(occ_tmp + kk, indocc + kk) = vju_c(ind)    &
                                    & + muj_tmp(occ_tmp + kk, indocc + kk)
                            muj_touch(occ_tmp + kk, indocc + kk) = 1
                        enddo
                    endif

                enddo

                indpar = indpar + nparamj_c(i)
                indocc = indocc + multj_c(i)

            elseif(ioptorbj_c(i).lt.900000) then

                if(nparamj_c(i).gt.1) then
                    if(rankrep.eq.0) write(6, *) 'Uncontracted Jas orbital with', nparamj_c(i), 'zeta!', rankcolrep
                    call error('read_fort10', ' It is supposed to be a single zeta orbital! ', 1, rank)
                endif

                jj = 1
                if(ioptorbcontr(ioptorbj_c(i), LBoxj, j).eq.200) then
                    ! the constant is always the constant independent of kion

                    do while(ioptorbcontr(ioptorbj_c(i), LBoxj, j).ne.ioptorbj(jj)&
                            &.and.jj.le.nshellj)
                        jj = jj + 1
                    enddo
                else
                    do while(((kionj_c(i).ne.kionj(jj).or.&
                            &ioptorbcontr(ioptorbj_c(i), LBoxj, j).ne.ioptorbj(jj).or.&
                            &vju_c(indpar + 1).ne.vjurold(vjutouch(jj))))&
                            &.and.jj.le.nshellj)
                        jj = jj + 1
                    enddo
                endif

                if(jj.gt.nshellj) then
                    nshellj = nshellj + 1
                    ioptorbj(nshellj) = ioptorbcontr(ioptorbj_c(i), LBoxj, j)
                    multj(nshellj) = multj_c(i)
                    kionj(nshellj) = kionj_c(i)
                    occshellj(nshellj) = occj
                    ! the constant does not have any parameter vju_c
                    if(ioptorbj(nshellj).ne.200) then
                        !otherwise it is zero (const
                        nparamj(nshellj) = 1
                        npar3body_fill = npar3body_fill + 1
                        vjurold(npar3body_fill) = vju_c(indpar + 1)
                        vjutouch(nshellj) = npar3body_fill
                        iesuptransj(npar3body_fill) = indpar + 1
                        iesuptransbj(indpar + 1) = npar3body_fill
                    endif
                    do kk = 1, multj(nshellj)
                        !                if(ioccj_c(indocc+kk).ne.0) then
                        ioccj(occj + kk) = 1
                        !                endif
                        muj_tmp(occj + kk, indocc + kk) = 1.d0
                        muj_touch(occj + kk, indocc + kk) = 1
                    enddo
                    occj = occj + multj(nshellj)
                else
                    if(ioptorbj(jj).ne.200)                           &
                            &                iesuptransbj(indpar + 1) = vjutouch(jj)
                    occ_tmp = occshellj(jj)
                    do kk = 1, multj(jj)
                        !                        if(ioccj_c(indocc+kk).ne.0) then
                        ioccj(occ_tmp + kk) = 1
                        !                        endif
                        muj_tmp(occ_tmp + kk, indocc + kk) = 1.d0
                        muj_touch(occ_tmp + kk, indocc + kk) = 1
                    enddo
                endif

                indpar = indpar + nparamj_c(i)
                indocc = indocc + multj_c(i)

            endif

        enddo

        nelorbj = 0
        do i = 1, occj
            if(ioccj(i).ne.0) nelorbj = nelorbj + 1
        enddo

        if(rank.eq.0) then
            write(6, *) 'Number of total Jas orbitals (root)', occj
            write(6, *) 'Number of occupied Jas orbitals (root)', nelorbj
        endif

        nelorbjh = nelorbj
        npar3body = npar3body_fill
        ! true value of nshell
        nshelljr = nshellj
        occtotj = occj

        do i = 1, npar3body
            vjur(i) = vjurold(i)
        enddo
        indpar = npar3body

        !      do i=1,nshell
        !      write(*,*) 'shell,iopt,mult,dup',i,ioptorb(i),mult(i),dupr(i)
        !      enddo
        !      stop

        npar3body_2 = npar3body

    else ! contractionj>0

        nelorbjh = nelorbj_c
        nelorbj = nelorbjh

        ALLOCATE(ioccj(max(occtotj, 1)))
        ioccj = 0
        do i = 1, occj_c
            ioccj(i) = ioccj_c(i)
        enddo

    endif

    allocate(ioptorbja(nelorbj))
    if(nelorbj.gt.0) ioptorbja = 0

    indorb = 0
    ind = 0
    do i = 1, nshellj
        do j = 1, multj(i)
            if(ioccj(ind + j).ne.0) then
                indorb = indorb + 1
                ioptorbja(indorb) = ioptorbj(i)
            endif
        enddo
        ind = ind + multj(i)
    enddo

    allocate(orbcostl(nelorbj))
    if(nelorbj.gt.0) orbcostl = .true.
    do i = 1, nelorbj
        if(ioptorbja(i).eq.200) then
            orbcostl(i) = .true.
        else
            orbcostl(i) = .false.
        endif
    enddo

    dimjas = 0
    do ii = 1, nshelljr
        if(ioptorbj(ii).ne.200) then
            dimjas = dimjas + 1
        endif
    enddo

    iscramax = nwm + 15 * npm
    if(itestr.eq.-5) then
        ! For bconstraint
        if(symmagp.and.ipc.eq.2) then
            iscramax = max(iscramax, iessw + 4 * nnozero_c)
        else
            iscramax = max(iscramax, iessw + ipc * nnozero_c)
        endif
    else
        if(symmagp.and.ipc.eq.2) then
            iscramax = max(iscramax, 4 * nnozero_c)
        else
            iscramax = max(iscramax, ipc * nnozero_c)
        endif
    endif
    iscramax = max(iscramax, 2 * nnozeroj_c)
    iscramax = max(iscramax, 2 * nelorbj)
    iscramax = max(iscramax, ipc * nelup_mat * nelup_mat + 4 * nelup_mat + (ipf - 1) * (ipc * nelup_mat * nelup_mat + 8 * ipc * nelup_mat + 4 * (ipc - 1) * nelup_mat)) ! needed for psireg and ipf=2 dsktrf and dgtsvx
    iscramax = max(iscramax, max(iesupr, ipc * iesupr_c) + iesupr) ! needed for transpsi

    if(.not.yesdft) then
        !        memory required for upscratch reduce at best
        !    if(iessz) then

        !       iscramax=max(iscramax,3*nel*nel+nelorbj*ip4*nel+(indt+2)*20&
        !            &,(indt+4+ip4)*nelorbj+(indt+2)*20&
        !            &,ipc*(indt+ip4)*nelorb+(indt+2)*20,ipc*nelup*nelup+nbdgetri*nelup)
        !    else
        !       iscramax=max(iscramax,2*nel*nel+nelorbj*ip4*nel+(indt+2)*20&
        !            &,(indt+ip4+4)*nelorbj,(indt+ip4)*nelorb*ipc+(indt+2)*20&
        !            &,ipc*nelup*nelup+nbdgetri*nelup)
        !    endif
        iscramax = max(iscramax, 20 * (indt + 1))  ! task4
        iscramax = max(iscramax, (nion * (nion - 1)) / 2) ! task3
        iscramax = max(iscramax, nparshellmax) ! task2 nparshellmax defined  after read_pseudo -->
        ! read_fort10 after read_pseudo

        !    task5 iscramax=max(iscramax,2*ipc*nelup*nmolfn) ! task5 after definition of nmolfn< nelorb_c
        iscramax = max(iscramax, 2 * ipc * nelup_mat * nelorb_c) ! task5 upper bound
        iscramax = max(iscramax, 2 * nelorbj_c * nel)   ! task5 Jastrow
        iscramax = max(iscramax, ipc * nelup_mat * nelup_mat + max(ipc * nbdgetri, 4 * nelup_mat)) ! task8 and psireg and scratchdet

        !    the true memory required in  task9
        if(iessz) then
            iscramax = max(iscramax, 2 * max(nel, indt + 5) + 2 * max(nelorbjh, 1)&
                    & + max(nelorbjh, 1) * (2 * indt + 12) + 27 * (indt + 1) * max(nshellj,nion))+nel+nion
        else
            iscramax = max(iscramax, 2 * max(nel, indt + 5) + max(ipj * nelorbjh, 1)&
                    & + max(nelorbjh, 1) * (2 * indt + 12) + 27 * (indt + 1) * max(nshellj,nion))+nel+nion
        endif
        iscramax = max(iscramax, ipc * (indt + 5) * nelorbh + nelorbh * (indt + 6) + 27 * (indt + 1) * max(nshell,nion) + nelorbh) ! task10
    endif

    dim_upwf = 27 * (indt + 1) * max(nshell,nion) + (indt + 5) * nelorbh
    dim_upwf = max(dim_upwf, 27 * (indt + 1) * max(nshellj,nion) + (indt + 5) * nelorbjh)
    dim_ratiovar = max(dim_upwf, ipc * 6 * nelup_mat)
    dim_uptabtot = max(dim_upwf, ipc * nmol * (5 + 2 * (1 + indt + ip4)))

    if(contraction.eq.0.and.yesfast.ne.0) then
        if(rank.eq.0) write(6, *) ' Warning without contracted orbitals yesfast=0, forced '
        yesfast = 0
    endif

    maxall = max(in1 * (nweight - ibinit), npm)
    maxall = max(maxall, npmn)
    maxall = max(maxall, 3 * ieskin)

    iscrapip = 1

    if(itestr.eq.-5) then
        lwork=npmn*(2+ILAENV( 1, 'DSYTRD', 'L', npmn, -1, -1, -1))  ! for dsyev
        lwork=max(lwork,4*npmn) ! for dgeev

        iscrapip =lwork+ max(npmn,npm)*max(3,2*ipc) ! when is complex yes_hessc=.true.
        iscrapip=max(iscrapip,6*max(npm,npmn)) ! for dgelscut
        iscrapip=max(iscrapip,nbin+3*max(npm,npmn)) ! for dgemv in reweight 
         iscramax = max(maxall + 4 * max(npm, npmn) + 13 * npmn, iscramax)
         iscramax = max(iscramax, ipc * iesup_c + ipc * iesupind, npar3body_c + iesmind) ! memory required for updating contracted orbitals
        if(itest.ne.2.and.nelorbj.gt.0) then ! memory required for uptabpip

            iscramax = max(iscramax, 2 * (indt + ip4 + 1) * nel)
        endif
        if(itest.ne.2) then ! memory required for upinvhop_fnf
            iscramax = max(iscramax, 5 * nelorb * ipc * ipf + 2 * (indt + ip4 + 1) * nelup_mat * ipc)
        endif
        iscramax = max(iscramax, ipc * iesup_c + ipc * iesupind, npar3body_c + iesmind) ! memory required for updating contracted orbitals

        iscramax = max(iscramax, nel * ipc + nel) ! memory required for updiag
!       iscramax = max(iscramax, 2 * npmn * npmn + 14 * npm + 10)
        iscramax = max(iscramax, nparshellmax) ! the correct memory required for pseudo
        !        iscramax=max(iscramax,2*npmn*npm)
        !        iscramax=max(iscramax,3*npmn*ieskin)
    endif

    nelcol = ipf * nelorbh + ndiff
    nelcol_c = nelorb_c + ndiff

    !        memory required for the scontraction
    !       definition nelorb_at
    if(contraction.ne.0) then
        nelorb_at = nelorb_c - molecular
        ndim_detmat = nelorb_c
    else
        nelorb_at = ipf * nelorbh
        ndim_detmat = nelorb_at
    endif

    if(molecular.gt.0.and.itestr.eq.-5) then
        iscramax = max(iscramax, ipc * nelorb_at * (nelorb_at + ndiff)) ! convertmol_c
        iscramax = max(iscramax, ipc * nnozero_c)  ! for symmetrize_agp in convertmol_c
        iscramax = max(iscramax, ipc * nelorb_at * nmolmat)  ! for projectmat in convertmol_c
        iscramax = max(iscramax, ipc * (indt + 4) * nelorbh + nelorbh * (indt + 4) + 20 * (indt + 1) + nelorbh) ! upnewwf
    endif


    !  in iscramax the memory required for convertmol_fast
    if(itestr.eq.-5) iscramax = max(iscramax, nelorb_c) ! gram-Scmidt orthogonalization

    if(itestr.eq.-5) iscramax = max(iscramax, 3 * np) ! memory needed for bootforce later
    !   iscramax=max(iscramax,np) !  for project
    iscramax = max(iscramax, in1) ! memory needed for branching later

!   iscramax = max(iscramax, iscrapip)
    nelorbp = ipj*max(nelorbjh,nelorbj_c) + 1
    nelorbpp = ipj*max(nelorbjh,nelorbj_c)+nelorbp
    if(yesdft) then
        iscramax = max(2 * nnozero_c * ipc + 1, 9 * nion + 1, iesupr + 1, 2 * nnozeroj_c + 1, nelorbjh * nelorbjh) ! replace a smaller allocation
    elseif(itest.ne.2) then
        iscramax = max(iscramax, nelorbpp + nelorbh * ipc * (indt + 5))
        iscramax = max(iscramax, nelorbpp + nel*(3*indt+14)) ! uptabpip
    endif

    !-------- GENERAL ALLOCATION --------!

    if(rank.eq.0) then
        write(6, *) ' Before allocation standard '
        write(6, *) ' nion = ', nion
        write(6, *) ' nel = ', nel
        write(6, *) ' indt=', indt
        write(6, *) ' npm= ', npm
        write(6, *) ' nwm= ', nwm
        write(6, *) ' nws= ', nws
        write(6, *) ' nelorbj= ', nelorbj
        write(6, *) ' iscrapip= ', iscrapip
        write(6, *) ' iscraipsip= ', iscraipsip
        write(6, *) ' nelorb= ', nelorb
        write(6, *) ' nelcol= ', nelcol
        write(6, *) ' nelup= ', nelup
        write(6, *) ' nelorb_c= ', nelorb_c
        write(6, *) ' nelorbj_c= ', nelorbj_c
        write(6, *) ' nelcol_c= ', nelcol_c
        write(6, *) ' nshell= ', nshell
        write(6, *) ' nshellj= ', nshellj
        write(6, *) ' npsamax= ', npsamax
        write(6, *) ' nintpseudo= ', nintpseudo
    endif


    nelorbjh2=2*nelorbj
    dim_jasmat=1
    if(yes_sparse) then
    dim_jasmat=nnozeroj
    elseif(contractionj.eq.0) then
    dim_jasmat=(ipj*nelorbj)**2
    endif
    dim_jasmat=max(dim_jasmat,1)
    if(yesdft) then
     
        ALLOCATE(iond(nion * nion)                                      &
                &, rmu(3 * nion * (max(indt, 1) + 1))               &
                &, r(nion * (max(indt, 1) + 1)), econf(in1 * npm), econfh(in1 * npm)       &
                &, wconfn(nwm), econfion(ieskindim * in1), factorsr(nws), vcut(nws)      &
                &, wsto(nws), zeta(nwm + 1), tabpip((indt + ip4) * nws), yescut(nws)       &
                &, table(indt * nws + nws), diffkin(3, nws)                &
                &, diag(nws), winv(nelorb * (indt + 5) * nws), ainv(nws)    &
                &, psiln(nws), psidetln(nws), psisn(nws), singdet(nws), vpot(nws), vpotreg(2, nws)   &
                &, enertrue(nws), diffuse(nws), winvdo((indt + ip4) * nws), tabler(indt * nws + nws)&
                &, winvup((indt + ip4) * nws), winvj(max(nelorbj, 1) * (indt + 5) * nws)   &
                &, winvbar(ipf * nelorbh * 2 * nws), winvjbar(max(ipj * nelorbjh, 1) * nws + nws)      &
                &, wint(nws), dist(nws * nion), tmu(max(indt, 1) * nws)    &
                &, diagfn(nws), enert(1, nws), berry_exp(2 * nws), wconfsav(nws)              &
                &, winvsj(max(nelorbj, 1) * (indt + 1 + ip4))                &
                &, psinew(nelorbh * (indt + 2 + ip4))&
                &, ainvup(nelorbh), ainvupb(nelorb), ainvs(2)&
                &, ainvdo(nelorbh)                                            &
                &, dists(nion), dists_kel(3, nion)            &
                &, jasmat(dim_jasmat), winvs(nelorb * (ip4 + 1) + 1)           &
                &, cnorm(nshell + nshellj), enerint(nws)  &
                &, iond_cart(3, nion * nion))

        ! initialize  ALL
        iond = 0.d0
        rmu = 0.d0
        r = 0.d0
        econf = 0.d0
        econfh = 0.d0
        wconfn = 0.d0
        econfion = 0.d0
        factorsr = 0.d0
        vcut = 0.d0
        wsto = 0.d0
        zeta = 0.d0
        tabpip = 0.d0
        yescut = .false.
        table = 0.d0
        diffkin = 0.d0
        diag = 0.d0
        winv = 0.d0
        ainv = 0.d0
        psiln = 0.d0
        psidetln = 0.d0
        psisn = 0
        singdet = .true.
        vpot = 0.d0
        vpotreg = 0.d0
        enertrue = 0.d0
        diffuse = 0.d0
        winvdo = 0.d0
        tabler = 0.d0
        winvup = 0.d0
        winvj = 0.d0
        winvbar = 0.d0
        winvjbar = 0.d0
        wint = 0.d0
        dist = 0.d0
        tmu = 0.d0
        diagfn = 0.d0
        enert = 0.d0
        berry_exp = 0.d0
        wconfsav = 0.d0
        winvsj = 0.d0
        psinew = 0.d0
        ainvup = 0.d0
        ainvupb = 0.d0
        ainvs = 0.d0
        ainvdo = 0.d0
        dists = 0.d0
        dists_kel = 0.d0
        jasmat = 0.d0
        winvs = 0.d0
        cnorm = 0.d0
        enerint = 0.d0
        iond_cart = 0.d0

        memtot = memtot + 2 * nion * nion + 3 * nion * max(indt, 1) + 1 + nion * (max(indt, 1) + 1) + &
                &2 * in1 * npm + nwm + ieskindim * in1 + 2 * nws + nws + nwm + 1 + (indt + ip4) * nws + nws + &
                &nwm + indt * nws + nws + 3 * nws + nws + nelorb * (indt + 12) * nws + (indt + ip4) * nws + &
                &indt * nws + nws + (indt + ip4) * nws + max(nelorbj, 1) * (indt + ip4 + 1) * nws + nelorbh * 2 * nws + &
                &2 * nws + nws * nion + (indt + ip4) * nws + 3 * nws + iscramax + nelorbj * (indt + ip4 + 1) + 1 + &
                &nelorb * (indt + ip4) + nelorbh + nelorb + 2 + nelorbh + 4 * nion + max(nelorbjh * nelorbjh, 1) + &
                &nelorb * (ip4 + 1) + 1 + nshell + nshellj + 2 * nws + indt + ip4 + 3 * nion * nion

        ! scratch vectors needed for LBox.eq.3.d0 option
        if(abs(LBox).eq.3.d0) then
            allocate(rmunew(3, 0:indt,nion), npip(0:indt,3, nion), rknew(0:indt,nion))
            npip = 0
            rmunew = 0.d0
            rknew = 0.d0
        endif

        if(rank.eq.0) write(6, *) ' Memory required DFT =', memtot * 8 / 1d9, 'Gbyte'
        if(rank.eq.0) write(6, *) ' Memory winv =', 8d-9 * nelorb * (indt + 5) * nws

    else
        ! in1 = nw/nproc = nws
        ! npm=6 by default (# of correlation functions to be computed)
        ALLOCATE(iond(nion * nion)                                      &
                &, rmu(3 * max(nel,nion) * nion * (max(indt, 1) + 1))               &
                &, r(nel * nion * (max(indt, 1) + 1))       &
                &, wconfn(nwm), econfion(ieskindim * in1), factorsr(nws), vcut(nws)      &
                &, zeta(nwm + 1), tabpip((indt + ip4) * nel * nws), yescut(nws)       &
                &, psiln(nws), psidetln(nws), vpot(nws), vpotreg(2, nel * nws)   &
                &, diffuse(nws), wint(nws)&
                &, winvj(max(nelorbj, 1) * nel * (indt4j + 1) * nws)   &
                &, winvjbar(max(ipj * nelorbjh, 1) * nel * nws + nws)      &
                &, dist(nel * nws * nion), tmu(max(indt, 1) * nel * nws)  &
                &, diagfn(nws), berry_exp(2 * nws), wconfsav(nws)              &
                &, winvsj(max(nelorbj, 1) * (indt + 1 + ip4))                &
                &, dists(nion), dists_kel(3, nion)            &
                &, jasmat(dim_jasmat)              &
                &, cnorm(nshell + nshellj), enerint(nws) &
                &, iond_cart(3, nion * nion), tabler(nel * indt * nws + nws))
        !------- variables necessary for complex algorithms -----!
        if(.not.yes_complex) then
            ! wave function
            ALLOCATE(table(nel * indt * nws + nws)                                            &
                    &, ainv(nelup_mat**2 * nws), winv(nelorb * nel * (indt4 + 1) * nws)                 &
                    &, winvs(nelorb * (indt4 + 1)), winvdo((indt + ip4) * neldomax * nws)          &
                    &, winvup((indt + ip4) * nelup * nws), winvbar(ipf * nelorbh * nel_mat * nws), &
                    & psisn(nws), psinew(ipf * nelorbh * (indt + 2 + ip4)), singdet(nws) &
                    &, ainvup(nelup * nelorbh), ainvupb(nelorb * nelup)            &
                    &, ainvdo(neldomax * nelorbh)                                     &
                    &, diffkin(3, nws), econf(in1 * npm), econfh(in1 * npm), diag(nws)              &
                    &, enertrue(nws), enert(1, nws), wsto(nws))
        else
            ! wave function evaluation
            ALLOCATE(table(2 * nel * indt * nws + 2 * nws)                                       &
                    &, ainv(2 * nelup_mat * nelup_mat * nws)        &
                    &, winvdo(2 * (indt + ip4) * neldomax * nws)     &
                    &, winv(2 * nelorb * nel * (indt4 + 1) * nws), winvup(2 * (indt + ip4) * nelup * nws)   &
                    &, winvs(2 * nelorb * (indt4 + 1)), winvbar(ipf * nelorbh * 2 * nel_mat * nws)           &
                    &, psisn(nws), singdet(nws), psinew(2 * ipf * nelorbh * (indt + 2 + ip4))                         &
                    &, ainvup(2 * nelup * nelorbh), ainvupb(2 * nelorb * nelup)     &
                    &, ainvdo(neldomax * nelorbh)&
                    &, diffkin(3, nws), econf(in1 * npm), econfh(in1 * npm), diag(nws)     &
                    &, enertrue(nws), enert(2, nws), wsto(nws))
        endif
        allocate(ainvs(max(ipc * nelup_mat * nelup_mat, 2 * ipc * nelup_mat)))
        if(epscuttype.gt.0) then 
        allocate(agp(ipc * nelup_mat, nelup_mat, nws),&
                & agpn(ipc * nelup_mat, nelup_mat))
        else
        allocate(agp(1, 1, 1),agpn(1,1))
        endif
        agp=0.d0
        agpn=0.d0
        ! scratch vectors needed for LBox.eq.3.d0 option
        if(abs(LBox).eq.3.d0) then
            allocate(rmunew(3, 0:indt,nion), npip(0:indt,3, nion), rknew(0:indt,nion))
            allocate(rmunewb(3, 0:indt,nion), rknewb(0:indt,nion))
            npip = 0
            rmunew = 0.d0
            rknew = 0.d0
            rmunewb = 0.d0
            rknewb = 0.d0
        endif

        !------- end complex allocation ------!

        iond = 0.d0
        rmu = 0.d0
        r = 0.d0
        econf = 0.d0
        econfh = 0.d0
        wconfn = 0.d0
        econfion = 0.d0
        factorsr = 0.d0
        vcut = 0.d0
        wsto = 0.d0
        zeta = 0.d0
        tabpip = 0.d0
        yescut = .false.
        table = 0.d0
        diffkin = 0.d0
        diag = 0.d0
        winv = 0.d0
        ainv = 0.d0
        if(allocated(agp)) then
            agp = 0.d0
            agpn = 0.d0
        endif
        psiln = 0.d0
        psidetln = 0.d0
        psisn = 0
        singdet = .true.
        vpot = 0.d0
        vpotreg = 0.d0
        enertrue = 0.d0
        diffuse = 0.d0
        winvdo = 0.d0
        tabler = 0.d0
        winvup = 0.d0
        winvj = 0.d0
        winvbar = 0.d0
        winvjbar = 0.d0
        wint = 0.d0
        dist = 0.d0
        tmu = 0.d0
        diagfn = 0.d0
        enert = 0.d0
        berry_exp = 0.d0
        wconfsav = 0.d0
        winvsj = 0.d0
        psinew = 0.d0
        ainvup = 0.d0
        ainvupb = 0.d0
        ainvs = 0.d0
        ainvdo = 0.d0
        dists = 0.d0
        dists_kel = 0.d0
        jasmat = 0.d0
        winvs = 0.d0
        cnorm = 0.d0
        enerint = 0.d0
        iond_cart = 0.d0

        memtot = memtot + 2 * nion * nion + 3 * nel * nion * max(indt, 1) + 1 + nel * nion * (max(indt, 1) + 1) + &
                &2 * in1 * npm + nwm + ieskindim * in1 + 2 * nws + nws + nwm + 1 + (indt + ip4) * nel * nws + nws + &
                &nwm + nel * indt * nws + nws + 3 * nws + nws + nelorb * nel * (indt4 + 2) * nws + nelup * nelup + &
                &(indt + ip4) * nws + indt * nws + nws + (indt + ip4) * neldo * nws + &
                &nel * max(nelorbj, 1) * (indt4j + 1) * nws + nelorbh * nelup * nws + nelorbjh * nel * nws + &
                &2 * nws + nws * nion + (indt + ip4) * nws + 3 * nws + iscramax + nelorbj * (indt + ip4 + 1) + 1 + &
                &nelorb * (indt + ip4) + nelup * nelorbh + nelup * nelorb + 2 + nelorbh + 4 * nion + &
                &max(nelorbjh * nelorbjh, 1) + &
                &nelorb * (ip4 + 1) + 1 + nshell + nshellj + 2 * nws + indt + 4 + 3 * nion * nion
        if(rank.eq.0) write(6, *) ' Memory required QMC =', memtot * 8 / 1d9, 'Gbyte'
        if(rank.eq.1) write(6, *) ' Memory required QMC slaves =', memtot * 8 / 1d9, 'Gbyte'
        if(rank.eq.0) write(6, *) ' Memory winv =', 8d-9 * nelorb * nel * (indt4 + 1) * nws

    endif

    !   Define yesdetmat & yesdetmatc
    if(contraction.eq.0) then
        if(yes_complex) then
            allocate(detmat(2 * ipf * nelorbh * nelcolh))
        else
            allocate(detmat(nelorbh * ipf * nelcolh))
        endif
        yesdetmatc = .true.
        if(contraction.eq.0) then
            yesdetmat = .true.
        else
            yesdetmat = .false.
        endif
    else
        if(yes_complex) then
            allocate(detmat(2 * ipf * nelorbh))
        else
            allocate(detmat(ipf * nelorbh))
        endif
        yesdetmat = .false.
        yesdetmatc = .false.
    endif

    !    memory required for scontract_mat
    if(contraction.gt.0.and.yesdetmatc) then
        if(yes_complex) then
!           iscrapip = max(iscrapip, iessw + 2 * nelcol * nelcol_c)
            iscramax = max(iscramax, iessw + 2 * nelcol * nelcol_c)
        else
!           iscrapip = max(iscrapip, nelcol * nelcol_c)
            iscramax = max(iscramax, nelcol * nelcol_c)
        endif
    endif

    if(contractionj.gt.0) then
!       iscrapip = max(iscrapip, nelorbj * nelorbj_c)
        iscramax = max(iscramax, nelorbj * nelorbj_c)
    endif
    iscramax = iscramax + 3 * nel ! in upscratch_global psip is called with address 3*nel+1
    if(rank.eq.0)  write(6, *) ' iscramax for master = ', iscramax
    if(rank.eq.1)  write(6, *) ' iscramax for slaves = ', iscramax


    !    iscramax=1000000
    allocate(psip(iscramax))

    psip = 0.d0
    !if(rank.eq.0)  write(6,*) ' iscramax for master = ',iscramax
    !if(rank.eq.1)  write(6,*) ' iscramax for slaves = ',iscramax

    jj = abs(iesinv) + iesm + iesd + abs(iesfree) + abs(iessw) + iesup

    !        if(jj.ne.0) then  ! only with minimization required
    !         nmaxder=(nelup+neldo)
    !   allocate(derleft(nmaxder*nelcol),derright(nmaxder*nelorb))
    !         if(contraction.eq.0) then
    !         allocate(derdet(nelorb*nelcol))
    !         else
    !         allocate(derdet(max(nelcol_c*nelcol,nelorb*nelcol)))
    !         endif
    !        allocate(derjas(nelorbj*nelorbj+1)&
    !    &,derdet_mu(nelorb*nelorb_c),derdet_c(nelorb_c*max(nmaxder,nelcol_c))&
    !    &,derjas_mu(nelorbj*nelorbj_c),derjas_c(nelorbj_c*nelorbj_c))
    !        endif

    diffuse = 0.d0
    winvjbar = 0.d0

    !************ CAFFAREL FORCES ALLOCATION **********

    !call dscalzero(iscramax,0.d0,psip,1)
    !psip=0.d0

    allocate(orbcostn(max(ipj * nelorbj_c, 1)))

    if(nelorbj_c.gt.0) orbcostn = .true.

    if(iessz) then
        ALLOCATE (orbcost(max(nelorbj_c, 1))&
                &, jasmatsz(max(nelorbjh * nelorbjh, 1)))
    else
        allocate(jasmatsz(1), orbcost(1))
    endif

    if(iessz.and..not.yesdft) then
        ALLOCATE (winvjbarsz(max(nelorbjh * nel * nws, 1)))
    else
        ALLOCATE (winvjbarsz(nws))
    endif

    winvjbarsz = 0.d0
    orbcost = .true.
    jasmatsz = 0.d0
    detmat = 0


    !    if(.not.yesdft) then
    !       if(itestrr.eq.-4) then
    !          if(yes_complex) then
    !             ALLOCATE(ainvdob(2*nelcol*nelup))
    !          else
    !             ALLOCATE(ainvdob(nelcol*nelup))
    !          endif
    !          ainvdob=0.d0
    !       endif
    !    endif

    ! read non zero elements of lambda and put in scale
    if(rankrep.eq.0) then
        read(ufort10, *, err = 142, end = 142)
        if(rank.eq.0) write(6, *) ' Reading det nnozero . . . . '
        do i = 1, nnozero_c + nnozero_eagp
            if(yes_complex) then
                read(ufort10, *, err = 142, end = 142) ix, iy, scale_c(2 * i - 1), scale_c(2 * i)
            else
                read(ufort10, *, err = 142, end = 142) ix, iy, scale_c(i)
            endif

            if(i.le.nnozero_c) then
                nozero_c(i) = (iy - 1) * nelorb_c + ix
                if(ix.gt.nelorb_c.or.iy.gt.nelcol_c) then
                    call error(' read_fort10 ', 'Raw/Col index too large in determinant!', 1, rank)
                endif
            else
                nozero_c(i) = (iy - 1) * ndiff + ix
                if(ix.gt.ndiff.or.iy.gt.ndiff) then
                    call error(' read_fort10 ', 'Raw index too large in Ghost det', 1, rank)
                endif
            endif
        enddo
    endif
#ifdef PARALLEL
     call bcast_integer(nozero_c,size(nozero_c),0,commrep_mpi)
     call bcast_real(scale_c,size(scale_c),0,commrep_mpi)
#endif



    ! read independent lambda
    if(rankrep.eq.0) then
        read(ufort10, *, err = 144, end = 144)
    endif
    nnozero_all = nnozero_c + nnozero_eagp
    jbradet = 0
    sjbradet = .false.
    nozeron = 0
    jbradetn = 0
    iessw0 = 0
    ind = 0
    indn = 0
    indnn = 0
    indp = 0
    indtot = 0
    iimax = 0
    if(rank.eq.0) write(*, *) ' Reading det nnozero symmetries ....'
    if(rankrep.eq.0) then
        do i = 1, iesswr
            read(ufort10, *, err = 146, end = 146) ii, (ipsip(j), j = 1, 2 * abs(ii))
            !#ifdef PARALLEL
            !        call bcast_integer(ii,1,0,commrep_mpi)
            !#ifdef UNREL
            !        !   For unreliable  networks.
            !        call mpi_barrier(commrep_mpi,ierr)
            !        !$omp barrier
            !#endif
            !        call bcast_integer(ipsip,2*abs(ii),0,commrep_mpi)
            !#ifdef UNREL
            !        !   For unreliable  networks.
            !        call mpi_barrier(commrep_mpi,ierr)
            !        !$omp barrier
            !#endif
            !#endif

            if(abs(ii).gt.iimax) iimax = abs(ii)

            if(ii.gt.0) then
                indp = indp + 1
                indtot = indtot + ii

                do j = 1, ii
                    indn = indn + 1
                    nozeron(indn) = nelorb_c * (abs(ipsip(2 * j)) - 1) + abs(ipsip(2 * j - 1))
                    if(ipsip(2 * j - 1).gt.0) then
                        jbradet(indn) = indp
                    elseif(ipsip(2 * j - 1).lt.0) then
                        jbradet(indn) = -indp
                    endif
                    if(ipsip(2 * j).lt.0) then
                        sjbradet(indn) = .true.
                    else
                        sjbradet(indn) = .false.
                    endif
                enddo
            else
                do j = 1, -ii
                    indn = indn + 1
                    nozeron(indn) = -(nelorb_c * (abs(ipsip(2 * j)) - 1) + abs(ipsip(2 * j - 1)))
                enddo
                iessw0 = iessw0 + 1
                indnn = indnn + 1
                indtot = indtot - ii
                jbradetn(indnn) = ii
                do j = 1, -2 * ii
                    indnn = indnn + 1
                    jbradetn(indnn) = ipsip(j)
                enddo
            endif
        enddo
        do i = 1, iesswr_eagp
            read(ufort10, *, err = 146, end = 146) ii, (ipsip(j), j = 1, 2 * abs(ii))
            if(abs(ii).gt.iimax) iimax = abs(ii)
            if(ii.gt.0) then
                indp = indp + 1
                indtot = indtot + ii
                do j = 1, ii
                    indn = indn + 1
                    nozeron(indn) = ndiff * (abs(ipsip(2 * j)) - 1) + abs(ipsip(2 * j - 1))
                    if(ipsip(1).gt.0) then
                        jbradet(indn) = indp
                    elseif(ipsip(1).lt.0) then
                        jbradet(indn) = -indp
                    endif
                enddo
            else
                do j = 1, -ii
                    indn = indn + 1
                    nozeron(indn) = -(ndiff * (abs(ipsip(2 * j)) - 1) + abs(ipsip(2 * j - 1)))
                enddo
                iessw0 = iessw0 + 1
                indnn = indnn + 1
                indtot = indtot - ii
                jbradetn(indnn) = ii
                do j = 1, -2 * ii
                    indnn = indnn + 1
                    jbradetn(indnn) = ipsip(j)
                enddo
            endif
        enddo
    endif  ! endif rankrep


#ifdef PARALLEL
        call bcast_integer(indp,1,0,commrep_mpi)
        call bcast_integer(indn,1,0,commrep_mpi)
        call bcast_integer(indtot,1,0,commrep_mpi)
        call bcast_integer(indnn,1,0,commrep_mpi)
        call bcast_integer(iessw0,1,0,commrep_mpi)
        call bcast_integer(iimax,1,0,commrep_mpi)
        call bcast_integer(nozeron,size(nozeron),0,commrep_mpi)
        call bcast_integer(jbradet,size(jbradet),0,commrep_mpi)
        call bcast_integer(jbradetn,size(jbradetn),0,commrep_mpi)
        call mpi_bcast(sjbradet,size(sjbradet),MPI_LOGICAL,0,commrep_mpi,ierr)
#endif

    iesswind = indp
    ind = indn

    if(iessw.ne.0.and.iesswind.ne.iessw) then
        if(rank.eq.0) write(6, *) ' Warning iessw changed to ', iesswind * ipc
        iessw = iesswind * ipc
        if(symmagp.and.ipc.eq.2.and.yes_correct) iessw = iessw * 2
    endif

    if(rank.eq.0) then
        write(6, *) 'Number of non zero geminal lambda for det ', indtot
        write(6, *) 'Number of non fixed geminal lambda for det', indn
    endif

    call checkmatrix(nnozero_c, indn - nnozero_eagp, nelorb_c, nozero_c, nozeron&
            &, ipsip, rank, nelcol_c)

    if(npar_eagp.gt.0) then
        call checkmatrix(nnozero_eagp, nnozero_eagp, ndiff, nozero_c(nnozero_c + 1)&
                &, nozeron(nnozero_c + 1), ipsip(nnozero_c + 1), rank, ndiff)
        do i = nnozero_c + 1, nnozero_all
            if(ipsip(i).gt.0) then
                ipsip(i) = ipsip(i) + nnozero_c
            elseif(ipsip(i).lt.0) then
                ipsip(i) = ipsip(i) - nnozero_c
            endif
        enddo
    endif


    !     if(rank.eq.0) then
    !     write(6,*) ' Mapping '
    !     do i=1,nnozero_c
    !     if(ipsip(i).ne.0) then
    !     write(6,*) i,ipsip(i),nozero_c(i),nozeron(abs(ipsip(i)))
    !     endif
    !     enddo
    !     endif

    ipsip(nnozero_all + 1:2 * nnozero_all) = jbradet(1:nnozero_all)

    !     The pointer to the optimizable parameters of the AGP-Pfaff matrix
    jbradet(1:nnozero_all) = 0
    do i = 1, nnozero_all
        if(ipsip(i).gt.0) then
            jbradet(i) = ipsip(nnozero_all + ipsip(i))
            !       if(rank.eq.0) write(6,*) i,jbradet(i)
        endif
    enddo
    !     the same sorting for sjbradet storing the original in ipsip
    do i = 1, nnozero_c
        if(sjbradet(i)) then
            ipsip(nnozero_all + i) = 1
        else
            ipsip(nnozero_all + i) = 0
        endif
    enddo
    sjbradet = .false.
    do i = 1, nnozero_c
        if(ipsip(i).gt.0) then
            if(ipsip(nnozero_all + ipsip(i)).eq.1) sjbradet(i) = .true.
        endif
    enddo
    !     nozeron(1:nnozero_c)=nozero_c(1:nnozero_c)
    !     psip(1:ipc*nnozero_c)=scale_c(1:ipc*nnozero_c)
    !     nozeron(1:indn)=abs(nozeron(1:indn))
    call checkrepeat(indn - nnozero_eagp, nelorb_c, nozeron, rank)
    if(npar_eagp.gt.0) then
        call checkrepeat(nnozero_eagp, ndiff, nozeron(indn - nnozero_eagp + 1), rank)
    endif
    !    do i=1,nnozero_c
    !       nozero_c(i)=nozeron(i)
    !       if(yes_complex) then
    !          scale_c(2*i)=psip(2*i)
    !          scale_c(2*i-1)=psip(2*i-1)
    !       else
    !          scale_c(i)=psip(i)
    !       endif
    !    enddo

    if(indtot.ne.nnozero_c + nnozero_eagp) then
        if(rank.eq.0) then
            write(6, *) ' Inconsistent input sym in Det. ', indtot, nnozero_c
            write(6, *) ' The # of written table entries =', nnozero_c
            write(6, *) ' The # of symmetric table entries =', indtot
        endif
        call error('read_fort10', ' Error checking matrices! ', 1, rank)
    elseif(indn.ne.nnozero_c + nnozero_eagp) then
        if(rank.eq.0) then
            write(6, *) ' The # of written table entries =', nnozero_c + nnozero_eagp
            write(6, *) ' The # of matches  =', indn
        endif
        call error('read_fort10', ' Error checking matrices! ', 1, rank)
    endif

    if(symmagp.and.ipc.eq.1.and.ipf.eq.1) then
        nmol = molecular - ndiff
    else
        !       if(ipf.ne.2) then
        nmol = (molecular - ndiff) / 2
        !       else
        !       nmol=2*(molecular-ndiff)
        !       endif
    endif

    occopt = .false.
    call update_kiontot
    if(npar_eagp.ne.0) then
        allocate(eagp_pfaff(ipc * ndiff, ndiff))
    else
        allocate(eagp_pfaff(1,1))
    endif
    eagp_pfaff = 0.d0
    if(contraction.ne.0) then
        ! default_allocation set yesdft to true
        ! in order to optimize memory allocation for DFT and tools
        if(yesdft) then
            if(yes_complex) then
                ALLOCATE(detmat_c(2 * nelorb_c * max(nelcol_c, nel)))
                memtot = memtot + 2 * nelorb_c * max(nelcol_c, nel)
            else
                ALLOCATE(detmat_c(nelorb_c * max(nelcol_c, nel)))
                memtot = memtot + nelorb_c * max(nelcol_c, nel)
            endif
        else
            ! allocate geminal matrix in atomic contracted basis
            if(.not.yes_complex) then
                allocate(detmat_c(nelorb_c * max(nelcol_c, nel))&
                        &, projm(nelorbh * ipf * nelcol_c))
            else
                allocate(detmat_c(2 * nelorb_c * max(nelcol_c, nel))&
                        &, projm(ipf * 2 * nelorbh * nelcol_c))
            endif
            !
            projm = 0.d0
            if(yes_complex) then
                memtot = memtot + 2 * nelorb_c * max(nelcol_c, nel)
            else
                memtot = memtot + nelorb_c * max(nelcol_c, nel)
            endif
        endif ! endif yesdft

        if(yes_complex) then

            call dscalzero(2 * nelorb_c * nelcol_c, 0.d0, detmat_c, 1)
            !

            do i = 1, nnozero_c
                call upsim_complex(detmat_c, nelorb_c, nozero_c(i), scale_c(2 * i - 1), symmagp, ipf)
            enddo
        else
            call dscalzero(nelorb_c * nelcol_c, 0.d0, detmat_c, 1)
            do i = 1, nnozero_c
                call upsim(detmat_c, nelorb_c, nozero_c(i), scale_c(i), symmagp, ipf)
            enddo
        endif
        if(npar_eagp.ne.0) then
            eagp_pfaff = 0.d0
            do i = 1, nnozero_eagp
                ind = nnozero_c + i
                iy = (nozero_c(ind) - 1) / ndiff + 1
                ix = nozero_c(ind) - (iy - 1) * ndiff
                if(ipc.eq.1) then
                    eagp_pfaff(ix, iy) = scale_c(ind)
                    eagp_pfaff(iy, ix) = -scale_c(ind)
                else
                    eagp_pfaff(2 * ix - 1, iy) = scale_c(2 * ind - 1)
                    eagp_pfaff(2 * ix, iy) = scale_c(2 * ind)
                    eagp_pfaff(2 * iy - 1, ix) = -scale_c(2 * ind - 1)
                    eagp_pfaff(2 * iy, ix) = -scale_c(2 * ind)
                endif
            enddo
        endif



        !          Here symmetrize agp for detecting the conventional phase
        !          if not molyes and manyfort10
        !          Here is the generic read in start.
        !            if(manyfort10.and.(.not.molyes.or.write_effective)) then
        !           call  symmetrizeagp_ref(nnozero_c,nozero_c,jbradet,sjbradet&
        !          &,jbradetn,dsw,iessw0,psip,ipsip,detmat_c,nelorb_c,nelorb_at&
        !          &,nelcol_c,symmagp,yes_hermite)
        !            endif

        if(allowed_averagek) then
            !          From the effective to the real one
            call attach_phase2det(.true., detmat_c)
            if(rank.eq.0) write(6, *) ' Passi qui from eff. to real I '
        endif
    else
        ! no contraction
        if(yes_complex) then
            ALLOCATE(detmat_c(2), mu_c(2 * ipf * nelorbh, 1), projm(2 * ipf * nelorbh))
        else
            ALLOCATE(detmat_c(1), mu_c(ipf * nelorbh, 1), projm(ipf * nelorbh))
        endif
        detmat_c = 0.d0
        mu_c = 0.d0
        projm = 0.d0
    endif

    !  to be optimized
    if(contraction.ne.0) then
        ! allocation of matrix of contracted coefficients mu_c
        if(allocated(mu_c)) deallocate(mu_c)
        !
        if(yes_complex) then
            ALLOCATE(mu_c(2 * ipf * nelorbh, nelorb_c))
            mu_c = 0.d0
        else
            ALLOCATE(mu_c(ipf * nelorbh, nelorb_c))
            mu_c = 0.d0
        endif
        !
        memtot = memtot + nelorbh * nelorb_c

        if(rankrep.eq.0) then
            if(allocated(addr_occ)) deallocate(addr_occ)
            allocate(addr_occ(occ_c))
            allocate(transpip_sav(maxshell))
            allocate(transpip_sav(1)%col(iesup_c))
            do i1 = 2, maxshell
                allocate(transpip_sav(i1)%col(iesup_atom))
            enddo

            transpip_sav = transpip
            ! shift the mu matrix according to ioccup
            ix = 0

            do i = 1, occ
                if(ioccup(i).ne.0) then
                    ix = ix + 1
                    iy = 0
                    do j = 1, occ_c
                        if(ioccup_c(j).ne.0) then
                            iy = iy + 1
                            if(mu_touch(i, j).ne.0) then
                                if(ipf.eq.2) then
                                    if(yes_complex) then
                                        mu_c(2 * ix - 1, iy) = mu_c(2 * ix - 1, iy) + mu_tmp(2 * i - 1, j)
                                        mu_c(2 * ix, iy) = mu_c(2 * ix, iy) + mu_tmp(2 * i, j)
                                        mu_c(2 * ix - 1 + 2 * nelorbh, iy + nelorb_at / 2) = &
                                                &       mu_c(2 * ix - 1 + 2 * nelorbh, iy + nelorb_at / 2) + mu_tmp(2 * i - 1, j)
                                        mu_c(2 * ix + 2 * nelorbh, iy + nelorb_at / 2) = &
                                                &       mu_c(2 * ix + 2 * nelorbh, iy + nelorb_at / 2) + mu_tmp(2 * i, j)
                                    else
                                        mu_c(ix, iy) = mu_c(ix, iy) + mu_tmp(i, j)
                                        mu_c(ix + nelorbh, iy + nelorb_at / 2) = mu_c(ix + nelorbh, iy + nelorb_at / 2) + mu_tmp(i, j)
                                    endif
                                else
                                    if(yes_complex) then
                                        mu_c(2 * ix - 1, iy) = mu_c(2 * ix - 1, iy) + mu_tmp(2 * i - 1, j)
                                        mu_c(2 * ix, iy) = mu_c(2 * ix, iy) + mu_tmp(2 * i, j)
                                    else
                                        mu_c(ix, iy) = mu_c(ix, iy) + mu_tmp(i, j)
                                    endif
                                endif
                            endif
                        endif
                    enddo
                endif
            enddo

            iy = 0
            addr_occ = 0
            do j = 1, occ_c
                if(ioccup_c(j).ne.0) then
                    iy = iy + 1
                    addr_occ(j) = iy
                endif
            enddo

            !            Simply put out of the loop

            if(ipf.eq.2) then
                deallocate(transpip)
                maxshell = 2 * maxshell
                allocate(transpip(maxshell))
                allocate(transpip(1)%col(iesup_c))
                transpip(1)%col = 0
                do i1 = 2, maxshell
                    allocate(transpip(i1)%col(iesup_atom))
                    transpip(i1)%col = 0
                enddo
                do kk = 1, iesup_c
                    do jj = 1, multranspip(kk)
                        ii = (transpip_sav(jj)%col(kk) - 1) / nelorbmax + 1
                        ix = transpip_sav(jj)%col(kk) - (ii - 1) * nelorbmax
                        iy = addr_occ(ii)
                        if(iy.eq.0) then
                            write(6, *) ' ERROR contraction coefficient !!! ', iy
                        else
                            transpip(jj)%col(kk) = (iy - 1) * ipf * nelorbh + ix
                            transpip(jj + multranspip(kk))%col(kk) = (iy + nelorb_at / 2 - 1) * ipf * nelorbh + ix + nelorbh
                        endif
                    enddo
                    multranspip(kk) = multranspip(kk) * 2
                enddo
            else
                do kk = 1, iesup_c
                    do jj = 1, multranspip(kk)
                        ii = (transpip_sav(jj)%col(kk) - 1) / nelorbmax + 1
                        ix = transpip_sav(jj)%col(kk) - (ii - 1) * nelorbmax
                        iy = addr_occ(ii)
                        if(iy.eq.0) then
                            write(6, *) ' ERROR contraction coefficient !!! ', iy
                        else
                            transpip(jj)%col(kk) = (iy - 1) * nelorbh + ix
                        endif
                    enddo
                enddo
            endif

            deallocate(addr_occ)
        endif ! endif rank=0

#ifdef PARALLEL
!  It has been changed  for ipf=2
     if(ipf.eq.2.and.rankrep.ne.0) then
     maxshell=2*maxshell
           deallocate(transpip)
           allocate(transpip(maxshell))
           allocate(transpip(1)%col(iesup_c))
           transpip(1)%col=0
           do i1=2,maxshell
           allocate(transpip(i1)%col(iesup_atom))
           transpip(i1)%col=0
           enddo
     endif
   do i1=1,maxshell
   call bcast_integer(transpip(i1)%col(1),size(transpip(i1)%col),0,commrep_mpi)
   enddo
   call bcast_real(mu_c,size(mu_c),0,commrep_mpi)
#endif




        if(molyesp) then
            indorb = 0
            indpar = 0
            ind = 0
            do i = 1, nshell_c
                if(ioptorb_c(i).ge.900000) then
                    ind = ind + 1
                    if(ioccup_c(ind).eq.1) then
                        indorb = indorb + 1
                        ! definition mu_c
                        do kk = 1, nparam_c(i) / 2
                            if(yes_complex) then
                                ix = nint(dup_c(2 * (indpar + kk) - 1))
                            else
                                ix = nint(dup_c(indpar + kk))
                            endif
                            if(ix.gt.nelorbh * ipf.or.ix.lt.1) then
                                if(ioptorb_c(i).eq.900000) then
                                    write(6, *) ' Error the hybrid atomic basis does not match with orbital # ', ix
                                    write(6, *) ' Minimum/Maximum atomic orbital # allowed ', 1, ipf * nelorbh
                                    call error(' read_fort10 ', 'Error setting hybrid orbitals.', 1, rank)
                                else
                                    write(6, *) ' Error the molecular basis does not match with orbital # ', ix
                                    write(6, *) ' Minimum/Maximum atomic orbital # allowed ', 1, ipf * nelorbh
                                    call error(' read_fort10 ', 'Error setting molecular orbitals.', 1, rank)
                                endif
                            endif
                            indn = indpar + kk + nparam_c(i) / 2
                            if(ioptorb_c(i).eq.1000000.and.ipf.eq.2) then
                                if(yes_complex) then
                                    mu_c(2 * ix - 1, indorb + nelorb_at / 2) = dup_c(2 * indn - 1)
                                    mu_c(2 * ix, indorb + nelorb_at / 2) = dup_c(2 * indn)
                                else
                                    mu_c(ix, indorb + nelorb_at / 2) = dup_c(indn)
                                endif
                            elseif(kk.le.nparam_c(i) / 4.or.ipf.eq.1) then
                                if(yes_complex) then
                                    mu_c(2 * ix - 1, indorb) = dup_c(2 * indn - 1)
                                    mu_c(2 * ix, indorb) = dup_c(2 * indn)
                                else
                                    mu_c(ix, indorb) = dup_c(indn)
                                endif
                            else
                                if(yes_complex) then
                                    mu_c(2 * ix - 1, indorb + nelorb_at / 2) = dup_c(2 * indn - 1)
                                    mu_c(2 * ix, indorb + nelorb_at / 2) = dup_c(2 * indn)
                                else
                                    mu_c(ix, indorb + nelorb_at / 2) = dup_c(indn)
                                endif
                            endif
                            !                multpointer(ix)=multpointer(ix)+1
                            !                mupointer(multpointer(ix),ix)=indorb
                            if(ioptorb_c(i).eq.1000000.and.ipf.eq.2) then
                                multranspip(indn) = 1
                                transpip(1)%col(indn) = nelorbh * ipf * (indorb + nelorb_at / 2 - 1) + ix
                            elseif(kk.le.nparam_c(i) / 4.or.ipf.eq.1) then
                                multranspip(indn) = 1
                                transpip(1)%col(indn) = nelorbh * ipf * (indorb - 1) + ix
                            else
                                multranspip(indn) = 1
                                transpip(1)%col(indn) = nelorbh * ipf * (indorb + nelorb_at / 2 - 1) + ix
                            endif
                        enddo
                    endif
                        indpar = indpar + nparam_c(i)

                else
                    do j = 1, mult_c(i)
                        ind = ind + 1
                        if(ioccup_c(ind).eq.1) indorb = indorb + 1  ! For ipf=2 the basis is doubled and indorb is the index for the up spin basis
                    enddo
                    indpar = indpar + nparam_c(i)
                endif
            enddo
        endif

        if(rankrep.eq.0) deallocate(mu_touch, mu_tmp, transpip_sav, occshell)


        ! multranspip: multeplicity of a linear coefficient in the contracted shell
        ! to the angular momentum degeneracy

        ! transpip(mult,ind): ind is the linear coefficient index in the contraction
        ! , mult is the multeplicity of that linear coefficient
        ! , transpip gives the coordinates of the coefficient in the mu_c matrix

        ! mu_c(i,j): i is the index for the root orbital
        ! , j is the index for the contracted orbital
        ! ,mu_c gives the linear coefficient of the contraction

        ! root orbital ix, given the root orbital ix and mult=1,multpointer(ix)

        ! ELIMINATED multpointer & mupointer A STUPID WAY TO DO CONTRACTION

        ! multpointer(ix): ix index of root orbitals
        ! multpointer number of contracted orbitals which include the orbital i
        ! expansion
        ! mupointer(mult,ix): is the index of a contracted orbital which include

        if(yesdetmatc) call scontract_mat_det(nelorbh, nelorbh, nelcolh&
                &, nelorb_c, nelcol_c, detmat, detmat_c, mu_c, psip)

        if(symmagp.or.ipf.eq.2) then
            if(ipf.eq.2) then
                indocc = (2 * nelorbh * (2 * nelorbh - 1)) / 2 + 2 * ndiff * nelorbh
            else
                indocc = (nelorbh * (nelorbh + 1)) / 2 + ndiff * nelorbh
            endif
        else
            indocc = nelorbh * nelorbh + ndiff * nelorbh
        endif
        nnozero = indocc

        if(rank.eq.0) write(6, *) 'Number of non zero lambda (root)', nnozero

        if(yesdetmatc) then
            if(yes_complex) then
                ALLOCATE(scale(2 * nnozero + 2 * nnozero_eagp), nozero(nnozero + nnozero_eagp))
            else
                ALLOCATE(scale(nnozero + nnozero_eagp), nozero(nnozero + nnozero_eagp))
            endif
        else
            ALLOCATE(scale(2), nozero(1))
        endif

        if(contraction.eq.0) then
            allocate(nozerodet(nnozero + nnozero_eagp))
        else
            allocate(nozerodet(1))
        endif
        nozerodet = 0
        scale = 0
        nozero = 0
        if(yesdetmatc) then
            if(symmagp.or.ipf.eq.2) then
                indocc = 0
                do ix = 1, nelorbh * ipf
                    do iy = ix + ipf - 1, nelorbh * ipf + ndiff
                        ind = (iy - 1) * nelorbh * ipf + ix
                        indocc = indocc + 1

                        if(yes_complex) then
                            call dcopy(2, detmat(2 * ind - 1), 1, scale(2 * indocc - 1), 1)
                        else
                            scale(indocc) = detmat(ind)
                        endif
                        if(contraction.eq.0) nozerodet(indocc) = (iy - 1) * nelorb * ipf + ix
                        nozero(indocc) = (iy - 1) * nelorbh * ipf + ix
                    enddo
                enddo
            else
                indocc = 0
                do ix = 1, nelorbh
                    do iy = 1, nelorbh + nelup - neldo
                        ind = (iy - 1) * nelorbh + ix
                        indocc = indocc + 1
                        if(yes_complex) then
                            call dcopy(2, detmat(2 * ind - 1), 1, scale(2 * indocc - 1), 1)
                        else
                            scale(indocc) = detmat(ind)
                        endif
                        if(contraction.eq.0) nozerodet(indocc) = (iy - 1) * nelorb + ix
                        nozero(indocc) = (iy - 1) * nelorbh * ipf + ix
                    enddo
                enddo
            endif  ! endif symmagp

        endif  ! endif yesdetmat

    else  ! contraction.gt.0

        scale = scale_c

        nozerodet = nozero_c
        nozero = nozero_c

        if(yes_complex) then
            call dscalzero(2 * ipf * nelorbh * nelcolh, 0.d0, detmat, 1)
            do i = 1, nnozero
                call upsim_complex(detmat, ipf * nelorbh, nozero(i), scale(2 * i - 1), symmagp, ipf)
            enddo

        else
            call dscalzero(ipf * nelorbh * nelcolh, 0.d0, detmat, 1)
            do i = 1, nnozero
                call upsim(detmat, ipf * nelorbh, nozero(i), scale(i), symmagp, ipf)
            enddo
        endif
        if(npar_eagp.ne.0) then
            eagp_pfaff = 0.d0
            do i = 1, nnozero_eagp
                ind = nnozero_c + i
                iy = (abs(nozero_c(ind)) - 1) / ndiff + 1
                ix = abs(nozero_c(ind)) - (iy - 1) * ndiff
                if(ipc.eq.1) then
                    eagp_pfaff(ix, iy) = scale_c(ind)
                    eagp_pfaff(iy, ix) = -scale_c(ind)
                else
                    eagp_pfaff(2 * ix - 1, iy) = scale_c(2 * ind - 1)
                    eagp_pfaff(2 * ix, iy) = scale_c(2 * ind)
                    eagp_pfaff(2 * iy - 1, ix) = -scale_c(2 * ind - 1)
                    eagp_pfaff(2 * iy, ix) = -scale_c(2 * ind)
                endif
            enddo
        endif

        !          Here symmetrize agp for detecting the conventional phase
        !          if not molyes and manyfort10

        !           if(manyfort10.and.(.not.molyes.or.write_effective)) then

        !           call  symmetrizeagp_ref(nnozero_c,nozero_c,jbradet,sjbradet&
        !          &,jbradetn,dsw,iessw0,psip,ipsip,detmat,nelorb_c,nelorb_at&
        !          &,nelcol_c,symmagp,yes_hermite)
        !           endif
        if(allowed_averagek) then
            !          From the effective to the real one
            call attach_phase2det(.true., detmat)
            if(rank.eq.0) write(6, *) ' Passi qui from eff. to real II '
        endif

    endif  ! endif contraction gt 0

    ! read lambda for jastrow and put in scalej
    if(rankrep.eq.0) then
        read(ufort10, *, err = 148, end = 148)
        if(rank.eq.0) write(*, *) ' Reading jas nnozero ....'
        if(yes_sparse) then
        do i = 1, nnozeroj_c
            read(ufort10, *, err = 148, end = 148) ix, iy, scalej_c(i)
            nozeroj_c(i) =ix
            nozeroj_c(i+nnozeroj_c)=iy
            if(ix.gt.ipj * nelorbj_c) then
                call error(' read_fort10 ', 'Raw index too large in Jastrow!', 1, rank)
            endif
        enddo
        else
        do i = 1, nnozeroj_c
            read(ufort10, *, err = 148, end = 148) ix, iy, scalej_c(i)
            nozeroj_c(i) = (iy - 1) * ipj * nelorbj_c + ix
            if(ix.gt.ipj * nelorbj_c) then
                call error(' read_fort10 ', 'Raw index too large in Jastrow!', 1, rank)
            endif
        enddo
        endif
    endif
#ifdef PARALLEL
#ifdef UNREL
     !   For unreliable  networks.
     call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
     call bcast_integer(nozeroj_c,size(nozeroj_c),0,commrep_mpi)
     call bcast_real(scalej_c,nnozeroj_c,0,commrep_mpi)
#ifdef UNREL
     call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
#endif
     if(iessz) then
    if(rankrep.eq.0) then
        write(*, *) ' Reading jas-sz nnozero ....'
        read(ufort10, *, err = 150, end = 150)
        do i = 1, nnozeroj_c
            iy = (nozeroj_c(i) - 1) / nelorbj_c + 1
            ix = nozeroj_c(i) - (iy - 1) * nelorbj_c
            read(ufort10, *, err = 150, end = 150) ixr, iyr, scalejsz_c(i)
            if(ixr.ne.ix.or.iyr.ne.iy) then
                call error(' read_fort10 ', 'Jastrowsz should be in the same order!', 1, rank)
            endif
        enddo
    endif

#ifdef PARALLEL
        call bcast_real(scalejsz_c,nnozeroj_c,0,commrep_mpi)
#ifdef UNREL
        !   For unreliable  networks.
        call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
#endif

    endif

! Calculation which orbitals are constant ones
! It is supposed that the contraction does not contain constant orbitals
ind = 0
indorb = 0
if(rank.eq.0) write(6, *) 'Location constant orbitals in Jastrow'
do i = 1, nshellj_c
    do j = 1, multj_c(i)
        ind = ind + 1
        if(ioccj_c(ind).eq.1) then
            indorb = indorb + 1

            if(ioptorbj_c(i).eq.200.or.ioptorbj_c(i).eq.199) then
                orbcostn(indorb) = .true.
                if(rank.eq.0) write(6, *) indorb, ' Constant orbital '
            else
                orbcostn(indorb) = .false.
            endif
            if(iessz) then
                if((ioptorbj_c(i).eq.200.or.ioptorbj_c(i).eq.199)&
                        &.and.iesgros) then
                    orbcost(indorb) = .true.
                else
                    orbcost(indorb) = .false.
                endif
            endif
        endif
    enddo
enddo

if(ipj.eq.2) orbcostn(indorb + 1:2 * indorb) = orbcostn(1:indorb)


! read the set of constraints  for the Jastrow Jastrowsz
if(rankrep.eq.0) then
    read(ufort10, *, err = 152, end = 152)
endif
    jbraj = 0
    nozeron = 0
    jbrajsz = 0

indn = 0
indnn = 0
indp = 0
indtot = 0
iesfreesz = 0
iijmax = 0
if(rank.eq.0) write(6, *) ' Reading jas nnozero symmetries ....', rankcolrep
if(rankrep.eq.0) then
    do i = 1, iesfreer
        read(ufort10, *, err = 154, end = 154) ii, (ipsip(j), j = 1, 2 * abs(ii))
        !#ifdef PARALLEL
        !        call bcast_integer(ii,1,0,commrep_mpi)
        !#ifdef UNREL
        !        !   For unreliable  networks.
        !        call mpi_barrier(commrep_mpi,ierr)
        !        !$omp barrier
        !#endif
        !        call bcast_integer(ipsip,2*abs(ii),0,commrep_mpi)
        !#ifdef UNREL
        !        !   For unreliable  networks.
        !        call mpi_barrier(commrep_mpi,ierr)
        !        !$omp barrier
        !#endif
        !#endif
        if(abs(ii).gt.iijmax) iijmax = abs(ii)

        if(ii.gt.0) then
            if(twobodyoff) then
                iescostd = .false.
                do j = 1, 2 * ii
                    if(orbcostn(abs(ipsip(j)))) iescostd = .true.
                enddo
                if(iescostd) indp = indp + 1
            else
                iescostd = .true.
                indp = indp + 1
            endif

            if((onebodysz.or.twobodyoff).and.iessz) then
                !   If there is constant orbital update the spin Jastrow Mat. el.

                iescost = .false.
                do j = 1, 2 * ii
                    if(orbcostn(abs(ipsip(j)))) iescost = .true.
                enddo
                if(iescost)  iesfreesz = iesfreesz + 1

            elseif(.not.onebodysz.and.iessz) then
                !   If there is no constant orbital update the spin Jastrow Mat. el.

                iescost = .true.
                do j = 1, 2 * ii
                    if(orbcost(abs(ipsip(j)))) iescost = .false.
                enddo
                if(iescost)  iesfreesz = iesfreesz + 1

            endif
            indtot = indtot + ii

            if(iescostd) then
                do j = 1, ii
                    indn = indn + 1
                    if(yes_sparse) then
                    nozeron(indn)=abs(ipsip(2 * j - 1))
                    nozeron(indn+nnozeroj)=abs(ipsip(2 * j))
                    else
                    nozeron(indn) = ipj * nelorbj_c * (abs(ipsip(2 * j)) - 1)              &
                            & + abs(ipsip(2 * j - 1))
                    endif
                    if(ipsip(2 * j - 1).gt.0) then
                        jbraj(indn) = indp
                    elseif(ipsip(2 * j - 1).lt.0) then
                        jbraj(indn) = -indp
                    endif

                    if(iescost) then
                        if(ipsip(2 * j - 1).gt.0) then
                            jbrajsz(indn) = iesfreesz
                        elseif(ipsip(2 * j - 1).lt.0) then
                            jbrajsz(indn) = -iesfreesz
                        endif
                    endif
                enddo
            else
                indnn = indnn + 1
                jbrajn(indnn) = ii
                do j = 1, 2 * abs(ii)
                    indnn = indnn + 1
                    jbrajn(indnn) = ipsip(j)
                enddo
                do j = 1, abs(ii)
                    indn = indn + 1
                    if(yes_sparse) then
                    nozeron(indn)=-abs(ipsip(2 * j - 1))
                    nozeron(indn+nnozeroj)=abs(ipsip(2 * j))
                    else 
                    nozeron(indn) = -(ipj * nelorbj_c * (abs(ipsip(2 * j)) - 1)&
                            & + abs(ipsip(2 * j - 1)))
                    endif
                enddo
            endif

        else
            indnn = indnn + 1
            jbrajn(indnn) = ii
            indtot = indtot - ii
            do j = 1, -2 * ii
                indnn = indnn + 1
                jbrajn(indnn) = ipsip(j)
            enddo
            do j = 1, abs(ii)
                indn = indn + 1
                if(yes_sparse) then
                nozeron(indn)=-abs(ipsip(2 * j - 1))
                nozeron(indn+nnozeroj)=abs(ipsip(2 * j))
                else
                nozeron(indn) = -(ipj * nelorbj_c * (abs(ipsip(2 * j)) - 1)&
                        & + abs(ipsip(2 * j - 1)))
                endif
            enddo
        endif
    enddo
endif

#ifdef PARALLEL
        call bcast_integer(indp,1,0,commrep_mpi)
        call bcast_integer(indn,1,0,commrep_mpi)
        call bcast_integer(indtot,1,0,commrep_mpi)
        call bcast_integer(indnn,1,0,commrep_mpi)
        call bcast_integer(iessw0,1,0,commrep_mpi)
        call bcast_integer(iijmax,1,0,commrep_mpi)
        call bcast_integer(iesfreesz,1,0,commrep_mpi)
        call bcast_integer(nozeron,size(nozeron),0,commrep_mpi)
        call bcast_integer(jbraj,size(jbraj),0,commrep_mpi)
        call bcast_integer(jbrajn,size(jbrajn),0,commrep_mpi)
        call bcast_integer(jbrajsz,size(jbrajsz),0,commrep_mpi)
        call mpi_bcast(iescost,1,MPI_LOGICAL,0,commrep_mpi,ierr)
        call mpi_bcast(iescostd,1,MPI_LOGICAL,0,commrep_mpi,ierr)
#endif

if(iesfree.ne.0.and.indp.ne.iesfree) then
    if(rank.eq.0) write(6, *) 'Warning iesfree changed to ', indp
    iesfree = indp
endif

if(iesinv.ne.0.and.iesinv.ne.iesfreesz) then
    if(rank.eq.0) write(6, *) 'Warning iesinv changed to ', iesfreesz
    iesinv = iesfreesz
endif

if(.not.iessz.or.(iesinv.eq.0.and..not.yesdft)) iesfreesz = 0
if(rank.eq.0) then
    write(6, *) 'Number of non zero geminal lambda for Jas', indtot
    write(6, *) 'Number of non fixed geminal lambda for Jas', indn
    write(6, *) ' Number of accepted nnozeron Jas Sz ', iesfreesz
endif

if(yes_sparse) then
call checkmatrix_sparse(nnozeroj_c, indn, ipj * nelorbj_c, nozeroj_c, nozeron&
        &, ipsip, rank, ipj * nelorbj_c)
else
call checkmatrix(nnozeroj_c, indn, ipj * nelorbj_c, nozeroj_c, nozeron&
        &, ipsip, rank, ipj * nelorbj_c)
endif

ipsip(nnozeroj_c + 1:2 * nnozeroj_c) = jbraj(1:nnozeroj_c)
jbraj = 0
do i = 1, nnozeroj_c
    if(ipsip(i).gt.0) then
        jbraj(i) = ipsip(nnozeroj_c + ipsip(i))
    endif
enddo
if(yes_sparse) then
   nnozeroj=0
   if(itestr.eq.-5) then ! in case of optimization the optimizable coeff can turn non zero
   do i=1,nnozeroj_c
     if(jbraj(i).ne.0.or.(jbraj(i).eq.0.and.scalej_c(i).ne.0.d0)) then
     nnozeroj=nnozeroj+1
     nozerojder(nnozeroj)=i
     endif
   enddo
   else
!  Only  the non zero coefficients
   do i=1,nnozeroj_c
     if(scalej_c(i).ne.0.d0) then
     nnozeroj=nnozeroj+1
     nozerojder(nnozeroj)=i
     endif
   enddo
   endif
   if(nnozeroj.gt.0) ipsip(nnozeroj_c+1:nnozeroj_c+nnozeroj)=nozerojder(1:nnozeroj)
   deallocate(jasmat,nozerojder,nozeroj)
!  save a bit of memory
   allocate(jasmat(max(nnozeroj,1)),nozerojder(max(nnozeroj,1)),nozeroj(max(2*nnozeroj,1)))
   jasmat=0.d0
   nozerojder=0
   nozeroj=0
   if(nnozeroj.gt.0) nozerojder(1:nnozeroj)=ipsip(nnozeroj_c+1:nnozeroj_c+nnozeroj)
   do i=1,nnozeroj
   nozeroj(i)=nozeroj_c(nozerojder(i))
   nozeroj(i+nnozeroj)=nozeroj_c(nozerojder(i)+nnozeroj_c)
   jasmat(i)=scalej_c(nozerojder(i))
   enddo
   if(rank.eq.0) write(6,*) ' SPARSE number of jasmat el.=',nnozeroj
endif

if(iessz) then
    ipsip(nnozeroj_c + 1:2 * nnozeroj_c) = jbrajsz(1:nnozeroj_c)
    do i = 1, nnozeroj_c
        if(ipsip(i).gt.0) then
            jbrajsz(i) = ipsip(nnozeroj_c + ipsip(i))
        endif
    enddo
endif

if(.not.yes_sparse) then
if(rank.eq.0) write(6, *) ' Check repeated in the symmetry table Jastrow  '
call checkrepeat(indn, nelorbj_c * ipj, nozeron, rank)
endif


!    do i=1,nnozeroj_c
!       nozeroj_c(i)=nozeron(i)
!       scalej_c(i)=psip(i)
!    enddo

!    if(iessz) then
!       do i=1,nnozeroj_c
!          scalejsz_c(i)=psip(i+nnozeroj_c)
!       enddo
!    endif

if(indtot.ne.nnozeroj_c) then
    if(rank.eq.0) then
        write(6, *) ' Inconsistent input sym in Jas. ', indtot, nnozeroj_c
        write(6, *) ' The # of written table entris =', nnozeroj_c
        write(6, *) ' The # of symmetric table entries =', indtot
    endif
    call error('read_fort10', ' Error checking matrices! ', 1, rank)
elseif(indn.ne.nnozeroj_c) then
    if(rank.eq.0) then
        write(6, *) ' The # of written table entries =', nnozeroj_c
        write(6, *) ' The # of matches  =', indn
    endif
    call error('read_fort10', ' Error checking matrices! ', 1, rank)
endif

if(contractionj.ne.0) then
    ALLOCATE(jasmat_c(max(ipj * ipj * nelorbj_c * nelorbj_c, 1)))
    if(iessz) then
        ALLOCATE(jasmatsz_c(max(nelorbj_c * nelorbj_c, 1)))
    else
        ALLOCATE(jasmatsz_c(1))
    endif

else
    ALLOCATE(jasmatsz_c(1), jasmat_c(1))
endif
jasmatsz_c = 0.d0
jasmat_c = 0.d0

if(contractionj.ne.0) then
    jasmat_c = 0.d0
    do i = 1, nnozeroj_c
        call upsim(jasmat_c, ipj * nelorbj_c, nozeroj_c(i), scalej_c(i), .true., 1)
    enddo

    if(iessz) then
        jasmatsz_c = 0.d0
        do i = 1, nnozeroj_c
            call upsim(jasmatsz_c, nelorbj_c, nozeroj_c(i), scalejsz_c(i), .true., 1)
        enddo
    endif

    ALLOCATE(muj_c(nelorbjh, nelorbj_c), transpipj_sav(maxshellj))
    do i1 = 1, maxshellj
        allocate(transpipj_sav(i1)%col(npar3body_c))
    enddo
    transpipj_sav = transpipj
    muj_c = 0.d0

    allocate(addr_occ(occj_c))
    ! shift the muj matrix according to ioccj
    ix = 0
    do i = 1, occj
        if(ioccj(i).ne.0) then
            ix = ix + 1
            iy = 0
            !              multpointerj(ix)=0
            do j = 1, occj_c
                if(ioccj_c(j).ne.0) then
                    iy = iy + 1
                    if(muj_touch(i, j).ne.0) then
                        muj_c(ix, iy) = muj_c(ix, iy) + muj_tmp(i, j)
                    endif
                endif
            enddo
        endif
    enddo

    iy = 0
    addr_occ = 0
    do j = 1, occj_c
        if(ioccj_c(j).ne.0) then
            iy = iy + 1
            addr_occ(j) = iy
        endif
    enddo


    !    Simply put  transpipj  out of the loop

    do kk = 1, npar3body_c
        do jj = 1, multranspipj(kk)
            ii = (transpipj_sav(jj)%col(kk) - 1) / nelorbmaxj + 1
            ix = transpipj_sav(jj)%col(kk) - (ii - 1) * nelorbmaxj
            iy = addr_occ(ii)
            if(iy.ne.0) then
                transpipj(jj)%col(kk) = (iy - 1) * nelorbjh + ix
            else
                if(rank.eq.0) write(6, *) ' ERROR input '
            endif
        enddo
    enddo

    deallocate(addr_occ)

    if(moljyes) then
        indorb = 0
        indpar = 0
        ind = 0

        do i = 1, nshellj_c
            if(ioptorbj_c(i).ge.900000) then
                ind = ind + 1
                if(ioccj_c(ind).eq.1) then
                    indorb = indorb + 1
                    !            definition mu_c
                    do kk = 1, nparamj_c(i) / 2
                        ix = nint(vju_c(indpar + kk))

                        indn = indpar + kk + nparamj_c(i) / 2
                        muj_c(ix, indorb) = vju_c(indn)
                        !                multpointerj(ix)=multpointerj(ix)+1
                        !                mupointerj(multpointerj(ix),ix)=indorb
                        multranspipj(indn) = 1
                        transpipj(1)%col(indn) = nelorbjh * (indorb - 1) + ix
                    enddo

                endif
                indpar = indpar + nparamj_c(i)

            else
                do j = 1, multj_c(i)
                    ind = ind + 1
                    if(ioccj_c(ind).eq.1) indorb = indorb + 1
                enddo
                indpar = indpar + nparamj_c(i)
            endif
        enddo

    endif

    deallocate(muj_touch, muj_tmp, transpipj_sav, occshellj           &
            &, vjutouch, vjurold)

!   if(ipj.eq.2) then
!       call scontract_genj(nelorbjh, nelorbj_c, jasmat, jasmat_c, muj_c, psip)
!   else

!       call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh, nelorbj_c     &
!               &, nelorbj_c, jasmat, jasmat_c, muj_c, psip)

!   endif

    if(iessz) call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh, nelorbj_c&
            &, nelorbj_c, jasmatsz, jasmatsz_c, muj_c, psip)

    ! fill scalej, scalejsz, nozeroj, and nnozeroj
    indocc = (ipj * nelorbjh * (ipj * nelorbjh + 1)) / 2
    nnozeroj = indocc

    if(rank.eq.0) write(6, *) 'Number of non zero Jas lambda (root)', nnozeroj
    if(iessz) then 
    allocate(scalejsz(max(nnozeroj,1)))
    else
    allocate(scalejsz(1))
    endif
    ALLOCATE(scalej(1),nozeroj(1),nozerojder(1))
    scalej = 0.d0
    scalejsz = 0.d0
    nozeroj = 0
    nozerojder = 0

!   indocc = 0
!   do ix = 1, ipj * nelorbjh
!       do iy = ix, ipj * nelorbjh
!           ind = (iy - 1) * nelorbjh * ipj + ix
!           indocc = indocc + 1
!           scalej(indocc) = jasmat(ind)
!           if(iessz) scalejsz(indocc) = jasmatsz(ind)
!           nozerojder(indocc) = (iy - 1) * nelorbj * ipj + ix
!           nozeroj(indocc) = (iy - 1) * nelorbjh * ipj + ix
!       enddo
!   enddo

else ! contractionj>0

    scalej = scalej_c
    scalejsz = scalejsz_c
    allocate(muj_c(1, 1))
    muj_c = 0.d0
    if(.not.yes_sparse) then
    nozeroj = nozeroj_c
    nozerojder = nozeroj_c
    jasmat = 0.d0
    do i = 1, nnozeroj
        call upsim(jasmat, ipj * nelorbjh, nozeroj(i), scalej(i), .true., 1)
    enddo
    endif
    if(iessz) then
        jasmatsz = 0.d0
        do i = 1, nnozeroj
            call upsim(jasmatsz, nelorbjh, nozeroj(i), scalejsz(i), .true., 1)
        enddo
    endif

    deallocate(vjurold)

endif


!   symmetry for det zeta
ALLOCATE(jbraiesup(iesupr_c)                               &
        &, jbraiesup_sav(iesup_c + iesupind))
jbraiesup = 0
jbraiesup_sav = 0

if(rank.eq.0) write(*, *) ' Reading Z-AGP symmetries ....'
if(rankrep.eq.0) then
    read(ufort10, *, err = 156, end = 156)
endif
if(iesupind.ne.0) then
    !            write(6,*) ' This loop '
    ind = 0
    iesuptouched = 0
    do i = 1, iesupind
        if(rankrep.eq.0) then
            read(ufort10, *, err = 158, end = 158) ii, (ipsip(j), j = 1, abs(ii))
        endif

#ifdef PARALLEL
           call bcast_integer(ii,1,0,commrep_mpi)
#ifdef UNREL
           !   For unreliable  networks.
           call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
           call bcast_integer(ipsip,2*abs(ii),0,commrep_mpi)
#ifdef UNREL
           !   For unreliable  networks.
           call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
#endif

        ind = ind+1
        jbraiesup_sav(ind) = ii
        do j = 1, abs(ii)
            jbraiesup_sav(ind + j) = ipsip(j)
        enddo
        ind = ind + abs(ii)

        if(ii.gt.0) iesuptouched = iesuptouched + 1

        !            write(6,*) ' ipsip read =',i,' ',(ipsip(j),j=1,ii)
        do j = 1, ii
            if(ipsip(j).gt.0) then
                jbraiesup(ipsip(j)) = iesuptouched
            elseif(ipsip(j).lt.0) then
                jbraiesup(abs(ipsip(j))) = -iesuptouched
            endif
        enddo
    enddo
endif


! same zeta with same angular momentum and ion belong to the same
! root shell and must be included in the same symmetry
if(contraction.ne.0) then

    maxiesup = 0
    do i = 1, iesup_c
        if(iesuptransb(i).ne.0) maxiesup = i
    enddo

    do i = 1, maxiesup
        ii = iesuptransb(i)
        if(ii.ne.0) then
            do j = 1, maxiesup
                jj = iesuptransb(j)
                if(ii.eq.jj) then
                    if(jbraiesup(i).ne.jbraiesup(j)) then
                        write(6, *) 'det zeta', i, 'and', j, 'must be joined with symmetry!'
                        write(6, *) dup_c(i), dup_c(j)
                        call error('read_fort10', ' Error in reading symmetries! ', 1, rank)
                    endif
                endif
            enddo
        endif
    enddo
endif

if(rank.eq.0) write(6, *) ' Touched det zeta par =', iesuptouched
if(iesup.ne.0.and.ipc * iesuptouched.ne.iesup) then
    if(rank.eq.0) write(6, *) ' Warning iesup changed to ', ipc * iesuptouched
    iesup = ipc * iesuptouched
endif

! symmetry for jas zeta
ALLOCATE(jbraiesm(max(npar3bodyr_c, 1))                       &
        &, jbraiesm_sav(max(npar3body_c + iesmind, 1)))
jbraiesm = 0
jbraiesm_sav = 0

if(rankrep.eq.0) then
    read(ufort10, *, err = 160, end = 160)
endif
if(iesmind.ne.0) then
    if(rank.eq.0) write(6, *) ' Reading Z-jas symmetries .... '
    ind = 0
    vjutouched = 0
    do i = 1, iesmind
        if(rankrep.eq.0)  then
            read(ufort10, *, err = 162, end = 162) ii, (ipsip(j), j = 1, abs(ii))
        endif
#ifdef PARALLEL
           call bcast_integer(ii,1,0,commrep_mpi)
#ifdef UNREL
           !   For unreliable  networks.
           call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
           call bcast_integer(ipsip,abs(ii),0,commrep_mpi)
#ifdef UNREL
           !   For unreliable  networks.
           call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
#endif
           ind = ind+1
        jbraiesm_sav(ind) = ii
        do j = 1, abs(ii)
            jbraiesm_sav(ind + j) = ipsip(j)
        enddo
        ind = ind + abs(ii)
        if(ii.gt.0) vjutouched = vjutouched + 1

        do j = 1, ii
            if(ipsip(j).gt.0) then
                jbraiesm(ipsip(j)) = vjutouched
            elseif(ipsip(j).lt.0) then
                jbraiesm(abs(ipsip(j))) = -vjutouched
            endif
        enddo
    enddo
endif

! same zeta with same angular momentum and ion belong to the same
! root shell and must be included in the same symmetry
if(contractionj.ne.0) then
    maxnpar = 0
    do i = 1, npar3body_c
        if(iesuptransbj(i).ne.0) maxnpar = i
    enddo
    do i = 1, maxnpar
        ii = iesuptransbj(i)
        if(ii.ne.0) then
            do j = 1, maxnpar
                jj = iesuptransbj(j)
                if(ii.eq.jj) then
                    if(jbraiesm(i).ne.jbraiesm(j)) then
                        write(6, *) 'Jas zeta', i, 'and', j, 'must be joined with symmetry!'
                        write(6, *) vju_c(i), vju_c(j)
                        call error('read_fort10', ' Error in reading symmetries! ', 1, rank)
                    endif
                endif
            enddo
        endif
    enddo
endif

if(rank.eq.0) write(6, *) ' Touched Jas zeta par =', vjutouched
if(iesm.ne.0.and.vjutouched.ne.iesm) then
    if(rank.eq.0) write(6, *) ' Warning iesm changed to ', vjutouched
    iesm = vjutouched
endif

if(contraction.gt.0) then

    if(.not.yesdft.and.contraction.ne.0) then
        if(yes_complex) then
            call zgemm('N', 'N', nelorbh * ipf, nelcol_c, nelorb_c, (1.d0, 0.d0)&
                    &, mu_c, nelorbh * ipf, detmat_c, nelorb_c, (0.d0, 0.d0), projm, nelorbh * ipf)
        else
            call dgemm('N', 'N', ipf * nelorbh, nelcol_c, nelorb_c, 1.d0&
                    &, mu_c, ipf * nelorbh, detmat_c, nelorb_c, 0.d0, projm, ipf * nelorbh)
        endif
    endif

endif

if(yesmin_read) then
    iesup_read = ipc * iesup_atom
else
    iesup_read = ipc * iesup_c
endif

! MMM reading parameters for complex wave functions with different k-points!
! To be optimized

! read the optimized parameters
if(ireadmin.gt.0) then

    if(itestr.ne.-5) &
            call error('read_fort10', ' You should run optimization first ! ', 1, rank)

    if(rankrep.eq.0) then

        read(ufort10, *, err = 164, end = 164)

        if(rank.eq.0) write(6, *) ' Reading new WF parameters ....'

        if(ieskint.ne.0.and.ieskint.eq.ieskinr_pos) then

            read(ufort10, *, err = 164, end = 164) (ddwsz(i), i = 1, iesinv), (vju(i), i = 1, iesm)         &
                    &, (vj(i), i = 1, abs(iesd)), (ddw(i), i = 1, iesfree)                       &
                    &, (scale_c(i), i = 1, ipc * (nnozero_c + nnozero_eagp))                          &
                    &, (dup_c(i), i = 1, iesup_read), ((rion(ii, jj), ii = 1, 3), jj = 1, nion)

        elseif(ieskint.ne.0.and.ieskint.gt.ieskinr_pos) then

            read(ufort10, *, err = 164, end = 164) (ddwsz(i), i = 1, iesinv), (vju(i), i = 1, iesm)&
                    &, (vj(i), i = 1, abs(iesd)), (ddw(i), i = 1, iesfree)                       &
                    &, (scale_c(i), i = 1, ipc * (nnozero_c + nnozero_eagp))                                       &
                    &, (dup_c(i), i = 1, iesup_read), ((rion(ii, jj), ii = 1, 3), jj = 1, nion)&
                    &, rs, celldm(2:3)

            omega = celldm(2) * celldm(3)

            celldm(1) = (PI * nel * 4.d0 / 3.d0 / omega)**(1.d0 / 3.d0) * rs
            celldm(4:6) = 90.d0 * PI / 180.d0
            if(rank.eq.0) write(6, *) ' rs read =', rs
#ifndef PARALLEL

            givens2r = .false.
            call InitCell(nion, nel, yes_complex)
            kappa = kappar / lmin
            if(yes_tilted.and.ksq.ne.0.5d0) kappa = kappa/cond_find(metric(1,2),metric(1,3),metric(2,3))
            rs = (omega / nel * 3.d0 / 4.d0 / pi)**(1.d0 / 3.d0)
#endif

        else

            read(ufort10, *, err = 164, end = 164) (ddwsz(i), i = 1, iesinv)&
                    &, (vju(i), i = 1, iesm), (vj(i), i = 1, abs(iesd)), (ddw(i), i = 1, iesfree)&
                    &, (scale_c(i), i = 1, ipc * (nnozero_c + nnozero_eagp) * min(iessw, 1))&
                    &, (dup_c(i), i = 1, iesup_read * min(iesup, 1))

        endif
    endif ! if rank=0


#ifdef PARALLEL
        call bcast_real(ddwsz,iesinv,0,commrep_mpi)
        call bcast_real(vju,iesm,0,commrep_mpi)
        call bcast_real(vj,abs(iesd),0,commrep_mpi)
        call bcast_real(ddw,iesfree,0,commrep_mpi)
        call bcast_real(scale_c,size(scale_c),0,commrep_mpi)
        call bcast_real(dup_c,size(dup_c),0,commrep_mpi)
        if(ieskint.ne.0) then
        call bcast_real(rion,3*nion,0,commrep_mpi)
        endif
        if(ieskint.gt.ieskinr_pos) then
        call bcast_real(celldm,6,0,commrep_mpi)
              givens2r=.false.
              call InitCell(nion,nel,yes_complex)
              
              kappa = kappar/lmin
              if(yes_tilted.and.ksq.ne.0.5d0) kappa = kappa/cond_find(metric(1,2),metric(1,3),metric(2,3))
              rs=(omega/nel*3.d0/4.d0/pi)**(1.d0/3.d0)
        endif
#ifdef UNREL
        !   For unreliable  networks.
        call mpi_barrier(commrep_mpi,ierr)
!$omp barrier
#endif
#endif

    !         call checkiflagerr(iflagerr,rank,' ERROR in read_fort10 check your input fort.10')

    !vju_c
    call bconstrbr(iesm, npar3bodyr_c, jbraiesm, vju_c, vju)

    !detmat_c
    if(contraction.ne.0) then
        if(yes_complex) then
            call dscalzero(2 * nelorb_c * nelcol_c, 0.d0, detmat_c, 1)
            do i = 1, nnozero_c
                call upsim_complex(detmat_c, nelorb_c, nozero_c(i), scale_c(2 * i - 1), symmagp, ipf)
            enddo
        else
            call dscalzero(nelorb_c * nelcol_c, 0.d0, detmat_c, 1)
            do i = 1, nnozero_c
                call upsim(detmat_c, nelorb_c, nozero_c(i), scale_c(i), symmagp, ipf)
            enddo
        endif
        if(npar_eagp.ne.0) then
            eagp_pfaff = 0.d0
            do i = 1, nnozero_eagp
                ind = nnozero_c + i
                iy = (nozero_c(ind) - 1) / ndiff + 1
                ix = nozero_c(ind) - (iy - 1) * ndiff
                if(ipc.eq.1) then
                    eagp_pfaff(ix, iy) = scale_c(ind)
                    eagp_pfaff(iy, ix) = -scale_c(ind)
                else
                    eagp_pfaff(2 * ix - 1, iy) = scale_c(2 * ind - 1)
                    eagp_pfaff(2 * ix, iy) = scale_c(2 * ind)
                    eagp_pfaff(2 * iy - 1, ix) = -scale_c(2 * ind - 1)
                    eagp_pfaff(2 * iy, ix) = -scale_c(2 * ind)
                endif
            enddo
        endif
        !           if(allowed_averagek) then
        !!          From the real one to the effective
        !            call attach_phase2det(.false.,detmat_c)
        !           endif
        !          Here symmetrize agp for detecting the conventional phase
        !          if not molyes and manyfort10
        !          Here is the read used to average parameters
        !            if(manyfort10.and.(.not.molyes.or.write_effective)) then
        !            call  symmetrizeagp_ref(nnozero_c,nozero_c,jbradet,sjbradet&
        !           &,jbradetn,dsw,iessw0,psip,ipsip,detmat_c,nelorb_c,nelorb_at&
        !           &,nelcol_c,symmagp,yes_hermite)
        !            endif
        if(allowed_averagek) then
            !          From the effective to the real one
            call attach_phase2det(.true., detmat_c)
            if(rank.eq.0) write(6, *) ' Passi qui from eff. to real III '
        endif
    else
        if(yes_complex) then
            call dscalzero(2 * ipf * nelorbh * nelcolh, 0.d0, detmat, 1)
            do i = 1, nnozero
                call upsim_complex(detmat, ipf * nelorbh, nozero(i), scale(2 * i - 1), symmagp, ipf)
            enddo
        else
            call dscalzero(ipf * nelorbh * nelcolh, 0.d0, detmat, 1)
            do i = 1, nnozero_c
                call upsim(detmat, ipf * nelorbh, nozero(i), scale_c(i), symmagp, ipf)
            enddo
        endif
        if(npar_eagp.ne.0) then
            eagp_pfaff = 0.d0
            do i = 1, nnozero_eagp
                ind = nnozero_c + i
                iy = (nozero(ind) - 1) / ndiff + 1
                ix = nozero(ind) - (iy - 1) * ndiff
                if(ipc.eq.1) then
                    eagp_pfaff(ix, iy) = scale_c(ind)
                    eagp_pfaff(iy, ix) = -scale_c(ind)
                else
                    eagp_pfaff(2 * ix - 1, iy) = scale_c(2 * ind - 1)
                    eagp_pfaff(2 * ix, iy) = scale_c(2 * ind)
                    eagp_pfaff(2 * iy - 1, ix) = -scale_c(2 * ind - 1)
                    eagp_pfaff(2 * iy, ix) = -scale_c(2 * ind)
                endif
            enddo
        endif
        !           if(allowed_averagek) then
        !          From the real one to the effective
        !           call attach_phase2det(.false.,detmat)
        !           endif
        !          Here symmetrize agp for detecting the conventional phase
        !          if not molyes and manyfort10
        !            if(manyfort10.and.(.not.molyes.or.write_effective)) then
        !            call  symmetrizeagp_ref(nnozero_c,nozero_c,jbradet,sjbradet&
        !           &,jbradetn,dsw,iessw0,psip,ipsip,detmat,nelorb_c,nelorb_at&
        !           &,nelcol_c,symmagp,yes_hermite)
        !            endif
        if(allowed_averagek) then
            !          From the effective to the real one
            call attach_phase2det(.true., detmat)
            if(rank.eq.0) write(6, *) ' Passi qui from eff. to real IV '
        endif
    endif

    !jasmat_c
    if(iesfree.ne.0) then
        if(contractionj.ne.0) then
            call bconstrbra(iesfree, nnozeroj_c, jbraj, nozeroj_c, jasmat_c&
                    &, ipj * nelorbj_c, ddw)
            do i = 1, nnozeroj_c
                scalej_c(i) = jasmat_c(nozeroj_c(i))
            enddo
        else
            if(yes_sparse) then
            call bconstrbra_sparse(iesfree, nnozeroj, jbraj, nozerojder,jasmat, ipj * nelorbjh, ddw)
            do i = 1, nnozeroj
            scalej_c(i) = jasmat(i)
            enddo
            else
            call bconstrbra(iesfree, nnozeroj, jbraj, nozeroj, jasmat, ipj * nelorbjh, ddw)
            do i = 1, nnozeroj_c
                scalej_c(i) = jasmat(nozeroj(i))
            enddo
            endif
        endif
    endif

    !jasmatsz_c
    if(iesinv.ne.0) then
        if(contractionj.ne.0) then
            call bconstrbra(iesinv, nnozeroj_c, jbrajsz, nozeroj_c, jasmatsz_c, nelorbj_c, ddwsz)
            do i = 1, nnozeroj_c
                scalejsz_c(i) = jasmatsz_c(nozeroj_c(i))
            enddo
        else
            call bconstrbra(iesinv, nnozeroj, jbrajsz, nozeroj, jasmatsz, nelorbjh, ddwsz)
            do i = 1, nnozeroj_c
                scalejsz_c(i) = jasmatsz(nozeroj(i))
            enddo
        endif
    endif

    if(contraction.eq.0) then
        if(yes_complex) then
            do jj = 1, iesupr_c
                dupr(jj) = dup_c(2 * jj - 1)
            enddo
        else
            dupr = dup_c
        endif
        scale = scale_c
    else

        ! mu_c
        allocate(mu_touch(ipf * nelorb, nelorb_c))
        mu_touch = 0
        do jj = 1, iesup_c
            do kk = 1, multranspip(jj)
                iy = (transpip(kk)%col(jj) - 1) / (ipf * nelorbh) + 1
                ix = transpip(kk)%col(jj) - (iy - 1) * (ipf * nelorbh)
                if(mu_touch(ix, iy).eq.0) then
                    if(yes_complex) then
                        mu_c(2 * ix - 1, iy) = dup_c(2 * jj - 1)
                        mu_c(2 * ix, iy) = dup_c(2 * jj)
                    else
                        mu_c(ix, iy) = dup_c(jj)
                    endif
                    mu_touch(ix, iy) = 1
                else
                    if(yes_complex) then
                        mu_c(2 * ix - 1, iy) = mu_c(2 * ix - 1, iy) + dup_c(2 * jj - 1)
                        mu_c(2 * ix, iy) = mu_c(2 * ix, iy) + dup_c(2 * jj)
                    else
                        mu_c(ix, iy) = mu_c(ix, iy) + dup_c(jj)
                    endif
                endif
            enddo
        enddo

        ! dupr
        if(yes_complex) then
            do jj = 1, iesupr_c
                if(iesuptransb(jj).ne.0) then
                    dupr(iesuptransb(jj)) = dup_c(2 * jj - 1)
                endif
            enddo
        else
            do jj = 1, iesupr_c
                if(iesuptransb(jj).ne.0) then
                    dupr(iesuptransb(jj)) = dup_c(jj)
                endif
            enddo
        endif

        ! detmat
        if(yesdetmatc) then
            call scontract_mat_det(nelorbh, nelorbh, nelcolh&
                    &, nelorb_c, nelcol_c, detmat, detmat_c, mu_c, psip)
        endif

        if(.not.yesdft.and.contraction.ne.0) then
            if(yes_complex) then
                call zgemm('N', 'N', ipf * nelorbh, nelcol_c, nelorb_c, (1.d0, 0.d0), mu_c, ipf * nelorbh&
                        &, detmat_c, nelorb_c, (0.d0, 0.d0), projm, ipf * nelorbh)
            else
                call dgemm('N', 'N', ipf * nelorbh, nelcol_c, nelorb_c, 1.d0, mu_c, ipf * nelorbh&
                        &, detmat_c, nelorb_c, 0.d0, projm, ipf * nelorbh)
            endif
        endif

        if(yesdetmatc) then

            scale = 0.d0
            indocc = 0

            if(symmagp.or.ipf.eq.2) then
                do ix = 1, ipf * nelorbh
                    do iy = ix + ipf - 1, ipf * nelorbh + ndiff
                        ind = (iy - 1) * nelorbh * ipf + ix
                        indocc = indocc + 1
                        if(yes_complex) then
                            scale(2 * indocc - 1) = detmat(2 * ind - 1)
                            scale(2 * indocc) = detmat(2 * ind)
                        else
                            scale(indocc) = detmat(ind)
                        endif
                    enddo
                enddo
            else
                do ix = 1, ipf * nelorbh
                    do iy = 1, nelorbh + nelup - neldo
                        ind = (iy - 1) * nelorbh + ix
                        indocc = indocc + 1
                        if(yes_complex) then
                            scale(2 * indocc - 1) = detmat(2 * ind - 1)
                            scale(2 * indocc) = detmat(2 * ind)
                        else
                            scale(indocc) = detmat(ind)
                        endif
                    enddo
                enddo
            endif
        endif
        deallocate(mu_touch)
    endif ! contraction.eq.0

    if(contractionj.eq.0) then
        vjur = vju_c

        scalej = scalej_c
        scalejsz = scalejsz_c
    else

        ! muj_c
        allocate(muj_touch(nelorbj, nelorbj_c))
        muj_touch = 0
        do jj = 1, npar3body_c
            do kk = 1, multranspipj(jj)
                iy = (transpipj(kk)%col(jj) - 1) / nelorbjh + 1
                ix = transpipj(kk)%col(jj) - (iy - 1) * nelorbjh
                if(muj_touch(ix, iy).eq.0) then
                    muj_touch(ix, iy) = 1
                    muj_c(ix, iy) = vju_c(jj)
                else
                    muj_c(ix, iy) = muj_c(ix, iy) + vju_c(jj)
                endif
            enddo
        enddo
        ! vjur
        do jj = 1, npar3bodyr_c
            if(iesuptransbj(jj).ne.0) then
                vjur(iesuptransbj(jj)) = vju_c(jj)
            endif
        enddo

        ! jasmat, jasmatsz
!       if(ipj.eq.2) then
!           call scontract_genj(nelorbjh, nelorbj_c, jasmat, jasmat_c, muj_c, psip)
!       else
!           call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh, nelorbj_c &
!                   &, nelorbj_c, jasmat, jasmat_c, muj_c, psip)
!       endif
        if(iessz) call scontract_mat_jas(nelorbjh, nelorbjh, nelorbjh, nelorbj_c&
                &, nelorbj_c, jasmatsz, jasmatsz_c, muj_c, psip)
        scalej = 0.d0
        scalejsz = 0.d0
!       indocc = 0
!       do ix = 1, nelorbjh * ipj
!           do iy = ix, nelorbjh * ipj
!               ind = (iy - 1) * nelorbjh * ipj + ix
!               indocc = indocc + 1
!               scalej(indocc) = jasmat(ind)
!               if(iessz) scalejsz(indocc) = jasmatsz(ind)
!           enddo
!       enddo
        deallocate(muj_touch)
    endif ! contractionj.eq.0

endif ! endif ireadmin -> finished reading the new WF parameters. End of the fort.10 file!
!
!
! computing address tables for each orbitals
!
kp_ion = np - ieskin
allocate(rpar(kp_ion), adrlambda(2, kp_ion), jas_invariant(4, kp_ion), typeorb(occj_c), orbps(kp_ion))
rpar = 0.d0
adrlambda = 0
jas_invariant = 0.d0
typeorb = 0
orbps = .false.

call write_type_orb(nshellj_c, multj_c, ioccj_c, typeorb)

allocate(indpar_tab(nshell + 1), indorb_tab(nshell + 1), indshell_tab(nshell + 1))
allocate(slaterorb_read(nshell+nshellj))
iflagpip = 2
indshell = 0
indorb = 0
indpar = 0
indpar_tab = 0
indshell_tab = 0
indorb_tab = 0
DO ii = 1, nshell
    !  update only indorb indpar indshell
    indpar_tab(ii) = indpar
    indorb_tab(ii) = indorb
    indshell_tab(ii) = indshell
    if(LBox.eq.1.d0.and..not.yes_crystal) then

        CALL MAKEFUN_PBC(ioptorb(ii), ioccup, 0, 1, 1, 0, &
                & 1, indpar, indorb, indshell, nelorb, winv, dupr, &
                & r, rmu, psip, iflagpip, &
                & cnorm(ii), rmucos, rmusin,sinphase,cosphase)

    else

        CALL MAKEFUN(ioptorb(ii), 0, 1, 1, 0, &
                & 1, indpar, indorb, indshell, nelorb, winv, dupr, zetar&
                &,r,rmu,psip,iflagpip, cnorm(ii))

    endif
    slaterorb_read(ii)=slaterorb(ioptorb(ii))
END DO
nshell_det=nshell
indpar_tab(nshell + 1) = indpar
indorb_tab(nshell + 1) = indorb
indshell_tab(nshell + 1) = indshell
allocate(adr_nion(nion + 1), ind_nion(nshell))
ind_nion = 0
ind = 0
adr_nion(1) = 1
do ii = 1, nion
    do jj = 1, nshell
        if(kion(jj).eq.ii) then
            ind = ind + 1
            ind_nion(ind) = jj
        endif
    enddo
    adr_nion(ii + 1) = ind + 1
enddo

if(nshellj.ne.0) then

    allocate(indparj_tab(nshellj + 1), indorbj_tab(nshellj + 1), indshellj_tab(nshellj + 1))

    indshell = 0
    indorb = 0
    indpar = 0
    indparj_tab = 0
    indshellj_tab = 0
    indorbj_tab = 0

    DO ii = 1, nshellj
        !  update only indorb indpar indshell
        indparj_tab(ii) = indpar
        indorbj_tab(ii) = indorb
        indshellj_tab(ii) = indshell
        if(LBoxj.eq.1.d0) then

            CALL MAKEFUN_PBC(ioptorbj(ii), ioccj, 0, 1, 1, 0, &
                    & 1, indpar, indorb, indshell, nelorbj, winvj, vjur, &
                    & r, rmu, psip, iflagpip, cnorm(ii + nshell),&
                    &  rmucos, rmusin,sinphase,cosphase)

        else

            CALL MAKEFUN(ioptorbj(ii), 0, 1, 1, 0, &
                    & 1, indpar, indorb, indshell, nelorbj, winvj, vjur, zetar&
                    &, r, rmu, psip &
                    &, iflagpip, cnorm(ii + nshell))

        endif
    slaterorb_read(ii+nshell_det)=slaterorb(ioptorbj(ii))
    END DO

    indparj_tab(nshellj + 1) = indpar
    indorbj_tab(nshellj + 1) = indorb
    indshellj_tab(nshellj + 1) = indshell

    !  To comunicate that this table refers to the Jastrow
    indparj_tab(1) = -1

    allocate(adrj_nion(nion + 1), indj_nion(nshellj))
    ind = 0
    indj_nion = 0
    adrj_nion(1) = 1
    do ii = 1, nion
        do jj = 1, nshellj
            if(kionj(jj).eq.ii) then
                ind = ind + 1
                indj_nion(ind) = jj
            endif
        enddo
        adrj_nion(ii + 1) = ind + 1
    enddo

else ! no Jastrow-3body present

    allocate(indparj_tab(1), indorbj_tab(1)&
            &, indshellj_tab(1), adrj_nion(1), indj_nion(1))
    indparj_tab(1) = -1
    indorbj_tab(1) = 0
    indshellj_tab(1) = 0
    adrj_nion(1) = 0
    indj_nion(1) = 0

endif

timeg=cclock()

call update_kgrid

if(rank.eq.0) write(6,*) ' Time spent in  update_kgrid=',cclock()-timeg

if(iespbc.and.yesbump.and.rank.eq.0) then
    write(6, *) ' Warning minimum Z allowed for bump Gaussian functions &
            &    check your fort.10 ', 36.d0 / Lmin**2
endif

if(iespbc.and..not.yesdft) then
    !     Unit Ry/a.u.^3
    pressclass = temp * dble(nion) / cellscale(1) / cellscale(2) / cellscale(3)
    dcellclass(1:3) = temp * dble(nion) / cellscale(1:3)
    if(rank.eq.0) then
        write(6, *) ' Warning contribution perfect gas &
                &    pressure a.u. NOT included ', pressclass / 2.d0
    endif
else
    pressclass = 0.d0
    dcellclass = 0.d0
endif

niont = 0
niong = 0
do ii = 1, nion
    if(atom_number(ii).gt.0) then
        niont = niont + 1
    else
        niong = niong + 1
    endif
enddo

if(rank.eq.0.and.niong.ne.0) then
    write(6, *) ' Number of ghost atoms =', niong
    write(6, *) ' Number of real  atoms =', niont
endif

if(iesking.eq.3 * niong.and.niong.eq.niont) then
    if(rank.eq.0) write(6, *) ' Warning allocated warp matrix '
    allocate(warpmat(niont, niong))
    warpmat = 0.d0
else
    allocate(warpmat(niont, 1))
    warpmat = 0.d0
endif

!    Initializing once for all pointvj
numvjpar = num_vjpar(iesdrr)

if(rank.eq.0) then
    if(npsar.gt.0) then
     if(niesd.ne.numvjpar) then
        write(6, *) ' Warning number niesd suggested (# par one/two body Jastrow) = ', numvjpar
     endif
    else
     if(niesd.ne.nmax_ion + numvjpar) then
        write(6, *) ' Warning number niesd suggested (# par one/two body Jastrow) = ', nmax_ion + numvjpar
        write(6, *) ' in such case one has an independent one body for each different  atomic specie'
     endif
    endif
endif

pointvj = 0
do jj = 1, nion
    pointvj(1, jj) = pointvjf(type_atom(jj), niesd, nmax_ion, numvjpar)
    if(costz(jj).eq.1.d0) then
        pointvj(2, jj) = 8
    else
        pointvj(2, jj) = iesdr1iesd(iesdr)
    endif
enddo

iesd_onebody = iesdr1iesd(iesdr)
iesd_twobody = iesdr2iesd(iesdr)
if(case_map.ne.0.and..not.chosen_map) then
    if((iesd_twobody.eq.4.or.iesd_twobody.eq.-4.or.iesd_twobody.eq.8.or.&
            &iesd_twobody.eq.5).and.(iesd_onebody.eq.4.or.iesd_onebody.eq.-4.&
        &or.iesd_onebody.eq.8.or.iesd_onebody.eq.5)) then
    case_map = 5
    if(rank.eq.0) write(6, *) ' Warning one&two body exponential (case_map=5) '
    endif
endif

if(n_body_on.ne.0) then
    !  The value of the exp for r=0
    scale_one_body = 0.d0
    if(.not.iespbc) then
        do j = 1, 3
            rion_ref0(j) = sum(rion(j, :)) / nion
        enddo
    else
        rion_ref0(:) = 0.d0
    endif
    do jj = 1, nion
        if(iespbc) then
            rc(1:3) = rion_ref0(1:3) - rion(1:3, jj)
            call CartesianToCrystal(rc, 1)
            do k = 1, 3
                rc(k) = costz(jj) * map(rc(k), cellscale(k))
            enddo
            r0 = norm_metric(rc, metric)
        else
            rc(:) = (rion_ref0(:) - rion(:, jj)) * costz(jj)
            r0 = dsqrt(sum(rc(:)**2))
        endif
        scale_one_body = scale_one_body - &
                &jastrow_ei(r0, vj(pointvj(1, jj)), pointvj(2, jj)) * costz3(jj)
    enddo

else
    scale_one_body = 0.d0
endif
if(rank.eq.0) write(6, *) ' scale one body =', scale_one_body

if(rank.eq.0) then
    write(6, *)
    write(6, *) ' END reading the wave function fort.10 '
    write(6, *)
endif
if(molecular.gt.0) then
    molyes = .true.
else
    molyes = .false.
endif
if(allocated(nozeron)) deallocate(nozeron)

! $$$$$$$$$$$$$$$$$$$$$$ END READING fort.10 $$$$$$$$$$$$$$$$$$$$$$$$$$$

return

! ERRORS section
100 call error(' read_fort10 ', 'fort.10 wrong or absent.', 1, rank)
102 call error(' read_fort10 ', 'Error reading celldm.', 1, rank)
104 call error(' read_fort10 ', 'Error reading the begin.', 1, rank)
106 call error(' read_fort10 ', 'Error reading the begin, next section.', 1, rank)
108 call error(' read_fort10 ', 'Error reading zeta and/or ion positions.', 1, rank)
110 call error(' read_fort10 ', 'Error reading forces constraints.', 1, rank)
112 call error(' read_fort10 ', 'Error reading forces constraints, next section.', 1, rank)
114 call error(' read_fort10 ', 'Error reading 2-body Jastrow.', 1, rank)
116 call error(' read_fort10 ', 'Error reading 2-body Jastrow, next section.', 1, rank)
118 call error(' read_fort10 ', 'Error reading 2-body Jastrow, next section.', 1, rank)
120 call error(' read_fort10 ', 'Error reading 2-body Jastrow, next section.', 1, rank)
122 call error(' read_fort10 ', 'Error reading AGP orbitals.', 1, rank)
124 call error(' read_fort10 ', 'Error reading AGP orbitals, next section.', 1, rank)
126 call error(' read_fort10 ', 'Error reading AGP orbitals, next next section.', 1, rank)
128 call error(' read_fort10 ', 'Error reading Jastrow orbitals.', 1, rank)
130 call error(' read_fort10 ', 'Error reading Jastrow orbitals, next section.', 1, rank)
132 call error(' read_fort10 ', 'Error reading Jastrow orbitals, next next section.', 1, rank)
134 call error(' read_fort10 ', 'Error reading determinant occupations.', 1, rank)
136 call error(' read_fort10 ', 'Error reading determinant occupations, next section.', 1, rank)
138 call error(' read_fort10 ', 'Error reading Jastrow occupations.', 1, rank)
140 call error(' read_fort10 ', 'Error reading Jastrow occupations, next section.', 1, rank)
142 call error(' read_fort10 ', 'Error reading non-zero AGP matrix elements.', 1, rank)
144 call error(' read_fort10 ', 'Error reading non-zero AGP matrix elements symmetries.', 1, rank)
146 call error(' read_fort10 ', 'Error reading non-zero AGP matrix elements symmetries, next section.', 1, rank)
148 call error(' read_fort10 ', 'Error reading non-zero Jastrow matrix elements.', 1, rank)
150 call error(' read_fort10 ', 'Error reading non-zero spin Jastrow matrix elements.', 1, rank)
152 call error(' read_fort10 ', 'Error reading non-zero Jastrow matrix elements symmetries.', 1, rank)
154 call error(' read_fort10 ', 'Error reading non-zero Jastrow matrix elements symmetries, next section.', 1, rank)
156 call error(' read_fort10 ', 'Error reading Z-AGP symmetries.', 1, rank)
158 call error(' read_fort10 ', 'Error reading Z-AGP symmetries, next section.', 1, rank)
160 call error(' read_fort10 ', 'Error reading Z-Jastrow symmetries.', 1, rank)
162 call error(' read_fort10 ', 'Error reading Z-Jastrow symmetries, next section.', 1, rank)
164 call error(' read_fort10 ', 'Error reading new wavefunction parameters.', 1, rank)
end subroutine read_fort10

        !--------------------------------------------------------------------------
        !  This subroutine has input iesinv,iesm,iesd,iesfree,iessw,iesup,ieskin
        !  and reads the input file fort.10 defining the wavefunction in TurboRVB.
        !  Each parameter of the input is modified if it is nonzero and remains
        !  with the same sign of input (case iesinv,iesfree,iessw).
        !  The absolute value of each output value corresponds to
        !  the number of parameters  that can be minimized
                !  consistently with fort.10.
                !--------------------------------------------------------------------------

subroutine read_fort10_fast
                use allio
                implicit none
                !     integer iesinv,iesm,iesd,iesfree,iessw,iesup,ieskin
                integer icq, icp, icd, icj, kionpar, ilaenv&
                &, dimpsip, dimipsip, indorb, numpaired,ind_nelorb
                integer i,j,k,ii,jj,ind,indpar,nshelljmax,ioptorbcontr,nelupr
                !     integer i,j,indpar,occj,occ,nelorbj,icq,icp,nelup,nel,nion,nshell
                !    1,nshellj,iesdrr,iesupr,npar3body,nnozero,nnozeroj,iesupind,iesmind
                !    1,iesfreer,iesswr,ieskinr,indorb,ind,icd,icj,kion,nparam,ioptorb
                !    1,niesd

                !     integer, dimension(:), allocatable:: ipsip,mult,multj,nparamj
                !    1,ioptorbj,ioccupj
                !     logical, dimension(:), allocatable:: orbcost
                character(2) chars
                character(3) checkpbc
        character(5) checkpbc_c
        !     real*8, dimension(:), allocatable:: psip
        logical iesonz, iestot, check1,check2, checks1, checks2, iespsip, iesipsip,yes_contraction
        logical, allocatable :: occupied(:)

                molyes = .false.
                ! it is assumed that the file fort.10 is already open
                rewind(10)
                yes_crystal = .false.
                yes_complex = .false.
                yes_tilted= .false.
                read(10, *, err = 100, end =100) chars, checkpbc_c
        if(checkpbc_c.eq.'PBC_T') then
        yes_tilted = .true.
        iespbc = .true.
                yes_crystal = .true.
                gamma_point = .true.
                elseif(checkpbc_c.eq.'PBC_C') then
                iespbc = .true.
                yes_crystal = .true.
                gamma_point = .true.
                else
                rewind(10)
                read(10, *, err = 100, end =100) chars, checkpbc
        if(checkpbc.eq.'PBC') then
        iespbc = .true.
        else
        iespbc =.false.
        gamma_point = .false.
        endif
        endif
                rewind(10)
                if(iespbc) then
                !
                celldm(4:6) = 90.d0*PI/180.d0
        read(10, *, err = 100, end =100)
                if(.not.yes_crystal) then
                read(10, *, err = 100, end =100) rs, celldm(2:3), phase(:)
        else
        if(yes_tilted) then
        read(10, *, err = 100, end =100) (s2r(:, i), i = 1, 3), phase(:), phase_down(:)
                else
                read(10, *, err = 100, end =100) rs, celldm(2:3), phase(:), phase_down(:)
                endif
                endif
                !
                opposite_phase = .true.
                same_phase = .true.

                do i = 1, 3
                if(abs(phase(i)-nint(phase(i))).ne.0.5d0.and.abs(phase(i)-nint(phase(i))).ne.0.d0) yes_complex = .true.
        ! treat gamma and other real boundaries, like (0.5,0.5,0.5), without complex w.f.
        ! unless complex algorithm is forced on input (nshell<0)
        if(yes_crystal) then
        if(abs(phase_down(i)-nint(phase_down(i))).ne.0.5d0.and.abs(phase_down(i)-nint(phase_down(i))).ne.0.d0) yes_complex = .true.
        if((abs(phase(i)-nint(phase(i))).ne.0.5.and.&
        &(phase_down(i)-nint(phase_down(i)).ne.nint(phase(i))-phase(i))).or.&
        &(abs(phase(i)-nint(phase(i))).eq.0.5.and.&
                &abs(phase_down(i)-nint(phase_down(i))).ne.0.5)) opposite_phase = .false.
        if((abs(phase(i)-nint(phase(i))).ne.0.5.and.&
        &(phase_down(i)-nint(phase_down(i)).ne.phase(i)-nint(phase(i)))).or.&
        &(abs(phase(i)-nint(phase(i))).eq.0.5.and.&
                &abs(phase_down(i)-nint(phase_down(i))).ne.0.5)) same_phase = .false.
        endif
        enddo
        if(opposite_phase.and.same_phase) same_phase = .false.
        else
        opposite_phase = .true.
                same_phase = .false.
                endif


                write(*, *) ' Reading the begin in fast. . . . '
                read(10, *, err = 101, end =101)
                read(10, *, err = 101, end =101) nelup, nel, nion

                if(nion.lt.0) then
                yesbump = .true.
                nion = -nion
        else
        yesbump = .false.
        endif

        if(nel.lt.0) then
        nel = -nel
        ipf = 2
                else
                ipf = 1
                endif

                if(nelup.le.0) then
        pfaffup = .true.
        nelup = -nelup
                else
                pfaffup = .false.
                endif



                if(mod(nelup, nel).ne.0) then
        numpaired = 2*(nelup/nel)
        nelupr = mod(nelup, nel)
                else
                numpaired = 2*(nelup/nel-1)
                nelupr = nel
                endif

                ndiff = nelupr-neldo
        npar_eagp = 0
        if(ipf.eq.2) then
        if(mod(nelupr+neldo, 2).ne.0) then
        ndiff = 1+numpaired
        else
        ndiff = numpaired
        endif
        npar_eagp =(ndiff*(ndiff-1))/2
        endif

        if(rs.lt.0.d0.and..not.yes_tilted) then
        yeslbox = .true.
        rs = (-3.d0/4.d0/Pi*rs**3*celldm(2)*celldm(3)/nel)**(1.d0/3.d0)
                else
                yeslbox = .false.
                endif


                read(10, *, err = 101, end =101)
        read(10, *, err = 101, end =101) nshell, nshellj

        ! we can use a complex w.f. for open systems too
        if(nshell.lt.0) then
        yes_complex = .true.
        nshell = -nshell
                endif
                yes_crystalj = .false.
                if(nshellj.lt.0.or.yes_tilted) then
                nshellj = abs(nshellj)
                yes_crystalj = .true.
                endif


                if(yes_complex) then
        ipc = 2
        else
        ipc =1
        endif

        read(10, *, err = 101, end =101)
                if(nelup.gt.nel) then
                read(10, *, err = 106, end =106) iesdrr, iesupr, npar3body, iesswr_eagp, nnozero_eagp
                else
                iesswr_eagp = 0
                nnozero_eagp = 0
                read(10, *, err = 106, end =106) iesdrr, iesupr, npar3body
        endif
        case_map=0
        if(abs(iesdrr).ge.100) then
        case_map=abs(iesdrr)/100
         if(iesdrr.gt.0) then
         iesdrr=iesdrr-100*case_map
         else
         iesdrr=iesdrr+100*case_map
         endif
        endif


        if(iesdrr.eq.-12.or.iesdrr.eq.-22.or.iesdrr.eq.-26) then
        ipj = 2
        else
        ipj =1
        endif

        ! Check the various cases for the spin Jastrow
        iesonz = .false.
        if(iesdrr.eq.-8.or.iesdrr.eq.-9.or.iesdrr.eq.-18.or.iesdrr.eq.-28.or.iesdrr.eq.-29&
        &.or.iesdrr.eq.-19.or.iesdrr.eq.-10.or.iesdrr.eq.-11.or.iesdrr.eq.-16) then
        iessz = .true.
        if(iesdrr.eq.-18.or.iesdrr.eq.-8.or.iesdrr.eq.-10.or.iesdrr.eq.-28) then
        iesonz = .true.
        elseif(onebodysz) then
        iesonz = .true.
        endif
        else
        iessz = .false.
        endif
        read(10, *, err = 101, end =101)
                read(10, *, err = 101, end =101)  nnozero, nnozeroj
        forcesymm = .false.
        if(nnozero.lt.0) then
        nnozero = -nnozero
        forcesymm = .true.
                opposite_phase = .false.
                same_phase = .true.
                endif
                read(10, *, err = 101, end =101)
        read(10, *, err = 101, end =101)  iesupind, iesmind
        read(10, *, err = 101, end =101)
                read(10, *, err = 101, end =101)  iesfreer, iesswr, ieskinr, ireadminr


        if(allocated(zetar_fast))&
        &deallocate(zetar_fast, rion_fast, atom_number_fast)
        allocate(zetar_fast(nion), rion_fast(3, nion), atom_number_fast(nion))

                write(*, *) ' Reading zeta and rion  in fast . . . . '


                read(10, *, err = 102, end =102)
                zmax = 0.d0
                npsar = 0
                core_pseudo=0.d0

                do j = 1, nion
                read(10, *, err = 102, end =102) zetar_fast(j), atom_number_fast(j),(rion_fast(i, j), i = 1, 3)
                if(zetar_fast(j).ne.atom_number_fast(j).and.atom_number_fast(j).gt.0) &
                & npsar = npsar+1
                  if(atom_number_fast(j).gt.zmax) then
                  zmax = atom_number_fast(j)
                  core_pseudo=zmax-zetar_fast(j)
                  endif
                enddo


        if(allocated(ipsip)) then
        dimipsip = size(ipsip)
                iesipsip = .true.
                deallocate(ipsip)
                else
                iesipsip = .false.
                endif

                allocate(ipsip(6*nion))


                write(*, *) ' Reading ieskin in fast. . . . '
        read(10, *, err = 103, end =103)
                check1 = .false.
                ieskinold = ieskin
                if(ieskin.eq.0.or.ieskin.gt.1) check1 = .true.
                ieskin = 0
                iesking =0
                do i = 1, ieskinr
                read(10, *, err = 103, end =103) icq,(ipsip(j), j = 1, 2*abs(icq))
                icj = abs(ipsip(1))
        if(icj.gt.0) then
        icp = atom_number_fast(icj)
                else
                icp = 0
                endif
                if(icq.gt.0.and.icp.gt.0) ieskin = ieskin+1
        if(icq.gt.0.and.icp.le.0) iesking = iesking+1
        enddo
        if(check1) then
        if(ieskinold.eq.0) then
        ieskin = 0
        iesking = 0
                elseif(ieskinold.gt.ieskin) then
                ieskin = ieskinold
                write(6, *) ' Warning chosen input ieskin> ieskinr ', ieskin
                !        elseif(ieskinold.le.3) then
                !        ieskin=ieskin+ieskinold
                !        write(6,*) ' Warning also cell forces computed ',ieskin
                elseif(ieskinold.le.ieskin) then
                ieskin = ieskinold
                iesking = 0
        endif
        endif
        ieskint = iesking+ieskin
                ieskinr_pos = ieskint

                deallocate(ipsip)

                if(allocated(psip)) then
        dimpsip = size(psip)
                deallocate(psip)
                iespsip = .true.
                else
                iespsip = .false.
        endif
        allocate(psip(10000))
                check1 = .false.
                if(iesd.eq.0) check1 = .true.
        write(*, *) ' Reading 2-body jastrow in fast . . . . '
        read(10, *, err = 104, end =104)
                symmagp = .true.
                if(abs(iesdrr).eq.0) then
        read(10, *, err = 104, end =104) niesd
        if(niesd.lt.0) symmagp = .false.
        iesd = 0
        else
                read(10, *, err = 104, end =104)  niesd,(psip(i), i = 1, abs(niesd))
                if(niesd.lt.0) symmagp = .false.
        iesd = abs(niesd)
                endif
                if(check1) iesd = 0

                deallocate(psip)

                write(*, *) ' Reading det shells in fast. . . . '
                read(10, *, err = 105, end =105)
        allocate(mult(nshell), nparam(nshell), ioptorb(nshell))
                if(yes_complex) then
        allocate(psip(2*iesupr))
                else
                allocate(psip(iesupr))
                endif
                allocate(occupied(iesupr))
        occupied = .false.
        yes_contraction = .false.
        iesup_atom= 0
        do i = 1, nshell
                read(10, *, err = 105, end =105)  mult(i), nparam(i), ioptorb(i)
        if(yes_complex.and.nparam(i).gt.1) then
        read(10, *, err = 105, end =105)  kionpar,(psip(j), j = 1, nparam(i)*3/2)
                else
                read(10, *, err = 105, end =105)  kionpar,(psip(j), j = 1, nparam(i))
        endif
        if(ioptorb(i).ge.900000) then
        do j = 1, nparam(i)/2
        occupied(nint(psip(j)))= .true.
        enddo
        endif
        if(ioptorb(i).eq.1000000) molyes = .true.
        if(ioptorb(i).ne.ioptorbcontr(ioptorb(i), LBox, 1)) yes_contraction = .true.
        if(ioptorb(i).ne.1000000) iesup_atom = iesup_atom+nparam(i)
                enddo
                deallocate(psip)

                nshelljmax = max(nshellj, 1)

                allocate(multj(nshelljmax), ioptorbj(nshelljmax), nparamj(nshelljmax))
        allocate(psip(max(npar3body, 1)))
                write(*, *) ' Reading jas shells in fast. . . . '

                read(10, *, err = 106, end =106)
                do i = 1, nshellj
                read(10, *, err = 106, end =106) multj(i), nparamj(i), ioptorbj(i)
                read(10, *, err = 106, end =106) kionpar,(psip(j), j = 1, nparamj(i))
                enddo
                deallocate(psip)

                occ = 0
                do i = 1,nshell
                occ = occ+mult(i)
                enddo

                write(*, *) ' Reading det occupation in fast. . . . ', occ
                read(10, *, err = 107, end =107)
                allocate(ioccup(occ))
                do i = 1, occ
                read(10, *, err = 107, end =107) ioccup(i)
                enddo

                !    Calculation of nelorb
                ind = 0
                !    Hybrid contribution to the basis. Hybrid are assumed to be occupied.
                nelorb= 0
                do i = 1, iesupr
                if(occupied(i)) nelorb = nelorb+1
                enddo
                nelorb =nelorb/ipf
             
                ind_nelorb=0
                do i = 1, nshell
                 do k=1,max(nparam(i)/2, 1)
                   do j= 1, mult(i)
                !        Normal contracted orbital or uncontracte (nparam=1)
                    if(ioccup(ind+j).ne.0.and.ioptorb(i).lt.900000) then
                    ind_nelorb=ind_nelorb+1
                     if(ind_nelorb.le.iesupr) then
                      if(.not.occupied(ind_nelorb)) nelorb=nelorb+1
                     else
                      nelorb=nelorb+1
                     endif
                    endif
                   enddo
                 enddo
                 ind = ind+mult(i)
                enddo
                deallocate(occupied)
      
                occj = 0
                do i = 1, nshellj
                occj =occj+multj(i)
                enddo

        allocate(ioccj(max(occj, 1)))
                write(*, *) ' Reading jas occupation in fast. . . . ', occj
                read(10, *, err = 108, end =108)
                nelorbj = 0
                do i = 1,occj
                read(10, *, err = 108, end =108) ioccj(i)
        if(ioccj(i).ne.0) nelorbj = nelorbj+1
        enddo


        allocate(orbcost(nelorbj))

                indorb = 0
                ind = 0
                orbcost= .false.

        do i = 1, nshellj
                do j = 1, multj(i)
                ind = ind+1
                if(ioccj(ind).ne.0) then
        indorb = indorb+1
        if(ioptorbj(i).eq.200) orbcost(indorb)= .true.
        endif
        enddo
        enddo


        write(*, *) ' Reading det nnozero in fast. . . . '
        read(10, *, err = 109, end =109)
                do i = 1, nnozero
                read(10, *, err = 109, end =109)
        enddo
        do i = 1, nnozero_eagp
                read(10, *, err = 109, end =109)
                enddo
                check1 = .false.
                checks1 = .false.
                if(iessw.eq.0) check1 = .true.
        if(iessw.lt.0) checks1 = .true.
        iessw = 0

        allocate(ipsip(2*(nnozero+nnozero_eagp)))
        write(*, *) ' Reading det nnozero symmetries in fast....'
        read(10, *, err = 110, end =110)
                do i = 1, iesswr
                read(10, *, err = 110, end =110) icd,(ipsip(j), j = 1, 2*abs(icd))
                if(icd.gt.0) then
        if(symmagp.and.ipc.eq.2.and.yes_correct) then
        iessw = iessw+4
        else
                iessw = iessw+ipc
                endif
                endif
                enddo
                do i = 1, iesswr_eagp
                read(10, *, err = 110, end =110) icd,(ipsip(j), j = 1, 2*abs(icd))
                if(icd.gt.0) then
        if(symmagp.and.ipc.eq.2.and.yes_correct) then
        iessw = iessw+4
        else
                iessw = iessw+ipc
                endif
                endif
                enddo
                deallocate(ipsip)
                if(check1) iessw = 0
                if(checks1) iessw = -iessw
        write(*, *) ' Reading jas nnozero in fast....'
        read(10, *, err = 111, end =111)
                do i = 1, nnozeroj
                !   ixj_1(i),iyj_1(i),jasmat_1(i)
                read(10, *, err = 111, end =111)
                enddo
                if(iessz) then
                write(*, *) ' Reading jas-sz nnozero in fast....'
        read(10, *, err = 112, end =112)
                do i = 1, nnozeroj
                ! ixj_1(i),iyj_1(i),jasmatsz_1(i)
                read(10, *, err = 112, end =112)
                enddo
                endif

                allocate(ipsip(2*nnozeroj))

                check1 = .false.
                check2 = .false.
                checks1= .false.
                checks2 = .false.
                if(iesfree.eq.0) check1 = .true.
                if(iesfree.lt.0) checks1 = .true.
        if(iesinv.eq.0) check2 = .true.
        if(iesinv.lt.0) checks2 = .true.
        iesfree = 0
        iesinv =0
        write(*, *) ' Reading jas nnozero symmetries in fast....'
        read(10, *, err = 113, end =113)
                do i = 1, iesfreer
                read(10, *, err = 113, end =113) icj,(ipsip(j), j = 1, 2*abs(icj))
                iestot = .false.
                !do j=1,2*abs(icj)
                !  if(orbcost(abs(ipsip(j)))) iestot=.true.
                !enddo
                do j= 1, 2*abs(icj)
                if (abs(ipsip(j))<=nelorbj) then
                if(orbcost(abs(ipsip(j)))) iestot = .true.
                else
                if(orbcost(abs(ipsip(j))-nelorbj)) iestot = .true.
        endif
        enddo

        if(icj.gt.0) then
        if((twobodyoff.and.iestot).or..not.twobodyoff) iesfree = iesfree+1
        if(onebodysz.or.twobodyoff) then
        if(.not.iesonz.or.(iesonz.and.iestot)) iesinv = iesinv+1
        else
        if(.not.iesonz.or.(iesonz.and..not.iestot)) iesinv = iesinv+1
        endif
        endif
        enddo


        if(check1) iesfree = 0
        if(checks1) iesfree = -iesfree
        if(check2) iesinv = 0
        if(checks2) iesinv = -iesinv
        if(.not.iessz) iesinv = 0 ! In any event if the Jastrow Sz is off iesinv=0

        write(*, *) ' Reading Z-det  symmetries in fast....'
        read(10, *, err = 114, end =114)
                check1 = .false.
                if(iesup.eq.0) check1 = .true.
                iesup = 0
                deallocate(ipsip)
                allocate(ipsip(iesupr))
                do i = 1, iesupind
                read(10, *, err = 114, end =114) icp,(ipsip(j), j = 1, abs(icp))
                        if(icp.gt.0) iesup = iesup+ipc
                enddo
                deallocate(ipsip)
                        if(check1) iesup = 0


        write(*, *) ' Reading Z-jas  symmetries in fast....'
        read(10, *, err = 115, end =115)
                check1 = .false.
                if(iesm.eq.0) check1 = .true.
                iesm = 0

                allocate(ipsip(max(npar3body, 1)))
                do i = 1, iesmind
                read(10, *, err = 115, end =115) icq,(ipsip(j), j = 1, abs(icq))
                if(icq.gt.0) iesm = iesm+1
                enddo

                if(check1) iesm = 0

                deallocate(mult, ioptorb, nparam, ioccup, multj,ioptorbj, nparamj, ioccj, orbcost, ipsip)
                if(iespsip) then
        allocate(psip(dimpsip))
                psip = 0.d0
                endif
                if(iesipsip) then
                allocate(ipsip(dimipsip))
                ipsip = 0
                endif
                rewind(10)
                return

                !********** ERRORS!! ***********

                100  call error(' read_fort10_fast ', 'Error reading PBC.', 1, rank)
                101  call error(' read_fort10_fast ', 'Error reading the begin.', 1, rank)
                        102  call error(' read_fort10_fast ', 'Error reading zeta/ionic positions.', 1, rank)
                        103  call error(' read_fort10_fast ', 'Error reading forces.', 1, rank)
                        104  call error(' read_fort10_fast ', 'Error reading Jastrow 2-body.', 1, rank)
                        105  call error(' read_fort10_fast ', 'Error reading determinant shells.', 1, rank)
                106  call error(' read_fort10_fast ', 'Error reading jastrow shells.', 1, rank)
                        107  call error(' read_fort10_fast ', 'Error reading determinant occupation.', 1, rank)
                        108  call error(' read_fort10_fast ', 'Error reading Jastrow occupation.', 1, rank)
                        109  call error(' read_fort10_fast ', 'Error reading AGP matrix elements.', 1, rank)
                        110  call error(' read_fort10_fast ', 'Error reading AGP matrix elements symmetries.', 1, rank)
                111  call error(' read_fort10_fast ', 'Error reading Jastrow matrix elements.', 1, rank)
                        112  call error(' read_fort10_fast ', 'Error reading spin Jastrow matrix elements.', 1, rank)
                        113  call error(' read_fort10_fast ', 'Error reading Jastrow matrix elements symmetries.', 1, rank)
                        114  call error(' read_fort10_fast ', 'Error reading determinant orbital symmetries.', 1, rank)
                        115  call error(' read_fort10_fast ', 'Error reading Jastrow orbital symmetries.', 1, rank)
                116  call error(' read_fort10_fast ', 'Error reading new wave function parameters.', 1, rank)

                        end subroutine read_fort10_fast

                        subroutine write_fort10(unit)
                        use allio
                        implicit none
                        integer unit
                        integer, dimension(:), allocatable :: indexv
                real*8, dimension(:), allocatable :: sortvect
                integer i, ii,j, iref, ind, indpar, numpaired
        logical, dimension(:), allocatable :: yespip
        rewind(unit)
                if(iespbc) then
        if(.not.yes_crystal.and..not.yes_complex) then ! for complex case use always PBC_C basis set
        write(unit, *) '# PBC rs, Ly/Lx, Lz/Lx phase '
                if(yeslbox) then
                write(unit, '(3f18.12,3f12.8)') -cellscale(1),(celldm(i), i = 2, 3),(phase(i), i = 1, 3)
                else
                write(unit, '(3f18.12,3f12.8)') rs,(celldm(i), i = 2, 3),(phase(i), i = 1, 3)
                endif
                else
                if(yeslbox) then
                write(unit, *) '# PBC_C rs, Ly/Lx, Lz/Lx phase_up  phase_down'
                write(unit, '(3f18.12,6f12.8)') -cellscale(1),(celldm(i), i = 2, 3),(phase(i), i = 1, 3),(phase_down(i), i = 1, 3)
                else
                if(yes_tilted) then
                write(unit, *) '# PBC_T a b c  phase_up  phase_down'
                write(unit, '(9f18.12,6f12.8)') (s2r(:, i), i = 1, 3),(phase(i), i = 1, 3),(phase_down(i), i = 1, 3)
                else
                write(unit, *) '# PBC_C rs, Ly/Lx, Lz/Lx phase_up  phase_down'
        write(unit, '(3f18.12,6f12.8)') rs,(celldm(i), i = 2, 3),(phase(i), i = 1, 3),(phase_down(i), i = 1, 3)
                endif
                endif
                endif
                endif
                write(unit, *) '# Nelup  #Nel  # Ion '
        if(ipf.eq.2) then
        numpaired = ndiff-mod(nel, 2)
                nelup = nelup+nel*numpaired/2
                nel = -nel
                endif
                if(pfaffup) nelup = -nelup
                if(yesbump) then
                write(unit, *) nelup, nel, -nion
        else
        write(unit, *) nelup, nel, nion
        endif
        if(pfaffup) nelup = -nelup
        if(ipf.eq.2) then
        nel = -nel
        nelup = nelup-nel*numpaired/2
        endif

        write(unit, *) '# Shell Det.   # Shell Jas. '
        if(yes_crystalj) nshellj_c = -nshellj_c
        if(forcecomplex) then
        write(unit, *) -nshell_c, nshellj_c
        else
        write(unit, *) nshell_c, nshellj_c
        endif
        if(yes_crystalj) nshellj_c = -nshellj_c

        if(chosen_map) then
         if(iesdrr.ge.0) then
         iesdrr=iesdrr+100*case_map
         else
         iesdrr=iesdrr-100*case_map
         endif
        endif
        if(npar_eagp.ne.0) then
        write(unit, *) '# Jas 2body  # Det   #  3 body atomic par.  iessw_eagp #eagp'
        write(unit, *) iesdrr, iesup_c, npar3body_c, iesswr_eagp, nnozero_eagp
                else
                write(unit, *) '# Jas 2body  # Det   #  3 body atomic par.  '
                write(unit, *) iesdrr, iesup_c, npar3body_c
                endif
                write(unit, *) '# Det mat. =/0  # Jas mat. =/0  '
        if(forcesymm) then
        write(unit, *)  -nnozero_c, nnozeroj_c
        else
        write(unit, *)  nnozero_c, nnozeroj_c
        endif
        write(unit, *) ' # Eq. Det atomic par.  # Eq. 3 body atomic. par.'
        write(unit, *)  iesupind, iesmind
        write(unit, *) '# unconstrained iesfree,iessw,ieskinr,I/O flag '
        write(unit, *)  iesfreer, iesswr, ieskinr, 0
        write(unit, *) '# Ion coordinates '

        do j = 1, nion
        write(unit, *) zetar(j), atom_number(j),(rion(i, j), i = 1, 3)
                enddo
                write(unit, *) '#  Constraints for forces: ion - coordinate'
                do i = 1, ieskinr
                ind= abs(ion_table(i)%mult)
                write(unit, *) ion_table(i)%mult,(ion_table(i)%ion(j), ion_table(i)%comp(j), j = 1, ind)
                enddo


                write(unit, *) '#          Parameters Jastrow two body'
                if(iesdr.eq.0) then
                if(symmagp) then
        write(unit, *) 0
        else
        write(unit, *) -1
        endif
        elseif(abs(iesdr).le.3) then
        write(unit, *) (vj(i), i = 1, abs(iesdr))
                else
                if(symmagp) then
        write(unit, *) abs(niesd),(vj(i), i = 1, abs(niesd))
                else
                write(unit, *) -abs(niesd),(vj(i), i = 1, abs(niesd))
                endif
                endif
                write(unit, *) '#          Parameters atomic wf'
                indpar = 0
                do i = 1, nshell_c
                write(unit, *) mult_c(i), nparam_c(i), ioptorb_c(i)
                if(ioptorb_c(i).eq.900000.or.ioptorb_c(i).eq.1000000) then  ! hybrids or molecular orbitals
                if(yes_complex) then
                write(unit, *) kion_c(i),(nint(dup_c(2*(indpar+j)-1)), j = 1, nparam_c(i)/2) &   ! writing only odd terms until
        ! nparam_c(i)/2 since they are real quantities
        &,(dup_c(2*(indpar+j)-1), dup_c(2*(indpar+j)), j = nparam_c(i)/2+1, nparam_c(i))  ! complex quantities
        else
        write(unit, *) kion_c(i),(nint(dup_c(indpar+j)), j = 1, nparam_c(i)/2) &
        &,(dup_c(indpar+j), j = nparam_c(i)/2+1, nparam_c(i))
                endif
                else  ! atomic orbitals
                if(yes_complex.and.nparam_c(i).gt.1) then
                write(unit, *) kion_c(i),(dup_c(2*(indpar+j)-1), j = 1, nparam_c(i)/2) &
        &,(dup_c(2*(indpar+j)-1), dup_c(2*(indpar+j)), j = nparam_c(i)/2+1, nparam_c(i))
                elseif(yes_complex.and.nparam_c(i).eq.1) then
                write(unit, *) kion_c(i), dup_c(2*indpar+1)
        else
        write(unit, *) kion_c(i),(dup_c(indpar+j), j = 1, nparam_c(i))
                endif
                endif
                indpar = indpar+nparam_c(i)
                enddo

                write(unit, *) '#  Parameters atomic Jastrow wf '
        indpar = 0
        do i = 1, nshellj_c
        write(unit, *) multj_c(i), nparamj_c(i), ioptorbj_c(i)
                if(ioptorbj_c(i).eq.900000.or.ioptorbj_c(i).eq.1000000) then
        write(unit, *) kionj_c(i),(nint(vju_c(indpar+j)), j = 1, nparamj_c(i)/2)&
        &,(vju_c(indpar+j), j = nparamj_c(i)/2+1, nparamj_c(i))
                else
                write(unit, *) kionj_c(i),(vju_c(indpar+j), j = 1, nparamj_c(i))
                endif
                indpar = indpar+nparamj_c(i)
                enddo

                write(unit, *) '#  Occupation atomic orbitals  '
                do i = 1, occ_c
                write(unit, *) ioccup_c(i)
                enddo

                write(unit, *) '#  Occupation atomic orbitals  Jastrow  '
                do i = 1, occj_c
                write(unit, *) ioccj_c(i)
                enddo

                write(unit, *) ' #          Nonzero values of  detmat '
        call update_kiontot  ! In case nelorb_c has changed in convertmol_c

        if(contraction.ne.0) then
        !       From the real one to the effective
        if(allowed_averagek) call attach_phase2det(.false., detmat_c)
                if(rank.eq.0) write(6, *) ' Passi qui from real to eff V '
        do i = 1, nnozero_c
        iy=(nozero_c(i)-1)/nelorb_c+1
        ix = nozero_c(i)-(iy-1)*nelorb_c
        if(yes_complex) then
        if(no_sjbra) then
        if(sjbradet(i)) then
        if(abs(detmat_c(2*nozero_c(i))).gt.eps8.or&
        &.abs(detmat_c(2*nozero_c(i)-1)).gt.eps8) write(6, *) ' Warning sjbradet &
        & mat  removed ', detmat_c(2*nozero_c(i)-1), detmat_c(2*nozero_c(i))
        write(unit, *) ix, iy, 0.d0, 0.d0
        elseif(real_agp) then
        if(abs(detmat_c(2*nozero_c(i))).gt.eps8) write(6, *) &
        &' Warning dirty imaginary part for real_agp removed ', detmat_c(2*nozero_c(i))
        write(unit, *) ix, iy, detmat_c(2*nozero_c(i)-1), 0.d0
        else
        write(unit, *) ix, iy, detmat_c(2*nozero_c(i)-1), detmat_c(2*nozero_c(i))
        endif
        else
        if(real_agp) then
        if(abs(detmat_c(2*nozero_c(i))).gt.eps8) write(6, *) &
        &' Warning dirty imaginary part for real_agp removed ', detmat_c(2*nozero_c(i))
                write(unit, *) ix, iy, detmat_c(2*nozero_c(i)-1), 0.d0
        else
        write(unit, *) ix, iy, detmat_c(2*nozero_c(i)-1), detmat_c(2*nozero_c(i))
                endif
                endif
                else
                write(unit, *) ix, iy, detmat_c(nozero_c(i))
                endif
                enddo
                !       From the effective to the real one
                if(allowed_averagek)&
                &call attach_phase2det(.true., detmat_c)
                if(rank.eq.0) write(6, *) ' Passi qui from eff to real  VI '
        else
        !       From the the real one to the effective
        if(allowed_averagek) call attach_phase2det(.false., detmat)
                if(rank.eq.0) write(6, *) ' Passi qui from real to eff VII '
                do i = 1, nnozero
                iy=(nozero(i)-1)/(ipf*nelorbh)+1
        ix = nozero(i)-(iy-1)*nelorbh*ipf
        if(yes_complex) then
        if(no_sjbra) then
        if(sjbradet(i)) then
        if(abs(detmat(2*nozero(i))).gt.eps8.or&
        &.abs(detmat(2*nozero(i)-1)).gt.eps8) write(6, *) &
        &' Warning sjbradet mat  removed ', detmat(2*nozero(i)-1), detmat(2*nozero(i))
        write(unit, *) ix, iy, 0.d0, 0.d0
        elseif(real_agp) then
        if(abs(detmat(2*nozero(i))).gt.eps8) write(6, *) &
        &' Warning dirty imaginary part for real_agp removed ', detmat(2*nozero(i))
        write(unit, *) ix, iy, detmat(2*nozero(i)-1), 0.d0
        else
        write(unit, *) ix, iy, detmat(2*nozero(i)-1), detmat(2*nozero(i))
        endif
        else
        if(real_agp) then
        if(abs(detmat(2*nozero(i))).gt.eps8) write(6, *) &
        &' Warning dirty imaginary part for real_agp removed ', detmat(2*nozero(i))
                write(unit, *) ix, iy, detmat(2*nozero(i)-1), 0.d0
        else
        write(unit, *) ix, iy, detmat(2*nozero(i)-1), detmat(2*nozero(i))
                endif
                endif
                else
                write(unit, *) ix, iy, detmat(nozero(i))
                endif
                enddo
                !       From the effective to real (to avoid that if called two times
                !       does not work)
                if(allowed_averagek) call attach_phase2det(.true., detmat)
                if(rank.eq.0) write(6, *) ' Passi qui from eff to real VIII '
        endif
        do i = 1, nnozero_eagp
        iy =(nozero_c(i+nnozero_c)-1)/ndiff+1
        ix = nozero_c(i+nnozero_c)-(iy-1)*ndiff
        if(ipc.eq.1) then
        write(unit, *) ix, iy, eagp_pfaff(ix, iy)
                else
                write(unit, *) ix, iy, eagp_pfaff(2*ix-1, iy), eagp_pfaff(2*ix, iy)
                endif
                enddo

                write(unit, *) '#   Grouped par.  in the chosen ordered basis'

        if(nnozero_c.gt.0) then
        allocate(sortvect(nnozero_c+1), indexv(nnozero_c))
                allocate(yespip(nnozero_c))
                sortvect = 0.d0
                sortvect(1:nnozero_c)= abs(jbradet(1:nnozero_c))
        call dsortx(sortvect, 1, nnozero_c, indexv)
                endif


                indnn = 0
                icount = 1
                do i = 1, iesswr

                do while(sortvect(icount).lt.i.and.icount.le.nnozero_c)
                icount = icount+1
                enddo

                ind= 0
        do while(sortvect(icount).eq.i.and.icount.le.nnozero_c)
                ind = ind+1
                j = indexv(icount)
        ipsip(ind)= j
        if(jbradet(j).lt.0) ipsip(ind)= -ipsip(ind)
                yespip(ind)= sjbradet(j)
                icount = icount+1
                enddo

                do j = 1, ind
                ii = abs(ipsip(j))
                ipsip(ind+2*j)=(nozero_c(ii)-1)/nelorb_c+1
                ipsip(ind+2*j-1)= nozero_c(ii)-(ipsip(ind+2*j)-1)*nelorb_c
        if(ipsip(j).lt.0) ipsip(ind+2*j-1)= -ipsip(ind+2*j-1)
        if(yespip(j)) ipsip(ind+2*j)= -ipsip(ind+2*j)
                enddo

                if(ind.ne.0) then
        !rimosso ipf che moltiplicava nelorb_at
        if(ireadminr.lt.0.and.abs(ipsip(ind+1)).le.nelorb_at) then
        write(unit, *) -ind,(ipsip(ind+2*j-1), ipsip(ind+2*j), j = 1, ind)
                else
                write(unit, *) ind,(ipsip(ind+2*j-1), ipsip(ind+2*j), j = 1, ind)
                endif
                else

                indnn = indnn+1
                ii = jbradetn(indnn)
        write(unit, *) ii,(jbradetn(indnn+j), j = 1, -2*ii)
                indnn = indnn-2*ii

                endif

                enddo
                if(nnozero_c.gt.0) deallocate(sortvect, indexv, yespip)
                if(npar_eagp.ne.0) then

        ind = nnozero_c
        allocate(sortvect(nnozero_eagp+1), indexv(nnozero_eagp))
                sortvect = 0.d0
                sortvect(1:nnozero_eagp)= abs(jbradet(nnozero_c+1:nnozero_c+nnozero_eagp))
                call dsortx(sortvect, 1, nnozero_eagp, indexv)
                icount = 1
                do while(sortvect(icount).eq.0.d0)
        icount = icount+1
        enddo
        iref= sortvect(icount)-1
        icount = 1


        do i =1, iesswr_eagp
        do while(sortvect(icount).lt.i+iref.and.icount.le.nnozero_eagp)
                icount = icount+1
                enddo

                ind= 0
        do while(sortvect(icount).eq.i+iref.and.icount.le.nnozero_eagp)
                ind = ind+1
                j = indexv(icount)
                ipsip(ind)= j
                if(jbradet(nnozero_c+j).lt.0) ipsip(ind)= -ipsip(ind)
                icount = icount+1
                enddo

                do j = 1, ind
        ii = abs(ipsip(j))
                ipsip(ind+2*j)=(nozero_c(nnozero_c+ii)-1)/ndiff+1
                ipsip(ind+2*j-1)= nozero_c(nnozero_c+ii)-(ipsip(ind+2*j)-1)*ndiff
        if(ipsip(j).lt.0) ipsip(ind+2*j-1)= -ipsip(ind+2*j-1)
                enddo
                if(ind.ne.0) then
        !rimosso ipf che moltiplicava nelorb_at
        if(ireadminr.lt.0) then
        write(unit, *) -ind,(ipsip(ind+2*j-1), ipsip(ind+2*j), j = 1, ind)
                else
                write(unit, *) ind,(ipsip(ind+2*j-1), ipsip(ind+2*j), j = 1, ind)
        endif
        else
        indnn = indnn+1
                ii = jbradetn(indnn)
                write(unit, *) ii,(jbradetn(indnn+j), j = 1, -2*ii)
                indnn = indnn-2*ii
                endif
                enddo
                deallocate(sortvect, indexv)
        endif


        write(unit, *) ' #          Nonzero values of  jasmat '

        if(contractionj.ne.0) then
         do i = 1, nnozeroj_c
         iy=(nozeroj_c(i)-1)/(ipj*nelorbj_c)+1
         ix = nozeroj_c(i)-(iy-1)*nelorbj_c*ipj
         write(unit, *) ix, iy, jasmat_c(nozeroj_c(i))
         enddo
        else
          if(yes_sparse) then
          psip(1:nnozeroj_c)=0.d0
          do i=1,nnozeroj
          psip(nozerojder(i))=jasmat(i)
          enddo
          endif
         do i = 1,nnozeroj_c
          if(yes_sparse) then
          iy = nozeroj_c(i+nnozeroj_c)
          ix = nozeroj_c(i)
          write(unit, *) ix, iy, psip(i)
          else
          iy =(nozeroj(i)-1)/(ipj*nelorbjh)+1
          ix = nozeroj(i)-(iy-1)*nelorbjh*ipj
          write(unit, *) ix, iy, jasmat(nozeroj(i))
          endif
         enddo
        endif


                if(iessz) then
        write(unit, *) ' #       Nonzero values of  jasmat Sz '
        !          if(contractionj.eq.0) then
        !          do i=1,nelorbj_c
        !           call dcopy(nelorbj_c,jasmatsz(nelorbj*(i-1)+1),1
        !    1,jasmatsz_c(nelorbj_c*(i-1)+1),1)
        !          enddo
        !          endif

        if(contractionj.ne.0) then
        do i = 1, nnozeroj_c
        iy=(nozeroj_c(i)-1)/nelorbj_c+1
        ix = nozeroj_c(i)-(iy-1)*nelorbj_c
        write(unit, *) ix, iy, jasmatsz_c(nozeroj_c(i))
                enddo
                else
                do i = 1,nnozeroj
                iy =(nozeroj(i)-1)/nelorbjh+1
                ix = nozeroj(i)-(iy-1)*nelorbjh
        write(unit, *) ix, iy, jasmatsz(nozeroj(i))
                enddo
                endif


                endif

                write(unit, *) '# Eq. par. in the 3-body Jastrow in the chosen basis '

                if(nnozeroj_c.gt.0) then
                allocate(sortvect(nnozeroj_c+1), indexv(nnozeroj_c))
        sortvect = 0.d0
        sortvect(1:nnozeroj_c)= abs(jbraj(1:nnozeroj_c))
                call dsortx(sortvect, 1, nnozeroj_c, indexv)
                endif
                indnn = 0
                icount = 1
                do i = 1, iesfreer


                !                do j=1,nnozeroj_c
                !                if(abs(jbraj(j)).eq.i) then
                !                ind=ind+1
                !                  ipsip(ind)=j
                !                  if(jbraj(j).lt.0) ipsip(ind)=-ipsip(ind)
                !                endif
                !                enddo

                do while(sortvect(icount).lt.i.and.icount.le.nnozeroj_c)
        icount = icount+1
        enddo

        ind= 0
        do while(sortvect(icount).eq.i.and.icount.le.nnozeroj_c)
                ind = ind+1
                j = indexv(icount)
                ipsip(ind)= j
                if(jbraj(j).lt.0) ipsip(ind)= -ipsip(ind)
                icount = icount+1
                enddo


                do j = 1, ind
                if(yes_sparse) then
                ipsip(ind+2*j-1)= nozeroj_c(abs(ipsip(j)))
                ipsip(ind+2*j)= nozeroj_c(abs(ipsip(j))+nnozeroj_c)
                else
                ipsip(ind+2*j)=(nozeroj_c(abs(ipsip(j)))-1)/(ipj*nelorbj_c)+1
                ipsip(ind+2*j-1)= nozeroj_c(abs(ipsip(j)))                     &
        &-(ipsip(ind+2*j)-1)*nelorbj_c*ipj
                endif
                if(ipsip(j).lt.0) ipsip(ind+2*j-1)= -ipsip(ind+2*j-1)
                enddo
        if(ind.ne.0) then
        write(unit, *) ind,(ipsip(ind+2*j-1), ipsip(ind+2*j), j = 1, ind)
                else
                indnn = indnn+1
                ii= jbrajn(indnn)
                write(unit, *) ii,(jbrajn(indnn+j), j = 1, 2*abs(ii))
                indnn = indnn+2*abs(ii)
                endif
                enddo
                if(nnozeroj_c.gt.0)  deallocate(sortvect, indexv)

                write(unit, *) '# Eq. par. in the atomic Det par.in the chosen basis '
        if(iesupind.ne.0) then
        ind = 0
        do i = 1, iesupind
        ind = ind+1
        ii= jbraiesup_sav(ind)
                write(unit, *) ii,(jbraiesup_sav(ind+j), j = 1, abs(ii))
                ind = ind+abs(ii)
                        enddo
                        endif
                        write(unit, *) '# Eq. par. in the atomic 3-body  par. in the chosen basis '
                if(iesmind.ne.0) then
                ind = 0
                do i = 1, iesmind
                ind = ind+1
                ii= jbraiesm_sav(ind)
                write(unit, *) ii,(jbraiesm_sav(ind+j), j = 1, abs(ii))
                        ind = ind+abs(ii)
                        enddo
                        endif
                        end subroutine write_fort10

                        subroutine update_kiontot
                use allio, only:kiontot, kiontotj, nelcol_c, nelorb_c, nelorbj_c, nshell_c, mult_c&
                &, multj_c, nshellj_c, kion_c, kionj_c, ioccup_c, ioccj_c, nelorb_at, ioptorb_c
                use constants, only:ipj, ipf
                implicit none
                integer ind, indj, i, j
                if(allocated(kiontot)) deallocate(kiontot)
                        if(allocated(kiontotj)) deallocate(kiontotj)
                allocate(kiontot(nelcol_c), kiontotj(ipj*nelorbj_c))
                        kiontot = 0
                        kiontotj = 0
                        ind= 0
                indj = 0
                do i =1, nshell_c
                do j = 1, mult_c(i)
                        ind = ind+1
                        if(ioccup_c(ind).eq.1) then
                indj = indj+1
                kiontot(indj)= kion_c(i)
                        if(ioptorb_c(i).eq.1000000) kiontot(indj)= 0
                        endif
                        enddo
                        enddo
                        if(indj.le.nelorb_c/2.and.ipf.eq.2) then
                        kiontot(nelorb_at/2+1:nelorb_at)= kiontot(1:nelorb_at/2)
                        endif

                        ind = 0
                        indj = 0
                        do i = 1, nshellj_c
                        do j= 1, multj_c(i)
                        ind = ind+1
                        if(ioccj_c(ind).eq.1) then
                        indj = indj+1
                        kiontotj(indj)= kionj_c(i)
                        endif
                        enddo
                        enddo
                        if(ipj.eq.2) kiontotj(indj+1:2*indj)= kiontotj(1:indj)
                        return
                        end subroutine update_kiontot
                        function pointvjf(j, niesd, nmax_ion, numvjpar)
        implicit none
        integer pointvjf, j, niesd, nmax_ion, numvjpar
        if(niesd.le.numvjpar+1) then
        !    Here the two body can be defined only by one parameter
        pointvjf = niesd
        else
        !  With this change the two body can be defined by an arbitrary number of parameters
        !    It is assumed only that the one body has only one parameter per atomic specie and the last one (nmax_ion) is at the position niesd
                pointvjf = min(j+numvjpar, niesd)
                !    Old version was assuming only one parameter for the two-body Jastrow
                !    pointvjf=j+1
                endif
                return
                end
                function num_vjpar(iesdrr)
                implicit none
                integer num_vjpar, iesdrr

                if(iesdrr.eq.0.or.iesdrr.eq.-10.or.iesdrr.eq.-11) then
        num_vjpar = 0
        elseif(iesdrr.eq.5.or.iesdrr.eq.3) then
        num_vjpar = 3
        elseif(iesdrr.eq.-27.or.iesdrr.eq.2.or.iesdrr.eq.9.or.iesdrr.eq.10&
        &.or.iesdrr.eq.6.or.iesdrr.eq.-2.or.iesdrr.eq.30.or.iesdrr.eq.31) then
        num_vjpar = 2
        else
        num_vjpar =1
endif

return
end


