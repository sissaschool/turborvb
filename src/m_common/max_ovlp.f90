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

!
! @file max_ovlp.f90
! @brief Module for calculating geminal matrix lambda with maximum overlap
!
! This module allows to calculate the geminal matrix lambda that has the larger overlap with a previous     !
! one in a different basis. The calculation is performed numerically to allow the imposition of constraints !
! on the imaginary part of the matrix and selecting the elements of the matrix to be taken in to account.   !
!                                                                                                          !
! Ref:                                                                                                      !
! - parbcs.pdf                                                                                              !
! - Numerical Recipes, Press et al.                                                                         !
!
! @author TurboRVB group
! @date 2022
! @version 1.0
!

subroutine max_ovlp (size1, size2, type_lambda, nnozero_c, nozero_c, jbradet&
        &, img, max_iter, optimize, L1, L1mod, SL2, SR2, SL12, SR12, L2&
        &, prec, Z, rank, nprocu, mpi_comm_world)
    use constants, only: ipf
    implicit none

    !> @brief Dimension of matrices in the previous basis (size1) and new basis (size2)
    integer, intent(in) :: size1, size2
    !> @brief Type of lambda matrix symmetry: 1=hermitian, 2=symmetric (also imaginary part), 3=non-symmetric
    integer, intent(in) :: type_lambda
    !> @brief Size of the vector nozero_c and jbradet
    integer, intent(in) :: nnozero_c
    !> @brief Vector containing information to obtain the whole matrix lambda from compact notation
    integer, dimension(nnozero_c), intent(in) :: nozero_c
    !> @brief Vector pointing to independent variable number associated with symmetry
    integer, dimension(nnozero_c), intent(in) :: jbradet
    !> @brief Input to understand if wave function is real (1) or complex (2)
    integer, intent(in) :: img
    !> @brief Maximum number of iterations for optimization
    integer, intent(in) :: max_iter
    !> @brief Information on which matrix elements to optimize (1=yes, 0=no)
    logical, dimension(size2*img, size2), intent(in) :: optimize
    !> @brief Lambda matrix in the previous basis (input)
    real(8), dimension(size1*img, size1), intent(in) :: L1
    !> @brief Module of the previous wave function (input)
    real(8), intent(in) :: L1mod
    !> @brief Overlap matrices of the new basis (spin-up=left, spin-down=right)
    real(8), dimension(size2*img, size2), intent(in) :: SL2, SR2
    !> @brief Overlap matrices calculated between the two different bases
    real(8), dimension(size1*img, size2), intent(in) :: SL12, SR12
    !> @brief Lambda matrix in the new basis (input/output)
    real(8), dimension(size2*img, size2), intent(inout) :: L2
    !> @brief Precision required for the overlap calculation
    real(8), intent(in) :: prec
    !> @brief Overlap between L1 and L2 (output)
    real(8), intent(out) :: Z
    !> @brief MPI rank, number of processes, and communicator
    integer, intent(in) :: rank, nprocu, mpi_comm_world

    !All the quantities referred to the previous basis have the label 1, the quantities of the new one with 2
    !(NB: in the notation of the parbcs S*1=S*, S*2=S*', L1=Lambda, L2=Lambda', S*12=\bar{S}*, Op=O' ,
    !Ob=\bar{O}, A=A)

    !Size is the dimension of the matrices, img are the input to understand if the wf is real==1 or complex==2,
    !rank, nprocu and mpicommworld are variables of mpi, nnozero_c is the size of the vector nozero_c and jbradet,
    !nsym is the total number of variables once applied the symmetries, ix and iy are auxiliary variables to indicate
    !a position of the matrix lambda, type_lambda contains the information on the symmetry of the lambda: 1=hermitian,
    !2=symmetric (also the imaginary part) and 3=non symmetric (nnozero_c=size2**2)
    integer :: i, j, size1, size2, img, max_iter, rank, nprocu, mpi_comm_world, nnozero_c, nsym, ix, iy, type_lambda
    !optimize contains the information on which element of the matrix have to be optimize (1==yes) (0==no)
    logical, dimension(size2*img, size2) :: optimize
    !nozero_c is the vector containing the information to obtain the whole matrix lambda starting from the compact
    !notation and containing only the upper diagonal of lambda, in a way that nozero_c(i)=ix+(iy-1)*nelorbc,
    !jbradet is the vector pointing at  the indipendent variable number associated to the symmetry (e.g. if the
    !elements i and j are simmetric jbradet(i)=jbradet(j). This is only to my further memory, so that I can understand
    !again in the future
    integer, dimension(nnozero_c) :: jbradet, nozero_c
    !Z is the overlap between L1 and L2, delta is the number multiplying the derivatives when applied to L2,
    !prec is the precision required for the overlap, A is a variable for the fast calculation, L*mod is the
    !module of the * wf, modder is the module of the derivative, delta0 is the initial delta for the optimization,
    !lambda is the factor that multiplies the direction, lambda is  lambda= (g h)/(h A h) and is the factor of
    ! the movement along the h direction, prevZ is the Z calculated size2**2 iteration ago
    real(8) :: Z, delta, prec, A, L1mod, L2mod, modder, delta0, ddL2, lambda, prevZ
    !SL and SR are the overlap matrices of the bases the parbs spinup==left and spindown==right, SR12 and SL12
    !are the overlap matrices calculated between the two different bases, L1 and L2 are the derivatives, Os and
    !Op are matrices necessary for the fast calculation, A is a variable for the fast calculation, dL2 is the
    !matrix of the components of the derivatives, L2eff is the L2 with the phase attached (L2 will be almost
    !always real)
    !prevH and prevG are the previous conjugate direction and gradient (NOW
    ! THE GRADIENT HAS THE DIMENSION OF THE INDIPENDENT VARIABLES)(Numerical
    !Recipes, Press, et al.)
    real(8), dimension(size2*img, size2) :: SL2, SR2, L2
    real(8), allocatable :: dL2(:, :), Op(:, :), Ob(:, :), L2eff(:, :), prevG(:), prevH(:)
    real(8), dimension(size1*img, size1) :: L1
    real(8), dimension(size1*img, size2) :: SL12, SR12
    !count_sym is the vector that tells the number of recurrency of every symmetry
    integer, allocatable :: count_sym(:)

    integer ierr
    !#ifdef PARALLEL
    !  include 'mpif.h'
    !#endif

    !CCC starting from here

    !Understanding which is the number of indipendent variables
    nsym = 0
    do i = 1, nnozero_c
        if (abs(jbradet(i)) .gt. nsym) nsym = jbradet(i)
    end do

    !  if (rank.eq.0) write(6,*) jbradet (1), jbradet(2), jbradet(3)
    if (rank .eq. 0) write (6, *) "Symmetry lambda =", type_lambda
    if (rank .eq. 0) write (6, *) "nsym =", nsym
    if (rank .eq. 0) write (6, *) "nnozero_c =", nnozero_c

    !Allocating the matrices and the vectors
    allocate (dL2(img*size2, size2), Op(img*size2, size2), Ob(img*size2, size2), L2eff(img*size2, size2))
    allocate (prevG(img*nsym), prevH(img*nsym), count_sym(nsym))

    !Initializing to the identity
    L2 = 0.d0
    do i = 1, size2
        if (img .eq. 1) then
            L2(i, i) = 1.d0
        else
            L2(2*i - 1, i) = 1.d0
        end if
    end do

    call initialize_symm(size2, img, nnozero_c, nsym, count_sym, jbradet, optimize, L2, nozero_c, type_lambda)

    !Attaching the phase to the initial L2
    L2eff = L2
    call attach_phase2det(.true., L2eff)

    prevG = 0.d0
    prevH = 0.d0

    !It calculates the module of the WF1 (L1mod = A = Tr[O1 L1^*] ) and of the WF2 (L2mod=Tr [ Op L2^*])
    !NB: O1 is always the same, while Op must be calculated every time that L2 changes
    !It calculates the Ob, must be calculated only at the beginning

    call calcOb(size1, size2, img, SL12, L1, SR12, Ob, rank, nprocu, mpi_comm_world)
    call calcOp(size2, img, SR2, L2eff, SL2, Op, rank, nprocu, mpi_comm_world)
    call calcA(size2, img, L2eff, L2mod, Op, rank, nprocu, mpi_comm_world)
    call calcZ(size1, size2, img, L2eff, Ob, L1mod, L2mod, Z, rank, nprocu, mpi_comm_world)

    if (rank .eq. 0) then
        write (*, *) "Starting overlap:"
        write (*, *) "L1=", L1mod, ",   L2=", L2mod
        write (*, *) "Z=", Z
        write (*, *) ""
        write (*, *) ""
    end if
    prevZ = 0.d0

    do i = 1, max_iter

        !It Calculates the new values with the new L2 for the overlap
        !if(rank.eq.0) write(6,*) ' after 1',i
        call calcOp(size2, img, SR2, L2eff, SL2, Op, rank, nprocu, mpi_comm_world)
        !if(rank.eq.0) write(6,*) ' after 2',i
        call calcA(size2, img, L2eff, L2mod, Op, rank, nprocu, mpi_comm_world)
        !if(rank.eq.0) write(6,*) ' after 3',i
        call calcZ(size1, size2, img, L2eff, Ob, L1mod, L2mod, Z, rank, nprocu, mpi_comm_world)
        !if(rank.eq.0) write(6,*) ' after 4',i
        call calcDeriv(size2, img, Ob, Op, L1mod, L2mod, L2eff, dL2, optimize, modder, rank, nprocu, mpi_comm_world)

        !Removing the phase to the derivatives
        call attach_phase2det(.false., dL2)
        !     if(rank.eq.0) write(6,*) ' Iteration =',i,Z,prevZ

        if (mod(i, 20) == 0) then
            if (rank .eq. 0) then
                write (*, *) "Iteration number", i
                write (*, *) "Overlap value Z=", Z
                write (*, *) ""
                write (*, *) ""
            end if
            !        call initialize_symm(size2, img, nnozero_c, nsym,  count_sym, jbradet, optimize, L2, nozero_c, type_lambda)
            if ((abs(Z - prevZ) < prec*0.1) .or. (Z < prevZ)) exit
            prevZ = Z
        end if

        !Minimization along the conjugate direction with the derivative bisection (it's done using the derivatives without the
        !phase
        !  if(rank.eq.0) write(6,*) ' after 5',i
        call minDB(size1, size2, type_lambda, nnozero_c, nozero_c, jbradet, nsym, count_sym, img, L1mod, L2mod, L2, Ob, SL2, &
                   SR2, A, dL2, optimize, prevH, prevG, prec, lambda, rank, nprocu, mpi_comm_world)

        !Calculating the new L2eff
        L2eff = L2
        call attach_phase2det(.true., L2eff)

    end do

    if (rank .eq. 0) then
        write (*, *) "Final iteration  number", i
        write (*, *) "Overlap value:"
        write (*, *) "Z=", Z
        write (*, *) "|dL2|=", modder
        !     write(*,*) "Lambda=", lambda
        write (*, *) ""
        write (*, *) ""
    end if
    L2 = L2eff
    deallocate (dL2, Op, Ob, prevG, prevH, count_sym, L2eff)

end subroutine max_ovlp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Conjugate gradient with line minimization performed
!with the improoved interpolation methond
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Conjugate gradient optimization with line minimization using improved interpolation method
!> @details This subroutine implements a conjugate gradient algorithm with line minimization
!> performed using an improved interpolation method. It calculates the gradient, updates
!> the conjugate direction using the Polak-Ribière formula, and performs line minimization
!> along the conjugate direction to find the optimal step size.
!>
!> @param[in] size1 Dimension of matrices in the previous basis
!> @param[in] size2 Dimension of matrices in the new basis  
!> @param[in] type_lambda Type of lambda matrix symmetry (1=hermitian, 2=symmetric, 3=non-symmetric)
!> @param[in] nnozero_c Size of the vector nozero_c and jbradet
!> @param[in] nozero_c Vector containing information to obtain the whole matrix lambda
!> @param[in] jbradet Vector pointing to independent variable number associated with symmetry
!> @param[in] nsym Total number of independent variables after applying symmetries
!> @param[in] count_sym Vector counting the recurrence of every symmetry
!> @param[in] img Input to understand if wave function is real (1) or complex (2)
!> @param[in] L1mod Module of the previous wave function
!> @param[in] L2mod Module of the new wave function
!> @param[in,out] L2 Lambda matrix in the new basis
!> @param[in] Ob Matrix Ob defined in parbcs
!> @param[in] SL2,SR2 Overlap matrices of the new basis
!> @param[in] A Variable for fast calculation
!> @param[in] dL2 Matrix of the components of the derivatives
!> @param[in] optimize Information on which matrix elements to optimize
!> @param[in,out] prevH Previous conjugate direction
!> @param[in,out] prevG Previous gradient
!> @param[in] prec Precision required for the overlap calculation
!> @param[out] lambda Factor that multiplies the direction
!> @param[in] rank MPI rank
!> @param[in] nprocu Number of MPI processes
!> @param[in] mpi_comm_world MPI communicator
!>
!> @note This subroutine uses the Polak-Ribière formula for updating the conjugate direction:
!>       gamma = (g_new - g_old) · g_new / (g_old · g_old)
!>       h_new = g_new + gamma * h_old
!>
!> @see lineminDB, symm_deriv
subroutine minDB(size1, size2, type_lambda, nnozero_c                    &
           &, nozero_c, jbradet, nsym, count_sym, img                    &
           &, L1mod, L2mod, L2, Ob, SL2, SR2, A, dL2                     &
           &, optimize, prevH, prevG, prec, lambda, rank                 &
           &, nprocu, mpi_comm_world)
    implicit none
    integer :: size2, img, size1, i, j, rank, nprocu, mpi_comm_world, nnozero_c, ix, iy, nsym, type_lambda, ierr
    !mac prec is the precision that we can expect from the zero of the derivative. If the single number precision is
    !10E-16 (to be safe) we have to consider that the zero of the derivative is the zero of all this noise also, so
    ! I expect a macprec=10E-16*size2
    real(8) :: L1mod, L2mod, delta, lambda, A, renorm, prec, Lauxmod, Z, macprec
    !It distinguishes between the real and the complex case
    real(8) :: gamma, gamman, gammad
    real(8), dimension(size2*img, size2) :: L2, Ob, SL2, SR2, dL2
    real(8), dimension(img*nsym) :: prevG, prevH
    integer, dimension(nsym) :: count_sym
    integer, dimension(nnozero_c) :: jbradet, nozero_c
    real(8), allocatable :: h(:), gaux(:), g(:), Op(:, :)
    logical, dimension(size2*img, size2) :: optimize
    allocate (h(img*nsym), gaux(nsym*img), g(nsym*img)) !, Op(size2*img,size2))

    g = 0.d0

    macprec = 10e-15*size2

    do i = 1, size2
        do j = 1, size2*img
            if (optimize(j, i) .eqv. .false.) dL2(j, i) = 0.d0
        end do
    end do

    !Taking only the elements reduced by symmetries of the derivatives
    call symm_deriv(img, size2, type_lambda, nnozero_c, nsym, dL2, g, nozero_c, jbradet, count_sym, rank)

    !  call numerical_g(img, size1, size2, type_lambda, nnozero_c, nsym, L2, g, nozero_c, jbradet, count_sym, rank, &
    !       L1mod, L2mod, Ob, SL2, SR2,   optimize,  nprocu, mpi_comm_world )

    gaux(:) = g(:) - prevG(:)
    gamman = 0.d0
    gammad = 0.d0

    !Calculation of the new conjugate direction
    do i = 1, img*nsym
        gamman = gamman + g(i)*Gaux(i)
        gammad = gammad + prevg(i)*prevg(i)
    end do
    gamma = gamman/gammad
    if (gammad .lt. 1e-14) gamma = 0.d0
    h(:) = g(:) + gamma*prevH(:)

    call lineminDB(size1, size2, type_lambda, nnozero_c, nozero_c, jbradet, nsym, count_sym, &
                   img, L1mod, L2mod, L2, Ob, SL2, SR2, g, h, lambda, prec, optimize, macprec, rank, nprocu, mpi_comm_world)
    prevH(:) = h(:)
    prevG(:) = g(:)

    deallocate (gaux, h, g)
end subroutine minDB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Maximization with the hybrid method along the direction h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Line minimization using hybrid method along the conjugate direction
!> @details This subroutine performs line minimization along the conjugate direction h
!> using a hybrid algorithm that combines linear interpolation and bisection methods.
!> It finds the optimal step size lambda that maximizes the overlap between the two
!> wave functions by finding the zero of the derivative along the search direction.
!>
!> The algorithm works as follows:
!> 1. Find two points where the derivative changes sign (bracketing)
!> 2. Use a hybrid method combining linear interpolation and bisection to find the zero
!> 3. Update the lambda matrix with the optimal step size
!>
!> @param[in] size1 Dimension of matrices in the previous basis
!> @param[in] size2 Dimension of matrices in the new basis
!> @param[in] type_lambda Type of lambda matrix symmetry (1=hermitian, 2=symmetric, 3=non-symmetric)
!> @param[in] nnozero_c Size of the vector nozero_c and jbradet
!> @param[in] nozero_c Vector containing information to obtain the whole matrix lambda
!> @param[in] jbradet Vector pointing to independent variable number associated with symmetry
!> @param[in] nsym Total number of independent variables after applying symmetries
!> @param[in] count_sym Vector counting the recurrence of every symmetry
!> @param[in] img Input to understand if wave function is real (1) or complex (2)
!> @param[in] L1mod Module of the previous wave function
!> @param[in,out] L2mod Module of the new wave function
!> @param[in,out] L2 Lambda matrix in the new basis
!> @param[in] Ob Matrix Ob defined in parbcs
!> @param[in] SL2,SR2 Overlap matrices of the new basis
!> @param[in] g Current gradient vector
!> @param[in] h Conjugate direction vector
!> @param[out] lambda Optimal step size along the direction h
!> @param[in] prec Precision required for the overlap calculation
!> @param[in] optimize Information on which matrix elements to optimize
!> @param[in] macprec Machine precision for derivative calculations
!> @param[in] rank MPI rank
!> @param[in] nprocu Number of MPI processes
!> @param[in] mpi_comm_world MPI communicator
!>
!> @note The hybrid algorithm uses linear interpolation when the new point is far from
!>       the boundaries, and switches to bisection when it gets too close to avoid
!>       slow convergence near the boundaries.
!>
!> @see translate_L, calcOp, calcA, calcZ, calcDeriv, symm_deriv
subroutine lineminDB(size1, size2, type_lambda, nnozero_c, nozero_c, jbradet, nsym, count_sym, img, L1mod, &
                     L2mod, L2, Ob, SL2, SR2, g, h, lambda, prec, optimize, macprec, rank, nprocu, mpi_comm_world)
    implicit none
    integer :: size1, size2, img, step, i, j, k, rank, nprocu, mpi_comm_world, ierr, debug, nsym, nnozero_c, ix, iy, type_lambda
    !lt1 and lt2 are the lambda trial
    logical, dimension(size2*img, size2) :: optimize
    real(8) :: L1mod, L2mod, Lauxmod, lambda, prec, dg1, dgnew, dg2, lt1, lt2, newlt, aux, distlt1, distlt2, macprec, Z, Zaux
    real(8), dimension(size2*img, size2) :: L2, Ob, SL2, SR2
    real(8), dimension(img*nnozero_c) :: g, h
    integer, dimension(nsym) :: count_sym
    integer, dimension(nnozero_c) :: jbradet, nozero_c
    real(8), allocatable :: g1(:), g2(:), newg(:), Op(:, :), L2aux(:, :), dL2aux(:, :)
    !Control is on whene no further minimization is needed
    logical :: control

    allocate (g1(img*nsym), g2(img*nsym), newg(img*nsym), Op(size2*img, size2) &
              , L2aux(size2*img, size2), dL2aux(size2*img, size2))
    control = .false.

    g1(:) = g(:)
    dg1 = 0.d0
    dgnew = 0.d0
    lt1 = 0.d0
    do i = 1, img*nsym
        dg1 = dg1 + g1(i)*h(i)
    end do
    if (dg1 .lt. macprec .and. rank .eq. 0) write (6, *) i, ' Warning dg1 too small  ', dg1

    Z = 0
    !Changing the direction of h if necessary (It shouldn't be necessary, only for check)
    if (dg1 .lt. 0 .and. rank .eq. 0) write (6, *) i, ' Warning h g <0 '

    !Initialization of the calculation, finding a point where g2 is negative
    do i = 1, 100
        lt2 = 5.d0**i

        !Calculating the value of the matrix at lt2
        call translate_L(img, size2, type_lambda, nnozero_c, nozero_c, jbradet, L2, L2aux, h, lt2)
        call attach_phase2det(.true., L2aux)
        call calcOp(size2, img, SR2, L2aux, SL2, Op, rank, nprocu, mpi_comm_world)
        call calcA(size2, img, L2aux, Lauxmod, Op, rank, nprocu, mpi_comm_world)
        call calcDeriv(size2, img, Ob, Op, L1mod, Lauxmod, L2aux, dL2aux, optimize, aux, rank, nprocu, mpi_comm_world)
        call calcZ(size1, size2, img, L2aux, Ob, L1mod, Lauxmod, Z, rank, nprocu, mpi_comm_world)
        !Removing the phase from the derivative
        call attach_phase2det(.false., dL2aux)

        dg2 = 0.d0
        !calculating g2

        call symm_deriv(img, size2, type_lambda, nnozero_c, nsym, dL2aux, g2, nozero_c, jbradet, count_sym, rank)
        !     call numerical_g(img, size1, size2, type_lambda, nnozero_c, nsym, L2, g, nozero_c, jbradet, count_sym, rank, &
        !          L1mod, L2mod, Ob, SL2, SR2,   optimize,  nprocu, mpi_comm_world )

        do k = 1, nsym*img
            dg2 = dg2 + g2(k)*h(k)
        end do
        !Checking for negative value of the derivative
        if (dg2 < 0) then
            exit
        else if (abs(dg2) .lt. macprec) then
            control = .true.
            exit
        else
            dg1 = dg2
            g1(:) = g2(:)
            lt1 = lt2
        end if
    end do

    !  if (rank.eq.0) write(6,*) "Passi qui 1", i, dg1, dg2, Z

    !Checking if a 0 has already been founded, in case exiting
    if (control) then
        lambda = lt2
        L2mod = Lauxmod
        aux = (L1mod/L2mod)**(0.5)
        call attach_phase2det(.false., L2aux)
        L2(:, :) = L2aux(:, :)*aux
    else

        !Once we have the two values, we can start the newton/bisection
        do i = 1, 1000
            !Cheching if newton is ok
            newlt = lt1 - dg1*(lt2 - lt1)/(dg2 - dg1)
            !Checking if the newlt and applying an hybrid algorithm between the
            !linear interpolation and the bisection (so much faster!)
            distlt1 = abs((lt1 - newlt)/(lt1 - lt2))
            distlt2 = 1.d0 - distlt1
            if (distlt1 .lt. 0.1d0) then
                newlt = lt1 + (lt2 - lt1)/5.d0
            else if (distlt2 .lt. 0.1d0) then
                newlt = lt1 + (lt2 - lt1)/1.25d0
            end if
            !       newlt=(lt1+lt2)/2
            !        newlt=lt1+(lt2-lt1)*i/1000

            !Calculating the value of the matrix at ltnew
            call translate_L(img, size2, type_lambda, nnozero_c, nozero_c, jbradet, L2, L2aux, h, newlt)

            Zaux = Z

            !Calculating the new derivative in ltnew
            call attach_phase2det(.true., L2aux)
            call calcOp(size2, img, SR2, L2aux, SL2, Op, rank, nprocu, mpi_comm_world)
            call calcA(size2, img, L2aux, Lauxmod, Op, rank, nprocu, mpi_comm_world)
            call calcZ(size1, size2, img, L2aux, Ob, L1mod, Lauxmod, Z, rank, nprocu, mpi_comm_world)
            call calcDeriv(size2, img, Ob, Op, L1mod, Lauxmod, L2aux, dL2aux, optimize, aux, rank, nprocu, mpi_comm_world)
            call attach_phase2det(.false., dL2aux)

            dgnew = 0.d0
            !calculating newg
            call symm_deriv(img, size2, type_lambda, nnozero_c, nsym, dL2aux, newg, nozero_c, jbradet, count_sym, rank)
            !call numerical_g(img, size1, size2, type_lambda, nnozero_c, nsym, L2, g, nozero_c, jbradet, count_sym, rank, &
            !     L1mod, L2mod, Ob, SL2, SR2,   optimize,  nprocu, mpi_comm_world )
            do k = 1, nsym*img
                dgnew = dgnew + newg(k)*h(k)
            end do
            !        if (rank.eq.0) write (6,*) i, newlt, Z, dgnew

            !Checking the point for the substitution
            if (dgnew > 0) then
                !Checking if concluding or not
                if (dgnew < macprec) then
                    lambda = newlt
                    L2mod = Lauxmod
                    aux = (L1mod/L2mod)**(0.5)
                    call attach_phase2det(.false., L2aux)
                    L2(:, :) = L2aux(:, :)*aux
                    exit
                else
                    lt1 = newlt
                    dg1 = dgnew
                end if
            else
                !Checking if concluding or not
                if (abs(dgnew) < macprec) then
                    lambda = newlt
                    L2mod = Lauxmod
                    aux = (L1mod/L2mod)**(0.5)
                    call attach_phase2det(.false., L2aux)
                    L2(:, :) = L2aux(:, :)*aux
                    exit
                else
                    lt2 = newlt
                    dg2 = dgnew
                end if
            end if
        end do
    end if

    !  if (rank.eq.0) write(6,*) "Passi qui 2", i, dg1, dg2, Z

    !  if(rank.eq.0) write(6,*) ' # Iteration lt2 =',i,lt1,lt2,dg1,dg2, dgnew, Z
    !  if(rank.eq.0) write(6,*) ' # Iteration min =',i
    deallocate (g1, g2, newg, Op, L2aux, dL2aux)
    !  if(i.gt.100) then
    !  call mpi_finalize(ierr)
    !  stop
    !  endif

end subroutine lineminDB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Numerical derivatives with symmetries (without translate_L)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Calculate numerical derivatives with symmetries (without translate_L)
!> @details This subroutine calculates numerical derivatives of the overlap with respect
!> to the lambda matrix elements using finite differences. It handles both real and complex
!> cases and respects the symmetry constraints of the lambda matrix.
!>
!> The numerical derivative is calculated as:
!> ∂Z/∂λᵢⱼ ≈ (Z(λᵢⱼ + ε) - Z(λᵢⱼ)) / ε
!>
!> where ε is a small step size (default: 0.00001).
!>
!> @param[in] img Input to understand if wave function is real (1) or complex (2)
!> @param[in] size1 Dimension of matrices in the previous basis
!> @param[in] size2 Dimension of matrices in the new basis
!> @param[in] type_lambda Type of lambda matrix symmetry (1=hermitian, 2=symmetric, 3=non-symmetric)
!> @param[in] nnozero_c Size of the vector nozero_c and jbradet
!> @param[in] nsym Total number of independent variables after applying symmetries
!> @param[in] L2 Lambda matrix in the new basis
!> @param[out] g Gradient vector (symmetry-reduced)
!> @param[in] nozero_c Vector containing information to obtain the whole matrix lambda
!> @param[in] jbradet Vector pointing to independent variable number associated with symmetry
!> @param[in] count_sym Vector counting the recurrence of every symmetry
!> @param[in] rank MPI rank
!> @param[in] L1mod Module of the previous wave function
!> @param[in] L2mod Module of the new wave function
!> @param[in] Ob Matrix Ob defined in parbcs
!> @param[in] SL2,SR2 Overlap matrices of the new basis
!> @param[in] optimize Information on which matrix elements to optimize
!> @param[in] nprocu Number of MPI processes
!> @param[in] mpi_comm_world MPI communicator
!>
!> @note This subroutine uses finite differences to calculate derivatives, which is
!>       more robust but computationally expensive compared to analytical derivatives.
!>       It is mainly used for testing and verification purposes.
!>
!> @see calcOp, calcA, calcZ, symm_deriv
subroutine numerical_g(img, size1, size2, type_lambda, nnozero_c, nsym, L2, g, nozero_c, jbradet, count_sym, rank, &
                       L1mod, L2mod, Ob, SL2, SR2, optimize, nprocu, mpi_comm_world)
    implicit none
    integer :: img, size1, size2, type_lambda, nnozero_c, nsym, rank, i, j, ix, iy, aux, nprocu, mpi_comm_world
    logical, dimension(size2*img, size2) :: optimize
    integer, dimension(nnozero_c) :: nozero_c, jbradet
    integer, dimension(nsym) :: count_sym
    real(8), dimension(nsym*img) :: g
    real(8), dimension(size2*img, size2) :: L2, Ob, SL2, SR2
    real(8) :: L1mod, L2mod, L2auxmod, Z, Zaux, step
    real(8), allocatable :: L2aux(:, :), Op(:, :)

    allocate (L2aux(size2*img, size2), Op(size2*img, size2))
    L2aux = L2
    step = 0.00001d0
    !Initial calculation of the Z
    call calcOp(size2, img, SR2, L2, SL2, Op, rank, nprocu, mpi_comm_world)
    call calcA(size2, img, L2, L2auxmod, Op, rank, nprocu, mpi_comm_world)
    call calcZ(size1, size2, img, L2, Ob, L1mod, L2auxmod, Z, rank, nprocu, mpi_comm_world)
    g = 0.d0

    !Preparing the Lambdaaux
    do i = 1, nsym
        L2aux(:, :) = L2(:, :)
        !Real Case
        if (img .eq. 1) then
            do j = 1, nnozero_c
                if (abs(jbradet(j)) .eq. i) then
                    iy = (nozero_c(j) - 1)/size2 + 1
                    ix = nozero_c(j) - (iy - 1)*size2
                    if (optimize(ix, iy)) then
                        L2aux(ix, iy) = L2(ix, iy) + step*sign(1, jbradet(j))
                        if ((type_lambda .ne. 3) .and. (ix .ne. iy)) L2aux(iy, ix) = L2(iy, ix) + step*sign(1, jbradet(j))
                    end if
                end if
            end do
            call calcOp(size2, img, SR2, L2aux, SL2, Op, rank, nprocu, mpi_comm_world)
            call calcA(size2, img, L2aux, L2auxmod, Op, rank, nprocu, mpi_comm_world)
            call calcZ(size1, size2, img, L2aux, Ob, L1mod, L2auxmod, Zaux, rank, nprocu, mpi_comm_world)
            g(i) = (Zaux - Z)/step
        end if
        !Complex case

        if (img .eq. 2) then
            !Real Part
            do j = 1, nnozero_c
                if (abs(jbradet(j)) .eq. i) then
                    iy = (nozero_c(j) - 1)/size2 + 1
                    ix = nozero_c(j) - (iy - 1)*size2
                    if (optimize(2*ix - 1, iy)) then
                        L2aux(2*ix - 1, iy) = L2(2*ix - 1, iy) + step*sign(1, jbradet(j))
                      if ((type_lambda .ne. 3)&
                         & .and. (ix .ne. iy))&
                         & L2aux(2*iy - 1, ix) = L2(2*iy - 1, ix) + step*sign(1, jbradet(j))
                    end if
                end if
            end do

            call calcOp(size2, img, SR2, L2aux, SL2, Op, rank, nprocu, mpi_comm_world)
            call calcA(size2, img, L2aux, L2auxmod, Op, rank, nprocu, mpi_comm_world)
            call calcZ(size1, size2, img, L2aux, Ob, L1mod, L2auxmod, Zaux, rank, nprocu, mpi_comm_world)
            g(2*i - 1) = (Zaux - Z)/step
            L2aux = L2
            !Immaginary Part
            do j = 1, nnozero_c
                if (abs(jbradet(j)) .eq. i) then
                    iy = (nozero_c(j) - 1)/size2 + 1
                    ix = nozero_c(j) - (iy - 1)*size2
                    if (optimize(2*ix, iy)) then
                        L2aux(2*ix, iy) = L2(2*ix, iy) + step*sign(1, jbradet(j))
                        if ((type_lambda .eq. 1) .and. (ix .ne. iy)) L2aux(2*iy, ix) = L2(2*iy, ix) - step*sign(1, jbradet(j))
                        if ((type_lambda .eq. 2) .and. (ix .ne. iy)) L2aux(2*iy, ix) = L2(2*iy, ix) + step*sign(1, jbradet(j))
                    end if
                end if
            end do
            call calcOp(size2, img, SR2, L2aux, SL2, Op, rank, nprocu, mpi_comm_world)
            call calcA(size2, img, L2aux, L2auxmod, Op, rank, nprocu, mpi_comm_world)
            call calcZ(size1, size2, img, L2aux, Ob, L1mod, L2auxmod, Zaux, rank, nprocu, mpi_comm_world)
            g(2*i) = (Zaux - Z)/step
        end if

    end do

    deallocate (L2aux, Op)

end subroutine numerical_g

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Numerical derivatives with symmetries (with translate_L)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Calculate numerical derivatives with symmetries using translate_L
!> @details This subroutine calculates numerical derivatives of the overlap with respect
!> to the lambda matrix elements using finite differences and the translate_L function.
!> It handles both real and complex cases and respects the symmetry constraints of the lambda matrix.
!>
!> The numerical derivative is calculated as:
!> ∂Z/∂λᵢⱼ ≈ (Z(λᵢⱼ + ε·h) - Z(λᵢⱼ)) / ε
!>
!> where ε is a small step size (default: 0.00001) and h is the direction vector.
!>
!> @param[in] img Input to understand if wave function is real (1) or complex (2)
!> @param[in] size1 Dimension of matrices in the previous basis
!> @param[in] size2 Dimension of matrices in the new basis
!> @param[in] type_lambda Type of lambda matrix symmetry (1=hermitian, 2=symmetric, 3=non-symmetric)
!> @param[in] nnozero_c Size of the vector nozero_c and jbradet
!> @param[in] nsym Total number of independent variables after applying symmetries
!> @param[in] L2 Lambda matrix in the new basis
!> @param[out] g Gradient vector (symmetry-reduced)
!> @param[in] nozero_c Vector containing information to obtain the whole matrix lambda
!> @param[in] jbradet Vector pointing to independent variable number associated with symmetry
!> @param[in] count_sym Vector counting the recurrence of every symmetry
!> @param[in] rank MPI rank
!> @param[in] L1mod Module of the previous wave function
!> @param[in] L2mod Module of the new wave function
!> @param[in] Ob Matrix Ob defined in parbcs
!> @param[in] SL2,SR2 Overlap matrices of the new basis
!> @param[in] optimize Information on which matrix elements to optimize
!> @param[in] nprocu Number of MPI processes
!> @param[in] mpi_comm_world MPI communicator
!>
!> @note This subroutine uses the translate_L function to generate perturbed matrices,
!>       which ensures that all symmetry constraints are properly maintained during
!>       the finite difference calculation.
!>
!> @see translate_L, calcOp, calcA, calcZ
subroutine numerical_g2(img, size1, size2, type_lambda, nnozero_c, nsym, L2, g, nozero_c, jbradet, count_sym, rank, &
                        L1mod, L2mod, Ob, SL2, SR2, optimize, nprocu, mpi_comm_world)
    implicit none
    integer :: img, size1, size2, type_lambda, nnozero_c, nsym, rank, i, j, ix, iy, aux, nprocu, mpi_comm_world
    logical, dimension(size2*img, size2) :: optimize
    integer, dimension(nnozero_c) :: nozero_c, jbradet
    integer, dimension(nsym) :: count_sym
    real(8), dimension(nsym*img) :: g
    real(8), dimension(size2*img, size2) :: L2, Ob, SL2, SR2
    real(8) :: L1mod, L2mod, L2auxmod, Z, Zaux, step, lt
    real(8), allocatable :: L2aux(:, :), Op(:, :), h(:)

    allocate (L2aux(size2*img, size2), Op(size2*img, size2), h(img*nsym))
    L2aux = L2
    step = 0.00001d0
    lt = 1.d0
    !Initial calculation of the Z
    call calcOp(size2, img, SR2, L2, SL2, Op, rank, nprocu, mpi_comm_world)
    call calcA(size2, img, L2, L2auxmod, Op, rank, nprocu, mpi_comm_world)
    call calcZ(size1, size2, img, L2, Ob, L1mod, L2auxmod, Z, rank, nprocu, mpi_comm_world)
    g = 0.d0

    !Preparing the Lambdaaux
    do i = 1, nsym
        h = 0.d0
        L2aux(:, :) = L2(:, :)
        !Real Case
        if (img .eq. 1) then
            h(i) = step
            call translate_L(img, size2, type_lambda, nnozero_c, nozero_c, jbradet, L2, L2aux, h, lt)
            call calcOp(size2, img, SR2, L2aux, SL2, Op, rank, nprocu, mpi_comm_world)
            call calcA(size2, img, L2aux, L2auxmod, Op, rank, nprocu, mpi_comm_world)
            call calcZ(size1, size2, img, L2aux, Ob, L1mod, L2auxmod, Zaux, rank, nprocu, mpi_comm_world)
            g(i) = (Zaux - Z)/step
        end if
        !Complex case

        if (img .eq. 2) then
            !Real Part
            h(2*i - 1) = step
            call translate_L(img, size2, type_lambda, nnozero_c, nozero_c, jbradet, L2, L2aux, h, lt)
            call calcOp(size2, img, SR2, L2aux, SL2, Op, rank, nprocu, mpi_comm_world)
            call calcA(size2, img, L2aux, L2auxmod, Op, rank, nprocu, mpi_comm_world)
            call calcZ(size1, size2, img, L2aux, Ob, L1mod, L2auxmod, Zaux, rank, nprocu, mpi_comm_world)
            g(2*i - 1) = (Zaux - Z)/step
            L2aux = L2
            h = 0.d0
            !Immaginary Part
            h(2*i) = step
            call translate_L(img, size2, type_lambda, nnozero_c, nozero_c, jbradet, L2, L2aux, h, lt)
            call calcOp(size2, img, SR2, L2aux, SL2, Op, rank, nprocu, mpi_comm_world)
            call calcA(size2, img, L2aux, L2auxmod, Op, rank, nprocu, mpi_comm_world)
            call calcZ(size1, size2, img, L2aux, Ob, L1mod, L2auxmod, Zaux, rank, nprocu, mpi_comm_world)
            g(2*i) = (Zaux - Z)/step
        end if

    end do

    deallocate (L2aux, Op, h)

end subroutine numerical_g2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculating the g vector made only by the irreducible elements, averaging the derivatives (maybe not necessary),
!using only the upperdiagonal matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Calculate symmetry-reduced gradient vector from full derivative matrix
!> @details This subroutine calculates the gradient vector g made only by the irreducible
!> elements, averaging the derivatives and using only the upper diagonal matrix.
!> It handles both real and complex cases, as well as different symmetry types
!> (hermitian, symmetric, non-symmetric) and particle-hole symmetry (ipf).
!>
!> The subroutine performs the following operations:
!> 1. Extracts derivatives for symmetry-equivalent elements
!> 2. Applies appropriate signs based on jbradet
!> 3. Handles particle-hole symmetry constraints
!> 4. Averages contributions from equivalent elements
!>
!> @param[in] img Input to understand if wave function is real (1) or complex (2)
!> @param[in] size2 Dimension of matrices in the new basis
!> @param[in] type_lambda Type of lambda matrix symmetry (1=hermitian, 2=symmetric, 3=non-symmetric)
!> @param[in] nnozero_c Size of the vector nozero_c and jbradet
!> @param[in] nsym Total number of independent variables after applying symmetries
!> @param[in] dL2 Full derivative matrix (size2*img × size2)
!> @param[out] g Symmetry-reduced gradient vector (img*nsym)
!> @param[in] nozero_c Vector containing information to obtain the whole matrix lambda
!> @param[in] jbradet Vector pointing to independent variable number associated with symmetry
!> @param[in] count_sym Vector counting the recurrence of every symmetry
!> @param[in] rank MPI rank
!>
!> @note The subroutine handles different cases:
!>       - Real vs complex wave functions (img = 1 or 2)
!>       - Particle-hole symmetry (ipf = 1 or 2)
!>       - Different matrix symmetries (type_lambda = 1, 2, or 3)
!>
!> @see translate_L, numerical_g, numerical_g2
subroutine symm_deriv(img, size2, type_lambda, nnozero_c, nsym, dL2, g, nozero_c, jbradet, count_sym, rank)
    use constants, only: ipf
    implicit none
    integer :: ix, iy, img, size2, nnozero_c, nsym, i, type_lambda, rank
    real(8), dimension(img*size2, size2) :: dL2
    real(8), dimension(img*nsym) :: g
    integer, dimension(nnozero_c) :: jbradet, nozero_c
    integer, dimension(nsym) :: count_sym

    !Distinguishing the imaginary and the real case and considering if the elements are symmetric or antysimmetric
    g = 0.d0
    if (img .eq. 1) then
        do i = 1, nnozero_c
            iy = (nozero_c(i) - 1)/size2 + 1
            ix = nozero_c(i) - (iy - 1)*size2
            if (jbradet(i) .ne. 0) then
                if (ipf .eq. 1) then
                    g(abs(jbradet(i))) = g(abs(jbradet(i))) + dL2(ix, iy)*sign(1, jbradet(i))
              if ((ix .ne. iy)&
                 & .and. (type_lambda .ne. 3)) g(abs(jbradet(i)))&
                 & = g(abs(jbradet(i))) + dL2(iy, ix)*sign(1, jbradet(i))
                else
                    g(abs(jbradet(i))) = g(abs(jbradet(i))) + dL2(ix, iy)*sign(1, jbradet(i))
                    !              g(abs(jbradet(i)))=g(abs(jbradet(i)))- dL2(iy,ix)*sign(1,jbradet(i))
                    if (type_lambda .ne. 3) then
                        if (iy .gt. size2/2) then
                            g(abs(jbradet(i))) = g(abs(jbradet(i))) - dL2(ix + size2/2, iy - size2/2)*sign(1, jbradet(i))
                        else
                            g(abs(jbradet(i))) = g(abs(jbradet(i))) - dL2(ix + size2/2, iy + size2/2)*sign(1, jbradet(i))
                        end if
                        !                 g(abs(jbradet(i)))=g(abs(jbradet(i)))+ dL2(iy-size2/2,ix+size2/2)*sign(1,jbradet(i))
                    end if
                end if
            end if
        end do
    else
        do i = 1, nnozero_c
            iy = (nozero_c(i) - 1)/size2 + 1
            ix = nozero_c(i) - (iy - 1)*size2
            if (jbradet(i) .ne. 0) then
                if (ipf .eq. 1) then
                    g(2*abs(jbradet(i)) - 1) = g(abs(jbradet(i))*2 - 1) + dL2(2*ix - 1, iy)*sign(1, jbradet(i))
                    g(2*abs(jbradet(i))) = g(2*abs(jbradet(i))) + dL2(2*ix, iy)*sign(1, jbradet(i))
                    if ((type_lambda .ne. 3) .and. (ix .ne. iy)) then
                        g(2*abs(jbradet(i)) - 1) = g(abs(jbradet(i))*2 - 1) + dL2(2*iy - 1, ix)*sign(1, jbradet(i))
                        if (type_lambda .eq. 1) g(2*abs(jbradet(i))) = g(2*abs(jbradet(i))) - dL2(2*iy, ix)*sign(1, jbradet(i))
                        if (type_lambda .eq. 2) g(2*abs(jbradet(i))) = g(2*abs(jbradet(i))) + dL2(2*iy, ix)*sign(1, jbradet(i))
                    end if
                else
                    g(2*abs(jbradet(i)) - 1) = g(abs(jbradet(i))*2 - 1) + dL2(2*ix - 1, iy)*sign(1, jbradet(i))
                    g(2*abs(jbradet(i))) = g(2*abs(jbradet(i))) + dL2(2*ix, iy)*sign(1, jbradet(i))
                    !              g(2*abs(jbradet(i))-1)=g(abs(jbradet(i))*2-1)- dL2(2*iy-1,ix)*sign(1,jbradet(i))
                    !              g(2*abs(jbradet(i)))=g(2*abs(jbradet(i)))- dL2(2*iy,ix)*sign(1,jbradet(i))

                    if (type_lambda .ne. 3) then
                        if (iy .gt. size2/2) then
                   g(2*abs(jbradet(i)) - 1) = g(abs(jbradet(i))*2 - 1)&
                      & + dL2(2*(iy - size2/2) - 1, ix + size2/2)*sign(1, jbradet(i))
                            if (type_lambda .eq. 1) then
                             g(2*abs(jbradet(i))) = g(2*abs(jbradet(i)*2))&
                                & - dL2(2*(iy - size2/2), ix + size2/2)*sign(1, jbradet(i))
                            else if (type_lambda .eq. 2) then
                             g(2*abs(jbradet(i))) = g(2*abs(jbradet(i)*2))&
                                & + dL2(2*(iy - size2/2), ix + size2/2)*sign(1, jbradet(i))
                            end if
                        else
                   g(2*abs(jbradet(i)) - 1) = g(abs(jbradet(i))*2 - 1)&
                      & - dL2(2*(iy + size2/2) - 1, ix + size2/2)*sign(1, jbradet(i))
                            if (type_lambda .eq. 1) then
                             g(2*abs(jbradet(i))) = g(2*abs(jbradet(i)*2))&
                                & + dL2(2*(iy + size2/2), ix + size2/2)*sign(1, jbradet(i))
                            else if (type_lambda .eq. 2) then
                             g(2*abs(jbradet(i))) = g(2*abs(jbradet(i)*2))&
                                & - dL2(2*(iy + size2/2), ix + size2/2)*sign(1, jbradet(i))
                            end if
                        end if
                    end if

                end if

            end if
        end do
    end if
end subroutine symm_deriv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculating the L2 in the new position given the vector of the direction in the symmetry reduced rappresentation
!and the movement in the direction, it uses also the given proprieties of lambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Transform lambda matrix along a given direction vector
!> @details This subroutine calculates the lambda matrix L2 in a new position given
!> the direction vector h in the symmetry-reduced representation and the movement
!> parameter lt. It uses the given properties of lambda (symmetry type) to ensure
!> the transformed matrix maintains the correct symmetry constraints.
!>
!> The transformation is performed as:
!> L2_new = L2_old + lt * h_transformed
!>
!> where h_transformed is the direction vector expanded to the full matrix space
!> according to the symmetry constraints.
!>
!> @param[in] img Input to understand if wave function is real (1) or complex (2)
!> @param[in] size2 Dimension of matrices in the new basis
!> @param[in] type_symm Type of lambda matrix symmetry (1=hermitian, 2=symmetric, 3=non-symmetric)
!> @param[in] nnozero_c Size of the vector nozero_c and jbradet
!> @param[in] nozero_c Vector containing information to obtain the whole matrix lambda
!> @param[in] jbradet Vector pointing to independent variable number associated with symmetry
!> @param[in] L2 Original lambda matrix (size2*img × size2)
!> @param[out] L2aux Transformed lambda matrix (size2*img × size2)
!> @param[in] h Direction vector in symmetry-reduced representation (img*nnozero_c)
!> @param[in] lt Step size along the direction
!>
!> @note The subroutine handles different cases:
!>       - Real vs complex wave functions (img = 1 or 2)
!>       - Particle-hole symmetry (ipf = 1 or 2)
!>       - Different matrix symmetries (type_symm = 1, 2, or 3)
!>       - Hermitian (type_symm = 1): L(i,j) = L*(j,i)
!>       - Symmetric (type_symm = 2): L(i,j) = L(j,i)
!>       - Non-symmetric (type_symm = 3): no constraints
!>
!> @see symm_deriv, numerical_g2
subroutine translate_L(img, size2, type_symm, nnozero_c, nozero_c, jbradet, L2, L2aux, h, lt)
    use constants, only: ipf
    implicit none
    integer :: img, size2, type_symm, nnozero_c, ix, iy, j
    integer, dimension(nnozero_c) :: jbradet, nozero_c
    real(8) :: lt
    real(8), dimension(size2*img, size2) :: L2, L2aux
    real(8), dimension(img*nnozero_c) :: h

    L2aux = 0.d0
    do j = 1, nnozero_c
        if (jbradet(j) .ne. 0) then
            iy = (nozero_c(j) - 1)/size2 + 1
            ix = nozero_c(j) - (iy - 1)*size2
            if (img .eq. 1) then
                if (ipf .eq. 1) then
                    L2aux(ix, iy) = L2(ix, iy) + lt*h(abs(jbradet(j)))*sign(1, jbradet(j))
                    !Calculating also the symmetric element for off diagonal terms
                    if ((ix .ne. iy) .and. (type_symm .ne. 3)) then
                        L2aux(iy, ix) = L2aux(ix, iy)
                    end if
                else
                    L2aux(ix, iy) = L2(ix, iy) + lt*h(abs(jbradet(j)))*sign(1, jbradet(j))
                    L2aux(iy, ix) = -L2aux(ix, iy)
                    !Calculating also the symmetric element for off diagonal terms
                    if ((ix .ne. iy) .and. (type_symm .ne. 3)) then
                        if (iy .gt. size2/2) then
                            L2aux(iy - size2/2, ix + size2/2) = L2aux(ix, iy)
                            L2aux(ix + size2/2, iy - size2/2) = -L2aux(ix, iy)
                        else
                            L2aux(ix + size2/2, iy + size2/2) = L2aux(ix, iy)
                            L2aux(iy + size2/2, ix + size2/2) = -L2aux(ix, iy)
                        end if
                    end if
                end if
            else
                if (ipf .eq. 1) then
                    L2aux(2*ix - 1, iy) = L2(2*ix - 1, iy) + lt*h(2*abs(jbradet(j)) - 1)*sign(1, jbradet(j))
                    L2aux(2*ix, iy) = L2(2*ix, iy) + lt*h(2*abs(jbradet(j)))*sign(1, jbradet(j))
                    !Calculating also the symmetric element for off diagonal terms
                    if ((type_symm .ne. 3) .and. (ix .ne. iy)) then
                        L2aux(2*iy - 1, ix) = L2aux(2*ix - 1, iy)
                        if (type_symm .eq. 1) then
                            L2aux(2*iy, ix) = -L2aux(2*ix, iy)
                        else
                            L2aux(2*iy, ix) = L2aux(2*ix, iy)
                        end if
                    end if
                else
                    L2aux(2*ix - 1, iy) = L2(2*ix - 1, iy) + lt*h(2*abs(jbradet(j)) - 1)*sign(1, jbradet(j))
                    L2aux(2*ix, iy) = L2(2*ix, iy) + lt*h(2*abs(jbradet(j)))*sign(1, jbradet(j))
                    L2aux(2*iy - 1, ix) = -L2aux(2*ix - 1, iy)
                    L2aux(2*iy, ix) = -L2aux(2*ix, iy)
                    !Calculating also the symmetric element for off diagonal terms
                    if ((type_symm .ne. 3) .and. (ix .ne. iy)) then
                        if (iy .gt. size2/2) then
                            L2aux(2*(iy - size2/2) - 1, ix + size2/2) = L2aux(2*ix - 1, iy)
                            if (type_symm .eq. 1) then
                                L2aux(2*(iy - size2/2), ix + size2/2) = -L2aux(2*ix, iy)
                            else
                                L2aux(2*(iy - size2/2), ix + size2/2) = L2aux(2*ix, iy)
                            end if
                        else
                            L2aux(2*(ix + size2/2) - 1, iy + size2/2) = L2aux(2*ix - 1, iy)
                            if (type_symm .eq. 1) then
                                L2aux(2*(ix + size2/2), iy + size2/2) = -L2aux(2*ix, iy)
                            else
                                L2aux(2*(ix + size2/2), iy + size2/2) = L2aux(2*ix, iy)
                            end if
                        end if
                    end if

                end if
            end if
        end if
    end do

end subroutine translate_L

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Initializing alla the symmetries and the element of the lambda (It is equivalent to the cleanfort.10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Initialize symmetries and elements of the lambda matrix
!> @details This subroutine initializes all the symmetries and elements of the lambda matrix.
!> It is equivalent to the cleanfort.10 operation and ensures that the lambda matrix
!> satisfies all symmetry constraints. The subroutine performs the following operations:
!>
!> 1. Sets up the optimize array to include all symmetry-equivalent elements
!> 2. Counts the recurrence of every symmetry
!> 3. Initializes the lambda matrix with averaged values for symmetry-equivalent elements
!> 4. Applies symmetry constraints to ensure consistency
!>
!> @param[in] size2 Dimension of matrices in the new basis
!> @param[in] img Input to understand if wave function is real (1) or complex (2)
!> @param[in] nnozero_c Size of the vector nozero_c and jbradet
!> @param[in] nsym Total number of independent variables after applying symmetries
!> @param[in] count_sym Vector counting the recurrence of every symmetry
!> @param[in] jbradet Vector pointing to independent variable number associated with symmetry
!> @param[in,out] optimize Information on which matrix elements to optimize
!> @param[in,out] L2 Lambda matrix in the new basis
!> @param[in] nozero_c Vector containing information to obtain the whole matrix lambda
!> @param[in] type_symm Type of lambda matrix symmetry (1=hermitian, 2=symmetric, 3=non-symmetric)
!>
!> @note The subroutine handles different cases:
!>       - Real vs complex wave functions (img = 1 or 2)
!>       - Particle-hole symmetry (ipf = 1 or 2)
!>       - Different matrix symmetries (type_symm = 1, 2, or 3)
!>       - Elements with jbradet = 0 are set to zero and not optimized
!>
!> @see translate_L, symm_deriv
subroutine initialize_symm(size2, img, nnozero_c, nsym, count_sym, jbradet, optimize, L2, nozero_c, type_symm)
    use constants, only: ipf
    implicit none
    integer :: i, j, ix, iy, size2, img, nnozero_c, type_symm, nsym
    integer, dimension(nsym) :: count_sym
    integer, dimension(nnozero_c) :: jbradet, nozero_c
    logical, dimension(size2*img, size2) :: optimize
    real(8), dimension(size2*img, size2) :: L2
    !slambda is the symmetry reduced version of the lambda for the cleanfort10 operation
    real(8), allocatable :: slambda(:)

    allocate (slambda(nsym*img))

    if (ipf .eq. 2) then
        do i = 1, nnozero_c
            iy = (nozero_c(i) - 1)/size2 + 1
            ix = nozero_c(i) - (iy - 1)*size2
            if (img .eq. 1) then
                if (optimize(ix, iy)) then
                    optimize(iy, ix) = .true.
                    if (type_symm .ne. 3) then
                        if (iy .gt. size2/2) then
                            optimize(iy - size2/2, ix + size2/2) = .true.
                            optimize(ix + size2/2, iy - size2/2) = .true.
                        else
                            optimize(iy + size2/2, ix + size2/2) = .true.
                            optimize(ix + size2/2, iy + size2/2) = .true.
                        end if
                    end if
                end if
            else
                if (optimize(2*ix - 1, iy)) then
                    optimize(2*iy - 1, ix) = .true.
                    if (type_symm .ne. 3) then
                        if (iy .gt. size2/2) then
                            optimize(2*(iy - size2/2) - 1, ix + size2/2) = .true.
                            optimize(2*(ix + size2/2) - 1, iy - size2/2) = .true.
                        else
                            optimize(2*(iy + size2/2) - 1, ix + size2/2) = .true.
                            optimize(2*(ix + size2/2) - 1, iy + size2/2) = .true.
                        end if
                    end if
                end if
                if (optimize(2*ix, iy)) then
                    optimize(2*iy, ix) = .true.
                    if (type_symm .ne. 3) then
                        if (iy .gt. size2/2) then
                            optimize(2*(iy - size2/2), ix + size2/2) = .true.
                            optimize(2*(ix + size2/2), iy - size2/2) = .true.
                        else
                            optimize(2*(iy + size2/2), ix + size2/2) = .true.
                            optimize(2*(ix + size2/2), iy + size2/2) = .true.
                        end if
                    end if
                end if

            end if
        end do
    end if

    !counting the recurrency of every symmetry (necessary for the average)
    count_sym = 0
    do i = 1, nnozero_c
        !Here we considerate the case of jbradet=0 that happens for elements =0 for symmetry,
        !if it is the case we do not optimize the relative element
        if (jbradet(i) .ne. 0) then
            count_sym(abs(jbradet(i))) = count_sym(abs(jbradet(i))) + 1
        else
            iy = (nozero_c(i) - 1)/size2 + 1
            ix = nozero_c(i) - (iy - 1)*size2
            if (img .eq. 1) then
                optimize(ix, iy) = .false.
            else
                optimize(2*ix - 1, iy) = .false.
                optimize(2*ix, iy) = .false.
            end if
        end if
    end do

    slambda = 0.d0
    !Initializing the values of the lambda for the lambda
    do i = 1, nnozero_c
        if (jbradet(i) .ne. 0) then
            iy = (nozero_c(i) - 1)/size2 + 1
            ix = nozero_c(i) - (iy - 1)*size2
            if (img .eq. 1) then
                slambda(abs(jbradet(i))) = slambda(abs(jbradet(i))) + L2(ix, iy)*sign(1, jbradet(i))/count_sym(abs(jbradet(i)))
            else
   slambda(2*abs(jbradet(i)) - 1) = slambda(2*abs(jbradet(i)) - 1)&
      & + L2(2*ix - 1, iy)*sign(1, jbradet(i))/count_sym(abs(jbradet(i)))
               slambda(2*abs(jbradet(i))) = slambda(2*abs(jbradet(i)))&
                  & + L2(2*ix, iy)*sign(1, jbradet(i))/count_sym(abs(jbradet(i)))
            end if
        end if
    end do

    !Assigning the right value to everly lambda element considering also the lambda symmetries
    do i = 1, nnozero_c
        if (jbradet(i) .ne. 0) then
            iy = (nozero_c(i) - 1)/size2 + 1
            ix = nozero_c(i) - (iy - 1)*size2
            !Considering the real case
            if (img .eq. 1) then
                if (ipf .eq. 1) then
                    L2(ix, iy) = slambda(abs(jbradet(i)))*sign(1, jbradet(i))
                    if (type_symm .ne. 3) then
                        L2(iy, ix) = slambda(abs(jbradet(i)))*sign(1, jbradet(i))
                    end if
                else
                    L2(ix, iy) = slambda(abs(jbradet(i)))*sign(1, jbradet(i))
                    L2(iy, ix) = -L2(ix, iy)
                    if (type_symm .ne. 3) then
                        if (iy .gt. size2/2) then
                            L2(iy - size2/2, ix + size2/2) = L2(ix, iy)
                            L2(ix + size2/2, iy - size2/2) = -L2(ix, iy)
                        else
                            L2(ix + size2/2, iy + size2/2) = L2(ix, iy)
                            L2(iy + size2/2, ix + size2/2) = -L2(ix, iy)
                        end if
                    end if
                end if
                !Considering the complex case
            else
                if (ipf .eq. 1) then
                    L2(2*ix - 1, iy) = slambda(2*abs(jbradet(i)) - 1)*sign(1, jbradet(i))
                    L2(2*ix, iy) = slambda(2*abs(jbradet(i)))*sign(1, jbradet(i))
                    if ((ix .ne. iy) .and. (type_symm .eq. 1)) then
                        L2(2*iy - 1, ix) = slambda(2*abs(jbradet(i)) - 1)*sign(1, jbradet(i))
                        L2(2*iy, ix) = -slambda(2*abs(jbradet(i)))*sign(1, jbradet(i))
                    else if ((ix .ne. iy) .and. (type_symm .eq. 2)) then
                        L2(2*iy - 1, ix) = slambda(2*abs(jbradet(i)) - 1)*sign(1, jbradet(i))
                        L2(2*iy, ix) = slambda(2*abs(jbradet(i)))*sign(1, jbradet(i))
                    end if
                else
                    L2(2*ix - 1, iy) = slambda(2*abs(jbradet(i)) - 1)*sign(1, jbradet(i))
                    L2(2*ix, iy) = slambda(2*abs(jbradet(i)))*sign(1, jbradet(i))
                    L2(2*iy - 1, ix) = -L2(2*ix - 1, iy)
                    L2(2*iy, ix) = -L2(2*ix, iy)
                    if (iy .gt. iy) then
                        if ((ix .ne. iy) .and. (type_symm .eq. 1)) then
                            L2(2*(iy - size2/2) - 1, ix + size2/2) = L2(2*ix - 1, iy)
                            L2(2*(iy - size2/2), ix + size2/2) = -L2(2*ix, iy)
                            L2(2*(ix + size2/2) - 1, iy - size2/2) = -L2(2*ix - 1, iy)
                            L2(2*(ix + size2/2), iy - size2/2) = L2(2*ix, iy)
                        else if ((ix .ne. iy) .and. (type_symm .eq. 2)) then
                            L2(2*(iy - size2/2) - 1, ix + size2/2) = L2(2*ix - 1, iy)
                            L2(2*(iy - size2/2), ix + size2/2) = L2(2*ix, iy)
                            L2(2*(ix + size2/2) - 1, iy - size2/2) = -L2(2*ix - 1, iy)
                            L2(2*(ix + size2/2), iy - size2/2) = -L2(2*ix, iy)
                        end if
                    else
                        if ((ix .ne. iy) .and. (type_symm .eq. 1)) then
                            L2(2*(ix + size2/2) - 1, iy + size2/2) = L2(2*ix - 1, iy)
                            L2(2*(ix + size2/2), iy + size2/2) = -L2(2*ix, iy)
                            L2(2*(iy + size2/2) - 1, ix + size2/2) = -L2(2*ix - 1, iy)
                            L2(2*(iy + size2/2), ix + size2/2) = L2(2*ix, iy)
                        else if ((ix .ne. iy) .and. (type_symm .eq. 2)) then
                            L2(2*(ix + size2/2) - 1, iy + size2/2) = L2(2*ix - 1, iy)
                            L2(2*(ix + size2/2), iy + size2/2) = L2(2*ix, iy)
                            L2(2*(iy + size2/2) - 1, ix + size2/2) = -L2(2*ix - 1, iy)
                            L2(2*(iy + size2/2), ix + size2/2) = -L2(2*ix, iy)
                        end if
                    end if

                end if
            end if
        end if
    end do

    !Setting to zero the required parameters
    do i = 1, size2*img
        do j = 1, size2
            if (optimize(i, j) .eqv. .false.) L2(i, j) = 0.d0
        end do
    end do

    deallocate (slambda)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!It calculates the derivatives of the overlap with respect of L2 elements
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Calculate analytical derivatives of the overlap with respect to L2 elements
!> @details This subroutine calculates the analytical derivatives of the overlap Z
!> with respect to the lambda matrix elements L2. The derivatives are calculated
!> using the analytical expression derived from the overlap formula.
!>
!> For complex wave functions (img = 2), the derivatives are:
!> ∂Z/∂L2ᵢⱼ = 2[(R_L12mod * Obⱼᵢ - I_L12mod * Obⱼᵢ) * L2mod - Opⱼᵢ * L12mod²] / (L1mod * L2mod²)
!>
!> For real wave functions (img = 1), the derivatives are:
!> ∂Z/∂L2ᵢⱼ = 2[R_L12mod * Obⱼᵢ * L2mod - Opⱼᵢ * L12mod²] / (L1mod * L2mod²)
!>
!> @param[in] size2 Dimension of matrices in the new basis
!> @param[in] img Input to understand if wave function is real (1) or complex (2)
!> @param[in] Ob Matrix Ob defined in parbcs
!> @param[in] Op Matrix Op defined in parbcs
!> @param[in] L1mod Module of the previous wave function
!> @param[in] L2mod Module of the new wave function
!> @param[in] L2 Lambda matrix in the new basis
!> @param[out] dL2 Matrix of the components of the derivatives (size2*img × size2)
!> @param[in] opL2 Information on which matrix elements to optimize
!> @param[out] modder Module of the derivative vector
!> @param[in] rank MPI rank
!> @param[in] nprocu Number of MPI processes
!> @param[in] mpi_comm_world MPI communicator
!>
!> @note The subroutine calculates derivatives only for elements where opL2(i,j) = .true.
!>       The module of the derivative vector is calculated as ||dL2|| = √(Σᵢⱼ |dL2ᵢⱼ|²)
!>
!> @see calcOp, calcA, calcZ, tracemm
subroutine calcDeriv(size2, img, Ob, Op, L1mod, L2mod, L2, dL2, opL2, modder, rank, nprocu, mpi_comm_world)
    implicit none
    integer :: size2, img, i, j, rank, nprocu, mpi_comm_world
    real(8), dimension(size2*img, size2) :: Ob, Op, L2, dL2
    real(8) :: L1mod, L2mod, modder, aux1, aux2, aux3
    !L12mod2 is the square module of L12mod
    real(8) :: R_L12mod, I_L12mod, L12mod2, denomin
    real(8), allocatable :: m_aux(:, :)
    logical, dimension(size2*img, size2) :: opL2
    allocate (m_aux(size2*img, size2))

    !It calculates the value of C1 (C2 is the cc)
    call tracemm("c", "t", size2, size2, size2, size2, size2, img, Ob, L2, R_L12mod, I_L12mod)
    if (img .eq. 1) I_L12mod = 0.d0

    modder = 0
    denomin = L1mod*(L2mod**2)
    L12mod2 = I_L12mod**2 + R_L12mod**2
    dL2 = 0.d0

    if (img == 2) then
        !It calculates the value of the derivative when required
        aux1 = 2*I_L12mod*L2mod/denomin
        aux2 = 2*R_L12mod*L2mod/denomin
        aux3 = 2*L12mod2/denomin
        do i = 1, size2
            do j = 1, size2
                if (opL2(2*i - 1, j) .eqv. .true.) then
                    dL2(2*i - 1, j) = aux2*Ob(2*j - 1, i) - aux1*Ob(2*j, i) - Op(2*j - 1, i)*aux3
                    !              dL2(2*i-1, j)= 2*(( - I_L12mod * Ob(2*j,i) + R_L12mod * Ob(2*j-1, i) ) &
                    !                   * L2mod -  Op(2*j-1, i) * L12mod2 )/denomin
                    modder = modder + dL2(2*i - 1, j)**2
                end if
                if (opL2(2*i, j) .eqv. .true.) then
                    dL2(2*i, j) = aux1*Ob(2*j - 1, i) + aux2*Ob(2*j, i) - Op(2*j, i)*aux3
                    !              dL2(2*i, j)= 2*((I_L12mod * Ob(2*j-1,i) + R_L12mod * Ob(2*j, i) ) &
                    !                   * L2mod - Op(2*j, i) * L12mod2 )/denomin
                    modder = modder + dL2(2*i, j)**2
                end if
            end do
        end do
    else
        aux1 = 2*R_L12mod*L2mod/denomin
        aux2 = 2*L12mod2/denomin
        do i = 1, size2
            do j = 1, size2
                if (opL2(i, j) .eqv. .true.) then
                    dL2(i, j) = aux1*Ob(j, i) - Op(j, i)*aux2
                    modder = modder + dL2(i, j)**2
                end if
            end do
        end do
    end if

    modder = modder**0.5
    deallocate (m_aux)
end subroutine calcDeriv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculation of Op defined in parbcs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Calculate matrix Op as defined in parbcs
!> @details This subroutine calculates the matrix Op defined in the parbcs formalism.
!> The matrix Op is calculated as:
!> Op = SR * L^T * SL^T
!>
!> where SR and SL are the overlap matrices of the new basis (spin-up=left, spin-down=right),
!> and L is the lambda matrix with phase attached.
!>
!> @param[in] size Dimension of the matrices
!> @param[in] img Input to understand if wave function is real (1) or complex (2)
!> @param[in] SR Right overlap matrix (size*img × size)
!> @param[in] L Lambda matrix with phase attached (size*img × size)
!> @param[in] SL Left overlap matrix (size*img × size)
!> @param[out] Op Calculated matrix Op (size*img × size)
!> @param[in] rank MPI rank
!> @param[in] nprocu Number of MPI processes
!> @param[in] mpi_comm_world MPI communicator
!>
!> @note The subroutine uses BLAS matrix multiplication routines:
!>       - For complex case: zgemm_my for complex matrix multiplication
!>       - For real case: dgemm_my for real matrix multiplication
!>
!> @see calcOb, calcA, calcZ
subroutine calcOp(size, img, SR, L, SL, Op, rank, nprocu, mpi_comm_world)
    implicit none
    integer :: size, img, rank, nprocu, mpi_comm_world
    real(8), dimension(size*img, size) :: SR, L, SL, Op, tmp1, Lcc

    if (img == 2) then
        call zgemm_my("n", "t", size, size, size, (1.d0, 0.d0), SR, size, L, size, (0.d0, 0.d0), &
                      tmp1, size, nprocu, rank, mpi_comm_world)
        call zgemm_my("n", "t", size, size, size, (1.d0, 0.d0), tmp1, size, SL, size, (0.d0, 0.d0), &
                      Op, size, nprocu, rank, mpi_comm_world)
    end if

    if (img == 1) then
        call dgemm_my("n", "t", size, size, size, 1.d0, SR, size, L, size, 0.d0, tmp1, size, &
                      nprocu, rank, mpi_comm_world)
        call dgemm_my("n", "t", size, size, size, 1.d0, tmp1, size, SL, size, 0.d0, Op, size, &
                      nprocu, rank, mpi_comm_world)
    end if

end subroutine calcOp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculation of Ob defined in parbcs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcOb(size1, size2, img, SL12, L1, SR12, Ob, rank, nprocu, mpi_comm_world)
    implicit none
    integer :: size1, size2, img, i, rank, nprocu, mpi_comm_world
    real(8), dimension(size1*img, size2) :: SR12, SL12
    real(8), allocatable :: SLcc(:, :), tmp1(:, :)
    real(8), dimension(size2*img, size2) :: Ob
    real(8), dimension(size1*img, size1) :: L1
    allocate (tmp1(size2*img, size1), SLcc(size1*img, size2))

    !SLcc is the complex conjugate of SL12
    if (img == 2) then
        do i = 1, size1
            SLcc(2*i - 1, :) = SL12(2*i - 1, :)
            SLcc(2*i, :) = -SL12(2*i, :)
        end do
    end if

    !It distinguishes the two cases to determine the LD of the matrices
    if (img == 2) then
        call zgemm_my("c", "t", size2, size1, size1, (1.d0, 0.d0), SR12, size1, L1, size1, &
                      (0.d0, 0.d0), tmp1, size2, nprocu, rank, mpi_comm_world)
        call zgemm_my("n", "n", size2, size2, size1, (1.d0, 0.d0), tmp1, size2, SLcc, size1, &
                      (0.d0, 0.d0), Ob, size2, nprocu, rank, mpi_comm_world)
    else

        call dgemm_my("c", "t", size2, size1, size1, 1.d0, SR12, size1, L1, size1, 0.d0, tmp1, &
                      size2, nprocu, rank, mpi_comm_world)
        call dgemm_my("n", "n", size2, size2, size1, 1.d0, tmp1, size2, SL12, size1, 0.d0, Ob, &
                      size2, nprocu, rank, mpi_comm_world)
    end if

    deallocate (tmp1, SLcc)

end subroutine calcOb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculation of A defined in parbcs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Calculate the module of the wave function (A) as defined in parbcs
!> @details This subroutine calculates the module of the wave function (A), which is
!> the trace of the product of the Op matrix and the lambda matrix (or its conjugate).
!> For complex wave functions, the conjugate of L is used. The result is used to
!> normalize the overlap and derivatives in the optimization process.
!>
!> The calculation is:
!>   A = Tr[Op * L] (real)
!>   A = Tr[Op * L^*] (complex)
!>
!> @param[in] size Dimension of the matrices
!> @param[in] img Input to understand if wave function is real (1) or complex (2)
!> @param[in] L Lambda matrix (size*img × size)
!> @param[out] A Calculated module of the wave function
!> @param[in] Op Matrix Op (size*img × size)
!> @param[in] rank MPI rank
!> @param[in] nprocu Number of MPI processes
!> @param[in] mpi_comm_world MPI communicator
!>
!> @note The subroutine uses the tracemm utility to compute the trace of the matrix product.
!>
!> @see calcOp, calcOb, calcZ, tracemm
subroutine calcA(size, img, L, A, Op, rank, nprocu, mpi_comm_world)
    implicit none
    integer :: size, img, i, rank, nprocu, mpi_comm_world
    real(8), dimension(size*img, size) :: SR, L, SL, Op
    real(8), allocatable :: tmp1(:, :), Lcc(:, :)
    real(8) :: A, A1
    real(8) :: aux, aux1

    allocate (tmp1(size*img, size), Lcc(size*img, size))

    aux = 0
    !L1cc is the complex conjugate of L1cc
    if (img == 2) then
        do i = 1, size
            Lcc(2*i - 1, :) = L(2*i - 1, :)
            Lcc(2*i, :) = -L(2*i, :)
        end do
    end if

    if (img .eq. 1) then
        call tracemm("n", "n", size, size, size, size, size, img, Op, L, A, aux)
    else
        call tracemm("n", "n", size, size, size, size, size, img, Op, Lcc, A, aux)
    end if

    deallocate (tmp1, Lcc)

end subroutine calcA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!It calculates the overlap between the two wf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Calculate the overlap between two wave functions
!> @details This subroutine calculates the overlap Z between two wave functions
!> represented by lambda matrices L1 and L2. The overlap is computed as:
!>   Z = |Tr[Ob * L2]|^2 / (L1mod * L2mod)
!> where Ob is the matrix constructed from the previous basis, and L1mod, L2mod
!> are the modules of the previous and new wave functions, respectively.
!>
!> @param[in] size1 Dimension of matrices in the previous basis
!> @param[in] size2 Dimension of matrices in the new basis
!> @param[in] img Input to understand if wave function is real (1) or complex (2)
!> @param[in] L2 Lambda matrix in the new basis (size2*img × size2)
!> @param[in] Ob Matrix Ob (size2*img × size2)
!> @param[in] L1mod Module of the previous wave function
!> @param[in] L2mod Module of the new wave function
!> @param[out] Z Calculated overlap value
!> @param[in] rank MPI rank
!> @param[in] nprocu Number of MPI processes
!> @param[in] mpi_comm_world MPI communicator
!>
!> @note The subroutine uses the tracemm utility to compute the trace of the matrix product.
!>
!> @see calcOp, calcOb, calcA, tracemm
subroutine calcZ(size1, size2, img, L2, Ob, L1mod, L2mod, Z, rank, nprocu, mpi_comm_world)
    implicit none

    integer :: size1, size2, img, rank, nprocu, mpi_comm_world
    real(8) :: Z, L1mod, L2mod, R_L12mod, I_L12mod, aux, aux1, aux2
    real(8), dimension(size2*img, size2) :: L2, Ob, m_aux

    call tracemm("c", "t", size2, size2, size2, size2, size2, img, Ob, L2, R_L12mod, I_L12mod)
    !  if (rank.eq.0) write (6,*) "Mine", aux1, aux2, "Old", R_L12mod, I_L12mod

    Z = (R_L12mod**2 + I_L12mod**2)/(L1mod*L2mod)

end subroutine calcZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!It calculates the trace of a matrix matrix multiplication
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Calculate the trace of a matrix-matrix multiplication
!> @details This subroutine calculates the trace of the product of two matrices A and B,
!> supporting both real and complex cases, and different transpose/conjugate options.
!> The result is split into real and imaginary parts as needed.
!>
!> The calculation is:
!>   R = Re[Tr(op(A) * op(B))]
!>   Im = Im[Tr(op(A) * op(B))]
!> where op(A) and op(B) are optionally transposed or conjugated versions of A and B.
!>
!> @param[in] TRANSA Character flag for operation on A ('n'=none, 't'=transpose, 'c'=conjugate)
!> @param[in] TRANSB Character flag for operation on B ('n'=none, 't'=transpose, 'c'=conjugate)
!> @param[in] m, n, k Matrix dimensions
!> @param[in] LDA, LDB Leading dimensions of A and B
!> @param[in] img Input to understand if wave function is real (1) or complex (2)
!> @param[in] A First matrix (LDA*img × *)
!> @param[in] B Second matrix (LDB*img × *)
!> @param[out] R Real part of the trace
!> @param[out] Im Imaginary part of the trace
!>
!> @note The subroutine handles all combinations of transpose/conjugate for real and complex matrices.
!>
!> @see traceR, traceC
subroutine tracemm(TRANSA, TRANSB, m, n, k, LDA, LDB, img, A, B, R, Im)
    implicit none
    !This are the dimensions of the matrices (like blas)
    integer :: LDA, LDB, LDC, img, m, n, k
    integer :: i, j, ii, ij
    real(8), dimension(LDA*img, *) :: A
    real(8), dimension(LDB*img, *) :: B
    real(8), allocatable :: Aaux(:, :), Baux(:, :)
    real(8) :: aux1, aux2, R, Im
    character :: TRANSA, TRANSB

    !Aaux will be the transposed matrix of op(A), while Baux is op(B) (look the blas documentation)
    allocate (Aaux(k*img, n), Baux(k*img, n))

    !Preparing A
    !IF TRANSA eq n
    if (TRANSA .eq. "n") then
        if (img .eq. 1) then
            do j = 1, k
                do i = 1, n
                    Aaux(j, i) = A(i, j)
                end do
            end do
        else
            do j = 1, k
                do i = 1, n
                    Aaux(2*j - 1, i) = A(2*i - 1, j)
                    Aaux(2*j, i) = A(2*i, j)
                end do
            end do
        end if
        !IF TRANSA eq t
    else if ((TRANSA .eq. "t") .or. ((TRANSA .eq. "c") .and. (img .eq. 1))) then
        Aaux(:, 1:n) = A(:, 1:n)
        !IF Trans eq h
    else if ((TRANSA .eq. "c") .and. (img .eq. 2)) then
        do j = 1, k
            do i = 1, n
                Aaux(2*i - 1, j) = A(2*i - 1, j)
                Aaux(2*i, j) = -A(2*i, j)
            end do
        end do
    end if

    !Preparing B
    !IF TRANSA eq n
    if (TRANSB .eq. "n") then
        Baux(:, 1:n) = B(:, 1:n)
        !IF TRANSA eq t
    else if ((TRANSB .eq. "t") .or. ((TRANSB .eq. "c") .and. (img .eq. 1))) then
        if (img .eq. 1) then
            do j = 1, k
                do i = 1, n
                    Baux(j, i) = B(i, j)
                end do
            end do
        else
            do j = 1, k
                do i = 1, n
                    Baux(2*j - 1, i) = B(2*i - 1, j)
                    Baux(2*j, i) = B(2*i, j)
                end do
            end do
        end if
        !IF Trans eq h
    else if ((TRANSB .eq. "c") .and. (img .eq. 2)) then
        do j = 1, k
            do i = 1, n
                Baux(2*j - 1, i) = B(2*i - 1, j)
                Baux(2*j, i) = -B(2*i, j)
            end do
        end do
    end if

    R = 0.d0
    Im = 0.d0
    !Real case
    if (img .eq. 1) then
        do i = 1, k
            aux1 = 0.d0
            do j = 1, k
                aux1 = aux1 + Aaux(j, i)*Baux(j, i)
            end do
            R = R + aux1
        end do
        !Complex case
    else
        do i = 1, k
            aux1 = 0.d0
            aux2 = 0.d0
            do j = 1, k
                aux1 = aux1 + Aaux(2*j - 1, i)*Baux(2*j - 1, i) - Aaux(2*j, i)*Baux(2*j, i)
                aux2 = aux2 + Aaux(2*j - 1, i)*Baux(2*j, i) + Aaux(2*j, i)*Baux(2*j - 1, i)
            end do
            R = R + aux1
            Im = Im + aux2
        end do
    end if
    deallocate (Aaux, Baux)

end subroutine tracemm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!It calculates the trace of a real matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Calculate the trace of a real square matrix
!> @details This subroutine calculates the trace of a real square matrix M of size (size, size).
!>
!> @param[in] size Dimension of the matrix
!> @param[in] M Real matrix (size × size)
!> @param[out] Tr Calculated trace value
!>
!> @see traceC, tracemm
subroutine traceR(size, M, Tr)
    implicit none
    integer :: size, i
    real(8), dimension(size, size) :: M
    real(8) :: Tr

    Tr = 0.d0
    do i = 1, size
        Tr = Tr + M(i, i)
    end do
end subroutine traceR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!It calculates the trace of a complex matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Calculate the trace of a complex square matrix
!> @details This subroutine calculates the trace of a complex square matrix M of size (2*size, size).
!> The real and imaginary parts of the trace are returned separately.
!>
!> @param[in] size Dimension of the matrix
!> @param[in] M Complex matrix (2*size × size), real and imaginary parts interleaved
!> @param[out] Tr Real part of the trace
!> @param[out] Im_Tr Imaginary part of the trace
!>
!> @see traceR, tracemm
subroutine traceC(size, M, Tr, Im_Tr)
    implicit none
    integer :: size, i
    real(8), dimension(2*size, size) :: M
    real(8) :: Tr, Im_Tr
    Tr = 0.d0
    Im_Tr = 0.d0
    do i = 1, size
        Tr = Tr + M(2*i - 1, i)
        Im_Tr = Im_Tr + M(2*i, i)
    end do
end subroutine traceC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!It writes the final wf in finalwf.dat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Write the final lambda matrix to file
!> @details This subroutine writes the final lambda matrix L2 to the file 'finalwf.dat'.
!> For complex wave functions, both real and imaginary parts are written.
!>
!> @param[in] size2 Dimension of the matrix
!> @param[in] img Input to understand if wave function is real (1) or complex (2)
!> @param[in] L2 Lambda matrix (size2*img × size2)
!>
!> @note The output format is:
!>   - For real: i, j, L2(i, j)
!>   - For complex: i, j, Re(L2), Im(L2)
!>
!> @see max_ovlp
subroutine finalwf(size2, img, L2)
    implicit none
    integer :: size2, img, i, j
    real(8), dimension(size2*img, size2) :: L2

    open (unit=1, file="finalwf.dat", status="replace")
    if (img == 2) then
        do i = 1, size2
            do j = 1, size2
                write (1, *) i, j, L2(2*i - 1, j), L2(2*i, j)
            end do
        end do
    else
        do i = 1, size2
            do j = 1, size2
                write (1, *) i, j, L2(i, j)
            end do
        end do
    end if
    close (1)
end subroutine finalwf

