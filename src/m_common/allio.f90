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

! =============================================
! FILE: allio.f90 - TurboRVB Input/Output and Control Module
! =============================================
!
! PURPOSE:
! This module serves as the central hub for all input/output operations,
! parameter management, and control variables in TurboRVB quantum Monte Carlo
! calculations. It defines the complete interface between user input and
! the computational engine.
!
! FILE STRUCTURE:
! ===============
!
! 1. MODULE DECLARATION AND IMPORTS (Lines 15-30)
!    - Module declaration and use statements for dependencies
!    - External module imports (constants, cell, Ewald, types, etc.)
!
! 2. VARIABLE DECLARATIONS (Lines 31-1446)
!    - Basic control and configuration variables
!    - System parameters (dimensions, electron counts, etc.)
!    - Optimization and convergence parameters
!    - Physical constants and thresholds
!    - Array declarations for matrices and vectors
!    - Workspace and memory management variables
!    - MPI and parallel computing variables
!    - File I/O and scratch management variables
!
! 3. NAMELIST DEFINITIONS (Lines 1447-1940)
!    - /simulation/     : Main simulation control parameters
!    - /pseudo/         : Pseudopotential settings
!    - /readio/         : I/O and file handling options
!    - /vmc/            : Variational Monte Carlo parameters
!    - /dmclrdmc/       : Diffusion Monte Carlo and LR-DMC settings
!    - /optimization/   : Wave function optimization parameters
!    - /parameters/     : Parameter control flags
!    - /fitpar/         : Parameter fitting settings
!    - /dynamic/        : Molecular dynamics parameters
!    - /unused/         : Legacy/unused parameters
!    - /molecul/        : Molecular system settings
!    - /link/           : Link atom parameters
!
! 4. SUBROUTINES AND FUNCTIONS (Lines 1959-2561)
!    - scontract_genj()     : Matrix contraction for Jastrow factors
!    - scontract_mat_jas()  : Jastrow matrix contraction operations
!    - scontract_mat_det()  : Determinant matrix contraction operations
!    - update_kgrid()       : K-point grid update for periodic systems
!    - norm_metric()        : Vector norm calculation using metric tensor
!
! 5. STANDALONE SUBROUTINES (Lines 2587-2614)
!    - prep_map()           : Lattice vector mapping for periodic systems
!
! KEY FEATURES:
! =============
! - Comprehensive parameter management for all TurboRVB calculations
! - Support for both open and periodic boundary conditions
! - Advanced optimization algorithms (SR, LM, ADAM)
! - Multi-level parallel computing support (MPI + OpenMP)
! - Flexible I/O system with scratch file management
! - Extensive debugging and monitoring capabilities
!
! USAGE:
! ======
! This module is automatically included in all TurboRVB executables.
! Users interact with it through namelist input files that define
! calculation parameters. The module handles parameter validation,
! default value assignment, and communication between different
! parts of the code.
!
! DEPENDENCIES:
! =============
! - constants: Physical constants and mathematical parameters
! - cell: Unit cell and periodic boundary condition handling
! - Ewald: Ewald summation for long-range interactions
! - types: Custom data type definitions
! - kpoints_mod: K-point sampling and Brillouin zone integration
! - extpot: External potential and QM/MM interface
! - van_der_waals: Van der Waals interaction handling
! - link_atoms: Link atom functionality for QM/MM
! - sub_comm: Sub-communicator management for parallel computing
! - mpiio: MPI I/O operations
! - dielectric: Dielectric continuum models
!
! =============================================

module allio
    use constants
    use cell
    use Ewald
    use types
    use kpoints_mod
    use io_m, only: lchlen
    ! by E. Coccia (22/11/10)
    use extpot, only: ext_pot, link_atom, mm_restr, write_rwalk
    ! by E. Coccia (23/12/10)
    use van_der_waals, only: vdw
    ! by E. Coccia (9/5/11)
    use link_atoms
    ! for sub communicator by Y. Luo (17/10/14)
    use sub_comm
    ! for mpiio by Y. Luo (24/2/15)
    use mpiio, only: file_obj
    use dielectric
    implicit none

    ! Basic control and configuration variables
    ! Main process identifier
integer :: i_main                                ! Main process identifier
integer :: nw                                    ! Number of walkers for quantum Monte Carlo
integer :: max_sparse_choice                     ! Maximum number of sparse matrix choices
integer :: iseed                                ! Random number seed for Monte Carlo
integer :: ngenc                                ! Number of generations for optimization
integer :: ngg, ngn, ngs                        ! Generation counters for different optimization stages
integer :: d                                    ! Dimension of the system
integer :: nmax_ion                             ! Maximum number of ions
integer :: ngen                                 ! Total number of generations
integer :: nscra                                ! Number of scratch files
integer :: Lz                                   ! Z-direction length
integer :: iout                                 ! Output control flag
integer :: indvic                               ! Vicinity index
integer :: nbra                                 ! Number of branches
integer :: ng                                   ! Number of groups
integer :: iopt                                 ! Optimization flag
integer :: nmatr                                ! Number of matrices
integer :: nshell_det                           ! Number of determinant shells
integer :: iscramax                             ! Maximum scratch size
integer :: iscrapip                             ! Scratch size for pip
integer :: iscraipsip                           ! Scratch size for ipsip
integer :: iscramax_c                           ! Maximum scratch size for contracted
integer :: testderiv                            ! Test derivative flag
integer :: nwr                                  ! Number of walkers for reading
integer :: iread                                ! Read flag
integer :: ireadr                               ! Read flag for restart
integer :: nwm                                  ! Number of walkers for memory
integer :: indt                                 ! Time step index
integer :: indtupt                              ! Update time step index
integer :: nbrar                                ! Number of branches for reading
integer :: npm                                  ! Number of parameters for memory
integer :: np                                   ! Number of parameters
integer :: npr                                  ! Number of parameters for reading
integer :: npbra                                ! Number of parameters per branch
integer :: npmn                                 ! Number of parameters for memory new
integer :: npbrar                               ! Number of parameters per branch for reading
integer :: np3                                  ! Number of parameters cubed
integer :: nmp                                  ! Number of matrix parameters
integer :: nmpr                                 ! Number of matrix parameters for reading
integer :: irank                                ! Rank index
integer :: nfatr                                ! Number of factors for reading
integer :: indcor                               ! Correlation index
integer :: indt4                                ! Time step index 4
integer :: indt4j                               ! Time step index 4 for Jastrow
integer :: npp                                  ! Number of parameters plus
integer :: nfat                                 ! Number of factors
integer :: itry                                 ! Try counter
integer :: ncore                                ! Number of cores
integer :: nbram                                ! Number of branches for memory
integer :: dimjas                               ! Dimension of Jastrow
integer :: iesupmax_c                           ! Maximum number of up electrons for contracted
integer :: nmat                                 ! Number of matrices
integer :: npsov                                ! Number of pseudo overlap
integer :: npsovr                               ! Number of pseudo overlap for reading
integer :: np3r                                 ! Number of parameters cubed for reading
integer :: ncg                                  ! Number of conjugate gradient steps
integer :: ncgdim                               ! Dimension of conjugate gradient
integer :: ncgread                              ! Number of conjugate gradient steps for reading
integer :: inext                                ! Next index
integer :: ibinit                               ! Initial branch index
integer :: iend                                 ! End index
integer :: iendr                                ! End index for reading
integer :: irec                                 ! Record index
integer :: writescratch                         ! Write scratch flag
integer :: countscra                            ! Scratch counter
integer :: nelm                                 ! Total number of electrons
integer :: nelup                                ! Number of up electrons
integer :: neldo                                ! Number of down electrons
integer :: nelup_mat                            ! Number of up electrons in matrix
integer :: nel                                  ! Number of electrons
integer :: nel_mat                              ! Number of electrons in matrix
integer :: indspin                              ! Spin index
integer :: ifz                                  ! Zero flag
integer :: jmax                                 ! Maximum J value
integer :: nwnel                                ! Number of walkers per electron
integer :: nwdw                                 ! Number of walkers for down
integer :: nwrep                                ! Number of walkers for rep
integer :: info                                 ! Information flag
integer :: nwfree                               ! Number of free walkers
integer :: nrep                                 ! Number of repetitions
integer :: nws                                  ! Number of walkers for scratch
integer :: iesfree                              ! Free electron flag
integer :: iesinv                               ! Inverse flag
integer :: ieskin                               ! Kinetic energy flag
integer :: iesking                              ! Kinetic energy gradient flag
integer :: nwkin                                ! Number of walkers for kinetic
integer :: iessw                                ! Switch flag
integer :: nwsw                                 ! Number of walkers for switch
integer :: indc                                 ! Index for contracted
integer :: indc0                                ! Initial index for contracted
integer :: iesup                                ! Up electron flag
integer :: nwup                                 ! Number of walkers for up
integer :: iesdr                                ! Derivative flag for reading
integer :: iesdrr                               ! Derivative flag for restart
integer :: iesfreer                             ! Free electron flag for reading
integer :: iesswr                               ! Switch flag for reading
integer :: ieskinr                              ! Kinetic energy flag for reading
integer :: ieskinrp                             ! Kinetic energy gradient flag for reading
integer :: iesupr                               ! Up electron flag for reading
integer :: nmats                                ! Number of matrices for scratch
integer :: iseedr                               ! Random seed for reading
integer :: npf                                  ! Number of parameters for force
integer :: iesupind                             ! Up electron index
integer :: iesmind                              ! Mind index
integer :: typedyncell                          ! Cell dynamics type
integer :: ieskint                              ! Kinetic energy type
integer :: nwking                               ! Number of walkers for kinetic gradient
integer :: iesswr_eagp                          ! Switch flag for EAGP
integer :: nwnp                                 ! Number of walkers for new parameters
integer :: np3m                                 ! Number of parameters cubed for memory
integer :: icount                               ! Counter
integer :: ireadmin                             ! Read admin flag
integer :: ireadminr                            ! Read admin flag for restart
integer :: iessw0                               ! Initial switch flag
integer :: iesup_read                           ! Up electron read flag
integer :: nelnw                                ! Number of electrons for new walkers
integer :: nelnw0                               ! Initial number of electrons for new walkers
integer :: jold                                 ! Old J value
integer :: niesd                                ! Number of iterations for energy and derivatives
integer :: ngentry                              ! Number of generation entries
integer :: numcost                              ! Number of cost functions
integer :: Ltab                                 ! Table length
integer :: Lztab                                ! Z-table length
integer :: Ltabb                                ! Table length for branches
integer :: indold                               ! Old index
integer :: indnn                                ! New index
integer :: js                                   ! Jastrow index
integer :: jtype2                               ! Jastrow type 2
integer :: nelorb                               ! Number of orbitals
integer :: nshell                               ! Number of shells
integer :: nshelldo                             ! Number of down shells
integer :: nrnel                                ! Number of real electrons
integer :: nelkel                               ! Number of kelvin electrons
integer :: dimee                                ! Dimension of electron-electron
integer :: dimei                                ! Dimension of electron-ion
integer :: iesconv                              ! Convergence flag
integer :: nelnion                              ! Number of electrons per ion
integer :: nion                                 ! Number of ions
integer :: indksij                              ! K-space index for ij
integer :: k_ion                                ! K-space index for ions
integer :: i_ion                                ! Ion index
integer :: nelorbj                              ! Number of Jastrow orbitals
integer :: indtj_nw                             ! Time step index for new walkers
integer :: indtabj_nw                           ! Table index for new walkers
integer :: indtabbj_nw                          ! Branch table index for new walkers
integer :: indkj_nw                             ! K-space index for new walkers
integer :: indksj_nw                            ! K-space shell index for new walkers
integer :: indksij_nw                           ! K-space ij index for new walkers
integer :: indwwj_nw                            ! Wave function index for new walkers
integer :: indwupj_nw                           ! Up wave function index for new walkers
integer :: indwdoj_nw                           ! Down wave function index for new walkers
integer :: indaupj_nw                           ! Up atomic index for new walkers
integer :: j_nw                                 ! J value for new walkers
integer :: j_nws                                ! J value for new walkers scratch
integer :: irej                                 ! Rejection index
integer :: nshellj                              ! Number of Jastrow shells
integer :: indwwjj                              ! Wave function index for Jastrow
integer :: indwwjj_nw                           ! Wave function index for Jastrow new walkers
integer :: nel2wtj                              ! Number of electrons for Jastrow
integer :: indbar                               ! Bar index
integer :: indjbar                              ! Jastrow bar index
integer :: indbar_nw                            ! Bar index for new walkers
integer :: indjbar_nw                           ! Jastrow bar index for new walkers
integer :: nelorbp                              ! Number of primitive orbitals
integer :: nelorbpp                             ! Number of primitive orbitals plus
integer :: indjbarsz                            ! Jastrow bar size index
integer :: indjbarsz_nw                         ! Jastrow bar size index for new walkers
integer :: npar3body                            ! Number of 3-body parameters
integer :: iesswind                             ! Switch index
integer :: nelcol                               ! Number of columns
integer :: nelcol_c                             ! Number of columns for contracted
integer :: inddsw                               ! Down switch index
integer :: indn                                 ! New index
integer :: iflag                                ! Flag
integer :: nelcolh                              ! Number of columns for half
integer :: irankdet                             ! Rank determinant
integer :: pippo                                ! Pippo parameter
integer :: cpippo                               ! Contracted pippo
integer :: indpippo                             ! Pippo index
integer :: dimpippo                             ! Pippo dimension
integer :: ccpippo                              ! Contracted contracted pippo
integer :: nnozeromax                           ! Maximum number of non-zero elements
integer :: skipforce                            ! Skip force flag
integer :: pippoc                               ! Pippo contracted
integer :: ndone                                ! Number of completed steps
integer :: lbin                                 ! Bin length
integer :: nbinr                                ! Number of bins for reading
integer :: nbin                                 ! Number of bins
integer :: nbinread                             ! Number of bins for reading
integer :: iesfreesz                            ! Free size flag
integer :: ixr                                  ! X index for reading
integer :: iyr                                  ! Y index for reading
integer :: nbinmax                              ! Maximum number of bins
integer :: stepcg                               ! Conjugate gradient step
integer :: nbindim                              ! Bin dimension
integer :: ieskindim                            ! Kinetic energy dimension
integer :: npdim                                ! Parameter dimension
integer :: ndimpdim                             ! Parameter dimension plus
integer :: lwork                                ! Work length
integer :: indfix                               ! Fixed index
integer :: repf                                 ! Repetition factor
integer :: kinf                                 ! Kinetic energy factor
integer :: idyn                                 ! Dynamics flag
integer :: iskipdyn                             ! Skip dynamics flag
integer :: iskipdynr                            ! Skip dynamics flag for reading
integer :: idynu                                ! Dynamics update flag
integer :: icost                                ! Cost flag
integer :: indmax                               ! Maximum index
integer :: indj                                 ! Jastrow index
integer :: kp_ion                               ! K-space ion index
integer :: indfree                              ! Free index
integer :: npar                                 ! Number of parameters
integer :: initpar                              ! Initial parameter
integer :: nparsw                               ! Number of parameters for switch
integer :: initparsw                            ! Initial parameter for switch
integer :: indsw                                ! Switch index
integer :: indinv                               ! Inverse index
integer :: endinv                               ! End inverse index
integer :: nparinv                              ! Number of parameters for inverse
integer :: initparinv                           ! Initial parameter for inverse
integer :: iflagk                               ! K-space flag
integer :: maxall                               ! Maximum all
integer :: molecular                            ! Molecular flag
integer :: molecularj                           ! Molecular Jastrow flag
integer :: nshell_max                           ! Maximum number of shells
integer :: nshellj_max                          ! Maximum number of Jastrow shells
integer :: maxiesup                             ! Maximum up electron flag
integer :: ieskingdim                           ! Kinetic energy gradient dimension
integer :: maxnpar                              ! Maximum number of parameters
integer :: developer                            ! Developer flag
integer :: molopt                               ! Molecular optimization flag
integer :: nparshellmax                         ! Maximum number of parameters per shell
integer :: optbra                               ! Branch optimization flag
integer :: unreliable                           ! Unreliable flag
integer :: iflagerrall                          ! Error flag for all
integer :: ieskinold                            ! Old kinetic energy flag
integer :: nel2wtfn                             ! Number of electrons for wave function
integer :: indwwjfn                             ! Wave function index
integer :: nel2barfn                            ! Number of electrons for bar wave function
integer :: indbarfn                             ! Bar wave function index
integer :: indbarfn_nw                          ! Bar wave function index for new walkers
integer :: indwwjfn_nw                          ! Wave function index for new walkers
integer :: powermin                             ! Minimum power
integer :: npower                               ! Number of powers
integer :: powerminsz                           ! Minimum power size
integer :: npowersz                             ! Number of power sizes
integer :: ieskinr_pos                          ! Positive kinetic energy flag for reading
integer :: indtjr                               ! Time step index for reading
integer :: Lztabr                               ! Z-table length for reading
integer :: indtjr_nw                            ! Time step index for new walkers
integer :: commopt_mpi                          ! MPI optimization communicator
integer :: commcolopt_mpi                       ! MPI column optimization communicator
integer :: rankopt                              ! Optimization rank
integer :: nprocopt                             ! Number of optimization processes
integer :: row_comm                             ! Row communicator
integer :: col_comm                             ! Column communicator
integer :: row_id                               ! Row ID
integer :: col_id                               ! Column ID
integer :: commcov_mpi                          ! MPI covariance communicator
integer :: nproccov                             ! Number of covariance processes
integer :: commsr_mpi                           ! MPI search communicator
integer :: nprocsr                              ! Number of search processes
integer :: ranksr                               ! Search rank
integer :: mcol                                 ! Column size
integer :: commrep_mpi                          ! MPI repetition communicator
integer :: nprocrep                             ! Number of repetition processes
integer :: rankrep                              ! Repetition rank
integer :: nrep_bead                            ! Number of repetition beads
integer :: mcol_rep                             ! Column size for repetition
integer :: rankcolrep                           ! Column repetition rank
integer :: commcolrep_mpi                       ! MPI column repetition communicator
integer :: nproccolrep                          ! Number of column repetition processes
integer :: commcolsr_mpi                        ! MPI column search communicator
integer :: nproc_diag                           ! Number of diagonalization processes
integer :: nprocu                               ! Number of update processes
integer :: ref_atom                             ! Reference atom
integer :: epscuttyper                          ! Epsilon cut type for reading
integer :: ndim_detmat                          ! Dimension of determinant matrix
integer :: delay_changeparr                     ! Delay for parameter change
integer :: iesdelay                             ! Delay flag
integer :: min_block                            ! Minimum block size
integer :: max_ortho                            ! Maximum orthogonalization
integer :: maxiter_changeparr                   ! Maximum iterations for parameter change
integer :: prep                                 ! Preparation flag
integer :: comm_col                             ! Column communicator
integer :: comm_raw                             ! Raw communicator
integer :: rankraw                              ! Raw rank
integer :: rankcol                              ! Column rank
integer :: nbra_cyrus                           ! Number of branches for Cyrus
integer :: dim_cyrus                            ! Dimension for Cyrus
integer :: nw_max                               ! Maximum number of walkers
integer :: npar_eagp                            ! Number of parameters for EAGP
integer :: dim_upwf                             ! Dimension of up wave function
integer :: dim_ratiovar                         ! Dimension of ratio variance
integer :: dim_uptabtot                         ! Dimension of up table total
integer :: nnozero_eagp                         ! Number of non-zero elements for EAGP
integer :: max_target                           ! Maximum target
integer :: max_targetsr                         ! Maximum target for search
integer :: nbra_cyrus_read                      ! Number of branches for Cyrus reading

type(mpi_sub_comm) :: sub_comm_diag             ! MPI sub-communicator for diagonalization

integer :: occ                                  ! Occupation number
integer :: nnozero                              ! Number of non-zero elements
integer :: nnozeroj                             ! Number of non-zero elements for j
integer :: occtot                               ! Total occupation
integer :: kl                                   ! k-point index
integer :: freqcheck                            ! Frequency check parameter
integer :: nmaxder                              ! Maximum derivative order
integer :: occj                                 ! Occupation number for j
integer :: occtotj                              ! Total occupation for j
integer :: indp                                 ! Index pointer
integer :: iflagnorm                            ! Normalization flag
integer :: nnozerojr                            ! Number of non-zero elements for j (real)
integer :: novar                                ! Number of variables
integer :: parcutg                              ! Parameter cutoff for g
integer :: equil_steps                          ! Equilibrium steps
integer :: firstmol                             ! First molecule index
integer :: nmolfn                               ! Number of molecular functions
integer :: yesfast                              ! Fast calculation flag
integer :: iesup_atom                           ! Atom index for iesup
integer :: dimtranspip                          ! Dimension of transpip
integer :: lastmol                              ! Last molecule index
integer :: nelorb_at                            ! Number of atomic orbitals
integer :: nel2up                               ! Number of up electrons (2)
integer :: nel2upt                              ! Number of up electrons (2) total
integer :: nel2dot                              ! Number of down electrons (2)
integer :: nel2wt                               ! Number of weighted electrons (2)
integer :: nel2bar                              ! Number of barred electrons (2)
integer :: indtot                               ! Total index
integer :: nel2jbar                             ! Number of barred electrons (2) for j
integer :: nel2jbarsz                           ! Size of barred electrons (2) for j
integer :: iesuptouched                         ! Flag for touched iesup
integer :: vjutouched                           ! Flag for touched vju
integer :: indtj                                ! Index for j
integer :: indtabj                              ! Table index for j
integer :: indtabbj                             ! Back table index for j
integer :: indkj                                ! k-index for j
integer :: indwupj                              ! Up weight index for j
integer :: indwdoj                              ! Down weight index for j
integer :: indaupj                              ! Up atomic index for j
integer :: iskip                                ! Skip flag
integer :: indwwj                               ! Weight index for j
integer :: iesupskip                            ! Skip flag for iesup
integer :: indksj                               ! k-shell index for j
integer :: nmatb                                ! Number of matrix b
integer :: nmatbr                               ! Number of matrix b (real)
integer :: nshellr                              ! Number of real shells
integer :: nshelljr                             ! Number of real shells for j
integer :: nelorbh                              ! Number of hybrid orbitals
integer :: nelorbjh                             ! Number of hybrid orbitals for j
integer :: nelorbjh2                            ! Number of hybrid orbitals for j (2)
integer :: irstart                              ! Initial restart index
integer :: nx                                   ! Number of points in x direction
integer :: ny                                   ! Number of points in y direction
integer :: nz                                   ! Number of points in z direction
integer :: nbufd                                ! Buffer dimension
integer :: nmol                                 ! Number of molecules
integer :: nmolmin                              ! Minimum number of molecules
integer :: nmolmax                              ! Maximum number of molecules
integer :: nel3                                 ! Number of electrons (3)
integer :: ndiff                                ! Number of differences
integer :: ndiffdim                             ! Dimension of differences
integer :: indndiff                             ! Index for differences
integer :: ifreqdump                            ! Frequency dump flag
integer :: iimax                                ! Maximum i index
integer :: iijmax                               ! Maximum ij index
integer :: nmolmaxw                             ! Maximum number of molecules (weighted)
integer :: typereg                              ! Type of regularization
integer :: niont                                ! Number of ions (total)
integer :: niong                                ! Number of ions (g)
integer :: epscuttype                           ! Type of epsilon cutoff
integer :: ncg_adr                              ! Number of CG addresses
integer :: kp0                                  ! Initial k-point
integer :: nbead                                ! Number of beads
integer :: neldomax                             ! Maximum number of electrons in domain

real*8 :: pot_aasunel                          ! Potential for aasunel
real*8 :: gamma                                ! Gamma parameter
real*8 :: etry                                 ! Trial energy
real*8 :: etryr                                ! Real trial energy
real*8 :: rata                                 ! Ratio parameter
real*8 :: time                                 ! Time parameter
real*8 :: wbra                                 ! Bra weight
real*8 :: lambda                               ! Lambda parameter
real*8 :: eps_umrigar                          ! Epsilon parameter for Umrigar
real*8 :: ener                                 ! Energy
real*8 :: ratio                                ! Ratio parameter
real*8 :: wbra1                                ! Bra weight 1
real*8 :: wbra2                                ! Bra weight 2
real*8 :: ratior(2)                            ! Ratio array (2 elements)
real*8 :: ratiodet                             ! Ratio determinant
real*8 :: sumdiff                              ! Sum of differences
real*8 :: timep                                ! Time parameter p
real*8 :: timepp                               ! Time parameter pp
real*8 :: timescra                             ! Scraping time
real*8 :: ratiorn(2)                           ! Ratio array n (2 elements)
real*8 :: psign                                ! Sign of psi
real*8 :: wbran                                ! Bra weight n
real*8 :: wcort                                ! Weight correction
real*8 :: timemc                               ! Monte Carlo time
real*8 :: timeopt                              ! Optimization time
real*8 :: timeinit                             ! Initialization time
real*8 :: veff                                 ! Effective potential
real*8 :: veffright                            ! Right effective potential
real*8 :: beta                                 ! Beta parameter
real*8 :: temp                                 ! Temperature
real*8 :: kappar                               ! Kappa parameter
real*8 :: cost                                 ! Cost function
real*8 :: cost1                                ! Cost function 1
real*8 :: cost2                                ! Cost function 2
real*8 :: parr                                 ! Parameter array
real*8 :: econtnew                             ! New energy contribution
real*8 :: jacold                               ! Old Jacobian
real*8 :: parcutr                              ! Parameter cutoff (real)
real*8 :: epst                                 ! Epsilon t
real*8 :: epsi                                 ! Epsilon i
real*8 :: psiav                                ! Average psi
real*8 :: psisav                               ! Saved psi
real*8 :: pdiag                                ! Diagonal parameter
real*8 :: ttry                                 ! Trial parameter
real*8 :: psirav                               ! Average psi (real)
real*8 :: countav                              ! Average count
real*8 :: countt                               ! Total count
real*8 :: avreweight                           ! Average reweight
real*8 :: tleft                                ! Left time
real*8 :: tbra                                 ! Bra time
real*8 :: tbrar                                ! Bra time (real)
real*8 :: counttot                             ! Total count
real*8 :: countreg                             ! Regular count
real*8 :: countcut                             ! Cut count
real*8 :: parcut                               ! Parameter cutoff
real*8 :: costpassed                           ! Passed cost
real*8 :: ttry1                                ! Trial parameter 1
real*8 :: ttry2                                ! Trial parameter 2
real*8 :: pdiag1                               ! Diagonal parameter 1
real*8 :: pdiag2                               ! Diagonal parameter 2
real*8 :: wcorwt                               ! Weight correction
real*8 :: epsdgel                              ! Epsilon for dgel
real*8 :: eps_dyn5                             ! Epsilon for dynamic 5
real*8 :: epstion                              ! Epsilon for ion
real*8 :: epsdgm                               ! Epsilon for dgm
real*8 :: memtot                               ! Total memory
real*8 :: wtotf(2)                             ! Total weight function (2 elements)
real*8 :: tpar                                 ! Parameter time
real*8 :: tparf                                ! Parameter time final
real*8 :: epstl                                ! Epsilon tl
real*8 :: epstlu                               ! Epsilon tlu
real*8 :: rsignr                               ! Real sign
real*8 :: parbest                              ! Best parameter
real*8 :: wtotf1                               ! Total weight function 1
real*8 :: wtotf2                               ! Total weight function 2
real*8 :: rweight                              ! Real weight
real*8 :: hopfraction                          ! Hop fraction
real*8 :: ration(2)                            ! Ratio array n (2 elements)
real*8 :: nontr                                ! Non-trivial parameter
real*8 :: tstep                                ! Time step
real*8 :: bcost                                ! Base cost
real*8 :: ccost                                ! Current cost
real*8 :: minz                                 ! Minimum z
real*8 :: maxz                                 ! Maximum z
real*8 :: minzj                                ! Minimum z for j
real*8 :: maxzj                                ! Maximum z for j
real*8 :: ukwald                               ! UK Wald parameter
real*8 :: f                                    ! Function value
real*8 :: fb                                   ! Function value b
real*8 :: vpotint                              ! Potential integral
real*8 :: diffint                              ! Difference integral
real*8 :: costexp                              ! Cost exponent
real*8 :: theta_reg                            ! Theta regularization
real*8 :: epscut                               ! Epsilon cutoff
real*8 :: epstlrat                             ! Epsilon tl ratio
real*8 :: enermin                              ! Minimum energy
real*8 :: varmin                               ! Minimum variance
real*8 :: spsi                                 ! Scaled psi
real*8 :: ratioreg                             ! Ratio regularization
real*8 :: epscutu                              ! Epsilon cutoff u
real*8 :: fmax                                 ! Maximum function value
real*8 :: parcutmin                            ! Minimum parameter cutoff
real*8 :: parcutpar                            ! Parameter cutoff parameter
real*8 :: parcute                              ! Parameter cutoff e
real*8 :: nacc                                 ! Number of acceptances
real*8 :: naccpseudo                           ! Number of pseudo acceptances
real*8 :: npow                                 ! Power parameter
real*8 :: normcorr                             ! Normalization correction
real*8 :: costwn                               ! Cost weight n
real*8 :: signflip                             ! Sign flip parameter
real*8 :: acclarge                             ! Large acceptance parameter
real*8 :: weightall                            ! All weights
real*8 :: tmes                                 ! Time measurement
real*8 :: rsr                                  ! Real space radius
real*8 :: rmax                                 ! Maximum radius
real*8 :: rmaxj                                ! Maximum radius for j
real*8 :: rmaxinv                              ! Inverse of maximum radius
real*8 :: epscutur                             ! Epsilon cutoff ur
real*8 :: epstlur                              ! Epsilon tl ur
real*8 :: time_main                            ! Main time
real*8 :: time_ratiovar                        ! Time for ratio variance
real*8 :: time_uptabtot                        ! Total update table time
real*8 :: timewf                               ! Wave function time
real*8 :: timepip                              ! Pipeline time
real*8 :: ax                                   ! x-axis parameter
real*8 :: ay                                   ! y-axis parameter
real*8 :: az                                   ! z-axis parameter
real*8 :: zmax                                 ! Maximum z
real*8 :: tstepfn                              ! Time step function
real*8 :: tion                                 ! Ion time
real*8 :: tcell                                ! Cell time
real*8 :: weight_loc                           ! Local weight
real*8 :: power                                ! Power parameter
real*8 :: Klrdmc                               ! KLR DMC parameter
real*8 :: epscutdmc                            ! Epsilon cutoff for DMC
real*8 :: epstldmc                             ! Epsilon tl for DMC
real*8 :: time_meas                            ! Measurement time
real*8 :: time_branch                          ! Branching time
real*8 :: cutreg                               ! Cut regularization
real*8 :: cutweight                            ! Cut weight
real*8 :: powerwarp                            ! Power warp
real*8 :: epsrem_contr                         ! Epsilon remaining contribution
real*8 :: tolcg                                ! Tolerance for CG
real*8 :: smoothcut                            ! Smooth cut
real*8 :: pressclass                           ! Pressure class
real*8 :: dcellclass(3)                        ! Cell class differences (3 elements)
real*8 :: scalepulay                           ! Pulay scale
real*8 :: shift                                ! Shift parameter
real*8 :: alat2v                               ! Lattice parameter 2v
real*8 :: maxtime                              ! Maximum time
real*8 :: smearing                             ! Smearing parameter
real*8 :: scale_mass                           ! Mass scale
real*8 :: scale_one_body                       ! One-body scale
real*8 :: scaleeloc                            ! Local energy scale
real*8 :: minjonetwobody                       ! Minimum Jastrow two-body
real*8 :: timings(11)                          ! Timing array (11 elements)
real*8 :: timingsb(11)                         ! Timing array b (11 elements)
real*8 :: parr_min                             ! Minimum parameter array
real*8 :: parr_max                             ! Maximum parameter array
real*8 :: epsvar                               ! Epsilon variance
real*8 :: tave_cyrus                           ! Average time for Cyrus
real*8 :: tcount_cyrus                         ! Count time for Cyrus
real*8 :: zmin                                 ! Minimum z
real*8 :: l0_kousuke                           ! L0 parameter for Kousuke
real*8 :: core_pseudo                          ! Core pseudo
real*8 :: count_zerowf                         ! Count of zero wave functions  
real*8 :: count_allwf                          ! Count of all wave functions

! Logical variables (added by Andrea Tirelli)
logical :: tpar_increased                      ! Flag for increased tpar
logical :: stop_increasing_tpar                ! Flag to stop increasing tpar
logical :: yes_cutweight                       ! Flag for cut weight
logical :: yes_scemama                         ! Flag for scemama
logical :: yes_scemama_open                    ! Flag for open scemama
logical :: yes_sparse                          ! Flag for sparse
logical :: yes_sparse_choose                   ! Flag for sparse choose
logical :: yes_dgelscut                        ! Flag for dgels cut

! Integer variables (added by Andrea Tirelli)
integer :: counter_very_unstable_tpar          ! Counter for very unstable tpar
integer :: tpar_unstble_stop                   ! Stop for unstable tpar
integer :: inc_tpar_frequency                  ! Frequency of tpar increase
integer :: counter_unstable_tpar               ! Counter for unstable tpar
integer :: counter_unstable_energy             ! Counter for unstable energy
integer :: counter_unstable_err                ! Counter for unstable error
integer :: len_tpar_stable_list                ! Length of stable tpar list

! Real variables (added by Andrea Tirelli)
real(8) :: min_running_ave                     ! Minimum running average
real(8) :: min_running_std                     ! Minimum running standard deviation
real(8) :: min_running_ave_energy              ! Minimum running average energy
real(8) :: min_running_std_energy              ! Minimum running standard deviation energy
real(8) :: tpar_max                           ! Maximum tpar
real(8) :: divide_tpar                        ! Division factor for tpar
real(8) :: multiply_tpar                      ! Multiplication factor for tpar
real(8) :: n_sigmas_tpar                      ! Number of sigmas for tpar

! Integer*8 variables
integer*8 :: handle                           ! CUDA handle for device operations

! Integer*4 variables  
integer*4 :: ldworkspace                      ! Dimension of double precision workspace
integer*4 :: lzworkspace                      ! Dimension of complex workspace
integer*4 :: dev_Info(1)                      ! Device information array (1 element)

! Real*8 allocatable arrays
real*8, allocatable, dimension(:) :: dev_dgetrf_workspace      ! Workspace for double precision GETRF
real*8, allocatable, dimension(:, :) :: dev_dgetri_workspace  ! Workspace for double precision GETRI

! Complex*16 allocatable arrays
complex*16, allocatable, dimension(:) :: dev_zgetrf_workspace  ! Workspace for complex GETRF
complex*16, allocatable, dimension(:, :) :: dev_zgetri_workspace ! Workspace for complex GETRI

! Real*8 variables
real*8 :: deriv1                              ! First derivative
real*8 :: deriv2                              ! Second derivative
real*8 :: deriv3                              ! Third derivative
real*8 :: pulay1                              ! First Pulay term
real*8 :: pulay2                              ! Second Pulay term
real*8 :: pulay3                              ! Third Pulay term
real*8 :: scalecell(3)                        ! Cell scaling factors (3 elements)
real*8 :: epsder                              ! Epsilon for derivatives
real*8 :: scale_grad                          ! Gradient scaling
real*8 :: norm_corr                           ! Normalization correction

! Integer variables
integer :: icdiff                             ! Index for difference calculation
integer :: itest                              ! Test index
integer :: itestr                             ! Real test index
integer :: jn                                 ! Index j
integer :: nrest                              ! Number of restarts
integer :: iesm                               ! Index for ESM
integer :: isfix                              ! Fixed index
integer :: nwfix                              ! Number of fixed weights
integer :: itestrr                            ! Real test index rr
integer :: itestr3                            ! Real test index 3
integer :: itestr4                            ! Real test index 4
integer :: iesd                               ! Index for ESD
integer :: iese                               ! Index for ESE
integer :: ieser                              ! Real index for ESE
integer :: npsamax                            ! Maximum PSA
integer :: nwm2                               ! Number of weights minus 2
integer :: nprest                             ! Number of prestarts
integer :: nindt                              ! Number of independent terms
integer :: nwdim                              ! Dimension of weights
integer :: nintpsa                            ! Number of PSA intervals
integer :: ix                                 ! x-index
integer :: iy                                 ! y-index
integer :: iboot                              ! Bootstrap index
integer :: xj                                 ! Index j
integer :: itestrfn                           ! Real test function index
integer :: npsar                              ! Real PSA
integer :: indberry                           ! Berry index
integer :: true_wagner                        ! True Wagner index
integer :: nelsquare                          ! Number of electrons squared

! Real*8 variables
real*8 :: wback                               ! Back weight
real*8 :: costw                               ! Cost weight
real*8 :: enerc                               ! Current energy
real*8 :: jacnew                              ! New Jacobian
real*8 :: gradold                             ! Old gradient
real*8 :: gradbarold                          ! Old barred gradient
real*8 :: friction                            ! Friction parameter
real*8 :: delta0                              ! Initial delta
real*8 :: delta0q                             ! Quadratic delta
real*8 :: delta0k                             ! Kinetic delta
real*8 :: dt                                  ! Time step
real*8 :: scalecov                            ! Covariance scale
real*8 :: maxdev_dyn                          ! Maximum deviation for dynamics
real*8 :: pressfixed                          ! Fixed pressure
real*8 :: versoralat(3, 12)                   ! Versor lattice (3x12)
real*8 :: epsbas                              ! Epsilon basis
real*8 :: lepsbas                             ! Length epsilon basis
real*8 :: weight_moroni                       ! Moroni weight

! Integer variables
integer :: ntry                               ! Number of tries
integer :: nwinv                              ! Number of inverse weights
integer :: nweight                            ! Number of weights
integer :: nweightr                           ! Number of real weights
integer :: nmore_force                        ! Number of additional forces

! Real*8 variables
real*8 :: alat                                ! Lattice parameter
real*8 :: alat2                               ! Squared lattice parameter
real*8 :: sigma_new                           ! New sigma
real*8 :: sigma_true(2)                       ! True sigma (2 elements)
real*8 :: costa                               ! Cost parameter
real*8 :: ener_true(2)                        ! True energy (2 elements)
real*8 :: plat(3)                             ! Lattice vectors (3)
real*8 :: dstep(3)                            ! Step size (3)

! Integer variables
integer :: ndim                               ! Dimension
integer :: ndimp                              ! Parameter dimension
integer :: ndimj                              ! J dimension
integer :: ndimjp                             ! JP dimension
integer :: ndims                              ! S dimension
integer :: ndimsp                             ! SP dimension
integer :: ndimiesup                          ! IESUP dimension
integer :: movedion                           ! Moved ion flag

! Allocatable real(8) arrays (2D)
real(8), dimension(:, :), allocatable :: rion                  ! Ion positions
real(8), dimension(:, :), allocatable :: rion_fast            ! Fast ion positions
real(8), dimension(:, :), allocatable :: rionsav              ! Saved ion positions
real(8), dimension(:, :), allocatable :: efenergy             ! Effective energy
real(8), dimension(:, :), allocatable :: gradpsi              ! Psi gradient
real(8), dimension(:, :), allocatable :: gradpsibar           ! Barred psi gradient
real(8), dimension(:, :), allocatable :: efpress              ! Effective pressure
real(8), dimension(:, :), allocatable :: gradpsiold           ! Old psi gradient
real(8), dimension(:, :), allocatable :: gradpsibarold        ! Old barred psi gradient
real(8), dimension(:, :), allocatable :: angle                ! Angle
real(8), dimension(:, :), allocatable :: fk                   ! Force constant
real(8), dimension(:, :), allocatable :: reduce               ! Reduction factor
real(8), dimension(:, :), allocatable :: reducel              ! Local reduction
real(8), dimension(:, :), allocatable :: rcne                 ! RCNE
real(8), dimension(:, :), allocatable :: warpmat              ! Warp matrix
real(8), dimension(:, :), allocatable :: projmat_c            ! Projection matrix C

! Allocatable real(8) arrays (4D)
real(8), dimension(:, :, :, :), allocatable :: jastrowall_ee    ! Jastrow all electron-electron
real(8), dimension(:, :, :, :), allocatable :: queue_cyrus      ! Cyrus queue

! Allocatable real(8) arrays (3D)
real(8), dimension(:, :, :), allocatable :: rmunew              ! New RMU
real(8), dimension(:, :, :), allocatable :: sov                 ! SOV
real(8), dimension(:, :, :), allocatable :: ef                  ! EF
real(8), dimension(:, :, :), allocatable :: efp                 ! EFP
real(8), dimension(:, :, :), allocatable :: ivic                ! IVIC
real(8), dimension(:, :, :), allocatable :: jastrowall_ei       ! Jastrow all electron-ion
real(8), dimension(:, :, :), allocatable :: rmunewb             ! New RMU B
real(8), dimension(:, :, :), allocatable :: agp                 ! AGP

! Allocatable integer arrays (3D)
integer, dimension(:, :, :), allocatable :: npip                ! NPIP

! Allocatable real(8) arrays (2D)
real(8), dimension(:, :), allocatable :: kel                    ! Kinetic energy
real(8), dimension(:, :), allocatable :: keln                   ! New kinetic energy
real(8), dimension(:, :), allocatable :: dists_kel              ! Kinetic energy distances
real(8), dimension(:, :), allocatable :: velion                 ! Ion velocity
real(8), dimension(:, :), allocatable :: diffkin                ! Kinetic difference
real(8), dimension(:, :), allocatable :: rknew                  ! New RK
real(8), dimension(:, :), allocatable :: rknewb                 ! New RK B

! Allocatable real(8) arrays (2D)
real(8), dimension(:, :), allocatable :: iond_cart              ! Ion distance Cartesian
real(8), dimension(:, :), allocatable :: vpotreg                ! Regularized potential
real(8), dimension(:, :), allocatable :: kdyn                   ! Kinetic dynamics
real(8), dimension(:, :), allocatable :: vpotsav_ee             ! Saved electron-electron potential

! winv vector for complex DFT
! variables required for complex total energy

! Complex variables
complex(8) :: enerc_c                                           ! Complex current energy
complex(8) :: ener_true_c                                       ! Complex true energy
complex(8) :: ener_c                                            ! Complex energy
complex(8) :: enermin_c                                         ! Complex minimum energy

! Allocatable real(8) arrays (1D)
real(8), dimension(:), allocatable :: duprold                   ! Old duplicate parameter
real(8), dimension(:), allocatable :: vjurold                   ! Old VJU
real(8), dimension(:), allocatable :: eig                       ! Eigenvalues
real(8), dimension(:), allocatable :: vj                        ! VJ parameter
real(8), dimension(:), allocatable :: dup                       ! Duplicate parameter
real(8), dimension(:), allocatable :: dsw                       ! DSW parameter
real(8), dimension(:), allocatable :: zetar                     ! Zeta real
real(8), dimension(:), allocatable :: zetar_fast                ! Fast zeta real
real(8), dimension(:), allocatable :: zetaq                     ! Zeta q
real(8), dimension(:), allocatable :: zetamin                   ! Minimum zeta
real(8), dimension(:), allocatable :: distmin                   ! Minimum distance
real(8), dimension(:), allocatable :: dek                       ! Kinetic energy difference
real(8), dimension(:), allocatable :: dekg                      ! Kinetic energy gradient difference
real(8), dimension(:), allocatable :: vju                       ! VJU parameter
real(8), dimension(:), allocatable :: ddw                       ! DDW parameter
real(8), dimension(:), allocatable :: ddwsz                     ! DDW size
real(8), dimension(:), allocatable :: vjur                      ! Real VJU
real(8), dimension(:), allocatable :: dupr                      ! Real duplicate parameter
real(8), dimension(:), allocatable :: winv                      ! Inverse wave function
real(8), dimension(:), allocatable :: ainv                      ! Inverse A
real(8), dimension(:), allocatable :: vpot                      ! Potential
real(8), dimension(:), allocatable :: enertrue                  ! True energy
real(8), dimension(:), allocatable :: diffuse                   ! Diffuse parameter
real(8), dimension(:), allocatable :: ainvup                    ! Up inverse A
real(8), dimension(:), allocatable :: ainvdo                    ! Down inverse A
real(8), dimension(:), allocatable :: winvup                    ! Up inverse wave function
real(8), dimension(:), allocatable :: winvdo                    ! Down inverse wave function
real(8), dimension(:), allocatable :: scale                     ! Scale parameter
real(8), dimension(:), allocatable :: scalpar                   ! Parameter scale
real(8), dimension(:), allocatable :: scalej                    ! Jastrow scale
real(8), dimension(:), allocatable :: scalejsz                  ! Jastrow scale size
real(8), dimension(:), allocatable :: winvj                     ! Jastrow inverse wave function
real(8), dimension(:), allocatable :: cnorm                     ! Normalization constant
real(8), dimension(:), allocatable :: psinew                    ! New psi
real(8), dimension(:), allocatable :: tmu                       ! TMU parameter
real(8), dimension(:), allocatable :: winvbar                   ! Barred inverse wave function
real(8), dimension(:), allocatable :: winvjbar                  ! Jastrow barred inverse wave function
real(8), dimension(:), allocatable :: detmat                    ! Determinant matrix
real(8), dimension(:), allocatable :: jasmat                    ! Jastrow matrix
real(8), dimension(:), allocatable :: jasmatsz                  ! Jastrow matrix size
real(8), dimension(:), allocatable :: cnorm_nw                  ! New walker normalization constant
real(8), dimension(:), allocatable :: psibar                    ! Barred psi
real(8), dimension(:), allocatable :: ainvupb                   ! Barred up inverse A
real(8), dimension(:), allocatable :: err                       ! Error
real(8), dimension(:), allocatable :: force                     ! Force
real(8), dimension(:), allocatable :: dist                      ! Distance
real(8), dimension(:), allocatable :: iond                      ! Ion distance
real(8), dimension(:), allocatable :: rmu                       ! RMU parameter
real(8), dimension(:), allocatable :: r                         ! Radius
real(8), dimension(:), allocatable :: winvs                     ! Scaled inverse wave function
real(8), dimension(:), allocatable :: winvsj                    ! Scaled Jastrow inverse wave function
real(8), dimension(:), allocatable :: ainvs                     ! Scaled inverse A
real(8), dimension(:), allocatable :: dists                     ! Scaled distances
real(8), dimension(:), allocatable :: wint                      ! Weight integral
real(8), dimension(:), allocatable :: alphavar                  ! Alpha variance
real(8), dimension(:), allocatable :: berry_exp                 ! Berry exponential
real(8), dimension(:), allocatable :: wcorw                     ! Weight correction
real(8), dimension(:), allocatable :: wintw                     ! Weight integral w
real(8), dimension(:), allocatable :: alphab                    ! Alpha b
real(8), dimension(:), allocatable :: diagfn                    ! Diagonal function
real(8), dimension(:), allocatable :: psip                      ! Psi prime
real(8), dimension(:), allocatable :: psip_reweight             ! Reweighted psi prime
real(8), dimension(:), allocatable :: econf                     ! Configuration energy
real(8), dimension(:), allocatable :: econfh                    ! Half configuration energy
real(8), dimension(:), allocatable :: wconfn                    ! New walker configuration weight
real(8), dimension(:), allocatable :: econfion                  ! Ion configuration energy
real(8), dimension(:), allocatable :: factorsr                  ! Factors r
real(8), dimension(:), allocatable :: vcut                      ! Potential cutoff
real(8), dimension(:), allocatable :: etot                      ! Total energy
real(8), dimension(:), allocatable :: wsto                      ! Weight storage
real(8), dimension(:), allocatable :: zeta                      ! Zeta
real(8), dimension(:), allocatable :: tabpip                    ! Table pip
real(8), dimension(:), allocatable :: table                     ! Table
real(8), dimension(:), allocatable :: tabler                    ! Table r
real(8), dimension(:), allocatable :: diag                      ! Diagonal
real(8), dimension(:), allocatable :: tcore                     ! Core time
real(8), dimension(:), allocatable :: tcost                     ! Cost time
real(8), dimension(:), allocatable :: gradtot                   ! Total gradient
real(8), dimension(:), allocatable :: gradtotbar                ! Barred total gradient
real(8), dimension(:), allocatable :: costz                     ! Cost z
real(8), dimension(:), allocatable :: costz3                    ! Cost z3
real(8), dimension(:), allocatable :: rcarto                    ! Old Cartesian radius
real(8), dimension(:), allocatable :: rcart                     ! Cartesian radius
real(8), dimension(:), allocatable :: dx_old                    ! Old dx
real(8), dimension(:), allocatable :: dx_new                    ! New dx
real(8), dimension(:), allocatable :: jasnew_ei                 ! New Jastrow electron-ion
real(8), dimension(:), allocatable :: fkav                      ! Average force constant
real(8), dimension(:), allocatable :: winvjbarsz                ! Jastrow barred inverse wave function size
real(8), dimension(:), allocatable :: okav                      ! Average ok
real(8), dimension(:), allocatable :: skdiag                    ! Scaled kinetic diagonal
real(8), dimension(:), allocatable :: cov                       ! Covariance
real(8), dimension(:), allocatable :: cov_old                   ! Old covariance
real(8), dimension(:), allocatable :: jasnew_ee                 ! New Jastrow electron-electron
real(8), dimension(:), allocatable :: atom_number               ! Atom number
real(8), dimension(:), allocatable :: atom_number_fast          ! Fast atom number
real(8), dimension(:), allocatable :: winvbarn                  ! Barred inverse wave function n
real(8), dimension(:), allocatable :: enerint                   ! Energy integral
real(8), dimension(:), allocatable :: enerintw                  ! Weighted energy integral
real(8), dimension(:), allocatable :: wconfsav                  ! Saved configuration weight
real(8), dimension(:), allocatable :: wtot                      ! Total weight
real(8), dimension(:), allocatable :: winvbarfn                 ! Barred inverse wave function fn
real(8), dimension(:), allocatable :: winvfn                    ! Inverse wave function fn
real(8), dimension(:), allocatable :: t_cyrus                   ! Cyrus time
real(8), dimension(:), allocatable :: first_moment              ! First moment
real(8), dimension(:), allocatable :: second_moment             ! Second moment

! Real variables
real :: maxoutput                                               ! Maximum output value

! Allocatable real(4) arrays (1D)
real(4), dimension(:), allocatable :: econfw                    ! Configuration energy with weight
real(4), dimension(:), allocatable :: wconfw                    ! Configuration weight

! Allocatable real(8) arrays (1D)
real(8), dimension(:), allocatable :: psilnw                    ! New wave function logarithm
real(8), dimension(:), allocatable :: bufscra                   ! Buffer scratch
real(8), dimension(:), allocatable :: v_adr                     ! V address (sensitive data)

! Allocatable integer arrays (1D)
integer, dimension(:), allocatable :: ipsip                     ! Psi prime index
integer, dimension(:), allocatable :: mult                      ! Multiplicity
integer, dimension(:), allocatable :: kion                      ! Ion index
integer, dimension(:), allocatable :: nparam                    ! Number of parameters
integer, dimension(:), allocatable :: multj                     ! Jastrow multiplicity
integer, dimension(:), allocatable :: nparamj                   ! Jastrow number of parameters
integer, dimension(:), allocatable :: kionj                     ! Jastrow ion index
integer, dimension(:), allocatable :: ioccj                     ! Jastrow occupation
integer, dimension(:), allocatable :: ioptorbj                  ! Jastrow orbital option
integer, dimension(:), allocatable :: nozero                    ! Non-zero elements
integer, dimension(:), allocatable :: nozerodet                 ! Non-zero determinant elements
integer, dimension(:), allocatable :: nozeroj                   ! Jastrow non-zero elements
integer, dimension(:), allocatable :: ioptorb                   ! Orbital option
integer, dimension(:), allocatable :: ioccup                    ! Occupation
integer, dimension(:), allocatable :: ioccdo                    ! Down occupation
integer, dimension(:), allocatable :: naccm                     ! Number of accepted moves
integer, dimension(:), allocatable :: jbra                      ! Bra index
integer, dimension(:), allocatable :: jbraw                     ! Weighted bra index
integer, dimension(:), allocatable :: jbraj                     ! Jastrow bra index
integer, dimension(:), allocatable :: jbrajsz                   ! Jastrow bra size
integer, dimension(:), allocatable :: jbradet                   ! Determinant bra index
integer, dimension(:), allocatable :: jbradetn                  ! New determinant bra index
integer, dimension(:), allocatable :: jbrajn                    ! New Jastrow bra index
integer, dimension(:), allocatable :: jbraiesup                 ! Up electron-ion bra index
integer, dimension(:), allocatable :: jbraiesm                  ! Down electron-ion bra index
integer, dimension(:), allocatable :: whereiesup                ! Up electron-ion position
integer, dimension(:), allocatable :: whereiesm                 ! Down electron-ion position
integer, dimension(:), allocatable :: vjutouch                  ! VJU touch
integer, dimension(:), allocatable :: itouch                    ! Touch index
integer, dimension(:), allocatable :: icore                     ! Core index
integer, dimension(:), allocatable :: jbraiesup_sav             ! Saved up electron-ion bra index
integer, dimension(:), allocatable :: jbraiesm_sav              ! Saved down electron-ion bra index
integer, dimension(:), allocatable :: nozerojder                ! Jastrow derivative non-zero elements
integer, dimension(:), allocatable :: ipip_adr                  ! Pip address
integer, dimension(:), allocatable :: first_cyrus               ! First Cyrus

! Allocatable logical arrays (1D)
logical, dimension(:), allocatable :: sjbradet                  ! Determinant bra flag
logical, dimension(:), allocatable :: slaterorb_read            ! Slater orbital read flag

! Allocatable ion component arrays (1D)
type(ion_comp), dimension(:), allocatable :: ion_table          ! Ion table

! Allocatable integer arrays (1D)
integer, allocatable :: indpar_tab(:)                           ! Parameter index table
integer, allocatable :: indorb_tab(:)                           ! Orbital index table
integer, allocatable :: indshell_tab(:)                         ! Shell index table
integer, allocatable :: indparj_tab(:)                          ! Jastrow parameter index table
integer, allocatable :: indorbj_tab(:)                          ! Jastrow orbital index table
integer, allocatable :: indshellj_tab(:)                        ! Jastrow shell index table
integer, allocatable :: adr_nion(:)                             ! Ion address
integer, allocatable :: ind_nion(:)                             ! Ion index
integer, allocatable :: adrj_nion(:)                            ! Jastrow ion address
integer, allocatable :: indj_nion(:)                            ! Jastrow ion index
integer, allocatable :: addr_occ(:)                             ! Occupation address
integer, allocatable :: type_atom(:)                            ! Atom type

! Allocatable integer arrays (2D)
integer, allocatable :: pointvj(:, :)                           ! VJ point

! Integer variables
integer :: npsa                                                 ! Number of pseudo atoms
integer :: lmax                                                 ! Maximum angular momentum
integer :: istart                                               ! Start index
integer :: indteff                                              ! Effective index
integer :: lzeff                                                ! Effective z

! Character variables
character(3) :: pseudoname                                      ! Pseudo name (3 characters)
character(60) :: pseudofile                                     ! Pseudo file name (60 characters)

! Allocatable integer arrays (2D)
integer, dimension(:, :), allocatable :: nparpshell             ! Number of parameters per shell
integer, dimension(:, :), allocatable :: jpseudo                ! Pseudo index
integer, dimension(:, :), allocatable :: indtm                  ! Index matrix

! Allocatable integer arrays (1D)
integer, dimension(:), allocatable :: kindion                   ! Ion kind
integer, dimension(:), allocatable :: pshell                    ! Pseudo shell

! Allocatable real(8) arrays (2D)
real(8), dimension(:, :), allocatable :: parshell               ! Shell parameters
real(8), dimension(:, :), allocatable :: legendre               ! Legendre polynomials
real(8), dimension(:, :), allocatable :: versor                 ! Versor
real(8), dimension(:, :), allocatable :: prefactor              ! Prefactor
real(8), dimension(:, :), allocatable :: enert                  ! Energy matrix
real(8), dimension(:, :), allocatable :: agpn                   ! AGP matrix

! Allocatable real(8) arrays (1D)
real(8), dimension(:), allocatable :: rcutoff                   ! Cutoff radius
real(8), dimension(:), allocatable :: wpseudo                   ! Pseudo weight
real(8), dimension(:), allocatable :: pseudolocal               ! Local pseudo
real(8), dimension(:), allocatable :: wintpseudo                ! Pseudo weight integral

! Integer variables
integer :: npseudopar                                           ! Number of pseudo parameters
integer :: npseudoparn                                          ! New number of pseudo parameters
integer :: ion                                                  ! Ion index
integer :: nintpseudo                                           ! Number of pseudo integrals

! Real(8) variables
real(8) :: coeff                                                ! Coefficient
real(8) :: beta_learning                                        ! Learning rate

! Variables for complex/quantum algorithms
real(8), dimension(:), allocatable :: psiln                     ! Wave function logarithm
real(8), dimension(:), allocatable :: psidetln                  ! Determinant logarithm
real(8), dimension(:), allocatable :: kdyn_eig                  ! Dynamic eigenvalues
real(8), dimension(:), allocatable :: psisn                     ! Wave function sign
complex(8), dimension(:), allocatable :: psip_c                 ! Complex wave function
complex(8) :: ratio_c                                           ! Complex ratio
complex(8) :: ratior_c                                          ! Complex real ratio
complex(8) :: ratiodet_c                                        ! Complex determinant ratio

! Real(8) variables
real(8) :: psioverpsi_new                                       ! New wave function ratio
real(8) :: psioverpsi_old                                       ! Old wave function ratio

! Real(8) variables
real(8) :: cutoff_p                                             ! Cutoff for bump orbitals (abs(Lbox) = 2), default = 9.d0

! Variables used only in scalapack (__SCALAPACK) but allocated always for compatibility
! (irrelevant memory overhead)

! Integer variables for BLACS (Basic Linear Algebra Communication Subprograms)
integer :: me_blacs = 0                                         ! BLACS processor index starting from 0
integer :: np_blacs = 1                                         ! BLACS number of processor
integer :: world_cntx = 0                                       ! BLACS context of all processor
integer :: ortho_cntx = 0                                       ! BLACS context for ortho_comm
integer :: me_ortho(2) = 0                                      ! Coordinates of the processors
integer :: me_ortho1 = 0                                        ! Task id for the ortho group
integer :: np_ortho(2) = 1                                      ! Size of the processor grid used in ortho
integer :: np_ortho1 = 1                                        ! Size of the ortho group
integer :: ortho_comm = 0                                       ! Communicator used for fast and memory saving ortho
integer :: ortho_comm_id = 0                                    ! ID of the ortho_comm
integer :: leg_ortho = 1                                        ! The distance in the father communicator of two neighbour processors in ortho_comm

! General integer variables
integer :: rank                                                 ! Process rank
integer :: nproc                                                ! Number of processes
integer :: nprocr                                               ! Number of processes in row
integer :: ist                                                  ! Start index
integer :: ien                                                  ! End index
integer :: id1                                                  ! ID 1
integer :: istm                                                 ! Start index m
integer :: in1                                                  ! Input 1
integer :: ierr                                                 ! Error code
integer :: skipreshuff                                          ! Skip reshuffle flag

! Character variables
character(14) :: chara                                          ! Character string (14 characters)
character(14) :: charaq                                         ! Character string q (14 characters)
character(60 + lchlen) :: wherescratch                         ! Scratch directory path (60 + lchlen characters)
character(lchlen) :: errmsg                                     ! Error message (lchlen characters)

! Integer variables
integer :: nshell_c                                             ! Number of shells
integer :: occ_c                                                ! Occupation
integer :: iesup_c                                              ! Up electron-ion
integer :: nnozero_c                                            ! Number of non-zero elements
integer :: nelorb_c                                             ! Number of orbitals
integer :: contraction                                          ! Contraction
integer :: ikshift                                              ! K shift
integer :: nshellj_c                                            ! Jastrow number of shells
integer :: occj_c                                               ! Jastrow occupation
integer :: npar3body_c                                          ! Number of 3-body parameters
integer :: nnozeroj_c                                           ! Jastrow number of non-zero elements
integer :: nelorbj_c                                            ! Jastrow number of orbitals
integer :: contractionj                                         ! Jastrow contraction
integer :: maxparam                                             ! Maximum number of parameters
integer :: maxioccmult                                          ! Maximum occupation multiplicity
integer :: nelorbmax                                            ! Maximum number of orbitals
integer :: indocc                                               ! Occupation index
integer :: maxshell                                             ! Maximum number of shells
integer :: maxparamj                                            ! Maximum number of Jastrow parameters
integer :: maxshellj                                            ! Maximum number of Jastrow shells
integer :: nelorbmaxj                                           ! Maximum number of Jastrow orbitals
integer :: ll                                                   ! Angular momentum l
integer :: mm                                                   ! Angular momentum m
integer :: iesupr_2                                             ! Up electron-ion 2
integer :: iesupr_c                                             ! Up electron-ion c
integer :: npar3bodyr_c                                         ! Real 3-body parameters
integer :: occ_tmp                                              ! Temporary occupation
integer :: npar3body_fill                                       ! Fill 3-body parameters
integer :: npar3body_2                                          ! 3-body parameters 2
integer :: nmolmat                                              ! Molecular matrix
integer :: nmolmatw                                             ! Weighted molecular matrix

! Allocatable integer arrays (1D)
integer, dimension(:), allocatable :: mult_c                    ! Multiplicity
integer, dimension(:), allocatable :: nparam_c                  ! Number of parameters
integer, dimension(:), allocatable :: ioptorb_c                 ! Orbital option
integer, dimension(:), allocatable :: kion_c                    ! Ion index
integer, dimension(:), allocatable :: nozero_c                  ! Non-zero elements
integer, dimension(:), allocatable :: ioccup_c                  ! Occupation
integer, dimension(:), allocatable :: occshell                  ! Shell occupation
integer, dimension(:), allocatable :: iesuptrans                ! Up electron-ion transform
integer, dimension(:), allocatable :: multranspip               ! Multi transform pip
integer, dimension(:), allocatable :: iesuptransb               ! Up electron-ion transform b
integer, dimension(:), allocatable :: multj_c                   ! Jastrow multiplicity
integer, dimension(:), allocatable :: nparamj_c                 ! Jastrow number of parameters
integer, dimension(:), allocatable :: ioptorbj_c                ! Jastrow orbital option
integer, dimension(:), allocatable :: kionj_c                   ! Jastrow ion index
integer, dimension(:), allocatable :: nozeroj_c                 ! Jastrow non-zero elements
integer, dimension(:), allocatable :: ioccj_c                   ! Jastrow occupation
integer, dimension(:), allocatable :: occshellj                 ! Jastrow shell occupation
integer, dimension(:), allocatable :: iesuptransj               ! Jastrow up electron-ion transform
integer, dimension(:), allocatable :: multranspipj              ! Jastrow multi transform pip
integer, dimension(:), allocatable :: iesuptransbj              ! Jastrow up electron-ion transform b
integer, dimension(:), allocatable :: kiontot                   ! Total ion index
integer, dimension(:), allocatable :: kiontotj                  ! Total Jastrow ion index
integer, dimension(:), allocatable :: ioptorbja                 ! Jastrow orbital option a
integer, dimension(:), allocatable :: typeorb                   ! Orbital type

! Allocatable array_int type arrays (1D)
type(array_int), allocatable :: transpip(:)                     ! Transform pip
type(array_int), allocatable :: transpip_sav(:)                ! Saved transform pip
type(array_int), allocatable :: transpipj(:)                    ! Jastrow transform pip
type(array_int), allocatable :: transpipj_sav(:)               ! Saved Jastrow transform pip

! Allocatable nkgrid type arrays (1D)
type(nkgrid), allocatable :: kgrid(:)                          ! K grid
type(nkgrid), allocatable :: kgrid_atom(:)                     ! Atom K grid

! Allocatable integer arrays (2D)
integer, dimension(:, :), allocatable :: mu_touch               ! Mu touch matrix
integer, dimension(:, :), allocatable :: muj_touch              ! Jastrow mu touch matrix
integer, dimension(:, :), allocatable :: adrlambda              ! Lambda address

! Allocatable real(8) arrays (1D)
real(8), dimension(:), allocatable :: dup_c                     ! Duplicate parameter
real(8), dimension(:), allocatable :: scale_c                   ! Scale parameter
real(8), dimension(:), allocatable :: detmat_c                  ! Determinant matrix
real(8), dimension(:), allocatable :: projm                     ! Projection matrix
real(8), dimension(:), allocatable :: detmat_proj               ! Projected determinant matrix
real(8), dimension(:), allocatable :: vju_c                     ! VJU parameter
real(8), dimension(:), allocatable :: scalej_c                  ! Jastrow scale
real(8), dimension(:), allocatable :: scalejsz_c                ! Jastrow scale size
real(8), dimension(:), allocatable :: jasmat_c                  ! Jastrow matrix
real(8), dimension(:), allocatable :: jasmatsz_c                ! Jastrow matrix size
real(8), dimension(:), allocatable :: rpar                      ! Real parameters

! Allocatable real(8) arrays (2D)
real(8), dimension(:, :), allocatable :: jas_invariant          ! Jastrow invariant
real(8), dimension(:, :), allocatable :: eagp_pfaff             ! EAGP Pfaffian
real(8), dimension(:, :), allocatable :: eagp_pfaffb            ! EAGP Pfaffian b
real(8), dimension(:, :), allocatable :: mu_tmp                 ! Temporary mu matrix
real(8), dimension(:, :), allocatable :: mu_c                   ! Mu matrix
real(8), dimension(:, :), allocatable :: muj_tmp                ! Temporary Jastrow mu matrix
real(8), dimension(:, :), allocatable :: muj_c                  ! Jastrow mu matrix
real(8), dimension(:, :), allocatable :: mat_adr                ! Matrix address

! Allocatable logical arrays (1D)
logical, dimension(:), allocatable :: orbcost                   ! Orbital cost flag
logical, dimension(:), allocatable :: orbcostn                  ! New orbital cost flag
logical, dimension(:), allocatable :: orbcostl                  ! Local orbital cost flag
logical, dimension(:), allocatable :: yescut                    ! Cutoff flag
logical, dimension(:), allocatable :: orbps                     ! Orbital psi flag
logical, dimension(:), allocatable :: singdet                   ! Single determinant flag
logical, dimension(:), allocatable :: allowed_par               ! Allowed parameters flag

! Allocatable logical arrays (2D)
logical, dimension(:, :), allocatable :: allowcontr             ! Allowed contraction flag
logical, dimension(:, :), allocatable :: allowcontrj            ! Allowed Jastrow contraction flag

! Logical variables for general settings
logical :: pseudologic                                          ! Pseudo logic flag
logical :: iessz                                                ! Electron spin flag
logical :: iesgros                                              ! Gross flag
logical :: iescost                                              ! Cost flag
logical :: fncont                                               ! Function contraction flag
logical :: iesbra                                               ! Bra flag
logical :: yeszj                                                ! ZJ flag
logical :: molyes                                               ! Molecular flag
logical :: moljyes                                              ! Jastrow molecular flag
logical :: iesrandoma                                           ! Random a flag
logical :: iesrandoml                                           ! Random l flag
logical :: pseudorandom                                         ! Pseudo random flag
logical :: iespbc                                               ! Periodic boundary conditions flag
logical :: fnloc                                                ! Local function flag
logical :: yesivic                                              ! IVIC flag
logical :: yesnleft                                             ! N left flag
logical :: rejweight                                            ! Rejection weight flag
logical :: orthoyes                                             ! Orthogonal flag
logical :: yeszagp                                              ! ZAGP flag
logical :: yesdetmat                                            ! Determinant matrix flag
logical :: symmagp                                              ! Symmetric AGP flag
logical :: yesdft                                               ! DFT flag
logical :: yesfmubar                                            ! FMU bar flag
logical :: membig                                               ! Memory big flag
logical :: membigcpu                                            ! CPU memory big flag
logical :: membigr                                              ! Real memory big flag
logical :: onebodysz                                            ! One body size flag
logical :: yespress                                             ! Pressure flag
logical :: iescostd                                             ! Cost d flag
logical :: twobodyoff                                           ! Two body off flag
logical :: printoverlap                                         ! Print overlap flag
logical :: yesdetmatc                                           ! Determinant matrix c flag
logical :: oldscra                                              ! Old scratch flag
logical :: symiesup                                             ! Symmetric up flag
logical :: iesdtwobodyoff                                       ! Two body off flag
logical :: iesdonebodyoff                                       ! One body off flag
logical :: warp                                                 ! Warp flag
logical :: yesbump                                              ! Bump flag
logical :: yeslbox                                              ! L box flag
logical :: yespulay                                             ! Pulay flag
logical :: stopdyn                                              ! Stop dynamics flag
logical :: allfit                                               ! All fit flag
logical :: yescutjas                                            ! Jastrow cutoff flag
logical :: yescutdet                                            ! Determinant cutoff flag
logical :: defparcutg                                           ! Default parameter cutoff flag
logical :: detc_proj                                            ! Determinant projection flag
logical :: yesread10                                            ! Read 10 flag
logical :: stepcg_recount                                       ! Step CG recount flag
logical :: write_cov                                            ! Write covariance flag
logical :: gramyes                                              ! Gram flag
logical :: fixpar                                               ! Fix parameter flag
logical :: symmetrize_agp                                       ! Symmetrize AGP flag
logical :: yesprimitive                                         ! Primitive flag

! Logical variables for quantum algorithm
logical :: yesQuantum                                           ! Quantum flag
logical :: yesavsr                                              ! Average SR flag
logical :: yesavcov                                             ! Average covariance flag
logical :: yesavopt                                             ! Average optimization flag
logical :: yeswrite10                                           ! Write 10 flag
logical :: yesperiodize                                         ! Periodize flag
logical :: yesturboq                                            ! Turbo Q flag
logical :: yes_complex                                          ! Complex flag

! Up electron flag
! yesupel=.true.    --> means dealing with spin up electrons, now one can choose a different phase for up/down spin
logical :: yesupel
! Second flag
logical :: yessecond

! Crystal flag
! yes_crystal=.true. --> use Crystal basis set defined as: \phi_k(r)=\sum_R \phi(r-R_a-R)*exp(ikR). Use the input variable
logical :: yes_crystal

! Jastrow crystal flag
!   yes_crystalj=.true. --> use Crystal basis  also for the Jastrow
!                        epsbas to set the cutoff on the sum over the direct lattice vectors.
logical :: yes_crystalj

! Test AAD flag
! test_aad=.true. --> used only when you are using the program testadc.x for testing AAD derivatives, otherwise set to .false.
logical :: test_aad

logical :: eqcellab                ! Force cell parameters (a,b) to be equal when using typedyncell>0 for cell relaxation
logical :: eqcellac                ! Force cell parameters (a,c) to be equal when using typedyncell>0 for cell relaxation  
logical :: eqcellbc                ! Force cell parameters (b,c) to be equal when using typedyncell>0 for cell relaxation
logical :: forcecomplex            ! Force complex flag
logical :: oldscaling              ! Old scaling flag
logical :: ldynsecond              ! Dynamic second flag
logical :: add_onebody2det         ! Add one body to determinant flag
logical :: yes_hermite             ! Hermite flag
logical :: allowed_averagek        ! Allowed average k flag
logical :: yes_correct             ! Correct flag
logical :: yes_real                ! Real flag
logical :: srcomplex               ! Complex SR flag
logical :: killcut                 ! Kill cutoff flag
logical :: change_epscut           ! Change epsilon cutoff flag
logical :: change_tstep            ! Change time step flag
logical :: better_dmc              ! Better DMC flag
logical :: yesalfe                 ! Alpha F flag
logical :: safelrdmc               ! Safe LR DMC flag
logical :: yesrootc                ! Root c flag
logical :: addrognoso              ! Add Rognoso flag
logical :: changelambda            ! Change lambda flag
logical :: cleanrognoso            ! Clean Rognoso flag
logical :: fixa                    ! Fix a flag
logical :: fixb                    ! Fix b flag
logical :: fixc                    ! Fix c flag
logical :: forcesymm               ! Force symmetry flag
logical :: signalnoise             ! Signal noise flag
logical :: real_contracted         ! Real contracted flag
logical :: gauge_fixing            ! Gauge fixing flag
logical :: yesmin_read             ! Minimum read flag
logical :: noopt_onebody           ! No optimization one body flag
logical :: real_agp                ! Real AGP flag
logical :: softcusp                ! Soft cusp flag
logical :: scalermax               ! Scale r max flag
logical :: yeswritebead            ! Write bead flag
logical :: yes_hessc               ! Hessian c flag
logical :: no_sjbra                ! No SJ bra flag
logical :: manyfort10              ! Many fort10 flag
logical :: shift_origin            ! Shift origin flag
logical :: shiftx                  ! Shift x flag
logical :: shifty                  ! Shift y flag
logical :: shiftz                  ! Shift z flag
logical :: double_mesh             ! Double mesh flag
logical :: change_parr             ! Change parameter flag
logical :: default_epsdgel         ! Default epsilon dgel flag
logical :: read_molecul            ! Read molecule flag
logical :: hybyes                  ! Hybrid flag
logical :: pfaffup                 ! Pfaffian up flag
logical :: k6gen                   ! K6 generation flag
logical :: noblocking              ! No blocking flag
logical :: add_diff                ! Add difference flag
logical :: lrdmc_der               ! LR DMC derivative flag
logical :: lrdmc_nonodes           ! LR DMC no nodes flag
logical :: nosingledet             ! No single determinant flag
logical :: enforce_detailb         ! Enforce detail b flag
logical :: nowrite12               ! No write 12 flag
logical :: yes_fastbranch          ! Fast branch flag
logical :: flush_write             ! Flush write flag
logical :: yes_adams               ! Adams flag
logical :: only_molecular          ! Only molecular flag
logical :: add_offmol              ! Add off molecular flag
logical :: novec_loop1             ! No vector loop1 flag

! Control IO behavior (Y. Luo, 15/2/15)
character(len=80) :: disk_io                ! Disk IO mode (80 characters)
integer :: io_level                         ! IO level

! MPI-IO file objects
type(file_obj) :: kelcont                  ! Kelvin container
type(file_obj) :: quantcont                ! Quantum container
type(file_obj) :: details_SP               ! Single precision details
type(file_obj) :: details_DP               ! Double precision details

! Parameters for automatic gutta cavat lapidem (to adjust tpar)
! Added by K. Nakano

! Parameters read from input file
logical :: change_tpar                      ! Change tpar flag
logical :: use_stable_tpar                  ! Use stable tpar flag

! Internal variables (not read from input file)
integer :: inc_counter_tpar                 ! Increment counter for tpar
integer :: dec_counter_tpar                 ! Decrement counter for tpar
integer :: tpar_buffer_len                  ! Tpar buffer length
integer :: len_shorter_buffer               ! Length of shorter buffer
integer :: times_tpar_decreased             ! Number of times tpar decreased
real(8) :: cut_sigma                       ! Cut sigma value

! Allocatable logical arrays (1D)
logical, dimension(:), allocatable :: tpar_buffer_filled  ! Tpar buffer filled flag

! Allocatable real(8) arrays (1D)
real(8), dimension(:), allocatable :: energy_list         ! Energy list
real(8), dimension(:), allocatable :: error_energy_list   ! Error energy list
real(8), dimension(:), allocatable :: tpar_stable_list    ! Stable tpar list

! ************ EWALD SUMS *******************!
! Integer variables
integer :: kmax2                            ! Maximum k squared
integer :: xi                               ! Xi parameter
integer :: n_body_on                        ! N-body on flag
integer :: yesmin                           ! Minimum flag
integer :: yesminr                          ! Real minimum flag

! Allocatable real(8) arrays (2D)
real(8), dimension(:, :), allocatable :: rmusin  ! Sine of rmu
real(8), dimension(:, :), allocatable :: rmucos  ! Cosine of rmu

! Real(8) variables
real(8) :: selfsum                          ! Self sum
real(8) :: LBox                             ! Box length
real(8) :: LBoxj                            ! Jastrow box length
real(8) :: Linv                             ! Inverse length
real(8) :: rs                               ! Radius
real(8) :: ris(5)                           ! Radius inverse (5 elements)
real(8) :: derEVp                           ! Derivative of energy with respect to volume (Pulay)
real(8) :: errEVp                           ! Error of derivative of energy with respect to volume (Pulay)
real(8) :: derEVnopulay                     ! Derivative of energy with respect to volume (no Pulay)
real(8) :: errEVnopulay                     ! Error of derivative of energy with respect to volume (no Pulay)

! Allocatable real(8) arrays (1D)
real(8), dimension(:), allocatable :: p_pulay  ! Pulay parameter
real(8), dimension(:), allocatable :: derEV    ! Derivative of energy with respect to volume

! Integer variables
integer :: add_pulay                        ! Add Pulay flag

! ************ CAFFAREL FORCES *************
! Note: These are comments describing the purpose of variables
! orbderiv = orbital derivatives respect to nuclei
! derpot   = derivatives of potential energy
! dercaf = laplacian of caffarel Q

! Default values
! rank = 0
! nproc = 1
! nw = 1

! =============================================
! /simulation/ namelist parameters
! =============================================
! Optimization Control
! ------------------
! itestr4        : Test type for optimization (4: VMC optimization)
! iopt           : Optimization type
! ngen           : Number of generations for optimization
! kappar         : Kappa parameter for optimization
! nbra           : Number of bra states

! Memory and Resource Management
! ----------------------------
! membig         : Memory size for big arrays
! membigcpu      : CPU memory size for big arrays
! nscra          : Number of scratch files
! maxtime        : Maximum computation time
! disk_io        : Disk I/O control flag

! Parallel Computing
! ----------------
! nproc_diag     : Number of processors for diagonalization
! ip_reshuff     : Reshuffling control parameter

! Walkers and States
! ----------------
! nw             : Number of walkers
! kSq            : Wave vector squared
! iseedr         : Random seed for reproducibility

! Checkpointing and Monitoring
! -------------------------
! freqcheck      : Frequency of checkpointing

! Developer and Performance Options
! -----------------------------
! developer      : Developer mode flag
! yesfast        : Fast mode flag
! novec_loop1    : Flag to disable vectorization in loop1

! Band Structure and Mesh
! ---------------------
! compute_bands  : Flag to compute energy bands
! double_mesh    : Flag for double mesh calculation

! Target and Block Control
! ---------------------
! min_block      : Minimum block size
! max_target     : Maximum target value
! max_targetsr   : Maximum target value for SR

! Dielectric Properties
! ------------------
! dielectric_ratio: Dielectric ratio parameter
! dielectric_length: Dielectric length parameter
! case_diel      : Dielectric case selector

! Neighbor and Sparse Matrix
! ------------------------
! neigh          : Neighbor list control
! yes_sparse     : Flag to enable sparse matrix operations
! yes_sparse_choose: Flag to enable sparse matrix selection
! max_sparse_choice: Maximum number of sparse matrix choices

namelist /simulation/ itestr4, iopt, ngen, nscra, nbra, iseedr, nw, kSq, &
    & kappar, freqcheck, membig, membigcpu, developer, yesfast, maxtime, &
    & nproc_diag, disk_io, ip_reshuff, compute_bands, double_mesh, &
    & min_block, max_target, max_targetsr, dielectric_ratio, &
    & dielectric_length, case_diel, neigh, novec_loop1, yes_sparse, &
    & yes_sparse_choose, max_sparse_choice

! =============================================
! /pseudo/ namelist parameters
! =============================================
! Pseudopotential Control
! ---------------------
! nintpsa        : Number of integration points for pseudopotential
! npsamax        : Maximum number of pseudopotential parameters
! pseudorandom   : Logical flag for pseudopotential randomization

namelist /pseudo/ nintpsa, npsamax, pseudorandom

! =============================================
! /readio/ namelist parameters
! =============================================
! Parallel Processing
! ----------------
! ncore          : Number of cores to use for parallel computation
! np3            : Number of processors for 3D parallelization
! np             : Number of processors for general parallelization

! I/O Control
! ----------
! iread          : Read control flag (0: normal, 1: restart, etc.)
! writescratch   : Logical flag to control writing to scratch files
! wherescratch   : Directory path for scratch files
! nowrite12      : Logical flag to control writing to file 12
! flush_write    : Logical flag to control buffer flushing behavior

! Data Management
! -------------
! unreliable     : Logical flag for handling unreliable data
! ifreqdump      : Frequency of data dumps during calculation

namelist /readio/ ncore, np3, np, iread, writescratch, wherescratch, unreliable, ifreqdump, nowrite12, flush_write

! =============================================
! /vmc/ namelist parameters
! =============================================
! Time Step and Dynamics
! --------------------
! tstep          : Time step for VMC (Variational Monte Carlo) simulation
! hopfraction    : Fraction of moves that are hops in VMC
! change_tstep   : Logical flag to allow changing tstep during run

! Energy Control
! -------------
! epscut         : Cutoff parameter for energy precision
! epstlrat       : Tolerance ratio for energy precision
! epscuttype     : Type of energy cutoff to use
! shift          : Energy shift parameter
! change_epscut  : Logical flag to allow changing epscut during run
! epsvar         : Variance threshold for energy

! Lattice Parameters
! ----------------
! alat2v         : Lattice parameter squared for VMC

! Optimization and Regularization
! ----------------------------
! theta_reg      : Regularization parameter
! typereg        : Type of regularization to use
! npow           : Power parameter for various calculations

! Algorithm Control
! --------------
! true_wagner    : Logical flag for Wagner algorithm
! cutweight      : Weight cutoff parameter
! nbra_cyrus     : Number of bra states for Cyrus algorithm

namelist /vmc/ tstep, hopfraction, epscut, epstlrat, epscuttype, alat2v, shift, change_epscut, change_tstep &
    &, epsvar, theta_reg, true_wagner, cutweight, nbra_cyrus, typereg, npow

! =============================================
! /dmclrdmc/ namelist parameters
! =============================================
! Energy and Energy Cutoff Parameters
! ----------------------------------
! etry           : Trial energy for DMC
! epscutdmc      : Energy cutoff for DMC
! epstldmc       : Tolerance for DMC
! gamma          : Damping parameter

! Lattice and Structure Parameters
! -------------------------------
! plat           : Lattice parameter
! alat2          : Lattice parameter squared
! alat           : Lattice parameter

! Time Step and Dynamics Parameters
! --------------------------------
! tstepfn        : Time step function parameter
! npow           : Power parameter for DMC
! tbra           : Bra state parameter
! Klrdmc         : K parameter for LR-DMC

! Optimization and Algorithm Control
! --------------------------------
! optbra         : Flag for optimizing bra states
! yesalfe        : Flag for alpha parameter optimization
! better_dmc     : Flag for improved DMC algorithm
! safelrdmc      : Flag for safe LR-DMC calculation
! changelambda    : Flag to allow lambda changes
! parcutg        : Parameter cutoff for gradient
! novar          : Flag to disable variance optimization
! typereg        : Type of regularization
! cutreg         : Cutoff for regularization

! Weight and Rejection Control
! --------------------------
! rejweight      : Rejection weight parameter
! cutweight      : Weight cutoff parameter
! weight_moroni  : Moroni weight parameter

! Walker and Population Control
! ---------------------------
! nw_max         : Maximum number of walkers
! nbra_cyrus     : Number of bra states for Cyrus algorithm
! noblocking     : Flag to disable blocking
! add_diff       : Flag to add diffusion term

! LR-DMC Specific Parameters
! ------------------------
! lrdmc_der      : Flag for LR-DMC derivatives
! lrdmc_nonodes  : Flag to disable nodes in LR-DMC
! enforce_detailb: Flag to enforce detailed balance
! iesrandoma     : Flag for random moves
! zmin           : Minimum z parameter
! yes_fastbranch : Flag for fast branching
! l0_kousuke     : Flag for Kousuke's l0 parameter
! true_wagner    : Flag for Wagner algorithm

namelist /dmclrdmc/ etry, npow, tbra, gamma, plat, alat2, alat        &
    &, tstepfn, Klrdmc, optbra, parcutg, novar, epscutdmc, typereg&
    &, epstldmc, rejweight, cutreg, cutweight, better_dmc, yesalfe&
    &, safelrdmc, changelambda, noblocking, add_diff, nbra_cyrus, lrdmc_der&
    &, lrdmc_nonodes, enforce_detailb, iesrandoma, zmin, yes_fastbranch&
    &, l0_kousuke, nw_max, true_wagner, weight_moroni


! =============================================
! /optimization/ namelist parameters
! =============================================
! Optimization Control Parameters
! -----------------------------
! tpar           : Optimization parameter
! nfat           : Number of force applications
! iboot          : Bootstrap parameter
! nweight        : Number of weights
! nmore_force    : Additional force parameter
! epsi           : Epsilon for optimization
! eps_dyn5       : Dynamic epsilon parameter
! epsdgel        : Epsilon for diagonal elements
! kl             : Kappa-lambda parameter
! idyn           : Dynamic optimization flag
! nbinr          : Number of bins for radial
! npbra          : Number of bra states for optimization
! ncg            : Number of conjugate gradient iterations

! Parameter Bounds and Limits
! -------------------------
! minz           : Minimum z parameter
! maxz           : Maximum z parameter
! minzj          : Minimum zj parameter
! maxzj          : Maximum zj parameter
! parr           : Parameter array
! parcute        : Parameter cutoff for energy
! parcut         : Parameter cutoff
! parcutmin      : Minimum parameter cutoff
! parcutpar      : Parameter cutoff parameter
! parr_max       : Maximum parameter array value
! parr_min       : Minimum parameter array value

! Molecular Dynamics and Cell Optimization
! --------------------------------------
! tion           : Temperature for ions
! tcell          : Temperature for cell
! molopt         : Molecular optimization flag
! epstion        : Epsilon for ion temperature

! One-body and Two-body Parameters
! ------------------------------
! onebodysz      : One-body size parameter
! twobodyoff     : Two-body offset parameter
! iesdtwobodyoff : Flag for two-body offset
! iesdonebodyoff : Flag for one-body offset
! minjonetwobody : Minimum one-body two-body parameter

! Symmetry and Cutting Parameters
! -----------------------------
! symiesup       : Symmetry flag
! yescutjas      : Flag to cut Jastrow parameters
! yescutdet      : Flag to cut determinant parameters
! fixpar         : Flag to fix parameters
! symmetrize_agp : Flag to symmetrize AGP

! Quantum and Bead Parameters
! -------------------------
! yesquantum     : Quantum optimization flag
! nbead          : Number of beads
! yeswritebead   : Flag to write bead information
! yesread10      : Flag to read from file 10

! Scaling and Complex Parameters
! ----------------------------
! oldscaling     : Old scaling flag
! srcomplex      : Stochastic reconfiguration complex flag
! scalermax      : Maximum scaling parameter

! Power and Signal Parameters
! -------------------------
! power          : Power parameter
! signalnoise    : Signal to noise ratio
! gauge_fixing   : Gauge fixing flag
! beta_learning  : Beta learning parameter

! Parameter Change Control
! ----------------------
! change_parr    : Flag to change parameter array
! delay_changeparr: Delay for parameter array change
! maxiter_changeparr: Maximum iterations for parameter array change
! change_tpar    : Flag to change tpar
! inc_tpar_frequency: Increment frequency for tpar
! use_stable_tpar: Flag to use stable tpar
! divide_tpar    : Flag to divide tpar
! multiply_tpar  : Flag to multiply tpar
! tpar_buffer_len: Buffer length for tpar
! tpar_max       : Maximum tpar value

! Advanced Optimization Parameters
! ------------------------------
! eps_umrigar    : Umrigar epsilon parameter
! yes_adams      : Adams method flag
! k6gen          : K6 generation parameter
! max_ortho      : Maximum orthogonality
! prep           : Preparation flag
! cut_sigma      : Sigma cutoff
! n_sigmas_tpar  : Number of sigmas for tpar
! len_tpar_stable_list: Length of stable tpar list
! yes_dgelscut   : Flag for DGELS cutoff
! tolcg          : Tolerance for conjugate gradient
! noopt_onebody  : Flag to disable one-body optimization

namelist /optimization/ tpar, nfat, iboot, nweight, nmore_force, epsi, eps_dyn5 &
    &, epsdgel, kl, idyn, nbinr, npbra, ncg, minz, maxz, minzj, maxzj, parr, parcute&
    &, parcut, parcutmin, parcutpar, tion, tcell, molopt, epstion&
    &, onebodysz, twobodyoff, iesdtwobodyoff, iesdonebodyoff, tolcg, minjonetwobody&
    &, symiesup, yescutjas, yescutdet, fixpar, symmetrize_agp&
    &, yesquantum, nbead, yeswrite10, oldscaling, srcomplex&
    &, power, signalnoise, gauge_fixing, beta_learning&
    &, noopt_onebody, scalermax, yeswritebead, yesread10, change_parr&
    &, parr_max, parr_min, delay_changeparr, maxiter_changeparr, k6gen, max_ortho, prep &
    &, change_tpar, inc_tpar_frequency, use_stable_tpar &
    &, eps_umrigar, yes_adams, divide_tpar, multiply_tpar, tpar_buffer_len, tpar_max &
    &, cut_sigma, n_sigmas_tpar, len_tpar_stable_list, yes_dgelscut
    

! =============================================
! /parameters/ namelist parameters
! =============================================
! Energy and Optimization Flags
! ---------------------------
! ieser          : Energy serial flag
! iesinv         : Energy inverse flag
! iesm           : Energy mean flag
! iesd           : Energy derivative flag
! isfix          : Fix flag
! iesfree        : Energy free flag
! iessw          : Energy switch flag
! iesup          : Energy update flag
! ieskin         : Kinetic energy flag

! Pressure and Cell Parameters
! --------------------------
! yespress       : Pressure flag
! warp           : Warp parameter
! powerwarp      : Power warp parameter
! pressfixed     : Fixed pressure flag
! fixa           : Fix lattice parameter a
! fixb           : Fix lattice parameter b
! fixc           : Fix lattice parameter c

! Pulay and Cell Dynamics
! ---------------------
! add_pulay      : Add Pulay correction flag
! yespulay       : Pulay correction flag
! typedyncell    : Type of dynamic cell
! scalepulay     : Scale Pulay parameter

! External Potential and Interactions
! --------------------------------
! ext_pot        : External potential flag
! vdw            : Van der Waals interaction flag
! link_atom      : Link atom flag
! mm_restr       : Molecular mechanics restraint flag

! Walk and Correction
! -----------------
! write_rwalk    : Write random walk flag
! yes_correct    : Correction flag

! Averaging and Statistics
! ----------------------
! yesavopt       : Average optimization flag
! yesavsr        : Average stochastic reconfiguration flag
! yesavcov       : Average covariance flag
! nrep_bead      : Number of replica beads

! Periodicity and K-points
! ----------------------
! yesperiodize   : Periodize flag
! yes_kpoints    : K-points flag

! Basis and Precision
! -----------------
! epsbas         : Epsilon for basis
! epsder         : Epsilon for derivatives

! Jastrow and AGP Parameters
! ------------------------
! yeszj          : Jastrow flag
! yeszagp        : AGP flag
! real_contracted: Real contracted flag
! real_agp       : Real AGP flag
! no_sjbra       : No SJ bra flag

! Run Control
! ----------
! decoupled_run  : Decoupled run flag
! scaleeloc      : Scale local energy
! cutoff_p       : Cutoff parameter
! read_molecul   : Read molecular flag

! Scemama Parameters
! ----------------
! yes_scemama    : Scemama flag
! yes_scemama_open: Open Scemama flag

namelist /parameters/ ieser, iesinv, iesm, iesd, isfix, iesfree &
    &, iessw, iesup, ieskin, yespress, warp, powerwarp&
    &, add_pulay, yespulay, typedyncell, scalepulay, ext_pot, vdw, link_atom&
    &, mm_restr, write_rwalk, yes_correct&
    &, yesavopt, yesavsr, yesavcov, nrep_bead, yesperiodize, yes_kpoints, epsbas&
    &, yeszj, yeszagp, decoupled_run, scaleeloc, cutoff_p, fixa, fixb, fixc&
    &, real_contracted, real_agp, no_sjbra, pressfixed, read_molecul, epsder, yes_scemama, yes_scemama_open


! =============================================
! /fitpar/ namelist parameters
! =============================================
! Parameter Fitting Control
! -----------------------
! nparinv        : Number of inverse parameters
! initparinv     : Initial inverse parameter
! rmaxinv        : Maximum radius for inverse
! npar           : Number of parameters
! initpar        : Initial parameter
! rmaxj          : Maximum radius for Jastrow
! nparsw         : Number of switch parameters
! initparsw      : Initial switch parameter
! rmax           : Maximum radius
! npower         : Number of powers
! powermin       : Minimum power
! npowersz       : Number of powers for size
! powerminsz     : Minimum power for size
! allfit         : Fit all parameters flag

namelist /fitpar/ nparinv, initparinv, rmaxinv, npar, initpar, rmaxj &
    &, nparsw, initparsw, rmax, npower, powermin, npowersz, powerminsz, allfit

! =============================================
! /dynamic/ namelist parameters
! =============================================
! Temperature and Friction
! ----------------------
! temp           : Temperature
! friction       : Friction coefficient

! Delta Parameters
! --------------
! delta0         : Initial delta parameter
! delta0q        : Quantum delta parameter
! delta0k        : Kinetic delta parameter

! Covariance and Scaling
! --------------------
! scalecov       : Scale covariance parameter
! scale_mass     : Scale mass parameter

! Dynamics Control
! --------------
! iskipdyn       : Skip dynamics flag
! maxdev_dyn     : Maximum deviation for dynamics
! stepcg_recount : Step conjugate gradient recount
! write_cov      : Write covariance flag
! normcorr       : Normalize correlation flag

! Turbo and Second Order
! --------------------
! yesturboq      : Turbo quantum flag
! yessecond      : Second order flag

! Cutting and Smoothing
! -------------------
! smoothcut      : Smooth cutoff flag
! killcut        : Kill cutoff flag

! Cell Equilibrium
! --------------
! eqcellab       : Equilibrium cell ab flag
! eqcellac       : Equilibrium cell ac flag
! eqcellbc       : Equilibrium cell bc flag

! Root and Rognoso
! --------------
! yesrootc       : Root c flag
! addrognoso     : Add Rognoso flag
! cleanrognoso   : Clean Rognoso flag

namelist /dynamic/ temp, friction, delta0, delta0q, delta0k, scalecov&
    &, iskipdyn, maxdev_dyn, stepcg_recount, write_cov, normcorr &
    &, yesturboq, yessecond, smoothcut, killcut, scale_mass, eqcellab &
    &, eqcellac, eqcellbc, yesrootc, addrognoso, cleanrognoso

! =============================================
! /unused/ namelist parameters
! =============================================
! Unused Parameters
! ----------------
! rsignr         : Random sign parameter
! beta           : Beta parameter
! testderiv      : Test derivative flag

namelist /unused/ rsignr, beta, testderiv

! =============================================
! /molecul/ namelist parameters
! =============================================
! Grid Parameters
! -------------
! nx             : Number of grid points in x
! ny             : Number of grid points in y
! nz             : Number of grid points in z
! nbufd          : Number of buffer dimensions

! Lattice Parameters
! ----------------
! ax             : Lattice parameter a in x
! ay             : Lattice parameter a in y
! az             : Lattice parameter a in z

! Molecular Control
! ---------------
! epsdgm         : Epsilon for molecular dynamics
! nmolmin        : Minimum number of molecules
! nmolmax        : Maximum number of molecules
! nmolmaxw       : Maximum number of molecules for walkers
! smearing       : Smearing parameter
! weight_loc     : Local weight parameter

! Orthogonality and Gram-Schmidt
! ----------------------------
! orthoyes       : Orthogonality flag
! epsrem_contr   : Epsilon for removing contractions
! gramyes        : Gram-Schmidt flag

! One-body and Origin
! -----------------
! add_onebody2det: Add one-body to determinant flag
! shift_origin   : Shift origin flag
! shiftx         : Shift in x direction
! shifty         : Shift in y direction
! shiftz         : Shift in z direction

namelist /molecul/ epsdgm, nx, ny, nz, nbufd, ax, ay, az, nmolmin, smearing&
    &, nmolmax, weight_loc, orthoyes, epsrem_contr, nmolmaxw, gramyes&
    &, add_onebody2det, shift_origin, shiftx, shifty, shiftz

! =============================================
! /link/ namelist parameters
! =============================================
! Link Parameters
! -------------
! calpha         : Alpha parameter for linking

namelist /link/ calpha

contains

    ! This module contains subroutines for handling I/O operations and data structures
    ! The original large allio module has been split into smaller files for faster compilation:
    ! - fort10_io.f90: Handles reading/writing of fort.10 files
    ! - fort11_io.f90: Handles reading/writing of fort.11 files 
    ! - read_datas.f90: Reads datas* control files
    ! - read_pseudo.f90: Reads pseudopotential files
    ! - memOP.f90: Default memory allocation/deallocation operations
    ! - writeoutput.f90: Output writing operations

    ! Contracts generalized Jastrow matrices
    ! Parameters:
    ! - nelorbh: Number of orbitals per spin channel
    ! - nelorb_c: Number of contracted orbitals
    ! - detmat: Output determinant matrix (2*nelorbh x 2*nelorbh)
    ! - detmat_c: Input contracted determinant matrix (2*nelorb_c x 2*nelorb_c)
    ! - mu_c: Contraction coefficients (nelorbh x nelorb_c)
    ! - psip: Temporary workspace (nelorbh x nelorb_c)
    subroutine scontract_genj(nelorbh, nelorb_c         &
         &, detmat, detmat_c, mu_c, psip)
        implicit none
        integer nelorbh, nelorb, nelorb_c, i
        real*8 detmat(2*nelorbh, 2*nelorbh), detmat_c(2*nelorb_c, 2*nelorb_c)&
             &, mu_c(nelorbh, *), psip(nelorbh, *)
!#ifdef __CASO
!    nprocu=nprocopt
!#else
!    nprocu=1
!#endif
        detmat = 0.d0
        ! Contract up-up block: detmat(1:nelorbh, 1:nelorbh) = mu_c * detmat_c(1:nelorb_c, 1:nelorb_c) * mu_c^T
        call dgemm_my('N', 'N', nelorbh, nelorb_c, nelorb_c, 1.d0, mu_c, nelorbh  &
             &, detmat_c, 2*nelorb_c, 0.d0, psip, nelorbh, nprocu, rankopt, commopt_mpi)
        call dgemm_my('N', 'T', nelorbh, nelorbh, nelorb_c, 1.d0, psip, nelorbh   &
             &, mu_c, nelorbh, 0.d0, detmat, 2*nelorbh, nprocu, rankopt, commopt_mpi)
!   down-down
        ! Contract down-down block: detmat(nelorbh+1:2*nelorbh, nelorbh+1:2*nelorbh) = mu_c * detmat_c(nelorb_c+1:2*nelorb_c, nelorb_c+1:2*nelorb_c) * mu_c^T
        call dgemm_my('N', 'N', nelorbh, nelorb_c, nelorb_c, 1.d0, mu_c, nelorbh  &
       &, detmat_c(nelorb_c + 1, nelorb_c + 1), 2*nelorb_c, 0.d0, psip, nelorbh, nprocu, rankopt, commopt_mpi)
        call dgemm_my('N', 'T', nelorbh, nelorbh, nelorb_c, 1.d0, psip, nelorbh   &
       &, mu_c, nelorbh, 0.d0, detmat(nelorbh + 1, nelorbh + 1), 2*nelorbh, nprocu, rankopt, commopt_mpi)

!   down-up
        ! Contract down-up block: detmat(nelorbh+1:2*nelorbh, 1:nelorbh) = mu_c * detmat_c(nelorb_c+1:2*nelorb_c, 1:nelorb_c) * mu_c^T
        call dgemm_my('N', 'N', nelorbh, nelorb_c, nelorb_c, 1.d0, mu_c, nelorbh  &
       &, detmat_c(nelorb_c + 1, 1), 2*nelorb_c, 0.d0, psip, nelorbh, nprocu, rankopt, commopt_mpi)
        call dgemm_my('N', 'T', nelorbh, nelorbh, nelorb_c, 1.d0, psip, nelorbh   &
       &, mu_c, nelorbh, 0.d0, detmat(nelorbh + 1, 1), 2*nelorbh, nprocu, rankopt, commopt_mpi)
!   up-down
        ! Contract up-down block: detmat(1:nelorbh, nelorbh+1:2*nelorbh) = mu_c * detmat_c(1:nelorb_c, nelorb_c+1:2*nelorb_c) * mu_c^T
        call dgemm_my('N', 'N', nelorbh, nelorb_c, nelorb_c, 1.d0, mu_c, nelorbh  &
       &, detmat_c(1, nelorb_c + 1), 2*nelorb_c, 0.d0, psip, nelorbh, nprocu, rankopt, commopt_mpi)
        call dgemm_my('N', 'T', nelorbh, nelorbh, nelorb_c, 1.d0, psip, nelorbh   &
       &, mu_c, nelorbh, 0.d0, detmat(1, nelorbh + 1), 2*nelorbh, nprocu, rankopt, commopt_mpi)
    end subroutine scontract_genj

    ! Contracts Jastrow matrices
    ! Parameters:
    ! - nelorbh: Number of orbitals per spin channel
    ! - nelorb: Total number of orbitals
    ! - nelcol: Number of columns
    ! - nelorb_c: Number of contracted orbitals
    ! - nelcol_c: Number of contracted columns
    ! - detmat: Output determinant matrix (nelorb x nelcol)
    ! - detmat_c: Input contracted determinant matrix (nelorb_c x nelcol_c)
    ! - mu_c: Contraction coefficients (nelorbh x nelorb_c)
    ! - psip: Temporary workspace (nelorbh x nelcol_c)
    subroutine scontract_mat_jas(nelorbh, nelorb, nelcol, nelorb_c         &
         &, nelcol_c, detmat, detmat_c, mu_c, psip)
        implicit none
        integer nelorbh, nelorb, nelorb_c, nelcol, nelcol_c, i
        real*8 detmat(nelorb, nelcol), detmat_c(nelorb_c, *)              &
             &, mu_c(nelorbh, *), psip(nelorbh, *)
! WARNING  IT REFERS ONLY TO jASTROW SO IPF AND UNPAIRED ARE NOT USEd
!#ifdef __CASO
!    nprocu=nprocopt
!#else
!    nprocu=1
!#endif
        detmat = 0.d0
        ! Perform matrix contractions for Jastrow part: detmat = mu_c * detmat_c * mu_c^T
        call dgemm_my('N', 'N', nelorbh, nelcol_c, nelorb_c, 1.d0, mu_c, nelorbh  &
             &, detmat_c, nelorb_c, 0.d0, psip, nelorbh, nprocu, rankopt, commopt_mpi)
        call dgemm_my('N', 'T', nelorbh, nelorbh, nelorb_c, 1.d0, psip, nelorbh   &
             &, mu_c, nelorbh, 0.d0, detmat, nelorb, nprocu, rankopt, commopt_mpi)
    end subroutine scontract_mat_jas

    ! Contracts determinant matrices
    ! Parameters:
    ! - nelorbh: Number of orbitals per spin channel
    ! - nelorb: Total number of orbitals
    ! - nelcol: Number of columns
    ! - nelorb_c: Number of contracted orbitals
    ! - nelcol_c: Number of contracted columns
    ! - detmat: Output determinant matrix (ipc*ipf*nelorb x nelcol)
    ! - detmat_c: Input contracted determinant matrix (ipc*nelorb_c x nelcol_c)
    ! - mu_c: Contraction coefficients (ipc*ipf*nelorbh x nelorb_c)
    ! - psip: Temporary workspace (ipf*ipc*nelorbh x nelcol_c)
    subroutine scontract_mat_det(nelorbh, nelorb, nelcol, nelorb_c         &
         &, nelcol_c, detmat, detmat_c, mu_c, psip)
        implicit none
        integer nelorbh, nelorb, nelorb_c, nelcol, nelcol_c, i
        real*8 detmat(ipc*ipf*nelorb, nelcol), detmat_c(ipc*nelorb_c, *)              &
             &, mu_c(ipc*ipf*nelorbh, *), psip(ipf*ipc*nelorbh, *)
!#ifdef __CASO
!    nprocu=nprocopt
!#else
!    nprocu=1
!#endif
        if (ipc .eq. 2) then
            ! Complex case: Use complex BLAS operations (zgemm_my)
            detmat = 0.d0
            call zgemm_my('N', 'N', ipf*nelorbh, nelcol_c, nelorb_c, (1.d0, 0.d0), mu_c, ipf*nelorbh  &
                 &, detmat_c, nelorb_c, (0.d0, 0.d0), psip, ipf*nelorbh, nprocu, rankopt, commopt_mpi)
            call zgemm_my('N', 'T', ipf*nelorbh, ipf*nelorbh, nelorb_c, (1.d0, 0.d0), psip, ipf*nelorbh   &
                 &, mu_c, ipf*nelorbh, (0.d0, 0.d0), detmat, ipf*nelorb, nprocu, rankopt, commopt_mpi)
            if (nelcol_c .gt. nelorb_c) then
                ! Copy remaining columns for complex case
                do i = nelorb_c + 1, nelcol_c
                    call zcopy(ipf*nelorbh, psip(1, i), 1, detmat(1, ipf*nelorb + i - nelorb_c), 1)
                end do
            end if
        else
            ! Real case: Use real BLAS operations (dgemm_my)
            detmat = 0.d0
            call dgemm_my('N', 'N', ipf*nelorbh, nelcol_c, nelorb_c, 1.d0, mu_c, ipf*nelorbh&
                 &, detmat_c, nelorb_c, 0.d0, psip, ipf*nelorbh, nprocu, rankopt, commopt_mpi)
            call dgemm_my('N', 'T', ipf*nelorbh, ipf*nelorbh, nelorb_c, 1.d0, psip, ipf*nelorbh&
                    &, mu_c, ipf*nelorbh, 0.d0, detmat, ipf*nelorb, nprocu, rankopt, commopt_mpi)
            if (nelcol_c .gt. nelorb_c) then
                ! Copy remaining columns for real case
                do i = nelorb_c + 1, nelcol_c
                    call dcopy(ipf*nelorbh, psip(1, i), 1, detmat(1, nelorb*ipf + i - nelorb_c), 1)
                end do
            end if
        end if
    end subroutine scontract_mat_det

    ! Updates the k-point grid for periodic calculations
    ! This subroutine:
    ! 1. Initializes direct lattice vectors for periodic basis set
    ! 2. Builds lattice vector map for basis set computation
    ! 3. Optimizes number of vectors in summation
    ! 4. Handles both real and complex wave functions
    ! 5. Supports both open and periodic boundary conditions
    subroutine update_kgrid
        implicit none
        integer i, j, ii, jj, kk, ll, count1, count2, count1j, count2j, indpar, indparp, kboundi&
       &, ind, kboundi_max, iii, jjj, kkk, k, maxdim
        integer, dimension(:, :), allocatable :: kshell_map
        real*8 kbound, map_tmp(3, 27), max_rejected, cost_tilted
        logical, external :: slaterorb
        logical not_found
        integer, dimension(:, :), allocatable :: kpip_sav

        ! Initialize direct lattice vectors for periodic basis set
        ! Uses same definition as Crystal DFT code
        ! Can be used for both real and complex wave functions
        ! For complex wave functions, can also handle open systems
        ! Use "PBC_C" keyword in wave function first line to enable
        if (allocated(kgrid)) deallocate (kgrid)
        if (allocated(kgrid_atom)) deallocate (kgrid_atom)
        ikshift = nshell
        if (abs(LBox) .eq. 3.d0) then
            ! Determine the maximum number of direct lattice vectors for each shell
            ! This is based on the basis set parameters and cutoff criteria
            if (yes_crystalj) then
                allocate (kshell_map(3, nshell + nshellj))
            else
                allocate (kshell_map(3, nshell))
            end if
            kshell_map = 0
            indpar = 0
            do i = 1, nshell
                indparp = indpar + 1
                if (lepsbas .gt. 0.d0) then
                    do j = 1, 3
                        if (cellscale(j) .ne. 0.d0) then
                            if (slaterorb(ioptorb(i))) then ! STO: Slater-type orbital
                                kbound = 0.5d0 + (lepsbas/dupr(indparp)/metric_min)/cellscale(j)
                            else ! GTO: Gaussian-type orbital
                                kbound = 0.5d0 + (dsqrt(lepsbas/dupr(indparp))/metric_min)/cellscale(j)
                            end if
                            kshell_map(j, i) = kbound + 1
                        end if
                    end do
                else
                    do j = 1, 3
                        kshell_map(j, i) = 0
                    end do
                end if
                indpar = indpar + nparam(i)
            end do
            if (yes_crystalj) then
                ! Handle Jastrow functions with periodic boundary conditions
                indpar = 0
                do i = 1, nshellj
                    indparp = indpar + 1
                    if (lepsbas .gt. 0.d0 .and. ioptorbj(i) .ne. 200) then
                        do j = 1, 3
                            if (cellscale(j) .ne. 0.d0) then
                                if (slaterorb(ioptorbj(i))) then ! STO for Jastrow
                                    kbound = 0.5d0 + (lepsbas/vjur(indparp)/metric_min)/cellscale(j)
                                else ! GTO for Jastrow
                                    kbound = 0.5d0 + (dsqrt(lepsbas/vjur(indparp))/metric_min)/cellscale(j)
                                end if
                                kshell_map(j, i + ikshift) = kbound + 1
                            end if
                        end do
                    else
                        do j = 1, 3
                            kshell_map(j, i + ikshift) = 0
                        end do
                    end if
                    indpar = indpar + nparamj(i)
                end do
            end if

            ! Build the lattice vector map needed to compute the basis set
            ! Allocate kgrid structures for shells and atoms
            if (yes_crystalj) then
                allocate (kgrid(nshell + nshellj))
                allocate (kgrid_atom(2*nion))
            else
                allocate (kgrid(nshell))
                allocate (kgrid_atom(nion))
            end if

            ! Handle 1D and 2D systems by setting z and y components to zero
            if (yes2d .or. yes1d) kshell_map(3, :) = 0
            if (yes1d) kshell_map(2, :) = 0

            if (iespbc) then
            ! Periodic boundary conditions: Build comprehensive grid for each ion
            do ii = 1, nion
                kboundi_max = 0
                do jj = adr_nion(ii), adr_nion(ii + 1) - 1
                    i = ind_nion(jj)
!    do i=1,nshell
                    kboundi = (2*kshell_map(1, i) + 1)*(2*kshell_map(2, i) + 1)*(2*kshell_map(3, i) + 1)
                    kboundi_max = max(kboundi_max, kboundi)
                    allocate (kgrid(i)%kpip(3, kboundi))
                    kgrid(i)%kpip(:, :) = 0
                end do
                allocate (kgrid_atom(ii)%kpip(3, kboundi_max))
                kgrid(ii)%kpip(:, :) = 0
            end do
            if (yes_crystalj) then
            ! Handle Jastrow functions for periodic boundary conditions
            do ii = 1, nion
                kboundi_max = 0
                do jj = adrj_nion(ii), adrj_nion(ii + 1) - 1
                    i = indj_nion(jj) + nshell
                    kboundi = (2*kshell_map(1, i) + 1)*(2*kshell_map(2, i) + 1)*(2*kshell_map(3, i) + 1)
                    kboundi_max = max(kboundi_max, kboundi)
                    allocate (kgrid(i)%kpip(3, kboundi))
                    kgrid(i)%kpip(:, :) = 0
                end do
                allocate (kgrid_atom(ii + nion)%kpip(3, kboundi_max))
                kgrid(ii + nion)%kpip(:, :) = 0
            end do
            end if
            else
            ! Open boundary conditions: Simple grid with only origin
            do i = 1, nshell
                allocate (kgrid(i)%kpip(3, 1))
                kgrid(i)%kpip(:, :) = 0
            end do
            end if

            max_rejected = lepsbas
            count1 = 0
            count2 = 0

            if (iespbc) then
                ! Optimize the number of vectors in the summation to be inside a sphere
                ! This reduces computational cost while maintaining accuracy
                indpar = 0
                do iii = 1, nion
                    kgrid_atom(iii)%dimshell = 0
                    do jjj = adr_nion(iii), adr_nion(iii + 1) - 1
                        i = ind_nion(jjj)
!       do i=1,nshell
                        indparp = indpar + 1
                        kgrid(i)%dimshell = 0
                        do kk = -kshell_map(3, i), kshell_map(3, i)
                            do jj = -kshell_map(2, i), kshell_map(2, i)
                                do ii = -kshell_map(1, i), kshell_map(1, i)
                                    ! Optimize the number of vectors in the summation
                                    ! to be inside a sphere based on metric cutoff
                                    count2 = count2 + 1
!In order to estimate a lower bound of exp(- Z_basis  |z|)
!for each ii,jj,kk it is computed  the minimum distance of a point z from origin
!    dist  = [ (x+ii Lx,y+jj Ly, z+kk Lz, metric (x+ii Lx, y+jj Ly, z+kk Lz)]
!  with the  condition |x|<Lx/2 |y|<Ly/2 |z|<Lz/2
! The minimum is typically at the boundaries x=+/-Lx/2,y=+/-Ly/2,z=+/-Lz/2
! apart cases when ii=0 or jj=0 or kk=0, when,  under some conditions
!  depending on the metric, x=0 or y=0 or z=0.

                                    map_tmp(1, 1) = cellscale(1)*ii
                                    map_tmp(2, 1) = cellscale(2)*jj
                                    map_tmp(3, 1) = cellscale(3)*kk

                                    call prep_map(map_tmp, cellscale)

! NB to avoid roundoff problems ii,jj,kk--> min(abs(ii),1) =0 ,1 only.

                                    cost_tilted = norm_metric(map_tmp(1, 1), metric)**2
                                    do ll = 2, 27
                                        cost_tilted = min(cost_tilted, norm_metric(map_tmp(1, ll), metric)**2)
                                    end do

                                    if (slaterorb(ioptorb(i))) then
                                        kbound = dupr(indparp)*dsqrt(cost_tilted)
                                    else
                                        kbound = dupr(indparp)*cost_tilted
                                    end if
                                    if (kbound .lt. lepsbas .or. lepsbas .le. 0.d0) then
                                        count1 = count1 + 1
                                        kgrid(i)%dimshell = kgrid(i)%dimshell + 1
                                        kgrid(i)%kpip(1, kgrid(i)%dimshell) = ii
                                        kgrid(i)%kpip(2, kgrid(i)%dimshell) = jj
                                        kgrid(i)%kpip(3, kgrid(i)%dimshell) = kk

                                        not_found = .true.
                                        do kkk = 1, kgrid_atom(iii)%dimshell
                                            if (kgrid_atom(iii)%kpip(1, kkk) .eq. ii .and.&
                                           & kgrid_atom(iii)%kpip(2, kkk) .eq. jj .and.&
                                           & kgrid_atom(iii)%kpip(3, kkk) .eq. kk) not_found = .false.
                                        end do
                                        if (not_found) then
                                            kgrid_atom(iii)%dimshell = kgrid_atom(iii)%dimshell + 1
                                            kgrid_atom(iii)%kpip(1, kgrid_atom(iii)%dimshell) = ii
                                            kgrid_atom(iii)%kpip(2, kgrid_atom(iii)%dimshell) = jj
                                            kgrid_atom(iii)%kpip(3, kgrid_atom(iii)%dimshell) = kk
                                        end if
                                    else
                                        if (kbound .lt. max_rejected .or. max_rejected .eq. lepsbas)&
                                        &max_rejected = kbound
                                    end if
                                end do
                            end do
                        end do
                        indpar = indpar + nparam(i)
                    end do ! end ii
                end do ! end ion
!       Shrink the memory allocated to optimize memory usage
                maxdim = 0
                do i = 1, nion
                    maxdim = max(maxdim, kgrid_atom(i)%dimshell)
                end do
                allocate (kpip_sav(3, maxdim))
                do i = 1, nion
                    do j = 1, kgrid_atom(i)%dimshell
                        kpip_sav(:, j) = kgrid_atom(i)%kpip(:, j)
                    end do
                    deallocate (kgrid_atom(i)%kpip)
                    allocate (kgrid_atom(i)%kpip(3, kgrid_atom(i)%dimshell))
                    do j = 1, kgrid_atom(i)%dimshell
                        kgrid_atom(i)%kpip(:, j) = kpip_sav(:, j)
                    end do
                end do
                do i = 1, nshell
                    do j = 1, kgrid(i)%dimshell
                        kpip_sav(:, j) = kgrid(i)%kpip(:, j)
                    end do
                    deallocate (kgrid(i)%kpip)
                    allocate (kgrid(i)%kpip(3, kgrid(i)%dimshell))
                    do j = 1, kgrid(i)%dimshell
                        kgrid(i)%kpip(:, j) = kpip_sav(:, j)
                    end do
                end do
                deallocate (kpip_sav)
                count1j = 0
                count2j = 0
                if (yes_crystalj) then
                    ! Handle Jastrow functions with periodic boundary conditions
                    indpar = 0
                    do iii = 1, nion
                        kgrid_atom(iii + nion)%dimshell = 0
                        do jjj = adrj_nion(iii), adrj_nion(iii + 1) - 1
!       do i=1,nshellj
                            i = indj_nion(jjj)
                            indparp = indpar + 1
                            kgrid(i + ikshift)%dimshell = 0
                            do kk = -kshell_map(3, i + ikshift), kshell_map(3, i + ikshift)
                                do jj = -kshell_map(2, i + ikshift), kshell_map(2, i + ikshift)
                                    do ii = -kshell_map(1, i + ikshift), kshell_map(1, i + ikshift)
                                        ! Optimize the number of vectors in the summation
                                        ! to be inside a sphere for Jastrow functions
                                        count2j = count2j + 1
!  for each ii,jj,kk compute  the minimum distance of a point z from origin
!    dist  = [ (x+ii Lx,y+jj Ly, z+kk Lz, metric (x+ii Lx, y+jj Ly, z+kk Lz)]
!  with the  condition |x|<Lx/2 |y|<Ly/2 |z|<Lz/2
! The minimum is typically at the boundaries x=+/-Lx/2,y=+/-Ly/2,z=+/-Lz/2
! apart cases where ii=0 or jj=0 or kk=0, when,  under some conditions
!  depending on the metric, x=0 or y=0 or z=0.

                                        map_tmp(1, 1) = cellscale(1)*ii
                                        map_tmp(2, 1) = cellscale(2)*jj
                                        map_tmp(3, 1) = cellscale(3)*kk

                                        call prep_map(map_tmp, cellscale)

                                        cost_tilted = norm_metric(map_tmp(1, 1), metric)**2
                                        do ll = 2, 27
                                            cost_tilted = min(cost_tilted, norm_metric(map_tmp(1, ll), metric)**2)
                                        end do

                                        if (ioptorbj(i) .ne. 200) then
                                        if (slaterorb(ioptorbj(i))) then
                                            kbound = vjur(indparp)*dsqrt(cost_tilted)
                                        else
                                            kbound = vjur(indparp)*cost_tilted
                                        end if
                                        else
                                        kbound = 0.d0
                                        end if
                                        if (kbound .lt. lepsbas .or. lepsbas .le. 0.d0) then
                                            count1j = count1j + 1
                                            kgrid(i + ikshift)%dimshell = kgrid(i + ikshift)%dimshell + 1
                                            kgrid(i + ikshift)%kpip(1, kgrid(i + ikshift)%dimshell) = ii
                                            kgrid(i + ikshift)%kpip(2, kgrid(i + ikshift)%dimshell) = jj
                                            kgrid(i + ikshift)%kpip(3, kgrid(i + ikshift)%dimshell) = kk
                                            not_found = .true.
                                            do kkk = 1, kgrid_atom(iii + nion)%dimshell
                                                if (kgrid_atom(iii + nion)%kpip(1, kkk) .eq. ii .and.&
                                               & kgrid_atom(iii + nion)%kpip(2, kkk) .eq. jj .and.&
                                         & kgrid_atom(iii + nion)%kpip(3, kkk) .eq. kk) not_found = .false.
                                            end do
                                            if (not_found) then
                                                kgrid_atom(iii + nion)%dimshell = kgrid_atom(iii + nion)%dimshell + 1
                                                kgrid_atom(iii + nion)%kpip(1, kgrid_atom(iii + nion)%dimshell) = ii
                                                kgrid_atom(iii + nion)%kpip(2, kgrid_atom(iii + nion)%dimshell) = jj
                                                kgrid_atom(iii + nion)%kpip(3, kgrid_atom(iii + nion)%dimshell) = kk
                                            end if
                                        else
                                            if (kbound .lt. max_rejected .or.&
                                                & max_rejected .eq. lepsbas .and.&
                                                & ioptorbj(i) .ne. 200) &
                                                &max_rejected = kbound
                                        end if
                                    end do
                                end do
                            end do
                            indpar = indpar + nparamj(i)
                        end do ! ii
                    end do ! nion
!       Shrink the memory allocated for Jastrow functions
                    maxdim = 0
                    do i = 1, nion
                        maxdim = max(maxdim, kgrid_atom(i + nion)%dimshell)
                    end do
                    allocate (kpip_sav(3, maxdim))
                    do i = 1, nion
                        do j = 1, kgrid_atom(i + nion)%dimshell
                            kpip_sav(:, j) = kgrid_atom(i + nion)%kpip(:, j)
                        end do
                        deallocate (kgrid_atom(i + nion)%kpip)
                        allocate (kgrid_atom(i + nion)%kpip(3, kgrid_atom(i + nion)%dimshell))
                        do j = 1, kgrid_atom(i + nion)%dimshell
                            kgrid_atom(i + nion)%kpip(:, j) = kpip_sav(:, j)
                        end do
                    end do
                    do i = 1, nshellj
                        do j = 1, kgrid(i + ikshift)%dimshell
                            kpip_sav(:, j) = kgrid(i + ikshift)%kpip(:, j)
                        end do
                        deallocate (kgrid(i + ikshift)%kpip)
                        allocate (kgrid(i + ikshift)%kpip(3, kgrid(i + ikshift)%dimshell))
                        do j = 1, kgrid(i + ikshift)%dimshell
                            kgrid(i + ikshift)%kpip(:, j) = kpip_sav(:, j)
                        end do
                    end do
                    deallocate (kpip_sav)
                end if
            else
                ! Open boundary conditions: Simple grid with only origin point
                count1 = 0
                do i = 1, nshell
                    kgrid(i)%dimshell = 1
                    kgrid(i)%kpip(1, kgrid(i)%dimshell) = 0
                    kgrid(i)%kpip(2, kgrid(i)%dimshell) = 0
                    kgrid(i)%kpip(3, kgrid(i)%dimshell) = 0
                end do
                cellscale = 0.d0 ! initialize to be sure the irrelevant vector cellscale
            end if
            deallocate (kshell_map)
            if (iespbc) then
            ! Set up tobedone flags for efficient computation
            do ii = 1, nion
                do jj = adr_nion(ii), adr_nion(ii + 1) - 1
                    i = ind_nion(jj)
                    allocate (kgrid(i)%tobedone(kgrid_atom(ii)%dimshell))
                    kgrid(i)%tobedone(:) = .false.
                end do
            end do
            if (yes_crystalj) then
            do ii = 1, nion
                do jj = adrj_nion(ii), adrj_nion(ii + 1) - 1
                    i = indj_nion(jj) + nshell
                    allocate (kgrid(i)%tobedone(kgrid_atom(ii + nion)%dimshell))
                    kgrid(i)%tobedone(:) = .false.
                end do
            end do
            end if
            do iii = 1, nion
                do jjj = adr_nion(iii), adr_nion(iii + 1) - 1
                    i = ind_nion(jjj)
                    do j = 1, kgrid_atom(iii)%dimshell
                        do k = 1, kgrid(i)%dimshell
                        if (kgrid_atom(iii)%kpip(1, j) .eq. kgrid(i)%kpip(1, k)&
                        &.and. kgrid_atom(iii)%kpip(2, j) .eq. kgrid(i)%kpip(2, k)&
                        &.and. kgrid_atom(iii)%kpip(3, j) .eq. kgrid(i)%kpip(3, k)) then
                            kgrid(i)%tobedone(j) = .true.
                        end if
                        end do
                    end do
                end do
            end do
            if (yes_crystalj) then
                do iii = 1, nion
                    do jjj = adrj_nion(iii), adrj_nion(iii + 1) - 1
                        i = indj_nion(jjj)
                        do j = 1, kgrid_atom(iii + nion)%dimshell
                            do k = 1, kgrid(i + ikshift)%dimshell
                                if (kgrid_atom(iii + nion)%kpip(1, j) .eq. kgrid(i + ikshift)%kpip(1, k)&
                               &.and. kgrid_atom(iii + nion)%kpip(2, j) .eq. kgrid(i + ikshift)%kpip(2, k)&
                               &.and. kgrid_atom(iii + nion)%kpip(3, j) .eq. kgrid(i + ikshift)%kpip(3, k)) then
                                    kgrid(i + ikshift)%tobedone(j) = .true.
                                end if
                            end do
                        end do
                    end do
                end do
            end if
            end if
            if (rank .eq. 0) then
                write (6, *) 'Warning: updated kgrid considering Det/Jastrow:', count1, count1j
                write (6, *) 'Warning: lowest wf discarded=', exp(-max_rejected)
#ifdef _DEBUG
                write (6, *) ' Grid considered /ion'
                do i = 1, nion
                    write (6, *) ' ion # ', i, 'dimension =', kgrid_atom(i)%dimshell
                    do j = 1, kgrid_atom(i)%dimshell
                        write (6, *) j, kgrid_atom(i)%kpip(1, j), kgrid_atom(i)%kpip(2, j), kgrid_atom(i)%kpip(3, j)
                    end do
                end do
                if (yes_crystalj) then
                    write (6, *) ' Grid considered Jastrow/ion '
                    do i = nion + 1, 2*nion
                        write (6, *) ' ion # ', i - nion, 'dimension =', kgrid_atom(i)%dimshell
                        do j = 1, kgrid_atom(i)%dimshell
                            write (6, *) j, kgrid_atom(i)%kpip(1, j), kgrid_atom(i)%kpip(2, j), kgrid_atom(i)%kpip(3, j)
                        end do
                    end do
                end if
                write (6, *) ' Grid considered '
                do i = 1, nshell
                    write (6, *) ' Shell # ', i, 'dimension =', kgrid(i)%dimshell
                    do j = 1, kgrid(i)%dimshell
                        write (6, *) j, kgrid(i)%kpip(1, j), kgrid(i)%kpip(2, j), kgrid(i)%kpip(3, j)
                    end do
                end do
                if (yes_crystalj) then
                    write (6, *) ' Grid considered Jastrow '
                    do i = 1, nshellj
                        write (6, *) ' Shell # ', i, 'dimension =', kgrid(i + ikshift)%dimshell
                        do j = 1, kgrid(i + ikshift)%dimshell
                            write (6, *) j, kgrid(i + ikshift)%kpip(1, j), &
                                kgrid(i + ikshift)%kpip(2, j), kgrid(i + ikshift)%kpip(3, j)
                        end do
                    end do
                end if
#endif
            end if
!    kgrid%kpip is no longer needed after optimization
            do i = 1, nshell
                if (allocated(kgrid(i)%kpip)) deallocate (kgrid(i)%kpip)
            end do
            if (yes_crystalj) then
            do i = nshell + 1, nshell + nshellj
                if (allocated(kgrid(i)%kpip)) deallocate (kgrid(i)%kpip)
            end do
            end if
        end if
    end subroutine update_kgrid

    ! =============================================
    ! norm_metric - Calculate Norm Using Metric Tensor
    ! =============================================
    ! Purpose: Calculates the norm of a 3D vector using a given metric tensor
    ! 
    ! Mathematical Definition:
    ! ||r||_metric = sqrt(r^T * metric * r)
    ! where r is a 3D vector and metric is a 3x3 symmetric matrix
    !
    ! Parameters:
    ! - r: 3D position vector (input)
    ! - metric: 3x3 metric tensor (input)
    ! Returns: Norm of the vector in the given metric (real*8)
    !
    ! Algorithm:
    ! 1. Compute quadratic form: r^T * metric * r
    ! 2. Take square root with protection from roundoff errors
    ! 3. Return the norm
    !
    ! Usage:
    ! This function is used in periodic boundary condition calculations
    ! to determine the distance between lattice points in a general metric
    function norm_metric(r, metric)
        implicit none
        real*8 norm_metric, r(3), metric(3, 3)
        ! Calculate norm using metric tensor: ||r|| = sqrt(r^T * metric * r)
        norm_metric = metric(1, 1)*r(1)*r(1) + metric(2, 2)*r(2)*r(2) + metric(3, 3)*r(3)*r(3)&
                & + 2.d0*(metric(1, 2)*r(1)*r(2) + metric(1, 3)*r(1)*r(3) + metric(2, 3)*r(2)*r(3))
        norm_metric = dsqrt(max(norm_metric, 0.d0)) ! Protection from roundoff errors
        return
    end

end module allio

! =============================================
! prep_map - Prepare Mapping of Lattice Vectors
! =============================================
! Purpose: Prepares mapping of lattice vectors for periodic boundary condition calculations
! 
! Algorithm:
! 1. For each coordinate (x, y, z), there are three possible minimum positions:
!    - Current position (when coordinate is zero)
!    - Lower boundary (-L/2)
!    - Upper boundary (+L/2)
! 2. This results in 27 possible combinations (3^3)
! 3. The function generates all these combinations to find the minimum distance
!
! Mathematical Background:
! When computing distances in periodic systems, the minimum distance
! between two points can occur at different positions within the unit cell.
! This subroutine explores all possible combinations to find the true minimum.
!
! Parameters:
! - map_tmp: Array to store mapped vectors (3 x 3 x 3 x 3)
! - cellscale: Cell scaling factors (3)
!
! Output:
! - map_tmp: Contains all 27 possible vector positions for distance calculation
subroutine prep_map(map_tmp, cellscale)
    implicit none
    real*8 map_tmp(3, 3, 3, 3), cellscale(3), cellhalf(3, 3)
    integer i, j, k
    ! Find possible minimum of the metric
    ! For each coordinate there are two possibilities:
    ! 1. Minimum at boundary +/- L/2
    ! 2. Minimum at current position (when zero)
    ! Results in 27 possibilities (3 for each coordinate)
    ! including the input map(:,1,1,1)
    cellhalf(:, 1) = 0.d0
    cellhalf(:, 2) = -cellscale(:)/2.d0
    cellhalf(:, 3) = cellscale(:)/2.d0

    ! Generate all possible combinations of boundary conditions
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                map_tmp(1, i, j, k) = map_tmp(1, 1, 1, 1) + cellhalf(1, i)
                map_tmp(2, i, j, k) = map_tmp(2, 1, 1, 1) + cellhalf(2, j)
                map_tmp(3, i, j, k) = map_tmp(3, 1, 1, 1) + cellhalf(3, k)
            end do
        end do
    end do
    return
end
