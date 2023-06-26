.. TurboRVB_manual documentation master file, created by
   sphinx-quickstart on Thu Jan 24 00:11:17 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Initializations of Wavefunctions
================================================

External quantum chemistry and DFT codes (Gaussian, GAMESS, and pySCF)
----------------------------------------------------------------------------
You can use external quantum chemistry and DFT codes such as Gaussian, GAMESS, and pySCF to generate a trail wavefunction. Please refer to the TurboGenius tutorial [https://github.com/kousuke-nakano/turbotutorials].

Built-in DFT code (prep.x)
------------------------------------------------
In this section, we describe the INPUT file that controls DFT calculation. The input file is built using Fortran namelists. Keywords are divided into different sections according to their meanings.

Simulation section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The sections and the main parameters are listed and shortly described in the following. Keywords that are not required can be omitted in the input file: they will assume a default value.

.. table:: Simulation section

   +-----------+----------+---------+------------------------------------------------------+
   | Parameter | Datatype | Default | Description                                          |
   +===========+==========+=========+======================================================+
   | itest4    | NA       | NA      |                                                      |
   |           |          |         | ``itest4 = -4``: Standard DFT run.                   |
   |           |          |         | ``itest4 = -8``: DFT calculation with twice larger   |
   |           |          |         | basis for the Hartree potential.                     |
   +-----------+----------+---------+------------------------------------------------------+
   | iopt      | NA       | NA      |                                                      |
   |           |          |         | ``iopt = 1``: Initialize with no potential (no       |
   |           |          |         | Hartree, xc, correlation).                           |
   |           |          |         | ``iopt = 0``: Continuation, starting from            |
   |           |          |         | wavefunctions read from ``fort.10_new`` and          |
   |           |          |         | occupation read from ``occupationlevels.dat`` (both  |
   |           |          |         | generated after ``iopt = 1`` run).                   |
   |           |          |         | ``iopt = 2``: Same as ``iopt = 0``, but write main   |
   |           |          |         | matrices (basis set/Hamiltonian overlaps, charge/spin|
   |           |          |         | density).                                            |
   +-----------+----------+---------+------------------------------------------------------+

Pseudo section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. table:: Pseudo section (summary)

    +------------------+------------+----------+---------------------------------------------------------------------+
    | Parameter        | Datatype   | Default  | Description                                                         |
    +==================+============+==========+=====================================================================+
    | ``nintpsa``      | NA         | NA       | Number of integer points for pseudopotential if present.            |
    +------------------+------------+----------+---------------------------------------------------------------------+
    | ``pseudorandom`` | NA         | NA       | Use a random integration mesh for pseudo with the algorithm for QMC |
    |                  |            |          | by R. Fahy.                                                         |
    +------------------+------------+----------+---------------------------------------------------------------------+
    | ``npsamax``      | NA         | NA       | Multiplication factor for the number of pseudo integration points.  |
    |                  |            |          | Note that, use ``npsmax > 2`` if the code terminates with the error |
    |                  |            |          | 'Increase npsamax'.                                                 |
    +------------------+------------+----------+---------------------------------------------------------------------+

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Optimization section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. table:: Optimization section

    +----------------+----------+---------+--------------------------------------------------+
    | Parameter      | Datatype | Default | Description                                      |
    +================+==========+=========+==================================================+
    | ``molopt``     | NA       | NA      | Do not change this value, as DFT works with      |
    |                |          |         | molecular orbitals only.                         |
    +----------------+----------+---------+--------------------------------------------------+

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Readio section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. table:: Readio section

   +----------------+----------+---------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | parameter name | datatype | default | description                                                                                                                                                                                                                                                                                                                                                                                                                             |
   +================+==========+=========+=========================================================================================================================================================================================================================================================================================================================================================================================================================================+
   | writescratch   | NA       | NA      | ``writescratch = 0``: Writes binary scratch files on disk to speed up continuation and allow non-self-consistent calculations and post-processing tools. The following files are written:                                                                                                                                                                                                                                               |
   |                |          |         |                                                                                                                                                                                                                                                                                                                                                                                                                                         |
   |                |          |         | - **tmp000xxx:** One for each processor. These files contain basis/Hamiltonian overlap matrix elements in the first record. In the second record, they contain charge/spin density distributed matrices.                                                                                                                                                                                                                                |
   |                |          |         |                                                                                                                                                                                                                                                                                                                                                                                                                                         |
   |                |          |         | - **total_densities.sav:** A single file containing the total (i.e., not distributed over the real space integration grid) charge density in the first record and the total spin density in the second record (in the case of LSDA calculations).                                                                                                                                                                                       |
   |                |          |         |                                                                                                                                                                                                                                                                                                                                                                                                                                         |
   |                |          |         | - **wavefunction.sav:** A single file containing the final Kohn-Sham eigenvectors for all the bands. If k-points are present, each record contains the eigenvectors for a single k-point in the order as they appear in the **occupationlevels.dat** file.                                                                                                                                                                              |
   |                |          |         |                                                                                                                                                                                                                                                                                                                                                                                                                                         |
   |                |          |         | ``writescratch = 1``: Does not write any scratch file on disk. Continuing from previous runs and non-self-consistent calculations will not be possible.                                                                                                                                                                                                                                                                                 |
   +----------------+----------+---------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Parameters section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. table:: Parameters section

    +------------------+----------+---------+-----------------------------------------------------------------------------------------------------------------------------+
    | Parameter        | Datatype | Default | Description                                                                                                                 |
    +==================+==========+=========+=============================================================================================================================+
    | ``decoupled_run``| NA       | NA      | If ``.true.`` the code starts a k-independent calculation. When k-points are activated, this option allows performing an    |
    |                  |          |         | independent self-consistent cycle for each k-point without performing the k-points average of electronic density. To be used|
    |                  |          |         | before a twist average calculation in QMC.                                                                                  |
    +------------------+----------+---------+-----------------------------------------------------------------------------------------------------------------------------+
    | ``yes_kpoints``  | NA       | NA      | Set it to ``.true.`` if you plan to do a calculation with k-points sampling. In this case, the phase of the wavefunction    |
    |                  |          |         | fort.10 is disregarded.                                                                                                     |
    +------------------+----------+---------+-----------------------------------------------------------------------------------------------------------------------------+
    | ``epsbas``       | NA       | NA      | Real space cutoff for the periodic basis set (keyword PBC_C in the first record of the wave function fort.10). Note that, if|
    |                  |          |         | DFT energy is different from what is expected, then try decreasing the cutoff.                                              |
    +------------------+----------+---------+-----------------------------------------------------------------------------------------------------------------------------+

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Molecule section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. table:: Molecule section

   +----------------+----------+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | parameter name | datatype | default | description                                                                                                                                                                                                                                                                                              |
   +================+==========+=========+==========================================================================================================================================================================================================================================================================================================+
   | nx             | NA       | NA      | Number of lattice points for the real space integration grid in the :math:`x` direction. By default, ``nx = ny = nz``                                                                                                                                                                                    |
   +----------------+----------+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ax             | NA       | NA      | Space (a.u.) for an open system in the :math:`x` direction; for periodic systems, it is chosen as the cell parameter in the same direction and need not be specified. By default, ``ax = ay = az``                                                                                                       |
   +----------------+----------+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | nbufd          | NA       | NA      | Input value for the buffer dimension. Note that, in the complex code, the buffer dimension is automatically doubled. In this case, consider decreasing the buffer dimension if you have a memory problem.                                                                                                |
   +----------------+----------+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Kpoints section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. table:: Kpoints section

   +-------------------+----------+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | parameter name    | datatype | default | description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
   +===================+==========+=========+============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================+
   | kp_type           | NA       | NA      | This integer specifies the type of k-points which will be chosen in the calculation.                                                                                                                                                                                                                                                                                                                                                                                                                                       |
   |                   |          |         |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
   |                   |          |         | - ``kp_type = 0``: Do not perform any k-points sampling and use the phase specified in the ``fort.10`` as the unique k-point.                                                                                                                                                                                                                                                                                                                                                                                              |
   |                   |          |         |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
   |                   |          |         | - ``kp_type = 1``: Use the Monkhorst-Pack algorithm to generate equally-spaced k-points in the first Brillouin zone. The size of the grid in the three Cartesian directions is determined by the integers ``nk1, nk2, nk3``. Note that ``nk1`` must be set to a value > 0. If ``nk2, nk3`` are not set, then they are taken to be equal to ``nk1``. Also, if skip_equivalence (see below) is set to ``.false.`` the number of k-points might be reduced.                                                                   |
   |                   |          |         |                    In this case, run the tool ``find_kpoints.x`` with the desired input in order to know how many processors must be allocated.                                                                                                                                                                                                                                                                                                                                                                            |
   |                   |          |         |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
   |                   |          |         | - ``kp_type = 2``: When this is used, k-points are set by the user and their number is specified by the integer ``nk1``. In this case, the section **KPOINTS** is needed (see below).                                                                                                                                                                                                                                                                                                                                      |
   |                   |          |         |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
   |                   |          |         | - ``kp_type = 3``: Generates k-points path along high-symmetry lines in the first Brillouin zone. The initial and final points of these segments are specified in the section **KPOINTS** (see below). The number of extremal points is specified by the integer ``nk1`` and the number of points in each segment is specified by the integer ``nk2``.                                                                                                                                                                     |
   |                   |          |         |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
   |                   |          |         | - ``kp_type = 4``: Generate k-points randomly within the first Brillouin zone. The number of k-points is specified by ``nk1``. **KPOINTS** section is not needed.                                                                                                                                                                                                                                                                                                                                                          |
   +-------------------+----------+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | nk1, nk2, nk3     | NA       | NA      | Meanings depend on the value of ``kp_type``, see above for a detailed explanation.                                                                                                                                                                                                                                                                                                                                                                                                                                         |
   +-------------------+----------+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | k1, k2, k3        | NA       | NA      | For ``kp_type = 1``, set ``k1, k2, k3`` equal to 1 in order to apply an offset to the k-point grid generated by the Monkhorst-Pack algorithm. In some cases, this can help to reach                                                                                                                                                                                                                                                                                                                                        |
   +-------------------+----------+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Additional information:

- Monkhorst-Pack mesh:

  .. code-block:: bash

    &kpoints
    kp_type=1
    nk1=4
    nk2=4
    nk3=4
    skip_equivalence=.false.
    double_kpgrid=.true.

- User-defined k-points are written in the following manner.
  ``wkp(i)`` denotes the weight corresponding to the kpoint
  ``xkp(:,i)`` if the total weight is different from one.:

  .. code-block:: bash

    wkp(i) is the weight corresponding to the the kpoint xkp(:,i).
    ! NB: if the total weight is different from on
    ! xkp(1,1) xkp(2,1) xkp(3,1) wkp(1)
    ! xkp(1,2) xkp(2,2) xkp(3,2) wkp(2)
    ! ......
    KPOINTS
    0.1667 0.1667 0.5000  0.5
    0.5000 0.5000 0.5000  0.5
    # blank line and after k-points for spin down electrons
    -0.1667 -0.1667 -0.5000  0.5
    -0.5000 -0.5000 -0.5000  0.5

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
DFT section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. table:: Parameter List

   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | Parameter Name    | Datatype | Default | Description                                                                                                                    |
   +===================+==========+=========+================================================================================================================================+
   | contracted_on     | NA       | NA      | If ``.true.`` it acts on the contracted basis (considerably faster).                                                           |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | maxit             | NA       | NA      | Maximum number of iterations in the self-consistent cycle.                                                                     |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | epsdft            | NA       | NA      | Tolerance in the convergence of total energy.                                                                                  |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | typeopt           | NA       | NA      | - ``typeopt = 0``: Use self consistency method with standard mixing.                                                           |
   |                   |          |         | - ``typeopt = 2``: Linear mixing scheme.                                                                                       |
   |                   |          |         | - ``typeopt = 3``: Conjugate gradients method with SR acceleration.                                                            |
   |                   |          |         | - ``typeopt = 4``: Anderson mixing scheme with Jacobian acceleration, no use of mixing is made; this method looks to be the    |
   |                   |          |         |   faster and therefore the preferred among the available ones. For information on the algorithm see doc/tex/parbcs.tex,        |
   |                   |          |         |   Ch. V                                                                                                                        |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | mixing            | NA       | NA      | Choose a small value for better convergence. If even in this way it does not converge, switch on the smearing technique        |
   |                   |          |         | setting ``optocc=1`` (suggested for open shell systems). Alternatively you can change iteration method with ``typeopt=3``      |
   |                   |          |         | (conjugate gradients) which will certainly converge for mixing small enough. In these cases mixing means just the maximum      |
   |                   |          |         | amplitude in the step.                                                                                                         |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | mixingder (abbr)  | NA       | NA      | - Case 1 (``typeopt = 3``): Used to evaluate numerically the first and second derivatives.                                     |
   |                   |          |         | - Case 2 (``typeopt = 4``): Used to be closer to the linear regime for the evaluation of the Jacobian (``mixingder`` << 1).    |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | tfcut             | NA       | NA      | Used only with ``typeopt = 0/2/4``. It is used for preconditioning to improve convergence of small q charge fluctuations.      |
   |                   |          |         | Suggested value of ``tfcut`` :math:`= \frac{1}{{xi_{TF}}^2}`  where :math:`xi` is the Thomas-Fermi length expressed in a.u.    |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | orthodiag         | NA       | NA      | ``.false.`` the Kohn-Sham eigenvectors are not orthogonalized after each Hamiltonian diagonalization.                          |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | maxold            | NA       | NA      | The number of previous iterations to be considered in the numerical evaluation of Jacobian with ``typeopt = 4``.               |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | typedft           | NA       | NA      | - ``typedft = 0``: DFT calculation with Hartree potential only.                                                                |
   |                   |          |         | - ``typedft = 1``: LDA (PZ 1981).                                                                                              |
   |                   |          |         | - ``typedft = 2``: LDA (OB 1994).                                                                                              |
   |                   |          |         | - ``typedft = -1,-2``: Same as the two above, but with the corresponding fit performed by imposing continuity in the           |
   |                   |          |         |   correlation energy at :math:`rs = 1`.                                                                                        |
   |                   |          |         | - ``typedft = 3``: KZK finite volume DFT: should be more accurate for finite volume.                                           |
   |                   |          |         | - ``typedft = -3``: Different fitting procedure, suitable for open systems. Could be used with periodic systems too, but it    |
   |                   |          |         |   is less stable.                                                                                                              |
   |                   |          |         | - ``typedft = 4``: Standard LSDA.                                                                                              |
   |                   |          |         | - ``typedft = -4``: Standard LSDA, but with the corresponding fit performed by imposing continuity in the correlation energy   |
   |                   |          |         |   at :math:`rs = 1`.                                                                                                           |
   |                   |          |         | - ``typedft = 5``: LSDA + KZK (not applied on spin) (``typedft = -5`` similar to ``typedft = -3``).                            |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | weightvh          | NA       | NA      | Weight of the hartree potential. Setting ``weightvh`` :math:`\neq 1` is used just for testing.                                 |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | weightxc          | NA       | NA      | Weight of the exchange energy. When ``weightxc`` :math:`= 0` no exchange and when :math:`weightcorr = 1` standard LDA is used. |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | weightcorr        | NA       | NA      | Weight of the correlation energy, e.g. ``weightcorr`` :math:`= 0` no correlation, ``weightcorr`` :math:`= 1` standard LDA.     |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | optocc            | NA       | NA      | - ``optocc = 0``: Use standard occupation of levels. It works well only for closed shell systems (insulators or special cases).|
   |                   |          |         |   Occupations are read from standard input (see below); the number of occupations read is chosen by the parameter ``nelocc``   |
   |                   |          |         |   (``neloccdo`` for down spin electrons) which must be specified in input. In this case we have ``occupations(1:nelocc) = 2``  |
   |                   |          |         |   for LDA; ``occupations(1:nelocc) = 1`` && ``occupationdo(1:neloccdo) = 1`` for LSDA.                                         |
   |                   |          |         | - ``optocc = 1``: Use a smeared Fermi distribution with a spread given by the parameter ``epsshell`` (see below). In this      |
   |                   |          |         |   case:  :math:`occupations(i) = \exp{ \frac{eig(i)-ef}{epsshell}+1}` where the Fermi energy :math:`ef` is determined by the   |
   |                   |          |         |   constraint, sum(occupations(1:bands)) = no. of electrons via bisection method . For LSDA (``|typedft| = 4, 5``) two Fermi    |
   |                   |          |         |   distributions are introduced for up and down electrons. In the case of k-points sampling, the Fermi energy is determined     |
   |                   |          |         |   by averaging over :math:`ef` computed for each k-point.                                                                      |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | epsshell          | NA       | NA      | Spread of the Fermi distribution used when ``optocc = 1``. The unit is Ha.                                                     |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | memlarge          | NA       | NA      | Optimize speed at the cost of much greater memory requirements. The whole basis set is saved on disk!                          |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | epsover           | NA       | NA      | Minimum tolerance for the lowest eigenvalues of the overlap matrix. If ``epsover`` < 0  no orthogonalization is implemented    |
   |                   |          |         | (faster but less stable).                                                                                                      |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | mincond           | NA       | NA      | Disregard the first ``mincond``-1 direction regardless of the condition number limited by ``epsover``. This is useful to have  |
   |                   |          |         | better cancellation errors as a function e.g. of pressure.                                                                     |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | maxcg             | NA       | NA      | With ``typeopt = 3`` (conjugate gradient) each ``maxcg`` steps restart the conj. grad. procedure. If ``maxcg = 0`` no          |
   |                   |          |         | restarting is performed (discouraged since numerically unstable).                                                              |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | bands             | NA       | NA      | The number of lowest eigenvalues of Khon-Sham equations to be evaluated. By default ``bands = nelup + 7`` where ``nelup``      |
   |                   |          |         | is the number of spin up electrons. This corresponds to assuming at most an 8-fold degenerancy in the last occupied shell.     |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | nelocc            | NA       | NA      | If ``nelocc`` > 0, it is the number of occupations that are read in the last record of the input file. Occupation values can   |
   |                   |          |         | only be (0,2] (paired orbitals), -1 (unpaired at the end) or 0 (unoccupied). If ``nxs`` :math:`\times` ``nys`` :math:`\times`  |
   |                   |          |         | ``nzs`` > 0 this record is just after the ones to read the input magnetization (see below).                                    |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | neloccdo          | NA       | NA      | It is similar to ``nelocc`` but for the spin down electrons which are assumed with no unpaired orbitals. Another record with   |
   |                   |          |         | ``neloccdo`` integers be written below. Note that, in this case occupations for up spin can take values 1 (occupied paired     |
   |                   |          |         | orbital), -1 (unpaired), 0 (unoccupied orbital). Instead occupations for down spin electrons can take values 1 (occupied       |
   |                   |          |         | paired orbital) and 0 (unoccupied orbital).                                                                                    |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | randspin          | NA       | NA      | Used for initializing magnetization. If ``randspin`` > 0, add random component to the orbitals. If ``randspin`` < 0            |
   |                   |          |         | initialize with maximum possible spin given density and grid (see below). If zero no action.                                   |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | jaccond           | NA       | NA      | Minimum threshold for the condition matrix in the self-consistent approach ``typeopt = 4``.                                    |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | nxs               | NA       | NA      | Dimension of the grid where the magnetization is defined along the :math:`x`, :math:`y`, :math:`z` direction. The format is    |
   | nys               |          |         | written below.                                                                                                                 |
   | nzs               |          |         |                                                                                                                                |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | h_field           | NA       | NA      | If ``h\_field`` > (<) 0 put a magnetic field increasing (decreasing) the magnetization with the staggering given by the table  |
   |                   |          |         | sxyz defined for ``randspin``.                                                                                                 |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | optimize_overs    | NA       | NA      | If ``.true.`` optimize the overlap matrices calculation if the phase for spin down electrons is equal or opposite to the       |
   |                   |          |         | phase for up spin electrons. Otherwise it is automatically set to ``.false.``.                                                 |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | write_den         | NA       | NA      | If ``.true.`` write the overlap matrix elements for effective Hamiltonian calculations to the disk.                            |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | zero_jas          | NA       | NA      | If ``.true.`` set the one-body Jastrow to zero at the end of the DFT calculation.                                              |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+
   | fix_density       | NA       | NA      | If the flag ``decoupled_run`` is set to ``.true.`` (&parameters card) as well as the ``yes_kpoints`` flag (&kpoints card),     |
   |                   |          |         | the k-points are evolved independently but using the averaged electronic density.                                              |
   +-------------------+----------+---------+--------------------------------------------------------------------------------------------------------------------------------+

Additional information:

- ``nxs`` ``nys`` ``nzs``

  Dimension of the grid where the magnetization is defined along the :math:`x`, :math:`y`, :math:`z` direction:
  The format is

  .. code-block:: bash

        !  After "/" or "occupation list"
        ! # empty line
        ! s111 s211 s311 ... snxs11
        ! s121 s221 s321 ... snxs11
        ! s1nys1 ... ... ... snxsnys1
        ! # empty line
        ! s112 s212 s312 ... snxs12
        ! ...
        ! # empty line
        ! s11nzs s21nzs s31nzs ... snxs1nzs
        ! ...
        ! s1nysnzs ... ... ... snxsnysnzs

..
  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  Band_structure section
  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  This card should be used if and only if the flag ``compute_bands`` is set to ``.true.`` in the **&simulation** card.
  ``task`` (default: `0`) Flag to specify which quantity to compute after a non self-consistent run. It is ignored if ``type_comp_dft = 0``.
      ``task = 0`` Do not compute anything.
      ``task = 1`` Band structure plot (use ``kp_type = 2/3`` in the k-points card to specify the path in the Brillouin zone).
      ``task = 2`` Density of States calculations using smearing parameter given by ``epsshell``. The integer ``optocc`` must be set to 1.
  ``emin`` min(eigenvalue): minimum value of the energy to be included in band structure or DOS plot.
  ``emax`` max(eigenvalue): maximum value of the energy to be included in band structure or DOS plot.
  ``deltaE`` Energy bin for computing the density of states (task = 2).

