.. TurboRVB_manual documentation master file, created by
   sphinx-quickstart on Thu Jan 24 00:11:17 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Wavefunction
===========================================

TurboRVB employs a many-body WF ansatz :math:`\Psi` which can be written as the product of two terms:

.. math::

    \Psi  =  \Phi _\text{AS} \times \exp J \,,

where the term :math:`\exp J`, conventionally dubbed Jastrow factor, is symmetric under electron exchange, and the term :math:`\Phi _\text{AS}`, also referred to as the determinant part of the WF, is antisymmetric.
The resulting WF :math:`\Psi` is antisymmetric, thus fermionic.

The Jastrow factor (:math:`\exp J`) plays an important role
in improving the correlation of the WF and in fulfilling Kato's cusp conditions.
TurboRVB implements the Jastrow term composed of one-body, two-body, and three/four-body factors (:math:`J = {J_1}+{J_2}+{J_{3/4}}`).

-------------------------------------------
Wavefunction file format (fort.10)
-------------------------------------------

In TurboRVB, ``fort.10`` is a fundamental file that contains all information about the nuclear positions and the wave function (WF) details and parameters. ::

    # fort.10 of the C2-dimer (the Pfaffian ansatz with the Filippi pseudo potential.)
    # Nelup  #Nel  # Ion
            4          -8           2
    # Shell Det.   # Shell Jas.
            50          43
    # Jas 2body  # Det   #  3 body atomic par.
            -22        1482          42
    # Det mat. =/0  # Jas mat. =/0
            120        8370
    # Eq. Det atomic par.  # Eq. 3 body atomic. par.
            741          21
    # unconstrained iesfree,iessw,ieskinr,I/O flag
            8370         120           6           0
    # Ion coordinates
    4.00000000000000        6.00000000000000       0.000000000000000E+000
    0.000000000000000E+000  -1.14999954166875
    4.00000000000000        6.00000000000000       0.000000000000000E+000
    0.000000000000000E+000   1.14999954166875
    #  Constraints for forces: ion - coordinate
            1           1           1
            1           1           2
            1           1           3
            1           2           1
            1           2           2
            1           2           3
    #          Parameters Jastrow two body
            -1  0.342214663461764
    ...

In the following each part of the file will be described. The file is divided in different sections and it contains numbers (in a
free format) as well as comment lines (which start with \#).

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Header
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This part accounts for the first twelve lines of the file.
It contains general information about the system (such as number
of electrons and nuclei), wavefunction general characteristic (such as
the number of parameters and symmetry information between them).
Here we report an header example and then we shortly describe the
keywords.::

     # Nelup  #Nel  # Ion
               2           4           1
     # Shell Det.   # Shell Jas.
               3           3
     # Jas 2body  # Det   #  3 body atomic par.
              -8          16           8
     # Det mat. =/0  # Jas mat. =/0
               6           6
     # Eq. Det atomic par.  # Eq. 3 body atomic. par.
              13           8
     # unconstrained iesfree,iessw,ieskinr,I/O flag
               4           4           0           0

``NELUP`` Number of spin up electrons in the system.

``NEL`` Total number of electrons in the system.

``ION`` Number of nuclei in the system.

``SHELL DET`` Number of shells  used to describe the determinantal term of the AGP wavefunction.

``SHELL JAS`` Number of shells used to describe the Jastrow term of the wavefunction.

``JAS 2BODY`` Type of the Jastrow 2B/1B term used to satisfy the electron-electron, electron-ion cusp conditions, respectively.

The one-body Jastrow factor :math:`J_1` is the sum of two parts, the homogeneous part (enforcing the electron-ion cusp condition) and the corresponding inhomogeneous parts.

The homogenenous part reads

.. math::

    J_1^h \left( \mathbf{r}_1,\ldots,\mathbf{r}_N \right) = \sum_{i=1}^N \sum_{a=1}^{N_\text{at}} \left( { { - {{\left( {2{Z_a}} \right)}^{3/4}}u_a\left( {(2{Z_a})^{1/4}\left| {{\mathbf{r}_i} - {{\mathbf{R}}_a}} \right|} \right)} } \right),

where the most common choice for :math:`u_a` in turborvb is:

.. math::

    u_a\left( r \right) = \frac{ 1 }{2 b_{\text{e}a}} \left( {1 - {e^{ - r b_{\text{e}a}}}} \right)

depending on a single variational parameter :math:`b_{\text{e}a}`, that may be optimized
independently for each atomic species, but we can choose several different forms as described below.

The two-body Jastrow factor is defined as:

.. math::

    {J_2}\left( {{{\mathbf{r}}_1}{\sigma _1}, \ldots, {{\mathbf{r}}_N}{\sigma _N}} \right) =  {\sum\limits_{i < j} {{v_{{\sigma _i},{\sigma _j}}}\left( {\left| {{{\mathbf{r}}_i} - {{\mathbf{r}}_j}} \right|} \right)} },

where :math:`v_{{\sigma _i},{\sigma _j}}` is another  simple bounded  function. There are several possible choices for :math:`v_{{\sigma _i},{\sigma _j}}` implemented in turborvb (all listed as below), and one of them is, for instance, the following spin-dependent form:

.. math::

    {v_{{\sigma _i},{\sigma _j}}}\left( {{r_{i,j}}} \right) =
    \begin{cases}
        \cfrac{{{r_{i,j}}}}{4} \cdot {\left( {1 + b_{\rm{ee}}^{\rm{para}} \cdot {{r_{i,j}}}} \right)^{ - 1}} & ({\sigma _i} = {\sigma _j}) \\
        \cfrac{{{r_{i,j}}}}{2} \cdot {\left( {1 + b_{\rm{ee}}^{\rm{anti}} \cdot {{r_{i,j}}}} \right)^{ - 1}} & ({\sigma _i} \neq {\sigma _j})
    \end{cases},

where :math:`{r_{i,j}} = \left| {{{\mathbf{r}}_i} - {{\mathbf{r}}_j}} \right|`, and :math:`b_{\rm{ee}}^{\rm{para}}` and :math:`b_{\rm{ee}}^{\rm{anti}}` are variational parameters.

Different functional form of the Jastrow one and two body terms are available and each one is identified with an integer number code iesdrr (the leftmost integer below the line::

      # Jas 2body  # Det   #  3 body atomic par.

, i.e. -8 in the above example)  that selects also the number of parameters :math:`p` used for the two body Jastrow part only. The input consists of one line below::

      #          Parameters Jastrow two body
      e.g.  2  1.0 1.0

There are two body Jastrow defined with a number :math:`p` of parameters larger than one, that are put in order from left to right in the input line according to the alphabetic order, :math:`a,b,c \cdots` . The absolute value of the integer in the record  indicates the number of parameters :math:`niesd=p+p_{obebody}`, where :math:`p_{onebody}` is the number of parameters used in the one body Jastrow. From left to right in the record the first :math:`p` parameters refer to the two-body Jastwow and the remaining ones to the one body Jastrow. The one body form is assumed to have at most one parameter for each atomic
species and its default form is given by a rescaled (-4, see below) in order to satisfy the electron-ion cusp.There are only three possible values allowed for :math:`niesd`:

	:math:`p_{onebody} = 0` the one body parameter is set equal to the last parameter of the record for all the atomic species.

	:math:`p_{onebody} = 1` the one body parameter is set equal to the last parameter of the record for all the atomic species, but is independent of the two body parameters.

	:math:`p_{onebody} = \#` different atoms as above, the one body parameters act independently on each atomic species, for instance in water there are two independent atomic species (two Hydrogen and one Oxygen) and the one-body is defined by two parameters. In this case the one body parameters are sorted in the record according to the atomic number, the leftmost corresponding to the lightest atomic number.

Be careful with this number, as its allowed value is not tested in the input (yet).
The following values of the two body Jastrow (iesdrr) are allowed:

.. table:: Jastrow ``iesdrr`` summary

   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | iesdrr  | description                                                                                                                                                                                                                                                                                                                                      |
   +=========+==================================================================================================================================================================================================================================================================================================================================================+
   | 0       | No two body and one body, 3B Jastrow may be on.                                                                                                                                                                                                                                                                                                  |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -5      | Two body :math:`\frac{r}{2(1+ar)}` one body rescaled same form.                                                                                                                                                                                                                                                                                  |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -1      | Two body :math:`\frac{r}{2(1+ar)}` for opposite spins and  :math:`\frac{r}{4(1+ar)}` for parallel spins, one body rescaled same form.                                                                                                                                                                                                            |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -4      | Two body :math:`\frac{1}{2a} (1 - e^{-ar})` one body rescaled. Not spin contaminated.                                                                                                                                                                                                                                                            |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -7, -6  | Two body :math:`\frac{1}{2a} (1 - e^{-ar})` one body rescaled, +cusp for parallel spins (divided by two).                                                                                                                                                                                                                                        |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | 2       | Two parameters Jastrow improved version of -1 with an independent parameter for the parellel spins, :math:`\frac{r}{4(1 + br)}` for (anti-)parallel spins, spin contaminated.                                                                                                                                                                    |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -8      | Two body :math:`\frac{1}{2a} (1 - e^{-ar})` one body rescaled, + cusp for (anti-)parallel spins + 3B Jastrow Sz.                                                                                                                                                                                                                                 |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | 8       | Two body :math:`\frac{1}{a} (1 - e^{-ar^3})` for pseudo soft.                                                                                                                                                                                                                                                                                    |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -9      | Two body :math:`-A \ln(1 + a{(1 - \frac{r}{b})}^2)` for RVB wavefunction, with :math:`A = \frac{b(1 + a)}{4a}` , to satisfy the cusp conditions for opposite spin electrons. Two parameters :math:`niesd \geq 2`.                                                                                                                                |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | 9       | Two parameters RVB two-body Jastrow. Two body :math:`-A \ln(1 + a{(1 - \frac{r}{b})}^2)` for RVB wavefunction, with :math:`A = \frac{r_0(1 + b)}{4b}` , to satisfy the cusp conditions for opposite spin electrons.                                                                                                                              |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | 10      | Two parameters RVB two-body Jastrow. Two body :math:`-A \ln(1 + a{(1 - \frac{r}{b})}^2)` for RVB wavefunction, with :math:`A = \frac{r_0(1 + b)}{4b}` , to satisfy the cusp conditions for opposite spin electrons. Rescaled :math:`r \to \frac{r}{2}` to satisfy the cusp condition for parallel spins.                                         |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | 5       | Three parameters :math:`\frac{1}{a + b*c} (1+c-\exp(-ar)-c\exp(-br))` improved version of -6. Warning! Implemented only for open systems.                                                                                                                                                                                                        |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | 6       | Two parameters, the second is used to rescale the electron-electron distance :math:`r_s = \frac{1-\exp(-br)}{r}` and the Jastrow is defined by :math:`J_{ee}=\frac{r_s}{2(1+ar_s)}` , no spin contamination and cusp condition for opposite spin electrons.                                                                                      |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -2      | Two parameter Jastrow :math:`r_z = \sqrt{a^2(x^2+y^2) + {(bz)}^2}` , and :math:`J_2 = \frac{1}{2} \frac{r}{1+r_z}` + cusp for (anti-)parallel spins for anisotropic phases. Warning! Implemented only for open systems.                                                                                                                          |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | 3       | Three parameters correction to :math:`-5` :math:`J_2 = \frac{r}{2}(\frac{1}{1+ar} + \frac{cr}{{(1+br)}^2})` + cusp for (anti-)parallel spins.                                                                                                                                                                                                    |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -9      | Two body :math:`\frac{1}{2b} (1 - e^{-br})` one body rescaled, + cusp for (anti-)parallel spins + 3B + 1B Jastrow Sz (for studying magnetic phases).                                                                                                                                                                                             |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -10     | No two body and one body, 3B Jastrow and Jastrow Sz is on.                                                                                                                                                                                                                                                                                       |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -11     | No two body and one body, 3B+1B Jastrow and Jastrow Sz are on (for studying magnetic phases).                                                                                                                                                                                                                                                    |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -12     | General spin-density Jastrow, one body and two body as -15, namely without spin dependent cusp condition.                                                                                                                                                                                                                                        |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -15     | Long range two body :math:`\frac{r}{2(1+br)}` ; short range one body :math:`\frac{1}{2b} (1-e^{-br})` .                                                                                                                                                                                                                                          |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -20     | Two parameters, spin dependent (as -7) long range two body :math:`\frac{r}{2(1+ar)}` ; short range one body :math:`\frac{1}{2b} (1-e^{-br})` .                                                                                                                                                                                                   |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -21     | Three parameters, first two same as Jastrow number 2; short range one body :math:`\frac{1}{2c} (1-e^{-cr})`.                                                                                                                                                                                                                                     |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -22     | General spin-density Jastrow one body and two body as -20, with spin dependent cusp condition, more appropriate in this case, as the spin contamination is already implied by the three and four body term.                                                                                                                                      |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -26     | General spin-density Jastrow one body and two body as -7, with spin dependent cusp condition, without long range power law tails.                                                                                                                                                                                                                |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -27     | General spin-density Jastrow one body and two body as -21, with spin dependent cusp condition.NB :math:`p=2` in this case, so one can put niesd=3 safely. Warning! If you put niesd>2, it is recommended to set niesd equal to 2 + # different atomic species, e.g. niesd=4 for benzene. In this way,                                            |
   |         | all different atomic species will have a different one-body term.                                                                                                                                                                                                                                                                                |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -30     | General spin-density Jastrow one body and two as 10, with spin dependent cusp condition. NB :math:`p=2` in this case, so one can put niesd=3 safely.                                                                                                                                                                                             |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -31     | General spin-density Jastrow one body and two body as 10, with spin dependent cusp condition. NB :math:`p=2` in this case, so one can put niesd=3 safely.                                                                                                                                                                                        |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -18     | same as :math:`iesdrr=-8` but with two body :math:`\frac{r}{2(1+br)}` .                                                                                                                                                                                                                                                                          |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -19     | same as :math:`iesdrr=-9` but with two body :math:`\frac{r}{2(1+br)}` as -7.                                                                                                                                                                                                                                                                     |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -28     | same as :math:`iesdrr=-8` but with two body/one body as -20.                                                                                                                                                                                                                                                                                     |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -29     | same as :math:`iesdrr=-9` but with two body/one body as -20.                                                                                                                                                                                                                                                                                     |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | -16     | same as :math:`iesdrr=-19` but with spin independent two body as -5.                                                                                                                                                                                                                                                                             |
   +---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


``DET`` Total number of atomic variational parameters used for the description of the determinant part of the JAGP. The atomic variational parameters are all parameters describing the shell used to expand the WF (i.e. gaussian coefficient and gaussian exponents). For example, if a shell is described by a single Gaussian function, :math:`\phi(r) \sim \exp(-Zr^2)` there is only one atomic parameter, the Z exponent of the Gaussian. If the shell is described by a contraacted orbital of two Gaussians, :math:`\psi = \alpha_1 \exp(-Z_1 r^2) + \alpha_2 \exp(-Z_2 r^2)` then there are four atomic parameters, :math:`Z_1` and :math:`Z_2` and the coefficients :math:`\alpha_1` and :math:`\alpha_2` of the linear combination.

``3 BODY ATOMIC PAR`` Total number of atomic parameters used to describe the Jastrow part of the JAGP. It is analogous to **DET** but for the Jastrow orbitals.

``DET MAT =/0`` Number of coefficients :math:`\{\lambda_{ij}\}` of the determinant that are different from zero. Note that the :math:`\Lambda` matrix is symmetric and only :math:`\lambda_{ij}` for :math:`i \geq j` are provided in the data file. This number corresponds to the number of lines read afterwards.

``JAS MAT =/0`` Number of coefficients :math:`\{\lambda_{ij}\}` in the 3B Jastrow that will be considered different from zero. This number corresponds to the number of lines that will be read afterwards.

``EQ DET ATOMIC PAR`` Number of symmetries involving the variational atomic parameters of the determinant. This number corresponds to the number of lines that will be read afterwards.

``EQ 3 BODY ATOMIC PAR`` Number of symmetries involving the variational atomic parameters of the 3B Jastrow. This number corresponds to the number of lines that will be read afterwards.

``UNCONSTRAINED IESFREE`` Number of independent coefficients :math:`\{\lambda\}` of the Jastrow. This number corresponds to the number of lines to be read when the symmetries for the Jastrow :math:`\{\Lambda\}` matrix are described.

``IESWW`` Number of independent coefficients :math:`\{\lambda\}` of the determinant. This number corresponds to the number of lines to be read when the symmetries for the determinant are described.

``IESKINR`` Number of cartesian nuclear components to be optimized. This number corresponds to the number of lines to be read in the relative section of the file.

``I\O FLAG`` It describes how to read the FORT.10 file. If '0' the file is read and all the information is used. If '1', at the end of the file, an extra part is expected where new parameters (for example averaged from previous simulations) are provided in free format.

In case the system has periodic boundary conditions (PBC), two additional lines appear as first lines at the beginning of the header. Here is an example::

   	#   PBC rs, Ly/Lx, Lz/Lx
	    1.3100  	   1.0000	1.0000		0.5	0.0	0.0

The first line is a comment line required to switch on the use of PBC and the second line lists the cell dimension in :math:`x` direction :math:`Lx`, the ratio between :math:`Ly` and :math:`Lx` and the ratio between :math:`Lz` and :math:`Lx` . The last three numbers correspond to the phase of the wave-function along the direction :math:`x` , :math:`y` , :math:`z` . Zero is used for a periodic wavefunction and 0.5 for an antiperiodic along a given direction.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
After the header the coordinates of the nuclei are provided
in the same line. Starting from the left the first number is the
number of valence electron in the atom (atomic number :math:`Z` - number of core electrons considered in the pseudopotential), whereas the second number
is the atomic number :math:`Z` of the
atom (:math:`N \ne Z` with pseudopotential calculation).
The data are free format.
The coordinates are in atomic units (BOHR). For example for :math:`n` nuclei::

     # Ion coordinates
      N1 Z1                x1     y1     z1
      N2 Z2                x2     y2     z2
        ..                ..     ..     ..
      Nn Zn                xn     yn     zn

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Ionic Forces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This part of FORT.10 lists the cartesian components of the nuclear forces that will be calculated and used for the structural optimazation or dynamics. The number of lines to be read is defined by the **N.FORCES** in the **HEADER**. If **N.FORCES=0** in the header, no line will be read. At the same time, it can eventually specify symmetries to be enforced on the nuclear coordinates. To identify a force component, two numbers have to be specified: the atom number (according to the ion coordinate list) and the cartesian component (1 for X, 2 for Y and 3 for Z). For example::

     # Constraints for forces: ion - coordinate
       		1      1      3

The first number specifies that there is only one cartesian component in this line. The component is therefore independent of others (no symmetry). The component corresponds to atom 1, :math:`z` (i.e. 3) coordinate. In the following example, symmetry is specified::

    # Constraints for forces: ion - coordinate
               2      1	     1      2      -3

The first number indicates that two components have to be read afterwards, forming a symmetry constraint. For each component as usual, two numbers are expected: the ion index followed by the kind of component (x, y or z). In the above case, the :math:`x` coordinate of nucleus number 1 and the :math:`z` coordinate of nucleus number 2 are set to have opposite values because the coordinate index for nucleus number 2 is negative (-3). If negative sign is not used, the two components would be set to be equal, i.e. with the following simpler input::

    # Constraints for forces: ion - coordinate
               2      1	     1      2      3

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The 2B Jastrow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameter(s) of the 2B Jastrow are listed in one line::

    #	   Parameters Jastrow two body
           1  0.549835086466315

The absolute value of the first number, dubbed as :math:`iesd` in the code, indicates how many parameters have to be read (they depend on the Jastrow type) for iesdrr different from zero. The subsequent numbers in the same record are the 2B Jastrow parameters. If the first integer is negative the AGP function is not assumed to be symmetric. If no one-two body Jastrow is used (iesdrr=0) the records::

    #		Parameters Jastrow two body
    		-1

means AGP is not symmetric whereas::

    #		Parameters Jastrow two body
    		0

would be the standard symmetric case (not spin contaminated). In other cases (iesdrr not equal to zero) the sign of the first integer number determines the AGP symmetry as before, and its absolute value determines the following possibilities:

      * If ``iesd = 1`` the one body and two body Jastrow share the same variational parameter.

      * If ``iesd = 2`` there are two independent variational parameters one for the one-body Jastrow and one for the two-body Jastrow.

      * If ``iesd > 2``, ``iesd`` should be equal to the number different atomic species in the system plus one (e.g. in water :math:`iesd = 3` because of two atomic species corresponding to H and O), because for each atomic species, we assume an independent variational parameter for the one-body Jastrow. The variational parameters are ordered from left to right in this record, in the order of increasing atomic number (e.g. in the water for example, the first one corresponds to the two-body term, the second to the Hydrogen one body parameter, and the third to the Oxygen one).


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Basis Set for Determinant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this section the Basis Set used for expanding the JAGP determinant is described. There are a fixed number of lines to be read (2*NSHELL). Each shell of the determinant is described by two lines. The first one contains the multiplicity, the number of variational parameters of the shell function and the code describing the function. The code numbers and the description of the corresponding shell are described in the file makefun.f90 of the source code. The multiplicity depends on the shell type: Shells S, P and D have the multiplicities of 1, 3 and 5 respectively. In the second line the index of the nucleus on which the shell is centered is first indicated. Then the parameter values are listed. Keep in mind that the number of parameters to be read is given in the first line.::

   #		Parameters atomic wf
   Shell_Multiplicity	   Number of par.		Shell code
   Ion index		   [par (1, NUMBER OF PAR.)]

   #   		Parameters atomic wf
   		1           1          16
		1  0.500000000000000
		3           1          36
		1   1.00000000000000
		1           1          16
		2  0.300000000000000
		1           1          16
		3  0.300000000000000
		1           1          16
		4  0.300000000000000
		1           1          16
		5  0.300000000000000

All primitive orbitals are written in the source file makefun.f90 (open boundary), makefun_pbc.f90 (pbc) and makefun_bump.f90 (finite range orbitals).

TurboRVB also implements standard contracted orbitals written  as a linear combination of :math:`p` primitive orbitals. The definitions are easily found (and can be easily implemented) in the fortran file: ioptorbcontr.f90. In this case, the number corresponding to "Number of par." is equal to :math:`2p`. In the next line, one writes these extra coefficients, :math:`C_i, i = 1,...2p:` the coefficient :math:`C_{i+p}` acts on the orbital number defined by the contracted orbital written in "Shell code", with exponent :math:`Z_i = C_i` (we omit the normalization, each orbital is assumed to be normalized), for instance a :math:`2s` contracted orbital:

.. math::

   \phi(r) = 3.231 \cdot \exp(-2.0 \cdot r^2) + 7.54 \cdot \exp(-1.0 \cdot r^2)

is written as::

   #	   Parameters atomic wf
        1	      4		300
        1	2.0   1.0  3.231  7.54

Shell code::

    16 -> s orbital
    36 -> p orbital
    68 -> d orbital
    48 -> f orbital
    51 -> g orbital
    72 -> h orbital
    73 -> i orbital

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Molecular orbital
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``convertfort10mol.x`` can add molecular orbitals to ``fort.10``.

.. math::

    \tilde \Phi_k=\sum_{i=1}^{N_{\rm basis}}c_{i,k} \cdot \phi_i \left({\bm r} \right)

In ``fort.10``, ``1000000`` indicates a molecular orbital.::

           #always 1,   the number of components, 100000
           #index of basis [1,2,....]
           #coefficients for basis [1,2,....]
            1         180     1000000
            1           1           2           3           4           5
            6           7           8           9          10          11
            12          13          14          15          16          17
            18          19          20          21          22          23
            24          25          26          27          28          29
            30          31          32          33          34          35
            36          37          38          39          40          41
            42          43          44          45          46          47
            48          49          50          51          52          53
            54          55          56          57          58          59
            60          61          62          63          64          65
            66          67          68          69          70          71
            72          73          74          75          76          77
            78          79          80          81          82          83
            84          85          86          87          88          89
            90  0.438271164894104      -4.608166217803955E-002
    0.189550578594208       7.299757003784180E-002 -0.129178702831268
    -0.241831779479980      -7.793867588043213E-002 -0.143670558929443
    -0.181271851062775      -0.265352427959442       0.374841809272766
    -5.072158575057983E-002 -0.286649286746979       0.421764492988586
    -0.147124171257019       0.281676769256592       0.136297583580017
    -6.065595149993896E-002 -0.442295849323273       6.872278451919556E-002
    -0.382538557052612       0.445110499858856       0.338095664978027
    -6.700992584228516E-002 -0.306204080581665      -0.206682741641998
    -0.368851244449615      -0.185477197170258       0.275709867477417
    4.901111125946045E-003 -0.484303355216980      -0.346608221530914
    0.117953181266785      -9.146583080291748E-002 -0.450106739997864
    -0.420913279056549      -5.812942981719971E-002  6.194984912872314E-002
    -0.185146868228912       9.111911058425903E-002 -0.494102835655212
    -0.187083423137665       0.336221218109131       0.465211153030396
    -0.328829050064087       0.109235763549805      -0.194452404975891
    0.369445860385895       0.259138166904449       0.417137503623962
    0.273341476917267       0.424752712249756       0.248025238513947
    -0.142549693584442       0.235680162906647       0.194104492664337
    -0.394946813583374       0.418918550014496       0.286247551441193
    -0.165200054645538      -0.198730885982513      -0.141955792903900
    0.255677580833435       0.398207664489746      -0.194129884243011
    -0.355173707008362      -0.489923000335693      -0.488654732704163
    0.106894016265869       6.377434730529785E-002 -0.169396936893463
    -7.532972097396851E-002  5.515247583389282E-002  0.329569637775421
    4.597079753875732E-002  0.160456597805023      -0.421718478202820
    0.326750636100769       0.339765250682831       0.246131539344788
    6.178361177444458E-002  0.332796216011047      -9.556728601455688E-002
    0.266949594020844      -5.606234073638916E-002 -0.166017174720764
    0.363827764987946       0.222376465797424       0.321450889110565
    0.293389737606049

Note that the paring function is represented by molecular orbitals as follows:

.. math::

    f\left( {{{\mathbf{r}}_i},{{\mathbf{r}}_j}} \right) = \sum\limits_{k=1}^{M} {{{\lambda}_{k}}{\tilde \Phi _{k}}\left( {{{\mathbf{r}}_i}} \right){\tilde \Phi _{k}}\left( {{{\mathbf{r}}_j}} \right)}

We recover the standard Slater determinant by using :math:`M = N_{\rm electron}/2`.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Hybrid orbital
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``convertfort10.x`` can add hybrid orbitals to ``fort.10``.

.. math::

    \tilde{\Phi}_{k}^{a}=\sum_{i}^{N_{basis}^{a}}c_{i,k}^{a} \cdot \phi_{i}^{a} \left({\bm r} \right)

In ``fort.10``, ``900000`` indicates a hybrid orbital.::

           #atom index,   the number of components, 900000
           #index of basis [1,2,....]
           #coefficients for basis [1,2,....]
           1          90      900000
           1           1           2           3           4           5
           6           7           8           9          10          11
          12          13          14          15          16          17
          18          19          20          21          22          23
          24          25          26          27          28          29
          30          31          32          33          34          35
          36          37          38          39          40          41
          42          43          44          45  4.415629124781804E-004
 -2.665471779107012E-003  6.462432564128138E-003 -2.835419348050690E-002
 -1.801631060924397E-003  0.000000000000000E+000  0.488607316300470
  -1.08438499459258        1.00000000000000      -8.862166918962697E-002
 -3.585186058053676E-006  1.336604741376926E-005 -1.034531513737405E-002
 -9.682848944243672E-005  7.780031865084339E-005 -2.330742724611005E-002
  2.406091695918141E-004 -1.314165271323025E-004 -1.383980745357883E-002
 -8.265782808695380E-005  6.484551172622954E-005 -6.461899158835088E-002
 -1.224707476287876E-003  3.784895470189310E-004 -0.150215939301643
  3.279641185411990E-003  4.146447536354074E-004 -5.370349262603082E-002
 -2.867011651206949E-003 -6.433419048360544E-004 -0.169339132045330
  0.000000000000000E+000  0.000000000000000E+000  5.388221636070670E-002
  3.865987502185416E-004  1.537139475756473E-004 -0.134079774055112
 -4.310821717930281E-004 -2.609351310060059E-004  2.739963165919702E-002
  2.693581696773797E-002 -3.679623706959143E-005  1.393726496455493E-004
 -1.728876043016380E-004  1.648834534792205E-004

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Basis Set for 3B Jastrow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section describes the Basis Set used to expand the 3B Jastrow. The data are organized as in the section **Basis Set for Determinant**. The code number for the shells used for expanding the 3B Jastrow is also described in makefun.f90.::

     #	     Parameters atomic Jastrow wf
     	 1		0      200
	     1
	     1		1	1000
	     1  0.993536719652206

Similar to the determinantal case, one can insert contracted orbitals in the Jastrow basis (see ioptorbcontr.f90 for the definition/implementation). The only difference is that, in this case, the orbitals are not normalized. Also, when working with periodic boundary conditions, they are periodized always without any twist, as it should. The :math:`2s` contraction defined above reads in the Jastrow basis section::

	#      Parameters atomic Jastrow wf
	       1           4          3000
           1      2.0  1.0 3.231  7.54


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Occupation Determinant Orbitals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This part provides the occupation state of the determinant orbitals. The number of lines is :math:`\sum_{i}^{NSHELL} shell\_multiplicity(i)`. If occupied the orbital takes value of one, and otherwise zero. The orbitals are numbered as in the **Basis Set for Determinant**. Keep in mind that a shell P counts for three orbitals  and a shell D counts five orbitals::

     #	  Occupation atomic orbitals
            1
            1
            1
            1
            1
            1
            1
            1


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Occupation Jastrow Orbitals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This part provides the occupation state of the Jastrow orbitals. See above::

     #	  Occupation atomic orbitals Jastrow
            1
            1


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Coefficient of the Determinant :math:`A` Matrix different from zero
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For simplicity, we consider a system with an even number :math:`N` of electrons here.
The WF, written in terms of pairing functions, is

.. math::

    \Phi_\text{AS} (\mathbf{1},\ldots,\mathbf{N}) = {\cal A} \left\{ g(\mathbf{1},\mathbf{2}) g(\mathbf{3},\mathbf{4}) \ldots g(\mathbf{N-1},\mathbf{N}) \right\},

where :math:`{\cal A}` is the antisymmetrization operator

.. math::

    {\cal A} \equiv {1\over N!} \sum_{P\in S_N} \epsilon_P \hat P,

:math:`S_N` the permutation group of :math:`N` elements, :math:`\hat P` the operator corresponding to the generic permutation :math:`P`, and :math:`\epsilon_P` its sign.


Let us define :math:`G` the :math:`N\times N` matrix with elements :math:`G_{i,j} = g(\mathbf{i},\mathbf{j})`.Notice that

.. math::

    g(\mathbf{i},\mathbf{j}) = -g(\mathbf{j},\mathbf{i}) \; (\text{and} \; G_{i,j} = -G_{j,i}),

as a consequence of the statistics of fermionic particles, thus :math:`G` is skew-symmetric ({\it i.e.}, :math:`G^T = -G`, being :math:`^T` the transpose operator), so the diagonal is zero and the number of independent entries is :math:`\sum_{n=1}^{N-1} n = N(N-1)/2`.

The most general representation of a many-body wavefunction in TurboRVB is
the Pfaffian of a :math:`N\times N` skew-symmetrix matrix :math:`G` is defined as

.. math::

    \Phi_\text{Pf} = \text{Pf}(G) \equiv {1\over 2^{N/2} (N/2)!} \sum_{P\in S_N} \epsilon_P G_{P(1),P(2)} \cdots G_{P(N-1),P(N)}

Notice that the :math:`\Phi_\text{Pf}` here defined allows the description of  any system with :math:`N_u` electrons with spin-up and :math:`N_d` electrons with spin-down, provided that :math:`N=N_u+N_d` is even. Indeed, with no loss of generality, we can assume that electrons :math:`i=1,\ldots,N_u` have :math:`\sigma_i=1/2` and electrons with :math:`i=N_u+1,\ldots,N` have :math:`\sigma_i=-1/2`.
Thus, the :math:`N\times N` skew-symmetric matrix :math:`G` is written as:

.. math::

    G = \left[\begin{array}{c|c} G_{uu} & G_{ud} \\ \hline
    G_{du} & G_{dd}\end{array}\right]

where,

- :math:`G_{uu}` is a :math:`N_u\times N_u` skew-symmetric matrix with elements :math:`g_{uu}(\mathbf{i},\mathbf{j})`,
- :math:`G_{dd}` is a :math:`N_d\times N_d` skew-symmetric matrix with elements :math:`g_{dd}(\mathbf{i},\mathbf{j})`,
- :math:`G_{ud}` is a :math:`N_u\times N_d` matrix with elements :math:`g_{ud}(\mathbf{i},\mathbf{j})`, and
- :math:`G_{du} = -{G_{ud}}^T`, i.e., :math:`g_{du}(\mathbf{i},\mathbf{j})=-g_{ud}(\mathbf{j},\mathbf{i})`.

:math:`g_{uu}` describes the pairing between a pair of electrons with spin-up:

.. math::
    g_{uu}(\mathbf{i},\mathbf{j}) = f_{uu}({\bf r}_i,{\bf r}_j) \left| \uparrow  \uparrow \right\rangle

where the function :math:`f_{uu}` describes the spatial dependence on the coordinates :math:`{\bf r}_i,{\bf r}_j` for :math:`i,j\le N_u`. The spin part :math:`\left| \uparrow  \uparrow \right\rangle` describes a system with unit total spin  and spin projection along the z-axis,  and will be indicated by :math:`\left| 1, +1 \right\rangle`.

Similarly, :math:`g_{dd}` describes the pairing between pairs of electrons with spin-down for :math:`i,j> N_u`:

.. math::

    g_{dd}(\mathbf{i},\mathbf{j}) = f_{dd}({\bf r}_i,{\bf r}_j) \left| \downarrow  \downarrow \right\rangle

with :math:`f_{dd}({\bf r}_j,{\bf r}_i) = - f_{dd}({\bf r}_i,{\bf r}_j) , and the spin part :math:`\left| \downarrow  \downarrow \right\rangle` describes a system with total unit spin  and negative spin projection along the z-axis, indicated with :math:`\left| 1, -1 \right\rangle`.

:math:`g_{ud}` describes the pairing between pairs of electrons with unlike spins.
Since two electrons with unlike spins can form a singlet
:math:`\left| 0,0 \right\rangle = { {\left| \uparrow  \downarrow \right\rangle - \left| \downarrow  \uparrow \right\rangle}\over \sqrt{2}}`
or a triplet :math:`\left| 1,0 \right\rangle = { {\left| \uparrow  \downarrow \right\rangle + \left| \downarrow  \uparrow \right\rangle}\over \sqrt{2}}`, in the general case
the pairing function :math:`g_{ud}` will be a linear combination of the the two components:

.. math::

    g_{ud}(\mathbf{i},\mathbf{j}) = f_{S}({\bf r}_i,{\bf r}_j) { {\left| \uparrow  \downarrow \right\rangle - \left| \downarrow  \uparrow \right\rangle}\over \sqrt{2}} + f_{T}({\bf r}_i,{\bf r}_j) { {\left| \uparrow  \downarrow \right\rangle + \left| \downarrow  \uparrow \right\rangle}\over \sqrt{2}}

where :math:`f_{S}({\bf r}_i,{\bf r}_j) = f_{S}({\bf r}_j,{\bf r}_i)` describes the spatial dependence of the singlet part of :math:`g_{ud}`, and :math:`f_{T}({\bf r}_i,{\bf r}_j) = -f_{T}({\bf r}_j,{\bf r}_i)` describes the spatial dependence of the triplet part.

Indeed, the generic pairing function :math:`g(\mathbf{i},\mathbf{j})` is the sum of  all the four components mentioned above, namely :

.. math::

    \begin{split}
    g\left( \mathbf{i},\mathbf{j} \right)
    &= f_{S}({\bf r}_i,{\bf r}_j) \left| 0,0 \right\rangle + f_{T}({\bf r}_i,{\bf r}_j) \left| 1,0 \right\rangle \\
    &+ f_{uu}({\bf r}_i,{\bf r}_j) \left| 1,+1 \right\rangle + f_{dd}({\bf r}_i,{\bf r}_j) \left| 1,-1 \right\rangle \,.
    \end{split}

The pairing functions :math:`f_{S}`, :math:`f_{T}`, :math:`f_{uu}`, and :math:`f_{dd}` are expanded over atomic orbitals Say, for a generic pairing function :math:`f` we have

.. math::

    f\left( {{{\mathbf{r}}_i},{{\mathbf{r}}_j}} \right) = \sum\limits_{l,m,a,b} {{{A}_{\left\{ {a,l} \right\},\left\{ {b,m} \right\}}}{\psi _{a,l}}\left( {{{\mathbf{r}}_i}} \right){\psi _{b,m}}\left( {{{\mathbf{r}}_j}} \right)},

where :math:`{\psi_{a,l}}` and :math:`{\psi_{b,m}}` are primitive or contracted atomic orbitals, their indices :math:`l` and :math:`m` indicate different orbitals centered on atoms :math:`a` and :math:`b`, while :math:`i` and :math:`j` label the electron coordinates.

Symmetries on the system, or properties of the underlying pairing function :math:`f` imply constraints on the coefficients. For instance, the coefficients of :math:`f_{S}` are such that :math:`A_{\left\{ {a,l} \right\},\left\{ {b,m} \right\}} = A_{\left\{ {b,m} \right\},\left\{ {a,l} \right\}}` because :math:`f_{S}({\bf r}_i,{\bf r}_j) = f_{S}({\bf r}_j,{\bf r}_i)`, whereas :math:`A_{\left\{ {a,l} \right\},\left\{ {b,m} \right\}} = -A_{\left\{ {b,m} \right\},\left\{ {a,l} \right\}}` for :math:`f_{T}`, :math:`f_{uu}`, and :math:`f_{dd}`.

In this part, all the values :math:`{A}_{\left\{ {a,l} \right\}}` different from zero are listed. The first two numbers indicate the orbital indices of the determinant :math:`A` matrix, the third number is the value for :math:`{A}_{\left\{ {a,l} \right\}}`. The number of elements different from zero is indicated in the **Header**: see :math:`DET \neq 0`.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
JsAGPs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simplest choice of considering only the case of a pairing function :math:`g(\mathbf{i},\mathbf{j})` that is a spin singlet (namely, :math:`f_{uu}`, :math:`f_{dd}` :math:`f_{T}` are set to zero, yielding :math:`g(\mathbf{i},\mathbf{j})=f_{S}({\bf r}_i,{\bf r}_j) \left| 0,0 \right\rangle`) then we obtain the singlet Antisymmetrized Geminal Power.

In this case, the matrices :math:`G_{uu}` and :math:`G_{dd}` defined are both zero matrices of size :math:`N/2\times N/2`, and the matrix :math:`G_{ud}` has only the contribution coming from the singlet, that we dub :math:`G_S` with :math:`G_S^T=G_S`.

The antisymmetrization operator implies the computation of

.. math::

     {\text{Pf}\left({\begin{array}{c|c} 0 & G_{S} \\ \hline
                    -G_{S}^T & 0\end{array}}\right)}
     = (-1)^{N/2\times (N/2-1)\over 2} \det(G_S)

where the equality follows from a property of the Pfaffian.
The overall sign is arbitrary for a WF; thus the antisymmetrized product of singlet pairs (geminals) is indeed equivalent to the computation of the determinant of the matrix :math:`G_S`:

.. math::

    \Phi_\text{AGPs} = \det(G_S) \,.

This is called ``JsAGPs``. Indeed, we consider only the singlet part of the paring function:

.. math::

    g_{ud}(\mathbf{i},\mathbf{j}) \equiv g_{s}(\mathbf{i},\mathbf{j}) = f_{S}({\bf r}_i,{\bf r}_j) { {\left| \uparrow  \downarrow \right\rangle - \left| \downarrow  \uparrow \right\rangle}\over \sqrt{2}},

where

.. math::

        f_{S}\left( {{{\mathbf{r}}_i},{{\mathbf{r}}_j}} \right) = \sum\limits_{l,m,a,b} {{{A}_{\left\{ {a,l} \right\},\left\{ {b,m} \right\}}}{\psi _{a,l}}\left( {{{\mathbf{r}}_i}} \right){\psi _{b,m}}\left( {{{\mathbf{r}}_j}} \right)}.

We have assumed the :math:`A` matrix is symmetric for the AGPs WF, only :math:`A_{ij}` for :math:`i \le j` are provided in this section::

   #	      Nonzero values of detmat
   	      1	      5	     9.421753101774391E-002
	      1	      6	     9.421753101774391E-002
	      1	      7	     9.421753101774391E-002

.. math::

    A =
    \begin{pmatrix}
    A_{11}          & A_{12} & \dots  & A_{1n} \\
                    & A_{22} & \dots  & A_{2n} \\
                    &        & \ddots & \vdots \\
                    &        &        & A_{nn}
    \end{pmatrix}

Please note that all :math:`A` parameters that are not explicitly declared in these lines are set to zero and are never optimized.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
JAGPu
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It should be noticed that it is not necessary that the matrix :math:`G_{ud}` is symmetric to reduce the Pfaffian to a single determinant evaluation. As long as the matrices :math:`G_{uu}` and :math:`G_{dd}` are zero, the Pfaffian is indeed equivalent to :math:`\det(G_{ud})` and describes an antisymmetric WF. However, if :math:`G_{ud}` is not symmetric the function

.. math::

    \Phi_\text{AGP} = \det(G_{ud})

is not an eigenstate of the spin. In other terms, there is a spin contamination, similarly to the case of unrestricted HF calculations. This is called ``JAGPu``. Indeed, we consider both singlet and triplet functions

.. math::

    g_{ud}(\mathbf{i},\mathbf{j}) = f_{S}({\bf r}_i,{\bf r}_j) { {\left| \uparrow  \downarrow \right\rangle - \left| \downarrow  \uparrow \right\rangle}\over \sqrt{2}} + f_{T}({\bf r}_i,{\bf r}_j) { {\left| \uparrow  \downarrow \right\rangle + \left| \downarrow  \uparrow \right\rangle}\over \sqrt{2}},

where

.. math::

        f_{X=S,T}\left( {{{\mathbf{r}}_i},{{\mathbf{r}}_j}} \right) = \sum\limits_{l,m,a,b} {{{A}_{\left\{ {a,l} \right\},\left\{ {b,m} \right\}}}{\psi _{a,l}}\left( {{{\mathbf{r}}_i}} \right){\psi _{b,m}}\left( {{{\mathbf{r}}_j}} \right)},

As written above, :math:`A` matrix is symmetric and skew-symmetric for the singlet part (:math:`f_{S}`) and the triplet parts (:math:`f_{T}`) respectively. So :math:`A_{ij}` for :math:`i \le j` and :math:`A_{ij}` for :math:`i > j` are the coefficients of the singlet and the triplet parts, respectively (i.e., the element :math:`i = j` of the skew-symmetric matrix should be zero)::

   #	      Nonzero values of detmat
            1       1       8.321544938822982E-001 <- singlet
            ....
            1       5       9.421753101774391E-002 <- singlet
            1       6       9.421753101774391E-002
            1       7       9.421753101774391E-002
            2       1       3.485892384239842E-003 <- triplet
            2       2       3.589529849283749E-001 <- singlet
            2       3       2.489548797987997E-002 <- singlet
            3       1       1.112333456889842E-003 <- triplet
            3       2       2.585777744345490E-001 <- triplet
            3       3       3.936485649473937E-002 <- singlet

.. math::

    A_S =
    \begin{pmatrix}
    A_{11}          & A_{12} & \dots  & A_{1n} \\
                    & A_{22} & \dots  & A_{2n} \\
                    &        & \ddots & \vdots \\
                    &        &        & A_{nn}
    \end{pmatrix}

.. math::

    A_T=
    \begin{pmatrix}
    0      & -A_{21} & \dots  & -A_{n1} \\
    A_{21} & 0       & \dots  & -A_{n2} \\
    \vdots & \vdots  & \ddots & \vdots \\
    A_{n1} & A_{n2}  & \dots  & 0
    \end{pmatrix}

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
JAGP (JPf)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most general case is the Pfaffian ansatz, which is called ``JAGP`` or ``JPf`` in TurboRVB.
:math:`A_{ij}` is::

   # This is a C2-dimer case
   # The number of basis set for each carbon is 4

   #       Nonzero values of detmat
                                                  <- (1,1) is always zero (G is skew-symmetric).
           1           2 -2.917621798712210E-002  <- A_{up,up}, triplet
           1           3  8.326474500954891E-003  <- A_{up,up}, triplet
           1           4 -0.228326252284219       <- A_{up,up}, triplet
           1           5  0.470855192339553       <- A_{up,up}, triplet
           1           6 -3.285258904186700E-002  <- A_{up,up}, triplet
           1           7 -5.097409720647310E-003  <- A_{up,up}, triplet
           1           8  5.679495868355650E-002  <- A_{up,up}, triplet
           1           9  0.684164602152446       <- A_{up,dn}, singlet
           ....
           1          16 -0.104285811627841       <- A_{up,dn}, singlet
           2           3 -2.076224374212450E-002  <- A_{up,up}, triplet
           ...
           2           8 -4.145465677435435E-003  <- A_{up,up}, triplet
           2           9  3.735515724267560E-003  <- A_{up,dn}, triplet
           2          10  0.520587530210701       <- A_{up,dn}, singlet
           ...
           2          16  4.428757569068110E-003  <- A_{up,dn}, singlet
           ...
           9          10  4.813787735439980E-003  <- A_{dn,dn}, triplet
           ...
          15          16  9.827312227017149E-003  <- A_{dn,dn}, triplet

.. math::

    \begin{align*}
        A = \left[\begin{array}{c|c} A_{\text{up}-\text{up}} & A_{\text{up}-\text{dn}} \\ \hline
        A_{\text{dn}-\text{up}} & A_{\text{dn}-\text{dn}}\end{array}\right]
    \end{align*}

where

.. math::

    A_{\text{up}-\text{up}}=
    \begin{pmatrix}
    0      & A_{1,2}  & \dots  & A_{1,8} \\
    -A_{1,2} & 0       & \dots  & A_{2,8} \\
    \vdots & \vdots  & \ddots & \vdots \\
    -A_{1,8} & -A_{2,8}  & \dots  & 0
    \end{pmatrix}

.. math::

    A_{\text{dn}-\text{dn}}=
    \begin{pmatrix}
    0      & A_{9,10}  & \dots  & A_{9,16} \\
    -A_{9,10} & 0       & \dots  & A_{10,16} \\
    \vdots & \vdots  & \ddots & \vdots \\
    -A_{9,16} & -A_{10,16}  & \dots  & 0
    \end{pmatrix}

and for the :math:`A_{\text{up}-\text{dn}}` part,

.. math::

    A_S =
    \begin{pmatrix}
    A_{1,9}         & A_{1,10} & \dots  & A_{1,16} \\
                    & A_{2,10} & \dots  & A_{2,16} \\
                    &        & \ddots & \vdots \\
                    &        &        & A_{8,16}
    \end{pmatrix}

.. math::

    A_T=
    \begin{pmatrix}
    0       & -A_{2,9} & \dots  & -A_{8,9} \\
    A_{2,9} & 0        & \dots  & -A_{8,10} \\
    \vdots & \vdots    & \ddots & \vdots \\
    A_{8,9} & A_{8,10}  & \dots  & 0
    \end{pmatrix}

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Symmetries of the Determinant :math:`A` matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This part lists all the symmetries involving the elements of the Determinant :math:`A` matrix. Each line indicates a symmetry (see **DET_L (iessw)**) in the **Header**. The first number sets the number of elements involved in the symmetry. Then the indices (i-j) of the :math:`A` matrix elements involved in the symmetry are listed::

     	  [NUMBER of Symmetric elements]     [(I-J)	  (L-K)	 [..]	     (N-M)]

In the first line of the following example, there are four symmetric elements indicated (1-5) (1-6) (1-7) (1-8); in the second line, three symmetric elements (2-2) (3-3) (4-4) and in the last line one symmetric element (1-1). This means that :math:`A_{15} = A_{16} = A_{17} = A_{18}, A_{22} = A_{33} = A_{44}.` If an element is not symmetric to others, the syntax will be as in the last line::

   #   	     Grouped par. in the chosen ordered basis
   	     4	     1	     5	 	1	6	1
	     7	     1	     8
	     3	     2	     2		3	3	4
	     4
	     1	     1	     1

The symmetries allow to reduce the dimension of the space of parameters and then speed up the optimization. It is also possible to freeze a set of parameters. In this case, these parameters will not be optimized. This can be done using the negative value for the number of the symmetric elements. For example, the matrix elements :math:`A_{15} = A_{16} = A_{17} = A_{18}` can be kept constant during optimization::

    #	       Grouped par. in the chosen ordered basis
    	       -4      1   5   1   6   1   7   1   8

Keep in mind that the elements listed in this part are the ones that are effectively optimized, so you must list all the :math:`\lambda` s different from zero. Finally, note that at least one element must be kept frozen during the optimization. For convenience, this element is :math:`\lambda_{11}`::

     	       -1     	1      1


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Coefficient of the Jastrow different from zero
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The three/four-body Jastrow factor reads:

.. math::

    J_{3/4}\left( {{{\mathbf{r}}_1}{\sigma _1}, \ldots, {{\mathbf{r}}_N}{\sigma _N}} \right) =
    \sum_{i < j}
    \left(
    \sum_{a,l} \sum_{b,m}
    M_{a,l,b,m}^{{\sigma _i},{\sigma _j}}
    \chi _{a,l}( \mathbf{r}_i )
    \chi _{b,m}( \mathbf{r}_j )
    \right),

where the indices :math:`l` and :math:`m` again indicate different orbitals centered on
corresponding atoms :math:`a` and :math:`b`,
and :math:`\{ M_{a,l,b,m}^{\sigma _i,\sigma_j} \}` are variational parameters.
If the three/four-body jastrow parts are spin-independent, which depends on ``JAS_2BODY``, the matrix elements for :math:`\sigma _i \neq \sigma_j` are set 0.

Sometimes it is convenient to set to zero part of the coefficients of the four-body Jastrow factor, namely those corresponding to :math:`a \ne b`, as they increase the overall variational space significantly and make the optimization more challenging, without being much more effective in improving the variational WF.

As the one-body Jastrow factor :math:`J_1` is the sum of two parts, the homogeneous part (enforcing the electron-ion cusp condition) and inhomogeneous parts. The inhomogeneous part reads:

.. math::

    {J_1^{inh}}\left( {{{\mathbf{r}}_1}{\sigma _1}, \ldots, {{\mathbf{r}}_N}{\sigma _N}} \right) =  \sum_{i=1}^N \sum_{a=1}^{N_\text{at}} \left( {\sum\limits_{l} {M_{a,l}^{{\sigma _i}} \chi_{a,l}\left( {{{\mathbf{r}}_i}} \right)} } \right) ,

where :math:`{{{\mathbf{r}}_i}}:math:` are the electron positions, :math:`{{{\mathbf{R}}_a}}` are the atomic positions with corresponding atomic number :math:`Z_a`, :math:`l` runs over atomic orbitals :math:`\chi _{a,l}` ({\it e.g.}, GTO) centered on the atom :math:`a`, :math:`{{\sigma _i}}` represents the electron spin (:math:`\uparrow` or :math:`\downarrow`), :math:`\{ M_{a,l}^{\sigma _i} \}` are variational parameters. The matrix elements of the inhomogeneous are also written in this section.

Similar to the section (**Coefficient of the determinant different from zero**), in this section, the elements of the Jastrow :math:`M` matrix that are different from zero are listed. The number of lines is provided in the **Header**: see **JAS\_L** :math:`\neq 0`.

If the three/four body jastrows are spin-independent::

  #          Nonzero values of  jasmat
           1           1  9.892797458899720E-003    <- up-up (dn-dn)
           1           2 -1.895931999217210E-002
           ....
           1          45 -2.042711544366610E-004
           1          91 -6.994272320227231E-004    <- inhomogeneous onebody M.

If the three/four body jastrows are spin-dependent::

  #          Nonzero values of  jasmat
           1           1  1.861662710209880E-003    <- up-up (dn-dn)
           1           2 -1.200055730317670E-002
           ...
           1          45 -9.409004718340540E-004    <- up-up (dn-dn)
           1          92  6.998162994255180E-003    <- up-dn (dn-up)
           1          93 -1.494759322340230E-002
           ...
           1         136 -4.067067721217150E-004
           1         182  3.207558596898210E-003    <- inhomogeneous onebody M.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Symmetries of the Jastrow :math:`M` matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here the symmetries involving the elements of the Jastrow :math:`M` matrix are listed. The syntax is as explained for the determinant. The number of symmetries read is given by **UNCONSTRAINED IESFREE** in the **Header**.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Symmetries on the Z coefficients in the Determinant Basis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this section, the set of symmetries involving the atomic variational parameters of the basis set used for expanding the determinant are listed. The number of these symmetries is provided in the **Header**: see **DET SYMM**. In each line, the first number indicates the elements involved in the symmetry, and then the elements are listed. The numbering depends upon the ordering of the basis set functions (see **The Basis Set for the Determinant**)::

   	#	 Eq. par. in the atomic Det par. in the chosen basis
		     1	     1
		     1	     2
		     4	     3		4	5	6

As described before, the elements that must be optimized must be listed and it is possible to freeze some parameters. In the latter case, the syntax is the same as in the other sections (use '-' sign).


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Symmetries on the Z coefficients in the Jastrow Basis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The number of symmetries is provided in the **Header: 3 BODY ATOMIC PAR SYMM (Eq. 3 body atomic. par.)**. The syntax is as described before.


---------------------------------------------
Generate a Wavefunction file (makefort10.x)
---------------------------------------------

``makefort10.x`` is a tool for generating a template JAGP WF(fort.10) from makefort10.input. Here, we show an example of ``makefort10.input``

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
system section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. table:: Parameter List

    +-----------------+----------+---------+----------------------------------------------------------------+
    | Parameter Name  | Datatype | Default | Description                                                    |
    +=================+==========+=========+================================================================+
    | ``posunits``    | str      | NA      | units for atomic positions (bohr| angstrom | crystal)          |
    +-----------------+----------+---------+----------------------------------------------------------------+
    | ``natoms``      | int      | NA      | Total number of electrons in the system.                       |
    +-----------------+----------+---------+----------------------------------------------------------------+
    | ``ntyp``        | int      | NA      | Total number of element types in the system.                   |
    +-----------------+----------+---------+----------------------------------------------------------------+
    |``complexfort10``| NA       | NA      | it generates a complex fort.10 if it is .true.                 |
    +-----------------+----------+---------+----------------------------------------------------------------+
    | ``pbcfort10``   | NA       | NA      | it generates a fort.10 for PBC if it is .true.                 |
    +-----------------+----------+---------+----------------------------------------------------------------+
    | ``yes_pfaff``   | NA       | NA      | it generates pfaffian WF it it is .true.                       |
    +-----------------+----------+---------+----------------------------------------------------------------+
    | ``celldm(1-6)`` | NA       | NA      | they specify lattice vectors following Quantum Espresso's      |
    |                 |          |         | convention.                                                    |
    +-----------------+----------+---------+----------------------------------------------------------------+
    | ``yes_tilted``  | NA       | NA      | non-orthorombic cell if it is .true. # specify celldm(4-6).    |
    +-----------------+----------+---------+----------------------------------------------------------------+
    | ``nxyz(1-3)``   | NA       | NA      | repetition of the cell in the three direction. Use this option |
    |                 |          |         | for exploiting translational symmetries.                       |
    +-----------------+----------+---------+----------------------------------------------------------------+
    | ``phase(1-3)``  | NA       | NA      | phase factors for up electrons                                 |
    +-----------------+----------+---------+----------------------------------------------------------------+
    | ``phasedo(1-3)``| NA       | NA      | phase factors for down electrons                               |
    +-----------------+----------+---------+----------------------------------------------------------------+

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
electron section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. table:: Parameter List

   +-----------------+----------+---------+-----------------------------------------------------+
   | Parameter       | Datatype | Default | Description                                         |
   +=================+==========+=========+=====================================================+
   | filling         | NA       | NA      |                                                     |
   |                 |          |         | ``filling = diagonal``: Fill the initial detmat(:,:)|
   |                 |          |         | and jasmat(:,:) with 1.d0 on the diagonal.          |
   |                 |          |         | ``filling = random``: Fill with random numbers      |
   |                 |          |         | between (0,1.d0).                                   |
   |                 |          |         | ``filling = semidiagonal``: Fill off-diagonal       |
   |                 |          |         | elements with random numbers between (0,0.1d0) and  |
   |                 |          |         | diagonal elements with 1.d0.                        |
   +-----------------+----------+---------+-----------------------------------------------------+
   | orbtype         | NA       | NA      |                                                     |
   |                 |          |         | ``orbtype = normal``: Use normal orbitals.          |
   |                 |          |         | ``orbtype = mixed``: Use mixed orbitals.            |
   |                 |          |         | ``orbtype = tempered``: Use tempered orbitals.      |
   |                 |          |         | The same applies for Jastrow orbitals with          |
   |                 |          |         | ``jorbtype``.                                       |
   +-----------------+----------+---------+-----------------------------------------------------+
   | twobody         | NA       | NA      | Type of the Jastrow 2B/1B term used to satisfy the  |
   |                 |          |         | electron-electron, electron-ion cusp conditions.    |
   +-----------------+----------+---------+-----------------------------------------------------+
   | twobodypar      | NA       | NA      | Twobody parameter, :math:`p`                        |
   +-----------------+----------+---------+-----------------------------------------------------+
   | onebodypar      | NA       | NA      | Onebody parameter, :math:`b`                        |
   +-----------------+----------+---------+-----------------------------------------------------+
   | yes_crystal     | NA       | NA      | Use the crystal basis for the determinant part.     |
   |                 |          |         | The default is ``.true.`` for a PBC case.           |
   +-----------------+----------+---------+-----------------------------------------------------+
   | yes_crystalj    | NA       | NA      | Use the crystal basis for the jastrow part.         |
   |                 |          |         | The default is ``.true.`` for a PBC case.           |
   +-----------------+----------+---------+-----------------------------------------------------+
   | no_4body_jas    | NA       | NA      | Does not use the 4-body jastrow factors when it is  |
   |                 |          |         | true.                                               |
   +-----------------+----------+---------+-----------------------------------------------------+
   | neldiff         | NA       | NA      | The difference in the number of up and down         |
   |                 |          |         | electrons.                                          |
   +-----------------+----------+---------+-----------------------------------------------------+


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
symmetry section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. table:: Parameter List

   +-------------+----------+---------+--------------------------------------------------------------+
   | Parameter   | Datatype | Default | Description                                                  |
   +=============+==========+=========+==============================================================+
   | nosym       | NA       | NA      |                                                              |
   |             |          |         | If ``nosym = .true.``, do not use symmetry, but only identity|
   |             |          |         | and inversion symmetries. The default is ``.false.``.        |
   +-------------+----------+---------+--------------------------------------------------------------+
   | eqatoms     | NA       | NA      |                                                              |
   |             |          |         | If ``eqatoms = .true.``, set the same value for all the      |
   |             |          |         | exponents of the atomic basis if acting on the same type of  |
   |             |          |         | atom. If ``eqatoms = .false.``, exponents corresponding to   |
   |             |          |         | different atomic positions are equal only if related by      |
   |             |          |         | spatial symmetries.                                          |
   +-------------+----------+---------+--------------------------------------------------------------+
   | rot_det     | NA       | NA      |                                                              |
   |             |          |         | This flag is used to exclude rotation symmetries from the    |
   |             |          |         | lambda matrix of the determinant. If ``rot_det = .false.``,  |
   |             |          |         | makefort10.x uses only translations and inversion symmetry,  |
   |             |          |         | if present. Note that the rotation symmetries are still used |
   |             |          |         | to determine the relation between the parameters of the      |
   |             |          |         | orbitals, so the result is slightly different from           |
   |             |          |         | ``nosym = .true.`` that does not use rotation for everything.|
   +-------------+----------+---------+--------------------------------------------------------------+
   | symmagp     | NA       | NA      |                                                              |
   |             |          |         | If ``symmagp = .false.``, create a fort.10 file with a lambda|
   |             |          |         | matrix that is not symmetric (e.g., JAGPu ansatz).           |
   +-------------+----------+---------+--------------------------------------------------------------+


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ATOMIC_POSITIONS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The unit is specified with ``posunits`` in the ``&system`` section.

ATOMIC_POSITIONS::

    4.0  6.0  0.31842955585522  0.63686011171043  0.00000000000000
    4.0  6.0  0.68157044414478  0.36313988828957  0.00000000000000

    # Ion coordinates
    N1  Z1                x1     y1     z1
    N2  Z2                x2     y2     z2
        ..                ..     ..     ..
    Nn  Zn                xn     yn     zn

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Basis set
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Basis Sets used for expanding the determinant and jastrow are described. Each shell of the determinant is described by two lines. The first one contains the multiplicity, the number of variational parameters of the shell function and the code describing the function. The code numbers and the description of the corresponding shell are described in the file makefun.f90 of the source code. The multiplicity depends on the shell type: Shells S, P and D have the multiplicities of 1, 3 and 5 respectively. In the second line the index of the nucleus on which the shell is centered is first indicated. Then the parameter values are listed. Keep in mind that the number of parameters to be read is given in the first line.::

    ATOM_6
    &shells
    nshelldet=18
    nshelljas=10
    !ndet_hyb=0
    /
    1   1   16
    1   13.073594000000
    1   1   16
    1   6.541187000000
    1   1   16
    1   3.272791000000
    1   1   16
    1   1.637494000000
    1   1   16
    1   0.819297000000
    1   1   16
    1   0.409924000000
    1   1   16
    1   0.205100000000
    1   1   16
    1   0.127852000000
    1   1   16
    1   0.102619000000
    3   1   36
    1   7.480076000000
    3   1   36
    1   3.741035000000
    3   1   36
    1   1.871016000000
    3   1   36
    1   0.935757000000
    3   1   36
    1   0.468003000000
    3   1   36
    1   0.234064000000
    3   1   36
    1   0.149161000000
    3   1   36
    1   0.117063000000
    5   1   68
    1   0.561160000000
    #  Parameters atomic Jastrow wf
    1   1   16
    1   1.637494000000
    1   1   16
    1   0.846879000000
    1   1   16
    1   0.409924000000
    1   1   16
    1   0.269659000000
    1   1   16
    1   0.109576000000
    3   1   36
    1   1.871016000000
    3   1   36
    1   0.935757000000
    3   1   36
    1   0.468003000000
    3   1   36
    1   0.117063000000
    5   1   68
    1   2.013760000000

All primitive orbitals are written in the source file makefun.f90 (open boundary), makefun_pbc.f90 (pbc) and makefun_bump.f90 (finite range orbitals). TurboRVB also implements standard contracted orbitals written  as a linear combination of :math:`p` primitive orbitals. The definitions are easily found (and can be easily implemented) in the fortran file: ioptorbcontr.f90. In this case, the number corresponding to "Number of par." is equal to :math:`2p`. In the next line, one writes these extra coefficients, :math:`C_i, i = 1,...2p:` the coefficient :math:`C_{i+p}` acts on the orbital number defined by the contracted orbital written in "Shell code", with exponent :math:`Z_i = C_i` (we omit the normalization, each orbital is assumed to be normalized), for instance a :math:`2s` contracted orbital:

.. math::

   \phi(r) = 3.231 \exp(-2r^2) + 7.54 \exp(-r^2)

is written as::

   #	   Parameters atomic wf
   	   1	      4		300
	   1	2.0   1.0  3.231  7.54

``ndet_hyb`` is the number of hybrid orbitals::

    ATOM_6
    &shells
    nshelldet=18
    nshelljas=10
    ndet_hyb=4
    /
    1   1   16
    1   13.073594000000
    1   1   16
    1   6.541187000000
    1   1   16
    1   3.272791000000
    ....

..
    -----------------------------------------
    How to choose basis set?
    -----------------------------------------

    There are four cases:

    - Open with pseudo potentials
    - Open with all-electrons
    - PBC with pseudo potentials
    - PBC with all-electrons

    Open with pseudo potentials
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    One can use a provided basis set as it is. Note that we recommend use **uncontracted** primitive basis instead of contracted basis this is because contraction coefficients are usually not suitable for QMC calculations.

    - `Energy-consistent pseudo potential <http://burkatzki.com/pseudos/index.2.html>`_
    - `Pseudopotential Library <https://pseudopotentiallibrary.org>`_

    Open with all electrons
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    - `Basis set exchange <https://www.basissetexchange.org>`_

    One can use a provided basis set, but it is better to cut several largest exponents, typically larger than :math:`8 \times Z^2`, where :math:`Z` is the atomic number. The WF in the vicinity of nuclei is compensated by the one-body Jastrow part in TurboRVB. Indeed, the electron-nuclei cusp conditions are exactly fulfilled for any basis ({\it i.e.}, Gaussian orbital) even within the DFT framework in TurboRVB. This is achieved by an appropriate modification of the standard basis sets commonly used ({\it e.g.}, ccpVTZ) for WF based calculations: the new basis
    is obtained by multiplying each element of the original basis by a suitably chosen one-body Jastrow factor introducing the correct cusps, namely:

    .. math::

        \tilde \phi _j^b\left( {{\mathbf{r}} - {{\mathbf{R}}_b}} \right) = \phi _j^b\left( {{\mathbf{r}} - {{\mathbf{R}}_b}} \right){{\tilde J}_1}\left( {\mathbf{r}} \right),

    where :math:`{{\tilde J}_1}\left( {\mathbf{r}} \right)` is the homogeneous one-body Jastrow part.


    PBC with pseudo potentials
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    TurboRVB uses the CRYSTAL periodic basis for PBC calculations:

    .. math::

        \psi _{l,m,I}^{{\text{PBC}}}\left( {{\mathbf{r}};\zeta } \right) = \sum\limits_{{{\mathbf{T}}_s}} {\psi _{l,m,I}^{}\left( {{\mathbf{r}} + {{\mathbf{T}}_s};\zeta } \right){e^{-i{{\mathbf{k}}_s} \cdot {{\mathbf{T}}_s}}}}

    where :math:`{{{\mathbf{k}}_{{s}}}}` is a twist vector (:math:`{{\mathbf{k}}_{{s}}} = \left( {k_s^x,k_s^y,k_s^z} \right)`), and :math:`{{\mathbf{T}}_s}` represents an arbitrary simulation cell vector, and :math:`{\psi _{l,m,I}^{}}` is a non-periodic real atomic orbital such as Gaussian.

    Unfortunately, the provided basis set is redundant for a periodic case, so we recommend that one should cut several smaller exponents, typically, smaller than 0.10.


    PBC with all-electrons
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    The same for all-electron cases. Basis sets provided for open systems such as `Basis set exchange <https://www.basissetexchange.org>`_ are usually redundant for a periodic case, so we recommend that one should cut several smaller exponents, typically, smaller than 0.10. Or, one can use all-electron basis set optimized for periodic cases such as ones used in the `CRYSTAL DFT code <https://www.crystal.unito.it/basis-sets.php>`_.
