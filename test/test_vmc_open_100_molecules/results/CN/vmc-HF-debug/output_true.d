 ----------------------------------------------------------------------------
 TurboRVB version 1.0.1  git rev. 6038e74f                                
 
   SISSA Quantum Monte Carlo Package
 
     Contacts: sorella@sissa.it (Prof. Sandro Sorella)
     Website: https://people.sissa.it/~sorella/TurboRVB_Manual/build/html
 
   When you publish a paper using TurboRVB, please cite the following paper.
 
     TurboRVB: a many-body toolkit for ab initio electronic simulations
     K. Nakano, C. Attaccalite, M. Barborini, L. Capriotti, M. Casula,
     E. Coccia, M. Dagrada, Y. Luo, G. Mazzola, A. Zen, and S. Sorella,
     J. Chem. Phys. 152, 204121 (2020), doi:10.1063/5.0005037
 
 ----------------------------------------------------------------------------
  Number of threads/mpi proc =           1
  Initial path : 
 /home/nkousuke/application/TurboRVB_kagayaki_test/test/test_vmc_open_100_molecu
 les/results/CN/vmc-HF-debug
  After  reading simulation  
  Warning: unrecognized disk_io, set as 'default' !
  Default value of pseudorandom = T
  After reading pseudo 
  Pseudopotential file name : pseudo.dat
  After reading vmc 
  Default value of parcut =   1.110223024625157E-014
  Default value of epsdgel =  1.000000000000000E-003
  Default value of nweight for forces nweight X           -1
  Default value for tolcg   =  1.000000000000000E-006
  Default value for minzj =  -1000.00000000000     
  Default value of minimum one-two body Jastrow =  5.000000000000000E-002
  After reading readio 
  Default value of  iread            0
  Default value of epsbas =  1.000000000000000E-008
  After reading parameters 
  Default value of add_pulay =           2
  Default value for ieser=           1
  Default value for isfix=           1
  Basis set cutoff chosen:  1.000000000000000E-008
  after parameters   
  Parameters: iesinv,iesm,iesd,iesfree,iessw,iesup,ieskin 
  Parameters before read            0           0           0           0
           0           0           0
  Reading the begin in fast. . . . 
  Reading zeta and rion  in fast . . . . 
  Reading ieskin in fast. . . . 
  Reading 2-body jastrow in fast . . . . 
  Reading det shells in fast. . . . 
  Reading jas shells in fast. . . . 
  Reading det occupation in fast. . . .          113
  Reading jas occupation in fast. . . .            0
  Reading det nnozero in fast. . . . 
  Reading det nnozero symmetries in fast....
  Reading jas nnozero in fast....
  Reading jas nnozero symmetries in fast....
  Reading Z-det  symmetries in fast....
  Reading Z-jas  symmetries in fast....
  Warning NO DMC cutoff on the local energy
  Default value of epscuttype=            0
  Read value of epscut=  0.000000000000000E+000
  Default value of tstep=   2.00000000000000     
  Default value for nintpsa           12
  Parameters after  read            1           0           0           0
           0           0           0           0
  Default value for nscra           20
  Default value for nbra           18
  kl read =          -7
  before dynamic 
  Warning srcomplex turned to false (real case) 
 
  tparf =  0.350000000000000     
 Default cutoff on the weight =  0.000000000000000E+000
 itestrfn =           2
  Scratch extension ./scratch
  iese_eff=           1
 #############################################
     EFFECTIVE CORE POTENTIAL CALCULATION     
 #############################################
  ************ nintpseudo read           12
 ------  parameters for pseudopotentials -------
 Max angular momentum pseudo            2
 # of quadrature points in the projector          12
 # of pseudo atoms           2
 # of gaussian for pseudo          14
 -----------------------------------------------
  Single mesh  =  0.000000000000000E+000
  scratch of determinant each           20
 sub_comm_diag uses           1 processors
  Number of corr functions written =           3
  iopt =           1
  initial iseed =           0   472827767
  initial random number =  3.715735906735063E-002
 
  START reading the wave function fort.10 
 
  Reading the begin . . . . 
  Number of different atomic species =           2
  Open Boundary Conditions 
  Reading ieskin . . . . 
 Number of moved ionic coordinates           2
 opposite phase = F
 same  phase = T
 yes_hermite before = F
  Reading det shells . . . . 
 Warn. Ions are not in ascending order in Det !
 Check shell det =          29
 USING           9 CONTRACTED DET ORBITALS
 USING           5 MOLECULAR  DET ORBITALS
  iesup dimension  atomic basis =          96        1816
  Reading jas shells . . . . 
  USING UNCONTRACTED JASTROW ORBITALS 
  Reading det occupation . . . . 
 Number of total det orbitals         113
 Number of occupied det orbitals           5
  Preliminary allocation   4.358400000000000E-005
  transpip =  1.033600000000000E-005
 Number of total det orbitals (root)         172
 Number of occupied det orbitals (root)         172
  Reading jas occupation . . . . 
 Number of total Jas orbitals           0
 Number of occupied Jas orbitals           0
  Before allocation standard 
  nion =            2
  nel =            9
  indt=          36
  npm=            6
  nwm=            1
  nws=            1
  nelorbj=            0
  iscrapip=            1
  iscraipsip=         1912
  nelorb=          172
  nelcol=          173
  nelup=            5
  nelorb_c=            5
  nelorbj_c=            0
  nelcol_c=            6
  nshell=           60
  nshellj=            0
  npsamax=            3
  nintpseudo=           12
  Memory required QMC =  1.286300000000000E-003 Gbyte
  Memory winv =  5.077440000000000E-004
  iscramax for master =        74415
  Reading det nnozero . . . . 
  Reading det nnozero symmetries ....
 Number of non zero geminal lambda for det            5
 Number of non fixed geminal lambda for det           5
 Number of non zero lambda (root)       15050
  Reading jas nnozero ....
 Location constant orbitals in Jastrow
  Reading jas nnozero symmetries ....           0
 Number of non zero geminal lambda for Jas           0
 Number of non fixed geminal lambda for Jas           0
  Number of accepted nnozeron Jas Sz            0
  Check repeated in the symmetry table Jastrow  
  Reading Z-AGP symmetries ....
  Touched det zeta par =          92
  Touched Jas zeta par =           0
  Time spent in  update_kgrid=  1.096725463867188E-005
  scale one body =  0.000000000000000E+000
 
  END reading the wave function fort.10 
 
  Warning TurboRVB needs   5.6000001E-08  Gygabyte RAM per processor 
  Default chosen yesfast            1
  Warning nmolfn, Speeding factor =             4   21.5000000000000     
 *********************************************
 *********************************************
 VARIATIONAL MONTE CARLO
 *********************************************
 *********************************************
  Passi qui I real-eff 
  Passi qui III eff-real 
 %%%%%%%%%%%%%%%%%% 2 BODY JASTROW %%%%%%%%%%%%%%%%
 Number 2 body Jastrow parameters            0
 2 body Jastrow not present
 Total independent parameters in Jastrow sector           0
 
  nwfix given =           1           2
  Number of parameters in SR =           2
  Hartree Atomic Units 
 ****************************************
 ****** INITIALIZATION ******************
 zeta, up and down            1   4.00000000000000                2           2
 zeta, up and down            2   5.00000000000000                3           2
  nmol_ip nelorbh_ip used =           4         172
  firstmol nmolfn after all            1           4
  Size arrays            1           1           1         172          45
        1032         860           1       74415          25        7052
          41       63468         369        7225           1          25
        1720          10           1           1         688         860
  Initialization OK 
           0 % progress starts!
              1 steps,  10 % done in       0.00 sec
              2 steps,  20 % done in       0.00 sec
              3 steps,  30 % done in       0.00 sec
              4 steps,  40 % done in       0.00 sec
              5 steps,  50 % done in       0.00 sec
              6 steps,  60 % done in       0.00 sec
              7 steps,  70 % done in       0.00 sec
              8 steps,  80 % done in       0.00 sec
              9 steps,  90 % done in       0.00 sec
             10 steps, 100 % done in       0.00 sec
  Before finalizeall 
  Length record unit 12  =           0
  All files written correctly ...
  #####   TurboRVB  profiling (sec.) #####   
  Time initialization =  6.959438323974609E-004
  Total time with no initialization =  3.944873809814453E-003
  Total time with no measures   9.412765502929688E-004
  Time measures =  3.537654876708984E-003
  Time main =  6.914138793945312E-006
  Tracing the qmc update  move 
  Time ratiovar=  7.328987121582031E-004
  Time uptabtot=  1.289844512939453E-004
  Tracing the main routines 
  Time uptabpip  in uptabtot=  0.000000000000000E+000
  Time upnewwf in ratiovar/uptabtot=  5.943775177001953E-004
  Time upscratch  =  3.440380096435547E-003
  accept. rate off diagonal moves =  0.586419753086420     
  Optimal nbra suggested =          68
  Average time for 1000 generations   0.394487380981445     
 Average inverse A  wf =  5.440961370598056E-003 +/-  5.165941976468545E-003
  Final tstep found    2.00000000000000     
