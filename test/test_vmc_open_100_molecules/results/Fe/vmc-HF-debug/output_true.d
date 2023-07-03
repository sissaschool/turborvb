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
 les/results/Fe/vmc-HF-debug
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
  Reading det occupation in fast. . . .           73
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
  Default value for nintpsa           18
  Parameters after  read            1           0           0           0
           0           0           0           0
  Default value for nscra           34
  Default value for nbra           32
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
  ************ nintpseudo read           18
 ------  parameters for pseudopotentials -------
 Max angular momentum pseudo            3
 # of quadrature points in the projector          18
 # of pseudo atoms           1
 # of gaussian for pseudo           9
 -----------------------------------------------
  Single mesh  =  0.000000000000000E+000
  scratch of determinant each           34
 sub_comm_diag uses           1 processors
  Number of corr functions written =           3
  iopt =           1
  initial iseed =           0   472827767
  initial random number =  3.715735906735063E-002
 
  START reading the wave function fort.10 
 
  Reading the begin . . . . 
  Number of different atomic species =           1
  Open Boundary Conditions 
  Reading ieskin . . . . 
 Number of moved ionic coordinates           0
 opposite phase = F
 same  phase = T
 yes_hermite before = F
  Reading det shells . . . . 
 USING          21 CONTRACTED DET ORBITALS
 USING          10 MOLECULAR  DET ORBITALS
  iesup dimension  atomic basis =         280        2880
  Reading jas shells . . . . 
  USING UNCONTRACTED JASTROW ORBITALS 
  Reading det occupation . . . . 
 Number of total det orbitals          73
 Number of occupied det orbitals          10
  Preliminary allocation   6.912000000000000E-005
  transpip =  2.048000000000000E-005
 Number of total det orbitals (root)         130
 Number of occupied det orbitals (root)         130
  Reading jas occupation . . . . 
 Number of total Jas orbitals           0
 Number of occupied Jas orbitals           0
  Before allocation standard 
  nion =            1
  nel =           16
  indt=          54
  npm=            6
  nwm=            1
  nws=            1
  nelorbj=            0
  iscrapip=            1
  iscraipsip=         3057
  nelorb=          130
  nelcol=          134
  nelup=           10
  nelorb_c=           10
  nelorbj_c=            0
  nelcol_c=           14
  nshell=           40
  nshellj=            0
  npsamax=            3
  nintpseudo=           18
  Memory required QMC =  1.841832000000000E-003 Gbyte
  Memory winv =  9.817599999999999E-004
  iscramax for master =        75048
  Reading det nnozero . . . . 
  Reading det nnozero symmetries ....
 Number of non zero geminal lambda for det           10
 Number of non fixed geminal lambda for det          10
 Number of non zero lambda (root)        9035
  Reading jas nnozero ....
 Location constant orbitals in Jastrow
  Reading jas nnozero symmetries ....           0
 Number of non zero geminal lambda for Jas           0
 Number of non fixed geminal lambda for Jas           0
  Number of accepted nnozeron Jas Sz            0
  Check repeated in the symmetry table Jastrow  
  Reading Z-AGP symmetries ....
  Touched det zeta par =         166
  Touched Jas zeta par =           0
  Time spent in  update_kgrid=  9.059906005859375E-006
  scale one body =  0.000000000000000E+000
 
  END reading the wave function fort.10 
 
  Warning TurboRVB needs   5.6000001E-08  Gygabyte RAM per processor 
  Default chosen yesfast            1
  Warning nmolfn, Speeding factor =             6   10.8333333333333     
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
 zeta, up and down            1   16.0000000000000               10           6
  nmol_ip nelorbh_ip used =           6         130
  firstmol nmolfn after all            1           6
  Size arrays            1           1           1         130         160
        1820        1300           1       75048         100        7670
          59      122720         944        7801           1         100
        2600          17           1           1         780        1300
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
  Time initialization =  7.629394531250000E-004
  Total time with no initialization =  5.521059036254883E-003
  Total time with no measures   1.073122024536133E-003
  Time measures =  5.097866058349609E-003
  Time main =  1.192092895507812E-005
  Tracing the qmc update  move 
  Time ratiovar=  8.039474487304688E-004
  Time uptabtot=  1.468658447265625E-004
  Tracing the main routines 
  Time uptabpip  in uptabtot=  0.000000000000000E+000
  Time upnewwf in ratiovar/uptabtot=  5.648136138916016E-004
  Time upscratch  =  5.006074905395508E-003
  accept. rate off diagonal moves =  0.472222222222222     
  Optimal nbra suggested =         152
  Average time for 1000 generations   0.552105903625488     
 Average inverse A  wf =  8.114306714636027E-003 +/-  7.112302517587491E-003
  Final tstep found    2.00000000000000     
