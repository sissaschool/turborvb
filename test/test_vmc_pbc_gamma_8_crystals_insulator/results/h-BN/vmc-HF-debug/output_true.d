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
 /home/nkousuke/application/TurboRVB_kagayaki_test/test/test_vmc_pbc_gamma_8_cry
 stals_insulator/results/h-BN/vmc-HF-debug
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
  Reading det occupation in fast. . . .           52
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
  Default value for nscra           34
  Default value for nbra           32
  kl read =          -7
  Default values for kSq (precision Ewald)   1.000000000000000E-008
  Default value for kappar=   8.00000000000000     
  Default value of neighbors in Ewald            1
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
 # of pseudo atoms           4
 # of gaussian for pseudo          28
 -----------------------------------------------
  Single mesh  =  0.000000000000000E+000
  scratch of determinant each           34
 sub_comm_diag uses           1 processors
  Number of corr functions written =           3
  iopt =           1
  initial iseed =           0   472827767
  initial random number =  3.715735906735063E-002
 
  START reading the wave function fort.10 
 
  Warning: using Crystal periodic basis set definition.  Complex wave function w
 ith phase is possible now! 
  Reading celldm . . . . 
  Warning new mapping used case_map=8 !!! 
  Reading the begin . . . . 
 Phase up spin  :     0.000000    0.000000    0.000000
 Phase down spin:     0.000000    0.000000    0.000000
  Gamma Point Caculation 
  Number of different atomic species =           2
  Periodic System 
  Rs =    1.53855171419025     
  Notice celldm(1) set to            1
  Unit of length a.u. 
  LMin =    4.73185532183400     
  Real Volume of the Cell =    244.087133568664     
  Rescaled tstep    2.00000000000000     
  celscale before    4.73185532183400        4.73185532183414     
   12.5878436694250     
  Warning increased kappa Ewald by =    1.15470053837924     
  Rescaling Distances . . . 
  Reading ieskin . . . . 
 Number of moved ionic coordinates          12
 opposite phase = F
 same  phase = T
 yes_hermite before = F
  Reading det shells . . . . 
 Warn. Ions are not in ascending order in Det !
 Check shell det =          17
 USING          16 CONTRACTED DET ORBITALS
 USING           8 MOLECULAR  DET ORBITALS
  iesup dimension  atomic basis =         128        2496
  Reading jas shells . . . . 
  USING UNCONTRACTED JASTROW ORBITALS 
  Reading det occupation . . . . 
 Number of total det orbitals          52
 Number of occupied det orbitals           8
  Preliminary allocation   5.990400000000000E-005
  transpip =  1.203200000000000E-005
 Number of total det orbitals (root)         148
 Number of occupied det orbitals (root)         148
  Reading jas occupation . . . . 
 Number of total Jas orbitals           0
 Number of occupied Jas orbitals           0
  Before allocation standard 
  nion =            4
  nel =           16
  indt=          36
  npm=            6
  nwm=            1
  nws=            1
  nelorbj=            0
  iscrapip=            1
  iscraipsip=         2590
  nelorb=          148
  nelcol=          148
  nelup=            8
  nelorb_c=            8
  nelorbj_c=            0
  nelcol_c=            8
  nshell=           68
  nshellj=            0
  npsamax=            3
  nintpseudo=           12
  Memory required QMC =  1.688800000000000E-003 Gbyte
  Memory winv =  7.767040000000000E-004
  iscramax for master =        80412
  Reading det nnozero . . . . 
  Reading det nnozero symmetries ....
 Number of non zero geminal lambda for det            8
 Number of non fixed geminal lambda for det           8
 Number of non zero lambda (root)       11026
  Reading jas nnozero ....
 Location constant orbitals in Jastrow
  Reading jas nnozero symmetries ....           0
 Number of non zero geminal lambda for Jas           0
 Number of non fixed geminal lambda for Jas           0
  Number of accepted nnozeron Jas Sz            0
  Check repeated in the symmetry table Jastrow  
  Reading Z-AGP symmetries ....
  Touched det zeta par =          86
  Touched Jas zeta par =           0
 Warning: updated kgrid considering Det/Jastrow:        1460           0
 Warning: lowest wf discarded=  9.200227180221947E-009
  Time spent in  update_kgrid=  2.674102783203125E-003
  Warning contribution perfect gas     pressure a.u. NOT included 
  0.000000000000000E+000
  scale one body =  0.000000000000000E+000
 
  END reading the wave function fort.10 
 
  Warning TurboRVB needs   5.6000001E-08  Gygabyte RAM per processor 
  Default chosen yesfast            1
  Warning nmolfn, Speeding factor =             8   9.25000000000000     
  Number of G vectors         9699
  Ewald Self Energy   -185.038561363001     
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
 zeta, up and down            1   3.00000000000000                1           2
 zeta, up and down            2   3.00000000000000                1           2
 zeta, up and down            3   5.00000000000000                3           2
 zeta, up and down            4   5.00000000000000                3           2
  nmol_ip nelorbh_ip used =           8         148
  firstmol nmolfn after all            1           8
  Size arrays            1           1           1         148         128
        1184        1184           1       80412          64        6068
          41       97088         656        6217           1          64
        2368          17           1           1        1184        1184
  Initialization OK 
           0 % progress starts!
              1 steps,  10 % done in       0.01 sec
              2 steps,  20 % done in       0.01 sec
              3 steps,  30 % done in       0.01 sec
              4 steps,  40 % done in       0.01 sec
              5 steps,  50 % done in       0.01 sec
              6 steps,  60 % done in       0.01 sec
              7 steps,  70 % done in       0.01 sec
              8 steps,  80 % done in       0.01 sec
              9 steps,  90 % done in       0.01 sec
             10 steps, 100 % done in       0.01 sec
  Before finalizeall 
  Length record unit 12  =           0
  All files written correctly ...
  #####   TurboRVB  profiling (sec.) #####   
  Time initialization =  6.812095642089844E-003
  Total time with no initialization =  6.688404083251953E-002
  Total time with no measures   1.473331451416016E-002
  Time measures =  5.879092216491699E-002
  Time main =  8.106231689453125E-006
  Tracing the qmc update  move 
  Time ratiovar=  1.428675651550293E-002
  Time uptabtot=  2.129077911376953E-004
  Tracing the main routines 
  Time uptabpip  in uptabtot=  0.000000000000000E+000
  Time upnewwf in ratiovar/uptabtot=  1.341915130615234E-002
  Time upscratch  =  5.870556831359863E-002
  accept. rate off diagonal moves =  0.569444444444444     
  Optimal nbra suggested =         128
  Average time for 1000 generations    6.68840408325195     
 Average inverse A  wf =  4.521737614368151E-003 +/-  2.730000204632972E-003
  Final tstep found    2.00000000000000     
