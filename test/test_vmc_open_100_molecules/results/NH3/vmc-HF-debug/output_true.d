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
 les/results/NH3/vmc-HF-debug
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
  Reading det occupation in fast. . . .          148
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
  Default value for nscra           18
  Default value for nbra           16
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
 # of pseudo atoms           4
 # of gaussian for pseudo          21
 -----------------------------------------------
  Single mesh  =  0.000000000000000E+000
  scratch of determinant each           18
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
 Number of moved ionic coordinates          10
 opposite phase = F
 same  phase = T
 yes_hermite before = F
  Reading det shells . . . . 
 Warn. Ions are not in ascending order in Det !
 Check shell det =          45
 USING           9 CONTRACTED DET ORBITALS
 USING           4 MOLECULAR  DET ORBITALS
  iesup dimension  atomic basis =         123        1699
  Reading jas shells . . . . 
  USING UNCONTRACTED JASTROW ORBITALS 
  Reading det occupation . . . . 
 Number of total det orbitals         148
 Number of occupied det orbitals           4
  Preliminary allocation   4.077600000000000E-005
  transpip =  1.073200000000000E-005
 Number of total det orbitals (root)         197
 Number of occupied det orbitals (root)         197
  Reading jas occupation . . . . 
 Number of total Jas orbitals           0
 Number of occupied Jas orbitals           0
  Before allocation standard 
  nion =            4
  nel =            8
  indt=          36
  npm=            6
  nwm=            1
  nws=            1
  nelorbj=            0
  iscrapip=            1
  iscraipsip=         1780
  nelorb=          197
  nelcol=          197
  nelup=            4
  nelorb_c=            4
  nelorbj_c=            0
  nelcol_c=            4
  nshell=           81
  nshellj=            0
  npsamax=            3
  nintpseudo=           12
  Memory required QMC =  1.501332000000000E-003 Gbyte
  Memory winv =  5.169280000000000E-004
  iscramax for master =        97491
  Reading det nnozero . . . . 
  Reading det nnozero symmetries ....
 Number of non zero geminal lambda for det            4
 Number of non fixed geminal lambda for det           4
 Number of non zero lambda (root)       19503
  Reading jas nnozero ....
 Location constant orbitals in Jastrow
  Reading jas nnozero symmetries ....           0
 Number of non zero geminal lambda for Jas           0
 Number of non fixed geminal lambda for Jas           0
  Number of accepted nnozeron Jas Sz            0
  Check repeated in the symmetry table Jastrow  
  Reading Z-AGP symmetries ....
  Touched det zeta par =          77
  Touched Jas zeta par =           0
  Time spent in  update_kgrid=  1.215934753417969E-005
  scale one body =  0.000000000000000E+000
 
  END reading the wave function fort.10 
 
  Warning TurboRVB needs   5.6000001E-08  Gygabyte RAM per processor 
  Default chosen yesfast            1
  Warning nmolfn, Speeding factor =             4   24.6250000000000     
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
 zeta, up and down            1   5.00000000000000                4           1
 zeta, up and down            2   1.00000000000000                0           1
 zeta, up and down            3   1.00000000000000                0           1
 zeta, up and down            4   1.00000000000000                0           1
  nmol_ip nelorbh_ip used =           4         197
  firstmol nmolfn after all            1           4
  Size arrays            1           1           1         197          32
         788         788           1       97491          16        8077
          41       64616         328        8275           1          16
        1576           9           1           1         788         788
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
  Time initialization =  4.980564117431641E-004
  Total time with no initialization =  3.127813339233398E-003
  Total time with no measures   7.798671722412109E-004
  Time measures =  2.742767333984375E-003
  Time main =  5.006790161132812E-006
  Tracing the qmc update  move 
  Time ratiovar=  6.115436553955078E-004
  Time uptabtot=  8.487701416015625E-005
  Tracing the main routines 
  Time uptabpip  in uptabtot=  0.000000000000000E+000
  Time upnewwf in ratiovar/uptabtot=  5.342960357666016E-004
  Time upscratch  =  2.694606781005859E-003
  accept. rate off diagonal moves =  0.472222222222222     
  Optimal nbra suggested =          56
  Average time for 1000 generations   0.312781333923340     
 Average inverse A  wf =  1.089901032499795E-002 +/-  5.703799527966010E-003
  Final tstep found    2.00000000000000     
