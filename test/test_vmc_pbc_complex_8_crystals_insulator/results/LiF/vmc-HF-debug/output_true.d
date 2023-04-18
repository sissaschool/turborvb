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
 /home/nkousuke/application/TurboRVB_kagayaki_test/test/test_vmc_pbc_complex_8_c
 rystals_insulator/results/LiF/vmc-HF-debug
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
  Reading det occupation in fast. . . .          124
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
  Parameters after  read            2           0           0           0
           0           0           0           0
  Default value for nscra           66
  Default value for nbra           64
  kl read =          -7
  Default values for kSq (precision Ewald)   1.000000000000000E-008
  Default value for kappar=   8.00000000000000     
  Default value of neighbors in Ewald            1
  before dynamic 
 
  tparf =  0.350000000000000     
 Default cutoff on the weight =  0.000000000000000E+000
 itestrfn =           2
  Scratch extension ./scratch
  iese_eff=           2
 #############################################
     EFFECTIVE CORE POTENTIAL CALCULATION     
 #############################################
  ************ nintpseudo read           12
 ------  parameters for pseudopotentials -------
 Max angular momentum pseudo            2
 # of quadrature points in the projector          12
 # of pseudo atoms           8
 # of gaussian for pseudo          40
 -----------------------------------------------
  Single mesh  =  0.000000000000000E+000
  scratch of determinant each           66
 sub_comm_diag uses           1 processors
  Number of corr functions written =           4
  iopt =           1
  initial iseed =           0   472827767
  initial random number =  3.715735906735063E-002
 
  START reading the wave function fort.10 
 
  Warning: using Crystal periodic basis set definition.  Complex wave function w
 ith phase is possible now! 
  Reading celldm . . . . 
  Reading the begin . . . . 
  Warning complex symmetric geminal wave function  !!!
 Phase up spin  :     0.250000    0.250000    0.250000
 Phase down spin:     0.250000    0.250000    0.250000
  Complex phase calculation 
  Number of different atomic species =           2
  Periodic System 
  Rs =    1.48338404281600     
  Notice celldm(1) set to            1
  Unit of length a.u. 
  LMin =    7.59159676532219     
  Real Volume of the Cell =    437.521496908889     
  Rescaled tstep    2.00000000000000     
  celscale before    7.59159676532219        7.59159676532219     
   7.59159676532219     
  Rescaling Distances . . . 
  Reading ieskin . . . . 
 Number of moved ionic coordinates           0
  Warning algorithm with effective lambda for k-average 
 opposite phase = F
 same  phase = T
 yes_hermite before = F
  Reading det shells . . . . 
 Warn. Ions are not in ascending order in Det !
 Check shell det =          37
 USING          48 CONTRACTED DET ORBITALS
 USING          32 MOLECULAR  DET ORBITALS
  iesup dimension  atomic basis =         284       20764
  Reading jas shells . . . . 
  USING UNCONTRACTED JASTROW ORBITALS 
  Reading det occupation . . . . 
 Number of total det orbitals         124
 Number of occupied det orbitals          32
  Preliminary allocation   4.983360000000000E-004
  transpip =  8.760000000000000E-005
 Number of total det orbitals (root)         320
 Number of occupied det orbitals (root)         320
  Reading jas occupation . . . . 
 Number of total Jas orbitals           0
 Number of occupied Jas orbitals           0
  Before allocation standard 
  nion =            8
  nel =           32
  indt=          48
  npm=            6
  nwm=            1
  nws=            1
  nelorbj=            0
  iscrapip=            1
  iscraipsip=        20868
  nelorb=          320
  nelcol=          320
  nelup=           16
  nelorb_c=           32
  nelorbj_c=            0
  nelcol_c=           32
  nshell=          152
  nshellj=            0
  npsamax=            4
  nintpseudo=           12
  Memory required QMC =  7.751536000000000E-003 Gbyte
  Memory winv =  4.341760000000000E-003
  iscramax for master =       252712
  Reading det nnozero . . . . 
  Reading det nnozero symmetries ....
 Number of non zero geminal lambda for det           16
 Number of non fixed geminal lambda for det          16
  Passi qui from eff. to real I 
 Number of non zero lambda (root)       51360
  Reading jas nnozero ....
 Location constant orbitals in Jastrow
  Reading jas nnozero symmetries ....           0
 Number of non zero geminal lambda for Jas           0
 Number of non fixed geminal lambda for Jas           0
  Number of accepted nnozeron Jas Sz            0
  Check repeated in the symmetry table Jastrow  
  Reading Z-AGP symmetries ....
  Touched det zeta par =          96
  Touched Jas zeta par =           0
 Warning: updated kgrid considering Det/Jastrow:        2416           0
 Warning: lowest wf discarded=  3.993762081822332E-009
  Time spent in  update_kgrid=  3.153085708618164E-003
  Warning contribution perfect gas     pressure a.u. NOT included 
  0.000000000000000E+000
  scale one body =  0.000000000000000E+000
 
  END reading the wave function fort.10 
 
  Warning TurboRVB needs   8.0000000E-08  Gygabyte RAM per processor 
  Default chosen yesfast            1
  Warning nmolfn, Speeding factor =            32   5.00000000000000     
  Number of G vectors         2724
  Ewald Self Energy   -275.867093428320     
 *********************************************
 *********************************************
 VARIATIONAL MONTE CARLO
 *********************************************
 *********************************************
  Warning allowed translation via phase flux attaching
  Passi qui I real-eff 
  Passi qui III eff-real 
 %%%%%%%%%%%%%%%%%% 2 BODY JASTROW %%%%%%%%%%%%%%%%
 Number 2 body Jastrow parameters            0
 2 body Jastrow not present
 Total independent parameters in Jastrow sector           0
 
  nwfix given =           2           3
  Number of parameters in SR =           3
  Hartree Atomic Units 
 ****************************************
 ****** INITIALIZATION ******************
 zeta, up and down            1   1.00000000000000                0           1
 zeta, up and down            2   1.00000000000000                0           1
 zeta, up and down            3   1.00000000000000                0           1
 zeta, up and down            4   1.00000000000000                0           1
 zeta, up and down            5   7.00000000000000                4           3
 zeta, up and down            6   7.00000000000000                4           3
 zeta, up and down            7   7.00000000000000                4           3
 zeta, up and down            8   7.00000000000000                4           3
  Warning ip_reshuff set to 2 
  nmol_ip nelorbh_ip used =          16         640
  firstmol nmolfn after all            1          32
  Size arrays            1           1           1         640        2048
       20480       20480           1      252712         512       33920
          53     1085440        1696       34561           1         512
       20480          33           2           2       10240       10240
  Initialization OK 
           0 % progress starts!
              1 steps,  10 % done in       0.03 sec
              2 steps,  20 % done in       0.03 sec
              3 steps,  30 % done in       0.03 sec
              4 steps,  40 % done in       0.03 sec
              5 steps,  50 % done in       0.03 sec
              6 steps,  60 % done in       0.03 sec
              7 steps,  70 % done in       0.03 sec
              8 steps,  80 % done in       0.03 sec
              9 steps,  90 % done in       0.03 sec
             10 steps, 100 % done in       0.03 sec
  Before finalizeall 
  Length record unit 12  =           0
  All files written correctly ...
  #####   TurboRVB  profiling (sec.) #####   
  Time initialization =  2.637195587158203E-002
  Total time with no initialization =  0.286484956741333     
  Total time with no measures   5.034589767456055E-002
  Time measures =  0.262180089950562     
  Time main =  3.004074096679688E-005
  Tracing the qmc update  move 
  Time ratiovar=  4.675555229187012E-002
  Time uptabtot=  2.936124801635742E-003
  Tracing the main routines 
  Time uptabpip  in uptabtot=  0.000000000000000E+000
  Time upnewwf in ratiovar/uptabtot=  4.325962066650391E-002
  Time upscratch  =  0.261975526809692     
  accept. rate off diagonal moves =  0.510416666666667     
  Optimal nbra suggested =         333
  Average time for 1000 generations    28.6484956741333     
 Average inverse A  wf =  1.652884624120742E-003 +/-  1.123008742489681E-003
  Final tstep found    2.00000000000000     
