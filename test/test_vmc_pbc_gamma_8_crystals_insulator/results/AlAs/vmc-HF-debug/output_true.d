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
 stals_insulator/results/AlAs/vmc-HF-debug
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
  Reading det occupation in fast. . . .          104
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
  Default value for nscra           66
  Default value for nbra           64
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
  ************ nintpseudo read           18
 ------  parameters for pseudopotentials -------
 Max angular momentum pseudo            4
 # of quadrature points in the projector          18
 # of pseudo atoms           8
 # of gaussian for pseudo          76
 -----------------------------------------------
  Single mesh  =  0.000000000000000E+000
  scratch of determinant each           66
 sub_comm_diag uses           1 processors
  Number of corr functions written =           3
  iopt =           1
  initial iseed =           0   472827767
  initial random number =  3.715735906735063E-002
 
  START reading the wave function fort.10 
 
  Warning: using Crystal periodic basis set definition.  Complex wave function w
 ith phase is possible now! 
  Reading celldm . . . . 
  Reading the begin . . . . 
 Phase up spin  :     0.000000    0.000000    0.000000
 Phase down spin:     0.000000    0.000000    0.000000
  Gamma Point Caculation 
  Number of different atomic species =           2
  Periodic System 
  Rs =    2.07517942912600     
  Notice celldm(1) set to            1
  Unit of length a.u. 
  LMin =    10.6202608272025     
  Real Volume of the Cell =    1197.85858188610     
  Rescaled tstep    2.00000000000000     
  celscale before    10.6202608272025        10.6202608272025     
   10.6202608272025     
  Rescaling Distances . . . 
  Reading ieskin . . . . 
 Number of moved ionic coordinates          12
 opposite phase = F
 same  phase = T
 yes_hermite before = F
  Reading det shells . . . . 
 Warn. Ions are not in ascending order in Det !
 Check shell det =          33
 USING          32 CONTRACTED DET ORBITALS
 USING          16 MOLECULAR  DET ORBITALS
  iesup dimension  atomic basis =         256        9472
  Reading jas shells . . . . 
  USING UNCONTRACTED JASTROW ORBITALS 
  Reading det occupation . . . . 
 Number of total det orbitals         104
 Number of occupied det orbitals          16
  Preliminary allocation   2.273280000000000E-004
  transpip =  4.198400000000000E-005
 Number of total det orbitals (root)         288
 Number of occupied det orbitals (root)         288
  Reading jas occupation . . . . 
 Number of total Jas orbitals           0
 Number of occupied Jas orbitals           0
  Before allocation standard 
  nion =            8
  nel =           32
  indt=          54
  npm=            6
  nwm=            1
  nws=            1
  nelorbj=            0
  iscrapip=            1
  iscraipsip=         9550
  nelorb=          288
  nelcol=          288
  nelup=           16
  nelorb_c=           16
  nelorbj_c=            0
  nelcol_c=           16
  nshell=          136
  nshellj=            0
  npsamax=            3
  nintpseudo=           18
  Memory required QMC =  7.347424000000000E-003 Gbyte
  Memory winv =  4.349952000000001E-003
  iscramax for master =       236616
  Reading det nnozero . . . . 
  Reading det nnozero symmetries ....
 Number of non zero geminal lambda for det           16
 Number of non fixed geminal lambda for det          16
 Number of non zero lambda (root)       41616
  Reading jas nnozero ....
 Location constant orbitals in Jastrow
  Reading jas nnozero symmetries ....           0
 Number of non zero geminal lambda for Jas           0
 Number of non fixed geminal lambda for Jas           0
  Number of accepted nnozeron Jas Sz            0
  Check repeated in the symmetry table Jastrow  
  Reading Z-AGP symmetries ....
  Touched det zeta par =          72
  Touched Jas zeta par =           0
 Warning: updated kgrid considering Det/Jastrow:        1240           0
 Warning: lowest wf discarded=  8.018177204421034E-009
  Time spent in  update_kgrid=  2.007007598876953E-003
  Warning contribution perfect gas     pressure a.u. NOT included 
  0.000000000000000E+000
  scale one body =  0.000000000000000E+000
 
  END reading the wave function fort.10 
 
  Warning TurboRVB needs   5.6000001E-08  Gygabyte RAM per processor 
  Default chosen yesfast            1
  Warning nmolfn, Speeding factor =            16   9.00000000000000     
  Number of G vectors         2724
  Ewald Self Energy   -142.797020266388     
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
 zeta, up and down            3   3.00000000000000                1           2
 zeta, up and down            4   3.00000000000000                1           2
 zeta, up and down            5   5.00000000000000                3           2
 zeta, up and down            6   5.00000000000000                3           2
 zeta, up and down            7   5.00000000000000                3           2
 zeta, up and down            8   5.00000000000000                3           2
  nmol_ip nelorbh_ip used =          16         288
  firstmol nmolfn after all            1          16
  Size arrays            1           1           1         288         512
        4608        4608           1      236616         256       16992
          59      543744        1888       17281           1         256
        9216          33           1           1        4608        4608
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
             10 steps, 100 % done in       0.02 sec
  Before finalizeall 
  Length record unit 12  =           0
  All files written correctly ...
  #####   TurboRVB  profiling (sec.) #####   
  Time initialization =  0.107069015502930     
  Total time with no initialization =  0.142822027206421     
  Total time with no measures   2.067375183105469E-002
  Time measures =  0.228939533233643     
  Time main =  2.908706665039062E-005
  Tracing the qmc update  move 
  Time ratiovar=  1.897072792053223E-002
  Time uptabtot=  1.067399978637695E-003
  Tracing the main routines 
  Time uptabpip  in uptabtot=  0.000000000000000E+000
  Time upnewwf in ratiovar/uptabtot=  1.677203178405762E-002
  Time upscratch  =  0.228733062744141     
  accept. rate off diagonal moves =  0.652777777777778     
  Optimal nbra suggested =         709
  Average time for 1000 generations    14.2822027206421     
 Average inverse A  wf =  2.929832345862917E-004 +/-  1.731860824241673E-004
  Final tstep found    2.00000000000000     
