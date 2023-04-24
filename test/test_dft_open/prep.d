&simulation   !
itestr4=-4   ! -4 standard,-8 use twice larger basis for the Hartree potential.
iopt=1       ! 1 initialize with no potential (no Hartree, xc, correlation)                   0 continue starting with wf read in fort.10_new and occupation                  read in occupationlevels.dat (both generated after iopt=1).                     NB
/
&pseudo
/
&vmc
/
&optimization
molopt=1   !  do not change (it works only with molecular orbitals)
/
&readio
writescratch=1 ! if writescratch=0 write all initial matrices on disc for faster                   continuation.
/

&parameters
yes_kpoints=.false.
/
&molecul        ! section molecular orbitals
ax=0.25
nx=16
nz=32
!nbufd=8192  !  input value for the buffer dimension,  Default=1000
/
&dft
contracted_on=.false. ! if .true. it acts in the contracted basis.                                      If .false (default) it acts in the uncontracted basis.
maxit=50       ! max number of iterations
epsdft=1d-5     ! Tollerance in convergence in total energy
!typeopt=4       !typeopt=0 use self consistency method with standard mixing,                     typeopt=2 Anderson mixing scheme looks much faster,                             typeopt=1 use steepest descent method with SR acceleration,
mixing=0.5d0   ! linear simple mixing in the density, choose a small value                       for convergence. If even this does not work, switch on                          optocc>0. You can also change iteration method with
typedft=1      ! typedft=0 Hartree, typedft=1 LDA (PZ 1981), typefit=2 LDA (OB                 1994), 1 and 2 differ in the interpolation of the correlation                   energy. typedft<0 (-1 and -2) corresponding fit by imposing
optocc=1        !  default optocc=-1 use standard occupation of levels.                            It works for closed shell (insulators or special cases).                        optocc>=0 ---> print levels at each iteration.
epsshell=1d-2    ! smearing of the Fermi distribution:                                       occupations(i) =2 (1 for LSDA)/( exp( (eig(i)-ef)/epsshell)+1)                  where ef is determined by sum(occupations(1:bands))=#electrons.                 Fo
memlarge=.false. ! memlarge=.true. optimize speed  with >> memory requirements.
!bands=13 !# el up+7 ! The number of lowest eigenvalues of Khon-Sham equations.                         Default bands=nelup (the number of spin up electrons)+7.
/
