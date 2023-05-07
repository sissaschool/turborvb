&simulation
itestr4=-6        ! LRDMC option 
ngen=10         ! Number of branching 
nscra=5         !Recomputing by scratch determinant each nscra accepted moves
nw=64             ! Number of walkers 
iseedr=536473883  ! Initial random number 
kappar=6.d0
ksq=0.5d0
!yesfast=0
max_target=0
/
&pseudo 
nintpsa=6        ! 6 points integration pseudo
npsamax=1       ! Max number of atoms closer than 2 rc (pseudo core radius).
/
&dmclrdmc
tbra=0.1d0      ! DMC time between consecutive branching 
etry=-5.42d0    ! Guess of the energy
plat(2)=1.2d0   ! eta in LRDMC eta = 1 + O(a^2) 
plat(1)=-0.16d0
alat=0.25d0     ! lattice space a 
parcutg=0       ! Standard LRDMC
alat2=2.64575131106459d0
iesrandoma=.false.
yes_fastbranch=.true.
/        
&readio 
wherescratch='old'
/
&parameters
epsbas=1d-7
ieser=1         ! Energy 
isfix=1         ! Variance 
/
