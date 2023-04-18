&simulation
itestr4=-6        ! LRDMC option 
ngen=10         ! Number of branching 
nscra=5         !Recomputing by scratch determinant each nscra accepted moves
nw=8             ! Number of walkers 
iseedr=536473883  ! Initial random number 
kappar=6.d0
ksq=0.5d0
membig=.false.
membigcpu=.false.
!yesfast=0
yes_sparse_choose=.false.
yes_sparse=.true.
/
&pseudo 
nintpsa=6        ! 6 points integration pseudo
!npsamax=1       ! Max number of atoms closer than 2 rc (pseudo core radius).
/
&dmclrdmc
tbra=0.01d0      ! DMC time between consecutive branching 
etry=-5.42d0    ! Guess of the energy
plat(1)=-0.16d0
plat(2)=1.2d0   ! eta in LRDMC eta = 1 + O(a^2) 
alat=0.25d0     ! lattice space a 
parcutg=0       ! Standard LRDMC
alat2=2.64575131106459d0
iesrandoma=.false.
yes_fastbranch=.false.
/        
&readio 
wherescratch='old'
/
&parameters
epsbas=1.d-6
ieser=1         ! Energy 
isfix=1         ! Variance 
yes_scemama=.false.
/
