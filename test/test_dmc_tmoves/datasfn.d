&simulation
itestr4=-3        ! LRDMC option 
ngen=10         ! Number of branching 
nscra=5         !Recomputing by scratch determinant each nscra accepted moves
nw=64             ! Number of walkers 
iseedr=536473883  ! Initial random number 
/
&pseudo 
nintpsa=6        ! 6 points integration pseudo
npsamax=1       ! Max number of atoms closer than 2 rc (pseudo core radius).
/
&dmclrdmc
tbra=0.01d0      ! DMC time between consecutive branching 
etry=-5.42d0    ! Guess of the energy
/        
&readio 
wherescratch='old'
/
&parameters
/
