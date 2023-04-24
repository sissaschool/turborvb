&simulation
!itestr4=2   ! vmc value  ! default
ngen=10   ! number of measures of the local energy per walker
iseedr=235113183 ! initial random seed
nbra=16
kappar=8
!nscra=20
ksq=0.5d0
yes_sparse_choose=.false.
yes_sparse=.true.
/
&pseudo
/
&vmc
/
&readio
!iread=3  ! 3 for correlated sampling (lot of space disc required)
/
&parameters
ieskin=1
typedyncell=2 ! a,b,c  derivatives
yespress=.true.
epsbas=1d-6
yes_scemama=.false.
/
