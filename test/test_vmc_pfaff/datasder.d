&simulation
ngen=10
nscra=32
nbra=32
kappar=7.d0
!ksq=0.5
iseedr=536473883
!yesfast=0
membig=.true.
membigcpu=.true.
!dielectric_ratio=4.d0
!dielectric_length=9.12d0
!case_diel=2
!neigh=27
yes_sparse=.true.
yes_sparse_choose=.false.
/
&pseudo
!nintpsa=32
!npsamax=3
/
&vmc
tstep=2.2
epscut=5d-5
/
&readio 
/
&parameters
!ieser=1
!isfix=1
ieskin=1
iesd=1
!iesfree=1
iesm=1
iesup=1
iessw=1
!yespress=.true.
epsder=1.d-5
epsbas=1.d-11
yes_scemama=.true.
!warp=.false.
!yes_scemama_open=.true.
/
