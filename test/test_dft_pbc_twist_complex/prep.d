&simulation  
itestr4=-4  
iopt=1     
ksq=0.5d0
/

&pseudo
/

&vmc
/

&optimization
molopt=1  
/

&readio
/

&parameters
!yes_kpoints=.true.
yes_scemama=.false.
epsbas=1d-8
/

&molecul 
nx=16
/

&dft
maxit=10
epsdft=1d-4
nelocc=1
mixing=1.
optocc=1
epsshell=0.005
mixingder=0.25d0
typedft=1
/
2 
