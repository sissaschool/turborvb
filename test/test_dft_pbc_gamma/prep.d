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
!epsbas=1.d-10
/

&molecul 
nx=32 
/

&dft
maxit=1000
memlarge=.false.
epsdft=1d-5
nelocc=1
mixing=0.2
optocc=1
epsshell=0.005
typedft=1
/
2
