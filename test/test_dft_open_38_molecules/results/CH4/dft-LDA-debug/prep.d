&simulation
    itestr4=-4
    iopt=1
    maxtime=172800
/

&pseudo
    npsamax=4
/

&vmc
/

&optimization
    molopt=1
/

&readio
    writescratch=1
/

&parameters
    yes_kpoints=.false.
/

&kpoints
/

&molecul
    ax=.18
    ay=.18
    az=.18
    nx=38
    ny=38
    nz=38
/

&dft
    contracted_on=.true.
    maxit=50
    epsdft=0.001
    mixing=0.5
    typedft=1
    optocc=0
    epsshell=0.0
    memlarge=.false.
    nelocc=4
/

2 2 2 2
