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
    ax=.30
    ay=.30
    az=.30
    nx=20
    ny=20
    nz=23
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
    nelocc=7
/

2 2 2 2 2 2 2
