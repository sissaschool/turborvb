&simulation
    itestr4=-4
    iopt=1
    maxtime=3600
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
    ax=0.25
    ay=0.25
    az=0.25
    nx=50
    ny=30
    nz=30
/

&dft
    contracted_on=.false.
    maxit=50
    epsdft=1e-05
    mixing=0.5
    typedft=1
    optocc=0
    epsshell=0.0
    memlarge=.false.
    nelocc=1
/

2
