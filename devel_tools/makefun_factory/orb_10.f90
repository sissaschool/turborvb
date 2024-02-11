  ! s orbital
  ! 
  ! - angmom = 0
  ! - type = Slater
  ! - normalized = yes
  ! - angtype = spherical
  ! - npar = 1
  ! - multiplicity = 1
  !
  ! = N * R
  !
  ! 3s single zeta
  ! and R is the radial part
  ! R(r) = r**2*exp(-z1*r)
  !

  indshellp=indshell+1
  indorbp=indorb+1
  dd1=dd(indpar+1)
  c=dd1**3.5d0*0.11894160774351807429d0

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)*r(i)**2
  end do

  if(typec.ne.1) then
    fun=(2.d0-dd1*r(0))*distp(0,1)
    fun2=(2.d0-4*dd1*r(0)+(dd1*r(0))**2)*distp(0,1)
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp
  indpar=indpar+1
  indshell=indshellp
