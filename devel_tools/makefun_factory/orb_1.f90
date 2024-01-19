  ! s orbital
  ! 
  ! - angmom = 0
  ! - type = Slater
  ! - normalized = yes
  ! - angtype = spherical
  ! - npar = 1
  ! - multiplicity = 1
  !

  indshellp=indshell+1
  dd1=dd(indpar+1)
  c=dd1*dsqrt(dd1)*0.56418958354775628695d0

  indorbp=indorb+1
  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)
  end do

  if(typec.ne.1) then
    fun=-dd1*distp(0,1)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)/r(0)
    end do

    z(indorbp,indt+4)=(-2.d0*dd1/r(0)+dd1**2)*distp(0,1)
  end if

  indorb=indorbp
  indpar=indpar+1
  indshell=indshellp
