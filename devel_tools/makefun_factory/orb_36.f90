  ! p orbital
  ! 
  ! - angmom = 1
  ! - type = Gaussian
  ! - normalized = yes
  ! - angtype = spherical
  ! - npar = 1
  ! - multiplicity = 3
  !

  dd1=dd(indpar+1)
  c=dd1**1.25d0*1.42541094070998d0

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  do ic=1,3
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,1)
    end do
  end do

  if(typec.ne.1) then
    fun0=distp(0,1)
    fun=-2.d0*dd1*distp(0,1)
    fun2=fun*(1.d0-2.d0*dd1*r(0)**2)

    do ic=1,3
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*fun
      end do
      z(indorbp,indt+ic)=z(indorbp,indt+ic)+fun0
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
    end do
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp
