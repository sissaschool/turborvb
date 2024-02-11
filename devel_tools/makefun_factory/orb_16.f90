  ! s orbital
  ! 
  ! - angmom = 0
  ! - type = Gaussian
  ! - normalized = yes
  ! - angtype = spherical
  ! - npar = 1
  ! - multiplicity = 1
  !
  ! = N * R
  !
  ! where N is the normalization constant
  ! N = (2*alpha/pi)**(3/4)
  !
  ! and R is the radial part
  ! R = exp(-alpha*r**2)
  !


  indshellp=indshell+1
  indorbp=indorb+1
  dd1=dd(indpar+1)

  if(dd1.ne.0.) then
    c=0.71270547035499016d0*dd1**0.75d0
  else
    c=1.d0
  end if

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)
  end do

  if(typec.ne.1) then
    ! the first derivative /r
    fun=-2.d0*dd1*distp(0,1)

    ! the second derivative
    fun2=fun*(1.d0-2.d0*dd1*r(0)*r(0))

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*fun+fun2

    if(typec.eq.2) then
      ! Backflow
      funb=(fun2-fun)/(r(0)*r(0))

      z(indorbp,indt+5)=funb*rmu(1,0)*rmu(1,0)+fun
      z(indorbp,indt+6)=funb*rmu(2,0)*rmu(2,0)+fun
      z(indorbp,indt+7)=funb*rmu(3,0)*rmu(3,0)+fun
      z(indorbp,indt+8)=funb*rmu(1,0)*rmu(2,0)
      z(indorbp,indt+9)=funb*rmu(1,0)*rmu(3,0)
      z(indorbp,indt+10)=funb*rmu(2,0)*rmu(3,0)

    end if
  end if

  indorb=indorbp
  indpar=indpar+1
  indshell=indshellp

