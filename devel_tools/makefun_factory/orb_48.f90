  ! f orbital
  ! 
  ! - angmom = 3
  ! - type = Gaussian
  ! - normalized = yes
  ! - angtype = spherical
  ! - npar = 1
  ! - multiplicity = 7
  !

  indparp=indpar+1
  dd1=dd(indparp)

  c=dd1**2.25d0*1.47215808929909374563d0

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  do i=indtmin,indtm
    distp(i,2)=cost1f*rmu(3,i)                                       &
        *(5.d0*rmu(3,i)**2-3.d0*r(i)**2)
    ! lz=0
    distp(i,3)=cost2f*rmu(1,i)                                       &
        *(5.d0*rmu(3,i)**2-r(i)**2)
    ! lz=+/-1
    distp(i,4)=cost2f*rmu(2,i)                                       &
        *(5.d0*rmu(3,i)**2-r(i)**2)
    ! lz=+/-1
    distp(i,5)=cost3f*rmu(3,i)                                       &
        *(rmu(1,i)**2-rmu(2,i)**2)
    ! lz=+/-2
    distp(i,6)=cost3f*2.d0*rmu(3,i)                                  &
        *rmu(1,i)*rmu(2,i)
    ! lz=+/-2
    distp(i,7)=cost4f*rmu(1,i)                                       &
        *(rmu(1,i)**2-3.d0*rmu(2,i)**2)
    ! lz=+/-3
    distp(i,8)=cost4f*rmu(2,i)                                       &
        *(3.d0*rmu(1,i)**2-rmu(2,i)**2)
    ! lz=+/-3
  end do

  do ic=1,7
    indorbp=indorb+ic
    do k=i0,indtm
      z(indorbp,k)=distp(k,1)*distp(k,1+ic)
    end do
  end do


  if(typec.ne.1) then

    dd1=dd(indparp)
    fun0=distp(0,1)
    fun=-2.d0*dd1*distp(0,1)
    fun2=fun*(1.d0-2.d0*dd1*r(0)**2)

    do ic=1,7
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)                     &
            *fun
      end do
      if(ic.eq.1) then
        z(indorbp,indt+1)=z(indorbp,indt+1)-                         &
            6.d0*cost1f*fun0*rmu(1,0)*rmu(3,0)
        z(indorbp,indt+2)=z(indorbp,indt+2)-                         &
            6.d0*cost1f*fun0*rmu(2,0)*rmu(3,0)
        z(indorbp,indt+3)=z(indorbp,indt+3)+                         &
            cost1f*fun0*(9.d0*rmu(3,0)**2-3.d0*r(0)**2)
      elseif(ic.eq.2) then
        z(indorbp,indt+1)=z(indorbp,indt+1)+                         &
            cost2f*fun0*(5.d0*rmu(3,0)**2-r(0)**2-2.d0*rmu(1,0)**2)
        z(indorbp,indt+2)=z(indorbp,indt+2)-                         &
            2.d0*cost2f*fun0*rmu(2,0)*rmu(1,0)
        z(indorbp,indt+3)=z(indorbp,indt+3)+                         &
            8.d0*cost2f*fun0*rmu(3,0)*rmu(1,0)
      elseif(ic.eq.3) then
        z(indorbp,indt+1)=z(indorbp,indt+1)-                         &
            2.d0*cost2f*fun0*rmu(1,0)*rmu(2,0)
        z(indorbp,indt+2)=z(indorbp,indt+2)+                         &
            cost2f*fun0*(5.d0*rmu(3,0)**2-r(0)**2-2.d0*rmu(2,0)**2)
        z(indorbp,indt+3)=z(indorbp,indt+3)+                         &
            8.d0*cost2f*fun0*rmu(3,0)*rmu(2,0)
      elseif(ic.eq.4) then
        z(indorbp,indt+1)=z(indorbp,indt+1)+                         &
            2.d0*cost3f*fun0*rmu(1,0)*rmu(3,0)
        z(indorbp,indt+2)=z(indorbp,indt+2)-                         &
            2.d0*cost3f*fun0*rmu(2,0)*rmu(3,0)
        z(indorbp,indt+3)=z(indorbp,indt+3)+                         &
            cost3f*fun0*(rmu(1,0)**2-rmu(2,0)**2)
      elseif(ic.eq.5) then
        z(indorbp,indt+1)=z(indorbp,indt+1)+                         &
            2.d0*cost3f*fun0*rmu(2,0)*rmu(3,0)
        z(indorbp,indt+2)=z(indorbp,indt+2)+                         &
            2.d0*cost3f*fun0*rmu(1,0)*rmu(3,0)
        z(indorbp,indt+3)=z(indorbp,indt+3)+                         &
            2.d0*cost3f*fun0*rmu(1,0)*rmu(2,0)
      elseif(ic.eq.6) then
        z(indorbp,indt+1)=z(indorbp,indt+1)+                         &
            3.d0*cost4f*fun0*(rmu(1,0)**2-rmu(2,0)**2)
        z(indorbp,indt+2)=z(indorbp,indt+2)-                         &
            6.d0*cost4f*fun0*rmu(1,0)*rmu(2,0)
      else
        z(indorbp,indt+1)=z(indorbp,indt+1)+                         &
            6.d0*cost4f*fun0*rmu(1,0)*rmu(2,0)
        z(indorbp,indt+2)=z(indorbp,indt+2)+                         &
            3.d0*cost4f*fun0*(rmu(1,0)**2-rmu(2,0)**2)
      end if
      z(indorbp,indt+4)=distp(0,1+ic)*(8.d0*fun+fun2)
    end do
  end if

  indpar=indpar+1
  indshell=indshell+7
  indorb=indorbp

