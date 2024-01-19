  !      3p without cusp condition
  !      r ( e^{-z2 r } + z1 e^{-z3 r } )



  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  dd3=dd(indpar+3)
  !        if(iflagnorm.gt.2) then
  c=1.d0/2.d0/dsqrt(240.d0*pi*(1.d0/(2.d0*dd1)**7                    &
      +2.d0*dd3/(dd1+dd2)**7+dd3**2/(2.d0*dd2)**7))
  !        endif
  !
  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=c*dexp(-dd2*r(k))
  end do
  !
  do i=indtmin,indtm
    distp(i,3)=r(i)*(distp(i,1)+dd3*distp(i,2))
  end do
  !
  !        indorbp=indorb
  !
  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,3)
    end do
    !           endif
  end do
  !
  !
  if(typec.ne.1) then
    fun0=distp(0,3)
    fun=(1.d0-dd1*r(0))*distp(0,1)                                   &
        +dd3*(1.d0-dd2*r(0))*distp(0,2)
    fun2=dd1*(dd1*r(0)-2.d0)*distp(0,1)                              &
        +dd3*dd2*(dd2*r(0)-2.d0)*distp(0,2)
    !
    !              indorbp=indorb
    !
    do ic=1,3
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun/r(0)
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*                                   &
          (4.d0*fun/r(0)+fun2)
      !
      !                endif
    end do
    !
    !
    !endif for indt
  end if
  !
  indpar=indpar+3
  indshell=indshell+3
  indorb=indorbp



  ! 4p single zeta
