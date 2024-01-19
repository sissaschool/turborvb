  !     2p single Z with no  cusp condition


  dd1=dd(indpar+1)
  !        if(iflagnorm.gt.2) then
  !        c=dsqrt((2.d0*dd1)**5/8.d0/pi)/2.d0
  c=dd1**2.5d0*0.5641895835477562d0
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
  end do

  !        indorbp=indorb
  !
  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,1)
    end do
    !           endif
  end do

  if(typec.ne.1) then
    fun0=distp(0,1)
    fun=-dd1*distp(0,1)
    fun2=dd1**2*distp(0,1)

    !              indorbp=indorb

    do ic=1,3
      !                if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun/r(0)
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun/r(0)+fun2)
      !                 endif
    end do
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp


  ! 2p double zeta
