  !     3d without cusp condition triple Z


  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd(indpar+3)
  dd3=dd(indpar+4)
  peff2=dd(indpar+5)

  !        if(iflagnorm.gt.2) then
  c=1/2.d0*dsqrt(5.d0/pi)                                            &
      /dsqrt(1/(2.d0*dd1)**7+2*peff/(dd1+dd2)**7                     &
      +peff**2/(2.d0*dd2)**7+2*peff2/(dd1+dd3)**7                    &
      +peff2**2/(2.d0*dd3)**7+2*peff*peff2/(dd2+dd3)**7)/dsqrt(720.d0)
  !        endif

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
    distp(k,2)=dexp(-dd2*r(k))
    distp(k,3)=dexp(-dd3*r(k))
  end do

  do i=indtmin,indtm
    distp(i,4)=c*(distp(i,1)+peff*distp(i,2)+peff2*distp(i,3))
    !lz=0
    distp(i,5)=(3.d0*rmu(3,i)**2-r(i)**2)*cost1d
    !lz=+/-2
    distp(i,6)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d
    ! lz=+/- 2
    distp(i,7)=rmu(1,i)*rmu(2,i)*cost3d
    ! lz=+/-1
    distp(i,8)=rmu(2,i)*rmu(3,i)*cost3d
    ! lz=+/-1
    distp(i,9)=rmu(1,i)*rmu(3,i)*cost3d
  end do

  !        indorbp=indorb

  do ic=1,5
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=distp(i,4+ic)*distp(i,4)
    end do
    !           endif
  end do

  if(typec.ne.1) then
    fun0=distp(0,4)
    fun=c*(-dd1*distp(0,1)-peff*dd2*distp(0,2)                       &
        -peff2*dd3*distp(0,3))
    fun2=c*(dd1**2*distp(0,1)+peff*dd2**2*distp(0,2)                 &
        +peff2*dd3**2*distp(0,3))

    !              indorbp=indorb

    do ic=1,5
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,4+ic)*rmu(i,0)*fun/r(0)
        if(ic.eq.1) then
          if(i.ne.3) then
            z(indorbp,indt+i)=z(indorbp,indt+i)-                     &
                2.d0*rmu(i,0)*fun0*cost1d
          else
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                4.d0*rmu(i,0)*fun0*cost1d
          end if
        elseif(ic.eq.2) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                2.d0*rmu(i,0)*fun0*cost2d
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)-                     &
                2.d0*rmu(i,0)*fun0*cost2d
          end if
        elseif(ic.eq.3) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                rmu(2,0)*fun0*cost3d
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                rmu(1,0)*fun0*cost3d
          end if
        elseif(ic.eq.4) then
          if(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                rmu(3,0)*fun0*cost3d
          elseif(i.eq.3) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                rmu(2,0)*fun0*cost3d
          end if
        elseif(ic.eq.5) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                rmu(3,0)*fun0*cost3d
          elseif(i.eq.3) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                rmu(1,0)*fun0*cost3d
            !endif for i
          end if
          !endif for ic
        end if
        !enddo for i
      end do
      z(indorbp,indt+4)=distp(0,4+ic)*(6.d0*fun/r(0)+fun2)
      !endif for iocc
      !                endif
      ! enddo fot ic
    end do


    !endif for indt
  end if

  indpar=indpar+5
  indshell=indshell+5
  indorb=indorbp

