  !     4d without cusp and one parmater  derivative of 33


  dd1=dd(indpar+1)
  !         if(iflagnorm.gt.2) then
  !     c=                                                           &
      !                                                              &1.d0/(2.d0**3*3.d0)/dsqrt(56.d0*pi)*(2.d0*dd1)**(9.d0/2.d0)
      c=dd1**4.5d0*0.0710812062076410d0
  !         endif

  c0=-c
  c1=4.5d0*c/dd1


  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
  end do

  do i=indtmin,indtm
    distp(i,3)=distp(i,1)*(c0*r(i)**2+c1*r(i))
    ! lz=0
    distp(i,4)=(3.d0*rmu(3,i)**2-r(i)**2)*cost1d
    ! lz=+/
    distp(i,5)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d
    ! lz=+/-2
    distp(i,6)=rmu(1,i)*rmu(2,i)*cost3d
    ! lz=+/-1
    distp(i,7)=rmu(2,i)*rmu(3,i)*cost3d
    ! lz=+/-1
    distp(i,8)=rmu(1,i)*rmu(3,i)*cost3d
  end do

  !         indorbp=indorb

  do ic=1,5
    !            if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=distp(i,3+ic)*distp(i,3)
    end do
    !            endif
  end do

  if(typec.ne.1) then
    fun0=distp(0,3)
    fun=-dd1*distp(0,3)+distp(0,1)*(2.d0*c0*r(0)+c1)
    fun2=dd1**2*distp(0,3)+distp(0,1)*                               &
        (-2.d0*dd1*(2.d0*c0*r(0)+c1)+2.d0*c0)
    !              indorbp=indorb

    do ic=1,5
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0)                     &
            *fun/r(0)
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
      z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0)+fun2)
      !endif for iocc
      !                endif
      ! enddo fot ic
    end do
  end if
  !
  indpar=indpar+1
  indshell=indshell+5
  indorb=indorbp


  ! derivative of 36 with respect zeta
