  ! d orbitals
  ! R(r)= c*exp(-z r^2)*(7/4/z-r^2)



  !        indorbp=indorb
  indparp=indpar+1
  dd1=dd(indparp)
  dd2=dsqrt(dd1)

  !        if(iflagnorm.gt.2) then
  ! overall normalization
  !        c=4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)*ratiocd
  c=dd1**1.75d0*ratiocd
  !        endif


  do k=indtmin,indtm
    cost=dd1*r(k)**2/(1.d0+dd2*r(k))
    distp(k,1)=c*dexp(-cost)
  end do




  do i=indtmin,indtm
    ! lz=0
    distp(i,2)=(3.d0*rmu(3,i)**2-r(i)**2)*cost1d
    ! lz=+/-2
    distp(i,3)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d
    ! lz=+/-2
    distp(i,4)=rmu(1,i)*rmu(2,i)*cost3d
    ! lz=+/-1
    distp(i,5)=rmu(2,i)*rmu(3,i)*cost3d
    ! lz=+/-1
    distp(i,6)=rmu(1,i)*rmu(3,i)*cost3d
  end do


  do ic=1,5
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do k=i0,indtm
      cost=(1.d0+0.5d0*dd2*r(k))/(1.d0+dd2*r(k))**2
      z(indorbp,k)=distp(k,1)*(7.d0/4.d0/dd1-r(k)**2*cost)*          &
          distp(k,1+ic)
    end do
    !           endif
  end do


  if(typec.ne.1) then
    rp1=dd1*r(0)**2
    rp2=dd2*r(0)
    rp3=(1.d0+rp2)**2

    cost=(1.d0+0.5d0*rp2)/rp3
    fun0=distp(0,1)*(7.d0/4.d0/dd1-r(0)**2*cost)

    fun=0.25d0*distp(0,1)*                                           &
        (-22.d0-49.d0*rp2-28.d0*rp1-rp1*rp2+2.d0*rp1**2)/rp3**2
    fun2=-0.25d0*distp(0,1)*                                         &
        (22.d0+54.d0*rp2-26.d0*rp1-158.d0*rp1*rp2-139.d0*rp1**2      &
        -42.d0*rp1**2*rp2+rp1**3+2.d0*rp1**3*rp2)/rp3**3

    !              indorbp=indorb
    do ic=1,5
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)                     &
            *fun
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
      z(indorbp,indt+4)=distp(0,1+ic)*(6.d0*fun+fun2)
      !endif for iocc
      !                endif
      ! enddo fot ic
    end do

    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshell+5
  indorb=indorbp

