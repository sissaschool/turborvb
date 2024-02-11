  ! d orbitals cartesian !!!
  ! R(r)= exp(-alpha r^2)
  ! each gaussian term is normalized



  !        indorbp=indorb
  indparp=indpar+1
  dd1=dd(indparp)

  !        if(iflagnorm.gt.2) then
  ! overall normalization
  c=dd1**1.75d0*1.64592278064948967213d0
  !        c=4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do


  do i=indtmin,indtm
    distp(i,2)=rmu(1,i)**2
    distp(i,3)=rmu(2,i)**2
    distp(i,4)=rmu(3,i)**2
    ! lz=+/-2
    distp(i,5)=rmu(1,i)*rmu(2,i)*cost3d
    ! lz=+/-1
    distp(i,6)=rmu(2,i)*rmu(3,i)*cost3d
    ! lz=+/-1
    distp(i,7)=rmu(1,i)*rmu(3,i)*cost3d
  end do


  do ic=1,6
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do k=i0,indtm
      z(indorbp,k)=distp(k,1)*distp(k,1+ic)
    end do
    !           endif
  end do


  if(typec.ne.1) then

    dd1=dd(indparp)
    fun0=distp(0,1)
    fun=-2.d0*dd1*distp(0,1)
    fun2=fun*(1.d0-2.d0*dd1*r(0)**2)

    !              indorbp=indorb
    do ic=1,6
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)                     &
            *fun
        if(ic.le.3) then
          if(i.eq.ic) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                2.d0*rmu(i,0)*fun0
          end if
        elseif(ic.eq.4) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                rmu(2,0)*fun0*cost3d
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                rmu(1,0)*fun0*cost3d
          end if
        elseif(ic.eq.5) then
          if(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                rmu(3,0)*fun0*cost3d
          elseif(i.eq.3) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                rmu(2,0)*fun0*cost3d
          end if
        elseif(ic.eq.6) then
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
      if(ic.le.3) then
        z(indorbp,indt+4)=distp(0,1+ic)*(6.d0*fun+fun2)+2.d0*distp(0,1)
      else
        z(indorbp,indt+4)=distp(0,1+ic)*(6.d0*fun+fun2)
      end if
      !endif for iocc
      !                endif
      ! enddo fot ic
    end do


    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshell+6
  indorb=indorbp

