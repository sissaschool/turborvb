  ! g single gaussian orbital
  ! derivative of 51
  ! R(r)= exp(-alpha r^2)
  ! normalized


  !        indorbp=indorb
  indparp=indpar+1

  dd1=dd(indparp)

  !        if(iflagnorm.gt.2) then
  ! overall normalization
  !        c=16.d0/dsqrt(105.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(11.d0/4.d0)
  c=dd1**2.75d0*1.11284691281640568826d0
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do


  do i=indtmin,indtm
    distp(i,2)=cost1g*(35.d0*rmu(3,i)**4                             &
        -30.d0*rmu(3,i)**2*r(i)**2+3.d0*r(i)**4)
    ! lz=0
    distp(i,3)=cost2g*rmu(1,i)*rmu(3,i)                              &
        *(7.d0*rmu(3,i)**2-3.d0*r(i)**2)
    ! lz=+/-1
    distp(i,4)=cost2g*rmu(2,i)*rmu(3,i)                              &
        *(7.d0*rmu(3,i)**2-3.d0*r(i)**2)
    ! lz=+/-1
    distp(i,5)=cost3g*(rmu(1,i)**2-rmu(2,i)**2)                      &
        *(7.d0*rmu(3,i)**2-r(i)**2)
    ! lz=+/-2
    distp(i,6)=cost3g*2.d0*rmu(1,i)*rmu(2,i)                         &
        *(7.d0*rmu(3,i)**2-r(i)**2)
    ! lz=+/-2
    distp(i,7)=cost4g*rmu(1,i)*rmu(3,i)                              &
        *(rmu(1,i)**2-3.0*rmu(2,i)**2)
    ! lz=+/-3
    distp(i,8)=cost4g*rmu(2,i)*rmu(3,i)                              &
        *(3.d0*rmu(1,i)**2-rmu(2,i)**2)
    ! lz=+/-3
    distp(i,9)=cost5g*(rmu(1,i)**4                                   &
        -6.d0*rmu(1,i)**2*rmu(2,i)**2+rmu(2,i)**4)
    ! lz=+/-4
    distp(i,10)=cost5g*4.d0*rmu(1,i)*rmu(2,i)                        &
        *(rmu(1,i)**2-rmu(2,i)**2)
    ! lz=+/-4
  end do


  do ic=1,9
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do k=i0,indtm
      z(indorbp,k)=distp(k,1)*(11.d0/4.d0/dd1-r(k)**2)*              &
          distp(k,1+ic)
    end do
    !           endif
  end do


  if(typec.ne.1) then

    dd1=dd(indparp)
    fun0=distp(0,1)*(11.d0/4.d0/dd1-r(0)**2)
    fun=distp(0,1)*(2.d0*dd1*r(0)**2-15.d0/2.d0)
    fun2=distp(0,1)*(-4.d0*dd1**2*r(0)**4                            &
        +21.d0*dd1*r(0)**2-15.d0/2.d0)


    !              indorbp=indorb
    do ic=1,9
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)*fun
        if(ic.eq.1) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost1g*fun0*(-60.d0*rmu(1,0)*rmu(3,0)**2+12.d0*rmu(1,0)*r(0)**2)
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost1g*fun0*(-60.d0*rmu(2,0)*rmu(3,0)**2+12.d0*rmu(2,0)*r(0)**2)
          else
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost1g*fun0*(80.d0*rmu(3,0)**3-48.d0*rmu(3,0)*r(0)**2)
          end if
        elseif(ic.eq.2) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost2g*fun0*(-9.d0*rmu(1,0)**2*rmu(3,0)-3.d0*rmu(2,0)**2*rmu(3,0)+4.d0*rmu(3,0)**3)
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost2g*fun0*(-6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0))
          else
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost2g*fun0*(-3.d0*rmu(1,0)*(rmu(1,0)**2+rmu(2,0)**2-4.d0*rmu(3,0)**2))
          end if
        elseif(ic.eq.3) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost2g*fun0*(-6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0))
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost2g*fun0*(-3.d0*rmu(1,0)**2*rmu(3,0)-9.d0*rmu(2,0)**2*rmu(3,0)+4.d0*rmu(3,0)**3)
          else
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost2g*fun0*(-3.d0*rmu(2,0)*(rmu(1,0)**2+rmu(2,0)**2-4.d0*rmu(3,0)**2))
          end if
        elseif(ic.eq.4) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost3g*fun0*(-4.d0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(3,0)**2))
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost3g*fun0*(4.d0*(rmu(2,0)**3-3.d0*rmu(2,0)*rmu(3,0)**2))
          else
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost3g*fun0*(12.d0*(rmu(1,0)**2-rmu(2,0)**2)*rmu(3,0))
          end if
        elseif(ic.eq.5) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost3g*fun0*(-2.d0*rmu(2,0)*(3.d0*rmu(1,0)**2+rmu(2,0)**2-6.d0*rmu(3,0)**2))
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost3g*fun0*(-2.d0*rmu(1,0)*(rmu(1,0)**2+3.d0*rmu(2,0)**2-6.d0*rmu(3,0)**2))
          else
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost3g*fun0*24.d0*rmu(1,0)*rmu(2,0)*rmu(3,0)
          end if
        elseif(ic.eq.6) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost4g*fun0*3.d0*(rmu(1,0)**2-rmu(2,0)**2)*rmu(3,0)
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                -cost4g*fun0*6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0)
          else
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost4g*fun0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(2,0)**2)
          end if
        elseif(ic.eq.7) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost4g*fun0*6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0)
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost4g*fun0*3.d0*(rmu(1,0)**2-rmu(2,0)**2)*rmu(3,0)
          else
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost4g*fun0*(3.d0*rmu(1,0)**2*rmu(2,0)-rmu(2,0)**3)
          end if
        elseif(ic.eq.8) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost5g*fun0*4.d0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(2,0)**2)
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost5g*fun0*4.d0*(rmu(2,0)**3-3.d0*rmu(1,0)**2*rmu(2,0))
          end if
        elseif(ic.eq.9) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost5g*fun0*4.d0*(3.d0*rmu(1,0)**2*rmu(2,0)-rmu(2,0)**3)
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)                      &
                +cost5g*fun0*4.d0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(2,0)**2)
          end if
        end if
        !enddo for i
      end do
      z(indorbp,indt+4)=distp(0,1+ic)*(10.d0*fun+fun2)
      !endif for iocc
      !                endif
      ! enddo fot ic
    end do


    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshell+9
  indorb=indorbp




