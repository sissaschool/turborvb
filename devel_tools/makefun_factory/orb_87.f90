  ! f orbitals
  ! R(r)= c*exp(-z r^2)*(9/4/z-r^2)



  !        indorbp=indorb
  indparp=indpar+1
  dd1=dd(indparp)
  dd2=dsqrt(dd1)

  !        if(iflagnorm.gt.2) then
  ! overall normalization
  !        c=8.d0/dsqrt(15.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(9.d0/4.d0)*ratiocf
  c=dd1**2.25d0*ratiocf
  !        endif

  do k=indtmin,indtm
    cost=dd1*r(k)**2/(1.d0+dd2*r(k))
    distp(k,1)=c*dexp(-cost)
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
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do k=i0,indtm
      cost=(1.d0+0.5d0*dd2*r(k))/(1.d0+dd2*r(k))**2
      z(indorbp,k)=distp(k,1)*(9.d0/4.d0/dd1-r(k)**2*cost)*          &
          distp(k,1+ic)
    end do
    !           endif
  end do


  if(typec.ne.1) then

    rp1=dd1*r(0)**2
    rp2=dd2*r(0)
    rp3=(1.d0+rp2)**2

    cost=(1.d0+0.5d0*rp2)/rp3
    fun0=distp(0,1)*(9.d0/4.d0/dd1-r(0)**2*cost)

    fun=0.25d0*distp(0,1)*                                           &
        (-26.d0-59.d0*rp2-36.d0*rp1-3.d0*rp1*rp2+2.d0*rp1**2)/rp3**2
    fun2=0.25d0*distp(0,1)*                                          &
        (-26.d0-66.d0*rp2+22.d0*rp1+178.d0*rp1*rp2+165.d0*rp1**2     &
        +54.d0*rp1**2*rp2+rp1**3-2.d0*rp1**3*rp2)/rp3**3


    !              indorbp=indorb
    do ic=1,7
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)                     &
            *fun
        if(ic.eq.1) then
          z(indorbp,indt+i)=z(indorbp,indt+i)-                       &
              6.d0*cost1f*fun0*rmu(i,0)*rmu(3,0)
          if(i.eq.3) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                cost1f*fun0*(15.d0*rmu(i,0)**2-3.d0*r(0)**2)
          end if
        elseif(ic.eq.2) then
          z(indorbp,indt+i)=z(indorbp,indt+i)-                       &
              2.d0*cost2f*fun0*rmu(i,0)*rmu(1,0)
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                cost2f*fun0*(5.d0*rmu(3,0)**2-r(0)**2)
          elseif(i.eq.3) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                10.d0*cost2f*fun0*rmu(i,0)*rmu(1,0)
          end if
        elseif(ic.eq.3) then
          z(indorbp,indt+i)=z(indorbp,indt+i)-                       &
              2.d0*cost2f*fun0*rmu(i,0)*rmu(2,0)
          if(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                cost2f*fun0*(5.d0*rmu(3,0)**2-r(0)**2)
          elseif(i.eq.3) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                10.d0*cost2f*fun0*rmu(i,0)*rmu(2,0)
          end if
        elseif(ic.eq.4) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                2.d0*cost3f*fun0*rmu(1,0)*rmu(3,0)
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)-                     &
                2.d0*cost3f*fun0*rmu(2,0)*rmu(3,0)
          else
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                cost3f*fun0*(rmu(1,0)**2-rmu(2,0)**2)
          end if
        elseif(ic.eq.5) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                2.d0*cost3f*fun0*rmu(2,0)*rmu(3,0)
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                2.d0*cost3f*fun0*rmu(1,0)*rmu(3,0)
          else
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                2.d0*cost3f*fun0*rmu(1,0)*rmu(2,0)
          end if
        elseif(ic.eq.6) then
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                3.d0*cost4f*fun0*(rmu(1,0)**2-rmu(2,0)**2)
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)-                     &
                6.d0*cost4f*fun0*rmu(1,0)*rmu(2,0)
          end if
        else
          if(i.eq.1) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                6.d0*cost4f*fun0*rmu(1,0)*rmu(2,0)
          elseif(i.eq.2) then
            z(indorbp,indt+i)=z(indorbp,indt+i)+                     &
                3.d0*cost4f*fun0*(rmu(1,0)**2-rmu(2,0)**2)
          end if
          !endif for ic
        end if
        !enddo for i
      end do
      z(indorbp,indt+4)=distp(0,1+ic)*(8.d0*fun+fun2)
      !endif for iocc
      !                endif
      ! enddo fot ic
    end do

    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshell+7
  indorb=indorbp

