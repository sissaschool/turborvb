  ! f single Slater orbital derivative of 70
  ! R(r)= (9.d0/2.0 1/dd1 - r) * exp(-alpha r)
  ! normalized


  !        indorbp=indorb
  indparp=indpar+1

  dd1=dd(indparp)

  !        if(iflagnorm.gt.2) then
  ! overall normalization
  ! l = 3
  ! \int d\omega Y*Y = 4 pi / (2 l + 1)
  ! \int dr r^{2 l + 2} Exp [- 2 dd1 r^2 ] =  7 * 5 * 3**2 / 2**2 / dd1**9
  !        c=1.d0/dsqrt(10.d0)*(2.d0/pi)**(1.d0/2.d0)*dd1**(9.d0/2.d0)/3.d0
  c=dd1**4.5d0*0.084104417400672d0
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
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
      z(indorbp,k)=distp(k,1)*distp(k,1+ic)*(9.d0/2.d0/dd1 - r(k))
    end do
    !           endif
  end do


  if(typec.ne.1) then

    dd1=dd(indparp)
    fun0=distp(0,1)*(9.d0/2.d0/dd1-r(0))
    fun=distp(0,1)*(dd1-11.d0/2.d0/r(0))
    fun2=dd1*distp(0,1)*(13.d0/2.d0-dd1*r(0))

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



  ! 3s -derivative of 34 with respect to dd1
