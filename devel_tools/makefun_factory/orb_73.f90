

  !        indorbp=indorb
  indparp=indpar+1

  dd1=dd(indparp)

  !        if(iflagnorm.gt.2) then
  ! overall normalization  obtained by Mathematica
  c=dd1**3.75d0*0.43985656185609913955d0
  !  C= dd1^15/4 /Sqrt[Integrate[x^14 Exp[-2 x^2],{x,0,Infinity}]]
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  do i=indtmin,indtm
    do k=1,6
      zv(k)=rmu(3,i)**k
      yv(k)=rmu(2,i)**k
      xv(k)=rmu(1,i)**k
    end do
    r2=xv(2)+yv(2)+zv(2)
    r4=r2*r2
    r6=r2*r4
    ! lz=0
    distp(i,2)=cost1i*(231.d0*zv(6)-315.d0*zv(4)*r2+105.d0*zv(2)*r4-5.d0*r6)

    cost=(33.d0*zv(5)-30.d0*zv(3)*r2+5.d0*zv(1)*r4)
    ! lz=+/-1
    distp(i,3)=cost2i*rmu(1,i)*cost
    ! lz=+/-1
    distp(i,4)=cost2i*rmu(2,i)*cost

    cost=33.d0*zv(4)-18.d0*zv(2)*r2+r4
    ! lz=+/-2
    distp(i,5)=cost3i*(xv(2)-yv(2))*cost
    ! lz=+/-2
    distp(i,6)=2.d0*cost3i*xv(1)*yv(1)*cost



    cost=11.d0*zv(3)-3.d0*zv(1)*r2
    ! lz=+/-3
    distp(i,7)=cost4i*(xv(3)-3.d0*xv(1)*yv(2))*cost
    ! lz=+/-3
    distp(i,8)=-cost4i*(yv(3)-3.d0*yv(1)*xv(2))*cost


    cost=11.d0*zv(2)-r2
    ! lz=+/-4
    distp(i,9)=cost5i*(xv(4)-6.d0*xv(2)*yv(2)+yv(4))*cost
    ! lz=+/-4
    distp(i,10)=cost5i*4.d0*(xv(3)*yv(1)-yv(3)*xv(1))*cost


    ! lz=+/-5
    distp(i,11)=cost6i*(xv(5)-10.d0*xv(3)*yv(2)+5.d0*xv(1)*yv(4))*zv(1)
    ! lz=+/-5
    distp(i,12)=-cost6i*(-5.d0*xv(4)*yv(1)+10.d0*xv(2)*yv(3)-yv(5))*zv(1)

    ! lz=+/-6
    distp(i,13)=cost7i*(xv(6)-15.d0*xv(4)*yv(2)+15.d0*xv(2)*yv(4)-yv(6))
    ! lz=+/-6
    distp(i,14)=-cost7i*(-6.d0*xv(5)*yv(1)+20.d0*xv(3)*yv(3)-6.d0*yv(5)*xv(1))



  end do


  do ic=1,13
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
    do k=1,6
      zv(k)=rmu(3,0)**k
      yv(k)=rmu(2,0)**k
      xv(k)=rmu(1,0)**k
    end do

    r2=xv(2)+yv(2)+zv(2)
    r4=r2*r2
    r6=r2*r4

    !              indorbp=indorb
    do ic=1,13
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)*fun
      end do
      if(ic.eq.1) then
        !                       if(i.eq.1) then
        !   lz =0
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost1i*fun0*(-30.d0*xv(5)-60.d0*xv(3)*yv(2)-30.d0*xv(1)*yv(4)&
            +360.d0*xv(3)*zv(2)+360.d0*xv(1)*yv(2)*zv(2)-240.d0*xv(1)*zv(4))
        !                      elseif(i.eq.2) then
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost1i*fun0*(-30.d0*xv(4)*yv(1)-60.d0*xv(2)*yv(3)-30.d0*yv(5)&
            +360.d0*xv(2)*yv(1)*zv(2)+360.d0*yv(3)*zv(2)-240.d0*yv(1)*zv(4))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost1i*fun0*(180.d0*xv(4)*zv(1)+360.d0*xv(2)*yv(2)*zv(1)+180.d0*yv(4)*zv(1)&
            -480.d0*xv(2)*zv(3)-480.d0*yv(2)*zv(3)+96.d0*zv(5))
      elseif(ic.eq.2) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost2i*fun0*(25.d0*xv(4)*zv(1)+30.d0*xv(2)*yv(2)*zv(1)+5.d0*yv(4)*zv(1)&
            -60.d0*xv(2)*zv(3)-20.d0*yv(2)*zv(3)+8.d0*zv(5))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost2i*fun0*(20.d0*xv(3)*yv(1)*zv(1)+20.d0*xv(1)*yv(3)*zv(1)&
            -40.d0*xv(1)*yv(1)*zv(3))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost2i*fun0*(5.d0*xv(5)+10.d0*xv(3)*yv(2)+5.d0*yv(4)*xv(1)&
            -60.d0*xv(3)*zv(2)-60.d0*xv(1)*yv(2)*zv(2)+40.d0*xv(1)*zv(4))
      elseif(ic.eq.3) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            -cost2i*fun0*(-20.d0*xv(3)*yv(1)*zv(1)-20.d0*xv(1)*yv(3)*zv(1)&
            +40.d0*xv(1)*yv(1)*zv(3))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            -cost2i*fun0*(-5.d0*xv(4)*zv(1)-30.d0*xv(2)*yv(2)*zv(1)-25.d0*yv(4)*zv(1)&
            +20.d0*xv(2)*zv(3)+60.d0*yv(2)*zv(3)-8.d0*zv(5))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            -cost2i*fun0*(-5.d0*xv(4)*yv(1)-10.d0*xv(2)*yv(3)-5.d0*yv(5)&
            +60.d0*xv(2)*yv(1)*zv(2)+60.d0*yv(3)*zv(2)-40.d0*yv(1)*zv(4))
      elseif(ic.eq.4) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost3i*fun0*(6.d0*xv(5)+4.d0*xv(3)*yv(2)-2.d0*xv(1)*yv(4)&
            -64.d0*xv(3)*zv(2)+32.d0*xv(1)*zv(4))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost3i*fun0*(2.d0*xv(4)*yv(1)-4.d0*xv(2)*yv(3)-6.d0*yv(5)&
            +64.d0*yv(3)*zv(2)-32.d0*yv(1)*zv(4))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost3i*fun0*(-32.d0*xv(4)*zv(1)+32.d0*yv(4)*zv(1)+64.d0*xv(2)*zv(3)&
            -64.d0*yv(2)*zv(3))

      elseif(ic.eq.5) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            -cost3i*fun0*(-10.d0*xv(4)*yv(1)-12.d0*xv(2)*yv(3)-2.d0*yv(5)&
            +96.d0*xv(2)*yv(1)*zv(2)+32.d0*yv(3)*zv(2)-32.d0*yv(1)*zv(4))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            -cost3i*fun0*(-2.d0*xv(5)-12.d0*xv(3)*yv(2)-10.d0*xv(1)*yv(4)&
            +32.d0*xv(3)*zv(2)+96.d0*xv(1)*yv(2)*zv(2)-32.d0*xv(1)*zv(4))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            -cost3i*fun0*(64.d0*xv(3)*yv(1)*zv(1)+64.d0*xv(1)*yv(3)*zv(1)-128.d0*xv(1)*yv(1)*zv(3))
      elseif(ic.eq.6) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost4i*fun0*(-15.d0*xv(4)*zv(1)+18.d0*xv(2)*yv(2)*zv(1)+9.d0*yv(4)*zv(1)&
            +24.d0*xv(2)*zv(3)-24.d0*yv(2)*zv(3))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost4i*fun0*(12.d0*xv(3)*yv(1)*zv(1)+36.d0*xv(1)*yv(3)*zv(1)-48.d0*xv(1)*yv(1)*zv(3))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost4i*fun0*(-3.d0*xv(5)+6.d0*xv(3)*yv(2)+9.d0*xv(1)*yv(4)+24.d0*xv(3)*zv(2)&
            -72.d0*xv(1)*yv(2)*zv(2))
      elseif(ic.eq.7) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            -cost4i*fun0*(36.d0*xv(3)*yv(1)*zv(1)+12.d0*xv(1)*yv(3)*zv(1)-48.d0*xv(1)*yv(1)*zv(3))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            -cost4i*fun0*(9.d0*xv(4)*zv(1)+18.d0*xv(2)*yv(2)*zv(1)-15.d0*yv(4)*zv(1)&
            -24.d0*xv(2)*zv(3)+24.d0*yv(2)*zv(3))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            -cost4i*fun0*(9.d0*xv(4)*yv(1)+6.d0*xv(2)*yv(3)-3.d0*yv(5)&
            -72.d0*xv(2)*yv(1)*zv(2)+24.d0*yv(3)*zv(2))
      elseif(ic.eq.8) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost5i*fun0*(-6.d0*xv(5)+20.d0*xv(3)*yv(2)+10.d0*xv(1)*yv(4)&
            +40.d0*xv(3)*zv(2)-120.d0*xv(1)*yv(2)*zv(2))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost5i*fun0*(10.d0*xv(4)*yv(1)+20.d0*xv(2)*yv(3)-6.d0*yv(5)&
            -120.d0*xv(2)*yv(1)*zv(2)+40.d0*yv(3)*zv(2))

        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost5i*fun0*(20.d0*xv(4)*zv(1)-120.d0*xv(2)*yv(2)*zv(1)+20.d0*yv(4)*zv(1))
      elseif(ic.eq.9) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            -cost5i*fun0*(20.d0*xv(4)*yv(1)-4.d0*yv(5)-120.d0*xv(2)*yv(1)*zv(2)&
            +40.d0*yv(3)*zv(2))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            -cost5i*fun0*(4.d0*xv(5)-20.d0*xv(1)*yv(4)-40.d0*xv(3)*zv(2)&
            +120.d0*xv(1)*yv(2)*zv(2))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            -cost5i*fun0*(-80.d0*xv(3)*yv(1)*zv(1)+80.d0*xv(1)*yv(3)*zv(1))
      elseif(ic.eq.10) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost6i*fun0*(5.d0*xv(4)*zv(1)-30.d0*xv(2)*yv(2)*zv(1)+5.d0*yv(4)*zv(1))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost6i*fun0*(-20.d0*xv(3)*yv(1)*zv(1)+20.d0*xv(1)*yv(3)*zv(1))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost6i*fun0*(xv(5)-10.d0*xv(3)*yv(2)+5.d0*xv(1)*yv(4))
      elseif(ic.eq.11) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            -cost6i*fun0*(-20.d0*xv(3)*yv(1)*zv(1)+20.d0*xv(1)*yv(3)*zv(1))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            -cost6i*fun0*(-5.d0*xv(4)*zv(1)+30.d0*xv(2)*yv(2)*zv(1)-5.d0*yv(4)*zv(1))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            -cost6i*fun0*(-5.d0*xv(4)*yv(1)+10.d0*xv(2)*yv(3)-yv(5))
      elseif(ic.eq.12) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost7i*fun0*(6.d0*xv(5)-60.d0*xv(3)*yv(2)+30.d0*xv(1)*yv(4))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost7i*fun0*(-30.d0*xv(4)*yv(1)+60.d0*xv(2)*yv(3)-6.d0*yv(5))
      elseif(ic.eq.13) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            -cost7i*fun0*(-30.d0*xv(4)*yv(1)+60.d0*xv(2)*yv(3)-6.d0*yv(5))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            -cost7i*fun0*(-6.d0*xv(5)+60.d0*xv(3)*yv(2)-30.d0*xv(1)*yv(4))
      end if
      z(indorbp,indt+4)=distp(0,1+ic)*(14.d0*fun+fun2)
      !endif for iocc
      !                endif
    end do ! enddo fot ic
    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshell+13
  indorb=indorbp



  ! 2s gaussian for pseudo
