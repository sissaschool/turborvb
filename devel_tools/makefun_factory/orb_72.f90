

  !        indorbp=indorb
  indparp=indpar+1

  dd1=dd(indparp)

  !        if(iflagnorm.gt.2) then
  ! overall normalization  obtained by Mathematica
  c=dd1**3.25d0*0.79296269381073167718d0
  !  C= dd1^13/4 /Sqrt[Integrate[x^12 Exp[-2 x^2],{x,0,Infinity}]]
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  do i=indtmin,indtm
    do k=1,5
      zv(k)=rmu(3,i)**k
      yv(k)=rmu(2,i)**k
      xv(k)=rmu(1,i)**k
    end do
    r2=xv(2)+yv(2)+zv(2)
    r4=r2*r2
    ! lz=0
    distp(i,2)=cost1h*(63.d0*zv(5)-70.d0*zv(3)*r2+15.d0*zv(1)*r4)

    cost=(21.d0*zv(4)-14.d0*zv(2)*r2+r4)
    ! lz=+/-1
    distp(i,3)=cost2h*rmu(1,i)*cost
    ! lz=+/-1
    distp(i,4)=cost2h*rmu(2,i)*cost

    cost=3.d0*zv(3)-zv(1)*r2
    ! lz=+/-2
    distp(i,5)=cost3h*(xv(2)-yv(2))*cost
    ! lz=+/-2
    distp(i,6)=2.d0*cost3h*xv(1)*yv(1)*cost

    cost=9.d0*zv(2)-r2
    ! lz=+/-3
    distp(i,7)=cost4h*(xv(3)-3.d0*xv(1)*yv(2))*cost
    ! lz=+/-3
    distp(i,8)=-cost4h*(yv(3)-3.d0*yv(1)*xv(2))*cost

    ! lz=+/-4
    distp(i,9)=cost5h*(xv(4)-6.d0*xv(2)*yv(2)+yv(4))*zv(1)
    ! lz=+/-4
    distp(i,10)=cost5h*4.d0*(xv(3)*yv(1)-yv(3)*xv(1))*zv(1)
    ! lz=+/-5
    distp(i,11)=cost6h*(xv(5)-10.d0*xv(3)*yv(2)+5.d0*xv(1)*yv(4))
    ! lz=+/-5
    distp(i,12)=-cost6h*(-5.d0*xv(4)*yv(1)+10.d0*xv(2)*yv(3)-yv(5))

  end do


  do ic=1,11
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
    do k=1,5
      zv(k)=rmu(3,0)**k
      yv(k)=rmu(2,0)**k
      xv(k)=rmu(1,0)**k
    end do

    r2=xv(2)+yv(2)+zv(2)
    r4=r2*r2


    !              indorbp=indorb
    do ic=1,11
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)*fun
      end do
      if(ic.eq.1) then
        !                       if(i.eq.1) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost1h*fun0*20.d0*xv(1)*zv(1)*(3.d0*xv(2)+3.d0*yv(2)-4.d0*zv(2))
        !                      elseif(i.eq.2) then
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost1h*fun0*20.d0*yv(1)*zv(1)*(3.d0*xv(2)+3.d0*yv(2)-4.d0*zv(2))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost1h*fun0*(175.d0*zv(4)-150.d0*zv(2)*r2+15.d0*r4)
      elseif(ic.eq.2) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost2h*fun0*(5.d0*xv(4)+6.d0*xv(2)*yv(2)+yv(4)-36.d0*xv(2)*zv(2)&
            -12.d0*yv(2)*zv(2)+8.d0*zv(4))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost2h*fun0*(4.d0*xv(3)*yv(1)+4.d0*xv(1)*yv(3)-24.d0*xv(1)*yv(1)*zv(2))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost2h*fun0*(-24.d0*xv(3)*zv(1)-24.d0*xv(1)*yv(2)*zv(1)+32.d0*zv(3)*xv(1))
      elseif(ic.eq.3) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            -cost2h*fun0*(-4.d0*xv(3)*yv(1)-4.d0*xv(1)*yv(3)+24.d0*xv(1)*yv(1)*zv(2))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost2h*fun0*(5.d0*yv(4)+6.d0*xv(2)*yv(2)+xv(4)-36.d0*yv(2)*zv(2)&
            -12.d0*xv(2)*zv(2)+8.d0*zv(4))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost2h*fun0*(-24.d0*yv(3)*zv(1)-24.d0*yv(1)*xv(2)*zv(1) &
            +32.d0*zv(3)*yv(1))
      elseif(ic.eq.4) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost3h*fun0*(-4.d0*xv(3)*zv(1)+4.d0*xv(1)*zv(3))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost3h*fun0*(4.d0*yv(3)*zv(1)-4.d0*yv(1)*zv(3))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost3h*fun0*(-xv(4)+yv(4)+6.d0*xv(2)*zv(2)-6.d0*yv(2)*zv(2))
      elseif(ic.eq.5) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            -cost3h*fun0*(6.d0*xv(2)*yv(1)*zv(1)+2.d0*yv(3)*zv(1)-4.d0*yv(1)*zv(3))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            -cost3h*fun0*(2.d0*xv(3)*zv(1)+6.d0*xv(1)*yv(2)*zv(1)-4.d0*xv(1)*zv(3))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            -cost3h*fun0*(2.d0*xv(3)*yv(1)+2.d0*xv(1)*yv(3)-12.d0*xv(1)*yv(1)*zv(2))
      elseif(ic.eq.6) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost4h*fun0*(-5.d0*xv(4)+6.d0*xv(2)*yv(2)+3.d0*yv(4)+24.d0*xv(2)*zv(2)-24.d0*yv(2)*zv(2))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost4h*fun0*(4.d0*xv(3)*yv(1)+12.d0*xv(1)*yv(3)-48.d0*xv(1)*yv(1)*zv(2))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost4h*fun0*(16.d0*xv(3)*zv(1)-48.d0*xv(1)*yv(2)*zv(1))
      elseif(ic.eq.7) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            -cost4h*fun0*(12.d0*xv(3)*yv(1)+4.d0*xv(1)*yv(3)-48.d0*xv(1)*yv(1)*zv(2))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            -cost4h*fun0*(3.d0*xv(4)+6.d0*xv(2)*yv(2)-5.d0*yv(4)-24.d0*xv(2)*zv(2)+24.d0*yv(2)*zv(2))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            -cost4h*fun0*(-48.d0*xv(2)*yv(1)*zv(1)+16.d0*yv(3)*zv(1))
      elseif(ic.eq.8) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost5h*fun0*(4.d0*xv(3)*zv(1)-12.d0*xv(1)*yv(2)*zv(1))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost5h*fun0*(-12.d0*xv(2)*yv(1)*zv(1)+4.d0*yv(3)*zv(1))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost5h*fun0*(xv(4)-6.d0*xv(2)*yv(2)+yv(4))
      elseif(ic.eq.9) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            -cost5h*fun0*(-12.d0*xv(2)*yv(1)*zv(1)+4.d0*yv(3)*zv(1))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            -cost5h*fun0*(-4.d0*xv(3)*zv(1)+12.d0*xv(1)*yv(2)*zv(1))
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            -cost5h*fun0*(-4.d0*xv(3)*yv(1)+4.d0*xv(1)*yv(3))
      elseif(ic.eq.10) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost6h*fun0*(5.d0*xv(4)-30.d0*xv(2)*yv(2)+5.d0*yv(4))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost6h*fun0*(-20.d0*xv(3)*yv(1)+20.d0*xv(1)*yv(3))
      elseif(ic.eq.11) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            -cost6h*fun0*(-20.d0*xv(3)*yv(1)+20.d0*xv(1)*yv(3))
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            -cost6h*fun0*(-5.d0*xv(4)+30.d0*xv(2)*yv(2)-5.d0*yv(4))
      end if
      z(indorbp,indt+4)=distp(0,1+ic)*(12.d0*fun+fun2)
      !endif for iocc
      !                endif
    end do ! enddo fot ic
    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshell+11
  indorb=indorbp



  ! 2s gaussian for pseudo

