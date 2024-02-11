!     d gaussian  r**(2*npower)*exp(-alpha*r**2)

npower=iopt-1200

!        indorbp=indorb

dd2=dd(indpar+1)
do k=indtmin,indtm
  distp(k,1)=r(k)**(2*npower)*dexp(-dd2*r(k)**2)
end do

do i=indtmin,indtm
  distp(i,2)=(3.d0*rmu(3,i)**2-r(i)**2)*cost1d
  distp(i,3)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d
  distp(i,4)=rmu(1,i)*rmu(2,i)*cost3d
  distp(i,5)=rmu(2,i)*rmu(3,i)*cost3d
  distp(i,6)=rmu(1,i)*rmu(3,i)*cost3d
end do

do ic=1,5
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=distp(i,1+ic)*distp(i,1)
  end do
  !           endif
end do

if(typec.ne.1) then


  rp1=r(0)**2
  fun0=distp(0,1)
  fun=(npower-dd2*rp1)*distp(0,1)*2.d0/rp1
  fun2=(npower*(2.d0*npower-1.d0)-                                   &
      (1.d0+4.d0*npower)*dd2*rp1+2.d0*(dd2*rp1)**2)*                 &
      distp(0,1)*2.d0/rp1


  !              indorbp=indorb
  do ic=1,5
    !                 if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)                       &
          *fun
      if(ic.eq.1) then
        if(i.ne.3) then
          z(indorbp,indt+i)=z(indorbp,indt+i)-                       &
              2.d0*rmu(i,0)*fun0*cost1d
        else
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              4.d0*rmu(i,0)*fun0*cost1d
        end if
      elseif(ic.eq.2) then
        if(i.eq.1) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              2.d0*rmu(i,0)*fun0*cost2d
        elseif(i.eq.2) then
          z(indorbp,indt+i)=z(indorbp,indt+i)-                       &
              2.d0*rmu(i,0)*fun0*cost2d
        end if
      elseif(ic.eq.3) then
        if(i.eq.1) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              rmu(2,0)*fun0*cost3d
        elseif(i.eq.2) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              rmu(1,0)*fun0*cost3d
        end if
      elseif(ic.eq.4) then
        if(i.eq.2) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              rmu(3,0)*fun0*cost3d
        elseif(i.eq.3) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              rmu(2,0)*fun0*cost3d
        end if
      elseif(ic.eq.5) then
        if(i.eq.1) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              rmu(3,0)*fun0*cost3d
        elseif(i.eq.3) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
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


