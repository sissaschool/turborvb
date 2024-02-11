! 3d single gaussian


dd1=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd1*r(k)**2)
end do

do i=indtmin,indtm
  distp(i,3)=distp(i,1)
  distp(i,4)=(3.d0*rmu(3,i)**2-r(i)**2)*cost1d
  distp(i,5)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d
  distp(i,6)=rmu(1,i)*rmu(2,i)*cost3d
  distp(i,7)=rmu(2,i)*rmu(3,i)*cost3d
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
  fun=-2.d0*dd1*distp(0,1)
  fun2=((2.d0*dd1*r(0))**2-2.d0*dd1)*distp(0,1)

  !              indorbp=indorb


  do ic=1,5
    !                 if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0)*fun
    end do
    if(ic.eq.1) then
      !                         if(i.ne.3) then
      z(indorbp,indt+1)=z(indorbp,indt+1)-                           &
          2.d0*rmu(1,0)*fun0*cost1d
      z(indorbp,indt+2)=z(indorbp,indt+2)-                           &
          2.d0*rmu(2,0)*fun0*cost1d
      !                         else
      z(indorbp,indt+3)=z(indorbp,indt+3)+                           &
          4.d0*rmu(3,0)*fun0*cost1d
      !                         endif
    elseif(ic.eq.2) then
      !                         if(i.eq.1) then
      z(indorbp,indt+1)=z(indorbp,indt+1)+                           &
          2.d0*rmu(1,0)*fun0*cost2d
      !                         elseif(i.eq.2) then
      z(indorbp,indt+2)=z(indorbp,indt+2)-                           &
          2.d0*rmu(2,0)*fun0*cost2d
      !                         endif
    elseif(ic.eq.3) then
      !                         if(i.eq.1) then
      z(indorbp,indt+1)=z(indorbp,indt+1)+                           &
          rmu(2,0)*fun0*cost3d
      !                         elseif(i.eq.2) then
      z(indorbp,indt+2)=z(indorbp,indt+2)+                           &
          rmu(1,0)*fun0*cost3d
      !                         endif
    elseif(ic.eq.4) then
      !                         if(i.eq.2) then
      z(indorbp,indt+2)=z(indorbp,indt+2)+                           &
          rmu(3,0)*fun0*cost3d
      !                         elseif(i.eq.3) then
      z(indorbp,indt+3)=z(indorbp,indt+3)+                           &
          rmu(2,0)*fun0*cost3d
      !                         endif
    elseif(ic.eq.5) then
      !                         if(i.eq.1) then
      z(indorbp,indt+1)=z(indorbp,indt+1)+                           &
          rmu(3,0)*fun0*cost3d
      !                         elseif(i.eq.3) then
      z(indorbp,indt+3)=z(indorbp,indt+3)+                           &
          rmu(1,0)*fun0*cost3d
      !endif for i
      !                         endif
      !endif for ic
    end if
    !enddo for i
    !                   enddo
    z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun+fun2)
    !endif for iocc
    !                endif
    ! enddo fot ic
  end do

  !endif for indt
end if
!
indpar=indpar+1
indshell=indshell+5
indorb=indorbp

