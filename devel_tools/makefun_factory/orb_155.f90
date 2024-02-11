! Jastrow single gaussian f orbital
! derivative of 154 with respect to z
! unnormalized f orbitals
! R(r)= -r^2*exp(-z r^2)



!        indorbp=indorb
indparp=indpar+1
dd1=dd(indparp)


do k=indtmin,indtm
  distp(k,1)=dexp(-dd1*r(k)**2)
end do


do i=indtmin,indtm
  distp(i,2)=cost1f*rmu(3,i)                                         &
      *(5.d0*rmu(3,i)**2-3.d0*r(i)**2)
  ! lz=0
  distp(i,3)=cost2f*rmu(1,i)                                         &
      *(5.d0*rmu(3,i)**2-r(i)**2)
  ! lz=+/-1
  distp(i,4)=cost2f*rmu(2,i)                                         &
      *(5.d0*rmu(3,i)**2-r(i)**2)
  ! lz=+/-1
  distp(i,5)=cost3f*rmu(3,i)                                         &
      *(rmu(1,i)**2-rmu(2,i)**2)
  ! lz=+/-2
  distp(i,6)=cost3f*2.d0*rmu(3,i)                                    &
      *rmu(1,i)*rmu(2,i)
  ! lz=+/-2
  distp(i,7)=cost4f*rmu(1,i)                                         &
      *(rmu(1,i)**2-3.d0*rmu(2,i)**2)
  ! lz=+/-3
  distp(i,8)=cost4f*rmu(2,i)                                         &
      *(3.d0*rmu(1,i)**2-rmu(2,i)**2)
  ! lz=+/-3
end do


do ic=1,7
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do k=i0,indtm
    z(indorbp,k)=-r(k)**2*distp(k,1)*distp(k,1+ic)
  end do
  !           endif
end do


if(typec.ne.1) then

  dd1=dd(indparp)
  fun0=-r(0)**2*distp(0,1)
  fun=2.d0*(dd1*r(0)**2-1.d0)*distp(0,1)
  fun2=-2.d0*(2.d0*dd1**2*r(0)**4+1.d0                               &
      -5.d0*dd1*r(0)**2)*distp(0,1)

  !              indorbp=indorb
  do ic=1,7
    !                 if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)                       &
          *fun
      if(ic.eq.1) then
        z(indorbp,indt+i)=z(indorbp,indt+i)-                         &
            6.d0*cost1f*fun0*rmu(i,0)*rmu(3,0)
        if(i.eq.3) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              cost1f*fun0*(15.d0*rmu(i,0)**2-3.d0*r(0)**2)
        end if
      elseif(ic.eq.2) then
        z(indorbp,indt+i)=z(indorbp,indt+i)-                         &
            2.d0*cost2f*fun0*rmu(i,0)*rmu(1,0)
        if(i.eq.1) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              cost2f*fun0*(5.d0*rmu(3,0)**2-r(0)**2)
        elseif(i.eq.3) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              10.d0*cost2f*fun0*rmu(i,0)*rmu(1,0)
        end if
      elseif(ic.eq.3) then
        z(indorbp,indt+i)=z(indorbp,indt+i)-                         &
            2.d0*cost2f*fun0*rmu(i,0)*rmu(2,0)
        if(i.eq.2) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              cost2f*fun0*(5.d0*rmu(3,0)**2-r(0)**2)
        elseif(i.eq.3) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              10.d0*cost2f*fun0*rmu(i,0)*rmu(2,0)
        end if
      elseif(ic.eq.4) then
        if(i.eq.1) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              2.d0*cost3f*fun0*rmu(1,0)*rmu(3,0)
        elseif(i.eq.2) then
          z(indorbp,indt+i)=z(indorbp,indt+i)-                       &
              2.d0*cost3f*fun0*rmu(2,0)*rmu(3,0)
        else
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              cost3f*fun0*(rmu(1,0)**2-rmu(2,0)**2)
        end if
      elseif(ic.eq.5) then
        if(i.eq.1) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              2.d0*cost3f*fun0*rmu(2,0)*rmu(3,0)
        elseif(i.eq.2) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              2.d0*cost3f*fun0*rmu(1,0)*rmu(3,0)
        else
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              2.d0*cost3f*fun0*rmu(1,0)*rmu(2,0)
        end if
      elseif(ic.eq.6) then
        if(i.eq.1) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              3.d0*cost4f*fun0*(rmu(1,0)**2-rmu(2,0)**2)
        elseif(i.eq.2) then
          z(indorbp,indt+i)=z(indorbp,indt+i)-                       &
              6.d0*cost4f*fun0*rmu(1,0)*rmu(2,0)
        end if
      else
        if(i.eq.1) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
              6.d0*cost4f*fun0*rmu(1,0)*rmu(2,0)
        elseif(i.eq.2) then
          z(indorbp,indt+i)=z(indorbp,indt+i)+                       &
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



