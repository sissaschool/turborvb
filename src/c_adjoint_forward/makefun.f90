!TL off
subroutine makefun(iopt,indt,i0,indtmin,indtm,typec,indpar               &
           &,indorb,indshell,nelskip,z,dd,zeta,r,rmu,distp               &
           &,iflagnorm_unused,cr)
    use constants
    implicit none
    integer iopt,indt,i,k,nelskip,indpar,indorbp                         &
           &,indorb,indshell,indshellp,ic,indtmin,i0                     &
           &,iflagnorm_unused,indparp,indtm,npower,typec                 &
           &,ii,jj,kk
    real*8 z(nelskip,i0:*),dd(*),zeta(*),rmu(3,0:indtm)                  &
           &,r(0:indtm)                                                  &
           &,distp(0:indtm,20),peff,fun,fun0,fun2                        &
           &,rp1,rp2,rp3,rp4,rp5,rp6,rp7,rp8                             &
           &,dd1,dd2,dd3,dd4,dd5,c,cr,funp,fun2p,funb                    &
           &,peff2,arg,c0,c1,cost,zv(6),yv(6),xv(6),r2,r4,r6 ! up to i

    integer :: count, multiplicity
    integer, parameter :: max_power = 20
    real*8 :: powers(3,0:max_power,0:indtm)
  !
  ! indorb are the number of orbitals occupied before calling
  ! this subroutine
  !
  ! indpar is the number of variational parameters used
  ! before calling this subroutine
  !
  ! indshell is the index of the last  occupied orbital
  ! in the shell, characterized by occupation number iocc(indshell)
  !
  ! z(i,indt+4)   contains the laplacian of the orbital i
  ! z(i,indt+mu) contains the gradient of the orbital i (mu=1,2,3)
  ! In the following given a radial part of the orbital f(r)
  ! fun=1/r  d f(r)/d r

  !print *,__FILE__
  !print *,'makefun: iopt=',iopt
  !print *,'makefun: i=',i0,' a=',indtmin,' b=',indtm
  !print *,'makefun: indpar=',indpar,' indorb=',indorb,' indshell=',indshell
  !print *,'makefun: nelskip=',nelskip
select case (iopt)
case (105)
!     2s double gaussian without  constant
!     (exp (-dd2 r^2)+dd4*exp(-dd5*r^2))



!        dd1=1.d0
dd2=dd(indpar+1)
!        dd3=dd(indpar+2)
!        dd4=dd(indpar+3)
!        dd5=dd(indpar+4)
dd4=dd(indpar+2)
dd5=dd(indpar+3)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k)*r(k))
  distp(k,2)=dexp(-dd5*r(k)*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)+dd4*distp(i,2)
!              write(6,*) ' function inside = ',z(indorbp,i)
end do
!           endif


if(typec.ne.1) then
  fun=-2.d0*(dd2*distp(0,1)+dd5*dd4*distp(0,2))
  fun2=r(0)**2

  !              write(6,*) ' fun inside = ',fun,fun2

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*(dd2*(-3.d0+2.d0*dd2*fun2)*                 &
      distp(0,1)+dd5*dd4*(-3.d0+2.d0*dd5*fun2)*distp(0,2))




  !endif for indt
end if

indpar=indpar+3
indshell=indshellp
indorb=indorbp

case (40)
  !      3p without cusp condition derivative of 20
  !      r e^{-z1 r }



  dd1=dd(indpar+1)
  !        if(iflagnorm.gt.2) then
  !        c=dsqrt((2.d0*dd1)**5/8.d0/pi)/2.d0
  c=dd1**2.5d0*0.5641895835477562d0
  !        endif

  c0=-c
  c1=2.5d0*c/dd1

  !
  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
    distp(k,2)=r(k)*distp(k,1)
  end do
  !
  !        indorbp=indorb
  !
  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*(c0*distp(i,2)+c1*distp(i,1))
    end do
    !           endif
  end do
  !
  !
  if(typec.ne.1) then
    fun0=c0*distp(0,2)+c1*distp(0,1)
    fun=(c0*(1.d0-dd1*r(0))-c1*dd1)*distp(0,1)
    fun2=(c0*dd1*(dd1*r(0)-2.d0)+c1*dd1**2)*distp(0,1)
    !
    !              indorbp=indorb
    !
    do ic=1,3
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun/r(0)
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*                                   &
          (4.d0*fun/r(0)+fun2)
      !
      !                endif
    end do
    !
    !
    !endif for indt
  end if
  !
  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp


  ! 4p single zeta
case (52)
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




case (31)
  !     3d without cusp condition double Z


  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd(indpar+3)
  !        if(iflagnorm.gt.2) then
  c=1/2.d0*dsqrt(5.d0/pi)                                            &
      /dsqrt(1/dd1**7/128.d0+2*peff/(dd1+dd2)**7                     &
      +peff**2/dd2**7/128.d0)/dsqrt(720.d0)
  !        endif

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
    distp(k,2)=dexp(-dd2*r(k))
  end do

  do i=indtmin,indtm
    distp(i,3)=c*(distp(i,1)+peff*distp(i,2))
    !lz=0
    distp(i,4)=(3.d0*rmu(3,i)**2-r(i)**2)*cost1d
    !lz=+/-2
    distp(i,5)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d
    ! lz=+/- 2
    distp(i,6)=rmu(1,i)*rmu(2,i)*cost3d
    ! lz=+/-1
    distp(i,7)=rmu(2,i)*rmu(3,i)*cost3d
    ! lz=+/-1
    distp(i,8)=rmu(1,i)*rmu(3,i)*cost3d
  end do

  !        indorbp=indorb

  do ic=1,5
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=distp(i,3+ic)*distp(i,3)
    end do
    !           endif
  end do

  if(typec.ne.1) then
    fun0=distp(0,3)
    fun=c*(-dd1*distp(0,1)-peff*dd2*distp(0,2))
    fun2=c*(dd1**2*distp(0,1)                                        &
        +peff*dd2**2*distp(0,2))

    !              indorbp=indorb

    do ic=1,5
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0)*fun/r(0)
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
      z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0)+fun2)
      !endif for iocc
      !                endif
      ! enddo fot ic
    end do


    !endif for indt
  end if

  indpar=indpar+3
  indshell=indshell+5
  indorb=indorbp


case (113)
!     2s without cusp condition
!     dd1*( dd3 +r^2/(1+dd2*r)^4)


dd2=dd(indpar+1)
dd3=dd(indpar+2)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=r(k)**2/(1.d0+dd2*r(k))**4
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)+dd3
end do
!           endif


if(typec.ne.1) then
  fun= (2.d0-2.d0*dd2*r(0))/(1+dd2*r(0))**5
  fun2=2.d0*(1.d0-6.d0*dd2*r(0)+3.d0*(dd2*r(0))**2)                  &
      /(1+dd2*r(0))**6
  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=fun2+2.d0*fun


  !endif for indt
end if

indpar=indpar+2
indshell=indshellp
indorb=indorbp

case (10000:11000)
  ! Reserved for dummy orbitals

  indpar=indpar+0
  indshell=indshell+(iopt - 10000)
  indorb=indorbp + (iopt - 10000)
case (107)
!     2p single  lorentian  parent of 103



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=1.d0/(1.d0+dd2*r(k)**2)
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)
  end do
  !           endif
end do


if(typec.ne.1) then
  fun0=distp(0,1)
  fun=-dd2*distp(0,1)**2*2.d0
  fun2=fun*distp(0,1)*(1.d0-3.d0*dd2*r(0)**2)

  !              indorbp=indorb

  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do


  !endif for indt
end if

indpar=indpar+1
indshell=indshell+3
indorb=indorbp

case (43)
  !     4d without cusp and one parmater  derivative of 33


  dd1=dd(indpar+1)
  !         if(iflagnorm.gt.2) then
  !     c=                                                           &
      !                                                              &1.d0/(2.d0**3*3.d0)/dsqrt(56.d0*pi)*(2.d0*dd1)**(9.d0/2.d0)
      c=dd1**4.5d0*0.0710812062076410d0
  !         endif

  c0=-c
  c1=4.5d0*c/dd1


  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
  end do

  do i=indtmin,indtm
    distp(i,3)=distp(i,1)*(c0*r(i)**2+c1*r(i))
    ! lz=0
    distp(i,4)=(3.d0*rmu(3,i)**2-r(i)**2)*cost1d
    ! lz=+/
    distp(i,5)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d
    ! lz=+/-2
    distp(i,6)=rmu(1,i)*rmu(2,i)*cost3d
    ! lz=+/-1
    distp(i,7)=rmu(2,i)*rmu(3,i)*cost3d
    ! lz=+/-1
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
    fun=-dd1*distp(0,3)+distp(0,1)*(2.d0*c0*r(0)+c1)
    fun2=dd1**2*distp(0,3)+distp(0,1)*                               &
        (-2.d0*dd1*(2.d0*c0*r(0)+c1)+2.d0*c0)
    !              indorbp=indorb

    do ic=1,5
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0)                     &
            *fun/r(0)
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
      z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0)+fun2)
      !endif for iocc
      !                endif
      ! enddo fot ic
    end do
  end if
  !
  indpar=indpar+1
  indshell=indshell+5
  indorb=indorbp


  ! derivative of 36 with respect zeta
case (6)
  ! normalized

  indshellp=indshell+1


  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd(indpar+3)

  !           if(iflagnorm.gt.2) then
  !            c=                                              WRONG
  !                                                                  &0.5*dsqrt((1.d0/dd1**5/32.d0+2.d0*peff/(dd1+dd2)**5
      !                                                              &+peff**2/dd2**5/32.d0)/(24.d0*pi))

  c=1.d0/dsqrt((3.d0*pi)*                                            &
      (1.d0/dd1**5+ 64.d0*peff/(dd1+dd2)**5+peff**2/dd2**5))

  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=c*dexp(-dd2*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=r(i)*(distp(i,1)+distp(i,2)*peff)
  end do

  if(typec.ne.1) then

    fun=distp(0,1)*(1.d0-dd1*r(0))                                   &
        +peff*distp(0,2)*(1.d0-dd2*r(0))
    fun2=distp(0,1)*(dd1**2*r(0)-2.d0*dd1)+peff*distp(0,2)           &
        *(dd2**2*r(0)-2.d0*dd2)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)/r(0)
    end do
    z(indorbp,indt+4)=(2.d0*fun/r(0)+fun2)

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+3
  indshell=indshellp

  ! 2s double Z NO CUSP
case (136)
!     2p single exponential  r^5 e^{-z r}  !



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)*r(i)**5
  end do
  !           endif
end do

if(typec.ne.1) then

  fun0=distp(0,1)*r(0)**5
  fun=distp(0,1)*(5.d0-dd2*r(0))*r(0)**3
  fun2=distp(0,1)*(20*r(0)**3-10*dd2*r(0)**4                         &
      +dd2**2*r(0)**5)
  !              indorbp=indorb
  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do
  !endif for indt
end if
indpar=indpar+1
indshell=indshell+3
indorb=indorbp

case (200)
!     THE  COSTANT

indorbp=indorb+1
indshellp=indshell+1

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=1.d0
end do
!           endif

if(typec.ne.1) then
  do i=1,3
    z(indorbp,indt+i)=0
  end do

  z(indorbp,indt+4)=0
  !endif for indt
end if

indshell=indshellp
indorb=indorbp


case (118)
!     2s double lorentian with constant  parent of 102
!     (dd1+  1/ (1 + Exp[  dd2 (r^2 - r_0^2) ] )   | dd3=r_0
!      Fermi distribution with r^2


dd1=dd(indpar+1)
dd2=dd(indpar+2)
dd3=-dd2*dd(indpar+3)**2

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  arg=dd2*r(k)**2+dd3
  if(arg.gt.200) then
    distp(k,1)=dexp(200.d0)
  else
    distp(k,1)=dexp(arg)
  end if
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=dd1+1.d0/(1.d0+distp(i,1))
!              write(6,*) ' function inside = ',z(indorbp,i)
end do
!           endif


if(typec.ne.1) then

  fun= -2.d0*dd2*distp(0,1)/(1.d0+distp(0,1))**2
  fun2=-2.d0*dd2*(-distp(0,1)*(-1.d0-2.d0*dd2*r(0)**2)               &
      +distp(0,1)**2*(1.d0-2.d0*dd2*r(0)**2))/(1.d0+distp(0,1))**3


  !              write(6,*) ' fun inside = ',fun,fun2

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=fun2+2.d0*fun

  !           write(6,*) ' lap 106 =',z(indorbp,indt+4)

  !endif for indt
end if

indpar=indpar+3
indshell=indshellp
indorb=indorbp

case (15)
  ! (r**2 + dd2*(1 + dd1*r))*exp(-dd1*r)  ! normalized


  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  dd2=dd(indpar+2)

  c=dsqrt(2.d0*dd1**7/pi/                                            &
      (45.d0+42.d0*dd1**2*dd2+14.d0*dd1**4*dd2**2))

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
  end do


  do i=i0,indtm
    z(indorbp,i)=(r(i)**2+dd2*(1.d0+dd1*r(i)))                       &
        *distp(i,1)
  end do

  if(typec.ne.1) then

    fun=distp(0,1)*r(0)*(2.d0-dd1**2*dd2-dd1*r(0))
    fun2=distp(0,1)*((1.d0-dd1*r(0))                                 &
        *(3.d0-dd1**2*dd2-dd1*r(0))-1.d0)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)/r(0)
    end do
    z(indorbp,indt+4)=2.d0*fun/r(0)+fun2

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+2
  indshell=indshellp

  ! 2s gaussian for pseudo
case (122)
!     2s with cusp condition
!     dd1*( dd3 +exp(-dd2*r)*(1+dd2*r))


dd2=dd(indpar+1)
dd3=dd(indpar+2)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)*(1.d0+dd2*r(i))+dd3
end do
!           endif


if(typec.ne.1) then
  fun=-dd2**2*distp(0,1)
  fun2=fun*(1.d0-dd2*r(0))
  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do
  z(indorbp,indt+4)=2.d0*fun+fun2
  !endif for indt
end if

indpar=indpar+2
indshell=indshellp
indorb=indorbp

case (128)
!     2s with cusp condition
!     ( r^2*exp(-dd2*r))  ! with no cusp condition


dd2=dd(indpar+1)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)*r(i)**2
end do
!           endif


if(typec.ne.1) then
  fun=(2.d0-dd2*r(0))*distp(0,1)
  fun2=(2.d0-4*dd2*r(0)+(dd2*r(0))**2)*distp(0,1)
  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*fun+fun2


  !endif for indt
end if

indpar=indpar+1
indshell=indshellp
indorb=indorbp

case (16)
  ! R(r)=exp(-z*r**2) single zeta

  indshellp=indshell+1
  indorbp=indorb+1
  dd1=dd(indpar+1)

  if(dd1.ne.0.) then
    ! c=(2.d0*dd1/pi)**(3.d0/4.d0)
    c=0.71270547035499016d0*dd1**0.75d0
  else
    c=1.d0
  end if

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)
  end do

  if(typec.ne.1) then
    ! the first derivative /r
    fun=-2.d0*dd1*distp(0,1)

    ! the second derivative
    fun2=fun*(1.d0-2.d0*dd1*r(0)*r(0))

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*fun+fun2

    if(typec.eq.2) then
      funb=(fun2-fun)/(r(0)*r(0))

      z(indorbp,indt+5)=funb*rmu(1,0)*rmu(1,0)+fun
      z(indorbp,indt+6)=funb*rmu(2,0)*rmu(2,0)+fun
      z(indorbp,indt+7)=funb*rmu(3,0)*rmu(3,0)+fun
      z(indorbp,indt+8)=funb*rmu(1,0)*rmu(2,0)
      z(indorbp,indt+9)=funb*rmu(1,0)*rmu(3,0)
      z(indorbp,indt+10)=funb*rmu(2,0)*rmu(3,0)

    end if


  end if

  indorb=indorbp

  indpar=indpar+1
  indshell=indshellp

case (2200:2299)
!     d gaussian  -r**(2*(npower+1))*exp(-alpha*r**2) derivative 1200

npower=iopt+1-2200

!        indorbp=indorb

dd2=dd(indpar+1)
do k=indtmin,indtm
  distp(k,1)=-r(k)**(2*npower)*dexp(-dd2*r(k)**2)
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


case (55)
  ! g single Slater orbital
  ! R(r)= exp(-alpha r)
  ! normalized


  !        indorbp=indorb
  indparp=indpar+1

  dd1=dd(indparp)

  !        if(iflagnorm.gt.2) then
  ! overall normalization
  ! l = 4
  ! \int d\omega Y*Y = 4 pi / (2 l + 1)
  ! \int dr r^{2 l + 2} Exp [- 2 dd1 r^2 ] =  7 * 5**2 * 3**4 / 2**3 / dd1**11
  !        c=1.d0/dsqrt(7.d0)*(2.d0/pi)**(1.d0/2.d0)*dd1**(11.d0/2.d0)/3.d0/5.d0
  c=dd1**5.5d0*.020104801169736915d0
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
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
      z(indorbp,k)=distp(k,1)*distp(k,1+ic)
    end do
    !           endif
  end do


  if(typec.ne.1) then

    dd1=dd(indparp)
    fun0=distp(0,1)
    fun=-dd1*distp(0,1)/r(0)
    fun2=dd1**2*distp(0,1)


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


case (2)
  !



  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1

  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=(zeta(1)-dd1)/(dd2-zeta(1))

  !           if(iflagnorm.gt.2) then
  c=1.d0/2.d0/dsqrt(2.d0*pi*(1.d0/(2.d0*dd1)**3                      &
      +2.d0*peff/(dd1+dd2)**3+peff**2/(2.d0*dd2)**3))
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=c*dexp(-dd2*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)+peff*distp(i,2)
  end do

  if(typec.ne.1) then
    fun=(-dd1*distp(0,1)-dd2*distp(0,2)*peff)/r(0)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=(-2.d0*dd1/r(0)+dd1**2)                        &
        *distp(0,1)+peff*(-2.d0*dd2/r(0)+dd2**2)                     &
        *distp(0,2)


  end if

  indorb=indorbp

  !        endif

  indpar=indpar+2
  indshell=indshellp


  ! 1s double Z NO CUSP
case (23)
  !      3p without cusp condition
  !      r ( e^{-z2 r } + z1 e^{-z3 r } )



  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  dd3=dd(indpar+3)
  !        if(iflagnorm.gt.2) then
  c=1.d0/2.d0/dsqrt(240.d0*pi*(1.d0/(2.d0*dd1)**7                    &
      +2.d0*dd3/(dd1+dd2)**7+dd3**2/(2.d0*dd2)**7))
  !        endif
  !
  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=c*dexp(-dd2*r(k))
  end do
  !
  do i=indtmin,indtm
    distp(i,3)=r(i)*(distp(i,1)+dd3*distp(i,2))
  end do
  !
  !        indorbp=indorb
  !
  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,3)
    end do
    !           endif
  end do
  !
  !
  if(typec.ne.1) then
    fun0=distp(0,3)
    fun=(1.d0-dd1*r(0))*distp(0,1)                                   &
        +dd3*(1.d0-dd2*r(0))*distp(0,2)
    fun2=dd1*(dd1*r(0)-2.d0)*distp(0,1)                              &
        +dd3*dd2*(dd2*r(0)-2.d0)*distp(0,2)
    !
    !              indorbp=indorb
    !
    do ic=1,3
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun/r(0)
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*                                   &
          (4.d0*fun/r(0)+fun2)
      !
      !                endif
    end do
    !
    !
    !endif for indt
  end if
  !
  indpar=indpar+3
  indshell=indshell+3
  indorb=indorbp



  ! 4p single zeta
case (80)
  ! R(r)=exp(-z*r**2) single zeta

  indshellp=indshell+1
  indorbp=indorb+1

  dd1=dd(indpar+1)
  dd2=dsqrt(dd1)

  !           if(iflagnorm.gt.2) then
  !           c=(2.d0*dd1/pi)**(3.d0/4.d0)*ratiocs
  !           ratiocs--> ratiocs*(2/pi)**3/4
  c=dd1**0.75d0*ratiocs
  !           endif

  do k=indtmin,indtm
    cost=dd1*r(k)**2/(1.d0+dd2*r(k))
    distp(k,1)=c*dexp(-cost)
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)
  end do

  if(typec.ne.1) then
    ! the first derivative /r
    rp1=dd1*r(0)**2
    rp2=dd2*r(0)
    rp3=(1.d0+rp2)**2

    fun=-dd1*distp(0,1)*(2.d0+rp2)/rp3

    ! the second derivative
    fun2=dd1*distp(0,1)*(-2.d0-2.d0*rp2+4.d0*rp1+4.d0*rp1*rp2+rp1**2)/rp3**2

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp
  indpar=indpar+1
  indshell=indshellp

case (17)
  ! 2s gaussian for pseudo
  ! R(r)=r**2*exp(-z*r**2) single zeta

  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)

  !           if(iflagnorm.gt.2) then
  !        c=4.d0*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)/dsqrt(15.d0)
  c=.73607904464954686606d0*dd1**1.75d0
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)*r(i)**2
  end do

  if(typec.ne.1) then
    rp1=r(0)**2
    !              the first derivative / r
    fun=2.d0*distp(0,1)*(1.d0-dd1*rp1)
    !              the second derivative
    fun2=2.d0*distp(0,1)*(1.d0-5.d0*dd1*rp1+2.d0*dd1**2*rp1**2)
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  !        endif

  indpar=indpar+1
  indshell=indshellp

  ! 2s gaussian for pseudo
case (10)
  ! R(r)=r**2*exp(-z1*r)



  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  !           if(iflagnorm.gt.2) then
  !           c=dsqrt((2*dd1)**7/720.d0/pi)/2.d0
  c=dd1**3.5d0*0.11894160774351807429d0
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)*r(i)**2
  end do

  if(typec.ne.1) then
    fun=(2.d0-dd1*r(0))*distp(0,1)
    fun2=(2.d0-4*dd1*r(0)+(dd1*r(0))**2)*distp(0,1)
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  !        endif
  indpar=indpar+1
  indshell=indshellp


  ! 3s double zeta
case (129)
!     2p single exponential  r e^{-z r}  ! parent of 121



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)*r(i)
  end do
  !           endif
end do


if(typec.ne.1) then

  fun0=distp(0,1)*r(0)
  fun=distp(0,1)*(1.d0-dd2*r(0))/r(0)
  fun2=dd2*(dd2*r(0)-2.d0)*distp(0,1)
  !              indorbp=indorb
  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do
  !endif for indt
end if
indpar=indpar+1
indshell=indshell+3
indorb=indorbp

case (110)
!     2s without cusp condition
!     dd1*( dd3 +1/(1+dd2*r^3))


dd2=dd(indpar+1)
dd3=dd(indpar+2)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=1.d0/(1.d0+dd2*r(k)**3)
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)+dd3
end do
!           endif


if(typec.ne.1) then
  fun=-dd2*distp(0,1)**2*3.d0*r(0)
  fun2=fun*distp(0,1)*(2.d0-4.d0*dd2*r(0)**3)

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=fun2+2.d0*fun


  !endif for indt
end if

indpar=indpar+2
indshell=indshellp
indorb=indorbp

case (46)
  ! R(r)=c*r**2*exp(-z*r**2)*(7/4/z-r**2)


  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)

  !           if(iflagnorm.gt.2) then
  c=4.d0*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)/dsqrt(15.d0)
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)*(7.d0/4.d0/dd1*r(i)**2                   &
        -r(i)**4)
  end do

  if(typec.ne.1) then
    rp1=r(0)**2
    !              the first derivative / r
    fun=distp(0,1)*(7.d0-15.d0*dd1*rp1                               &
        +4.d0*(dd1*rp1)**2)/2.d0/dd1
    !              the second derivative
    fun2=distp(0,1)*(7.d0-59*dd1*rp1+50*(dd1*rp1)**2                 &
        -8*(dd1*rp1)**3)/2.d0/dd1
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  !        endif

  indpar=indpar+1
  indshell=indshellp

  ! 5s single zeta derivative of 12
case (143)
!     4d  one parmater der of 133


dd1=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd1*r(k))
end do

do i=indtmin,indtm
  distp(i,3)=distp(i,1)*r(i)**2
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
    z(indorbp,i)=-distp(i,3+ic)*distp(i,3)
  end do
  !            endif
end do

if(typec.ne.1) then
  fun0=-distp(0,3)
  fun=-(-2.d0+dd1*r(0))*distp(0,1)
  fun2=((dd1*r(0))**2 -4.d0*r(0)*dd1+2.d0)*distp(0,1)

  !              indorbp=indorb

  do ic=1,5
    !                 if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0)                       &
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

case (7)
  ! normalized IS WRONG!!!

  indshellp=indshell+1


  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd(indpar+3)

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
    distp(k,2)=dexp(-dd2*r(k))
  end do
  !           if(iflagnorm.gt.2) then
  c=                                                                 &
      1/dsqrt(1/(3.D0/4.D0/dd1**5+peff**2/dd2**3/4+12*peff/          &
      (dd1+dd2)**4))*1.d0/dsqrt(4.0*pi)
  !           endif

  do i=i0,indtm
    z(indorbp,i)=c*(distp(i,1)+r(i)*distp(i,2)*peff)
  end do

  if(typec.ne.1) then

    fun=-dd1*distp(0,1)+peff*distp(0,2)*(1.d0-dd2*r(0))
    fun2=distp(0,1)*dd1**2                                           &
        +peff*distp(0,2)*(dd2**2*r(0)-2.d0*dd2)

    do i=1,3
      z(indorbp,indt+i)=fun*c*rmu(i,0)/r(0)
    end do
    z(indorbp,indt+4)=c*(2.d0*fun/r(0)+fun2)

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+3
  indshell=indshellp



  ! 2s double Z WITH CUSP
case (36)



  dd1=dd(indpar+1)

  !\print *, "i0, indtmin, indtm"
  !\print *,i0, indtmin, indtm


  !        if(iflagnorm.gt.2) then
  !        c=dsqrt(2.d0)*pi**(-0.75d0)*(2.d0*dd1)**1.25d0
  c=dd1**1.25d0*1.42541094070998d0
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  !        indorbp=indorb
  !
  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,1)
    end do
    !           endif
  end do

  if(typec.ne.1) then
    fun0=distp(0,1)
    fun=-2.d0*dd1*distp(0,1)
    fun2=fun*(1.d0-2.d0*dd1*r(0)**2)

    !              indorbp=indorb

    do ic=1,3
      !                if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun
      end do
      z(indorbp,indt+ic)=z(indorbp,indt+ic)+fun0
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
      !                 endif
    end do
    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp

case (29)
  ! derivative of (28)


  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1

  dd1=dd(indpar+1)

  !           if(iflagnorm.gt.2) then
  !             if(dd1.gt.0.) then
  c=cost1s*dd1**1.5d0
  !             else
  !               c=1.d0
  !             endif
  !           endif
  !             if(dd1.gt.0.) then
  c1=1.5d0/dd1
  !             else
  !               c1=0.d0
  !             endif

  do i=indtmin,indtm
    distp(i,1)=c*dexp(-dd1*r(i))
  end do

  do i=i0,indtm
    !              rp1=(b1s*r(i))**4*dd1**3
    !              rp4=rp1*dd1
    !              rp5=dd1*r(i)
    !              z(indorbp,i)=distp(i,1)*                          &
        !                                                            &(c1*rp4/(1+rp4)-rp1*(-4+rp5+rp4*rp5)/(1+rp4)**2)
        rp4=(b1s*dd1*r(i))**4
    rp5=dd1*r(i)
    z(indorbp,i)=distp(i,1)*rp4/(1+rp4)*                             &
        (c1 - (1.d0/dd1)*(-4+rp5+rp4*rp5)/(1+rp4))
  end do

  if(typec.ne.1) then
    rp1=dd1*b1s*r(0)
    rp2=rp1**2
    rp4=rp2**2
    rp5=rp4*rp1
    rp8=rp4*rp4

    fun=distp(0,1)* (dd1*rp2*(4*b1s**2*(11-5*rp4) +2*(rp1+rp5)**2    &
        -b1s*rp1*(21+26*rp4+5*rp8)))/(2.*(1+rp4)**3)

    fun2=distp(0,1)*(dd1*rp2*(b1s*(31 + 7*rp4)*(rp1 + rp5)**2        &
        - 2*(rp1 + rp5)**3 + 64*b1s**2*rp1*(-2 - rp4 + rp8) +        &
        4*b1s**3*(33 - 134*rp4 + 25*rp8)))/(2.*b1s*(1 + rp4)**4)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*fun+fun2

    !           endif

    indorb=indorbp

  end if

  indpar=indpar+1
  indshell=indshellp





case (44)
  ! R(r)=x*exp(-z*r**2)*(5/4/z-r**2)


  dd1=dd(indpar+1)
  !        if(iflagnorm.gt.2) then
  !        c=dsqrt(2.d0)*pi**(-0.75d0)*(2.d0*dd1)**1.25d0
  c=dd1**1.25d0*1.42541094070998d0
  !        endif
  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  !        indorbp=indorb
  !
  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,1)*                             &
          (5.d0/4.d0/dd1-r(i)**2)
    end do
    !           endif
  end do

  if(typec.ne.1) then
    fun0=distp(0,1)*(5.d0/4.d0/dd1-r(0)**2)
    fun=distp(0,1)*(2.d0*dd1*r(0)**2-9.d0/2.d0)
    fun2=distp(0,1)*(-4.d0*dd1**2*r(0)**4                            &
        +15.d0*dd1*r(0)**2-9.d0/2.d0)

    !              indorbp=indorb

    do ic=1,3
      !                if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
      !                 endif
    end do


    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp


  ! derivative of 37 with respect to z
case (64)
  ! d orbitals
  ! R(r)= r exp(-alpha r^2)
  ! each gaussian term is normalized



  !        indorbp=indorb
  indparp=indpar+1

  dd1=dd(indparp)

  !        if(iflagnorm.gt.2) then
  !        c=8.d0/dsqrt(21.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(9.d0/4.d0)
  c=dd1**2.25d0*1.24420067280413253d0
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
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
      z(indorbp,k)=distp(k,1)*distp(k,1+ic)*r(k)
    end do
    !           endif
  end do


  if(typec.ne.1) then

    dd1=dd(indparp)

    rp1=2.d0*dd1*r(0)
    rp2=rp1*r(0)
    fun0=distp(0,1)*r(0)
    fun=(1.d0-rp2)*distp(0,1)/r(0)
    fun2=distp(0,1)*rp1*(rp2-3.d0)


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
  end if

  indpar=indpar+1
  indshell=indshell+5
  indorb=indorbp


case (106)
!     2s without cusp condition
!     dd1*( dd3 +1/(1+dd2*r^2))


dd2=dd(indpar+1)
dd3=dd(indpar+2)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=1.d0/(1.d0+dd2*r(k)*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)+dd3
end do
!           endif


if(typec.ne.1) then
  fun=-dd2*distp(0,1)**2*2.d0
  fun2=fun*distp(0,1)*(1.-3.d0*dd2*r(0)**2)

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=fun2+2.d0*fun


  !endif for indt
end if

indpar=indpar+2
indshell=indshellp
indorb=indorbp

case (71)
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
case (14)
  ! (1.d0 + dd1 r) * exp(-dd1 * r)   ! normalized


  indshellp=indshell+1


  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
  end do

  !           if(iflagnorm.gt.2) then
  !           c=dsqrt(dd1**3.d0/7.d0/pi)
  c=dd1**1.5d0*0.213243618622923d0
  !           endif

  do i=i0,indtm
    z(indorbp,i)=c*(1.d0+dd1*r(i))*distp(i,1)
  end do

  if(typec.ne.1) then
    fun=-distp(0,1)*dd1**2*r(0)
    fun2=-distp(0,1)*dd1**2*(1.d0-dd1*r(0))
    do i=1,3
      z(indorbp,indt+i)=c*fun*rmu(i,0)/r(0)
    end do
    z(indorbp,indt+4)=c*2.d0*fun/r(0)+c*fun2
  end if

  indorb=indorbp

  !        endif
  indpar=indpar+1
  indshell=indshellp



  ! 1s single Z pseudo
case (60)
  ! R(r)=r**3*exp(-z*r**2) single zeta



  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)

  !           if(iflagnorm.gt.2) then
  !     c=2.d0/pi**(3.d0/4.d0)*(2.d0*dd1)**(9.d0/4.d0)*dsqrt(2.d0/105.d0)
  c=dd1**2.25d0*.55642345640820284397d0
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)*r(k)
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)*r(i)**2
  end do

  if(typec.ne.1) then
    rp1=r(0)**2*dd1
    !              the first derivative / r
    fun=distp(0,1)*(3.d0-2.d0*rp1)
    !              the second derivative
    fun2=distp(0,1)*(6.d0-14.d0*rp1+4.d0*rp1**2)
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  !        endif

  indpar=indpar+1
  indshell=indshellp

  ! 3s -derivative of 60 with respect to dd1
case (19)
  ! R(r)=c*exp(-z*r**2)*(3/4/z-r**2)




  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1

  dd1=dd(indpar+1)

  !         if(iflagnorm.gt.2) then
  !          if(dd1.ne.0.) then
  !            c=(2.d0*dd1/pi)**(3.d0/4.d0)
  c=0.71270547035499016d0*dd1**0.75d0
  !          else
  !          c=1.d0
  !          endif
  !         endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)*(3.d0/4.d0/dd1-r(i)**2)
  end do

  if(typec.ne.1) then
    !              the first derivative /r
    fun=distp(0,1)*(2.d0*dd1*r(0)**2-7.d0/2.d0)

    !              the second derivative
    fun2=distp(0,1)*(-4.d0*dd1**2*r(0)**4                            &
        +13.d0*dd1*r(0)**2-7.d0/2.d0)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*fun+fun2

  end if

  indorb=indorbp

  !        endif

  indpar=indpar+1
  indshell=indshellp



  ! 2p single zeta
case (51)
  ! g single gaussian orbital
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
    do ic=1,9
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)*fun
      end do
      if(ic.eq.1) then
        !                      if(i.eq.1) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost1g*fun0*(-60.d0*rmu(1,0)*rmu(3,0)**2+12.d0*rmu(1,0)*r(0)**2)
        !                      elseif(i.eq.2) then
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost1g*fun0*(-60.d0*rmu(2,0)*rmu(3,0)**2+12.d0*rmu(2,0)*r(0)**2)
        !                      else
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost1g*fun0*(80.d0*rmu(3,0)**3-48.d0*rmu(3,0)*r(0)**2)
        !                      endif
      elseif(ic.eq.2) then
        !                      if(i.eq.1) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost2g*fun0*(-9.d0*rmu(1,0)**2*rmu(3,0)-3.d0*rmu(2,0)**2*rmu(3,0)+4.d0*rmu(3,0)**3)
        !                      elseif(i.eq.2) then
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost2g*fun0*(-6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0))
        !                      else
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost2g*fun0*(-3.d0*rmu(1,0)*(rmu(1,0)**2+rmu(2,0)**2-4.d0*rmu(3,0)**2))
        !                      endif
      elseif(ic.eq.3) then
        !                      if(i.eq.1) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost2g*fun0*(-6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0))
        !                      elseif(i.eq.2) then
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost2g*fun0*(-3.d0*rmu(1,0)**2*rmu(3,0)-9.d0*rmu(2,0)**2*rmu(3,0)+4.d0*rmu(3,0)**3)
        !                      else
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost2g*fun0*(-3.d0*rmu(2,0)*(rmu(1,0)**2+rmu(2,0)**2-4.d0*rmu(3,0)**2))
        !                      endif
      elseif(ic.eq.4) then
        !                      if(i.eq.1) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost3g*fun0*(-4.d0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(3,0)**2))
        !                      elseif(i.eq.2) then
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost3g*fun0*(4.d0*(rmu(2,0)**3-3.d0*rmu(2,0)*rmu(3,0)**2))
        !                      else
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost3g*fun0*(12.d0*(rmu(1,0)**2-rmu(2,0)**2)*rmu(3,0))
        !                      endif
      elseif(ic.eq.5) then
        !                      if(i.eq.1) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost3g*fun0*(-2.d0*rmu(2,0)*(3.d0*rmu(1,0)**2+rmu(2,0)**2-6.d0*rmu(3,0)**2))
        !                      elseif(i.eq.2) then
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost3g*fun0*(-2.d0*rmu(1,0)*(rmu(1,0)**2+3.d0*rmu(2,0)**2-6.d0*rmu(3,0)**2))
        !                      else
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost3g*fun0*24.d0*rmu(1,0)*rmu(2,0)*rmu(3,0)
        !                      endif
      elseif(ic.eq.6) then
        !                      if(i.eq.1) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost4g*fun0*3.d0*(rmu(1,0)**2-rmu(2,0)**2)*rmu(3,0)
        !                      elseif(i.eq.2) then
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            -cost4g*fun0*6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0)
        !                      else
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost4g*fun0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(2,0)**2)
        !                      endif
      elseif(ic.eq.7) then
        !                      if(i.eq.1) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost4g*fun0*6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0)
        !                      elseif(i.eq.2) then
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost4g*fun0*3.d0*(rmu(1,0)**2-rmu(2,0)**2)*rmu(3,0)
        !                      else
        z(indorbp,indt+3)=z(indorbp,indt+3)                          &
            +cost4g*fun0*(3.d0*rmu(1,0)**2*rmu(2,0)-rmu(2,0)**3)
        !                      endif
      elseif(ic.eq.8) then
        !                      if(i.eq.1) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost5g*fun0*4.d0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(2,0)**2)
        !                      elseif(i.eq.2) then
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost5g*fun0*4.d0*(rmu(2,0)**3-3.d0*rmu(1,0)**2*rmu(2,0))
        !                      endif
      elseif(ic.eq.9) then
        !                      if(i.eq.1) then
        z(indorbp,indt+1)=z(indorbp,indt+1)                          &
            +cost5g*fun0*4.d0*(3.d0*rmu(1,0)**2*rmu(2,0)-rmu(2,0)**3)
        !                      elseif(i.eq.2) then
        z(indorbp,indt+2)=z(indorbp,indt+2)                          &
            +cost5g*fun0*4.d0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(2,0)**2)
        !                      endif
      end if
      !enddo for i
      !                   enddo
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


case (142)
!     4d  one parmater


dd1=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd1*r(k))
end do

do i=indtmin,indtm
  distp(i,3)=distp(i,1)*r(i)
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
    z(indorbp,i)=-distp(i,3+ic)*distp(i,3)
  end do
  !            endif
end do

if(typec.ne.1) then
  fun0=-distp(0,3)
  fun=-(1.d0-dd1*r(0))*distp(0,1)
  fun2=-dd1*(dd1*r(0)-2.d0)*distp(0,1)

  !              indorbp=indorb

  do ic=1,5
    !                 if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0)                       &
          *fun/r(0)
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
    z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0)+fun2)
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

case (33)
  !     4d without cusp and one parmater


  dd1=dd(indpar+1)
  !         if(iflagnorm.gt.2) then
  !     c=
  !                                                                  &1/(2.d0**3*3.d0)*dsqrt(1.d0/pi)*(2.d0*dd1)**(7.d0/2.d0)
      !     c=                                                       &
      !                                                              &1.d0/(2.d0**3*3.d0)/dsqrt(56.d0*pi)*(2.d0*dd1)**(9.d0/2.d0)
      c=dd1**4.5d0*0.0710812062076410d0
  !         endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
  end do

  do i=indtmin,indtm
    distp(i,3)=distp(i,1)*r(i)
    ! lz=0
    distp(i,4)=(3.d0*rmu(3,i)**2-r(i)**2)*cost1d
    ! lz=+/
    distp(i,5)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d
    ! lz=+/-2
    distp(i,6)=rmu(1,i)*rmu(2,i)*cost3d
    ! lz=+/-1
    distp(i,7)=rmu(2,i)*rmu(3,i)*cost3d
    ! lz=+/-1
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
    fun=-dd1*distp(0,3)+distp(0,1)
    fun2=dd1**2*distp(0,3)-2.d0*dd1*distp(0,1)

    !              indorbp=indorb

    do ic=1,5
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0)                     &
            *fun/r(0)
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
      z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0)+fun2)
      !endif for iocc
      !                endif
      ! enddo fot ic
    end do
  end if
  !
  indpar=indpar+1
  indshell=indshell+5
  indorb=indorbp


  ! 2s single  Z WITH CUSP zero
case (154)
! Jastrow single gaussian f orbital
! R(r)= exp(-alpha r^2)
! unnormalized


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
    z(indorbp,k)=distp(k,1)*distp(k,1+ic)
  end do
  !           endif
end do


if(typec.ne.1) then

  !            dd1=dd(indparp)
  fun0=distp(0,1)
  fun=-2.d0*dd1*distp(0,1)
  fun2=fun*(1.d0-2.d0*dd1*r(0)**2)


  !              indorbp=indorb


  do ic=1,7
    !                 if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)                       &
          *fun
    end do
    if(ic.eq.1) then
      z(indorbp,indt+1)=z(indorbp,indt+1)-                           &
          6.d0*cost1f*fun0*rmu(1,0)*rmu(3,0)
      z(indorbp,indt+2)=z(indorbp,indt+2)-                           &
          6.d0*cost1f*fun0*rmu(2,0)*rmu(3,0)
      z(indorbp,indt+3)=z(indorbp,indt+3)+                           &
          cost1f*fun0*(9.d0*rmu(3,0)**2-3.d0*r(0)**2)
    elseif(ic.eq.2) then
      z(indorbp,indt+1)=z(indorbp,indt+1)+                           &
          cost2f*fun0*(5.d0*rmu(3,0)**2-r(0)**2-2.d0*rmu(1,0)**2)
      z(indorbp,indt+2)=z(indorbp,indt+2)-                           &
          2.d0*cost2f*fun0*rmu(2,0)*rmu(1,0)
      z(indorbp,indt+3)=z(indorbp,indt+3)+                           &
          8.d0*cost2f*fun0*rmu(3,0)*rmu(1,0)
    elseif(ic.eq.3) then
      z(indorbp,indt+1)=z(indorbp,indt+1)-                           &
          2.d0*cost2f*fun0*rmu(1,0)*rmu(2,0)
      z(indorbp,indt+2)=z(indorbp,indt+2)+                           &
          cost2f*fun0*(5.d0*rmu(3,0)**2-r(0)**2-2.d0*rmu(2,0)**2)
      z(indorbp,indt+3)=z(indorbp,indt+3)+                           &
          8.d0*cost2f*fun0*rmu(3,0)*rmu(2,0)
    elseif(ic.eq.4) then
      z(indorbp,indt+1)=z(indorbp,indt+1)+                           &
          2.d0*cost3f*fun0*rmu(1,0)*rmu(3,0)
      z(indorbp,indt+2)=z(indorbp,indt+2)-                           &
          2.d0*cost3f*fun0*rmu(2,0)*rmu(3,0)
      z(indorbp,indt+3)=z(indorbp,indt+3)+                           &
          cost3f*fun0*(rmu(1,0)**2-rmu(2,0)**2)
    elseif(ic.eq.5) then
      z(indorbp,indt+1)=z(indorbp,indt+1)+                           &
          2.d0*cost3f*fun0*rmu(2,0)*rmu(3,0)
      z(indorbp,indt+2)=z(indorbp,indt+2)+                           &
          2.d0*cost3f*fun0*rmu(1,0)*rmu(3,0)
      z(indorbp,indt+3)=z(indorbp,indt+3)+                           &
          2.d0*cost3f*fun0*rmu(1,0)*rmu(2,0)
    elseif(ic.eq.6) then
      z(indorbp,indt+1)=z(indorbp,indt+1)+                           &
          3.d0*cost4f*fun0*(rmu(1,0)**2-rmu(2,0)**2)
      z(indorbp,indt+2)=z(indorbp,indt+2)-                           &
          6.d0*cost4f*fun0*rmu(1,0)*rmu(2,0)
    else
      z(indorbp,indt+1)=z(indorbp,indt+1)+                           &
          6.d0*cost4f*fun0*rmu(1,0)*rmu(2,0)
      z(indorbp,indt+2)=z(indorbp,indt+2)+                           &
          3.d0*cost4f*fun0*(rmu(1,0)**2-rmu(2,0)**2)
    end if
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



case (34)
  ! normalized
  ! exp(-dd1*r) + dd1*r*exp(-dd1*r)


  indshellp=indshell+1


  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  !           peff=dd1


  !           if(iflagnorm.gt.2) then
  !              c=1.d0/dsqrt(1.d0/4.d0/dd1**3+12.d0*dd1/(2.d0*dd1)**4+&
      !                                                              &3.d0*dd1**2/4.d0/dd1**5)/dsqrt(4.d0*pi)
      !             c=dd1*dsqrt(dd1)/dsqrt(7.d0*pi)
  c=dd1*dsqrt(dd1)*.2132436186229231d0
  !           endif

  do i=indtmin,indtm
    distp(i,1)=c*dexp(-dd1*r(i))
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)*(1.d0+r(i)*dd1)
  end do

  if(typec.ne.1) then

    fun=-dd1**2*distp(0,1)
    fun2=fun*(1.d0-dd1*r(0))

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=(2.d0*fun+fun2)

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+1
  indshell=indshellp




  ! 2s single  Z WITH CUSP
case (18)
  ! R(r)=r**4*exp(-z*r**2) single zeta



  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)

  !           if(iflagnorm.gt.2) then
  !           c=(2.d0*dd1**11/pi)**(1.d0/4.d0)*(512.d0/945.d0/pi)
  c=dd1**2.75d0*0.1540487967684377d0
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  do i=i0,indtm
    z(indorbp,i)=r(i)**4*distp(i,1)
  end do

  if(typec.ne.1) then
    rp1=r(0)**2

    !              the first derivative
    fun=distp(0,1)*rp1*(4.d0-2.d0*dd1*rp1)

    !              the second derivative
    fun2=distp(0,1)*rp1*(12.d0-18.d0*dd1*rp1                         &
        +4.d0*dd1**2*rp1**2)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*fun+fun2

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+1
  indshell=indshellp

  ! derivative of 16 with respect to z
case (41)
  !c     4p without cusp condition derivative of 22
  !c      r^2   e^{-z1 r }


  dd1=dd(indpar+1)
  !        if(iflagnorm.gt.2) then
  !        c=dsqrt((2.d0*dd1)**7/240.d0/pi)/2.d0
  c=dd1**3.5d0*0.2060129077457011d0
  !        endif
  c0=-c

  c1=3.5d0*c/dd1

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
  end do

  do i=indtmin,indtm
    distp(i,3)=r(i)**2*distp(i,1)
  end do

  !        indorbp=indorb

  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*(c0*distp(i,3)+c1*r(i)*distp(i,1))
    end do
    !           endif
  end do


  if(typec.ne.1) then
    !           fun=(1.d0-dd1*r(0))*distp(0,1)
    !           fun2=dd1*(dd1*r(0)-2.d0)*distp(0,1)

    fun0=c0*distp(0,3)+c1*r(0)*distp(0,1)
    fun=(c0*(2.d0-dd1*r(0))*r(0)                                     &
        +c1*(1.d0-dd1*r(0)))*distp(0,1)
    fun2=(c0*((dd1*r(0))**2+2.d0-4.d0*dd1*r(0))                      &
        +c1*dd1*(dd1*r(0)-2.d0))*distp(0,1)

    !              indorbp=indorb
    do ic=1,3
      !                if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun/r(0)
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun/r(0)+fun2)
      !                 endif
    end do
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp

case (125)
!     2s with cusp condition
!     dd1*( dd3 +exp(-dd2*r))  ! with no cusp condition


dd2=dd(indpar+1)
dd3=dd(indpar+2)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)+dd3
end do
!           endif


if(typec.ne.1) then
  fun=-dd2*distp(0,1)/r(0)
  fun2=dd2**2*distp(0,1)

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*fun+fun2





  !endif for indt
end if

indpar=indpar+2
indshell=indshellp
indorb=indorbp

case (116)
!     2p  double  Lorentian
!       dd1 * x_mu  (L^3(dd2 r)+dd3 r * L(dd4*r)^4) ; L(x)=1/(1+x)



dd2=dd(indpar+1)
dd3=dd(indpar+2)
dd4=dd(indpar+3)

do k=indtmin,indtm
  distp(k,1)=1.d0/(1.d0+dd2*r(k))**3
  distp(k,2)=r(k)/(1.d0+dd4*r(k))**4
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*(distp(i,1)+dd3*distp(i,2))
  end do
  !           endif
end do


if(typec.ne.1) then


  fun0=distp(0,1)+dd3*distp(0,2)
  fun=-3.d0*dd2*distp(0,1)/(r(0)*(1.d0+dd2*r(0)))                    &
      +dd3*distp(0,2)/r(0)**2*(1.d0-3*dd4*r(0))                      &
      /(1.d0+dd4*r(0))
  fun2=12.d0*dd2**2/(1.+dd2*r(0))**5                                 &
      +dd3*4.d0*dd4*(-2.d0+3.d0*dd4*r(0))/(1.+dd4*r(0))**6

  !           fun0=distp(0,1)+dd3*distp(0,2)
  !           fun=2.d0*(-dd2*distp(0,1)**2-dd4*dd3*distp(0,2)**2)

  !       fun2=2*dd2*distp(0,1)**3*(-1.d0+3.d0*dd2*r(0)**2)
  !    1+2*dd3*dd4*distp(0,2)**3*(-1.d0+3.d0*dd4*r(0)**2)

  !              indorbp=indorb

  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do



  !endif for indt
end if

indpar=indpar+3
indshell=indshell+3
indorb=indorbp

case (48)
  ! f single gaussian orbital
  ! R(r)= exp(-alpha r^2)
  ! normalized


  !        indorbp=indorb
  indparp=indpar+1

  dd1=dd(indparp)

  !        if(iflagnorm.gt.2) then
  ! overall normalization
  !        c=8.d0/dsqrt(15.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(9.d0/4.d0)
  c=dd1**2.25d0*1.47215808929909374563d0
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
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
    do ic=1,7
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)                     &
            *fun
      end do
      if(ic.eq.1) then
        z(indorbp,indt+1)=z(indorbp,indt+1)-                         &
            6.d0*cost1f*fun0*rmu(1,0)*rmu(3,0)
        z(indorbp,indt+2)=z(indorbp,indt+2)-                         &
            6.d0*cost1f*fun0*rmu(2,0)*rmu(3,0)
        z(indorbp,indt+3)=z(indorbp,indt+3)+                         &
            cost1f*fun0*(9.d0*rmu(3,0)**2-3.d0*r(0)**2)
      elseif(ic.eq.2) then
        z(indorbp,indt+1)=z(indorbp,indt+1)+                         &
            cost2f*fun0*(5.d0*rmu(3,0)**2-r(0)**2-2.d0*rmu(1,0)**2)
        z(indorbp,indt+2)=z(indorbp,indt+2)-                         &
            2.d0*cost2f*fun0*rmu(2,0)*rmu(1,0)
        z(indorbp,indt+3)=z(indorbp,indt+3)+                         &
            8.d0*cost2f*fun0*rmu(3,0)*rmu(1,0)
      elseif(ic.eq.3) then
        z(indorbp,indt+1)=z(indorbp,indt+1)-                         &
            2.d0*cost2f*fun0*rmu(1,0)*rmu(2,0)
        z(indorbp,indt+2)=z(indorbp,indt+2)+                         &
            cost2f*fun0*(5.d0*rmu(3,0)**2-r(0)**2-2.d0*rmu(2,0)**2)
        z(indorbp,indt+3)=z(indorbp,indt+3)+                         &
            8.d0*cost2f*fun0*rmu(3,0)*rmu(2,0)
      elseif(ic.eq.4) then
        z(indorbp,indt+1)=z(indorbp,indt+1)+                         &
            2.d0*cost3f*fun0*rmu(1,0)*rmu(3,0)
        z(indorbp,indt+2)=z(indorbp,indt+2)-                         &
            2.d0*cost3f*fun0*rmu(2,0)*rmu(3,0)
        z(indorbp,indt+3)=z(indorbp,indt+3)+                         &
            cost3f*fun0*(rmu(1,0)**2-rmu(2,0)**2)
      elseif(ic.eq.5) then
        z(indorbp,indt+1)=z(indorbp,indt+1)+                         &
            2.d0*cost3f*fun0*rmu(2,0)*rmu(3,0)
        z(indorbp,indt+2)=z(indorbp,indt+2)+                         &
            2.d0*cost3f*fun0*rmu(1,0)*rmu(3,0)
        z(indorbp,indt+3)=z(indorbp,indt+3)+                         &
            2.d0*cost3f*fun0*rmu(1,0)*rmu(2,0)
      elseif(ic.eq.6) then
        z(indorbp,indt+1)=z(indorbp,indt+1)+                         &
            3.d0*cost4f*fun0*(rmu(1,0)**2-rmu(2,0)**2)
        z(indorbp,indt+2)=z(indorbp,indt+2)-                         &
            6.d0*cost4f*fun0*rmu(1,0)*rmu(2,0)
      else
        z(indorbp,indt+1)=z(indorbp,indt+1)+                         &
            6.d0*cost4f*fun0*rmu(1,0)*rmu(2,0)
        z(indorbp,indt+2)=z(indorbp,indt+2)+                         &
            3.d0*cost4f*fun0*(rmu(1,0)**2-rmu(2,0)**2)
      end if
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


  ! derivative of 48 with respect to z
case (102)
  !     2s double gaussian with constant
  !     (dd3+ exp (-dd2 r^2)+dd4*exp(-dd5*r^2))



  dd2=dd(indpar+1)
  dd3=dd(indpar+2)
  dd4=dd(indpar+3)
  dd5=dd(indpar+4)

  indorbp=indorb+1
  indshellp=indshell+1
  do k=indtmin,indtm
    distp(k,1)=dexp(-dd2*r(k)*r(k))
    distp(k,2)=dexp(-dd5*r(k)*r(k))
  end do

  !           if(iocc(indshellp).eq.1) then
  do i=i0,indtm
    z(indorbp,i)=(distp(i,1)+dd3+dd4*distp(i,2))
!              write(6,*) ' function inside = ',z(indorbp,i)
end do
!           endif


if(typec.ne.1) then
  fun=-2.d0*(dd2*distp(0,1)+dd5*dd4*distp(0,2))
  fun2=r(0)**2

  !              write(6,*) ' fun inside = ',fun,fun2

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*(dd2*(-3.d0+2.d0*dd2*fun2)*                 &
      distp(0,1)+dd5*dd4*(-3.d0+2.d0*dd5*fun2)*distp(0,2))

  !           write(6,*) ' lap 106 =',z(indorbp,indt+4)


  !           stop

  !endif for indt
end if

indpar=indpar+4
indshell=indshellp
indorb=indorbp



case (35)
  ! normalized
  ! exp(-dd1*r) + dd1* r * exp(-dd2*r)


  indshellp=indshell+1


  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd1

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
    distp(k,2)=dexp(-dd2*r(k))
  end do

  !           if(iflagnorm.gt.2) then
  c=1.d0/dsqrt(1/4.d0/dd1**3+12*peff/(dd1+dd2)**4+                   &
      3*peff**2/4/dd2**5)/dsqrt(4.0*pi)
  !           endif

  do i=i0,indtm
    z(indorbp,i)=c*(distp(i,1)+r(i)*distp(i,2)*peff)
  end do

  if(typec.ne.1) then

    fun=-dd1*distp(0,1)+peff*distp(0,2)*(1.d0-dd2*r(0))
    fun2=distp(0,1)*dd1**2                                           &
        +peff*distp(0,2)*(dd2**2*r(0)-2.d0*dd2)

    do i=1,3
      z(indorbp,indt+i)=fun*c*rmu(i,0)/r(0)
    end do
    z(indorbp,indt+4)=c*(2.d0*fun/r(0)+fun2)

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+2
  indshell=indshellp

  ! single gaussian p orbitals
case (103)
!     2p single gaussian



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k)**2)
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)
  end do
  !           endif
end do


if(typec.ne.1) then
  fun0=distp(0,1)
  fun=-dd2*distp(0,1)*2.d0
  fun2=2.d0*dd2*(-1.d0+2.d0*dd2*r(0)**2)*                            &
      distp(0,1)

  !              indorbp=indorb
  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
    end do
    z(indorbp,indt+ic)=z(indorbp,indt+ic)+fun0
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do

  !endif for indt
end if

indpar=indpar+1
indshell=indshell+3
indorb=indorbp

case (120)
!     2p  double  cubic
!       dd1 * x_mu  (L^3(dd2 r)+dd3 L(dd4*r)^3) ; L(x)=1/(1+x)



dd2=dd(indpar+1)
dd3=dd(indpar+2)
dd4=dd(indpar+3)

do k=indtmin,indtm
  distp(k,1)=1.d0/(1.d0+dd2*r(k))**3
  distp(k,2)=1.d0/(1.d0+dd4*r(k))**3
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*(distp(i,1)+dd3*distp(i,2))
  end do
  !           endif
end do


if(typec.ne.1) then


  fun0=distp(0,1)+dd3*distp(0,2)
  fun=-3.d0*dd2*distp(0,1)/(r(0)*(1.d0+dd2*r(0)))                    &
      -3.d0*dd4*dd3*distp(0,2)/(r(0)*(1.d0+dd4*r(0)))
  fun2=12.d0*dd2**2/(1.+dd2*r(0))**5                                 &
      +12.d0*dd3*dd4**2/(1.+dd4*r(0))**5

  !           fun0=distp(0,1)+dd3*distp(0,2)
  !           fun=2.d0*(-dd2*distp(0,1)**2-dd4*dd3*distp(0,2)**2)

  !       fun2=2*dd2*distp(0,1)**3*(-1.d0+3.d0*dd2*r(0)**2)
  !    1+2*dd3*dd4*distp(0,2)**3*(-1.d0+3.d0*dd4*r(0)**2)

  !              indorbp=indorb

  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do

  !endif for indt
end if

indpar=indpar+3
indshell=indshell+3
indorb=indorbp

case (135)
!     2p single exponential  r^4 e^{-z r}  !



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)*r(i)**4
  end do
  !           endif
end do

if(typec.ne.1) then

  fun0=distp(0,1)*r(0)**4
  fun=distp(0,1)*(4.d0-dd2*r(0))*r(0)**2
  fun2=distp(0,1)*(12*r(0)**2-8*dd2*r(0)**3                          &
      +dd2**2*r(0)**4)
  !              indorbp=indorb
  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do
  !endif for indt
end if
indpar=indpar+1
indshell=indshell+3
indorb=indorbp


case (114)
!     2s without cusp condition
!     dd1*( dd3 +r^2/(1+dd2*r)^3)


dd2=dd(indpar+1)
dd3=dd(indpar+2)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=r(k)**2/(1.d0+dd2*r(k))**3
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)+dd3
end do
!           endif


if(typec.ne.1) then
  fun= (2.d0-dd2*r(0))/(1+dd2*r(0))**4
  fun2=2.d0*(1.d0-4.d0*dd2*r(0)+(dd2*r(0))**2)                       &
      /(1+dd2*r(0))**5
  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=fun2+2.d0*fun


  !endif for indt
end if

indpar=indpar+2
indshell=indshellp
indorb=indorbp

case (63)
  ! R(r)=c x*exp(-z*r**2)*r (c1 - r^2)


  dd1=dd(indpar+1)
  !        if(iflagnorm.gt.2) then
  c=dd1**1.75d0*1.2749263037197753d0
  !        c=2.d0/pi**0.75d0*(2.d0*dd1)**1.75d0/dsqrt(5.d0)
  !        endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do
  c1=1.75d0/dd1

  !        indorbp=indorb
  !
  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,1)*                             &
          (c1-r(i)**2)*r(i)
    end do
    !           endif
  end do

  if(typec.ne.1) then



    rp1=dd1*r(0)**2
    cost=2.d0*rp1

    fun0=distp(0,1)*r(0)*(c1-r(0)**2)
    fun=distp(0,1)*(c1*(1.d0-cost)/r(0)+                             &
        (-3.d0+cost)*r(0))
    !           My bug !!!
    !           fun2=distp(0,1)*(c1*2.d0*dd1*fun0*(cost-3.d0)
    !                                                                &-2.d0*r(0)*(3.d0-7.d0*rp1+2.d0*rp1**2))
        fun2=-2.d0*distp(0,1)*r(0)*                                  &
        (3.d0-7.d0*rp1+2.d0*rp1**2+c1*dd1*(3.d0-cost))


    !              indorbp=indorb
    do ic=1,3
      !                if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
      !                 endif
    end do
    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp

case (148)
!  derivative of 147 with respect to dd1


dd1=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd1*r(k)**2)
end do

do i=indtmin,indtm
  distp(i,3)=-r(i)**2*distp(i,1)
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
  fun=2.d0*(dd1*r(0)**2-1.d0)*r(0)*distp(0,1)
  fun2=-2.d0*(2.d0*dd1**2*r(0)**4+1.d0                               &
      -5.d0*dd1*r(0)**2)*distp(0,1)

  !              indorbp=indorb

  do ic=1,5
    !                 if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0)                       &
          *fun/r(0)
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
    z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0)+fun2)
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

case (12)
  ! R(r)=r**3*exp(-z1*r)
  !
  indshellp=indshell+1



  !       if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  !           if(iflagnorm.gt.2) then
  !           c=dsqrt((2*dd1)**9/40320.d0/pi)/2.d0
  c=dd1**4.5d0*.03178848180059307346d0
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)*r(i)**3
  end do

  if(typec.ne.1) then
    rp1=r(0)**3
    rp2=r(0)**2
    !
    !c              the first derivative
    fun=distp(0,1)*(3.d0*rp2-dd1*rp1)
    !c
    !c              the second derivative
    fun2=distp(0,1)*(6.d0*r(0)-6.d0*dd1*rp2+dd1**2*rp1)
    !c
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)/r(0)
    end do

    z(indorbp,indt+4)=2.d0*fun/r(0)+fun2

  end if
  !
  indorb=indorbp
  !
  !        endif
  indpar=indpar+1
  indshell=indshellp
  !


  ! 4s double zeta
case (1000:1099)
!     s gaussian  r**(2*npower)*exp(-alpha*r**2)

npower=iopt-1000

indorbp=indorb+1
indshellp=indshell+1

dd2=dd(indpar+1)
do k=indtmin,indtm
  distp(k,1)=r(k)**(2*npower)*dexp(-dd2*r(k)**2)
end do

!        if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)
end do
!        endif


if(typec.ne.1) then


  rp1=r(0)**2
  fun=(npower-dd2*rp1)*distp(0,1)*2.d0/rp1
  fun2=(npower*(2.d0*npower-1.d0)-                                   &
      (1.d0+4.d0*npower)*dd2*rp1+2.d0*(dd2*rp1)**2)*                 &
      distp(0,1)*2.d0/rp1

  !           if(iocc(indshellp).eq.1) then
  do i=1,3
    z(indorbp,indt+i)=rmu(i,0)*fun
  end do
  z(indorbp,indt+4)=2.d0*fun+fun2
  !           endif


  !endif for indt
end if

indpar=indpar+1
indshell=indshell+1
indorb=indorbp

case (144)
!     2p single exponential  -r^3 e^{-z r}  ! derivative of  130



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=-dexp(-dd2*r(k))
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)*r(i)**3
  end do
  !           endif
end do

if(typec.ne.1) then

  fun0=distp(0,1)*r(0)**3
  fun=distp(0,1)*(3.d0-dd2*r(0))*r(0)
  !           fun= derivative of fun0 respect to r divided dy r
  fun2=distp(0,1)*(dd2**2*r(0)**3-6*dd2*r(0)**2                      &
      +6*r(0))
  !           fun2= second derivative of fun0 respect to r
  !              indorbp=indorb
  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do
  !endif for indt
end if
indpar=indpar+1
indshell=indshell+3
indorb=indorbp

case (70)
  ! f single Slater orbital
  ! R(r)= exp(-alpha r)
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
      z(indorbp,k)=distp(k,1)*distp(k,1+ic)
    end do
    !           endif
  end do


  if(typec.ne.1) then

    dd1=dd(indparp)
    fun0=distp(0,1)
    fun=-dd1*distp(0,1)/r(0)
    fun2=dd1**2*distp(0,1)

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


case (100)
  !     2s single gaussian
  !     exp(-dd2*r^2)


  dd2=dd(indpar+1)


  indorbp=indorb+1
  indshellp=indshell+1
  do k=indtmin,indtm
    distp(k,1)=dexp(-dd2*r(k)*r(k))
  end do

  !           if(iocc(indshellp).eq.1) then
  do i=i0,indtm
    z(indorbp,i)=distp(i,1)
  end do
  !           endif


  if(typec.ne.1) then
    fun=-dd2*distp(0,1)*2.d0
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*dd2*(-3.d0+2.d0*dd2*r(0)**2)*             &
        distp(0,1)
    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshellp
  indorb=indorbp



case (138)
!     2s with cusp condition
!     ( -dd2*r^2*exp(-dd2*r))  ! with no cusp condition der of 137


dd2=dd(indpar+1)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=-dd2*dexp(-dd2*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)*r(i)**2
end do
!           endif


if(typec.ne.1) then
  fun=(2.d0-dd2*r(0))*distp(0,1)
  fun2=(2.d0-4*dd2*r(0)+(dd2*r(0))**2)*distp(0,1)
  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*fun+fun2


  !endif for indt
end if

indpar=indpar+1
indshell=indshellp
indorb=indorbp

case (56)
  ! g single Slater orbital derivative of 55
  ! R(r)= (11.d0/2.0 1/dd1 - r) * exp(-alpha r)
  ! normalized


  !        indorbp=indorb
  indparp=indpar+1

  dd1=dd(indparp)

  !        if(iflagnorm.gt.2) then
  ! overall normalization
  ! l = 4
  ! \int d\omega Y*Y = 4 pi / (2 l + 1)
  ! \int dr r^{2 l + 2} Exp [- 2 dd1 r^2 ] =  7 * 5**2 * 3**4 / 2**3 / dd1**11
  c=dd1**5.5d0*.020104801169736915d0
  !        c=1.d0/dsqrt(7.d0)*(2.d0/pi)**(1.d0/2.d0)*dd1**(11.d0/2.d0)/3.d0/5.d0
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
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
      z(indorbp,k)=distp(k,1)*distp(k,1+ic)*(11.d0/2.d0/dd1 - r(k))
    end do
    !           endif
  end do


  if(typec.ne.1) then

    dd1=dd(indparp)
    fun0=distp(0,1)*(11.d0/2.d0/dd1-r(0))
    fun=distp(0,1)*(dd1-13.d0/2.d0/r(0))
    fun2=dd1*distp(0,1)*(15.d0/2.d0-dd1*r(0))


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


case (1)



  indshellp=indshell+1


  !        if(iocc(indshellp).eq.1) then

  dd1=dd(indpar+1)
  !           if(iflagnorm.gt.2) then
  !           c=dd1*dsqrt(dd1)/dsqrt(pi)
  c=dd1*dsqrt(dd1)*0.56418958354775628695d0
  !           endif

  indorbp=indorb+1
  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)
  end do

  if(typec.ne.1) then
    fun=-dd1*distp(0,1)


    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)/r(0)
    end do

    z(indorbp,indt+4)=(-2.d0*dd1/r(0)+dd1**2)                        &
        *distp(0,1)


  end if

  indorb=indorbp

  !        endif

  indpar=indpar+1
  indshell=indshellp


  ! 1s double Z with cusp cond
case (49)
  ! f orbitals
  ! R(r)= c*exp(-z r^2)*(9/4/z-r^2)



  !        indorbp=indorb
  indparp=indpar+1
  dd1=dd(indparp)

  !        if(iflagnorm.gt.2) then
  ! overall normalization
  !        c=8.d0/dsqrt(15.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(9.d0/4.d0)
  c=dd1**2.25d0*1.47215808929909374563d0
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
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
      z(indorbp,k)=distp(k,1)*(9.d0/4.d0/dd1-r(k)**2)*               &
          distp(k,1+ic)
    end do
    !           endif
  end do


  if(typec.ne.1) then

    dd1=dd(indparp)
    fun0=distp(0,1)*(9.d0/4.d0/dd1-r(0)**2)
    fun=distp(0,1)*(2.d0*dd1*r(0)**2-13.d0/2.d0)
    fun2=distp(0,1)*(-4.d0*dd1**2*r(0)**4                            &
        +19.d0*dd1*r(0)**2-13.d0/2.d0)


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



case (141)
!     2p single exponential  r^2  e^{-z r}  ! parent of 121



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=-dexp(-dd2*r(k))
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)*r(i)**2
  end do
  !           endif
end do


if(typec.ne.1) then

  fun0=distp(0,1)*r(0)**2
  fun=distp(0,1)*(2.d0-dd2*r(0))
  fun2=(2.d0-4.d0*dd2*r(0)+(dd2*r(0))**2)*distp(0,1)
  !              indorbp=indorb
  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do
  !endif for indt
end if
indpar=indpar+1
indshell=indshell+3
indorb=indorbp

! der of 127
case (26)
  !     2p without cusp condition



  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd(indpar+3)
  dd3=dd(indpar+4)
  peff2=dd(indpar+5)

  !        if(iflagnorm.gt.2) then
  c=1.d0/2.d0/dsqrt(8.d0*pi*(1.d0/(2.d0*dd1)**5                      &
      +2.d0*peff/(dd1+dd2)**5+peff**2/(2.d0*dd2)**5                  &
      +2.d0*peff2/(dd1+dd3)**5+peff2**2/(2.d0*dd3)**5                &
      +2.d0*peff2*peff/(dd2+dd3)**5))
  !        endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=c*dexp(-dd2*r(k))
    distp(k,3)=c*dexp(-dd3*r(k))
  end do

  do i=indtmin,indtm
    distp(i,4)=distp(i,1)+peff*distp(i,2)+peff2*distp(i,3)
  end do

  !        indorbp=indorb

  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,4)
    end do
    !           endif
  end do


  if(typec.ne.1) then
    fun=(-dd1*distp(0,1)-dd2*peff*distp(0,2)                         &
        -dd3*peff2*distp(0,3))/r(0)
    fun2=dd1**2*distp(0,1)+peff*dd2**2*distp(0,2)                    &
        +peff2*dd3**2*distp(0,3)

    !              indorbp=indorb

    do ic=1,3
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*fun
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+distp(0,4)
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
      !                 endif
    end do


    !endif for indt
  end if

  indpar=indpar+5
  indshell=indshell+3
  indorb=indorbp

  ! 3p triple zeta
case (86)
  ! f single gaussian orbital
  ! R(r)= exp(-alpha r^2)
  ! normalized


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
      z(indorbp,k)=distp(k,1)*distp(k,1+ic)
    end do
    !           endif
  end do

  if(typec.ne.1) then

    fun0=distp(0,1)
    rp1=dd1*r(0)**2
    rp2=dd2*r(0)
    rp3=(1.d0+rp2)**2

    fun=-dd1*distp(0,1)*(2.d0+rp2)/rp3


    !              the second derivative
    fun2=dd1*distp(0,1)*(-2.d0-2.d0*rp2+4.d0*rp1+4.d0*rp1*rp2+rp1**2)/rp3**2



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


  ! derivative of 48 with respect to z
case (101)
  !     2s without cusp condition
  !     dd1*( dd3 +exp(-dd2*r^2))


  dd2=dd(indpar+1)
  dd3=dd(indpar+2)

  indorbp=indorb+1
  indshellp=indshell+1
  do k=indtmin,indtm
    distp(k,1)=dexp(-dd2*r(k)*r(k))
  end do

  !           if(iocc(indshellp).eq.1) then
  do i=i0,indtm
    z(indorbp,i)=distp(i,1)+dd3
  end do
  !           endif


  if(typec.ne.1) then
    fun=-dd2*distp(0,1)*2.d0

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*dd2*(-3.d0+2.d0*dd2*r(0)**2)*             &
        distp(0,1)


    !endif for indt
  end if

  indpar=indpar+2
  indshell=indshellp
  indorb=indorbp


case (150)
!     2p single exponential  r e^{-z r^2}



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k)**2)
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)*r(i)
  end do
  !           endif
end do


if(typec.ne.1) then

  fun0=distp(0,1)*r(0)
  cost=2.d0*dd2*r(0)**2
  fun=distp(0,1)*(1.d0-cost)/r(0)
  fun2=2.d0*dd2*fun0*(cost-3.d0)
  !              indorbp=indorb
  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do
  !endif for indt
end if
indpar=indpar+1
indshell=indshell+3
indorb=indorbp

case (155)
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



case (83)
  ! R(r)=x*exp(-z*r**2)*(5/4/z-r**2)

  dd1=dd(indpar+1)
  dd2=dsqrt(dd1)
  c=dd1**1.25d0*ratiocp

  do k=indtmin,indtm
    cost=dd1*r(k)**2/(1.d0+dd2*r(k))
    distp(k,1)=c*dexp(-cost)
  end do
  !        indorbp=indorb
  !
  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      cost=(1.d0+0.5d0*dd2*r(i))/(1.d0+dd2*r(i))**2
      z(indorbp,i)=rmu(ic,i)*distp(i,1)*(1.25d0/dd1-r(i)**2*cost)
    end do
    !           endif
  end do

  if(typec.ne.1) then
    rp1=dd1*r(0)**2
    rp2=dd2*r(0)
    rp3=(1.d0+rp2)**2
    cost=(1.d0+0.5d0*rp2)/rp3
    fun0=distp(0,1)*(1.25d0/dd1-r(0)**2*cost)
    fun=0.25d0*distp(0,1)*                                           &
        (-18.d0-39.d0*rp2-20.d0*rp1+rp1*rp2+2.d0*rp1**2)/rp3**2
    fun2=0.25d0*distp(0,1)*                                          &
        (-18.d0-42.d0*rp2+30.d0*rp1+138.d0*rp1*rp2+113.d0*rp1**2     &
        +30.d0*rp1**2*rp2-3.d0*rp1**3-2.d0*rp1**3*rp2)/rp3**3

    !              indorbp=indorb

    do ic=1,3
      !                if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
      !                 endif
    end do
    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp



case (81)
  ! R(r)=c*exp(-z*r**2)*(3/4/z-r**2)

  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1

  dd1=dd(indpar+1)
  dd2=dsqrt(dd1)

  !         c=(2.d0*dd1/pi)**(3.d0/4.d0)*ratiocs
  c=dd1**0.75d0*ratiocs

  do k=indtmin,indtm
    cost=dd1*r(k)**2/(1.d0+dd2*r(k))
    distp(k,1)=c*dexp(-cost)
  end do

  do i=i0,indtm
    cost=(1.d0+0.5d0*dd2*r(i))/(1.d0+dd2*r(i))**2
    z(indorbp,i)=distp(i,1)*(3.d0/4.d0/dd1-r(i)**2*cost)
  end do

  if(typec.ne.1) then
    !              the first derivative /r

    rp1=dd1*r(0)**2
    rp2=dd2*r(0)
    rp3=(1.d0+rp2)**2

    fun=0.25d0*distp(0,1)*                                           &
        (-14.d0-29.d0*rp2-12.d0*rp1+3.d0*rp1*rp2+2.d0*rp1**2)/rp3**2
    !              the second derivative
    fun2=0.25d0*distp(0,1)*                                          &
        (-14.d0-30.d0*rp2+34.d0*rp1+118.d0*rp1*rp2+87.d0*rp1**2      &
        +18.d0*rp1**2*rp2-5.d0*rp1**3-2.d0*rp1**3*rp2)/rp3**3

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  indpar=indpar+1
  indshell=indshellp

case (130)
!     2p single exponential  r^2  e^{-z r}  ! parent of 121



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)*r(i)**2
  end do
  !           endif
end do


if(typec.ne.1) then

  fun0=distp(0,1)*r(0)**2
  fun=distp(0,1)*(2.d0-dd2*r(0))
  fun2=(2.d0-4.d0*dd2*r(0)+(dd2*r(0))**2)*distp(0,1)
  !              indorbp=indorb
  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do
  !endif for indt
end if
indpar=indpar+1
indshell=indshell+3
indorb=indorbp

case (89)
  ! g single gaussian orbital
  ! derivative of 51
  ! R(r)= exp(-alpha r^2)
  ! normalized


  !        indorbp=indorb
  indparp=indpar+1

  dd1=dd(indparp)
  dd2=dsqrt(dd1)

  !     if(iflagnorm.gt.2) then
  ! overall normalization
  !     c=16.d0/dsqrt(105.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(11.d0/4.d0)*ratiocg
  c=dd1**2.75d0*ratiocg
  !     endif

  do k=indtmin,indtm
    cost=dd1*r(k)**2/(1.d0+dd2*r(k))
    distp(k,1)=c*dexp(-cost)
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
      cost=(1.d0+0.5d0*dd2*r(k))/(1.d0+dd2*r(k))**2
      z(indorbp,k)=distp(k,1)*(11.d0/4.d0/dd1-r(k)**2*cost)*         &
          distp(k,1+ic)
    end do
    !           endif
  end do


  if(typec.ne.1) then

    rp1=dd1*r(0)**2
    rp2=dd2*r(0)
    rp3=(1.d0+rp2)**2

    cost=(1.d0+0.5d0*rp2)/rp3
    fun0=distp(0,1)*(11.d0/4.d0/dd1-r(0)**2*cost)

    fun=0.25d0*distp(0,1)*                                           &
        (-30.d0-69.d0*rp2-44.d0*rp1-5.d0*rp1*rp2+2.d0*rp1**2)/rp3**2
    fun2=0.25d0*distp(0,1)*                                          &
        (-30.d0-78.d0*rp2+18.d0*rp1+198.d0*rp1*rp2+191.d0*rp1**2     &
        +66.d0*rp1**2*rp2+3.d0*rp1**3-2.d0*rp1**3*rp2)/rp3**3


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

  ! WARNING IN DFT it is assumed that UMRIGAR orbitals could be extended
  ! up to number 99, so i,h,... are possible extensions.



  ! 1s single Z NO CUSP!
case (1100:1199)
!     p gaussian  r**(2*npower)*exp(-alpha*r**2)

npower=iopt-1100

!        indorbp=indorb

dd2=dd(indpar+1)
do k=indtmin,indtm
  distp(k,1)=r(k)**(2*npower)*dexp(-dd2*r(k)**2)
end do

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)
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

  !           indorbp=indorb
  do ic=1,3
    !              if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
    !              endif
  end do



  !endif for indt
end if

indpar=indpar+1
indshell=indshell+3
indorb=indorbp


case (119)
!     2p single   r_mu/(1+b r^2)^(3/2)   parent of 103



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=1.d0/(1.d0+dd2*r(k)**2)**1.5d0
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)
  end do
  !           endif
end do



if(typec.ne.1) then

  fun0=distp(0,1)
  fun=-3.d0*dd2*distp(0,1)/(1.d0+dd2*r(0)**2)
  fun2=3.d0*dd2*(-1.d0+4.d0*dd2*r(0)**2)                             &
      /(1.d0+dd2*r(0)**2)**3.5d0

  !              indorbp=indorb

  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do

  !endif for indt
end if

indpar=indpar+1
indshell=indshell+3
indorb=indorbp

case (27)
  !     2p without cusp condition



  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd(indpar+3)
  dd3=dd(indpar+4)
  peff2=dd(indpar+5)

  !        if(iflagnorm.gt.2) then
  c=1.d0/2.d0/dsqrt(240.d0*pi*(1.d0/(2.d0*dd1)**7                    &
      +2.d0*peff/(dd1+dd2)**7+peff**2/(2.d0*dd2)**7                  &
      +2.d0*peff2/(dd1+dd3)**7+peff2**2/(2.d0*dd3)**7                &
      +2.d0*peff2*peff/(dd2+dd3)**7))
  !        endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=c*dexp(-dd2*r(k))
    distp(k,3)=c*dexp(-dd3*r(k))
  end do

  do i=indtmin,indtm
    distp(i,4)=r(i)*(distp(i,1)+peff*distp(i,2)                      &
        +peff2*distp(i,3))
  end do

  !        indorbp=indorb

  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,4)
    end do
    !           endif
  end do


  if(typec.ne.1) then
    fun0=distp(0,4)
    fun=(1.d0-dd1*r(0))*distp(0,1)                                   &
        +peff*(1.d0-dd2*r(0))*distp(0,2)                             &
        +peff2*(1.d0-dd3*r(0))*distp(0,3)
    fun2=dd1*(dd1*r(0)-2.d0)*distp(0,1)                              &
        +peff*dd2*(dd2*r(0)-2.d0)*distp(0,2)                         &
        +peff2*dd3*(dd3*r(0)-2.d0)*distp(0,3)

    !              indorbp=indorb

    do ic=1,3
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun/r(0)
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*                                   &
          (4.d0*fun/r(0)+fun2)
      !                 endif
    end do


    !endif for indt
  end if

  indpar=indpar+5
  indshell=indshell+3
  indorb=indorbp


case (85)
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

case (115)
!     2s double lorentian with constant  parent of 102
!     (dd3+ r^2/(1+dd2*r)^3+dd4*r^3/(1+dd5*r)^4;



dd2=dd(indpar+1)
dd3=dd(indpar+2)
dd4=dd(indpar+3)
dd5=dd(indpar+4)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=r(k)**2/(1.d0+dd2*r(k))**3
  distp(k,2)=r(k)**3/(1.d0+dd5*r(k))**4
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=(distp(i,1)+dd3+dd4*distp(i,2))
!              write(6,*) ' function inside = ',z(indorbp,i)
end do
!           endif


if(typec.ne.1) then

  fun= (2.d0-dd2*r(0))/(1+dd2*r(0))**4                               &
      -dd4*r(0)*(-3.d0+dd5*r(0))/(1.d0+dd5*r(0))**5
  fun2=2.d0*(1.d0-4.d0*dd2*r(0)+(dd2*r(0))**2)                       &
      /(1+dd2*r(0))**5                                               &
      +dd4*2.d0*r(0)*(3.d0-6.d0*dd5*r(0)+(dd5*r(0))**2)              &
      /(1.d0+dd5*r(0))**6


  !              write(6,*) ' fun inside = ',fun,fun2

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=fun2+2.d0*fun

  !           write(6,*) ' lap 106 =',z(indorbp,indt+4)

  !endif for indt
end if

indpar=indpar+4
indshell=indshellp
indorb=indorbp

case (22)
  !      3p without cusp condition
  !      r e^{-z1 r }



  dd1=dd(indpar+1)
  !        c=dsqrt((2.d0*dd1)**7/240.d0/pi)/2.d0
  c=dd1**3.5d0*0.2060129077457011d0
  !
  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=r(k)*distp(k,1)
  end do
  !
  !        indorbp=indorb
  !
  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,2)
    end do
    !           endif
  end do
  !
  !
  if(typec.ne.1) then
    fun0=distp(0,2)
    fun=(1.d0-dd1*r(0))*distp(0,1)
    fun2=dd1*(dd1*r(0)-2.d0)*distp(0,1)
    !
    !              indorbp=indorb
    !
    do ic=1,3
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun/r(0)
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*                                   &
          (4.d0*fun/r(0)+fun2)
      !
      !                endif
    end do
    !
    !
    !endif for indt
  end if
  !
  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp


  ! 3p double zeta
case (13)
  ! R(r)=r**3*(exp(-z1*r)+z3*exp(-z2*r))
  !
  indshellp=indshell+1

  !
  !
  !        if(iocc(indshellp).eq.1) then
  !
  indorbp=indorb+1
  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  dd3=dd(indpar+3)
  !           if(iflagnorm.gt.2) then
  c=1.d0/2.d0/dsqrt(pi*40320.d0*(1.d0/(2.d0*dd1)**9+                 &
      2.d0*dd3/(dd1+dd2)**9+dd3**2/(2.d0*dd2)**9))
  !           endif

  !
  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=c*dexp(-dd2*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=(distp(i,1)+dd3*distp(i,2))*r(i)**3
  end do
  !
  if(typec.ne.1) then
    rp1=r(0)**3
    rp2=r(0)**2
    !
    !c              the first derivative
    fun=distp(0,1)*(3.d0*rp2-dd1*rp1)                                &
        +dd3*distp(0,2)*(3.d0*rp2-dd2*rp1)
    !c
    !              the second derivative
    fun2=distp(0,1)*(6.d0*r(0)-6.d0*dd1*rp2+dd1**2*rp1)              &
        +dd3*distp(0,2)*(6.d0*r(0)-6.d0*dd2*rp2+dd2**2*rp1)
    !c
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)/r(0)
    end do
    !
    z(indorbp,indt+4)=2.d0*fun/r(0)+fun2
    !
  end if

  indorb=indorbp

  !        endif
  indpar=indpar+3
  indshell=indshellp

  ! 1s single Z pseudo
case (37,68)
  ! d orbitals
  ! R(r)= exp(-alpha r^2)
  ! each gaussian term is normalized


  !        indorbp=indorb
  indparp=indpar+1

  dd1=dd(indparp)

  !        if(iflagnorm.gt.2) then
  ! overall normalization
  !        c=4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)
  c=dd1**1.75d0*1.64592278064948967213d0
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
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
    do ic=1,5
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)                     &
            *fun
      end do
      if(ic.eq.1) then
        !                         if(i.ne.3) then
        z(indorbp,indt+1)=z(indorbp,indt+1)-                         &
            2.d0*rmu(1,0)*fun0*cost1d
        z(indorbp,indt+2)=z(indorbp,indt+2)-                         &
            2.d0*rmu(2,0)*fun0*cost1d
        !                         else
        z(indorbp,indt+3)=z(indorbp,indt+3)+                         &
            4.d0*rmu(3,0)*fun0*cost1d
        !                         endif
      elseif(ic.eq.2) then
        !                         if(i.eq.1) then
        z(indorbp,indt+1)=z(indorbp,indt+1)+                         &
            2.d0*rmu(1,0)*fun0*cost2d
        !                         elseif(i.eq.2) then
        z(indorbp,indt+2)=z(indorbp,indt+2)-                         &
            2.d0*rmu(2,0)*fun0*cost2d
        !                         endif
      elseif(ic.eq.3) then
        !                         if(i.eq.1) then
        z(indorbp,indt+1)=z(indorbp,indt+1)+                         &
            rmu(2,0)*fun0*cost3d
        !                         elseif(i.eq.2) then
        z(indorbp,indt+2)=z(indorbp,indt+2)+                         &
            rmu(1,0)*fun0*cost3d
        !                         endif
      elseif(ic.eq.4) then
        !                         if(i.eq.2) then
        z(indorbp,indt+2)=z(indorbp,indt+2)+                         &
            rmu(3,0)*fun0*cost3d
        !                         elseif(i.eq.3) then
        z(indorbp,indt+3)=z(indorbp,indt+3)+                         &
            rmu(2,0)*fun0*cost3d
        !                         endif
      elseif(ic.eq.5) then
        !                         if(i.eq.1) then
        z(indorbp,indt+1)=z(indorbp,indt+1)+                         &
            rmu(3,0)*fun0*cost3d
        !                         elseif(i.eq.3) then
        z(indorbp,indt+3)=z(indorbp,indt+3)+                         &
            rmu(1,0)*fun0*cost3d
        !endif for i
        !                         endif
        !endif for ic
      end if
      !enddo for i
      !                   enddo
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

case (140)
!     2p single exponential  -r e^{-z r}  ! der  of 121



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=-dexp(-dd2*r(k))
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)*r(i)
  end do
  !           endif
end do


if(typec.ne.1) then

  fun0=distp(0,1)*r(0)
  fun=distp(0,1)*(1.d0-dd2*r(0))/r(0)
  fun2=dd2*(dd2*r(0)-2.d0)*distp(0,1)
  !              indorbp=indorb
  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do
  !endif for indt
end if
indpar=indpar+1
indshell=indshell+3
indorb=indorbp

case (2000:2099)
!     s gaussian  -r**(2*(npower+1))*exp(-alpha*r**2) derivative of 1000

npower=iopt+1-2000

indorbp=indorb+1
indshellp=indshell+1


dd2=dd(indpar+1)
do k=indtmin,indtm
  distp(k,1)=-r(k)**(2*npower)*dexp(-dd2*r(k)**2)
end do

!        if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)
end do
!        endif


if(typec.ne.1) then


  rp1=r(0)**2
  fun0=distp(0,1)
  fun=(npower-dd2*rp1)*distp(0,1)*2.d0/rp1
  fun2=(npower*(2.d0*npower-1.d0)-                                   &
      (1.d0+4.d0*npower)*dd2*rp1+2.d0*(dd2*rp1)**2)*                 &
      distp(0,1)*2.d0/rp1

  !           if(iocc(indshellp).eq.1) then
  do i=1,3
    z(indorbp,indt+i)=rmu(i,0)*fun
  end do
  z(indorbp,indt+4)=2.d0*fun+fun2
  !           endif


  !endif for indt
end if

indpar=indpar+1
indshell=indshell+1
indorb=indorbp


case (61)
  ! R(r)=c (-r**5*exp(-z1*r**2)+c1 r**3 exp(-z1*r**2))


  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)

  !          if(iflagnorm.gt.2) then
  !     c=2.d0/pi**(3.d0/4.d0)*(2.d0*dd1)**(9.d0/4.d0)*dsqrt(2.d0/105.d0)
  c=dd1**2.25d0*.55642345640820284397d0
  !           endif

  c1=2.25d0/dd1

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)*r(k)
  end do

  do i=i0,indtm
    z(indorbp,i)=(-r(i)**4+c1*r(i)**2)*distp(i,1)
  end do


  if(typec.ne.1) then
    rp1=r(0)**2
    rp2=rp1*dd1

    fun=c1*distp(0,1)*(3.d0-2.d0*rp2)                                &
        +distp(0,1)*rp1*(-5.d0+2.d0*rp2)
    !              the second derivative
    fun2=c1*distp(0,1)*(6.d0-14.d0*rp2+4.d0*rp2**2)                  &
        +distp(0,1)*rp1*(-20.d0+22.d0*rp2-4.d0*rp2**2)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  !        endif
  indpar=indpar+1
  indshell=indshellp

  ! single gaussianx r  p orbitals
case (20)
  !     2p single Z with no  cusp condition


  dd1=dd(indpar+1)
  !        if(iflagnorm.gt.2) then
  !        c=dsqrt((2.d0*dd1)**5/8.d0/pi)/2.d0
  c=dd1**2.5d0*0.5641895835477562d0
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
  end do

  !        indorbp=indorb
  !
  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,1)
    end do
    !           endif
  end do

  if(typec.ne.1) then
    fun0=distp(0,1)
    fun=-dd1*distp(0,1)
    fun2=dd1**2*distp(0,1)

    !              indorbp=indorb

    do ic=1,3
      !                if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun/r(0)
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun/r(0)+fun2)
      !                 endif
    end do
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp


  ! 2p double zeta
case (38)
  ! R(r)=r**2*exp(-z1*r)



  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  !           if(iflagnorm.gt.2) then
  !              c=1.d0/dsqrt(1.d0/4.d0/dd1**3+12.d0*dd1/(2.d0*dd1)**4+&
      !                                                              &3.d0*dd1**2/4.d0/dd1**5)/dsqrt(4.d0*pi)
      !             c=dd1*dsqrt(dd1)/dsqrt(7.d0*pi)
  c=dd1*dsqrt(dd1)*0.21324361862292308211d0
  !           endif

  c0=-c*dd1

  c1=1.5d0*c/dd1



  do i=indtmin,indtm
    distp(i,1)=dexp(-dd1*r(i))
  end do

  do i=i0,indtm
    z(indorbp,i)=(c0*r(i)**2+c1*(1.d0+dd1*r(i)))                     &
        *distp(i,1)
  end do

  c1=c1*dd1**2

  if(typec.ne.1) then
    fun=(c0*(2.d0-dd1*r(0))-c1)*distp(0,1)
    fun2=(c0*(2.d0-4*dd1*r(0)+(dd1*r(0))**2)                         &
        +c1*(dd1*r(0)-1.d0))*distp(0,1)
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  !        endif
  indpar=indpar+1
  indshell=indshellp

  ! 4s single zeta derivative of 10
case (84)
  ! d orbitals
  ! R(r)= exp(-alpha r^2)
  ! each gaussian term is normalized



  !        indorbp=indorb
  indparp=indpar+1

  dd1=dd(indparp)
  dd2=dsqrt(dd1)

  !        if(iflagnorm.gt.2) then
  ! overall normalization
  !        c=4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)*ratiocd
  c=ratiocd*dd1**1.75d0
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
      z(indorbp,k)=distp(k,1)*distp(k,1+ic)
    end do
    !           endif
  end do


  if(typec.ne.1) then

    fun0=distp(0,1)
    rp1=dd1*r(0)**2
    rp2=dd2*r(0)
    rp3=(1.d0+rp2)**2

    fun=-dd1*distp(0,1)*(2.d0+rp2)/rp3


    !              the second derivative
    fun2=dd1*distp(0,1)*(-2.d0-2.d0*rp2+4.d0*rp1+4.d0*rp1*rp2+rp1**2)/rp3**2



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

  ! derivative of 37 with respect to z
case (24)
  !c     4p without cusp condition
  !c      r^2   e^{-z1 r }


  dd1=dd(indpar+1)
  !        if(iflagnorm.gt.2) then
  !        c=dsqrt((2.d0*dd1)**9/120960.d0/pi)/2.d0
  c=dd1**4.5d0*0.01835308852470193d0
  !        endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
  end do

  do i=indtmin,indtm
    distp(i,3)=r(i)**2*distp(i,1)
  end do

  !        indorbp=indorb

  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,3)
    end do
    !           endif
  end do


  if(typec.ne.1) then
    fun0=distp(0,3)
    fun=(2.d0*r(0)-dd1*r(0)**2)*distp(0,1)
    fun2=((dd1*r(0))**2+2.d0-4.d0*dd1*r(0))*distp(0,1)
    !              indorbp=indorb
    do ic=1,3
      !                if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun/r(0)
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun/r(0)+fun2)
      !                 endif
    end do

    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp


  ! 4p double zeta
case (5)
  ! normalized

  indshellp=indshell+1


  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
  end do

  !           if(iflagnorm.gt.2) then
  !           c=dd1**2.5d0/dsqrt(3.d0*pi)
  c=dd1**2.5d0*0.32573500793527994772d0
  !           endif

  do i=i0,indtm
    z(indorbp,i)=c*r(i)*distp(i,1)
  end do

  if(typec.ne.1) then

    fun=distp(0,1)*(1.d0-dd1*r(0))
    fun2=distp(0,1)*(dd1**2*r(0)-2.d0*dd1)

    do i=1,3
      z(indorbp,indt+i)=c*fun*rmu(i,0)/r(0)
    end do
    z(indorbp,indt+4)=c*2.d0*fun/r(0)+c*fun2

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+1
  indshell=indshellp


  ! 2s double  Z NO CUSP
case (88)
  ! g single gaussian orbital
  ! R(r)= exp(-alpha r^2)
  ! normalized

  !        indorbp=indorb
  indparp=indpar+1

  dd1=dd(indparp)
  dd2=dsqrt(dd1)

  !     if(iflagnorm.gt.2) then
  ! overall normalization
  !     c=16.d0/dsqrt(105.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(11.d0/4.d0)*ratiocg
  c=dd1**2.75d0*ratiocg
  !     endif

  do k=indtmin,indtm
    cost=dd1*r(k)**2/(1.d0+dd2*r(k))
    distp(k,1)=c*dexp(-cost)
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
      z(indorbp,k)=distp(k,1)*distp(k,1+ic)
    end do
    !           endif
  end do


  if(typec.ne.1) then


    fun0=distp(0,1)
    rp1=dd1*r(0)**2
    rp2=dd2*r(0)
    rp3=(1.d0+rp2)**2

    fun=-dd1*distp(0,1)*(2.d0+rp2)/rp3


    !              the second derivative
    fun2=dd1*distp(0,1)*(-2.d0-2.d0*rp2+4.d0*rp1+4.d0*rp1*rp2+rp1**2)/rp3**2




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


case (2100:2199)
!     p gaussian  -r**(2*(npower+1))*exp(-alpha*r**2) derivative 1100

npower=iopt+1-2100

!        indorbp=indorb

dd2=dd(indpar+1)
do k=indtmin,indtm
  distp(k,1)=-r(k)**(2*npower)*dexp(-dd2*r(k)**2)
end do

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)
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

  !           indorbp=indorb
  do ic=1,3
    !              if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
    !              endif
  end do



  !endif for indt
end if

indpar=indpar+1
indshell=indshell+3
indorb=indorbp


case (72)


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

case (152)
!     2s with cusp condition
!     ( r^3*exp(-dd2*r^2))  ! with no cusp condition

dd2=dd(indpar+1)
indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k)**2)*r(k)
end do
!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)*r(i)**2
end do
!           endif
if(typec.ne.1) then
  rp1=r(0)**2*dd2
  fun=(3.d0-2.d0*rp1)*distp(0,1)
  fun2=(6.d0-14.d0*rp1+4.d0*rp1**2)*distp(0,1)
  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do
  z(indorbp,indt+4)=2.d0*fun+fun2
  !endif for indt
end if
indpar=indpar+1
indshell=indshellp
indorb=indorbp
case (126)
!     2s double exp  with constant
!     (dd3+ exp (-dd2 r)+dd4*exp(-dd5*r))



dd2=dd(indpar+1)
dd3=dd(indpar+2)
dd4=dd(indpar+3)
dd5=dd(indpar+4)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))
  distp(k,2)=dexp(-dd5*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)+dd3+dd4*distp(i,2)
!              write(6,*) ' function inside = ',z(indorbp,i)
end do
!           endif


if(typec.ne.1) then


  fun=-(dd2*distp(0,1)+dd5*dd4*distp(0,2))/r(0)
  fun2=dd2**2*distp(0,1)+dd4*dd5**2*distp(0,2)



  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*fun+fun2



  !endif for indt
end if

indpar=indpar+4
indshell=indshellp
indorb=indorbp

case (153)
!     2s with cusp condition
!     (-r^5*exp(-dd2*r^2))  ! derivative of 152

dd2=dd(indpar+1)
indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k)**2)*r(k)**3
end do
!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=-distp(i,1)*r(i)**2
end do
!           endif
if(typec.ne.1) then
  rp1=dd2*r(0)**2
  fun=(-5.d0+2.d0*rp1)*distp(0,1)
  fun2=(-20.d0+22.d0*rp1-4.d0*rp1**2)*distp(0,1)
  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do
  z(indorbp,indt+4)=2.d0*fun+fun2
  !endif for indt
end if
indpar=indpar+1
indshell=indshellp
indorb=indorbp


case (121)
!     2p single exponential



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)
  end do
  !           endif
end do


if(typec.ne.1) then


  fun0=distp(0,1)
  fun=-dd2*distp(0,1)/r(0)
  fun2=dd2**2*distp(0,1)

  !              indorbp=indorb
  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do



  !endif for indt
end if

indpar=indpar+1
indshell=indshell+3
indorb=indorbp

case (149)
! derivative of 131 with respect z_1
! - r^4 exp(-z_1 r^2)



dd2=dd(indpar+1)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k)*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=-distp(i,1)*r(i)**4
end do
!           endif


if(typec.ne.1) then
  fun0=dd2*r(0)**2
  fun=-2.d0*r(0)**2*distp(0,1)*(2.d0-fun0)
  fun2=-2.d0*r(0)**2*distp(0,1)*(6.d0-9.d0*fun0+2.d0*fun0**2)

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*fun+fun2


  !endif for indt
end if

indpar=indpar+1
indshell=indshellp
indorb=indorbp

case (147)
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

case (134)
!     2p single exponential  r^3 e^{-z r}  !



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)*r(i)**3
  end do
  !           endif
end do

if(typec.ne.1) then

  fun0=distp(0,1)*r(0)**3
  fun=distp(0,1)*(3.d0-dd2*r(0))*r(0)
  !           fun= derivative of fun0 respect to r divided dy r
  fun2=distp(0,1)*(dd2**2*r(0)**3-6*dd2*r(0)**2                      &
      +6*r(0))
  !           fun2= second derivative of fun0 respect to r
  !              indorbp=indorb
  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do
  !endif for indt
end if
indpar=indpar+1
indshell=indshell+3
indorb=indorbp


case (146)
!     2p single exponential  -r^2  e^{-z r^2}  ! derivative  of 103



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k)*r(k))
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=-rmu(ic,i)*distp(i,1)*r(i)*r(i)
  end do
  !           endif
end do


if(typec.ne.1) then
  rp2=dd2*r(0)*r(0)
  fun0=-distp(0,1)*r(0)*r(0)
  fun=distp(0,1)*(-2.d0+2.d0*rp2)
  fun2=(-2.d0+10.d0*rp2-4.d0*rp2*rp2)*distp(0,1)
  !              indorbp=indorb
  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do
  !endif for indt
end if
indpar=indpar+1
indshell=indshell+3
indorb=indorbp

case (25)
  !     4p without cusp condition
  !      r^2  ( e^{-z2 r } + z1 e^{-z3 r } )



  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  dd3=dd(indpar+3)
  !        if(iflagnorm.gt.2) then
  c=1.d0/2.d0/dsqrt(120960.d0*pi*(1.d0/(2.d0*dd1)**9                 &
      +2.d0*dd3/(dd1+dd2)**9+dd3**2/(2.d0*dd2)**9))
  !        endif
  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=c*dexp(-dd2*r(k))
  end do

  do i=indtmin,indtm
    distp(i,3)=r(i)**2*(distp(i,1)+dd3*distp(i,2))
  end do

  !        indorbp=indorb

  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,3)
    end do
    !           endif
  end do


  if(typec.ne.1) then
    fun0=distp(0,3)
    fun=(2.d0*r(0)-dd1*r(0)**2)*distp(0,1)                           &
        +dd3*(2.d0*r(0)-dd2*r(0)**2)*distp(0,2)
    fun2=((dd1*r(0))**2+2.d0-4.d0*dd1*r(0))*distp(0,1)               &
        +dd3*((dd2*r(0))**2+2.d0-4.d0*dd2*r(0))*distp(0,2)
    !              indorbp=indorb
    do ic=1,3
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun/r(0)
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun/r(0)+fun2)
      !                 endif
    end do

    !endif for indt
  end if

  indpar=indpar+3
  indshell=indshell+3
  indorb=indorbp


  ! 2p triple zeta
case (32)
  !     3d without cusp condition triple Z


  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd(indpar+3)
  dd3=dd(indpar+4)
  peff2=dd(indpar+5)

  !        if(iflagnorm.gt.2) then
  c=1/2.d0*dsqrt(5.d0/pi)                                            &
      /dsqrt(1/(2.d0*dd1)**7+2*peff/(dd1+dd2)**7                     &
      +peff**2/(2.d0*dd2)**7+2*peff2/(dd1+dd3)**7                    &
      +peff2**2/(2.d0*dd3)**7+2*peff*peff2/(dd2+dd3)**7)/dsqrt(720.d0)
  !        endif

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
    distp(k,2)=dexp(-dd2*r(k))
    distp(k,3)=dexp(-dd3*r(k))
  end do

  do i=indtmin,indtm
    distp(i,4)=c*(distp(i,1)+peff*distp(i,2)+peff2*distp(i,3))
    !lz=0
    distp(i,5)=(3.d0*rmu(3,i)**2-r(i)**2)*cost1d
    !lz=+/-2
    distp(i,6)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d
    ! lz=+/- 2
    distp(i,7)=rmu(1,i)*rmu(2,i)*cost3d
    ! lz=+/-1
    distp(i,8)=rmu(2,i)*rmu(3,i)*cost3d
    ! lz=+/-1
    distp(i,9)=rmu(1,i)*rmu(3,i)*cost3d
  end do

  !        indorbp=indorb

  do ic=1,5
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=distp(i,4+ic)*distp(i,4)
    end do
    !           endif
  end do

  if(typec.ne.1) then
    fun0=distp(0,4)
    fun=c*(-dd1*distp(0,1)-peff*dd2*distp(0,2)                       &
        -peff2*dd3*distp(0,3))
    fun2=c*(dd1**2*distp(0,1)+peff*dd2**2*distp(0,2)                 &
        +peff2*dd3**2*distp(0,3))

    !              indorbp=indorb

    do ic=1,5
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,4+ic)*rmu(i,0)*fun/r(0)
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
      z(indorbp,indt+4)=distp(0,4+ic)*(6.d0*fun/r(0)+fun2)
      !endif for iocc
      !                endif
      ! enddo fot ic
    end do


    !endif for indt
  end if

  indpar=indpar+5
  indshell=indshell+5
  indorb=indorbp

case (145)
!     2s without cusp condition  !derivative 100
!     -(r^2*exp(-dd2*r^2))


dd2=dd(indpar+1)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k)*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=-distp(i,1)*r(i)**2
end do
!           endif


if(typec.ne.1) then
  fun0=dd2*r(0)**2
  fun=-2.d0*distp(0,1)*(1.d0-fun0)
  fun2=-2.d0*distp(0,1)*(1.d0-5.d0*fun0+2.d0*fun0**2)

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*fun+fun2


  !endif for indt
end if

indpar=indpar+1
indshell=indshellp
indorb=indorbp

case (21)
  !     2p without cusp condition



  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd(indpar+3)
  c=0.5d0/dsqrt(8.d0*pi*(1.d0/(2.d0*dd1)**5                          &
      +2.d0*peff/(dd1+dd2)**5+peff**2/(2.d0*dd2)**5))
  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=c*dexp(-dd2*r(k))
  end do

  do i=indtmin,indtm
    distp(i,3)=distp(i,1)+peff*distp(i,2)
  end do

  !        indorbp=indorb

  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,3)
    end do
    !           endif
  end do


  if(typec.ne.1) then
    fun=(-dd1*distp(0,1)-dd2*peff*distp(0,2))/r(0)
    fun2=dd1**2*distp(0,1)+peff*dd2**2*distp(0,2)

    !              indorbp=indorb

    do ic=1,3
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*fun
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+distp(0,3)
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
      !                 endif
    end do


    !endif for indt
  end if

  indpar=indpar+3
  indshell=indshell+3
  indorb=indorbp


  ! 3p single zeta
case (108)
!     2s double lorentian with constant  parent of 102
!     (dd3+ L(dd2 r^2)+dd4*L(dd5*r^2)) ;  L(x)=1/1+x^2



dd2=dd(indpar+1)
dd3=dd(indpar+2)
dd4=dd(indpar+3)
dd5=dd(indpar+4)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=1.d0/(1.d0+dd2*r(k)*r(k))
  distp(k,2)=1.d0/(1.d0+dd5*r(k)*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=(distp(i,1)+dd3+dd4*distp(i,2))
!              write(6,*) ' function inside = ',z(indorbp,i)
end do
!           endif


if(typec.ne.1) then
  fun=-2.d0*(dd2*distp(0,1)**2+dd5*dd4*distp(0,2)**2)
  fun2=2.d0*dd2*distp(0,1)**3*(-1.d0+3.d0*dd2*r(0)**2)               &
      +2.d0*dd5*dd4*distp(0,2)**3*(-1.d0+3.d0*dd5*r(0)**2)

  !              write(6,*) ' fun inside = ',fun,fun2

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=fun2+2.d0*fun

  !           write(6,*) ' lap 106 =',z(indorbp,indt+4)

  !endif for indt
end if

indpar=indpar+4
indshell=indshellp
indorb=indorbp


case (131)
!     2s without cusp condition
!     dd1*(r^2*exp(-dd2*r^2))


dd2=dd(indpar+1)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k)*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)*r(i)**2
end do
!           endif


if(typec.ne.1) then
  fun0=dd2*r(0)**2
  fun=2.d0*distp(0,1)*(1.d0-fun0)
  fun2=2.d0*distp(0,1)*(1.d0-5.d0*fun0+2.d0*fun0**2)

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*fun+fun2


  !endif for indt
end if

indpar=indpar+1
indshell=indshellp
indorb=indorbp

case (133)
!     4d  one parmater


dd1=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd1*r(k))
end do

do i=indtmin,indtm
  distp(i,3)=distp(i,1)*r(i)
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
  fun=(1.d0-dd1*r(0))*distp(0,1)
  fun2=dd1*(dd1*r(0)-2.d0)*distp(0,1)

  !              indorbp=indorb

  do ic=1,5
    !                 if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0)                       &
          *fun/r(0)
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
    z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0)+fun2)
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

case (66)
  !     derivative of 57 (orbital 1s STO regolarized for r->0)
  !  dR(r)/dz    with    R(r)= C(z) * P(z*r) * exp(-z*r)
  !     P(x) = (x+a)^n/(1+(x+a)^n)    with a so that P(0)==dP(0)/dx
  !     C(z) = const * z^(3/2)       normalization
  ! the following definitions are in module constants
  !     n         -> costSTO1s_n = 4
  !     a         -> costSTO1s_a = 1.2263393530877080588
  !     const(n)  -> costSTO1s_c = 0.58542132302621750732
  !


  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1

  dd1=dd(indpar+1)

  !           if(iflagnorm.gt.2) then
  !             if(dd1.gt.0.) then
  c=costSTO1s_c*dd1**1.5d0
  !             else
  !               c=1.d0
  !             endif
  !           endif

  do i=indtmin,indtm
    distp(i,1)=c*dexp(-dd1*r(i))
  end do

  do i=i0,indtm
    rp1=dd1*r(i)+costSTO1s_a
    rp4=rp1**costSTO1s_n
    z(indorbp,i)=distp(i,1)*rp4/(1.d0+rp4)*                          &
        (1.5d0/dd1 + r(i)*                                           &
        (-1.d0+(costSTO1s_n)/(rp1*(1.d0+rp4))))
  end do

  if(typec.ne.1) then
    rp1=dd1*r(0)+costSTO1s_a
    rp2=rp1**2
    rp4=rp1**costSTO1s_n
    rp6=rp4**2
    !              the first derivative /r
    fun=distp(0,1)*(dd1*rp4*(-2.d0*costSTO1s_a*(costSTO1s_n**2*      &
        (-1.d0+rp4)+costSTO1s_n*(1.d0+2.d0*rp1)*(1.d0+rp4)-rp2*      &
        (1.d0+rp4)**2) +rp1*(2*costSTO1s_n**2*(-1+rp4)+costSTO1s_n   &
        *(-3.d0+4.d0*rp1)*(1.d0+rp4)-rp1*(-5.d0+2.d0*rp1)*(1.d0+     &
        rp4)**2)))/(2.d0*rp2*(costSTO1s_a-rp1)*(1.d0+rp4)**3)

    !              fun=+distp(0,1)*rp4/(1.d0+rp4)*(1.5d0/dd1 + r(0)* &
        !                                                            &(-1.d0+(costSTO1s_n)/(rp1*(1.d0+rp4))))               &
        !                                                            &*(dd1*(2.d0 + 5.d0/(-dd1*r(0)) +                 &
        !                                                            &(4.d0*costSTO1s_n**2)/(rp2*(1.d0+rp4)**2) +           &
        !                                                            &(costSTO1s_n*((3.d0 - 2.d0*costSTO1s_n - 4.d0*rp1)*rp1&
        !                                                            &+ 2.d0*costSTO1s_a*(1.d0+costSTO1s_n+2.d0*rp1)))/    &
        !                                                            &(rp2*(dd1*r(0))*(1 + rp4))))/2.d0

    !              the second derivative derivative
    fun2=-distp(0,1)*(dd1*rp4*(rp1*(-(costSTO1s_n*(-3.d0-8.d0*rp1+   &
        6.d0*rp2)*(1.d0+rp4)**2)+rp2*(-7.d0+2.d0*rp1)*(1.d0+rp4)**3- &
        costSTO1s_n**2*(-1.d0+6.d0*rp1)*(-1.d0+rp6)-2*costSTO1s_n**3*&
        (1.d0+rp4*(-4.d0+rp4))) + 2.d0*costSTO1s_a*(-(rp1*rp2*(1.d0 +&
        rp4)**3) + 3.d0*costSTO1s_n**2*(1.d0+rp1)*(-1.d0+rp6)+       &
        costSTO1s_n*(1.d0+rp4)**2*(2.d0+3.d0*rp1*(1.d0+rp1)) +       &
        costSTO1s_n**3*(1.d0+rp4*(-4.d0+rp4)))))/                    &
        (2.d0*rp1*rp2*(1+rp4)**4)

    !              fun2=-distp(0,1)*rp4/(1.d0+rp4)*(1.5d0/dd1 + r(0)*&
        !                                                            &(-1.d0+(costSTO1s_n)/(rp1*(1.d0+rp4))))               &
        !                                                            &*(dd1*(rp1*(-(costSTO1s_n*(-3 - 8*rp1 &
        !                                                            &+ 6*rp2)*(1 + rp4)**2) + rp2*(-7 + 2*rp1)*(1 + rp4)**3 -     &
        !                                                            &costSTO1s_n**2*(-1 + 6*rp1)*(-1 + rp6) - 2*costSTO1s_n**3*  &
        !                                                            &(1 + rp4*(-4 + rp4))) - 2*costSTO1s_a*(-(rp1*rp2*(1 + rp4)**3)&
        !                                                            &+ 3*costSTO1s_n**2*(1 + rp1)*(-1 + rp6) + costSTO1s_n*(1 +   &
        !                                                            &rp4)**2 *(2 + 3*rp1*(1 +rp1)) + costSTO1s_n**3*(1 + rp4*(-4 +&
        !                                                            &rp4)))))/(2.*rp1*rp2*(1 + rp4)**3)

    !    gradient: dR(r)/dr_i=r_i*fun
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    !    laplacian = 2*fun+fun2
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  !        endif

  indpar=indpar+1
  indshell=indshellp








case (57)
  !     orbital 1s (no cusp) - STO regolarized for r->0
  ! R(r)= C(z) * P(z*r) * exp(-z*r)
  !     P(x) = (x+a)^n/(1+(x+a)^n)    with a so that P(0)==dP(0)/dx
  !     C(z) = const * z^(3/2)       normalization
  ! the following definitions are in module constants
  !     n         -> costSTO1s_n = 4
  !     a         -> costSTO1s_a = 1.2263393530877080588
  !     const(n)  -> costSTO1s_c = 0.58542132302621750732
  !
  !


  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1

  dd1=dd(indpar+1)

  !           if(iflagnorm.gt.2) then
  !             if(dd1.gt.0.) then
  c=costSTO1s_c*dd1**1.5d0
  !             else
  !               c=1.d0
  !             endif
  !           endif

  do i=indtmin,indtm
    distp(i,1)=c*dexp(-dd1*r(i))
  end do

  do i=i0,indtm
    rp4=(dd1*r(i)+costSTO1s_a)**costSTO1s_n
    z(indorbp,i)=distp(i,1)*rp4/(1.d0+rp4)
  end do

  if(typec.ne.1) then
    rp1=dd1*r(0)+costSTO1s_a
    rp2=rp1**2
    rp4=rp1**costSTO1s_n
    rp6=rp4**2
    !              the first derivative /r
    !fun=-z(indorbp,0)*((dd1**2*(-costSTO1s_n+rp1+rp1*rp4))/         &
        !                                                            &(rp1*(-costSTO1s_a+rp1)*(1.d0+rp4)))
        fun=-distp(0,1)*rp4*                                         &
        ((dd1**2*(-costSTO1s_n+rp1+rp1*rp4))/                        &
        (rp1*(-costSTO1s_a+rp1)*(1.d0+rp4)**2))
    !              the second derivative derivative
    fun2=+distp(0,1)*rp4*(dd1**2*(-(costSTO1s_n**2*                  &
        (-1.d0+rp4))-costSTO1s_n*(1.d0+2.d0*rp1)*(1.d0+rp4)          &
        +rp2*(1.d0+rp4)**2)) / (rp2*(1.d0+rp4)**3)
    !    gradient: dR(r)/dr_i=r_i*fun
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    !    laplacian = 2*fun+fun2
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  !        endif

  indpar=indpar+1
  indshell=indshellp






case (123)
!     2p  double  exp
!       dd1 * x_mu  (exp(-dd2 r)+dd3 * exp(-dd4*r))


dd2=dd(indpar+1)
dd3=dd(indpar+2)
dd4=dd(indpar+3)


do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))
  distp(k,2)=dexp(-dd4*r(k))
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*(distp(i,1)+dd3*distp(i,2))
  end do
  !           endif
end do


if(typec.ne.1) then

  fun0=distp(0,1)+dd3*distp(0,2)
  fun=-(dd2*distp(0,1)+dd3*dd4*distp(0,2))/r(0)
  fun2=dd2**2*distp(0,1)+dd3*dd4**2*distp(0,2)


  !              indorbp=indorb

  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do



  !endif for indt
end if

indpar=indpar+3
indshell=indshell+3
indorb=indorbp

case (87)
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

case (47)
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

case (104)
!     2p  double gaussian
!       dd1 * x_mu  (exp(-dd2 r^2)+dd3 * exp(-dd4*r^2))


dd2=dd(indpar+1)
dd3=dd(indpar+2)
dd4=dd(indpar+3)


do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k)**2)
  distp(k,2)=dexp(-dd4*r(k)**2)
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*(distp(i,1)+dd3*distp(i,2))
  end do
  !           endif
end do


if(typec.ne.1) then
  fun0=(distp(0,1)+dd3*distp(0,2))
  fun=2.d0*(-dd2*distp(0,1)                                          &
      -dd4*dd3*distp(0,2))
  fun2=2.d0*(dd2*(-1.d0+2.d0*dd2*r(0)**2)*distp(0,1)                 &
      +dd4*dd3*(-1.d0+2.d0*dd4*r(0)**2)*distp(0,2))


  !              indorbp=indorb

  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do

  !endif for indt
end if

indpar=indpar+3
indshell=indshell+3
indorb=indorbp

case (199)
! derivative of 200 LA COSTANTE

indorbp=indorb+1
indshellp=indshell+1

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=0.d0
end do
!           endif

if(typec.ne.1) then
  do i=1,3
    z(indorbp,indt+i)=0.d0
  end do

  z(indorbp,indt+4)=0
  !endif for indt
end if

indshell=indshellp
indorb=indorbp
case (11)
  ! R(r)=r**2*(exp(-z1*r)+p*exp(-z2*r))



  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd(indpar+3)
  !           if(iflagnorm.gt.2) then
  c=1.d0/2.d0/dsqrt(pi*720.d0*(1.d0/(2.d0*dd1)**7+                   &
      2.d0*peff/(dd1+dd2)**7+peff**2/(2.d0*dd2)**7))
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=c*dexp(-dd2*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=(distp(i,1)+peff*distp(i,2))*r(i)**2
  end do

  if(typec.ne.1) then
    rp1=r(0)**2
    !              the first derivative
    fun=distp(0,1)*(2.d0*r(0)-dd1*rp1)                               &
        +peff*distp(0,2)*(2.d0*r(0)-dd2*rp1)
    !
    !              the second derivative
    fun2=distp(0,1)*(2.d0-4.d0*dd1*r(0)+dd1**2*rp1)                  &
        +peff*distp(0,2)*(2.d0-4.d0*dd2*r(0)+dd2**2*rp1)
    !
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)/r(0)
    end do

    z(indorbp,indt+4)=2.d0*fun/r(0)+fun2

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+3
  indshell=indshellp


  ! 4s single zeta
case (39)
  ! R(r)=r**3*exp(-z1*r)
  !
  indshellp=indshell+1



  !       if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  !           if(iflagnorm.gt.2) then
  !           c=dsqrt((2*dd1)**9/40320.d0/pi)/2.d0
  !           c=dsqrt((2*dd1)**7/720.d0/pi)/2.d0
  c=dd1**3.5d0*0.11894160774351807429d0
  !           c=-c
  !           endif

  c0=-c
  c1=3.5d0*c/dd1

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=(c0*r(i)**3+c1*r(i)**2)*distp(i,1)
  end do

  if(typec.ne.1) then
    rp1=r(0)**3
    rp2=r(0)**2

    !              fun=(2.d0-dd1*r(0))*distp(0,1)
    !              fun2=(2.d0-4*dd1*r(0)+(dd1*r(0))**2)*distp(0,1)
    !
    !c              the first derivative/r
    fun=distp(0,1)*(c0*(3.d0*r(0)-dd1*rp2)                           &
        +c1*(2.d0-dd1*r(0)))

    !c

    !c              the second derivative
    fun2=distp(0,1)*                                                 &
        (c0*(6.d0*r(0)-6.d0*dd1*rp2+dd1**2*rp1)                      &
        +c1*(2.d0-4*dd1*r(0)+(dd1*r(0))**2))
    !c
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*fun+fun2

  end if
  !
  indorb=indorbp
  !
  !        endif
  indpar=indpar+1
  indshell=indshellp
  !
  ! 3p single zeta
case (132)
!     2s with cusp condition
!     ( r^3*exp(-dd2*r))  ! with no cusp condition


dd2=dd(indpar+1)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))*r(k)
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)*r(i)**2
end do
!           endif


if(typec.ne.1) then
  fun=(3.d0-dd2*r(0))*distp(0,1)
  fun2=(6.d0-6*dd2*r(0)+(dd2*r(0))**2)*distp(0,1)
  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*fun+fun2


  !endif for indt
end if

indpar=indpar+1
indshell=indshellp
indorb=indorbp


case (30)
  !     3d without cusp and one parmater


  dd1=dd(indpar+1)
  !         if(iflagnorm.gt.2) then
  !     c=                                                           &
      !                                                              &1.d0/(2.d0**3*3.d0)*dsqrt(1.d0/pi)*(2.d0*dd1)**(7.d0/2.d0)
      c=dd1**3.5d0*0.26596152026762178d0
  !         endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
  end do

  do i=indtmin,indtm
    distp(i,3)=distp(i,1)
    ! lz=0
    distp(i,4)=(3.d0*rmu(3,i)**2-r(i)**2)*cost1d
    ! lz=+/-2
    distp(i,5)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d
    ! lz=+/-2
    distp(i,6)=rmu(1,i)*rmu(2,i)*cost3d
    ! lz=+/-1
    distp(i,7)=rmu(2,i)*rmu(3,i)*cost3d
    ! lz=+/-1
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
    fun=-dd1*distp(0,1)
    fun2=dd1**2*distp(0,1)

    !              indorbp=indorb

    do ic=1,5
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0)                     &
            *fun/r(0)
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
      z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0)+fun2)
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

case (73)


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
case (65)
  ! d orbitals
  ! R(r)= c0 r^3  exp(-alpha r^2)+ c1 r exp(-alpha r^2)
  ! each gaussian term is normalized



  !        indorbp=indorb
  indparp=indpar+1

  dd1=dd(indparp)

  !        if(iflagnorm.gt.2) then
  ! overall normalization to be done
  !        c=8.d0/dsqrt(21.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(9.d0/4.d0)
  c=dd1**2.25d0*1.24420067280413253d0
  !        endif

  c0=-c

  c1=2.25d0*c/dd1

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
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
      z(indorbp,k)=distp(k,1)*(c0*distp(k,1+ic)*r(k)**3+             &
          c1*r(k))
    end do
    !           endif
  end do


  if(typec.ne.1) then

    dd1=dd(indparp)

    rp1=2.d0*dd1*r(0)
    rp2=rp1*r(0)
    fun0=distp(0,1)*(c1*r(0)+c0*r(0)**3)
    fun=(c1*(1.d0-rp2)+c0*r(0)**2*(3.d0-rp2))                        &
        *distp(0,1)/r(0)

    fun2=distp(0,1)*(c1*rp1*(rp2-3.d0)+c0*r(0)                       &
        *(3.d0-3.5d0*rp2+0.5d0*rp2**2))


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




  ! ******************* END GAUSSIAN BASIS ************************

  ! ** ** ** ** ** ** **  JASTROW ORBITALS ** ** ** ** ** ** ** ** *
case (82)

  dd1=dd(indpar+1)
  dd2=dsqrt(dd1)

  !        if(iflagnorm.gt.2) then
  !        ratiocp--> ratiocp*dsqrt(2.d0)*pi**(-0.75d0)*2**1.25
  !        c=dsqrt(2.d0)*pi**(-0.75d0)*(2.d0*dd1)**1.25d0*ratiocp
  c=dd1**1.25d0*ratiocp
  !        endif

  do k=indtmin,indtm
    cost=dd1*r(k)**2/(1.d0+dd2*r(k))
    distp(k,1)=c*dexp(-cost)
  end do
  ! indorbp=indorb
  !
  do ic=1,3
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,1)
    end do
  end do

  if(typec.ne.1) then
    fun0=distp(0,1)
    ! fun=-2.d0*dd1*distp(0,1)
    ! fun2=fun*(1.d0-2.d0*dd1*r(0)**2)
    rp1=dd1*r(0)**2
    rp2=dd2*r(0)
    rp3=(1.d0+rp2)**2

    fun=-dd1*distp(0,1)*(2.d0+rp2)/rp3

    ! the second derivative
    fun2=dd1*distp(0,1)*(-2.d0-2.d0*rp2+4.d0*rp1+4.d0*rp1*rp2+rp1**2)/rp3**2
    ! indorbp=indorb
    do ic=1,3
      ! if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
    end do
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp

case (111)
!     2p single   r_mu/(1+b r^3)   parent of 103



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=1.d0/(1.d0+dd2*r(k)**3)
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)
  end do
  !           endif
end do


if(typec.ne.1) then
  fun0=distp(0,1)
  fun=-dd2*distp(0,1)**2*3.d0*r(0)
  fun2=fun*distp(0,1)*(2.d0-4.d0*dd2*r(0)**3)

  !              indorbp=indorb

  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do


  !endif for indt
end if

indpar=indpar+1
indshell=indshell+3
indorb=indorbp
case (62)



  dd1=dd(indpar+1)


  !        if(iflagnorm.gt.2) then
  !        c=2.d0/pi**0.75d0*(2.d0*dd1)**1.75d0/dsqrt(5.d0)
  c=dd1**1.75d0*1.2749263037197753d0
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  !        indorbp=indorb
  !
  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,1)*r(i)
    end do
    !           endif
  end do

  if(typec.ne.1) then
    fun0=distp(0,1)*r(0)
    cost=2.d0*dd1*r(0)**2
    fun=distp(0,1)*(1.d0-cost)/r(0)
    fun2=2.d0*dd1*fun0*(cost-3.d0)
    !              indorbp=indorb

    do ic=1,3
      !                if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
      !                 endif
    end do

    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp

  ! derivative of 62 with respect zeta
case (42)
  !     4d without cusp and one parmater derivative of 30


  dd1=dd(indpar+1)
  !         if(iflagnorm.gt.2) then
  !     c=                                                           &
      !                                                              &1.d0/(2.d0**3*3.d0)*dsqrt(1.d0/pi)*(2.d0*dd1)**(7.d0/2.d0)
      c=dd1**3.5d0*0.26596152026762178d0
  !     c=
  !                                                                  &1.d0/(2.d0**3*3.d0)/dsqrt(56.d0*pi)*(2.d0*dd1)**(9.d0/2.d0)
      !         endif

  c0=-c
  c1=3.5d0*c/dd1

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
  end do

  do i=indtmin,indtm
    distp(i,3)=distp(i,1)*(c0*r(i)+c1)
    ! lz=0
    distp(i,4)=(3.d0*rmu(3,i)**2-r(i)**2)*cost1d
    ! lz=+/
    distp(i,5)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d
    ! lz=+/-2
    distp(i,6)=rmu(1,i)*rmu(2,i)*cost3d
    ! lz=+/-1
    distp(i,7)=rmu(2,i)*rmu(3,i)*cost3d
    ! lz=+/-1
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
    fun=-dd1*distp(0,3)+c0*distp(0,1)
    fun2=dd1**2*distp(0,3)-2.d0*dd1*c0*distp(0,1)

    !              indorbp=indorb

    do ic=1,5
      !                 if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0)                     &
            *fun/r(0)
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
      z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0)+fun2)
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

case (4)
  ! normalized

  indshellp=indshell+1


  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  dd2=dd(indpar+2)

  !           if(iflagnorm.gt.2) then
  c=dd1**2.5d0/dsqrt(3.d0*pi*(1.d0+dd2**2/3.d0))
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
  end do


  do i=i0,indtm
    z(indorbp,i)=(r(i)+dd2*rmu(3,i))*distp(i,1)
  end do

  if(typec.ne.1) then

    fun=distp(0,1)*(1.d0-dd1*r(0))
    funp=-dd2*dd1*distp(0,1)*rmu(3,0)
    fun2=distp(0,1)*(dd1**2*r(0)-2.d0*dd1)
    fun2p=dd1**2*dd2*distp(0,1)*rmu(3,0)

    do i=1,3
      z(indorbp,indt+i)=(fun+funp)*rmu(i,0)/r(0)
    end do
    z(indorbp,indt+3)=z(indorbp,indt+3)+dd2*distp(0,1)
    z(indorbp,indt+4)=(2.d0*fun+4.d0*funp)/r(0)                      &
        +(fun2+fun2p)

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+2
  indshell=indshellp


  ! 2s single Z NO CUSP
case (137)
!     2s with cusp condition
!     dd1*(exp(-dd2*r)*(1+dd2*r))


dd2=dd(indpar+1)

!           if(iflagnorm.gt.2) then
!           c=1.d0/dsqrt(1/4.d0/dd2**3+12*dd2/(2.d0*dd2)**4+
!                                                                    &3*dd2**2/4/dd2**5)/dsqrt(4.0*pi)
    !           endif

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)*(1.d0+dd2*r(i))
end do
!           endif


if(typec.ne.1) then
  fun=-dd2**2*distp(0,1)
  fun2=fun*(1.d0-dd2*r(0))


  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*fun+fun2


  !endif for indt
end if

indpar=indpar+1
indshell=indshellp
indorb=indorbp

case (8)
  ! normalized
  ! exp(-dd1*r) + (dd1-zeta) * r * exp(-dd2*r)


  indshellp=indshell+1


  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd1-zeta(1)

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
    distp(k,2)=dexp(-dd2*r(k))
  end do

  !           if(iflagnorm.gt.2) then
  c=1.d0/dsqrt(1/4.d0/dd1**3+12*peff/(dd1+dd2)**4+                   &
      3*peff**2/4/dd2**5)/dsqrt(4.0*pi)
  !           endif

  do i=i0,indtm
    z(indorbp,i)=c*(distp(i,1)+r(i)*distp(i,2)*peff)
  end do

  if(typec.ne.1) then

    fun=-dd1*distp(0,1)+peff*distp(0,2)*(1.d0-dd2*r(0))
    fun2=distp(0,1)*dd1**2                                           &
        +peff*distp(0,2)*(dd2**2*r(0)-2.d0*dd2)

    do i=1,3
      z(indorbp,indt+i)=fun*c*rmu(i,0)/r(0)
    end do
    z(indorbp,indt+4)=c*(2.d0*fun/r(0)+fun2)

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+2
  indshell=indshellp

  ! 3s single zeta
case (109)
!     2p  double  Lorentian
!       dd1 * x_mu  (L(dd2 r^2)+dd3 * L(dd4*r^2)) ; L(x)=1/(1+x^2)



dd2=dd(indpar+1)
dd3=dd(indpar+2)
dd4=dd(indpar+3)

do k=indtmin,indtm
  distp(k,1)=1.d0/(1.d0+dd2*r(k)**2)
  distp(k,2)=1.d0/(1.d0+dd4*r(k)**2)
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*(distp(i,1)+dd3*distp(i,2))
  end do
  !           endif
end do


if(typec.ne.1) then
  fun0=distp(0,1)+dd3*distp(0,2)

  fun=2.d0*(-dd2*distp(0,1)**2-dd4*dd3*distp(0,2)**2)
  !     fun2=2.d0*(dd2*(-1.d0+2.d0*dd2*r(0)**2)*distp(0,1)
  !    1+dd4*dd3*(-1.d0+2.d0*dd4*r(0)**2)*distp(0,2))

  fun2=2*dd2*distp(0,1)**3*(-1.d0+3.d0*dd2*r(0)**2)                  &
      +2*dd3*dd4*distp(0,2)**3*(-1.d0+3.d0*dd4*r(0)**2)

  !              indorbp=indorb

  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do

  !endif for indt
end if

indpar=indpar+3
indshell=indshell+3
indorb=indorbp

case (112)
!     2p single   r_mu/(1+b r)^3   parent of 103



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=1.d0/(1.d0+dd2*r(k))**3
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)
  end do
  !           endif
end do



if(typec.ne.1) then

  fun0=distp(0,1)
  fun=-3.d0*dd2*distp(0,1)/(r(0)*(1.d0+dd2*r(0)))
  fun2=12.d0*dd2**2/(1.+dd2*r(0))**5

  !              indorbp=indorb

  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do

  !endif for indt
end if

indpar=indpar+1
indshell=indshell+3
indorb=indorbp

case (151)
!     2p single exponential  -r^3  e^{-z r^2}  ! parent of 150



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k)**2)
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=-rmu(ic,i)*distp(i,1)*r(i)**3
  end do
  !           endif
end do


if(typec.ne.1) then

  fun0=-distp(0,1)*r(0)**3
  cost=dd2*r(0)**2
  fun=distp(0,1)*(-3.d0+2.d0*cost)*r(0)
  fun2=-2.d0*distp(0,1)*r(0)*(3.d0-7.d0*cost+2.d0*cost**2)
  !              indorbp=indorb
  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do
  !endif for indt
end if
indpar=indpar+1
indshell=indshell+3
indorb=indorbp



case (127)
!     3d without cusp and one parmater


dd1=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd1*r(k))
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
  fun=-dd1*distp(0,1)
  fun2=dd1**2*distp(0,1)

  !              indorbp=indorb

  do ic=1,5
    !                 if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=distp(0,3+ic)*rmu(i,0)                       &
          *fun/r(0)
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
    z(indorbp,indt+4)=distp(0,3+ic)*(6.d0*fun/r(0)+fun2)
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

case (139)
!     2s with cusp condition
!     ( r^3*exp(-dd2*r))  !  der of 128


dd2=dd(indpar+1)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=-dexp(-dd2*r(k))*r(k)
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)*r(i)**2
end do
!           endif


if(typec.ne.1) then
  fun=(3.d0-dd2*r(0))*distp(0,1)
  fun2=(6.d0-6*dd2*r(0)+(dd2*r(0))**2)*distp(0,1)
  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*fun+fun2


  !endif for indt
end if

indpar=indpar+1
indshell=indshellp
indorb=indorbp


case (45,69)
  ! d orbitals
  ! R(r)= c*exp(-z r^2)*(7/4/z-r^2)



  !        indorbp=indorb
  indparp=indpar+1
  dd1=dd(indparp)

  !        if(iflagnorm.gt.2) then
  ! overall normalization
  !        c=4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)
  c=dd1**1.75d0*1.64592278064948967213d0
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
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
      z(indorbp,k)=distp(k,1)*(7.d0/4.d0/dd1-r(k)**2)*               &
          distp(k,1+ic)
    end do
    !           endif
  end do


  if(typec.ne.1) then

    dd1=dd(indparp)
    fun0=distp(0,1)*(7.d0/4.d0/dd1-r(0)**2)
    fun=distp(0,1)*(2.d0*dd1*r(0)**2-11.d0/2.d0)
    fun2=distp(0,1)*(-4.d0*dd1**2*r(0)**4                            &
        +17.d0*dd1*r(0)**2-11.d0/2.d0)


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

  ! derivative of 17 with respect to z
case (1200:1299)
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


case (90:99)

    indshellp=indshell+1
    indorbp=indorb+1
    dd1=dd(indpar+1)

    multiplicity = (iopt - 90 + 2) * (iopt - 90 + 1) / 2

    powers(:,0,:) = 1.0d0

    do i = 1, max_power
        do k = indtmin, indtm
            powers(1, i, k) = powers(1, i-1, k) * rmu(1, k)
            powers(2, i, k) = powers(2, i-1, k) * rmu(2, k)
            powers(3, i, k) = powers(3, i-1, k) * rmu(3, k)
        end do
    end do


    c = 0.712705470354990_8 * dd1 ** 0.75_8! * 2.829
    if (iopt - 90 .ne. 0) then
        c = c * (8_4 * dd1) ** ((iopt - 90)/2.0_8)
    end if
    do k = i0, indtm
        distp(k,1) = dexp(-1.0_8 * dd1 * r(k) * r(k)) * c
    end do
    do k = i0, indtm
        count = 0
        do ii = (iopt - 90), 0, -1
            do jj = (iopt - 90) - ii, 0, -1
                kk = (iopt - 90) - ii - jj
                z(indorbp + count, k) = 1.0_8
                rp1 = 1.0_8
                do i = ii + 1, 2 * ii
                    rp1 = rp1 * i
                end do
                z(indorbp + count, k) = z(indorbp + count, k) / dsqrt(rp1)
                rp1 = 1.0_8
                do i = jj + 1, 2 * jj
                    rp1 = rp1 * i
                end do
                z(indorbp + count, k) = z(indorbp + count, k) / dsqrt(rp1)
                rp1 = 1.0_8
                do i = kk + 1, 2 * kk
                    rp1 = rp1 * i
                end do
                z(indorbp + count, k) = z(indorbp + count, k) / dsqrt(rp1)
                count = count + 1
            end do
        end do
    end do

    ! We need to calculate it again for derivatives, it could not be done in previous loop because of case if i0 /= indtmin
    if (typec .ne. 1) then
        count = 0
        do ii = (iopt - 90), 0, -1
            do jj = (iopt - 90) - ii, 0, -1
                kk = (iopt - 90) - ii - jj
                z(indorbp + count, indt + 1) = 1.0_8
                z(indorbp + count, indt + 2) = 1.0_8
                z(indorbp + count, indt + 3) = 1.0_8
                z(indorbp + count, indt + 4) = 1.0_8
                rp1 = 1.0_8
                do i = ii + 1, 2 * ii
                    rp1 = rp1 * i
                end do
                z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) / dsqrt(rp1)
                z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) / dsqrt(rp1)
                z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) / dsqrt(rp1)
                z(indorbp + count, indt + 4) = z(indorbp + count, indt + 4) / dsqrt(rp1)
                rp1 = 1.0_8
                do i = jj + 1, 2 * jj
                    rp1 = rp1 * i
                end do
                z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) / dsqrt(rp1)
                z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) / dsqrt(rp1)
                z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) / dsqrt(rp1)
                z(indorbp + count, indt + 4) = z(indorbp + count, indt + 4) / dsqrt(rp1)
                rp1 = 1.0_8
                do i = kk + 1, 2 * kk
                    rp1 = rp1 * i
                end do
                z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) / dsqrt(rp1)
                z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) / dsqrt(rp1)
                z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) / dsqrt(rp1)
                z(indorbp + count, indt + 4) = z(indorbp + count, indt + 4) / dsqrt(rp1)
                count = count + 1
            end do
        end do
    end if

    ! Initialize gradients and laplacians (radial part)

    if (typec .ne. 1) then
        distp(indt + 1, 1) = -2.0d0 * dd1 * rmu(1, 0) * distp(0, 1)
        distp(indt + 2, 1) = -2.0d0 * dd1 * rmu(2, 0) * distp(0, 1)
        distp(indt + 3, 1) = -2.0d0 * dd1 * rmu(3, 0) * distp(0, 1)
        distp(indt + 4, 1) = dd1 * (4.0d0 * dd1 * (r(0) * r(0)) - 6.0d0) * distp(0, 1)
    end if

    do k = i0, indtm
        count = 0
        do ii = (iopt - 90), 0, -1
            do jj = (iopt - 90) - ii, 0, -1
                kk = (iopt - 90) - ii - jj
                z(indorbp + count, k) = z(indorbp + count, k) * powers(1, ii, k)
                z(indorbp + count, k) = z(indorbp + count, k) * powers(2, jj, k)
                z(indorbp + count, k) = z(indorbp + count, k) * powers(3, kk, k)
                count = count + 1
            end do
        end do
    end do

    if (typec .ne. 1) then
        ! Solve ang_mom = 0, 1 separately
        if (iopt - 90 .eq. 0) then
            z(indorbp, indt + 1) = distp(indt + 1, 1) 
            z(indorbp, indt + 2) = distp(indt + 2, 1)
            z(indorbp, indt + 3) = distp(indt + 3, 1)
            z(indorbp, indt + 4) = distp(indt + 4, 1)
        else if (iopt - 90 .eq. 1) then
            rp1 = dsqrt(2.0_8)
            z(indorbp    , indt + 1) = (distp(indt + 1, 1) * rmu(1, indtmin) + distp(0, 1)) / rp1
            z(indorbp    , indt + 2) = (distp(indt + 2, 1) * rmu(1, indtmin)) / rp1
            z(indorbp    , indt + 3) = (distp(indt + 3, 1) * rmu(1, indtmin)) / rp1

            z(indorbp + 1, indt + 1) = (distp(indt + 1, 1) * rmu(2, indtmin)) / rp1
            z(indorbp + 1, indt + 2) = (distp(indt + 2, 1) * rmu(2, indtmin) + distp(0, 1)) / rp1
            z(indorbp + 1, indt + 3) = (distp(indt + 3, 1) * rmu(2, indtmin)) / rp1

            z(indorbp + 2, indt + 1) = (distp(indt + 1, 1) * rmu(3, indtmin)) / rp1
            z(indorbp + 2, indt + 2) = (distp(indt + 2, 1) * rmu(3, indtmin)) / rp1
            z(indorbp + 2, indt + 3) = (distp(indt + 3, 1) * rmu(3, indtmin) + distp(0, 1)) / rp1

            z(indorbp    , indt + 4) = (distp(indt + 4, 1) * rmu(1, indtmin) + 2.0d0 * distp(indt + 1, 1)) / rp1
            z(indorbp + 1, indt + 4) = (distp(indt + 4, 1) * rmu(2, indtmin) + 2.0d0 * distp(indt + 2, 1)) / rp1
            z(indorbp + 2, indt + 4) = (distp(indt + 4, 1) * rmu(3, indtmin) + 2.0d0 * distp(indt + 3, 1)) / rp1
        else if (iopt - 90 .eq. 2) then
            rp1 = 2.0_8
            rp2 = dsqrt(12.0_8)
            z(indorbp    , indt + 1) = (distp(indt + 1, 1) * rmu(1, indtmin) * rmu(1, indtmin)  + 2 * rmu(1, indtmin) * distp(0, 1)) / rp2
            z(indorbp    , indt + 2) = (distp(indt + 2, 1) * rmu(1, indtmin) * rmu(1, indtmin)) / rp2
            z(indorbp    , indt + 3) = (distp(indt + 3, 1) * rmu(1, indtmin) * rmu(1, indtmin)) / rp2

            z(indorbp + 1, indt + 1) = (distp(indt + 1, 1) * rmu(1, indtmin) * rmu(2, indtmin)  + rmu(2, indtmin) * distp(0, 1)) / rp1
            z(indorbp + 1, indt + 2) = (distp(indt + 2, 1) * rmu(1, indtmin) * rmu(2, indtmin)  + rmu(1, indtmin) * distp(0, 1)) / rp1
            z(indorbp + 1, indt + 3) = (distp(indt + 3, 1) * rmu(1, indtmin) * rmu(2, indtmin)) / rp1

            z(indorbp + 2, indt + 1) = (distp(indt + 1, 1) * rmu(1, indtmin) * rmu(3, indtmin)  + rmu(3, indtmin) * distp(0, 1)) / rp1
            z(indorbp + 2, indt + 2) = (distp(indt + 2, 1) * rmu(1, indtmin) * rmu(3, indtmin)) / rp1
            z(indorbp + 2, indt + 3) = (distp(indt + 3, 1) * rmu(1, indtmin) * rmu(3, indtmin)) + rmu(1, indtmin) * distp(0, 1)/ rp1

            z(indorbp + 3, indt + 1) = (distp(indt + 1, 1) * rmu(2, indtmin) * rmu(2, indtmin)) / rp2
            z(indorbp + 3, indt + 2) = (distp(indt + 2, 1) * rmu(2, indtmin) * rmu(2, indtmin)  + 2 * rmu(2, indtmin) * distp(0, 1)) / rp2
            z(indorbp + 3, indt + 3) = (distp(indt + 3, 1) * rmu(2, indtmin) * rmu(2, indtmin)) / rp2

            z(indorbp + 4, indt + 1) = (distp(indt + 1, 1) * rmu(2, indtmin) * rmu(3, indtmin)) / rp1
            z(indorbp + 4, indt + 2) = (distp(indt + 2, 1) * rmu(2, indtmin) * rmu(3, indtmin)  + rmu(3, indtmin) * distp(0, 1)) / rp1
            z(indorbp + 4, indt + 3) = (distp(indt + 3, 1) * rmu(2, indtmin) * rmu(3, indtmin)  + rmu(2, indtmin) * distp(0, 1)) / rp1

            z(indorbp + 5, indt + 1) = (distp(indt + 1, 1) * rmu(3, indtmin) * rmu(3, indtmin)) / rp2
            z(indorbp + 5, indt + 2) = (distp(indt + 2, 1) * rmu(3, indtmin) * rmu(3, indtmin)) / rp2
            z(indorbp + 5, indt + 3) = (distp(indt + 3, 1) * rmu(3, indtmin) * rmu(3, indtmin)  + 2 * rmu(3, indtmin) * distp(0, 1)) / rp2

            z(indorbp    , indt + 4) = (1.0d0 * distp(indt + 4, 1) * rmu(1, indtmin) * rmu(1, indtmin)&
                                   & +  4.0d0 * distp(indt + 1, 1) * rmu(1, indtmin)&
                                   & +  2.0d0 * distp(0, 1)) / rp2
            z(indorbp + 1, indt + 4) = (1.0d0 * distp(indt + 4, 1) * rmu(1, indtmin) * rmu(2, indtmin)&
                                   & +  2.0d0 * distp(indt + 2, 1) * rmu(1, indtmin)&
                                   & +  2.0d0 * distp(indt + 1, 1) * rmu(2, indtmin)) / rp1
            z(indorbp + 2, indt + 4) = (1.0d0 * distp(indt + 4, 1) * rmu(1, indtmin) * rmu(3, indtmin)&
                                   & +  2.0d0 * distp(indt + 3, 1) * rmu(1, indtmin)&
                                   & +  2.0d0 * distp(indt + 1, 1) * rmu(3, indtmin)) / rp1
            z(indorbp + 3, indt + 4) = (1.0d0 * distp(indt + 4, 1) * rmu(2, indtmin) * rmu(2, indtmin)&
                                   & +  4.0d0 * distp(indt + 2, 1) * rmu(2, indtmin)&
                                   & +  2.0d0 * distp(0, 1)) / rp2
            z(indorbp + 4, indt + 4) = (1.0d0 * distp(indt + 4, 1) * rmu(2, indtmin) * rmu(3, indtmin)&
                                   & +  2.0d0 * distp(indt + 3, 1) * rmu(2, indtmin)&
                                   & +  2.0d0 * distp(indt + 2, 1) * rmu(3, indtmin)) / rp1
            z(indorbp + 5, indt + 4) = (1.0d0 * distp(indt + 4, 1) * rmu(3, indtmin) * rmu(3, indtmin)&
                                   & +  4.0d0 * distp(indt + 3, 1) * rmu(3, indtmin)&
                                   & +  2.0d0 * distp(0, 1)) / rp2
        else
            count = 0
            do ii = (iopt - 90), 0, -1
                do jj = (iopt - 90) - ii, 0, -1
                    kk = (iopt - 90) - ii - jj

                    ! First store polynomial part into respective places
                    ! Then solve full laplacian using using lower derivatives
                    ! Then do the same thing for gradients
                    ! Then finally the values

                    z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) * powers(1, ii-1, 0)
                    z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) * powers(2, jj, 0)
                    z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) * powers(3, kk, 0)
                    z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) * ii

                    z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) * powers(1, ii, 0)
                    z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) * powers(2, jj-1, 0)
                    z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) * powers(3, kk, 0)
                    z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) * jj

                    z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) * powers(1, ii, 0)
                    z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) * powers(2, jj, 0)
                    z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) * powers(3, kk-1, 0)
                    z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) * kk

                    z(indorbp + count, indt + 4) = z(indorbp + count, indt + 4) &
                                               & * (powers(1, ii-2, 0) * powers(2, jj, 0) * powers(3, kk, 0) * ii * (ii-1)&
                                               & +  powers(1, ii, 0) * powers(2, jj-2, 0) * powers(3, kk, 0) * jj * (jj-1)&
                                               & +  powers(1, ii, 0) * powers(2, jj, 0) * powers(3, kk-2, 0) * kk * (kk-1))

                     
                    ! All polynomial parts are now stored
                    ! Now solve laplacian
                    z(indorbp + count, indt + 4) =         z(indorbp + count, indt + 4) * distp(0, 1) &
                                               & + 2.0_8 * z(indorbp + count, indt + 1) * distp(indt + 1, 1) &
                                               & + 2.0_8 * z(indorbp + count, indt + 2) * distp(indt + 2, 1) &
                                               & + 2.0_8 * z(indorbp + count, indt + 3) * distp(indt + 3, 1) &
                                               & +         z(indorbp + count, indtmin)  * distp(indt + 4, 1)

                    ! Now solve gradients
                    z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) * distp(0, 1) &
                                               & + z(indorbp + count, indtmin)  * distp(indt + 1, 1)
                    z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) * distp(0, 1) &
                                               & + z(indorbp + count, indtmin)  * distp(indt + 2, 1)
                    z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) * distp(0, 1) &
                                               & + z(indorbp + count, indtmin)  * distp(indt + 3, 1)
                    count = count + 1
                end do
            end do
        end if

    end if

    ! Multiply by radial part for values
    do ii = 1, multiplicity
        do kk = i0, indtm
            z(indorbp + ii - 1, kk) = z(indorbp + ii - 1, kk) * distp(kk, 1)
        end do
    end do

    indpar=indpar + 1
    indshell=indshell + multiplicity
    indorb=indorb + multiplicity
case (117)
!     2s double lorentian with constant  parent of 102
!     (dd3+r^3/(1+dd5*r)^4;



dd3=dd(indpar+1)
dd5=dd(indpar+2)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=r(k)**3/(1.d0+dd5*r(k))**4
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=dd3+distp(i,1)
!              write(6,*) ' function inside = ',z(indorbp,i)
end do
!           endif


if(typec.ne.1) then

  fun=                                                               &
      -r(0)*(-3.d0+dd5*r(0))/(1.d0+dd5*r(0))**5
  fun2=                                                              &
      +2.d0*r(0)*(3.d0-6.d0*dd5*r(0)+(dd5*r(0))**2)                  &
      /(1.d0+dd5*r(0))**6


  !              write(6,*) ' fun inside = ',fun,fun2

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=fun2+2.d0*fun

  !           write(6,*) ' lap 106 =',z(indorbp,indt+4)

  !endif for indt
end if

indpar=indpar+2
indshell=indshellp
indorb=indorbp

case (50)
  ! R(r)=(c0*r**4+c1*r**3)*exp(-z1*r)
  !
  indshellp=indshell+1



  !       if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  !           if(iflagnorm.gt.2) then
  c=dsqrt((2*dd1)**9/40320.d0/pi)/2.d0
  !           endif

  c0=-c
  c1=4.5d0*c/dd1

  do k=indtmin,indtm
    distp(k,1)=r(k)*dexp(-dd1*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=(c0*r(i)**3+c1*r(i)**2)*distp(i,1)
  end do

  if(typec.ne.1) then
    rp1=r(0)*dd1
    rp2=rp1*rp1

    !c              the first derivative/r
    fun=-distp(0,1)*(c0*r(0)*(rp1-4.d0)+c1*(rp1-3.d0))

    !c

    !c              the second derivative
    fun2=distp(0,1)*                                                 &
        (c0*r(0)*(12.d0-8.d0*rp1+rp2)+c1*(6.d0-6*rp1+rp2))
    !c
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*fun+fun2

  end if
  !
  indorb=indorbp
  !
  !        endif
  indpar=indpar+1
  indshell=indshellp
  !



case (3)
  !



  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1

  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd(indpar+3)

  !           if(iflagnorm.gt.2) then
  c=1.d0/2.d0/dsqrt(2.d0*pi*(1.d0/(2.d0*dd1)**3                      &
      +2.d0*peff/(dd1+dd2)**3+peff**2/(2.d0*dd2)**3))
  !           endif

  do i=indpar+1,indpar+2
    do k=indtmin,indtm
      distp(k,i-indpar)=c*dexp(-dd(i)*r(k))
    end do
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)+peff*distp(i,2)
  end do

  if(typec.ne.1) then
    fun=-dd1*distp(0,1)-peff*dd2*distp(0,2)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)/r(0)
    end do

    z(indorbp,indt+4)=(-2.d0*dd1/r(0)+dd1**2)                        &
        *distp(0,1)+peff*(-2.d0*dd2/r(0)+dd2**2)                     &
        *distp(0,2)


  end if

  indorb=indorbp

  !        endif

  indpar=indpar+3
  indshell=indshellp

  ! 2s 2pz Hybryd single Z
case (124)
!     2s double exp  with constant and cusp cond.
!     (dd3+ exp (-dd2 r)*(1+dd2*r)+dd4*exp(-dd5*r)*(1+dd5*r))



dd2=dd(indpar+1)
dd3=dd(indpar+2)
dd4=dd(indpar+3)
dd5=dd(indpar+4)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,3)=dexp(-dd2*r(k))
  distp(k,4)=dexp(-dd5*r(k))
  distp(k,1)=distp(k,3)*(1.d0+dd2*r(k))
  distp(k,2)=distp(k,4)*(1.d0+dd5*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)+dd3+dd4*distp(i,2)
!              write(6,*) ' function inside = ',z(indorbp,i)
end do
!           endif


if(typec.ne.1) then

  fun=-dd2**2*distp(0,3)-dd5**2*dd4*distp(0,4)
  fun2=-dd2**2*distp(0,3)*(1.d0-dd2*r(0))                            &
      -dd4*dd5**2*distp(0,4)*(1.d0-dd5*r(0))



  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*fun+fun2


  !endif for indt
end if

indpar=indpar+4
indshell=indshellp
indorb=indorbp

case (28)
  ! R(r)=(Z*b*r)^4/(1+(Z*b*r)^4)*exp(-z*r)     orbital 1s (no cusp)
  ! d -> b1s (defined in module constants)
  ! normadization: cost1s, depends on b1s


  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1

  dd1=dd(indpar+1)

  !           if(iflagnorm.gt.2) then
  !             if(dd1.gt.0.) then
  c=cost1s*dd1**1.5d0
  !             else
  !               c=1.d0
  !             endif
  !           endif

  do i=indtmin,indtm
    distp(i,1)=c*dexp(-dd1*r(i))
  end do

  do i=i0,indtm
    rp4=(dd1*b1s*r(i))**4
    z(indorbp,i)=distp(i,1)*rp4/(1.d0+rp4)
  end do

  if(typec.ne.1) then
    rp1=dd1*b1s*r(0)
    rp2=rp1**2
    rp4=rp2**2
    rp5=r(0)*dd1
    rp6=(b1s*dd1)**2*rp2
    !              the first derivative /r
    fun=-distp(0,1)*rp6*(-4.d0+rp5+rp4*rp5)/(1.d0+rp4)**2
    !              the second derivative derivative
    fun2=distp(0,1)*rp6*(12.d0-8*rp5+rp5**2-20*rp4-                  &
        8*rp4*rp5+2*rp4*rp5**2+(rp4*rp5)**2)/(1.d0+rp4)**3
    !    gradient: dR(r)/dr_i=r_i*fun
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    !    laplacian = 2*fun+fun2
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  !        endif

  indpar=indpar+1
  indshell=indshellp

case default
write(6,*) 'WARNING makefun: orbital',iopt,'not found'

iflagerr=1

end select
! ** ** ** ** ** ** ** END OF JASTROW ORBITALS ** ** ** ** ** ** ** ** *

return
end
