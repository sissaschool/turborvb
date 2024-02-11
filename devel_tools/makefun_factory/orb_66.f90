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








