!TL off
SUBROUTINE NGIVEJ_B(rcart1, rcart1b, rcart2, rcart2b, lbox, ngivejb, cellscaleb, s2rb)
    use Cell, only : cellscale, s2r, car2cry
    IMPLICIT NONE
    REAL*8 :: dxyz(3), ngivej, tempb, cellscaleb(3), s2rb(3, 3)
    REAL*8 :: dxyzb(3), ngivejb
    REAL*8 :: rcart1(3), rcart2(3), lbox
    REAL*8 :: rcart1b(3), rcart2b(3)
    real*8 :: npip(3)
    dxyz(:) = rcart1(:) - rcart2(:)
    IF (lbox .GT. 0.d0) THEN
        npip(:) = dxyz(:)
!       Call CartesianToCrystal(npip, 1)
        npip(:)=car2cry(:,1)*npip(1)+car2cry(:,2)*npip(2)+car2cry(:,3)*npip(3)
        npip(1) = anint(npip(1) / cellscale(1))
        npip(2) = anint(npip(2) / cellscale(2))
        npip(3) = anint(npip(3) / cellscale(3))
        call dgemv('N', 3, 3, -1.d0, s2r, 3, npip, 1, 1.d0, dxyz, 1)
        !      dxyz(1)=dxyz(1)-cellscale(1)*anint(dxyz(1)/cellscale(1))
        !      dxyz(2)=dxyz(2)-cellscale(2)*anint(dxyz(2)/cellscale(2))
        !      dxyz(3)=dxyz(3)-cellscale(3)*anint(dxyz(3)/cellscale(3))
    END IF
    ngivej = DSQRT(dxyz(1)**2 + dxyz(2)**2 + dxyz(3)**2)
    IF (ngivej .LT. 1d-9) ngivej = 0.0_8
    dxyzb(:) = 0.0_8
    IF (ngivej.EQ.0.0) THEN
        tempb = 0.0
    ELSE
        !   tempb = ngivejb/(DSQRT(dxyz(1)**2+dxyz(2)**2+dxyz(3)**2))
        tempb = ngivejb / ngivej
    END IF
    dxyzb(1) = dxyz(1) * tempb
    dxyzb(2) = dxyz(2) * tempb
    dxyzb(3) = dxyz(3) * tempb
    rcart1b(:) = rcart1b(:) + dxyzb(:)
    rcart2b(:) = rcart2b(:) - dxyzb(:)
    if(Lbox.gt.0) then
        !  reverse of      call dgemv('N',3,3,-1.d0,s2r,3,npip,1,1.d0,dxyz,1)
        !  npip fixed
        call dger(3, 3, -1.d0, dxyzb, 1, npip, 1, s2rb, 3)
    endif
END SUBROUTINE NGIVEJ_B
SUBROUTINE DGIVEJ_B(rcart1, rcart1b, rcart2, rcart2b, lbox, dxyzb, cellscaleb, s2rb)
    use Cell, only : cellscale, s2r, car2cry
    IMPLICIT NONE
    REAL*8 :: dxyz(3),  tempb, cellscaleb(3), s2rb(3, 3)
    REAL*8 :: dxyzb(3)
    REAL*8 :: rcart1(3), rcart2(3), lbox
    REAL*8 :: rcart1b(3), rcart2b(3)
    real*8 :: npip(3)
    dxyz(:) = rcart1(:) - rcart2(:)
    IF (lbox .GT. 0.d0) THEN
        npip(:) = dxyz(:)
!       Call CartesianToCrystal(npip, 1)
        npip(:)=car2cry(:,1)*npip(1)+car2cry(:,2)*npip(2)+car2cry(:,3)*npip(3)
        npip(1) = anint(npip(1) / cellscale(1))
        npip(2) = anint(npip(2) / cellscale(2))
        npip(3) = anint(npip(3) / cellscale(3))
        call dgemv('N', 3, 3, -1.d0, s2r, 3, npip, 1, 1.d0, dxyz, 1)
        !      dxyz(1)=dxyz(1)-cellscale(1)*anint(dxyz(1)/cellscale(1))
        !      dxyz(2)=dxyz(2)-cellscale(2)*anint(dxyz(2)/cellscale(2))
        !      dxyz(3)=dxyz(3)-cellscale(3)*anint(dxyz(3)/cellscale(3))
    END IF
    rcart1b(:) = rcart1b(:) + dxyzb(:)
    rcart2b(:) = rcart2b(:) - dxyzb(:)
    if(Lbox.gt.0) then
        !  reverse of      call dgemv('N',3,3,-1.d0,s2r,3,npip,1,1.d0,dxyz,1)
        !  npip fixed
        call dger(3, 3, -1.d0, dxyzb, 1, npip, 1, s2rb, 3)
    endif
    dxyzb=0.d0
END SUBROUTINE DGIVEJ_B
