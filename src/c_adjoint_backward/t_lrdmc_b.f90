!TL off
subroutine lrdmc(kel, kel_b, nion, rion, rion_b, indvic, alat, ivic, plat, cellscale, cellscale_b, s2r, s2rb, iespbc, t_lrdmc_b)
    use cell, only : yes_tilted, CartesianToCrystal
    use allio, only : zmin, zetar
    implicit none
    integer nion, nel, i
    integer indvic, imu, jmu
    real*8  t_lrdmc_b, alat, plat(*), ivic(3, *), distmu, distmin, rion(3, *), fun     &
            &, rmu(3), srmu(3), cellscale(3), cellscale_b(3), rion_b(3, *)
    real*8 kel(3), kel_b(3), fun_b, distmu_b, rmu_b(3), s2r(3, 3), s2rb(3, 3)
    logical iespbc
    ! the kinetic part on the lattice is required (first 12 ivic)
    !      Compute the minimum distance between the electron
    !      and the nuclei's with maximum Z
    if(plat(1).ne.0.d0) then
        if(indvic.le.6) then
            !        imu=mod(indvic-1,3)+1
            !        t_lrdmc=plat(2)*fun/sum(ivic(:,indvic)**2)
            fun_b = t_lrdmc_b * plat(2) / sum(ivic(:, indvic)**2)
        else
            !        imu=mod(indvic-1,3)+1
            !        t_lrdmc=plat(3)*(1.d0-fun)/sum(ivic(:,indvic)**2)
            fun_b = -t_lrdmc_b * plat(3) / sum(ivic(:, indvic)**2)
        endif
        t_lrdmc_b = 0.d0
        !      Recomputation distmu
        rmu(:) = kel(:) + 0.5d0 * ivic(:, indvic) - rion(:, 1)
        if(iespbc) then
            if(yes_tilted) then
                srmu(:) = rmu(:)
                call CartesianToCrystal(srmu, 1)
                srmu(1) = anint(srmu(1) / cellscale(1))
                srmu(2) = anint(srmu(2) / cellscale(2))
                srmu(3) = anint(srmu(3) / cellscale(3))
                call dgemv('N', 3, 3, -1.d0, s2r, 3, srmu, 1, 1.d0, rmu, 1)
            else
                rmu(:) = rmu(:) - cellscale(:) * anint(rmu(:) / cellscale(:))
            endif
        endif
        distmin = 0.d0
        if(zetar(1).ge.zmin) distmin = sum(rmu(:)**2)
        imu = 1
        !     Calculate min el-ion  distance --> distmu
        do i = 2, nion
            if(zetar(i).ge.zmin) then
                rmu(:) = kel(:) + 0.5d0 * ivic(:, indvic) - rion(:, i)
                if(iespbc) then
                    if(yes_tilted) then
                        srmu(:) = rmu(:)
                        call CartesianToCrystal(srmu, 1)
                        srmu(1) = anint(srmu(1) / cellscale(1))
                        srmu(2) = anint(srmu(2) / cellscale(2))
                        srmu(3) = anint(srmu(3) / cellscale(3))
                        call dgemv('N', 3, 3, -1.d0, s2r, 3, srmu, 1, 1.d0, rmu, 1)
                    else
                        rmu(:) = rmu(:) - cellscale(:) * anint(rmu(:) / cellscale(:))
                    endif
                endif
                distmu = sum(rmu(:)**2)
                if(distmu.lt.distmin.or.distmin.eq.0.d0) then
                    distmin = distmu
                    imu = i
                endif
            endif
        enddo
        distmu = distmin
        if(plat(1).gt.0.d0) then
            !      fun=1.d0/(1.d0+plat(1)*distmu)
            distmu_b = -(1.d0 + plat(1) * distmu)**2 * plat(1) * fun_b
        else
            fun = exp(distmu / plat(1))
            distmu_b = fun / plat(1) * fun_b
        endif
        !        distmu=sum(rmu(:)**2)
        rmu(:) = kel(:) + 0.5d0 * ivic(:, indvic) - rion(:, imu)
        rmu_b(:) = 2.d0 * distmu_b * rmu(:)
        !        if(iespbc) rmu(:)=rmu(:)-cellscale(:)*anint(rmu(:)/cellscale(:))
        if(iespbc) then
            if(yes_tilted) then
                srmu(:) = rmu(:)
                call CartesianToCrystal(srmu, 1)
                srmu(1) = anint(srmu(1) / cellscale(1))
                srmu(2) = anint(srmu(2) / cellscale(2))
                srmu(3) = anint(srmu(3) / cellscale(3))
                !  reverse of         call dgemv('N',3,3,-1.d0,s2r,3,srmu,1,1.d0,rmu,1)
                call dger(3, 3, -1.d0, rmu_b, 1, srmu, 1, s2rb, 3)
            else
                cellscale_b(:) = cellscale_b(:) - rmu_b(:) * anint(rmu(:) / cellscale(:))
            endif
        endif
        !      rmu(:)=kel(:)+0.5d0*ivic(:,indvic)-rion(:,i)
        kel_b(:) = kel_b(:) + rmu_b(:)
        rion_b(:, imu) = rion_b(:, imu) - rmu_b(:)
    else
        t_lrdmc_b = 0.d0
    endif

    return
END
