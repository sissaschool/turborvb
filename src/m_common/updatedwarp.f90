!TL off
!> @brief Update coordinate warping and calculate forces for quantum Monte Carlo
!>
!> This subroutine performs coordinate warping calculations and updates forces
!> for quantum Monte Carlo simulations. It handles both isolated and periodic
!> boundary conditions, with special treatment for electron-ion and ion-ion
!> interactions using warping functions.
!>
!> Parameters
!> ----------
!> NEL : integer, in
!>     Number of electrons.
!> NION : integer, in
!>     Number of ions.
!> KEL : real*8 array, in
!>     Electron positions (3 × NEL).
!> RION : real*8 array, in
!>     Ion positions (3 × NION).
!> RMU : real*8 array, inout
!>     Relative position vectors and their derivatives.
!> R : real*8 array, out
!>     Distances between electrons and ions (NEL × NION).
!> KELELOCB : real*8 array, in
!>     Electron local forces (3 × NEL).
!> KELLOGB : real*8 array, in
!>     Electron log derivatives (3 × NEL).
!> RIONELOCB : real*8 array, in
!>     Ion local forces (3 × NION).
!> RIONLOGB : real*8 array, in
!>     Ion log derivatives (3 × NION).
!> CELLB : real*8 array, inout
!>     Cell forces for pressure calculation.
!> CELLLB : real*8 array, inout
!>     Cell log derivatives for pressure calculation.
!> FORCE : real*8 array, out
!>     Total forces on ions (3 × NION).
!> PULAY : real*8 array, out
!>     Pulay forces on ions (3 × NION).
!> IESPBC : logical, in
!>     Flag for periodic boundary conditions.
!> WARP : logical, in
!>     Flag for coordinate warping.
!> POWER : real*8, in
!>     Power in warping function (1/r^power).
!> ATOM_NUMBER : integer array, in
!>     Atomic numbers of ions.
!> WARPMAT : real*8 array, out
!>     Warping matrix for ghost atoms.
!> NIONG : integer, in
!>     Number of ghost ions.
!>
!> Notes
!> -----
!> - Handles both isolated and periodic boundary conditions.
!> - Uses coordinate warping for improved sampling efficiency.
!> - Calculates forces including Pulay corrections.
!> - Supports pressure calculations for periodic systems.
!> - Distinguishes between real atoms (ATOM_NUMBER > 0) and ghost atoms.
!> - Uses metric tensor for periodic boundary conditions.
!>
!> Algorithm
!> ---------
!> 1. Calculate relative positions (electron-ion and ion-ion)
!> 2. Apply periodic boundary conditions if needed
!> 3. Calculate distances using appropriate metric
!> 4. Compute warping functions and derivatives
!> 5. Update forces with warping contributions
!> 6. Calculate pressure contributions for periodic systems
!>
!> Example
!> -------
!> Used in variational Monte Carlo for force calculations and coordinate warping.
subroutine updatedwarp(nel, nion, kel, rion, rmu, r, kelelocb, kellogb&
        &, rionelocb, rionlogb, cellb, celllb, force, pulay&
        &, iespbc, warp, power, atom_number, warpmat, niong)
    use allio, only: norm_metric
    use Cell
    implicit none
    logical iespbc, warp
    integer nel, nion, i, j, k, n, niong, irefg, ireft
    real*8 kel(3, nel),rion(3, nion),rmu(3,max(nel,nion),nion), r(nel,nion)&
            &, kelelocb(3, nel), kellogb(3, nel), rionelocb(3, nion), rionlogb(3, nion)&
            &, force(3, nion), pulay(3, nion), sumw, sumdw, wfunc, wderiv, cellb(3), celllb(3)&
            &, derpul, wder, xmu(3), atom_number(*), warpmat(nion - niong, *)
    real*8 power, rc

!> @brief Calculate relative positions for periodic boundary conditions
    if(iespbc) then
        do k = 1, nion
            if(atom_number(k).gt.0) then
!> @brief Electron-ion relative positions
                do i = 1, nel
                    rmu(1, i, k) = kel(1, i) - rion(1, k)
                    rmu(2, i, k) = kel(2, i) - rion(2, k)
                    rmu(3, i, k) = kel(3, i) - rion(3, k)
                enddo
                call CartesianToCrystal(rmu(1, 1, k), nel)
                do i = 1, nel
                    rmu(1, i, k) = map(rmu(1, i, k), cellscale(1))
                    rmu(2, i, k) = map(rmu(2, i, k), cellscale(2))
                    rmu(3, i, k) = map(rmu(3, i, k), cellscale(3))
                enddo
            else
!> @brief Ion-ion relative positions for ghost atoms
                do i = 1, nion
                    rmu(1, i, k) = rion(1, k) - rion(1, i)
                    rmu(2, i, k) = rion(2, k) - rion(2, i)
                    rmu(3, i, k) = rion(3, k) - rion(3, i)
                enddo
                call CartesianToCrystal(rmu(1, 1, k), nion)
                do i = 1, nion
                    rmu(1, i, k) = map(rmu(1, i, k), cellscale(1))
                    rmu(2, i, k) = map(rmu(2, i, k), cellscale(2))
                    rmu(3, i, k) = map(rmu(3, i, k), cellscale(3))
                enddo
            endif
        enddo
    else
!> @brief Calculate relative positions for isolated systems
        do k = 1, nion
            if(atom_number(k).gt.0) then
                do i = 1, nel
                    rmu(1, i, k) = kel(1, i) - rion(1, k)
                    rmu(2, i, k) = kel(2, i) - rion(2, k)
                    rmu(3, i, k) = kel(3, i) - rion(3, k)
                enddo
            else
                do i = 1, nion
                    rmu(1, i, k) = rion(1, k) - rion(1, i)
                    rmu(2, i, k) = rion(2, k) - rion(2, i)
                    rmu(3, i, k) = rion(3, k) - rion(3, i)
                enddo
            endif
        enddo
    endif

!> @brief Calculate distances using appropriate metric
    if(iespbc) then
        do k = 1, nion
            do i = 1, nel
                r(i, k) = max(norm_metric(rmu(1, i, k), metric), 1d-9)
            enddo
        enddo
    else
        do k = 1, nion
            do i = 1, nel
                r(i, k) = max(dsqrt(sum(rmu(:, i, k)**2)), 1d-9)
            enddo
        enddo
    endif

!> @brief Calculate derivatives for periodic boundary conditions
    if(iespbc) then
        do k = 1, nion
            if(atom_number(k).gt.0) then
                do i = 1, nel
                    rmu(1, i, k) = kel(1, i) - rion(1, k)
                    rmu(2, i, k) = kel(2, i) - rion(2, k)
                    rmu(3, i, k) = kel(3, i) - rion(3, k)
                enddo
                call CartesianToCrystal(rmu(1, 1, k), nel)
                do i = 1, nel
                xmu(1)=map(rmu(1,i,k),cellscale(1))
                xmu(2)=map(rmu(2,i,k),cellscale(2))
                xmu(3)=map(rmu(3,i,k),cellscale(3))
   rmu(1, i, k)=(metric(1,1)*xmu(1)+metric(1,2)*xmu(2)+metric(1,3)*xmu(3))&
  &*dmap(rmu(1, i, k), cellscale(1))
   rmu(2, i, k) =(metric(2,1)*xmu(1)+metric(2,2)*xmu(2)+metric(2,3)*xmu(3))&
  &*dmap(rmu(2, i, k), cellscale(2))
   rmu(3, i, k) =(metric(3,1)*xmu(1)+metric(3,2)*xmu(2)+metric(3,3)*xmu(3))&
  &*dmap(rmu(3, i, k), cellscale(3))
!> @brief Apply chain rule for coordinate transformation
  xmu(:)=rmu(:,i,k)
  rmu(1,i,k)=xmu(1)*car2cry(1,1)+xmu(2)*car2cry(2,1)+xmu(3)*car2cry(3,1)
  rmu(2,i,k)=xmu(1)*car2cry(1,2)+xmu(2)*car2cry(2,2)+xmu(3)*car2cry(3,2)
  rmu(3,i,k)=xmu(1)*car2cry(1,3)+xmu(2)*car2cry(2,3)+xmu(3)*car2cry(3,3)
                enddo
            else
                do i = 1, nion
                    rmu(1, i, k) = rion(1, k) - rion(1, i)
                    rmu(2, i, k) = rion(2, k) - rion(2, i)
                    rmu(3, i, k) = rion(3, k) - rion(3, i)
                enddo
                call CartesianToCrystal(rmu(1, 1, k), nion)
                do i = 1, nion
                xmu(1)=map(rmu(1,i,k),cellscale(1))
                xmu(2)=map(rmu(2,i,k),cellscale(2))
                xmu(3)=map(rmu(3,i,k),cellscale(3))
   rmu(1, i, k)=(metric(1,1)*xmu(1)+metric(1,2)*xmu(2)+metric(1,3)*xmu(3))&
  &*dmap(rmu(1, i, k), cellscale(1))
   rmu(2, i, k) =(metric(2,1)*xmu(1)+metric(2,2)*xmu(2)+metric(2,3)*xmu(3))&
  &*dmap(rmu(2, i, k), cellscale(2))
   rmu(3, i, k) =(metric(3,1)*xmu(1)+metric(3,2)*xmu(2)+metric(3,3)*xmu(3))&
  &*dmap(rmu(3, i, k), cellscale(3))
!> @brief Apply chain rule for coordinate transformation
  xmu(:)=rmu(:,i,k)
  rmu(1,i,k)=xmu(1)*car2cry(1,1)+xmu(2)*car2cry(2,1)+xmu(3)*car2cry(3,1)
  rmu(2,i,k)=xmu(1)*car2cry(1,2)+xmu(2)*car2cry(2,2)+xmu(3)*car2cry(3,2)
  rmu(3,i,k)=xmu(1)*car2cry(1,3)+xmu(2)*car2cry(2,3)+xmu(3)*car2cry(3,3)
                enddo
            endif
        enddo
    endif

!> @brief Initialize forces from local contributions
    force = rionelocb
    pulay = rionlogb

!> @brief Apply coordinate warping if enabled
    if(warp) then
        do n = 1, nel
            sumw = 0.d0
            if(power.ne.0.d0) then
                do k = 1, nion
                    if(atom_number(k).gt.0.d0) then
                        wfunc = 1.d0 / r(n, k)**power
                        sumw = sumw + wfunc
                    endif
                enddo
            endif
            do j = 1, 3
                sumdw = 0.d0
                if(power.ne.0.d0) then
                    do k = 1, nion
                        if(atom_number(k).gt.0.d0) then
                            wderiv = -power / 2.d0 / r(n, k)**(power + 2) * rmu(j, n, k)
                            sumdw = sumdw + wderiv
                        endif
                    enddo
                endif
                do i = 1, nion
                    if(atom_number(i).gt.0.d0) then
                        if(power.eq.0.d0) then
                            wder = 1.d0
                            derpul = 0.d0
                        else
                            wderiv = -(power / 2.d0) / r(n, i)**(power + 2) * rmu(j, n, i)
                            wfunc = 1.d0 / r(n, i)**power
                            derpul = wderiv / sumw - wfunc / sumw**2 * sumdw
                            wder = wfunc / sumw
                        endif
!> @brief Update forces with warping contributions
                        pulay(j, i) = pulay(j, i) + wder * kellogb(j, n) + derpul
                        force(j, i) = force(j, i) + wder * kelelocb(j, n)
                    endif
                enddo
            enddo
        enddo
!> @brief Handle ghost atoms for warping
        irefg = 0
        do n = 1, nion
            if(atom_number(n).le.0) then
                irefg = irefg + 1
                sumw = 0.d0
                if(power.ne.0.d0) then
                    do k = 1, nion
                        if(atom_number(k).gt.0.d0) then
                            wfunc = 1.d0 / r(k, n)**power
                            sumw = sumw + wfunc
                        endif
                    enddo
                endif
                do j = 1, 3
                    ireft = 0
                    do i = 1, nion
                        if(atom_number(i).gt.0.d0) then
                            ireft = ireft + 1
                            if(power.ne.0.d0) then
                                wfunc = 1.d0 / r(i, n)**power
                                wder = wfunc / sumw
                            else
                                wder = 1.d0
                            endif
!> @brief Update forces for ghost atoms
                            pulay(j, i) = pulay(j, i) + wder * rionlogb(j, n)
                            force(j, i) = force(j, i) + wder * rionelocb(j, n)
                            warpmat(ireft, irefg) = wder
                        endif
                    enddo
                enddo
            endif
        enddo
    endif
!> @brief Calculate pressure contributions for periodic systems
    if(iespbc) then
        rmu(:,1:nion,1)=rion(:,1:nion)
        call CartesianToCrystal(rmu,nion)

        do k = 1, nion
            cellb(1) = cellb(1) + rionelocb(1, k) * rmu(1, k,1) / cellscale(1)
            celllb(1) = celllb(1) + rionlogb(1, k) * rmu(1, k,1) / cellscale(1)

            cellb(2) = cellb(2) + rionelocb(2, k) * rmu(2, k,1) / cellscale(2)
            celllb(2) = celllb(2) + rionlogb(2, k) * rmu(2, k,1) / cellscale(2)

            cellb(3) = cellb(3) + rionelocb(3, k) * rmu(3, k,1) / cellscale(3)
            celllb(3) = celllb(3) + rionlogb(3, k) * rmu(3, k,1) / cellscale(3)
        enddo
        rmu(:,1:nel,1)=kel(:,1:nel)
        call CartesianToCrystal(rmu,nel)

        do k = 1, nel
            cellb(1) = cellb(1) + kelelocb(1, k) * rmu(1, k,1) / cellscale(1)
            celllb(1) = celllb(1) + kellogb(1, k) * rmu(1, k,1) / cellscale(1)

            cellb(2) = cellb(2) + kelelocb(2, k) * rmu(2, k,1) / cellscale(2)
            celllb(2) = celllb(2) + kellogb(2, k) * rmu(2, k,1) / cellscale(2)

            cellb(3) = cellb(3) + kelelocb(3, k) * rmu(3, k,1) / cellscale(3)
            celllb(3) = celllb(3) + kellogb(3, k) * rmu(3, k,1) / cellscale(3)
        enddo
    endif

    return
end
