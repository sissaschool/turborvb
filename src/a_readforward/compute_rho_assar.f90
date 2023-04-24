! Copyright (C) 2022 TurboRVB group
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.

module Assar_module

    logical :: ifrho_assar !(Matteo) flag for Assaraf density

    real*8, allocatable :: am(:, :), bm(:), Gx(:)
    real*8 :: assar_cut, assar_parr, assar_parr2, assar_parr3, Gn
    real :: dxas, dyas, dzas
    integer :: nx, ny, nz, coefx, coefy, coefz, bins_cut(3)
    real*8 :: kswitch
    integer :: nswitch
    real*8, allocatable :: dist_grid_ion(:, :)

contains

    subroutine initialize_assaraf(vell, ngrid_l)
        use grid_module, only: out_grid, grid_points, ext_grid
        use allio, only: nion, rion
        implicit none

        real*8 :: vell(3)
        integer :: ngrid_l(3)

        integer i, k

        allocate (dist_grid_ion(grid_points, nion), am(grid_points, nion), bm(grid_points), Gx(grid_points))

        assar_parr2 = assar_parr**2
        assar_parr3 = assar_parr**3

        !bins_cut=int(assar_cut/abs(out_grid(1,1)-out_grid(1,2)))
        bm(:) = 1

        ! am(:,:)=1
        ! Gx(:)=1

        do i = 1, grid_points
            Gx(i) = 0.d0
            do k = 1, nion
                dist_grid_ion(i, k) = sqrt((rion(1, k) - out_grid(1, i))**2 &
                                           + (rion(2, k) - out_grid(2, i))**2 + (rion(3, k) - out_grid(3, i))**2)
                am(i, k) = kswitch**nswitch/(kswitch**nswitch + dist_grid_ion(i, k)**nswitch)
                bm(i) = bm(i)*(dist_grid_ion(i, k)**nswitch)/(kswitch**nswitch + dist_grid_ion(i, k)**nswitch)
                if (dist_grid_ion(i, k) .le. kswitch) then
                    goto 117
                end if
                Gx(i) = Gx(i) + 1./dist_grid_ion(i, k)
            end do
            Gx(i) = Gx(i)/nion
            goto 118
117         Gx(i) = 0.d0
118         dxas = 0
        end do
        if (.not. ext_grid) then
            dxas = vell(1)
            dyas = vell(2)
            dzas = vell(3)
            !dxas=abs(out_grid(1,1)-out_grid(1,2))
            !dyas=abs(out_grid(2,1)-out_grid(2,2))
            !dzas=abs(out_grid(3,1)-out_grid(3,2))
            !dyas=dxas
            !dzas=dxas
            !nx=int(grid_points**(1./3))
            !ny=nx
            !nz=nx
            nx = ngrid_l(1)
            ny = ngrid_l(2)
            nz = ngrid_l(3)
            coefx = int(abs(out_grid(1, 1))/dxas)
            coefy = int(abs(out_grid(2, 1))/dyas)
            coefz = int(abs(out_grid(3, 1))/dzas)
            !coefy=coefx
            !coefz=coefx
            bins_cut(1) = int(assar_cut/dxas)
            bins_cut(2) = int(assar_cut/dyas)
            bins_cut(3) = int(assar_cut/dzas)
        end if

    end subroutine initialize_assaraf

    subroutine compute_rho_assar(dist, rion, nion, kel, nelup, neldo, indt, tabpip, winvup, winvdo, density, zetaq)
        use grid_module
        use constants, only: pi
        implicit none
        integer i, j, k, jp, indt, nelup, neldo, nion !nelup, neldo - number spinup spindown, nion - # of atoms
        !Coordinates and distances------------------------------------------------------------------------------------
        real*8 rion(3, nion) !coord'tes of ions, x,y,z and nion - # of atoms
        real*8 kel(3, *) !all the electron cooord computed in VMC run
        real*8 dist(nelup + neldo, nion) !distance between an ion and an electron (VMC)
        real*8 :: comp_grid_ele(3), dist_grid_ele !coordinates and distance between a grid point and an electron (VMC)
        !Functions-----------------------------------------------------------------------------------------------------
        real*8 density(*) !vector of density
        real*8 :: Fx
        real*8 :: gradF(3), nablaF !Operators on F grad laplas
        real*8 :: wave_laplacian, wave_gradient(3)
        real*8 :: winvup(nelup, *), winvdo(neldo, *), tabpip(nelup + neldo, *) !* - as long as num of electrons
        !Parameters----------------------------------------------------------------------------------------------------
        real*8 :: zetaq(nion) !, assar_cut                     !Zq-the atomic number,assar_cut - cutoff
        real*8 :: exp_grid !exp(-lambda|ri-r|), distance |ri-r| cut off
        integer ix, iy, iz, ax, ay, az, bx, by, bz, jx, jy, jz
        !--------------------------------------------------------------------------------------------------------------
        do i = 1, grid_points
            density(i) = 0.d0
        end do

        if (ext_grid) then
            !Per gli elettroni UP
            do j = 1, nelup
                do i = 1, grid_points
                    comp_grid_ele(:) = -out_grid(:, i) + kel(:, j)
                    dist_grid_ele = sqrt(sum(comp_grid_ele(:)**2))
                    exp_grid = exp(-assar_parr*dist_grid_ele)

                    !Calculating f, grad(f) and nabla(f)------------------------------------------
                    Fx = 1
                    gradF(:) = 0.d0
                    nablaF = 0.d0
                    do k = 1, nion
                        Fx = Fx + am(i, k)*2*zetaq(k)*(dist(j, k) - dist_grid_ion(i, k))
                        gradF(:) = gradF(:) + am(i, k)*2*zetaq(k)*(kel(:, j) - rion(:, k))/dist(j, k)
                        nablaF = nablaF + am(i, k)*4*zetaq(k)/dist(j, k)
                    end do
                    Fx = Fx - bm(i) + bm(i)*(1 + assar_parr*dist_grid_ele)*exp_grid
                    gradF(:) = gradF(:) - bm(i)*assar_parr2*comp_grid_ele(:)*exp_grid
                    nablaF = nablaF + bm(i)*(assar_parr3*dist_grid_ele - 3*assar_parr2)*exp_grid
                    !------------------------------------------------------------------------------
                    wave_gradient(:) = (winvup(j, indt + 1:indt + 3) + tabpip(j, indt + 1:indt + 3))
                    wave_laplacian = winvup(j, indt + 4) + tabpip(j, indt + 4) &
                                     + 2*sum(winvup(j, indt + 1:indt + 3)*tabpip(j, indt + 1:indt + 3)) &
                                     + sum(tabpip(j, indt + 1:indt + 3)*tabpip(j, indt + 1:indt + 3))
                    density(i) = density(i) &
                                 - 1/(4*pi)*(1/dist_grid_ele - bm(i)*Gx(i)) &
                                   &*(nablaF + 2*Fx*(wave_laplacian + sum(wave_gradient(:)**2)) &
                                   &+ 4*sum(gradF(:)*wave_gradient(:)))
                end do
            end do
            !Per gli elettroni DOWN
            do j = 1, neldo
                jp = j + nelup
                do i = 1, grid_points
                    comp_grid_ele(:) = -out_grid(:, i) + kel(:, jp)
                    dist_grid_ele = sqrt(sum(comp_grid_ele(:)**2))
                    exp_grid = exp(-assar_parr*dist_grid_ele)
                    !Calculating f, grad(f) and nabla(f)------------------------------------------
                    Fx = 1
                    gradF(:) = 0.d0
                    nablaF = 0.d0
                    do k = 1, nion
                        Fx = Fx + am(i, k)*2*zetaq(k)*(dist(jp, k) - dist_grid_ion(i, k))
                        gradF(:) = gradF(:) + am(i, k)*2*zetaq(k)*(kel(:, jp) - rion(:, k))/dist(jp, k)
                        nablaF = nablaF + am(i, k)*4*zetaq(k)/dist(jp, k)
                    end do
                    Fx = Fx - bm(i) + bm(i)*(1 + assar_parr*dist_grid_ele)*exp_grid
                    gradF(:) = gradF(:) - bm(i)*assar_parr2*comp_grid_ele(:)*exp_grid
                    nablaF = nablaF + bm(i)*(assar_parr3*dist_grid_ele - 3*assar_parr2)*exp_grid
                    !------------------------------------------------------------------------------
                    wave_gradient(:) = winvdo(j, indt + 1:indt + 3) + tabpip(jp, indt + 1:indt + 3)
                    wave_laplacian = winvdo(j, indt + 4) + tabpip(jp, indt + 4) &
                                    & + 2*sum(winvdo(j, indt + 1:indt + 3)*tabpip(jp, indt + 1:indt + 3))&
                                    & + sum(tabpip(jp, indt + 1:indt + 3)*tabpip(jp, indt + 1:indt + 3))
                    density(i) = density(i) &
                                 - 1/(4*pi)*(1/dist_grid_ele - bm(i)*Gx(i)) &
                                 *(nablaF + 2*Fx*(wave_laplacian + sum(wave_gradient(:)**2)) &
                                   + 4*sum(gradF(:)*wave_gradient(:)))
                end do
            end do
        else
            do j = 1, nelup
                !------------------------------------------------------
                ix = int(kel(1, j)/dxas) + coefx
                ax = max(1, ix - bins_cut(1))
                bx = min(nx, ix + bins_cut(1))
                !-----------------------------------------------------
                iy = int(kel(2, j)/dyas) + coefy
                ay = max(1, iy - bins_cut(2))
                by = min(ny, iy + bins_cut(2))
                !-----------------------------------------------------
                iz = int(kel(3, j)/dzas) + coefz
                az = max(1, iz - bins_cut(3))
                bz = min(nz, iz + bins_cut(3))
                !-----------------------------------------------------
                do jz = az, bz
                    do jy = ay, by
                        do jx = ax, bx
                            i = jx + (jy - 1)*nx + (jz - 1)*nx*ny
                            comp_grid_ele(:) = -out_grid(:, i) + kel(:, j)
                            dist_grid_ele = sqrt(sum(comp_grid_ele(:)**2))
                            exp_grid = exp(-assar_parr*dist_grid_ele)
                            !Calculating f, grad(f) and nabla(f)------------------------------------------
                            Fx = 1
                            gradF(:) = 0.d0
                            nablaF = 0.d0
                            do k = 1, nion
                                Fx = Fx + am(i, k)*2*zetaq(k)*(dist(j, k) - dist_grid_ion(i, k))
                                gradF(:) = gradF(:) + am(i, k)*2*zetaq(k)*(kel(:, j) - rion(:, k))/dist(j, k)
                                nablaF = nablaF + am(i, k)*4*zetaq(k)/dist(j, k)
                            end do
                            Fx = Fx - bm(i) + bm(i)*(1 + assar_parr*dist_grid_ele)*exp_grid
                            gradF(:) = gradF(:) - bm(i)*assar_parr2*comp_grid_ele(:)*exp_grid
                            nablaF = nablaF + bm(i)*(assar_parr3*dist_grid_ele - 3*assar_parr2)*exp_grid
                            !------------------------------------------------------------------------------
                            wave_gradient(:) = (winvup(j, indt + 1:indt + 3) + tabpip(j, indt + 1:indt + 3))
                            wave_laplacian = winvup(j, indt + 4) + tabpip(j, indt + 4) &
                                             + 2*sum(winvup(j, indt + 1:indt + 3)*tabpip(j, indt + 1:indt + 3)) &
                                             + sum(tabpip(j, indt + 1:indt + 3)*tabpip(j, indt + 1:indt + 3))
                            density(i) = density(i) &
                                         - 1/(4*pi)*(1/dist_grid_ele - bm(i)*Gx(i)) &
                                         *(nablaF + 2*Fx*(wave_laplacian + &
                                                          sum(wave_gradient(:)**2)) + 4*sum(gradF(:)*wave_gradient(:)))
                        end do
                    end do
                end do
            end do
            do j = 1, neldo
                jp = j + nelup
                !------------------------------------------------------
                ix = int(kel(1, jp)/dxas) + coefx
                ax = max(1, ix - bins_cut(1))
                bx = min(nx, ix + bins_cut(1))
                !-----------------------------------------------------
                iy = int(kel(2, jp)/dyas) + coefy
                ay = max(1, iy - bins_cut(2))
                by = min(ny, iy + bins_cut(2))
                !------------------------------------------------------
                iz = int(kel(3, jp)/dzas) + coefz
                az = max(1, iz - bins_cut(3))
                bz = min(nz, iz + bins_cut(3))
                !-----------------------------------------------------
                do jz = az, bz
                    do jy = ay, by
                        do jx = ax, bx
                            i = jx + (jy - 1)*nx + (jz - 1)*nx*ny

                            comp_grid_ele(:) = -out_grid(:, i) + kel(:, jp)
                            dist_grid_ele = sqrt(sum(comp_grid_ele(:)**2))
                            exp_grid = exp(-assar_parr*dist_grid_ele)
                            !Calculating f, grad(f) and nabla(f)------------------------------------------
                            Fx = 1
                            gradF(:) = 0.d0
                            nablaF = 0.d0
                            do k = 1, nion
                                Fx = Fx + am(i, k)*2*zetaq(k)*(dist(jp, k) - dist_grid_ion(i, k))
                                gradF(:) = gradF(:) + am(i, k)*2*zetaq(k)*(kel(:, jp) - rion(:, k))/dist(jp, k)
                                nablaF = nablaF + am(i, k)*4*zetaq(k)/dist(jp, k)
                            end do
                            Fx = Fx - bm(i) + bm(i)*(1 + assar_parr*dist_grid_ele)*exp_grid
                            gradF(:) = gradF(:) - bm(i)*assar_parr2*comp_grid_ele(:)*exp_grid
                            nablaF = nablaF + bm(i)*(assar_parr3*dist_grid_ele - 3*assar_parr2)*exp_grid
                            !------------------------------------------------------------------------------
                            wave_gradient(:) = winvdo(j, indt + 1:indt + 3) + tabpip(jp, indt + 1:indt + 3)
                            wave_laplacian = winvdo(j, indt + 4) + tabpip(jp, indt + 4) &
                                             + 2*sum(winvdo(j, indt + 1:indt + 3)*tabpip(jp, indt + 1:indt + 3)) &
                                             + sum(tabpip(jp, indt + 1:indt + 3)*tabpip(jp, indt + 1:indt + 3))
                            density(i) = density(i) - 1/(4*pi)*(1/dist_grid_ele - bm(i)*Gx(i))*(nablaF + 2*Fx*(wave_laplacian&
                                    & + sum(wave_gradient(:)**2)) + 4*sum(gradF(:)*wave_gradient(:)))
                        end do
                    end do
                end do
            end do
        end if

    end subroutine compute_rho_assar

end module Assar_module
