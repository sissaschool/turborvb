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

subroutine ratio_psi(Rold, Rnew, parcut, ratioreg, reweight, epstl)
    implicit none
    real*8 Rold, Rnew, parcut, ratioreg, reweight, epstl
    real*8 Roldreg, Rnewreg

    !      This function computes serveral values related to regularization
    !      which is defined as \psi_g(x)=R_\epsilon(x)*[\psi_T(x)/R(x)].
    !      R_\epsilon(x) is explicitly expressed here.

    !      Input  : Rold,Rnew,parcut,epstl
    !               Rold and Rnew : old and new (a proposed move) value
    !                 of R(x) function which is computed in ratiovar()
    !                 as ratiodetr
    !               parcut : regularization parameter
    !               epstl : precision control parameter (if zero not used)

    !      Output : ratioreg,reweight
    !               ratioreg : new R_\epsilon(x) / old R_\epsilon(x)
    !               reweight : [new \psi_T(x) / new \psi_g(x)]^2

    if (Rnew .ne. 0.d0) then
        if (Rnew .gt. epstl .or. epstl .eq. 0.d0) then
            !   OLD method  C.  Attaccalite & S. Sorella 2008
            !         if(Rnew.le.parcut.and.Rold.gt.parcut) then
            !           Rnewreg=parcut*(Rnew/parcut)**(Rnew/parcut)
            !           reweight=(Rnew/Rnewreg)**2
            !           ratioreg=Rnewreg/Rold
            !         elseif(Rnew.gt.parcut.and.Rold.le.parcut) then
            !           Roldreg=parcut*(Rold/parcut)**(Rold/parcut)
            !           reweight=1.d0
            !           ratioreg=Rnew/Roldreg
            !         elseif(Rnew.le.parcut.and.Rold.le.parcut) then
            !           Rnewreg=parcut*(Rnew/parcut)**(Rnew/parcut)
            !           Roldreg=parcut*(Rold/parcut)**(Rold/parcut)
            !           reweight=(Rnew/Rnewreg)**2
            !           ratioreg=Rnewreg/Roldreg
            !         elseif(Rnew.gt.parcut.and.Rold.gt.parcut) then
            !           ratioreg=Rnew/Rold
            !           reweight=1.d0
            !         endif
            ! HERE R_epsilon = Max( \epsilon, R(x))
            if (Rnew .le. parcut .and. Rold .gt. parcut) then
                Rnewreg = parcut
                reweight = (Rnew/Rnewreg)**2
                ratioreg = Rnewreg/Rold
            elseif (Rnew .gt. parcut .and. Rold .le. parcut) then
                Roldreg = parcut
                reweight = 1.d0
                ratioreg = Rnew/Roldreg
            elseif (Rnew .le. parcut .and. Rold .le. parcut) then
                Rnewreg = parcut
                Roldreg = parcut
                reweight = (Rnew/Rnewreg)**2
                ratioreg = 1.d0
            elseif (Rnew .gt. parcut .and. Rold .gt. parcut) then
                ratioreg = Rnew/Rold
                reweight = 1.d0
            end if
        else
            ratioreg = 0.d0
            reweight = 0.d0
        end if
    else
        ratioreg = 0.d0
        reweight = 0.d0
    end if

    return
end
