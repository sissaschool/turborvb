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

subroutine pseudoset(jel, kel, ivic, prefactor, pseudolocal &
        &, nel, dist, rion, wpseudo, ncore, skip &
        &, indtm, iopt, angle, psip, Lbox, mindist, iflagpseudo, vec)
    use allio, only: rcutoff, versor, nintpseudo, kindion, pshell, nparpshell &
            &, parshell, lmax, wintpseudo, jpseudo, enforce_detailb, norm_metric
    use Cell, only: cellscale, s2r, car2cry, Lmin, ApplyPBC, yes_tilted, metric
    implicit none
    integer jel, core, ncore, lang, jj, kk, nel, skip, indtm, npseudo &
            &, iflagpseudo, kimages, maximages
    real*8 ivic(3, *), prefactor(skip, *), dist(*), rion(3, *), distmin &
            &, wpseudo(*), pseudofun, pseudolocal, scal, legfun &
            &, rver(3), r(3), rsav(3), angle(3, 3), psip(*), Lbox, mindist, vec(3)

    !      input kel,rion,dist
    !      output prefactor (passed ) ,legendre,wpseudo, ivic (passed)
    !      pseudolocal (passed)
    logical iopt
    real*8 kel(3)

    if (iflagpseudo .ne. 0 .or. nintpseudo .le. 0) return

    ! if iopt.eq.true generate a new rotation matrix
    if (iopt) then
        if (enforce_detailb) then
            call fillmatrix_vec(angle, vec)
        else
            call fillmatrix(angle)
        end if
    end if

    do jj = 1, skip
        do kk = 1, 3
            ivic(kk, jj) = 0.d0
        end do
        prefactor(jj, jel) = 0.d0
    end do
    pseudolocal = 0.d0

    npseudo = 0
    do core = 1, ncore
        distmin = dist(kindion(core))

        if (distmin .le. rcutoff(core)) then

            if (distmin .le. mindist) distmin = mindist

            do kk = 1, 3
                rsav(kk) = kel(kk) - rion(kk, kindion(core))
            end do
            if (Lbox .gt. 0.d0) then
                !   Find the image ion rion(:,kindion(core))-rdiff closest to kel
                !   NB I am not using PBC on electron positions
                !   here so it works also for APBC or whatever boundary condition.
                !          call ApplyPBC(rsav,1)  ! the  output will be rsav-rion(:,..)+rdiff
                rsav(:) = car2cry(:, 1)*rsav(1) + car2cry(:, 2)*rsav(2) + car2cry(:, 3)*rsav(3)
                rsav(1:3) = rsav(1:3) - cellscale(1:3)*anint(rsav(1:3)/cellscale(1:3))
                if (rcutoff(core) .gt. Lmin/2.d0) then
                    maximages = 4
                else
                    maximages = 1
                end if
            else
                maximages = 1
            end if

            do kimages = 1, maximages
                !           redefine all for the images
                r = rsav
                !           if(kimages.eq.1) then
                !           d=dsav
                if (kimages .eq. 2) then ! below PBC is certain
                    !              first image over x
                    if (r(1) .gt. 0.d0) then
                        r(1) = rsav(1) - cellscale(1)
                    else
                        r(1) = rsav(1) + cellscale(1)
                    end if
                    distmin = norm_metric(r, metric)
                elseif (kimages .eq. 3) then
                    !              second image over y
                    if (r(2) .gt. 0.d0) then
                        r(2) = rsav(2) - cellscale(2)
                    else
                        r(2) = rsav(2) + cellscale(2)
                    end if
                    distmin = norm_metric(r, metric)
                elseif (kimages .eq. 4) then
                    !              third image over z
                    if (r(3) .gt. 0.d0) then
                        r(3) = rsav(3) - cellscale(3)
                    else
                        r(3) = rsav(3) + cellscale(3)
                    end if
                    distmin = norm_metric(r, metric)
                end if

                if (LBox .gt. 0.d0) then
                    !          go back in Cartesian Coordinates
                    r(1:3) = r(1:3)/cellscale(1:3)
                    r(:) = s2r(:, 1)*r(1) + s2r(:, 2)*r(2) + s2r(:, 3)*r(3)
                end if

                if (distmin .le. rcutoff(core)) then

                    indtm = indtm + nintpseudo

                    npseudo = npseudo + 1

                    if (npseudo .gt. skip/nintpseudo) then
                        write (6, *) ' ERROR nuclei too close for pseudo !!! '
                        write (6, *) 'increase max # of nuclei, npsamax>', skip/nintpseudo
                        iflagpseudo = 1
                        return
                    end if

                    do lang = 1, pshell(core) - 1
                        wpseudo(lang) = pseudofun(nparpshell(lang, core)         &
                                &, distmin, parshell(1, jpseudo(lang, core)), psip)
                    end do
                    pseudolocal = pseudolocal                                     &
                            & + pseudofun(nparpshell(pshell(core), core)              &
                                    &, distmin, parshell(1, jpseudo(pshell(core), core)), psip)

                    do jj = (npseudo - 1)*nintpseudo + 1, npseudo*nintpseudo

                        !              call rotation(angle
                        !    &         ,versor(1,jj-(npseudo-1)*nintpseudo),rver)
                        call dgemv('N', 3, 3, 1.d0, angle, 3                          &
                                &, versor(1, jj - (npseudo - 1)*nintpseudo), 1, 0.d0, rver, 1)

                        scal = 0.d0
                        do kk = 1, 3
                            scal = scal + rver(kk)*r(kk)
                            !              This has to be referred to the chosen Image.
                            !              rver(kk)=distmin*rver(kk)+rion(kk,kindion(core))-rdiff(kk)
                            !              ivic(kk,jj)=rver(kk)-kel(kk)
                            !             N.B.  rion-rdiff-kel=-r  by the above definition

                            !              ivic is the versor to be added to kel to have the new
                            !              mesh electronic point in the pseudo integration.
                            ivic(kk, jj) = rver(kk)*distmin - r(kk)
                        end do

                        scal = scal/distmin

                        do lang = 1, pshell(core) - 1
                            prefactor(jj, jel) = prefactor(jj, jel) + &
                                    &            (2.d0*lang - 1.d0)*wpseudo(lang)*legfun(lang - 1, scal)
                        end do
                        prefactor(jj, jel) = prefactor(jj, jel)                      &
                                            & *2.d0*wintpseudo(jj - (npseudo - 1)*nintpseudo)

                    end do

                end if ! endif check distance images

            end do ! enddo kimages

        end if ! endif check distance

    end do ! enddo core

    return
end subroutine pseudoset
