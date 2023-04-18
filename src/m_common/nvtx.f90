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

module nvtx

    use iso_c_binding
    implicit none

#ifdef _NVTX
    integer, private :: col(7) = [Z'0000ff00', Z'000000ff', Z'00ffff00', Z'00ff00ff', Z'0000ffff', Z'00ff0000', Z'00ffffff']
    character, private, target :: tempName(256)

    type, bind(C) :: nvtxEventAttributes
        integer(c_int16_t) :: version = 1
        integer(c_int16_t) :: size = 48 !
        integer(c_int) :: category = 0
        integer(c_int) :: colorType = 1 ! NVTX_COLOR_ARGB = 1
        integer(c_int) :: color
        integer(c_int) :: payloadType = 0 ! NVTX_PAYLOAD_UNKNOWN = 0
        integer(c_int) :: reserved0
        integer(c_int64_t) :: payload ! union uint,int,double
        integer(c_int) :: messageType = 1 ! NVTX_MESSAGE_TYPE_ASCII     = 1
        type(c_ptr) :: message ! ascii char
    end type

    interface nvtxRangePush
        ! push range with custom label and standard color
        subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
            use iso_c_binding
            character(kind=c_char) :: name(256)
        end subroutine

        ! push range with custom label and custom color
        subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
            use iso_c_binding
            import :: nvtxEventAttributes
            type(nvtxEventAttributes) :: event
        end subroutine
    end interface

    interface nvtxRangePop
        subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
        end subroutine
    end interface

contains

    subroutine nvtxStartRange(name, id)
        character(kind=c_char, len=*) :: name
        integer, optional :: id
        type(nvtxEventAttributes) :: event
        character(kind=c_char, len=256) :: trimmed_name
        integer :: i

        trimmed_name = trim(name)//c_null_char

        ! move scalar trimmed_name into character array tempName
        do i = 1, len(trim(name)) + 1
            tempName(i) = trimmed_name(i:i)
        end do

        if (.not. present(id)) then
            call nvtxRangePush(tempName)
        else
            event%color = col(mod(id, 7) + 1)
            event%message = c_loc(tempName)
            call nvtxRangePushEx(event)
        end if
    end subroutine

    subroutine nvtxEndRange
        call nvtxRangePop
    end subroutine

#endif

end module nvtx
