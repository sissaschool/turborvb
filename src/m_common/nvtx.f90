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

/**
 * @brief NVTX profiling module for performance analysis
 *
 * This module provides an interface to NVIDIA Tools Extension (NVTX) for
 * performance profiling and analysis. It allows marking code regions and
 * events for visualization in NVIDIA profiling tools like Nsight Systems
 * and Nsight Compute.
 *
 * @details
 * The module provides the following functionality:
 * - Range markers for profiling code sections
 * - Custom colors for visual distinction
 * - Event attributes for detailed profiling
 * - C-compatible string handling for NVTX interface
 * - Conditional compilation for NVTX support
 *
 * Key features:
 * - Automatic color cycling for different regions
 * - Support for custom event attributes
 * - Integration with NVIDIA profiling tools
 * - No performance impact when NVTX is disabled
 * - Thread-safe profiling markers
 *
 * @note
 * - Only available when _NVTX is defined
 * - Requires NVIDIA GPU and compatible drivers
 * - Used for performance analysis and optimization
 * - No functionality when compiled without NVTX support
 *
 * @see nvtxStartRange(), nvtxEndRange()
 *
 * @author TurboRVB group
 * @date 2022
 */
module nvtx

    use iso_c_binding
    implicit none

#ifdef _NVTX
    /**
     * @brief Color array for NVTX range markers
     *
     * Array of predefined colors for NVTX range markers, providing
     * visual distinction between different code regions in profiling tools.
     * Colors are specified in ARGB format.
     */
    integer, private :: col(7) = [Z'0000ff00', Z'000000ff', Z'00ffff00', Z'00ff00ff', Z'0000ffff', Z'00ff0000', Z'00ffffff']
    
    /**
     * @brief Temporary character array for NVTX string handling
     *
     * Character array used for converting Fortran strings to C-style
     * null-terminated strings for NVTX interface calls.
     */
    character, private, target :: tempName(256)

    /**
     * @brief NVTX event attributes structure
     *
     * C-compatible structure for NVTX event attributes, used to specify
     * custom properties for range markers including color, message, and
     * payload information.
     */
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

    /**
     * @brief Interface for NVTX range push operations
     *
     * Interface providing two methods for pushing NVTX range markers:
     * - nvtxRangePushA: Simple range with custom label and standard color
     * - nvtxRangePushEx: Range with custom label and custom color/attributes
     */
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

    /**
     * @brief Interface for NVTX range pop operations
     *
     * Interface for popping NVTX range markers, ending the current
     * profiling range.
     */
    interface nvtxRangePop
        subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
        end subroutine
    end interface

contains

    /**
     * @brief Start an NVTX profiling range
     *
     * This subroutine starts an NVTX profiling range with the specified name.
     * It can optionally use a custom color based on an ID parameter.
     *
     * @param[in] name Name of the profiling range (will be displayed in profiling tools)
     * @param[in] id Optional ID for color selection (if not provided, uses standard color)
     *
     * @details
     * The subroutine:
     * 1. Converts the Fortran string to C-style null-terminated string
     * 2. If ID is provided, creates custom event attributes with color
     * 3. Calls appropriate NVTX function to start the range
     * 4. Colors are cycled through a predefined array based on ID
     *
     * @note
     * - Only available when _NVTX is defined
     * - Colors cycle through 7 predefined colors
     * - String is automatically null-terminated
     * - No effect when NVTX is not available
     *
     * @see nvtxEndRange()
     */
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

    /**
     * @brief End the current NVTX profiling range
     *
     * This subroutine ends the current NVTX profiling range, marking
     * the end of a code section for profiling analysis.
     *
     * @details
     * The subroutine calls nvtxRangePop to end the current range
     * and return to the previous profiling context.
     *
     * @note
     * - Only available when _NVTX is defined
     * - Must be called to match each nvtxStartRange call
     * - No effect when NVTX is not available
     * - No parameters required
     *
     * @see nvtxStartRange()
     */
    subroutine nvtxEndRange
        call nvtxRangePop
    end subroutine

#endif

end module nvtx
