module logger_io
  use stdlib_logger
  use m_tostring, only: tostring
  use, intrinsic :: iso_fortran_env, only : error_unit, input_unit, output_unit
  implicit none

  private
  character(len=4096) :: log_msg_  ! buffer
  logical :: export_ = .true.
  integer,parameter :: info_level = information_level

  public :: log_error, log_warning, log_info, log_debug, logger_config
  public :: all_level, debug_level, info_level, warning_level, error_level, none_level

contains
  subroutine log_error(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    class(*),intent(in),optional :: arg1, arg2, arg3, arg4, arg5
    class(*),intent(in),optional :: arg6, arg7, arg8, arg9
    if (export_) then
       log_msg_ = tostring(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
       call global_logger%log_error(message=log_msg_)
    end if
  end subroutine log_error

  subroutine log_warning(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    class(*),intent(in),optional :: arg1, arg2, arg3, arg4, arg5
    class(*),intent(in),optional :: arg6, arg7, arg8, arg9
    if (export_) then
       log_msg_ = tostring(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
       call global_logger%log_warning(message=log_msg_)
    end if
  end subroutine log_warning

  subroutine log_info(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    class(*),intent(in),optional :: arg1, arg2, arg3, arg4, arg5
    class(*),intent(in),optional :: arg6, arg7, arg8, arg9
    if (export_) then
       log_msg_ = tostring(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
       call global_logger%log_information(message=log_msg_)
    end if
  end subroutine log_info

  subroutine log_debug(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    class(*),intent(in),optional :: arg1, arg2, arg3, arg4, arg5
    class(*),intent(in),optional :: arg6, arg7, arg8, arg9
    if (export_) then
       log_msg_ = tostring(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
       call global_logger%log_debug(message=log_msg_)
    end if
  end subroutine log_debug

  subroutine logger_config(level, rank, force_output, filename)
    integer,intent(in),optional :: level
    integer,intent(in),optional :: rank
    logical,intent(in),optional :: force_output
    character(*),intent(in),optional :: filename

    integer :: rank_info = 0
    character(len=256) :: filename_buf
    integer :: stat
    
    if (present(level)) then
       call global_logger%configure(level=level)
    end if
    if (present(rank)) then
       rank_info = rank
       export_ = rank .eq. 0
    end if
    if (present(force_output)) then
       export_ = force_output
    end if
    if (present(filename).and.export_) then
       write(filename_buf, '(A,i0)') trim(filename) // '.', rank_info
       call global_logger%add_log_file(trim(filename_buf), position='APPEND', stat=stat)
       if (stat .ne. success) then
          error stop 'unable to open output file.'
       end if
       if (rank_info .eq. 0) then
          call global_logger%add_log_unit(output_unit)
       end if
    end if
  end subroutine logger_config
  
end module logger_io
