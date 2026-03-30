! https://fortranwiki.org/fortran/show/tostring
module M_tostring
  implicit none
  private
  public tostring
contains
  function tostring(generic1, generic2, generic3, generic4, generic5, generic6, generic7, generic8, generic9)
    implicit none
    ! convert up to nine scalar intrinsic values to a string
    class(*),intent(in),optional  :: generic1 ,generic2 ,generic3 ,generic4, generic5
    class(*),intent(in),optional  :: generic6 ,generic7 ,generic8 ,generic9
    character(len=:), allocatable :: tostring
    character(len=4096)        :: line
    integer                    :: istart
    line=''
    istart=1
    if(present(generic1))call print_generic(generic1)
    if(present(generic2))call print_generic(generic2)
    if(present(generic3))call print_generic(generic3)
    if(present(generic4))call print_generic(generic4)
    if(present(generic5))call print_generic(generic5)
    if(present(generic6))call print_generic(generic6)
    if(present(generic7))call print_generic(generic7)
    if(present(generic8))call print_generic(generic8)
    if(present(generic9))call print_generic(generic9)
    tostring=trim(line)
  contains
    subroutine print_generic(generic)
      use,intrinsic :: iso_fortran_env, only : int8, int16, int32, int64, real32, real64, real128
      class(*),intent(in),optional :: generic
      select type(generic)
      type is (integer(kind=int8));     write(line(istart:),'(i0)') generic
      type is (integer(kind=int16));    write(line(istart:),'(i0)') generic
      type is (integer(kind=int32));    write(line(istart:),'(i0)') generic
      type is (integer(kind=int64));    write(line(istart:),'(i0)') generic
      type is (real(kind=real32));      write(line(istart:),'(1pg0)') generic
      type is (real(kind=real64));      write(line(istart:),'(1pg0)') generic
      type is (real(kind=real128));     write(line(istart:),'(1pg0)') generic
      type is (logical);                write(line(istart:),'(1l)') generic
      type is (character(len=*));       write(line(istart:),'(a)') generic
      type is (complex);                write(line(istart:),'("(",1pg0,",",1pg0,")")') generic
      end select
      istart=len_trim(line)+2
    end subroutine print_generic
  end function tostring
end module M_tostring
