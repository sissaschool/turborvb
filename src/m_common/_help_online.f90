subroutine help_online(name_tool)
    implicit none
    integer ncount
    character(100) :: name_dir, name_tool, name_ext
    character(200) :: name_file
    character(256) :: linedata

    name_dir = 'NAME_DIR'
    !commented out by K.Nakano on 13th Sep.
    !Readmes will be integrated into the TurboRVB manual.
    !name_ext='.README'
    !name_file=TRIM(name_dir)//'/doc/readme/'//TRIM(name_tool)//TRIM(name_ext)
    !
    !write(*,*)name_dir
    !write(*,*)name_tool
    !write(*,*)name_file
    !write(*,*)name_ext
    !
    !open(unit=22,file=name_file,status='unknown',form='formatted')
    !
    !
    !ncount=0
    !do while(.true.)
    ! read(22,'(a)', end=100) linedata
    ! ncount=ncount+1
    ! write(*,*)TRIM(linedata)
    ! if(mod(ncount,22).eq.0) then
    ! read(5,*)
    ! endif
    !end do
    !100   continue
    !if(ncount.eq.0) then
    !write(6,*) ' Sorry no on-line documentation exists for this tool '
    !close(22,status='DELETE') ! delete this empty file to avoid confusion
    !endif
    !
    !if(ncount.ne.0) close(22)

    name_ext = '.input'
    name_file = trim(name_tool)//trim(name_ext)
    open (unit=23, file=name_file, status='unknown', form='formatted')
!     check if the output is empty
    ncount = 0
    do while (.true.)
        read (23, '(a)', end=110) linedata
        ncount = ncount + 1
    end do
110 continue
    if (ncount .eq. 0) then ! only if the output file is empty
!     for ibm
        close (23)
        open (unit=23, file=name_file, status='unknown', form='formatted')
        name_file = trim(name_dir)//'/template/'//trim(name_tool)//trim(name_ext)
        !write(6,*) name_file
        open (unit=22, file=name_file, status='unknown', form='formatted')
        ncount = 0
        do while (.true.)
            read (22, '(a)', end=120) linedata
            ncount = ncount + 1
            write (23, '(a)') trim(linedata)
        end do
120     continue

        if (ncount .gt. 0) then
            write (6, *)
            write (6, *) ' Warning a sample input is given  ', trim(name_tool)//trim(name_ext)
            close (22)
            close (23)
        else
            close (22, status='DELETE')
            close (23, status='DELETE')
        end if
    else
        close (23) ! close safely the existing file
    end if
    return
end

