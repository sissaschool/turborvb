!> @file _help_online.f90
!> @brief Subroutine for providing online help and sample input files in TurboRVB
!> @author TurboRVB group
!> @date 2022
!> @version 1.0
!> @details
!> This subroutine provides online help functionality by checking for existing
!> input files and generating sample input files from templates when needed.
!> It was originally designed to display README files but has been modified
!> to focus on input file management and template generation.
!>
!> The subroutine performs the following operations:
!> - Checks if an input file exists for the specified tool
!> - If the input file is empty, copies a template from the template directory
!> - Provides user feedback about the generated sample input
!> - Handles file operations safely with proper cleanup
!>
!> @note README functionality has been commented out and integrated into TurboRVB manual
!> @warning File operations use hardcoded unit numbers (22, 23)

!> @brief Provide online help and sample input files for TurboRVB tools
!> @details This subroutine implements online help functionality by managing
!> input files and templates. It first checks if an input file exists for the
!> specified tool. If the file is empty or doesn't exist, it copies a template
!> from the template directory to provide a sample input file for the user.
!>
!> The algorithm works as follows:
!> 1. Check if input file exists and has content
!> 2. If empty, look for template in template directory
!> 3. Copy template to create sample input file
!> 4. Provide user feedback about the generated file
!> 5. Clean up file handles and temporary files
!>
!> Originally, this subroutine also displayed README files, but that functionality
!> has been commented out as READMEs are now integrated into the TurboRVB manual.
!>
!> @param[in] name_tool Name of the TurboRVB tool for which help is requested
!> @note The tool name should match the corresponding template file name
!> @warning File operations use hardcoded unit numbers that may conflict with other routines
!> @see Template files in template/ directory for available sample inputs
subroutine help_online(name_tool)
    implicit none
    
    !> @var ncount Counter for lines read from files
    integer ncount
    
    !> @var name_dir Directory name for templates (hardcoded as 'NAME_DIR')
    character(100) :: name_dir
    
    !> @var name_tool Name of the tool for which help is requested
    character(100) :: name_tool
    
    !> @var name_ext File extension for input files ('.input')
    character(100) :: name_ext
    
    !> @var name_file Full path to the input or template file
    character(200) :: name_file
    
    !> @var linedata Buffer for reading file lines
    character(256) :: linedata

    ! Initialize directory name for templates
    name_dir = 'NAME_DIR'
    
    ! Commented out README functionality - now integrated into TurboRVB manual
    ! Original code for displaying README files has been disabled
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

    ! Set up input file handling
    name_ext = '.input'
    name_file = trim(name_tool)//trim(name_ext)
    
    ! Open the input file to check if it exists and has content
    open (unit=23, file=name_file, status='unknown', form='formatted')
    
    ! Count lines in the existing input file to check if it's empty
    ncount = 0
    do while (.true.)
        read (23, '(a)', end=110) linedata
        ncount = ncount + 1
    end do
110 continue
    
    ! If the input file is empty, generate a sample from template
    if (ncount .eq. 0) then ! only if the output file is empty
        ! Close and reopen the file for writing
        close (23)
        open (unit=23, file=name_file, status='unknown', form='formatted')
        
        ! Construct path to template file
        name_file = trim(name_dir)//'/template/'//trim(name_tool)//trim(name_ext)
        !write(6,*) name_file
        
        ! Open template file and copy its contents
        open (unit=22, file=name_file, status='unknown', form='formatted')
        ncount = 0
        do while (.true.)
            read (22, '(a)', end=120) linedata
            ncount = ncount + 1
            write (23, '(a)') trim(linedata)
        end do
120     continue

        ! Provide feedback to user about the generated sample input
        if (ncount .gt. 0) then
            write (6, *)
            write (6, *) ' Warning a sample input is given  ', trim(name_tool)//trim(name_ext)
            close (22)
            close (23)
        else
            ! If template doesn't exist, clean up empty files
            close (22, status='DELETE')
            close (23, status='DELETE')
        end if
    else
        ! If input file already has content, just close it safely
        close (23) ! close safely the existing file
    end if
    
    return
end

