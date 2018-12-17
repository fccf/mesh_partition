!> author: FANGCHAO
!  data: 2017/7/12
!
!## description
! this module provide some command line operation.
! contains: `get_option`,`have_option`
! the option usually is string, like `--help`,
! the arguments can be integer, real, logical, string


module system

  implicit none

  public :: get_option, have_option, mkdir, cwd, cd

  private


contains
  !=============================================================================
  subroutine mkdir(dir)
    !< make directory
    character(*), intent(in):: dir

    call system('mkdir '//dir)

  end subroutine mkdir
  !=============================================================================
  function cwd() result(path)
    !< get current work directory
    character(:),allocatable :: path
    character(256) :: path_

    integer :: ierr

    ierr = getcwd(path_)
    path = trim(adjustl(path_))

  end function cwd
  !=============================================================================
  subroutine cd(dir)
    !< change directory
    character(*),intent(in) :: dir
    integer :: ierr

    ierr = chdir(dir)

  end subroutine cd
  !=============================================================================
  function have_option(option) result(found)
    !< confirm if appear option
    character(*),intent(in) :: option
    logical   :: found

    integer :: argc
    character(32),allocatable :: argv(:)
    integer :: i
    found = .false.

    argc = command_argument_count()
    if(argc /=0) then
      allocate(argv(argc))
    end if

    do i = 1,argc
      call get_command_argument(i,argv(i))
    end do

    do i = 1,argc
      if(option == trim(adjustl(argv(i)))) then
        found = .true.
        exit
      end if
    enddo

  end function have_option
  !=============================================================================
  subroutine get_option(option,arg,found)
    !< get option for integer, real, logical, character(*)
    character(*),intent(in)      :: option
    class(*),intent(out)         :: arg
    logical,intent(out),optional :: found

    integer :: argc
    character(32),allocatable :: argv(:)
    integer :: i
    logical :: found_ = .false.

    ! arg  = 0
    argc = command_argument_count()
    if(argc /=0) then
      allocate(argv(argc))
    end if


    do i = 1,argc
      call get_command_argument(i,argv(i))
    end do

    do i = 1,argc
      if(option == trim(adjustl(argv(i)))) then
        select type (arg)
        type is(integer)
          read(argv(i+1),*) arg
        type is(real)
          read(argv(i+1),*) arg
        type is(logical)
          read(argv(i+1),*) arg
        type is(character(*))
          read(argv(i+1),*) arg
        class default
            error stop "unknow option type"
        end select

        found_ = .true.
        exit
      end if
    enddo

    if(present(found)) found = found_

  end subroutine get_option


end module system
