program main
use constants
use structs
use init_final
use correlate_pivot
!use assign_eigenchann
use stdio
implicit none 
character(len=longlong_str) :: DIR

call cmdline_parser(DIR)
call step1(DIR)
!call step2

contains 
  subroutine cmdline_parser(DIR)
    implicit none
    character(len = longlong_str) :: arg
    integer :: i, ie
    character(len=longlong_str), intent(out) :: DIR
    !if(iargc() == 0) stop "No arguments! you must specify:&
    !  &-p $some_where_to_find_the_smat"

    DIR ='./'
    i = 1
    do while( i <= iargc() )
       call getarg( i, arg)
       if((trim(arg) == '-p') .or. (trim(arg) == '-path')) then
          i = i + 1
          call getarg( i, arg)
          read(arg, '(A)' ) DIR
          if(index(DIR,'/') == 0) write(DIR,"(A,'/')") trim(DIR)
       end if
       i = i + 1
   end do
  end subroutine cmdline_parser
end program main
