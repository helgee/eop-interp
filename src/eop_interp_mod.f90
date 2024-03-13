module eop_interp_mod

    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    integer, parameter :: dp = kind(0.d0)

contains

    ! Subroutine: stop_error
    !   Print an error string to stderr and stop the program with exitcode 1.
    !
    ! Parameter:
    !   str - Error string
    subroutine stop_error(str, file, line)
        character(len=*), intent(in) :: str
        character(len=*), intent(in), optional :: file
        integer, intent(in), optional :: line
        if (present(file)) write (error_unit, *) "File:", file
        if (present(line)) write (error_unit, *) "Line:", line
        write (error_unit, *) str
        stop 1
    end subroutine stop_error

end module eop_interp_mod
