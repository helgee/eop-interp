program eop_interp

    use eop_interp_mod
    use interp_mod

    implicit none

    real(dp) :: mjd0

    integer :: i
    integer :: iostat
    integer :: n_args
    integer :: n_lines = 0
    integer :: u
    character(len=1024) :: filename
    character(len=18) :: buf
    real(dp), dimension(:), allocatable :: mjd
    real(dp), dimension(:), allocatable :: x
    real(dp), dimension(:), allocatable :: y
    character(len=10), dimension(:), allocatable :: typ
    real(dp), dimension(:), allocatable :: delta_ut1_utc
    real(dp) :: x_int
    real(dp) :: y_int
    real(dp) :: ut1_int
    ! character(len=128) :: tmp

    n_args = command_argument_count()

    if (n_args /= 2) then
        call stop_error("Needs two arguments")
    end if

    call get_command_argument(1, filename)
    open (newunit=u, file=filename, status="old", action="read")

    do
        read (u, *, iostat=iostat)
        if (iostat /= 0) exit
        n_lines = n_lines + 1
    end do

    rewind (u)

    ! Skip header
    n_lines = n_lines - 1
    read (u, *)

    allocate (mjd(n_lines))
    allocate (x(n_lines))
    allocate (y(n_lines))
    allocate (typ(n_lines))
    allocate (delta_ut1_utc(n_lines))

    do i = 1, n_lines
        read (u, *) mjd(i), x(i), y(i), typ(i), delta_ut1_utc(i)
    end do

    close (u)

    call get_command_argument(2, buf)
    read (buf, *) mjd0

    call interp(mjd, x, y, delta_ut1_utc, n_lines, mjd0, x_int, y_int, ut1_int)

    print *, "Data: ", trim(filename)
    print *, "Date:", mjd0
    print *, "Interpolated pole x: ", x_int
    print *, "Interpolated pole y: ", y_int
    print *, "Interpolated UT1-UTC:", ut1_int

end program eop_interp
