module interp_mod
!----------------------------------------------------------------

! Code converted using TO_F90 by Alan Miller
! Date: 2024-03-13  Time: 09:28:36

contains

    SUBROUTINE interp(rjd, x, y, ut1, n, rjd_int, x_int, y_int, ut1_int)

!     This subroutine takes a series of x, y, and UT1-UTC values
!     and interpolates them to an epoch of choice. This routine
!     assumes that the values of x and y are in seconds of
!     arc and that UT1-UTC is in seconds of time. At least
!     one point before and one point after the epoch of the
!     interpolation point are necessary in order for the
!     interpolation scheme to work.

!     parameters are :
!     RJD     - array of the epochs of data (given in mjd)
!     X       - array of x polar motion (arcsec)
!     Y       - array of y polar motion (arcsec)
!     UT1     - array of UT1-UTC (sec)
!     n       - number of points in arrays
!     rjd_int - epoch for the interpolated value
!     x_int   - interpolated value of x
!     y_int   - interpolated value of y
!     ut1_int - interpolated value of ut1-utc

!     CALLED SUBROUTINE : LAGINT (Lagrange interpolation)
!                         PMUT1_OCEANS (Diurnal and semidiurnal oceanic effects)
!                         PM_GRAVI (Diurnal and semidiurnal lunisolar effects)

!      coded by Ch. BIZOUARD (Observatoire de Paris) : November 2002
!                                          Corrected : September 2007
        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN)         :: rjd(n)
        DOUBLE PRECISION, INTENT(IN)         :: x(n)
        DOUBLE PRECISION, INTENT(IN)         :: y(n)
        DOUBLE PRECISION, INTENT(IN)         :: ut1(n)
        INTEGER, INTENT(IN)                  :: n
        DOUBLE PRECISION, INTENT(IN)         :: rjd_int
        DOUBLE PRECISION, INTENT(OUT)            :: x_int
        DOUBLE PRECISION, INTENT(OUT)            :: y_int
        DOUBLE PRECISION, INTENT(OUT)            :: ut1_int
        DOUBLE PRECISION :: cor_x, cor_y, cor_ut1, cor_lod

        CALL lagint(rjd, x, n, rjd_int, x_int)

        CALL lagint(rjd, y, n, rjd_int, y_int)

        CALL lagint(rjd, ut1, n, rjd_int, ut1_int)

! --------------
! Oceanic effect
! --------------

        CALL pmut1_oceans(rjd_int, cor_x, cor_y, cor_ut1, cor_lod)

        x_int = x_int + cor_x
        y_int = y_int + cor_y
        ut1_int = ut1_int + cor_ut1

! Lunisolar effect
        CALL pm_gravi(rjd_int, cor_x, cor_y)

        x_int = x_int + cor_x
        y_int = y_int + cor_y

        RETURN

    END SUBROUTINE interp

!----------------------------------------------------------------

    SUBROUTINE lagint(x, y, n, xint, yout)

!     This subroutine performs lagrangian interpolation
!     within a set of (X,Y) pairs to give the y
!     value corresponding to xint. This program uses a
!     window of 4 data points to perform the interpolation.
!     if the window size needs to be changed, this can be
!     done by changing the indices in the do loops for
!     variables m and j.

!     PARAMETERS ARE :
!     X     - array of values of the independent variable
!     Y     - array of function values corresponding to x
!     n     - number of points
!     xint  - the x-value for which estimate of y is desired
!     yout  - the y value returned to caller
        IMPLICIT NONE

        REAL(8), INTENT(IN)                       :: x(n)
        REAL(8), INTENT(IN)                       :: y(n)
        INTEGER, INTENT(IN)                      :: n
        REAL(8), INTENT(IN)                       :: xint
        REAL(8), INTENT(OUT)                      :: yout
        REAL(8) :: term
        INTEGER :: i, j, k, m

        yout = 0.d0
        DO i = 1, n - 1
            IF (xint >= x(i) .AND. xint < x(i + 1)) k = i
        END DO

        IF (k < 2) k = 2
        IF (k > n - 2) k = n - 2

        DO m = k - 1, k + 2
            term = y(m)
            DO j = k - 1, k + 2
                IF (m /= j) THEN
                    term = term*(xint - x(j))/(x(m) - x(j))
                END IF
            END DO
            yout = yout + term
        END DO

        RETURN
    END SUBROUTINE lagint

!----------------------------------------------------------------

    SUBROUTINE pmut1_oceans(rjd, cor_x, cor_y, cor_ut1, cor_lod)

!    This subroutine provides, in time domain, the diurnal/subdiurnal
!    tidal effets on polar motion ("), UT1 (s) and LOD (s). The tidal terms,
!    listed in the program above, have been extracted from the procedure
!    ortho_eop.f coed by Eanes in 1997.

!    N.B.:  The fundamental lunisolar arguments are those of Simon et al.

!    These corrections should be added to "average"
!    EOP values to get estimates of the instantaneous values.

!     PARAMETERS ARE :
!     rjd      - epoch of interest given in mjd
!     cor_x    - tidal correction in x (sec. of arc)
!     cor_y    - tidal correction in y (sec. of arc)
!     cor_ut1  - tidal correction in UT1-UTC (sec. of time)
!     cor_lod  - tidal correction in length of day (sec. of time)

!     coded by Ch. Bizouard (2002), initially coded by McCarthy and
!     D.Gambis(1997) for the 8 prominent tidal waves.
        IMPLICIT NONE

        REAL(8), INTENT(IN)                   :: rjd
        REAL(8), INTENT(OUT)                      :: cor_x
        REAL(8), INTENT(OUT)                      :: cor_y
        REAL(8), INTENT(OUT)                      :: cor_ut1
        REAL(8), INTENT(OUT)                      :: cor_lod

        INTEGER, PARAMETER :: nlines = 71
        DOUBLE PRECISION :: arg(6)   ! Array of the tidal arguments
        DOUBLE PRECISION :: darg(6) ! Array of their time derivative

        REAL(4) :: xcos(nlines), xsin(nlines), ycos(nlines), ysin(nlines), utcos(nlines), utsin(nlines)

        REAL(8) :: t, ag, dag, halfpi, secrad
        INTEGER :: narg(nlines, 6), i, j

        halfpi = 1.5707963267948966D0
        secrad = 2.d0*halfpi/(180.d0*3600.d0)

!  Oceanic tidal terms present in x (microas),y(microas),ut1(microseconds)
!  NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments.

        DATA(narg(j, 1), narg(j, 2), narg(j, 3), narg(j, 4), narg(j, 5), narg(j, 6), &
             xsin(j), xcos(j), ysin(j), ycos(j), utsin(j), utcos(j), j=1, nlines)/ &
            1, -1, 0, -2, -2, -2, -0.05, 0.94, -0.94, -0.05, 0.396, -0.078, &
            1, -2, 0, -2, 0, -1, 0.06, 0.64, -0.64, 0.06, 0.195, -0.059, &
            1, -2, 0, -2, 0, -2, 0.30, 3.42, -3.42, 0.30, 1.034, -0.314, &
            1, 0, 0, -2, -2, -1, 0.08, 0.78, -0.78, 0.08, 0.224, -0.073, &
            1, 0, 0, -2, -2, -2, 0.46, 4.15, -4.15, 0.45, 1.187, -0.387, &
            1, -1, 0, -2, 0, -1, 1.19, 4.96, -4.96, 1.19, 0.966, -0.474, &
            1, -1, 0, -2, 0, -2, 6.24, 26.31, -26.31, 6.23, 5.118, -2.499, &
            1, 1, 0, -2, -2, -1, 0.24, 0.94, -0.94, 0.24, 0.172, -0.090, &
            1, 1, 0, -2, -2, -2, 1.28, 4.99, -4.99, 1.28, 0.911, -0.475, &
            1, 0, 0, -2, 0, 0, -0.28, -0.77, 0.77, -0.28, -0.093, 0.070, &
            1, 0, 0, -2, 0, -1, 9.22, 25.06, -25.06, 9.22, 3.025, -2.280, &
            1, 0, 0, -2, 0, -2, 48.82, 132.91, -132.90, 48.82, 16.020, -12.069, &
            1, -2, 0, 0, 0, 0, -0.32, -0.86, 0.86, -0.32, -0.103, 0.078, &
            1, 0, 0, 0, -2, 0, -0.66, -1.72, 1.72, -0.66, -0.194, 0.154, &
            1, -1, 0, -2, 2, -2, -0.42, -0.92, 0.92, -0.42, -0.083, 0.074, &
            1, 1, 0, -2, 0, -1, -0.30, -0.64, 0.64, -0.30, -0.057, 0.050, &
            1, 1, 0, -2, 0, -2, -1.61, -3.46, 3.46, -1.61, -0.308, 0.271, &
            1, -1, 0, 0, 0, 0, -4.48, -9.61, 9.61, -4.48, -0.856, 0.751, &
            1, -1, 0, 0, 0, -1, -0.90, -1.93, 1.93, -0.90, -0.172, 0.151, &
            1, 1, 0, 0, -2, 0, -0.86, -1.81, 1.81, -0.86, -0.161, 0.137, &
            1, 0, -1, -2, 2, -2, 1.54, 3.03, -3.03, 1.54, 0.315, -0.189, &
            1, 0, 0, -2, 2, -1, -0.29, -0.58, 0.58, -0.29, -0.062, 0.035, &
            1, 0, 0, -2, 2, -2, 26.13, 51.25, -51.25, 26.13, 5.512, -3.095, &
            1, 0, 1, -2, 2, -2, -0.22, -0.42, 0.42, -0.22, -0.047, 0.025, &
            1, 0, -1, 0, 0, 0, -0.61, -1.20, 1.20, -0.61, -0.134, 0.070, &
            1, 0, 0, 0, 0, 1, 1.54, 3.00, -3.00, 1.54, 0.348, -0.171, &
            1, 0, 0, 0, 0, 0, -77.48, -151.74, 151.74, -77.48, -17.620, 8.548, &
            1, 0, 0, 0, 0, -1, -10.52, -20.56, 20.56, -10.52, -2.392, 1.159, &
            1, 0, 0, 0, 0, -2, 0.23, 0.44, -0.44, 0.23, 0.052, -0.025, &
            1, 0, 1, 0, 0, 0, -0.61, -1.19, 1.19, -0.61, -0.144, 0.065, &
            1, 0, 0, 2, -2, 2, -1.09, -2.11, 2.11, -1.09, -0.267, 0.111, &
            1, -1, 0, 0, 2, 0, -0.69, -1.43, 1.43, -0.69, -0.288, 0.043, &
            1, 1, 0, 0, 0, 0, -3.46, -7.28, 7.28, -3.46, -1.610, 0.187, &
            1, 1, 0, 0, 0, -1, -0.69, -1.44, 1.44, -0.69, -0.320, 0.037, &
            1, 0, 0, 0, 2, 0, -0.37, -1.06, 1.06, -0.37, -0.407, -0.005, &
            1, 2, 0, 0, 0, 0, -0.17, -0.51, 0.51, -0.17, -0.213, -0.005, &
            1, 0, 0, 2, 0, 2, -1.10, -3.42, 3.42, -1.09, -1.436, -0.037, &
            1, 0, 0, 2, 0, 1, -0.70, -2.19, 2.19, -0.70, -0.921, -0.023, &
            1, 0, 0, 2, 0, 0, -0.15, -0.46, 0.46, -0.15, -0.193, -0.005, &
            1, 1, 0, 2, 0, 2, -0.03, -0.59, 0.59, -0.03, -0.396, -0.024, &
            1, 1, 0, 2, 0, 1, -0.02, -0.38, 0.38, -0.02, -0.253, -0.015, &
            2, -3, 0, -2, 0, -2, -0.49, -0.04, 0.63, 0.24, -0.089, -0.011, &
            2, -1, 0, -2, -2, -2, -1.33, -0.17, 1.53, 0.68, -0.224, -0.032, &
            2, -2, 0, -2, 0, -2, -6.08, -1.61, 3.13, 3.35, -0.637, -0.177, &
            2, 0, 0, -2, -2, -2, -7.59, -2.05, 3.44, 4.23, -0.745, -0.222, &
            2, 0, 1, -2, -2, -2, -0.52, -0.14, 0.22, 0.29, -0.049, -0.015, &
            2, -1, -1, -2, 0, -2, 0.47, 0.11, -0.10, -0.27, 0.033, 0.013, &
            2, -1, 0, -2, 0, -1, 2.12, 0.49, -0.41, -1.23, 0.141, 0.058, &
            2, -1, 0, -2, 0, -2, -56.87, -12.93, 11.15, 32.88, -3.795, -1.556, &
            2, -1, 1, -2, 0, -2, -0.54, -0.12, 0.10, 0.31, -0.035, -0.015, &
            2, 1, 0, -2, -2, -2, -11.01, -2.40, 1.89, 6.41, -0.698, -0.298, &
            2, 1, 1, -2, -2, -2, -0.51, -0.11, 0.08, 0.30, -0.032, -0.014, &
            2, -2, 0, -2, 2, -2, 0.98, 0.11, -0.11, -0.58, 0.050, 0.022, &
            2, 0, -1, -2, 0, -2, 1.13, 0.11, -0.13, -0.67, 0.056, 0.025, &
            2, 0, 0, -2, 0, -1, 12.32, 1.00, -1.41, -7.31, 0.605, 0.266, &
            2, 0, 0, -2, 0, -2, -330.15, -26.96, 37.58, 195.92, -16.195, -7.140, &
            2, 0, 1, -2, 0, -2, -1.01, -0.07, 0.11, 0.60, -0.049, -0.021, &
            2, -1, 0, -2, 2, -2, 2.47, -0.28, -0.44, -1.48, 0.111, 0.034, &
            2, 1, 0, -2, 0, -2, 9.40, -1.44, -1.88, -5.65, 0.425, 0.117, &
            2, -1, 0, 0, 0, 0, -2.35, 0.37, 0.47, 1.41, -0.106, -0.029, &
            2, -1, 0, 0, 0, -1, -1.04, 0.17, 0.21, 0.62, -0.047, -0.013, &
            2, 0, -1, -2, 2, -2, -8.51, 3.50, 3.29, 5.11, -0.437, -0.019, &
            2, 0, 0, -2, 2, -2, -144.13, 63.56, 59.23, 86.56, -7.547, -0.159, &
            2, 0, 1, -2, 2, -2, 1.19, -0.56, -0.52, -0.72, 0.064, 0.000, &
            2, 0, 0, 0, 0, 1, 0.49, -0.25, -0.23, -0.29, 0.027, -0.001, &
            2, 0, 0, 0, 0, 0, -38.48, 19.14, 17.72, 23.11, -2.104, 0.041, &
            2, 0, 0, 0, 0, -1, -11.44, 5.75, 5.32, 6.87, -0.627, 0.015, &
            2, 0, 0, 0, 0, -2, -1.24, 0.63, 0.58, 0.75, -0.068, 0.002, &
            2, 1, 0, 0, 0, 0, -1.77, 1.79, 1.71, 1.04, -0.146, 0.037, &
            2, 1, 0, 0, 0, -1, -0.77, 0.78, 0.75, 0.45, -0.064, 0.017, &
            2, 0, 0, 2, 0, 2, -0.33, 0.62, 0.65, 0.19, -0.049, 0.018/

        t = (rjd - 51544.5D0)/36525.0D0  ! julian century

! Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
! et leur derivee temporelle

        arg(1) = (67310.54841D0 + (876600D0*3600D0 + 8640184.812866D0)*t + &
                  0.093104D0*t**2 - 6.2D-6*t**3)*15.0D0 + 648000.0D0
        arg(1) = DMOD(arg(1), 1296000D0)*secrad

        darg(1) = (876600D0*3600D0 + 8640184.812866D0 &
                   + 2.d0*0.093104D0*t - 3.d0*6.2D-6*t**2)*15.d0
        darg(1) = darg(1)*secrad/36525.0D0   ! rad/day

        arg(2) = -0.00024470D0*t**4 + 0.051635D0*t**3 + 31.8792D0*t**2 &
                 + 1717915923.2178D0*t + 485868.249036D0
        arg(2) = DMOD(arg(2), 1296000D0)*secrad

        darg(2) = -4.d0*0.00024470D0*t**3 + 3.d0*0.051635D0*t**2 &
                  + 2.d0*31.8792D0*t + 1717915923.2178D0
        darg(2) = darg(2)*secrad/36525.0D0   ! rad/day

        arg(3) = -0.00001149D0*t**4 + 0.000136D0*t**3 - 0.5532D0*t**2 &
                 + 129596581.0481D0*t + 1287104.79305D0
        arg(3) = DMOD(arg(3), 1296000D0)*secrad

        darg(3) = -4.d0*0.00001149D0*t**3 - 3.d0*0.000136D0*t**2 &
                  - 2.d0*0.5532D0*t + 129596581.0481D0
        darg(3) = darg(3)*secrad/36525.0D0   ! rad/day

        arg(4) = 0.00000417D0*t**4 - 0.001037D0*t**3 - 12.7512D0*t**2 &
                 + 1739527262.8478D0*t + 335779.526232D0
        arg(4) = DMOD(arg(4), 1296000D0)*secrad

        darg(4) = 4.d0*0.00000417D0*t**3 - 3.d0*0.001037D0*t**2 &
                  - 2.d0*12.7512D0*t + 1739527262.8478D0
        darg(4) = darg(4)*secrad/36525.0D0   ! rad/day

        arg(5) = -0.00003169D0*t**4 + 0.006593D0*t**3 - 6.3706D0*t**2 &
                 + 1602961601.2090D0*t + 1072260.70369D0
        arg(5) = DMOD(arg(5), 1296000D0)*secrad

        darg(5) = -4.d0*0.00003169D0*t**3 + 3.d0*0.006593D0*t**2 &
                  - 2.d0*6.3706D0*t + 1602961601.2090D0
        darg(5) = darg(5)*secrad/36525.0D0   ! rad/day

        arg(6) = -0.00005939D0*t**4 + 0.007702D0*t**3 + 7.4722D0*t**2 &
                 - 6962890.2665D0*t + 450160.398036D0
        arg(6) = DMOD(arg(6), 1296000D0)*secrad

        darg(6) = -4.d0*0.00005939D0*t**3 + 3.d0*0.007702D0*t**2 &
                  + 2.d0*7.4722D0*t - 6962890.2665D0
        darg(6) = darg(6)*secrad/36525.0D0   ! rad/day

! CORRECTIONS

        cor_x = 0.d0
        cor_y = 0.d0
        cor_ut1 = 0.d0
        cor_lod = 0.d0

        DO j = 1, nlines

            ag = 0.d0
            dag = 0.d0
            DO i = 1, 6
                ag = ag + DBLE(narg(j, i))*arg(i)
                dag = dag + DBLE(narg(j, i))*darg(i)
            END DO
            ag = DMOD(ag, 4.d0*halfpi)

            cor_x = cor_x + DBLE(xcos(j))*DCOS(ag) + DBLE(xsin(j))*DSIN(ag)
            cor_y = cor_y + DBLE(ycos(j))*DCOS(ag) + DBLE(ysin(j))*DSIN(ag)
            cor_ut1 = cor_ut1 + DBLE(utcos(j))*DCOS(ag) + DBLE(utsin(j))*DSIN(ag)
            cor_lod = cor_lod - (-DBLE(utcos(j))*DSIN(ag) &
                                 + DBLE(utsin(j))*DCOS(ag))*dag

        END DO

        cor_x = cor_x*1.0D-6   ! arcsecond (")
        cor_y = cor_y*1.0D-6   ! arcsecond (")
        cor_ut1 = cor_ut1*1.0D-6 ! second (s)
        cor_lod = cor_lod*1.0D-6 ! second (s)

        RETURN
    END SUBROUTINE pmut1_oceans

!----------------------------------------------------------------

    SUBROUTINE pm_gravi(rjd, cor_x, cor_y)

!    This subroutine provides, in time domain, the diurnal
!    lunisolar effet on polar motion (")

!    N.B.:  The fundamental lunisolar arguments are those of Simon et al.

!    These corrections should be added to "average"
!    EOP values to get estimates of the instantaneous values.

!     PARAMETERS ARE :
!     rjd      - epoch of interest given in mjd
!     cor_x    - tidal correction in x (sec. of arc)
!     cor_y    - tidal correction in y (sec. of arc)

!     coded by Ch. Bizouard (2002)
        IMPLICIT DOUBLE PRECISION(a - h, o - z)

        DOUBLE PRECISION, INTENT(IN)         :: rjd
        DOUBLE PRECISION, INTENT(OUT)            :: cor_x
        DOUBLE PRECISION, INTENT(OUT)            :: cor_y

        INTEGER, PARAMETER :: nlines = 10
        DOUBLE PRECISION :: arg(6)    ! Array of the tidal arguments
        REAL(4) :: xcos(nlines), xsin(nlines), ycos(nlines), ysin(nlines)
        INTEGER :: narg(nlines, 6)

        halfpi = 1.5707963267948966D0
        secrad = 2.d0*halfpi/(180.d0*3600.d0)

!  Diurnal lunisolar tidal terms present in x (microas),y(microas)
!  NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments.

        DATA(narg(j, 1), narg(j, 2), narg(j, 3), narg(j, 4), narg(j, 5), narg(j, 6), &
             xsin(j), xcos(j), ysin(j), ycos(j), j=1, nlines)/ &
            1, -1, 0, -2, 0, -1, -.44, .25, -.25, -.44, &
            1, -1, 0, -2, 0, -2, -2.31, 1.32, -1.32, -2.31, &
            1, 1, 0, -2, -2, -2, -.44, .25, -.25, -.44, &
            1, 0, 0, -2, 0, -1, -2.14, 1.23, -1.23, -2.14, &
            1, 0, 0, -2, 0, -2, -11.36, 6.52, -6.52, -11.36, &
            1, -1, 0, 0, 0, 0, .84, -.48, .48, .84, &
            1, 0, 0, -2, 2, -2, -4.76, 2.73, -2.73, -4.76, &
            1, 0, 0, 0, 0, 0, 14.27, -8.19, 8.19, 14.27, &
            1, 0, 0, 0, 0, -1, 1.93, -1.11, 1.11, 1.93, &
            1, 1, 0, 0, 0, 0, .76, -.43, .43, .76/

        t = (rjd - 51544.5D0)/36525.0D0  ! julian century

! Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
! et leur derivee temporelle

        arg(1) = (67310.54841D0 + (876600D0*3600D0 + 8640184.812866D0)*t + &
                  0.093104D0*t**2 - 6.2D-6*t**3)*15.0D0 + 648000.0D0
        arg(1) = DMOD(arg(1), 1296000D0)*secrad

        arg(2) = -0.00024470D0*t**4 + 0.051635D0*t**3 + 31.8792D0*t**2 &
                 + 1717915923.2178D0*t + 485868.249036D0
        arg(2) = DMOD(arg(2), 1296000D0)*secrad

        arg(3) = -0.00001149D0*t**4 - 0.000136D0*t**3 - 0.5532D0*t**2 &
                 + 129596581.0481D0*t + 1287104.79305D0
        arg(3) = DMOD(arg(3), 1296000D0)*secrad

        arg(4) = 0.00000417D0*t**4 - 0.001037D0*t**3 - 12.7512D0*t**2 &
                 + 1739527262.8478D0*t + 335779.526232D0
        arg(4) = DMOD(arg(4), 1296000D0)*secrad

        arg(5) = -0.00003169D0*t**4 + 0.006593D0*t**3 - 6.3706D0*t**2 &
                 + 1602961601.2090D0*t + 1072260.70369D0
        arg(5) = DMOD(arg(5), 1296000D0)*secrad

        arg(6) = -0.00005939D0*t**4 + 0.007702D0*t**3 + 7.4722D0*t**2 &
                 - 6962890.2665D0*t + 450160.398036D0
        arg(6) = DMOD(arg(6), 1296000D0)*secrad

! CORRECTIONS

        cor_x = 0.d0
        cor_y = 0.d0

        DO j = 1, nlines

            ag = 0.d0
            DO i = 1, 6
                ag = ag + DBLE(narg(j, i))*arg(i)
            END DO
            ag = DMOD(ag, 4.d0*halfpi)

            cor_x = cor_x + DBLE(xcos(j))*DCOS(ag) + DBLE(xsin(j))*DSIN(ag)
            cor_y = cor_y + DBLE(ycos(j))*DCOS(ag) + DBLE(ysin(j))*DSIN(ag)

        END DO

        cor_x = cor_x*1.0D-6   ! arcsecond (")
        cor_y = cor_y*1.0D-6   ! arcsecond (")

        RETURN

    END SUBROUTINE pm_gravi

end module interp_mod
