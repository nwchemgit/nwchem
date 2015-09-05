Module err_func

    Public erf

Contains

Pure Function erf(x) Result(fn_val)
!-----------------------------------------------------------------------
!             EVALUATION OF THE REAL ERROR FUNCTION
! Based upon a Fortran 66 routine in the Naval Surface Warfare Center's
! Mathematics Library (1993 version).
! Adapted by Alan.Miller @ vic.cmis.csiro.au
! Adapted for DIM by Seth Morton
!-----------------------------------------------------------------------
!    use constants
    Implicit None
#include "dimqm_constants.fh"
    Real (KINDR), Intent(In) :: x
    Real (KINDR)             :: fn_val

!   Local variables
    Real (KINDR), Parameter :: c = .564189583547756_KINDR
    Real (KINDR), Parameter ::  &
        a(5) = (/ .771058495001320E-04_KINDR, -.133733772997339E-02_KINDR,    &
                  .323076579225834E-01_KINDR,  .479137145607681E-01_KINDR,    &
                  .128379167095513E+00_KINDR /),                              &
        b(3) = (/ .301048631703895E-02_KINDR,  .538971687740286E-01_KINDR,    &
                  .375795757275549E+00_KINDR /),                              &
        p(8) = (/ -1.36864857382717E-07_KINDR, 5.64195517478974E-01_KINDR,    &
                   7.21175825088309E+00_KINDR, 4.31622272220567E+01_KINDR,    &
                   1.52989285046940E+02_KINDR, 3.39320816734344E+02_KINDR,    &
                   4.51918953711873E+02_KINDR, 3.00459261020162E+02_KINDR /), &
        q(8) = (/  1.00000000000000E+00_KINDR, 1.27827273196294E+01_KINDR,    &
                   7.70001529352295E+01_KINDR, 2.77585444743988E+02_KINDR,    &
                   6.38980264465631E+02_KINDR, 9.31354094850610E+02_KINDR,    &
                   7.90950925327898E+02_KINDR, 3.00459260956983E+02_KINDR /), &
        r(5) = (/  2.10144126479064E+00_KINDR, 2.62370141675169E+01_KINDR,    &
                   2.13688200555087E+01_KINDR, 4.65807828718470E+00_KINDR,    &
                   2.82094791773523E-01_KINDR /),                             &
        s(4) = (/  9.41537750555460E+01_KINDR, 1.87114811799590E+02_KINDR,    &
                   9.90191814623914E+01_KINDR, 1.80124575948747E+01_KINDR /)
    Real (KINDR) :: ax, bot, t, top, x2

    ax = ABS(x)

    if (ax <= HALF) then
        t = x*x
        top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + ONE
        bot = ((b(1)*t + b(2))*t + b(3))*t + ONE
        fn_val = x*(top/bot)
        return
    end if

    if (ax <= FOUR) then
        top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax  &
            + p(6))*ax + p(7))*ax + p(8)
        bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax  &
            + q(6))*ax + q(7))*ax + q(8)
        fn_val = HALF + (HALF - EXP(-x*x)*top/bot)
        if (x < ZERO) fn_val = -fn_val
        return
    end if

    if (ax < 5.8_KINDR) then
        x2 = x*x
        t = ONE / x2
        top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
        bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + ONE
        fn_val = (c - top/(x2*bot)) / ax
        fn_val = HALF + (HALF - EXP(-x2)*fn_val)
        if (x < ZERO) fn_val = -fn_val
        return
    end if

    fn_val = SIGN(ONE, x)
    return
End Function erf

End Module err_func
