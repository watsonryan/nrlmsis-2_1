! Author: Watosn
module msis21_fortran_utils
  use iso_c_binding, only: c_double
  implicit none
contains
  function msis21_fortran_alt2gph(lat, alt) bind(C, name="msis21_fortran_alt2gph") result(gph)
    implicit none
    real(c_double), value :: lat
    real(c_double), value :: alt
    real(c_double) :: gph

    real(c_double), parameter  :: deg2rad = 0.017453292519943295d0
    real(c_double), parameter  :: a = 6378.1370d0 * 1d3
    real(c_double), parameter  :: finv = 298.257223563d0
    real(c_double), parameter  :: w = 7292115d-11
    real(c_double), parameter  :: GM = 398600.4418 * 1d9

    real(c_double), parameter  :: asq = a*a
    real(c_double), parameter  :: wsq = w*w
    real(c_double), parameter  :: f = 1.0d0 / finv
    real(c_double), parameter  :: esq = 2*f - f*f
    real(c_double), parameter  :: e = sqrt(esq)
    real(c_double), parameter  :: Elin = a*e
    real(c_double), parameter  :: Elinsq = Elin*Elin
    real(c_double), parameter  :: epr = e / (1-f)
    real(c_double), parameter  :: q0 = ((1.0d0 + 3.0d0/(epr*epr))*atan(epr) - 3.0d0/epr)/2.0d0
    real(c_double), parameter  :: U0 = -GM*atan(epr)/Elin - wsq*asq/3d0
    real(c_double), parameter  :: g0 = 9.80665d0
    real(c_double), parameter  :: GMdivElin = GM / Elin
    real(c_double), parameter  :: x0sq = 2d7**2
    real(c_double), parameter  :: Hsq = 1.2d7**2

    real(c_double) :: altm, sinsqlat, v, xsq, zsq
    real(c_double) :: rsqminElinsq, usq, cossqdelta, epru, atanepru, q, U, Vc

    altm = alt * 1000.0d0
    sinsqlat = sin(lat*deg2rad)**2
    v = a / sqrt(1-esq*sinsqlat)
    xsq = (v + altm)**2 * (1 - sinsqlat)
    zsq = (v*(1-esq) + altm)**2 * sinsqlat
    rsqminElinsq = xsq + zsq - Elinsq
    usq = rsqminElinsq/2.0d0 + sqrt(rsqminElinsq**2 / 4.0d0 + Elinsq*zsq)
    cossqdelta = zsq / usq

    epru = Elin / sqrt(usq)
    atanepru = atan(epru)
    q = ((1+3.0d0/(epru*epru))*atanepru - 3.0d0/epru)/2.0d0
    U = -GMdivElin * atanepru - wsq * ( asq * q * (cossqdelta - 1/3.0d0) / q0 ) / 2.0d0

    if (xsq .le. x0sq) then
      Vc = (wsq/2.0d0) * xsq
    else
      Vc = (wsq/2.0d0) * (Hsq*tanh((xsq-x0sq)/Hsq) + x0sq)
    endif
    U = U - Vc

    gph = (U - U0) / g0 / 1000.0d0
  end function msis21_fortran_alt2gph
end module msis21_fortran_utils
