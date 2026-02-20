! Author: Watosn
module msis21_fortran_bridge
  use iso_c_binding, only: c_int, c_float, c_double
  use msis_constants, only: rp
  use msis_init, only: msisinit, initflag
  use msis_calc, only: msiscalc
  implicit none
contains
  subroutine msis21_fortran_eval(iyd, sec, alt, glat, glon, f107a, f107, apd, tn, dn, status) &
      bind(C, name="msis21_fortran_eval")
    integer(c_int), value :: iyd
    real(c_float), value :: sec, alt, glat, glon, f107a, f107, apd
    real(c_double) :: tn
    real(c_double) :: dn(10)
    integer(c_int) :: status

    real(kind=rp) :: day_rp, sec_rp, alt_rp, glat_rp, glon_rp, f107a_rp, f107_rp, ap_rp(7), tn_rp, dn_rp(10)

    if (.not. initflag) then
#ifdef MSIS21_PARM_FILE
      call msisinit(parmfile=MSIS21_PARM_FILE)
#else
      call msisinit()
#endif
    endif

    day_rp = mod(iyd, 1000)
    sec_rp = sec
    alt_rp = alt
    glat_rp = glat
    glon_rp = glon
    f107a_rp = f107a
    f107_rp = f107
    ap_rp = apd

    call msiscalc(day_rp, sec_rp, alt_rp, glat_rp, glon_rp, f107a_rp, f107_rp, ap_rp, tn_rp, dn_rp)

    tn = tn_rp
    dn = dn_rp
    status = 0_c_int
  end subroutine msis21_fortran_eval
end module msis21_fortran_bridge
