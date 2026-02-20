! Author: Watosn
module msis21_ref_bridge
  use iso_c_binding
  use msis_init, only : msisinit
  implicit none
contains

  subroutine msis21_ref_init(parm_dir_c, status) bind(c, name="msis21_ref_init")
    character(kind=c_char), intent(in) :: parm_dir_c(*)
    integer(c_int), intent(out) :: status

    character(len=512) :: parm_dir
    integer :: i

    parm_dir = ''
    i = 1
    do while (i <= len(parm_dir))
      if (parm_dir_c(i) == c_null_char) exit
      parm_dir(i:i) = parm_dir_c(i)
      i = i + 1
    end do

    call msisinit(parmpath=trim(parm_dir), parmfile='msis21.parm')
    status = 0_c_int
  end subroutine msis21_ref_init

  subroutine msis21_ref_gtd8d(iyd, sec, alt, glat, glong, stl, f107a, f107, ap_daily, d, t, status) &
      bind(c, name="msis21_ref_gtd8d")
    use iso_c_binding
    integer(c_int), value :: iyd
    real(c_float), value :: sec, alt, glat, glong, stl, f107a, f107, ap_daily
    real(c_float), intent(out) :: d(10), t(2)
    integer(c_int), intent(out) :: status

    real(c_float) :: ap(7)
    integer(c_int) :: mass

    ap = ap_daily
    mass = 48
    call gtd8d(iyd, sec, alt, glat, glong, stl, f107a, f107, ap, mass, d, t)
    status = 0_c_int
  end subroutine msis21_ref_gtd8d

end module msis21_ref_bridge
