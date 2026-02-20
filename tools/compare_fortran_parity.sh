#!/usr/bin/env bash
# Author: Watosn
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FORTRAN_DIR="${1:-/Users/rmw/Documents/code/old_NRL}"
BUILD_DIR="${ROOT_DIR}/build/macos-clang-debug"
TMP_DIR="${TMPDIR:-/tmp}/msis_parity.$$"
mkdir -p "${TMP_DIR}"
trap 'rm -rf "${TMP_DIR}"' EXIT

if ! command -v gfortran >/dev/null 2>&1; then
  echo "gfortran not found"
  exit 1
fi

if [[ ! -d "${FORTRAN_DIR}" ]]; then
  echo "Fortran source dir not found: ${FORTRAN_DIR}"
  exit 1
fi

cat > "${TMP_DIR}/msis_ref_msiscalc.F90" <<'EOF'
program msis_ref_msiscalc
  use msis_constants, only : rp, dmissing
  use msis_init, only : msisinit
  use msis_calc, only : msiscalc
  implicit none

  integer :: iyd
  real(4) :: sec4, alt4, glat4, glong4, stl4, f107a4, f1074, apd4
  real(kind=rp) :: day, sec, alt, glat, glong, f107a, f107, ap(7), tn, dn(10)

  call msisinit()
  read(*,*) iyd, sec4, alt4, glat4, glong4, stl4, f107a4, f1074, apd4
  day = mod(iyd,1000)
  sec = sec4
  alt = alt4
  glat = glat4
  glong = glong4
  f107a = f107a4
  f107 = f1074
  ap = apd4

  call msiscalc(day,sec,alt,glat,glong,f107a,f107,ap,tn,dn)

  where (dn .ne. dmissing) dn = dn*1.0e-6_rp
  if (dn(1) .ne. dmissing) dn(1) = dn(1)*1.0e3_rp
  write(*,'(11(1X,ES24.16E3))') dn(5),dn(4),dn(2),dn(3),dn(7),dn(1),dn(6),dn(8),dn(9),dn(10),tn
end program msis_ref_msiscalc
EOF

gfortran -cpp -DDBLE -O0 -o "${TMP_DIR}/msis_ref_msiscalc" \
  "${FORTRAN_DIR}/msis_constants.F90" \
  "${FORTRAN_DIR}/msis_init.F90" \
  "${FORTRAN_DIR}/msis_utils.F90" \
  "${FORTRAN_DIR}/msis_gfn.F90" \
  "${FORTRAN_DIR}/msis_tfn.F90" \
  "${FORTRAN_DIR}/msis_dfn.F90" \
  "${FORTRAN_DIR}/msis_calc.F90" \
  "${TMP_DIR}/msis_ref_msiscalc.F90"

cp -f "${FORTRAN_DIR}/msis21.parm" "${TMP_DIR}/msis21.parm"

cmake --build --preset macos-clang-debug -j >/dev/null

INPUT_FILE="${ROOT_DIR}/data/msis2.1_test_in.txt"
FORTRAN_OUT="${TMP_DIR}/fortran.txt"
CPP_OUT="${TMP_DIR}/cpp.txt"

tail -n +2 "${INPUT_FILE}" | awk 'NF>=9{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | while read -r a b c d e f g h i; do
  echo "$a $b $c $d $e $f $g $h $i" | (cd "${TMP_DIR}" && ./msis_ref_msiscalc) >> "${FORTRAN_OUT}"
  "${BUILD_DIR}/msis21_cli" "$a" "$b" "$c" "$d" "$e" "$f" "$g" "$h" "$i" \
    | tail -n 1 | sed -E 's/[a-z0-9_]+=//g' >> "${CPP_OUT}"
done

paste "${FORTRAN_OUT}" "${CPP_OUT}" | awk '
BEGIN{
  tiny=1e-30;
  for(i=1;i<=11;i++){maxr[i]=0;maxa[i]=0;maxrg[i]=0}
}
{
  for(i=1;i<=11;i++){
    r=$i; o=$(i+11); d=o-r; ad=(d<0?-d:d); ar=(r!=0?ad/((r<0?-r:r)):ad);
    rg=ad/(((r<0?-r:r)>tiny)?(r<0?-r:r):tiny);
    if(ar>maxr[i])maxr[i]=ar;
    if(rg>maxrg[i])maxrg[i]=rg;
    if(ad>maxa[i])maxa[i]=ad;
  }
}
END{
  for(i=1;i<=11;i++) {
    printf("col%d max_rel=%.6e max_rel_guarded=%.6e max_abs=%.6e\n",i,maxr[i],maxrg[i],maxa[i]);
  }
}' 
