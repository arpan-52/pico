#!/usr/bin/env bash
# Build pico_image. Requires cmake>=3.25, gcc>=7 (C++17), cfitsio.
# Honours CC / CXX if set; otherwise picks the first gcc/g++ on PATH.
# CFITSIO_INCLUDE_DIR / CFITSIO_LIB override the cfitsio location.
set -euo pipefail
: "${CC:=$(command -v gcc)}"
: "${CXX:=$(command -v g++)}"
export CC CXX
echo "Using CC=$CC"
echo "Using CXX=$CXX"

# Auto-detect cfitsio. Prefer a matching header+library pair in $SOFT_ROOT
# (paramrudra cluster layout) before falling back to system /usr.
SOFT_ROOT_DEFAULT=/lustre_archive/apps/astrosoft/test_imaging/soft
: "${SOFT_ROOT:=$SOFT_ROOT_DEFAULT}"

CFITSIO_ARGS=()
if [[ -z "${CFITSIO_INCLUDE_DIR:-}" || -z "${CFITSIO_LIB:-}" ]]; then
  for prefix in "$SOFT_ROOT" /usr/local /usr; do
    if [[ -f "$prefix/include/fitsio.h" ]]; then
      for libdir in lib64 lib; do
        for ext in so a; do
          if [[ -f "$prefix/$libdir/libcfitsio.$ext" ]]; then
            : "${CFITSIO_INCLUDE_DIR:=$prefix/include}"
            : "${CFITSIO_LIB:=$prefix/$libdir/libcfitsio.$ext}"
            break 3
          fi
        done
      done
    fi
  done
fi

if [[ -n "${CFITSIO_INCLUDE_DIR:-}" && -n "${CFITSIO_LIB:-}" ]]; then
  echo "Using cfitsio header: $CFITSIO_INCLUDE_DIR/fitsio.h"
  echo "Using cfitsio lib   : $CFITSIO_LIB"
  CFITSIO_ARGS=(
    -DCFITSIO_INCLUDE_DIR="$CFITSIO_INCLUDE_DIR"
    -DCFITSIO_LIB="$CFITSIO_LIB"
  )
else
  echo "cfitsio: letting cmake find_library auto-detect"
fi

# The cluster conda env exports -ffast-math in CFLAGS/CXXFLAGS. That breaks
# IEEE NaN/inf semantics, which both our NaN guards AND SuperNOVAS's NaN-sentinel
# precession/nutation caches depend on — with fast-math the uvw J2000 rotation
# silently returns an all-zero matrix and the dirty image fills with NaN.
# Append the negations LAST (gcc honours the final flag) so they win over any
# inherited -ffast-math without us having to know the full inherited flag set.
SAFE_MATH="-fno-fast-math -fno-finite-math-only -fsigned-zeros -ftrapping-math"

cd "$(dirname "$0")/pico_image"
cmake -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" \
  -DCMAKE_C_FLAGS="${CFLAGS:-} $SAFE_MATH" \
  -DCMAKE_CXX_FLAGS="${CXXFLAGS:-} $SAFE_MATH" \
  "${CFITSIO_ARGS[@]}" "$@"
cmake --build build -j"$(nproc)"
echo
echo "Built: $(pwd)/build/pico_image"
