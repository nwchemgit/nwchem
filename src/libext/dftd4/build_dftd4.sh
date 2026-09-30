#!/usr/bin/env bash
# DFTD4 interface build for NWChem.
# Developed with assistance from OpenAI Codex
# integration and numerical validation were
# reviewed and tested by MPH.

set -euo pipefail

version=${1:-4.2.0}
prefix=${2:-"$PWD/install"}
archive="dftd4-${version}.tar.gz"
source_dir="dftd4-${version}"
url="https://github.com/dftd4/dftd4/archive/refs/tags/v${version}.tar.gz"

if [[ -f "$archive" ]] && gzip -t "$archive" 2>/dev/null; then
    echo "using existing $archive"
else
    echo "downloading DFTD4 $version"
    curl -fL --retry 3 "$url" -o "$archive"
    gzip -t "$archive"
fi

mkdir -p "$source_dir"
tar -xzf "$archive" -C "$source_dir" --strip-components=1
ln -sfn "$source_dir" dftd4

cmake_command=${CMAKE:-cmake}
blas_size=${DFTD4_BLAS_SIZE:-4}
blas_libraries=${DFTD4_BLAS_LIBRARIES:-}
lapack_libraries=${DFTD4_LAPACK_LIBRARIES:-}
if [[ "$blas_size" == 8 ]]; then
    ilp64=ON
else
    ilp64=OFF
fi
"$cmake_command" -S dftd4 -B dftd4/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="${CC:-cc}" \
    -DCMAKE_Fortran_COMPILER="${FC:-gfortran}" \
    -DCMAKE_INSTALL_PREFIX="$prefix" \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DBUILD_SHARED_LIBS=OFF \
    -DWITH_API=ON \
    -DWITH_ILP64="$ilp64" \
    -DBLAS_LIBRARIES="$blas_libraries" \
    -DLAPACK_LIBRARIES="$lapack_libraries" \
    -DWITH_OpenMP=OFF \
    -Ddftd4-dependency-method=fetch
"$cmake_command" --build dftd4/build --parallel "${MAKE_JOBS:-4}"
"$cmake_command" --install dftd4/build
