cmake_minimum_required(VERSION 3.25)
project(pico_image LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Bake the FFTW / conda env library paths into the binary so it runs without
# LD_LIBRARY_PATH on cluster nodes.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)
if(DEFINED ENV{CONDA_PREFIX})
  list(APPEND CMAKE_BUILD_RPATH   "$ENV{CONDA_PREFIX}/lib")
  list(APPEND CMAKE_INSTALL_RPATH "$ENV{CONDA_PREFIX}/lib")
endif()
list(APPEND CMAKE_BUILD_RPATH   "/usr/local/lib" "/usr/local/lib64")
list(APPEND CMAKE_INSTALL_RPATH "/usr/local/lib" "/usr/local/lib64")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# ---- Dependencies ------------------------------------------------------------

find_package(OpenMP REQUIRED)

# cfitsio: provided system-wide (libcfitsio.so) on Linux. Fall back to manual.
find_library(CFITSIO_LIB cfitsio REQUIRED)
find_path(CFITSIO_INCLUDE_DIR fitsio.h)

# FINUFFT and SuperNOVAS as sibling source trees (../finufft, ../SuperNOVAS).
# Allow overrides via -DFINUFFT_DIR=... / -DSUPERNOVAS_DIR=...
set(FINUFFT_DIR    "${CMAKE_CURRENT_SOURCE_DIR}/../finufft"     CACHE PATH "")
set(SUPERNOVAS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../SuperNOVAS"  CACHE PATH "")

# Disable optional finufft features so the subbuild is fast.
set(FINUFFT_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(FINUFFT_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
set(FINUFFT_BUILD_DEVEL    OFF CACHE BOOL "" FORCE)
add_subdirectory(${FINUFFT_DIR}    ${CMAKE_BINARY_DIR}/_finufft)

# SuperNOVAS: core C99 lib only (no ephemeris plugins, no C++ extras, no tests,
# no cURL — we don't fetch EOP at runtime; reduced-accuracy precession+nutation
# is what we use for J2000 rotation).
set(BUILD_TESTING   OFF CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES  OFF CACHE BOOL "" FORCE)
set(BUILD_BENCHMARK OFF CACHE BOOL "" FORCE)
set(WITHOUT_CURL    ON  CACHE BOOL "" FORCE)
add_subdirectory(${SUPERNOVAS_DIR} ${CMAKE_BINARY_DIR}/_supernovas)

# ---- pico_image binary -------------------------------------------------------

add_executable(pico_image
  src/config.cpp
  src/antsamp.cpp
  src/raw_io.cpp
  src/uvw.cpp
  src/dm.cpp
  src/calib.cpp
  src/grid.cpp
  src/clean.cpp
  src/fits_out.cpp
  src/pipeline.cpp
  src/main.cpp
)

target_include_directories(pico_image PRIVATE
  include
  ${CFITSIO_INCLUDE_DIR}
)

target_link_libraries(pico_image PRIVATE
  finufft::finufft
  supernovas::core
  ${CFITSIO_LIB}
  OpenMP::OpenMP_CXX
)

target_compile_options(pico_image PRIVATE -Wall -Wno-unused-parameter)

install(TARGETS pico_image RUNTIME DESTINATION bin)
