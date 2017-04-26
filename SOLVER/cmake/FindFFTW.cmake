# - Find the FFTW library
#
# Usage:
#   find_package(FFTW [REQUIRED] [QUIET] )
#     
# It sets the following variables:
#   FFTW_FOUND               ... true if fftw is found on the system
#   FFTW_LIBRARIES           ... full path to fftw library
#   FFTW_INCLUDE_DIR         ... fftw include directory

if(NOT FFTW_INCLUDE_DIR OR NOT FFTW_LIBRARIES)
    
    #find libs
    find_library(
        FFTW_LIB
        NAMES fftw3
        HINTS 
        ${FFTW_ROOT}
        $ENV{FFTW_ROOT}
        PATH_SUFFIXES lib
    )
    
    if(NOT FFTW_LIB)
        message(STATUS "Double-precision FFTW library is not found. Install FFTW both with and without --enable-float.")
    endif(NOT FFTW_LIB)

    find_library(
        FFTWF_LIB
        NAMES fftw3f
        HINTS 
        ${FFTW_ROOT}
        $ENV{FFTW_ROOT}
        PATH_SUFFIXES lib
    )
    
    if(NOT FFTWF_LIB)
        message(STATUS "Single-precision FFTW library is not found. Install FFTW both with and without --enable-float.")
    endif(NOT FFTWF_LIB)

    #find includes
    find_path(
        FFTW_INCLUDE_DIR
        NAMES fftw3.h
        HINTS 
        ${FFTW_ROOT}
        $ENV{FFTW_ROOT}
        PATH_SUFFIXES include
    )
    
endif()

set(FFTW_LIBRARIES ${FFTW_LIB} ${FFTWF_LIB})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG
FFTW_INCLUDE_DIR FFTW_LIBRARIES)

mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARIES FFTW_LIB FFTWF_LIB)

