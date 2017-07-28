# - Find the NETCDF library
#
# Usage:
#   find_package(NETCDF [REQUIRED] [QUIET] )
#     
# It sets the following variables:
#   NETCDF_FOUND               ... true if NETCDF is found on the system
#   NETCDF_LIBRARIES           ... full path to NETCDF library
#   NETCDF_INCLUDE_DIR         ... NETCDF include directory

if(NOT NETCDF_INCLUDE_DIR OR NOT NETCDF_LIBRARIES)
    
    #find libs
    find_library(
        NETCDF_LIB
        NAMES netcdf
        HINTS 
        ${NETCDF_ROOT}
        $ENV{NETCDF_ROOT}
        PATH_SUFFIXES lib
    )
    
    #find includes
    find_path(
        NETCDF_INCLUDE_DIR
        NAMES netcdf.h
        HINTS 
        ${NETCDF_ROOT}
        $ENV{NETCDF_ROOT}
        PATH_SUFFIXES include
    )
    
endif()

set(PARALLEL_NETCDF_FOUND TRUE)
if(USE_PARALLEL_NETCDF)
    if(NOT EXISTS ${NETCDF_INCLUDE_DIR}/netcdf_par.h)
        message(STATUS "Parallel NetCDF header netcdf_par.h is not found.")
        set(PARALLEL_NETCDF_FOUND FALSE)
    endif()
endif()

set(NETCDF_LIBRARIES ${NETCDF_LIB})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NETCDF DEFAULT_MSG
NETCDF_INCLUDE_DIR NETCDF_LIBRARIES PARALLEL_NETCDF_FOUND)

mark_as_advanced(NETCDF_INCLUDE_DIR NETCDF_LIBRARIES)

