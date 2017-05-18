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

set(NETCDF_LIBRARIES ${NETCDF_LIB})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NETCDF DEFAULT_MSG
NETCDF_INCLUDE_DIR NETCDF_LIBRARIES)

mark_as_advanced(NETCDF_INCLUDE_DIR NETCDF_LIBRARIES)

