# - Try to find METIS
# Once done this will define
#
#  METIS_FOUND        - system has METIS
#  METIS_INCLUDE_DIR  - include directories for METIS
#  METIS_LIBRARIES    - libraries for METIS

# check metis.h
macro(_metis_check_32bit)
    SET(METIS_32BIT TRUE)
    file(STRINGS "${METIS_INCLUDE_DIR}/metis.h" metis_file)
    foreach(line ${metis_file})
        string(FIND "${line}" "#define" pos_def)
        if (${pos_def} EQUAL "-1")
            continue()
        endif()
        
        string(FIND "${line}" "64" pos_val)
        if (${pos_val} EQUAL "-1")
            continue()
        endif()
        
        string(FIND "${line}" "IDXTYPEWIDTH" pos_key)
        if (NOT ${pos_key} EQUAL "-1") 
            SET(METIS_32BIT FALSE)
            break()
        endif()
        
        string(FIND "${line}" "REALTYPEWIDTH" pos_key)
        if (NOT ${pos_key} EQUAL "-1") 
            SET(METIS_32BIT FALSE)
            break()
        endif()
    endforeach()
    
    if (NOT METIS_32BIT) 
        message(STATUS "Incompatible METIS build found: ${METIS_INCLUDE_DIR}")
        message(STATUS "Please re-install METIS with 32-bit configuration, and/or edit "
            "METIS_ROOT in CMakeLists.txt to locate the 32-bit build.")
    endif()
endmacro(_metis_check_32bit)

# lib
if (NOT METIS_LIBRARIES)
    find_library(METIS_LIBRARIES
        NAMES metis
        HINTS 
        ${METIS_ROOT}
        $ENV{METIS_ROOT}
        PATH_SUFFIXES lib
    )
endif()    

# include
if(METIS_INCLUDE_DIR)
    _metis_check_32bit()
    set(METIS_FOUND ${METIS_32BIT})
else()
    find_path(METIS_INCLUDE_DIR 
        NAMES metis.h
        HINTS
        ${METIS_ROOT}
        $ENV{METIS_ROOT}
        PATH_SUFFIXES include)
    
    if (METIS_INCLUDE_DIR)
        _metis_check_32bit()
    endif()
endif()

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS DEFAULT_MSG
METIS_INCLUDE_DIR METIS_LIBRARIES METIS_32BIT)

mark_as_advanced(METIS_INCLUDE_DIR METIS_LIBRARIES)



