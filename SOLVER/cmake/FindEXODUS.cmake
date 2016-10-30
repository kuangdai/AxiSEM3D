#
# Find the Exodus finite element data model library 
#
#  EXODUS_FOUND - System has Exodus
#  EXODUS_INCLUDE_DIR - The LibXml2 include directory
#  EXODUS_LIBRARIES - The libraries needed to use LibXml2

if (NOT EXODUS_INCLUDE_DIR)
    find_path(EXODUS_INCLUDE_DIR NAMES exodusII.h
        HINTS
        ${EXODUS_ROOT}
        $ENV{EXODUS_ROOT}
        PATH_SUFFIXES include
    )
endif ()

if (NOT EXODUS_LIBRARIES)
    find_library(EXODUS_LIBRARIES NAMES exodus
        HINTS
        ${EXODUS_ROOT}
        $ENV{EXODUS_ROOT}
        PATH_SUFFIXES lib
    )
endif ()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EXODUS DEFAULT_MSG
EXODUS_INCLUDE_DIR EXODUS_LIBRARIES)

mark_as_advanced(EXODUS_INCLUDE_DIR EXODUS_LIBRARIES)
