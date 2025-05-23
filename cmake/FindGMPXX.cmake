# - Find GMPXX library
# Once done, this will define
#  GMPXX_FOUND
#  GMPXX_INCLUDE_DIR
#  GMPXX_LIBRARIES
include(FindPackageHandleStandardArgs)
find_package(PkgConfig)

if (PKG_CONFIG_FOUND)
  pkg_check_modules(GREENX_GMPXX IMPORTED_TARGET GLOBAL gmpxx)
  pkg_check_modules(GREENX_GMP IMPORTED_TARGET GLOBAL gmp)
endif()
if(NOT GREENX_GMP_FOUND)
  find_path(GREENX_GMPXX_INCLUDE_DIR NAMES gmpxx.h)
  find_library(GREENX_GMPXX_LIBRARY NAMES gmpxx)
  find_library(GREENX_GMP_LIBRARY NAMES gmp)
else()
  set(GREENX_GMP_LIBRARY ${GREENX_GMP_LINK_LIBRARIES})
  set(GREENX_GMPXX_LIBRARY ${GREENX_GMPXX_LINK_LIBRARIES})
  set(GREENX_GMPXX_INCLUDE_DIR ${GREENX_GMP_INCLUDE_DIRS})
endif()

find_package_handle_standard_args(GMPXX DEFAULT_MSG GREENX_GMPXX_INCLUDE_DIRS GREENX_GMPXX_LIBRARY GREENX_GMP_LIBRARY)

set(GREENX_GMPXX_LIBRARIES ${GREENX_GMPXX_LIBRARY} ${GREENX_GMP_LIBRARY})

if (NOT TARGET greenX::gmpxx)
  add_library(greenX::gmpxx INTERFACE IMPORTED)
  set_target_properties(greenX::gmpxx
                        PROPERTIES
                        INTERFACE_LINK_LIBRARIES "${GREENX_GMPXX_LIBRARIES}")
  if (GREENX_GMPXX_INCLUDE_DIRS)
    set_target_properties(greenX::gmpxx
                          PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES ${GREENX_GMPXX_INCLUDE_DIRS})
  endif()
endif()

mark_as_advanced(GMPXX_INCLUDE_DIR GMPXX_LIBRARY)
