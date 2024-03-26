# - Find GMPXX library
# Once done, this will define
#  GMPXX_FOUND
#  GMPXX_INCLUDE_DIR
#  GMPXX_LIBRARIES

find_path(GMPXX_INCLUDE_DIR NAMES gmpxx.h)
find_library(GMPXX_LIBRARY NAMES gmpxx)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMPXX DEFAULT_MSG GMPXX_INCLUDE_DIR GMPXX_LIBRARY)

if(GMPXX_FOUND)
  set(GMPXX_INCLUDE_DIRS ${GMPXX_INCLUDE_DIR})
  set(GMPXX_LIBRARIES ${GMPXX_LIBRARY})
else()
  set(GMPXX_INCLUDE_DIRS)
  set(GMPXX_LIBRARIES)
endif()

mark_as_advanced(GMPXX_INCLUDE_DIR GMPXX_LIBRARY)
