include(CMakeFindDependencyMacro)

# If greenX depends on other libraries, do:
find_dependency(BLAS REQUIRED)
find_dependency(LAPACK REQUIRED)
# Finally, pull in the targets you exported:
include("${CMAKE_CURRENT_LIST_DIR}/greenXTargets.cmake")
