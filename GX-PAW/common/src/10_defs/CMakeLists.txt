# Obtain files under the folder
file(GLOB ABINIT_COMMON_SRC_FOLDER
    "*.F90"
    "*.f90"
    "*.c"
    "*.cpp")

set(ABINIT_COMMON_SRC ${ABINIT_COMMON_SRC} ${ABINIT_COMMON_SRC_FOLDER} PARENT_SCOPE)
