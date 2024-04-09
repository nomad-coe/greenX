file(REMOVE_RECURSE
  "libabinit_common.a"
  "libabinit_common.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C Fortran)
  include(CMakeFiles/abinit_common.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
