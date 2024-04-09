file(REMOVE_RECURSE
  "libpaw.pdb"
  "libpaw.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang C Fortran)
  include(CMakeFiles/paw.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
