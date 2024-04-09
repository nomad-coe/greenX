# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/mrm/progs/greenX/GX-PAW"
  "/home/mrm/progs/greenX/GX-PAW/src/LIBPAW-build"
  "/home/mrm/progs/greenX/GX-PAW"
  "/home/mrm/progs/greenX/GX-PAW/tmp"
  "/home/mrm/progs/greenX/GX-PAW/src/LIBPAW-stamp"
  "/home/mrm/progs/greenX/GX-PAW/src"
  "/home/mrm/progs/greenX/GX-PAW/src/LIBPAW-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/mrm/progs/greenX/GX-PAW/src/LIBPAW-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/mrm/progs/greenX/GX-PAW/src/LIBPAW-stamp${cfgdir}") # cfgdir has leading slash
endif()
