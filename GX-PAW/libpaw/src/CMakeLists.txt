# Because of the MACROS setting includes and use we need to set
# the ordering by hand.
add_compile_definitions(HAVE_LIBPAW_ABINIT)
set(SRCLIBPAWCORE libpaw_libxc.c
                  m_libpaw_mpi.F90
                  m_libpaw_tools.F90
                  m_libpaw_libxc.F90
                  m_paral_atom.F90
                  m_paw_numeric.F90
                  m_paw_sphharm.F90
                  m_paw_io.F90
                  m_pawang.F90
                  m_pawrad.F90
                  m_pawtab.F90
                  m_paw_an.F90
                  m_paw_ij.F90
                  m_pawfgrtab.F90
                  m_paw_finegrid.F90
                  m_paw_onsite.F90
                  m_pawcprj.F90
                  m_pawrhoij.F90
                  m_pawdij.F90
                  m_pawxc.F90
                  m_paw_atom.F90
                  m_paw_gaussfit.F90
                  m_pawxmlps.F90
                  m_pawpsp.F90)

# Files under paw
add_library(libpaw_objects OBJECT ${SRCLIBPAWCORE})
target_include_directories(libpaw_objects PUBLIC ${PROJECT_SOURCE_DIR}/common/src/incs/)
target_include_directories(libpaw_objects PUBLIC ${MOD_DIR})
add_library(paw SHARED $<TARGET_OBJECTS:libpaw_objects>)
target_link_libraries(paw PUBLIC abinit_common)

