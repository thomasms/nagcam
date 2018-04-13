project (EXAMPLELIB)
enable_language (Fortran)

set (EXE_NAME nagcamexample)
set (NAGCAM_SOURCES solver_m.F90
                    types_m.F90
                    RHS_m.F90
                    CFL_m.F90
                    IO_m.F90)

set (EXE_SOURCE fd1d_heat_explicit.f90 ${NAGCAM_SOURCES})
set_source_files_properties(${EXE_SOURCE} PROPERTIES COMPILE_FLAGS "${F90_FLAGS}")

add_executable (${EXE_NAME} ${EXE_SOURCE})

set_target_properties(${EXE_NAME}
                      PROPERTIES OUTPUT_NAME ${EXE_NAME}
                      RUNTIME_OUTPUT_DIRECTORY ${PROJ_BUILD_DIR}/bin
                      LIBRARY_OUTPUT_DIRECTORY ${PROJ_BUILD_DIR}/lib
                      ARCHIVE_OUTPUT_DIRECTORY ${PROJ_BUILD_DIR}/lib
                      Fortran_MODULE_DIRECTORY ${PROJ_MODULE_DIR})

