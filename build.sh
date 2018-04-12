#!bin/bash

rm *.mod *.o test h_test01.txt
nagfor src/types_m.F90 src/RHS_m.F90 src/CFL_m.F90 src/IO_m.F90 src/solver_m.F90 src/fd1d_heat_explicit.f90 -o test
./test
diff ref/h_test01.txt_bak h_test01.txt
