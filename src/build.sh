#!bin/bash

rm *.mod *.o test h_test01.txt
nagfor types_m.F90 fd1d_heat_explicit.f90 -o test
./test
diff h_test01.txt_bak h_test01.txt
