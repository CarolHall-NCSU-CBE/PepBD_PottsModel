ifort -c -O2 constant.f90
ifort -c -O2 input.f90
ifort -c -O2 utilities.f90
ifort -c -O2 pdbfile.f90
ifort -c -O2 database.f90
ifort -c -O2 surface_area.f90
ifort -c -O2 energy_calculation.f90
ifort -c -O2 advanced_function.f90
ifort -c -O2 PepBD_Potts.f90
ifort -o -O2 PepBD_Potts *.o
rm *.o *.mod
