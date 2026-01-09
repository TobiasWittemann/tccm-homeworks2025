## Installation

### test
```bash
# compile module
gfortran -c src/basis_set_parser.f90 -J build/ -o build/basis_set_parser.o
# run basis set parsing test
gfortran tests/test_basis_parsing.f90 build/basis_set_parser.o -I build/ -o build/run_basis_test
./build/run_basis_test
# run matrix generation test
gfortran tests/test_matrix_generation.f90 build/basis_set_parser.o -I build/ -o build/run_matrix_test
./build/run_matrix_test
```

### Compilation
```bash
gfortran -c src/basis_set_parser.f90 -J build/ -o build/basis_set_parser.o
gfortran src/MacMurchie_Davidson.f90 build/basis_set_parser.o -I build/ -o build/MacMurchie_Davidson
./build/MacMurchie_Davidson
gfortran src/Spherical2Cartesian.f90 -o build/Spherical2Cartesian
./build/Spherical2Cartesian 
```
