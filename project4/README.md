## Usage
# 3.1 SphericalBasisExpressions
This program prints the mathematical expressions for the basis functions of a specified element. As inputs, the program takes the atomic number and the provided file 6-31G.1.dalton. 

# 3.2 Spherical2Cartesian
This program calculates the transformation matrix from Spherical GTOs to Cartesian GTOs as well as the self-overlap of the
spherical / cartesian GTOs provided the angular momentum l and the exponent a of the primitive Gaussian as inputs.

# 3.3 MacMurchieDavidson
This program calculates the overlap matrix for an arbitrary molecule in a cartesian Gaussian basis set.
Using the MacMurchie-Davidson recursion relations, it allows the treatment of basis functions with arbitrary angular momentum.
As inputs, the geometry of the molecule has to be provided as .xyz file and the basis set needs to be provided in .dalton format.

#  Testing of the programs
Tests for all programs are providedin the /tests/ directory. The .f90 can be compiled according to the instructions given in INSTALL.md. The .sh files can
be exectued to test programs 3.2 and 3.3 after compilation.


## Directory structure
```text
project4
├── bulid/
│   ├── basis_set_parser.o
│   ├── basisset_module.mod
│   ├── run_basis_test
│   ├── run_matrix_test
│   ├── spherical_basis_expressions
│   ├── MacMurchie_Davidson
│   └── Spherical2Cartesian
├── src/
│   ├── basis_set_parser.f90
│   ├── MacMurchie_Davidson.f90
│   └── Spherical2Cartesian.f90
└── tests/
    ├── 6-31g.1.dalton              # basis set file
    ├── CO.xyz                      # example geometry
    ├── orca_test_co_6-31G.out      # orca output for comparison
    ├── test_basis_parsing_C.f90.   # test programs
    ├── test_basis_parsing.f90
    ├── test_macmurchie_davidson.sh
    ├── test_matrix_generation.f90
    └── test_spherical2cartesian.sh 

