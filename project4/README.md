## Usage
# MacMurchie-Davidson
This program calculates the overlap matrix for an arbitrary molecule in a cartesian Gaussian basis set.
Using the MacMurchie-Davidson recursion relations, it allows the treatment of basis functions with arbitrary angular momentum.
As inputs, the geometry of the molecule has to be provided as .xyz file and the basis set needs to be provided in .dalton format.

# Spherical2Cartesian
This program calculates the transformation matrix from Spherical GTOs 


## Directory structure
```text
project4/
├── build/                 # Compiled binaries and intermediate files
│   ├── basis_set_parser.o # Object file for the basis set module
│   ├── basisset_module.mod# Fortran module definition file
│   ├── run_basis_test     # Executable for testing basis set parsing
│   └── run_matrix_test    # Executable for testing matrix generation
├── src/                   # Core source code
│   └── basis_set_parser.f90  # Main module containing basis set parsing logic
└── tests/                 # Unit and integration tests
    ├── test_basis_parsing.f90      # Tests for parsing logic
    └── test_matrix_generation.f90  # Tests for generating intermediate matricies
