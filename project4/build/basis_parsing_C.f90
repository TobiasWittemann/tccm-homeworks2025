! Test program for parsing basis set for Carbon (Z=6)
program test_basis
    use BasisSet_Module
    implicit none

    type(AtomicBasis) :: basis_C
    integer :: unit_idx, nb, mc, target_Z
    
    target_Z = 6 ! Carbon
    
    open(newunit=unit_idx, file='6-31g.1.dalton', status='old', action='read')
    call ParseElement(unit_idx, target_Z, basis_C)
    close(unit_idx)
    
    nb = GetNBasis(basis_C)
    mc = GetMaxCont(basis_C)
    
    print *, "--- Basis Set Test for Carbon (Z=6) ---"
    print *, "Expected Nbasis: 9"
    print *, "Calculated Nbasis:", nb
    print *, ""
    print *, "Expected MaxCont: 6"
    print *, "Calculated MaxCont:", mc
    print *, "---------------------------------------"
    
    print *, "Checking Spherical Basis Expressions:"
    call PrintSphericalBasis(basis_C)

end program test_basis
