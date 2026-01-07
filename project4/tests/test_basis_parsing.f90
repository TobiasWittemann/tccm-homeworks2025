! Test program for basis set parsing get nbasis
program BasisSet_Test
    use BasisSet_Module
    implicit none

    integer :: unit_nb
    type(AtomicBasis) :: basis
    integer :: n_basis
    logical :: success

    success = .true.
    unit_nb = 10

    ! Test 1 Hydrogen (Z=1): expect 2 basis functions
    open(unit=unit_nb, file='6-31g.1.dalton', status='old', action='read')
    call ParseElement(unit_nb, 1, basis)
    
    n_basis = GetNBasis(basis)
    if (n_basis /= 2) then
        print *, "FAILED: Hydrogen (Z=1) Nbasis should be 2, got", n_basis
        print *, "       >>> check s orbital parsing (L=0)"
    else
        print *, "PASSED: Hydrogen (Z=1) Nbasis is 2"
    end if
    close(unit_nb)

    ! Test 2 Carbon (Z=6): expect 9 basis functions
    open(unit=unit_nb, file='6-31g.1.dalton', status='old', action='read')
    call ParseElement(unit_nb, 6, basis)
    n_basis = GetNBasis(basis)
    if (n_basis /= 9) then
        print *, "FAILED: Carbon (Z=6) Nbasis should be 9, got", n_basis
        print *, "       >>> check p orbital parsing (L=1)"
        success = .false.
    else
        print *, "PASSED: Carbon (Z=6) Nbasis is 9"
    end if
    close(unit_nb)

    ! Test 3 Zinc (Z=30): expect 28 basis functions
    open(unit=unit_nb, file='6-31g.1.dalton', status='old', action='read')
    call ParseElement(unit_nb, 30, basis)
    
    n_basis = GetNBasis(basis)
    if (n_basis /= 29) then
        print *, "FAILED: Zinc (Z=30) Nbasis should be 29, got", n_basis
        print *, "       >>> check d orbital parsing (L=2)"
        success = .false.
    else
        print *, "PASSED: Zinc (Z=30) Nbasis is 29"
    end if
    close(unit_nb)

    if (success) then
        print *, "This parser works for atoms with s, p, and d orbitals."
    else
        stop "Some tests failed!"
    end if

end program BasisSet_Test