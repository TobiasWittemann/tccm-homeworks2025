program test_matrix_generation
    use BasisSet_Module
    implicit none

    integer :: unit_nb
    type(AtomicBasis) :: basis
    double precision, allocatable :: basis_centers(:,:), primitive_exponents(:,:), &
                                     cartesian_exponents(:,:), contraction_coeffs(:,:)
    integer :: Nbasis, maxcont, offset
    double precision :: pos_A(3)
    logical :: success

    success = .true.
    unit_nb = 10
    
    ! TEST CASE: one Carbon atom (Z=6)
    open(unit=unit_nb, file='6-31g.1.dalton', status='old', action='read')
    call ParseElement(unit_nb, 6, basis)
    close(unit_nb)

    Nbasis  = GetNBasis(basis)
    maxcont = GetMaxCont(basis)
    
    allocate(basis_centers(Nbasis, 3))
    allocate(primitive_exponents(Nbasis, maxcont))
    allocate(cartesian_exponents(Nbasis, 3))
    allocate(contraction_coeffs(Nbasis, maxcont))
    
    pos_A = [0.0d0, 0.0d0, 0.0d0]
    offset = 0
    
    call GenerateMDMatrices(basis, pos_A, Nbasis, maxcont, offset, &
                                basis_centers, primitive_exponents, &
                                cartesian_exponents, contraction_coeffs)
    
    ! Check final offset
    if (offset /= Nbasis) then
        print *, "[FAIL] Final offset should match Nbasis. Expected:", Nbasis, " Got:", offset
        success = .false.
    end if

    ! Check Cartesian exponents for p-orbitals (Rows 4-6)
    ! Expecting (1,0,0), (0,1,0), (0,0,1) for the first p-shell
    if (int(cartesian_exponents(4,1)) /= 1 .or. int(cartesian_exponents(5,2)) /= 1) then
        print *, "[FAIL] Cartesian exponents mapping is incorrect."
        success = .false.
    else
        print *, "[PASS] Cartesian mapping (l,m,n) verified."
    end if

    ! Check primitive padding
    ! Carbon p-orbitals only have 4 primitives. Columns 5-10 should be 0.
    if (abs(primitive_exponents(4, 5)) > 1.d-10) then
        print *, "[FAIL] Padding logic failed. Column 5 of p-orbital should be 0."
        success = .false.
    else
        print *, "[PASS] Primitive padding verified."
    end if

    if (success) then
        print *, "The matrices were generated correctly."
    else
        stop "RESULT: Generator Tests FAILED!"
    end if

end program test_matrix_generation