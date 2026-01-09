! print the spherical basis expressions for a given atomic number
program spherical_basis_expressions
    use BasisSet_Module
    implicit none

    type(AtomicBasis) :: basis_C
    integer :: unit_idx, nb, mc, target_Z
    character(len=20) :: input_buffer
    integer :: ios

    do
        write(*, '(A)', advance='no') "Enter the Atomic Number (Z) for the basis (1-36): "
        read(*, '(A)') input_buffer
        
        ! Try to read an integer from the string buffer
        read(input_buffer, *, iostat=ios) target_Z
        
        if (ios /= 0) then
            write(*, *) "Error: Invalid format. Please enter a whole integer."
            cycle ! Restart the loop
        end if
        
        if (target_Z < 1 .or. target_Z > 36) then
            write(*, *) "Error: Atomic number must be between 1 and 36."
            cycle ! Restart the loop
        end if
        exit
    end do

    open(newunit=unit_idx, file='tests/6-31g.1.dalton', status='old', action='read')
    call ParseElement(unit_idx, target_Z, basis_C)
    close(unit_idx)
    
    nb = GetNBasis(basis_C)
    mc = GetMaxCont(basis_C)
    
    print *, "Selected atom Z =", target_Z
    print *, "Calculated Nbasis:", nb
    print *, ""
    print *, "Calculated MaxCont:", mc
    
    print *, "Spherical Basis Expressions:"
    call PrintSphericalBasis(basis_C)

end program spherical_basis_expressions
