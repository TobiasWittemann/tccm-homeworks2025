module BasisSet_Module
    implicit none
    private
    public :: Primitive, Shell, AtomicBasis, ParseElement, PrintSphericalBasis, &
              GenerateMDMatrices, GetNBasis, GetMaxCont

    ! Represents a single Gaussian primitive
    type Primitive
        real(8) :: alpha              ! Gaussian exponent
        real(8), allocatable :: d(:)   ! Contraction coefficients for different AOs
    end type Primitive

    ! Represents a shell
    type Shell
        integer :: L                   ! Angular momentum: 0=s, 1=p, 2=d
        integer :: n_prim              ! Number of primitives in this shell
        integer :: n_ao                ! Number of contracted atomic orbitals (AOs)
        type(Primitive), allocatable :: prims(:)
    end type Shell

    ! Represents the full basis set for a specific atom
    type AtomicBasis
        integer :: Z                   ! Atomic number
        type(Shell), allocatable :: shells(:)
    end type AtomicBasis

contains

    ! Calculate primitive normalization constant
    function calc_prim_norm(alpha, L) result(N)
        real(8), intent(in) :: alpha
        integer, intent(in) :: L
        real(8) :: N
        real(8), parameter :: PI = acos(-1.0d0)
        integer :: i
        real(8) :: dfact

        dfact = 1.0d0
        if (L > 0) then
            ! Compute double factorial (2L-1)!!
            do i = 2*L-1, 1, -2
                dfact = dfact * i
            end do
        end if

        N = ((2.0d0 * alpha / PI)**0.75d0) * &
            (sqrt((4.0d0 * alpha)**L / dfact))
    end function calc_prim_norm

    ! Calculate contracted normalization factor
    function get_contracted_norm(sh, ao_index, L_value) result(N_tot)
        type(Shell), intent(in) :: sh
        integer, intent(in) :: ao_index
        integer, intent(in) :: L_value
        real(8) :: N_tot
        
        integer :: i, j
        real(8) :: c_i, c_j, alpha_i, alpha_j, overlap, sum_overlap 
        real(8) :: power_factor

        sum_overlap = 0.0d0
        ! Overlap power factor: L + 1.5
        power_factor = real(L_value, 8) + 1.5d0

        do i = 1, sh%n_prim
            do j = 1, sh%n_prim
                alpha_i = sh%prims(i)%alpha
                alpha_j = sh%prims(j)%alpha
                c_i = sh%prims(i)%d(ao_index) 
                c_j = sh%prims(j)%d(ao_index) 

                ! Overlap S_ij between two normalized primitives
                overlap = ( (2.0d0 * sqrt(alpha_i * alpha_j)) / (alpha_i + alpha_j) )**power_factor
                
                sum_overlap = sum_overlap + (c_i * c_j * overlap)
            end do
        end do
        
        if (sum_overlap < 1.0d-20) then
             N_tot = 0.0d0 ! Prevent division by zero
        else
             N_tot = 1.0d0 / sqrt(sum_overlap)
        end if
    end function get_contracted_norm

    ! Parser for dalton basis files
    subroutine ParseElement(unit_nb, target_Z, basis)
        integer, intent(in) :: unit_nb
        integer, intent(in) :: target_Z
        type(AtomicBasis), intent(out) :: basis
        character(len=256) :: line
        integer :: ios, current_Z, s_idx, i, n_prim, n_ao, n_shells
        integer :: current_L  ! Track L based on comments 
        character(len=2) :: sym

        rewind(unit_nb)
        basis%Z = target_Z

        ! Locate target atom and count shells 
        find_z_1: do
            read(unit_nb, '(A)', iostat=ios) line
            if (ios < 0) stop "Error: Atomic number not found."
            if (line(1:1) == 'a' .or. line(1:1) == 'A') then
                read(line(2:), *, iostat=ios) current_Z
                if (ios == 0 .and. current_Z == target_Z) exit find_z_1
            end if
        end do find_z_1

        n_shells = 0
        do 
            read(unit_nb, '(A)', iostat=ios) line
            if (ios < 0) exit

            if ((line(1:1) == 'a' .or. line(1:1) == 'A') .and. &
                len_trim(line) > 1) then
                exit 
            end if

            if (len_trim(line) > 0 .and. line(1:1) /= '!' .and. line(1:1) /= ' ') then
                n_shells = n_shells + 1 
                read(line, *) sym, n_prim, n_ao
                do i = 1, n_prim; read(unit_nb, *); end do
            end if
        end do
        
        allocate(basis%shells(n_shells))
        
        ! Read detailed data and assign L
        rewind(unit_nb) 
        find_z_2: do
            read(unit_nb, '(A)', iostat=ios) line
            if (line(1:1) == 'a' .or. line(1:1) == 'A') then
                read(line(2:), *, iostat=ios) current_Z
                if (ios == 0 .and. current_Z == target_Z) exit find_z_2
            end if
        end do find_z_2

        current_L = 0  ! Default to s-type
        s_idx = 0
        
        do while (s_idx < n_shells)
        read(unit_nb, '(A)', iostat=ios) line
        if (ios /= 0) exit

        if ((line(1:1) == 'a' .or. line(1:1) == 'A') .and. len_trim(line) > 1) then
             exit 
        end if
            
            ! Check comment lines for L-type keywords
            if (line(1:1) == '!') then
                if (index(line, 's functions') > 0) current_L = 0
                if (index(line, 'p functions') > 0) current_L = 1
                if (index(line, 'd functions') > 0) current_L = 2
                cycle ! Skip comment line
            end if
            
            ! Read data headers (starting with 'H')
            if (len_trim(line) > 0 .and. line(1:1) /= ' ') then
                s_idx = s_idx + 1
                read(line, *) sym, n_prim, n_ao
                
                basis%shells(s_idx)%L = current_L
                basis%shells(s_idx)%n_prim = n_prim
                basis%shells(s_idx)%n_ao = n_ao
                
                allocate(basis%shells(s_idx)%prims(n_prim))
                
                do i = 1, n_prim 
                    allocate(basis%shells(s_idx)%prims(i)%d(n_ao))
                    read(unit_nb, *) basis%shells(s_idx)%prims(i)%alpha, &
                                     basis%shells(s_idx)%prims(i)%d(1:n_ao)
                end do 
            end if
        end do
    end subroutine ParseElement

    ! Print basis info in spherical notation
    subroutine PrintSphericalBasis(basis)
        type(AtomicBasis), intent(in) :: basis
        integer :: s, a, p, m, current_L
        real(8) :: alpha, raw_coeff
        character(len=10) :: ao_name

        do s = 1, size(basis%shells)
            current_L = basis%shells(s)%L

            do a = 1, basis%shells(s)%n_ao
                ! Generate AO name (1s, 2s, 1p, etc.) 
                if (current_L == 0) then
                    write(ao_name, '(I1, "s")') a
                else
                    write(ao_name, '(I1, "p")') a
                end if

                ! Print s-type (L=0)
                if (current_L == 0) then
                    write(*, '(/,"chi ", A, " =")') trim(ao_name)
                    do p = 1, basis%shells(s)%n_prim
                        alpha = basis%shells(s)%prims(p)%alpha
                        raw_coeff = basis%shells(s)%prims(p)%d(a)
                        if (abs(raw_coeff) > 1.d-10) then
                            write(*, '(F15.8, " * exp(", F15.8, " r^2)")') raw_coeff, -alpha
                        end if 
                    end do
                else
                    ! Print p-type (m = -1, 0, 1)
                    do m = -1, 1
                        write(*, '(/,"chi ", A, " m=", I2, " =")') trim(ao_name), m
                        do p = 1, basis%shells(s)%n_prim
                            alpha = basis%shells(s)%prims(p)%alpha
                            raw_coeff = basis%shells(s)%prims(p)%d(a)
                            if (abs(raw_coeff) > 1.d-10) then
                                write(*, '(F15.8, " * S", I1, I2, "(r) * exp(", F15.8, " r^2)")') &
                                raw_coeff, current_L, m, -alpha 
                            end if
                        end do
                    end do
                end if
            end do
        end do
    end subroutine PrintSphericalBasis

    ! Flatten nested basis set into MD matrices
    subroutine GenerateMDMatrices(basis, atom_pos, n_basis, max_prim, &
                                 offset, basis_centers, primitive_exponents, &
                                 cartesian_exponents, contraction_coeffs)
        type(AtomicBasis), intent(in) :: basis
        real(8), intent(in)           :: atom_pos(3)
        integer, intent(in)           :: n_basis
        integer, intent(in)           :: max_prim
        integer, intent(inout)        :: offset
        
        double precision, intent(inout) :: basis_centers(n_basis, 3)
        double precision, intent(inout) :: primitive_exponents(n_basis, max_prim)
        double precision, intent(inout) :: cartesian_exponents(n_basis, 3)
        double precision, intent(inout) :: contraction_coeffs(n_basis, max_prim)

        integer :: s, a, p, L_val, b_idx, m_count, n_prim_in_shell
        
        b_idx = offset 

        do s = 1, size(basis%shells)
            L_val = basis%shells(s)%L
            n_prim_in_shell = basis%shells(s)%n_prim 
            
            ! Loop over each Shell
            do a = 1, basis%shells(s)%n_ao
                
                ! Loop over each spherical harmonic component of the contracted orbital (2L + 1)
                ! S orbitals loop 1 time; P orbitals loop 3 times
                do m_count = 1, (2 * L_val + 1)
                    
                    b_idx = b_idx + 1
                    
                    ! atom coordinates for the current basis function
                    basis_centers(b_idx, :) = atom_pos
                    
                    ! Cartesian exponents (l, m, n)
                    ! for s and p only
                    if (L_val == 0) then
                        cartesian_exponents(b_idx, :) = [0.0d0, 0.0d0, 0.0d0]
                    else if (L_val == 1) then
                        if (m_count == 1) cartesian_exponents(b_idx, :) = [1.0d0, 0.0d0, 0.0d0] ! px
                        if (m_count == 2) cartesian_exponents(b_idx, :) = [0.0d0, 1.0d0, 0.0d0] ! py
                        if (m_count == 3) cartesian_exponents(b_idx, :) = [0.0d0, 0.0d0, 1.0d0] ! pz
                    end if

                    ! Populate primitive exponents and contraction coefficients
                    do p = 1, n_prim_in_shell
                        primitive_exponents(b_idx, p) = basis%shells(s)%prims(p)%alpha
                        
                        contraction_coeffs(b_idx, p) = basis%shells(s)%prims(p)%d(a)
                    end do
                    
                end do
            end do
        end do
        
        offset = b_idx 
    end subroutine GenerateMDMatrices

    ! Calculate total Cartesian basis function count (Nbasis)
    function GetNBasis(basis) result(n_total)
        type(AtomicBasis), intent(in) :: basis
        integer :: n_total, s, a, L_val
        
        n_total = 0
        do s = 1, size(basis%shells)
            L_val = basis%shells(s)%L
            ! Each AO contributes (L+1)(L+2)/2 Cartesian components
            do a = 1, basis%shells(s)%n_ao
                n_total = n_total + (L_val + 1) * (L_val + 2) / 2
            end do
        end do
    end function GetNBasis

    ! Calculate maximum contractions (maxcont)
    function GetMaxCont(basis) result(maxcont)
        type(AtomicBasis), intent(in) :: basis
        integer :: maxcont, s
        
        maxcont = 0
        do s = 1, size(basis%shells)
            ! debug
            ! print *, "Shell index:", s, " L=", basis%shells(s)%L, " n_prim=", basis%shells(s)%n_prim
            if (basis%shells(s)%n_prim > maxcont) then
                maxcont = basis%shells(s)%n_prim
            end if
        end do
    end function GetMaxCont

end module BasisSet_Module
