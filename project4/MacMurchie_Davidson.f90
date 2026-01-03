program MacMurchieDavidson
implicit none
integer :: Nbasis, MaxCont, i, j, t, t_max, k, l, i_max, j_max, d
double precision, allocatable, dimension(:,:) :: basis_centers, primitive_exponents, cartesian_exponents, contraction_coeffs
double precision, allocatable, dimension(:,:,:) ::  E_tij
double precision, dimension(3) :: rA, rB, rAB, rPA, rPB, rP
double precision :: p, a, b, mu


! Read basis set data from provided files
Nbasis = 18
maxcont = 6

allocate(basis_centers(Nbasis,3), primitive_exponents(Nbasis, maxcont), cartesian_exponents(Nbasis,maxcont), contraction_coeffs(Nbasis,maxcont))

call read_array('basis_centers.txt', Nbasis, 3, basis_centers)
call read_array('primitive_exponents.txt', Nbasis, maxcont, primitive_exponents)
call read_array('cartesian_exponents.txt', Nbasis, 3, cartesian_exponents)
call read_array('contraction_coeffs.txt', Nbasis, maxcont, contraction_coeffs)

!write(*,*) "basis centers"
!call write_array(basis_centers, Nbasis, 3)
!write(*,*) "primitive exponets"
!call write_array(primitive_exponents, Nbasis, maxcont)
!write(*,*) "cartesian exponents"
!call write_array(cartesian_exponents, Nbasis, 3)
!write(*,*) "contraction coeffs"
!call write_array(contraction_coeffs, Nbasis, maxcont)

! Use 2pz orbitals of C and O to test Hermitian recurrence relations

! Apply Gaussian Product theorem
a = primitive_exponents(9,1)
b = primitive_exponents(18,1)
p = a + b
mu = a*b/p
rA = basis_centers(9,:)
rB = basis_centers(18,:)
rP = (a*rA + b*rB)/p
rAB = rA - rB
rPA = rP-rA
rPB = rP-rB
! Consider z direction for MacMurchie-Davidson
d = 3

i_max = cartesian_exponents(9,d)
j_max = cartesian_exponents(18,d)
t_max = i_max + j_max

write(*,*) "i_max, j_max, t_max", i_max, j_max, t_max

allocate(E_tij(t_max+1, i_max+2, j_max+2))
E_tij = 0.0

! Calculate MacMurchie-Davidson expansion coefficients
do i=0,i_max
write(*,*) "i", i
  do j=0,j_max
    write(*,*) "j", j
    do t=0,i+j+1
      write(*,*) "t", t
      if (t > t_max) then
        write(*,*) "Case 1"
      else if (t == t_max) then
        write(*,*) "Case 2"
        E_tij(t,i+1,j) = 1/(2*p)*E_tij(t-1,i,j) + rPA(d)*E_tij(t,i,j)
        E_tij(t,i,j+1) = 1/(2*p)*E_tij(t-1,i,j) + rPB(d)*E_tij(t,i,j) 
      else
        write(*,*) "Case 3"
        E_tij(t,i+1,j) = 1/(2*p)*E_tij(t-1,i,j) + rPA(d)*E_tij(t,i,j) + (t+1)*E_tij(t+1,i,j)
        E_tij(t,i+1,j) = 1/(2*p)*E_tij(t-1,i,j) + rPB(d)*E_tij(t,i,j) + (t+1)*E_tij(t+1,i,j) 
      end if
    end do
  end do
end do

end program MacMurchieDavidson



subroutine write_array(arr, m, n)
integer, intent(in) :: m, n
double precision, intent(in), dimension(m,n) :: arr
integer :: i
do i=1,m
  write(*,*) arr(i,:)
end do
end subroutine write_array



subroutine read_array(file, m, n, arr)
implicit none
character(len=*) :: file
integer, intent(in) :: m, n
double precision, intent(out), dimension(m,n) :: arr
integer :: i

open(1, file=file, status="old", action='read')
do i=1,m
   read(1,*) arr(i,:)
end do
end subroutine read_array
