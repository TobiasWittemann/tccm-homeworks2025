program MacMurchieDavidson
implicit none
integer :: Nbasis, MaxCont, i
double precision, allocatable, dimension(:,:) :: basis_centers, primitive_exponents, cartesian_exponents, contraction_coeffs

Nbasis = 18
MaxCont = 6

allocate(basis_centers(Nbasis,3))

open(1, file='basis_centers.txt', status="old", action='read')
do i=1,Nbasis
   read(1,*) basis_centers(i,:)
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

