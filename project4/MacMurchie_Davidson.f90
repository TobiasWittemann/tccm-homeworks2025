program MacMurchieDavidson
implicit none
integer :: Nbasis, MaxCont, i, j, t, t_max, k, l, i_max, j_max, d, m, n
double precision, allocatable, dimension(:,:) :: basis_centers, primitive_exponents, cartesian_exponents, contraction_coeffs
double precision, allocatable, dimension(:,:,:) ::  E_tij
double precision, dimension(3) :: rA, rB, rAB, rPA, rPB, rP
double precision :: p, a, b, mu, xAB, xPA, xPB, Nn, Nm, double_factorial, S_ij, dm, dn, S_prim, S_mm, S_nn
double precision,dimension(:,:), allocatable :: S
double precision, parameter :: pi = acos(-1.0d0)

! Read basis set data from provided files
Nbasis = 10
maxcont = 3

allocate(basis_centers(Nbasis,3), primitive_exponents(Nbasis, maxcont), cartesian_exponents(Nbasis,maxcont), contraction_coeffs(Nbasis,maxcont))

call read_array('basis_centers_sto3g.txt', Nbasis, 3, basis_centers)
call read_array('primitive_exponents_sto3g.txt', Nbasis, maxcont, primitive_exponents)
call read_array('cartesian_exponents_sto3g.txt', Nbasis, 3, cartesian_exponents)
call read_array('contraction_coeffs_sto3g.txt', Nbasis, maxcont, contraction_coeffs)

! Initialize overlap matrix
allocate(S(NBasis,NBasis))
S = 0.0d0

! Iterate over contracted Gaussians
do j=1,NBasis
  rB = basis_centers(j,:)
do i=1,j
  rA = basis_centers(i,:)
  rAB = rA - rB
  S_ij = 0.0d0
  ! Iterate over primitive Gaussians
  do n=1,maxcont
    dn = contraction_coeffs(j,n)
    if (dn == 0.0d0) exit
    b = primitive_exponents(j,n) 
  do m=1,maxcont
    dm = contraction_coeffs(i,m)
    if (dm == 0.0d0) exit
    a = primitive_exponents(i,m)
    p = a + b
    mu = a*b/p
    rP = (a*rA + b*rB)/p
    rPA = rP-rA
    rPB = rP-rB

  S_prim = dm*dn
  ! Iterate over x,y,z
  do d=1,3
    i_max = cartesian_exponents(i,d)
    j_max = cartesian_exponents(j,d)
    t_max = i_max + j_max
    xPA = rPA(d)
    xPB = rPB(d)
    xAB = rAB(d)
    allocate(E_tij(0:t_max, 0:i_max+1, 0:j_max+1))
    call calc_E_tij(i_max, j_max, p, mu, xAB, xPA, xPB, E_tij)     
    Nm = 1/SQRT(double_factorial(REAL(2*i_max-1,kind=8))/((4*a)**i_max)*SQRT(pi/(2*a)))
    Nn = 1/SQRT(double_factorial(REAL(2*j_max-1,kind=8))/((4*b)**j_max)*SQRT(pi/(2*b)))
    S_prim = S_prim*Nm*Nn*E_tij(0,i_max,j_max)*SQRT(pi/p)
    deallocate(E_tij)
  end do
  S_ij = S_ij + S_prim
  end do
  end do
  S(i,j) = S_ij
  S(j,i) = S_ij
end do
end do

call write_array(S,Nbasis,Nbasis)

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


subroutine calc_E_tij(i_max, j_max, p, mu, xAB, xPA, xPB, E_tij)
implicit none
integer, intent(in) :: i_max, j_max
double precision, intent(in) :: p, mu, xAB, xPA, xPB
double precision, dimension(0:i_max+j_max, 0:i_max+1, 0:j_max+1), intent(out) :: E_tij
integer :: i,j,t,t_max
E_tij = 0.0
E_tij(0,0,0) = EXP(-mu*xAB**2)
t_max = i_max + j_max
do i=0,i_max
  do j=0,j_max
    do t=0,i+j+1
      if (t > t_max) then
      else if (t == t_max) then
        E_tij(t,i+1,j) = 1/(2*p)*E_tij(t-1,i,j) + xPA*E_tij(t,i,j)
        E_tij(t,i,j+1) = 1/(2*p)*E_tij(t-1,i,j) + xPB*E_tij(t,i,j)
      else
        E_tij(t,i+1,j) = 1/(2*p)*E_tij(t-1,i,j) + xPA*E_tij(t,i,j) + (t+1)*E_tij(t+1,i,j)
        E_tij(t,i,j+1) = 1/(2*p)*E_tij(t-1,i,j) + xPB*E_tij(t,i,j) + (t+1)*E_tij(t+1,i,j)
      end if
    end do
  end do
end do
end subroutine calc_E_tij


function double_factorial(n_in) result(res_out)
implicit none
double precision, intent(in) :: n_in
integer :: res, i, n
double precision :: res_out
n = int(n_in)
if (n==0) then
  res = 1
else if ((n>0).and.(mod(n,2)==0)) then  
  res = 1
  do i=2,n,+2
   res = res*i
  end do
  res_out = real(res, kind=8)
else if ((n>0).and.(mod(n,2)==1)) then
  res = 1
  do i=1,n,2
   res = res*i
  end do
  res_out = real(res, kind=8)
else if ((n<0).and.(mod(n,2)==-1)) then
  res = 1
  do i=n+2,1,2
   res = res*i
  end do
   res_out = 1/real(res, kind=8)
end if
end function double_factorial


