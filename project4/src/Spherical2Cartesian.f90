program main
implicit none
integer :: l, m, x_exp, y_exp, z_exp, i, j
integer, dimension(:,:), allocatable :: index_list
double precision, dimension(:,:), allocatable :: M_trans
double precision :: a, S_self, S_self_cart
integer,dimension(3) :: exponents
character(len=8) :: x_exp_char, y_exp_char, z_exp_char

! Get angular momentum input from user
write(*,*) "Please select angular momentum l"
read(*,*) l


! Calculate transformation matrix
allocate(index_list((l+1)*(l+2)/2,3))
index_list = 0
allocate(M_trans((l+1)*(l+2)/2,-l:+l))
M_trans = 0.0d0


! Use hard-coded matrix from table provided in the exercise if l < 3
if (l < 3) then
  call calc_M_trans_hardcoded(l,M_trans)
! calculate transformation matrix if l > 2
else
  call calc_M_trans(l,M_trans) 
end if

! Print calculated transformation matrix in terminal
call get_index_list(l,index_list)
write(*,'(A,I3)') "Transformation matrix for l = ", l
write(*,'(/)')
write(*,'(10X,A,5X)',advance='no') "m :"
do m = -l, l
 write(*,'(I3,9X)', advance='no') m  ! I3 = width 3, 1X = space, no newline
end do
write(*,'(/)') 
do i=1,(l+1)*(l+2)/2
 write(x_exp_char,'(I0)') index_list(i,1)
 write(y_exp_char,'(I0)') index_list(i,2)
 write(z_exp_char,'(I0)') index_list(i,3)
 write(*,'(A)',advance='no') "x^"//trim(x_exp_char)//" y^"//trim(y_exp_char)//" z^"//trim(z_exp_char)//" : "
 do m = -l, l
   write(*,'(F12.6)',advance='no') M_trans(i,m)
 end  do
write(*,*)
end do




!call get_index_list(l,index_list)
!a = 0.5
! Iterate over m = -l,...,+l
!do m = -l,l
!S_self = 0.0d0
!write(*,*) "m", m
  ! Iterate over cartesian functions i and j
!  do i = 1,(l+1)*(l+2)/2
!    do j = 1,(l+1)*(l+2)/2
!      ! The exponents of the cartesians 
!      x_exp = index_list(i,1) + index_list(j,1)
!      y_exp = index_list(i,2) + index_list(j,2)
!      z_exp = index_list(i,3) + index_list(j,3)
!      !write(*,*) "exponents", x_exp, y_exp, z_exp
!      if (mod(x_exp,2)==0.and.mod(y_exp,2)==0.and.mod(z_exp,2)==0) then
!        !write(*,*) "Non-vanishing contribution to integral!"
!        call calc_S_self_cart(x_exp/2,y_exp/2,z_exp/2,a,S_self_cart)
!        S_self = S_self + M_trans(i,m)*M_trans(j,m)*S_self_cart
!        !write(*,*) "S_self", S_self
!      end if
!    end do
!  end do
!write(*,*) "Self overlap for m = ", m, ": ", S_self
!end do

end program main

!Calculation of 1D self-overlap function from Eq. (15)
double precision function S1D_self(i,a) result(S)
  implicit none
  integer, intent(in) :: i
  double precision, intent(in) :: a
  double precision, parameter :: pi = acos(-1.d0)
  double precision :: double_factorial

  ! Eq. (15): <G_i|G_i> = (2i-1)!! / ( (4a)^i ) * sqrt(pi/(2a))
  S = double_factorial(2.d0*real(i,8) - 1.d0) / ((4.d0*a)**i) * sqrt(pi/(2.d0*a))
end function S1D_self


! Calculate 3D Cartesian self-overlap using product of three 1D overlaps
subroutine calc_S_self_cart_ijk(i,j,k,a,S)
  implicit none
  integer, intent(in) :: i,j,k
  double precision, intent(in) :: a
  double precision, intent(out) :: S
  double precision :: S1D_self

  S = S1D_self(i,a) * S1D_self(j,a) * S1D_self(k,a)
end subroutine calc_S_self_cart_ijk


subroutine write_array(arr, m, n)
integer, intent(in) :: m, n
integer, intent(in), dimension(m,n) :: arr
integer :: i
do i=1,m
  write(*,*) arr(i,:)
end do
end subroutine write_array

function factorial(n) result(res)
implicit none
integer, intent(in) :: n
integer :: i
double precision :: res
res = 1.0d0
do i = 1,n
  res = res*real(i,kind=8)
end do
end function factorial

! This function calculates the double factorial of a given double-precision real
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

function delta(a,b) result(res)
implicit none
integer, intent(in) :: a,b
integer :: res
if (a == b) then
res = 1
else
  res = 0
end if
end function delta

function nCk(n,k) result(res)
implicit none
integer,intent(in) :: n,k
double precision :: res,factorial
res=factorial(n)/(factorial(k)*factorial(n-k))
end function nCk

! This subroutine creates a sorted list of the cartesian exponents (a,b,c) in  x^a y^b z^c for a given l = a+b+c.
! The list is ordered in such a way that functions with contributions from one cartesian component come first, contributions from two components come second and contributions from three components come after that.
! E.g., for l = 3, the desired order is x^3, y^3, z^3, x^2y, x^2z, xy^2, y^2z, xz^2, yz^2, xyz and the output is 
subroutine get_index_list(l,index_list_fin)
implicit none
integer, intent(in) :: l
integer, dimension(3 + 3*l + 2*(l**2),3) :: index_list
integer, dimension((l+1)*(l+2)/2,3),intent(out) :: index_list_fin
integer :: i, j, idx, a, b, d
logical :: found
index_list(:,:) = 0
! x,y,z with exponent l
do i = 1,3
  index_list(i,i) = l
end do
idx = 4
do a = l-1,ceiling(real(l)/2),-1
  index_list(idx,:) = [a, l-a, 0]
  index_list(idx+1,:) = [a, 0, l-a] 
  index_list(idx+2,:) = [l-a, a, 0]
  index_list(idx+3,:) = [0, a, l-a]
  index_list(idx+4,:) = [l-a, 0, a]
  index_list(idx+5,:) = [0, l-a, a]
  idx = idx + 6
end do

do a = l-2,ceiling(real(l)/3),-1
  do b = 1,min(l-a-1,a),+1
    d = l-a-b
    index_list(idx,:) = [a,b,d]
    index_list(idx+1,:) = [a,d,b]
    index_list(idx+2,:) = [b,a,d]
    index_list(idx+3,:) = [d,a,b]
    index_list(idx+4,:) = [d,b,a]
    index_list(idx+5,:) = [b,d,a]
    idx = idx + 6
  end do
end do
!write(*,*) "Index list pre"
!call write_array(index_list,3 + 3*l + 2*(l**2),3)
!write(*,*) "Final index list"
index_list_fin = 0
idx = 1
do i = 1,(l+1)*(l+2)
  found = .false.
  do j = 1,(l+1)*(l+2)/2
    if (all(index_list_fin(j,:)==index_list(i,:))) then
      found=.true.
    end if
  end do
  if ( .not. found .and. .not. all(index_list(i,:) == [0,0,0]) ) then
    index_list_fin(idx,:) = index_list(i,:)
    idx = idx + 1
  end if
end do

end subroutine get_index_list

subroutine get_index(l,index_list,exponents,i)
integer,intent(in) :: l
integer,dimension((l+1)*(l+2)/2,3),intent(in) :: index_list
integer,dimension(3),intent(in) :: exponents
integer,intent(out) :: i
do i=1,(l+1)*(l+2)/2
  if (all(index_list(i,:) == exponents)) then
    exit
  end if
end do
end subroutine get_index

subroutine calc_M_trans(l,M_trans)
implicit none
integer, intent(in) :: l
integer, dimension((l+1)*(l+2)/2,3) :: index_list
double precision, dimension((l+1)*(l+2)/2,-l:+l), intent(out) :: M_trans
integer :: delta, m, t, u, v_count, x_exp, y_exp, z_exp, i
double precision :: vm, v, C, Nlm, factorial, nCk
integer, dimension(3) :: exponents
call get_index_list(l,index_list)
do m = -l,l
  Nlm = 1/(2.0d0**abs(m)*factorial(l))*SQRT(real(2*factorial(l+abs(m))*factorial(l-abs(m))/(2**delta(0,m)),kind=8))
  if (m.ge.0) then
    vm = 0.0d0
  else
    vm = 0.5d0
  end if
  do t = 0, floor(real((l-abs(m)),kind=8)/2)
    do u = 0, t
      do v_count = 0, floor(real(abs(m),kind=8)/2-vm)
        v = v_count + vm
        C = (-1)**(t+v-vm)*(0.25d0**t)*nCk(l,t)*nCk(l-t,abs(m)+t)*nCk(t,u)*nCk(abs(m),int(2*v))
        x_exp = 2*t+abs(m)-int(2*(u+v))
        y_exp = int(2*(u+v))
        z_exp = l-2*t-abs(m)
        exponents = (/ x_exp, y_exp, z_exp /)
        call get_index(l,index_list,exponents,i)
        M_trans(i,m) = M_trans(i,m) + Nlm*c
      end do
    end do
  end do
end do
end subroutine calc_M_trans

subroutine calc_M_trans_hardcoded(l,M_trans)
integer, intent(in) :: l
double precision, dimension((l+1)*(l+2)/2,-l:+l), intent(out) :: M_trans
if (l==0) then
  M_trans(1,0) = 1.0d0
else if (l==1) then
  M_trans = 0.0d0
  M_trans(1,1) = 1.0d0
  M_trans(2,-1) = 1.0d0
  M_trans(3,0) = 1.0d0
else if (l==2) then
 M_trans = 0.0d0
 ! m = 2
 M_trans(1,2) = SQRT(3.0d0)/2
 M_trans(2,2) = -SQRT(3.0d0)/2
 ! m = 1
 M_trans(5,1) = SQRT(3.0d0)
 ! m = 0
 M_trans(1,0) = -0.5
 M_trans(2,0) = -0.5
 M_trans(3,0) = 1
 ! m = -1
 M_trans(6,-1) = SQRT(3.0d0)
 ! m = -2
 M_trans(4,-2) = SQRT(3.0d0)
end if
end subroutine calc_M_trans_hardcoded

subroutine calc_S_self_cart(i,j,k,a,S)
implicit none
integer, intent(in) :: i,j,k
double precision, intent(in) :: a
double precision, intent(out) :: S
double precision, parameter :: pi=acos(-1.d0)
double precision :: double_factorial

S = double_factorial(2*real(i,8)-1)*double_factorial(2*real(j,8)-1)*double_factorial(2*real(k,8)-1)/((4*a)**(i+j+k))*(pi/(2*a))**(3/2)

end subroutine calc_S_self_cart






