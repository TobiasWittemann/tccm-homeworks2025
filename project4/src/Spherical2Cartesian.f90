program main
implicit none
integer :: l, m, t, u, v_count, factorial, delta,nCk, x_exp, y_exp, z_exp, i, j
integer, dimension(:,:), allocatable :: index_list
double precision :: vm, v, C
double precision :: Nlm
double precision, dimension(:,:), allocatable :: M_trans
logical :: found
integer,dimension(3) :: exponents
character(len=1) :: x_exp_char, y_exp_char, z_exp_char

! Get angular momentum input from user
write(*,*) "Please select angular momentum l"
read(*,*) l


! Calculate transformation matrix
allocate(index_list((l+1)*(l+2)/2,3))
index_list = 0
allocate(M_trans((l+1)*(l+2)/2,-l:+l))
M_trans = 0.0d0
call get_index_list(l,index_list)
call calc_M_trans(l,index_list,M_trans) 


! Print calculated transformation matrix in terminal
write(*,'(A,I3)') "Transformation matrix for l = ", l
write(*,'(/)')
write(*,'(10X,A,5X)',advance='no') "m :"
do m = -l, l
 write(*,'(I3,7X)', advance='no') m  ! I3 = width 3, 1X = space, no newline
end do
write(*,'(/)') 
do i=1,(l+1)*(l+2)/2
 write(x_exp_char,'(I0)') index_list(i,1)
 write(y_exp_char,'(I0)') index_list(i,2)
 write(z_exp_char,'(I0)') index_list(i,3)
 write(*,'(A)',advance='no') "x^"//x_exp_char//" y^"//y_exp_char//" z^"//z_exp_char//" : "
 do m = -l, l
   write(*,'(F10.6)',advance='no') M_trans(i,m)
 end  do
write(*,*)
end do


end program main


subroutine write_array(arr, m, n)
integer, intent(in) :: m, n
double precision, intent(in), dimension(m,n) :: arr
integer :: i
do i=1,m
  write(*,*) arr(i,:)
end do
end subroutine write_array

function factorial(n) result(res)
implicit none
integer, intent(in) :: n
integer :: res, i
res = 1
do i = 1,n
  res = res*i
end do
end function factorial

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
integer :: res,factorial
res=factorial(n)/(factorial(k)*factorial(n-k))
end function nCk

! This subroutine creates a sorted list of the cartesian exponents (a,b,c) in  x^a y^b z^c for a given l = a+b+c.
! The list is ordered in such a way that functions with contributions from one cartesian component come first, contributions from two components come second and contributions from three components come after that.
! E.g., for l = 3, the desired order is x^3, y^3, z^3, x^2y, x^2z, xy^2, y^2z, xz^2, yz^2, xyz and the output is 
subroutine get_index_list(l,index_list_fin)
implicit none
integer, intent(in) :: l
integer, dimension((l+1)*(l+2),3) :: index_list
integer, dimension((l+1)*(l+2)/2,3),intent(out) :: index_list_fin
integer :: i, j, idx, a, b, d
logical :: found
idx = 1
! x,y,z with exponent l
do i = 1,3
  index_list(i,i) = l
  idx = idx + 1
end do
do a = l-1,ceiling(real(l)/2),-1
  do i = 1,3
    if (i==1) then
      index_list(idx,i) = a
      index_list(idx,2) = l-a
      idx = idx + 1
      index_list(idx,i) = a
      index_list(idx,3) = l-a
      idx = idx + 1
    else if (i==2) then
      index_list(idx,i) = a
      index_list(idx,1) = l-a
      idx = idx + 1
      index_list(idx,i) = a
      index_list(idx,3) = l-a
      idx = idx + 1
    else if (i==3) then
      index_list(idx,i) = a
      index_list(idx,1) = l-a
      idx = idx + 1
      index_list(idx,i) = a
      index_list(idx,2) = l-a
      idx = idx + 1
    end if
  end do
end do

do a = l-2,ceiling(real(l-1)/2),-1
  do b = 1,floor(real(l-1)/2)
    d = l-a-b
    index_list(idx,1) = a
    index_list(idx,2) = b
    index_list(idx,3) = d
    idx = idx + 1
    index_list(idx,1) = a
    index_list(idx,2) = d
    index_list(idx,3) = b
    idx = idx + 1
    index_list(idx,1) = b
    index_list(idx,2) = a
    index_list(idx,3) = d
    idx = idx + 1
    index_list(idx,1) = d
    index_list(idx,2) = a
    index_list(idx,3) = b
    idx = idx + 1
    index_list(idx,2) = b
    index_list(idx,1) = d
    index_list(idx,3) = a
    idx = idx + 1
    index_list(idx,2) = d
    index_list(idx,1) = b
    index_list(idx,3) = a
    idx = idx + 1
  end do
end do

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

subroutine calc_M_trans(l,index_list,M_trans)
implicit none
integer, intent(in) :: l
integer, dimension((l+1)*(l+2)/2,3),intent(in) :: index_list
double precision, dimension((l+1)*(l+2)/2,-l:+l), intent(out) :: M_trans
integer :: factorial, nCk, delta, m, t, u, v_count, x_exp, y_exp, z_exp, i
double precision :: vm, v, C, Nlm
integer, dimension(3) :: exponents
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

