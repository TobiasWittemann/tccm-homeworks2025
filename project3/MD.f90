program MD
implicit none
integer :: Natoms, i, j
double precision :: dx2, dy2, dz2, m, dist, epsilon, sigma, lj_potential, V
double precision, allocatable :: coord(:,:), mass(:)
double precision, allocatable :: distance(:,:), acceleration(:,:)

! Open and read input file
open(1, file='inp.txt', status="old", action='read')
read(1,*) Natoms    ! read number of atoms
read(1,*)           ! read comment line
write(*,*) "Number of atoms", Natoms
allocate(coord(Natoms, NAtoms), mass(NAtoms))
do i=1,Natoms
  read(1,*) coord(i,1), coord(i,2), coord(i,3), mass(i)    ! read coordinates and mass line by line
end do
close(1)

allocate(distance(NAtoms,NAtoms))

! Compute distances
call compute_distances(NAtoms, coord, distance)
do i=1,NAtoms
  write(*,*) distance(i,:)
end do

! Compute acceleration
epsilon = 0.997 ! kJ/mol
sigma = 3.405   ! Angstrom
allocate(acceleration(NAtoms,NAtoms))
call compute_acc(Natoms, epsilon, sigma, coord, mass, distance, acceleration)
write(*,*) acceleration

! Compute Lennard-Jones potential
lj_potential = V(epsilon, sigma, Natoms, distance)
write(*,*) 'Lennard-Jones Potential (kJ/mol): ', lj_potential
end program MD


subroutine compute_distances(NAtoms, coord, distance)
integer,intent(in) :: NAtoms
integer :: i,j
double precision :: dx2, dy2, dz2, m, dist
double precision, intent(in) :: coord(NAtoms,NAtoms)
double precision, intent(out) :: distance(NAtoms,NAtoms)
do j=1,NAtoms
  do i=1,j
    write(*,*) 'i, j', i, j
    dx2 = (coord(i,1)-coord(j,1))**2
    dy2 = (coord(i,2)-coord(j,2))**2
    dz2 = (coord(i,3)-coord(j,3))**2
    dist = SQRT(dx2 + dy2 + dz2)
    distance(i,j) = dist
    distance(j,i) = dist
  end do
end do
end subroutine compute_distances


subroutine compute_acc(Natoms, epsilon, sigma, coord, mass, distance, acceleration)
implicit none
integer,intent(in) :: NAtoms
integer :: i, j
double precision, intent(in) :: coord(NAtoms,NAtoms), distance(NAtoms,NAtoms), mass(NAtoms), epsilon, sigma
double precision, intent(out) :: acceleration(NAtoms,3)
double precision, dimension(3) :: F, r, U_ij
do i=1,NAtoms
  write(*,*) 'i', i
  F = 0
  do j=1,NAtoms
    if (i.ne.j) then
      r = distance(i,j)
      U_ij = 4*epsilon*((sigma/r)**12 - (sigma/r)**6)
      F = F + U_ij*(coord(i,:)-coord(j,:))/r
    end if
  end do
  acceleration(i,:) = -F/mass(i)
end do
end subroutine compute_acc

double precision function V(epsilon, sigma, Natoms, distance)
    implicit none
    double precision, intent(in) :: epsilon, sigma
    integer, intent(in) :: Natoms
    double precision, intent(in) :: distance(Natoms, Natoms)

    integer :: i, j
    double precision :: r, sr6, sr12, pair_potential
    V = 0.0d0

    ! sum over i and j>i to avoid double counting
    do j = 2, Natoms
        do i = 1, j-1 ! Fortran is column major, i loop inside
            r = distance(i, j)

            ! V_ij = 4*epsilon*((sigma/r)^12 - (sigma/r)^6)
            if (r > 0.0d0) then
                sr6 = (sigma / r)**6
                sr12 = sr6 * sr6
                pair_potential = 4.0d0 * epsilon * (sr12 - sr6)
                V = V + pair_potential
            end if
        end do
    end do
end function V