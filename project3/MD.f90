program MD
implicit none
integer :: Natoms, i, j
double precision :: dx2, dy2, dz2, m, dist, epsilon, sigma, lj_potential, V
double precision :: m, dist, epsilon, sigma, epsilon_ext, sigma_ext
double precision, allocatable :: coord(:,:), coord_ext(:,:), mass(:)
double precision, allocatable :: distance(:,:), acceleration(:,:)
character(len=100) :: file

! Get input filename from command line argument
call get_command_argument(1, file)

! Read the input file and allocate arrays
call read_NAtoms(file,NAtoms)
allocate(coord(NAtoms,3), coord_ext(NAtoms,3), mass(NAtoms))
allocate(distance(NAtoms,NAtoms),acceleration(NAtoms,3))
call read_molecule(file, NAtoms, epsilon, sigma, coord, coord_ext, mass)
 
! Compute distances
call compute_distances(NAtoms, coord, distance)
write(*,*) "distance"
call write_array(distance, NAtoms, NAtoms)

! Compute acceleration
call compute_acc(Natoms, epsilon, sigma, coord, mass, distance, acceleration)
write(*,*) acceleration

! Compute Lennard-Jones potential
lj_potential = V(epsilon, sigma, Natoms, distance)
write(*,*) 'Lennard-Jones Potential (kJ/mol): ', lj_potential
end program MD


subroutine write_array(arr, m, n)
integer, intent(in) :: m, n
double precision, intent(in), dimension(m,n) :: arr
integer :: i
do i=1,m
  write(*,*) arr(i,:)
end do
end subroutine write_array

subroutine read_NAtoms(file, NAtoms)
implicit none
character(len=100), intent(in) :: file
integer, intent(out) :: NAtoms
open(1, file=file, status="old", action='read')
read(1,*) Natoms 
close(1)        ! read number of atoms
end subroutine read_NAtoms

subroutine read_molecule(file, NAtoms, epsilon, sigma, coord, coord_ext, mass)
implicit none
character(len=100), intent(in) :: file
integer, intent(in) :: Natoms
integer :: i,j, NAtomsDummy
double precision :: m, dist, epsilon_ext, sigma_ext
double precision, intent(out) :: epsilon, sigma 
double precision, intent(out), dimension(NAtoms,3) :: coord, coord_ext
double precision, intent(out), dimension(NAtoms) :: mass
! Open and read input file
open(1, file='inp.txt', status="old", action='read')
read(1,*) NAtomsDummy         ! read number of atoms
read(1,*) epsilon_ext, sigma_ext ! read epsilon (kJ/mol) and sigma (A) from comment line
! Convert epsilon and sigma from externally used units into the units used in the program
epsilon = epsilon_ext ! 1 kJ/mol = 1 gnm^2/ps^2
sigma = 0.1*sigma_ext ! 1 A = 0.1 nm
! Writing to the output file
write(*,*) "Number of atoms  : ", Natoms
write(*,*) "Epsilon / kJ/mol : ", epsilon_ext
write(*,*) "Sigma / A        : ", sigma_ext
! Allocate Coordinate array
do i=1,NAtoms
  write(*,*) "i", i
  read(1,*) coord_ext(i,1), coord_ext(i,2), coord_ext(i,3), mass(i)    ! read coordinates and mass line by line
end do
coord = 0.1*coord_ext ! Convert coordinates from A into pm
close(1)
end subroutine read_molecule

subroutine compute_distances(NAtoms, coord, distance)
integer,intent(in) :: NAtoms
integer :: i,j
double precision :: dx2, dy2, dz2, m, dist
double precision, intent(in) :: coord(NAtoms,NAtoms)
double precision, intent(out) :: distance(NAtoms,NAtoms)
do j=1,NAtoms
  do i=1,j
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