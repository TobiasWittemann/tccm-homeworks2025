program MD
implicit none
integer :: NAtoms, i, j, nsteps
double precision :: m, dist, epsilon, sigma, epsilon_ext, sigma_ext, lj_potential, V, dt
double precision, allocatable :: coord(:,:), coord_ext(:,:), mass(:)
double precision, allocatable :: distance(:,:), acceleration(:,:), velocity(:,:)
character(len=100) :: file

! Get input filename from command line argument
! call get_command_argument(1, file)
call get_user_input(file, nsteps, dt)

! Read the input file and allocate arrays
call read_NAtoms(file,NAtoms)
allocate(coord(NAtoms,3), coord_ext(NAtoms,3), mass(NAtoms))
allocate(distance(NAtoms,NAtoms),acceleration(NAtoms,3))
allocate(velocity(NAtoms,3))
call read_molecule(file, NAtoms, epsilon, sigma, coord, coord_ext, mass)
 
! Compute distances
call compute_distances(NAtoms, coord, distance)
write(*,*) "distance"
call write_array(distance, NAtoms, NAtoms)

! Compute acceleration
call compute_acc(NAtoms, epsilon, sigma, coord, mass, distance, acceleration)
write(*,*) "Acceleration"
call write_array(acceleration, NAtoms,3)

! Compute Lennard-Jones potential
lj_potential = V(epsilon, sigma, NAtoms, distance)
write(*,*) 'Lennard-Jones Potential (kJ/mol): ', lj_potential

! Run MD simulation

! Input the number of MD steps and the time step in ps
! do
!     write(*,*) 'Enter the number of MD steps ():'
!     read(*,*, iostat=io_status) resolution
!     if (io_status /= 0 .or. resolution <= 0 .or. resolution > max_line)  then
!         write(*,*) 'Invalid input. Resolution must be positive and less than ', max_line
!         cycle
!     else
!         exit
!     end if
! end do


nsteps = 1000 ! number of MD steps
dt = 0.02d0 ! time step in ps
velocity = 0.0d0 ! initialize velocities to zero
call run_md(NAtoms, nsteps, dt, distance, mass, epsilon, sigma, coord, velocity, acceleration)
write(*,*) 'Simulation completed.'
write(*,*) 'Final Coordinates after MD (nm):'
call write_array(coord, NAtoms, 3)
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
read(1,*) NAtoms 
close(1)        ! read number of atoms
end subroutine read_NAtoms

subroutine read_molecule(file, NAtoms, epsilon, sigma, coord, coord_ext, mass)
implicit none
character(len=100), intent(in) :: file
integer, intent(in) :: NAtoms
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
write(*,*) "Number of atoms  : ", NAtoms
write(*,*) "Epsilon / kJ/mol : ", epsilon_ext
write(*,*) "Sigma / A        : ", sigma_ext
! Allocate Coordinate array
do i=1,NAtoms
  write(*,*) "i", i
  read(1,*) coord_ext(i,1), coord_ext(i,2), coord_ext(i,3), mass(i)    ! read coordinates and mass line by line
end do
coord = 0.1*coord_ext ! Convert coordinates from A into nm
close(1)
end subroutine read_molecule

subroutine compute_distances(NAtoms, coord, distance)
integer,intent(in) :: NAtoms
integer :: i,j
double precision :: dx2, dy2, dz2, m, dist
double precision, intent(in) :: coord(NAtoms,3)
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


subroutine compute_acc(NAtoms, epsilon, sigma, coord, mass, distance, acceleration)
implicit none
integer,intent(in) :: NAtoms
integer :: i, j
double precision, intent(in) :: coord(NAtoms,3), distance(NAtoms,NAtoms), mass(NAtoms), epsilon, sigma
double precision, intent(out) :: acceleration(NAtoms,3)
double precision, dimension(3) :: F, r, U_ij
do i=1,NAtoms
  F = 0
  do j=1,NAtoms
    if (i.ne.j) then
      r = distance(i,j)
      U_ij = 24*epsilon/r*((sigma/r)**6 - 2*(sigma/r)**12)
      F = F + U_ij*(coord(i,:)-coord(j,:))/r
    end if
  end do
  acceleration(i,:) = -F/mass(i)
end do
end subroutine compute_acc

double precision function V(epsilon, sigma, NAtoms, distance)
    implicit none
    double precision, intent(in) :: epsilon, sigma
    integer, intent(in) :: NAtoms
    double precision, intent(in) :: distance(NAtoms, NAtoms)

    integer :: i, j
    double precision :: r, sr6, sr12, pair_potential
    V = 0.0d0

    ! sum over i and j>i to avoid double counting
    do j = 2, NAtoms
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

subroutine run_md(NAtoms, nsteps, dt, distance, mass, epsilon, sigma, coord, velocity, acceleration)
    implicit none
    integer, intent(in) :: NAtoms, nsteps
    double precision, intent(in) :: dt, mass(NAtoms), epsilon, sigma
    double precision, intent(inout), dimension(NAtoms,3) :: coord, velocity, acceleration
    double precision, intent(inout) :: distance(NAtoms, NAtoms)
    integer :: step, i, unit_traj
    double precision :: lj_potential, V
    ! double precision :: ke, pe, total_e, initial_e
    ! logical :: first_output = .true. 
    
    unit_traj = 20
    open(unit_traj, file='traj.xyz', status='replace', action='write')
    
    call compute_distances(NAtoms, coord, distance)
    call compute_acc(NAtoms, epsilon, sigma, coord, mass, distance, acceleration)
    do step = 1, nsteps
        coord = coord + velocity*dt + 0.5d0*acceleration*dt**2 ! update position r(n+1)
        call compute_distances(NAtoms, coord, distance) ! update distances r(n+1)
        velocity = velocity + 0.5d0*acceleration*dt ! update half-step velocity v(n+1/2)
        call compute_acc(NAtoms, epsilon, sigma, coord, mass, distance, acceleration) ! update acceleration a(n+1)
        velocity = velocity + 0.5d0*acceleration*dt ! update final step velocity v(n+1)
    
        ! Write trajectory every 10 steps
        if (mod(step, 10) == 0) then
            lj_potential = V(epsilon, sigma, NAtoms, distance)
            ! TODO: compute kinetic energy T, total_energy and check energy conservation
            ! ke = 0.0d0
            ! do i = 1, NAtoms
            !     ! 速度的平方和: v_x^2 + v_y^2 + v_z^2
            !     ke = ke + 0.5d0 * mass(i) * sum(velocity(i,:)**2)
            ! end do

            ! ! 2. 计算势能 (PE)
            ! ! 确保此时的 distance 矩阵是最新的
            ! pe = V(epsilon, sigma, NAtoms, distance)

            ! ! 3. 计算总能量
            ! total_e = ke + pe

            ! ! 4. 检查是否收敛/守恒
            ! if (first_output) then
            !     initial_e = total_e
            !     first_output = .false.
            !     write(*, '(A)') "  Step      Kinetic        Potential          Total        Relative ΔE"
            ! end if

            ! ! 计算相对误差 (Relative Drift)
            ! ! 如果值很小（如 10^-4 以下），说明算法收敛且稳定
            ! write(*, '(I6, 3F15.6, E15.4)') &
            !     step, ke, pe, total_e, abs((total_e - initial_e)/initial_e)

            ! Write coordinates to XYZ trajectory file
            write(unit_traj,*) NAtoms
            ! TODO: write comment line with energy info
            do i = 1, NAtoms
                write(unit_traj,'(A,3F12.6)') 'Ar', coord(i, :) * 10.0d0 ! write coordinates in Angstroms
            end do
        end if
    end do
    close(unit_traj)
end subroutine run_md


subroutine get_user_input(file_name, nsteps, dt)
    implicit none
    character(len=100), intent(out) :: file_name
    integer, intent(out) :: nsteps
    double precision, intent(out) :: dt
    
    character(len=100) :: input_str
    integer :: io_status

    print *, "=========================================="
    print *, "      Molecular Dynamics Configuration      "
    print *, "=========================================="

    ! get input filename
    do
        write(*, '(A)', advance='no') '>> Input filename (default "inp.txt"): '
        read(*, '(A)') input_str
        if (len_trim(input_str) == 0) then
            file_name = 'inp.txt'
        else
            file_name = trim(input_str)
        end if
        
        open(unit=10, file=trim(file_name), status='old', iostat=io_status)
        if (io_status == 0) then
            close(10)    ! check if file exists
            exit
        else
            write(*,*) 'Error: File not found. Please try again.'
        end if
    end do

    ! get number of steps
    do
        write(*, '(A)', advance='no') '>> Number of steps (default 1000): '
        read(*, '(A)') input_str
        if (len_trim(input_str) == 0) then
            nsteps = 1000
            exit
        else
            read(input_str, *, iostat=io_status) nsteps
            if (io_status == 0 .and. nsteps > 0) exit
            write(*,*) 'Invalid input. Enter a positive integer.'
        end if
    end do

    ! get time step
    do
        write(*, '(A)') '>> Time step dt in ps (default 0.02): '
        read(*, '(A)') input_str
        if (len_trim(input_str) == 0) then
            dt = 0.02d0
            exit
        else
            read(input_str, *, iostat=io_status) dt
            if (io_status == 0 .and. dt > 0.0d0) exit
            write(*,*) 'Invalid input. Enter a positive number.'
        end if
    end do
    
    print *, "------------------------------------------"
end subroutine get_user_input