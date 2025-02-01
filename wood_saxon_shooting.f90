program wood_saxon
    use potentials
    use ieee_arithmetic
    use, intrinsic :: iso_fortran_env
    implicit none
    character(len=100) :: line
    character(len=50) :: filename
    real(dp), parameter :: h = 1.e-3_dp ! fm, mesh
    real(dp), parameter :: R_max = 10.0_dp ! Maximum
    INTEGER(8) :: i, ii, num, &
    l, N, Z, isospin, unit_number, &
    radial_quantum_number, node_count, &
    num_entries, io_status
    
    REAL(dp), allocatable :: E_array(:)
    INTEGER(8), ALLOCATABLE :: A_array(:), N_array(:), Z_array(:)
    REAL(dp) :: j, vx, vx_next, vx_back, E_up, E_down, E

    radial_quantum_number = 1 ! n = k + 1 (number of nodes)
    N = 8
    Z = 8
    l = 2
    j = 1.0_dp/2.0_dp
    isospin = 1

    
    filename = "stable_nuclei.csv"  ! Change if needed

    ! First, count the number of rows
    num_entries = 0
    open(unit=10, file=filename, status="old", action="read", iostat=io_status)
    if (io_status /= 0) then
        print *, "Error opening file!"
        stop
    end if
    ! Skip first two lines
    do i = 1, 12
        read(10, "(A)", iostat=io_status) line
    enddo
    
    do
        read(10, "(A)", iostat=io_status) line
        if (io_status /= 0) exit
        num_entries = num_entries + 1
    end do
    close(10)

    ! Print number of entries found
    print *, "Number of stable nuclei found: ", num_entries

    ! Allocate arrays dynamically based on count
    allocate(Z_array(num_entries), N_array(num_entries),&
     A_array(num_entries), E_array(num_entries))
    
    ! Read data into arrays
    
    open(unit=10, file=filename, status="old", action="read", iostat=io_status)
    ! Skip first two lines
    do i = 1, 12
        read(10, "(A)", iostat=io_status) line
    enddo
    if (io_status /= 0) then
        print *, "Error reopening file!"
        stop
    end if

    i = 1
    do
        read(10, "(A)", iostat=io_status) line
        if (io_status /= 0) exit
        read(line, *) Z_array(i), N_array(i), A_array(i)  ! Read and parse CSV
        i = i + 1
    end do
    close(10)

    ! Print some values to verify
    print *, "First 5 entries (Z, N):"
    do i = 1, min(5, num_entries)
        print *, "Z = ", Z_array(i), ", N =", N_array(i)
    end do
    ! N = 0, 1, 2
    ! l = N, N-2,..., 1 or 0
    ! k = (N-l)/2
    radial_quantum_number = 1 ! n = k + 1 (number of nodes)
    N = 8
    Z = 8
    l = 0
    j = 1.0_dp/2.0_dp
    isospin = 1
    E_up = 100.0_dp
    E_down = -1.e6_dp
    ! CALL WS_shooting(h, R_max, &
    !     radial_quantum_number, l, j, isospin, &
    !     N, Z, E_up, E_down, E)
    do i = 1, num_entries, 10
        E_up = 100.0_dp
        E_down = -1.e6_dp
        CALL WS_shooting(h, R_max, &
        radial_quantum_number, l, j, isospin, &
        N_array(i), Z_array(i), E_up, E_down, E)
        E_array(i) = E
    
    enddo
    unit_number = 20
    open(unit_number, file="k=0,l=0,j=0.5.txt",&
     status="unknown", action="write")

    ! Write header to the file
    write(unit_number, '(A)') "r [fm]    Potential [MeV]"

    ! Write data into two columns
    do i = 1, num_entries
        write(unit_number, '(I5, F15.7)') &
        A_array(i), E_array(i)
    end do

    ! Close the file
    close(unit_number)
    DEALLOCATE(Z_array, N_array, A_array, E_array)
end program

