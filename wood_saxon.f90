program wood_saxon
    use potentials
    use ieee_arithmetic
    use, intrinsic :: iso_fortran_env
    implicit none
    ! Ground state
    ! n = 0, l = 0, j = 1/2 # 2 nucleons
    ! n = 1, l = 1, j = 3/2 # 4 nucleons
    ! n = 1, l = 1, j = 1/2 # 2 nucleons
    ! Schrodinger Eq.
    ! -hbar^2/2m \psi" + V_WS/r + hbar^2/2m * l(l+1)/r^2 + V_SO + V_Coul
    ! B.C. u(0), u(R_max) = 0.0
    real(dp), parameter :: h = 1.e-2_dp ! fm, mesh
    real(dp), parameter :: R_max = 10.0_dp ! Maximum
    INTEGER :: i, ii, iiii,k, num, &
    l, N, Z, isospin, info, lda, lwork, unit_number
    REAL(dp), allocatable :: psi(:), r_array(:)
    REAL(dp) :: j, vx
    REAL(dp) :: dummy(1, 1)
    ! print *, hbar ** 2 / (2 * effective_mass)
    num = int(R_max/h) + 1


    allocate(psi(num), r_array(num))
    print *, e2

    do i = 2, num, 1
        r_array(i) = DBLE(i-1) * h
    enddo
    r_array(1) = 1.e-8
    ! psi(1) = 0.0_dp
    ! psi(2) = 1.e-10_dp

    N = 8
    Z = 8
    l = 2
    j = 1.0_dp/2.0_dp
    isospin = 1

    unit_number = 20
    open(unit_number, file="potential.txt",&
     status="unknown", action="write")

    ! Write header to the file
    write(unit_number, '(A)') "r [fm]    Potential [MeV]"

    ! Write data into two columns
    do i = 2, num-1
        vx = V_SO(r_array(i), l, j, N, Z, isospin) + &
            V_centrifugal(r_array(i), l) + V_WS(r_array(i), N, Z, isospin) + &
            V_C(r_array(i) , N, Z, isospin)
        write(unit_number, '(F10.5, F15.12, F10.5)') vx, r_array(i)
    end do

    ! Close the file
    close(unit_number)
end program

