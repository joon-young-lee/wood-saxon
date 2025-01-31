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
    real(dp), parameter :: R_max = 5.0_dp ! Maximum
    INTEGER :: i, ii, k, num, l, N, Z, isospin, info, lda, lwork, unit_number
    REAL(dp), allocatable :: psi(:), r_array(:),&
     Mat(:, :), WR(:), WI(:), WORK(:)
    REAL(dp) :: j, vx, vx_back, vx_next, step, &
      V1, VN
    REAL(dp) :: dummy(1, 1)
    ! print *, hbar ** 2 / (2 * effective_mass)
    num = int(R_max/h) + 1
    lda = num
    lwork = 4 * num
    allocate(psi(num), r_array(num), &
    Mat(num,num), WR(num), WI(num), WORK(lwork))
    print *, sqrt(e2)
    do i = 1, num, 1
        do k = 1, num, 1
            Mat(i, k) = 0.0_dp
        enddo
    enddo


    do i = 1, num, 1
        r_array(i) = DBLE(i-1) * h
    enddo
    r_array(1) = 1.e-7
    ! psi(1) = 0.0_dp
    ! psi(2) = 1.e-10_dp

    N = 8
    Z = 8
    l = 1
    j = 1.0_dp/2.0_dp
    isospin = -1
    print *, (N+Z) ** (1/3)
    ! print *, r_array(1)
    ! print *, r_array(num-2)
        
    V1 = V_SO(r_array(1), l, j, N, Z, isospin) + &
        V_centrifugal(r_array(1), l) + V_WS(r_array(1), N, Z, isospin) + &
        V_C(r_array(1) , N, Z, isospin)
    Mat(1, 1) = 1/h**2 ! -2.0_dp / h ** 2 + V1
    ! Mat(1, 2) = 0.0_dp ! 1.0_dp / h**2
    ! Mat(num, num-1) = 0.0_dp ! 1.0_dp / h**2
    ! VN = V_SO(r_array(num), l, j, N, Z, isospin) + &
    !     V_centrifugal(r_array(num), l) + V_WS(r_array(num), N, Z, isospin) + &
    !     V_C(r_array(num) , N, Z, isospin)
    Mat(num, num) = 1/h**2
    
    do i = 2, num-1, 1
            vx = V_SO(r_array(i), l, j, N, Z, isospin) + &
            V_centrifugal(r_array(i), l) + V_WS(r_array(i), N, Z, isospin) + &
            V_C(r_array(i), N, Z, isospin) 
            ! print *, vx
            Mat(i, i - 1) = 1.0_dp / h**2
            Mat(i, i) = -2.0_dp / h**2 + vx
            Mat(i, i + 1) = 1.0_dp / h**2
        ! enddo
        ! psi(i+1) = (psi(i) * a2 - psi(i-1) * a3) / a1
        if (ieee_is_nan(vx)) then
            print *, "vx is NaN"
            exit
        endif
    enddo


    CALL DGEEV('N', 'N', lda, Mat, lda, WR, WI, &
            dummy, 1, dummy, 1, WORK, lwork, info)
    ! Check if computation was successful
    if (info /= 0) then
        print *, "Error: DGEEV failed with INFO =", info
        stop
    end if

    ! Print eigenvalues
    print *, "Eigenvalues (Real, Imaginary):"
    do i = 1, lda
        print *, WR(i) * hbar ** 2 / (2*m_n), WI(i)

        ! Free workspace
    ! deallocate(WORK)
        
        ! print *, "Eigenvalue (Energy): ", eigenEnergy
        ! print *, "Boundary Value: ", psi(n)
    
    
    enddo

    unit_number = 20
    open(unit_number, file="potential.txt",&
     status="unknown", action="write")

    ! Write header to the file
    write(unit_number, '(A)') "x    y"

    ! Write data into two columns
    do i = 2, num-1
        vx = V_SO(r_array(i), l, j, N, Z, isospin) + &
            V_centrifugal(r_array(i), l) + V_WS(r_array(i), N, Z, isospin) + &
            V_C(r_array(i) , N, Z, isospin)
        write(unit_number, '(F10.5, F15.12)') vx, r_array(i)
    end do

    ! Close the file
    close(unit_number)
end program

