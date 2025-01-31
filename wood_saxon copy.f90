program wood_saxon
    use potentials
    use ieee_arithmetic
    implicit none
    ! Ground state
    ! n = 0, l = 0, j = 1/2 # 2 nucleons
    ! n = 1, l = 1, j = 3/2 # 4 nucleons
    ! n = 1, l = 1, j = 1/2 # 2 nucleons
    ! Schrodinger Eq.
    ! -hbar^2/2m \psi" + V_WS/r + hbar^2/2m * l(l+1)/r^2 + V_SO + V_Coul
    ! B.C. u(0), u(R_max) = 0.0
    real(dp), parameter :: h = 1.e-5_dp ! fm, mesh
    real(dp), parameter :: R_max = 7.0_dp ! Maximum
    INTEGER :: i, ii, k, num, l, N, Z, isospin
    REAL(dp), allocatable :: psi(:), r_array(:), A(:, :)
    REAL(dp) :: j, vx, vx_back, vx_next, step, &
     E_up, E_down, E_trial, a1, a2, a3
    
    ! print *, hbar ** 2 / (2 * effective_mass)
    num = int(R_max/h) + 1
    allocate(psi(num), r_array(num), A(num,num))
    do i = 1, num, 1
        do k = 1, num, 1
            A(i, j) = 0.0_dp
        enddo
    enddo


    do i = 1, num, 1
        r_array(i) = DBLE(i-1) * h
    enddo
    r_array(1) = 1.e-20
    psi(1) = 0.0_dp
    psi(2) = 1.e-10_dp
    E_up = 0.0_dp
    E_down = -1000.0_dp
    N = 8
    Z = 8
    l = 0
    j = 1/2
    isospin = 1
    ! print *, r_array(1)
    ! print *, r_array(num-2)
    do ii = 1, 100
      E_trial = (E_up + E_down)/2.0_dp  
        do i = 2, num-1, 1
            do k = 1, 3, 1
                vx = V_SO(r_array(i), l, j, N, Z) + &
                V_centrifugal(r_array(i), l) + V_WS(r_array(i), N, Z) + &
                V_C(r_array(i) , N, Z, isospin) ! + &
                ! 2 * effective_mass / hbar ** 2 * E_trial
                ! print *, vx
                vx_next = V_SO(r_array(i+1), l, j, N, Z) + &
                V_centrifugal(r_array(i+1), l) + &
                V_WS(r_array(i+1), N, Z) + V_C(r_array(i+1) , N, Z, isospin) ! + &
                ! 2 * effective_mass / hbar ** 2 * E_trial
                
                vx_back = V_SO(r_array(i-1), l, j, N, Z) + &
                V_centrifugal(r_array(i-1), l) + &
                V_WS(r_array(i-1), N, Z) + V_C(r_array(i-1) , N, Z, isospin) ! + &
                ! 2 * effective_mass / hbar ** 2 * E_trial
                

                a1 = 1 + h**2 / 12 * vx_next
                a2 = 2 * (1 - 5/12 * vx * h**2)
                a3 = 1 + h**2/12 * vx_back
            
            enddo
            
            ! psi(i+1) = (psi(i) * a2 - psi(i-1) * a3) / a1
            if (ieee_is_nan(vx)) then
                print *, "vx is NaN"
                exit
            elseif (ieee_is_nan(vx_back)) then
                print *, "vx_back is NaN"
                exit
            elseif (ieee_is_nan(vx_next)) then
                print *, "vx_next is NaN"
                exit
            endif
        enddo
        print *, "Eigenvalue (Energy): ", E_trial
        print *, "Boundary Value: ", psi(n)
        if (abs(E_up - E_down) < 1.e-5) then
            exit
        elseif (psi(n) < 0.0_dp) then
            E_up = E_trial
        elseif (psi(n) > 0.0_dp) then
            E_down = E_trial
        endif
    
    
    enddo



end program

