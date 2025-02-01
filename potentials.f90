module Potentials
    implicit none
    integer(16), parameter :: dp = selected_real_kind(31)
    real(dp), parameter :: hbar = 6.58211899_dp * 1.e-22_dp ! MeV * s
    real(dp), parameter :: c = 299792458 * 1.e15 ! fm/s
    real(dp), parameter :: m_n = 938.272013 / c ** 2 ! MeV/c^2
    real(dp), parameter :: alpha = 1.0_dp/137.036_dp ! hbar ** 2 / 2m
    real(dp), parameter :: a = 0.67_dp ! fm
    real(dp), parameter :: pi = 3.141592_dp
    real(dp), parameter :: r0 = 1.27_dp ! fm
    real(dp), parameter :: e2 = hbar * c * alpha
    real(dp), parameter :: m_star = 1.0_dp * m_n
    real(dp), parameter :: effective_mass = m_n
    real(dp), parameter :: pp = hbar ** 2 / (2 * effective_mass)
     
    ! print *, V_centrifugal(1._dp, 1)



contains

real(dp) function V_WS(r, N, Z, isospin) ! , effective_mass)
    implicit none
    REAL(dp), INTENT(in) :: r
    INTEGER(8), INTENT(in) :: N, Z, isospin
    REAL(dp) :: R_q, t, V0
    R_q = r0 * DBLE(N+Z) ** (1.0_dp/3.0_dp)
    t = (r - R_q) / a
    if (isospin == 1) then
        V0 = (51 + 33 * (N-Z)/(N+Z))
    elseif (isospin == -1) then
        V0 = (51 - 33 * (N-Z)/(N+Z))
    endif

    V_WS = -V0 / (1 + exp(t)) ! * 0.0_dp
end function V_WS


real(dp) function V_centrifugal(r, l) ! , effective_mass)
    REAL(dp), INTENT(in) :: r ! , effective_mass
    INTEGER(8), INTENT(in) :: l

    V_centrifugal = pp * (DBLE(l) * DBLE(l+1)) / (r**2.0_dp)
    ! hbar ** 2 / (2 * effective_mass) *

end function V_centrifugal

real(dp) function V_SO(r, l, j, N, Z, isospin)! effective_mass, N, Z)
    REAL(dp), INTENT(in) :: r, j
    INTEGER(8), INTENT(in) :: l, N, Z, isospin
    REAL(dp) :: t, R_q, V0
    R_q = r0 * DBLE(N+Z) ** (1.0_dp/3.0_dp)
    t = (r - R_q)/a
    if (isospin == 1) then
        V0 = 0.44 * (51 + 33 * (N-Z)/(N+Z))
    elseif (isospin == -1) then
        V0 = 0.44 * (51 - 33 * (N-Z)/(N+Z))
    endif
    
    V_SO = -(r0 ** 2.0_dp) * exp(t) / a /&
     ((1 + exp(t))) ** 2 / r * &
    (j * (j+1) - DBLE(l) * DBLE(l+1) - 3.0_dp/4.0_dp) * V0 / 2.0_dp

end FUNCTION V_SO

REAL(dp) FUNCTION V_C(r , N, Z, isospin)
    REAL(dp), INTENT(in) :: r
    INTEGER(8), INTENT(in) :: N, Z, isospin 
    REAL(dp) :: R_p
    R_p = r0 * DBLE(N+Z) ** (1.0_dp/3.0_dp)

    if (isospin == -1) then
        V_C = 0.0_dp
    
    elseif (isospin == 1) then
        if (r <= r_p) then
            V_C = Z * e2 * (3 - (r/R_p) ** 2.0_dp) / (2.0_dp * R_p)
        
        elseif (r > R_p) then
            V_C = Z * e2 / r
        endif
    else
        print*, "Invalid Isosplin!"
    
    endif
    end FUNCTION V_C


! length in fm, energy in MeV
subroutine WS_shooting(h, R_max, radial_quantum_number, l, j, isospin, &
                        N, Z, E_up, E_down, E)
    INTEGER(8), INTENT(in) :: l, radial_quantum_number, isospin, N, Z
    REAL(dp), INTENT(in) :: h, R_max, j
    REAL(DP), INTENT(inout) :: E_up, E_down
    INTEGER(8) :: num, i, ii, Node_count
    REAL(dp) :: vx, vx_next, vx_back, E_trial,  a1, a2, a3
    REAL(dp), ALLOCATABLE :: r_array(:), u(:)
    REAL(dp), INTENT(out) :: E

    num = int(R_max/h) + 1

    ALLOCATE(r_array(num), u(num))

    do i = 2, num, 1
        r_array(i) = DBLE(i-1) * h
    enddo
    r_array(1) = 1.e-7_dp
    u(1) = 0.0_dp
    u(2) = 1.e-10_dp

    do ii = 1, 1000
        
        E_trial = (E_up + E_down) / 2.0_dp
        
        do i = 2, num-1
            
            ! v(x)
            vx = E_trial - (V_SO(r_array(i), l, j, N, Z, isospin) + &
            V_centrifugal(r_array(i), l) + &
            V_WS(r_array(i), N, Z, isospin) + &
            V_C(r_array(i) , N, Z, isospin))
            vx = vx / pp
            
            ! v(x+h)
            vx_next = E_trial - (V_SO(r_array(i+1), l, j, N, Z, isospin) + &
            V_centrifugal(r_array(i+1), l) + &
            V_WS(r_array(i+1), N, Z, isospin) + &
            V_C(r_array(i+1) , N, Z, isospin))
            vx_next = vx_next / pp
            
            ! v(x-h)
            vx_back = E_trial - (V_SO(r_array(i-1), l, j, N, Z, isospin) + &
            V_centrifugal(r_array(i-1), l) + &
            V_WS(r_array(i-1), N, Z, isospin) + &
            V_C(r_array(i-1) , N, Z, isospin))
            vx_back = vx_back / pp
            
            a1 = 1.0_dp + h**2 / 12.0_dp * vx_next
            a2 = 2.0_dp * (1.0_dp - 5.0_dp * h**2 / 12.0_dp * vx)
            a3 = 1.0_dp + h**2 / 12.0_dp * vx_back
            u(i + 1) = (a2 * u(i) - a3 * u(i-1)) / a1
        
        enddo


!--------Node Count---------
        Node_count = 0

        do i = 1, num-1
        
            if (u(i) * u(i+1) < 0.0) then
                Node_count = Node_count + 1
            endif
        
        enddo
        
        write(*, "(A, i2)") "Total Nodes for this iteration: ", Node_count
!------Adjust E_trial--------
        if (abs(E_up - E_down) < 1.e-6) then
            exit
        elseif (Node_count < radial_quantum_number) then
            E_down = E_trial
        
        elseif (Node_count >= radial_quantum_number) then
            E_up = E_trial
        
        endif

        print *, ' '
        write (*, "(A, i4)") "Trial: ", ii
        write (*, "(A, i2)") "Total Nodes: ", Node_count
        ! print *, "PSI(N)", psi(n)
        write (*, "(f25.10, A)") ABS(E_up - E_down), " Up Down Difference"
        write (*, "(A, f15.10)") "E_trial = ", E_trial
    
    enddo
    
    write (*, "(A)") "----------------------------"
    write (*, "(A)") "State"
    write (*, "(A, i3)") "A: ", N+Z
    write (*, "(A, i3)") "N: ", N
    write (*, "(A, i3)") "Z: ", Z
    write (*, "(A, i1)") "Radial Quantum Number: ", radial_quantum_number-1
    write (*, "(A, i1)") "L: ", l
    write (*, "(A, f3.1)") "J: ", j
    E = E_trial
    write (*, "(A, f10.5, A)") "State Energy: ", E, "MeV"
    write (*, "(A)") "----------------------------"
    DEALLOCATE(r_array, u)
    


end subroutine WS_shooting
end module

