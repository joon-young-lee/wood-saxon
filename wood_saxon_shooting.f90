program wood_saxon
    use potentials
    use ieee_arithmetic
    use, intrinsic :: iso_fortran_env
    implicit none

    real(dp), parameter :: h = 1.e-3_dp ! fm, mesh
    real(dp), parameter :: R_max = 10.0_dp ! Maximum
    INTEGER :: i, ii, num, &
    l, N, Z, isospin, unit_number, &
    radial_quntum_number, node_count
    REAL(dp), allocatable :: psi(:), r_array(:)
    REAL(dp) :: j, vx, vx_next, vx_back,&
    E_up, E_down, E_trial,&
    a1, a2, a3
    
    num = int(R_max/h) + 1

    allocate(psi(num), r_array(num))

    do i = 2, num, 1
        r_array(i) = DBLE(i-1) * h
    enddo
    r_array(1) = 1.e-7
    psi(1) = 0.0_dp
    psi(2) = 1.e-10_dp

    radial_quntum_number = 1 ! n = k + 1 (number of nodes)
    N = 8
    Z = 8
    l = 2
    j = 5.0_dp/2.0_dp
    isospin = 1

    E_up = 100.0_dp
    E_down = -1.e6_dp

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
            psi(i + 1) = (a2 * psi(i) - a3 * psi(i-1)) / a1
        
        enddo


!--------Node Count---------
        Node_count = 0

        do i = 1, num-1
        
            if (psi(i) * psi(i+1) < 0.0) then
                Node_count = Node_count + 1
            endif
        
        enddo
        
        print *, "Total Nodes for this iteration: ", Node_count
!------Adjust E_trial--------
        if (abs(E_up - E_down) < 1.e-8) then
            exit
        elseif (Node_count < radial_quntum_number) then
            E_down = E_trial
        
        elseif (Node_count >= radial_quntum_number) then
            E_up = E_trial
        
        endif

    
    
    print *, ' '
    print *, "Trial: ", ii
    print *, "Total Nodes: ", Node_count
    ! print *, "PSI(N)", psi(n)
    print *, ABS(E_up - E_down), "Up Down Difference"
    print *, "E_trial = ", E_trial
    ! print '(F10.2)',"| n = ", radial_quntum_number-1," l = ", l, "j = ", j,">"  
    
    enddo


end program