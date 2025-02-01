program Shrodinger
    implicit none
    integer(8), parameter :: dp = selected_real_kind(15)
! Infinite potential well
! Solving the Shrodinger equation using Numerov method
! -hbar/2m * \psi" = E\psi
    real(dp), parameter :: hbar = 6.58211899 * 1.e-22 ! MeV * s
    real(dp), parameter :: c = 299792458 * 1.e15 ! fm/s
    real(dp), parameter :: m_n = 938.272013 / c ** 2 ! MeV/c^2
    real(dp), parameter :: alpha = 20.75_dp ! hbar ** 2 / 2m
    real(dp), parameter :: a = 6.0_dp ! fm
    real(dp), parameter :: h = 1.e-3_dp ! fm
    real(dp), parameter :: pi = 3.141592 ! 

! Numerov method
! (1 + h**2/12 * v(x+h)) * f(x+h) = 2 * (1-5*h**2/12 * v(x))*f(x)
! - (1+h**2/12 * v(x-h)) * f(x-h)
    integer(8) :: n, i, ii, Node, Node_count, unit_number
    real(dp), allocatable :: psi(:)
    real(dp) :: v, E_trial, E_up, E_down, a1, a2, a3, SUM, rand_num
    n = int(a/h) + 1
    allocate(psi(n))
    psi(1) = 0.0_dp
    psi(2) = 1.e-10_dp
    E_down = -3000.0
    E_up = 0.0
    E_trial = (E_up + E_down) / 2
    v = E_trial/alpha
    Node = 1
    
    print *, hbar ** 2 / (2 * m_n)
    do ii = 1, 1000
        E_trial = (E_up + E_down) / 2.
        v = E_trial * (2 * m_n) / hbar ** 2
        a1 = 1. + h**2 / 12. * v
        a2 = 2. * (1. - 5. * h**2 / 12. * v)
        a3 = 1. + h**2 / 12. * v
        do i = 2, n-1
            
            psi(i + 1) = (a2 * psi(i) - a3 * psi(i-1)) / a1
        enddo


!--------Node Count---------
        Node_count = 0

        do i = 1, n-1
        
            if (psi(i) * psi(i+1) < 0.0) then
                Node_count = Node_count + 1
            endif
        
        enddo
        
        print *, "Total Nodes for this iteration: ", Node_count
!------Adjust E_trial--------
        if (abs(E_up - E_down) < 1.e-10) then
            exit
        elseif (Node_count < Node) then
            E_down = E_trial
        
        elseif (Node_count >= Node) then
            E_up = E_trial
        
        endif

    
    
    print *, ' '
    print *, "Trial: ", ii
    print *, "Total Nodes: ", Node_count
    ! print *, "PSI(N)", psi(n)
    print *, ABS(E_up - E_down), "Up Down Difference"
    print *, "E_trial = ", E_trial
    print *, "E = ", (Node * pi * hbar / a) ** 2 / (2 * m_n)
    
    
    enddo


!----Renormalization of Psi----
    SUM = 0.0_dp
    do iii = 1, n-1
        SUM = h / 2. * (psi(iii+1) ** 2 + psi(iii) ** 2)
    enddo

    ! psi = psi / sqrt(SUM)

    ! SUM = 0.0_dp
    ! do iii = 1, n-1
    !     SUM = h / 2. * (psi(iii+1) ** 2 + psi(iii) ** 2)
    ! enddo
    ! print *, SUM


unit_number = 20
    open(unit_number, file="output1.txt", status="unknown", action="write")

    ! Write header to the file
    write(unit_number, '(A)') "x    y"

    ! Write data into two columns
    do i = 1, n
        write(unit_number, '(F10.5, F15.12)') h * (i - 1), psi(i)
    end do

    ! Close the file
    close(unit_number)

end program Shrodinger