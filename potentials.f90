module Potentials
    implicit none
    integer(8), parameter :: dp = selected_real_kind(15)
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
    ! print *, V_centrifugal(1._dp, 1)



contains

real(dp) function V_WS(r, N, Z, isospin) ! , effective_mass)
    implicit none
    REAL(dp), INTENT(in) :: r
    INTEGER, INTENT(in) :: N, Z, isospin
    REAL(dp) :: R_q, t, V0
    R_q = r0 * DBLE(N+Z) ** (1.0_dp/3.0_dp)
    t = (r - R_q) / a
    if (isospin == 1) then
        V0 = (51 + 33 * (N-Z)/(N+Z))
    elseif (isospin == -1) then
        V0 = (51 - 33 * (N-Z)/(N+Z))
    endif

    V_WS = -2 * effective_mass / hbar ** 2 * (-1) * &
            V0 / (1 + exp(t))
end function V_WS


real(dp) function V_centrifugal(r, l) ! , effective_mass)
    REAL(dp), INTENT(in) :: r ! , effective_mass
    INTEGER, INTENT(in) :: l

    V_centrifugal = hbar ** 2 * (DBLE(l) * DBLE(l+1)) / r**2.0_dp 
    ! hbar ** 2 / (2 * effective_mass) *

end function V_centrifugal

real(dp) function V_SO(r, l, j, N, Z, isospin)! effective_mass, N, Z)
    REAL(dp), INTENT(in) :: r, j
    INTEGER, INTENT(in) :: l, N, Z, isospin
    REAL(dp) :: t, R_q, V0
    R_q = r0 * DBLE(N+Z) ** (1.0_dp/3.0_dp)
    t = (r - R_q)/a
    if (isospin == 1) then
        V0 = 0.44 * (51 + 33 * (N-Z)/(N+Z))
    elseif (isospin == -1) then
        V0 = 0.44 * (51 - 33 * (N-Z)/(N+Z))
    endif
    
    V_SO = -2 * effective_mass * r0 ** 2 * &
    (-1) * exp(t) / a / (1 + exp(t)) ** 2 / r * &
    (j * (j+1) - l * (l+1) - 3.0_dp/4.0_dp) * V0 / 2.0_dp ! Constant

end FUNCTION V_SO

REAL(dp) FUNCTION V_C(r , N, Z, isospin)
    REAL(dp), INTENT(in) :: r
    INTEGER, INTENT(in) :: N, Z, isospin 
    REAL(dp) :: R_p
    R_p = r0 * DBLE(N+Z) ** (1.0_dp/3.0_dp)

    if (isospin == -1) then
        V_C = 0.0_dp
    
    elseif (isospin == 1) then
        if (r <= r_p) then
            V_C = Z * e2 * (3 - (r/R_p) ** 2.0_dp) / (2.0_dp * R_p)&
            * (-2 * effective_mass / hbar ** 2.0_dp)
        
        elseif (r > R_p) then
            V_C = Z * e2 / r * (-2 * effective_mass / hbar ** 2.0_dp)
        endif
    else
        print*, "Invalid Isosplin!"
    
    endif
    end FUNCTION V_C



end module