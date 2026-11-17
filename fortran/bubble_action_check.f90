! Independent radial finite-difference check of the O(3) bubble action.
! Build: gfortran -O2 -std=f2008 -Wall -Wextra -fcheck=all \
!        bubble_action_check.f90 -o bubble_action_check.exe
! Run:   bubble_action_check.exe quartic 0.001 0.025
!        bubble_action_check.exe sextic 0.005 0.025
!
! The second argument is the POTENTIAL parameter, not the exact gap.
! The program computes Delta W from the true vacuum before comparing
! S_3/B_TW - 1 with the analytical first-order coefficient.

module bubble_action_model
  use, intrinsic :: iso_fortran_env, only : real64
  implicit none
  private
  integer, parameter, public :: dp = real64
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  public :: potential, potential_prime, potential_second
  public :: true_vacuum, barrier_value, wall_alpha, wall_tension

contains

  pure real(dp) function wall_alpha() result(alpha)
    alpha = log(9.0_dp) / (log(3.0_dp + 2.0_dp * sqrt(3.0_dp)) / sqrt(2.0_dp))
  end function wall_alpha

  pure real(dp) function potential(model, x, parameter) result(w)
    character(len=*), intent(in) :: model
    real(dp), intent(in) :: x, parameter
    real(dp) :: alpha
    select case (trim(model))
    case ('quartic')
      w = 0.25_dp * x**2 * (x**2 - 4.0_dp * (1.0_dp + parameter) * x + 4.0_dp)
    case ('sextic')
      alpha = wall_alpha()
      w = alpha**2 * x**2 * (1.0_dp - x**2)**2 / 8.0_dp - parameter * x**3
    case default
      w = huge(1.0_dp)
    end select
  end function potential

  pure real(dp) function potential_prime(model, x, parameter) result(wp)
    character(len=*), intent(in) :: model
    real(dp), intent(in) :: x, parameter
    real(dp) :: alpha
    select case (trim(model))
    case ('quartic')
      wp = x * (x**2 - 3.0_dp * (1.0_dp + parameter) * x + 2.0_dp)
    case ('sextic')
      alpha = wall_alpha()
      wp = alpha**2 * (x - 4.0_dp*x**3 + 3.0_dp*x**5) / 4.0_dp &
           - 3.0_dp * parameter * x**2
    case default
      wp = huge(1.0_dp)
    end select
  end function potential_prime

  pure real(dp) function potential_second(model, x, parameter) result(wpp)
    character(len=*), intent(in) :: model
    real(dp), intent(in) :: x, parameter
    real(dp) :: alpha
    select case (trim(model))
    case ('quartic')
      wpp = 3.0_dp*x**2 - 6.0_dp*(1.0_dp + parameter)*x + 2.0_dp
    case ('sextic')
      alpha = wall_alpha()
      wpp = alpha**2 * (1.0_dp - 12.0_dp*x**2 + 15.0_dp*x**4) / 4.0_dp &
            - 6.0_dp * parameter * x
    case default
      wpp = huge(1.0_dp)
    end select
  end function potential_second

  real(dp) function true_vacuum(model, parameter) result(tv)
    character(len=*), intent(in) :: model
    real(dp), intent(in) :: parameter
    real(dp) :: low, high, middle
    integer :: j
    if (trim(model) == 'quartic') then
      tv = (3.0_dp*(1.0_dp + parameter) &
           + sqrt(9.0_dp*(1.0_dp + parameter)**2 - 8.0_dp)) / 2.0_dp
      return
    end if
    low = 1.0_dp
    high = 2.0_dp
    do j = 1, 100
      middle = (low + high) / 2.0_dp
      if (potential_prime(model, middle, parameter) < 0.0_dp) then
        low = middle
      else
        high = middle
      end if
    end do
    tv = (low + high) / 2.0_dp
  end function true_vacuum

  real(dp) function barrier_value(model, parameter) result(b)
    character(len=*), intent(in) :: model
    real(dp), intent(in) :: parameter
    real(dp) :: low, high, middle
    integer :: j
    if (trim(model) == 'quartic') then
      b = (3.0_dp*(1.0_dp + parameter) &
           - sqrt(9.0_dp*(1.0_dp + parameter)**2 - 8.0_dp)) / 2.0_dp
      return
    end if
    low = 0.1_dp
    high = 0.8_dp
    do j = 1, 100
      middle = (low + high) / 2.0_dp
      if (potential_prime(model, middle, parameter) > 0.0_dp) then
        low = middle
      else
        high = middle
      end if
    end do
    b = (low + high) / 2.0_dp
  end function barrier_value

  pure real(dp) function wall_tension(model) result(sigma)
    character(len=*), intent(in) :: model
    if (trim(model) == 'quartic') then
      sigma = 4.0_dp / (3.0_dp * sqrt(2.0_dp))
    else
      sigma = wall_alpha() / 8.0_dp
    end if
  end function wall_tension

end module bubble_action_model

module radial_solver
  use bubble_action_model
  implicit none
  private
  public :: solve_bubble, compute_action, half_height_radius

contains

  subroutine equations(model, parameter, h, phi, f, lower, diagonal, upper)
    character(len=*), intent(in) :: model
    real(dp), intent(in) :: parameter, h, phi(0:)
    real(dp), intent(out) :: f(0:), lower(0:), diagonal(0:), upper(0:)
    real(dp) :: r, lo, hi
    integer :: i, n
    n = ubound(phi, 1)
    lower = 0.0_dp
    diagonal = 0.0_dp
    upper = 0.0_dp
    f(0) = 6.0_dp*(phi(1)-phi(0))/h**2 &
           - potential_prime(model, phi(0), parameter)
    diagonal(0) = -6.0_dp/h**2 - potential_second(model, phi(0), parameter)
    upper(0) = 6.0_dp/h**2
    do i = 1, n-1
      r = real(i, dp)*h
      lo = 1.0_dp/h**2 - 1.0_dp/(r*h)
      hi = 1.0_dp/h**2 + 1.0_dp/(r*h)
      f(i) = lo*phi(i-1) - 2.0_dp*phi(i)/h**2 + hi*phi(i+1) &
             - potential_prime(model, phi(i), parameter)
      lower(i) = lo
      diagonal(i) = -2.0_dp/h**2 - potential_second(model, phi(i), parameter)
      upper(i) = hi
    end do
    f(n) = phi(n)
    diagonal(n) = 1.0_dp
  end subroutine equations

  subroutine tridiagonal(lower, diagonal, upper, rhs, x)
    real(dp), intent(in) :: lower(0:), diagonal(0:), upper(0:), rhs(0:)
    real(dp), intent(out) :: x(0:)
    real(dp), allocatable :: c(:), d(:)
    real(dp) :: pivot
    integer :: i, n
    n = ubound(rhs, 1)
    allocate(c(0:n), d(0:n))
    pivot = diagonal(0)
    if (abs(pivot) < 1.0e-14_dp) error stop 'zero tridiagonal pivot'
    c(0) = upper(0)/pivot
    d(0) = rhs(0)/pivot
    do i = 1, n
      pivot = diagonal(i) - lower(i)*c(i-1)
      if (abs(pivot) < 1.0e-14_dp) error stop 'zero tridiagonal pivot'
      c(i) = upper(i)/pivot
      d(i) = (rhs(i) - lower(i)*d(i-1))/pivot
    end do
    x(n) = d(n)
    do i = n-1, 0, -1
      x(i) = d(i) - c(i)*x(i+1)
    end do
  end subroutine tridiagonal

  subroutine solve_bubble(model, parameter, h, phi, n, rmax, tv, gap, &
                          iterations, final_norm)
    character(len=*), intent(in) :: model
    real(dp), intent(in) :: parameter, h
    real(dp), allocatable, intent(out) :: phi(:)
    integer, intent(out) :: n, iterations
    real(dp), intent(out) :: rmax, tv, gap, final_norm
    real(dp), allocatable :: trial(:), f(:), lower(:), diagonal(:), upper(:), step(:)
    real(dp) :: sigma, rtw, r, exponent, factor, norm_old, norm_new, barrier
    integer :: i, j
    tv = true_vacuum(model, parameter)
    barrier = barrier_value(model, parameter)
    gap = -potential(model, tv, parameter)
    if (gap <= 0.0_dp) error stop 'nonpositive vacuum gap'
    sigma = wall_tension(model)
    rtw = 2.0_dp*sigma/gap
    n = ceiling((rtw + 25.0_dp)/h)
    if (mod(n, 2) /= 0) n = n + 1
    rmax = real(n, dp)*h
    allocate(phi(0:n), trial(0:n), f(0:n), lower(0:n), &
             diagonal(0:n), upper(0:n), step(0:n))
    do i = 0, n
      r = real(i, dp)*h
      if (trim(model) == 'quartic') then
        exponent = sqrt(2.0_dp)*(r-rtw-1.5_dp)
        phi(i) = tv/(1.0_dp + (sqrt(2.0_dp)-1.0_dp)*exp(exponent))
      else
        exponent = wall_alpha()*(r-rtw-1.3_dp)
        phi(i) = tv/sqrt(1.0_dp + exp(exponent))
      end if
    end do
    phi(n) = 0.0_dp
    do j = 1, 90
      call equations(model, parameter, h, phi, f, lower, diagonal, upper)
      norm_old = maxval(abs(f))
      if (norm_old < 2.0e-9_dp) exit
      call tridiagonal(lower, diagonal, upper, -f, step)
      factor = 1.0_dp
      do while (factor >= 1.0_dp/8192.0_dp)
        trial = phi + factor*step
        if (trial(0) <= barrier .or. any(trial < -1.0e-10_dp) .or. &
            any(trial > tv*1.1_dp)) then
          factor = factor/2.0_dp
          cycle
        end if
        call equations(model, parameter, h, trial, f, lower, diagonal, upper)
        norm_new = maxval(abs(f))
        if (norm_new < norm_old) exit
        factor = factor/2.0_dp
      end do
      if (factor < 1.0_dp/8192.0_dp) error stop 'Newton line search stalled'
      phi = trial
    end do
    call equations(model, parameter, h, phi, f, lower, diagonal, upper)
    final_norm = maxval(abs(f))
    iterations = min(j-1, 90)
    if (final_norm >= 2.0e-9_dp) then
      write(*,'(a,i0,a,es12.4)') 'Newton steps=', iterations, ', residual=', final_norm
      error stop 'Newton did not converge'
    end if
    if (any(phi(1:n) > phi(0:n-1) + 1.0e-9_dp)) &
      error stop 'computed profile is not decreasing'
  end subroutine solve_bubble

  real(dp) function half_height_radius(phi, h) result(radius)
    real(dp), intent(in) :: phi(0:), h
    real(dp) :: target, fraction
    integer :: i, n
    n = ubound(phi, 1)
    target = phi(0)/sqrt(2.0_dp)
    do i = 1, n
      if (phi(i) <= target) then
        fraction = (phi(i-1)-target)/(phi(i-1)-phi(i))
        radius = (real(i-1, dp)+fraction)*h
        return
      end if
    end do
    error stop 'half-height not crossed'
  end function half_height_radius

  subroutine compute_action(model, parameter, phi, h, kinetic, potential_part)
    character(len=*), intent(in) :: model
    real(dp), intent(in) :: parameter, phi(0:), h
    real(dp), intent(out) :: kinetic, potential_part
    real(dp), allocatable :: derivative(:)
    real(dp) :: r, weight
    integer :: i, n
    n = ubound(phi, 1)
    allocate(derivative(0:n))
    derivative(0) = 0.0_dp
    ! The outer endpoint lies 25 radial units beyond R_TW; its kinetic
    ! contribution is exponentially small at the tested parameters.
    derivative(n) = 0.0_dp
    do i = 1, n-1
      derivative(i) = (phi(i+1)-phi(i-1))/(2.0_dp*h)
    end do
    kinetic = 0.0_dp
    potential_part = 0.0_dp
    do i = 0, n
      r = real(i, dp)*h
      if (i == 0 .or. i == n) then
        weight = 1.0_dp
      else if (mod(i, 2) == 0) then
        weight = 2.0_dp
      else
        weight = 4.0_dp
      end if
      kinetic = kinetic + weight*r**2*derivative(i)**2/2.0_dp
      potential_part = potential_part + weight*r**2*potential(model, phi(i), parameter)
    end do
    kinetic = 4.0_dp*pi*h*kinetic/3.0_dp
    potential_part = 4.0_dp*pi*h*potential_part/3.0_dp
  end subroutine compute_action

end module radial_solver

program bubble_action_check
  use bubble_action_model
  use radial_solver
  implicit none
  character(len=128) :: model, arg, filename
  real(dp), allocatable :: phi(:)
  real(dp) :: parameter, h, rmax, tv, gap, radius, sigma, btw, &
              kinetic, potential_part, action, action_derrick, &
              relative, prediction, discrepancy, derrick, solver_norm
  integer :: argc, ios, n, iterations, i, unit
  argc = command_argument_count()
  if (argc < 3 .or. argc > 4) then
    print '(a)', 'Usage: bubble_action_check MODEL PARAMETER H [profile.csv]'
    print '(a)', 'MODEL = quartic or sextic; PARAMETER is epsilon or sextic delta.'
    error stop 'invalid command line'
  end if
  call get_command_argument(1, model)
  if (trim(model) /= 'quartic' .and. trim(model) /= 'sextic') &
    error stop 'MODEL must be quartic or sextic'
  call get_command_argument(2, arg)
  read(arg, *, iostat=ios) parameter
  if (ios /= 0) error stop 'invalid PARAMETER'
  if (parameter <= 0.0_dp) error stop 'invalid PARAMETER'
  call get_command_argument(3, arg)
  read(arg, *, iostat=ios) h
  if (ios /= 0) error stop 'invalid H'
  if (h <= 0.0_dp) error stop 'invalid H'

  call solve_bubble(trim(model), parameter, h, phi, n, rmax, tv, gap, &
                    iterations, solver_norm)
  call compute_action(trim(model), parameter, phi, h, kinetic, potential_part)
  radius = half_height_radius(phi, h)
  sigma = wall_tension(trim(model))
  btw = 16.0_dp*pi*sigma**3/(3.0_dp*gap**2)
  action = kinetic + potential_part
  action_derrick = 2.0_dp*kinetic/3.0_dp
  relative = action/btw - 1.0_dp
  if (trim(model) == 'quartic') then
    prediction = 27.0_dp*gap/8.0_dp
  else
    prediction = 24.0_dp*(3.0_dp-2.0_dp*log(2.0_dp)) &
                 *gap/wall_alpha()**2
  end if
  discrepancy = relative - prediction
  derrick = abs(kinetic+3.0_dp*potential_part) &
            /(abs(kinetic)+3.0_dp*abs(potential_part))
  write(*,'(a,a)') 'model = ', trim(model)
  write(*,'(a,es24.15)') 'potential parameter = ', parameter
  write(*,'(a,es24.15)') 'grid spacing h = ', h
  write(*,'(a,i0,a,es24.15)') 'grid intervals = ', n, ', rmax = ', rmax
  write(*,'(a,i0)') 'Newton steps = ', iterations
  write(*,'(a,es24.15)') 'maximum equation residual = ', solver_norm
  write(*,'(a,es24.15)') 'true vacuum = ', tv
  write(*,'(a,es24.15)') 'exact gap Delta W = ', gap
  write(*,'(a,es24.15)') 'half-height radius R = ', radius
  write(*,'(a,es24.15)') 'kinetic T = ', kinetic
  write(*,'(a,es24.15)') 'potential V = ', potential_part
  write(*,'(a,es24.15)') 'action T+V = ', action
  write(*,'(a,es24.15)') 'action 2T/3 = ', action_derrick
  write(*,'(a,es24.15)') 'thin-wall B_TW = ', btw
  write(*,'(a,es24.15)') 'observed B/B_TW - 1 = ', relative
  write(*,'(a,es24.15)') 'first-order prediction = ', prediction
  write(*,'(a,es24.15)') 'difference = ', discrepancy
  write(*,'(a,es24.15)') 'difference / (Delta W)^2 = ', discrepancy/gap**2
  write(*,'(a,es24.15)') 'Derrick residual = ', derrick

  if (argc == 4) then
    call get_command_argument(4, filename)
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'cannot open profile output'
    write(unit,'(a)') 'r,phi'
    do i = 0, n
      write(unit,'(es24.15,a,es24.15)') real(i, dp)*h, ',', phi(i)
    end do
    close(unit)
  end if
end program bubble_action_check
