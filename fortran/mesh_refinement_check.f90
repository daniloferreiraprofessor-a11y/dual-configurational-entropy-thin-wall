program mesh_refinement_check
  use, intrinsic :: iso_fortran_env, only: real64, error_unit
  implicit none

  ! Run this program in the directory containing bubble_action_check.exe.
  ! Each run writes its complete stdout to a distinct .log file.
  call check_case('quartic', '0.00025')
  call check_case('sextic',  '0.00125')

contains

  subroutine check_case(model, potential_parameter)
    character(*), intent(in) :: model, potential_parameter
    character(len=9), parameter :: h_arg(3) = &
      [character(len=9) :: '0.0125', '0.00625', '0.003125']
    character(len=6), parameter :: h_tag(3) = &
      [character(len=6) :: '0125', '00625', '003125']
    character(len=256) :: log_file, command
    real(real64) :: q(3), derrick(3), equation_residual(3)
    real(real64) :: gap, first_order, r12, r23, ratio
    integer :: i, exit_status, command_status

    do i = 1, 3
      log_file = 'mesh_' // trim(model) // '_h' // trim(h_tag(i)) // '.log'
      command = 'bubble_action_check.exe ' // trim(model) // ' ' // &
        trim(potential_parameter) // ' ' // trim(h_arg(i)) // &
        ' > "' // trim(log_file) // '"'
      call execute_command_line(trim(command), exitstat=exit_status, &
        cmdstat=command_status)
      if (command_status /= 0 .or. exit_status /= 0) then
        write(error_unit, '(a)') 'Failed: ' // trim(command)
        error stop 1
      end if
      call read_scalar(log_file, 'observed B/B_TW - 1', q(i))
      call read_scalar(log_file, 'Derrick residual', derrick(i))
      call read_scalar(log_file, 'maximum equation residual', equation_residual(i))
    end do

    ! These physical quantities do not depend on grid spacing.
    call read_scalar(log_file, 'exact gap Delta W', gap)
    call read_scalar(log_file, 'first-order prediction', first_order)

    r12 = (4.0_real64 * q(2) - q(1)) / 3.0_real64
    r23 = (4.0_real64 * q(3) - q(2)) / 3.0_real64
    ratio = (q(2) - q(1)) / (q(3) - q(2))

    write(*, '(a)') repeat('=', 72)
    write(*, '(a)') trim(model) // ', potential parameter = ' // trim(potential_parameter)
    write(*, '(a,es23.15)') 'Delta W                         = ', gap
    write(*, '(a,es23.15)') 'First-order Q                   = ', first_order
    do i = 1, 3
      write(*, '(a,a,a,es23.15,a,es11.3,a,es11.3)') &
        'Q(h=', trim(h_arg(i)), ') = ', q(i), &
        '  Derrick=', derrick(i), '  EqRes=', equation_residual(i)
    end do
    write(*, '(a,f14.9)')  'Successive-change ratio        = ', ratio
    write(*, '(a,es23.15)') 'Richardson h1,h2             = ', r12
    write(*, '(a,es23.15)') 'Richardson h2,h3             = ', r23
    write(*, '(a,es23.15)') 'Difference of extrapolations = ', r23 - r12
    write(*, '(a,es23.15)') 'Quadratic quotient h1,h2     = ', &
      (r12 - first_order) / gap**2
    write(*, '(a,es23.15)') 'Quadratic quotient h2,h3     = ', &
      (r23 - first_order) / gap**2
  end subroutine check_case

  subroutine read_scalar(filename, key, value)
    character(*), intent(in) :: filename, key
    real(real64), intent(out) :: value
    character(len=512) :: line
    integer :: file_unit, io_status, equals_position

    open(newunit=file_unit, file=trim(filename), status='old', &
      action='read', iostat=io_status)
    if (io_status /= 0) error stop 'Cannot open a run log.'
    do
      read(file_unit, '(a)', iostat=io_status) line
      if (io_status /= 0) exit
      if (index(line, trim(key)) == 0) cycle
      equals_position = index(line, '=')
      if (equals_position == 0) exit
      read(line(equals_position + 1:), *, iostat=io_status) value
      close(file_unit)
      if (io_status /= 0) error stop 'Cannot parse a value in a run log.'
      return
    end do
    close(file_unit)
    write(error_unit, '(a)') 'Missing key: ' // trim(key) // ' in ' // trim(filename)
    error stop 1
  end subroutine read_scalar
end program mesh_refinement_check
