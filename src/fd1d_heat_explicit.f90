program fd1d_heat_explicit_prb
    use types_m
    use CFL_m
    use IO_m
    use solver_m

    implicit none

    integer, parameter :: T_NUM = 201
    integer, parameter :: X_NUM = 21

    integer :: stat

    real (kind=dp) :: cfl
    real (kind=dp) :: dt

    real (kind=dp), dimension(:), allocatable   :: h
    real (kind=dp), dimension(:), allocatable   :: h_new
    real (kind=dp), dimension(:,:), allocatable :: hmat
    real (kind=dp), dimension(:), allocatable   :: t
    real (kind=dp), dimension(:), allocatable   :: x

! the "matrix" stores all x-values for all t-values
! remember Fortran is column major, meaning that rows are contiguous
    integer :: i
    integer :: j
    real (kind=dp) :: k

    real (kind=dp) :: t_max
    real (kind=dp) :: t_min
    real (kind=dp) :: x_max
    real (kind=dp) :: x_min

    write (*, '(a)') ' '
    write (*, '(a)') 'FD1D_HEAT_EXPLICIT_PRB:'
    write (*, '(a)') '  FORTRAN77 version.'
    write (*, '(a)') '  Test the FD1D_HEAT_EXPLICIT library.'

    write (*, '(a)') ' '
    write (*, '(a)') 'FD1D_HEAT_EXPLICIT_PRB:'
    write (*, '(a)') '  Normal end of execution.'
    write (*, '(a)') ' '

    write (*, '(a)') ' '
    write (*, '(a)') 'FD1D_HEAT_EXPLICIT_TEST01:'
    write (*, '(a)') '  Compute an approximate solution to the time-dependent'
    write (*, '(a)') '  one dimensional heat equation:'
    write (*, '(a)') ' '
    write (*, '(a)') '    dH/dt - K * d2H/dx2 = f(x,t)'
    write (*, '(a)') ' '
    write (*, '(a)') '  Run a simple test case.'

    ! Set the status to zero by default
    stat = 0

    ! do the allocation for all global arrays
    allocate(h(1:X_NUM), h_new(1:X_NUM), hmat(1:X_NUM, 1:T_NUM), x(1:X_NUM), t(1:T_NUM), stat = stat)
    if(stat /= 0)then
        write (*, '(a)') 'Failed to allocate arrays!'
        stop 1
    end if

! heat coefficient
    k = 0.002e+00_dp

! the x-range values
    x_min = 0.0e+00_dp
    x_max = 1.0e+00_dp
! x_num is the number of intervals in the x-direction
    call r8vec_linspace(x_min, x_max, x)

! the t-range values. integrate from t_min to t_max
    t_min = 0.0e+00_dp
    t_max = 80.0e+00_dp

! t_num is the number of intervals in the t-direction
    dt = (t_max-t_min)/real(t_num-1, kind=dp)
    call r8vec_linspace(t_min, t_max, t)

! get the CFL coefficient
    call fd1d_heat_explicit_cfl(k, T_NUM, t_min, t_max, X_NUM, x_min, x_max, &
      cfl)

    if (0.5e+00_dp<=cfl) then
        write (*, '(a)') ' '
        write (*, '(a)') 'FD1D_HEAT_EXPLICIT_CFL - Fatal error!'
        write (*, '(a)') '  CFL condition failed.'
        write (*, '(a)') '  0.5 <= K * dT / dX / dX = CFL.'
        stop
    end if

! set the initial condition
    do j = 1, x_num
        h(j) = 50.0e+00_dp
    end do

! set the bounday condition
    h(1) = 90.0e+00_dp
    h(x_num) = 70.0e+00_dp

! initialise the matrix to the initial condition
    do i = 1, x_num
        hmat(i, 1) = h(i)
    end do

! the main time integration loop 
    do j = 2, t_num
        call fd1d_heat_explicit(x, t(j-1), dt, cfl, h, h_new)

        do i = 1, x_num
            hmat(i, j) = h_new(i)
            h(i) = h_new(i)
        end do
    end do

! write data to files
    call r8mat_write('h_test01.txt', hmat)
    call r8vec_write('t_test01.txt', t)
    call r8vec_write('x_test01.txt', x)

    deallocate(h, h_new, hmat, x, t)

end program
