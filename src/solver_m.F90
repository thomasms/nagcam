module solver_m
    use types_m
    use RHS_m
    implicit none
    private

    public :: fd1d_heat_explicit

contains

    subroutine fd1d_heat_explicit(x, t, dt, cfl, h, h_new)
        implicit none

        real (kind=dp), dimension(:), contiguous, intent(in) :: x
        real (kind=dp), intent(in)   :: t
        real (kind=dp), intent(in)   :: dt
        real (kind=dp), intent(in)   :: cfl
        real (kind=dp), dimension(:), contiguous, intent(in)   :: h
        real (kind=dp), dimension(:), contiguous, intent(inout):: h_new

        integer :: j
        real (kind=dp) :: f

        h_new(1) = 0.0e+00_dp

        do j = 2, ubound( x, 1 ) - 1
            f = func(j, x)
            h_new(j) = h(j) + dt*f + cfl*(h(j-1)-2.0e+00_dp*h(j)+h(j+1))
        end do

        ! set the boundary conditions again
        h_new(1) = 90.0e+00_dp
        h_new(ubound( x, 1 )) = 70.0e+00_dp
    end subroutine

end module solver_m
