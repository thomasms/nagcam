module RHS_m
    use types_m
    implicit none
    private

    public :: func

contains

    function func(j, x) result (d)
        implicit none

        integer, intent(in)                       :: j
        real (kind=dp), dimension(:), intent(in)  :: x
        real (kind=dp)  :: d

        d = 0.0e+00_dp
    end function

end module RHS_m
