module IO_m
    use types_m
    implicit none
    private

    public :: r8vec_linspace
    public :: r8vec_write
    public :: r8mat_write

contains

    subroutine r8vec_linspace( a_first, a_last, a)

        implicit none

        real (kind=dp), intent(in) :: a_first
        real (kind=dp), intent(in) :: a_last
        real (kind=dp), dimension(:), intent(inout) :: a
        integer :: i

        do i = lbound(a, 1), ubound(a, 1)
            a(i) = (real(ubound(a, 1)-i,kind=dp)*a_first+real(i-1,kind=dp)*a_last)/ &
              real(ubound(a, 1)-1, kind=dp)
        end do
    end subroutine

    subroutine r8vec_write(output_filename, x)

        implicit none

        character (len=*), intent(in) :: output_filename
        real (kind=dp), dimension(:), intent(inout) :: x

        integer :: m
        integer :: j
        integer :: output_unit_id

        output_unit_id = 11
        open (unit=output_unit_id, file=output_filename, status='replace')

        do j = lbound(x, 1), ubound(x, 1)
            write (output_unit_id, '(2x,g24.16)') x(j)
        end do

        close (unit=output_unit_id)
    end subroutine


    subroutine r8mat_write(output_filename, table)
        implicit none

        character (len=*), intent(in) :: output_filename
        real (kind=dp), dimension(:, :), intent(inout) :: table

        integer :: j, m, n
        integer :: output_unit_id
        character (len=30) :: string

        output_unit_id = 10
        open (unit=output_unit_id, file=output_filename, status='replace')

        m = size( table(:, :), 1 )
        n = size( table(:, :), 2 )
        write (string, '(a1,i8,a1,i8,a1,i8,a1)') '(', m, 'g', 24, '.', 16, ')'

        do j = 1, n
            write (output_unit_id, string) table(1:m, j)
        end do

        close (unit=output_unit_id)
    end subroutine

end module IO_m
