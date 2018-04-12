module types_m
    use, intrinsic :: iso_fortran_env
    implicit none
    private

    !> Single precision - [32 bit] / [4 bytes]
    integer, parameter, public :: sp = real32

    !> Double precision - [64 bit] / [8 bytes]
    integer, parameter, public :: dp = real64

    !> Quad precision - [128 bit] / [16 bytes]
    integer, parameter, public :: qp = real128

end module types_m
