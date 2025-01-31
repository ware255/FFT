program main
    use, intrinsic :: ieee_arithmetic
    implicit none
    complex(8) :: x(8) = [(1,0), (1,0), (1,0), (1,0), (0,0), (0,0), (0,0), (0,0)]
    integer(8) i

    call fft(x, .false.)
    do i = 1, 8
        print '(F11.6, F11.6, A)', real(x(i)), aimag(x(i)), "i"
    end do

    call fft(x, .true.)
    print *, ""
    do i = 1, 8
        print '(F11.6, F11.6, A)', real(x(i)), aimag(x(i)), "i"
    end do

contains
    recursive subroutine fft(x, invert)
        implicit none
        complex(8), intent(inout) :: x(:)
        logical, intent(in) :: invert
        complex(8), allocatable :: even(:), odd(:)
        real(8), parameter :: PI = 4.0 * atan(1.0)
        complex(8) w, t, j
        integer(8) i, n
        n = ubound(x, 1)
        if (n <= 1) return

        allocate(even(n / 2))
        allocate(odd(n / 2))

        do i = 1, n / 2
            even(i) = x(2 * i - 1)
            odd(i) = x(2 * i)
        end do

        call fft(even, invert)
        call fft(odd, invert)

        j = cmplx(0.0d0, 1.0d0, 8)
        do i = 1, n / 2
            w = exp(sum(pack([2.0d0], invert, [-2.0d0])) * PI * j * (i - 1) / n)
            t = w * odd(i)
            x(i) = even(i) + t
            x(i + n/2) = even(i) - t
            if (invert) then
                x(i) = x(i) / 2.0d0
                x(i + n/2) = x(i + n/2) / 2.0d0
            end if
        end do

        deallocate(odd)
        deallocate(even)
    end subroutine fft
end program main
