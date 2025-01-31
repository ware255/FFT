program main
    implicit none
    complex(8) :: x(8) = [(1,0), (1,0), (1,0), (1,0), (0,0), (0,0), (0,0), (0,0)]
    integer(8) i
    call fft(x, .false.)

    do i = 1, 8
        print '(F8.3, F8.3, A)', real(x(i)), aimag(x(i)), "i"
    end do

contains
    recursive subroutine fft(x, invert)
        implicit none
        complex(8), intent(inout) :: x(:)
        logical, intent(in) :: invert
        complex(8), allocatable :: even(:), odd(:)
        real(8), parameter :: PI = 4 * atan(1.0)
        complex(8) w, wn
        real(8) ang
        integer(8) i, n
        n = ubound(x, 1)
        if (n <= 1) return

        allocate(even(n/2))
        allocate(odd(n/2))

        do i = 1, n / 2
            even(i) = x(i * 2)
            odd(i) = x(i * 2 - 1)
        end do
        call fft(even, invert)
        call fft(odd, invert)

        if (invert) then
            ang = -2 * PI / n
        else
            ang = 2 * PI / n
        end if
        w = (1.0, 0.0)
        wn = cmplx(cos(ang), sin(ang), 8)
        do i = 1, n / 2
            x(i) = even(i) + w * odd(i)
            x(i + n/2) = even(i) - w * odd(i)
            if (invert) then
                x(i) = x(i) / 2
                x(i + n/2) = x(i + n/2) / 2
            end if
            w = w * wn
        end do

        deallocate(odd)
        deallocate(even)
    end subroutine fft
end program main
