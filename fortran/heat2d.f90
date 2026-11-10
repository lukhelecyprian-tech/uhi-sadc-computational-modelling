
program heat2d
    implicit none

    integer, parameter :: nx = 51, ny = 51
    integer, parameter :: nt = 200
    integer :: i, j, n

    real(8), parameter :: alpha = 1.0d-5
    real(8), parameter :: dx = 1.0d0
    real(8), parameter :: dy = 1.0d0
    real(8), parameter :: dt = 1000.0d0

    real(8) :: T(nx,ny), Tnew(nx,ny)
    real(8) :: rx, ry

    ! Numerical stability coefficients
    rx = alpha * dt / dx**2
    ry = alpha * dt / dy**2

    if (rx + ry > 0.5d0) then
        print *, "ERROR: Unstable time step."
        stop 1
    end if

    ! Initial temperature: 20 degrees Celsius
    T = 20.0d0

    ! Hot central square: 40 degrees Celsius
    T(21:31,21:31) = 40.0d0

    ! Save initial temperature field
    call save_field("results/initial.csv", T)

    ! Time integration
    do n = 1, nt
        Tnew = T

        do j = 2, ny-1
            do i = 2, nx-1
                Tnew(i,j) = T(i,j) &
                    + rx*(T(i+1,j)-2*T(i,j)+T(i-1,j)) &
                    + ry*(T(i,j+1)-2*T(i,j)+T(i,j-1))
            end do
        end do

        ! Fixed boundary temperature: 20 C
        Tnew(1,:)  = 20.0d0
        Tnew(nx,:) = 20.0d0
        Tnew(:,1)  = 20.0d0
        Tnew(:,ny) = 20.0d0

        T = Tnew
    end do

    ! Save final temperature field
    call save_field("results/final.csv", T)

    print *, "Simulation completed successfully."
    print *, "Stability sum =", rx + ry
    print *, "Initial maximum temperature = 40 C"
    print *, "Final maximum temperature =", maxval(T)

contains

    subroutine save_field(filename, field)
        character(len=*), intent(in) :: filename
        real(8), intent(in) :: field(nx,ny)
        integer :: unit, row, col

        open(newunit=unit, file=filename, &
             status="replace", action="write")

        do row = 1, ny
            do col = 1, nx
                if (col < nx) then
                    write(unit,'(F12.6,A)',advance="no") &
                        field(col,row), ","
                else
                    write(unit,'(F12.6)') field(col,row)
                end if
            end do
        end do

        close(unit)
    end subroutine save_field

end program heat2d
