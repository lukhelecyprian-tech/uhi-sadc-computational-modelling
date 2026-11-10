
program verify_heat2d
    implicit none

    integer, parameter :: nx = 51, ny = 51
    integer, parameter :: nt = 200
    integer :: i, j, n

    real(8), parameter :: alpha = 1.0d-5
    real(8), parameter :: L = 50.0d0
    real(8), parameter :: dt = 1000.0d0
    real(8), parameter :: pi = acos(-1.0d0)

    real(8), parameter :: dx = L / (nx-1)
    real(8), parameter :: dy = L / (ny-1)

    real(8) :: T(nx,ny), Tnew(nx,ny)
    real(8) :: exact(nx,ny)
    real(8) :: x, y, time
    real(8) :: rx, ry, rmse, max_error

    rx = alpha * dt / dx**2
    ry = alpha * dt / dy**2

    if (rx + ry > 0.5d0) then
        print *, "ERROR: Unstable time step."
        stop 1
    end if

    ! Initial condition with known analytical solution
    do j = 1, ny
        y = (j-1)*dy
        do i = 1, nx
            x = (i-1)*dx

            T(i,j) = 20.0d0 + 10.0d0 * &
                sin(pi*x/L)*sin(pi*y/L)
        end do
    end do

    ! Explicit finite-difference time stepping
    do n = 1, nt
        Tnew = T

        do j = 2, ny-1
            do i = 2, nx-1
                Tnew(i,j) = T(i,j) &
                    + rx*(T(i+1,j)-2*T(i,j)+T(i-1,j)) &
                    + ry*(T(i,j+1)-2*T(i,j)+T(i,j-1))
            end do
        end do

        ! Dirichlet boundary conditions
        Tnew(1,:) = 20.0d0
        Tnew(nx,:) = 20.0d0
        Tnew(:,1) = 20.0d0
        Tnew(:,ny) = 20.0d0

        T = Tnew
    end do

    ! Exact solution at the final simulated time
    time = nt*dt

    do j = 1, ny
        y = (j-1)*dy
        do i = 1, nx
            x = (i-1)*dx

            exact(i,j) = 20.0d0 + 10.0d0 * &
                sin(pi*x/L)*sin(pi*y/L) * &
                exp(-2.0d0*alpha*pi**2*time/L**2)
        end do
    end do

    ! Numerical error measurements
    rmse = sqrt(sum((T-exact)**2)/real(nx*ny,8))
    max_error = maxval(abs(T-exact))

    print *, "--------------------------------"
    print *, "2D HEAT EQUATION VERIFICATION"
    print *, "--------------------------------"
    print *, "Grid:", nx, "x", ny
    print *, "Final time (s):", time
    print *, "Stability sum:", rx+ry
    print *, "Numerical maximum (C):", maxval(T)
    print *, "Analytical maximum (C):", maxval(exact)
    print *, "RMSE (C):", rmse
    print *, "Maximum absolute error (C):", max_error

    call save_field("results/verification_numerical.csv", T)
    call save_field("results/verification_exact.csv", exact)

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
                    write(unit,'(F14.8,A)',advance="no") &
                        field(col,row), ","
                else
                    write(unit,'(F14.8)') field(col,row)
                end if
            end do
        end do

        close(unit)
    end subroutine save_field

end program verify_heat2d

