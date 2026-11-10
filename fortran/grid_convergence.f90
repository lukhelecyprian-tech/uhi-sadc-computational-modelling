
program grid_convergence
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: ncases = 3

    integer :: grids(ncases)
    integer :: k, n, nx, ny, nt
    integer :: i, j, step

    real(dp), parameter :: alpha = 1.0d-5
    real(dp), parameter :: L = 50.0d0
    real(dp), parameter :: final_time = 200000.0d0
    real(dp), parameter :: pi = acos(-1.0d0)

    real(dp) :: dx, dy, dt, rx, ry
    real(dp) :: x, y, decay
    real(dp) :: rmse(ncases)
    real(dp) :: max_error(ncases)
    real(dp) :: order

    real(dp), allocatable :: T(:,:), Tnew(:,:), exact(:,:)

    grids = [26, 51, 101]

    print *, "======================================"
    print *, "2D HEAT EQUATION: GRID CONVERGENCE"
    print *, "======================================"

    do k = 1, ncases

        n = grids(k)
        nx = n
        ny = n

        dx = L / real(nx-1, dp)
        dy = L / real(ny-1, dp)

        ! Choose dt proportional to dx squared.
        ! For the selected grids this gives:
        ! 400 s, 100 s and 25 s.
        dt = 100.0d0 * dx**2

        nt = nint(final_time / dt)

        rx = alpha * dt / dx**2
        ry = alpha * dt / dy**2

        if (rx + ry > 0.5d0) then
            print *, "ERROR: Unstable time step."
            stop 1
        end if

        allocate(T(nx,ny), Tnew(nx,ny), exact(nx,ny))

        ! Initialise the analytical sine-mode profile
        do j = 1, ny
            y = real(j-1,dp)*dy

            do i = 1, nx
                x = real(i-1,dp)*dx

                T(i,j) = 20.0d0 + 10.0d0 * &
                    sin(pi*x/L)*sin(pi*y/L)
            end do
        end do

        ! Explicit finite-difference time integration
        do step = 1, nt

            Tnew = T

            do j = 2, ny-1
                do i = 2, nx-1

                    Tnew(i,j) = T(i,j) &
                        + rx*(T(i+1,j)-2.0d0*T(i,j) &
                        +T(i-1,j)) &
                        + ry*(T(i,j+1)-2.0d0*T(i,j) &
                        +T(i,j-1))

                end do
            end do

            ! Fixed temperature at boundaries
            Tnew(1,:) = 20.0d0
            Tnew(nx,:) = 20.0d0
            Tnew(:,1) = 20.0d0
            Tnew(:,ny) = 20.0d0

            T = Tnew

        end do

        ! Analytical solution
        decay = exp(-2.0d0*alpha*pi**2 * &
                    final_time/L**2)

        do j = 1, ny
            y = real(j-1,dp)*dy

            do i = 1, nx
                x = real(i-1,dp)*dx

                exact(i,j) = 20.0d0 + 10.0d0 * &
                    sin(pi*x/L)*sin(pi*y/L)*decay

            end do
        end do

        rmse(k) = sqrt(sum((T-exact)**2) / &
                  real(nx*ny,dp))

        max_error(k) = maxval(abs(T-exact))

        print *
        print *, "Grid:", nx, "x", ny
        print *, "dx =", dx, "m"
        print *, "dt =", dt, "s"
        print *, "Time steps =", nt
        print *, "RMSE =", rmse(k)
        print *, "Maximum error =", max_error(k)

        deallocate(T, Tnew, exact)

    end do

    print *
    print *, "======================================"
    print *, "OBSERVED CONVERGENCE ORDER"
    print *, "======================================"

    do k = 1, ncases-1

        order = log(rmse(k)/rmse(k+1)) / log(2.0d0)

        print *, "Grid", grids(k), "to", grids(k+1)
        print *, "Estimated order =", order

    end do

    ! Save numerical errors for Python plotting
    open(unit=10, file="results/convergence.csv", &
         status="replace", action="write")

    write(10,'(A)') "grid,dx,rmse,max_error"

    do k = 1, ncases

        dx = L / real(grids(k)-1,dp)

        write(10,'(I0,",",ES16.8,",",ES16.8,",",ES16.8)') &
            grids(k), dx, rmse(k), max_error(k)

    end do

    close(10)

    print *
    print *, "Results saved to results/convergence.csv"

end program grid_convergence

