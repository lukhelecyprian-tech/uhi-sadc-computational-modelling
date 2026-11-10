
program time_convergence
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: nx = 51, ny = 51
    integer, parameter :: ncases = 4

    real(dp), parameter :: alpha = 1.0d-5
    real(dp), parameter :: L = 50.0d0
    real(dp), parameter :: final_time = 200000.0d0
    real(dp), parameter :: pi = acos(-1.0d0)

    real(dp), parameter :: dx = L / real(nx-1, dp)
    real(dp), parameter :: dy = L / real(ny-1, dp)

    real(dp) :: dt_values(ncases)
    real(dp) :: rmse(ncases), max_error(ncases)
    real(dp) :: T(nx,ny), Tnew(nx,ny)
    real(dp) :: reference(nx,ny)

    real(dp) :: x, y, dt, rx, ry
    real(dp) :: lambda_h, amplitude
    real(dp) :: order

    integer :: i, j, k, step, nt, unit

    ! Four time steps for temporal convergence
    dt_values = [1000.0d0, 500.0d0, &
                 250.0d0, 125.0d0]

    ! Discrete Laplacian eigenvalue
    lambda_h = -4.0d0 * ( &
        sin(pi*dx/(2.0d0*L))**2 / dx**2 + &
        sin(pi*dy/(2.0d0*L))**2 / dy**2 )

    ! Exact-in-time reference for the spatial grid
    amplitude = 10.0d0 * &
        exp(alpha*lambda_h*final_time)

    do j = 1, ny
        y = real(j-1,dp)*dy

        do i = 1, nx
            x = real(i-1,dp)*dx

            reference(i,j) = 20.0d0 + amplitude * &
                sin(pi*x/L)*sin(pi*y/L)
        end do
    end do

    print *, "===================================="
    print *, "TEMPORAL CONVERGENCE VERIFICATION"
    print *, "===================================="
    print *, "Grid:", nx, "x", ny
    print *, "Final time (s):", final_time

    ! Run each time-step experiment
    do k = 1, ncases

        dt = dt_values(k)
        nt = nint(final_time/dt)

        rx = alpha*dt/dx**2
        ry = alpha*dt/dy**2

        ! Stability check for explicit scheme
        if (rx + ry > 0.5d0) then
            print *, "ERROR: Unstable time step."
            stop 1
        end if

        ! Initial temperature distribution
        do j = 1, ny
            y = real(j-1,dp)*dy

            do i = 1, nx
                x = real(i-1,dp)*dx

                T(i,j) = 20.0d0 + 10.0d0 * &
                    sin(pi*x/L)*sin(pi*y/L)
            end do
        end do

        ! Forward Euler time integration
        do step = 1, nt

            Tnew = T

            do j = 2, ny-1
                do i = 2, nx-1

                    Tnew(i,j) = T(i,j) &
                        + rx*(T(i+1,j) - &
                        2.0d0*T(i,j) + T(i-1,j)) &
                        + ry*(T(i,j+1) - &
                        2.0d0*T(i,j) + T(i,j-1))

                end do
            end do

            ! Fixed boundary temperature
            Tnew(1,:) = 20.0d0
            Tnew(nx,:) = 20.0d0
            Tnew(:,1) = 20.0d0
            Tnew(:,ny) = 20.0d0

            T = Tnew

        end do

        ! Compare with exact-in-time reference
        rmse(k) = sqrt( &
            sum((T-reference)**2) / &
            real(nx*ny,dp) )

        max_error(k) = maxval(abs(T-reference))

        print *
        print *, "Time step (s):", dt
        print *, "Number of steps:", nt
        print *, "Stability sum:", rx+ry
        print *, "RMSE (C):", rmse(k)
        print *, "Maximum error (C):", max_error(k)

    end do

    print *
    print *, "===================================="
    print *, "OBSERVED TEMPORAL CONVERGENCE ORDER"
    print *, "===================================="

    do k = 1, ncases-1

        order = log(rmse(k)/rmse(k+1)) / &
                log(2.0d0)

        print *, "dt", dt_values(k), &
                 "to", dt_values(k+1)
        print *, "Estimated order:", order

    end do

    ! Save results to CSV for Python plotting
    open(newunit=unit, &
         file="results/time_convergence.csv", &
         status="replace", action="write")

    write(unit,'(A)') "dt,rmse,max_error"

    do k = 1, ncases

        write(unit, &
            '(ES16.8,",",ES16.8,",",ES16.8)') &
            dt_values(k), rmse(k), max_error(k)

    end do

    close(unit)

    print *
    print *, "Results saved to:"
    print *, "results/time_convergence.csv"

end program time_convergence

