
program mbabane_coupled_spatial_diagnostic

    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: ngrids = 4
    integer, parameter :: sizes(ngrids) = [20, 40, 80, 160]

    real(dp), parameter :: pi = acos(-1.0d0)
    real(dp), parameter :: length = 1200.0d0
    real(dp), parameter :: t0 = 293.15d0
    real(dp), parameter :: amplitude = 5.0d0
    real(dp), parameter :: tau = 21600.0d0

    real(dp), parameter :: ca = 2.0d6
    real(dp), parameter :: kappa = 1.0d0
    real(dp), parameter :: sigma = 5.670374419d-8
    real(dp), parameter :: emissivity = 0.95d0

    real(dp), parameter :: solar = 500.0d0
    real(dp), parameter :: albedo = 0.12d0
    real(dp), parameter :: tair = 293.15d0
    real(dp), parameter :: tsky = 283.15d0
    real(dp), parameter :: h = 10.0d0
    real(dp), parameter :: latent = 50.0d0

    real(dp), parameter :: final_time = 21600.0d0
    real(dp), parameter :: dt_target = 0.1d0

    integer :: g, nx, ny, i, j, n, nt, ios
    real(dp) :: dx, dy, dt, x, y, time
    real(dp) :: laplacian, net_flux, source
    real(dp) :: exact_value, difference
    real(dp) :: rms_error, max_error
    real(dp) :: errors(ngrids), order
    real(dp) :: stability

    real(dp), allocatable :: temp(:,:), temp_new(:,:)

    open(20, file="results/mbabane_coupled_spatial_diagnostic_dt01.csv", &
         status="replace", action="write", iostat=ios)

    if (ios /= 0) error stop "Cannot create CSV file"

    write(20,'(A)') "grid_n,dx_m,dt_s,rmse_C,max_error_C"

    print *
    print *, "=========================================="
    print *, "COUPLED MODEL SPATIAL CONVERGENCE TEST"
    print *, "=========================================="
    print *, "Manufactured exact solution"
    print *, "Zero-flux boundaries"
    print *
    print *, "Grid       dx(m)        RMSE(C)        Max(C)"
    print *, "-----------------------------------------------"

    do g = 1, ngrids

        nx = sizes(g)
        ny = nx

        dx = length / real(nx,dp)
        dy = dx

        allocate(temp(nx,ny), temp_new(nx,ny))

        nt = ceiling(final_time/dt_target)
        dt = final_time / real(nt,dp)

        stability = kappa*dt*(1d0/dx**2+1d0/dy**2)

        if (stability > 0.5d0) &
            error stop "Explicit diffusion stability failure"

        ! ---------------------------------------------
        ! Initialise from manufactured exact solution
        ! ---------------------------------------------

        do j = 1, ny
            y = (real(j,dp)-0.5d0)*dy

            do i = 1, nx
                x = (real(i,dp)-0.5d0)*dx

                temp(i,j) = exact_temperature(x,y,0d0)
            end do
        end do

        ! ---------------------------------------------
        ! Coupled numerical integration
        ! ---------------------------------------------

        do n = 1, nt

            time = real(n-1,dp)*dt
            temp_new = temp

            do j = 1, ny
                y = (real(j,dp)-0.5d0)*dy

                do i = 1, nx
                    x = (real(i,dp)-0.5d0)*dx

                    laplacian = 0d0

                    if (i > 1) then
                        laplacian = laplacian + &
                            (temp(i-1,j)-temp(i,j))/dx**2
                    end if

                    if (i < nx) then
                        laplacian = laplacian + &
                            (temp(i+1,j)-temp(i,j))/dx**2
                    end if

                    if (j > 1) then
                        laplacian = laplacian + &
                            (temp(i,j-1)-temp(i,j))/dy**2
                    end if

                    if (j < ny) then
                        laplacian = laplacian + &
                            (temp(i,j+1)-temp(i,j))/dy**2
                    end if

                    net_flux = (1d0-albedo)*solar &
                        - emissivity*sigma * &
                          (temp(i,j)**4-tsky**4) &
                        - h*(temp(i,j)-tair) &
                        - latent

                    source = manufactured_source(x,y,time)

                    temp_new(i,j) = temp(i,j) + dt * &
                        (kappa*laplacian + &
                         (net_flux+source)/ca)

                end do
            end do

            temp = temp_new

        end do

        ! ---------------------------------------------
        ! Compare against exact final temperature
        ! ---------------------------------------------

        rms_error = 0d0
        max_error = 0d0

        do j = 1, ny
            y = (real(j,dp)-0.5d0)*dy

            do i = 1, nx
                x = (real(i,dp)-0.5d0)*dx

                exact_value = exact_temperature(x,y,final_time)
                difference = temp(i,j)-exact_value

                rms_error = rms_error + difference**2
                max_error = max(max_error,abs(difference))
            end do
        end do

        rms_error = sqrt(rms_error/real(nx*ny,dp))
        errors(g) = rms_error

        print '(I5,4X,F10.4,4X,ES13.5,4X,ES13.5)', &
            nx, dx, rms_error, max_error

        write(20,'(I6,A,F12.6,A,F10.4,A,ES18.10,A,ES18.10)') &
            nx,",",dx,",",dt,",",rms_error,",",max_error

        deallocate(temp,temp_new)

    end do

    close(20)

    print *
    print *, "Estimated convergence orders"
    print *, "------------------------------------------"

    do g = 1, ngrids-1

        if (errors(g) > 0d0 .and. errors(g+1) > 0d0) then

            order = log(errors(g)/errors(g+1))/log(2d0)

            print '(I5,A,I5,A,F12.6)', &
                sizes(g)," -> ",sizes(g+1),": ",order

        else
            print *, "Cannot calculate convergence order"
        end if

    end do

    print *
    print *, "Saved spatial convergence CSV."
    print *, "Test completed."

contains

    function exact_temperature(x,y,t) result(value)

        real(dp), intent(in) :: x,y,t
        real(dp) :: value

        value = t0 + amplitude*exp(-t/tau) * &
            cos(pi*x/length)*cos(pi*y/length)

    end function exact_temperature

    function manufactured_source(x,y,t) result(value)

        real(dp), intent(in) :: x,y,t
        real(dp) :: value
        real(dp) :: texact, theta, qexact
        real(dp) :: temporal, spatial

        texact = exact_temperature(x,y,t)
        theta = texact-t0

        temporal = -theta/tau

        spatial = -2d0*pi*pi*theta/length**2

        qexact = (1d0-albedo)*solar &
            - emissivity*sigma*(texact**4-tsky**4) &
            - h*(texact-tair) &
            - latent

        value = ca*(temporal-kappa*spatial)-qexact

    end function manufactured_source

end program mbabane_coupled_spatial_diagnostic
