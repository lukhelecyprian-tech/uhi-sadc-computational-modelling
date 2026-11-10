
program mbabane_coupled_time_convergence

    use, intrinsic :: ieee_arithmetic
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: nruns = 4

    integer :: nx, ny, i, j, k, n, nt, ios
    integer :: active_count

    integer, allocatable :: mask(:,:)

    real(dp), allocatable :: initial(:,:), temp(:,:)
    real(dp), allocatable :: temp_new(:,:), albedo(:,:)
    real(dp), allocatable :: solutions(:,:,:)

    real(dp) :: dx, dy, dt, duration
    real(dp) :: sigma, emissivity, ca, kappa
    real(dp) :: solar, tair, tsky, h, latent
    real(dp) :: laplacian, net_flux, stability
    real(dp) :: mean_initial, mean_final
    real(dp) :: rmse, diff_norm(3), order
    real(dp) :: stepsizes(nruns)

    character(len=120) :: filename

    ! ========================================================
    ! 1. GRID AND INPUT FILES
    ! ========================================================

    open(10, file="data/mbabane_energy_grid.txt", &
         status="old", action="read", iostat=ios)

    if (ios /= 0) error stop "Cannot open grid file"

    read(10, *, iostat=ios) nx, ny
    if (ios /= 0) error stop "Cannot read dimensions"

    read(10, *, iostat=ios) dx, dy
    if (ios /= 0) error stop "Cannot read grid spacing"

    close(10)

    if (nx < 2 .or. ny < 2) error stop "Invalid dimensions"
    if (dx <= 0d0 .or. dy <= 0d0) error stop "Invalid spacing"

    allocate(initial(nx,ny), temp(nx,ny))
    allocate(temp_new(nx,ny), albedo(nx,ny))
    allocate(mask(nx,ny))
    allocate(solutions(nx,ny,nruns))

    open(11, file="data/mbabane_energy_temperature_K.txt", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Cannot open temperature"

    open(12, file="data/mbabane_energy_albedo.txt", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Cannot open albedo"

    open(13, file="data/mbabane_energy_mask.txt", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Cannot open mask"

    do j = 1, ny

        read(11, *, iostat=ios) (initial(i,j), i=1,nx)
        if (ios /= 0) error stop "Temperature read error"

        read(12, *, iostat=ios) (albedo(i,j), i=1,nx)
        if (ios /= 0) error stop "Albedo read error"

        read(13, *, iostat=ios) (mask(i,j), i=1,nx)
        if (ios /= 0) error stop "Mask read error"

    end do

    close(11)
    close(12)
    close(13)

    if (any(mask /= 0 .and. mask /= 1)) &
        error stop "Invalid mask values"

    active_count = count(mask == 1)
    if (active_count == 0) error stop "No active cells"

    if (any((.not. ieee_is_finite(initial)) .and. mask == 1)) &
        error stop "Non-finite initial temperatures"

    if (any((.not. ieee_is_finite(albedo)) .and. mask == 1)) &
        error stop "Non-finite albedo values"

    if (any((albedo < 0d0 .or. albedo > 1d0) .and. mask == 1)) &
        error stop "Albedo outside physical range"

    ! ========================================================
    ! 2. SAME PARAMETERS AS COUPLED BASELINE
    ! ========================================================

    sigma = 5.670374419d-8
    emissivity = 0.95d0

    ca = 2.0d6
    kappa = 1.0d-5

    solar = 500.0d0
    tair = 293.15d0
    tsky = 283.15d0
    h = 10.0d0
    latent = 50.0d0

    duration = 21600.0d0

    stepsizes = [120.0d0, 60.0d0, 30.0d0, 15.0d0]

    mean_initial = sum(initial, mask=mask==1) &
                   / real(active_count,dp)

    print *
    print *, "============================================"
    print *, "MBABANE COUPLED MODEL TIME CONVERGENCE"
    print *, "============================================"
    print *, "Grid:", nx, "x", ny
    print *, "Active pixels:", active_count
    print *, "Initial mean (C):", mean_initial-273.15d0
    print *
    print *, "dt(s)   steps    final mean (C)"
    print *, "--------------------------------------------"

    ! ========================================================
    ! 3. REPEAT SIMULATION FOR FOUR TIME STEPS
    ! ========================================================

    do k = 1, nruns

        dt = stepsizes(k)
        nt = nint(duration/dt)

        stability = kappa*dt * &
                    (1d0/dx**2 + 1d0/dy**2)

        if (stability > 0.5d0) &
            error stop "Diffusion stability condition failed"

        ! Reset to identical initial satellite temperature

        temp = initial

        do n = 1, nt

            temp_new = temp

            do j = 1, ny
                do i = 1, nx

                    if (mask(i,j) /= 1) cycle

                    laplacian = 0d0

                    ! West
                    if (i > 1) then
                        if (mask(i-1,j) == 1) then
                            laplacian = laplacian + &
                                (temp(i-1,j)-temp(i,j))/dx**2
                        end if
                    end if

                    ! East
                    if (i < nx) then
                        if (mask(i+1,j) == 1) then
                            laplacian = laplacian + &
                                (temp(i+1,j)-temp(i,j))/dx**2
                        end if
                    end if

                    ! South / previous raster row
                    if (j > 1) then
                        if (mask(i,j-1) == 1) then
                            laplacian = laplacian + &
                                (temp(i,j-1)-temp(i,j))/dy**2
                        end if
                    end if

                    ! North / next raster row
                    if (j < ny) then
                        if (mask(i,j+1) == 1) then
                            laplacian = laplacian + &
                                (temp(i,j+1)-temp(i,j))/dy**2
                        end if
                    end if

                    net_flux = (1d0-albedo(i,j))*solar &
                         - emissivity*sigma * &
                           (temp(i,j)**4-tsky**4) &
                         - h*(temp(i,j)-tair) &
                         - latent

                    temp_new(i,j) = temp(i,j) + dt * &
                         (kappa*laplacian + net_flux/ca)

                end do
            end do

            temp = temp_new

            if (any((.not. ieee_is_finite(temp)) .and. mask == 1)) &
                error stop "Non-finite temperature"

        end do

        solutions(:,:,k) = temp

        mean_final = sum(temp, mask=mask==1) &
                     / real(active_count,dp)

        print '(F7.1,2X,I6,2X,F18.10)', &
            dt, nt, mean_final-273.15d0

    end do

    ! ========================================================
    ! 4. RMSE AGAINST FINEST TIME STEP (15 SECONDS)
    ! ========================================================

    print *
    print *, "RMSE against dt = 15 s"
    print *, "--------------------------------------------"

    do k = 1, nruns-1

        rmse = sqrt( &
            sum((solutions(:,:,k)-solutions(:,:,nruns))**2, &
                mask=mask==1) / real(active_count,dp))

        print '(A,F7.1,A,ES14.6,A)', &
            "dt = ",stepsizes(k), &
            " s   RMSE = ",rmse," C"

    end do

    ! ========================================================
    ! 5. SELF-CONVERGENCE AND ESTIMATED ORDER
    ! ========================================================

    do k = 1, nruns-1

        diff_norm(k) = sqrt( &
            sum((solutions(:,:,k)-solutions(:,:,k+1))**2, &
                mask=mask==1) / real(active_count,dp))

    end do

    print *
    print *, "Successive solution differences (RMSE)"
    print *, "--------------------------------------------"

    do k = 1, nruns-1

        print '(F7.1,A,F7.1,A,ES14.6)', &
            stepsizes(k)," vs ",stepsizes(k+1), &
            " s: ",diff_norm(k)

    end do

    print *
    print *, "Estimated temporal convergence orders"
    print *, "--------------------------------------------"

    do k = 1, nruns-2

        if (diff_norm(k+1) > 0d0 .and. &
            diff_norm(k) > 0d0) then

            order = log(diff_norm(k)/diff_norm(k+1)) &
                    / log(2d0)

            print '(A,I1,A,F12.6)', &
                "Order estimate ",k,": ",order

        else
            print *, "Cannot estimate order: zero difference"
        end if

    end do

    ! ========================================================
    ! 6. SAVE RESULTS TO CSV
    ! ========================================================

    filename = "results/mbabane_coupled_time_convergence.csv"

    open(20, file=trim(filename), status="replace", &
         action="write", iostat=ios)

    if (ios /= 0) error stop "Cannot create CSV file"

    write(20,'(A)') "dt_seconds,steps,final_mean_C,rmse_vs_15s_C"

    do k = 1, nruns

        mean_final = sum(solutions(:,:,k),mask=mask==1) &
                     / real(active_count,dp)

        rmse = sqrt( &
            sum((solutions(:,:,k)-solutions(:,:,nruns))**2, &
                mask=mask==1) / real(active_count,dp))

        write(20,'(F8.2,A,I8,A,F18.10,A,ES16.8)') &
            stepsizes(k),",",nint(duration/stepsizes(k)), &
            ",",mean_final-273.15d0,",",rmse

    end do

    close(20)

    print *
    print *, "Saved: ",trim(filename)
    print *, "Time convergence test completed."

end program mbabane_coupled_time_convergence
