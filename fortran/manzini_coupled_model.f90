
program manzini_coupled_model

    use, intrinsic :: ieee_arithmetic
    implicit none

    integer :: nx, ny, i, j, n, nt, ios
    integer :: active_count
    integer, allocatable :: mask(:,:)

    real(8), allocatable :: temp(:,:), temp_new(:,:)
    real(8), allocatable :: initial(:,:), albedo(:,:)
    real(8), allocatable :: ndvi(:,:), reflectance(:,:)

    real(8) :: dx, dy, dt, duration
    real(8) :: sigma, emissivity, ca, kappa, ka
    real(8) :: solar, tair, tsky, h, latent
    real(8) :: absorbed, longwave, sensible
    real(8) :: net_flux, laplacian, change
    real(8) :: initial_mean, final_mean
    real(8) :: max_change, stability

    ! ========================================================
    ! 1. READ GRID INFORMATION
    ! ========================================================

    open(10, file="data/manzini_energy_grid.txt", &
         status="old", action="read", iostat=ios)

    if (ios /= 0) error stop "Cannot open grid file"

    read(10, *, iostat=ios) nx, ny
    if (ios /= 0) error stop "Invalid grid dimensions"

    read(10, *, iostat=ios) dx, dy
    if (ios /= 0) error stop "Invalid grid spacing"

    close(10)

    if (nx < 2 .or. ny < 2) error stop "Invalid grid size"
    if (dx <= 0d0 .or. dy <= 0d0) error stop "Invalid spacing"

    allocate(temp(nx,ny), temp_new(nx,ny))
    allocate(initial(nx,ny), albedo(nx,ny))
    allocate(ndvi(nx,ny), reflectance(nx,ny))
    allocate(mask(nx,ny))

    ! ========================================================
    ! 2. READ SATELLITE INPUTS
    ! ========================================================

    open(11, file="data/manzini_energy_temperature_K.txt", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Cannot open temperature"

    open(12, file="data/manzini_energy_ndvi.txt", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Cannot open NDVI"

    open(13, file="data/manzini_energy_albedo.txt", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Cannot open albedo"

    open(14, file="data/manzini_energy_reflectance.txt", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Cannot open reflectance"

    open(15, file="data/manzini_energy_mask.txt", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Cannot open mask"

    do j = 1, ny

        read(11, *, iostat=ios) (temp(i,j), i=1,nx)
        if (ios /= 0) error stop "Temperature read error"

        read(12, *, iostat=ios) (ndvi(i,j), i=1,nx)
        if (ios /= 0) error stop "NDVI read error"

        read(13, *, iostat=ios) (albedo(i,j), i=1,nx)
        if (ios /= 0) error stop "Albedo read error"

        read(14, *, iostat=ios) (reflectance(i,j), i=1,nx)
        if (ios /= 0) error stop "Reflectance read error"

        read(15, *, iostat=ios) (mask(i,j), i=1,nx)
        if (ios /= 0) error stop "Mask read error"

    end do

    close(11)
    close(12)
    close(13)
    close(14)
    close(15)

    if (any(mask /= 0 .and. mask /= 1)) then
        error stop "Mask must contain only 0 and 1"
    end if

    active_count = count(mask == 1)
    if (active_count == 0) error stop "No active cells"

    if (any((.not. ieee_is_finite(temp)) .and. mask == 1)) &
        error stop "Non-finite initial temperature"

    if (any((.not. ieee_is_finite(albedo)) .and. mask == 1)) &
        error stop "Non-finite albedo"

    if (any((albedo < 0d0 .or. albedo > 1d0) .and. mask == 1)) &
        error stop "Albedo outside 0 to 1"

    initial = temp

    ! ========================================================
    ! 3. PHYSICAL PARAMETERS
    ! Illustrative values; not yet calibrated
    ! ========================================================

    sigma = 5.670374419d-8
    emissivity = 0.95d0

    ca = 2.0d6
    kappa = 1.0d-5
    ka = ca * kappa

    solar = 500.0d0
    tair = 293.15d0
    tsky = 283.15d0
    h = 10.0d0
    latent = 50.0d0

    dt = 60.0d0
    duration = 6.0d0 * 3600.0d0
    nt = nint(duration / dt)

    stability = kappa * dt * &
        (1.0d0 / dx**2 + 1.0d0 / dy**2)

    if (stability > 0.5d0) then
        error stop "Explicit diffusion stability violated"
    end if

    ! ========================================================
    ! 4. INITIAL STATISTICS
    ! ========================================================

    initial_mean = sum(temp, mask=mask==1) / dble(active_count)

    print *
    print *, "==========================================="
    print *, "MANZINI COUPLED SURFACE ENERGY MODEL"
    print *, "==========================================="
    print *, "Grid:", nx, "x", ny
    print *, "Active pixels:", active_count
    print *, "Duration (hours):", duration / 3600d0
    print *, "Diffusivity (m2/s):", kappa
    print *, "Effective lateral conductance (W/K):", ka
    print *, "Diffusion stability sum:", stability

    ! ========================================================
    ! 5. COUPLED EXPLICIT EULER SOLVER
    ! ========================================================

    do n = 1, nt

        temp_new = temp

        do j = 1, ny
            do i = 1, nx

                if (mask(i,j) /= 1) cycle

                ! --------------------------------------------
                ! MODEL A: LATERAL HEAT DIFFUSION
                ! Zero flux across masked boundaries
                ! --------------------------------------------

                laplacian = 0.0d0

                if (i > 1) then
                    if (mask(i-1,j) == 1) then
                        laplacian = laplacian + &
                            (temp(i-1,j)-temp(i,j))/dx**2
                    end if
                end if

                if (i < nx) then
                    if (mask(i+1,j) == 1) then
                        laplacian = laplacian + &
                            (temp(i+1,j)-temp(i,j))/dx**2
                    end if
                end if

                if (j > 1) then
                    if (mask(i,j-1) == 1) then
                        laplacian = laplacian + &
                            (temp(i,j-1)-temp(i,j))/dy**2
                    end if
                end if

                if (j < ny) then
                    if (mask(i,j+1) == 1) then
                        laplacian = laplacian + &
                            (temp(i,j+1)-temp(i,j))/dy**2
                    end if
                end if

                ! --------------------------------------------
                ! MODEL B: SURFACE ENERGY BALANCE
                ! --------------------------------------------

                absorbed = (1d0-albedo(i,j)) * solar

                longwave = emissivity * sigma * &
                    (temp(i,j)**4 - tsky**4)

                sensible = h * (temp(i,j)-tair)

                net_flux = absorbed - longwave &
                         - sensible - latent

                ! --------------------------------------------
                ! COUPLED TEMPERATURE UPDATE
                ! --------------------------------------------

                change = dt * ( &
                    kappa * laplacian + net_flux / ca )

                temp_new(i,j) = temp(i,j) + change

            end do
        end do

        temp = temp_new

        if (any((.not. ieee_is_finite(temp)) .and. mask == 1)) &
            error stop "Non-finite temperature during simulation"

    end do

    ! ========================================================
    ! 6. FINAL STATISTICS
    ! ========================================================

    final_mean = sum(temp, mask=mask==1) / dble(active_count)

    max_change = maxval(abs(temp-initial), mask=mask==1)

    print *
    print *, "Initial mean temperature (C):", &
        initial_mean - 273.15d0

    print *, "Final mean temperature (C):", &
        final_mean - 273.15d0

    print *, "Mean temperature change (C):", &
        final_mean - initial_mean

    print *, "Maximum absolute cell change (C):", max_change

    ! ========================================================
    ! 7. EXPORT RESULTS
    ! ========================================================

    open(20, file="results/manzini_coupled_final_K.txt", &
         status="replace", action="write", iostat=ios)

    if (ios /= 0) error stop "Cannot open output file"

    do j = 1, ny
        write(20, '(*(ES18.10,1X))', iostat=ios) &
            (temp(i,j), i=1,nx)

        if (ios /= 0) error stop "Output write failed"
    end do

    close(20)

    print *
    print *, "Saved: results/manzini_coupled_final_K.txt"
    print *, "Coupled simulation completed successfully."

end program manzini_coupled_model
