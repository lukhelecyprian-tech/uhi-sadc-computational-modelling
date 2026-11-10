
program mbabane_energy_balance

    use, intrinsic :: ieee_arithmetic
    implicit none

    ! ========================================================
    ! MSc COMPUTATIONAL PHYSICS - UHI THESIS
    !
    ! MBABANE SPATIAL SURFACE ENERGY BALANCE MODEL
    !
    ! Numerical method: Explicit Euler
    ! Simulation duration: 6 hours
    ! Spatial resolution: 30 metres
    !
    ! Governing equation for this baseline:
    !
    ! C_A dT/dt =
    !     (1 - albedo) * S
    !   - epsilon * sigma * (T**4 - T_sky**4)
    !   - h * (T - T_air)
    !   - LE
    !
    ! No lateral diffusion is included in this first test.
    !
    ! Meteorological and thermal parameters are illustrative,
    ! not calibrated to Mbabane observations.
    ! ========================================================

    integer :: nx, ny
    integer :: i, j, n, nt
    integer :: active_count
    integer :: ios

    real(8) :: dx, dy, dt
    real(8) :: sigma, epsilon
    real(8) :: heat_capacity, h
    real(8) :: solar, air_temp, sky_temp
    real(8) :: latent_flux

    real(8) :: absorbed
    real(8) :: longwave
    real(8) :: sensible
    real(8) :: net_flux

    real(8) :: initial_mean, final_mean
    real(8) :: initial_min, initial_max
    real(8) :: final_min, final_max
    real(8) :: max_change

    real(8), allocatable :: temp(:,:)
    real(8), allocatable :: temp_initial(:,:)
    real(8), allocatable :: albedo(:,:)
    real(8), allocatable :: ndvi(:,:)
    real(8), allocatable :: reflectance(:,:)

    integer, allocatable :: mask(:,:)

    ! ========================================================
    ! 1. READ GRID INFORMATION
    ! ========================================================

    open(10, file="data/mbabane_energy_grid.txt", &
         status="old", action="read", iostat=ios)

    if (ios /= 0) error stop "Cannot open grid file"

    read(10, *, iostat=ios) nx, ny
    if (ios /= 0) error stop "Cannot read grid dimensions"

    read(10, *, iostat=ios) dx, dy
    if (ios /= 0) error stop "Cannot read grid spacing"

    close(10)

    if (nx <= 0 .or. ny <= 0) then
        error stop "Invalid grid dimensions"
    end if

    if (dx <= 0.0d0 .or. dy <= 0.0d0) then
        error stop "Invalid grid spacing"
    end if

    allocate(temp(nx,ny))
    allocate(temp_initial(nx,ny))
    allocate(albedo(nx,ny))
    allocate(ndvi(nx,ny))
    allocate(reflectance(nx,ny))
    allocate(mask(nx,ny))

    ! ========================================================
    ! 2. READ SATELLITE DATA
    !
    ! Python exports ny rows and nx columns.
    ! Fortran stores arrays as (column, row).
    ! ========================================================

    open(11, file="data/mbabane_energy_temperature_K.txt", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Cannot open temperature file"

    open(12, file="data/mbabane_energy_ndvi.txt", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Cannot open NDVI file"

    open(13, file="data/mbabane_energy_albedo.txt", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Cannot open albedo file"

    open(14, file="data/mbabane_energy_reflectance.txt", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Cannot open reflectance file"

    open(15, file="data/mbabane_energy_mask.txt", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Cannot open mask file"

    do j = 1, ny

        read(11, *, iostat=ios) (temp(i,j), i=1,nx)
        if (ios /= 0) error stop "Temperature read failed"

        read(12, *, iostat=ios) (ndvi(i,j), i=1,nx)
        if (ios /= 0) error stop "NDVI read failed"

        read(13, *, iostat=ios) (albedo(i,j), i=1,nx)
        if (ios /= 0) error stop "Albedo read failed"

        read(14, *, iostat=ios) (reflectance(i,j), i=1,nx)
        if (ios /= 0) error stop "Reflectance read failed"

        read(15, *, iostat=ios) (mask(i,j), i=1,nx)
        if (ios /= 0) error stop "Mask read failed"

    end do

    close(11)
    close(12)
    close(13)
    close(14)
    close(15)

    ! ========================================================
    ! 3. VALIDATE INPUTS
    ! ========================================================

    active_count = count(mask == 1)

    if (active_count == 0) then
        error stop "No active pixels"
    end if

    if (any(mask /= 0 .and. mask /= 1)) then
        error stop "Mask must contain only 0 and 1"
    end if

    if (any((.not. ieee_is_finite(temp)) .and. &
            (mask == 1))) then
        error stop "Invalid initial temperature"
    end if

    if (any((.not. ieee_is_finite(albedo)) .and. &
            (mask == 1))) then
        error stop "Invalid albedo"
    end if

    if (any((.not. ieee_is_finite(ndvi)) .and. &
            (mask == 1))) then
        error stop "Invalid NDVI"
    end if

    if (any((.not. ieee_is_finite(reflectance)) .and. &
            (mask == 1))) then
        error stop "Invalid reflectance"
    end if

    if (any((albedo < 0.0d0 .or. albedo > 1.0d0) &
            .and. (mask == 1))) then
        error stop "Albedo outside physical range"
    end if

    if (any((temp <= 0.0d0) .and. (mask == 1))) then
        error stop "Temperature must be in Kelvin"
    end if

    temp_initial = temp

    ! ========================================================
    ! 4. DEFINE ILLUSTRATIVE PHYSICAL PARAMETERS
    ! ========================================================

    sigma = 5.670374419d-8  ! W m^-2 K^-4
    epsilon = 0.95d0         ! Dimensionless emissivity

    heat_capacity = 2.0d6    ! J m^-2 K^-1
    h = 10.0d0               ! W m^-2 K^-1

    solar = 500.0d0          ! W m^-2
    air_temp = 293.15d0     ! K
    sky_temp = 283.15d0     ! K

    latent_flux = 50.0d0    ! W m^-2

    dt = 60.0d0             ! Seconds
    nt = 360                ! Six hours

    ! ========================================================
    ! 5. INITIAL TEMPERATURE STATISTICS
    ! ========================================================

    initial_mean = &
        sum(temp, mask=mask==1) / dble(active_count)

    initial_min = minval(temp, mask=mask==1)
    initial_max = maxval(temp, mask=mask==1)

    ! ========================================================
    ! 6. EXPLICIT EULER ENERGY-BALANCE SIMULATION
    ! ========================================================

    do n = 1, nt

        do j = 1, ny

            do i = 1, nx

                if (mask(i,j) /= 1) cycle

                absorbed = &
                    (1.0d0 - albedo(i,j)) * solar

                longwave = epsilon * sigma * &
                    (temp(i,j)**4 - sky_temp**4)

                sensible = &
                    h * (temp(i,j) - air_temp)

                net_flux = absorbed &
                         - longwave &
                         - sensible &
                         - latent_flux

                temp(i,j) = temp(i,j) + &
                    dt * net_flux / heat_capacity

            end do

        end do

    end do

    ! ========================================================
    ! 7. FINAL TEMPERATURE STATISTICS
    ! ========================================================

    if (any((.not. ieee_is_finite(temp)) .and. &
            (mask == 1))) then
        error stop "Non-finite final temperature"
    end if

    final_mean = &
        sum(temp, mask=mask==1) / dble(active_count)

    final_min = minval(temp, mask=mask==1)
    final_max = maxval(temp, mask=mask==1)

    max_change = maxval( &
        abs(temp - temp_initial), mask=mask==1)

    ! ========================================================
    ! 8. PRINT RESULTS
    ! ========================================================

    print *
    print *, "============================================"
    print *, "MBABANE SPATIAL SURFACE ENERGY BALANCE"
    print *, "============================================"

    print *, "Grid:", nx, "x", ny
    print *, "Pixel size (m):", dx, dy
    print *, "Active pixels:", active_count

    print *
    print *, "Simulation duration (hours):", &
             dble(nt) * dt / 3600.0d0

    print *, "Time steps:", nt
    print *, "Time step (seconds):", dt

    print *
    print *, "Initial mean temperature (C):", &
             initial_mean - 273.15d0

    print *, "Final mean temperature (C):", &
             final_mean - 273.15d0

    print *, "Mean temperature change (C):", &
             final_mean - initial_mean

    print *
    print *, "Initial minimum (C):", &
             initial_min - 273.15d0

    print *, "Final minimum (C):", &
             final_min - 273.15d0

    print *, "Initial maximum (C):", &
             initial_max - 273.15d0

    print *, "Final maximum (C):", &
             final_max - 273.15d0

    print *
    print *, "Maximum absolute cell change (C):", &
             max_change

    ! ========================================================
    ! 9. EXPORT FINAL TEMPERATURE GRID
    !
    ! Active cells: simulated temperature in Kelvin.
    ! Inactive cells: zero placeholder, not a measurement.
    ! ========================================================

    open(20, file="results/mbabane_energy_baseline_K.txt", &
         status="replace", action="write", iostat=ios)

    if (ios /= 0) error stop "Cannot open output file"

    do j = 1, ny

        write(20, '(*(ES18.10,1X))', iostat=ios) &
            (temp(i,j), i=1,nx)

        if (ios /= 0) error stop "Output write failed"

    end do

    close(20)

    print *
    print *, "Saved: results/mbabane_energy_baseline_K.txt"
    print *, "Baseline simulation completed successfully."

end program mbabane_energy_balance

