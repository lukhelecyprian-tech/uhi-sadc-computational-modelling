
program mbabane_coupled_conservation

    use, intrinsic :: ieee_arithmetic
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    integer :: nx, ny, i, j, n, nt, ios
    integer :: active_count
    integer, allocatable :: mask(:,:)

    real(dp), allocatable :: temp(:,:), temp_new(:,:)
    real(dp), allocatable :: initial(:,:), albedo(:,:)

    real(dp) :: dx, dy, dt, duration, area
    real(dp) :: sigma, emissivity, ca, kappa
    real(dp) :: solar, tair, tsky, h, latent
    real(dp) :: laplacian, net_flux, diff_term
    real(dp) :: total_flux, total_diff
    real(dp) :: energy_initial, energy_final
    real(dp) :: energy_input, energy_change
    real(dp) :: residual, relative_residual
    real(dp) :: diffusion_energy, stability
    real(dp) :: initial_mean, final_mean

    ! =====================================================
    ! 1. READ GRID
    ! =====================================================

    open(10, file="data/mbabane_energy_grid.txt", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Cannot open grid file"

    read(10, *, iostat=ios) nx, ny
    if (ios /= 0) error stop "Cannot read dimensions"

    read(10, *, iostat=ios) dx, dy
    if (ios /= 0) error stop "Cannot read spacing"

    close(10)

    if (nx < 2 .or. ny < 2) error stop "Invalid grid"
    if (dx <= 0d0 .or. dy <= 0d0) error stop "Invalid spacing"

    area = dx * dy

    allocate(temp(nx,ny), temp_new(nx,ny))
    allocate(initial(nx,ny), albedo(nx,ny))
    allocate(mask(nx,ny))

    ! =====================================================
    ! 2. READ INPUTS
    ! =====================================================

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
        read(11, *, iostat=ios) (temp(i,j), i=1,nx)
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
        error stop "Mask must contain only 0 and 1"

    active_count = count(mask == 1)
    if (active_count == 0) error stop "No active pixels"

    if (any((.not. ieee_is_finite(temp)) .and. mask == 1)) &
        error stop "Non-finite initial temperature"

    if (any((.not. ieee_is_finite(albedo)) .and. mask == 1)) &
        error stop "Non-finite albedo"

    if (any((albedo < 0d0 .or. albedo > 1d0) .and. mask == 1)) &
        error stop "Invalid albedo"

    initial = temp

    ! =====================================================
    ! 3. PARAMETERS
    ! =====================================================

    sigma = 5.670374419d-8
    emissivity = 0.95d0

    ca = 2.0d6
    kappa = 1.0d-5

    solar = 500.0d0
    tair = 293.15d0
    tsky = 283.15d0
    h = 10.0d0
    latent = 50.0d0

    dt = 60.0d0
    duration = 21600.0d0
    nt = nint(duration / dt)

    stability = kappa * dt * &
        (1d0/dx**2 + 1d0/dy**2)

    if (stability > 0.5d0) &
        error stop "Diffusion stability condition failed"

    ! =====================================================
    ! 4. INITIAL ENERGY
    ! =====================================================

    energy_initial = ca * area * &
        sum(temp, mask=mask==1)

    energy_input = 0d0
    diffusion_energy = 0d0

    initial_mean = sum(temp, mask=mask==1) &
                   / real(active_count,dp)

    ! =====================================================
    ! 5. SOLVE AND ACCUMULATE ENERGY CONTRIBUTIONS
    ! =====================================================

    do n = 1, nt

        temp_new = temp

        total_flux = 0d0
        total_diff = 0d0

        do j = 1, ny
            do i = 1, nx

                if (mask(i,j) /= 1) cycle

                laplacian = 0d0

                if (i > 1) then
                    if (mask(i-1,j) == 1) &
                        laplacian = laplacian + &
                        (temp(i-1,j)-temp(i,j))/dx**2
                end if

                if (i < nx) then
                    if (mask(i+1,j) == 1) &
                        laplacian = laplacian + &
                        (temp(i+1,j)-temp(i,j))/dx**2
                end if

                if (j > 1) then
                    if (mask(i,j-1) == 1) &
                        laplacian = laplacian + &
                        (temp(i,j-1)-temp(i,j))/dy**2
                end if

                if (j < ny) then
                    if (mask(i,j+1) == 1) &
                        laplacian = laplacian + &
                        (temp(i,j+1)-temp(i,j))/dy**2
                end if

                net_flux = (1d0-albedo(i,j))*solar &
                    - emissivity*sigma * &
                      (temp(i,j)**4-tsky**4) &
                    - h*(temp(i,j)-tair) &
                    - latent

                diff_term = kappa * laplacian

                temp_new(i,j) = temp(i,j) + dt * &
                    (diff_term + net_flux/ca)

                total_flux = total_flux + net_flux
                total_diff = total_diff + diff_term

            end do
        end do

        energy_input = energy_input + &
            total_flux * area * dt

        diffusion_energy = diffusion_energy + &
            total_diff * ca * area * dt

        temp = temp_new

        if (any((.not. ieee_is_finite(temp)) .and. mask == 1)) &
            error stop "Non-finite temperature"

    end do

    ! =====================================================
    ! 6. ENERGY BUDGET
    ! =====================================================

    energy_final = ca * area * &
        sum(temp, mask=mask==1)

    energy_change = energy_final - energy_initial

    residual = energy_change - energy_input

    relative_residual = abs(residual) / &
        max(abs(energy_change), abs(energy_input), 1d0)

    final_mean = sum(temp, mask=mask==1) &
                 / real(active_count,dp)

    print *
    print *, "============================================"
    print *, "MBABANE COUPLED MODEL CONSERVATION TEST"
    print *, "============================================"
    print *, "Active pixels:", active_count
    print *, "Initial mean (C):", initial_mean-273.15d0
    print *, "Final mean (C):", final_mean-273.15d0
    print *
    print *, "Initial energy (J):", energy_initial
    print *, "Final energy (J):", energy_final
    print *, "Stored energy change (J):", energy_change
    print *, "Surface energy input (J):", energy_input
    print *, "Diffusion contribution (J):", diffusion_energy
    print *, "Energy balance residual (J):", residual
    print *, "Relative residual:", relative_residual

    if (relative_residual < 1d-8) then
        print *
        print *, "PASS: Discrete energy budget consistent."
    else
        print *
        print *, "CHECK: Energy budget residual needs review."
    end if

    print *
    print *, "Conservation test completed."

end program mbabane_coupled_conservation
