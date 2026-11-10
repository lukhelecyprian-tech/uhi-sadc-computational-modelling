
program manzini_physical_sensitivity
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: ntests = 15
    real(dp), parameter :: solar_values(3) = &
        [400.0d0, 500.0d0, 600.0d0]
    real(dp), parameter :: ca_values(3) = &
        [1.5d6, 2.0d6, 2.5d6]
    real(dp), parameter :: emissivity_values(3) = &
        [0.90d0, 0.95d0, 0.98d0]
    real(dp), parameter :: h_values(3) = &
        [5.0d0, 10.0d0, 15.0d0]
    real(dp), parameter :: kappa_values(3) = &
        [5.0d-6, 1.0d-5, 2.0d-5]

    real(dp), parameter :: sigma = 5.670374419d-8
    real(dp), parameter :: emissivity_ref = 0.95d0
    real(dp) :: emissivity = emissivity_ref
    real(dp), parameter :: ca_ref = 2.0d6
    real(dp) :: ca = ca_ref
    real(dp), parameter :: kappa_ref = 1.0d-5
    real(dp) :: kappa = kappa_ref
    real(dp), parameter :: solar_ref = 500.0d0
    real(dp) :: solar = solar_ref
    real(dp), parameter :: tair = 293.15d0
    real(dp), parameter :: tsky = 283.15d0
    real(dp), parameter :: h_ref = 10.0d0
    real(dp) :: h = h_ref
    real(dp), parameter :: latent_base = 50.0d0
    real(dp), parameter :: duration = 21600.0d0
    real(dp), parameter :: dt = 60.0d0
    real(dp), parameter :: max_albedo = 0.80d0

    real(dp), parameter :: thresholds(3) = &
        [0.20d0, 0.30d0, 0.40d0]
    real(dp), parameter :: albedo_changes(3) = &
        [0.05d0, 0.10d0, 0.15d0]
    real(dp), parameter :: latent_changes(3) = &
        [15.0d0, 30.0d0, 45.0d0]

    integer :: nx, ny, nt, i, j, n, test, scenario, ios
    integer :: active_count, eligible_count
    integer, allocatable :: mask(:,:)
    real(dp) :: dx, dy, stability
    real(dp) :: threshold, dalbedo, dlatent
    real(dp) :: lap, flux, alb, latent
    real(dp) :: final_mean, cooling_domain
    real(dp) :: s0_mean, parameter_value
    character(len=32) :: parameter_name
    real(dp) :: cooling_eligible, max_cooling
    real(dp) :: reference_difference
    real(dp), allocatable :: initial(:,:), ndvi(:,:)
    real(dp), allocatable :: albedo(:,:), reflectance(:,:)
    real(dp), allocatable :: baseline(:,:)
    real(dp), allocatable :: original_s3(:,:)
    real(dp), allocatable :: temp(:,:), next_temp(:,:)
    real(dp), allocatable :: simulated_s0(:,:), simulated_s3(:,:)
    logical, allocatable :: eligible(:,:)

    ! Read original grid geometry
    open(10, file="data/manzini_energy_grid.txt", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Cannot open grid file"

    read(10,*,iostat=ios) nx, ny
    if (ios /= 0) error stop "Invalid grid dimensions"

    read(10,*,iostat=ios) dx, dy
    if (ios /= 0) error stop "Invalid pixel spacing"
    close(10)

    if (nx < 2 .or. ny < 2) error stop "Invalid grid"
    if (dx <= 0d0 .or. dy <= 0d0) &
        error stop "Invalid pixel size"

    allocate(initial(nx,ny), ndvi(nx,ny))
    allocate(albedo(nx,ny), reflectance(nx,ny))
    allocate(baseline(nx,ny), original_s3(nx,ny))
    allocate(temp(nx,ny), next_temp(nx,ny))
    allocate(simulated_s0(nx,ny), simulated_s3(nx,ny))
    allocate(mask(nx,ny), eligible(nx,ny))

    call read_real_grid( &
        "data/manzini_energy_temperature_K.txt", initial)
    call read_real_grid( &
        "data/manzini_energy_ndvi.txt", ndvi)
    call read_real_grid( &
        "data/manzini_energy_albedo.txt", albedo)
    call read_real_grid( &
        "data/manzini_energy_reflectance.txt", reflectance)

    call read_integer_grid( &
        "data/manzini_energy_mask.txt", mask)

    ! Previously completed simulations, used for comparison.
    call read_real_grid( &
        "results/manzini_mitigation_S0_K.txt", baseline)
    call read_real_grid( &
        "results/manzini_mitigation_S3_K.txt", original_s3)

    if (any(mask /= 0 .and. mask /= 1)) &
        error stop "Invalid mask"

    active_count = count(mask == 1)
    if (active_count == 0) error stop "Empty mask"

    do j = 1, ny
        do i = 1, nx
            if (mask(i,j) /= 1) cycle

            if (.not.ieee_is_finite(initial(i,j)) .or. &
                .not.ieee_is_finite(ndvi(i,j)) .or. &
                .not.ieee_is_finite(albedo(i,j)) .or. &
                .not.ieee_is_finite(reflectance(i,j)) .or. &
                .not.ieee_is_finite(baseline(i,j)) .or. &
                .not.ieee_is_finite(original_s3(i,j))) &
                error stop "Non-finite input"

            if (initial(i,j) <= 0d0) &
                error stop "Invalid initial temperature"

            if (albedo(i,j) < 0d0 .or. &
                albedo(i,j) > 1d0) &
                error stop "Invalid albedo"
        end do
    end do

    nt = nint(duration/dt)
    if (abs(real(nt,dp)*dt-duration) > 1d-8) &
        error stop "Invalid time interval"


    open(30, &
      file="results/manzini_physical_sensitivity.csv", &
      status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "Cannot create CSV"

    write(30,'(A)') &
      "test,parameter,parameter_value," // &
      "solar_Wm2,heat_capacity_Jm2K,emissivity," // &
      "heat_transfer_Wm2K,diffusivity_m2s," // &
      "eligible_pixels,S0_mean_C,S3_mean_C," // &
      "domain_mean_cooling_C," // &
      "eligible_mean_cooling_C,max_cooling_C"

    print *
    print *, "======================================"
    print *, "MANZINI PHYSICAL SENSITIVITY ANALYSIS"
    print *, "======================================"
    print *, "Grid:", nx, "x", ny
    print *, "Active pixels:", active_count
    print *, "Simulation hours:", duration/3600d0
    print *, "Number of tests:", ntests
    print *

    do test = 1, ntests

        ! Reset all physical parameters to reference values
        solar = solar_ref
        ca = ca_ref
        emissivity = emissivity_ref
        h = h_ref
        kappa = kappa_ref

        ! Keep mitigation settings fixed
        threshold = 0.30d0
        dalbedo = 0.10d0
        dlatent = 30.0d0

        ! One-at-a-time physical parameter variations
        if (test <= 3) then
            solar = solar_values(test)
            parameter_name = "solar"
            parameter_value = solar
        else if (test <= 6) then
            ca = ca_values(test-3)
            parameter_name = "heat_capacity"
            parameter_value = ca
        else if (test <= 9) then
            emissivity = emissivity_values(test-6)
            parameter_name = "emissivity"
            parameter_value = emissivity
        else if (test <= 12) then
            h = h_values(test-9)
            parameter_name = "heat_transfer"
            parameter_value = h
        else
            kappa = kappa_values(test-12)
            parameter_name = "diffusivity"
            parameter_value = kappa
        end if

        ! Check stability using this experiment's diffusivity
        stability = kappa*dt*(1d0/dx**2+1d0/dy**2)
        if (stability > 0.5d0) &
            error stop "Diffusion stability condition failed"

        eligible = (mask == 1 .and. ndvi < threshold)
        eligible_count = count(eligible)

        ! Run independent baseline and mitigation simulations
        do scenario = 0, 1

            ! Both scenarios start from the same Landsat field
            temp = initial

            do n = 1, nt
            next_temp = temp

            do j = 1, ny
                do i = 1, nx
                    if (mask(i,j) /= 1) cycle

                    ! Zero-flux boundaries
                    lap = 0d0

                    if (i > 1) then
                        if (mask(i-1,j) == 1) &
                            lap = lap + &
                            (temp(i-1,j)-temp(i,j))/dx**2
                    end if

                    if (i < nx) then
                        if (mask(i+1,j) == 1) &
                            lap = lap + &
                            (temp(i+1,j)-temp(i,j))/dx**2
                    end if

                    if (j > 1) then
                        if (mask(i,j-1) == 1) &
                            lap = lap + &
                            (temp(i,j-1)-temp(i,j))/dy**2
                    end if

                    if (j < ny) then
                        if (mask(i,j+1) == 1) &
                            lap = lap + &
                            (temp(i,j+1)-temp(i,j))/dy**2
                    end if

                    alb = albedo(i,j)
                    latent = latent_base

                    ! Mitigation applies only to S3
                    if (scenario == 1 .and. eligible(i,j)) then
                        alb = max(alb, min(max_albedo, alb+dalbedo))
                        latent = latent_base+dlatent
                    end if

                    flux = (1d0-alb)*solar &
                        - emissivity*sigma* &
                          (temp(i,j)**4-tsky**4) &
                        - h*(temp(i,j)-tair) &
                        - latent

                    next_temp(i,j) = temp(i,j) + dt * &
                        (kappa*lap + flux/ca)

                end do
            end do

            temp = next_temp
            end do

            if (.not.all(ieee_is_finite(temp))) &
                error stop "Non-finite scenario output"

            if (scenario == 0) then
                simulated_s0 = temp
            else
                simulated_s3 = temp
            end if

        end do

        if (.not.all(ieee_is_finite(temp))) &
            error stop "Non-finite output"

        s0_mean = sum(simulated_s0, mask=mask==1) / &
            real(active_count,dp) - 273.15d0

        final_mean = sum(simulated_s3, mask=mask==1) / &
            real(active_count,dp) - 273.15d0

        cooling_domain = &
            sum(simulated_s0-simulated_s3, mask=mask==1) / &
            real(active_count,dp)

        cooling_eligible = 0d0
        if (eligible_count > 0) then
            cooling_eligible = &
                sum(simulated_s0-simulated_s3, mask=eligible) / &
                real(eligible_count,dp)
        end if

        max_cooling = maxval( &
            simulated_s0-simulated_s3, mask=mask==1)

        print *, "Test:", test
        print *, " Parameter:", trim(parameter_name)
        print *, " Parameter value:", parameter_value
        print *, " Eligible pixels:", eligible_count
        print *, " S0 mean (C):", s0_mean
        print *, " S3 mean (C):", final_mean
        print *, " Domain cooling (C):", cooling_domain
        print *, " Eligible cooling (C):", cooling_eligible
        print *, " Maximum cooling (C):", max_cooling

        ! Verify both scenarios at the reference settings
        if (mod(test-2,3) == 0) then

            reference_difference = maxval( &
                abs(simulated_s0-baseline), mask=mask==1)

            print *, " S0 reference maximum difference (K):", &
                reference_difference

            if (reference_difference > 1d-8) &
                error stop "Original S0 reproduction failed"

            reference_difference = maxval( &
                abs(simulated_s3-original_s3), mask=mask==1)

            print *, " S3 reference maximum difference (K):", &
                reference_difference

            if (reference_difference > 1d-8) &
                error stop "Original S3 reproduction failed"

        end if

        ! Write the completed experiment only after verification
        call save_row()

        print *

    end do

    close(30)

    print *
    print *, "Sensitivity analysis completed."
    print *, "CSV: results/manzini_physical_sensitivity.csv"
    print *, "Original scenario files were not modified."

contains

    subroutine save_row()
        write(30, &
          '(I0,",",A,12(",",ES24.15))') &
          test, trim(parameter_name), parameter_value, &
          solar, ca, emissivity, h, kappa, &
          real(eligible_count,dp), s0_mean, final_mean, &
          cooling_domain, cooling_eligible, max_cooling
    end subroutine save_row

    subroutine read_real_grid(filename, arr)
        character(len=*), intent(in) :: filename
        real(dp), intent(out) :: arr(:,:)
        integer :: unit, stat, ii, jj

        open(newunit=unit, file=filename, &
             status="old", action="read", iostat=stat)
        if (stat /= 0) then
            print *, "Cannot open:", trim(filename)
            error stop
        end if

        do jj = 1, size(arr,2)
            read(unit,*,iostat=stat) &
                (arr(ii,jj), ii=1,size(arr,1))
            if (stat /= 0) then
                print *, "Invalid data:", trim(filename)
                error stop
            end if
        end do

        close(unit)
    end subroutine read_real_grid

    subroutine read_integer_grid(filename, arr)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: arr(:,:)
        integer :: unit, stat, ii, jj

        open(newunit=unit, file=filename, &
             status="old", action="read", iostat=stat)
        if (stat /= 0) then
            print *, "Cannot open:", trim(filename)
            error stop
        end if

        do jj = 1, size(arr,2)
            read(unit,*,iostat=stat) &
                (arr(ii,jj), ii=1,size(arr,1))
            if (stat /= 0) then
                print *, "Invalid mask:", trim(filename)
                error stop
            end if
        end do

        close(unit)
    end subroutine read_integer_grid

end program manzini_physical_sensitivity
