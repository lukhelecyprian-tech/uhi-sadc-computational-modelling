
program johannesburg_sensitivity_analysis
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: ntests = 9

    real(dp), parameter :: sigma = 5.670374419d-8
    real(dp), parameter :: emissivity = 0.95d0
    real(dp), parameter :: ca = 2.0d6
    real(dp), parameter :: kappa = 1.0d-5
    real(dp), parameter :: solar = 500.0d0
    real(dp), parameter :: tair = 293.15d0
    real(dp), parameter :: tsky = 283.15d0
    real(dp), parameter :: h = 10.0d0
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

    integer :: nx, ny, nt, i, j, n, test, ios
    integer :: active_count, eligible_count
    integer, allocatable :: mask(:,:)
    real(dp) :: dx, dy, stability
    real(dp) :: threshold, dalbedo, dlatent
    real(dp) :: lap, flux, alb, latent
    real(dp) :: final_mean, cooling_domain
    real(dp) :: cooling_eligible, max_cooling
    real(dp) :: reference_difference
    real(dp), allocatable :: initial(:,:), ndvi(:,:)
    real(dp), allocatable :: albedo(:,:), reflectance(:,:)
    real(dp), allocatable :: baseline(:,:)
    real(dp), allocatable :: original_s3(:,:)
    real(dp), allocatable :: temp(:,:), next_temp(:,:)
    logical, allocatable :: eligible(:,:)

    ! Read original grid geometry
    open(10, file="data/johannesburg_energy_grid.txt", &
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
    allocate(mask(nx,ny), eligible(nx,ny))

    call read_real_grid( &
        "data/johannesburg_energy_temperature_K.txt", initial)
    call read_real_grid( &
        "data/johannesburg_energy_ndvi.txt", ndvi)
    call read_real_grid( &
        "data/johannesburg_energy_albedo.txt", albedo)
    call read_real_grid( &
        "data/johannesburg_energy_reflectance.txt", reflectance)

    call read_integer_grid( &
        "data/johannesburg_energy_mask.txt", mask)

    ! Previously completed simulations, used for comparison.
    call read_real_grid( &
        "results/johannesburg_mitigation_S0_K.txt", baseline)
    call read_real_grid( &
        "results/johannesburg_mitigation_S3_K.txt", original_s3)

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

    stability = kappa*dt*(1d0/dx**2+1d0/dy**2)
    if (stability > 0.5d0) &
        error stop "Diffusion stability condition failed"

    open(30, &
      file="results/johannesburg_sensitivity_spatial.csv", &
      status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "Cannot create CSV"

    write(30,'(A)') &
      "test,parameter,ndvi_threshold,delta_albedo," // &
      "delta_latent_Wm2,eligible_pixels," // &
      "final_mean_C,domain_mean_cooling_C," // &
      "eligible_mean_cooling_C,max_cooling_C"

    print *
    print *, "======================================"
    print *, "JOHANNESBURG PARAMETER SENSITIVITY ANALYSIS"
    print *, "======================================"
    print *, "Grid:", nx, "x", ny
    print *, "Active pixels:", active_count
    print *, "Simulation hours:", duration/3600d0
    print *, "Number of tests:", ntests
    print *

    do test = 1, ntests

        ! One-at-a-time parameter variations
        threshold = 0.30d0
        dalbedo = 0.10d0
        dlatent = 30.0d0

        if (test <= 3) then
            threshold = thresholds(test)
        else if (test <= 6) then
            dalbedo = albedo_changes(test-3)
        else
            dlatent = latent_changes(test-6)
        end if

        eligible = (mask == 1 .and. ndvi < threshold)
        eligible_count = count(eligible)

        temp = initial

        ! Same numerical solver as the original S3 run
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

                    if (eligible(i,j)) then
                        alb = min(max_albedo, alb+dalbedo)
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
            error stop "Non-finite output"

        final_mean = sum(temp, mask=mask==1) / &
            real(active_count,dp) - 273.15d0

        cooling_domain = &
            sum(baseline-temp, mask=mask==1) / &
            real(active_count,dp)

        cooling_eligible = 0d0
        if (eligible_count > 0) then
            cooling_eligible = &
                sum(baseline-temp, mask=eligible) / &
                real(eligible_count,dp)
        end if

        max_cooling = maxval( &
            baseline-temp, mask=mask==1)

        print *, "Test:", test
        print *, " NDVI threshold:", threshold
        print *, " Albedo increase:", dalbedo
        print *, " Latent increase (W/m2):", dlatent
        print *, " Eligible pixels:", eligible_count
        print *, " Final mean (C):", final_mean
        print *, " Domain cooling (C):", cooling_domain
        print *, " Eligible cooling (C):", cooling_eligible
        print *, " Maximum cooling (C):", max_cooling
        print *

        if (test <= 3) then
            call save_row("NDVI_threshold")
        else if (test <= 6) then
            call save_row("albedo_increase")
        else
            call save_row("latent_increase")
        end if

        ! Save spatial cooling for NDVI sensitivity tests.
        if (test <= 3) then
            call save_ndvi_cooling(test)
        end if

        ! Three centre-value tests must reproduce original S3
        if (test == 2 .or. test == 5 .or. test == 8) then
            reference_difference = maxval( &
                abs(temp-original_s3), mask=mask==1)

            print *, " S3 reference maximum difference (K):", &
                reference_difference

            if (reference_difference > 1d-8) &
                error stop "Original S3 reproduction failed"
        end if

    end do

    close(30)

    print *
    print *, "Sensitivity analysis completed."
    print *, "CSV: results/johannesburg_sensitivity_spatial.csv"
    print *, "Original scenario files were not modified."

contains

    subroutine save_ndvi_cooling(test_number)
        integer, intent(in) :: test_number
        integer :: unit_number, status_code
        integer :: row, col
        character(len=120) :: filename
        character(len=3) :: threshold_label

        select case (test_number)
        case (1)
            threshold_label = "020"
        case (2)
            threshold_label = "030"
        case (3)
            threshold_label = "040"
        case default
            error stop "Invalid NDVI sensitivity test"
        end select

        filename = "results/johannesburg_ndvi" // &
                   threshold_label // "_cooling.txt"

        open(newunit=unit_number, file=trim(filename), &
             status="new", action="write", &
             iostat=status_code)

        if (status_code /= 0) &
            error stop "Cannot create spatial cooling file"

        do row = 1, ny
            write(unit_number, '(*(ES24.15,1X))') &
                (baseline(col,row)-temp(col,row), col=1,nx)
        end do

        close(unit_number)

        print *, " Spatial cooling saved: ", trim(filename)

    end subroutine save_ndvi_cooling


    subroutine save_row(parameter_name)
        character(len=*), intent(in) :: parameter_name

        write(30, &
          '(I0,",",A,",",F8.4,",",F8.4,",",' // &
          'F10.4,",",I0,",",ES24.15,",",' // &
          'ES24.15,",",ES24.15,",",ES24.15)') &
          test, trim(parameter_name), threshold, &
          dalbedo, dlatent, eligible_count, &
          final_mean, cooling_domain, &
          cooling_eligible, max_cooling
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

end program johannesburg_sensitivity_analysis
