
program manzini_sensitivity_ndvi020

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: ns = 4

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

    real(dp), parameter :: ndvi_threshold = 0.20d0
    real(dp), parameter :: delta_albedo = 0.10d0
    real(dp), parameter :: max_albedo = 0.80d0
    real(dp), parameter :: delta_latent = 30.0d0

    real(dp), parameter :: kelvin_offset = 273.15d0

    integer :: nx, ny, nt, np
    integer :: i, j, n, s, ios
    integer :: active_count, eligible_count
    integer, allocatable :: mask(:,:)

    real(dp) :: dx, dy, stability
    real(dp) :: lap, flux, alb, latent
    real(dp) :: sum_initial, sum_final
    real(dp) :: mean_initial, mean_final
    real(dp) :: mean_change, max_change

    real(dp), allocatable :: initial(:,:), ndvi(:,:)
    real(dp), allocatable :: albedo(:,:), reflectance(:,:)
    real(dp), allocatable :: temp(:,:), next_temp(:,:)
    real(dp), allocatable :: final_fields(:,:,:)
    logical, allocatable :: eligible(:,:)

    character(len=80), parameter :: output_files(ns) = [ &
        character(len=80) :: &
        "results/manzini_ndvi020_S0_K.txt", &
        "results/manzini_ndvi020_S1_K.txt", &
        "results/manzini_ndvi020_S2_K.txt", &
        "results/manzini_ndvi020_S3_K.txt" ]

    ! Read dimensions and pixel spacing
    open(10,file="data/manzini_energy_grid.txt", &
         status="old",action="read",iostat=ios)
    if (ios /= 0) error stop "Cannot open grid file"

    read(10,*,iostat=ios) nx, ny
    if (ios /= 0) error stop "Invalid grid dimensions"

    read(10,*,iostat=ios) dx, dy
    if (ios /= 0) error stop "Invalid pixel spacing"
    close(10)

    if (nx < 2 .or. ny < 2) error stop "Invalid grid size"
    if (dx <= 0d0 .or. dy <= 0d0) &
        error stop "Invalid pixel size"

    allocate(initial(nx,ny),ndvi(nx,ny))
    allocate(albedo(nx,ny),reflectance(nx,ny))
    allocate(mask(nx,ny),eligible(nx,ny))
    allocate(temp(nx,ny),next_temp(nx,ny))
    allocate(final_fields(nx,ny,ns))

    call read_real_grid( &
        "data/manzini_energy_temperature_K.txt",initial)

    call read_real_grid( &
        "data/manzini_energy_ndvi.txt",ndvi)

    call read_real_grid( &
        "data/manzini_energy_albedo.txt",albedo)

    call read_real_grid( &
        "data/manzini_energy_reflectance.txt",reflectance)

    call read_integer_grid( &
        "data/manzini_energy_mask.txt",mask)

    if (any(mask /= 0 .and. mask /= 1)) &
        error stop "Mask must contain only 0 or 1"

    active_count = count(mask == 1)

    if (active_count == 0) &
        error stop "No active pixels"

    do j = 1, ny
        do i = 1, nx
            if (mask(i,j) == 1) then
                if (.not.ieee_is_finite(initial(i,j)) .or. &
                    .not.ieee_is_finite(ndvi(i,j)) .or. &
                    .not.ieee_is_finite(albedo(i,j)) .or. &
                    .not.ieee_is_finite(reflectance(i,j))) &
                    error stop "Non-finite active input"

                if (initial(i,j) <= 0d0) &
                    error stop "Invalid Kelvin temperature"

                if (albedo(i,j) < 0d0 .or. &
                    albedo(i,j) > 1d0) &
                    error stop "Albedo outside physical range"
            end if
        end do
    end do

    eligible = (mask == 1 .and. ndvi < ndvi_threshold)
    eligible_count = count(eligible)

    nt = nint(duration/dt)

    if (abs(real(nt,dp)*dt-duration) > 1d-8) &
        error stop "Duration not divisible by dt"

    stability = kappa*dt*(1d0/dx**2+1d0/dy**2)

    if (stability > 0.5d0) &
        error stop "Diffusion stability condition failed"

    np = active_count

    sum_initial = sum(initial,mask=mask==1)
    mean_initial = sum_initial/real(np,dp)-kelvin_offset

    open(30,file="results/manzini_ndvi020_summary.csv", &
         status="replace",action="write",iostat=ios)

    if (ios /= 0) error stop "Cannot create summary CSV"

    write(30,'(A)') &
      "scenario,active_pixels,eligible_pixels," // &
      "initial_mean_C,final_mean_C,mean_change_C," // &
      "max_abs_cell_change_C"

    print *
    print *, "=========================================="
    print *, "MANZINI COUPLED MITIGATION SCENARIOS"
    print *, "=========================================="
    print *, "Grid:",nx,"x",ny
    print *, "Active pixels:",active_count
    print *, "Eligible pixels:",eligible_count
    print *, "Duration (hours):",duration/3600d0
    print *, "Time step (seconds):",dt
    print *, "Initial mean (C):",mean_initial
    print *

    ! Run all scenarios from identical initial conditions
    do s = 1, ns

        temp = initial

        do n = 1, nt

            next_temp = temp

            do j = 1, ny
                do i = 1, nx

                    if (mask(i,j) /= 1) cycle

                    ! Zero-flux across grid and mask boundaries
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

                        if (s == 2 .or. s == 4) then
                            alb = max(alb, min(max_albedo, &
                                      alb+delta_albedo))
                        end if

                        if (s == 3 .or. s == 4) then
                            latent = latent_base+delta_latent
                        end if

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
            error stop "Non-finite simulated temperature"

        final_fields(:,:,s) = temp

        sum_final = sum(temp,mask=mask==1)

        mean_final = &
            sum_final/real(np,dp)-kelvin_offset

        mean_change = mean_final-mean_initial

        max_change = maxval( &
            abs(temp-initial),mask=mask==1)

        print *, "Scenario S",s-1
        print *, " Final mean (C):",mean_final
        print *, " Mean warming (C):",mean_change
        print *, " Max absolute cell change (C):",max_change
        print *

        write(30,'(I0,A,I0,A,I0,A,ES24.15,A,' // &
                  'ES24.15,A,ES24.15,A,ES24.15)') &
            s-1,",",active_count,",",eligible_count,",", &
            mean_initial,",",mean_final,",", &
            mean_change,",",max_change

        call write_real_grid(trim(output_files(s)),temp)

    end do

    close(30)

    print *, "------------------------------------------"
    print *, "Cooling relative to baseline (C)"
    print *, "------------------------------------------"

    do s = 2, ns

        print '(A,I1,A,F15.9)', &
            "Scenario S",s-1,": ", &
            sum(final_fields(:,:,1)-final_fields(:,:,s), &
                mask=mask==1)/real(np,dp)

    end do

    print *
    print *, "Saved scenario outputs and summary CSV."
    print *, "Mitigation simulations completed."

contains

    subroutine read_real_grid(filename,arr)

        character(len=*),intent(in) :: filename
        real(dp),intent(out) :: arr(:,:)
        integer :: unit,stat,ii,jj

        open(newunit=unit,file=filename, &
             status="old",action="read",iostat=stat)
        if (stat /= 0) then
            print *, "Cannot open:",trim(filename)
            error stop
        end if

        ! Row-major input: nx values per raster row
        do jj = 1, size(arr,2)
            read(unit,*,iostat=stat) &
                (arr(ii,jj),ii=1,size(arr,1))
            if (stat /= 0) then
                print *, "Invalid data:",trim(filename)
                error stop
            end if
        end do

        close(unit)

    end subroutine read_real_grid

    subroutine read_integer_grid(filename,arr)

        character(len=*),intent(in) :: filename
        integer,intent(out) :: arr(:,:)
        integer :: unit,stat,ii,jj

        open(newunit=unit,file=filename, &
             status="old",action="read",iostat=stat)
        if (stat /= 0) then
            print *, "Cannot open:",trim(filename)
            error stop
        end if

        do jj = 1, size(arr,2)
            read(unit,*,iostat=stat) &
                (arr(ii,jj),ii=1,size(arr,1))
            if (stat /= 0) then
                print *, "Invalid mask:",trim(filename)
                error stop
            end if
        end do

        close(unit)

    end subroutine read_integer_grid

    subroutine write_real_grid(filename,arr)

        character(len=*),intent(in) :: filename
        real(dp),intent(in) :: arr(:,:)
        integer :: unit,stat,ii,jj

        open(newunit=unit,file=filename, &
             status="replace",action="write",iostat=stat)
        if (stat /= 0) then
            print *, "Cannot write:",trim(filename)
            error stop
        end if

        do jj = 1, size(arr,2)
            write(unit,'(*(ES24.15,1X))') &
                (arr(ii,jj),ii=1,size(arr,1))
        end do

        close(unit)

    end subroutine write_real_grid

end program manzini_sensitivity_ndvi020
