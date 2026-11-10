
program mbabane_heat2d
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: alpha = 1.0d-5
    real(dp), parameter :: dt = 1000.0d0
    integer, parameter :: nt = 200

    integer :: nx, ny, i, j, step, unit
    integer :: nactive
    integer, allocatable :: mask(:,:)

    real(dp) :: dx, dy, rx, ry
    real(dp) :: flux_x, flux_y
    real(dp) :: mean_initial, mean_final
    real(dp) :: max_change

    real(dp), allocatable :: T(:,:), Tnew(:,:)
    real(dp), allocatable :: Tinitial(:,:)

    ! Read raster dimensions and grid spacing
    open(newunit=unit, file="data/mbabane_grid.txt", &
         status="old", action="read")
    read(unit,*) nx, ny
    read(unit,*) dx, dy
    close(unit)

    allocate(T(nx,ny), Tnew(nx,ny))
    allocate(Tinitial(nx,ny), mask(nx,ny))

    ! Read temperatures in raster row order
    open(newunit=unit, &
         file="data/mbabane_temperature.txt", &
         status="old", action="read")

    do j = 1, ny
        read(unit,*) (T(i,j), i=1,nx)
    end do
    close(unit)

    ! Read active-domain mask
    open(newunit=unit, &
         file="data/mbabane_domain_mask.txt", &
         status="old", action="read")

    do j = 1, ny
        read(unit,*) (mask(i,j), i=1,nx)
    end do
    close(unit)

    Tinitial = T

    rx = alpha*dt/(dx*dx)
    ry = alpha*dt/(dy*dy)

    if (rx+ry > 0.5d0) then
        print *, "ERROR: Explicit scheme unstable."
        stop 1
    end if

    nactive = count(mask == 1)

    if (nactive == 0) then
        print *, "ERROR: No active cells."
        stop 1
    end if

    mean_initial = sum(T, mask=mask==1) / &
                   real(nactive,dp)

    print *, "===================================="
    print *, "MBABANE LANDSAT HEAT DIFFUSION"
    print *, "===================================="
    print *, "Grid:", nx, "x", ny
    print *, "Spacing (m):", dx, dy
    print *, "Active cells:", nactive
    print *, "Diffusivity:", alpha
    print *, "Time step (s):", dt
    print *, "Number of steps:", nt
    print *, "Stability sum:", rx+ry
    print *, "Initial mean (C):", mean_initial

    ! Explicit finite-difference time stepping
    do step = 1, nt

        Tnew = T

        do j = 1, ny
            do i = 1, nx

                if (mask(i,j) /= 1) cycle

                flux_x = 0.0d0
                flux_y = 0.0d0

                ! East neighbour
                if (i < nx) then
                    if (mask(i+1,j) == 1) then
                        flux_x = flux_x + T(i+1,j)-T(i,j)
                    end if
                end if

                ! West neighbour
                if (i > 1) then
                    if (mask(i-1,j) == 1) then
                        flux_x = flux_x + T(i-1,j)-T(i,j)
                    end if
                end if

                ! South neighbour
                if (j < ny) then
                    if (mask(i,j+1) == 1) then
                        flux_y = flux_y + T(i,j+1)-T(i,j)
                    end if
                end if

                ! North neighbour
                if (j > 1) then
                    if (mask(i,j-1) == 1) then
                        flux_y = flux_y + T(i,j-1)-T(i,j)
                    end if
                end if

                ! No flux through inactive neighbours
                Tnew(i,j) = T(i,j) + &
                            rx*flux_x + ry*flux_y

            end do
        end do

        T = Tnew

    end do

    mean_final = sum(T, mask=mask==1) / &
                 real(nactive,dp)

    max_change = maxval(abs(T-Tinitial), &
                        mask=mask==1)

    print *
    print *, "===================================="
    print *, "SIMULATION RESULTS"
    print *, "===================================="
    print *, "Final mean (C):", mean_final
    print *, "Mean difference (C):", &
             mean_final-mean_initial
    print *, "Maximum absolute change (C):", &
             max_change
    print *, "Initial maximum (C):", &
             maxval(Tinitial, mask=mask==1)
    print *, "Final maximum (C):", &
             maxval(T, mask=mask==1)
    print *, "Initial minimum (C):", &
             minval(Tinitial, mask=mask==1)
    print *, "Final minimum (C):", &
             minval(T, mask=mask==1)

    ! Save final temperature in raster row order
    open(newunit=unit, &
         file="results/mbabane_diffusion_final.txt", &
         status="replace", action="write")

    do j = 1, ny
        write(unit,'(*(ES18.10,1X))') &
             (T(i,j), i=1,nx)
    end do

    close(unit)

    print *
    print *, "Output saved to:"
    print *, "results/mbabane_diffusion_final.txt"

end program mbabane_heat2d

