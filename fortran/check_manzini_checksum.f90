
program check_manzini_checksum
    use iso_fortran_env, only: real64
    implicit none

    integer, parameter :: nx = 442, ny = 378

    real(real64), allocatable :: temperature(:,:)
    integer, allocatable :: mask(:,:)

    real(real64) :: total, mean_temp, sum_squares
    integer :: i, j, nvalid, ios

    allocate(temperature(nx,ny))
    allocate(mask(nx,ny))

    ! Read temperature data row by row.
    open(10, file="data/manzini_energy_temperature_K.txt", &
         status="old", action="read", iostat=ios)

    if (ios /= 0) stop "Cannot open temperature file"

    do j = 1, ny
        read(10,*,iostat=ios) (temperature(i,j), i=1,nx)
        if (ios /= 0) stop "Temperature reading error"
    end do

    close(10)

    ! Read the corresponding common mask.
    open(11, file="data/manzini_energy_mask.txt", &
         status="old", action="read", iostat=ios)

    if (ios /= 0) stop "Cannot open mask file"

    do j = 1, ny
        read(11,*,iostat=ios) (mask(i,j), i=1,nx)
        if (ios /= 0) stop "Mask reading error"
    end do

    close(11)

    total = 0.0_real64
    sum_squares = 0.0_real64
    nvalid = 0

    do j = 1, ny
        do i = 1, nx
            if (mask(i,j) == 1) then
                nvalid = nvalid + 1

                total = total + temperature(i,j)

                sum_squares = sum_squares + &
                    temperature(i,j)**2
            end if
        end do
    end do

    if (nvalid == 0) stop "No valid cells found"

    mean_temp = total / real(nvalid, real64)

    print *, "MANZINI FORTRAN CHECKSUM"
    print *, "-----------------------"
    print *, "Valid pixels:", nvalid
    print *, "Temperature sum (K):", total
    print *, "Mean temperature (K):", mean_temp
    print *, "Sum of temperature squares:", sum_squares

    deallocate(temperature, mask)

end program check_manzini_checksum
