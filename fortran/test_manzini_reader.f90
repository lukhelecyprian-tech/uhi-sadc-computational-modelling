
program test_manzini_reader
    implicit none

    integer, parameter :: nx = 442
    integer, parameter :: ny = 378

    real(8), allocatable :: temperature(:,:)
    integer :: i, j, ios

    allocate(temperature(nx,ny))

    open(unit=10, &
         file="data/manzini_energy_temperature_K.txt", &
         status="old", action="read", iostat=ios)

    if (ios /= 0) stop "Cannot open temperature file"

    ! Read each geographical row, one value at a time.
    do j = 1, ny
        read(10,*,iostat=ios) (temperature(i,j), i=1,nx)
        if (ios /= 0) stop "Error reading temperature row"
    end do

    close(10)

    print *, "FORTRAN MANZINI READER TEST"
    print *, "--------------------------"
    print *, "Grid columns:", nx
    print *, "Grid rows:", ny

    print *, "First cell:", temperature(1,1)
    print *, "Cell (100,100):", temperature(100,100)
    print *, "Cell (200,150):", temperature(200,150)

    deallocate(temperature)

end program test_manzini_reader
