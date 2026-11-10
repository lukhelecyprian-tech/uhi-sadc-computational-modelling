
program verify_coupled_spatial_operator

    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: ngrids = 5

    integer, parameter :: sizes(ngrids) = &
        [20, 40, 80, 160, 320]

    real(dp), parameter :: pi = acos(-1.0d0)
    real(dp), parameter :: domain_length = 1200.0d0
    real(dp), parameter :: t0 = 293.15d0
    real(dp), parameter :: amplitude = 5.0d0

    integer :: g, nx, ny, i, j, ios

    real(dp) :: dx, dy, x, y
    real(dp) :: lap_num, lap_exact
    real(dp) :: error, rmse, max_error
    real(dp) :: errors(ngrids), order
    real(dp) :: left_t, right_t
    real(dp) :: bottom_t, top_t

    real(dp), allocatable :: temp(:,:)

    ! =====================================================
    ! 1. OPEN OUTPUT CSV
    ! =====================================================

    open(20, &
        file="results/coupled_spatial_operator_verification.csv", &
        status="replace", action="write", iostat=ios)

    if (ios /= 0) error stop "Cannot create CSV"

    write(20,'(A)') &
        "grid_n,dx_m,rmse_K_per_m2,max_error_K_per_m2"

    print *
    print *, "=============================================="
    print *, "COUPLED MODEL SPATIAL OPERATOR VERIFICATION"
    print *, "=============================================="
    print *, "Testing five-point Laplacian"
    print *, "Cell-centred grid, zero-flux boundaries"
    print *
    print *, "Grid       dx(m)          RMSE           Max"
    print *, "------------------------------------------------"

    ! =====================================================
    ! 2. GRID REFINEMENT
    ! =====================================================

    do g = 1, ngrids

        nx = sizes(g)
        ny = nx

        dx = domain_length / real(nx,dp)
        dy = dx

        allocate(temp(nx,ny))

        ! -------------------------------------------------
        ! Initialise analytical temperature field
        ! -------------------------------------------------

        do j = 1, ny

            y = (real(j,dp)-0.5d0)*dy

            do i = 1, nx

                x = (real(i,dp)-0.5d0)*dx

                temp(i,j) = exact_temperature(x,y)

            end do
        end do

        rmse = 0d0
        max_error = 0d0

        ! -------------------------------------------------
        ! Calculate numerical and exact Laplacians
        ! -------------------------------------------------

        do j = 1, ny

            y = (real(j,dp)-0.5d0)*dy

            do i = 1, nx

                x = (real(i,dp)-0.5d0)*dx

                ! At physical boundaries, use the cell
                ! temperature itself as the ghost value.
                ! This implements zero normal flux.

                left_t = temp(i,j)
                right_t = temp(i,j)
                bottom_t = temp(i,j)
                top_t = temp(i,j)

                if (i > 1) left_t = temp(i-1,j)
                if (i < nx) right_t = temp(i+1,j)

                if (j > 1) bottom_t = temp(i,j-1)
                if (j < ny) top_t = temp(i,j+1)

                lap_num = &
                    (left_t - 2d0*temp(i,j) + right_t)/dx**2 &
                    + &
                    (bottom_t - 2d0*temp(i,j) + top_t)/dy**2

                lap_exact = exact_laplacian(x,y)

                error = lap_num - lap_exact

                rmse = rmse + error**2
                max_error = max(max_error,abs(error))

            end do
        end do

        rmse = sqrt(rmse/real(nx*ny,dp))
        errors(g) = rmse

        print '(I5,4X,F10.4,4X,ES13.5,4X,ES13.5)', &
            nx, dx, rmse, max_error

        write(20,'(I6,A,F12.6,A,ES18.10,A,ES18.10)') &
            nx,",",dx,",",rmse,",",max_error

        deallocate(temp)

    end do

    ! =====================================================
    ! 3. CALCULATE SPATIAL CONVERGENCE ORDER
    ! =====================================================

    print *
    print *, "Estimated spatial convergence orders"
    print *, "----------------------------------------------"

    do g = 1, ngrids-1

        if (errors(g) > 0d0 .and. errors(g+1) > 0d0) then

            order = log(errors(g)/errors(g+1))/log(2d0)

            print '(I5,A,I5,A,F12.6)', &
                sizes(g)," -> ",sizes(g+1),": ",order

        else

            print *, "Unable to estimate order"

        end if

    end do

    print *
    print *, "Expected spatial convergence order: 2"
    print *, "Results saved to:"
    print *, "results/coupled_spatial_operator_verification.csv"
    print *
    print *, "Verification completed."

contains

    ! =====================================================
    ! EXACT TEMPERATURE FIELD
    ! =====================================================

    function exact_temperature(x,y) result(value)

        real(dp), intent(in) :: x,y
        real(dp) :: value

        value = t0 + amplitude * &
            cos(pi*x/domain_length) * &
            cos(pi*y/domain_length)

    end function exact_temperature

    ! =====================================================
    ! EXACT CONTINUOUS LAPLACIAN
    ! =====================================================

    function exact_laplacian(x,y) result(value)

        real(dp), intent(in) :: x,y
        real(dp) :: value

        value = -2d0*pi*pi*amplitude / &
            domain_length**2 * &
            cos(pi*x/domain_length) * &
            cos(pi*y/domain_length)

    end function exact_laplacian

end program verify_coupled_spatial_operator
