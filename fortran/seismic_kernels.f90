module seismic_kernels
    use, intrinsic :: iso_c_binding
    implicit none

contains

    !==============================================================================
    ! CPML 2D Anisotropic Elastic Time Step
    !==============================================================================
    subroutine cpml_anisotropic_step_c(nx, ny, deltax, deltay, dt, it, &
                                       isource, jsource, force_x, force_y, &
                                       c11, c12, c22, c33, rho, &
                                       vx, vy, sigmaxx, sigmayy, sigmaxy, &
                                       d_x, d_x_half, d_y, d_y_half, &
                                       k_x, k_x_half, k_y, k_y_half, &
                                       alpha_x, alpha_x_half, alpha_y, alpha_y_half, &
                                       a_x, a_x_half, a_y, a_y_half, &
                                       b_x, b_x_half, b_y, b_y_half, &
                                       memory_dvx_dx, memory_dvx_dy, &
                                       memory_dvy_dx, memory_dvy_dy, &
                                       memory_dsigmaxx_dx, memory_dsigmayy_dy, &
                                       memory_dsigmaxy_dx, memory_dsigmaxy_dy) &
                                       bind(c, name="cpml_anisotropic_step_c")
        
        integer(c_int), value, intent(in) :: nx, ny, it, isource, jsource
        real(c_double), value, intent(in) :: deltax, deltay, dt, force_x, force_y
        
        ! Material properties (2D arrays in column-major order)
        real(c_double), intent(in) :: c11(nx,ny), c12(nx,ny), c22(nx,ny), c33(nx,ny), rho(nx,ny)
        
        ! Wave fields (2D arrays in column-major order)
        real(c_double), intent(inout) :: vx(nx,ny), vy(nx,ny)
        real(c_double), intent(inout) :: sigmaxx(nx,ny), sigmayy(nx,ny), sigmaxy(nx,ny)
        
        ! PML damping profiles (1D arrays)
        real(c_double), intent(in) :: d_x(nx), d_x_half(nx), d_y(ny), d_y_half(ny)
        real(c_double), intent(in) :: k_x(nx), k_x_half(nx), k_y(ny), k_y_half(ny)
        real(c_double), intent(in) :: alpha_x(nx), alpha_x_half(nx), alpha_y(ny), alpha_y_half(ny)
        real(c_double), intent(in) :: a_x(nx), a_x_half(nx), a_y(ny), a_y_half(ny)
        real(c_double), intent(in) :: b_x(nx), b_x_half(nx), b_y(ny), b_y_half(ny)
        
        ! PML memory variables (2D arrays)
        real(c_double), intent(inout) :: memory_dvx_dx(nx,ny), memory_dvx_dy(nx,ny)
        real(c_double), intent(inout) :: memory_dvy_dx(nx,ny), memory_dvy_dy(nx,ny)
        real(c_double), intent(inout) :: memory_dsigmaxx_dx(nx,ny), memory_dsigmayy_dy(nx,ny)
        real(c_double), intent(inout) :: memory_dsigmaxy_dx(nx,ny), memory_dsigmaxy_dy(nx,ny)
        
        ! Local variables
        integer :: i, j
        real(c_double) :: value_dvx_dx, value_dvx_dy, value_dvy_dx, value_dvy_dy
        real(c_double) :: value_dsigmaxx_dx, value_dsigmayy_dy, value_dsigmaxy_dx, value_dsigmaxy_dy
        real(c_double) :: DELTAT_over_rho, DELTAT_over_DELTAX, DELTAT_over_DELTAY
        
        DELTAT_over_DELTAX = dt / deltax
        DELTAT_over_DELTAY = dt / deltay
        
        !------------------------------------------------------------
        ! compute stress sigma and update memory variables for C-PML
        !------------------------------------------------------------
        
        do j = 2, ny
            do i = 1, nx-1
                
                value_dvx_dx = (vx(i+1,j) - vx(i,j)) / deltax
                value_dvy_dy = (vy(i,j) - vy(i,j-1)) / deltay
                
                memory_dvx_dx(i,j) = b_x_half(i) * memory_dvx_dx(i,j) + a_x_half(i) * value_dvx_dx
                memory_dvy_dy(i,j) = b_y(j) * memory_dvy_dy(i,j) + a_y(j) * value_dvy_dy
                
                value_dvx_dx = value_dvx_dx / k_x_half(i) + memory_dvx_dx(i,j)
                value_dvy_dy = value_dvy_dy / k_y(j) + memory_dvy_dy(i,j)
                
                sigmaxx(i,j) = sigmaxx(i,j) + (c11(i,j) * value_dvx_dx + c12(i,j) * value_dvy_dy) * dt
                sigmayy(i,j) = sigmayy(i,j) + (c12(i,j) * value_dvx_dx + c22(i,j) * value_dvy_dy) * dt
                
            enddo
        enddo
        
        do j = 1, ny-1
            do i = 2, nx
                
                value_dvy_dx = (vy(i,j) - vy(i-1,j)) / deltax
                value_dvx_dy = (vx(i,j+1) - vx(i,j)) / deltay
                
                memory_dvy_dx(i,j) = b_x(i) * memory_dvy_dx(i,j) + a_x(i) * value_dvy_dx
                memory_dvx_dy(i,j) = b_y_half(j) * memory_dvx_dy(i,j) + a_y_half(j) * value_dvx_dy
                
                value_dvy_dx = value_dvy_dx / k_x(i) + memory_dvy_dx(i,j)
                value_dvx_dy = value_dvx_dy / k_y_half(j) + memory_dvx_dy(i,j)
                
                sigmaxy(i,j) = sigmaxy(i,j) + c33(i,j) * (value_dvy_dx + value_dvx_dy) * dt
                
            enddo
        enddo
        
        !--------------------------------------------------------
        ! compute velocity and update memory variables for C-PML
        !--------------------------------------------------------
        
        do j = 2, ny
            do i = 2, nx
                
                value_dsigmaxx_dx = (sigmaxx(i,j) - sigmaxx(i-1,j)) / deltax
                value_dsigmaxy_dy = (sigmaxy(i,j) - sigmaxy(i,j-1)) / deltay
                
                memory_dsigmaxx_dx(i,j) = b_x(i) * memory_dsigmaxx_dx(i,j) + a_x(i) * value_dsigmaxx_dx
                memory_dsigmaxy_dy(i,j) = b_y(j) * memory_dsigmaxy_dy(i,j) + a_y(j) * value_dsigmaxy_dy
                
                value_dsigmaxx_dx = value_dsigmaxx_dx / k_x(i) + memory_dsigmaxx_dx(i,j)
                value_dsigmaxy_dy = value_dsigmaxy_dy / k_y(j) + memory_dsigmaxy_dy(i,j)
                
                DELTAT_over_rho = dt / rho(i,j)
                vx(i,j) = vx(i,j) + (value_dsigmaxx_dx + value_dsigmaxy_dy) * DELTAT_over_rho
                
            enddo
        enddo
        
        do j = 1, ny-1
            do i = 1, nx-1
                
                value_dsigmaxy_dx = (sigmaxy(i+1,j) - sigmaxy(i,j)) / deltax
                value_dsigmayy_dy = (sigmayy(i,j+1) - sigmayy(i,j)) / deltay
                
                memory_dsigmaxy_dx(i,j) = b_x_half(i) * memory_dsigmaxy_dx(i,j) + a_x_half(i) * value_dsigmaxy_dx
                memory_dsigmayy_dy(i,j) = b_y_half(j) * memory_dsigmayy_dy(i,j) + a_y_half(j) * value_dsigmayy_dy
                
                value_dsigmaxy_dx = value_dsigmaxy_dx / k_x_half(i) + memory_dsigmaxy_dx(i,j)
                value_dsigmayy_dy = value_dsigmayy_dy / k_y_half(j) + memory_dsigmayy_dy(i,j)
                
                DELTAT_over_rho = dt / rho(i,j)
                vy(i,j) = vy(i,j) + (value_dsigmaxy_dx + value_dsigmayy_dy) * DELTAT_over_rho
                
            enddo
        enddo
        
        ! Add the source (force vector located at a given grid point)
        if (isource >= 1 .and. isource <= nx .and. jsource >= 1 .and. jsource <= ny) then
            vx(isource,jsource) = vx(isource,jsource) + force_x * dt / rho(isource,jsource)
            vy(isource,jsource) = vy(isource,jsource) + force_y * dt / rho(isource,jsource)
        endif
        
    end subroutine cpml_anisotropic_step_c

    !==============================================================================
    ! ADE-PML 2D Elastic RK4 8th Order Time Step  
    !==============================================================================
    subroutine adepml_rk4_step_c(nx, ny, deltax, deltay, dt, rk_substep, &
                                  isource, jsource, force_x, force_y, &
                                  lambda, mu, rho, &
                                  vx, vy, sigmaxx, sigmayy, sigmaxy, &
                                  dvx, dvy, dsigmaxx, dsigmayy, dsigmaxy, &
                                  d_x, d_x_half, d_y, d_y_half, &
                                  k_x, k_x_half, k_y, k_y_half, &
                                  alpha_x, alpha_x_half, alpha_y, alpha_y_half, &
                                  a_x_rk4, a_x_half_rk4, a_y_rk4, a_y_half_rk4, &
                                  b_x_rk4, b_x_half_rk4, b_y_rk4, b_y_half_rk4, &
                                  memory_dvx_dx, memory_dvx_dy, &
                                  memory_dvy_dx, memory_dvy_dy, &
                                  memory_dsigmaxx_dx, memory_dsigmayy_dy, &
                                  memory_dsigmaxy_dx, memory_dsigmaxy_dy) &
                                  bind(c, name="adepml_rk4_step_c")
        
        integer(c_int), value, intent(in) :: nx, ny, rk_substep, isource, jsource
        real(c_double), value, intent(in) :: deltax, deltay, dt, force_x, force_y
        
        ! Material properties (extended grid: -4 to nx+4, -4 to ny+4)
        real(c_double), intent(in) :: lambda(-4:nx+4,-4:ny+4), mu(-4:nx+4,-4:ny+4), rho(-4:nx+4,-4:ny+4)
        
        ! Wave fields (extended grid)
        real(c_double), intent(inout) :: vx(-4:nx+4,-4:ny+4), vy(-4:nx+4,-4:ny+4)
        real(c_double), intent(inout) :: sigmaxx(-4:nx+4,-4:ny+4), sigmayy(-4:nx+4,-4:ny+4), sigmaxy(-4:nx+4,-4:ny+4)
        
        ! RK4 intermediate arrays
        real(c_double), intent(inout) :: dvx(3,-4:nx+4,-4:ny+4), dvy(3,-4:nx+4,-4:ny+4)
        real(c_double), intent(inout) :: dsigmaxx(3,-4:nx+4,-4:ny+4), dsigmayy(3,-4:nx+4,-4:ny+4), dsigmaxy(3,-4:nx+4,-4:ny+4)
        
        ! PML damping profiles
        real(c_double), intent(in) :: d_x(-4:nx+4), d_x_half(-4:nx+4), d_y(-4:ny+4), d_y_half(-4:ny+4)
        real(c_double), intent(in) :: k_x(-4:nx+4), k_x_half(-4:nx+4), k_y(-4:ny+4), k_y_half(-4:ny+4)
        real(c_double), intent(in) :: alpha_x(-4:nx+4), alpha_x_half(-4:nx+4), alpha_y(-4:ny+4), alpha_y_half(-4:ny+4)
        
        ! ADE-PML RK4 coefficients (2D arrays: [4 substeps][spatial points])
        real(c_double), intent(in) :: a_x_rk4(4,-4:nx+4), a_x_half_rk4(4,-4:nx+4), a_y_rk4(4,-4:ny+4), a_y_half_rk4(4,-4:ny+4)
        real(c_double), intent(in) :: b_x_rk4(4,-4:nx+4), b_x_half_rk4(4,-4:nx+4), b_y_rk4(4,-4:ny+4), b_y_half_rk4(4,-4:ny+4)
        
        ! PML memory variables
        real(c_double), intent(inout) :: memory_dvx_dx(-4:nx+4,-4:ny+4), memory_dvx_dy(-4:nx+4,-4:ny+4)
        real(c_double), intent(inout) :: memory_dvy_dx(-4:nx+4,-4:ny+4), memory_dvy_dy(-4:nx+4,-4:ny+4)
        real(c_double), intent(inout) :: memory_dsigmaxx_dx(-4:nx+4,-4:ny+4), memory_dsigmayy_dy(-4:nx+4,-4:ny+4)
        real(c_double), intent(inout) :: memory_dsigmaxy_dx(-4:nx+4,-4:ny+4), memory_dsigmaxy_dy(-4:nx+4,-4:ny+4)
        
        ! Local variables
        integer :: i, j
        real(c_double) :: value_dvx_dx, value_dvx_dy, value_dvy_dx, value_dvy_dy
        real(c_double) :: value_dsigmaxx_dx, value_dsigmayy_dy, value_dsigmaxy_dx, value_dsigmaxy_dy
        real(c_double) :: lambda_half_x, mu_half_x, lambda_plus_two_mu_half_x, mu_half_y, rho_half_x_half_y
        
        ! Holberg coefficients for 8th order spatial derivatives
        real(c_double), parameter :: c1 = 1.231666d0
        real(c_double), parameter :: c2 = -1.041182d-1
        real(c_double), parameter :: c3 = 2.063707d-2
        real(c_double), parameter :: c4 = -3.570998d-3
        
        ! RK4 coefficients
        real(c_double), parameter :: rk41(4) = [0.0d0, 0.5d0, 0.5d0, 1.0d0]
        real(c_double), parameter :: rk42(4) = [1.0d0/6.0d0, 2.0d0/6.0d0, 2.0d0/6.0d0, 1.0d0/6.0d0]
        
        real(c_double) :: dt_rk
        
        dt_rk = dt * rk42(rk_substep)
        
        !------------------------------------------------------------
        ! compute stress sigma and update memory variables for ADE-PML
        ! using 8th-order Holberg spatial differences
        !------------------------------------------------------------
        
        do j = 1, ny
            do i = 1, nx-1
                
                ! 8th-order spatial derivative using Holberg coefficients
                value_dvx_dx = (c1*(vx(i+1,j) - vx(i,j)) + &
                               c2*(vx(i+2,j) - vx(i-1,j)) + &
                               c3*(vx(i+3,j) - vx(i-2,j)) + &
                               c4*(vx(i+4,j) - vx(i-3,j))) / deltax
                
                value_dvy_dy = (c1*(vy(i,j+1) - vy(i,j)) + &
                               c2*(vy(i,j+2) - vy(i,j-1)) + &
                               c3*(vy(i,j+3) - vy(i,j-2)) + &
                               c4*(vy(i,j+4) - vy(i,j-3))) / deltay
                
                ! ADE-PML memory variable updates using RK4-specific coefficients
                memory_dvx_dx(i,j) = b_x_half_rk4(rk_substep,i) * memory_dvx_dx(i,j) + &
                                    a_x_half_rk4(rk_substep,i) * value_dvx_dx
                memory_dvy_dy(i,j) = b_y_rk4(rk_substep,j) * memory_dvy_dy(i,j) + &
                                    a_y_rk4(rk_substep,j) * value_dvy_dy
                
                value_dvx_dx = value_dvx_dx / k_x_half(i) + memory_dvx_dx(i,j)
                value_dvy_dy = value_dvy_dy / k_y(j) + memory_dvy_dy(i,j)
                
                ! Interpolate material parameters
                lambda_half_x = 0.5d0 * (lambda(i+1,j) + lambda(i,j))
                mu_half_x = 0.5d0 * (mu(i+1,j) + mu(i,j))
                lambda_plus_two_mu_half_x = lambda_half_x + 2.0d0 * mu_half_x
                
                ! Update stress components
                dsigmaxx(1,i,j) = lambda_plus_two_mu_half_x * value_dvx_dx + lambda_half_x * value_dvy_dy
                dsigmayy(1,i,j) = lambda_half_x * value_dvx_dx + lambda_plus_two_mu_half_x * value_dvy_dy
                
            enddo
        enddo
        
        do j = 1, ny-1
            do i = 1, nx
                
                ! 8th-order spatial derivatives
                value_dvy_dx = (c1*(vy(i,j) - vy(i-1,j)) + &
                               c2*(vy(i+1,j) - vy(i-2,j)) + &
                               c3*(vy(i+2,j) - vy(i-3,j)) + &
                               c4*(vy(i+3,j) - vy(i-4,j))) / deltax
                
                value_dvx_dy = (c1*(vx(i,j+1) - vx(i,j)) + &
                               c2*(vx(i,j+2) - vx(i,j-1)) + &
                               c3*(vx(i,j+3) - vx(i,j-2)) + &
                               c4*(vx(i,j+4) - vx(i,j-3))) / deltay
                
                ! ADE-PML memory variable updates using RK4-specific coefficients
                memory_dvy_dx(i,j) = b_x_rk4(rk_substep,i) * memory_dvy_dx(i,j) + &
                                    a_x_rk4(rk_substep,i) * value_dvy_dx
                memory_dvx_dy(i,j) = b_y_half_rk4(rk_substep,j) * memory_dvx_dy(i,j) + &
                                    a_y_half_rk4(rk_substep,j) * value_dvx_dy
                
                value_dvy_dx = value_dvy_dx / k_x(i) + memory_dvy_dx(i,j)
                value_dvx_dy = value_dvx_dy / k_y_half(j) + memory_dvx_dy(i,j)
                
                ! Interpolate material parameters
                mu_half_y = 0.5d0 * (mu(i,j+1) + mu(i,j))
                
                ! Update shear stress
                dsigmaxy(1,i,j) = mu_half_y * (value_dvy_dx + value_dvx_dy)
                
            enddo
        enddo
        
        !--------------------------------------------------------
        ! compute velocity and update memory variables for ADE-PML
        !--------------------------------------------------------
        
        do j = 1, ny
            do i = 1, nx
                
                ! 8th-order spatial derivatives for stress
                value_dsigmaxx_dx = (c1*(sigmaxx(i,j) - sigmaxx(i-1,j)) + &
                                    c2*(sigmaxx(i+1,j) - sigmaxx(i-2,j)) + &
                                    c3*(sigmaxx(i+2,j) - sigmaxx(i-3,j)) + &
                                    c4*(sigmaxx(i+3,j) - sigmaxx(i-4,j))) / deltax
                
                value_dsigmaxy_dy = (c1*(sigmaxy(i,j) - sigmaxy(i,j-1)) + &
                                    c2*(sigmaxy(i,j+1) - sigmaxy(i,j-2)) + &
                                    c3*(sigmaxy(i,j+2) - sigmaxy(i,j-3)) + &
                                    c4*(sigmaxy(i,j+3) - sigmaxy(i,j-4))) / deltay
                
                ! ADE-PML memory variable updates using RK4-specific coefficients
                memory_dsigmaxx_dx(i,j) = b_x_rk4(rk_substep,i) * memory_dsigmaxx_dx(i,j) + &
                                         a_x_rk4(rk_substep,i) * value_dsigmaxx_dx
                memory_dsigmaxy_dy(i,j) = b_y_rk4(rk_substep,j) * memory_dsigmaxy_dy(i,j) + &
                                         a_y_rk4(rk_substep,j) * value_dsigmaxy_dy
                
                value_dsigmaxx_dx = value_dsigmaxx_dx / k_x(i) + memory_dsigmaxx_dx(i,j)
                value_dsigmaxy_dy = value_dsigmaxy_dy / k_y(j) + memory_dsigmaxy_dy(i,j)
                
                ! Update x-velocity
                dvx(1,i,j) = (value_dsigmaxx_dx + value_dsigmaxy_dy) / rho(i,j)
                
            enddo
        enddo
        
        do j = 1, ny
            do i = 1, nx
                
                ! 8th-order spatial derivatives
                value_dsigmaxy_dx = (c1*(sigmaxy(i+1,j) - sigmaxy(i,j)) + &
                                    c2*(sigmaxy(i+2,j) - sigmaxy(i-1,j)) + &
                                    c3*(sigmaxy(i+3,j) - sigmaxy(i-2,j)) + &
                                    c4*(sigmaxy(i+4,j) - sigmaxy(i-3,j))) / deltax
                
                value_dsigmayy_dy = (c1*(sigmayy(i,j+1) - sigmayy(i,j)) + &
                                    c2*(sigmayy(i,j+2) - sigmayy(i,j-1)) + &
                                    c3*(sigmayy(i,j+3) - sigmayy(i,j-2)) + &
                                    c4*(sigmayy(i,j+4) - sigmayy(i,j-3))) / deltay
                
                ! ADE-PML memory variable updates using RK4-specific coefficients
                memory_dsigmaxy_dx(i,j) = b_x_half_rk4(rk_substep,i) * memory_dsigmaxy_dx(i,j) + &
                                          a_x_half_rk4(rk_substep,i) * value_dsigmaxy_dx
                memory_dsigmayy_dy(i,j) = b_y_half_rk4(rk_substep,j) * memory_dsigmayy_dy(i,j) + &
                                          a_y_half_rk4(rk_substep,j) * value_dsigmayy_dy
                
                value_dsigmaxy_dx = value_dsigmaxy_dx / k_x_half(i) + memory_dsigmaxy_dx(i,j)
                value_dsigmayy_dy = value_dsigmayy_dy / k_y_half(j) + memory_dsigmayy_dy(i,j)
                
                ! Update y-velocity
                dvy(1,i,j) = (value_dsigmaxy_dx + value_dsigmayy_dy) / rho(i,j)
                
            enddo
        enddo
        
        ! Apply RK4 updates to the wave fields
        ! Store computed derivatives in dvx(1,:,:), dvy(1,:,:), etc.
        ! These will be used by C++ to perform the RK4 accumulation
        ! DO NOT update vx, vy, sigmaxx, sigmayy, sigmaxy here!
        do j = 1, ny
            do i = 1, nx
                ! Store the computed derivatives for RK4 accumulation in C++
                ! dvx(1,i,j), dvy(1,i,j), etc. already contain the computed derivatives
                ! from the main computation loop above
            enddo
        enddo
        
        ! Note: Source addition is now handled in C++ RK4Integrator
        ! to maintain proper RK4 time integration
        
    end subroutine adepml_rk4_step_c

end module seismic_kernels
