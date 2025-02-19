!>
!! @file m_time_steppers.f90
!! @brief Contains module m_time_steppers

#:include 'macros.fpp'

!> @brief The following module features a variety of time-stepping schemes.
!!              Currently, it includes the following Runge-Kutta (RK) algorithms:
!!                   1) 1st Order TVD RK
!!                   2) 2nd Order TVD RK
!!                   3) 3rd Order TVD RK
!!              where TVD designates a total-variation-diminishing time-stepper.
module m_time_steppers

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_rhs                  !< Right-hane-side (RHS) evaluation procedures

    use m_data_output          !< Run-time info & solution data output procedures

    use m_bubbles_EE           !< Ensemble-averaged bubble dynamics routines

    use m_bubbles_EL           !< Lagrange bubble dynamics routines

    use m_ibm

    use m_hyperelastic

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_boundary_conditions

    use m_helper

    use m_sim_helpers

    use m_fftw

    use m_nvtx

    use m_thermochem, only: num_species

    use m_body_forces

    implicit none

    type(vector_field), allocatable, dimension(:) :: q_cons_ts !<
    !! Cell-average conservative variables at each time-stage (TS)

    type(scalar_field), allocatable, dimension(:) :: q_prim_vf !<
    !! Cell-average primitive variables at the current time-stage

    type(scalar_field), allocatable, dimension(:) :: rhs_vf !<
    !! Cell-average RHS variables at the current time-stage

    type(vector_field), allocatable, dimension(:) :: rhs_ts_rkck
    !! Cell-average RHS variables at each time-stage (TS)
    !! Adaptive 4th/5th order Runge—Kutta–Cash–Karp (RKCK) time stepper

    type(vector_field), allocatable, dimension(:) :: q_prim_ts !<
    !! Cell-average primitive variables at consecutive TIMESTEPS

    real(wp), allocatable, dimension(:, :, :, :, :) :: rhs_pb

    type(scalar_field) :: q_T_sf !<
    !! Cell-average temperature variables at the current time-stage

    real(wp), allocatable, dimension(:, :, :, :, :) :: rhs_mv

    real(wp), allocatable, dimension(:, :, :) :: max_dt

    integer, private :: num_ts !<
    !! Number of time stages in the time-stepping scheme

    ! drag calculation variables 
    type(scalar_field), allocatable, dimension(:) :: rhs_rhouu
    type(vector_field), allocatable, dimension(:) :: du_dxyz ! dudx, dudy, dudz 
    real(wp), allocatable, dimension(:) :: F_D_vi
    real(wp), allocatable, dimension(:) :: F_D_si
    real(wp) :: x_surf, y_surf, z_surf, pressure_surf
    real(wp), allocatable, dimension(:) :: dudx_surf, dudy_surf, dudz_surf, F_D_mat
    real(wp), allocatable, dimension(:, :) :: stress_tensor

    ! periodic forcing variables
    type(scalar_field), allocatable, dimension(:) :: q_bar
    type(scalar_field), allocatable, dimension(:) :: q_periodic_force
    real(wp), allocatable, dimension(:) :: q_spatial_avg, q_spatial_avg_glb
    real(wp) :: N_x_total_glb, volfrac_phi, rho_inf, u_inf, T_inf

    !$acc declare create(q_cons_ts, q_prim_vf, q_T_sf, rhs_vf, rhs_ts_rkck, q_prim_ts, rhs_mv, rhs_pb, max_dt)
    !$acc declare create(rhs_rhouu, du_dxyz, F_D_vi, F_D_si)
    !$acc declare create(x_surf, y_surf, z_surf, pressure_surf, dudx_surf, dudy_surf, dudz_surf, stress_tensor)
    !$acc declare create(q_bar, q_periodic_force, q_spatial_avg, q_spatial_avg_glb, N_x_total_glb, volfrac_phi, rho_inf, u_inf, T_inf)

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_time_steppers_module

        integer :: i, j !< Generic loop iterators

        ! Setting number of time-stages for selected time-stepping scheme
        if (time_stepper == 1) then
            num_ts = 1
        elseif (any(time_stepper == (/2, 3, 4/))) then
            num_ts = 2
        end if

        ! Allocating the cell-average conservative variables
        @:ALLOCATE(q_cons_ts(1:num_ts))

        do i = 1, num_ts
            @:ALLOCATE(q_cons_ts(i)%vf(1:sys_size))
        end do

        do i = 1, num_ts
            do j = 1, sys_size
                @:ALLOCATE(q_cons_ts(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
            end do
            @:ACC_SETUP_VFs(q_cons_ts(i))
        end do

        ! Allocating the cell-average primitive ts variables
        if (probe_wrt) then
            @:ALLOCATE(q_prim_ts(0:3))

            do i = 0, 3
                @:ALLOCATE(q_prim_ts(i)%vf(1:sys_size))
            end do

            do i = 0, 3
                do j = 1, sys_size
                    @:ALLOCATE(q_prim_ts(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                end do
            end do

            do i = 0, 3
                @:ACC_SETUP_VFs(q_prim_ts(i))
            end do
        end if

        ! Allocating the cell-average primitive variables
        @:ALLOCATE(q_prim_vf(1:sys_size))

        do i = 1, adv_idx%end
            @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(q_prim_vf(i))
        end do

        if (bubbles_euler) then
            do i = bub_idx%beg, bub_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(i))
            end do
            if (adv_n) then
                @:ALLOCATE(q_prim_vf(n_idx)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(n_idx))
            end if
        end if

        if (elasticity) then
            do i = stress_idx%beg, stress_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(i))
            end do
        end if

        if (hyperelasticity) then
            do i = xibeg, xiend + 1
                @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(i))
            end do
        end if

        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(i))
            end do
        end if

        if (surface_tension) then
            @:ALLOCATE(q_prim_vf(c_idx)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(q_prim_vf(c_idx))
        end if

        if (chemistry) then
            do i = chemxb, chemxe
                @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(i))
            end do

            @:ALLOCATE(q_T_sf%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(q_T_sf)
        end if

        @:ALLOCATE(pb_ts(1:2))
        !Initialize bubble variables pb and mv at all quadrature nodes for all R0 bins
        if (qbmm .and. (.not. polytropic)) then
            @:ALLOCATE(pb_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(1))

            @:ALLOCATE(pb_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(2))

            @:ALLOCATE(rhs_pb(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))
        else if (qbmm .and. polytropic) then
            @:ALLOCATE(pb_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(1))

            @:ALLOCATE(pb_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(2))

            @:ALLOCATE(rhs_pb(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
        end if

        @:ALLOCATE(mv_ts(1:2))

        if (qbmm .and. (.not. polytropic)) then
            @:ALLOCATE(mv_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(1))

            @:ALLOCATE(mv_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(2))

            @:ALLOCATE(rhs_mv(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))

        else if (qbmm .and. polytropic) then
            @:ALLOCATE(mv_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(1))

            @:ALLOCATE(mv_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(2))

            @:ALLOCATE(rhs_mv(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
        end if

        ! Allocating the cell-average RHS time-stages for adaptive RKCK stepper
        if (bubbles_lagrange .and. time_stepper == 4) then
            @:ALLOCATE(rhs_ts_rkck(1:num_ts_rkck))
            do i = 1, num_ts_rkck
                @:ALLOCATE(rhs_ts_rkck(i)%vf(1:sys_size))
            end do
            do i = 1, num_ts_rkck
                do j = 1, sys_size
                    @:ALLOCATE(rhs_ts_rkck(i)%vf(j)%sf(0:m, 0:n, 0:p))
                end do
                @:ACC_SETUP_VFs(rhs_ts_rkck(i))
            end do
        else
            ! Allocating the cell-average RHS variables
            @:ALLOCATE(rhs_vf(1:sys_size))

            do i = 1, sys_size
                @:ALLOCATE(rhs_vf(i)%sf(0:m, 0:n, 0:p))
                @:ACC_SETUP_SFs(rhs_vf(i))
            end do
        end if

        ! allocating drag calculation variables
        if (compute_CD_vi) then
            @:ALLOCATE(rhs_rhouu(1:sys_size))

            do i = 1, sys_size
                @:ALLOCATE(rhs_rhouu(i)%sf(0:m, 0:n, 0:p))
                @:ACC_SETUP_SFs(rhs_rhouu(i))
            end do

            @:ALLOCATE(F_D_vi(num_ibs))
        end if

        if (compute_CD_si) then
            @:ALLOCATE(du_dxyz(1:3))

            do i = 1, 3
                @:ALLOCATE(du_dxyz(i)%vf(1:3))
            end do

            do i = 1, 3
                do j = 1, 3
                    @:ALLOCATE(du_dxyz(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                end do
                @:ACC_SETUP_VFs(du_dxyz(i))
            end do

            @:ALLOCATE(F_D_si(num_ibs))

            @:ALLOCATE(dudx_surf(1:3))
            @:ALLOCATE(dudy_surf(1:3))
            @:ALLOCATE(dudz_surf(1:3))

            @:ALLOCATE(stress_tensor(1:3, 1:3))
            @:ALLOCATE(F_D_mat(1:3))
        end if

        ! allocating periodic forcing variables
        if (periodic_forcing) then
            @:ALLOCATE(q_bar(1:5))
            @:ALLOCATE(q_periodic_force(1:8))

            do i = 1, 5
                @:ALLOCATE(q_bar(i)%sf(0:m, 0:n, 0:p))
                @:ACC_SETUP_SFs(q_bar(i))
            end do

            do i = 1, 8
                @:ALLOCATE(q_periodic_force(i)%sf(0:m, 0:n, 0:p))
                @:ACC_SETUP_SFs(q_periodic_force(i))
            end do 
            
            @:ALLOCATE(q_spatial_avg(1:5))
            @:ALLOCATE(q_spatial_avg_glb(1:5))

            N_x_total_glb = (m_glb + 1) * (n_glb + 1) * (p_glb + 1)
            !$acc update device(N_x_total_glb)

            volfrac_phi = num_ibs * 4._wp/3._wp * pi * patch_ib(1)%radius**3.0 / ((x_domain%end - x_domain%beg)*(y_domain%end - y_domain%beg)*(z_domain%end - z_domain%beg))
            !$acc update device(volfrac_phi)
        end if

        ! Opening and writing the header of the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_open_run_time_information_file()
        end if

        if (cfl_dt) then
            @:ALLOCATE(max_dt(0:m, 0:n, 0:p))
        end if

        ! clear drag output file for use in m_time_steppers
        open(unit=100, file='FD_vi.bin', status='replace', form='unformatted', access='stream')
        open(unit=101, file='FD_si.bin', status='replace', form='unformatted', access='stream')
        open(unit=102, file='xmom_spatialavg.bin', status='replace', form='unformatted', access='stream')

    end subroutine s_initialize_time_steppers_module

    !> 1st order TVD RK time-stepping algorithm
        !! @param t_step Current time step
    subroutine s_1st_order_tvd_rk(t_step, time_avg)

        integer, intent(in) :: t_step
        real(wp), intent(inout) :: time_avg

        integer :: i, j, k, l, q !< Generic loop iterator

        ! Stage 1 of 1
        call nvtxStartRange("TIMESTEP")

        call s_compute_rhs(q_cons_ts(1)%vf, q_T_sf, q_prim_vf, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, &
        rhs_rhouu, du_dxyz)

#ifdef DEBUG
        print *, 'got rhs'
#endif

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

#ifdef DEBUG
        print *, 'wrote runtime info'
#endif

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (cfl_dt) then
            if (mytime >= t_stop) return
        else
            if (t_step == t_step_stop) return
        end if

        if (bubbles_lagrange) then
            call s_compute_EL_coupled_solver(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, stage=1)
            call s_update_lagrange_tdv_rk(stage=1)
        end if

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        !Evolve pb and mv for non-polytropic qbmm
        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    pb_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_pb(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    mv_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_mv(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (bodyForces) call s_apply_bodyforces(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, dt)

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(1)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        call nvtxEndRange

    end subroutine s_1st_order_tvd_rk

    !> 2nd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_2nd_order_tvd_rk(t_step, time_avg)

        integer, intent(in) :: t_step
        real(wp), intent(inout) :: time_avg

        integer :: i, j, k, l, q!< Generic loop iterator
        real(wp) :: start, finish

        ! Stage 1 of 2

        call cpu_time(start)

        call nvtxStartRange("TIMESTEP")

        call s_compute_rhs(q_cons_ts(1)%vf, q_T_sf, q_prim_vf, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, &
        rhs_rhouu, du_dxyz)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (cfl_dt) then
            if (mytime >= t_stop) return
        else
            if (t_step == t_step_stop) return
        end if

        if (bubbles_lagrange) then
            call s_compute_EL_coupled_solver(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, stage=1)
            call s_update_lagrange_tdv_rk(stage=1)
        end if

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        !Evolve pb and mv for non-polytropic qbmm
        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    pb_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_pb(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    mv_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_mv(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (bodyForces) call s_apply_bodyforces(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, dt)

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(2)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf, pb_ts(2)%sf, mv_ts(2)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf)
            end if
        end if

        ! Stage 2 of 2

        call s_compute_rhs(q_cons_ts(2)%vf, q_T_sf, q_prim_vf, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, t_step, time_avg, &
        rhs_rhouu, du_dxyz)

        if (bubbles_lagrange) then
            call s_compute_EL_coupled_solver(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, stage=2)
            call s_update_lagrange_tdv_rk(stage=2)
        end if

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + dt*rhs_vf(i)%sf(j, k, l))/2._wp
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    (pb_ts(1)%sf(j, k, l, q, i) &
                                     + pb_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_pb(j, k, l, q, i))/2._wp
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    (mv_ts(1)%sf(j, k, l, q, i) &
                                     + mv_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_mv(j, k, l, q, i))/2._wp
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (bodyForces) call s_apply_bodyforces(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, 2._wp*dt/3._wp)

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(1)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        call nvtxEndRange

        call cpu_time(finish)

    end subroutine s_2nd_order_tvd_rk

    !> 3rd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_3rd_order_tvd_rk(t_step, time_avg)

        integer, intent(IN) :: t_step
        real(wp), intent(INOUT) :: time_avg

        integer :: i, j, k, l, q !< Generic loop iterator

        real(wp) :: start, finish

        ! Stage 1 of 3

        if (.not. adap_dt) then
            call cpu_time(start)
            call nvtxStartRange("TIMESTEP")
        end if

        if (periodic_forcing) then 
            call s_compute_phase_average(q_cons_ts(1)%vf, q_T_sf, q_prim_vf, q_bar, q_spatial_avg, q_spatial_avg_glb, t_step+1)
            call s_compute_periodic_forcing(q_cons_ts(1)%vf, q_T_sf, q_prim_vf, q_bar, q_periodic_force)
        end if

        call s_compute_rhs(q_cons_ts(1)%vf, q_T_sf, q_prim_vf, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, &
        rhs_rhouu, du_dxyz)
        
        if (compute_CD_vi) then
            call s_compute_dragforce_vi(rhs_rhouu, q_prim_vf)
        end if
        if (compute_CD_si) then
            call s_compute_dragforce_si(q_prim_vf, du_dxyz)
        end if

        if (periodic_forcing) then
            call s_add_periodic_forcing(rhs_vf, q_periodic_force)
        end if

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (cfl_dt) then
            if (mytime >= t_stop) return
        else
            if (t_step == t_step_stop) return
        end if

        if (bubbles_lagrange) then
            call s_compute_EL_coupled_solver(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, stage=1)
            call s_update_lagrange_tdv_rk(stage=1)
        end if

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        !Evolve pb and mv for non-polytropic qbmm
        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    pb_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_pb(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    mv_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_mv(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (bodyForces) call s_apply_bodyforces(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, dt)

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(2)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf, pb_ts(2)%sf, mv_ts(2)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf)
            end if
        end if

        ! Stage 2 of 3

        call s_compute_rhs(q_cons_ts(2)%vf, q_T_sf, q_prim_vf, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, t_step, time_avg, &
        rhs_rhouu, du_dxyz)

        if (periodic_forcing) then
            call s_add_periodic_forcing(rhs_vf, q_periodic_force)
        end if

        if (bubbles_lagrange) then
            call s_compute_EL_coupled_solver(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, stage=2)
            call s_update_lagrange_tdv_rk(stage=2)
        end if

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            (3._wp*q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + dt*rhs_vf(i)%sf(j, k, l))/4._wp
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    (3._wp*pb_ts(1)%sf(j, k, l, q, i) &
                                     + pb_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_pb(j, k, l, q, i))/4._wp
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    (3._wp*mv_ts(1)%sf(j, k, l, q, i) &
                                     + mv_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_mv(j, k, l, q, i))/4._wp
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (bodyForces) call s_apply_bodyforces(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, dt/4._wp)

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(2)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf, pb_ts(2)%sf, mv_ts(2)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf)
            end if
        end if

        ! Stage 3 of 3
        call s_compute_rhs(q_cons_ts(2)%vf, q_T_sf, q_prim_vf, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, t_step, time_avg, &
        rhs_rhouu, du_dxyz)

        if (periodic_forcing) then
            call s_add_periodic_forcing(rhs_vf, q_periodic_force)
        end if

        if (bubbles_lagrange) then
            call s_compute_EL_coupled_solver(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, stage=3)
            call s_update_lagrange_tdv_rk(stage=3)
        end if

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + 2._wp*q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + 2._wp*dt*rhs_vf(i)%sf(j, k, l))/3._wp
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    (pb_ts(1)%sf(j, k, l, q, i) &
                                     + 2._wp*pb_ts(2)%sf(j, k, l, q, i) &
                                     + 2._wp*dt*rhs_pb(j, k, l, q, i))/3._wp
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    (mv_ts(1)%sf(j, k, l, q, i) &
                                     + 2._wp*mv_ts(2)%sf(j, k, l, q, i) &
                                     + 2._wp*dt*rhs_mv(j, k, l, q, i))/3._wp
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (bodyForces) call s_apply_bodyforces(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, 2._wp*dt/3._wp)

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3 .and. (.not. relax)) then
            call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
        end if

        call nvtxStartRange("RHS-ELASTIC")
        if (hyperelasticity) call s_hyperelastic_rmt_stress_update(q_cons_ts(1)%vf, q_prim_vf)
        call nvtxEndRange

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(1)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        if (.not. adap_dt) then
            call nvtxEndRange
            call cpu_time(finish)

            time = time + (finish - start)
        end if
    end subroutine s_3rd_order_tvd_rk

    ! computes the drag force for every immersed boundary using the volume integral formulation
    subroutine s_compute_dragforce_vi(rhs_rhouu, q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(in) :: rhs_rhouu
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp), dimension(1:num_ibs) :: F_D_global
        real(wp) :: C_D
        integer :: i, j, k, pos

        ! initialize F_D to zero
        !$acc parallel loop gang vector default(present)
        do i = 1, num_ibs
            F_D_vi(i) = 0._wp
        end do

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n
                do k = 0, p 
                    if (ib_markers%sf(i, j, k) /= 0) then 
                        
                        !$acc atomic 
                        F_D_vi(ib_markers%sf(i, j, k)) = F_D_vi(ib_markers%sf(i, j, k)) + rhs_rhouu(2)%sf(i, j, k) * dx(i) * dy(j) * dz(k)
                    
                    end if
                end do 
            end do 
        end do

        !$acc update host(F_D_vi)

        ! reduce drag force to one value
        do i = 1, num_ibs
            call s_mpi_allreduce_sum(F_D_vi(i), F_D_global(i))
        end do

        ! calculate C_D, C_D = F_D/(1/2 * rho * Uinf^2 * A)
        C_D = F_D_global(1) / (0.5 * rho_inf * (u_inf**2.0) * pi * (patch_ib(1)%radius**2.0))

        !print *, 'C_D (vi): ', C_D
        write(100) F_D_global  

    end subroutine s_compute_dragforce_vi

    ! computes the drag force for every immersed boundary using the surface integral formulation
    subroutine s_compute_dragforce_si(q_prim_vf, du_dxyz)
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(vector_field), dimension(1:3), intent(in) :: du_dxyz
        real(wp) :: mu_visc, C_D, divergence_u

        real(wp), dimension(1:num_ibs) :: F_D_global
    
        integer :: i, j, k, l, q, i_ibs, i_sphere_markers, pos

        mu_visc = 0.03334881107326017 ! M=1.2, Re=1500

        do i = 1, num_ibs
            F_D_si(i) = 0._wp
        end do

        do i = 1, 3
            do j = 1, 3
                !$acc update host(du_dxyz(i)%vf(j)%sf)
            end do
        end do

        !$acc update host(q_prim_vf(5)%sf)

        do i_ibs = 1, num_ibs
            do i_sphere_markers = 1, num_sphere_markers(i_ibs)

                i = sphere_markers_loc(i_ibs)%sf(1, i_sphere_markers, 1)
                j = sphere_markers_loc(i_ibs)%sf(1, i_sphere_markers, 2)
                k = sphere_markers_loc(i_ibs)%sf(1, i_sphere_markers, 3)

                ! x-dir ============================================================
                if (abs(levelset_norm%sf(i, j, k, i_ibs, 1)) >= abs(levelset_norm%sf(i, j, k, i_ibs, 2)) .and. abs(levelset_norm%sf(i, j, k, i_ibs, 1)) >= abs(levelset_norm%sf(i, j, k, i_ibs, 3))) then 
                    if (x_cc(i) < patch_ib(i_ibs)%x_centroid) then
                        x_surf = patch_ib(i_ibs)%x_centroid & 
                            - sqrt(patch_ib(i_ibs)%radius**2.0 & 
                            - (y_cc(j) - patch_ib(i_ibs)%y_centroid)**2.0 & 
                            - (z_cc(k) - patch_ib(i_ibs)%z_centroid)**2.0)
                        
                        do l = 1, 3
                            dudx_surf(l) = du_dxyz(l)%vf(1)%sf(i-1, j, k) + (x_surf - x_cc(i-1))/(x_cc(i) - x_cc(i-1))*(du_dxyz(l)%vf(1)%sf(i, j, k) - du_dxyz(l)%vf(1)%sf(i-1, j, k))
                            dudy_surf(l) = du_dxyz(l)%vf(2)%sf(i-1, j, k) + (x_surf - x_cc(i-1))/(x_cc(i) - x_cc(i-1))*(du_dxyz(l)%vf(2)%sf(i, j, k) - du_dxyz(l)%vf(2)%sf(i-1, j, k))
                            dudz_surf(l) = du_dxyz(l)%vf(3)%sf(i-1, j, k) + (x_surf - x_cc(i-1))/(x_cc(i) - x_cc(i-1))*(du_dxyz(l)%vf(3)%sf(i, j, k) - du_dxyz(l)%vf(3)%sf(i-1, j, k))

                        end do
                        pressure_surf = q_prim_vf(5)%sf(i-1, j, k) + (x_surf - x_cc(i-1))/(x_cc(i) - x_cc(i-1))*(q_prim_vf(5)%sf(i, j, k) - q_prim_vf(5)%sf(i-1, j, k))

                    else if (x_cc(i) > patch_ib(i_ibs)%x_centroid) then 
                        x_surf = patch_ib(i_ibs)%x_centroid & 
                            + sqrt(patch_ib(i_ibs)%radius**2.0 & 
                            - (y_cc(j) - patch_ib(i_ibs)%y_centroid)**2.0 & 
                            - (z_cc(k) - patch_ib(i_ibs)%z_centroid)**2.0)

                        do l = 1, 3
                            dudx_surf(l) = du_dxyz(l)%vf(1)%sf(i+1, j, k) + (x_surf - x_cc(i+1))/(x_cc(i) - x_cc(i+1))*(du_dxyz(l)%vf(1)%sf(i, j, k) - du_dxyz(l)%vf(1)%sf(i+1, j, k))
                            dudy_surf(l) = du_dxyz(l)%vf(2)%sf(i+1, j, k) + (x_surf - x_cc(i+1))/(x_cc(i) - x_cc(i+1))*(du_dxyz(l)%vf(2)%sf(i, j, k) - du_dxyz(l)%vf(2)%sf(i+1, j, k))
                            dudz_surf(l) = du_dxyz(l)%vf(3)%sf(i+1, j, k) + (x_surf - x_cc(i+1))/(x_cc(i) - x_cc(i+1))*(du_dxyz(l)%vf(3)%sf(i, j, k) - du_dxyz(l)%vf(3)%sf(i+1, j, k))

                        end do
                        pressure_surf = q_prim_vf(5)%sf(i+1, j, k) + (x_surf - x_cc(i+1))/(x_cc(i) - x_cc(i+1))*(q_prim_vf(5)%sf(i, j, k) - q_prim_vf(5)%sf(i+1, j, k))
                    
                    end if

                ! y-dir ============================================================
                else if (abs(levelset_norm%sf(i, j, k, i_ibs, 2)) >= abs(levelset_norm%sf(i, j, k, i_ibs, 1)) .and. abs(levelset_norm%sf(i, j, k, i_ibs, 2)) >= abs(levelset_norm%sf(i, j, k, i_ibs, 3))) then 
                    if (y_cc(j) < patch_ib(i_ibs)%y_centroid) then 
                        y_surf = patch_ib(i_ibs)%y_centroid & 
                            - sqrt(patch_ib(i_ibs)%radius**2.0 & 
                            - (x_cc(i) - patch_ib(i_ibs)%x_centroid)**2.0 & 
                            - (z_cc(k) - patch_ib(i_ibs)%z_centroid)**2.0)

                        do l = 1, 3
                            dudx_surf(l) = du_dxyz(l)%vf(1)%sf(i, j-1, k) + (y_surf - y_cc(j-1))/(y_cc(j) - y_cc(j-1))*(du_dxyz(l)%vf(1)%sf(i, j, k) - du_dxyz(l)%vf(1)%sf(i, j-1, k))
                            dudy_surf(l) = du_dxyz(l)%vf(2)%sf(i, j-1, k) + (y_surf - y_cc(j-1))/(y_cc(j) - y_cc(j-1))*(du_dxyz(l)%vf(2)%sf(i, j, k) - du_dxyz(l)%vf(2)%sf(i, j-1, k))
                            dudz_surf(l) = du_dxyz(l)%vf(3)%sf(i, j-1, k) + (y_surf - y_cc(j-1))/(y_cc(j) - y_cc(j-1))*(du_dxyz(l)%vf(3)%sf(i, j, k) - du_dxyz(l)%vf(3)%sf(i, j-1, k))

                        end do
                        pressure_surf = q_prim_vf(5)%sf(i, j-1, k) + (y_surf - y_cc(j-1))/(y_cc(j) - y_cc(j-1))*(q_prim_vf(5)%sf(i, j, k) - q_prim_vf(5)%sf(i, j-1, k))

                    else if (y_cc(j) > patch_ib(i_ibs)%y_centroid) then 
                        y_surf = patch_ib(i_ibs)%y_centroid & 
                            + sqrt(patch_ib(i_ibs)%radius**2.0 & 
                            - (x_cc(i) - patch_ib(i_ibs)%x_centroid)**2.0 & 
                            - (z_cc(k) - patch_ib(i_ibs)%z_centroid)**2.0)

                        do l = 1, 3
                            dudx_surf(l) = du_dxyz(l)%vf(1)%sf(i, j+1, k) + (y_surf - y_cc(j+1))/(y_cc(j) - y_cc(j+1))*(du_dxyz(l)%vf(1)%sf(i, j, k) - du_dxyz(l)%vf(1)%sf(i, j+1, k))
                            dudy_surf(l) = du_dxyz(l)%vf(2)%sf(i, j+1, k) + (y_surf - y_cc(j+1))/(y_cc(j) - y_cc(j+1))*(du_dxyz(l)%vf(2)%sf(i, j, k) - du_dxyz(l)%vf(2)%sf(i, j+1, k))
                            dudz_surf(l) = du_dxyz(l)%vf(3)%sf(i, j+1, k) + (y_surf - y_cc(j+1))/(y_cc(j) - y_cc(j+1))*(du_dxyz(l)%vf(3)%sf(i, j, k) - du_dxyz(l)%vf(3)%sf(i, j+1, k))

                        end do
                        pressure_surf = q_prim_vf(5)%sf(i, j+1, k) + (y_surf - y_cc(j+1))/(y_cc(j) - y_cc(j+1))*(q_prim_vf(5)%sf(i, j, k) - q_prim_vf(5)%sf(i, j+1, k))
                    
                    end if

                ! z-dir ============================================================
                else if (abs(levelset_norm%sf(i, j, k, i_ibs, 3)) >= abs(levelset_norm%sf(i, j, k, i_ibs, 1)) .and. abs(levelset_norm%sf(i, j, k, i_ibs, 3)) >= abs(levelset_norm%sf(i, j, k, i_ibs, 2))) then 
                    if (z_cc(k) < patch_ib(i_ibs)%z_centroid) then
                        z_surf = patch_ib(i_ibs)%z_centroid & 
                            - sqrt(patch_ib(i_ibs)%radius**2.0 & 
                            - (x_cc(i) - patch_ib(i_ibs)%x_centroid)**2.0 & 
                            - (y_cc(j) - patch_ib(i_ibs)%y_centroid)**2.0)

                        do l = 1, 3
                            dudx_surf(l) = du_dxyz(l)%vf(1)%sf(i, j, k-1) + (z_surf - z_cc(k-1))/(z_cc(k) - z_cc(k-1))*(du_dxyz(l)%vf(1)%sf(i, j, k) - du_dxyz(l)%vf(1)%sf(i, j, k-1))
                            dudy_surf(l) = du_dxyz(l)%vf(2)%sf(i, j, k-1) + (z_surf - z_cc(k-1))/(z_cc(k) - z_cc(k-1))*(du_dxyz(l)%vf(2)%sf(i, j, k) - du_dxyz(l)%vf(2)%sf(i, j, k-1))
                            dudz_surf(l) = du_dxyz(l)%vf(3)%sf(i, j, k-1) + (z_surf - z_cc(k-1))/(z_cc(k) - z_cc(k-1))*(du_dxyz(l)%vf(3)%sf(i, j, k) - du_dxyz(l)%vf(3)%sf(i, j, k-1))

                        end do
                        pressure_surf = q_prim_vf(5)%sf(i, j, k-1) + (z_surf - z_cc(k-1))/(z_cc(k) - z_cc(k-1))*(q_prim_vf(5)%sf(i, j, k) - q_prim_vf(5)%sf(i, j, k-1))

                    else if (z_cc(k) > patch_ib(i_ibs)%z_centroid) then 
                        z_surf = patch_ib(i_ibs)%z_centroid & 
                            + sqrt(patch_ib(i_ibs)%radius**2.0 & 
                            - (x_cc(i) - patch_ib(i_ibs)%x_centroid)**2.0 & 
                            - (y_cc(j) - patch_ib(i_ibs)%y_centroid)**2.0)

                        do l = 1, 3
                            dudx_surf(l) = du_dxyz(l)%vf(1)%sf(i, j, k+1) + (z_surf - z_cc(k+1))/(z_cc(k) - z_cc(k+1))*(du_dxyz(l)%vf(1)%sf(i, j, k) - du_dxyz(l)%vf(1)%sf(i, j, k+1))
                            dudy_surf(l) = du_dxyz(l)%vf(2)%sf(i, j, k+1) + (z_surf - z_cc(k+1))/(z_cc(k) - z_cc(k+1))*(du_dxyz(l)%vf(2)%sf(i, j, k) - du_dxyz(l)%vf(2)%sf(i, j, k+1))
                            dudz_surf(l) = du_dxyz(l)%vf(3)%sf(i, j, k+1) + (z_surf - z_cc(k+1))/(z_cc(k) - z_cc(k+1))*(du_dxyz(l)%vf(3)%sf(i, j, k) - du_dxyz(l)%vf(3)%sf(i, j, k+1))

                        end do
                        pressure_surf = q_prim_vf(5)%sf(i, j, k+1) + (z_surf - z_cc(k+1))/(z_cc(k) - z_cc(k+1))*(q_prim_vf(5)%sf(i, j, k) - q_prim_vf(5)%sf(i, j, k+1))

                    end if

                else 
                    print *, 'FAILURE TO INTERPOLATE: SURFACE NORM'
                end if ! end interpolation to surface

                ! drag force calculation ============================================================
                do l = 1, 3
                    F_D_mat(l) = 0._wp
                end do

                do l = 1, 3
                    stress_tensor(l, 1) = dudx_surf(l)
                    stress_tensor(l, 2) = dudy_surf(l)
                    stress_tensor(l, 3) = dudz_surf(l)
                end do 

                do l = 1, 3
                    stress_tensor(1, l) = stress_tensor(1, l) + dudx_surf(l)
                    stress_tensor(2, l) = stress_tensor(2, l) + dudy_surf(l)
                    stress_tensor(3, l) = stress_tensor(3, l) + dudz_surf(l)
                end do

                divergence_u = dudx_surf(1) + dudy_surf(2) + dudz_surf(3) ! du1dx + du2dy + du3dz

                do l = 1, 3
                    stress_tensor(l, l) = stress_tensor(l, l) - (2._wp/3._wp * divergence_u)
                end do

                stress_tensor = -mu_visc * stress_tensor

                do l = 1, 3
                    stress_tensor(l, l) = stress_tensor(l,l) + pressure_surf
                end do

                do l = 1, 3
                    do q = 1, 3
                        F_D_mat(l) = F_D_mat(l) + stress_tensor(l, q) * levelset_norm%sf(i, j, k, i_ibs, q)
                    end do
                end do

                F_D_mat = F_D_mat * data_plane_area(i_ibs)%sf(1, 1, i_sphere_markers)

                F_D_si(i_ibs) = F_D_si(i_ibs) + F_D_mat(1)

            end do ! marker loop
        end do ! ib loop
        
        do i = 1, num_ibs
            call s_mpi_allreduce_sum(F_D_si(i), F_D_global(i))
        end do

        C_D = -F_D_global(4) / (0.5 * rho_inf * (u_inf**2.0) * pi * (patch_ib(1)%radius**2.0))

        print *, 'C_D (si): ', C_D

        write(101) F_D_global

    end subroutine s_compute_dragforce_si 

    ! computes the phase average (time and space) of rho*u, rho, and T
    subroutine s_compute_phase_average(q_cons_vf, q_T_sf, q_prim_vf, q_bar, q_spatial_avg, q_spatial_avg_glb, t_step)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), intent(inout) :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), dimension(1:5), intent(inout) :: q_bar ! 1:3 rho*u, 4 rho, 5 T
        real(wp), dimension(1:5), intent(inout) :: q_spatial_avg ! 1:3 rho*u, 4 rho, 5 T
        real(wp), dimension(1:5), intent(inout) :: q_spatial_avg_glb
        integer :: i, j, k, t_step, pos

        ! initialize
        !$acc parallel loop gang vector default(present)
        do i = 1, 5
            q_spatial_avg(i) = 0._wp
        end do

        ! spatial average
        !$acc parallel loop collapse(3) gang vector default(present) reduction(+:q_spatial_avg(1), q_spatial_avg(2), q_spatial_avg(3), q_spatial_avg(4), q_spatial_avg(5))
        do i = 0, m 
            do j = 0, n 
                do k = 0, p 
                    if (ib_markers%sf(i, j, k) == 0) then
                        q_spatial_avg(4) = q_spatial_avg(4) + q_cons_vf(1)%sf(i, j, k)
                        q_spatial_avg(5) = q_spatial_avg(5) + 1.4 * (q_cons_vf(5)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k) & 
                                           - 0.5 * ((q_cons_vf(2)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k))**2 & 
                                           + (q_cons_vf(3)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k))**2 & 
                                           + (q_cons_vf(4)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k))**2)) 
                                           

                        q_spatial_avg(1) = q_spatial_avg(1) + q_cons_vf(2)%sf(i, j, k)
                        q_spatial_avg(2) = q_spatial_avg(2) + q_cons_vf(3)%sf(i, j, k)
                        q_spatial_avg(3) = q_spatial_avg(3) + q_cons_vf(4)%sf(i, j, k)
                    end if
                end do
            end do
        end do

        !$acc update host(q_spatial_avg(1), q_spatial_avg(2), q_spatial_avg(3), q_spatial_avg(4), q_spatial_avg(5))

        call s_mpi_allreduce_sum(q_spatial_avg(4), q_spatial_avg_glb(4))
        call s_mpi_allreduce_sum(q_spatial_avg(5), q_spatial_avg_glb(5))
        call s_mpi_allreduce_sum(q_spatial_avg(1), q_spatial_avg_glb(1))
        call s_mpi_allreduce_sum(q_spatial_avg(2), q_spatial_avg_glb(2))
        call s_mpi_allreduce_sum(q_spatial_avg(3), q_spatial_avg_glb(3))

        q_spatial_avg_glb(4) = q_spatial_avg_glb(4) / N_x_total_glb
        q_spatial_avg_glb(5) = q_spatial_avg_glb(5) / N_x_total_glb
        q_spatial_avg_glb(1) = q_spatial_avg_glb(1) / N_x_total_glb
        q_spatial_avg_glb(2) = q_spatial_avg_glb(2) / N_x_total_glb
        q_spatial_avg_glb(3) = q_spatial_avg_glb(3) / N_x_total_glb

        ! write the spatial avg of the x-mom 
        write(102) q_spatial_avg_glb(1) 
        print *, 'T', q_spatial_avg_glb(5)

        ! set reference quantities
        if ((t_step - 1) == 0) then
            u_inf = (q_spatial_avg_glb(1) / q_spatial_avg_glb(4)) / (1._wp - volfrac_phi)
            rho_inf = q_spatial_avg_glb(4) / (1._wp - volfrac_phi)
            T_inf = q_spatial_avg_glb(5) / (1._wp - volfrac_phi)
            !$acc update device(u_inf, rho_inf, T_inf)
        end if
        
        !$acc update device(q_spatial_avg_glb(1), q_spatial_avg_glb(2), q_spatial_avg_glb(3), q_spatial_avg_glb(4), q_spatial_avg_glb(5))

        ! time average
        !$acc parallel loop collapse(3) gang vector default(present) copyin(t_step)
        do i = 0, m 
            do j = 0, n
                do k = 0, p 
                    ! rho*u bar
                    q_bar(1)%sf(i, j, k) = ( (q_spatial_avg_glb(1) + (t_step - 1._wp)*q_bar(1)%sf(i, j, k)) / t_step ) 
                    q_bar(2)%sf(i, j, k) = ( (q_spatial_avg_glb(2) + (t_step - 1._wp)*q_bar(2)%sf(i, j, k)) / t_step ) 
                    q_bar(3)%sf(i, j, k) = ( (q_spatial_avg_glb(3) + (t_step - 1._wp)*q_bar(3)%sf(i, j, k)) / t_step ) 
                    
                    ! rho bar
                    q_bar(4)%sf(i, j, k) = ( (q_spatial_avg_glb(4) + (t_step - 1._wp)*q_bar(4)%sf(i, j, k)) / t_step ) 
                    ! T bar
                    q_bar(5)%sf(i, j, k) = ( (q_spatial_avg_glb(5) + (t_step - 1._wp)*q_bar(5)%sf(i, j, k)) / t_step ) 
                end do 
            end do 
        end do

    end subroutine s_compute_phase_average

    ! computes the periodic forcing terms described in Khalloufi and Capecelatro
    subroutine s_compute_periodic_forcing(q_cons_vf, q_T_sf, q_prim_vf, q_bar, q_periodic_force)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), intent(inout) :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), dimension(1:5), intent(inout) :: q_bar
        type(scalar_field), dimension(1:8), intent(inout) :: q_periodic_force

        integer :: i, j, k, l

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    ! f_u
                    q_periodic_force(1)%sf(i, j, k) = (rho_inf*u_inf - q_bar(1)%sf(i, j, k)/(1._wp - volfrac_phi)) / dt
                    q_periodic_force(2)%sf(i, j, k) = (rho_inf*u_inf - q_bar(2)%sf(i, j, k)/(1._wp - volfrac_phi)) / dt
                    q_periodic_force(3)%sf(i, j, k) = (rho_inf*u_inf - q_bar(3)%sf(i, j, k)/(1._wp - volfrac_phi)) / dt

                    ! u*f_u
                    q_periodic_force(4)%sf(i, j, k) = q_cons_vf(2)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k) * q_periodic_force(1)%sf(i, j, k)
                    q_periodic_force(5)%sf(i, j, k) = q_cons_vf(3)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k) * q_periodic_force(2)%sf(i, j, k)
                    q_periodic_force(6)%sf(i, j, k) = q_cons_vf(4)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k) * q_periodic_force(3)%sf(i, j, k)

                    ! f_rho
                    q_periodic_force(7)%sf(i, j, k) = (rho_inf - q_bar(4)%sf(i, j, k)/(1._wp - volfrac_phi)) / dt

                    ! f_T
                    q_periodic_force(8)%sf(i, j, k) = (q_cons_vf(1)%sf(i, j, k) / 1.4) * (T_inf - q_bar(5)%sf(i, j, k)/(1._wp - volfrac_phi)) / dt
                end do 
            end do
        end do

    end subroutine s_compute_periodic_forcing

    ! adds periodic forcing terms to RHS, as detailed in Khalloufi and Capecelatro
    subroutine s_add_periodic_forcing(rhs_vf, q_periodic_force)
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), dimension(1:8), intent(inout) :: q_periodic_force

        integer :: i, j, k

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    !if (ib_markers%sf(i, j, k) == 0) then
                    rhs_vf(1)%sf(i, j, k) = rhs_vf(1)%sf(i, j, k) + q_periodic_force(7)%sf(i, j, k) ! continuity
                    rhs_vf(2)%sf(i, j, k) = rhs_vf(2)%sf(i, j, k) + q_periodic_force(1)%sf(i, j, k) ! x momentum
                    rhs_vf(5)%sf(i, j, k) = rhs_vf(5)%sf(i, j, k) + q_periodic_force(4)%sf(i, j, k) + q_periodic_force(8)%sf(i, j, k) ! energy
                    !end if 
                end do
            end do
        end do

    end subroutine s_add_periodic_forcing

    !> Strang splitting scheme with 3rd order TVD RK time-stepping algorithm for
        !!      the flux term and adaptive time stepping algorithm for
        !!      the source term
        !! @param t_step Current time-step
    subroutine s_strang_splitting(t_step, time_avg)

        integer, intent(in) :: t_step
        real(wp), intent(inout) :: time_avg

        real(wp) :: start, finish

        call cpu_time(start)

        call nvtxStartRange("TIMESTEP")

        ! Stage 1 of 3
        call s_adaptive_dt_bubble(t_step)

        ! Stage 2 of 3
        call s_3rd_order_tvd_rk(t_step, time_avg)

        ! Stage 3 of 3
        call s_adaptive_dt_bubble(t_step)

        call nvtxEndRange

        call cpu_time(finish)

        time = time + (finish - start)

    end subroutine s_strang_splitting

    !> Bubble source part in Strang operator splitting scheme
        !! @param t_step Current time-step
    subroutine s_adaptive_dt_bubble(t_step)

        integer, intent(in) :: t_step

        type(vector_field) :: gm_alpha_qp

        call s_convert_conservative_to_primitive_variables( &
            q_cons_ts(1)%vf, &
            q_T_sf, &
            q_prim_vf, &
            idwint, &
            gm_alpha_qp%vf)

        call s_compute_bubble_EE_source(q_cons_ts(1)%vf, q_prim_vf, t_step, rhs_vf)

        call s_comp_alpha_from_n(q_cons_ts(1)%vf)

    end subroutine s_adaptive_dt_bubble

    subroutine s_compute_dt()

        real(wp) :: rho        !< Cell-avg. density
        real(wp), dimension(num_dims) :: vel        !< Cell-avg. velocity
        real(wp) :: vel_sum    !< Cell-avg. velocity sum
        real(wp) :: pres       !< Cell-avg. pressure
        real(wp), dimension(num_fluids) :: alpha      !< Cell-avg. volume fraction
        real(wp) :: gamma      !< Cell-avg. sp. heat ratio
        real(wp) :: pi_inf     !< Cell-avg. liquid stiffness function
        real(wp) :: c          !< Cell-avg. sound speed
        real(wp) :: H          !< Cell-avg. enthalpy
        real(wp), dimension(2) :: Re         !< Cell-avg. Reynolds numbers
        type(vector_field) :: gm_alpha_qp

        real(wp) :: dt_local
        integer :: j, k, l !< Generic loop iterators

        call s_convert_conservative_to_primitive_variables( &
            q_cons_ts(1)%vf, &
            q_T_sf, &
            q_prim_vf, &
            idwint, &
            gm_alpha_qp%vf)

        !$acc parallel loop collapse(3) gang vector default(present) private(vel, alpha, Re)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    call s_compute_enthalpy(q_prim_vf, pres, rho, gamma, pi_inf, Re, H, alpha, vel, vel_sum, j, k, l)

                    ! Compute mixture sound speed
                    call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, alpha, vel_sum, 0._wp, c)

                    call s_compute_dt_from_cfl(vel, c, max_dt, rho, Re, j, k, l)
                end do
            end do
        end do

        !$acc kernels
        dt_local = minval(max_dt)
        !$acc end kernels

        if (num_procs == 1) then
            dt = dt_local
        else
            call s_mpi_allreduce_min(dt_local, dt)
        end if

        !$acc update device(dt)

    end subroutine s_compute_dt

    !> This subroutine applies the body forces source term at each
        !! Runge-Kutta stage
    subroutine s_apply_bodyforces(q_cons_vf, q_prim_vf, rhs_vf, ldt)

        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(1:sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(1:sys_size), intent(inout) :: rhs_vf

        real(wp), intent(in) :: ldt !< local dt

        integer :: i, j, k, l

        call nvtxStartRange("RHS-BODYFORCES")
        call s_compute_body_forces_rhs(q_prim_vf, q_cons_vf, rhs_vf)

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = momxb, E_idx
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l) + &
                                                   ldt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        call nvtxEndRange

    end subroutine s_apply_bodyforces

    !> This subroutine saves the temporary q_prim_vf vector
        !!      into the q_prim_ts vector that is then used in p_main
        !! @param t_step current time-step
    subroutine s_time_step_cycling(t_step)

        integer, intent(in) :: t_step

        integer :: i !< Generic loop iterator

        do i = 1, sys_size
            !$acc update host(q_prim_vf(i)%sf)
        end do

        if (t_step == t_step_start) then
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 1) then
            do i = 1, sys_size
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 2) then
            do i = 1, sys_size
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 3) then
            do i = 1, sys_size
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        else ! All other timesteps
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_ts(2)%vf(i)%sf(:, :, :)
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_ts(1)%vf(i)%sf(:, :, :)
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_ts(0)%vf(i)%sf(:, :, :)
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        end if

    end subroutine s_time_step_cycling

    !> (Adaptive) 4th/5th order Runge—Kutta–Cash–Karp (RKCK) time-stepping algorithm (Cash J. and Karp A., 1990)
        !!      Method for initial value problems with rapidly varying RHS. A maximum error between the 4th and 5th
        !!      order Runge-Kutta-Cash-Karp solutions for the same time step size is calculated. If the error is
        !!      smaller than a tolerance, then the algorithm employs the 5th order solution, while if not, both
        !!      eulerian/lagrangian variables are re-calculated with a smaller time step size.
        !! @param t_step Current time-step
        !! @param hdid Advanced time increment (adaptive time stepping)
    subroutine s_4th_5th_order_rkck(t_step, time_avg)

        integer, intent(in) :: t_step
        real(wp), intent(out) :: time_avg

        logical :: restart_rkck_step, start_rkck_step
        real(wp) :: lag_largestep, rkck_errmax, dt_did
        integer :: RKstep

        mytime = mytime - dt

        start_rkck_step = .true.
        restart_rkck_step = .false.

        do while (start_rkck_step .or. restart_rkck_step)

            start_rkck_step = .false.
            restart_rkck_step = .false.

            ! FIRST TIME-STAGE
            RKstep = 1
            rkck_time_tmp = mytime + rkck_c1*dt
!$acc update device (rkck_time_tmp)

#ifdef DEBUG
            if (proc_rank == 0) print *, 'RKCK 1st time-stage at', rkck_time_tmp
#endif
            call s_compute_rhs(q_cons_ts(1)%vf, q_T_sf, q_prim_vf, rhs_ts_rkck(1)%vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, &
            rhs_rhouu, du_dxyz)
            call s_compute_EL_coupled_solver(q_cons_ts(1)%vf, q_prim_vf, rhs_ts_rkck(1)%vf, RKstep)
            call s_update_tmp_rkck(RKstep, q_cons_ts, rhs_ts_rkck, lag_largestep)
            if (lag_largestep > 0._wp) call s_compute_rkck_dt(lag_largestep, restart_rkck_step)
            if (restart_rkck_step) cycle

            ! SECOND TIME-STAGE
            RKstep = 2
            rkck_time_tmp = mytime + rkck_c2*dt
!$acc update device (rkck_time_tmp)

#ifdef DEBUG
            if (proc_rank == 0) print *, 'RKCK 2nd time-stage at', rkck_time_tmp
#endif
            call s_compute_rhs(q_cons_ts(2)%vf, q_T_sf, q_prim_vf, rhs_ts_rkck(2)%vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, &
            rhs_rhouu, du_dxyz)
            call s_compute_EL_coupled_solver(q_cons_ts(2)%vf, q_prim_vf, rhs_ts_rkck(2)%vf, RKstep)
            call s_update_tmp_rkck(RKstep, q_cons_ts, rhs_ts_rkck, lag_largestep)
            if (lag_largestep > 0._wp) call s_compute_rkck_dt(lag_largestep, restart_rkck_step)
            if (restart_rkck_step) cycle

            ! THIRD TIME-STAGE
            RKstep = 3
            rkck_time_tmp = mytime + rkck_c3*dt
!$acc update device (rkck_time_tmp)

#ifdef DEBUG
            if (proc_rank == 0) print *, 'RKCK 3rd time-stage at', rkck_time_tmp
#endif
            call s_compute_rhs(q_cons_ts(2)%vf, q_T_sf, q_prim_vf, rhs_ts_rkck(3)%vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, &
            rhs_rhouu, du_dxyz)
            call s_compute_EL_coupled_solver(q_cons_ts(2)%vf, q_prim_vf, rhs_ts_rkck(3)%vf, RKstep)
            call s_update_tmp_rkck(RKstep, q_cons_ts, rhs_ts_rkck, lag_largestep)
            if (lag_largestep > 0._wp) call s_compute_rkck_dt(lag_largestep, restart_rkck_step)
            if (restart_rkck_step) cycle

            ! FOURTH TIME-STAGE
            RKstep = 4
            rkck_time_tmp = mytime + rkck_c4*dt
!$acc update device (rkck_time_tmp)

#ifdef DEBUG
            if (proc_rank == 0) print *, 'RKCK 4th time-stage at', rkck_time_tmp
#endif
            call s_compute_rhs(q_cons_ts(2)%vf, q_T_sf, q_prim_vf, rhs_ts_rkck(4)%vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, &
            rhs_rhouu, du_dxyz)
            call s_compute_EL_coupled_solver(q_cons_ts(2)%vf, q_prim_vf, rhs_ts_rkck(4)%vf, RKstep)
            call s_update_tmp_rkck(RKstep, q_cons_ts, rhs_ts_rkck, lag_largestep)
            if (lag_largestep > 0._wp) call s_compute_rkck_dt(lag_largestep, restart_rkck_step)
            if (restart_rkck_step) cycle

            ! FIFTH TIME-STAGE
            RKstep = 5
            rkck_time_tmp = mytime + rkck_c5*dt
!$acc update device (rkck_time_tmp)

#ifdef DEBUG
            if (proc_rank == 0) print *, 'RKCK 5th time-stage at', rkck_time_tmp
#endif
            call s_compute_rhs(q_cons_ts(2)%vf, q_T_sf, q_prim_vf, rhs_ts_rkck(5)%vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, &
            rhs_rhouu, du_dxyz)
            call s_compute_EL_coupled_solver(q_cons_ts(2)%vf, q_prim_vf, rhs_ts_rkck(5)%vf, 5)
            call s_update_tmp_rkck(5, q_cons_ts, rhs_ts_rkck, lag_largestep)
            if (lag_largestep > 0._wp) call s_compute_rkck_dt(lag_largestep, restart_rkck_step)
            if (restart_rkck_step) cycle

            ! SIXTH TIME-STAGE
            RKstep = 6
            rkck_time_tmp = mytime + rkck_c6*dt
!$acc update device (rkck_time_tmp)

#ifdef DEBUG
            if (proc_rank == 0) print *, 'RKCK 6th time-stage at', rkck_time_tmp
#endif
            call s_compute_rhs(q_cons_ts(2)%vf, q_T_sf, q_prim_vf, rhs_ts_rkck(6)%vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, &
            rhs_rhouu, du_dxyz)
            call s_compute_EL_coupled_solver(q_cons_ts(2)%vf, q_prim_vf, rhs_ts_rkck(6)%vf, 6)
            call s_update_tmp_rkck(6, q_cons_ts, rhs_ts_rkck, lag_largestep)
            if (lag_largestep > 0._wp) call s_compute_rkck_dt(lag_largestep, restart_rkck_step)
            if (restart_rkck_step) cycle

            dt_did = dt

            if (rkck_adap_dt) then
                ! TRUNCATION ERROR
#ifdef DEBUG
                if (proc_rank == 0) print *, 'Computing truncation error (4th/5th RKCK)'
#endif
                call s_calculate_rkck_truncation_error(rkck_errmax)
                call s_compute_rkck_dt(lag_largestep, restart_rkck_step, rkck_errmax)
                if (restart_rkck_step) cycle
            end if

        end do

        !> Update values
        mytime = mytime + dt_did
        call s_update_rkck(q_cons_ts)

        call s_write_void_evol(mytime)
        if (lag_params%write_bubbles_stats) call s_calculate_lag_bubble_stats()

        if (lag_params%write_bubbles) then
            !$acc update host(gas_p, gas_mv, intfc_rad, intfc_vel)
            call s_write_lag_particles(mytime)
        end if

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

    end subroutine s_4th_5th_order_rkck

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_time_steppers_module

        integer :: i, j !< Generic loop iterators

        ! Deallocating the cell-average conservative variables
        do i = 1, num_ts

            do j = 1, sys_size
                @:DEALLOCATE(q_cons_ts(i)%vf(j)%sf)
            end do

            @:DEALLOCATE(q_cons_ts(i)%vf)

        end do

        @:DEALLOCATE(q_cons_ts)

        ! Deallocating the cell-average primitive ts variables
        if (probe_wrt) then
            do i = 0, 3
                do j = 1, sys_size
                    @:DEALLOCATE(q_prim_ts(i)%vf(j)%sf)
                end do
                @:DEALLOCATE(q_prim_ts(i)%vf)
            end do
            @:DEALLOCATE(q_prim_ts)
        end if

        ! Deallocating the cell-average primitive variables
        do i = 1, adv_idx%end
            @:DEALLOCATE(q_prim_vf(i)%sf)
        end do

        if (elasticity) then
            do i = stress_idx%beg, stress_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        if (hyperelasticity) then
            do i = xibeg, xiend + 1
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        if (bubbles_euler) then
            do i = bub_idx%beg, bub_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        @:DEALLOCATE(q_prim_vf)

        ! Deallocating the cell-average RHS variable for adaptive method, Lagrangian solver
        if (bubbles_lagrange .and. time_stepper == 4) then ! RKCK stepper
            do i = 1, num_ts_rkck
                do j = 1, sys_size
                    @:DEALLOCATE(rhs_ts_rkck(i)%vf(j)%sf)
                end do
                @:DEALLOCATE(rhs_ts_rkck(i)%vf)
            end do
            @:DEALLOCATE(rhs_ts_rkck)
        else
            ! Deallocating the cell-average RHS variables
            do i = 1, sys_size
                @:DEALLOCATE(rhs_vf(i)%sf)
            end do

            @:DEALLOCATE(rhs_vf)
        end if

        do i = 1, sys_size
            @:DEALLOCATE(rhs_rhouu(i)%sf)
        end do 
        @:DEALLOCATE(rhs_rhouu)

        do i = 1, 3
            do j = 1, 3
                @:DEALLOCATE(du_dxyz(i)%vf(j)%sf)
            end do
            @:DEALLOCATE(du_dxyz(i)%vf)
        end do
        @:DEALLOCATE(du_dxyz)

        do i = 1, 5
            @:DEALLOCATE(q_bar(i)%sf)
        end do
        @:DEALLOCATE(q_bar)

        @:DEALLOCATE(F_D_vi)
        @:DEALLOCATE(F_D_si)

        @:DEALLOCATE(dudx_surf)
        @:DEALLOCATE(dudy_surf)
        @:DEALLOCATE(dudz_surf)

        @:DEALLOCATE(stress_tensor)
        @:DEALLOCATE(F_D_mat)
        
        do i = 1, 8
            @:DEALLOCATE(q_periodic_force(i)%sf)
        end do
        @:DEALLOCATE(q_periodic_force)

        @:DEALLOCATE(q_spatial_avg)
        @:DEALLOCATE(q_spatial_avg_glb)

        close(100)
        close(101)
        close(102)

        ! Writing the footer of and closing the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_close_run_time_information_file()
        end if

    end subroutine s_finalize_time_steppers_module

end module m_time_steppers
