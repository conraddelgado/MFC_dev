!>
!! @file m_fftw.f90
!! @brief Contains module m_fftw

#:include 'macros.fpp'

!> @brief The module contains the subroutines for the FFT routines
module m_fftw

    use, intrinsic :: iso_c_binding

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_ibm

#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif

#if defined(MFC_OpenACC) && defined(__PGI)
    use cufft
#elif defined(MFC_OpenACC)
    use hipfort
    use hipfort_check
    use hipfort_hipfft
#endif

    implicit none

    private; public :: s_initialize_fftw_module, &
 s_apply_fourier_filter, &
 s_finalize_fftw_module, & 
 s_initialize_fftw_explicit_filter_module, &
 s_apply_fftw_filter_cons, & 
 s_initialize_filtering_kernel, s_initialize_filtered_fluid_indicator_function, & 
 s_finalize_fftw_explicit_filter_module, & 
 s_apply_fftw_filter_tensor, s_apply_fftw_filter_scalarfield, & 
 s_transpose_z2y_mpi, s_transpose_y2x_mpi, s_transpose_x2y_mpi, s_transpose_y2z_mpi, & 
 s_mpi_perform_transpose_forward, s_mpi_perform_transpose_backward

#if !defined(MFC_OpenACC)
    include 'fftw3.f03'
#endif

    type(c_ptr) :: fwd_plan, bwd_plan
    type(c_ptr) :: fftw_real_data, fftw_cmplx_data, fftw_fltr_cmplx_data
    integer :: real_size, cmplx_size, x_size, batch_size, Nfq

    real(c_double), pointer :: data_real(:) !< Real data

    complex(c_double_complex), pointer :: data_cmplx(:) !<
    !! Complex data in Fourier space

    complex(c_double_complex), pointer :: data_fltr_cmplx(:) !<
    !! Filtered complex data in Fourier space

#if defined(MFC_OpenACC)
    !$acc declare create(real_size, cmplx_size, x_size, batch_size, Nfq)

    real(dp), allocatable, target :: data_real_gpu(:)
    complex(dp), allocatable, target :: data_cmplx_gpu(:)
    complex(dp), allocatable, target :: data_fltr_cmplx_gpu(:)
!$acc declare create(data_real_gpu, data_cmplx_gpu, data_fltr_cmplx_gpu)

#if defined(__PGI)
    integer :: fwd_plan_gpu, bwd_plan_gpu
#else
    type(c_ptr) :: fwd_plan_gpu, bwd_plan_gpu
#endif

    integer, allocatable :: gpu_fft_size(:), iembed(:), oembed(:)

    integer :: istride, ostride, idist, odist, rank
#endif

    integer :: ierr

    ! new for explicit filtering of data
    type(scalar_field), public :: fluid_indicator_function_I
    !$acc declare create(fluid_indicator_function_I)

#if defined(MFC_OpenACC)
    real(dp), allocatable, target :: real_data_gpu(:, :, :)
    complex(dp), allocatable, target :: complex_data_gpu(:, :, :)

    complex(dp), allocatable, target :: complex_y_gpu(:, :, :)
    complex(dp), allocatable, target :: complex_x_gpu(:, :, :)

    real(dp), allocatable, target :: real_kernelG_gpu(:, :, :)
    complex(dp), allocatable, target :: complex_kernelG_gpu(:, :, :)

    complex(dp), allocatable, target :: complex_kernelG_gpu_mpi(:, :, :)
    complex(dp), allocatable, target :: complex_data_gpu_mpi(:, :, :)

    integer :: forward_plan_gpu
    integer :: backward_plan_gpu

    integer :: forward_plan_kernelG_gpu

    integer :: plan_y_gpu 
    integer :: plan_x_gpu

    !$acc declare create(real_data_gpu, complex_data_gpu, complex_y_gpu, complex_x_gpu, real_kernelG_gpu, complex_kernelG_gpu, complex_kernelG_gpu_mpi, complex_data_gpu_mpi)
#endif
    type(c_ptr) :: plan_forward
    type(c_ptr) :: plan_backward

    real(c_double), pointer :: array_in(:, :, :)
    complex(c_double_complex), pointer :: array_out(:, :, :)

    type(c_ptr) :: p_real
    type(c_ptr) :: p_complex

    ! for filtering kernel
    type(c_ptr) :: plan_kernelG_forward

    real(c_double), pointer :: array_kernelG_in(:, :, :)
    complex(c_double_complex), pointer :: array_kernelG_out(:, :, :)

    type(c_ptr) :: p_kernelG_real
    type(c_ptr) :: p_kernelG_complex

contains

    !>  The purpose of this subroutine is to create the fftw plan
        !!      that will be used in the forward and backward DFTs when
        !!      applying the Fourier filter in the azimuthal direction.
    subroutine s_initialize_fftw_module

        ! Size of input array going into DFT
        real_size = p + 1
        ! Size of output array coming out of DFT
        cmplx_size = (p + 1)/2 + 1

        x_size = m + 1

        batch_size = x_size*sys_size

#if defined(MFC_OpenACC)
        rank = 1; istride = 1; ostride = 1

        allocate (gpu_fft_size(1:rank), iembed(1:rank), oembed(1:rank))

        gpu_fft_size(1) = real_size; 
        iembed(1) = 0
        oembed(1) = 0
        !$acc enter data copyin(real_size, cmplx_size, x_size, sys_size, batch_size, Nfq)
        !$acc update device(real_size, cmplx_size, x_size, sys_size, batch_size)
#else
        ! Allocate input and output DFT data sizes
        fftw_real_data = fftw_alloc_real(int(real_size, c_size_t))
        fftw_cmplx_data = fftw_alloc_complex(int(cmplx_size, c_size_t))
        fftw_fltr_cmplx_data = fftw_alloc_complex(int(cmplx_size, c_size_t))
        ! Associate input and output data pointers with allocated memory
        call c_f_pointer(fftw_real_data, data_real, [real_size])
        call c_f_pointer(fftw_cmplx_data, data_cmplx, [cmplx_size])
        call c_f_pointer(fftw_fltr_cmplx_data, data_fltr_cmplx, [cmplx_size])

        ! Generate plans for forward and backward DFTs
        fwd_plan = fftw_plan_dft_r2c_1d(real_size, data_real, data_cmplx, FFTW_ESTIMATE)
        bwd_plan = fftw_plan_dft_c2r_1d(real_size, data_fltr_cmplx, data_real, FFTW_ESTIMATE)
#endif

#if defined(MFC_OpenACC)
        @:ALLOCATE(data_real_gpu(1:real_size*x_size*sys_size))
        @:ALLOCATE(data_cmplx_gpu(1:cmplx_size*x_size*sys_size))
        @:ALLOCATE(data_fltr_cmplx_gpu(1:cmplx_size*x_size*sys_size))

#if defined(__PGI)
        ierr = cufftPlanMany(fwd_plan_gpu, rank, gpu_fft_size, iembed, istride, real_size, oembed, ostride, cmplx_size, CUFFT_D2Z, batch_size)
        ierr = cufftPlanMany(bwd_plan_gpu, rank, gpu_fft_size, iembed, istride, cmplx_size, oembed, ostride, real_size, CUFFT_Z2D, batch_size)
#else
        ierr = hipfftPlanMany(fwd_plan_gpu, rank, gpu_fft_size, iembed, istride, real_size, oembed, ostride, cmplx_size, HIPFFT_D2Z, batch_size)
        ierr = hipfftPlanMany(bwd_plan_gpu, rank, gpu_fft_size, iembed, istride, cmplx_size, oembed, ostride, real_size, HIPFFT_Z2D, batch_size)
#endif

#endif

    end subroutine s_initialize_fftw_module

    ! create fftw plan to be used for explicit filtering of data -> volume filtering
    subroutine s_initialize_fftw_explicit_filter_module
        integer :: fft_rank
        integer :: nfft(1), inembed(1), onembed(1)
        integer :: batch, istride, idist, ostride, odist

        print *, 'FFTW SETUP...'

        @:ALLOCATE(fluid_indicator_function_I%sf(0:m, 0:n, 0:p))
        @:ACC_SETUP_SFs(fluid_indicator_function_I)

#if defined(MFC_OpenACC)
        ! gpu data setup
        @:ALLOCATE(real_data_gpu(m+1, n+1, p+1))
        @:ALLOCATE(complex_data_gpu(m+1, n+1, (p+1)/2+1))

        ! MPI transpose arrays
        @:ALLOCATE(complex_y_gpu(m+1, (p+1)/2+1, n+1))
        @:ALLOCATE(complex_x_gpu(n+1, (p+1)/2+1, m+1))

        ! gpu kernel arrays
        @:ALLOCATE(real_kernelG_gpu(m+1, n+1, p+1))
        @:ALLOCATE(complex_kernelG_gpu(m+1, n+1, (p+1)/2+1))

        ! MPI transposed complex kernel array
        @:ALLOCATE(complex_kernelG_gpu_mpi(n+1, (p+1)/2+1, m+1))

        ! MPI transposed complex data array
        @:ALLOCATE(complex_data_gpu_mpi(n+1, (p+1)/2+1, m+1))

        ! gpu plan creation
        ierr = cufftPlan3d(forward_plan_gpu, m+1, n+1, p+1, CUFFT_D2Z)
        ierr = cufftPlan3d(backward_plan_gpu, m+1, n+1, p+1, CUFFT_Z2D)

        ierr = cufftPlan3d(forward_plan_kernelG_gpu, m+1, n+1, p+1, CUFFT_D2Z)

        fft_rank = 1

        ! MPI transpose y plan
        nfft(1) = n+1
        inembed(1) = n+1
        onembed(1) = n+1

        istride = 1
        idist = n+1
        ostride = 1
        odist = n+1

        batch = (m+1) * ((p+1)/2+1)

        ierr = cufftPlanMany(plan_y_gpu, fft_rank, nfft, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, batch)

        ! MPI transpose x plan
        nfft(1) = m+1
        inembed(1) = m+1
        onembed(1) = m+1
        
        istride = 1
        idist = m+1
        ostride = 1
        odist = m+1
        
        batch = (n+1) * ((p+1)/2+1)
        
        ierr = cufftPlanMany(plan_x_gpu, fft_rank, nfft, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, batch)                     

#else
        ! CPU
        ! data setup 
        p_real = fftw_alloc_real(int((m+1)*(n+1)*(p+1), C_SIZE_T))
        p_complex = fftw_alloc_complex(int(((m+1)/2+1)*(n+1)*(p+1), C_SIZE_T))

        call c_f_pointer(p_real, array_in, [m+1, n+1, p+1])
        call c_f_pointer(p_complex, array_out, [(m+1)/2+1, n+1, p+1])

        plan_forward = fftw_plan_dft_r2c_3d(p+1, n+1, m+1, array_in, array_out, FFTW_MEASURE)
        plan_backward = fftw_plan_dft_c2r_3d(p+1, n+1, m+1, array_out, array_in, FFTW_MEASURE)

        ! kernel setup
        p_kernelG_real = fftw_alloc_real(int((m+1)*(n+1)*(p+1), C_SIZE_T))
        p_kernelG_complex = fftw_alloc_complex(int(((m+1)/2+1)*(n+1)*(p+1), C_SIZE_T))

        call c_f_pointer(p_kernelG_real, array_kernelG_in, [m+1, n+1, p+1])
        call c_f_pointer(p_kernelG_complex, array_kernelG_out, [(m+1)/2+1, n+1, p+1])

        plan_kernelG_forward = fftw_plan_dft_r2c_3d(p+1, n+1, m+1, array_kernelG_in, array_kernelG_out, FFTW_MEASURE)
#endif

    end subroutine s_initialize_fftw_explicit_filter_module

    !< initialize the gaussian filtering kernel in real space and then compute its DFT
    subroutine s_initialize_filtering_kernel
        real(dp) :: sigma_stddev
        real(dp) :: Lx, Ly, Lz
        real(dp) :: x_r, y_r, z_r  
        real(dp) :: r2
        real(dp) :: G_norm_int, G_norm_int_glb
        integer :: i, j, k

        ! gaussian filter
        sigma_stddev = 3.0_dp * 0.05_dp

        Lx = x_domain_end_glb - x_domain_beg_glb
        Ly = y_domain_end_glb - y_domain_beg_glb  
        Lz = z_domain_end_glb - z_domain_beg_glb    
        
        G_norm_int = 0._dp

        print *, 'x, y, z local index sizes', m, n, p, num_procs

#if defined(MFC_OpenACC)    
        !$acc parallel loop collapse(3) gang vector default(present) reduction(+:G_norm_int) copyin(Lx, Ly, Lz, sigma_stddev) private(x_r, y_r, z_r, r2)
        do i = 0, m 
            do j = 0, n 
                do k = 0, p 
                    x_r = min(abs(x_cc(i) - x_domain_beg_glb), Lx - abs(x_cc(i) - x_domain_beg_glb))
                    y_r = min(abs(y_cc(j) - y_domain_beg_glb), Ly - abs(y_cc(j) - y_domain_beg_glb))
                    z_r = min(abs(z_cc(k) - z_domain_beg_glb), Lz - abs(z_cc(k) - z_domain_beg_glb))

                    r2 = x_r**2 + y_r**2 + z_r**2

                    real_kernelG_gpu(i+1, j+1, k+1) = exp(-r2 / (2.0_dp*sigma_stddev**2))

                    G_norm_int = G_norm_int + real_kernelG_gpu(i+1, j+1, k+1)*dx(i)*dy(j)*dz(k)
                end do 
            end do
        end do

        call s_mpi_allreduce_sum(G_norm_int, G_norm_int_glb) 

        ! normalize the gaussian kernel
        !$acc parallel loop collapse(3) gang vector default(present) copyin(G_norm_int_glb)
        do i = 1, m+1 
            do j = 1, n+1 
                do k = 1, p+1
                    real_kernelG_gpu(i, j, k) = real_kernelG_gpu(i, j, k) / G_norm_int_glb
                end do 
            end do 
        end do

        ! perform fourier transform on gaussian kernel
        ierr = cufftExecD2Z(forward_plan_kernelG_gpu, real_kernelG_gpu, complex_kernelG_gpu)

        ! transpose DFT data, only necessary with MPI domain decomposition and more than one rank
        if (num_procs > 1) then
            call s_transpose_z2y_mpi(complex_kernelG_gpu) ! -> gives complex_y_gpu

            ierr = cufftExecZ2Z(plan_y_gpu, complex_y_gpu, complex_y_gpu, CUFFT_FORWARD) 

            call s_transpose_y2x_mpi ! gives complex_x_gpu

            ierr = cufftExecZ2Z(plan_x_gpu, complex_x_gpu, complex_kernelG_gpu_mpi, CUFFT_FORWARD)

            ! normalize filtering kernel after fourier transform
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, n+1
                do j = 1, (p+1)/2+1
                    do k = 1, m+1
                        complex_kernelG_gpu_mpi(i, j, k) = complex_kernelG_gpu_mpi(i, j, k) / (real(m_glb+1, dp)*real(n_glb+1, dp)*real(p_glb+1, dp))
                    end do 
                end do 
            end do
        
        else 
            ! normalize transformed kernel (divide by Nx*Ny*Nz)
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, m+1
                do j = 1, n+1
                    do k = 1, (p+1)/2+1
                        complex_kernelG_gpu(i, j, k) = complex_kernelG_gpu(i, j, k) / (real(m_glb+1, dp)*real(n_glb+1, dp)*real(p_glb+1, dp))
                    end do 
                end do 
            end do
        end if

#else
        ! CPU
        do i = 0, m 
            do j = 0, n 
                do k = 0, p 
                    x_r = min(abs(x_cc(i) - x_domain_beg_glb), Lx - abs(x_cc(i) - x_domain_beg_glb))
                    y_r = min(abs(y_cc(j) - y_domain_beg_glb), Ly - abs(y_cc(j) - y_domain_beg_glb))
                    z_r = min(abs(z_cc(k) - z_domain_beg_glb), Lz - abs(z_cc(k) - z_domain_beg_glb))

                    r2 = x_r**2 + y_r**2 + z_r**2

                    array_kernelG_in(i+1, j+1, k+1) = exp(-r2 / (2.0_dp*sigma_stddev**2))

                    G_norm_int = G_norm_int + array_kernelG_in(i+1, j+1, k+1)*dx(i)*dy(j)*dz(k)
                end do 
            end do
        end do

        call s_mpi_allreduce_sum(G_norm_int, G_norm_int_glb) 

        array_kernelG_in = array_kernelG_in / G_norm_int_glb ! normalize gaussian, integrate to unity over domain

        call fftw_execute_dft_r2c(plan_kernelG_forward, array_kernelG_in, array_kernelG_out)
        
        array_kernelG_out = array_kernelG_out / (real(m_glb+1, dp)*real(n_glb+1, dp)*real(p_glb+1, dp)) ! normalize DFT
#endif

    end subroutine s_initialize_filtering_kernel

    !< initialize the fluid indicator function and compute its filtered counterpart
    subroutine s_initialize_filtered_fluid_indicator_function(filtered_fluid_indicator_function)
        type(scalar_field) :: filtered_fluid_indicator_function

        integer :: i, j, k

        ! define fluid indicator function
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n 
                do k = 0, p
                    if (ib_markers%sf(i, j, k) == 0) then 
                        fluid_indicator_function_I%sf(i, j, k) = 1.0_dp
                    else 
                        fluid_indicator_function_I%sf(i, j, k) = 0.0_dp
                    end if
                end do
            end do
        end do

        ! filter fluid indicator function -> stored in q_cons_vf(advxb)
#if defined(MFC_OpenACC)
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m 
            do j = 0, n 
                do k = 0, p
                    real_data_gpu(i+1, j+1, k+1) = fluid_indicator_function_I%sf(i, j, k)
                end do 
            end do 
        end do 

        ierr = cufftExecD2Z(forward_plan_gpu, real_data_gpu, complex_data_gpu)

        if (num_procs > 1) then
            call s_mpi_perform_transpose_forward(complex_data_gpu) !< gives complex_data_gpu_mpi

            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, n+1
                do j = 1, (p+1)/2+1
                    do k = 1, m+1
                        complex_data_gpu_mpi(i, j, k) = complex_data_gpu_mpi(i, j, k) * complex_kernelG_gpu_mpi(i, j, k)
                    end do 
                end do 
            end do

            call s_mpi_perform_transpose_backward(complex_data_gpu)

        else ! 1 rank
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, m+1
                do j = 1, n+1 
                    do k = 1, (p+1)/2+1
                        complex_data_gpu(i, j, k) = complex_data_gpu(i, j, k) * complex_kernelG_gpu(i, j, k)
                    end do 
                end do 
            end do
        end if

        ierr = cufftExecZ2D(backward_plan_gpu, complex_data_gpu, real_data_gpu)

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m 
            do j = 0, n 
                do k = 0, p 
                    filtered_fluid_indicator_function%sf(i, j, k) = real_data_gpu(i+1, j+1, k+1) / (real(m_glb+1, dp)*real(n_glb+1, dp)*real(p_glb+1, dp))
                end do
            end do
        end do

#else
        ! CPU
        array_in(1:m+1, 1:n+1, 1:p+1) = fluid_indicator_function_I%sf(0:m, 0:n, 0:p)

        call fftw_execute_dft_r2c(plan_forward, array_in, array_out)

        array_out(:, :, :) = array_out(:, :, :) * array_kernelG_out(:, :, :)

        call fftw_execute_dft_c2r(plan_backward, array_out, array_in)

        filtered_fluid_indicator_function%sf(0:m, 0:n, 0:p) = array_in(1:m+1, 1:n+1, 1:p+1) / (real(m_glb+1, dp)*real(n_glb+1, dp)*real(p_glb+1, dp))
#endif

    end subroutine s_initialize_filtered_fluid_indicator_function

    !< apply the gaussian filter to the conservative variables and compute their filtered components
    subroutine s_apply_fftw_filter_cons(q_cons_vf, q_cons_filtered)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_filtered

        integer :: i, j, k, l, q

        do l = 1, sys_size-1
            call s_apply_fftw_filter_scalarfield(q_cons_filtered(advxb), .true., q_cons_vf(l), q_cons_filtered(l))
        end do 

    end subroutine s_apply_fftw_filter_cons

    !< applies the gaussian filter to an arbitrary scalar field
    subroutine s_apply_fftw_filter_scalarfield(filtered_fluid_indicator_function, fluid_quantity, q_temp_in, q_temp_out)
        type(scalar_field), intent(inout) :: q_temp_in
        type(scalar_field), intent(inout), optional :: q_temp_out
        type(scalar_field), intent(in) :: filtered_fluid_indicator_function

        logical :: fluid_quantity !< whether or not convolution integral is over V_f or V_p^(i) - integral over fluid volume or particle volume

        integer :: i, j, k

#if defined(MFC_OpenACC)
        if (fluid_quantity) then 
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m 
                do j = 0, n 
                    do k = 0, p 
                        real_data_gpu(i+1, j+1, k+1) = q_temp_in%sf(i, j, k) * fluid_indicator_function_I%sf(i, j, k)
                    end do
                end do
            end do
        else 
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m 
                do j = 0, n 
                    do k = 0, p 
                        real_data_gpu(i+1, j+1, k+1) = q_temp_in%sf(i, j, k) * (1.0_dp - fluid_indicator_function_I%sf(i, j, k))
                    end do
                end do
            end do
        end if
            
        ierr = cufftExecD2Z(forward_plan_gpu, real_data_gpu, complex_data_gpu)

        if (num_procs > 1) then 
            call s_mpi_perform_transpose_forward(complex_data_gpu)

            do i = 1, n+1
                do j = 1, (p+1)/2+1
                    do k = 1, m+1
                        complex_data_gpu_mpi(i, j, k) = complex_data_gpu_mpi(i, j, k) * complex_kernelG_gpu_mpi(i, j, k)
                    end do 
                end do 
            end do

            call s_mpi_perform_transpose_backward(complex_data_gpu)

        else ! 1 rank
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 1, m+1
                do j = 1, n+1
                    do k = 1, (p+1)/2+1
                        complex_data_gpu(i, j, k) = complex_data_gpu(i, j, k) * complex_kernelG_gpu(i, j, k)
                    end do
                end do
            end do
        end if

        ierr = cufftExecZ2D(backward_plan_gpu, complex_data_gpu, real_data_gpu)
        
        if (present(q_temp_out)) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m 
                do j = 0, n 
                    do k = 0, p 
                        q_temp_out%sf(i, j, k) = real_data_gpu(i+1, j+1, k+1) / (real(m_glb+1, dp)*real(n_glb+1, dp)*real(p_glb+1, dp) * filtered_fluid_indicator_function%sf(i, j, k))
                    end do
                end do
            end do
        else 
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m 
                do j = 0, n 
                    do k = 0, p 
                        q_temp_in%sf(i, j, k) = real_data_gpu(i+1, j+1, k+1) / (real(m_glb+1, dp)*real(n_glb+1, dp)*real(p_glb+1, dp) * filtered_fluid_indicator_function%sf(i, j, k))
                    end do
                end do
            end do
        end if 

#else 
        ! CPU
        array_in(1:m+1, 1:n+1, 1:p+1) = q_temp_in%sf(0:m, 0:n, 0:p) * fluid_indicator_function_I%sf(0:m, 0:n, 0:p)

        call fftw_execute_dft_r2c(plan_forward, array_in, array_out)

        array_out(:, :, :) = array_out(:, :, :) * array_kernelG_out(:, :, :)

        call fftw_execute_dft_c2r(plan_backward, array_out, array_in)

        if (present(q_temp_out)) then 
            q_temp_out%sf(0:m, 0:n, 0:p) = array_in(1:m+1, 1:n+1, 1:p+1) / (real(m_glb+1, dp)*real(n_glb+1, dp)*real(p_glb+1, dp) * filtered_fluid_indicator_function%sf(0:m, 0:n, 0:p)) 
        else 
            q_temp_in%sf(0:m, 0:n, 0:p) = array_in(1:m+1, 1:n+1, 1:p+1) / (real(m_glb+1, dp)*real(n_glb+1, dp)*real(p_glb+1, dp) * filtered_fluid_indicator_function%sf(0:m, 0:n, 0:p)) 
        end if
#endif

    end subroutine s_apply_fftw_filter_scalarfield

    !< apply the gaussian filter to the requisite tensors to compute unclosed terms of interest
    subroutine s_apply_fftw_filter_tensor(pt_Re_stress, R_mu, q_cons_filtered, rhs_rhouu, pImT_filtered)
        type(vector_field), dimension(1:num_dims), intent(inout) :: pt_Re_stress
        type(vector_field), dimension(1:num_dims), intent(inout) :: R_mu
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_filtered
        type(scalar_field), dimension(momxb:momxe), intent(inout) :: rhs_rhouu
        type(scalar_field), dimension(1:num_dims), intent(inout) :: pImT_filtered

        integer :: i, j, k, l, q

        ! pseudo turbulent reynolds stress
        do l = 1, num_dims 
            do q = 1, num_dims
                call s_apply_fftw_filter_scalarfield(q_cons_filtered(advxb), .true., pt_Re_stress(l)%vf(q))
            end do
        end do 

        ! effective viscosity
        do l = 1, num_dims 
            do q = 1, num_dims
                call s_apply_fftw_filter_scalarfield(q_cons_filtered(advxb), .true., R_mu(l)%vf(q))
            end do
        end do 

        ! interphase momentum exchange
        do l = 1, num_dims
            call s_apply_fftw_filter_scalarfield(q_cons_filtered(advxb), .false., rhs_rhouu(momxb-1+l), pImT_filtered(l))
        end do 

    end subroutine s_apply_fftw_filter_tensor

    !< wrapper for performing MPI data transpose steps to compute 3D DFT, FORWARD (only for flow quantities)
    subroutine s_mpi_perform_transpose_forward(complex_in_gpu_mpi)
        complex(dp), intent(inout) :: complex_in_gpu_mpi(m+1, n+1, (p+1)/2+1)

#if defined(MFC_OpenACC)
        call s_transpose_z2y_mpi(complex_in_gpu_mpi) ! -> gives complex_y_gpu

        ierr = cufftExecZ2Z(plan_y_gpu, complex_y_gpu, complex_y_gpu, CUFFT_FORWARD) 

        call s_transpose_y2x_mpi ! gives complex_x_gpu

        ierr = cufftExecZ2Z(plan_x_gpu, complex_x_gpu, complex_data_gpu_mpi, CUFFT_FORWARD) 
#endif

    end subroutine s_mpi_perform_transpose_forward
    
    !< wrapper for performing MPI data transpose steps to compute 3D DFT, BACKWARD
    subroutine s_mpi_perform_transpose_backward(complex_in_gpu_mpi)
        complex(dp), intent(inout) :: complex_in_gpu_mpi(m+1, n+1, (p+1)/2+1)

#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_x_gpu, complex_data_gpu_mpi, complex_x_gpu, CUFFT_INVERSE)

        call s_transpose_x2y_mpi !< gives complex_y_gpu

        ierr = cufftExecZ2Z(plan_y_gpu, complex_y_gpu, complex_y_gpu, CUFFT_INVERSE) 

        call s_transpose_y2z_mpi(complex_in_gpu_mpi) 
#endif

    end subroutine s_mpi_perform_transpose_backward

    subroutine s_transpose_z2y_mpi(complex_z_gpu)
        complex(dp), intent(inout) :: complex_z_gpu(m+1, n+1, (p+1)/2+1)
        complex(dp), allocatable :: sendbuf(:), recvbuf(:)
        integer :: sendcounts(num_procs), recvcounts(num_procs)
        integer :: sdispls(num_procs), rdispls(num_procs)
        integer :: total_size, offset
        integer :: i, j, k, r, s

#if defined(MFC_OpenACC)
        total_size = (m+1) * (n+1) * ((p+1)/2+1)

        allocate(sendbuf(total_size))
        allocate(recvbuf(total_size))

        !$acc parallel loop collapse(3) gang vector default(present) copy(sendbuf)
        do i = 1, m+1
            do j = 1, n+1
                do k = 1, (p+1)/2+1
                    offset = (i-1)*(n+1)*((p+1)/2+1) + (j-1)*((p+1)/2+1) + k
                    sendbuf(offset) = complex_z_gpu(i, j, k)
                end do
            end do
        end do

        do r = 1, num_procs
            sendcounts(r) = total_size / num_procs
            recvcounts(r) = total_size / num_procs
            sdispls(r) = (r-1) * sendcounts(r) + 1 
            rdispls(r) = (r-1) * recvcounts(r) + 1
        end do

        call MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_DOUBLE_COMPLEX, &
                           recvbuf, recvcounts, rdispls, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

        !$acc parallel loop collapse(3) gang vector default(present) copy(recvbuf)
        do i = 1, m+1
            do j = 1, (p+1)/2+1
                do k = 1, n+1
                    offset = (i-1)*((p+1)/2+1)*(n+1) + (j-1)*(n+1) + k
                    complex_y_gpu(i, j, k) = recvbuf(offset)
                end do
            end do
        end do

        deallocate(sendbuf, recvbuf)
#endif
    end subroutine s_transpose_z2y_mpi

    subroutine s_transpose_y2x_mpi
        complex(dp), allocatable :: sendbuf(:), recvbuf(:)
        integer :: sendcounts(num_procs), recvcounts(num_procs)
        integer :: sdispls(num_procs), rdispls(num_procs)
        integer :: total_size, offset
        integer :: i, j, k, s, r

#if defined(MFC_OpenACC)
        total_size = (m+1) * (n+1) * ((p+1)/2+1)

        allocate(sendbuf(total_size))
        allocate(recvbuf(total_size))

        !$acc parallel loop collapse(3) gang vector default(present) copy(sendbuf)
        do i = 1, m+1
            do j = 1, (p+1)/2+1
                do k = 1, n+1
                    offset = (i-1)*((p+1)/2+1)*(n+1) + (j-1)*(n+1) + k
                    sendbuf(offset) = complex_y_gpu(i, j, k)
                end do
            end do
        end do

        do r = 1, num_procs
            sendcounts(r) = total_size / num_procs
            recvcounts(r) = total_size / num_procs
            sdispls(r) = (r-1) * sendcounts(r) + 1
            rdispls(r) = (r-1) * recvcounts(r) + 1
        end do

        call MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_DOUBLE_COMPLEX, &
                           recvbuf, recvcounts, rdispls, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

        !$acc parallel loop collapse(3) gang vector default(present) copy(recvbuf)
        do i = 1, n+1
            do j = 1, (p+1)/2+1
                do k = 1, m+1
                    offset = (i-1)*((p+1)/2+1)*(m+1) + (j-1)*(m+1) + k 
                    complex_x_gpu(i, j, k) = recvbuf(offset)
                end do
            end do
        end do

        deallocate(sendbuf, recvbuf)
#endif
    end subroutine s_transpose_y2x_mpi

    subroutine s_transpose_x2y_mpi      
        complex(dp), allocatable :: sendbuf(:), recvbuf(:)
        integer :: sendcounts(num_procs), recvcounts(num_procs)
        integer :: sdispls(num_procs), rdispls(num_procs)
        integer :: total_size, offset
        integer :: i, j, k, r
      
#if defined(MFC_OpenACC)
        total_size = (m+1) * (n+1) * ((p+1)/2+1)

        allocate(sendbuf(total_size))
        allocate(recvbuf(total_size))
      
        offset = 0
        !$acc parallel loop collapse(3) gang vector default(present) copy(offset, sendbuf)
        do i = 1, n+1
            do j = 1, (p+1)/2+1
                do k = 1, m+1
                    offset = offset + 1
                    sendbuf(offset) = complex_x_gpu(i, j, k)
                end do
            end do
        end do
      
        do r = 1, num_procs
            sendcounts(r) = total_size / num_procs
            recvcounts(r) = total_size / num_procs
            sdispls(r) = (r-1) * sendcounts(r) + 1
            rdispls(r) = (r-1) * recvcounts(r) + 1
        end do
      
        call MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_DOUBLE_COMPLEX, &
                           recvbuf, recvcounts, rdispls, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
      
        offset = 0
        !$acc parallel loop collapse(3) gang vector default(present) copy(offset, recvbuf)
        do i = 1, m+1
            do j = 1, (p+1)/2+1
                do k = 1, n+1
                    offset = offset + 1
                    complex_y_gpu(i, j, k) = recvbuf(offset)
                end do
            end do
        end do
      
        deallocate(sendbuf, recvbuf)
#endif
    end subroutine s_transpose_x2y_mpi

    subroutine s_transpose_y2z_mpi(complex_z_gpu)
        complex(dp), intent(inout) :: complex_z_gpu(m+1, n+1, (p+1)/2+1)
        complex(dp), allocatable :: sendbuf(:), recvbuf(:)
        integer :: sendcounts(num_procs), recvcounts(num_procs)
        integer :: sdispls(num_procs), rdispls(num_procs)
        integer :: total_size, offset
        integer :: i, j, k, r
      
#if defined(MFC_OpenACC)
        total_size = (m+1) * (n+1) * ((p+1)/2+1)

        allocate(sendbuf(total_size))
        allocate(recvbuf(total_size))
      
        offset = 0
        !$acc parallel loop collapse(3) gang vector default(present) copy(offset, sendbuf)
        do i = 1, m+1
            do j = 1, (p+1)/2+1
                do k = 1, n+1
                    offset = offset + 1
                    sendbuf(offset) = complex_y_gpu(i, j, k)
                end do
            end do
        end do
      
        do r = 1, num_procs
            sendcounts(r) = total_size / num_procs
            recvcounts(r) = total_size / num_procs
            sdispls(r) = (r-1) * sendcounts(r) + 1
            rdispls(r) = (r-1) * recvcounts(r) + 1
        end do
      
        call MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_DOUBLE_COMPLEX, &
                           recvbuf, recvcounts, rdispls, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
      
        offset = 0
        !$acc parallel loop collapse(3) gang vector default(present) copy(offset, recvbuf) 
        do i = 1, m+1
            do j = 1, n+1
                do k = 1, (p+1)/2+1
                    offset = offset + 1
                    complex_z_gpu(i, j, k) = recvbuf(offset)
                end do
            end do
        end do
      
        deallocate(sendbuf, recvbuf)
#endif
    end subroutine s_transpose_y2z_mpi

    !>  The purpose of this subroutine is to apply a Fourier low-
        !!      pass filter to the flow variables in the azimuthal direction
        !!      to remove the high-frequency content. This alleviates the
        !!      restrictive CFL condition arising from cells near the axis.
        !! @param q_cons_vf Conservative variables
    subroutine s_apply_fourier_filter(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(c_double), pointer :: p_real(:)
        complex(c_double_complex), pointer :: p_cmplx(:), p_fltr_cmplx(:)
        integer :: i, j, k, l !< Generic loop iterators

        ! Restrict filter to processors that have cells adjacent to axis
        if (bc_y%beg >= 0) return
#if defined(MFC_OpenACC)

        !$acc parallel loop collapse(3) gang vector default(present)
        do k = 1, sys_size
            do j = 0, m
                do l = 1, cmplx_size
                    data_fltr_cmplx_gpu(l + j*cmplx_size + (k - 1)*cmplx_size*x_size) = (0_dp, 0_dp)
                end do
            end do
        end do

        !$acc parallel loop collapse(3) gang vector default(present)
        do k = 1, sys_size
            do j = 0, m
                do l = 0, p
                    data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size) = q_cons_vf(k)%sf(j, 0, l)
                end do
            end do
        end do

        p_real => data_real_gpu
        p_cmplx => data_cmplx_gpu
        p_fltr_cmplx => data_fltr_cmplx_gpu

!$acc data attach(p_real, p_cmplx, p_fltr_cmplx)
!$acc host_data use_device(p_real, p_cmplx, p_fltr_cmplx)
#if defined(__PGI)
        ierr = cufftExecD2Z(fwd_plan_gpu, data_real_gpu, data_cmplx_gpu)
#else
        ierr = hipfftExecD2Z(fwd_plan_gpu, c_loc(p_real), c_loc(p_cmplx))
        call hipCheck(hipDeviceSynchronize())
#endif
        !$acc end host_data
        Nfq = 3
        !$acc update device(Nfq)

        !$acc parallel loop collapse(3) gang vector default(present)
        do k = 1, sys_size
            do j = 0, m
                do l = 1, Nfq
                    data_fltr_cmplx_gpu(l + j*cmplx_size + (k - 1)*cmplx_size*x_size) = data_cmplx_gpu(l + j*cmplx_size + (k - 1)*cmplx_size*x_size)
                end do
            end do
        end do

!$acc host_data use_device(p_real, p_fltr_cmplx)
#if defined(__PGI)
        ierr = cufftExecZ2D(bwd_plan_gpu, data_fltr_cmplx_gpu, data_real_gpu)
#else
        ierr = hipfftExecZ2D(bwd_plan_gpu, c_loc(p_fltr_cmplx), c_loc(p_real))
        call hipCheck(hipDeviceSynchronize())
#endif
        !$acc end host_data

        !$acc parallel loop collapse(3) gang vector default(present)
        do k = 1, sys_size
            do j = 0, m
                do l = 0, p
                    data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size) = data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size)/real(real_size, dp)
                    q_cons_vf(k)%sf(j, 0, l) = data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size)
                end do
            end do
        end do

        do i = 1, fourier_rings

            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 1, sys_size
                do j = 0, m
                    do l = 1, cmplx_size
                        data_fltr_cmplx_gpu(l + j*cmplx_size + (k - 1)*cmplx_size*x_size) = (0_dp, 0_dp)
                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present) firstprivate(i)
            do k = 1, sys_size
                do j = 0, m
                    do l = 0, p
                        data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size) = q_cons_vf(k)%sf(j, i, l)
                    end do
                end do
            end do

!$acc host_data use_device(p_real, p_cmplx)
#if defined(__PGI)
            ierr = cufftExecD2Z(fwd_plan_gpu, data_real_gpu, data_cmplx_gpu)
#else
            ierr = hipfftExecD2Z(fwd_plan_gpu, c_loc(p_real), c_loc(p_cmplx))
            call hipCheck(hipDeviceSynchronize())
#endif
            !$acc end host_data

            Nfq = min(floor(2_dp*real(i, dp)*pi), cmplx_size)
            !$acc update device(Nfq)

            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 1, sys_size
                do j = 0, m
                    do l = 1, Nfq
                        data_fltr_cmplx_gpu(l + j*cmplx_size + (k - 1)*cmplx_size*x_size) = data_cmplx_gpu(l + j*cmplx_size + (k - 1)*cmplx_size*x_size)
                    end do
                end do
            end do

!$acc host_data use_device(p_real, p_fltr_cmplx)
#if defined(__PGI)
            ierr = cufftExecZ2D(bwd_plan_gpu, data_fltr_cmplx_gpu, data_real_gpu)
#else
            ierr = hipfftExecZ2D(bwd_plan_gpu, c_loc(p_fltr_cmplx), c_loc(p_real))
            call hipCheck(hipDeviceSynchronize())
#endif
            !$acc end host_data

            !$acc parallel loop collapse(3) gang vector default(present) firstprivate(i)
            do k = 1, sys_size
                do j = 0, m
                    do l = 0, p
                        data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size) = data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size)/real(real_size, dp)
                        q_cons_vf(k)%sf(j, i, l) = data_real_gpu(l + j*real_size + 1 + (k - 1)*real_size*x_size)
                    end do
                end do
            end do

        end do

#else
        Nfq = 3
        do j = 0, m
            do k = 1, sys_size
                data_fltr_cmplx(:) = (0_dp, 0_dp)
                data_real(1:p + 1) = q_cons_vf(k)%sf(j, 0, 0:p)
                call fftw_execute_dft_r2c(fwd_plan, data_real, data_cmplx)
                data_fltr_cmplx(1:Nfq) = data_cmplx(1:Nfq)
                call fftw_execute_dft_c2r(bwd_plan, data_fltr_cmplx, data_real)
                data_real(:) = data_real(:)/real(real_size, dp)
                q_cons_vf(k)%sf(j, 0, 0:p) = data_real(1:p + 1)
            end do
        end do

        ! Apply Fourier filter to additional rings
        do i = 1, fourier_rings
            Nfq = min(floor(2_dp*real(i, dp)*pi), cmplx_size)
            do j = 0, m
                do k = 1, sys_size
                    data_fltr_cmplx(:) = (0_dp, 0_dp)
                    data_real(1:p + 1) = q_cons_vf(k)%sf(j, i, 0:p)
                    call fftw_execute_dft_r2c(fwd_plan, data_real, data_cmplx)
                    data_fltr_cmplx(1:Nfq) = data_cmplx(1:Nfq)
                    call fftw_execute_dft_c2r(bwd_plan, data_fltr_cmplx, data_real)
                    data_real(:) = data_real(:)/real(real_size, dp)
                    q_cons_vf(k)%sf(j, i, 0:p) = data_real(1:p + 1)
                end do
            end do
        end do
#endif
!$acc end data
    end subroutine s_apply_fourier_filter

    !>  The purpose of this subroutine is to destroy the fftw plan
        !!      that will be used in the forward and backward DFTs when
        !!      applying the Fourier filter in the azimuthal direction.
    subroutine s_finalize_fftw_module

#if defined(MFC_OpenACC)
        @:DEALLOCATE(data_real_gpu, data_fltr_cmplx_gpu, data_cmplx_gpu)
#if defined(__PGI)

        ierr = cufftDestroy(fwd_plan_gpu)
        ierr = cufftDestroy(bwd_plan_gpu)
#else
        ierr = hipfftDestroy(fwd_plan_gpu)
        ierr = hipfftDestroy(bwd_plan_gpu)
#endif
#else
        call fftw_free(fftw_real_data)
        call fftw_free(fftw_cmplx_data)
        call fftw_free(fftw_fltr_cmplx_data)

        call fftw_destroy_plan(fwd_plan)
        call fftw_destroy_plan(bwd_plan)
#endif

    end subroutine s_finalize_fftw_module

    subroutine s_finalize_fftw_explicit_filter_module
        @:DEALLOCATE(fluid_indicator_function_I%sf)

#if defined(MFC_OpenACC)
        @:DEALLOCATE(real_data_gpu, complex_data_gpu, real_kernelG_gpu, complex_kernelG_gpu)
        @:DEALLOCATE(complex_data_gpu_mpi, complex_kernelG_gpu_mpi, complex_y_gpu, complex_x_gpu)

        ierr = cufftDestroy(forward_plan_gpu)
        ierr = cufftDestroy(backward_plan_gpu)

        ierr = cufftDestroy(forward_plan_kernelG_gpu)
#else
        call fftw_free(p_real)
        call fftw_free(p_complex)
        call fftw_free(p_kernelG_real)
        call fftw_free(p_kernelG_complex)

        call fftw_destroy_plan(plan_forward)
        call fftw_destroy_plan(plan_backward)
        call fftw_destroy_plan(plan_kernelG_forward)
#endif

    end subroutine s_finalize_fftw_explicit_filter_module

end module m_fftw
