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
 s_apply_fftw_filter_tensor, s_apply_fftw_filter_scalarfield

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

    ! ============ NEW ADDITIONS ============

    ! fluid indicator function (1 = fluid, 0 = otherwise)
    type(scalar_field), public :: fluid_indicator_function_I

    !$acc declare create(fluid_indicator_function_I)

#if defined(MFC_OpenACC)
    ! GPU plans
    integer :: plan_x_fwd_gpu, plan_x_bwd_gpu, plan_y_gpu, plan_z_gpu
#else
    ! CPU plans
    type(c_ptr) :: plan_x_r2c_fwd, plan_x_c2r_bwd
    type(c_ptr) :: plan_y_c2c_fwd, plan_y_c2c_bwd 
    type(c_ptr) :: plan_z_c2c_fwd, plan_z_c2c_bwd
    type(c_ptr) :: plan_x_r2c_kernelG, plan_y_c2c_kernelG, plan_z_c2c_kernelG
#endif

    ! domain size information
    integer :: Nx, Ny, Nz, NxC, Nyloc, Nzloc

    ! 1D real and complex vectors for FFT routines
    real(c_double), allocatable :: data_real_in1d(:) 
    complex(c_double_complex), allocatable :: data_cmplx_out1d(:)
    complex(c_double_complex), allocatable :: data_cmplx_out1dy(:)

    ! 3D arrays for slab transposes
    complex(c_double_complex), allocatable :: data_cmplx_slabz(:, :, :), data_cmplx_slaby(:, :, :)

    ! input array for FFT
    real(c_double), allocatable :: data_real_3D_slabz(:, :, :)

    ! filtering kernel in physical space
    real(c_double), allocatable :: real_kernelG_in(:, :, :)

    ! FFT of filtering kernel
    complex(c_double_complex), allocatable :: cmplx_kernelG1d(:)

    !$acc declare create(Nx, Ny, Nz, NxC, Nyloc, Nzloc)
    !$acc declare create(data_real_in1d, data_cmplx_out1d, data_cmplx_out1dy, data_cmplx_slabz, data_cmplx_slaby, data_real_3D_slabz, real_kernelG_in, cmplx_kernelG1d)

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

    !< create fft plans to be used for explicit filtering of data 
    subroutine s_initialize_fftw_explicit_filter_module
        integer :: size_n(1), inembed(1), onembed(1)

        print *, 'FFTW SETUP...'

        @:ALLOCATE(fluid_indicator_function_I%sf(0:m, 0:n, 0:p))
        @:ACC_SETUP_SFs(fluid_indicator_function_I)

        !< global sizes 
        Nx = m_glb + 1
        Ny = n_glb + 1
        Nz = p_glb + 1

        !< complex size
        NxC = Nx/2 + 1

        !< local sizes on each processor
        Nyloc = Ny / num_procs
        Nzloc = p + 1

        !$acc update device(Nx, Ny, Nz, NxC, Nyloc, Nzloc)

        @:ALLOCATE(data_real_in1d(Nx*Ny*Nzloc))
        @:ALLOCATE(data_cmplx_out1d(NxC*Ny*Nz/num_procs))
        @:ALLOCATE(data_cmplx_out1dy(NxC*Ny*Nz/num_procs))
        @:ALLOCATE(cmplx_kernelG1d(NxC*Nyloc*Nz))
        @:ALLOCATE(real_kernelG_in(Nx, Ny, Nzloc))
        @:ALLOCATE(data_real_3D_slabz(Nx, Ny, Nzloc))
        @:ALLOCATE(data_cmplx_slabz(NxC, Ny, Nzloc))
        @:ALLOCATE(data_cmplx_slaby(NxC, Nyloc, Nz))

#if defined(MFC_OpenACC)
        !< GPU FFT plans
        !< X - plans
        size_n(1) = Nx
        inembed(1) = Nx
        onembed(1) = NxC
        ierr = cufftPlanMany(plan_x_fwd_gpu, 1, size_n, inembed, 1, Nx, onembed, 1, NxC, CUFFT_D2Z, Ny*Nzloc)
        size_n(1) = Nx
        inembed(1) = NxC
        onembed(1) = Nx  
        ierr = cufftPlanMany(plan_x_bwd_gpu, 1, size_n, inembed, 1, NxC, onembed, 1, Nx, CUFFT_Z2D, Ny*Nzloc)
        !< Y - plans
        size_n(1) = Ny
        inembed(1) = Ny
        onembed(1) = Ny
        ierr = cufftPlanMany(plan_y_gpu, 1, size_n, inembed, 1, Ny, onembed, 1, Ny, CUFFT_Z2Z, NxC*Nzloc)
        !< Z - plans
        size_n(1) = Nz 
        inembed(1) = Nz 
        onembed(1) = Nz 
        ierr = cufftPlanMany(plan_z_gpu, 1, size_n, inembed, 1, Nz, onembed, 1, Nz, CUFFT_Z2Z, NxC*Nyloc)
#else
        !< CPU FFT plans
        !< X - direction plans
        size_n(1) = Nx
        inembed(1) = Nx
        onembed(1) = NxC
        plan_x_r2c_fwd = fftw_plan_many_dft_r2c(1, size_n, Ny*Nzloc, &                  ! rank, n, howmany
                                                data_real_in1d, inembed, 1, Nx, &       ! in, inembed, istride, idist
                                                data_cmplx_out1d, onembed, 1, NxC, &    ! out, onembed, ostride, odist
                                                FFTW_MEASURE)                           ! sign, flags
        size_n(1) = Nx
        inembed(1) = NxC
        onembed(1) = Nx                                                         
        plan_x_c2r_bwd = fftw_plan_many_dft_c2r(1, size_n, Ny*Nzloc, & 
                                                data_cmplx_out1d, inembed, 1, NxC, & 
                                                data_real_in1d, onembed, 1, Nx, & 
                                                FFTW_MEASURE)
        !< Y - direction plans
        size_n(1) = Ny
        inembed(1) = Ny
        onembed(1) = Ny
        plan_y_c2c_fwd = fftw_plan_many_dft(1, size_n, NxC*Nzloc, & 
                                            data_cmplx_out1dy, inembed, 1, Ny, & 
                                            data_cmplx_out1dy, onembed, 1, Ny, & 
                                            FFTW_FORWARD, FFTW_MEASURE)
        plan_y_c2c_bwd = fftw_plan_many_dft(1, size_n, NxC*Nzloc, & 
                                            data_cmplx_out1dy, inembed, 1, Ny, & 
                                            data_cmplx_out1dy, onembed, 1, Ny, & 
                                            FFTW_BACKWARD, FFTW_MEASURE)
        !< Z - direction plans
        size_n(1) = Nz 
        inembed(1) = Nz 
        onembed(1) = Nz 
        plan_z_c2c_fwd = fftw_plan_many_dft(1, size_n, NxC*Nyloc, & 
                                            data_cmplx_out1d, inembed, 1, Nz, & 
                                            data_cmplx_out1d, onembed, 1, Nz, & 
                                            FFTW_FORWARD, FFTW_MEASURE)
        plan_z_c2c_bwd = fftw_plan_many_dft(1, size_n, NxC*Nyloc, & 
                                            data_cmplx_out1d, inembed, 1, Nz, &
                                            data_cmplx_out1d, onembed, 1, Nz, & 
                                            FFTW_BACKWARD, FFTW_MEASURE)
        ! forward plans for filtering kernel
        ! X kernel plan
        size_n(1) = Nx
        inembed(1) = Nx
        onembed(1) = NxC
        plan_x_r2c_kernelG = fftw_plan_many_dft_r2c(1, size_n, Ny*Nzloc, &                    
                                                    data_real_in1d, inembed, 1, Nx, &        
                                                    cmplx_kernelG1d, onembed, 1, NxC, &    
                                                    FFTW_MEASURE)          
        ! Y kernel plan                  
        size_n(1) = Ny
        inembed(1) = Ny
        onembed(1) = Ny
        plan_y_c2c_kernelG = fftw_plan_many_dft(1, size_n, NxC*Nzloc, & 
                                                data_cmplx_out1dy, inembed, 1, Ny, & 
                                                data_cmplx_out1dy, onembed, 1, Ny, & 
                                                FFTW_FORWARD, FFTW_MEASURE)
        ! Z kernel plan
        size_n(1) = Nz 
        inembed(1) = Nz 
        onembed(1) = Nz 
        plan_z_c2c_kernelG = fftw_plan_many_dft(1, size_n, NxC*Nyloc, & 
                                                cmplx_kernelG1d, inembed, 1, Nz, & 
                                                cmplx_kernelG1d, onembed, 1, Nz, & 
                                                FFTW_FORWARD, FFTW_MEASURE)
#endif
    end subroutine s_initialize_fftw_explicit_filter_module

    !< initialize the gaussian filtering kernel in real space and then compute its DFT
    subroutine s_initialize_filtering_kernel
        real(dp) :: sigma_stddev
        real(dp) :: Lx, Ly, Lz
        real(dp) :: x_r, y_r, z_r  
        real(dp) :: r2
        real(dp) :: G_norm_int, G_norm_int_glb
        integer :: i, j, k, idx

        ! gaussian filter
        sigma_stddev = 3.0_dp * 0.05_dp

        Lx = x_domain_end_glb - x_domain_beg_glb
        Ly = y_domain_end_glb - y_domain_beg_glb  
        Lz = z_domain_end_glb - z_domain_beg_glb    
        
        G_norm_int = 0.0_dp
   
        !$acc parallel loop collapse(3) gang vector default(present) reduction(+:G_norm_int) copyin(Lx, Ly, Lz, sigma_stddev) private(x_r, y_r, z_r, r2)
        do i = 0, m 
            do j = 0, n 
                do k = 0, p 
                    x_r = min(abs(x_cc(i) - x_domain_beg_glb), Lx - abs(x_cc(i) - x_domain_beg_glb))
                    y_r = min(abs(y_cc(j) - y_domain_beg_glb), Ly - abs(y_cc(j) - y_domain_beg_glb))
                    z_r = min(abs(z_cc(k) - z_domain_beg_glb), Lz - abs(z_cc(k) - z_domain_beg_glb))

                    r2 = x_r**2 + y_r**2 + z_r**2

                    real_kernelG_in(i+1, j+1, k+1) = exp(-r2 / (2.0_dp*sigma_stddev**2))

                    G_norm_int = G_norm_int + real_kernelG_in(i+1, j+1, k+1)*dx(i)*dy(j)*dz(k)
                end do 
            end do
        end do

        call s_mpi_allreduce_sum(G_norm_int, G_norm_int_glb) 

        ! FFT of kernel
        ! normalize kernel
        !$acc parallel loop collapse(3) gang vector default(present) copyin(G_norm_int_glb)
        do i = 1, Nx 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_real_3D_slabz(i, j, k) = real_kernelG_in(i, j, k) / G_norm_int_glb
                end do 
            end do 
        end do 

        ! 3D z-slab -> 1D x, y, z
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, Nx 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_real_in1d(i + (j-1)*Nx + (k-1)*Nx*Ny) = data_real_3D_slabz(i, j, k)
                end do 
            end do 
        end do

        ! X FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecD2Z(plan_x_fwd_gpu, data_real_in1d, cmplx_kernelG1d)
#else
        call fftw_execute_dft_r2c(plan_x_r2c_kernelG, data_real_in1d, cmplx_kernelG1d)
#endif

        ! 1D x, y, z -> 1D y, x, z (CMPLX)
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_cmplx_out1dy(j + (i-1)*Ny + (k-1)*Ny*NxC) = cmplx_kernelG1d(i + (j-1)*NxC + (k-1)*NxC*Ny)
                end do 
            end do 
        end do

        ! Y FFT 
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_FORWARD)
#else
        call fftw_execute_dft(plan_y_c2c_kernelG, data_cmplx_out1dy, data_cmplx_out1dy)
#endif

        ! 1D y, x, z -> 3D z-slab
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_cmplx_slabz(i, j, k) = data_cmplx_out1dy(j + (i-1)*Ny + (k-1)*Ny*NxC)
                end do 
            end do 
        end do 

        ! transpose z-slab to y-slab
        call s_mpi_transpose_slabZ2Y 

        ! 3D y-slab -> 1D z, x, y
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Nyloc 
                do k = 1, Nz
                    cmplx_kernelG1d(k + (i-1)*Nz + (j-1)*Nz*NxC) = data_cmplx_slaby(i, j, k)
                end do 
            end do 
        end do

        ! Z FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_z_gpu, cmplx_kernelG1d, cmplx_kernelG1d, CUFFT_FORWARD)
#else
        call fftw_execute_dft(plan_z_c2c_kernelG, cmplx_kernelG1d, cmplx_kernelG1d)
#endif

        ! normalize FFT 
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Nyloc 
                do k = 1, Nz
                    cmplx_kernelG1d(k + (i-1)*Nz + (j-1)*Nz*NxC) = cmplx_kernelG1d(k + (i-1)*Nz + (j-1)*Nz*NxC) / (real(Nx*Ny*Nz, dp))
                end do 
            end do 
        end do

        ! return cmplx_kernelG1d: 1D z, x, y
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
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, Nx 
            do j = 1, Ny 
                do k = 1, Nzloc 
                    data_real_3D_slabz(i, j, k) = fluid_indicator_function_I%sf(i-1, j-1, k-1)
                end do 
            end do 
        end do 

        call s_mpi_FFT_fwd 

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Nyloc 
                do k = 1, Nz
                    data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) = data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) * cmplx_kernelG1d(k + (i-1)*Nz + (j-1)*Nz*NxC)
                end do
            end do 
        end do

        call s_mpi_FFT_bwd

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, Nx 
            do j = 1, Ny 
                do k = 1, Nzloc
                    filtered_fluid_indicator_function%sf(i-1, j-1, k-1) = data_real_3D_slabz(i, j, k) / (real(Nx*Ny*Nz, dp))
                end do 
            end do
        end do

    end subroutine s_initialize_filtered_fluid_indicator_function

    !< apply the gaussian filter to the conservative variables and compute their filtered components
    subroutine s_apply_fftw_filter_cons(q_cons_vf, q_cons_filtered)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_filtered

        integer :: l

        do l = 1, sys_size-1
            call s_apply_fftw_filter_scalarfield(q_cons_filtered(advxb), .true., q_cons_vf(l), q_cons_filtered(l))
        end do 

    end subroutine s_apply_fftw_filter_cons

    !< applies the gaussian filter to an arbitrary scalar field
    subroutine s_apply_fftw_filter_scalarfield(filtered_fluid_indicator_function, fluid_quantity, q_temp_in, q_temp_out)
        type(scalar_field), intent(in) :: filtered_fluid_indicator_function
        type(scalar_field), intent(inout) :: q_temp_in
        type(scalar_field), intent(inout), optional :: q_temp_out

        logical, intent(in) :: fluid_quantity !< whether or not convolution integral is over V_f or V_p^(i) - integral over fluid volume or particle volume

        integer :: i, j, k

        if (fluid_quantity) then 
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m 
                do j = 0, n 
                    do k = 0, p 
                        data_real_3D_slabz(i+1, j+1, k+1) = q_temp_in%sf(i, j, k) * fluid_indicator_function_I%sf(i, j, k)
                    end do 
                end do 
            end do
        else 
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m 
                do j = 0, n 
                    do k = 0, p 
                        data_real_3D_slabz(i+1, j+1, k+1) = q_temp_in%sf(i, j, k) * (1.0_dp - fluid_indicator_function_I%sf(i, j, k))
                    end do 
                end do 
            end do
        end if

        call s_mpi_FFT_fwd 

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Nyloc 
                do k = 1, Nz 
                    data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) = data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) * cmplx_kernelG1d(k + (i-1)*Nz + (j-1)*Nz*NxC)
                end do 
            end do 
        end do

        call s_mpi_FFT_bwd

        if (present(q_temp_out)) then 
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m
                do j = 0, n
                    do k = 0, p
                        q_temp_out%sf(i, j, k) = data_real_3D_slabz(i+1, j+1, k+1) / (real(Nx*Ny*Nz, dp) * filtered_fluid_indicator_function%sf(i, j, k))
                    end do 
                end do 
            end do
        else 
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m
                do j = 0, n 
                    do k = 0, p 
                        q_temp_in%sf(i, j, k) = data_real_3D_slabz(i+1, j+1, k+1) / (real(Nx*Ny*Nz, dp) * filtered_fluid_indicator_function%sf(i, j, k))      
                    end do 
                end do 
            end do
        end if

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

    !< transpose domain from z-slabs to y-slabs on each processor
    subroutine s_mpi_transpose_slabZ2Y
        complex(c_double_complex), allocatable :: sendbuf(:), recvbuf(:)
        integer :: dest_rank, src_rank
        integer :: i, j, k

        allocate(sendbuf(NxC*Nyloc*Nzloc*num_procs))
        allocate(recvbuf(NxC*Nyloc*Nzloc*num_procs))

        !$acc parallel loop gang vector default(present) copy(sendbuf) private(dest_rank)
        do dest_rank = 0, num_procs-1
            !$acc loop collapse(3) gang vector default(present) 
            do k = 1, Nzloc 
                do j = dest_rank*Nyloc+1, (dest_rank+1)*Nyloc
                    do i = 1, NxC
                        sendbuf(i + (j-(dest_rank*Nyloc+1))*NxC + (k-1)*NxC*Nyloc + dest_rank*NxC*Nyloc*Nzloc) = data_cmplx_slabz(i, j, k)
                    end do 
                end do
            end do
        end do

        call MPI_Alltoall(sendbuf, NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, & 
                          recvbuf, NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

        !$acc parallel loop gang vector default(present) copy(recvbuf) private(src_rank)
        do src_rank = 0, num_procs-1
            !$acc loop collapse(3) gang vector default(present) 
            do k = src_rank*Nzloc+1, (src_rank+1)*Nzloc
                do j = 1, Nyloc
                    do i = 1, NxC
                        data_cmplx_slaby(i, j, k) = recvbuf(i + (j-1)*NxC + (k-(src_rank*Nzloc+1))*NxC*Nyloc + src_rank*NxC*Nyloc*Nzloc)
                    end do 
                end do
            end do 
        end do

        deallocate(sendbuf, recvbuf)
    end subroutine s_mpi_transpose_slabZ2Y

    !< transpose domain from y-slabs to z-slabs on each processor
    subroutine s_mpi_transpose_slabY2Z 
        complex(c_double_complex), allocatable :: sendbuf(:), recvbuf(:)
        integer :: dest_rank, src_rank
        integer :: i, j, k

        allocate(sendbuf(NxC*Nyloc*Nzloc*num_procs))
        allocate(recvbuf(NxC*Nyloc*Nzloc*num_procs))

        !$acc parallel loop gang vector default(present) copy(sendbuf) private(dest_rank)
        do dest_rank = 0, num_procs-1
            !$acc parallel loop collapse(3) gang vector default(present) 
            do k = dest_rank*Nzloc+1, (dest_rank+1)*Nzloc
                do j = 1, Nyloc 
                    do i = 1, NxC 
                        sendbuf(i + (j-1)*NxC + (k-(dest_rank*Nzloc+1))*NxC*Nyloc + dest_rank*NxC*Nyloc*Nzloc) = data_cmplx_slaby(i, j, k)
                    end do 
                end do 
            end do 
        end do

        call MPI_Alltoall(sendbuf, NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, & 
                          recvbuf, NxC*Nyloc*Nzloc, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

        !$acc parallel loop gang vector default(present) copy(recvbuf) private(src_rank)
        do src_rank = 0, num_procs-1
            !$acc parallel loop collapse(3) gang vector default(present)
            do k = 1, Nzloc
                do j = src_rank*Nyloc+1, (src_rank+1)*Nyloc
                    do i = 1, NxC 
                        data_cmplx_slabz(i, j, k) = recvbuf(i + (j-(src_rank*Nyloc+1))*NxC + (k-1)*NxC*Nyloc + src_rank*NxC*Nyloc*Nzloc)
                    end do 
                end do
            end do 
        end do
        
        deallocate(sendbuf, recvbuf)
    end subroutine s_mpi_transpose_slabY2Z

    !< compute forward FFT, input: data_real_3D_slabz, output: data_cmplx_out1d
    subroutine s_mpi_FFT_fwd
        integer :: i, j, k

        ! 3D z-slab -> 1D x, y, z
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, Nx 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_real_in1d(i + (j-1)*Nx + (k-1)*Nx*Ny) = data_real_3D_slabz(i, j, k)
                end do 
            end do 
        end do

        ! X FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecD2Z(plan_x_fwd_gpu, data_real_in1d, data_cmplx_out1d)
#else
        call fftw_execute_dft_r2c(plan_x_r2c_fwd, data_real_in1d, data_cmplx_out1d)
#endif

        ! 1D x, y, z -> 1D y, x, z (CMPLX)
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_cmplx_out1dy(j + (i-1)*Ny + (k-1)*Ny*NxC) = data_cmplx_out1d(i + (j-1)*NxC + (k-1)*NxC*Ny)
                end do 
            end do 
        end do

        ! Y FFT 
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_FORWARD)
#else
        call fftw_execute_dft(plan_y_c2c_fwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif 

        ! 1D y, x, z -> 3D z-slab
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_cmplx_slabz(i, j, k) = data_cmplx_out1dy(j + (i-1)*Ny + (k-1)*Ny*NxC)
                end do 
            end do 
        end do 

        ! transpose z-slab to y-slab
        call s_mpi_transpose_slabZ2Y 

        ! 3D y-slab -> 1D z, x, y
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Nyloc 
                do k = 1, Nz
                    data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC) = data_cmplx_slaby(i, j, k)
                end do 
            end do 
        end do

        ! Z FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_z_gpu, data_cmplx_out1d, data_cmplx_out1d, CUFFT_FORWARD)
#else
        call fftw_execute_dft(plan_z_c2c_fwd, data_cmplx_out1d, data_cmplx_out1d)
#endif

        ! return data_cmplx_out1d: 1D z, x, y
    end subroutine s_mpi_FFT_fwd

    !< compute inverse FFT, input: data_cmplx_out1d, output: data_real_3D_slabz
    subroutine s_mpi_FFT_bwd
        integer :: i, j, k

        ! Z inv FFT 
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_z_gpu, data_cmplx_out1d, data_cmplx_out1d, CUFFT_INVERSE)
#else
        call fftw_execute_dft(plan_z_c2c_bwd, data_cmplx_out1d, data_cmplx_out1d)
#endif

        ! 1D z, x, y -> 3D y-slab
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Nyloc 
                do k = 1, Nz 
                    data_cmplx_slaby(i, j, k) = data_cmplx_out1d(k + (i-1)*Nz + (j-1)*Nz*NxC)
                end do 
            end do 
        end do

        ! transpose y-slab to z-slab
        call s_mpi_transpose_slabY2Z

        ! 3D z-slab -> 1D y, x, z
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_cmplx_out1dy(j + (i-1)*Ny + (k-1)*Ny*NxC) = data_cmplx_slabz(i, j, k)
                end do 
            end do 
        end do

        ! Y inv FFT 
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2Z(plan_y_gpu, data_cmplx_out1dy, data_cmplx_out1dy, CUFFT_INVERSE)
#else
        call fftw_execute_dft(plan_y_c2c_bwd, data_cmplx_out1dy, data_cmplx_out1dy)
#endif

        ! 1D y, x, z -> 1D x, y, z 
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, NxC 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_cmplx_out1d(i + (j-1)*NxC + (k-1)*NxC*Ny) = data_cmplx_out1dy(j + (i-1)*Ny + (k-1)*Ny*NxC)
                end do 
            end do 
        end do

        ! X inv FFT
#if defined(MFC_OpenACC)
        ierr = cufftExecZ2D(plan_x_bwd_gpu, data_cmplx_out1d, data_real_in1d)
#else
        call fftw_execute_dft_c2r(plan_x_c2r_bwd, data_cmplx_out1d, data_real_in1d)
#endif

        ! 1D x, y, z -> 3D z-slab
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 1, Nx 
            do j = 1, Ny 
                do k = 1, Nzloc
                    data_real_3D_slabz(i, j, k) = data_real_in1d(i + (j-1)*Nx + (k-1)*Nx*Ny)
                end do 
            end do 
        end do

    end subroutine s_mpi_FFT_bwd
 
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

        @:DEALLOCATE(data_real_in1d, data_cmplx_out1d, data_cmplx_out1dy)
        @:DEALLOCATE(cmplx_kernelG1d, real_kernelG_in)
        @:DEALLOCATE(data_real_3D_slabz, data_cmplx_slabz, data_cmplx_slaby)

#if defined(MFC_OpenACC)
        ierr = cufftDestroy(plan_x_fwd_gpu)
        ierr = cufftDestroy(plan_x_bwd_gpu) 
        ierr = cufftDestroy(plan_y_gpu)
        ierr = cufftDestroy(plan_z_gpu)
#else
        call fftw_destroy_plan(plan_x_r2c_fwd)
        call fftw_destroy_plan(plan_x_c2r_bwd)
        call fftw_destroy_plan(plan_y_c2c_fwd) 
        call fftw_destroy_plan(plan_y_c2c_bwd) 
        call fftw_destroy_plan(plan_z_c2c_fwd) 
        call fftw_destroy_plan(plan_z_c2c_bwd) 
        call fftw_destroy_plan(plan_x_r2c_kernelG)
        call fftw_destroy_plan(plan_y_c2c_kernelG)
        call fftw_destroy_plan(plan_z_c2c_kernelG)
#endif

    end subroutine s_finalize_fftw_explicit_filter_module

end module m_fftw
