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
 s_initialize_fftw_explicit_filter, &
 s_apply_fftw_explicit_filter, & 
 s_initialize_fftw_kernelG

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
    integer :: ierr

    integer, allocatable :: gpu_fft_size(:), iembed(:), oembed(:)

    integer :: istride, ostride, idist, odist, rank
#endif
    ! new
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
    subroutine s_initialize_fftw_explicit_filter
        print *, 'FFTW SETUP...'

        ! data setup 
        p_real = fftw_alloc_real(int((m+1)*(n+1)*(p+1), C_SIZE_T))
        p_complex = fftw_alloc_complex(int(((m+1)/2+1)*(n+1)*(p+1), C_SIZE_T))

        call c_f_pointer(p_real, array_in, [m+1, n+1, p+1])
        call c_f_pointer(p_complex, array_out, [(m+1)/2+1, n+1, p+1])

        plan_forward = fftw_plan_dft_r2c_3d(p+1, n+1, m+1, array_in, array_out, FFTW_ESTIMATE)
        plan_backward = fftw_plan_dft_c2r_3d(p+1, n+1, m+1, array_out, array_in, FFTW_ESTIMATE)

        ! kernel setup
        p_kernelG_real = fftw_alloc_real(int((m+1)*(n+1)*(p+1), C_SIZE_T))
        p_kernelG_complex = fftw_alloc_complex(int(((m+1)/2+1)*(n+1)*(p+1), C_SIZE_T))

        call c_f_pointer(p_kernelG_real, array_kernelG_in, [m+1, n+1, p+1])
        call c_f_pointer(p_kernelG_complex, array_kernelG_out, [(m+1)/2+1, n+1, p+1])

        plan_kernelG_forward = fftw_plan_dft_r2c_3d(p+1, n+1, m+1, array_kernelG_in, array_kernelG_out, FFTW_ESTIMATE)

    end subroutine s_initialize_fftw_explicit_filter

    subroutine s_initialize_fftw_kernelG

        real(dp) :: r, delta
        integer :: i, j, k

        delta = 4*patch_ib(1)%radius

        do i = 1, m+1
            do j = 1, n+1 
                do k = 1, p+1
                    r = sqrt(x_cc(i)**2 + y_cc(i)**2 + z_cc(i)**2)
                    if (r >= 0 .and. r <= delta) then
                        array_kernelG_in(i, j, k) = 21_dp/(2_dp*pi*delta**3)*(4*r/delta + 1)*(1 - r/delta)**4
                    else 
                        array_kernelG_in(i, j, k) = 0 
                    end if 
                end do
            end do
        end do
        print *, array_kernelG_in(int(m/2), int(n/2), int(p/2))

    end subroutine s_initialize_fftw_kernelG

    subroutine s_apply_fftw_explicit_filter(q_cons_vf, q_filtered, volfrac_phi)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_filtered

        real(dp) :: volfrac_phi

        integer :: i, j, k, l

        do l = 1, sys_size
            do i = 0, m
                do j = 0, n 
                    do k = 0, p
                        if (ib_markers%sf(i, j, k) == 0) then
                            array_in(i+1, j+1, k+1) = q_cons_vf(l)%sf(i, j, k)
                        else
                            array_in(i+1, j+1, k+1) = 0_dp
                        end if
                    end do 
                end do
            end do

            call fftw_execute_dft_r2c(plan_forward, array_in, array_out)
            call fftw_execute_dft_r2c(plan_kernelG_forward, array_kernelG_in, array_kernelG_out)

            array_out(:, :, :) = array_out(:, :, :) * array_kernelG_out(:, :, :)

            call fftw_execute_dft_c2r(plan_backward, array_out, array_in)

            q_filtered(l)%sf(0:m, 0:n, 0:p) = array_in(1:m+1, 1:n+1, 1:p+1) / ((m+1)*(n+1)*(p+1) * (1._wp - volfrac_phi)) ! unnormalized DFT

            !print *, q_cons_vf(l)%sf(10, 10, 10)
            !print *, q_filtered(l)%sf(10, 10, 10)
        end do 


    end subroutine s_apply_fftw_explicit_filter

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

end module m_fftw
