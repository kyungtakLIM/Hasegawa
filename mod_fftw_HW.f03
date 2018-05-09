MODULE FFTW3

  use, intrinsic :: iso_c_binding
  use GLOBAL
  IMPLICIT NONE
  include 'fftw3.f03'

  type(C_PTR) :: plan_x_f !FFT PLAN
  type(C_PTR) :: plan_x_b !IFFT PLAN
  type(C_PTR) :: plan_y_f !FFT PLAN
  type(C_PTR) :: plan_y_b !IFFT PLAN
  type(C_PTR) :: plan_2d_f !FFT PLAN
  type(C_PTR) :: plan_2d_b !IFFT PLAN

  DOUBLE PRECISION, DIMENSION(NX) :: in_x_f, out_x_b
  COMPLEX(kind = 8), DIMENSION (NX/2 + 1) :: out_x_f, in_x_b

  DOUBLE PRECISION, DIMENSION(NY) :: in_y_f, out_y_b
  COMPLEX(kind = 8), DIMENSION (NY/2 + 1) :: out_y_f, in_y_b

  DOUBLE PRECISION, DIMENSION(NX,NY) :: in_2d_f, out_2d_b
  COMPLEX(kind = 8), DIMENSION (NX,NY/2 + 1) :: out_2d_f, in_2d_b
  !  TYPE :: fftw3Needs
  !    type(C_PTR) :: plan_f, plan_b
  ! input of FFT and output (result) of inverse FFT
  !    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: in_f, out_b
  ! output (result) of FFT and input of inverse FFT
  !    COMPLEX(kind = 8), ALLOCATABLE, DIMENSION (:) :: out_f, in_b
  !  END TYPE
  !  TYPE (fftw3Needs) :: fftw3

CONTAINS
  SUBROUTINE init_fft
  plan_x_f  = fftw_plan_dft_r2c_1d(NX, in_x_f, out_x_f, FFTW_ESTIMATE)
  plan_x_b  = fftw_plan_dft_c2r_1d(NX, in_x_b, out_x_b, FFTW_ESTIMATE)
  plan_y_f  = fftw_plan_dft_r2c_1d(NY, in_y_f, out_y_f, FFTW_ESTIMATE)
  plan_y_b  = fftw_plan_dft_c2r_1d(NY, in_y_b, out_y_b, FFTW_ESTIMATE)
  plan_2d_f = fftw_plan_dft_r2c_2d(NX,NY,in_2d_f,out_2D_f,FFTW_ESTIMATE)
  plan_2d_b = fftw_plan_dft_c2r_2d(NX,NY,in_2d_b,out_2D_b,FFTW_ESTIMATE)
  END SUBROUTINE

  SUBROUTINE destroy_fft
    CALL dfftw_destroy_plan(plan_x_f)
    CALL dfftw_destroy_plan(plan_x_b)
    CALL dfftw_destroy_plan(plan_y_f)
    CALL dfftw_destroy_plan(plan_y_b)
    CALL dfftw_destroy_plan(plan_2d_f)
    CALL dfftw_destroy_plan(plan_2d_b)
  END SUBROUTINE


END MODULE
