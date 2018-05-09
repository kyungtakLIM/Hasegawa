subroutine der2x(q1,q2)
  use FFTW3
  implicit none
  DOUBLE PRECISION, DIMENSION(NX) :: q1,q2
  INTEGER :: i

  complex ( kind = 8 ) temp_comp(NX/2 + 1)
  !       FORWARD DFT
  call fftw_execute_dft_r2c(plan_x_f, q1, temp_comp)

  !TAKES THE INTEGRAL OF COSINE AND SINE
  temp_comp(1)  = DCMPLX(0.d0, 0.d0)
  do i = 2, (nx/2)
    temp_comp(i) = -temp_comp(i) * zk2x(i-1)
  end do
  temp_comp(nx/2+1) = DCMPLX(0.d0,0.d0)

  !       BACKWARD DFT
  call fftw_execute_dft_c2r(plan_x_b, temp_comp, q2)

  return
  end

