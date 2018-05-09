subroutine der2y(q1,q2)
  use FFTW3
  implicit none
  DOUBLE PRECISION, DIMENSION(NY) :: q1,q2
  INTEGER :: i

  complex ( kind = 8 ) temp_comp(NY/2 + 1)
  !       FORWARD DFT
  call fftw_execute_dft_r2c(plan_y_f, q1, temp_comp)

  !TAKES THE INTEGRAL OF COSINE AND SINE
  temp_comp(1)  = DCMPLX(0.d0, 0.d0)
  do i = 2, (ny/2)
    temp_comp(i) = -temp_comp(i) * zk2y(i-1)
  end do
  temp_comp(ny/2+1) = DCMPLX(0.d0,0.d0)

  !       BACKWARD DFT
  call fftw_execute_dft_c2r(plan_y_b, temp_comp, q2)

  return
  end
