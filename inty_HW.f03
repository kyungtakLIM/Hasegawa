subroutine inty(q1, q2, c0y)
  use FFTW3
  implicit none

  !       integral function on Y
  DOUBLE PRECISION, DIMENSION(NY) :: q1,q2
  DOUBLE PRECISION :: c0y, A, B, NY2
  INTEGER :: i
  complex ( kind = 8 ) temp_comp(NY/2 + 1)
  !       FORWARD DFT
  call fftw_execute_dft_r2c(plan_y_f, q1, temp_comp)

  c0y = REALPART(temp_comp(1)) / NY
  !TAKES THE INTEGRAL OF COSINE AND SINE
  temp_comp(1)  = DCMPLX(0.d0, 0.d0)
  NY2 = NY * NY
  do i = 2, (NY/2)
    A = IMAGPART(temp_comp(i)) / (zk1y(i-1) * NY2)
    B = -REALPART(temp_comp(i)) / (zk1y(i-1) * NY2)
    temp_comp(i) = DCMPLX(A,B)
  end do
  temp_comp(NY/2+1) = DCMPLX(0.d0,0.d0)

  !       BACKWARD DFT
  call fftw_execute_dft_c2r(plan_y_b, temp_comp, q2)

  return
  end
