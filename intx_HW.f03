subroutine intx(q1, q2, c0x)
  use FFTW3
  implicit none

  !       integral function on X
  DOUBLE PRECISION, DIMENSION(NX) :: q1,q2
  DOUBLE PRECISION :: c0x, A, B, NX2
  INTEGER :: i
  complex ( kind = 8 ) temp_comp(NX/2 + 1)
  !       FORWARD DFT
  call fftw_execute_dft_r2c(plan_x_f, q1, temp_comp)

  c0x = REALPART(temp_comp(1)) / nx
  !TAKES THE INTEGRAL OF COSINE AND SINE
  temp_comp(1)  = DCMPLX(0.d0, 0.d0)
  NX2 = NX* NX
  do i = 2, (NX/2)
    A = IMAGPART(temp_comp(i)) / (zk1x(i-1) * NX2)
    B = -REALPART(temp_comp(i)) / (zk1x(i-1) * NX2)
    temp_comp(i) = DCMPLX(A,B)
  end do
  temp_comp(nx/2+1) = DCMPLX(0.d0,0.d0)

  !       BACKWARD DFT
  call fftw_execute_dft_c2r(plan_x_b, temp_comp, q2)

  return
  end
