subroutine rcstr(gmean,g)
  use FFTW3
  implicit none
  !  include 'par.inc'

  DOUBLE PRECISION, DIMENSION(NX, NY) :: gmean, g, g_3
  DOUBLE PRECISION, DIMENSION(NX) :: g_1, g_2
  DOUBLE PRECISION, DIMENSION(NY) :: g_4, g_5
  DOUBLE PRECISION :: zi, A, B
  INTEGER :: i,j

  complex ( kind = 8 ) t_c_x(NX/2 + 1)
  complex ( kind = 8 ) t_c_y(NY/2 + 1)

  do j = 1,ny
    g_1(:) = gmean(:,j)

    call fftw_execute_dft_r2c(plan_x_f, g_1, t_c_x)
    t_c_x(1)  = DCMPLX(REALPART(t_c_x(1)/nx), 0.d0)
    do i = 2,nx/2
      zi = (i-1)*gradx*dx/2.0d0
      A = ( REALPART(t_c_x(i))*cos(zi) + IMAGPART(t_c_x(i))*sin(zi)) &
        *zi/sin(zi)/nx
      B = (-REALPART(t_c_x(i))*sin(zi) + IMAGPART(t_c_x(i))*cos(zi)) &
        *zi/sin(zi)/nx
      t_c_x(i) = DCMPLX(A,B)
    enddo
    t_c_x(nx/2+1) = DCMPLX(0.d0,0.d0)
    call fftw_execute_dft_c2r(plan_x_b, t_c_x,g_2)
    g_3(:,j) = g_2(:)
  enddo

  do i = 1,nx
    g_4(:) = g_3(i,:)

    call fftw_execute_dft_r2c(plan_y_f, g_4, t_c_y)
    t_c_y(1)  = DCMPLX(REALPART(t_c_y(1)/ny), 0.d0)
    do j = 2,ny/2
      zi = (j-1)*grady*dy/2.0d0
      A = ( REALPART(t_c_y(j))*cos(zi) + IMAGPART(t_c_y(j))*sin(zi)) &
        *zi/sin(zi)/ny
      B = (-REALPART(t_c_y(j))*sin(zi) + IMAGPART(t_c_y(j))*cos(zi)) &
        *zi/sin(zi)/ny
      t_c_y(j) = DCMPLX(A,B)
    enddo
    t_c_y(ny/2+1) = DCMPLX(0.d0,0.d0)
    call fftw_execute_dft_c2r(plan_y_b,t_c_y,g_5)
    g(i,:) = g_5(:)
  enddo

  return
  end
