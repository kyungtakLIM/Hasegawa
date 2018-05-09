subroutine filter(g)
  USE FFTW3
  implicit none

  DOUBLE PRECISION, DIMENSION(NX,NY) :: g
  DOUBLE PRECISION, DIMENSION(NX) :: wkx
  DOUBLE PRECISION, DIMENSION(NY) :: wky
  DOUBLE PRECISION :: A,B
  INTEGER :: j, i
  complex ( kind = 8 ) t_c_x(NX/2 + 1)
  complex ( kind = 8 ) t_c_y(NY/2 + 1)

  ! X-Filtering
  do j = 1, ny
    wkx(:) = g(:,j)

    call fftw_execute_dft_r2c(plan_x_f, wkx, t_c_x)
    do i = 2, (nx/2)
      B = IMAGPART(t_c_x(i)) *work3(i-1)
      A = REALPART(t_c_x(i)) *work3(i-1)
      t_c_x(i) = DCMPLX(A,B)
    enddo
    t_c_x(nx/2+1) = DCMPLX(0.d0,0.d0)
    t_c_x = t_c_x / nx

    call fftw_execute_dft_c2r(plan_x_b, t_c_x, wkx)

    g(:,j) = wkx(:)
  enddo

  ! Y-Filtering
  do i = 1, nx
    wky(:) = g(i,:)

    call fftw_execute_dft_r2c(plan_y_f, wky, t_c_y)
    do j = 2, (ny/2)
      B = IMAGPART(t_c_y(j)) *work4(j-1)
      A = REALPART(t_c_y(j)) *work4(j-1)
      t_c_y(j) = DCMPLX(A,B)
    enddo
    t_c_y(ny/2+1) = DCMPLX(0.d0,0.d0)
    t_c_y = t_c_y / ny

    call fftw_execute_dft_c2r(plan_y_b, t_c_y, wky)

    g(i,:) = wky(:)
  enddo


  return
  end 
