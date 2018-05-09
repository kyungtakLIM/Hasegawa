subroutine grid_init(xl,yl,omega)
  use GLOBAL
  !  This routine defines the grid
  implicit none

  DOUBLE PRECISION :: w, xl, yl, work1, work2, zeta
  DOUBLE PRECISION :: alpha, beta, delta, gamma, omega
  DOUBLE PRECISION :: pi

  DOUBLE PRECISION, DIMENSION(3*NX/2) :: trigsx
  DOUBLE PRECISION, DIMENSION(3*NY/2) :: trigsy
  INTEGER, DIMENSION(10):: ifaxx,ifaxy
  INTEGER :: ix, kk, iy, i

  dx = xl / nx
  dy = yl / ny
  dxy = dx * dy

  pi = dACOS(-1.0d0)

  if(nx.LE.2) dx = 1.0d0
  if(ny.LE.2) dy = 1.0d0

  do ix = 1, nx
    x(ix) = - xl / 2.0 + (ix-1) * dx
  enddo

  do iy = 1, ny
    y(iy) = - yl / 2.0 + (iy-1) * dy
  enddo

  do i = 1, nx/2-1
    zk1x(i) = gradx * i
    zk2x(i) = zk1x(i) * zk1x(i) / nx
    zk1x(i) = zk1x(i) / nx
  end do

  do i = 1, ny/2-1
    zk1y(i) = grady * i
    zk2y(i) = zk1y(i) * zk1y(i) / ny
    zk1y(i) = zk1y(i) / ny
  end do

  !-------------------------------------------------------------

  write(*,*) 'x_0, y_0, x_n, y_n'
  write(*,'(F9.6,F9.6,F9.6,F9.6)') x(1), y(1), x(nx), y(ny)

  !  de = electron skin depth
  !  rhos = sqrt(T_e/m_i)/omega_{ci}
  !  (normalized to Lx = pi = xl/2.)
  de2 = de * de
  rhos2 = rhos * rhos
  rhoi2 = rhoi * rhoi

  !   new FFT parameters
  call drffti(nx, wsavex)
  call drffti(ny, wsavey)

  !   FFT parameters
  CALL FFTRIG(TRIGSX,nx,3)
  CALL FAX(IFAXX,nx,3)
  CALL FFTRIG(TRIGSY,ny,3)
  CALL FAX(IFAXY,ny,3)

  !   filter parameters
  zeta  = (3.d0 - 2.d0*omega)/10.d0
  alpha = (2.d0 + 3.d0*omega)/4.d0
  beta  = (6.d0 + 7.d0*omega)/8.d0
  gamma = (6.d0 + omega)/20.d0
  delta = (2.d0 - 3.d0*omega)/40.d0

  write(*,*) 'alpha,beta,gamma,delta: '
  write(*,'(F9.6,F9.6,F9.6,F9.6)')  alpha, beta, gamma, delta
  write(*,*) 'omega, zeta: '
  write(*,'(F9.6,F9.6)') omega, zeta

  do kk = 1, nx/2
    w = 2.d0 * pi * gradx * kk / nx

    work1 = alpha + beta*dcos(w) + &
      gamma*dcos(2.d0*w) + delta*dcos(3.d0*w)

    work2 = 1.d0+2.d0*omega*dcos(w)+2.d0*zeta*dcos(2.d0*w)

    work3(kk) = work1/work2
  enddo


  do kk = 1, ny/2
    w = 2.d0 * pi * grady * kk / ny

    work1 = alpha + beta*dcos(w) + &
      gamma*dcos(2.d0*w) + delta*dcos(3.d0*w)

    work2 = 1.d0+2.d0*omega*dcos(w)+2.d0*zeta*dcos(2.d0*w)

    work4(kk) = work1/work2
  enddo

  return
  end
