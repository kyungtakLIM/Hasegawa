subroutine cond_init(gp,gf,hp,hf,x0a, y0a, x0b, y0b, sigma_a, sigma_b, pso)
  use GLOBAL
  implicit none

  DOUBLE PRECISION :: x0a, y0a, x0b, y0b, radius_b2, radius_a2
  DOUBLE PRECISION :: sigma_a, sigma_b, RRa, Rrb, pso
  ! DOUBLE PRECISION ::  m,n,n1
  ! OUBLE PRECISION :: preangle_a2, preangle_b1, preangle_b2
  ! DOUBLE PRECISION :: shape_a, shape_b, preangle_a1
  ! DOUBLE PRECISION :: theta_a, theta_b


  DOUBLE PRECISION, DIMENSION(NX,NY) :: gp,gf,hp,hf
  DOUBLE PRECISION, DIMENSION(NX,NY) :: gpsi,gphi,gpsiL
  INTEGER :: iy, ix

  ! Calling random array for initial condition of phi, psi

  do iy = 1, ny
    do ix = 1, nx
      radius_a2 =  (x(ix)-x0a)**2 + (y(iy)-y0a)**2
      radius_b2 =  (x(ix)-x0b)**2 + (y(iy)-y0b)**2
      RRa = radius_a2 / sigma_a
      RRb = radius_b2 / sigma_b
      gpsi(ix,iy) = exp(-RRa)  +  exp(-RRb)
      !gphi(ix,iy) = exp(-RRa)  +  exp(-RRb)
      gphi(ix,iy) = 0
    end do
  end do

  call rcstr(gpsi,psi)
  call rcstr(gphi,phi)
  call laplace(gpsi,gpsiL)


  do iy = 1, ny
    do ix = 1, nx
      gp(ix,iy) = gpsi(ix,iy) - de*de * gpsiL(ix,iy)
      gf(ix,iy) = 0.d0
    end do
  end do

  call rcstr(gp,hp)
  call rcstr(gf,hf)

  return
  end
