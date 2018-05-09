subroutine flux(g, Fx, gx, Fy, gy)
  use FFTW3
  implicit none

  ! gx, gy are calculated flux
  
  DOUBLE PRECISION, DIMENSION(NX,NY) :: DD, g, Fx, gx, Fy, gy, fmx, fmy
  DOUBLE PRECISION, DIMENSION(NY) :: ptx, qtx
  DOUBLE PRECISION, DIMENSION(NX) :: pty, qty
  DOUBLE PRECISION :: c0y, c0x
  INTEGER :: i,j

  ! COMPUTE X-Flux
  !     compute DD = G * Fx
  ! And remember, for efficency in Fortran
  ! the first index varies fastest

  fmx = g * Fx
  fmy = g * Fy
  do j = 1, ny
    do i = 1, nx - 1
      DD(i,j) = (fmx(i+1,j) - fmx(i,j)) / dxy
    enddo
  enddo
  DD(nx,:) = (fmx(1,:) - fmx(nx,:)) / dxy

  do i = 1, nx
    ptx(:) = DD(i,:)
    call inty(ptx, qtx, c0y)
    do j = 1, ny-1
      gx(i,j) = qtx(j+1) - qtx(j) + c0y * dy
    end do
    gx(i,ny) = qtx(1) - qtx(ny) + c0y * dy
  end do

  ! COMPUTE Y-FLUX
  !     compute DD = G * Fx
  do j = 1, ny-1
    do i = 1, nx
      DD(i,j) = (fmy(i,j+1) - fmy(i,j)) / dxy
    enddo
  enddo
  DD(:,ny) = (fmy(:,1) - fmy(:,ny)) / dxy

  do j = 1, ny
    pty(:) = DD(:,j)
    call intx(pty, qty, c0x)
    do i = 1, nx-1
      gy(i,j) = qty(i+1) - qty(i) + c0x * dx
    end do
    gy(nx,j)= qty(1) - qty(nx) + c0x * dx
  end do

  return
  end
