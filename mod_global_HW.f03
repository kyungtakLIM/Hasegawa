MODULE GLOBAL

  IMPLICIT NONE

  INTEGER, PARAMETER :: nx=512, ny=512

  common/boite1/dx,dy,dxy
  common/boitehalf/gradx,grady
  common/boite2/time
  common/boite3/de,de2,rhos,rhos2,rhoi,rhoi2
  common/Fields1/x(nx),y(ny)
  common/Fields2/work3,work4
  common/Fields3/phi(nx,ny),psi(nx,ny)
  common/Fields4/wsavex(2*nx+15), wsavey(2*ny+15)
  common/Fourier1/zk1x(nx/2-1),zk1y(ny/2-1)
  common/Fourier2/zk2x(nx/2-1),zk2y(ny/2-1)

  DOUBLE PRECISION :: x, y, phi, psi
  DOUBLE PRECISION :: wsavex, wsavey
  DOUBLE PRECISION, DIMENSION(NX/2) :: work3
  DOUBLE PRECISION, DIMENSION(NY/2) :: work4
  DOUBLE PRECISION :: zk1x, zk1y, zk2x, zk2y
  DOUBLE PRECISION de,de2,rhos,rhos2,rhoi,rhoi2
  DOUBLE PRECISION time, grady, gradx
  DOUBLE PRECISION dx,dy,dxy

  !CONTAINS
  !  SUBROUTINE common_thing
  !  END SUBROUTINE

  !  SUBROUTINE common_thing2
  !  END SUBROUTINE
  END MODULE
