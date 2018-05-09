subroutine der(Fxp, Fyp, Fxf, Fyf)
  use GLOBAL
  implicit none
  !   derivo (rispetto a x) psi e phi

  DOUBLE PRECISION, DIMENSION(NX,NY) :: Fxp,Fyp,Fxf,Fyf
  DOUBLE PRECISION, DIMENSION(NX) :: d1,d2,f1,f2
  DOUBLE PRECISION, DIMENSION(NY) :: d3,d4,f3,f4
  INTEGER :: iy, ix

  do iy = 1,ny
    f1(:) = psi(:,iy)
    f2(:) = phi(:,iy)
    CALL der1x(f1,d1)
    CALL der1x(f2,d2)
    Fyp(:,iy) = d1(:)
    Fyf(:,iy) = d2(:)
  enddo

  ! Fyp = dérivé de psi en y, Fyf = dérivé de phi en y
  !   derivo (rispetto a y) psi e phi

  do ix = 1,nx
    f3(:) = psi(ix,:)
    f4(:) = phi(ix,:)
    CALL der1y(f3,d3)
    CALL der1y(f4,d4)
    Fxp(ix,:) = - d3(:)
    Fxf(ix,:) = - d4(:)
  enddo

  ! Fxp = dérivé de psi en x, Fxf = dérivé de phi en x

  return
  end
