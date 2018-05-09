subroutine laplace(q1,q2)
  use GLOBAL
  implicit none

  DOUBLE PRECISION, DIMENSION(NX,NY) :: q1,q2,q1x,q1y
  DOUBLE PRECISION, DIMENSION(NX) :: fx,f2x
  DOUBLE PRECISION, DIMENSION(NY) :: fy,f2y

  INTEGER :: iy, ix

  do iy = 1, ny
    fx(:) = q1(:,iy)
    CALL der2x(fx, f2x)
    q1x(:,iy) = f2x(:)
  end do

  do ix = 1, nx
    fy(:) = q1(ix,:)
    CALL der2y(fy, f2y)
    q1y(ix,:) = f2y(:)          
  end do 

  q2 = q1x + q1y 

  return
  end      
