subroutine conserv(gp,gf, F_sum, U_sum)
  use GLOBAL
  implicit none

  DOUBLE PRECISION :: F_sum, U_sum
  DOUBLE PRECISION, DIMENSION(NX,NY) :: gp,gf
  INTEGER :: ix, iy


  F_sum = sum(gp)*dxy
  U_sum = sum(gf)*dxy

  return
  end
