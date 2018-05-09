subroutine outfield(field1, file_num)
  use GLOBAL
  implicit none

  DOUBLE PRECISION, DIMENSION(NX,NY) :: field1
  INTEGER :: ix, iy, file_num

  write(file_num) time
  write(file_num) ((field1(ix,iy), ix = 1,nx), iy = 1,ny)

  ! ordering of numbering
  ! (ix=1, iy=1), (ix=2, iy=1) .... (ix=nx, iy=1)
  ! (ix=1, iy=2), (ix=2, iy=2) .... (ix=nx, iy=2)
  ! ...
  ! ...
  ! (ix=1, iy=ny), (ix=2, iy=ny) .... (ix=nx, iy=ny)

  return
  end
  
