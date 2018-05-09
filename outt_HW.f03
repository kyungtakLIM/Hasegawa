subroutine outt(Emxy,Eke,Eki,Epe,Er,Etot,F_sum,U_sum)
  use GLOBAL
  implicit none

  DOUBLE PRECISION :: Emxy, Eke, Eki, Epe, Er, Etot
  DOUBLE PRECISION :: F_sum, U_sum

  !       unit=17 --> 'Energy.dat'
  write(17) time, Emxy, Eke, Eki, Epe, Er, Etot

  return
  end
