subroutine helm1(hp) 

  !       based on version (FLR PPCF march 2011)
  !       calculates psi from F through the F=psi-d_e^2nabla^2psi
  !       inverted in the fourier space
  !       input:  hp (F)
  !       output: psi (via par.inc)
  use GLOBAL
  implicit none

  INTEGER :: ix, iy, m, n
  DOUBLE PRECISION :: aky, akk, akx, elambda 
  DOUBLE PRECISION, DIMENSION(NX) :: f1 
  DOUBLE PRECISION, DIMENSION(NY) :: f2
  DOUBLE PRECISION, DIMENSION(NX,NY) :: d, dd, hp 

  IF (de < 1.0d-10) THEN
    psi=hp
  ELSE

    elambda = -1.d0/de2

    dd = - hp/de2

    !dd fourier transfo in x and y of F/(de^2)

    do iy = 1, ny
      f1(:) = dd(:,iy)
      call drfftf(nx, f1, wsavex)
      dd(:,iy) = f1(:)
    enddo

    do ix = 1, nx
      f2(:) = dd(ix,:)
      call drfftf(ny, f2, wsavey)
      dd(ix,:) = f2(:)
    enddo           

    !          formula inversion: dd becomes the fourier transfo of psi

    !	   m = 1, n = 1

    akk = - elambda
    if (abs(akk) < 1.0d-14) then
      d(1,1) = 0.0
    else 
      akk = akk * nx * ny
      d(1,1) = - dd(1,1) / akk
    endif 

    !	   m = 1

    do n = 2, ny-1, 2
      aky = grady * (n/2)
      akk = aky * aky - elambda
      akk = akk * nx * ny
      d(1,n)     = - dd(1,n)   / akk
      d(1,n+1)   = - dd(1,n+1) / akk
    enddo

    !	   n = 1

    do m = 2, nx-1, 2
      akx = gradx * (m/2)
      akk = akx * akx - elambda
      akk = akk * nx * ny
      d(m,1)     = - dd(m,1)   / akk
      d(m+1,1)   = - dd(m+1,1) / akk
    enddo

    do n = 2, ny-1, 2
      aky = grady * (n/2)
      do m = 2, nx-1, 2
        akx = gradx * (m/2)
        akk = akx * akx + aky * aky - elambda
        akk = akk * nx * ny
        d(m,n)     = - dd(m,n)     / akk
        d(m,n+1)   = - dd(m,n+1)   / akk
        d(m+1,n)   = - dd(m+1,n)   / akk
        d(m+1,n+1) = - dd(m+1,n+1) / akk
      enddo
    enddo 

    d(:,ny) = 0.0
    d(nx,:) = 0.0

    !Fourier transfo inverse: dd becomes the result psi 

    do iy = 1, ny
      f1(:) = d(:,iy)
      call drfftb(nx, f1, wsavex)  
      d(:,iy) = f1(:)
    enddo

    do ix = 1, nx
      f2(:) = d(ix,:)
      call drfftb(ny, f2, wsavey)
      d(ix,:) = f2(:)
    enddo  

    psi = d 

  ENDIF 

  return        
  end  

