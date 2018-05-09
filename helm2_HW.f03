subroutine helm2(hf) 

  !  risolve Helmholtz: nabla^2 phi  = g2 (psi output)
  !  g2 = - hf
  use FFTW3
  implicit none

  DOUBLE PRECISION, DIMENSION(NX) :: f1
  DOUBLE PRECISION, DIMENSION(NY) :: f2
  DOUBLE PRECISION, DIMENSION(NX,NY) :: d,dd,hf
  DOUBLE PRECISION :: akx, aky, akk
  INTEGER :: n, m, ix, iy

  dd =  hf

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

  !	m = 1, n = 1
  d(1,1) = 0.0

  !	m = 1
  do n = 2, ny-1, 2
    aky = grady * n/2
    akk = aky * aky * nx * ny
    d(1,n)     = - dd(1,n)   / akk
    d(1,n+1)   = - dd(1,n+1) / akk
  enddo

  !	n = 1

  do m = 2, nx-1, 2
    akx = gradx * m/2
    akk = akx * akx * nx * ny
    d(m,1)     = - dd(m,1)   / akk
    d(m+1,1)   = - dd(m+1,1) / akk
  enddo

  do n = 2, ny-1, 2
    aky = grady * n/2
    do m = 2, nx-1, 2
      akx = gradx * m/2
      akk = (akx * akx + aky * aky) * nx * ny
      d(m,n)     = - dd(m,n)     / akk
      d(m,n+1)   = - dd(m,n+1)   / akk
      d(m+1,n)   = - dd(m+1,n)   / akk
      d(m+1,n+1) = - dd(m+1,n+1) / akk
    enddo
  enddo 

  d(:,ny) = 0.0
  d(nx,:) = 0.0

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

  phi = d 

  return        
  end  


