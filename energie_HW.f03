subroutine energie(hf, Emxy, Eke, Eki, Epe, Er, Etot)
  use GLOBAL
  implicit none

  DOUBLE PRECISION, DIMENSION(NX,NY) :: zmagn,ddx,ddy,hf,psiL
  DOUBLE PRECISION, DIMENSION(NX) :: fx,f1x,k2x

  DOUBLE PRECISION, DIMENSION(NY) :: fy,f1y,k2y
  DOUBLE PRECISION :: Emxy, Eke, Eki, Epe, Er, Etot
  INTEGER :: iy, ix

  !-------------------------------------------------------
  ! Emxy, Energy of the magnetic field in the plane
  !          (B_x + B_y)^2 = (\nabla \psi)^2
  !-------------------------------------------------------
  do iy = 1, ny
    fx(:) = psi(:,iy)
    CALL der1x(fx, f1x)
    zmagn(:,iy) = f1x(:) * f1x(:)
  enddo

  do ix = 1, nx
    fy(:) = psi(ix,:)
    CALL der1y(fy, f1y)
    zmagn(ix,:) = zmagn(ix,:) + f1y(:) * f1y(:)
  enddo

  Emxy = 0.0d0

  do iy = 1, ny
    do ix = 1, nx
      Emxy = zmagn(ix,iy) + Emxy
    end do
  end do

  Emxy = Emxy * dxy


  !-------------------------------------------------------
  ! Eki, Kinetic Energy of Ions
  !        (\nabla \phi)^2
  !-------------------------------------------------------
  do iy = 1, ny
    fx(:) = phi(:,iy)
    CALL der1x(fx, f1x)
    zmagn(:,iy) = f1x(:) * f1x(:)
  enddo

  do ix = 1, nx
    fy(:) = phi(ix,:)
    CALL der1y(fy, f1y)
    zmagn(ix,:) = zmagn(ix,:) + f1y(:) * f1y(:)
  enddo

  Eki = 0.0d0

  do iy = 1, ny
    do ix = 1, nx
      Eki = zmagn(ix,iy) + Eki
    end do
  end do

  Eki = Eki * dxy


  !-------------------------------------------------------
  ! Eke, kinetic energy of the electrons
  ! SQUARE OF CURRENT DENSITY (d_e \nabla^2 \psi)^2 = J^2
  !-------------------------------------------------------

  call laplace(psi,psiL)

  Eke = 0.0d0
  do iy = 1, ny
    do ix = 1, nx
      zmagn(ix,iy) = psiL(ix,iy) * psiL(ix,iy)
      Eke = de*de*zmagn(ix,iy) + Eke
    enddo
  enddo

  Eke = Eke * dxy

  !-------------------------------------------------------
  !                 ENSTROPHY IN Z DIRECTION
  !        rho_s^2 * U^2 = rho_s^2 * (\nabla^2 \phi)^2
  !-------------------------------------------------------

  Epe = 0.0d0

  do iy = 1, ny
    do ix = 1, nx
      zmagn(ix,iy) = rhos*rhos * hf(ix,iy) * hf(ix,iy)
      Epe = zmagn(ix,iy) + Epe
    end do
  end do

  Epe = Epe * dxy

  !-------------------------------------------------------
  !c energia potenziale degli elettroni (correzione dovuta a effetti FLR):
  !c      -  U * phi
  ! N.B.: NEL LIMITE rhoi=0 DIVENTA Eki = (grad phi)^2

  ddx = 0.d0
  ddy = 0.d0

  do ix = 2,nx-1,2
    k2x(ix) = zk2x(ix/2)*nx
    k2x(ix+1) = zk2x(ix/2)*nx
  enddo

  k2x(1) = 0.0d0
  k2x(nx) = 0.0d0

  !------
  do iy = 2,ny-1,2
    k2y(iy) = zk2y(iy/2)*ny
    k2y(iy+1) = zk2y(iy/2)*ny
  enddo

  k2y(1) = 0.0d0
  k2y(ny) = 0.0d0

  !------
  do iy = 1,ny
    fx(:) = phi(:,iy)
    call drfftf(nx, fx, wsavex)
    ddx(:,iy) = fx(:)
  enddo

  do ix = 1,nx
    fy(:) = ddx(ix,:)
    call drfftf(ny, fy, wsavey)
    ddx(ix,:) = fy(:)
  enddo

  do iy = 1,ny
    do ix = 1, nx
      ddx(ix,iy) = ddx(ix,iy)*ddx(ix,iy)
    enddo
  enddo

  do iy = 1, ny
    do ix = 1, nx
      zmagn(ix,iy) = k2x(ix) + k2y(iy)
      zmagn(ix,iy) =  zmagn(ix,iy) &
        /(1.0d0 + rhoi2 * zmagn(ix,iy))
    enddo
  enddo

  !=====================================

  do iy = 2, ny-1, 2
    do ix = 2, nx-1, 2
      ddy(ix,iy) = (ddx(ix,iy) + ddx(ix+1,iy) &
        + ddx(ix,iy+1) + ddx(ix+1,iy+1) )
    enddo
  enddo

  do iy = 2, ny-1, 2
    do ix = 2, nx-1, 2
      ddy(ix,iy) = zmagn(ix,iy) * ddy(ix,iy)
    enddo
  enddo

  Er = 0.0d0

  do iy = 1, ny
    do ix = 1, nx
      Er = ddy(ix,iy)/(nx*ny) + Er
    end do
  end do

  Er = 4.* Er * dxy

  !---------------------------------------------------------
  !  L' energia dissipata dalla resistivita' e':
  !          Ediss = Dres * J^2 =
  !                = Dres/de2 * (Ekxy + Ekz)
  !---------------------------------------------------------

  !TOTAL ENERGY
  Etot = Emxy + Eke + Epe + Er

  return
  end
