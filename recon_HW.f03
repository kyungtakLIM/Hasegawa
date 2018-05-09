!MAIN FILE FOR DELSARTOS 2D CODE
program recon
  USE FFTW3 ! Originally USE FFTW3
  implicit none

  DOUBLE PRECISION, DIMENSION(NX,NY) :: flxf1,flxf2,flxf3
  DOUBLE PRECISION, DIMENSION(NX,NY) :: flxp1,flxp2,flxp3
  DOUBLE PRECISION, DIMENSION(NX,NY) :: gxp, gyp, gxf, gyf
  DOUBLE PRECISION, DIMENSION(NX,NY) :: gxf2, gyf2, gp, gf, hp, hf
  DOUBLE PRECISION, DIMENSION(NX,NY) :: Fxp, Fyp, Fxf, Fyf
  DOUBLE PRECISION, DIMENSION(NX,NY) :: psiL, gxrh, gyrh
  DOUBLE PRECISION ab1,ab2,ab3
  INTEGER :: it,ir ! Counters
  REAL :: start, finish
  DOUBLE PRECISION :: ioutt, ioutf, pi

  ! Quantities to print at each noutt
  DOUBLE PRECISION :: Emxy, Eke, Eki, Epe, Er, Etot, Hc
  DOUBLE PRECISION :: U_sum, F_sum

  ! Initial Parameters
  INTEGER :: nstep, noutt, noutf
  DOUBLE PRECISION :: xl, yl, dt, omega, m,n,n1,pso
  DOUBLE PRECISION :: x0a, y0a, x0b, y0b
  DOUBLE PRECISION :: sigma_a, sigma_b


  !Read in Parameters from command line or in.com file
  call read_arg(nstep, noutt, noutf, xl, yl, dt, omega, x0a, y0a, x0b, y0b, sigma_a, sigma_b,m,n,n1,pso)

  open(unit=12,status='unknown',file='Field.dat', &
    ACCESS='STREAM')
  open(unit=13,status='unknown',file='G.dat', &
    ACCESS='STREAM')
  open(unit=17,status='unknown',file='Energy.dat', &
    ACCESS='STREAM')
  ! open(unit=18,status='unknown',file='Helicity.dat', &
  !   ACCESS='STREAM') \ not necessary in HW-type equation

  !The reasoning here is legit, don't change.
  pi = dACOS(-1.0d0)

  !Set up spatial dimensions
  gradx = 1.0d0 / xl
  grady = 1.0d0 / yl

  xl = 2.0d0 * pi * xl
  yl = 2.0d0 * pi * yl

  !Set up timers
  time = 0.
  ioutt = 1
  ioutf = 1

  !ADAMS-BASHFORTH COEFFICIENTS
  ab1   = dt * 5.d0 / 12.d0
  ab2   = - dt * 4.d0 / 3.d0
  ab3   = dt * 23.d0 / 12.d0

  !PREPARE THE FFT PLANS
  write(*,*) ''
  write(*,*) 'INIT CALLED...'
  call init_fft
  call grid_init(xl,yl,omega)
  call cond_init(gp, gf, hp, hf, x0a, y0a, x0b, y0b, sigma_a, sigma_b, pso)

  write(*,*) ''
  write(*,*) 'SIMULATION PARAMETERS'
  write(*,*) 'de, rhos, rhoi:'
  write(*,'(F9.6,F9.6,F9.6)') de, rhos, rhoi
  write(*,*) 'nstep, noutt, noutf'
  write(*,*) nstep, noutt, noutf
  write(*,*) 'Nx, Ny:'
  write(*,*) nx, ny
  write(*,*) 'dt, dx, dy:'
  write(*,'(F9.6,F9.6,F9.6)') dt, dx, dy
  write(*,*) 'dxy'
  write(*,'(F9.7)') dxy

  call conserv(gp,gf, F_sum, U_sum)
  call energie(hf, Emxy, Eke, Eki, Epe, Er, Etot)
  ! call helicity(Hc) \ not necessary in HW-type equation

  call outt(Emxy,Eke, Eki, Epe, Er, Etot,Hc, F_sum, U_sum)
  !Write Outfield Header
  !write(12) time, nx, ny, xl, yl, nstep, noutf
  ! LIM : I erased Outfield Header on the purposer of simplicity of result.
  call outfield(phi,12)
  call outfield(hf, 13)

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Naming  Ff and Fpsi the divergence-free vectors of components
  ! Fphi = (Fxf, Fyf)  and  Fpsi = (Fxp, Fyp)
  ! and by using  (Fphi\cdot\nabla) g = div(Fphi g),
  ! we see that once integrating over a grid cell between x_i, x_i+i
  ! and y_j, y_j+1, the divergence term becomes a flux term.
  ! The contributions deriving from these flux integrals are
  ! evaluated by means of the subroutine 'flux' whereas
  ! the components of Ff and Fpsi are evaluated by 'der'
  !
  ! So we evaluate the rhs terms of the time
  ! evolution of the average values (i.e. integrated over a cell
  ! volume) of F (i.e. hp) and U (i.e. hf).
  !
  ! The subroutine 'rcnst' evaluates F and U on the grid points
  ! once the 'average values' of F and U are known after the time
  ! advancement.
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  !Make 2 normal Euler steps to seed for the Adams-Bashforth Algorithm
  time = dt
  call laplace(psi,psiL)
  call der(Fxp, Fyp, Fxf, Fyf)
  call flux(hp,Fxf,gxp, Fyf, gyp)
  call flux(hf,Fxf,gxf, Fyf, gyf)
  call flux(hf,Fxp,gxrh, Fyp, gyrh)
  call flux(psiL,Fxp,gxf2, Fyp, gyf2)

  flxp1 = - gxp - gyp - rhos2 * (gxrh + gyrh)
  flxf1 = - gxf - gyf + gxf2 + gyf2
  gp = gp + flxp1 * dt
  gf = gf + flxf1 * dt

  call filter(gp)
  call filter(gf)
  call rcstr(gp,hp)
  call rcstr(gf,hf)
  call helm1(hp)
  call helm2(hf)
  call laplace(psi,psiL)

  !------------------------------------------------------
  time = dt + time
  call der(Fxp, Fyp, Fxf, Fyf)
  call flux(hp,Fxf,gxp, Fyf, gyp)
  call flux(hf,Fxf,gxf, Fyf, gyf)
  call flux(hf,Fxp,gxrh, Fyp, gyrh)
  call flux(psiL,Fxp,gxf2, Fyp, gyf2)

  flxp2 = - gxp - gyp - rhos2 * (gxrh + gyrh)
  flxf2 = - gxf - gyf + gxf2 + gyf2
  gp = gp + flxp2 * dt
  gf = gf + flxf2 * dt

  call filter(gp)
  call filter(gf)
  call rcstr(gp,hp)
  call rcstr(gf,hf)
  call helm1(hp)
  call helm2(hf)
  call laplace(psi,psiL)

  !-------------------------- ----------------------------
  !      Begin Adams-Bashforth Time-Step Loop
  !-----------------------------------------------------
  CALL CPU_TIME(start)
  ir = 1

  do it = 1, nstep
    time = time + dt

    call der(Fxp, Fyp, Fxf, Fyf)
    call flux(hp,Fxf,gxp, Fyf, gyp)
    call flux(hf,Fxf,gxf, Fyf, gyf)
    call flux(hf,Fxp,gxrh, Fyp, gyrh)
    call flux(psiL,Fxp,gxf2, Fyp, gyf2)

    if(ir.EQ.1) then
      flxp3 = - gxp - gyp - rhos2*(gxrh + gyrh)
      flxf3 = - gxf - gyf + gxf2 + gyf2
      gp = gp + ab1*flxp1 + ab2*flxp2 + ab3*flxp3
      gf = gf + ab1*flxf1 + ab2*flxf2 + ab3*flxf3

    elseif(ir.EQ.2) then
      flxp1 = -gxp - gyp - rhos2*(gxrh + gyrh)
      flxf1 = -gxf - gyf + gxf2 + gyf2
      gp = gp + ab1*flxp2 + ab2*flxp3 + ab3*flxp1
      gf = gf + ab1*flxf2 + ab2*flxf3 + ab3*flxf1

    else
      flxp2 = - gxp - gyp - rhos2 * (gxrh + gyrh)
      flxf2 = - gxf - gyf + gxf2 + gyf2
      gp = gp + ab1*flxp3 + ab2*flxp1 + ab3*flxp2
      gf = gf + ab1*flxf3 + ab2*flxf1 + ab3*flxf2

      ir = 0
    endif
    ir = ir + 1

    call filter(gp)
    call filter(gf)
    call rcstr(gp,hp)
    call rcstr(gf,hf)
    call helm1(hp)
    call helm2(hf)
    call laplace(psi,psiL)

    !------------------------------------------------------------
    if(ioutt.ge.noutt) then
      call conserv(gp,gf, F_sum, U_sum)
      call energie(hf,Emxy, Eke, Eki, Epe, Er, Etot)
      ! call helicity(Hc) / not necessary in HW-type equation
      call outt(Emxy,Eke, Eki, Epe, Er, Etot,Hc, F_sum, U_sum)
      ioutt = 0
    endif

    if(ioutf.ge.noutf) then
      call outfield(phi,12)
      call outfield(hf,13)
      ioutf = 0
    endif

    ioutt = ioutt + 1
    ioutf = ioutf + 1

    !TIME LOOP ENDS
  end do


  write(*,*)  't, steps=', time
  write(*,*)  'ioutt, ioutf=', ioutt, ioutf
  write(*,*)  ''
  write(*,*)  'DING! ALL DONE'

  !DESTORY THE FFT PLANS
  call destroy_fft
  call cpu_time(finish)
  write (*,*) 'Calculation Time = ',(finish-start)

  stop
  end
