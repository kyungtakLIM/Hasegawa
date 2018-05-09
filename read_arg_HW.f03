subroutine read_arg(nstep, noutt, noutf, xl, yl, dt, omega, x0a, y0a, x0b, y0b, sigma_a, sigma_b, m,n,n1,pso)

  use GLOBAL
  implicit none

  !This Routine Reads in simulation parameters
  !either from the command line, or the 'in.com' file
  !if command line arguements are given
  INTEGER :: narg !#of arg & counter of arg
  INTEGER :: nstep, noutt, noutf
  DOUBLE PRECISION :: xl, yl, dt, omega, m,n,n1,pso
  DOUBLE PRECISION :: x0a, y0a, x0b, y0b
  DOUBLE PRECISION :: sigma_a, sigma_b
  CHARACTER(len=20) :: name !Arg name
  CHARACTER(len=20) :: COND

  !Check if any arguments are found
  narg=command_argument_count()
  !Loop over the arguments
  IF (narg.EQ.17) THEN
    call get_command_argument(1,name)
    READ(name,*)dt
    call get_command_argument(2,name)
    READ(name,*)nstep
    call get_command_argument(3,name)
    READ(name,*)noutt
    call get_command_argument(4,name)
    READ(name,*)noutf
    call get_command_argument(5,name)
    READ(name,*)xl
    call get_command_argument(6,name)
    READ(name,*)yl
    call get_command_argument(7,name)
    READ(name,*)x0a
    call get_command_argument(8,name)
    READ(name,*)y0a
    call get_command_argument(9,name)
    READ(name,*)x0b
    call get_command_argument(10,name)
    READ(name,*)y0b
    call get_command_argument(11,name)
    READ(name,*)sigma_a
    call get_command_argument(12,name)
    READ(name,*)sigma_b
    call get_command_argument(13,name)
    READ(name,*)de
    call get_command_argument(14,name)
    READ(name,*)rhos
    call get_command_argument(15,name)
    READ(name,*)rhoi
    call get_command_argument(16,name)
    READ(name,*)omega
    call get_command_argument(17,name)
    READ(name,*)COND
    call get_command_argument(18,name)
    READ(name,*)m
    call get_command_argument(19,name)
    READ(name,*)n
    call get_command_argument(20,name)
    READ(name,*)n1
    call get_command_argument(21,name)
    READ(name,*)pso
  ELSE
    !Read data from file='in.com'
    print *,'Reading Data from in.com ....'
    open(unit=7,status='old',file='in.com')
    read(7,*)
    read(7,*)
    read(7,*) dt, nstep, noutt, noutf
    read(7,*)
    read(7,*) xl, yl
    read(7,*)
    read(7,*) x0a, y0a, x0b, y0b
    ! centers of the theta-pinches
    read(7,*)
    read(7,*) sigma_a, sigma_b
    ! square half-width of gaussians
    read(7,*)
    read(7,*) de, rhos, rhoi
    read(7,*)
    read(7,*) omega
    read(7,*)
    read(7,*) pso
  END IF

  return
  end
