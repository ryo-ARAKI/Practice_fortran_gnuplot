!http://fluid.web.nitech.ac.jp/Gotoh_Home_page/Edu/Under_Graduate_course/Computational_Phys_I/Diffusion/diffusion_cn.f90
!Diffusion equation in 1D by simple Euler

program main

  implicit none
  real(8) x0, x, diff, sigmai
  integer i, n
  real(8) :: t = 0.0d0
  real(8), parameter :: dt = 0.004  !Time increments
  integer, parameter :: nx = 50     !Number of particion
  !=0: Dirichlet  =1:Neumann
  integer, parameter :: ibs = 0     !Boundary condition
  integer, parameter :: tmax = 100  !Maximum time loop
  integer, parameter :: gdraw = 10  !Output interval
  real(8), dimension(0:nx)   :: u     !Value
  real(8), dimension(0:nx)   :: uold  !Value in the previous time
  real(8), dimension(0:nx)   :: alpha !Coefficient for Crank-Nicolson
  real(8), dimension(0:nx)   :: beta  !Coefficient for Crank-Nicolson
  real(8), dimension(0:tmax) :: q     !Coefficient for Crank-Nicolson

 !Make parameters
  real(8), parameter :: kappa = 0.1
  real(8), parameter :: dx = 1.0d0/dble(nx)
  real(8), parameter :: dx2i = 1.0d0/(dx**2)
  real(8), parameter :: bi = 2.0d0 * dx**2 /(kappa *dt)
  real(8), parameter :: pi = 2.0d0 * acos(0.0d0)
  real(8), parameter :: amp = 1.0d-1
  real(8), parameter :: t_init = 0.001

 !Initial condition
  x0 = dx * (nx/2.0d0)
  do i = 1, nx
    x = i*dx
    u(i) = amp/sqrt(2.0d0*pi*t_init) *exp(-(x-x0)**2/(4.0d0*kappa*t_init))
    if( u(i) .lt. 1.0d-20) u(i) = 0.0d0
  end do

 !Boundary condition
  if(ibs .eq. 0) then
    write(*,*) "Dirichlet condition"
    open(10, file="Diffusion_1D_Crank-Nicolson_Dirichlet.dat")
    u(0) = 0.0d0
    u(nx) = 0.0d0
  else
    write(*,*) "Neumann condition"
    open(10, file="Diffusion_1D_Crank-Nicolson_Neumann.dat")
    u(0) = u(1)
    u(nx) = u(nx-1)
  endif

 !Time advancing
  do n=1, tmax
   !Enforce boundary conditions at x=0
    if(ibs .eq. 0) then
      alpha(0) = 0.0d0
      beta(0) = 0.0d0
    else
      alpha(0) = 2.0d0/(1.0d0 +bi)
      beta(0) = (u(1) -(2.0d0 -bi) *u(0) +u(1))/(2.0d0 +bi)
    endif

   !Forward elimination
    do i=1, nx-1
      sigmai = 1.0d0/(2.d0 +bi -alpha(i-1))
      alpha(i) = sigmai
      diff = u(i+1) -(2.d0 -bi) *u(i) +u(i-1)
      beta(i) = sigmai *(diff +beta(i-1))
    end do

   !Enforce Boundary conditions at x=1
    if(ibs .eq. 0) then
      u(nx) = 0.0d0
    else
      alpha(nx) = 2.d0/(2.d0 +bi)
      beta(nx)=(u(nx-1) -(2.d0 -bi) *u(nx) +u(nx-1))/(2.d0 +bi)
      u(nx) = alpha(nx) *u(nx-1) +beta(nx)
    endif

   !Backward elimination
    do i = nx, 1, -1
      u(i-1) = alpha(i-1) * u(i) + beta(i-1)
    end do

   !Output data
    if(mod(n,gdraw) .eq. 0) then
      do i = 0, nx
        write(10, '(2f12.5)') i*dx, u(i)
      end do
      write(10, *) ""
      write(10, *) ""
    endif

  end do

end program main
