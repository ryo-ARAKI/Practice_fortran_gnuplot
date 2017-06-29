!http://fluid.web.nitech.ac.jp/Gotoh_Home_page/Edu/Under_Graduate_course/Computational_Phys_I/Diffusion/diffusion_e.f90
!Diffusion equation in 1D by simple Euler

program main

  implicit none
  real(8) x0, x, diff
  integer i, n
  real(8) :: t = 0.0d0
  real(8), parameter :: dt = 0.001  !Time increments
  integer, parameter :: nx = 50     !Number of particion
  !=0: Dirichlet  =1:Neumann
  integer, parameter :: ibs = 0     !Boundary condition
  integer, parameter :: tmax = 400  !Maximum time loop
  integer, parameter :: gdraw = 40  !Output interval
  real(8), dimension(0:nx) :: u     !Value
  real(8), dimension(0:nx) :: uold  !Value in the previous time

 !Make parameters
 real(8), parameter :: kappa = 0.1
  real(8), parameter :: dx = 1.0d0/dble(nx)
  real(8), parameter :: dx2i = 1.0d0/(dx**2)
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
    open(10, file="Diffusion_1D_Euler_Dirichlet.dat")
    u(0) = 0.0d0
    u(nx) = 0.0d0
  else
    write(*,*) "Neumann condition"
    open(10, file="Diffusion_1D_Euler_Neumann.dat")
    u(0) = u(1)
    u(nx) = u(nx-1)
  endif

 !Time advancing
  do n=1, tmax
    do i = 0, nx
      uold(i) = u(i)
    end do

    do i = 1, nx-1
      diff = kappa *(uold(i+1) -2.d0 *uold(i) +uold(i-1)) *dx2i
      u(i) = uold(i) + diff*dt
    end do
   !Set boundary condition again
    if(ibs .eq. 0) then
      u(0) = 0.0d0
      u(nx) = 0.0d0
    else
      u(0) = u(1)
      u(nx) = u(nx-1)
    endif

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
