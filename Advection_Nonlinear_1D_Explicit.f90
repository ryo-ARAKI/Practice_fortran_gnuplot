!Nonlinear first-order advection equation
!http://www.rain.hyarc.nagoya-u.ac.jp/~satoki/main/calc/Fortran/PDE/index.html
!http://www.rain.hyarc.nagoya-u.ac.jp/~satoki/main/calc/Fortran/PDE/non_linear_advection_explicit.f
program main

  implicit none
 !Set parameters
  integer, parameter :: nmax = 5000 !Total time step
  integer, parameter :: mmax = 80   !Total particions in space
  real(8) u(0:nmax, 0:mmax)         !Function
  real(8) x(0:mmax)                 !x-axis position
  real(8) y(0:mmax)                 !Copy u(i+1,:) a value of function in certaion time step
!
  real(8) :: xmin = 0.0d0   !Minimum x range
  real(8) :: xmax = 1.0d0   !Maximum x range
!
  real(8) dx, dt
  integer i, j
 !
  dx = (xmax-xmin)/mmax   !x increments
  dt = 1.0d0/nmax         !Time increments

 !Open dat file
  open(10, file="Advection_Nonlinear_1D_Explicit.dat")

 !Set space coordinates
  x(0) = xmin
  x(mmax) = xmax
  do j = 1, mmax-1
    x(j) = xmin + dx*dble(j)
  end do

 !Set initial condition
  do i = 0, mmax
    u(0,i) = exp(-20.0d0 *(x(i)-0.5d0)**2)
  end do

 !Core calculation
  do i = 0, nmax-1    !Time development
    do j = 1, mmax-1  !Space development
      u(i+1,j) = u(i,j) -0.5d0 *u(i,j) *(dt/dx) *(u(i,j)-u(i,j-1))  !Up-wind difference
      y(j) = u(i+1,j)
    end do
   !Set boundary condition
    u(i+1,0) = u(i+1, mmax-1) !Not sure
    u(i+1,mmax) = u(i+1, 1)   !Not sure
    y(0) = u(i+1, 0)          !Because do-loop does not includes j=0
    y(mmax) = u(i+1, mmax)    !Because do-loop does not includes j=mmax

   !Output results for each 20 time step
    if(mod(i,20) .eq. 0) then
      do j = 0, mmax
        write(10,'(2f12.5)') x(j), u(i+1,j)
      end do
      write(10,*) ''
      write(10,*) ''
    endif



  end do

  close(10)

end program main
