!Nonlinear first-order advection equation
!http://www.rain.hyarc.nagoya-u.ac.jp/~satoki/main/calc/Fortran/PDE/index.html
!http://www.rain.hyarc.nagoya-u.ac.jp/~satoki/main/calc/Fortran/PDE/NL_advection_periodic_2d.f
program main

  implicit none
 !Set parameters
  integer, parameter :: nmax = 1000   !Maximum number of time iteration
  integer, parameter :: mmax = 100    !Maximum number of space particion
  !
  real(8) u(nmax+1, mmax+1, mmax+1)   !Value of function on whole time
  real(8) v(mmax+1, mmax+1)           !Value of function on the previous time
  real(8) x(mmax+1)                   !x-axis coordinate
  real(8) y(mmax+1)                   !y-axis coordinate
  !
  real(8) :: xmin = 0.0d0             !Minimum x
  real(8) :: xmax = 1.0d0             !Maximum x
  real(8) :: ymin = 0.0d0             !Minimum y
  real(8) :: ymax = 1.0d0             !Maximum x
  !
  real(8) dx, dy, dt                  !Increments
  integer i, j, k                     !Dummy integer for iteration

 !Determine space increments
  dx = (xmax-xmin)/dble(mmax)
  dy = (ymax-ymin)/dble(mmax)

 !Determine time increment
  dt = 1.0d0/dble(nmax)

 !Set space coordinates
  do i = 1, mmax+1
    x(i) = xmin +dx *(i-1)
    y(i) = ymin +dy *(i-1)
  end do

 !Set initial condition
  do i = 1, mmax+1
    do j = 1, mmax+1
      u(1,i,j) = exp(-2.0d0 *((x(i)-0.5d0)**2 +(y(j)-0.5d0)**2))
      v(i,j)   = u(1,i,j)
    end do
  end do

  !Open dat file
   open(10, file="Advection_Nonlinear_2D_Explicit.dat")

 !Time iteration
  do i = 1, nmax
    do j = 2, mmax
      do k = 2, mmax
        !Upwind difference scheme
        u(i+1,j,k) = u(i,j,k) +u(i,j,k) *(dt/dx) *(u(i,j-1,k) -2.0d0*u(i,j,k) +u(i,j,k-1))
      end do
    end do

   !Set boundary condition
    do j = 2, mmax
      u(i+1, j, 1     ) = u(i+1, j, 2)
      u(i+1, j, mmax+1) = u(i+1, j, mmax   )
    end do
    do j = 2, mmax
      u(i+1, 1     , j) = u(i+1, 2 , j)
      u(i+1, mmax+1, j) = u(i+1, mmax , j)
    end do

   !Set boundary conditions in the four corners
    !Averaging two (calculated) data next to the one in the corner
    u(i+1, 1     , 1     ) = 0.5d0 * (u(i+1, 2  , 1     ) +u(i+1, 1     , 2  ))
    u(i+1, 1     , mmax+1) = 0.5d0 * (u(i+1, 2  , mmax+1) +u(i+1, 1     , mmax     ))
    u(i+1, mmax+1, 1     ) = 0.5d0 * (u(i+1, mmax     , 1     ) +u(i+1, mmax+1, 2  ))
    u(i+1, mmax+1, mmax+1) = 0.5d0 * (u(i+1, mmax, mmax+1     ) +u(i+1, mmax+1     , mmax))

   !Copy the result to v
    do k = 1, mmax+1
      do j = 1, mmax+1
        v(j,k) = u(i+1, j, k)
      end do
    end do

   !Output for each 20 time iteration
    !Output evary three data in the space
    if(mod(i,20) .eq. 0) then
      do k = 1, mmax+1, 3
        do j = 1, mmax+1, 3
          write(10,'(3f12.5)') x(j), y(k), v(j,k)
        end do
        write(10,*) ''
      end do
      write(10,*) ''
    endif

  end do

  close(10)



end program main
