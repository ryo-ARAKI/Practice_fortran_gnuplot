!http://www.geo.titech.ac.jp/lab/ida/numexe/manual/8/diffusion.html
program main

  implicit none
  integer, parameter :: jmax = 10   !Number of particion of the rod
  integer, parameter :: nmax = 250  !Maximum number of time development
  integer, parameter :: nint = 25   !Output time interval

  real(8), dimension(0:jmax) :: u
  real(8), dimension(jmax-1) :: a, b, c, d, x, w

  real(8) dx, dt, r

  integer j, n

  character(len=33) :: fname = 'output_Diffusion_Crank-Nicolson.d'

 !Define time and space increments
  dx = 1.0d0/dble(jmax)
  dt = 0.004d0
  r  = dt/dx**2

 !Set initial condition
  u(0) = 0.0d0
  do j = 1, jmax-1
    u(j) = 1.0d0
  enddo
  u(jmax) = 0.0d0

!Time development
  do n = 0, nmax

   !Output results for certain time increments
    if(mod(n,nint)==0) then
      open(10, file=fname, position='append')
      do j = 0, jmax
        write(10, '(2f12.5)') dble(j)*dx, u(j)
      enddo
      write(10,*) ''
      write(10,*) ''
      close(10)
    endif

   !Set coefficients
    do j = 1, jmax-1
      a(j) = -r
      b(j) = 2.0d0*(r+1.0d0)
      c(j) = -r
    enddo

   !Set right-hand side of simultaneous linear equation
    !d(1) and d(jmax-1) are needed for different boundary conditions
    !Difference in coefficients of d(1) and d(jmax-1) stands for boundary condition
    d(1) = u(2)*r -u(1)*(r-1.0d0)*2.0d0 +u(0)*r*2.0d0
    do j = 2, jmax-2
      d(j) = u(j+1)*r -u(j)*(r-1.0d0)*2.0d0 +u(j-1)*r
    end do
    d(jmax-1) = u(jmax)*r*2.0d0 -u(jmax-1)*(r-1.0d0)*2.0d0 +u(jmax-2)*r

   !Solve the simultaneous linear equation in subroutine
    call tridag(a,b,c,d,x,w,jmax-1)
    !Copy the solution
    do j = 1, jmax-1
      u(j) = x(j)
    enddo

  enddo

end program main

!Solve simultaneous linear equation in which coefficient matrix is tridiagonal
subroutine tridag(a,b,c,d,x,w,n)
  implicit none

  integer n
  real(8), dimension(n) :: a, b, c, d, x, w

  real(8) bet
  integer j

  bet = b(1)
  x(1) = d(1)/bet
  do j = 2, n
    w(j) = c(j-1)/bet
    bet = b(j) -a(j)*w(j)
    x(j) = (d(j) -a(j)*x(j-1))/bet
  enddo

  do j = n-1, 1, -1
    x(j) = x(j) -w(j+1)*x(j+1)
  enddo

  return

end subroutine tridag
