module subprog_Laplace_time

  implicit none
contains
 !Subroutine to allocate all matricies
  !To determine grid particions
  subroutine allocate_values(phi, phi2, x, n1, n2, d1, d2)
    real(8), allocatable :: phi(:,:)
    real(8), allocatable :: phi2(:,:)
    real(8), allocatable :: x(:,:,:)
    integer n1, n2
    real(8) d1, d2
    real(8) delta_x, delta_y, delta_t
    real(8) :: alpha = 5.0d-1
    integer i, j

    write(*,*) "Enter number of point in x-axis, n1"
    read(*,*) n1
    write(*,*) "Enter number of point in y-axis, n2"
    read(*,*) n2
    allocate(phi(n1,n2))
    allocate(phi2(n1,n2))
    allocate(x(2,n1,n2))

    write(*,*) "Enter time-span delta_t"
    read(*,*) delta_t

    !Set coordinate increment
    delta_x = 1.0d0/(dble(n1)-1.0d0)
    delta_y = 1.0d0/(dble(n2)-1.0d0)

    !Set coefficient (See textbook)
    d1 = alpha * delta_t/(delta_x**2)
    d2 = alpha * delta_t/(delta_y**2)
    if(d1 + d2 > 0.5d0) stop 'We need smaller delta_t'
    ! write(*,*) d1, d2

    !Insert coordinates information in matrix x
    do j = 1, n2
      do i = 1, n1
        x(1,i,j) = (dble(i)-1.0d0) * delta_x
        x(2,i,j) = (dble(j)-1.0d0) * delta_y
      end do
    end do
  end subroutine allocate_values

 !Subroutine to set Dirichlet boundary conditions
  subroutine set_dbc(phi, phi2, x, n1, n2, d1, d2)
    integer, intent(in) :: n1, n2
    real(8) phi(n1,n2)
    real(8) phi2(n1,n2)
    real(8) x(2,n1,n2)
    real(8) d1, d2
    real(8) delta_x
    integer i, j
    real(8) :: pi = 2.0d0 *acos(0.0d0)

    delta_x = 1.0d0/(dble(n1)-1.0d0)

    !Boundary condition
    !For phin and phi2(next time step)
    do i = 1, n1
      phi(i,1) = sin(pi*(dble(i)-1.0d0)*delta_x)
      phi2(i,1) = sin(pi*(dble(i)-1.0d0)*delta_x)
    end do
    phi(:,n2) = 0.0d0
    phi(1,:)  = 0.0d0
    phi(n1,:) = 0.0d0
    phi2(:,n2) = 0.0d0
    phi2(1,:)  = 0.0d0
    phi2(n1,:) = 0.0d0


  end subroutine set_dbc


 !Subroutine to calculate (and return) error
  function check_steady(phi, phi2, n1, n2) result(er)
    integer, intent(in) :: n1, n2
    real(8), intent(in) :: phi(n1,n2)
    real(8), intent(in) :: phi2(n1,n2)
    real(8) er
    integer i, j

    er = 0.0d0
    do j = 2, n2 -1
      do i = 2, n1 -1
        er = er + (phi2(i,j) - phi(i,j))**2
      end do
    end do
  end function check_steady

 !Subroutine to output results
  subroutine output(phi, x, n1, n2,istep)
    integer, intent(in) :: n1, n2, istep
    real(8), intent(in) :: phi(n1,n2)
    real(8), intent(in) :: x(2,n1,n2)
    integer i, j
    integer :: fo = 11
    !Open new file for istep ==1
    !Open the same file for istep>1
    !It is not working well
    if(istep == 1) then
      open(fo, file = 'Diffusion_explicit.d', status='new', action = 'write')
    else
      open(fo, file = 'Diffusion_explicit.d', position='append', action = 'write')
    endif

    do j = 1, n2
      do i = 1, n1
        ! write(*,'(100f7.3)') phi(i,:)
        write(fo,'(3e12.4)') x(:,i,j), phi(i,j)
      enddo
      write(fo,*) ''
    enddo

    !Insert two blank lines
    write(fo,*)''
    ! write(fo,*)''

    close(fo)
  end subroutine output

 !Subroutine to output final result
  ! subroutine output_final(phi, x, n1, n2)
  !   integer, intent(in) :: n1, n2
  !   real(8), intent(in) :: phi(n1,n2)
  !   real(8), intent(in) :: x(2,n1,n2)
  !   integer i, j
  !   integer :: fo = 11
  !   !Open new file for istep ==1
  !   !Open the same file for istep>1
  !     open(fo, file = 'output_ex6_19_final.d', status='replace', action = 'write')
  !
  !   do j = 1, n2
  !     do i = 1, n1
  !       ! write(*,'(100f7.3)') phi(i,:)
  !       write(fo,'(3e12.4)') x(:,i,j), phi(i,j)
  !     enddo
  !     write(fo,*) ''
  !   enddo
  !
  !   close(fo)
  ! end subroutine output_final

end module subprog_Laplace_time




program main
  use subprog_Laplace_time
  implicit none
  real(8), allocatable :: phi(:,:)  !Fuction
  real(8), allocatable :: phi2(:,:) !Fuction
  real(8), allocatable :: x(:,:,:)  !Information of grid point
  integer n1, n2                    !Number of particions
  integer i, j                      !Dummy integer
  real(8) d1, d2                    !Coefficients
  real(8) :: er0 = 1.0d-5           !Error criteria
  real(8) er                        !error
  integer istep                     !Number of loop
  integer :: nstep = 100000         !Maximum number of loop

  call allocate_values(phi, phi2, x, n1, n2, d1, d2)
  !Set initial condition
  phi(2:n1-1, 2:n2-1) = 0.0d0
  call set_dbc(phi, phi2, x, n1, n2, d1, d2)

  do istep = 1, nstep !Time development
    do j = 2, n2-1
      do i = 1, n1-1
        phi2(i,j) = phi(i,j) &
                    + d1 * (phi(i-1,j  ) -2.0d0 * phi(i,j) + phi(i+1,j  )) &
                    + d2 * (phi(i  ,j-1) -2.0d0 * phi(i,j) + phi(i  ,j+1))
      end do
    end do
  er = check_steady(phi, phi2, n1, n2)  !Check time steadyness
  phi(:,:) = phi2(:,:)
  !Output all information for each 10 time step
  if(mod(istep,10)==0) then
    call output(phi, x, n1, n2,istep)
  endif
  if( er<er0 ) exit
  end do

  write(*,*) "Loop number is: ", istep
  ! call output_final(phi, x, n1, n2)


end program main
