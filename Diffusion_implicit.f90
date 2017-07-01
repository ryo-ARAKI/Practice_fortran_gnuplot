!陰的解法による拡散方程式の時間発展
!連立一次方程式がうまく解けていない(右辺が0ではない)
!→教科書P177
module subprog_Laplace_time

  implicit none
contains
 !Subroutine to allocate all matricies
  !Subroutine to determine grid particions
  subroutine allocate_values(phi, phi2, x, n1, n2, d1, d2, f, e, gamma, delta_t)
    real(8), allocatable :: phi(:,:)
    real(8), allocatable :: phi2(:,:)
    real(8), allocatable :: x(:,:,:)
    integer n1, n2
    real(8) d1, d2
    real(8) f, e, gamma
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
    ! if(d1 + d2 > 0.5d0) stop 'We need smaller delta_t'
    ! write(*,*) d1, d2

    gamma = 1.0d0/(1.0d0 +2.0d0*d1 +2.0d0*d2)
    f = -gamma*d2
    e = -gamma*d1

    !Insert coordinates information in matrix x
    do j = 1, n2
      do i = 1, n1
        x(1,i,j) = (dble(i)-1.0d0) * delta_x
        x(2,i,j) = (dble(j)-1.0d0) * delta_y
      end do
    end do
  end subroutine allocate_values

 !Subroutine to set Dirichlet boundary conditions
  subroutine set_dbc(phi, phi2, x, n1, n2)
    integer, intent(in) :: n1, n2
    real(8), intent(out):: phi(n1,n2)
    real(8), intent(out):: phi2(n1,n2)
    real(8) x(2,n1,n2)
    real(8) delta_x
    integer i, j, k
    real(8) :: pi = 2.0d0 *acos(0.0d0)

    delta_x = 1.0d0/(dble(n1)-1.0d0)

    !Boundary condition
    !For phin and phi2(next time step)
    do i = 1, n1
      phi(i,1) = sin(pi*(dble(i)-1.0d0)*delta_x)
      phi2(i,1) = sin(pi*(dble(i)-1.0d0)*delta_x)
      ! write(*,*) phi2(i,1)
    end do
    phi(:,n2) = 0.0d0
    phi(1,:)  = 0.0d0
    phi(n1,:) = 0.0d0
    phi2(:,n2) = 0.0d0
    phi2(1,:)  = 0.0d0
    phi2(n1,:) = 0.0d0

  end subroutine set_dbc


 !Subroutine to calculate (and return) error
  !In the specific time loop
  !NEED TO BE REVIEWED
  function check_error(phi, phi2, f, e, gamma, n1, n2) result(er)
    integer, intent(in) :: n1, n2
    real(8), intent(in) :: phi(n1,n2)
    real(8), intent(in) :: phi2(n1,n2)
    real(8), intent(in) :: f, e, gamma
    real(8) er, calc
    integer i, j

    er = 0.0d0
    do j = 2, n2 -1
      do i = 2, n1 -1
        calc = f * (phi2(i,  j-1) + phi2(i,  j+1)) &
              +e * (phi2(i-1,j  ) + phi2(i+1,j  )) &
              -gamma * phi(i,j)
        er = er + calc**2
      end do
    end do
  end function check_error


 !Subroutine to calculate (and return) error
  !For time steadyness
  function check_steady(phi, phi2, n1, n2, delta_t) result(er)
    integer, intent(in) :: n1, n2
    real(8), intent(in) :: phi(n1,n2)
    real(8), intent(in) :: phi2(n1,n2)
    real(8), intent(in) :: delta_t
    real(8) er, er_calc
    integer i, j

    er_calc = 0.0d0
    do j = 2, n2 -1
      do i = 2, n1 -1
        er_calc = er_calc + (phi2(i,j) - phi(i,j)) ** 2
      end do
    end do
    er = er_calc/delta_t
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
       open(fo, file = 'Diffusion_implicit.d', status='replace', action = 'write')
     else
       open(fo, file = 'Diffusion_implicit.d', position='append', action = 'write')
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
end module subprog_Laplace_time




program main
  use subprog_Laplace_time
  implicit none
  real(8), allocatable :: phi(:,:)  !Fuction
  real(8), allocatable :: phi2(:,:) !Fuction
  real(8), allocatable :: phi3(:,:) !Fuction
  real(8), allocatable :: x(:,:,:)  !Information of grid point
  integer n1, n2                    !Number of particions
  integer i, j, k                   !Dummy integer
  real(8) d1, d2                    !Coefficients
  real(8) f, e, gamma               !Coefficients
  real(8) rhs                       !For calculation
  real(8) :: er0 = 1.0d-3           !Error criteria
  real(8) :: er0_time = 1.0d-5      !Error criteria for time development
  real(8) er                        !error
  real(8) er_time                   !error for time development
  integer istep                     !Number of loop for time development
  integer :: nstep = 100            !Maximum number of loop for time development
  integer itr                       !Number of loop for simultaneous linear equation
  integer :: itrmax = 50            !Maximum number of loop for simultaneous linear equation
  real(8) delta_t                   !time span
  real(8) :: omg = 1.5d0            !SOR coefficients


  call allocate_values(phi, phi2, x, n1, n2, d1, d2, f, e, gamma, delta_t)
    ! write(*,*) "f=", f
    ! write(*,*) "e=", e
    ! write(*,*) "gamma=", gamma

  !Set Dirichlet condition
  call set_dbc(phi, phi2, x, n1, n2)


  do istep = 1, nstep   !Time development
    do itr = 1, itrmax  !Solve simultaneous linear equations
      do j = 2, n2-1
        do i = 2, n1-1
          !Cannot solve simultaneous linear equations
            !I want to solve
            !f\phi^{n+1}_{i,j-1} +e\phi^{n+1}_{i-1,j} +\phi^{n+1}_{i,j} +e\phi^{n+1}_{i+1,j} +f\phi^{n+1}_{i,j+1} = \gamma \phi^{n}_{i,j}
            !for \phi^{n+1}_{i,j} using SOR method
          rhs = -e * (phi2(i-1,  j) + phi2(i+1,  j)) &
                -f * (phi2(i  ,j-1) + phi2(i  ,j+1))
          ! phi2(i,j) = gamma *phi(i,j) + omg * (rhs - phi(i,j))
          phi2(i,j) = gamma *phi(i,j) + rhs
        end do
      end do

      er = check_error(phi, phi2, f, e, gamma, n1, n2)
      if(mod(itr,10)==0) then
        write(*,*) 'Errror is: ', er
      endif

      !Check--------------------------
      if( itr == 10 ) then
           write(*,*) 'phi2 is:'
           do k = 1, n1
           write(*,'(100f7.3)') phi2(k,:)
           enddo
           write(*,*) ''
      endif
      !Check--------------------------

      if( er<er0 ) exit
    enddo

    write(*,*) "Number of loop is: ",itr, "Error is: ",er
    er_time = check_steady(phi, phi2, n1, n2, delta_t)
    phi(:,:) = phi2(:,:)
    if( mod(istep,3) == 0 ) call output(phi, x, n1, n2,istep)

    if( er_time<er0_time ) exit
  end do

  write(*,*) "Number of loop for time development is: ", istep


end program main

! !Check--------------------------
! do k = 1, n1
!   ! write(*,'(100f7.3)') phi(k,:)
!   write(*,'(100f7.3)') phi2(k,:)
! enddo
! write(*,*) ''
! !Check--------------------------
!gnuplot
! cd 'C:\Users\R_Araki\Documents\F90_bluebook'
! set pm3d
! set palette rgbformulae 33,13,10
! splot 'output_ex6_22.d' with pm3d
