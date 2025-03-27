module commondata
        use mpi
        implicit none

        ! Grid parameters
        integer(kind=4), parameter :: nx=513, ny=nx
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1

        real(kind=8), parameter :: Ra=1e5          !Rayleigh number
        real(kind=8), parameter :: Ma=0.1d0         !Mach number
        real(kind=8), parameter :: Pr=0.71d0        !Prandtl number
        real(kind=8), parameter :: Thot=1.0d0, Tcold=0.0d0, Tref=0.0d0
        real(kind=8) :: kappa                        !Thermal expansion coefficient
        real(kind=8) :: gbeta                        !Volume expansion coefficient * Gravitational acceleration

        ! Parameters in 0 system
        real(kind=8), parameter :: rho0 = 1.0d0
        real(kind=8) :: length0                   ! Characteristic length
        real(kind=8) :: viscosity0                ! Kinematic viscosity
        real(kind=8) :: dt0                       ! Time step
        real(kind=8) :: dx0                       ! Grid spacing

        ! Parameters in LB system
        real(kind=8) :: length_LB                 ! Characteristic length
        real(kind=8) :: viscosity_LB              ! Kinematic viscosity
        real(kind=8) :: tauf, taut                ! Relaxation time
        real(kind=8) :: dt                       ! Time step

        ! Iteration control
        integer(kind=4) :: itc = 0            ! Current iteration count
        integer(kind=4), parameter :: itc_max = 20000  ! Maximum iterations

        ! Convergence criteria
        real(kind=8) :: errorU, errorT              ! Current error
        real(kind=8), parameter :: epsU=1e-6, epsT=1e-6  ! Convergence threshold

        ! Grid coordinates
        real(8) :: x(0:nx+1),y(0:ny+1)
        real(8), allocatable:: xGrid(:), yGrid(:)

        ! Flow variables
        real(kind=8), allocatable :: u(:,:), v(:,:)      ! Velocity components
        real(kind=8), allocatable :: temp(:,:)           ! temperature
        real(kind=8), allocatable :: rho(:,:)            ! Density
        real(kind=8), allocatable :: up(:,:), vp(:,:)    ! Previous velocity components for error checking
        real(kind=8), allocatable :: utemp(:,:)          ! Previous temperature for error checking
        real(kind=8), allocatable :: Fx(:,:), Fy(:,:)    ! force components
        real(8), allocatable :: u_all(:,:), v_all(:,:), rho_all(:,:), T_all(:,:)

        ! Distribution functions
        real(kind=8), allocatable :: f(:,:,:), f_post(:,:,:)  ! Current and post-collision distributions
        real(kind=8),allocatable :: g(:,:,:), g_post(:,:,:)

        ! MRT relaxation parameters
        real(kind=8) :: omega_U(0:8), omega_T(0:4)  ! Relaxation rates for MRT

        ! Lattice directions
        integer(kind=4) :: ex(0:8), ey(0:8)  ! Lattice velocity directions
        data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
        data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/

        ! Additional MRT parameters
        real(kind=8) :: Snu, Sq, sig_k

        integer(kind=4), allocatable :: inter_x(:,:), inter_y(:,:)
!----------------------------MPI---------------------------------------------------
        integer(kind=4) :: NPROC, COMM2D, IERR, ISTATUS(MPI_STATUS_SIZE), MYID, MY_COORD(2)
        integer(kind=4), allocatable :: countX(:), countY(:)
        integer(kind=4), allocatable :: i_start(:), j_start(:)
        integer(kind=4) :: top, down, left, right
        integer(kind=4) :: dims(2)
        logical :: periods(2)
        data periods/2*.FALSE./
        integer(kind=4) :: nxLocal, nyLocal
        integer(kind=4), allocatable :: displX(:), displY(:)
end module commondata

subroutine mesh()
    use commondata
    implicit none

    integer(kind=4) :: i, j
    real(kind=8) :: dx(nx+1), dy(ny+1)
    real(kind=8) :: constA
    allocate(xGrid(0:nxLocal+1))
    allocate(yGrid(0:nyLocal+1))


    if(MYID == 0)then
        ! Compute grid coordinates
        constA = 3.2d0
        do i = 0, nx+1
            x(i) = 0.5d0 * (erf(constA  * (dble(i) / dble(nx+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
        end do
        do j = 0, ny+1
            y(j) = 0.5d0 * (erf(constA  * (dble(j) / dble(ny+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
        end do

        ! Compute grid spacing using array slicing
        dx(1:nx+1) = x(1:nx+1) - x(0:nx)
        dy(1:ny+1) = y(1:ny+1) - y(0:ny)
        dx0 = dx(1)
        dt0 = dx0
        write(*,*) "nx =", nx, ", ny =", ny
        write(*,*) "Ra =", real(Ra)
        write(*,*) "CPU in each dimension is: ", dims(1), dims(2)
        write(*,*) "  "
        write(*,*) "---in 0 system---"
        write(*,*) "deltaX = ", dx0
        write(*,*) "deltaT = ", dt0

        length0 = x(nx+1)-x(1)
        write(*,*) "length0 = ", real(length0)
        write(*,*) "  "

        ! Calculate viscosity based on LB unit
        length_LB = 1.0d0 / dx(1)
        dt = dt0 * length_LB

        x(1:nx) = x(1:nx)-dx0/2.0d0
        y(1:ny) = y(1:ny)-dx0/2.0d0
        x=x*length_LB
        y=y*length_LB
    endif

    call MPI_SCATTERV(x(1),countX,displX,MPI_REAL8,xGrid(1),countX(MYID),MPI_REAL8,0,COMM2D,IERR)
    call MPI_SCATTERV(y(1),countY,displY,MPI_REAL8,yGrid(1),countY(MYID),MPI_REAL8,0,COMM2D,IERR)

    call MPI_SENDRECV(xGrid(nxLocal),1,MPI_REAL8,right,0,xGrid(0),1,MPI_REAL8,left,0,COMM2D,ISTATUS,IERR)
    call MPI_SENDRECV(xGrid(1),1,MPI_REAL8,left,1,xGrid(nxLocal+1),1,MPI_REAL8,right,1,COMM2D,ISTATUS,IERR)

    call MPI_SENDRECV(yGrid(nyLocal),1,MPI_REAL8,top,0,yGrid(0),1,MPI_REAL8,down,0,COMM2D,ISTATUS,IERR)
    call MPI_SENDRECV(yGrid(1),1,MPI_REAL8,down,1,yGrid(nyLocal+1),1,MPI_REAL8,top,1,COMM2D,ISTATUS,IERR)

    call MPI_BCAST(dt,1,MPI_REAL8,0,COMM2D,IERR)
    call MPI_BCAST(length_LB,1,MPI_REAL8,0,COMM2D,IERR)

    return
end subroutine

subroutine initial()
    use commondata
    implicit none

    integer(kind=4) :: i, j, alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: velocitySquare
    allocate(inter_x(nxLocal,3))
    allocate(inter_y(nyLocal,3))
    ! Initialize iteration count and error
    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0

    ! Calculate viscosity based on LB unit
    viscosity_LB = (Ma*(length_LB-1.0d0)*dsqrt(Pr))/dsqrt(3.0d0*Ra)
    kappa = viscosity_LB/Pr
    gbeta = Ra*viscosity_LB*kappa/((length_LB-1.0d0)**3)

    ! Calculate relaxation time
    tauf = viscosity_LB * 3.0d0 + 0.5d0

    ! Calculate MRT relaxation parameters
    Snu = 1.0d0/tauf
    Sq = 8.0d0*(2.0d0-Snu)/(8.0d0-Snu)
    sig_k = 1.0d0/(0.5d0+4.0d0*(1.0d0/Snu-0.5d0)/(3.0d0*Pr))

    if(MYID == 0)then
         !Output initial parameters for verification
        write(*,*) "---in system LB"
        write(*,*) "deltaX = ", x(1)-x(0)
        write(*,*) "deltaT = ", dt
        write(*,*) "characteristic length   =", real(length_LB), "l.u."
        write(*,*) "viscosity_LB =", real(viscosity_LB), "l.u.^2/t.s."
        write(*,*) "timeStep ratio for (uniform) / (non-uniform) : ", real(length_LB / dble(nx))
        write(*,*) "tauf =", real(tauf)
        write(*,*) "    "
    end if

    !Allocate flow variables
    allocate (u(nxLocal,nyLocal))
    allocate (v(nxLocal,nyLocal))
    allocate (rho(nxLocal,nyLocal))
    allocate (up(nxLocal,nyLocal))
    allocate (vp(nxLocal,nyLocal))
    allocate (temp(nxLocal,nyLocal))
    allocate (utemp(nxLocal,nyLocal))
    allocate (Fx(nxLocal,nyLocal))
    allocate (Fy(nxLocal,nyLocal))

    if(MYID == 0)then
        allocate (u_all(nx,ny))
        allocate (v_all(nx,ny))
        allocate (rho_all(nx,ny))
        allocate (T_all(nx,ny))
    end if

    allocate (f(0:8,nxLocal,nyLocal))
    allocate (f_post(0:8,0:nxLocal+1,0:nyLocal+1))
    allocate (g(0:4,nxLocal,nyLocal))
    allocate (g_post(0:4,0:nxLocal+1,0:nyLocal+1))

    !Initialize flow variables
    rho = rho0
    temp = 0.0d0
    utemp=0.0d0
    u = 0.0d0
    v = 0.0d0
    up = 0.0d0
    vp = 0.0d0

    do j=1,nyLocal
        do i=1,nxLocal
            temp(i,j) = dble(i_start(MYID)+i-2)/dble(nx-1)*(Tcold-Thot)+Thot
        enddo
    enddo

    omega_U(0) = 4.0d0/9.0d0
    do alpha=1,4
        omega_U(alpha) = 1.0d0/9.0d0
    enddo
    do alpha=5,8
        omega_U(alpha) = 1.0d0/36.0d0
    enddo

    omega_T(0) = 1.0d0/2.0d0
    do alpha=1,4
        omega_T(alpha) = 1.0d0/8.0d0
    enddo

    do j = 1, nyLocal
        do i = 1, nxLocal
            velocitySquare = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            do alpha = 0,8
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                f(alpha, i, j) = rho(i,j)*omega_U(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*velocitySquare)
            enddo

            do alpha=0,4
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                g(alpha,i,j)=omega_T(alpha)*temp(i,j)*(1.0d0+4.0d0*un(alpha))
            end do
        enddo
    enddo

    if(MY_COORD(1) == 0)then
        do i = 1, nxLocal
            if(i == 1)then
                inter_x(i,:) = (/i+1, i, i+2/)

            else
                inter_x(i,:) = (/i-1, i, i+1/)
            end if
        end do

    elseif(MY_COORD(1) == dims(1)-1)then
        do i = 1, nxLocal
            if(i == nxLocal)then
                inter_x(i,:) = (/i-1, i, i-2/)

            else
                inter_x(i,:) = (/i-1, i, i+1/)
            end if
        end do

    else
        do i = 1, nxLocal
            inter_x(i,:) = (/i-1, i, i+1/)
        end do
    end if

    if(MY_COORD(2) == 0)then
        do j = 1, nyLocal
            if(j == 1)then
                inter_y(j,:) = (/j+1, j, j+2/)

            else
                inter_y(j,:) = (/j-1, j, j+1/)
            end if
        end do

    elseif(MY_COORD(2) == dims(2)-1)then
        do j = 1, nyLocal
            if(j == nyLocal)then
                inter_y(j,:) = (/j-1, j, j-2/)

            else
                inter_y(j,:) = (/j-1, j, j+1/)
            end if
        end do

    else
        do j = 1, nyLocal
            inter_y(j,:) = (/j-1, j, j+1/)
        end do
    end if

    return
end subroutine initial

program main
    use commondata
    implicit none
    real(kind=8) :: timestart, timeEnd
    !---------------------------------------------------------------------------------------
    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERR)
    allocate(countX(0:NPROC-1))
    allocate(countY(0:NPROC-1))
    allocate(i_start(0:NPROC-1))
    allocate(j_start(0:NPROC-1))
    allocate(displX(0:NPROC-1))
    allocate(displY(0:NPROC-1))

    dims = 0
    call MPI_DIMS_CREATE(NPROC, 2, dims, IERR)
    call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods, .TRUE.,COMM2D,IERR)
    call MPI_COMM_RANK(COMM2D,MYID,IERR)
    call MPI_CART_COORDS(COMM2D,MYID,2,MY_COORD,IERR)
    call MPI_CART_SHIFT(COMM2D,0,1,left,right,IERR)
    call MPI_CART_SHIFT(COMM2D,1,1,down,top,IERR)

    call StartEnd(1,nx,1,ny,MY_COORD)
    nxLocal = countX(MYID)
    nyLocal = countY(MYID)

   ! --------------------------------------------------------------------------------------------
    call MPI_BARRIER(MPI_COMM_WORLD, IERR)
    timeStart = MPI_WTIME()

    call mesh()

    call initial()

    do while(((errorU > epsU).or.(errorT > epsT)).AND.(itc < itc_max))

        itc = itc+1

        call collision_U()

        call collision_T()

        call update()

        call interpolate()

        call bounceback_u()

        call bounceback_T()

        call macro_u()

        call macro_t()

        if(MOD(itc,2000).EQ.0) then
            call check()
        endif
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD, IERR)
    timeEnd = MPI_WTIME()

    if(MYID == 0)then
        write(*,*)"Time=", timeEnd-timestart, "s"
        write(*,*) "MLUPS = ", real(dble(nx*ny)/1e6*dble(itc)/(timeEnd-timeStart))
    end if

    call gather_data()

    if(MYID == 0)then
        call output_ASCII()
        write(*,*) "Successfully: DNS completed!"
    end if

    deallocate(u)
    deallocate(v)
    deallocate(rho)
    deallocate(up)
    deallocate(vp)
    deallocate(f)
    deallocate(f_post)
    deallocate(g)
    deallocate(g_post)
    deallocate(temp)
    deallocate(utemp)
    deallocate(xGrid)
    deallocate(yGrid)

    if(MYID == 0)then
        deallocate(u_all)
        deallocate(v_all)
        deallocate(rho_all)
        deallocate(T_all)
    end if

    deallocate(countX)
    deallocate(countY)
    deallocate(i_start)
    deallocate(j_start)
    deallocate(displX)
    deallocate(displY)

    call MPI_FINALIZE(IERR)

    stop
end program main

subroutine StartEnd(iSx1, iSx2, iSy1, iSy2, coords)
    use commondata
    implicit none
    integer(kind=4) :: coords(1:2), iBlock(1:2)
    integer(kind=4) :: leng(1:2), ir(1:2)
    integer(kind=4) :: iSx1, iSx2, iSy1, iSy2
    integer(kind=4) :: startX, startY, localCountX, localCountY
    integer(kind=4) :: i, j
    integer(kind=4) :: startPoint(1:2), endPoint(1:2)
    integer(kind=4) :: dir

    startPoint(1) = iSx1
    endPoint(1) = iSx2

    startPoint(2) = iSy1
    endPoint(2) = iSy2

    do dir = 1, 2
        leng(dir) = endPoint(dir)-startPoint(dir)+1
        iBlock(dir) = leng(dir)/dims(dir)
        ir(dir) = leng(dir)-iBlock(dir)*dims(dir)
    enddo

    if(coords(1).LT.ir(1)) then
        localCountX = iBlock(1)+1
        startX = startPoint(1)+coords(1)*(iBlock(1)+1)
    elseif(coords(1).GE.ir(1)) then
        localCountX = iBlock(1)
        startX = startPoint(1)+coords(1)*iBlock(1)+ir(1)
    endif

    if(coords(2).LT.ir(2)) then
        localCountY = iBlock(2)+1
        startY = startPoint(2)+coords(2)*(iBlock(2)+1)
    elseif(coords(2).GE.ir(2)) then
        localCountY = iBlock(2)
        startY = startPoint(2)+coords(2)*iBlock(2)+ir(2)
    endif


    call MPI_ALLGATHER(localCountX, 1, MPI_INTEGER, countX, 1, MPI_INTEGER, COMM2D, IERR)
    call MPI_ALLGATHER(localCountY, 1, MPI_INTEGER, countY, 1, MPI_INTEGER, COMM2D, IERR)

    call MPI_ALLGATHER(startX, 1, MPI_INTEGER, i_start, 1, MPI_INTEGER, COMM2D, IERR)
    call MPI_ALLGATHER(startY, 1, MPI_INTEGER, j_start, 1, MPI_INTEGER, COMM2D, IERR)

    do i=0,NPROC-1
       displX(i)=i_start(i)-1
       displY(i)=j_start(i)-1
    end do

    return
end subroutine StartEnd

subroutine collision_U()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: s(0:8)
    real(kind=8) :: m(0:8)
    real(kind=8) :: m_post(0:8)
    real(kind=8) :: meq(0:8)
    real(kind=8) :: fSource(0:8)

    do j=1,nyLocal
        do i=1,nxLocal

    m(0) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(1) = -4.0d0*f(0,i,j)-f(1,i,j)-f(2,i,j)-f(3,i,j)-f(4,i,j)+2.0d0*(f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j))
    m(2) = 4.0d0*f(0,i,j)-2.0d0*(f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j))+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(3) = f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
    m(4) = -2.0d0*f(1,i,j)+2.0d0*f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
    m(5) = f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
    m(6) = -2.0d0*f(2,i,j)+2.0d0*f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
    m(7) = f(1,i,j)-f(2,i,j)+f(3,i,j)-f(4,i,j)
    m(8) = f(5,i,j)-f(6,i,j)+f(7,i,j)-f(8,i,j)

            meq(0) = rho(i,j)
            meq(1) = rho(i,j)*( -2.0d0+3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(2) = rho(i,j)*( 1.0d0-3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(3) = rho(i,j)*u(i,j)
            meq(4) = -rho(i,j)*u(i,j)
            meq(5) = rho(i,j)*v(i,j)
            meq(6) = -rho(i,j)*v(i,j)
            meq(7) = rho(i,j)*(u(i,j)*u(i,j)-v(i,j)*v(i,j))
            meq(8) = rho(i,j)*(u(i,j)*v(i,j))

            s(0) = 0.0d0
            s(1) = Snu
            s(2) = Snu
            s(3) = 0.0d0
            s(4) = Sq
            s(5) = 0.0d0
            s(6) = Sq
            s(7) = Snu
            s(8) = Snu

            Fx(i,j) = 0.0d0
            Fy(i,j) = gbeta*temp(i,j)*rho(i,j)

            fSource(0) = 0.0d0
            fSource(1) = (6.0d0-3.0d0*s(1))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(2) = -(6.0d0-3.0d0*s(2))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(3) = Fx(i,j)
            fSource(4) = -(1.0d0-0.5d0*s(4))*Fx(i,j)
            fSource(5) = Fy(i,j)
            fSource(6) = -(1.0d0-0.5d0*s(6))*Fy(i,j)
            fSource(7) = (2.0d0-s(7))*(u(i,j)*Fx(i,j)-v(i,j)*Fy(i,j))
            fSource(8) = (1-0.5d0*s(8))*(v(i,j)*Fx(i,j)+u(i,j)*Fy(i,j))


            do alpha=0,8
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)*dt
            enddo

    f_post(0,i,j) = ( m_post(0)-m_post(1)+m_post(2) )/9.0d0
    f_post(1,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                    +m_post(7)*0.25d0
    f_post(2,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)*0.25d0
    f_post(3,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                    +m_post(7)*0.25d0
    f_post(4,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)*0.25d0
    f_post(5,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)*0.25d0
    f_post(6,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)*0.25d0
    f_post(7,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)*0.25d0
    f_post(8,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)*0.25d0
        enddo
    enddo

    do j = 1, nyLocal
        do i = 1, nxLocal
            f(0,i,j) = f_post(0,i,j)
        enddo
    enddo

    return
end subroutine collision_U

subroutine collision_T()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: Q(0:4)

    do j=1,nyLocal
        do i=1,nxLocal
            n(0) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
            n(1) = g(1,i,j)-g(3,i,j)
            n(2) = g(2,i,j)-g(4,i,j)
            n(3) = g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
            n(4) = g(1,i,j)-g(2,i,j)+g(3,i,j)-g(4,i,j)

            Q(0) = 1.0d0
            Q(1) = sig_k
            Q(2) = sig_k
            Q(3) = 1.9d0
            Q(4) = 1.9d0

            neq(0) = temp(i,j)
            neq(1) = temp(i,j)*u(i,j)
            neq(2) = temp(i,j)*v(i,j)
            neq(3) = temp(i,j)*0.5d0
            neq(4) = 0.0d0

            do alpha=0,4
                n_post(alpha)=n(alpha)-Q(alpha)*(n(alpha)-neq(alpha))
            enddo

            g_post(0,i,j) = n_post(0)-n_post(3)
            g_post(1,i,j) = n_post(1)/2.0d0+n_post(3)/4.0d0+n_post(4)/4.0d0
            g_post(2,i,j) = n_post(2)/2.0d0+n_post(3)/4.0d0-n_post(4)/4.0d0
            g_post(3,i,j) = -n_post(1)/2.0d0+n_post(3)/4.0d0+n_post(4)/4.0d0
            g_post(4,i,j) = -n_post(2)/2.0d0+n_post(3)/4.0d0-n_post(4)/4.0d0

        enddo
    enddo

    do j = 1, nyLocal
        do i = 1, nxLocal
            g(0,i,j) = g_post(0,i,j)
        enddo
    enddo

    return
end subroutine collision_T

subroutine update()
    use commondata
    implicit none
    integer :: i, j, alpha
    real(kind=8) :: sendr_data_f(8*(nyLocal+2)), sendl_data_f(8*(nyLocal+2))
    real(kind=8) :: recvr_data_f(8*(nyLocal+2)), recvl_data_f(8*(nyLocal+2))
    real(kind=8) :: sendr_data_g(4*(nyLocal+2)), sendl_data_g(4*(nyLocal+2))
    real(kind=8) :: recvr_data_g(4*(nyLocal+2)), recvl_data_g(4*(nyLocal+2))

    !Update (top/down)
    !speed distribution function
    call MPI_SENDRECV(f_post(0,0,1),9*(nxLocal+2),MPI_REAL8,down,0,f_post(0,0,nyLocal+1),9*(nxLocal+2),MPI_REAL8,top,0&
    ,COMM2D,ISTATUS,IERR)

    call MPI_SENDRECV(f_post(0,0,nyLocal),9*(nxLocal+2),MPI_REAL8,top,1,f_post(0,0,0),9*(nxLocal+2),MPI_REAL8,down,1&
    ,COMM2D,ISTATUS,IERR)

    !temperature distribution function
    call MPI_SENDRECV(g_post(0,0,1),5*(nxLocal+2),MPI_REAL8,down,0,g_post(0,0,nyLocal+1),5*(nxLocal+2),MPI_REAL8,top,0&
    ,COMM2D,ISTATUS,IERR)

    call MPI_SENDRECV(g_post(0,0,nyLocal),5*(nxLocal+2),MPI_REAL8,top,1,g_post(0,0,0),5*(nxLocal+2),MPI_REAL8,down,1&
    ,COMM2D,ISTATUS,IERR)

    !Update (left/right)
    do j=0,nyLocal+1
        !speed distribution function
        do alpha=1,8
            sendr_data_f(8*j+alpha)=f_post(alpha,nxLocal,j)
            sendl_data_f(8*j+alpha)=f_post(alpha,1,j)
        end do
        !temperature distribution function
        do alpha=1,4
            sendr_data_g(4*j+alpha)=g_post(alpha,nxLocal,j)
            sendl_data_g(4*j+alpha)=g_post(alpha,1,j)
        end do
    end do

    call MPI_SENDRECV(sendr_data_f,8*(nyLocal+2),MPI_REAL8,right,2,recvl_data_f,8*(nyLocal+2),MPI_REAL8,left,2,COMM2D,ISTATUS,IERR)
    call MPI_SENDRECV(sendl_data_f,8*(nyLocal+2),MPI_REAL8,left,3,recvr_data_f,8*(nyLocal+2),MPI_REAL8,right,3,COMM2D,ISTATUS,IERR)

    call MPI_SENDRECV(sendr_data_g,4*(nyLocal+2),MPI_REAL8,right,2,recvl_data_g,4*(nyLocal+2),MPI_REAL8,left,2,COMM2D,ISTATUS,IERR)
    call MPI_SENDRECV(sendl_data_g,4*(nyLocal+2),MPI_REAL8,left,3,recvr_data_g,4*(nyLocal+2),MPI_REAL8,right,3,COMM2D,ISTATUS,IERR)

    do j=0,nyLocal+1
        do alpha=1,8
            f_post(alpha,0,j)=recvl_data_f(8*j+alpha)
            f_post(alpha,nxLocal+1,j)=recvr_data_f(8*j+alpha)
        end do

        do alpha=1,4
            g_post(alpha,0,j)=recvl_data_g(4*j+alpha)
            g_post(alpha,nxLocal+1,j)=recvr_data_g(4*j+alpha)
        end do
    end do
    return
end subroutine

subroutine interpolate()
    use commondata
    implicit none
    real(kind=8) :: interpolateF, delta_x, delta_y
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: f0, f1, f2, g0, g1, g2

        do j = 1, nyLocal
            do i = 1, nxLocal
                do alpha = 1, 8
                    delta_x=dble(ex(alpha))*dt
                    delta_y=dble(ey(alpha))*dt

            f0 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), f_post(alpha,inter_x(i,1), inter_y(j,1)), f_post(alpha,inter_x(i,1), inter_y(j,2))&
                , f_post(alpha,inter_x(i,1), inter_y(j,3)))

            f1 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), f_post(alpha,inter_x(i,2), inter_y(j,1)), f_post(alpha,inter_x(i,2), inter_y(j,2))&
                , f_post(alpha,inter_x(i,2), inter_y(j,3)))

            f2 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), f_post(alpha,inter_x(i,3), inter_y(j,1)), f_post(alpha,inter_x(i,3), inter_y(j,2))&
                , f_post(alpha,inter_x(i,3), inter_y(j,3)))

            f(alpha, i, j) = interpolateF(xGrid(inter_x(i,1))+delta_x, xGrid(inter_x(i,2))+delta_x, &
                            xGrid(inter_x(i,3))+delta_x, xGrid(inter_x(i,2)), f0, f1, f2)

                end do
            enddo
        enddo

        do j = 1, nyLocal
            do i = 1, nxLocal
                do alpha = 1, 4
                    delta_x=dble(ex(alpha))*dt
                    delta_y=dble(ey(alpha))*dt

            g0 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), g_post(alpha,inter_x(i,1), inter_y(j,1)), g_post(alpha,inter_x(i,1), inter_y(j,2))&
                , g_post(alpha,inter_x(i,1), inter_y(j,3)))

            g1 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), g_post(alpha,inter_x(i,2), inter_y(j,1)), g_post(alpha,inter_x(i,2), inter_y(j,2))&
                , g_post(alpha,inter_x(i,2), inter_y(j,3)))

            g2 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), g_post(alpha,inter_x(i,3), inter_y(j,1)), g_post(alpha,inter_x(i,3), inter_y(j,2))&
                , g_post(alpha,inter_x(i,3), inter_y(j,3)))

            g(alpha, i, j) = interpolateF(xGrid(inter_x(i,1))+delta_x, xGrid(inter_x(i,2))+delta_x, &
                            xGrid(inter_x(i,3))+delta_x, xGrid(inter_x(i,2)), g0, g1, g2)
                enddo
            end do
        end do
end subroutine

!!NOTE: consider using compiler-specific directives to suggest inlining if necessary.
pure function interpolateF(x0, x1, x2, x, f0, f1, f2) result(f_interp)
    implicit none
    real(kind=8), intent(in) :: x0, x1, x2, x, f0, f1, f2
    real(kind=8) :: f_interp

    ! Interpolation formula
    f_interp = ((x - x1) * (x - x2)) / ((x0 - x1) * (x0 - x2)) * f0 + &
               ((x - x0) * (x - x2)) / ((x1 - x0) * (x1 - x2)) * f1 + &
               ((x - x0) * (x - x1)) / ((x2 - x0) * (x2 - x1)) * f2

    return
end function interpolateF

subroutine bounceback_u()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    !Left side
    if(MY_COORD(1)==0)then
        do j=1,nyLocal
            f(1,1,j) = f_post(3,1,j)
            f(5,1,j) = f_post(7,1,j)
            f(8,1,j) = f_post(6,1,j)
        end do

    !Right side
    elseif(MY_COORD(1)==dims(1)-1)then
        do j=1,nyLocal
            f(3,nxLocal,j) = f_post(1,nxLocal,j)
            f(6,nxLocal,j) = f_post(8,nxLocal,j)
            f(7,nxLocal,j) = f_post(5,nxLocal,j)
        end do
    end if

    !Top side
    if(MY_COORD(2)==dims(2)-1)then
        do i=1,nxLocal
            f(4,i,nyLocal) = f_post(2,i,nyLocal)
            f(7,i,nyLocal) = f_post(5,i,nyLocal)
            f(8,i,nyLocal) = f_post(6,i,nyLocal)
        end do
    !Bottom side
    elseif(MY_COORD(2)==0)then
        do i=1,nxLocal
            f(2,i,1) = f_post(4,i,1)
            f(5,i,1) = f_post(7,i,1)
            f(6,i,1) = f_post(8,i,1)
        end do
    end if
return
end subroutine bounceback_u

subroutine bounceback_T()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !Left side
    if(MY_COORD(1)==0)then
        do j=1,nyLocal
            g(1,1,j) = -g_post(3,1,j) + 2.0d0 * omega_T(1) * Thot
        end do

    !Right side
    elseif(MY_COORD(1)==dims(1)-1)then
        do j=1,nyLocal
            g(3,nxLocal,j) = -g_post(1,nxLocal,j) + 2.0d0 * omega_T(3) * Tcold
        end do
    end if

    !Top side
    if(MY_COORD(2)==dims(2)-1)then
        do i=1,nxLocal
            g(4,i,nyLocal) = g_post(2,i,nyLocal)
        end do
    !Bottom side
    elseif(MY_COORD(2)==0)then
        do i=1,nxLocal
            g(2,i,1) = g_post(4,i,1)
        end do
    end if

    return
end subroutine

subroutine macro_u()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    do j=1, nyLocal
        do i=1, nxLocal
            rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
            u(i,j) = ( f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)+0.5d0*dt*Fx(i,j))/rho(i,j)
            v(i,j) = ( f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)+0.5d0*dt*Fy(i,j))/rho(i,j)
        enddo
    enddo
    return
end subroutine macro_u

subroutine macro_t()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    do j=1, nyLocal
        do i=1, nxLocal
            temp(i,j) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
        end do
    end do

    return
end subroutine macro_t

subroutine check()
    use commondata
    implicit none
    integer :: i, j
    real(kind=8) :: error1, error2, error1_all, error2_all
    real(kind=8) :: error3, error4, error3_all, error4_all

    error1 = 0.0d0
    error2 = 0.0d0
    error1_all = 0.0d0
    error2_all = 0.0d0

    error3 = 0.0d0
    error4 = 0.0d0
    error3_all = 0.0d0
    error4_all = 0.0d0
    !-------------------------errorU------------------------------------------------
    do j=1,nyLocal
        do i=1,nxLocal
            error1 = error1+(u(i,j)-up(i,j))**2+(v(i,j)-vp(i,j))**2
            error2 = error2+u(i,j)*u(i,j)+v(i,j)*v(i,j)

            up(i,j) = u(i,j)
            vp(i,j) = v(i,j)
        enddo
    enddo

    call MPI_ALLREDUCE(error1,error1_all,1,MPI_REAL8,MPI_SUM,COMM2D,IERR)
    call MPI_ALLREDUCE(error2,error2_all,1,MPI_REAL8,MPI_SUM,COMM2D,IERR)

    errorU = sqrt(error1_all)/sqrt(error2_all)

    !--------------------------errorT------------------------------------------------
    do j=1,nyLocal
        do i=1,nxLocal
            error3 = error3+(temp(i,j)-utemp(i,j))**2
            error4 = error4+temp(i,j)**2

            utemp(i,j) = temp(i,j)
        enddo
    enddo

    call MPI_ALLREDUCE(error3,error3_all,1,MPI_REAL8,MPI_SUM,COMM2D,IERR)
    call MPI_ALLREDUCE(error4,error4_all,1,MPI_REAL8,MPI_SUM,COMM2D,IERR)

    errorT = sqrt(error3_all)/sqrt(error4_all)

    if(MYID==0)then
       write(*,*) itc,' ',errorU,' ',errorT
    end if

    return
end subroutine check

subroutine gather_data()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: iSrc
    real(kind=8), allocatable :: u_temp(:,:), v_temp(:,:), rho_temp(:,:), T_temp(:,:)

    if(MYID==0) then
        do j = 1, nyLocal
            do i = 1, nxLocal
                u_all(i+i_start(0)-1, j+j_start(0)-1) = u(i, j)
                v_all(i+i_start(0)-1, j+j_start(0)-1) = v(i, j)
                rho_all(i+i_start(0)-1, j+j_start(0)-1) = rho(i, j)
                T_all(i+i_start(0)-1, j+j_start(0)-1) = temp(i, j)
            enddo
        enddo

        do iSrc = 1, NPROC-1
            allocate(u_temp(countX(iSrc),countY(iSrc)))
            allocate(v_temp(countX(iSrc),countY(iSrc)))
            allocate(rho_temp(countX(iSrc),countY(iSrc)))
            allocate(T_temp(countX(iSrc),countY(iSrc)))

            call MPI_RECV(u_temp, countX(iSrc)*countY(iSrc), MPI_REAL8, iSrc, 10, COMM2D, ISTATUS, IERR)
            call MPI_RECV(v_temp, countX(iSrc)*countY(iSrc), MPI_REAL8, iSrc, 11, COMM2D, ISTATUS, IERR)
            call MPI_RECV(rho_temp, countX(iSrc)*countY(iSrc), MPI_REAL8, iSrc, 12, COMM2D, ISTATUS, IERR)
            call MPI_RECV(T_temp, countX(iSrc)*countY(iSrc), MPI_REAL8, iSrc, 13, COMM2D, ISTATUS, IERR)

            do j = 1, countY(iSrc)
                do i = 1, countX(iSrc)
                    u_all(i+i_start(iSrc)-1, j+j_start(iSrc)-1) = u_temp(i, j)
                    v_all(i+i_start(iSrc)-1, j+j_start(iSrc)-1) = v_temp(i, j)
                    rho_all(i+i_start(iSrc)-1, j+j_start(iSrc)-1) = rho_temp(i, j)
                    T_all(i+i_start(iSrc)-1, j+j_start(iSrc)-1) = T_temp(i, j)
                enddo
            enddo

            deallocate(u_temp)
            deallocate(v_temp)
            deallocate(rho_temp)
            deallocate(T_temp)
        enddo

    else
        call MPI_SEND(u,countX(MYID)*countY(MYID), MPI_REAL8, 0, 10, COMM2D, IERR)
        call MPI_SEND(v,countX(MYID)*countY(MYID), MPI_REAL8, 0, 11, COMM2D, IERR)
        call MPI_SEND(rho,countX(MYID)*countY(MYID), MPI_REAL8, 0, 12, COMM2D, IERR)
        call MPI_SEND(temp,countX(MYID)*countY(MYID), MPI_REAL8, 0, 13, COMM2D, IERR)
    endif

    return
end subroutine

subroutine output_ASCII()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    character(len=100) :: filename

    write(filename,*) itc
    filename = adjustl(filename)

    open(unit=02,file='MRTcavity-'//trim(filename)//'.dat',status='unknown')

    write(02,*) 'TITLE="thermal convective flows"'
    write(02,*) 'VARIABLES="X" "Y" "U" "V" "T" "rho"'
    write(02,101) nx, ny
    do j=1,ny
        do i=1,nx
            write(02,100) x(i), y(j), u_all(i,j), v_all(i,j), T_all(i,j), rho_all(i,j)
        enddo
    enddo
100 format(1x,2(e11.4,' '),10(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(02)

    return
end subroutine output_ASCII
