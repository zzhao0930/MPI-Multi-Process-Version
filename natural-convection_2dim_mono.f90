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
        real(kind=8) :: xGrid(0:nx+1), y(0:ny+1)
        real(kind=8), allocatable :: yGrid(:)

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
        integer(kind=4) :: NPROC, MYID, IERR,ISTATUS(MPI_STATUS_SIZE)
        integer(kind=4) :: nyLocal
        integer(kind=4) :: upid, downid
        integer(kind=4), allocatable :: start1d(:), end1d(:), count1d(:), displ1d(:)
        integer(kind=4), allocatable :: start2d(:), end2d(:), count2d(:), displ2d(:)
end module commondata

subroutine mesh()
    use commondata
    implicit none

    integer(kind=4) :: i, j
    real(kind=8) :: dx(nx+1), dy(ny+1)
    real(kind=8) :: constA
    allocate(yGrid(0:nyLocal+1))

    if(MYID == 0)then
        ! Compute grid coordinates
        constA = 3.2d0
        do i = 0, nx+1
            xGrid(i) = 0.5d0 * (erf(constA  * (dble(i) / dble(nx+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
        end do
        do j = 0, ny+1
            y(j) = 0.5d0 * (erf(constA  * (dble(j) / dble(ny+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
        end do

        ! Compute grid spacing using array slicing
        dx(1:nx+1) = xGrid(1:nx+1) - xGrid(0:nx)
        dy(1:ny+1) = y(1:ny+1) - y(0:ny)
        dx0 = dx(1)
        dt0 = dx0
        write(*,*) "nx =", nx, ", ny =", ny
        write(*,*) "Ra =", real(Ra)
        write(*,*) "  "
        write(*,*) "---in 0 system---"
        write(*,*) "deltaX = ", dx0
        write(*,*) "deltaT = ", dt0

        length0 = xGrid(nx+1)-xGrid(1)
        write(*,*) "length0 = ", real(length0)
        write(*,*) "  "

        ! Calculate viscosity based on LB unit
        length_LB = 1.0d0 / dx(1)
        dt = dt0 * length_LB

        xGrid(1:nx) = xGrid(1:nx)-dx0/2.0d0
        y(1:ny) = y(1:ny)-dx0/2.0d0
        xGrid=xGrid*length_LB
        y=y*length_LB
    endif

    call MPI_BCAST(dt,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    call MPI_BCAST(length_LB,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    call MPI_BCAST(xGrid(1),nx,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    call MPI_SCATTERV(y(1),count1d,displ1d,MPI_REAL8,yGrid(1),count1d(MYID),MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    call MPI_SENDRECV(yGrid(1),1,MPI_REAL8,downid,0,yGrid(nyLocal+1),1,MPI_REAL8,upid,0&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    call MPI_SENDRECV(yGrid(nyLocal),1,MPI_REAL8,upid,1,yGrid(0),1,MPI_REAL8,downid,1&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    return
end subroutine

subroutine initial()
    use commondata
    implicit none

    integer(kind=4) :: i, j, alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: velocitySquare
    allocate(inter_x(nx,3))
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
        ! Output initial parameters for verification
        write(*,*) "---in system LB"
        write(*,*) "deltaX = ", xGrid(1)-xGrid(0)
        write(*,*) "deltaT = ", dt
        write(*,*) "characteristic length   =", real(length_LB), "l.u."
        write(*,*) "viscosity_LB =", real(viscosity_LB), "l.u.^2/t.s."
        write(*,*) "timeStep ratio for (uniform) / (non-uniform) : ", real(length_LB / dble(nx))
        write(*,*) "tauf =", real(tauf)
        write(*,*) "    "
    end if

    ! Allocate flow variables
    allocate (u(nx,nyLocal))
    allocate (v(nx,nyLocal))
    allocate (rho(nx,nyLocal))
    allocate (up(nx,nyLocal))
    allocate (vp(nx,nyLocal))
    allocate (temp(nx,nyLocal))
    allocate (utemp(nx,nyLocal))
    allocate (Fx(nx,nyLocal))
    allocate (Fy(nx,nyLocal))

    if(MYID == 0)then
        allocate (u_all(nx,ny))
        allocate (v_all(nx,ny))
        allocate (rho_all(nx,ny))
        allocate (T_all(nx,ny))
    end if

    allocate (f(0:8,nx,nyLocal))
    allocate (f_post(0:8,nx,0:nyLocal+1))
    allocate (g(0:4,nx,nyLocal))
    allocate (g_post(0:4,nx,0:nyLocal+1))

    ! Initialize flow variables
    rho = rho0
    temp = 0.0d0
    utemp=0.0d0
    u = 0.0d0
    v = 0.0d0
    up = 0.0d0
    vp = 0.0d0

    do j=1,nyLocal
        do i=1,nx
            temp(i,j) = dble(i-1)/dble(nx-1)*(Tcold-Thot)+Thot
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
        do i = 1, nx
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

    do i = 1, nx
        if(i == 1)then
            inter_x(i,:) = (/i+1, i, i+2/)

        elseif(i == nx)then
            inter_x(i,:) = (/i-1, i, i-2/)

        else
            inter_x(i,:) = (/i-1, i, i+1/)
        end if
    enddo

    if(MYID == 0)then
        do j = 1, nyLocal
            if(j == 1)then
                inter_y(j,:) = (/j+1, j, j+2/)

            else
                inter_y(j,:) = (/j-1, j, j+1/)
            end if
        end do

    elseif(MYID == NPROC-1)then
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
    real(8) :: timestart, timeEnd
    !----------------------------------------------------------
    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)

    upid=MYID+1
    downid=MYID-1
    if(MYID == 0)downid=MPI_PROC_NULL
    if(MYID == NPROC-1)upid=MPI_PROC_NULL
    !----------------------------------------------------------
    call StartEnd(1, ny)
    nyLocal = count1d(MYID)

    call MPI_BARRIER(MPI_COMM_WORLD, IERR)
    timestart = MPI_WTIME()

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
    deallocate(yGrid)
    deallocate(inter_y)

    if(MYID == 0)then
        deallocate(u_all)
        deallocate(v_all)
        deallocate(rho_all)
        deallocate(T_all)
    end if

    deallocate(start1d)
    deallocate(end1d)
    deallocate(count1d)
    deallocate(displ1d)
    deallocate(start2d)
    deallocate(end2d)
    deallocate(count2d)
    deallocate(displ2d)

    call MPI_FINALIZE(IERR)

    stop
end program main

subroutine StartEnd(iS1, iS2)
    use commondata
    implicit none
    integer(kind=4):: leng, iBlock
    integer(kind=4) :: ir
    integer(kind=4) :: iS1, iS2
    integer(kind=4) :: i
    allocate (start1d(0:NPROC-1))
    allocate (end1d(0:NPROC-1))
    allocate (count1d(0:NPROC-1))
    allocate (displ1d(0:NPROC-1))

    allocate (start2d(0:NPROC-1))
    allocate (end2d(0:NPROC-1))
    allocate (count2d(0:NPROC-1))
    allocate (displ2d(0:NPROC-1))

    leng = iS2-iS1+1
    iBlock = leng/NPROC
    ir= leng-iBlock*NPROC

    do i = 0, NPROC-1

        if(i.LT.ir) then
            count1d(i) = iBlock+1
            start1d(i) = iS1+i*(iBlock+1)
            end1d(i) = start1d(i)+count1d(i)-1
            !-----------------------------------------------------------
            count2d(i) = (iBlock+1)*nx
            start2d(i) = iS1+i*(iBlock+1)*nx
            end2d(i) = start2d(i)+count2d(i)-1

        else
            count1d(i) = iBlock
            start1d(i) = iS1+i*iBlock+ir
            end1d(i) = start1d(i)+count1d(i)-1
            !-----------------------------------------------------------
            count2d(i) = iBlock*nx
            start2d(i) = iS1+i*iBlock*nx+ir*nx
            end2d(i) = start2d(i)+count2d(i)-1
        endif

        displ1d(i) = start1d(i)-iS1
        displ2d(i) = start2d(i)-iS1

    enddo
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
        do i=1,nx

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
        do i = 1, nx
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
        do i=1,nx
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
        do i = 1, nx
            g(0,i,j) = g_post(0,i,j)
        enddo
    enddo

    return
end subroutine collision_T


subroutine update()
    use commondata
    implicit none

    call MPI_SENDRECV(f_post(0,1,1),9*nx,MPI_REAL8,downid,0,f_post(0,1,nyLocal+1),9*nx,MPI_REAL8,upid,0&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    call MPI_SENDRECV(f_post(0,1,nyLocal),9*nx,MPI_REAL8,upid,1,f_post(0,1,0),9*nx,MPI_REAL8,downid,1&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    call MPI_SENDRECV(g_post(0,1,1),5*nx,MPI_REAL8,downid,0,g_post(0,1,nyLocal+1),5*nx,MPI_REAL8,upid,0&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    call MPI_SENDRECV(g_post(0,1,nyLocal),5*nx,MPI_REAL8,upid,1,g_post(0,1,0),5*nx,MPI_REAL8,downid,1&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    return
end subroutine

subroutine interpolate()
    use commondata
    implicit none
    real(kind=8) :: interpolateF, delta_x, delta_y
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: f0, f1, f2, g0, g1, g2

        do j = 1, nyLocal
            do i = 1, nx
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
            do i = 1, nx
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

    do j=1,nyLocal
        !Left side
        f(1,1,j) = f_post(3,1,j)
        f(5,1,j) = f_post(7,1,j)
        f(8,1,j) = f_post(6,1,j)

        !Right side
        f(3,nx,j) = f_post(1,nx,j)
        f(6,nx,j) = f_post(8,nx,j)
        f(7,nx,j) = f_post(5,nx,j)
    enddo

    if(MYID==NPROC-1)then
        do i=1,nx
            !Top side
            f(4,i,nyLocal) = f_post(2,i,nyLocal)
            f(7,i,nyLocal) = f_post(5,i,nyLocal)
            f(8,i,nyLocal) = f_post(6,i,nyLocal)
        end do

    elseif(MYID==0)then
        do i=1,nx
            !Bottom side
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

    do j=1, nyLocal
        !left
        g(1,1,j) = -g_post(3,1,j) + 2.0d0 * omega_T(1) * Thot
        !!right
        g(3,nx,j) = -g_post(1,nx,j) + 2.0d0 * omega_T(3) * Tcold
    end do

    if(MYID==NPROC-1)then
        do i=1,nx
            !Top side
            g(4,i,nyLocal) = g_post(2,i,nyLocal)
        end do

    elseif(MYID==0)then
        do i=1,nx
            !Bottom side
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
        do i=1, nx
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
        do i=1, nx
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

    do j=1,nyLocal
        do i=1,nx
            error1 = error1+(u(i,j)-up(i,j))*(u(i,j)-up(i,j))+(v(i,j)-vp(i,j))*(v(i,j)-vp(i,j))
            error2 = error2+u(i,j)*u(i,j)+v(i,j)*v(i,j)

            up(i,j) = u(i,j)
            vp(i,j) = v(i,j)
        enddo
    enddo

    call MPI_ALLREDUCE(error1,error1_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
    call MPI_ALLREDUCE(error2,error2_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)

    errorU = sqrt(error1_all)/sqrt(error2_all)

    do j=1,nyLocal
        do i=1,nx
            error3 = error3+(temp(i,j)-utemp(i,j))**2
            error4 = error4+temp(i,j)**2

            utemp(i,j) = temp(i,j)
        enddo
    enddo

    call MPI_ALLREDUCE(error3,error3_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
    call MPI_ALLREDUCE(error4,error4_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)

    errorT = sqrt(error3_all)/sqrt(error4_all)

    if(MYID==0)then
       write(*,*) itc,' ',errorU,' ',errorT
    end if

    return
end subroutine check

subroutine gather_data()
    use commondata
    implicit none

    call MPI_GATHERV(u,count2d(MYID),MPI_REAL8,u_all,count2d,displ2d,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    call MPI_GATHERV(v,count2d(MYID),MPI_REAL8,v_all,count2d,displ2d,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    call MPI_GATHERV(rho,count2d(MYID),MPI_REAL8,rho_all,count2d,displ2d,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    call MPI_GATHERV(temp,count2d(MYID),MPI_REAL8,T_all,count2d,displ2d,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
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
            write(02,100) xGrid(i), y(j), u_all(i,j), v_all(i,j), T_all(i,j), rho_all(i,j)
        enddo
    enddo
100 format(1x,2(e11.4,' '),10(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(02)

    return
end subroutine output_ASCII
