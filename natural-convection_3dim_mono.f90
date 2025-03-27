module commondata
        use mpi
        implicit none

        ! Grid parameters
        integer(kind=4), parameter :: nx=257, ny=nx, nz=nx
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1, nzHalf=(nz-1)/2+1

        real(kind=8), parameter :: Ra=1e6          !Rayleigh number
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
        real(kind=8) :: tauf                     ! Relaxation time
        real(kind=8) :: dt                       ! Time step

        ! Iteration control
        integer(kind=4) :: itc = 0            ! Current iteration count
        integer(kind=4), parameter :: itc_max = 2000 ! Maximum iterations

        ! Convergence criteria
        real(kind=8) :: errorU, errorT              ! Current error
        real(kind=8), parameter :: epsU=1e-6, epsT=1e-6  ! Convergence threshold

        ! Grid coordinates
        real(kind=8) :: xGrid(0:nx+1), yGrid(0:ny+1), z(0:nz+1)
        real(kind=8), allocatable :: zGrid(:)

        ! Flow variables
        real(kind=8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)      ! Velocity components
        real(kind=8), allocatable :: temp(:,:,:)           ! temperature
        real(kind=8), allocatable :: rho(:,:,:)            ! Density
        real(kind=8), allocatable :: u_all(:,:,:), v_all(:,:,:), w_all(:,:,:)      ! Velocity components
        real(kind=8), allocatable :: T_all(:,:,:)           ! temperature
        real(kind=8), allocatable :: rho_all(:,:,:)            ! Density
        real(kind=8), allocatable :: up(:,:,:), vp(:,:,:), wp(:,:,:)    ! Previous velocity components for error checking
        real(kind=8), allocatable :: utemp(:,:,:)          ! Previous temperature for error checking
        real(kind=8), allocatable :: Fz(:,:,:)    ! force components

        ! Distribution functions
        real(kind=8), allocatable :: f(:,:,:,:), f_post(:,:,:,:)  ! Current and post-collision distributions
        real(kind=8),allocatable :: g(:,:,:,:), g_post(:,:,:,:)

        ! MRT relaxation parameters
        real(kind=8) :: omega_U(0:18), omega_T(0:6)  ! Relaxation rates for MRT

        ! Lattice directions
        integer(kind=4) :: ex(0:18), ey(0:18), ez(0:18)  ! Lattice velocity directions
        data ex/ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0/
        data ey/ 0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1/
        data ez/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1, -1,-1, 1, 1,-1,-1/

        ! Additional MRT parameters
        real(kind=8) :: Snu, Sq, sig_k

        integer(kind=4), allocatable :: inter_x(:,:), inter_y(:,:), inter_z(:,:)

    !----------------------------MPI---------------------------------------------------
        integer(kind=4) :: NPROC, MYID, IERR,ISTATUS(MPI_STATUS_SIZE)
        integer(kind=4) :: nzLocal
        integer(kind=4) :: upid, downid
        integer(kind=4), allocatable :: start1d(:), count1d(:), displ1d(:)
        integer(kind=4), allocatable :: start2d(:), count2d(:), displ2d(:)

end module commondata

subroutine StartEnd(iS1, iS2)
    use commondata
    implicit none
    integer(kind=4) :: leng, iBlock
    integer(kind=4) :: ir
    integer(kind=4) :: iS1, iS2
    integer(kind=4) :: i
    allocate (start1d(0:NPROC-1))
    allocate (count1d(0:NPROC-1))
    allocate (displ1d(0:NPROC-1))

    allocate (start2d(0:NPROC-1))
    allocate (count2d(0:NPROC-1))
    allocate (displ2d(0:NPROC-1))

    leng = iS2-iS1+1
    iBlock = leng/NPROC
    ir= leng-iBlock*NPROC

    do i=0,NPROC-1

        if(i < ir) then
            count1d(i) = iBlock+1
            start1d(i) = iS1+i*(iBlock+1)
            !-----------------------------------------------------------
            count2d(i) = (iBlock+1)*nx*ny
            start2d(i) = iS1+i*(iBlock+1)*nx*ny

        else
            count1d(i) = iBlock
            start1d(i) = iS1+i*iBlock+ir
            !-----------------------------------------------------------
            count2d(i) = iBlock*nx*ny
            start2d(i) = iS1+i*iBlock*nx*ny+ir*nx*ny
        endif

        displ1d(i) = start1d(i)-iS1
        displ2d(i) = start2d(i)-iS1

    enddo
    return
end subroutine StartEnd

subroutine mesh()
    use commondata
    implicit none

    integer(kind=4) :: i, j, k
    real(kind=8) :: dx(nx+1)
    real(kind=8) :: constA
    allocate(zGrid(0:nzLocal+1))

    if(MYID == 0)then
        ! Compute grid coordinates
        constA = 3.2d0
        do i = 0, nx+1
            xGrid(i) = 0.5d0 * (erf(constA  * (dble(i) / dble(nx+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
        end do

        do j = 0, ny+1
            yGrid(j) = 0.5d0 * (erf(constA  * (dble(j) / dble(ny+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
        end do

        do k = 0, nz+1
            z(k) = 0.5d0 * (erf(constA  * (dble(k) / dble(nz+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
        end do
        ! Compute grid spacing using array slicing
        dx(1:nx+1) = xGrid(1:nx+1) - xGrid(0:nx)
        dx0 = dx(1)
        dt0 = dx0
        write(*,*) "nx =", nx, ", ny =", ny, ", nz =", nz
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
        yGrid(1:ny) = yGrid(1:ny)-dx0/2.0d0
        z(1:nz) = z(1:nz)-dx0/2.0d0
        xGrid=xGrid*length_LB
        yGrid=yGrid*length_LB
        z=z*length_LB
    endif

    call MPI_BCAST(dt,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    call MPI_BCAST(length_LB,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    call MPI_BCAST(xGrid(1),nx,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    call MPI_BCAST(yGrid(1),ny,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    call MPI_SCATTERV(z(1),count1d,displ1d,MPI_REAL8,zGrid(1),count1d(MYID),MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    call MPI_SENDRECV(zGrid(1),1,MPI_REAL8,downid,0,zGrid(nzLocal+1),1,MPI_REAL8,upid,0&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    call MPI_SENDRECV(zGrid(nzLocal),1,MPI_REAL8,upid,1,zGrid(0),1,MPI_REAL8,downid,1&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    return
end subroutine

subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j, K, alpha
    real(kind=8) :: un(0:18)
    real(kind=8) :: us2
    allocate(inter_x(nx,3))
    allocate(inter_y(ny,3))
    allocate(inter_z(nzLocal,3))

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

    if(MYID == 0)then
        write(*,*) "---in LB unit---"
        write(*,*) "characteristic length   =", real(length_LB), "l.u."
        write(*,*) "viscosity_LB =", real(viscosity_LB), "l.u.^2/t.s."
        write(*,*) "timeStep ratio for (uniform) / (non-uniform) : ", real(length_LB / dble(nx))
        write(*,*) "    "

        write(*,*) "tauf =", real(tauf)
    end if

    ! Calculate MRT relaxation parameters
    Snu = 1.0d0/tauf

    Sq = 8.0d0*(2.0d0-Snu)/(8.0d0-Snu)

    sig_k = 1.0d0/(0.5d0+(1.0d0/Snu-0.5d0)/(0.5d0*Pr))

    ! Allocate flow variables
    allocate (u(nx,ny,nzLocal))
    allocate (v(nx,ny,nzLocal))
    allocate (w(nx,ny,nzLocal))
    allocate (rho(nx,ny,nzLocal))
    allocate (up(nx,ny,nzLocal))
    allocate (vp(nx,ny,nzLocal))
    allocate (wp(nx,ny,nzLocal))
    allocate (temp(nx,ny,nzLocal))
    allocate (utemp(nx,ny,nzLocal))
    allocate (Fz(nx,ny,nzLocal))

    if(MYID == 0)then
        allocate (u_all(nx,ny,nz))
        allocate (v_all(nx,ny,nz))
        allocate (w_all(nx,ny,nz))
        allocate (rho_all(nx,ny,nz))
        allocate (T_all(nx,ny,nz))
    end if

    allocate (f(0:18,nx,ny,nzLocal))
    allocate (f_post(0:18,nx,ny,0:nzLocal+1))
    allocate (g(0:6,nx,ny,nzLocal))
    allocate (g_post(0:6,nx,ny,0:nzLocal+1))

    ! Initialize flow variables
    rho = rho0
    temp = 0.0d0
    utemp=0.0d0
    u = 0.0d0
    v = 0.0d0
    up = 0.0d0
    vp = 0.0d0

    do k=1,nzLocal
        do j=1,ny
            do i=1,nx
                temp(i,j,k) = dble(j-1)/dble(ny-1)*(Tcold-Thot)+Thot
            enddo
        enddo
    end do

    omega_U(0)=1.0d0/3.0d0
    do alpha=1,6
        omega_U(alpha)=1.0d0/18.0d0
    end do
    do alpha=7,18
        omega_U(alpha)=1.0d0/36.0d0
    end do

    omega_T(0) = 1.0d0/2.0d0
    do alpha=1,6
        omega_T(alpha) = 1.0d0/12.0d0
    enddo

    do k=1,nzLocal
        do j=1,ny
            do i=1,nx
                us2 = u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)
                do alpha=0,18
                    un(alpha) = u(i,j,k)*ex(alpha)+v(i,j,k)*ey(alpha)+w(i,j,k)*ez(alpha)
                    f(alpha,i,j,k) = omega_U(alpha)*rho(i,j,k)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
                enddo

                do alpha=0,6
                    un(alpha) = u(i,j,k)*ex(alpha)+v(i,j,k)*ey(alpha)+w(i,j,k)*ez(alpha)
                    g(alpha,i,j,k)=omega_T(alpha)*temp(i,j,k)*(1.0d0+6.0d0*un(alpha))
                end do
            enddo
        enddo
    end do

    do i = 1, nx
        if(i == 1)then
            inter_x(i,:) = (/i+1, i, i+2/)
            inter_y(i,:) = (/i+1, i, i+2/)

        elseif(i == nx)then
            inter_x(i,:) = (/i-1, i, i-2/)
            inter_y(i,:) = (/i-1, i, i-2/)

        else
            inter_x(i,:) = (/i-1, i, i+1/)
            inter_y(i,:) = (/i-1, i, i+1/)
        end if
    enddo

    if(MYID == 0)then
        do k = 1, nzLocal
            if(k == 1)then
                inter_z(k,:) = (/k+1, k, k+2/)

            else
                inter_z(k,:) = (/k-1, k, k+1/)
            end if
        end do

    elseif(MYID == NPROC-1)then
        do k = 1, nzLocal
            if(k == nzLocal)then
                inter_z(k,:) = (/k-1, k, k-2/)

            else
                inter_z(k,:) = (/k-1, k, k+1/)
            end if
        end do

    else
        do k = 1, nzLocal
            inter_z(k,:) = (/k-1, k, k+1/)
        end do
    end if

    return
end subroutine initial

program main
    use commondata
    implicit none

    real(8) :: timestart, timeEnd
!--------------------------------------------------------------------
    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)

    upid=MYID+1
    downid=MYID-1
    if(MYID == 0)downid=MPI_PROC_NULL
    if(MYID == NPROC-1)upid=MPI_PROC_NULL

    call StartEnd(1, nz)
    nzLocal=count1d(MYID)
!--------------------------------------------------------------------

    call MPI_BARRIER(MPI_COMM_WORLD, IERR)
    timestart = MPI_WTIME()

    call mesh()

    call initial()

    do while(((errorU > epsU).or.(errorT > epsT)).AND.(itc < itc_max))

        itc = itc+1

        call collision_U()

        call collision_T()

        call update()

        call interploate()

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
        write(*,*) "Time (CPU) =", real(timeEnd-timestart),"s"
        write(*,*) "MLUPS = ", real(dble(nx*ny*nz)/1e6*dble(itc)/(timeEnd-timeStart))
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
    deallocate(Fz)
    deallocate(zGrid)
    if(MYID == 0)then
        deallocate (u_all)
        deallocate (v_all)
        deallocate (w_all)
        deallocate (rho_all)
        deallocate (T_all)
    end if
    deallocate (start1d)
    deallocate (count1d)
    deallocate (displ1d)

    deallocate (start2d)
    deallocate (count2d)
    deallocate (displ2d)

    call MPI_FINALIZE(IERR)
    stop
end program main

subroutine collision_U()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
    integer(kind=4) :: alpha
    real(kind=8) :: m(0:18), m_post(0:18), meq(0:18)
    real(kind=8) :: s(0:18)
    real(kind=8) :: fSource(0:18)
    real(kind=8) :: Vsquare

    do k=1,nzLocal
        do j=1,ny
            do i=1,nx
                Vsquare=u(i,j,k)**2.0d0+v(i,j,k)**2.0d0+w(i,j,k)**2.0d0

                m(0)=f(0,i,j,k)+f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+&
                    &f(6,i,j,k)+f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+&
                    &f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)&
                    &+f(18,i,j,k)

                m(1)=-30.0d0*f(0,i,j,k)-11.0d0*(f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+&
                    &f(5,i,j,k)+f(6,i,j,k))+8.0d0*(f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+&
                    &f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+&
                    &f(17,i,j,k)+f(18,i,j,k))

                m(2)=12.0d0*f(0,i,j,k)-4.0d0*(f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+&
                    &f(5,i,j,k)+f(6,i,j,k))+(f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+&
                    &f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+&
                    &f(17,i,j,k)+f(18,i,j,k))

                m(3)=f(1,i,j,k)-f(2,i,j,k)+f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k)+f(11,i,j,k)&
                    &-f(12,i,j,k)+f(13,i,j,k)-f(14,i,j,k)

                m(4)=4.0d0*(-f(1,i,j,k)+f(2,i,j,k))+f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k)&
                    &+f(11,i,j,k)-f(12,i,j,k)+f(13,i,j,k)-f(14,i,j,k)

                m(5)=f(3,i,j,k)-f(4,i,j,k)+f(7,i,j,k)+f(8,i,j,k)-f(9,i,j,k)-f(10,i,j,k)+f(15,i,j,k)&
                    &-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k)

                m(6)=-4.0d0*(f(3,i,j,k)-f(4,i,j,k))+f(7,i,j,k)+f(8,i,j,k)-f(9,i,j,k)-f(10,i,j,k)+f(15,i,j,k)&
                    &-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k)

                m(7)=f(5,i,j,k)-f(6,i,j,k)+f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k)&
                    &+f(15,i,j,k)+f(16,i,j,k)-f(17,i,j,k)-f(18,i,j,k)

                m(8)=-4.0d0*(f(5,i,j,k)-f(6,i,j,k))+f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k)&
                    &+f(15,i,j,k)+f(16,i,j,k)-f(17,i,j,k)-f(18,i,j,k)

                m(9)=2.0d0*(f(1,i,j,k)+f(2,i,j,k))-(f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k))&
                    &+f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)&
                    &+f(14,i,j,k)-2.0d0*(f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k))

                m(10)=-4.0d0*(f(1,i,j,k)+f(2,i,j,k))+2.0d0*(f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k))&
                    &+f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)&
                    &+f(14,i,j,k)-2.0d0*(f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k))

                m(11)=f(3,i,j,k)+f(4,i,j,k)-f(5,i,j,k)-f(6,i,j,k)+f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)&
                    &+f(10,i,j,k)-f(11,i,j,k)-f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k)

                m(12)=-2.0d0*(f(3,i,j,k)+f(4,i,j,k)-f(5,i,j,k)-f(6,i,j,k))+f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)&
                     &+f(10,i,j,k)-f(11,i,j,k)-f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k)

                m(13)=f(7,i,j,k)-f(8,i,j,k)-f(9,i,j,k)+f(10,i,j,k)

                m(14)=f(15,i,j,k)-f(16,i,j,k)-f(17,i,j,k)+f(18,i,j,k)

                m(15)=f(11,i,j,k)-f(12,i,j,k)-f(13,i,j,k)+f(14,i,j,k)

                m(16)=f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k)&
                    &-f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)+f(14,i,j,k)

                m(17)=-f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)&
                    &+f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k)

                m(18)=f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k)&
                    &-f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k)

                meq(0)=rho(i,j,k)
                meq(1)=-11.0d0*rho(i,j,k)+19.0d0*rho(i,j,k)*Vsquare
                meq(2)=3.0d0*rho(i,j,k)-5.5d0*rho(i,j,k)*Vsquare
                meq(3)=rho(i,j,k)*u(i,j,k)
                meq(4)=-2.0d0*rho(i,j,k)*u(i,j,k)/3.0d0
                meq(5)=rho(i,j,k)*v(i,j,k)
                meq(6)=-2.0d0*rho(i,j,k)*v(i,j,k)/3.0d0
                meq(7)=rho(i,j,k)*w(i,j,k)
                meq(8)=-2.0d0*rho(i,j,k)*w(i,j,k)/3.0d0
                meq(9)=3.0d0*rho(i,j,k)*u(i,j,k)**2.0d0-rho(i,j,k)*Vsquare
                meq(10)=-0.5d0*(3.0d0*rho(i,j,k)*u(i,j,k)**2.0d0-rho(i,j,k)*Vsquare)
                meq(11)=rho(i,j,k)*(v(i,j,k)**2.0d0-w(i,j,k)**2.0d0)
                meq(12)=-0.5d0*rho(i,j,k)*(v(i,j,k)**2.0d0-w(i,j,k)**2.0d0)
                meq(13)=rho(i,j,k)*u(i,j,k)*v(i,j,k)
                meq(14)=rho(i,j,k)*v(i,j,k)*w(i,j,k)
                meq(15)=rho(i,j,k)*u(i,j,k)*w(i,j,k)
                meq(16)=0.0d0
                meq(17)=0.0d0
                meq(18)=0.0d0

                s(0)=0.0d0 !-----------S_rho
                s(1)=Snu !------------S_e
                s(2)=Snu !------------S_vapersilon
                s(3)=0.0d0 !------------S_j
                s(4)=Sq !---------------S_q
                s(5)=0.0d0 !----------------S_j
                s(6)=Sq !-----------------S_q
                s(7)=0.0d0 !----------------S_j
                s(8)=Sq !-----------------S_q
                s(9)=Snu !---------------S_nu
                s(10)=Snu !--------------S_Pi
                s(11)=Snu !--------------S_nu
                s(12)=Snu !--------------S_Pi
                s(13)=Snu !--------------S_nu
                s(14)=Snu !--------------S_nu
                s(15)=Snu !--------------S_nu
                s(16)=Sq !----------------S_m
                s(17)=Sq !----------------S_m
                s(18)=Sq !----------------S_m


            Fz(i,j,k)=gbeta*temp(i,j,k)*rho(i,j,k)

            fSource(0) = 0.0d0
            fSource(1) = (38.0d0-19.0d0*s(1))*w(i,j,k)*Fz(i,j,k)
            fSource(2) = -(11.0d0-5.5d0*s(2))*w(i,j,k)*Fz(i,j,k)
            fSource(3) = 0.0d0
            fSource(4) = 0.0d0
            fSource(5) = 0.0d0
            fSource(6) = 0.0d0
            fSource(7) = (1.0d0-0.5d0*s(7))*Fz(i,j,k)
            fSource(8) = -2.0d0*(1.0d0-0.5d0*s(8))*Fz(i,j,k)/3.0d0
            fSource(9) = (-2.0d0+s(9))*w(i,j,k)*Fz(i,j,k)
            fSource(10) = (1.0d0-0.5d0*s(10))*w(i,j,k)*Fz(i,j,k)
            fSource(11) = -(2.0d0-s(11))*w(i,j,k)*Fz(i,j,k)
            fSource(12) = (1.0d0-0.5d0*s(12))*w(i,j,k)*Fz(i,j,k)
            fSource(13) = 0.0d0
            fSource(14) = (1.0d0-0.5d0*s(14))*v(i,j,k)*Fz(i,j,k)
            fSource(15) = (1.0d0-0.5d0*s(15))*u(i,j,k)*Fz(i,j,k)
            fSource(16) = 0.0d0
            fSource(17) = 0.0d0
            fSource(18) = 0.0d0

            do alpha=0,18
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)*dt
            enddo

                f_post(0,i,j,k)=m_post(0)/19.0d0-5.0d0*m_post(1)/399.0d0+m_post(2)/21.0d0

                f_post(1,i,j,k)=m_post(0)/19.0d0-11.0d0*m_post(1)/2394.0d0-m_post(2)/63.0d0+m_post(3)/10.0d0&
                    &-m_post(4)/10.0d0+m_post(9)/18.0d0-m_post(10)/18.0d0

                f_post(2,i,j,k)=m_post(0)/19.0d0-11.0d0*m_post(1)/2394.0d0-m_post(2)/63.0d0-m_post(3)/10.0d0&
                    &+m_post(4)/10.0d0+m_post(9)/18.0d0-m_post(10)/18.0d0

                f_post(3,i,j,k)=m_post(0)/19.0d0-11.0d0*m_post(1)/2394.0d0-m_post(2)/63.0d0+m_post(5)/10.0d0&
                    &-m_post(6)/10.0d0-m_post(9)/36.0d0+m_post(10)/36.0d0+m_post(11)/12.0d0-m_post(12)/12.0d0

                f_post(4,i,j,k)=m_post(0)/19.0d0-11.0d0*m_post(1)/2394.0d0-m_post(2)/63.0d0-m_post(5)/10.0d0&
                    &+m_post(6)/10.0d0-m_post(9)/36.0d0+m_post(10)/36.0d0+m_post(11)/12.0d0-m_post(12)/12.0d0

                f_post(5,i,j,k)=m_post(0)/19.0d0-11.0d0*m_post(1)/2394.0d0-m_post(2)/63.0d0+m_post(7)/10.0d0&
                    &-m_post(8)/10.0d0-m_post(9)/36.0d0+m_post(10)/36.0d0-m_post(11)/12.0d0+m_post(12)/12.0d0

                f_post(6,i,j,k)=m_post(0)/19.0d0-11.0d0*m_post(1)/2394.0d0-m_post(2)/63.0d0-m_post(7)/10.0d0&
                    &+m_post(8)/10.0d0-m_post(9)/36.0d0+m_post(10)/36.0d0-m_post(11)/12.0d0+m_post(12)/12.0d0

                f_post(7,i,j,k)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0+m_post(3)/10.0d0&
                    &+m_post(4)/40.0d0+m_post(5)/10.0d0+m_post(6)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &+m_post(11)/12.0d0+m_post(12)/24.0d0+m_post(13)/4.0d0+m_post(16)/8.0d0-m_post(17)/8.0d0

                f_post(8,i,j,k)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0-m_post(3)/10.0d0&
                    &-m_post(4)/40.0d0+m_post(5)/10.0d0+m_post(6)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &+m_post(11)/12.0d0+m_post(12)/24.0d0-m_post(13)/4.0d0-m_post(16)/8.0d0-m_post(17)/8.0d0

                f_post(9,i,j,k)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0+m_post(3)/10.0d0&
                    &+m_post(4)/40.0d0-m_post(5)/10.0d0-m_post(6)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &+m_post(11)/12.0d0+m_post(12)/24.0d0-m_post(13)/4.0d0+m_post(16)/8.0d0+m_post(17)/8.0d0

                f_post(10,i,j,k)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0-m_post(3)/10.0d0&
                    &-m_post(4)/40.0d0-m_post(5)/10.0d0-m_post(6)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &+m_post(11)/12.0d0+m_post(12)/24.0d0+m_post(13)/4.0d0-m_post(16)/8.0d0+m_post(17)/8.0d0

                f_post(11,i,j,k)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0+m_post(3)/10.0d0&
                    &+m_post(4)/40.0d0+m_post(7)/10.0d0+m_post(8)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &-m_post(11)/12.0d0-m_post(12)/24.0d0+m_post(15)/4.0d0-m_post(16)/8.0d0+m_post(18)/8.0d0

                f_post(12,i,j,k)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0-m_post(3)/10.0d0&
                    &-m_post(4)/40.0d0+m_post(7)/10.0d0+m_post(8)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &-m_post(11)/12.0d0-m_post(12)/24.0d0-m_post(15)/4.0d0+m_post(16)/8.0d0+m_post(18)/8.0d0

                f_post(13,i,j,k)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0+m_post(3)/10.0d0&
                    &+m_post(4)/40.0d0-m_post(7)/10.0d0-m_post(8)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &-m_post(11)/12.0d0-m_post(12)/24.0d0-m_post(15)/4.0d0-m_post(16)/8.0d0-m_post(18)/8.0d0

                f_post(14,i,j,k)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0-m_post(3)/10.0d0&
                    &-m_post(4)/40.0d0-m_post(7)/10.0d0-m_post(8)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &-m_post(11)/12.0d0-m_post(12)/24.0d0+m_post(15)/4.0d0+m_post(16)/8.0d0-m_post(18)/8.0d0

                f_post(15,i,j,k)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0+m_post(5)/10.0d0&
                    &+m_post(6)/40.0d0+m_post(7)/10.0d0+m_post(8)/40.0d0-m_post(9)/18.0d0-m_post(10)/36.0d0&
                    &+m_post(14)/4.0d0+m_post(17)/8.0d0-m_post(18)/8.0d0

                f_post(16,i,j,k)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0-m_post(5)/10.0d0&
                    &-m_post(6)/40.0d0+m_post(7)/10.0d0+m_post(8)/40.0d0-m_post(9)/18.0d0-m_post(10)/36.0d0&
                    &-m_post(14)/4.0d0-m_post(17)/8.0d0-m_post(18)/8.0d0

                f_post(17,i,j,k)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0+m_post(5)/10.0d0&
                    &+m_post(6)/40.0d0-m_post(7)/10.0d0-m_post(8)/40.0d0-m_post(9)/18.0d0-m_post(10)/36.0d0&
                    &-m_post(14)/4.0d0+m_post(17)/8.0d0+m_post(18)/8.0d0

                f_post(18,i,j,k)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0-m_post(5)/10.0d0&
                    &-m_post(6)/40.0d0-m_post(7)/10.0d0-m_post(8)/40.0d0-m_post(9)/18.0d0-m_post(10)/36.0d0&
                    &+m_post(14)/4.0d0-m_post(17)/8.0d0+m_post(18)/8.0d0
            end do
        end do
    end do

    do k = 1, nzLocal
       do j = 1, ny
            do i = 1, nx
                f(0,i,j,k)=f_post(0,i,j,k)
            enddo
        enddo
    end do

    return
end subroutine collision_U

subroutine collision_T()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k, alpha
    real(kind=8) :: n(0:6), n_post(0:6), neq(0:6)
    real(kind=8) :: Q(0:6)

     do k=1,nzLocal
        do j=1,ny
            do i=1,nx
                n(0) = g(0,i,j,k)+g(1,i,j,k)+g(2,i,j,k)+g(3,i,j,k)+g(4,i,j,k)+g(5,i,j,k)+g(6,i,j,k)
                n(1) = g(1,i,j,k)-g(2,i,j,k)
                n(2) = g(3,i,j,k)-g(4,i,j,k)
                n(3) = g(5,i,j,k)-g(6,i,j,k)
                n(4) = g(1,i,j,k)+g(2,i,j,k)+g(3,i,j,k)+g(4,i,j,k)+g(5,i,j,k)+g(6,i,j,k)
                n(5) = g(1,i,j,k)+g(2,i,j,k)-g(3,i,j,k)-g(4,i,j,k)
                n(6) = g(1,i,j,k)+g(2,i,j,k)-g(5,i,j,k)-g(6,i,j,k)

                Q(0) = 1.0d0
                Q(1) = sig_k
                Q(2) = sig_k
                Q(3) = sig_k
                Q(4) = 1.2d0
                Q(5) = 1.2d0
                Q(6) = 1.2d0

                neq(0) = temp(i,j,k)
                neq(1) = temp(i,j,k)*u(i,j,k)
                neq(2) = temp(i,j,k)*v(i,j,k)
                neq(3) = temp(i,j,k)*w(i,j,k)
                neq(4) = temp(i,j,k)*0.5d0
                neq(5) = 0.0d0
                neq(6) = 0.0d0

                do alpha=0,6
                    n_post(alpha)=n(alpha)-Q(alpha)*(n(alpha)-neq(alpha))
                enddo

                g_post(0,i,j,k) = n_post(0)-n_post(4)
                g_post(1,i,j,k) = n_post(1)/2.0d0+n_post(4)/6.0d0+n_post(5)/6.0d0+n_post(6)/6.0d0
                g_post(2,i,j,k) = -n_post(1)/2.0d0+n_post(4)/6.0d0+n_post(5)/6.0d0+n_post(6)/6.0d0
                g_post(3,i,j,k) = n_post(2)/2.0d0+n_post(4)/6.0d0-n_post(5)/3.0d0+n_post(6)/6.0d0
                g_post(4,i,j,k) = -n_post(2)/2.0d0+n_post(4)/6.0d0-n_post(5)/3.0d0+n_post(6)/6.0d0
                g_post(5,i,j,k) = n_post(3)/2.0d0+n_post(4)/6.0d0+n_post(5)/6.0d0-n_post(6)/3.0d0
                g_post(6,i,j,k) = -n_post(3)/2.0d0+n_post(4)/6.0d0+n_post(5)/6.0d0-n_post(6)/3.0d0

            enddo
        enddo
    end do

    do k = 1, nzLocal
       do j = 1, ny
            do i = 1, nx
                g(0,i,j,k)=g_post(0,i,j,k)
            enddo
        enddo
    end do

    return
end subroutine collision_T

subroutine update()
    use commondata
    implicit none

    call MPI_SENDRECV(f_post(0,1,1,1),19*nx*ny,MPI_REAL8,downid,0,f_post(0,1,1,nzLocal+1),19*nx*ny,MPI_REAL8,upid,0&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    call MPI_SENDRECV(f_post(0,1,1,nzLocal),19*nx*ny,MPI_REAL8,upid,1,f_post(0,1,1,0),19*nx*ny,MPI_REAL8,downid,1&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    call MPI_SENDRECV(g_post(0,1,1,1),7*nx*ny,MPI_REAL8,downid,2,g_post(0,1,1,nzLocal+1),7*nx*ny,MPI_REAL8,upid,2&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    call MPI_SENDRECV(g_post(0,1,1,nzLocal),7*nx*ny,MPI_REAL8,upid,3,g_post(0,1,1,0),7*nx*ny,MPI_REAL8,downid,3&
    ,MPI_COMM_WORLD,ISTATUS,IERR)

    return
end subroutine

subroutine interploate()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k, alpha
    real(kind=8) :: delta_x, delta_y, delta_z
    real(kind=8) :: interpolateF
    real(kind=8) :: f0,f1,f2
!--------------------------------------------------------------------------------------------------------------------
    do k=1, nzLocal
        do j=1, ny
            do i=1, nx
                do alpha=1, 18
                    delta_x=dble(ex(alpha))*dt
                    delta_y=dble(ey(alpha))*dt
                    delta_z=dble(ez(alpha))*dt
!------------------------------------------yoz------------------------------------------------------------------------
    if(alpha==5 .or. alpha==6 .or. alpha==15 .or. alpha==16 .or. alpha==17 .or. alpha==18)then
        f0 = interpolateF(zGrid(inter_z(k,1))+delta_z, zGrid(inter_z(k,2))+delta_z, zGrid(inter_z(k,3))+delta_z, &
                        zGrid(inter_z(k,2)), f_post(alpha,inter_x(i,2),inter_y(j,1),inter_z(k,1)), &
                        f_post(alpha,inter_x(i,2),inter_y(j,1),inter_z(k,2)), &
                        f_post(alpha,inter_x(i,2),inter_y(j,1),inter_z(k,3)))

        f1 = interpolateF(zGrid(inter_z(k,1))+delta_z, zGrid(inter_z(k,2))+delta_z, zGrid(inter_z(k,3))+delta_z, &
                        zGrid(inter_z(k,2)), f_post(alpha,inter_x(i,2),inter_y(j,2),inter_z(k,1)), &
                        f_post(alpha,inter_x(i,2),inter_y(j,2),inter_z(k,2)), &
                        f_post(alpha,inter_x(i,2),inter_y(j,2),inter_z(k,3)))

        f2 = interpolateF(zGrid(inter_z(k,1))+delta_z, zGrid(inter_z(k,2))+delta_z, zGrid(inter_z(k,3))+delta_z, &
                        zGrid(inter_z(k,2)), f_post(alpha,inter_x(i,2),inter_y(j,3),inter_z(k,1)), &
                        f_post(alpha,inter_x(i,2),inter_y(j,3),inter_z(k,2)), &
                        f_post(alpha,inter_x(i,2),inter_y(j,3),inter_z(k,3)))

        f(alpha,i,j,k) = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                            , yGrid(inter_y(j,2)), f0, f1, f2)

!-----------------------------------------xoz---------------------------------------------------------------------------
    else if(alpha==1 .or. alpha==2 .or. alpha==11 .or. alpha==12 .or. alpha==13 .or. alpha==14)then
        f0 = interpolateF(zGrid(inter_z(k,1))+delta_z, zGrid(inter_z(k,2))+delta_z, zGrid(inter_z(k,3))+delta_z, &
                        zGrid(inter_z(k,2)), f_post(alpha,inter_x(i,1),inter_y(j,2),inter_z(k,1)), &
                        f_post(alpha,inter_x(i,1),inter_y(j,2),inter_z(k,2)), &
                        f_post(alpha,inter_x(i,1),inter_y(j,2),inter_z(k,3)))

        f1 = interpolateF(zGrid(inter_z(k,1))+delta_z, zGrid(inter_z(k,2))+delta_z, zGrid(inter_z(k,3))+delta_z,&
                        zGrid(inter_z(k,2)), f_post(alpha,inter_x(i,2),inter_y(j,2),inter_z(k,1)), &
                        f_post(alpha,inter_x(i,2),inter_y(j,2),inter_z(k,2)), &
                        f_post(alpha,inter_x(i,2),inter_y(j,2),inter_z(k,3)))

        f2 = interpolateF(zGrid(inter_z(k,1))+delta_z, zGrid(inter_z(k,2))+delta_z, zGrid(inter_z(k,3))+delta_z, &
                        zGrid(inter_z(k,2)), f_post(alpha,inter_x(i,3),inter_y(j,2),inter_z(k,1)), &
                        f_post(alpha,inter_x(i,3),inter_y(j,2),inter_z(k,2)), &
                        f_post(alpha,inter_x(i,3),inter_y(j,2),inter_z(k,3)))

        f(alpha,i,j,k) = interpolateF(xGrid(inter_x(i,1))+delta_x, xGrid(inter_x(i,2))+delta_x, xGrid(inter_x(i,3))+delta_x&
                        , xGrid(inter_x(i,2)), f0, f1, f2)

!-----------------------------------------xoy---------------------------------------------------------------------------
    else if(alpha==3 .or. alpha==4 .or. alpha==7 .or. alpha==8 .or. alpha==9 .or. alpha==10)then
        f0 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y, &
                        yGrid(inter_y(j,2)), f_post(alpha,inter_x(i,1),inter_y(j,1),inter_z(k,2)), &
                        f_post(alpha,inter_x(i,1),inter_y(j,2),inter_z(k,2)), &
                        f_post(alpha,inter_x(i,1),inter_y(j,3),inter_z(k,2)))

        f1 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y, &
                        yGrid(inter_y(j,2)), f_post(alpha,inter_x(i,2),inter_y(j,1),inter_z(k,2)), &
                        f_post(alpha,inter_x(i,2),inter_y(j,2),inter_z(k,2)), &
                        f_post(alpha,inter_x(i,2),inter_y(j,3),inter_z(k,2)))

        f2 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y, &
                        yGrid(inter_y(j,2)), f_post(alpha,inter_x(i,3),inter_y(j,1),inter_z(k,2)), &
                        f_post(alpha,inter_x(i,3),inter_y(j,2),inter_z(k,2)), &
                        f_post(alpha,inter_x(i,3),inter_y(j,3),inter_z(k,2)))

        f(alpha,i,j,k) = interpolateF(xGrid(inter_x(i,1))+delta_x, xGrid(inter_x(i,2))+delta_x, xGrid(inter_x(i,3))+delta_x&
                        , xGrid(inter_x(i,2)), f0, f1, f2)
    end if
                end do
            end do
        end do
    end do

    do k=1, nzLocal
        do j=1, ny
            do i=1, nx
                do alpha=1, 6
                    delta_x=dble(ex(alpha))*dt
                    delta_y=dble(ey(alpha))*dt
                    delta_z=dble(ez(alpha))*dt
    if(alpha==3 .or. alpha==4)then
        g(alpha,i,j,k)=interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                                    , yGrid(inter_y(j,2)),g_post(alpha,inter_x(i,2),inter_y(j,1),inter_z(k,2))&
                                    , g_post(alpha,inter_x(i,2),inter_y(j,2),inter_z(k,2))&
                                    , g_post(alpha,inter_x(i,2),inter_y(j,3),inter_z(k,2)))

    elseif(alpha==5 .or. alpha==6)then
        g(alpha,i,j,k)=interpolateF(zGrid(inter_z(k,1))+delta_z, zGrid(inter_z(k,2))+delta_z, zGrid(inter_z(k,3))+delta_z&
                                    , zGrid(inter_z(k,2)),g_post(alpha,inter_x(i,2),inter_y(j,2),inter_z(k,1))&
                                    , g_post(alpha,inter_x(i,2),inter_y(j,2),inter_z(k,2))&
                                    , g_post(alpha,inter_x(i,2),inter_y(j,2),inter_z(k,3)))

    elseif(alpha==1 .or. alpha==2)then
        g(alpha,i,j,k)=interpolateF(xGrid(inter_x(i,1))+delta_x, xGrid(inter_x(i,2))+delta_x, xGrid(inter_x(i,3))+delta_x&
                                    , xGrid(inter_x(i,2)),g_post(alpha,inter_x(i,1),inter_y(j,2),inter_z(k,2))&
                                    , g_post(alpha,inter_x(i,2),inter_y(j,2),inter_z(k,2))&
                                    , g_post(alpha,inter_x(i,3),inter_y(j,2),inter_z(k,2)))
    end if
                end do
            end do
        end do
    end do

return
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
    integer :: i,j,k

    do k=1,nzLocal
        do i=1,nx
            !left side
            f(3,i,1,k)=f_post(4,i,1,k)
            f(7,i,1,k)=f_post(10,i,1,k)
            f(8,i,1,k)=f_post(9,i,1,k)
            f(15,i,1,k)=f_post(18,i,1,k)
            f(17,i,1,k)=f_post(16,i,1,k)

            !right side
            f(4,i,ny,k)=f_post(3,i,ny,k)
            f(9,i,ny,k)=f_post(8,i,ny,k)
            f(10,i,ny,k)=f_post(7,i,ny,k)
            f(16,i,ny,k)=f_post(17,i,ny,k)
            f(18,i,ny,k)=f_post(15,i,ny,k)
        end do
    end do

    do k=1,nzLocal
        do j=1,ny
            !front side
            f(2,nx,j,k)=f_post(1,nx,j,k)
            f(8,nx,j,k)=f_post(9,nx,j,k)
            f(10,nx,j,k)=f_post(7,nx,j,k)
            f(12,nx,j,k)=f_post(13,nx,j,k)
            f(14,nx,j,k)=f_post(11,nx,j,k)

            !back side
            f(1,1,j,k)=f_post(2,1,j,k)
            f(9,1,j,k)=f_post(8,1,j,k)
            f(7,1,j,k)=f_post(10,1,j,k)
            f(13,1,j,k)=f_post(12,1,j,k)
            f(11,1,j,k)=f_post(14,1,j,k)
        end do
    end do

    if(MYID==0)then
        do j=1,ny
            do i=1,nx
                f(5,i,j,1)=f_post(6,i,j,1)
                f(11,i,j,1)=f_post(14,i,j,1)
                f(12,i,j,1)=f_post(13,i,j,1)
                f(15,i,j,1)=f_post(18,i,j,1)
                f(16,i,j,1)=f_post(17,i,j,1)
            end do
        end do

    elseif(MYID==NPROC-1)then
        do j=1,ny
            do i=1,nx
                f(6,i,j,nzLocal)=f_post(5,i,j,nzLocal)
                f(13,i,j,nzLocal)=f_post(12,i,j,nzLocal)
                f(14,i,j,nzLocal)=f_post(11,i,j,nzLocal)
                f(17,i,j,nzLocal)=f_post(16,i,j,nzLocal)
                f(18,i,j,nzLocal)=f_post(15,i,j,nzLocal)
            end do
        end do
    end if

    return
end subroutine bounceback_u

subroutine bounceback_T()
    use commondata
    implicit none
    integer :: i,j,k

    do k=1,nzLocal
        do i=1,nx
            !left side
            g(3,i,1,k) = -g_post(4,i,1,k) + 2.0d0*omega_T(3)*Thot

            !right side
            g(4,i,ny,k) = -g_post(3,i,ny,k) + 2.0d0*omega_T(4)*Tcold
        end do
    end do

    do k=1,nzLocal
        do j=1,ny
            !front side
            g(2,nx,j,k) = g_post(1,nx,j,k)

            !back side
            g(1,1,j,k) = g_post(2,1,j,k)
        end do
    end do

    if(MYID==0)then
        do j=1,ny
            do i=1,nx
                g(5,i,j,1) = g_post(6,i,j,1)
            end do
        end do

    elseif(MYID==NPROC-1)then
        do j=1,ny
            do i=1,nx
                g(6,i,j,nzLocal) = g_post(5,i,j,nzLocal)
            end do
        end do
    end if

    return
end subroutine bounceback_T

subroutine macro_u()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k

    do k=1,nzLocal
        do j=1,ny
            do i=1,nx
                rho(i,j,k)=f(0,i,j,k)+f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)&
                        &+f(6,i,j,k)+f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)&
                        &+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k)

                u(i,j,k)=(f(1,i,j,k)-f(2,i,j,k)+f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k)&
                        &+f(11,i,j,k)-f(12,i,j,k)+f(13,i,j,k)-f(14,i,j,k))/rho(i,j,k)
                v(i,j,k)=(f(3,i,j,k)-f(4,i,j,k)+f(7,i,j,k)+f(8,i,j,k)-f(9,i,j,k)-f(10,i,j,k)&
                        &+f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k))/rho(i,j,k)
                w(i,j,k)=(f(5,i,j,k)-f(6,i,j,k)+f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k)&
                        &+f(15,i,j,k)+f(16,i,j,k)-f(17,i,j,k)-f(18,i,j,k)+0.5d0*dt*Fz(i,j,k))/rho(i,j,k)
            end do
        end do
    end do
    return
end subroutine macro_u

subroutine macro_t()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k

    do k=1,nzLocal
        do j=1,ny
            do i=1,nx
                temp(i,j,k) = g(0,i,j,k)+g(1,i,j,k)+g(2,i,j,k)+g(3,i,j,k)+g(4,i,j,k)+g(5,i,j,k)+g(6,i,j,k)
            end do
        end do
    end do
    return
end subroutine macro_t

subroutine check()
    use commondata
    implicit none
    integer(kind=4) :: i , j, k
    real(kind=8) :: error1, error2, error3, error4
    real(kind=8) :: error1_all, error2_all,error3_all, error4_all

    error1 = 0.0d0
    error2 = 0.0d0
    error1_all = 0.0d0
    error2_all = 0.0d0
    error3 = 0.0d0
    error4 = 0.0d0
    error3_all = 0.0d0
    error4_all = 0.0d0

    do k=1,nzLocal
        do j=1,ny
            do i=1,nx
                error1 = error1+(u(i,j,k)-up(i,j,k))**2+(v(i,j,k)-vp(i,j,k))**2+(w(i,j,k)-wp(i,j,k))**2
                error2 = error2+u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)

                up(i,j,k) = u(i,j,k)
                vp(i,j,k) = v(i,j,k)
                wp(i,j,k) = w(i,j,k)
            end do
        end do
    end do

    call MPI_ALLREDUCE(error1,error1_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
    call MPI_ALLREDUCE(error2,error2_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)

    errorU = sqrt(error1_all)/sqrt(error2_all)

    do k=1,nzLocal
        do j=1,ny
            do i=1,nx
                error3 = error3+(temp(i,j,k)-utemp(i,j,k))**2
                error4 = error4+temp(i,j,k)*temp(i,j,k)

                utemp(i,j,k) = temp(i,j,k)
            end do
        end do
    end do

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
    call MPI_GATHERV(w,count2d(MYID),MPI_REAL8,w_all,count2d,displ2d,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    call MPI_GATHERV(rho,count2d(MYID),MPI_REAL8,rho_all,count2d,displ2d,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    call MPI_GATHERV(temp,count2d(MYID),MPI_REAL8,T_all,count2d,displ2d,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    return
end subroutine

subroutine output_ASCII
    use commondata
    implicit none
    integer :: i, j, k
    character(len=100) :: filename

    write(filename,*) int(nx)
    filename = adjustl(filename)

    open(unit=02,file="side-heated convection"//trim(filename)//'.dat',status='unknown')
    write(02,*) 'TITLE="side-heated convection "'
    write(02,*) 'VARIABLES="X" "Y" "Z" "U" "V" "W" "RHO" "T" '
    write(02,101) nx, ny ,nz
    do k=1,nz
        do j=1,ny
            do i=1,nx
                write(02,100) xGrid(i), yGrid(j), z(k),u_all(i,j,k), v_all(i,j,k), w_all(i,j,k), rho_all(i,j,k), T_all(i,j,k)
            enddo
        enddo
    end do
100 format(1x,2(e11.4,' '),10(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'K=',1x,i5,1x,'F=POINT')
    close(02)
    return
end subroutine output_ASCII
