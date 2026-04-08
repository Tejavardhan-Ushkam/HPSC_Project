!======================================================================
! dem_physics.f90
! All heavy-compute Fortran 90 routines for the DEM solver.
! Every subroutine is OpenMP-parallelised where safe.
! Called from C++ via ISO C binding (bind(C) + trailing underscore).
!
! Modules provided:
!   1. set_omp_threads / get_omp_threads
!   2. init_particles
!   3. zero_forces
!   4. add_gravity
!   5. compute_particle_contacts
!   6. compute_wall_contacts
!   7. integrate_particles
!   8. compute_kinetic_energy
!======================================================================
module dem_physics
    use omp_lib
    implicit none
    integer, save :: g_nthreads = 1
contains

!----------------------------------------------------------------------
! 1. set_omp_threads / get_omp_threads
!----------------------------------------------------------------------
subroutine set_omp_threads(nthreads) bind(C,name='set_omp_threads_')
    use iso_c_binding, only: c_int
    integer(c_int), intent(in) :: nthreads
    g_nthreads = nthreads
    call omp_set_num_threads(g_nthreads)
end subroutine

function get_omp_threads() result(nt) bind(C,name='get_omp_threads_')
    use iso_c_binding, only: c_int
    integer(c_int) :: nt
    nt = omp_get_max_threads()
end function

!======================================================================
! 2. init_particles
! mode 1 : random cloud  (multi-particle default)
! mode 2 : single-particle free-fall  (z0=0.9*Lz, v=0)
! mode 3 : single-particle bounce     (z0=0.8*Lz, v=0)
! mode 4 : constant-velocity test     (z0=0.1*Lz, vz=2.0, g=0)
!======================================================================
subroutine init_particles(x,y,z,vx,vy,vz,mass,radius, &
                          n,mode,lx,ly,lz,seed) &
    bind(C,name='init_particles_')
    use iso_c_binding, only: c_int,c_double
    integer(c_int),  intent(in)    :: n,mode
    real(c_double),  intent(in)    :: lx,ly,lz
    integer(c_int),  intent(inout) :: seed
    real(c_double),  intent(out)   :: x(n),y(n),z(n)
    real(c_double),  intent(out)   :: vx(n),vy(n),vz(n)
    real(c_double),  intent(out)   :: mass(n),radius(n)

    integer :: i, iseed
    real(c_double) :: rv, r0, m0

    r0    = 0.02d0 * lx    ! radius = 2% of domain
    m0    = 1.0d0
    iseed = seed

    select case(mode)

    case(2)   !── free-fall ──────────────────────────────────────
        x(1)=0.5d0*lx; y(1)=0.5d0*ly; z(1)=0.9d0*lz
        vx(1)=0d0; vy(1)=0d0; vz(1)=0d0
        mass(1)=m0; radius(1)=r0

    case(3)   !── bounce ────────────────────────────────────────
        x(1)=0.5d0*lx; y(1)=0.5d0*ly; z(1)=0.8d0*lz
        vx(1)=0d0; vy(1)=0d0; vz(1)=0d0
        mass(1)=m0; radius(1)=r0

    case(4)   !── constant velocity (caller must set gz=0) ──────
        x(1)=0.5d0*lx; y(1)=0.5d0*ly; z(1)=0.1d0*lz
        vx(1)=0.8d0; vy(1)=-0.35d0; vz(1)=0.5d0   ! non-zero 3D velocity
        mass(1)=m0; radius(1)=r0

    case default  !── random cloud ───────────────────────────────
        do i=1,n
            iseed=mod(iseed*1664525+1013904223,2147483647)
            rv=real(mod(abs(iseed),10000),8)/10000d0
            x(i)=r0+rv*(lx-2d0*r0)

            iseed=mod(iseed*1664525+1013904223,2147483647)
            rv=real(mod(abs(iseed),10000),8)/10000d0
            y(i)=r0+rv*(ly-2d0*r0)

            iseed=mod(iseed*1664525+1013904223,2147483647)
            rv=real(mod(abs(iseed),10000),8)/10000d0
            z(i)=r0+rv*(lz-2d0*r0)

            vx(i)=0d0; vy(i)=0d0; vz(i)=0d0
            mass(i)=m0; radius(i)=r0
        end do
    end select
end subroutine

!======================================================================
! 3. zero_forces  – OpenMP static parallel
!======================================================================
subroutine zero_forces(fx,fy,fz,n) bind(C,name='zero_forces_')
    use iso_c_binding, only: c_int,c_double
    integer(c_int), intent(in)  :: n
    real(c_double), intent(out) :: fx(n),fy(n),fz(n)
    integer :: i
    !$omp parallel do schedule(static)
    do i=1,n
        fx(i)=0d0; fy(i)=0d0; fz(i)=0d0
    end do
    !$omp end parallel do
end subroutine

!======================================================================
! 4. add_gravity  – OpenMP static parallel
!======================================================================
subroutine add_gravity(fx,fy,fz,mass,n,gx,gy,gz) &
    bind(C,name='add_gravity_')
    use iso_c_binding, only: c_int,c_double
    integer(c_int), intent(in)    :: n
    real(c_double), intent(inout) :: fx(n),fy(n),fz(n)
    real(c_double), intent(in)    :: mass(n),gx,gy,gz
    integer :: i
    !$omp parallel do schedule(static)
    do i=1,n
        fx(i)=fx(i)+mass(i)*gx
        fy(i)=fy(i)+mass(i)*gy
        fz(i)=fz(i)+mass(i)*gz
    end do
    !$omp end parallel do
end subroutine

!======================================================================
! 5. compute_particle_contacts
! O(N^2) all-pairs; OpenMP dynamic over outer loop;
! atomic updates for Newton-3 force accumulation.
!======================================================================
subroutine compute_particle_contacts(x,y,z,vx,vy,vz, &
                                     fx,fy,fz,radius,  &
                                     n,kn,gamma_n,ncontacts) &
    bind(C,name='compute_particle_contacts_')
    use iso_c_binding, only: c_int,c_double
    integer(c_int), intent(in)    :: n
    real(c_double), intent(in)    :: x(n),y(n),z(n)
    real(c_double), intent(in)    :: vx(n),vy(n),vz(n)
    real(c_double), intent(inout) :: fx(n),fy(n),fz(n)
    real(c_double), intent(in)    :: radius(n),kn,gamma_n
    integer(c_int), intent(out)   :: ncontacts

    integer :: i,j
    real(c_double) :: rx,ry,rz,d,delta,nx,ny,nz,dvx,dvy,dvz,vn,fn,fcx,fcy,fcz

    ncontacts = 0

    !$omp parallel do schedule(dynamic,16) &
    !$omp   private(j,rx,ry,rz,d,delta,nx,ny,nz, &
    !$omp           dvx,dvy,dvz,vn,fn,fcx,fcy,fcz) &
    !$omp   reduction(+:ncontacts)
    do i=1,n-1
        do j=i+1,n
            rx=x(j)-x(i); ry=y(j)-y(i); rz=z(j)-z(i)
            d=sqrt(rx*rx+ry*ry+rz*rz)
            if(d<1d-14) cycle
            delta=radius(i)+radius(j)-d
            if(delta<=0d0) cycle

            nx=rx/d; ny=ry/d; nz=rz/d
            dvx=vx(j)-vx(i); dvy=vy(j)-vy(i); dvz=vz(j)-vz(i)
            vn=dvx*nx+dvy*ny+dvz*nz
            fn=max(0d0, kn*delta-gamma_n*vn)

            fcx=fn*nx; fcy=fn*ny; fcz=fn*nz

            !$omp atomic
            fx(i)=fx(i)-fcx
            !$omp atomic
            fy(i)=fy(i)-fcy
            !$omp atomic
            fz(i)=fz(i)-fcz

            !$omp atomic
            fx(j)=fx(j)+fcx
            !$omp atomic
            fy(j)=fy(j)+fcy
            !$omp atomic
            fz(j)=fz(j)+fcz

            ncontacts=ncontacts+1
        end do
    end do
    !$omp end parallel do
end subroutine

!======================================================================
! 6. compute_wall_contacts
! Six walls; also returns max_wall_overlap for verification test.
! OpenMP static parallel.
!======================================================================
subroutine compute_wall_contacts(x,y,z,vx,vy,vz, &
                                  fx,fy,fz,radius,  &
                                  n,kn,gamma_n,      &
                                  lx,ly,lz,max_wall_overlap) &
    bind(C,name='compute_wall_contacts_')
    use iso_c_binding, only: c_int,c_double
    integer(c_int), intent(in)    :: n
    real(c_double), intent(in)    :: x(n),y(n),z(n)
    real(c_double), intent(in)    :: vx(n),vy(n),vz(n)
    real(c_double), intent(inout) :: fx(n),fy(n),fz(n)
    real(c_double), intent(in)    :: radius(n),kn,gamma_n,lx,ly,lz
    real(c_double), intent(out)   :: max_wall_overlap

    integer :: i
    real(c_double) :: delta,fn,local_max

    max_wall_overlap = 0d0

    !$omp parallel do schedule(static) &
    !$omp   private(delta,fn,local_max) &
    !$omp   reduction(max:max_wall_overlap)
    do i=1,n
        local_max = 0d0

        ! z=0 floor
        delta=radius(i)-z(i)
        if(delta>0d0)then
            fn=max(0d0, kn*delta+gamma_n*vz(i))   ! vz towards floor is negative, sign correct
            fz(i)=fz(i)+fn
            if(delta>local_max) local_max=delta
        end if

        ! z=Lz ceiling
        delta=radius(i)-(lz-z(i))
        if(delta>0d0)then
            fn=max(0d0, kn*delta-gamma_n*vz(i))
            fz(i)=fz(i)-fn
            if(delta>local_max) local_max=delta
        end if

        ! x=0
        delta=radius(i)-x(i)
        if(delta>0d0)then
            fn=max(0d0, kn*delta+gamma_n*vx(i))
            fx(i)=fx(i)+fn
            if(delta>local_max) local_max=delta
        end if

        ! x=Lx
        delta=radius(i)-(lx-x(i))
        if(delta>0d0)then
            fn=max(0d0, kn*delta-gamma_n*vx(i))
            fx(i)=fx(i)-fn
            if(delta>local_max) local_max=delta
        end if

        ! y=0
        delta=radius(i)-y(i)
        if(delta>0d0)then
            fn=max(0d0, kn*delta+gamma_n*vy(i))
            fy(i)=fy(i)+fn
            if(delta>local_max) local_max=delta
        end if

        ! y=Ly
        delta=radius(i)-(ly-y(i))
        if(delta>0d0)then
            fn=max(0d0, kn*delta-gamma_n*vy(i))
            fy(i)=fy(i)-fn
            if(delta>local_max) local_max=delta
        end if

        if(local_max>max_wall_overlap) max_wall_overlap=local_max
    end do
    !$omp end parallel do
end subroutine

!======================================================================
! 7. integrate_particles  – semi-implicit Euler, OpenMP static
!======================================================================
subroutine integrate_particles(x,y,z,vx,vy,vz, &
                                fx,fy,fz,mass,n,dt) &
    bind(C,name='integrate_particles_')
    use iso_c_binding, only: c_int,c_double
    integer(c_int), intent(in)    :: n
    real(c_double), intent(inout) :: x(n),y(n),z(n)
    real(c_double), intent(inout) :: vx(n),vy(n),vz(n)
    real(c_double), intent(in)    :: fx(n),fy(n),fz(n),mass(n),dt
    integer :: i
    !$omp parallel do schedule(static)
    do i=1,n
        vx(i)=vx(i)+fx(i)/mass(i)*dt
        vy(i)=vy(i)+fy(i)/mass(i)*dt
        vz(i)=vz(i)+fz(i)/mass(i)*dt
        x(i)=x(i)+vx(i)*dt
        y(i)=y(i)+vy(i)*dt
        z(i)=z(i)+vz(i)*dt
    end do
    !$omp end parallel do
end subroutine

!======================================================================
! 8. compute_kinetic_energy  – OpenMP reduction
!======================================================================
subroutine compute_kinetic_energy(vx,vy,vz,mass,n,ke) &
    bind(C,name='compute_kinetic_energy_')
    use iso_c_binding, only: c_int,c_double
    integer(c_int), intent(in)  :: n
    real(c_double), intent(in)  :: vx(n),vy(n),vz(n),mass(n)
    real(c_double), intent(out) :: ke
    integer :: i
    ke=0d0
    !$omp parallel do reduction(+:ke) schedule(static)
    do i=1,n
        ke=ke+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
    end do
    !$omp end parallel do
end subroutine

end module dem_physics
