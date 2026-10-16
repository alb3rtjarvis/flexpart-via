!*****************************************************************************
!                                                                            *
! L. Bakels 2022: This module contains all interpolation subroutines         *
!                 Code has been organised into subroutines                   *
!                 Vertical logarithmic interpolation is optional (par_mod)   *
!                                                                            *
!*****************************************************************************

module interpol_mod
   use par_mod
   use com_mod
   use windfields_mod
   use particle_mod

   implicit none

   real, allocatable, dimension(:, :) :: &
      uprof, vprof, wprof, &
      usigprof, vsigprof, wsigprof, &
      rhoprof, rhogradprof
   logical, allocatable, dimension(:, :) :: &
      indzindicator, depoindicator

   real :: u, v, w, usig, vsig, wsig

   real :: p1, p2, p3, p4, ddx, ddy, rddx, rddy, dtt, dt1, dt2
   real :: xtn, ytn
   real :: dz1out, dz2out
   integer :: nix, njy
   integer :: ix, jy, ixp, jyp, ngrid, indz, indzp
   integer :: induv, indpuv
   logical :: lbounds(2) ! marking particles below or above bounds
#ifdef ETA
   real, allocatable, dimension(:, :) :: wprofeta ! ,wsigprofeta
   real :: ueta, veta, weta, wsigeta
   integer :: indzeta, indzpeta
   logical :: lbounds_w(2), lbounds_uv(2) ! marking particles below or above bounds
#endif

#ifdef ETA
   private :: interpol_wind_eta, stdev_eta, interpol_partoutput_val_eta
#else
   private :: interpol_wind_meter, stdev_meter, interpol_partoutput_val_meter
#endif

#ifdef VIA

   abstract interface
      subroutine interp_via_signature(itime, jpart, yint_sr, yint_rh, ytint, icmv)
         implicit none
         integer, intent(in) :: icmv, itime, jpart
         real, intent(out) :: yint_rh, yint_sr, ytint
      end subroutine interp_via_signature
   end interface

   procedure(interp_via_signature), pointer :: interp_via_calc => null()


   abstract interface
      subroutine decay_signature(jpart, sr_interp, rh_interp, decay_rate)
         implicit none
         integer, intent(in) :: jpart
         real, intent(in) :: sr_interp, rh_interp
         real, intent(out) :: decay_rate
      end subroutine decay_signature
   end interface

   procedure(decay_signature), pointer :: decay_calc => null()
#endif

   interface hor_interpol
      procedure hor_interpol_4d, hor_interpol_2d
   end interface hor_interpol

   interface hor_interpol_nest
      procedure hor_interpol_4d_nest, hor_interpol_2d_nest
   end interface hor_interpol_nest

   interface find_ngrid
      procedure find_ngrid_dp, find_ngrid_sp
   end interface find_ngrid

! uprof,vprof,wprof,usigprof,vsigprof,wsigprof,indzindicator, wsigprofeta,
! rhoprof,rhogradprof,wprofeta,depoindicator,
#ifdef ETA
!$OMP THREADPRIVATE( &
!$OMP u,v,w,usig,vsig,wsig, &
!$OMP p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2,ix,jy,ixp,jyp, &
!$OMP ngrid,indz,indzp, &
!$OMP induv,indpuv,lbounds,lbounds_w,lbounds_uv, &
!$OMP indzeta,indzpeta,ueta,veta,weta,wsigeta, &
!$OMP xtn,ytn,nix,njy,dz1out,dz2out)
#else
!$OMP THREADPRIVATE( &
!$OMP u,v,w,usig,vsig,wsig, &
!$OMP p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2,ix,jy,ixp,jyp, &
!$OMP ngrid,indz,indzp, &
!$OMP induv,indpuv,lbounds,xtn,ytn,nix,njy,dz1out,dz2out)
#endif

contains

   subroutine alloc_interpol ! wsigprofeta(nzmax,numthreads),
      implicit none
      integer :: stat

      allocate (uprof(nzmax, numthreads), vprof(nzmax, numthreads), wprof(nzmax, numthreads), &
         usigprof(nzmax, numthreads), vsigprof(nzmax, numthreads), wsigprof(nzmax, numthreads), &
         rhoprof(nzmax, numthreads), rhogradprof(nzmax, numthreads), &
         indzindicator(nzmax, numthreads), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate interpol prof arrays"
#ifdef ETA
      allocate (wprofeta(nzmax, numthreads), stat=stat)
      if (stat .ne. 0) error stop "Could not allocate wprofeta"
#endif
      if (DRYDEP) then
         allocate (depoindicator(maxspec, numthreads), stat=stat)
         if (stat .ne. 0) error stop "Could not allocate depoindicator"
      end if
   end subroutine alloc_interpol

   subroutine dealloc_interpol ! wsigprofeta,
      deallocate (uprof, vprof, wprof, &
         usigprof, vsigprof, wsigprof, &
         rhoprof, rhogradprof, indzindicator)
#ifdef ETA
      deallocate (wprofeta)
#endif
      if (DRYDEP) deallocate (depoindicator)
   end subroutine dealloc_interpol

   subroutine init_interpol(itime, xt, yt, zt, zteta)

      ! This routine initialises all important values used in the interpol module
      ! This includes:
      ! - The current grid number in which the particle is positioned
      ! - The interpolation fractions of the grid (x,y,z) and of time

      integer, intent(in) :: itime             ! time step
      real, intent(in)    :: xt, yt             ! particle positions
      real, intent(in)    :: zt                ! height in meters
      real, intent(in)    :: zteta             ! height in eta coordinates

      call find_ngrid(xt, yt)
      call find_grid_indices(xt, yt)
      call find_grid_distances(xt, yt)
      call find_time_vars(itime)
      call find_z_level(zt, zteta)

   end subroutine init_interpol

   subroutine find_grid_indices(xt, yt)

      real, intent(in) :: xt, yt                 ! particle positions

      if (ngrid .gt. 0) then ! Nest
         xtn = (xt - xln(ngrid))*xresoln(ngrid)
         ytn = (yt - yln(ngrid))*yresoln(ngrid)
         ! ix=int(xtn)
         ! jy=int(ytn)
         ! nix=nint(xtn)
         ! njy=nint(ytn)
         nix = max(min(int(xtn), nxn(ngrid) - 1), 0)
         njy = max(min(int(ytn), nyn(ngrid) - 1), 0)
         ix = nix
         jy = njy
         ixp = ix + 1
         jyp = jy + 1
         return
      else
         ix = int(xt)
         jy = int(yt)
         nix = ix!nint(xt)
         njy = jy!nint(yt)
         ixp = ix + 1
         jyp = jy + 1
      end if

      ! eso: Temporary fix for particle exactly at north pole
      if (jyp .ge. nymax) then
         write (*, *) 'WARNING: interpol_mod.f90 jyp >= nymax. xt,yt:', xt, yt
         jyp = jyp - 1
      end if

      if (ixp .ge. nxmax) then
         write (*, *) 'WARNING: interpol_mod.f90 ixp >= nxmax. xt,yt:', xt, yt
         ixp = ixp - nxmax
      end if

   end subroutine find_grid_indices

   subroutine find_grid_distances(xt, yt)

      implicit none

      real, intent(in) :: xt, yt                 ! particle positions

      if (ngrid .le. 0) then
         ddx = xt - real(ix)
         ddy = yt - real(jy)
      else ! Nest
         ddx = xtn - real(ix)
         ddy = ytn - real(jy)
      end if
      rddx = 1.-ddx
      rddy = 1.-ddy
      p1 = rddx*rddy
      p2 = ddx*rddy
      p3 = rddx*ddy
      p4 = ddx*ddy

   end subroutine find_grid_distances

   subroutine find_time_vars(itime)

      integer, intent(in) :: itime             ! time step

      dt1 = real(itime - memtime(1))
      dt2 = real(memtime(2) - itime)
      dtt = 1./(dt1 + dt2)

   end subroutine find_time_vars

   subroutine find_z_level(zt, zteta)

      real, intent(in)     :: &
         zt, & ! height in meters
         zteta                   ! height in eta

#ifdef ETA
      call find_z_level_meters(zt)
      call find_z_level_eta(zteta)
#else
      call find_z_level_meters(zt)
#endif

   end subroutine find_z_level

   subroutine find_z_level_meters(zt)

      real, intent(in)     :: zt       ! height in meters
      integer              :: i

      indz = nz - 1
      indzp = nz
      if (zt .le. height(1)) then
         lbounds(1) = .true.
         lbounds(2) = .false.
         indz = 1
         indzp = 2
      else if (zt .ge. height(nz)) then
         lbounds(1) = .false.
         lbounds(2) = .true.
      else
         lbounds(1) = .false.
         lbounds(2) = .false.
         do i = 2, nz
            if (height(i) .gt. zt) then
               indz = i - 1
               indzp = i
               exit
            end if
         end do
      end if

   end subroutine find_z_level_meters

   subroutine find_z_level_meters_safe(zt, kz)

      real, intent(in)     :: zt       ! height in meters
      integer, intent(out) :: kz
      integer              :: i

      ! indz = nz - 1
      ! indzp = nz
      kz = nz - 1
      if (zt .le. height(1)) then
         ! lbounds(1) = .true.
         ! lbounds(2) = .false.
         ! indz = 1
         ! indzp = 2
         kz = 1
      else if (zt .ge. height(nz)) then
         ! lbounds(1) = .false.
         ! lbounds(2) = .true.
      else
         ! lbounds(1) = .false.
         ! lbounds(2) = .false.
         do i = 2, nz
            if (height(i) .gt. zt) then
               ! indz = i - 1
               ! indzp = i
               kz = i - 1
               exit
            end if
         end do
      end if

   end subroutine find_z_level_meters_safe

#ifdef ETA
   subroutine find_z_level_eta(zteta)

      real, intent(in)       :: zteta    ! height in eta coordinates

      call find_z_level_eta_w(zteta)

      call find_z_level_eta_uv(zteta)

   end subroutine find_z_level_eta

   subroutine find_z_level_eta_w(zteta)

      real, intent(in)       :: zteta    ! height in eta coordinates
      integer                :: i        ! loop variable

      indzeta = nz - 1
      indzpeta = nz
      ! Flag particles that are above or below bounds
      if (zteta .ge. wheight(1)) then
         lbounds_w(1) = .true.
         lbounds_w(2) = .false.
         indzeta = 1
         indzpeta = 2
      else if (zteta .le. wheight(nz)) then
         lbounds_w(1) = .false.
         lbounds_w(2) = .true.
      else
         lbounds_w(1) = .false.
         lbounds_w(2) = .false.
         do i = 2, nz
            if (wheight(i) .lt. zteta) then
               indzeta = i - 1
               indzpeta = i
               exit
            end if
         end do
      end if

   end subroutine find_z_level_eta_w

   subroutine find_z_level_eta_uv(zteta)

      real, intent(in)       :: zteta    ! height in eta coordinates
      integer                :: i        ! loop variable

      induv = nz - 1
      indpuv = nz
      if (zteta .gt. uvheight(1)) then
         lbounds_uv(1) = .true.
         lbounds_uv(2) = .false.
         induv = 1
         indpuv = 2
      else if (zteta .lt. uvheight(nz)) then
         lbounds_uv(1) = .false.
         lbounds_uv(2) = .true.
      else
         lbounds_uv(1) = .false.
         lbounds_uv(2) = .false.
         do i = 2, nz
            if (uvheight(i) .lt. zteta) then
               induv = i - 1
               indpuv = i
               exit
            end if
         end do
      end if

   end subroutine find_z_level_eta_uv

   subroutine find_z_level_eta_uv_safe(zteta, kz)

      real, intent(in)       :: zteta    ! height in eta coordinates
      integer, intent(out)   :: kz
      integer                :: i        ! loop variable

      kz = nz - 1
      ! indpuv = nz
      if (zteta .gt. uvheight(1)) then
         ! lbounds_uv(1) = .true.
         ! lbounds_uv(2) = .false.
         ! induv = 1
         ! indpuv = 2
         kz = 1
      else if (zteta .lt. uvheight(nz)) then
         ! lbounds_uv(1) = .false.
         ! lbounds_uv(2) = .true.
      else
         ! lbounds_uv(1) = .false.
         ! lbounds_uv(2) = .false.
         do i = 2, nz
            if (uvheight(i) .lt. zteta) then
               ! induv = i - 1
               ! indpuv = i
               kz = i - 1
               exit
            end if
         end do
      end if

   end subroutine find_z_level_eta_uv_safe
#endif

   subroutine find_vert_vars(vertlevels, zpos, zlevel, dz1, dz2, bounds, wlevel)

      !*****************************************************************************
      !                                                                            *
      ! This subroutine computes the vertical interpolation variables              *
      ! logarithmically, unless log_interpol=.false. in the par_mod                *
      !                                                                            *
      ! Author: L. Bakels                                                          *
      !*****************************************************************************

      real, intent(in)    :: vertlevels(:)    ! vertical levels in coordinate system
      real, intent(in)    :: zpos             ! verticle particle position
      integer, intent(in) :: zlevel           ! vertical level of interest
      logical, intent(in) :: bounds(2), wlevel ! flag marking if particles are
      ! outside bounds
      real, intent(inout) :: dz1, dz2          ! fractional distance to point 1
      ! (closer to ground) and 2
      real                :: dz, dh1, dh, pfact
      real                :: psint1(2), psint, pr1, pr2, pr_test
      ! pressure of encompassing levels
      integer             :: m

      ! Only do logarithmic interpolation when using ETA coordinates, since the
      ! levels are following pressure, while METER levels are linear.
      !##############################################################
      if (.not. log_interpol) then
         call find_vert_vars_lin(vertlevels, zpos, zlevel, dz1, dz2, bounds)
         return
      end if

      ! To check if taking the logarithm is safe
      if (wlevel) then
         pr_test = akm(zlevel + 1) + bkm(zlevel + 1)
      else
         pr_test = akz(zlevel + 1) + bkz(zlevel + 1)
      end if

      ! If the particle is below bounds (bounds(1)==.true.):
      if (bounds(1)) then
         dz1 = 0.
         dz2 = 1.
         ! If above bounds (bounds(2)==.true.):
      else if (bounds(2)) then
         dz1 = 1.
         dz2 = 0.

         ! Instead of the linear z variables, we need the ones that correspond to
         ! the pressure of the height of the particle in relation to the model levels
         !***************************************************************************
      else if (pr_test .eq. 0) then
         dz = 1./(vertlevels(zlevel + 1) - vertlevels(zlevel))
         dz1 = (zpos - vertlevels(zlevel))*dz
         dz2 = (vertlevels(zlevel + 1) - zpos)*dz
      else
         if (ngrid .le. 0) then
            do m = 1, 2
               call hor_interpol(ps, psint1(m), 1, memind(m), 1)
            end do
         else
            do m = 1, 2
               call hor_interpol_nest(psn, psint1(m), 1, memind(m), 1)
            end do
         end if
         call temporal_interpolation(psint1(1), psint1(2), psint)
         dh = vertlevels(zlevel + 1) - vertlevels(zlevel)
         dh1 = zpos - vertlevels(zlevel)
         if (wlevel) then
            pr1 = akm(zlevel) + bkm(zlevel)*psint
            pr2 = akm(zlevel + 1) + bkm(zlevel + 1)*psint
         else
            pr1 = akz(zlevel) + bkz(zlevel)*psint
            pr2 = akz(zlevel + 1) + bkz(zlevel + 1)*psint
         end if
         pfact = log(pr2/pr1)*dh1/dh
         dz = 1./(pr2 - pr1)
         dz1 = pr1*(exp(pfact) - 1.)*dz
         dz2 = 1.-dz1
      end if
      ! else if ((vertlevels(zlevel).eq.0).or.(vertlevels(zlevel+1).eq.0)) then
      !   ! Linear interpolation for bottom or top layer is zero
      !   dz=1./(vertlevels(zlevel+1)-vertlevels(zlevel))
      !   dz1=(zpos-vertlevels(zlevel))*dz
      !   dz2=(vertlevels(zlevel+1)-zpos)*dz
      ! else
      !   ! Logaritmic interpolation
      !   dz=1./(log(vertlevels(zlevel+1))-log(vertlevels(zlevel)))
      !   dz1=(log(zpos)-log(vertlevels(zlevel)))*dz
      !   dz2=(log(vertlevels(zlevel+1))-log(zpos))*dz
      ! endif
   end subroutine find_vert_vars

   subroutine find_vert_vars_lin(vertlevels, zpos, zlevel, dz1, dz2, bounds)

      real, intent(in)    :: vertlevels(:)   ! vertical levels in coordinate system
      real, intent(in)    :: zpos            ! verticle particle position
      integer, intent(in) :: zlevel          ! vertical level of interest
      logical, intent(in) :: bounds(2)       ! flag marking if particles are outside
      ! bounds
      real, intent(inout) :: dz1, dz2         ! fractional distance to point 1
      ! (closer to ground) and 2
      real                :: dz

      ! If the particle is below bounds (bounds(1)==.true.):
      if (bounds(1)) then
         dz1 = 0.
         dz2 = 1.
         ! If above bounds (bounds(2)==.true.):
      else if (bounds(2)) then
         dz1 = 1.
         dz2 = 0.
      else
         dz = 1./(vertlevels(zlevel + 1) - vertlevels(zlevel))
         dz1 = (zpos - vertlevels(zlevel))*dz
         dz2 = (vertlevels(zlevel + 1) - zpos)*dz
      end if
   end subroutine find_vert_vars_lin

   subroutine find_ngrid_dp(xt, yt)

      real eps
      real(kind=dp), intent(in) :: xt, yt ! particle positions on grid
      integer :: j

      eps = nxmax/3.e5
      if (nglobal .and. (real(yt) .gt. switchnorthg)) then
         ngrid = -1
      else if (sglobal .and. (real(yt) .lt. switchsouthg)) then
         ngrid = -2
      else
         ngrid = 0
         ! Temporary fix for nested layer edges: replaced eps with dxn and dyn (LB)
         do j = numbnests, 1, -1
            if (real(xt) .gt. xln(j) + dxn(j) .and. real(xt) .lt. xrn(j) - dxn(j) .and. &
               real(yt) .gt. yln(j) + dyn(j) .and. real(yt) .lt. yrn(j) - dyn(j)) then
               ngrid = j
               exit
            end if
         end do
      end if

   end subroutine find_ngrid_dp

   subroutine find_ngrid_sp(xt, yt)

      real :: eps
      real, intent(in) :: xt, yt ! particle positions on grid
      integer :: j

      eps = nxmax/3.e5
      if (nglobal .and. yt .gt. switchnorthg) then
         ngrid = -1
      else if (sglobal .and. yt .lt. switchsouthg) then
         ngrid = -2
      else
         ngrid = 0
         ! Temporary fix for nested layer edges: replaced eps with dxn and dyn (LB)
         do j = numbnests, 1, -1
            if (xt .gt. xln(j) + dxn(j) .and. xt .lt. xrn(j) - dxn(j) .and. &
               yt .gt. yln(j) + dyn(j) .and. yt .lt. yrn(j) - dyn(j)) then
               ngrid = j
               exit
            end if
         end do
      end if
   end subroutine find_ngrid_sp

   subroutine hor_interpol_4d(field, output, zlevel, indexh, ztot)

      integer, intent(in) :: zlevel, ztot, indexh   ! interpolation z level, z
      real, intent(in)    :: field(0:nxmax - 1, 0:nymax - 1, ztot, numwfmem)
      ! input field to interpolate
      real, intent(inout) :: output  ! interpolated values

      output = p1*field(ix, jy, zlevel, indexh) &
         + p2*field(ixp, jy, zlevel, indexh) &
         + p3*field(ix, jyp, zlevel, indexh) &
         + p4*field(ixp, jyp, zlevel, indexh)
   end subroutine hor_interpol_4d

   subroutine hor_interpol_2d(field, output)
      implicit none

      real, intent(in)    :: field(0:nxmax - 1, 0:nymax - 1)       ! 2D imput field
      real, intent(inout) :: output                           ! Interpolated value

      output = p1*field(ix, jy) &
         + p2*field(ixp, jy) &
         + p3*field(ix, jyp) &
         + p4*field(ixp, jyp)
   end subroutine hor_interpol_2d

   subroutine hor_interpol_4d_nest(field, output, zlevel, indexh, ztot)

      integer, intent(in) :: zlevel, ztot, indexh   ! interpolation z level, z
      real, intent(in)    :: field(0:nxmaxn - 1, 0:nymaxn - 1, ztot, numwfmem, numbnests)
      ! input field to interpolate
      real, intent(inout) :: output               ! interpolated values

      output = p1*field(ix, jy, zlevel, indexh, ngrid) &
         + p2*field(ixp, jy, zlevel, indexh, ngrid) &
         + p3*field(ix, jyp, zlevel, indexh, ngrid) &
         + p4*field(ixp, jyp, zlevel, indexh, ngrid)
   end subroutine hor_interpol_4d_nest

   subroutine hor_interpol_2d_nest(field, output)

      real, intent(in)    :: field(0:nxmaxn - 1, 0:nymaxn - 1, numbnests)
      ! input field to interpolate
      real, intent(inout) :: output           ! interpolated values

      output = p1*field(ix, jy, ngrid) &
         + p2*field(ixp, jy, ngrid) &
         + p3*field(ix, jyp, ngrid) &
         + p4*field(ixp, jyp, ngrid)
   end subroutine hor_interpol_2d_nest

   subroutine temporal_interpolation(time1, time2, output)

      real, intent(in)    :: time1, time2     ! input data at two timesteps
      real, intent(inout) :: output          ! interpolated data

      output = (time1*dt2 + time2*dt1)*dtt
   end subroutine temporal_interpolation

   subroutine vert_interpol(input1, input2, dz1, dz2, output)

      real, intent(in)    :: input1, input2   ! input data at two vertical levels,
      ! 1 being closer to ground
      real, intent(in)    :: dz1, dz2         ! logarithmic interpolation values
      real, intent(inout) :: output          ! interpolated data

      output = input1*dz2 + input2*dz1 ! input1**dz2 * input2**dz1
   end subroutine vert_interpol

   subroutine bilin_spatial_interpol(field, output, zlevel, dz1, dz2, ztot)

      integer, intent(in) :: zlevel, ztot        ! interpolation z level
      real, intent(in)    :: field(0:nxmax - 1, 0:nymax - 1, ztot, numwfmem)
      ! input field to interpolate
      real, intent(in)    :: dz1, dz2
      real, intent(inout) :: output(2)          ! interpolated values
      integer             :: m, n, indzh
      real                :: output1(2)

      do m = 1, 2

         do n = 1, 2
            indzh = zlevel + n - 1
            call hor_interpol_4d(field, output1(n), indzh, memind(m), ztot)
         end do

         !**********************************
         ! 2.) Linear vertical interpolation on logarithmic scale
         !**********************************
         call vert_interpol(output1(1), output1(2), dz1, dz2, output(m))

      end do

   end subroutine bilin_spatial_interpol

   subroutine bilin_spatial_interpol_nest(field, output, zlevel, dz1, dz2, ztot)

      integer, intent(in) :: zlevel, ztot       ! interpolation z level
      real, intent(in)    :: field(0:nxmaxn - 1, 0:nymaxn - 1, ztot, numwfmem, numbnests)
      ! input field to interpolate
      real, intent(in)    :: dz1, dz2
      real, intent(inout) :: output(2)         ! interpolated values
      integer             :: m, n, indzh
      real                :: output1(2)

      do m = 1, 2

         do n = 1, 2
            indzh = zlevel + n - 1
            call hor_interpol_4d_nest(field, output1(n), indzh, memind(m), ztot)
         end do

         !**********************************
         ! 2.) Linear vertical interpolation on logarithmic scale
         !**********************************
         call vert_interpol(output1(1), output1(2), dz1, dz2, output(m))

      end do

   end subroutine bilin_spatial_interpol_nest

   subroutine compute_sl_sq(field, sl, sq, zlevel, indexh, ztot)

      integer, intent(in) :: zlevel, ztot, indexh   ! interpolation z levels
      real, intent(in)    :: field(0:nxmax - 1, 0:nymax - 1, ztot, numwfmem)
      ! input field to interpolate
      real, intent(inout) :: sl, sq                ! standard deviation

      sl = sl + field(ix, jy, zlevel, indexh) + field(ixp, jy, zlevel, indexh) &
         + field(ix, jyp, zlevel, indexh) + field(ixp, jyp, zlevel, indexh)

      sq = sq + field(ix, jy, zlevel, indexh)*field(ix, jy, zlevel, indexh) + &
         field(ixp, jy, zlevel, indexh)*field(ixp, jy, zlevel, indexh) + &
         field(ix, jyp, zlevel, indexh)*field(ix, jyp, zlevel, indexh) + &
         field(ixp, jyp, zlevel, indexh)*field(ixp, jyp, zlevel, indexh)

   end subroutine compute_sl_sq

   subroutine compute_sl_sq_nest(field, sl, sq, zlevel, indexh, ztot)
      integer, intent(in) :: zlevel, ztot, indexh   ! interpolation z levels
      real, intent(in)    :: field(0:nxmaxn - 1, 0:nymaxn - 1, ztot, numwfmem, numbnests)
      ! input field to interpolate
      real, intent(inout) :: sl, sq                ! standard deviation

      sl = sl + field(ix, jy, zlevel, indexh, ngrid) + field(ixp, jy, zlevel, indexh, ngrid) &
         + field(ix, jyp, zlevel, indexh, ngrid) + field(ixp, jyp, zlevel, indexh, ngrid)

      sq = sq + field(ix, jy, zlevel, indexh, ngrid)*field(ix, jy, zlevel, indexh, ngrid) + &
         field(ixp, jy, zlevel, indexh, ngrid)*field(ixp, jy, zlevel, indexh, ngrid) + &
         field(ix, jyp, zlevel, indexh, ngrid)*field(ix, jyp, zlevel, indexh, ngrid) + &
         field(ixp, jyp, zlevel, indexh, ngrid)*field(ixp, jyp, zlevel, indexh, ngrid)
   end subroutine compute_sl_sq_nest

   subroutine stdev(sl, sq, divisor, output)

      real, intent(in) :: sl, sq, divisor
      real, intent(out) :: output
      real :: xaux
      real, parameter      :: eps = 1.0e-30

      xaux = sq - sl*sl/divisor

      if (xaux .lt. eps) then
         output = 0.
      else
         output = sqrt(xaux/(divisor - 1.))
      end if

   end subroutine stdev

! Interpolation functions
!************************

   subroutine interpol_pbl(itime, xt, yt, zt, zteta, ithread)
      !                          i   i  i  i
      !*****************************************************************************
      !                                                                            *
      !  This subroutine interpolates everything that is needed for calculating the*
      !  dispersion.                                                               *
      !                                                                            *
      !    Author: A. Stohl                                                        *
      !                                                                            *
      !    16 December 1997                                                        *
      !                                                                            *
      !  Revision March 2005 by AST : all output variables in common block cal-    *
      !                               culation of standard deviation done in this  *
      !                               routine rather than subroutine call in order *
      !                               to save computation time                     *
      !                                                                            *
      !*****************************************************************************
      !                                                                            *
      ! Variables:                                                                 *
      ! itime [s]          current temporal position                               *
      ! memtime(3) [s]     times of the wind fields in memory                      *
      ! xt,yt,zt           coordinates position for which wind data shall be       *
      !                    culated                                                 *
      !                                                                            *
      ! Constants:                                                                 *
      !                                                                            *
      !*****************************************************************************

      use turbulence_mod

      integer, intent(in)  :: ithread ! If OMP, number of the thread, otherwise 1
      integer, intent(in) :: itime
      real, intent(in)    :: xt, yt, zt, zteta
      integer             :: m, n, indexh
      integer             :: iw(2), iweta(2)
      real                :: uh1(2), vh1(2), wh1(2), wetah1(2), rho1(2), rhograd1(2)
      real                :: dz1weta, dz2weta
      real, parameter      :: eps = 1.0e-30

      ! Auxiliary variables needed for interpolation
      real :: ust1(2), wst1(2), oli1(2), oliaux

      !********************************************
      ! Multilinear interpolation in time and space
      !********************************************

      ! ngrid and grid coordinates have already been definded, and are included
      ! in the input (for nested: xtn,ytn; for not nested: xts,yts)

      !************************************************************************
      ! Determine the lower left corner and its distance to the current position
      !*************************************************************************

      call find_grid_distances(xt, yt)

      ! Calculate variables for time interpolation
      !*******************************************

      call find_time_vars(itime)

      !********************************************************
      ! 1. Interpolate u*, w* and Obukhov length for turbulence
      !********************************************************

      ! a) Bilinear horizontal interpolation
      if (ngrid .le. 0) then ! No nest
         do m = 1, 2
            indexh = memind(m)
            call hor_interpol(ustar, ust1(m), 1, memind(m), 1)
            call hor_interpol(wstar, wst1(m), 1, memind(m), 1)
            call hor_interpol(oli, oli1(m), 1, memind(m), 1)
         end do
      else ! Nest
         do m = 1, 2
            indexh = memind(m)
            call hor_interpol_nest(ustarn, ust1(m), 1, memind(m), 1)
            call hor_interpol_nest(wstarn, wst1(m), 1, memind(m), 1)
            call hor_interpol_nest(olin, oli1(m), 1, memind(m), 1)
         end do
      end if

      ! b) Temporal interpolation
      call temporal_interpolation(ust1(1), ust1(2), ust)
      call temporal_interpolation(wst1(1), wst1(2), wst)
      call temporal_interpolation(oli1(1), oli1(2), oliaux)

      if (oliaux .ne. 0.) then
         ol = 1./oliaux
      else
         ol = 99999.
      end if

      ! Within the PBL, only METER coordinates are used
      ! with the exception of mesoscale turbulence,
      ! which uses wsigeta computed in interpol_mesoscale
      !**************************************************

      ! Determine the level below the current position
      !***********************************************
      call find_z_level_meters(zt)

      iw(:) = (/indz, indzp/)

      ! w(eta) velocities are necessary for the Petterssen correction
      !**************************************************************
#ifdef ETA
      call find_z_level_eta(zteta)
      iweta(:) = (/indzeta, indzpeta/)
#endif

      !**************************************
      ! 1.) Bilinear horizontal interpolation
      ! 2.) Temporal interpolation (linear)
      !**************************************

      ! Loop over 2 time steps and indz levels
      !***************************************
      if (ngrid .le. 0) then ! No nest
         do n = 1, 2
            do m = 1, 2
               call hor_interpol(ww, wh1(m), iw(n), memind(m), nzmax)
#ifdef ETA
               call hor_interpol(wweta, wetah1(m), iweta(n), memind(m), nzmax)
#endif
               call hor_interpol(rho, rho1(m), iw(n), memind(m), nzmax)
               call hor_interpol(drhodz, rhograd1(m), iw(n), memind(m), nzmax)
               if (ngrid .lt. 0) then
                  call hor_interpol(uupol, uh1(m), iw(n), memind(m), nzmax)
                  call hor_interpol(vvpol, vh1(m), iw(n), memind(m), nzmax)
               else
                  call hor_interpol(uu, uh1(m), iw(n), memind(m), nzmax)
                  call hor_interpol(vv, vh1(m), iw(n), memind(m), nzmax)
               end if
            end do
            call temporal_interpolation(wh1(1), wh1(2), wprof(iw(n), ithread))
#ifdef ETA
            call temporal_interpolation(wetah1(1), wetah1(2), wprofeta(iweta(n), ithread))
#endif
            call temporal_interpolation(uh1(1), uh1(2), uprof(iw(n), ithread))
            call temporal_interpolation(vh1(1), vh1(2), vprof(iw(n), ithread))
            call temporal_interpolation(rho1(1), rho1(2), rhoprof(iw(n), ithread))
            call temporal_interpolation(rhograd1(1), rhograd1(2), rhogradprof(iw(n), ithread))
         end do
      else ! Nest
         do n = 1, 2
            do m = 1, 2
               call hor_interpol_nest(wwn, wh1(m), iw(n), memind(m), nzmax)
#ifdef ETA
               call hor_interpol_nest(wwetan, wetah1(m), iweta(n), memind(m), nzmax)
#endif
               call hor_interpol_nest(uun, uh1(m), iw(n), memind(m), nzmax)
               call hor_interpol_nest(vvn, vh1(m), iw(n), memind(m), nzmax)
               call hor_interpol_nest(rhon, rho1(m), iw(n), memind(m), nzmax)
               call hor_interpol_nest(drhodzn, rhograd1(m), iw(n), memind(m), nzmax)
            end do
            call temporal_interpolation(wh1(1), wh1(2), wprof(iw(n), ithread))
#ifdef ETA
            call temporal_interpolation(wetah1(1), wetah1(2), wprofeta(iweta(n), ithread))
#endif
            call temporal_interpolation(uh1(1), uh1(2), uprof(iw(n), ithread))
            call temporal_interpolation(vh1(1), vh1(2), vprof(iw(n), ithread))
            call temporal_interpolation(rho1(1), rho1(2), rhoprof(iw(n), ithread))
            call temporal_interpolation(rhograd1(1), rhograd1(2), rhogradprof(iw(n), ithread))

            indzindicator(iw(n), ithread) = .false.
         end do
      end if

      ! Only necessary for the Petterssen correction
#ifdef ETA
      call find_vert_vars(wheight, zteta, indzeta, &
         dz1weta, dz2weta, lbounds_w, .true.)
      call vert_interpol(wprofeta(indzeta, ithread), wprofeta(indzpeta, ithread), &
         dz1weta, dz2weta, weta)
#endif

   end subroutine interpol_pbl

   subroutine interpol_pbl_misslev(ithread)
      !*****************************************************************************
      !                                                                            *
      !  This subroutine interpolates u,v,w, density and density gradients.        *
      !                                                                            *
      !    Author: A. Stohl                                                        *
      !                                                                            *
      !    16 December 1997                                                        *
      !    Update: 2 March 1999                                                    *
      !                                                                            *
      !  Revision March 2005 by AST : all output variables in common block cal-    *
      !                               culation of standard deviation done in this  *
      !                               routine rather than subroutine call in order *
      !                               to save computation time                     *
      !                                                                            *
      !*****************************************************************************
      !                                                                            *
      ! Variables:                                                                 *
      ! n                  level                                                   *
      !                                                                            *
      ! Constants:                                                                 *
      !                                                                            *
      !*****************************************************************************

      integer, intent(in)  :: ithread ! number of OMP thread starting at 1
      real                :: uh1(2), vh1(2), wh1(2), rho1(2), rhograd1(2)
      integer             :: m, n, iw(2)

      ! Within the PBL, only METER coordinates are used
      ! with the exception of mesoscale turbulence,
      ! which uses wsigeta computed in interpol_mesoscale
      !**************************************************

      !********************************************
      ! Multilinear interpolation in time and space
      !********************************************

      iw(:) = (/indz, indzp/)
      do n = 1, 2
         if (indzindicator(iw(n), ithread)) then
            if (ngrid .le. 0) then ! No nest
               do m = 1, 2
                  call hor_interpol(ww, wh1(m), iw(n), memind(m), nzmax)
                  call hor_interpol(rho, rho1(m), iw(n), memind(m), nzmax)
                  call hor_interpol(drhodz, rhograd1(m), iw(n), memind(m), nzmax)
                  if (ngrid .lt. 0) then
                     call hor_interpol(uupol, uh1(m), iw(n), memind(m), nzmax)
                     call hor_interpol(vvpol, vh1(m), iw(n), memind(m), nzmax)
                  else
                     call hor_interpol(uu, uh1(m), iw(n), memind(m), nzmax)
                     call hor_interpol(vv, vh1(m), iw(n), memind(m), nzmax)
                  end if
               end do
            else ! Nest
               do m = 1, 2
                  call hor_interpol_nest(wwn, wh1(m), iw(n), memind(m), nzmax)
                  call hor_interpol_nest(uun, uh1(m), iw(n), memind(m), nzmax)
                  call hor_interpol_nest(vvn, vh1(m), iw(n), memind(m), nzmax)
                  call hor_interpol_nest(rhon, rho1(m), iw(n), memind(m), nzmax)
                  call hor_interpol_nest(drhodzn, rhograd1(m), iw(n), memind(m), nzmax)
               end do
            end if
            call temporal_interpolation(wh1(1), wh1(2), wprof(iw(n), ithread))
            call temporal_interpolation(uh1(1), uh1(2), uprof(iw(n), ithread))
            call temporal_interpolation(vh1(1), vh1(2), vprof(iw(n), ithread))
            call temporal_interpolation(rho1(1), rho1(2), rhoprof(iw(n), ithread))
            call temporal_interpolation(rhograd1(1), rhograd1(2), rhogradprof(iw(n), ithread))

            indzindicator(iw(n), ithread) = .false.
         end if
      end do
   end subroutine interpol_pbl_misslev

   subroutine interpol_pbl_short(zt, rhoa, rhograd, ithread)
      implicit none

      integer, intent(in)  :: ithread ! number of OMP thread starting at 1
      real, intent(in)    :: zt
      real, intent(inout) :: rhoa, rhograd
      real                :: dz1, dz2

      call find_vert_vars(height, zt, indz, dz1, dz2, lbounds, .false.)

      call vert_interpol(wprof(indz, ithread), wprof(indzp, ithread), dz1, dz2, w)
      call vert_interpol(uprof(indz, ithread), uprof(indzp, ithread), dz1, dz2, u)
      call vert_interpol(vprof(indz, ithread), vprof(indzp, ithread), dz1, dz2, v)
      call vert_interpol(rhoprof(indz, ithread), rhoprof(indzp, ithread), dz1, dz2, rhoa)
      call vert_interpol(rhogradprof(indz, ithread), rhogradprof(indzp, ithread), dz1, dz2, rhograd)
   end subroutine interpol_pbl_short

   subroutine interpol_mesoscale(xt, yt, zt, zteta)

      use turbulence_mod

      real, intent(in)    :: xt, yt, zt, zteta
#ifdef ETA
      integer             :: iuv(2), iweta(2)
#endif
      integer             :: iw(2)

      ! Where in the grid? Stereographic (ngrid<0) or nested (ngrid>0)
      !***************************************************************
      call find_ngrid(xt, yt)

      call find_grid_indices(xt, yt)

      ! Determine the level below the current position
      !***********************************************
      call find_z_level_meters(zt)
      iw(:) = (/indz, indzp/)

#ifdef ETA
      call find_z_level_eta(zteta)
      iuv(:) = (/induv, indpuv/)
      iweta(:) = (/indzeta, indzpeta/)
      call stdev_eta(iw, iuv, iweta)
#else
      iw(:) = (/indz, indzp/)
      call stdev_meter(iw)
#endif

   end subroutine interpol_mesoscale

   subroutine interpol_wind(itime, xt, yt, zt, zteta)
      !                           i   i  i  i

      !*****************************************************************************
      !                                                                            *
      !  This subroutine interpolates the wind data to current trajectory position.*
      !                                                                            *
      !    Author: A. Stohl                                                        *
      !                                                                            *
      !    16 December 1997                                                        *
      !                                                                            *
      !  Revision March 2005 by AST : all output variables in common block cal-    *
      !                               culation of standard deviation done in this  *
      !                               routine rather than subroutine call in order *
      !                               to save computation time                     *
      !                                                                            *
      !*****************************************************************************
      !                                                                            *
      ! Variables:                                                                 *
      ! u,v,w              wind components                                         *
      ! itime [s]          current temporal position                               *
      ! memtime(3) [s]     times of the wind fields in memory                      *
      ! xt,yt,zt           coordinates position for which wind data shall be       *
      !                    calculated                                              *
      !                                                                            *
      ! Constants:                                                                 *
      !                                                                            *
      !*****************************************************************************

      integer, intent(in) :: itime
      real, intent(in)    :: xt, yt, zt, zteta
#ifdef ETA
      integer             :: iuv(2), iweta(2)
#else
      integer             :: iw(2)
#endif

      ! Where in the grid? Stereographic (ngrid<0) or nested (ngrid>0)
      !***************************************************************
      call find_ngrid(xt, yt)

      call find_grid_indices(xt, yt)
      ! ! Multilinear interpolation in time and space
      ! !********************************************

      ! Determine the lower left corner and its distance to the current position
      !*************************************************************************
      call find_grid_distances(xt, yt)

      ! Calculate variables for time interpolation
      !*******************************************
      call find_time_vars(itime)

      ! Interpolate over the windfields depending on the prefered
      ! coordinate system
      !**********************************************************
#ifdef ETA
      ! Same for eta coordinates
      !*************************
      call find_z_level_eta(zteta)

      iuv(:) = (/induv, indpuv/)
      iweta(:) = (/indzeta, indzpeta/)
      call interpol_wind_eta(zteta, iuv, iweta)
      !call stdev_wind_eta(iw,iuv,iweta)
#else
      ! Determine the level below the current position for u,v
      !*******************************************************
      call find_z_level_meters(zt)

      iw(:) = (/indz, indzp/)
      call interpol_wind_meter(zt, iw)
      !call stdev_wind_meter(iw)
#endif

   end subroutine interpol_wind

   subroutine interpol_wind_short(itime, xt, yt, zt, zteta)
!                                i   i  i  i  i

      !*****************************************************************************
      !                                                                            *
      !  This subroutine interpolates the wind data to current trajectory position.*
      !                                                                            *
      !    Author: A. Stohl                                                        *
      !                                                                            *
      !    16 December 1997                                                        *
      !                                                                            *
      !  Revision March 2005 by AST : all output variables in common block         *
      !                                                                            *
      !*****************************************************************************
      !                                                                            *
      ! Variables:                                                                 *
      ! u,v,w              wind components                                         *
      ! itime [s]          current temporal position                               *
      ! memtime(3) [s]     times of the wind fields in memory                      *
      ! xt,yt,zt           coordinates position for which wind data shall be       *
      !                    calculated                                              *
      !                                                                            *
      ! Constants:                                                                 *
      !                                                                            *
      !*****************************************************************************

      integer, intent(in) :: itime
      real, intent(in) :: xt, yt, zt, zteta
#ifdef ETA
      integer             :: iuv(2), iweta(2)
#else
      integer             :: iw(2)
#endif

      !********************************************
      ! Multilinear interpolation in time and space
      !********************************************

      ! Where in the grid? Stereographic (ngrid<0) or nested (ngrid>0)
      !***************************************************************
      call find_ngrid(xt, yt)
      call find_grid_indices(xt, yt)
      call find_grid_distances(xt, yt)

      ! Calculate variables for time interpolation
      !*******************************************
      call find_time_vars(itime)

      ! Interpolate over the windfields depending on the prefered
      ! coordinate system
      !**********************************************************
#ifdef ETA
      ! Determine the level below the current position for eta coordinates
      !*******************************************************************
      call find_z_level_eta(zteta)

      iuv(:) = (/induv, indpuv/)
      iweta(:) = (/indzeta, indzpeta/)
      ! Interpolate the u, v, weta windfields
      !**************************************
      call interpol_wind_eta(zteta, iuv, iweta)
#else

      ! Determine the level below the current position for u,v
      !*******************************************************
      call find_z_level_meters(zt)

      iw(:) = (/indz, indzp/)
      call interpol_wind_meter(zt, iw)
#endif

   end subroutine interpol_wind_short

   subroutine interpol_partoutput_val(fieldname, output, j)
      integer, intent(in)         :: j          ! particle number
      character(2), intent(in)    :: fieldname  ! input field to interpolate over
      real, intent(inout)         :: output

      ! Interpolate over the windfields depending on the prefered
      ! coordinate system
      !**********************************************************
#ifdef ETA
      call interpol_partoutput_val_eta(fieldname, output, j)
#else
      call interpol_partoutput_val_meter(fieldname, output, j)
#endif
   end subroutine interpol_partoutput_val

   subroutine interpol_htropo_hmix(tropop, h)

      real, intent(inout) :: tropop   ! height of troposphere
      real, intent(inout) :: h        ! mixing height
      real    :: h1(2)                ! mixing height of 2 timesteps
      integer :: mind                 ! windfield index
      integer :: i, j, k, m              ! loop variables

      h = 0.
      if (ngrid .le. 0) then
         if (interpolhmix) then
            do m = 1, 2
               call hor_interpol(hmix, h1(m), 1, memind(m), 1)
            end do
         else
            do k = 1, 2
               mind = memind(k) ! eso: compatibility with 3-field version
               do j = jy, jyp
                  do i = ix, ixp
                     if (hmix(i, j, 1, mind) .gt. h) h = hmix(i, j, 1, mind)
                  end do
               end do
            end do
         end if
         tropop = tropopause(nix, njy, 1, memind(1))
      else
         do k = 1, 2
            mind = memind(k)
            do j = jy, jyp
               do i = ix, ixp
                  if (hmixn(i, j, 1, mind, ngrid) .gt. h) h = hmixn(i, j, 1, mind, ngrid)
               end do
            end do
         end do
         tropop = tropopausen(nix, njy, 1, memind(1), ngrid)
      end if

      if (interpolhmix) h = (h1(1)*dt2 + h1(2)*dt1)*dtt

   end subroutine interpol_htropo_hmix

   subroutine interpol_density(itime, ipart, output)

      integer, intent(in) :: itime, ipart  ! time and particle index
      real, intent(inout) :: output ! output density (rhoi)
      integer :: ind
      real :: dz1, dz2, rhoprof(2)

      ! Where in the grid? Stereographic (ngrid<0) or nested (ngrid>0)
      !***************************************************************
      call find_ngrid(part(ipart)%xlon, part(ipart)%ylat)
      call find_grid_indices(real(part(ipart)%xlon), real(part(ipart)%ylat))
      call find_grid_distances(real(part(ipart)%xlon), real(part(ipart)%ylat))
      call find_time_vars(itime)

      ! Take density from 2nd wind field in memory
      !(accurate enough, no time interpolation needed)
      !***********************************************

#ifdef ETA
      call find_z_level_eta(real(part(ipart)%zeta))
      call find_vert_vars(uvheight, real(part(ipart)%zeta), induv, &
         dz1, dz2, lbounds_uv, .false.)
      if (ngrid .le. 0) then
         do ind = induv, indpuv
            call hor_interpol(rhoeta, rhoprof(ind - induv + 1), ind, memind(2), nzmax)
         end do
      else
         do ind = induv, indpuv
            call hor_interpol_nest(rhoetan, rhoprof(ind - induv + 1), ind, memind(2), &
               nzmax)
         end do
      end if
#else
      call find_z_level_meters(real(part(ipart)%z))
      call find_vert_vars(height, real(part(ipart)%z), indz, &
         dz1, dz2, lbounds, .false.)
      if (ngrid .le. 0) then
         do ind = indz, indzp
            call hor_interpol(rho, rhoprof(ind - indz + 1), ind, memind(2), nzmax)
         end do
      else
         do ind = indz, indzp
            call hor_interpol_nest(rhon, rhoprof(ind - indz + 1), ind, memind(2), nzmax)
         end do
      end if
#endif
      call vert_interpol(rhoprof(1), rhoprof(2), dz1, dz2, output)
   end subroutine interpol_density

   subroutine interpol_rain(itime, kz, yint1, yint2, yint3, ytint, yint4, intiy1, intiy2, icmv)
      !                       i     i   o     o     o     o     o     o      o     i
      !****************************************************************************
      !                                                                           *
      !  Interpolation of meteorological fields on 2-d model layers.              *
      !  In horizontal direction bilinear interpolation is used.                  *
      !  Temporally a linear interpolation is used.                               *
      !  Seven fields are interpolated at the same time.                          *
      !                                                                           *
      !  This is a special version of levlininterpol to save CPU time.            *
      !                                                                           *
      !  1 first time                                                             *
      !  2 second time                                                            *
      !                                                                           *
      !                                                                           *
      !     Author: A. Stohl                                                      *
      !                                                                           *
      !     30 August 1996                                                        *
      !                                                                           *
      !                                                                           *
      ! PS, AP 04/2019, 11/2020:                                                  *
      !                 put back temporal interpolation of rain, from v10.01      *
      !      and cloud bottom / thickness interpolation                           *
      ! PS, AP 01/2021:                                                           *
      !      interpolate particle temperature and cloud total water               *
      ! PS, AP 02/2021:                                                           *
      !      interpolation of precipitation using two additional fields           *
      !      which are temporally equidistant between the main fields             *
      !                                                                           *
      !****************************************************************************
      !                                                                           *
      ! Variables:                                                                *
      !                                                                           *
      ! dt1,dt2              time differences between fields and current position *
      ! dz1,dz2              z distance between levels and current position       *
      ! height(nzmax)        heights of the model levels                          *
      ! mm                   help variable                                        *
      ! indz                 the level closest to the current trajectory position *
      ! indzh                help variable                                        *
      ! itime                current time                                         *
      ! ix,jy                x,y coordinates of lower left subgrid point          *
      ! level                level at which interpolation shall be done 2d, =1    *
      ! kz                   level at which interpolation shall be done           *
      ! memind(3)            points to the places of the wind fields              *
      ! nx,ny                actual field dimensions in x,y and z direction       *
      ! nxmax,nymax,nzmax    maximum field dimensions in x,y and z direction      *
      ! xt                   current x coordinate                                 *
      ! yint                 the final interpolated value                         *
      ! yt                   current y coordinate                                 *
      ! yy?(0:nxmax,0:nymax,nzx2d,3) meteorological field used for interpolation  *
      ! yyt(0:nxmax,0:nymax,nzmax,3) tt field                                     *
      ! iy1,iy2(0:nxmax,0:nymax,3) cloud bottom, thickness fields (integer)       *
      ! zt                   current z coordinate                                 *
      !                                                                           *
      !****************************************************************************
      use par_mod, only: numwfmem, numpf
      use com_mod, only: lcw

      implicit none

      integer, intent(in) :: kz, icmv, itime
      integer, intent(out) :: intiy1, intiy2
      real, intent(out) :: yint1, yint2, yint3, ytint, yint4
      integer :: m
      integer :: mm
      real :: ip1, ip2, ip3, ip4, ipsum
      real :: dt, dtp1, dtp2, rt
      real :: y1(2), y2(2), y3(2), y4(2), yi1(2), yi2(2), ytt(2) ! interpolated values

      integer, dimension(2) :: ip, mp

      !**********************************************************************
      ! 1.) Bilinear horizontal interpolation
      ! This has to be done separately for 2 fields (Temporal)
      !*******************************************************
      ! Loop over 2 time steps
      !***********************

      !-------------------------------------------------------------------------
      ! PS, AT new interpolation of precip with 2 additional fields
      ! therefore, we need a special treatment of lsp,cp which are in yy1,yy2
      !-------------------------------------------------------------------------
      !
      !      1.1 1.2 1.3               ip(1).mp(1)
      !      1.2 1.3 2.1               ip(2).mp(2)
      !
      !   ||___|___|___||___|___|___||
      !
      ! ip  1   2   3    1   2   3    1
      ! m        1            2
      !
      !-------------------------------------------------------------------------

      dt1 = real(itime - memtime(1))
      dt2 = real(memtime(2) - itime)
      dt = real(memtime(2) - memtime(1))
      if (dt .eq. 0. .and. dt1 .eq. 0.) then ! Fix if last last timestep and memtime(2)=memtime(1)
         dt = 1.
         dt2 = 1.
      end if
      dtt = dt/3.
      if (numpf .eq. 1) then
         mp(1) = 1
         mp(2) = 2
         ip(1) = 1
         ip(2) = 1
         dtp1 = dt1
         dtp2 = dt2
      else
         rt = abs(dt1/dt)
         if (0 .le. rt .and. rt .lt. 1./3.) then
            mp(1) = 1
            mp(2) = 1
            ip(1) = 1
            ip(2) = 2
            dtp1 = dt1
            dtp2 = dt2 - 2.*dtt
         elseif (1./3. .le. rt .and. rt .lt. 2./3.) then
            mp(1) = 1
            mp(2) = 1
            ip(1) = 2
            ip(2) = 3
            dtp1 = dt1 - dtt
            dtp2 = dt2 - dtt
         elseif (2./3. .le. rt .and. rt .lt. 1.) then
            mp(1) = 1
            mp(2) = 2
            ip(1) = 3
            ip(2) = 1
            dtp1 = dt1 - 2.*dtt
            dtp2 = dt2
         end if
      end if

      if (ngrid .le. 0) then ! No nest
         do m = 1, 2
            mm = memind(mp(m))
            y1(m) = p1*lsprec(ix, jy, 1, ip(m), mm) &
               + p2*lsprec(ixp, jy, 1, ip(m), mm) &
               + p3*lsprec(ix, jyp, 1, ip(m), mm) &
               + p4*lsprec(ixp, jyp, 1, ip(m), mm)
            y2(m) = p1*convprec(ix, jy, 1, ip(m), mm) &
               + p2*convprec(ixp, jy, 1, ip(m), mm) &
               + p3*convprec(ix, jyp, 1, ip(m), mm) &
               + p4*convprec(ixp, jyp, 1, ip(m), mm)

            mm = memind(m)
            y3(m) = p1*tcc(ix, jy, 1, mm) &
               + p2*tcc(ixp, jy, 1, mm) &
               + p3*tcc(ix, jyp, 1, mm) &
               + p4*tcc(ixp, jyp, 1, mm)
#ifdef ETA
            ytt(m) = p1*tteta(ix, jy, kz, mm) &
               + p2*tteta(ixp, jy, kz, mm) &
               + p3*tteta(ix, jyp, kz, mm) &
               + p4*tteta(ixp, jyp, kz, mm)
#else
            ytt(m) = p1*tt(ix, jy, kz, mm) &
               + p2*tt(ixp, jy, kz, mm) &
               + p3*tt(ix, jyp, kz, mm) &
               + p4*tt(ixp, jyp, kz, mm)
#endif
            if (lcw) &
               y4(m) = p1*ctwc(ix, jy, mm) &
               + p2*ctwc(ixp, jy, mm) &
               + p3*ctwc(ix, jyp, mm) &
               + p4*ctwc(ixp, jyp, mm)

            !PS clouds:
            ip1 = 1.
            ip2 = 1.
            ip3 = 1.
            ip4 = 1.
            ipsum = 1.
            if (icloudbot(ix, jy, mm) .eq. icmv) then
               ip1 = 0.
               ipsum = ipsum - p1
            end if
            if (icloudbot(ixp, jy, mm) .eq. icmv) then
               ip2 = 0.
               ipsum = ipsum - p2
            end if
            if (icloudbot(ix, jyp, mm) .eq. icmv) then
               ip3 = 0.
               ipsum = ipsum - p3
            end if
            if (icloudbot(ixp, jyp, mm) .eq. icmv) then
               ip4 = 0.
               ipsum = ipsum - p4
            end if
            if (ipsum .eq. 0.) then
               yi1(m) = icmv
            else
               yi1(m) = (ip1*p1*icloudbot(ix, jy, mm) &
                  + ip2*p2*icloudbot(ixp, jy, mm) &
                  + ip3*p3*icloudbot(ix, jyp, mm) &
                  + ip4*p4*icloudbot(ixp, jyp, mm))/ipsum
               ! AP test output
               !      if (yi1(m) .lt. 1.) then
               !        write(*,*) ip1,ip2,ip3,ip4,ipsum
               !        write(*,*) p1,p2,p3,p4
               !        write(*,*) iy1(ix ,jy ,mm), iy1(ixp,jy ,mm),  iy1(ix ,jyp,mm), iy1(ixp,jyp,mm)
               !      endif

            end if

            ip1 = 1.
            ip2 = 1.
            ip3 = 1.
            ip4 = 1.
            ipsum = 1.
            if (icloudtop(ix, jy, mm) .eq. icmv) then
               ip1 = 0.
               ipsum = ipsum - p1
            end if
            if (icloudtop(ixp, jy, mm) .eq. icmv) then
               ip2 = 0.
               ipsum = ipsum - p2
            end if
            if (icloudtop(ix, jyp, mm) .eq. icmv) then
               ip3 = 0.
               ipsum = ipsum - p3
            end if
            if (icloudtop(ixp, jyp, mm) .eq. icmv) then
               ip4 = 0.
               ipsum = ipsum - p4
            end if
            if (ipsum .eq. 0.) then
               yi2(m) = icmv
            else
               yi2(m) = (ip1*p1*icloudtop(ix, jy, mm) &
                  + ip2*p2*icloudtop(ixp, jy, mm) &
                  + ip3*p3*icloudtop(ix, jyp, mm) &
                  + ip4*p4*icloudtop(ixp, jyp, mm))/ipsum
            end if
            !PS end clouds
         end do
      else ! Nest
         do m = 1, 2
            mm = memind(mp(m))
            y1(m) = p1*lsprecn(ix, jy, 1, ip(m), mm, ngrid) &
               + p2*lsprecn(ixp, jy, 1, ip(m), mm, ngrid) &
               + p3*lsprecn(ix, jyp, 1, ip(m), mm, ngrid) &
               + p4*lsprecn(ixp, jyp, 1, ip(m), mm, ngrid)
            y2(m) = p1*convprecn(ix, jy, 1, ip(m), mm, ngrid) &
               + p2*convprecn(ixp, jy, 1, ip(m), mm, ngrid) &
               + p3*convprecn(ix, jyp, 1, ip(m), mm, ngrid) &
               + p4*convprecn(ixp, jyp, 1, ip(m), mm, ngrid)

            mm = memind(m)
            y3(m) = p1*tccn(ix, jy, 1, mm, ngrid) &
               + p2*tccn(ixp, jy, 1, mm, ngrid) &
               + p3*tccn(ix, jyp, 1, mm, ngrid) &
               + p4*tccn(ixp, jyp, 1, mm, ngrid)
#ifdef ETA
            ytt(m) = p1*ttetan(ix, jy, kz, mm, ngrid) &
               + p2*ttetan(ixp, jy, kz, mm, ngrid) &
               + p3*ttetan(ix, jyp, kz, mm, ngrid) &
               + p4*ttetan(ixp, jyp, kz, mm, ngrid)
#else
            ytt(m) = p1*ttn(ix, jy, kz, mm, ngrid) &
               + p2*ttn(ixp, jy, kz, mm, ngrid) &
               + p3*ttn(ix, jyp, kz, mm, ngrid) &
               + p4*ttn(ixp, jyp, kz, mm, ngrid)
#endif
            if (lcw_nest(ngrid)) &
               y4(m) = p1*ctwcn(ix, jy, mm, ngrid) &
               + p2*ctwcn(ixp, jy, mm, ngrid) &
               + p3*ctwcn(ix, jyp, mm, ngrid) &
               + p4*ctwcn(ixp, jyp, mm, ngrid)

            !PS clouds:
            ip1 = 1.
            ip2 = 1.
            ip3 = 1.
            ip4 = 1.
            ipsum = 1.
            if (icloudbotn(ix, jy, mm, ngrid) .eq. icmv) then
               ip1 = 0.
               ipsum = ipsum - p1
            end if
            if (icloudbotn(ixp, jy, mm, ngrid) .eq. icmv) then
               ip2 = 0.
               ipsum = ipsum - p2
            end if
            if (icloudbotn(ix, jyp, mm, ngrid) .eq. icmv) then
               ip3 = 0.
               ipsum = ipsum - p3
            end if
            if (icloudbotn(ixp, jyp, mm, ngrid) .eq. icmv) then
               ip4 = 0.
               ipsum = ipsum - p4
            end if
            if (ipsum .eq. 0.) then
               yi1(m) = icmv
            else
               yi1(m) = (ip1*p1*icloudbotn(ix, jy, mm, ngrid) &
                  + ip2*p2*icloudbotn(ixp, jy, mm, ngrid) &
                  + ip3*p3*icloudbotn(ix, jyp, mm, ngrid) &
                  + ip4*p4*icloudbotn(ixp, jyp, mm, ngrid))/ipsum
            end if

            ip1 = 1.
            ip2 = 1.
            ip3 = 1.
            ip4 = 1.
            ipsum = 1.
            if (icloudtopn(ix, jy, mm, ngrid) .eq. icmv) then
               ip1 = 0.
               ipsum = ipsum - p1
            end if
            if (icloudtopn(ixp, jy, mm, ngrid) .eq. icmv) then
               ip2 = 0.
               ipsum = ipsum - p2
            end if
            if (icloudtopn(ix, jyp, mm, ngrid) .eq. icmv) then
               ip3 = 0.
               ipsum = ipsum - p3
            end if
            if (icloudtopn(ixp, jyp, mm, ngrid) .eq. icmv) then
               ip4 = 0.
               ipsum = ipsum - p4
            end if
            if (ipsum .eq. 0.) then
               yi2(m) = icmv
            else
               yi2(m) = (ip1*p1*icloudtopn(ix, jy, mm, ngrid) &
                  + ip2*p2*icloudtopn(ixp, jy, mm, ngrid) &
                  + ip3*p3*icloudtopn(ix, jyp, mm, ngrid) &
                  + ip4*p4*icloudtopn(ixp, jyp, mm, ngrid))/ipsum
            end if
            !PS end clouds
         end do
      end if

      !************************************
      ! 2.) Temporal interpolation (linear)
      !************************************

      yint1 = (y1(1)*dtp2 + y1(2)*dtp1)/dtt ! lsp
      yint2 = (y2(1)*dtp2 + y2(2)*dtp1)/dtt ! cp

      yint3 = (y3(1)*dt2 + y3(2)*dt1)/dt
      yint4 = (y4(1)*dt2 + y4(2)*dt1)/dt
      ytint = (ytt(1)*dt2 + ytt(2)*dt1)/dt

!PS clouds:450.

      ! write (*, *) yi1(1), yi1(2), yi2(1), yi2(2), dt, dt1, dt2
      intiy1 = int((yi1(1)*dt2 + yi1(2)*dt1)/dt)
      if (int(yi1(1)) .eq. icmv) intiy1 = int(yi1(2))
      if (int(yi1(2)) .eq. icmv) intiy1 = int(yi1(1))

      intiy2 = int((yi2(1)*dt2 + yi2(2)*dt1)/dt)
      if (int(yi2(1)) .eq. icmv) intiy2 = int(yi2(2))
      if (int(yi2(2)) .eq. icmv) intiy2 = int(yi2(1))

!  write(*,*) 'before cbot: ', intiy1, ' cthick: ', intiy2
      if (intiy1 .ne. icmv .and. intiy2 .ne. icmv) then
         intiy2 = intiy2 !intiy1 + intiy2 ! convert cloud thickness to cloud top
      else
         intiy1 = icmv
         intiy2 = icmv
      end if
      if (intiy2 .ne. icmv .and. intiy2 .lt. 0) then
         write (*, *) itime, memind(1), memind(2)
         write (*, *) yi1(1), yi1(2), yi2(1), yi2(2), dt, dt1, dt2
         write (*, *) 'final cbot: ', intiy1, ' ctop: ', intiy2
         stop 'intiy2 (cloud top) negative'
      end if
!PS end clouds

   end subroutine interpol_rain

!!! VIABILITY FUNCTIONS
#ifdef VIA

   subroutine init_viability_model(viability_fit, interp_via_method)
      ! This routine assigns the pointer ONCE at startup.
      implicit none
      integer, intent(in) :: viability_fit, interp_via_method

      select case (viability_fit)
       case (0)
         decay_calc => calculate_decay_constant_exponential
       case (1)
         decay_calc => calculate_decay_isard
       case (2)
         decay_calc => calculate_decay_two_stage
       case (3)
         decay_calc => calculate_decay_two_stage_shoulder
       case default
         print *, "FLEXPART ERROR: Unknown viability_fit", viability_fit
         stop
      end select

      select case (interp_via_method)
       case(0)
         interp_via_calc => null()
       case(1)
         interp_via_calc => interpol_viability_isard
       case(2)
         interp_via_calc => interpol_viability_agpads
       case(3)
         interp_via_calc => interpol_viability_full
       case default
         print *, "Constant exponential per/hr decay rate with value", real(interp_via_method) / 1000.0
      end select

   end subroutine init_viability_model

   subroutine interpol_viability_full(itime, jpart, yint_sr, yint_rh, ytint, icmv)
      !                             i      i       o        o       o      i
      !****************************************************************************
      !                                                                           *
      !  Interpolation of meteorological fields for particle viability.           *
      !  In horizontal direction bilinear interpolation is used.                  *
      !  Temporally a linear interpolation is used.                               *
      !                                                                           *
      !  Special version adapted from levlininterpol / interpol_rain for particle  *
      !  viability calculations. Interpolates solar radiation (attenuated through *
      !  cloud layers), relative humidity, and temperature.                       *
      !                                                                           *
      !     Author: A. Stohl                                                      *
      !                                                                           *
      !     30 August 1996                                                        *
      !                                                                           *
      ! PS, AP 04/2019, 11/2020:                                                  *
      !                 put back temporal interpolation of rain, from v10.01      *
      !      and cloud bottom / thickness interpolation                           *
      ! PS, AP 01/2021:                                                           *
      !      interpolate particle temperature and cloud total water               *
      ! PS, AP 02/2021:                                                           *
      !      interpolation of precipitation using two additional fields           *
      !      which are temporally equidistant between the main fields             *
      !                                                                           *
      ! AJ 07/2026 - Viability interpolation updates:                                          *
      !      adapted for full viability model: interpolates solar radiation with  *
      !      cloud layer attenuation, relative humidity, and air temperature      *
      !                                                                           *
      !****************************************************************************
      !                                                                           *
      ! Variables:                                                                *
      !                                                                           *
      ! itime                current simulation time [s]                          *
      ! jpart                index of the particle to be interpolated             *
      ! yint_sr              interpolated solar radiation at particle [W/m2]      *
      ! yint_rh              interpolated relative humidity at particle [%]       *
      ! ytint                interpolated temperature at particle [K]             *
      ! icmv                 viability model control / component flag             *
      !                                                                           *
      ! dt, dtp1, dtp2       time differences / weights for temporal interpol.    *
      ! ip1,ip2,ip3,ip4      bilinear interpolation weights                       *
      ! kz                   vertical model level of the particle                 *
      ! memind(3)            points to the places of the wind/met fields          *
      ! xts, yts             particle coordinates in grid units                   *
      ! yrh(2)               relative humidity at bounding time steps             *
      ! ysr(2)               solar radiation at bounding time steps               *
      ! ytt(2)               temperature at bounding time steps                   *
      ! interp_ssr           surface solar radiation                              *
      ! toa_rad              top-of-atmosphere solar radiation                    *
      ! eta_part             vertical eta coordinate of the particle              *
      ! eta_toc, eta_boc     eta coordinates of cloud top and cloud bottom        *
      ! ictop, icbot         vertical grid indices for cloud top and bottom       *
      ! fdepth               fractional cloud optical/physical depth              *
      ! Kt                   clearness index / cloud transmission factor          *
      !                                                                           *
      !****************************************************************************
      use par_mod, only: numwfmem, numpf
      use com_mod, only: lcw
      use qvsat_mod
      use date_mod, only: caldate

      implicit none

      integer, intent(in) :: icmv, itime, jpart
      real, intent(out) :: yint_rh, yint_sr, ytint
      integer :: m
      integer :: mm
      integer :: ictop, icbot, jjjjmmdd, hhmmss, kz
      real :: ip1, ip2, ip3, ip4, ipsum
      real :: dt, dtp1, dtp2, rt
      real :: yrh(2), ysr(2), yi1(2), yi2(2), ytt(2) ! interpolated values
      real :: interp_ssr, toa_rad, lon_deg, lat_deg, eta_part, eta_toc, eta_boc, fdepth, Kt, xts, yts
      real(kind=dp) :: jul
      logical :: below_top, above_bot

      integer, dimension(2) :: ip, mp

      !**********************************************************************
      ! 1.) Bilinear horizontal interpolation
      ! This has to be done separately for 2 fields (Temporal)
      !*******************************************************
      ! Loop over 2 time steps
      !***********************

      !-------------------------------------------------------------------------
      ! PS, AT new interpolation of precip with 2 additional fields
      ! therefore, we need a special treatment of lsp,cp which are in yy1,yy2
      !-------------------------------------------------------------------------
      !
      !      1.1 1.2 1.3               ip(1).mp(1)
      !      1.2 1.3 2.1               ip(2).mp(2)
      !
      !   ||___|___|___||___|___|___||
      !
      ! ip  1   2   3    1   2   3    1
      ! m        1            2
      !
      !-------------------------------------------------------------------------

#ifdef ETA
      call find_z_level_eta_uv_safe(real(part(jpart)%zeta), kz)
!      kz = induv
#else
      call find_z_level_meters_safe(real(part(jpart)%z), kz)
!      kz = indz
#endif
      xts = real(part(jpart)%xlon)
      yts = real(part(jpart)%ylat)


      if (xts .ge. real(nx - 1)) xts = real(nx - 1) - 0.00001
      if (yts .ge. real(ny - 1)) yts = real(ny - 1) - 0.00001

      call find_ngrid(xts, yts)
      call find_grid_indices(xts, yts)
      call find_grid_distances(xts, yts)
      dt1 = real(itime - memtime(1))
      dt2 = real(memtime(2) - itime)
      dt = real(memtime(2) - memtime(1))
      ! Wrap Longitude for the 360-degree line
      if (ixp .ge. nx) ixp = ixp - nx

      ! Clamp Latitude at the poles
      if (jyp .ge. ny) jyp = ny - 1

      ! Clamp Vertical to prevent top-of-atmosphere buffer read
      kz = max(1, min(kz, nz))

      dt1 = real(itime - memtime(1))
      dt2 = real(memtime(2) - itime)
      dt = real(memtime(2) - memtime(1))
      if (dt .eq. 0. .and. dt1 .eq. 0.) then ! Fix if last last timestep and memtime(2)=memtime(1)
         dt = 1.
         dt2 = 1.
      end if
      dtt = dt/3.
      if (numpf .eq. 1) then
         mp(1) = 1
         mp(2) = 2
         ip(1) = 1
         ip(2) = 1
         dtp1 = dt1
         dtp2 = dt2
      else
         rt = abs(dt1/dt)
         if (0 .le. rt .and. rt .lt. 1./3.) then
            mp(1) = 1
            mp(2) = 1
            ip(1) = 1
            ip(2) = 2
            dtp1 = dt1
            dtp2 = dt2 - 2.*dtt
         elseif (1./3. .le. rt .and. rt .lt. 2./3.) then
            mp(1) = 1
            mp(2) = 1
            ip(1) = 2
            ip(2) = 3
            dtp1 = dt1 - dtt
            dtp2 = dt2 - dtt
         elseif (2./3. .le. rt .and. rt .lt. 1.) then
            mp(1) = 1
            mp(2) = 2
            ip(1) = 3
            ip(2) = 1
            dtp1 = dt1 - 2.*dtt
            dtp2 = dt2
         end if
      end if

      if (ngrid .le. 0) then ! No nest
         do m = 1, 2
            mm = memind(m)
            yrh(m) = p1*calc_rh(ix, jy, kz, mm) &
               + p2*calc_rh(ixp, jy, kz, mm) &
               + p3*calc_rh(ix, jyp, kz, mm) &
               + p4*calc_rh(ixp, jyp, kz, mm)
#ifdef ETA
            ytt(m) = p1*tteta(ix, jy, kz, mm) &
               + p2*tteta(ixp, jy, kz, mm) &
               + p3*tteta(ix, jyp, kz, mm) &
               + p4*tteta(ixp, jyp, kz, mm)
#else
            ytt(m) = p1*tt(ix, jy, kz, mm) &
               + p2*tt(ixp, jy, kz, mm) &
               + p3*tt(ix, jyp, kz, mm) &
               + p4*tt(ixp, jyp, kz, mm)
#endif
            ip1 = 1.
            ip2 = 1.
            ip3 = 1.
            ip4 = 1.
            ipsum = 1.

            !PS clouds:
            if (icloudbot(ix, jy, mm) .eq. icmv) then
               ip1 = 0.
               ipsum = ipsum - p1
            end if
            if (icloudbot(ixp, jy, mm) .eq. icmv) then
               ip2 = 0.
               ipsum = ipsum - p2
            end if
            if (icloudbot(ix, jyp, mm) .eq. icmv) then
               ip3 = 0.
               ipsum = ipsum - p3
            end if
            if (icloudbot(ixp, jyp, mm) .eq. icmv) then
               ip4 = 0.
               ipsum = ipsum - p4
            end if
            if (ipsum .eq. 0.) then
               yi1(m) = icmv
            else
               yi1(m) = (ip1*p1*icloudbot(ix, jy, mm) &
                  + ip2*p2*icloudbot(ixp, jy, mm) &
                  + ip3*p3*icloudbot(ix, jyp, mm) &
                  + ip4*p4*icloudbot(ixp, jyp, mm))/ipsum
               ! AP test output
               !      if (yi1(m) .lt. 1.) then
               !        write(*,*) ip1,ip2,ip3,ip4,ipsum
               !        write(*,*) p1,p2,p3,p4
               !        write(*,*) iy1(ix ,jy ,mm), iy1(ixp,jy ,mm),  iy1(ix ,jyp,mm), iy1(ixp,jyp,mm)
               !      endif

            end if

            ip1 = 1.
            ip2 = 1.
            ip3 = 1.
            ip4 = 1.
            ipsum = 1.
            if (icloudtop(ix, jy, mm) .eq. icmv) then
               ip1 = 0.
               ipsum = ipsum - p1
            end if
            if (icloudtop(ixp, jy, mm) .eq. icmv) then
               ip2 = 0.
               ipsum = ipsum - p2
            end if
            if (icloudtop(ix, jyp, mm) .eq. icmv) then
               ip3 = 0.
               ipsum = ipsum - p3
            end if
            if (icloudtop(ixp, jyp, mm) .eq. icmv) then
               ip4 = 0.
               ipsum = ipsum - p4
            end if
            if (ipsum .eq. 0.) then
               yi2(m) = icmv
            else
               yi2(m) = (ip1*p1*icloudtop(ix, jy, mm) &
                  + ip2*p2*icloudtop(ixp, jy, mm) &
                  + ip3*p3*icloudtop(ix, jyp, mm) &
                  + ip4*p4*icloudtop(ixp, jyp, mm))/ipsum
            end if
            !PS end clouds
         end do
      else ! Nest
         do m = 1, 2
            mm = memind(m)
#ifdef ETA
            ytt(m) = p1*ttetan(ix, jy, kz, mm, ngrid) &
               + p2*ttetan(ixp, jy, kz, mm, ngrid) &
               + p3*ttetan(ix, jyp, kz, mm, ngrid) &
               + p4*ttetan(ixp, jyp, kz, mm, ngrid)
#else
            ytt(m) = p1*ttn(ix, jy, kz, mm, ngrid) &
               + p2*ttn(ixp, jy, kz, mm, ngrid) &
               + p3*ttn(ix, jyp, kz, mm, ngrid) &
               + p4*ttn(ixp, jyp, kz, mm, ngrid)
#endif
            ! if (lcw_nest(ngrid)) &
            ! y4(m) = p1*ctwcn(ix, jy, mm, ngrid) &
            !         + p2*ctwcn(ixp, jy, mm, ngrid) &
            !         + p3*ctwcn(ix, jyp, mm, ngrid) &
            !         + p4*ctwcn(ixp, jyp, mm, ngrid)
            ip1 = 1.
            ip2 = 1.
            ip3 = 1.
            ip4 = 1.
            ipsum = 1.
            if (icloudbotn(ix, jy, mm, ngrid) .eq. icmv) then
               ip1 = 0.
               ipsum = ipsum - p1
            end if
            if (icloudbotn(ixp, jy, mm, ngrid) .eq. icmv) then
               ip2 = 0.
               ipsum = ipsum - p2
            end if
            if (icloudbotn(ix, jyp, mm, ngrid) .eq. icmv) then
               ip3 = 0.
               ipsum = ipsum - p3
            end if
            if (icloudbotn(ixp, jyp, mm, ngrid) .eq. icmv) then
               ip4 = 0.
               ipsum = ipsum - p4
            end if
            if (ipsum .eq. 0.) then
               yi1(m) = icmv
            else
               yi1(m) = (ip1*p1*icloudbotn(ix, jy, mm, ngrid) &
                  + ip2*p2*icloudbotn(ixp, jy, mm, ngrid) &
                  + ip3*p3*icloudbotn(ix, jyp, mm, ngrid) &
                  + ip4*p4*icloudbotn(ixp, jyp, mm, ngrid))/ipsum
            end if
            !PS clouds:
            ip1 = 1.
            ip2 = 1.
            ip3 = 1.
            ip4 = 1.
            ipsum = 1.
            if (icloudtopn(ix, jy, mm, ngrid) .eq. icmv) then
               ip1 = 0.
               ipsum = ipsum - p1
            end if
            if (icloudtopn(ixp, jy, mm, ngrid) .eq. icmv) then
               ip2 = 0.
               ipsum = ipsum - p2
            end if
            if (icloudtopn(ix, jyp, mm, ngrid) .eq. icmv) then
               ip3 = 0.
               ipsum = ipsum - p3
            end if
            if (icloudtopn(ixp, jyp, mm, ngrid) .eq. icmv) then
               ip4 = 0.
               ipsum = ipsum - p4
            end if
            if (ipsum .eq. 0.) then
               yi2(m) = icmv
            else
               yi2(m) = (ip1*p1*icloudtopn(ix, jy, mm, ngrid) &
                  + ip2*p2*icloudtopn(ixp, jy, mm, ngrid) &
                  + ip3*p3*icloudtopn(ix, jyp, mm, ngrid) &
                  + ip4*p4*icloudtopn(ixp, jyp, mm, ngrid))/ipsum
            end if
            !PS end clouds
         end do
      end if

      !************************************
      ! 2.) Temporal interpolation (linear)
      !************************************

      ! Why do they use different dtt and dtp2 here compared to whats below?
      ! yint1 = (y1(1)*dtp2 + y1(2)*dtp1)/dtt ! lsp
      ! yint2 = (y2(1)*dtp2 + y2(2)*dtp1)/dtt ! cp

      yint_rh = (yrh(1)*dt2 + yrh(2)*dt1)/dt
      ytint = (ytt(1)*dt2 + ytt(2)*dt1)/dt

!PS clouds:450.

!  write(*,*) yi1(1),yi1(2),yi2(1),yi2(2),dt,dt1,dt2
      icbot = int((yi1(1)*dt2 + yi1(2)*dt1)/dt)
      if (int(yi1(1)) .eq. icmv) icbot = int(yi1(2))
      if (int(yi1(2)) .eq. icmv) icbot = int(yi1(1))

      ictop = int((yi2(1)*dt2 + yi2(2)*dt1)/dt)
      if (int(yi2(1)) .eq. icmv) ictop = int(yi2(2))
      if (int(yi2(2)) .eq. icmv) ictop = int(yi2(1))

      ! Calculate unified particle coordinate (0 = TOA, 1 = Surface)
#ifdef ETA
      eta_part = part(jpart)%zeta
#else
      eta_part = exp(-part(jpart)%z / 8400.0)
#endif

      ! Safely evaluate cloud boundaries
      if (ictop .ne. icmv .and. icbot .ne. icmv) then
#ifdef ETA
         eta_toc = real(ictop) / eta_convert
         eta_boc = real(icbot) / eta_convert
#else
         eta_toc = exp(-real(ictop) / 8400.0)
         eta_boc = exp(-real(icbot) / 8400.0)
#endif
         below_top = eta_part .gt. eta_toc
         above_bot = eta_part .lt. eta_boc
      else
         ! No cloud present, force clear-sky logic
         below_top = .false.
         above_bot = .false.
      end if

      ! Compute Top-of-Atmosphere Radiation
      jul = bdate + real(itime, kind=dp)/86400._dp
      call caldate(jul, jjjjmmdd, hhmmss)
      lon_deg = xlon0 + real(part(jpart)%xlon)*dx
      lat_deg = ylat0 + real(part(jpart)%ylat)*dy
      toa_rad = 1361.0 * max(cos(zenithangle(lat_deg, lon_deg, jjjjmmdd, hhmmss)*pi180), 0.0)

      ! Early exit for nighttime to avoid 0**0 floating point traps
      if (toa_rad .le. 1e-6) then
         yint_sr = 0.0
      else
         ! Evaluate Radiation Regimes
         if (below_top) then
            Kt = 0.75  ! Default clear-sky attenuation factor

            if (above_bot) then  ! IN CLOUD
               do m = 1, 2
                  mm = memind(m)
                  ysr(m) = p1*ssr(ix, jy, 1, mm) &
                     + p2*ssr(ixp, jy, 1, mm) &
                     + p3*ssr(ix, jyp, 1, mm) &
                     + p4*ssr(ixp, jyp, 1, mm)
               end do
               interp_ssr = max((ysr(1)*dtp2 + ysr(2)*dtp1)/dtt, 0.0)

               ! Boundary radiation calculations, reuse variables
               toa_rad = toa_rad * (Kt**eta_toc)
               interp_ssr = interp_ssr / (Kt**(1.0 - eta_boc))

               ! Fractional depth in cloud (0 = top, 1 = bottom)
               fdepth = (eta_part - eta_toc) / max(1e-6, eta_boc - eta_toc)
               fdepth = max(0.0, min(1.0, fdepth))

               ! Linear diffusive profile inside cloud
               yint_sr = toa_rad * (1.0 - fdepth) + (interp_ssr * fdepth)

            else  ! BELOW CLOUD
               do m = 1, 2
                  mm = memind(m)
                  ysr(m) = p1*ssr(ix, jy, 1, mm) &
                     + p2*ssr(ixp, jy, 1, mm) &
                     + p3*ssr(ix, jyp, 1, mm) &
                     + p4*ssr(ixp, jyp, 1, mm)
               end do
               interp_ssr = max((ysr(1)*dtp2 + ysr(2)*dtp1)/dtt, 0.0)

               ! Calculate radiation at bottom of cloud and scale to particle height
               interp_ssr = interp_ssr / (Kt**(1.0 - eta_boc))
               yint_sr = interp_ssr * (Kt**(eta_part - eta_boc))
            end if

         else  ! ABOVE CLOUD OR NO CLOUD
            if (ictop .eq. icmv) then  ! Clear sky (no cloud)
               do m = 1, 2
                  mm = memind(m)
                  ysr(m) = p1*ssr(ix, jy, 1, mm) &
                     + p2*ssr(ixp, jy, 1, mm) &
                     + p3*ssr(ix, jyp, 1, mm) &
                     + p4*ssr(ixp, jyp, 1, mm)
               end do
               interp_ssr = max((ysr(1)*dtp2 + ysr(2)*dtp1)/dtt, 0.0)
               ! Calculate clearness index and clamp between 0 and 1
               Kt = max(0.0, min(1.0, interp_ssr / toa_rad))
            else  ! Above cloud
               Kt = 0.75
            end if

            ! Scale exponentially from TOA down to particle height
            yint_sr = toa_rad * (Kt**eta_part)
         end if
      end if
      jul = bdate + real(itime, kind=dp)/86400._dp
      call caldate(jul, jjjjmmdd, hhmmss)
      !!! Uncomment for debugging !
      ! if (jpart .eq. 1) then
      !    write (*, *) '--- VIABILITY DEBUG STEP (Full) ---'
      !    write (*, *) 'Cloud status (0 = below, 1 = in, 2 = above): ', &
      !       merge(1, 0, above_bot) + merge(1, 0,.not. below_top)
      !    write (*, *) 'Time (hr): ', itime/3600.0
      !    write (* ,*) 'Date: jjjjmmdd', jjjjmmdd
      !    write (* ,*) 'Date: hhmmss', hhmmss
      !    write (*, *) 'Part Z (m): ', part(jpart)%z
      !    write (*, *) 'Part Lon (deg): ', xlon0+real(part(jpart)%xlon)*dx
      !    write (*, *) 'Part Lat (deg): ', ylat0+real(part(jpart)%ylat)*dy
      !    write (*, *) 'Part eta: ', part(jpart)%zeta
      !    write (*, *) 'Part etaupdate?: ', part(jpart)%etaupdate
      !    write (*, *) 'Part meterupdate?: ', part(jpart)%meterupdate
      !    write (*, *) 'Cloud Top (m): ', ictop, ' | Cloud Bot (m): ', icbot
      !    write (*, *) 'TOA Rad (W/m2): ', toa_rad
      !    write (*, *) 'Interp SSR (W/m2): ', interp_ssr
      !    write (*, *) 'FINAL yint_sr (W/m2): ', yint_sr
      !    write (*, *) 'FINAL yint_t (K): ', ytint
      !    write (*, *) 'FINAL yint_rh : ', yint_rh
      !    write (*, *) 'Particle Alive?: ', part(jpart)%alive
      !    write (*, *) 'Particle accumulated solar radiation (MJ): ', part(jpart)%sr_acc
      !    write (*, *) '----------------------------'
      ! end if

   end subroutine interpol_viability_full

   subroutine interpol_viability_agpads(itime, jpart, yint_sr, yint_rh, ytint, icmv)
      !                               i      i       o        o       o      i
      !****************************************************************************
      !                                                                           *
      !  Interpolation of meteorological fields for particle viability.           *
      !  In horizontal direction bilinear interpolation is used.                  *
      !  Temporally a linear interpolation is used.                               *
      !                                                                           *
      !  Special version adapted from levlininterpol / interpol_rain for the       *
      !  AgPads viability model. Interpolates solar radiation, relative           *
      !  humidity, and temperature at the particle position.                      *
      !                                                                           *
      !     Author: A. Stohl                                                      *
      !                                                                           *
      !     30 August 1996                                                        *
      !                                                                           *
      ! PS, AP 04/2019, 11/2020:                                                  *
      !                 put back temporal interpolation of rain, from v10.01      *
      !      and cloud bottom / thickness interpolation                           *
      ! PS, AP 01/2021:                                                           *
      !      interpolate particle temperature and cloud total water               *
      ! PS, AP 02/2021:                                                           *
      !      interpolation of precipitation using two additional fields           *
      !      which are temporally equidistant between the main fields             *
      !                                                                           *
      ! AJ 07/2026 - Viability interpolation updates:                                          *
      !      adapted for AgPads viability model: interpolates solar radiation,    *
      !      relative humidity, and temperature                                   *
      !                                                                           *
      !****************************************************************************
      !                                                                           *
      ! Variables:                                                                *
      !                                                                           *
      ! itime                current simulation time [s]                          *
      ! jpart                index of the particle to be interpolated             *
      ! yint_sr              interpolated solar radiation at particle [W/m2]      *
      ! yint_rh              interpolated relative humidity at particle [%]       *
      ! ytint                interpolated temperature at particle [K]             *
      ! icmv                 viability model control / component flag             *
      !                                                                           *
      ! dt, dtp1, dtp2       time differences / weights for temporal interpol.    *
      ! ip1,ip2,ip3,ip4      bilinear interpolation weights                       *
      ! kz                   vertical model level of the particle                 *
      ! memind(3)            points to the places of the wind/met fields          *
      ! xts, yts             particle coordinates in grid units                   *
      ! yrh(2)               relative humidity at bounding time steps             *
      ! ysr(2)               solar radiation at bounding time steps               *
      ! ytt(2)               temperature at bounding time steps                   *
      !                                                                           *
      !****************************************************************************
      use par_mod, only: numwfmem, numpf
      use com_mod, only: lcw
      use qvsat_mod

      implicit none

      integer, intent(in) :: itime, jpart, icmv
      real, intent(out) :: yint_rh, yint_sr, ytint
      integer :: m
      integer :: mm, kz
      real :: ip1, ip2, ip3, ip4, ipsum
      real :: dt, dtp1, dtp2, rt
      real :: yrh(2), ysr(2), ytt(2), xts, yts ! interpolated values
      integer, dimension(2) :: ip, mp

      !**********************************************************************
      ! 1.) Bilinear horizontal interpolation
      ! This has to be done separately for 2 fields (Temporal)
      !*******************************************************
      ! Loop over 2 time steps
      !***********************

      !-------------------------------------------------------------------------
      ! PS, AT new interpolation of precip with 2 additional fields
      ! therefore, we need a special treatment of lsp,cp which are in yy1,yy2
      !-------------------------------------------------------------------------
      !
      !      1.1 1.2 1.3               ip(1).mp(1)
      !      1.2 1.3 2.1               ip(2).mp(2)
      !
      !   ||___|___|___||___|___|___||
      !
      ! ip  1   2   3    1   2   3    1
      ! m        1            2
      !
      !-------------------------------------------------------------------------
#ifdef ETA
      call find_z_level_eta_uv_safe(real(part(jpart)%zeta), kz)
!      kz = induv
#else
      call find_z_level_meters_safe(real(part(jpart)%z), kz)
!      kz = indz
#endif
      xts = real(part(jpart)%xlon)
      yts = real(part(jpart)%ylat)


      if (xts .ge. real(nx - 1)) xts = real(nx - 1) - 0.00001
      if (yts .ge. real(ny - 1)) yts = real(ny - 1) - 0.00001

      call find_ngrid(xts, yts)
      call find_grid_indices(xts, yts)
      call find_grid_distances(xts, yts)
      dt1 = real(itime - memtime(1))
      dt2 = real(memtime(2) - itime)
      dt = real(memtime(2) - memtime(1))
      ! Wrap Longitude for the 360-degree line
      if (ixp .ge. nx) ixp = ixp - nx

      ! Clamp Latitude at the poles
      if (jyp .ge. ny) jyp = ny - 1

      ! Clamp Vertical to prevent top-of-atmosphere buffer read
      kz = max(1, min(kz, nz))


      dt1 = real(itime - memtime(1))
      dt2 = real(memtime(2) - itime)
      dt = real(memtime(2) - memtime(1))
      if (dt .eq. 0. .and. dt1 .eq. 0.) then ! Fix if last last timestep and memtime(2)=memtime(1)
         dt = 1.
         dt2 = 1.
      end if
      dtt = dt/3.
      if (numpf .eq. 1) then
         mp(1) = 1
         mp(2) = 2
         ip(1) = 1
         ip(2) = 1
         dtp1 = dt1
         dtp2 = dt2
      else
         rt = abs(dt1/dt)
         if (0 .le. rt .and. rt .lt. 1./3.) then
            mp(1) = 1
            mp(2) = 1
            ip(1) = 1
            ip(2) = 2
            dtp1 = dt1
            dtp2 = dt2 - 2.*dtt
         elseif (1./3. .le. rt .and. rt .lt. 2./3.) then
            mp(1) = 1
            mp(2) = 1
            ip(1) = 2
            ip(2) = 3
            dtp1 = dt1 - dtt
            dtp2 = dt2 - dtt
         elseif (2./3. .le. rt .and. rt .lt. 1.) then
            mp(1) = 1
            mp(2) = 2
            ip(1) = 3
            ip(2) = 1
            dtp1 = dt1 - 2.*dtt
            dtp2 = dt2
         end if
      end if

      if (ngrid .le. 0) then ! No nest
         do m = 1, 2
            mm = memind(m)
            yrh(m) = p1*calc_rh(ix, jy, kz, mm) &
               + p2*calc_rh(ixp, jy, kz, mm) &
               + p3*calc_rh(ix, jyp, kz, mm) &
               + p4*calc_rh(ixp, jyp, kz, mm)
            ! Temperature
         end do
      else ! Nest
         do m = 1, 2
            mm = memind(m)
#ifdef ETA
            ytt(m) = p1*ttetan(ix, jy, kz, mm, ngrid) &
               + p2*ttetan(ixp, jy, kz, mm, ngrid) &
               + p3*ttetan(ix, jyp, kz, mm, ngrid) &
               + p4*ttetan(ixp, jyp, kz, mm, ngrid)
#else
            ytt(m) = p1*ttn(ix, jy, kz, mm, ngrid) &
               + p2*ttn(ixp, jy, kz, mm, ngrid) &
               + p3*ttn(ix, jyp, kz, mm, ngrid) &
               + p4*ttn(ixp, jyp, kz, mm, ngrid)
#endif
            ! if (lcw_nest(ngrid)) &
            ! y4(m) = p1*ctwcn(ix, jy, mm, ngrid) &
            !         + p2*ctwcn(ixp, jy, mm, ngrid) &
            !         + p3*ctwcn(ix, jyp, mm, ngrid) &
            !         + p4*ctwcn(ixp, jyp, mm, ngrid)
         end do
      end if

      !************************************
      ! 2.) Temporal interpolation (linear)
      !************************************

      ! Why do they use different dtt and dtp2 here compared to whats below?
      ! yint1 = (y1(1)*dtp2 + y1(2)*dtp1)/dtt ! lsp
      ! yint2 = (y2(1)*dtp2 + y2(2)*dtp1)/dtt ! cp

      yint_rh = (yrh(1)*dt2 + yrh(2)*dt1)/dt

      do m = 1, 2
         mm = memind(m)

         ysr(m) = p1*ssr(ix, jy, 1, mm) &
            + p2*ssr(ixp, jy, 1, mm) &
            + p3*ssr(ix, jyp, 1, mm) &
            + p4*ssr(ixp, jyp, 1, mm)
      end do
      yint_sr = max((ysr(1)*dtp2 + ysr(2)*dtp1)/dtt, 0.0)
      ! yint_sr = max((ysr(1)*dt2 + ysr(2)*dt1)/dt, 0.0) * 1.25
      yint_sr = yint_sr + yint_sr*part(jpart)%z / 1000.0
      yint_sr = min(yint_sr, 1361.0)

      !!! Uncomment for debugging
      ! if (jpart .eq. 1) then
      !    write (*, *) '--- VIABILITY DEBUG STEP (AgPaDS) ---'
      !    write (*, *) 'Time (hr): ', itime/3600.0
      !    write (*, *) 'Part Z (m): ', part(jpart)%z
      !    write (*, *) 'Part eta: ', part(jpart)%zeta
      !    write (*, *) 'Part etaupdate?: ', part(jpart)%etaupdate
      !    write (*, *) 'Part meterupdate?: ', part(jpart)%meterupdate
      !    ! write (*, *) 'Interp SSR (W/m2): ', interp_ssr
      !    write (*, *) 'FINAL yint_sr (W/m2): ', yint_sr
      !    write (*, *) 'FINAL yint_rh : ', yint_rh
      !    write (*, *) 'Particle Alive?: ', part(jpart)%alive
      !    write (*, *) 'Particle accumulated solar radiation (MJ): ', part(jpart)%sr_acc
      !    write (*, *) '----------------------------'
      ! end if

   end subroutine interpol_viability_agpads

   subroutine interpol_viability_isard(itime, jpart, yint_sr, yint_rh, ytint, icmv)
      !                               i      i       o        o       o      i
      !****************************************************************************
      !                                                                           *
      !  Interpolation of meteorological fields for particle viability.           *
      !  In horizontal direction bilinear interpolation is used.                  *
      !  Temporally a linear interpolation is used.                               *
      !                                                                           *
      !  Special version adapted from levlininterpol / interpol_rain for the       *
      !  Isard viability model. Interpolates solar radiation, relative            *
      !  humidity, and temperature at the particle position.                      *
      !                                                                           *
      !     Author: A. Stohl                                                      *
      !                                                                           *
      !     30 August 1996                                                        *
      !                                                                           *
      ! PS, AP 04/2019, 11/2020:                                                  *
      !                 put back temporal interpolation of rain, from v10.01      *
      !      and cloud bottom / thickness interpolation                           *
      ! PS, AP 01/2021:                                                           *
      !      interpolate particle temperature and cloud total water               *
      ! PS, AP 02/2021:                                                           *
      !      interpolation of precipitation using two additional fields           *
      !      which are temporally equidistant between the main fields             *
      !                                                                           *
      ! Viability interpolation updates:                                          *
      !      adapted for Isard viability model: interpolates solar radiation,     *
      !      relative humidity, and temperature                                   *
      !                                                                           *
      !****************************************************************************
      !                                                                           *
      ! Variables:                                                                *
      !                                                                           *
      ! itime                current simulation time [s]                          *
      ! jpart                index of the particle to be interpolated             *
      ! yint_sr              interpolated solar radiation at particle [W/m2]      *
      ! yint_rh              interpolated relative humidity at particle [%]       *
      ! ytint                interpolated temperature at particle [K]             *
      ! icmv                 viability model control / component flag             *
      !                                                                           *
      ! dt, dtp1, dtp2       time differences / weights for temporal interpol.    *
      ! ip1,ip2,ip3,ip4      bilinear interpolation weights                       *
      ! kz                   vertical model level of the particle                 *
      ! memind(3)            points to the places of the wind/met fields          *
      ! xts, yts             particle coordinates in grid units                   *
      ! yrh(2)               relative humidity at bounding time steps             *
      ! ysr(2)               solar radiation at bounding time steps               *
      ! ytt(2)               temperature at bounding time steps                   *
      !                                                                           *
      !****************************************************************************
      use par_mod, only: numwfmem, numpf
      use com_mod, only: lcw
      use qvsat_mod

      implicit none

      integer, intent(in) :: itime, jpart, icmv
      real, intent(out) :: yint_rh, yint_sr, ytint
      integer :: m
      integer :: mm, kz
      real :: ip1, ip2, ip3, ip4, ipsum
      real :: dt, dtp1, dtp2, rt
      real :: yrh(2), ysr(2), ytt(2), xts, yts ! interpolated values
      integer, dimension(2) :: ip, mp

      !**********************************************************************
      ! 1.) Bilinear horizontal interpolation
      ! This has to be done separately for 2 fields (Temporal)
      !*******************************************************
      ! Loop over 2 time steps
      !***********************

      !-------------------------------------------------------------------------
      ! PS, AT new interpolation of precip with 2 additional fields
      ! therefore, we need a special treatment of lsp,cp which are in yy1,yy2
      !-------------------------------------------------------------------------
      !
      !      1.1 1.2 1.3               ip(1).mp(1)
      !      1.2 1.3 2.1               ip(2).mp(2)
      !
      !   ||___|___|___||___|___|___||
      !
      ! ip  1   2   3    1   2   3    1
      ! m        1            2
      !
      !-------------------------------------------------------------------------


      !xts = real(part(jpart)%xlon)
      !yts = real(part(jpart)%ylat)

      !if (xts .ge. real(nx - 1)) xts = real(nx - 1) - 0.00001
      !if (yts .ge. real(ny - 1)) yts = real(ny - 1) - 0.00001
      !call init_interpol(itime, xts, yts, real(part(jpart)%z), real(part(jpart)%zeta))
#ifdef ETA
      call find_z_level_eta_uv_safe(real(part(jpart)%zeta), kz)
!      kz = induv
#else
      call find_z_level_meters_safe(real(part(jpart)%z), kz)
!      kz = indz
#endif
      xts = real(part(jpart)%xlon)
      yts = real(part(jpart)%ylat)


      if (xts .ge. real(nx - 1)) xts = real(nx - 1) - 0.00001
      if (yts .ge. real(ny - 1)) yts = real(ny - 1) - 0.00001

      call find_ngrid(xts, yts)
      call find_grid_indices(xts, yts)
      call find_grid_distances(xts, yts)
      dt1 = real(itime - memtime(1))
      dt2 = real(memtime(2) - itime)
      dt = real(memtime(2) - memtime(1))
      ! Wrap Longitude for the 360-degree line
      if (ixp .ge. nx) ixp = ixp - nx

      ! Clamp Latitude at the poles
      if (jyp .ge. ny) jyp = ny - 1

      ! Clamp Vertical to prevent top-of-atmosphere buffer read
      kz = max(1, min(kz, nz))


      if (dt .eq. 0. .and. dt1 .eq. 0.) then ! Fix if last last timestep and memtime(2)=memtime(1)
         dt = 1.
         dt2 = 1.
      end if
      dtt = dt/3.
      if (numpf .eq. 1) then
         mp(1) = 1
         mp(2) = 2
         ip(1) = 1
         ip(2) = 1
         dtp1 = dt1
         dtp2 = dt2
      else
         rt = abs(dt1/dt)
         if (0 .le. rt .and. rt .lt. 1./3.) then
            mp(1) = 1
            mp(2) = 1
            ip(1) = 1
            ip(2) = 2
            dtp1 = dt1
            dtp2 = dt2 - 2.*dtt
         elseif (1./3. .le. rt .and. rt .lt. 2./3.) then
            mp(1) = 1
            mp(2) = 1
            ip(1) = 2
            ip(2) = 3
            dtp1 = dt1 - dtt
            dtp2 = dt2 - dtt
         elseif (2./3. .le. rt .and. rt .lt. 1.) then
            mp(1) = 1
            mp(2) = 2
            ip(1) = 3
            ip(2) = 1
            dtp1 = dt1 - 2.*dtt
            dtp2 = dt2
         end if
      end if

      if (ngrid .le. 0) then ! No nest
         do m = 1, 2
            mm = memind(m)
            yrh(m) = p1*calc_rh(ix, jy, kz, mm) &
               + p2*calc_rh(ixp, jy, kz, mm) &
               + p3*calc_rh(ix, jyp, kz, mm) &
               + p4*calc_rh(ixp, jyp, kz, mm)
            ! Temperature
         end do
      else ! Nest
         do m = 1, 2
            mm = memind(m)
#ifdef ETA
            ytt(m) = p1*ttetan(ix, jy, kz, mm, ngrid) &
               + p2*ttetan(ixp, jy, kz, mm, ngrid) &
               + p3*ttetan(ix, jyp, kz, mm, ngrid) &
               + p4*ttetan(ixp, jyp, kz, mm, ngrid)
#else
            ytt(m) = p1*ttn(ix, jy, kz, mm, ngrid) &
               + p2*ttn(ixp, jy, kz, mm, ngrid) &
               + p3*ttn(ix, jyp, kz, mm, ngrid) &
               + p4*ttn(ixp, jyp, kz, mm, ngrid)
#endif
            ! if (lcw_nest(ngrid)) &
            ! y4(m) = p1*ctwcn(ix, jy, mm, ngrid) &
            !         + p2*ctwcn(ixp, jy, mm, ngrid) &
            !         + p3*ctwcn(ix, jyp, mm, ngrid) &
            !         + p4*ctwcn(ixp, jyp, mm, ngrid)
         end do
      end if

      !************************************
      ! 2.) Temporal interpolation (linear)
      !************************************

      ! Why do they use different dtt and dtp2 here compared to whats below?
      ! yint1 = (y1(1)*dtp2 + y1(2)*dtp1)/dtt ! lsp
      ! yint2 = (y2(1)*dtp2 + y2(2)*dtp1)/dtt ! cp

      yint_rh = (yrh(1)*dt2 + yrh(2)*dt1)/dt

      do m = 1, 2
         mm = memind(m)

         ysr(m) = p1*ssr(ix, jy, 1, mm) &
            + p2*ssr(ixp, jy, 1, mm) &
            + p3*ssr(ix, jyp, 1, mm) &
            + p4*ssr(ixp, jyp, 1, mm)
      end do
      yint_sr = max((ysr(1)*dtp2 + ysr(2)*dtp1)/dtt, 0.0)
      ! yint_sr = max((ysr(1)*dt2 + ysr(2)*dt1)/dt, 0.0) * 1.25

      !!! Uncomment for debugging
      ! if (jpart .eq. 1) then
      !    write (*, *) '--- VIABILITY DEBUG STEP (SSR Method) ---'
      !    write (*, *) 'Time (hr): ', itime/3600.0
      !    write (*, *) 'Part Z (m): ', part(jpart)%z
      !    write (*, *) 'Part eta: ', part(jpart)%zeta
      !    write (*, *) 'Part etaupdate?: ', part(jpart)%etaupdate
      !    write (*, *) 'Part meterupdate?: ', part(jpart)%meterupdate
      !    ! write (*, *) 'Interp SSR (W/m2): ', interp_ssr
      !    write (*, *) 'FINAL yint_sr (W/m2): ', yint_sr
      !    write (*, *) 'FINAL yint_rh : ', yint_rh
      !    write (*, *) 'Particle Alive?: ', part(jpart)%alive
      !    write (*, *) 'Particle accumulated solar radiation (MJ): ', part(jpart)%sr_acc
      !    write (*, *) '----------------------------'
      ! end if

   end subroutine interpol_viability_isard


   function calc_rh(ix, jy, kz, n) result(rh)
      implicit none
      integer, intent(in) :: ix, jy, kz, n
      real :: rh, pressure

      pressure = rho(ix, jy, kz, n)*r_air*tt(ix, jy, kz, n)
      rh = min(1.0, qv(ix, jy, kz, n)/min(f_qvsat(pressure, tt(ix, jy, kz, n)), 1e-8))

   end function calc_rh

   function calc_surf_rh(ix, jy, n) result(rh)
      implicit none
      integer, intent(in) :: ix, jy, n
      real :: rh, pressure
      rh = ew(td2(ix, jy, 1, n), ps(ix, jy, 1, n))/ew(tt2(ix, jy, 1, n), ps(ix, jy, 1, n))

   end function calc_surf_rh

   function zenithangle(ylat, xlon, jjjjmmdd, hhmmss) result(sza)

      !*****************************************************************************
      !                                                                            *
      !    This function returns the sinus of solar elevation as a function        *
      !    of geographic longitde, latitude and GMT-Time.                          *
      !                                                                            *
      !    Author: G. Wotawa, Nov 1993                                             *
      !                                                                            *
      !*****************************************************************************
      !                                                                            *
      !     INPUT:                                                                 *
      !                                                                            *
      !            ylat          geographical latitude  [DEG]                      *
      !            xlon          geographical longitude [DEG]                      *
      !            jjjj          Year                                              *
      !            mm            Month                                             *
      !            dd            Day                                               *
      !            hh            Hour                                              *
      !            minute        Minute                                            *
      !                                                                            *
      !*****************************************************************************

      implicit none

      integer, intent(in)  :: jjjjmmdd, hhmmss
      real, intent(in)  :: ylat, xlon
      real    :: sza
      integer :: jjjj, mm, id, iu, minute
      integer :: ndaynum
      real :: sinsol, solelev
      real :: rnum, rylat, ttime, dekl, rdekl, eq
      real, parameter :: pi = 3.1415927

      jjjj = jjjjmmdd/10000
      mm = jjjjmmdd/100 - jjjj*100
      id = jjjjmmdd - jjjj*10000 - mm*100
      iu = hhmmss/10000
      minute = hhmmss/100 - 100*iu

      ndaynum = 31*(mm - 1) + id
      if (mm .gt. 2) ndaynum = ndaynum - int(0.4*mm + 2.3)
      if ((mm .gt. 2) .and. (jjjj/4*4 .eq. jjjj)) ndaynum = ndaynum + 1

      rnum = 2.*pi*ndaynum/365.
      rylat = pi*ylat/180.
      ttime = real(iu) + real(minute)/60.

      dekl = 0.396 + 3.631*sin(rnum) + 0.038*sin(2.*rnum) + 0.077*sin(3.*rnum) - &
         22.97*cos(rnum) - 0.389*cos(2.*rnum) - 0.158*cos(3.*rnum)
      rdekl = pi*dekl/180.

      eq = (0.003 - 7.343*sin(rnum) - 9.47*sin(2.*rnum) - &
         0.329*sin(3.*rnum) - 0.196*sin(4.*rnum) + &
         0.552*cos(rnum) - 3.020*cos(2.*rnum) - &
         0.076*cos(3.*rnum) - 0.125*cos(4.*rnum))/60.
      sinsol = sin(rylat)*sin(rdekl) + cos(rylat)*cos(rdekl)* &
         cos((ttime - 12.+xlon/15.+eq)*pi/12.)
      ! Calculate the maximum solar elevation on that day
      !sinsol=sin(rylat)*sin(rdekl)+cos(rylat)*cos(rdekl)*
      !    &       cos((eq)*pi/12.)
      solelev = asin(sinsol)*180./pi
      sza = 90.-solelev

   end function zenithangle

   subroutine calculate_logistic_decay_sr(jpart, sr_interp, rh_interp, decay_rate)

      ! Fits a decaying logistic function to calculate viable fraction due to solar radiation.
      ! returns the rate of decay at the current time instance
      implicit none
      integer, intent(in) :: jpart
      real, intent(in) :: sr_interp, rh_interp
      real, intent(out) :: decay_rate
      real :: V_sr
      ! real, parameter :: k_sr = 0.55, sr_half = 12.0

      V_sr = 1.0/(1.0 + exp(k_sr*(part(jpart)%sr_acc) - sr_half))
      decay_rate = k_sr*(1.0 - V_sr)*(sr_interp*1.0E-6)

   end subroutine calculate_logistic_decay_sr

   subroutine calculate_logistic_decay_sr_rh(jpart, sr_interp, rh_interp, decay_rate)

      ! Fits a decaying logistic function to calculate viable fraction due to solar radiation.
      ! returns the rate of decay at the current time instance
      implicit none
      integer, intent(in) :: jpart
      real, intent(in) :: sr_interp, rh_interp
      real, intent(out) :: decay_rate
      real :: V_sr, s
      ! real, parameter :: k_sr = 0.55, sr_half = 12.0
      if (rh_interp >= 0.9) then
         s = 2.0
      else
         s = 1.0
      end if

      V_sr = 1.0/(1.0 + exp(k_sr*(part(jpart)%sr_acc) - sr_half))
      decay_rate = k_sr*(1.0 - V_sr)*(s*sr_interp*1.0E-6)

   end subroutine calculate_logistic_decay_sr_rh

   subroutine calculate_decay_isard(jpart, sr_interp, rh_interp, decay_rate)

      ! Fits a decaying logistic function to calculate viable fraction due to solar radiation
      ! based on Isard 2005 paper. Returns the rate of decay at the current time instance.
      implicit none
      integer, intent(in) :: jpart
      real, intent(in) :: sr_interp, rh_interp
      real, intent(out) :: decay_rate
      ! real, parameter :: k_isard = 14.0

      decay_rate = (sr_interp*1.0E-6) / k_isard

   end subroutine calculate_decay_isard

   subroutine calculate_decay_uv_rh_agpads(jpart, sr_interp, rh_interp, decay_rate)

      ! Fits a decaying logistic function to calculate viable fraction due to solar radiation.
      ! returns the rate of decay at the current time instance
      implicit none
      integer, intent(in) :: jpart
      real, intent(in) :: sr_interp, rh_interp
      real, intent(out) :: decay_rate
      real :: s
      ! real, parameter :: k = 14.0
      if (rh_interp >= 0.9) then
         s = 2.0
      else
         s = 1.0
      end if

      decay_rate = s*(sr_interp*1.0E-6) / k_isard
   end subroutine calculate_decay_uv_rh_agpads


   subroutine calculate_decay_two_stage(jpart, sr_interp, rh_interp, decay_rate)

      ! Fits a decaying logistic function to calculate viable fraction due to solar radiation.
      ! returns the rate of decay at the current time instance
      implicit none
      integer, intent(in) :: jpart
      real, intent(in) :: sr_interp, rh_interp
      real, intent(out) :: decay_rate
      real :: f_fast, f_slow, D_sr, denom, lambda_dose
      ! real, parameter :: strong_frac = 0.05
      real :: k1_rec, k2_rec
      k1_rec = 1.0 / k_isard
      k2_rec = 1.0 / k2_slow

      D_sr = part(jpart)%sr_acc

      f_fast = (1.0 - strong_frac) * exp(-D_sr * k1_rec)
      f_slow = strong_frac * exp(-D_sr * k2_rec)
      denom = f_fast + f_slow
      if (denom < 1e-10) then
         lambda_dose = k2_rec
      else
         lambda_dose = (k1_rec * f_fast + k2_rec * f_slow) / denom
      end if

      decay_rate = lambda_dose * (sr_interp * 1.0E-6)

   end subroutine calculate_decay_two_stage

   subroutine calculate_decay_two_stage_shoulder(jpart, sr_interp, rh_interp, decay_rate)

      ! Fits a decaying logistic function to calculate viable fraction due to solar radiation.
      ! returns the rate of decay at the current time instance
      implicit none
      integer, intent(in) :: jpart
      real, intent(in) :: sr_interp, rh_interp
      real, intent(out) :: decay_rate
      real :: exp_fast, exp_slow, D_sr, denom, lambda_dose, shoulder_term
      ! real, parameter :: strong_frac = 0.05
      ! integer, parameter :: n1_s = 4
      real :: k1_rec, k2_rec
      k1_rec = 1.0 / k_isard
      k2_rec = 1.0 / k2_slow

      D_sr = part(jpart)%sr_acc

      exp_fast = exp(-D_sr * k1_rec)
      exp_slow = exp(-D_sr * k2_rec)
      shoulder_term = max(1.0 - exp_fast, 0.0)
      denom = (1.0 - strong_frac) * (1.0 - shoulder_term**n1_s) + strong_frac * exp_slow
      if (denom < 1e-10) then
         lambda_dose = k2_rec
      else
         lambda_dose = ((1.0-strong_frac)*real(n1_s)*k1_rec*exp_fast*shoulder_term**(n1_s - 1)+strong_frac*k2_rec*exp_slow) / denom
      end if

      decay_rate = lambda_dose * (sr_interp * 1.0E-6)

   end subroutine calculate_decay_two_stage_shoulder

   subroutine calculate_decay_constant_exponential(jpart, sr_interp, rh_interp, decay_rate)

      ! Fits a decaying logistic function to calculate viable fraction due to solar radiation.
      ! returns the rate of decay at the current time instance
      implicit none
      integer, intent(in) :: jpart
      real, intent(in) :: sr_interp, rh_interp
      real, intent(out) :: decay_rate
      real, parameter :: scale_decay = 1.0 / (1000.0 * 3600.0)

      ! since this version does not use interp_via_method, we use interp_via_method to store the
      ! decay_rate, in per/hr but since interp_via_method is an integer, we use 1000 * decay_rate
      ! as interp_via_method
      decay_rate = real(interp_via_method) * scale_decay
   end subroutine calculate_decay_constant_exponential

#endif

!*********************
!* PRIVATE FUNCTIONS *
!*********************
! Interpolation of wind fields
!*****************************
#ifdef ETA
   subroutine interpol_wind_eta(zteta, iuv, iweta)

!* PRIVATE FUNCTION *

      real, intent(in)    :: zteta
      integer, intent(in)  :: iuv(2), iweta(2)
      integer             :: n, m
      real                :: uh(2), vh(2), wetah(2), uh1(2), vh1(2), wetah1(2)
      real                :: dz1uv, dz2uv, dz1weta, dz2weta

      !**********************************************************************
      ! 1.) Bilinear horizontal interpolation
      ! This has to be done separately for 6 fields (Temporal(2)*Vertical(3))
      !**********************************************************************

      ! Vertical distance to the level below and above current position
      !****************************************************************
      call find_vert_vars(uvheight, zteta, induv, dz1uv, dz2uv, lbounds_uv, .false.)
      call find_vert_vars(wheight, zteta, indzeta, dz1weta, dz2weta, lbounds_w, .true.)

      ! Loop over 2 time steps and 2 levels
      !************************************
      if (ngrid .le. 0) then ! No nest
         do m = 1, 2
            do n = 1, 2
               call hor_interpol(wweta, wetah1(n), iweta(n), memind(m), nzmax)
               if (ngrid .lt. 0) then
                  call hor_interpol(uupoleta, uh1(n), iuv(n), memind(m), nzmax)
                  call hor_interpol(vvpoleta, vh1(n), iuv(n), memind(m), nzmax)
               else
                  call hor_interpol(uueta, uh1(n), iuv(n), memind(m), nzmax)
                  call hor_interpol(vveta, vh1(n), iuv(n), memind(m), nzmax)
               end if
            end do
            call vert_interpol(uh1(1), uh1(2), dz1uv, dz2uv, uh(m))
            call vert_interpol(vh1(1), vh1(2), dz1uv, dz2uv, vh(m))
            call vert_interpol(wetah1(1), wetah1(2), dz1weta, dz2weta, wetah(m))
         end do
      else ! Nest
         do m = 1, 2
            do n = 1, 2

               ! wetah1(n) = p1*wwetan(ix ,jy ,iweta(n),memind(m),ngrid) &
               !           + p2*wwetan(ixp,jy ,iweta(n),memind(m),ngrid) &
               !           + p3*wwetan(ix ,jyp,iweta(n),memind(m),ngrid) &
               !           + p4*wwetan(ixp,jyp,iweta(n),memind(m),ngrid)
               call hor_interpol_nest(wwetan, wetah1(n), iweta(n), memind(m), nzmax)
               call hor_interpol_nest(uuetan, uh1(n), iuv(n), memind(m), nzmax)
               call hor_interpol_nest(vvetan, vh1(n), iuv(n), memind(m), nzmax)
            end do
            call vert_interpol(uh1(1), uh1(2), dz1uv, dz2uv, uh(m))
            call vert_interpol(vh1(1), vh1(2), dz1uv, dz2uv, vh(m))
            call vert_interpol(wetah1(1), wetah1(2), dz1weta, dz2weta, wetah(m))
         end do
      end if
      call temporal_interpolation(uh(1), uh(2), u)
      call temporal_interpolation(vh(1), vh(2), v)
      call temporal_interpolation(wetah(1), wetah(2), weta)
   end subroutine interpol_wind_eta
#else

   subroutine interpol_wind_meter(zt, iw)

!* PRIVATE FUNCTION *

      real, intent(in)    :: zt
      integer, intent(in)  :: iw(2)
      integer             :: n, m
      real                :: uh(2), vh(2), wh(2), uh1(2), vh1(2), wh1(2)
      real                :: dz1w, dz2w

      !**********************************************************************
      ! 1.) Bilinear horizontal interpolation
      ! This has to be done separately for 6 fields (Temporal(2)*Vertical(3))
      !**********************************************************************

      ! Vertical distance to the level below and above current position
      !****************************************************************
      call find_vert_vars(height, zt, indz, dz1w, dz2w, lbounds, .false.)

      ! Loop over 2 time steps and 2 levels
      !************************************
      if (ngrid .le. 0) then ! No nest
         do m = 1, 2
            do n = 1, 2
               call hor_interpol(ww, wh1(n), iw(n), memind(m), nzmax)
               if (ngrid .lt. 0) then
                  call hor_interpol(uupol, uh1(n), iw(n), memind(m), nzmax)
                  call hor_interpol(vvpol, vh1(n), iw(n), memind(m), nzmax)
               else
                  call hor_interpol(uu, uh1(n), iw(n), memind(m), nzmax)
                  call hor_interpol(vv, vh1(n), iw(n), memind(m), nzmax)
               end if
            end do
            call vert_interpol(wh1(1), wh1(2), dz1w, dz2w, wh(m))
            call vert_interpol(uh1(1), uh1(2), dz1w, dz2w, uh(m))
            call vert_interpol(vh1(1), vh1(2), dz1w, dz2w, vh(m))
         end do
      else ! Nest
         do m = 1, 2
            do n = 1, 2
               call hor_interpol_nest(wwn, wh1(n), iw(n), memind(m), nzmax)
               call hor_interpol_nest(uun, uh1(n), iw(n), memind(m), nzmax)
               call hor_interpol_nest(vvn, vh1(n), iw(n), memind(m), nzmax)
            end do
            call vert_interpol(wh1(1), wh1(2), dz1w, dz2w, wh(m))
            call vert_interpol(uh1(1), uh1(2), dz1w, dz2w, uh(m))
            call vert_interpol(vh1(1), vh1(2), dz1w, dz2w, vh(m))
         end do
      end if
      call temporal_interpolation(wh(1), wh(2), w)
      call temporal_interpolation(uh(1), uh(2), u)
      call temporal_interpolation(vh(1), vh(2), v)
   end subroutine interpol_wind_meter
#endif

#ifdef ETA
   subroutine interpol_partoutput_val_eta(fieldname, output, j)
      implicit none

      integer, intent(in)         :: j          ! particle number
      character(2), intent(in)    :: fieldname  ! input field to interpolate over
      real, intent(inout)         :: output
      real                        :: field1(2)

      if (int(dz1out) .eq. -1) then
         call find_z_level_eta(real(part(j)%zeta))
         call find_vert_vars(uvheight, real(part(j)%zeta), induv, dz1out, dz2out, &
            lbounds_uv, .false.)
      end if

      select case (fieldname)
       case ('PR', 'pr')
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(prseta, field1, induv, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(prsetan, field1, induv, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
       case ('PV', 'pv')
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(pveta, field1, induv, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(pvetan, field1, induv, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
       case ('QV', 'qv')
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(qv, field1, induv, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(qvn, field1, induv, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
       case ('TT', 'tt')
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(tteta, field1, induv, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(ttetan, field1, induv, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
       case ('UU', 'uu')
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(uueta, field1, induv, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(uuetan, field1, induv, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
       case ('VV', 'vv')
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(vveta, field1, induv, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(vvetan, field1, induv, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
       case ('WW', 'ww')
         call find_z_level_meters(real(part(j)%z))
         call find_vert_vars(height, real(part(j)%z), indz, dz1out, dz2out, lbounds, .false.)
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(ww, field1, induv, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(wwn, field1, induv, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
         dz1out = -1
       case ('RH', 'rh')
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(rhoeta, field1, induv, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(rhoetan, field1, induv, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
      end select
   end subroutine interpol_partoutput_val_eta
#else

   subroutine interpol_partoutput_val_meter(fieldname, output, j)
      implicit none

!* PRIVATE FUNCTION *

      integer, intent(in)         :: j          ! particle number
      character(2), intent(in)    :: fieldname  ! input field to interpolate over
      real, intent(inout)         :: output
      real                        :: field1(2)

      if (int(dz1out) .eq. -1) then
         call find_z_level_meters(real(part(j)%z))
         call find_vert_vars(height, real(part(j)%z), indz, dz1out, dz2out, lbounds, .false.)
      end if

      select case (fieldname)
       case ('PR', 'pr')
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(prs, field1, indz, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(prsn, field1, indz, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
       case ('PV', 'pv')
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(pv, field1, indz, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(pvn, field1, indz, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
       case ('QV', 'qv')
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(qv, field1, indz, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(qvn, field1, indz, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
       case ('TT', 'tt')
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(tt, field1, indz, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(ttn, field1, indz, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
       case ('UU', 'uu')
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(uu, field1, indz, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(uun, field1, indz, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
       case ('VV', 'vv')
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(vv, field1, indz, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(vvn, field1, indz, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
       case ('WW', 'ww')
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(ww, field1, indz, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(wwn, field1, indz, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
       case ('RH', 'rh')
         if (ngrid .le. 0) then
            call bilin_spatial_interpol(rho, field1, indz, dz1out, dz2out, nzmax)
         else
            call bilin_spatial_interpol_nest(rhon, field1, indz, dz1out, dz2out, nzmax)
         end if
         call temporal_interpolation(field1(1), field1(2), output)
      end select
   end subroutine interpol_partoutput_val_meter
#endif
! #ifdef ETA
! subroutine interpol_pbl_eta(zt,zteta,rhoa,rhograd,ithread)

!   integer,intent(in)  :: ithread
!   real, intent(in)    :: zt,zteta
!   real, intent(inout) :: rhoa,rhograd
!   real                :: dz1w,dz2w,dz1uv,dz2uv,dz1weta,dz2weta

!   call find_vert_vars(height,zt,indz,dz1w,dz2w,lbounds,.false.)
!   call find_vert_vars(uvheight,zteta,induv,dz1uv,dz2uv,lbounds_uv,.false.)
!   call find_vert_vars(wheight,zteta,indzeta,dz1weta,dz2weta,lbounds_w,.true.)

!   call vert_interpol(wprof(indz,ithread),wprof(indzp,ithread),dz1w,dz2w,w)
!   call vert_interpol(uprof(induv,ithread),uprof(indpuv,ithread),dz1uv,dz2uv,u)
!   call vert_interpol(vprof(induv,ithread),vprof(indpuv,ithread),dz1uv,dz2uv,v)
!   call vert_interpol(rhoprof(induv,ithread),rhoprof(indpuv,ithread),dz1uv,dz2uv,rhoa)
!   call vert_interpol(rhogradprof(induv,ithread),rhogradprof(indpuv,ithread),dz1uv,dz2uv,rhograd)
!   call vert_interpol(wprofeta(indzeta,ithread),wprofeta(indzpeta,ithread),dz1weta,dz2weta,weta)
! end subroutine interpol_pbl_eta
! #endif

#ifdef ETA
   subroutine stdev_eta(iw, iuv, iweta)

!* PRIVATE FUNCTION *

      ! Standard deviation of surrounding grid points
      ! Only used in mesoscale turbulence calculations
      !***********************************************

      integer, intent(in)  :: iw(2), iuv(2), iweta(2)
      real :: wsl, wsq, usl, usq, vsl, vsq, wetasl, wetasq
      integer             :: n, m
      real, parameter      :: eps = 1.0e-30

      ! Standard deviations
      !********************
      wsl = 0.
      wsq = 0.
      usl = 0.
      usq = 0.
      vsl = 0.
      vsq = 0.
      wetasl = 0.
      wetasq = 0.

      if (ngrid .le. 0) then ! No nest
         do m = 1, 2
            do n = 1, 2
               call compute_sl_sq(ww, wsl, wsq, iw(n), memind(m), nzmax)
               call compute_sl_sq(wweta, wetasl, wetasq, iweta(n), memind(m), nzmax)
               if (ngrid .lt. 0) then
                  call compute_sl_sq(uupoleta, usl, usq, iuv(n), memind(m), nzmax)
                  call compute_sl_sq(vvpoleta, vsl, vsq, iuv(n), memind(m), nzmax)
               else
                  call compute_sl_sq(uueta, usl, usq, iuv(n), memind(m), nzmax)
                  call compute_sl_sq(vveta, vsl, vsq, iuv(n), memind(m), nzmax)
               end if
            end do
         end do
      else ! Nest
         do m = 1, 2
            do n = 1, 2
               call compute_sl_sq_nest(wwn, wsl, wsq, iw(n), memind(m), nzmax)
               call compute_sl_sq_nest(wwetan, wetasl, wetasq, iweta(n), memind(m), nzmax)
               call compute_sl_sq_nest(uuetan, usl, usq, iuv(n), memind(m), nzmax)
               call compute_sl_sq_nest(vvetan, vsl, vsq, iuv(n), memind(m), nzmax)
            end do
         end do
      end if

      call stdev(wsl, wsq, 16., wsig)
      call stdev(usl, usq, 16., usig)
      call stdev(vsl, vsq, 16., vsig)
      call stdev(wetasl, wetasq, 16., wsigeta)

   end subroutine stdev_eta
#else

   subroutine stdev_meter(iw)

!* PRIVATE FUNCTION *

      ! Standard deviation of surrounding grid points
      ! Only used in mesoscale turbulence calculations
      !***********************************************

      integer, intent(in)  :: iw(2)
      real                :: wsl, wsq, usl, usq, vsl, vsq
      integer             :: n, m
      real, parameter      :: eps = 1.0e-30

      ! Standard deviations
      !********************
      wsl = 0.
      wsq = 0.
      usl = 0.
      usq = 0.
      vsl = 0.
      vsq = 0.

      if (ngrid .le. 0) then ! No nest
         do m = 1, 2
            do n = 1, 2
               call compute_sl_sq(ww, wsl, wsq, iw(n), memind(m), nzmax)
               if (ngrid .lt. 0) then
                  call compute_sl_sq(uupol, usl, usq, iw(n), memind(m), nzmax)
                  call compute_sl_sq(vvpol, vsl, vsq, iw(n), memind(m), nzmax)
               else
                  call compute_sl_sq(uu, usl, usq, iw(n), memind(m), nzmax)
                  call compute_sl_sq(vv, vsl, vsq, iw(n), memind(m), nzmax)
               end if
            end do
         end do
      else ! Nest
         do m = 1, 2
            do n = 1, 2
               call compute_sl_sq_nest(wwn, wsl, wsq, iw(n), memind(m), nzmax)
               call compute_sl_sq_nest(uun, usl, usq, iw(n), memind(m), nzmax)
               call compute_sl_sq_nest(vvn, vsl, vsq, iw(n), memind(m), nzmax)
            end do
         end do
      end if

      call stdev(wsl, wsq, 16., wsig)
      call stdev(usl, usq, 16., usig)
      call stdev(vsl, vsq, 16., vsig)

   end subroutine stdev_meter
#endif
end module interpol_mod
