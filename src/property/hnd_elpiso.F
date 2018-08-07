C     This routine calculates the electrostatic potential
C     on a given electronic density isosurface.
C
C     Mark J. Williamson Spring 2016
C     mw529@cam.ac.uk
C
      subroutine hnd_elpiso(rtdb,basis,geom)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "bas.fh"
      integer rtdb           ! [Input] rtdb
      integer basis          ! [Input] Basis set
      integer geom           ! [Input] Geometry

      double precision iso   ! Electron density iso surface to work on
      double precision tol   ! Tolerance of the isovalue

      integer nat            ! Number of atoms in system
      character*16 at_tag    ! Geom centre tag
      integer g_dens(3)      ! Handle to density
      integer ndens          ! Number of active density handles (RHF=1, UHF=3)
      integer nclosed(2)     ! Occupation info
      integer nopen(2)       ! Occupation info
      integer nvirt(2)       ! Occupation info
      character*3 scftyp     ! Type of wave function
      character*30 theory    ! Input theory

      integer h_xyzpt        ! Memhandle to atoms coordinate array
      integer i_xyzpt        ! Index to atoms coordinate array
      integer h_zanpt        ! Memhandle to atoms charge array
      integer i_zanpt        ! Index to atoms charge array

      integer sub_nprp       ! No. of points on the electron density isosurface
      integer sub_h_prp_c    ! Memhandle points on the electron density iso
      integer sub_i_prp_c    ! Index to points on the electron density cube

      integer h_epot         ! Memhandle to electronic contribution of the electronic integral
      integer i_epot         ! Index to electronic contribution of the electronic integral

      integer h_tepot        ! Memhandle to total electric field array
      integer i_tepot        ! Index to total electric field array

      integer iat            ! Atom counter variable
      integer i              ! General counter variable
      double precision xp    ! Cartesian component of point on electron density isosurface
      double precision yp    ! Cartesian component of point on electron density isosurface
      double precision zp    ! Cartesian component of point on electron density isosurface
      double precision xn    ! Cartesian component of atom 
      double precision yn    ! Cartesian component of atom 
      double precision zn    ! Cartesian component of atom 

      double precision zan   ! Charge on an atom
      double precision elpotn! Electrostatic potential
      double precision rr    ! Distance

      character*255 esp_name ! File name for unformatted output cube

C set inputs, if none, fall back to a default
      if(.not.rtdb_get(rtdb,'prop:espiso:iso',MT_DBL,1,iso)) then
         iso = 0.002D0
      end if
      if(.not.rtdb_get(rtdb,'prop:espiso:tol',MT_DBL,1,tol)) then
         tol = 0.00003D0
      end if

C     Initialize integrals
C
      call int_init(rtdb,1, basis)
      call schwarz_init(geom, basis)
C
C     Get density matrix
C
      if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
     + call errquit('task: no task input for theory?',0, RTDB_ERR)

      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)

C     Determine how many atoms we have
      if (.not.geom_ncent(geom,nat))
     &    call errquit('hnd_elpiso: geom_ncent',911,GEOM_ERR)

C     Allocate arrays for atomic positions and charges
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',h_xyzpt,i_xyzpt))
     &    call errquit('hnd_elpiso: h_xyzpt ma failed',911,MA_ERR)

      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',h_zanpt,i_zanpt))
     &    call errquit('hnd_elpiso: h_zanpt ma failed',911,MA_ERR)

C     Populate arrays with atom locations and their charge
      do 30 iat=1,nat
        if(.not.geom_cent_get(geom,iat,at_tag,dbl_mb(i_xyzpt+3*(iat-1)),
     &     dbl_mb(i_zanpt+iat-1))) call
     &     errquit('hnd_elpiso: geom_cent_get',911,GEOM_ERR)
   30 continue

C     Calculate the electron density and extract a subset of points
      call get_points_from_iso(rtdb,basis,geom,iso,tol,
     & sub_nprp,sub_i_prp_c)

C     Allocate total electric field array
      if (.not. ma_push_get(mt_dbl,sub_nprp,'tot epot',h_tepot,i_tepot))
     &    call errquit('hnd_elpiso: h_tepot ma failed',911,MA_ERR)

C     
      if (.not. ma_push_get(mt_dbl,sub_nprp,'epot pnt',h_epot,i_epot))
     &    call errquit('hnd_elpiso: h_epot ma failed',911,MA_ERR)


C     This routine calculates the electronic contribution of the
C     electronic integral defined by nder for a given density at
C     the grid points defined in points.
C     nder=0; compute the electrostatic potential
      call hnd_elfcon(basis,geom,g_dens(ndens),
     &       dbl_mb(sub_i_prp_c),sub_nprp,dbl_mb(i_epot),0)


C Main loop
      if (ga_nodeid().gt.0) goto 300
C
C     ----- calculate electrostatic potential -----
C
      do 230  iat=1,sub_nprp
         xp = dbl_mb(sub_i_prp_c  +3*(iat-1))
         yp = dbl_mb(sub_i_prp_c+1+3*(iat-1))
         zp = dbl_mb(sub_i_prp_c+2+3*(iat-1))

         elpotn = -dbl_mb(i_epot+iat-1)
C
C     ----- add nuclear contribution -----
C
         do 210 i = 1,nat
            xn  = dbl_mb(i_xyzpt  +3*(i-1)) - xp
            yn  = dbl_mb(i_xyzpt+1+3*(i-1)) - yp
            zn  = dbl_mb(i_xyzpt+2+3*(i-1)) - zp
            zan = dbl_mb(i_zanpt+i-1)

            rr =  sqrt(xn*xn + yn*yn + zn*zn)
            if (rr.lt.1.0d-3) go to 210
            elpotn = elpotn + zan/rr
  210    continue

         dbl_mb(i_tepot+iat-1) = elpotn

  230    continue ! Assembling and printing next atom



  300 call ga_sync()

      if (.not. rtdb_cget(rtdb, "prop:grid:output",1,esp_name))
     >    call util_file_prefix("cube",esp_name)

      if(ga_nodeid().eq.0) then

C     Write to an unformatted cube file
        call prop_grid_write_cube_unformatted(geom,sub_nprp,
     >      dbl_mb(i_tepot), sub_i_prp_c, esp_name)
      end if


C
C Clean up phase
C
      if (.not.ma_pop_stack(h_epot)) call errquit
     &   ('hnd_elpiso, ma_pop_stack of h_epot failed',911,MA_ERR)

      if (.not.ma_pop_stack(h_tepot)) call errquit
     &   ('hnd_elpiso, ma_pop_stack of h_tepot failed',911,MA_ERR)

      if (.not.ma_pop_stack(h_zanpt)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of h_zanpt failed',911,MA_ERR)

      if (.not.ma_pop_stack(h_xyzpt)) call errquit
     &   ('hnd_elpiso, ma_pop_stack of h_xyzpt failed',911,MA_ERR)

C
C     Terminate integrals
C
      call schwarz_tidy()
      call int_terminate()

      end subroutine













C Returns the coordinates of grid points that lie within a TOLERANCE,
C on electron density isosurface of value ISOVALUE.
C Also calculates the number of voxels enclosed within the 
C isosurface, storing this in the rtdb as prop:epsiso:nenc.
C
C ISOVALUE == 0.002 ~ Vdw radius
C
C Mark J. Williamson Summer 2015
C Mark D. Driver Spring 2018
C
      subroutine get_points_from_iso(rtdb, basis, geom, isovalue,
     > tolerance, sub_nprp, i_sub_prp_c)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "geom.fh"
#include "rtdb.fh"
      integer rtdb                ! [Input] rtdb        
      integer basis               ! [Input] Basis set
      integer geom                ! [Input] Geometry
      double precision isovalue   ! [Input] Electron density iso surface to work on
      double precision tolerance  ! [Input] Tolerance of the isovalue
      integer sub_nprp            ! [Output] Number of points lying on this surface
      integer h_sub_prp_c         ! Memhandle to points lying on the iso surface
      integer i_sub_prp_c         ! [Output] Index to points lying on this surface
      integer sub_nprp_enc   ! [Output] Number of points lying within surface
      integer i,j,n,nat,iat
      integer g_dens(3)      ! Handle to density
      integer ndens          ! Number of active density handles (RHF=1, UHF=3)
      integer nclosed(2)     ! Occupation info
      integer nopen(2)       ! Occupation info
      integer nvirt(2)       ! Occupation info
      character*3 scftyp     ! Type of wave function
      character*16 at_tag    ! Geom centre tag

      integer nprp           ! Number of points on the entire density grid
      integer i_prp_c        ! Index to points on the entire density grid
      integer ngrid(3)       ! Grid size in x,y,z basis

      integer h_xyzpt        ! Memhandle to atoms coordinate array
      integer i_xyzpt        ! Index to atoms coordinate array
      integer h_zanpt        ! Memhandle to atoms charge array
      integer i_zanpt        ! Index to atoms charge array

      integer h_eden         ! Memhandle to electron density array
      integer i_eden         ! Index to electron density array
      double precision xp    ! Cartesian component of point on electron density isosurface
      double precision yp    ! Cartesian component of point on electron density isosurface
      double precision zp    ! Cartesian component of point on electron density isosurface

C     Determine how many atoms we have
      if (.not.geom_ncent(geom,nat)) call
     &    errquit('hnd_elpiso: geom_ncent',911,GEOM_ERR)

C     Allocate arrays for atomic positions and charges
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',h_xyzpt,i_xyzpt))
     &    call errquit('hnd_elpiso: h_xyzpt ma failed',911,MA_ERR)

      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',h_zanpt,i_zanpt))
     &    call errquit('hnd_elpiso: h_zanpt ma failed',911,MA_ERR)

C     Populate arrays with atom locations and their charge
      do 30 iat=1,nat
        if(.not.geom_cent_get(geom,iat,at_tag,dbl_mb(i_xyzpt+3*(iat-1)),
     &     dbl_mb(i_zanpt+iat-1))) call
     &     errquit('hnd_elpiso: geom_cent_get',911,GEOM_ERR)
   30 continue

C Set the grid up for the initial electron density calculation
      call prop_grid_initialize(rtdb,nat,dbl_mb(i_xyzpt))

C Get the number of points and their coordinate point index, i_prp_c
      call prop_grid_get_r_ptr(nprp,i_prp_c,ngrid)

C allocate memory and an index, i_eden, to the density values
      if (.not. ma_push_get(mt_dbl,nprp,'eden',h_eden,i_eden))
     &    call errquit('hnd_elpiso: eden h_eden ma failed',911,MA_ERR)

C Calculates density matrix
      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)

C Calculates the electronic contribution of the
C electronic integral defined by nder for a given density at
C the grid points defined in points.
C nder=-1; compute the 1-electron density
      call hnd_elfcon(basis,geom,g_dens(ndens),dbl_mb(i_prp_c),nprp, 
     & dbl_mb(i_eden),-1)

C Work out how many points lie within the given iso surface (sub_nprp_enc).
C Scan ahead and work out how many points (sub_nprp) lie on the isosurface
      sub_nprp = 0
      sub_nprp_enc = 0
      do i=1, nprp
        if (dbl_mb(i_eden + (i-1)) > isovalue) then
          sub_nprp_enc = sub_nprp_enc + 1
        end if
        if ( abs(dbl_mb(i_eden + (i-1)) - isovalue) < tolerance ) then
          sub_nprp = sub_nprp + 1
        end if
      end do

C Store number of points enclosed by surface in rtdb.
      if (.not. rtdb_put(rtdb, 'prop:epsiso:nenc', mt_int, 1,
     &                   sub_nprp_enc))
     &      call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
C Allocate array for storing the location of these points
      if(.not.ma_alloc_get(MT_DBL, 3*sub_nprp, 'hnd_elpiso grid array',
     &      h_sub_prp_c, i_sub_prp_c) ) call errquit(
     &      'grid_initialize: unable to allocate heap space',
     &      0, MA_ERR)

      sub_nprp = 0
      do i=1, nprp
        if ( abs(dbl_mb(i_eden + (i-1)) - isovalue) < tolerance ) then
C Coordinates of grid point
           xp = dbl_mb(i_prp_c  +3*(i-1))
           yp = dbl_mb(i_prp_c+1+3*(i-1))
           zp = dbl_mb(i_prp_c+2+3*(i-1))

C Copy matching coordinates to new array
           dbl_mb(i_sub_prp_c  +3*(sub_nprp)) = xp
           dbl_mb(i_sub_prp_c+1+3*(sub_nprp)) = yp
           dbl_mb(i_sub_prp_c+2+3*(sub_nprp)) = zp

           sub_nprp = sub_nprp + 1
        end if
      end do
 

C Clean up phase
      call prop_grid_destroy()

      if (.not.ma_pop_stack(h_eden)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of h_eden failed',911,MA_ERR)

      if (.not.ma_pop_stack(h_zanpt)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of h_zanpt failed',911,MA_ERR)

      if (.not.ma_pop_stack(h_xyzpt)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of h_xyzpt failed',911,MA_ERR)
    
1000  format(4F13.8)
      end subroutine

