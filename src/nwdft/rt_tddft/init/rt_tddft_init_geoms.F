C
C     Get geometry info from rtdb and load into struct.
C
C     Each geometry (fragment) has a "map" which is nbf_ao x nbf_ao, where nbf_ao
C     is the number of basis functions in the "active" geometry (ie set
C     geometry in input deck).  You can then project the full nbf_ao x nbf_ao
C     dipole integrals, etc onto each individual fragment via
C     ga_elem_multiply ().
C
C     NOTE: params%geom_active_handle refers to the nwchem geom
C     subsystem handle for "geometry", while params%active_geom refers
C     to the params%geom(:) array of geom structs.
C
      subroutine rt_tddft_init_geoms (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "inp.fh"
#include "geomP.fh"
#include "geom.fh"
#include "bas.fh"
#include "sym.fh"
#include "rt_tddft.fh"

C
C     Note: geomP.fh supplies ngeom_rtdb, names_rtdb, max_geom_rtdb, ndipole
C
      

C     == In/out ==
      type(rt_params_t) params  !geom params stored in here


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_init_geoms: "
      character(*), parameter :: geom_tmp_name = "rt_tmp_geom"
      double precision, parameter :: dist_thresh = 1d-8

      
C     == Variables ==
      integer me
      integer ig
      
      integer geom_tmp
      character*16 geom_rtdb_name, geom_excite_name
      character*16 active_geom_name
      integer geom_tmp_ncent

      integer ibf, icen, jcen
      character*16 icen_tag, jcen_tag
      double precision icen_loc(3), jcen_loc(3)
      double precision icen_charge, jcen_charge
      double precision dist
      logical bf_in_geom
      integer nbf_geom
      integer per_bf, per_cen
      double precision geom_charge_nuc
      logical got_excite
      integer is

C      character*25 role_string
      character*2 role_string

      logical found_active

      double precision ndip(3)


      me = ga_nodeid ()


C
C     Check that rt hardcoded max num geoms is less than = nwchem hardlimit
C
      if (rt_max_geoms .lt. nw_max_geoms) !hardcoded in nwc_const.fh
     $     call errquit (pname//"rt_max_geoms < nw_max_geoms",0,0)


C
C     Get geom info from rtdb
C
      if (.not.rtdb_get (params%rtdb, "geometry:ngeom",
     $     mt_int, 1, ngeom_rtdb))
     $     call errquit (pname//"failed to get ngeom from rtdb",0,0)
      
      if (ngeom_rtdb .lt. 1)
     $     call errquit (pname//"invalid ngeom_rtdb < 1",0,0)

c$$$      if (ngeom_rtdb .gt. nw_max_geoms)
c$$$     $     call errquit (pname//"invalid ngeom_rtdb > nw_max_geoms",0,0)

      if (ngeom_rtdb .gt. max_geom_rtdb)
     $     call errquit(pname//"invalid ngeom_rtdb > max_geom_rtdb",0,0)
 
      if (.not. rtdb_cget (params%rtdb,'geometry:names',
     $     ngeom_rtdb,names_rtdb))
     $     call errquit (pname//"failed to read names from rtdb",0,0)


C
C     Get internal index for active geometry from the rtdb name for the
C     geometry.
C

      if (.not. rtdb_cget (params%rtdb,'geometry',1,active_geom_name))
     $     call errquit (pname//"failed to read active geom from rtdb",
     $     0,0)

      found_active = .false.
      
      do ig = 1, ngeom_rtdb
         if ( trim(names_rtdb(ig)) .eq. trim(active_geom_name) ) then
            if (found_active) then
               call errquit (pname//"found multiple active geoms",0,0)
            else
               params%geom_active = ig
               found_active = .true.
            endif
         endif
      enddo

      if (.not. found_active)
     $     call errquit (pname//"failed to find active geom",0,0)


C
C     Loop over all geoms found, create a temporary geometry object for
C     each, calculate some properties, store in the rt struct, and
C     delete the geom object (leaving geoms untouched).
C
      do ig = 1, ngeom_rtdb

C     (get from rtdb)
         if (.not.geom_create (geom_tmp, geom_tmp_name))
     $        call errquit (pname//"geom_create failed",0,0)
         
         if (.not. geom_check_handle (geom_tmp,pname))
     $        call errquit (pname//"not a valid geom",0,0)

         if (.not. geom_rtdb_load (params%rtdb,
     $        geom_tmp, names_rtdb(ig)))
     $        call errquit (pname//"failed to load geom",0,0)

C     (extract info from geom)
         if (.not. geom_ncent (geom_tmp, geom_tmp_ncent))
     $        call errquit (pname//"geom_ncent failed",0,0)

         if (.not. geom_nuc_dipole (geom_tmp, ndip))
     $        call errquit (pname//"geom_nucl_dipole failed",0,0)

         geom_rtdb_name = names_rtdb(ig)

C     
C     Check that geometry fragment is not larger than active geometry.
C
         if (geom_tmp_ncent .gt. params%natoms) call errquit (pname//
     $        'geom "'//trim(geom_rtdb_name)//

     $        '" has more atoms than active geometry',
     $        geom_tmp_ncent, GEOM_ERR)

         
C     (alloc and compute mask GA)
         if (.not.ga_create (mt_dbl, params%nbf_ao, params%nbf_ao,
     $        trim(geom_rtdb_name)//"_mask", 0, 0,
     $        params%geom(ig)%g_mask))
     $        call errquit (pname//"alloc mask failed", 0, GA_ERR)


C
C     Mask starts as all 1's, then we 0 any row/column which refers to a
C     basis function not in this geom.
C
         call ga_fill (params%geom(ig)%g_mask, 1d0)
         

C
C     Loop over basis functions in overall/active basis set.  For each
C     basis function, find the "owning" atom, and then check if that
C     atom is part of the "ig" geometry.
C
         nbf_geom = 0
         do ibf = 1, params%nbf_ao
            if (.not. bas_bf2ce (params%ao_bas_han, ibf, icen))
     $           call errquit (pname//"bas_bf2ce failed", 0, 0)
            
C     (note this acts on full active geom, specified by the handle
C     stored in params)
            if (.not. geom_cent_get (params%geom_active_handle, icen,
     $           icen_tag, icen_loc, icen_charge))
     $           call errquit (pname//"geom_cent_get active failed",0,0)

            bf_in_geom = .false.
            do jcen = 1, geom_tmp_ncent  !increment over atoms in "ig" geom

               if (.not. geom_cent_get (geom_tmp, jcen,
     $              jcen_tag, jcen_loc, jcen_charge))
     $              call errquit(pname//"geom_cent_get jcen failed",0,0)
               
C     (icen_loc = location of ibf; jcen_loc = location of jcen'th atom/center in geom)
               dist = dsqrt (
     $              ( icen_loc(1) - jcen_loc(1) )**2 +
     $              ( icen_loc(2) - jcen_loc(2) )**2 +
     $              ( icen_loc(3) - jcen_loc(3) )**2 )

               if (dist < dist_thresh) bf_in_geom = .true.
               
            enddo ! jcen loop


C     (this basis func not in this geom, zero corresponding rows+columns)
            if (.not. bf_in_geom) then
               call ga_fill_patch (params%geom(ig)%g_mask,
     $              1, params%nbf_ao, ibf, ibf, 0d0)
               call ga_fill_patch (params%geom(ig)%g_mask,
     $              ibf, ibf, 1, params%nbf_ao, 0d0)
            else ! this bf in this geom
               nbf_geom = nbf_geom + 1
            endif
         enddo ! ibf loop

         
C     (compute total nuclear charge on this geom)
         geom_charge_nuc = 0d0
         do jcen = 1, geom_tmp_ncent

            if (.not. geom_cent_get (geom_tmp, jcen,
     $           jcen_tag, jcen_loc, jcen_charge))
     $           call errquit(pname//"geom_cent_get jcen failed",0,0)
            
            geom_charge_nuc =  geom_charge_nuc + jcen_charge
         enddo


C     (store misc params)
         params%ngeoms = ngeom_rtdb
         params%geom(ig)%name  = trim (geom_rtdb_name)
         params%geom(ig)%ncent = geom_tmp_ncent
         params%geom(ig)%nbf = nbf_geom
         params%geom(ig)%charge_nuc = geom_charge_nuc
         params%geom(ig)%ndip%x = ndip(1)
         params%geom(ig)%ndip%y = ndip(2)
         params%geom(ig)%ndip%z = ndip(3)
         params%geom(ig)%measure = .true. !XXX HARDCODED MEASURE ALL

         
C     (remove tmp geom object, rtdb untouched)
         if (.not. geom_destroy (geom_tmp))
     $        call errquit (pname//"geom_destroy failed", 0, 0)

      enddo ! ig loop

      
C
C     Print geom info to stdout
C
      if (me .eq.0) then
         write (luout, *) " "
         write (luout, *) " "
         
         write (luout, *) "Geometry               "//
     $        "Atoms      Basis func.   "//
     $        "Nuc. charge          Nuc. dip. mom."

         write (luout, *) "-----------------"//
     $        "------------------------------------"//
     $        "------------------------------------"

      endif

      do ig = 1, params%ngeoms

         per_bf = 100 * params%geom(ig)%nbf / params%nbf_ao
         per_cen = 100 * params%geom(ig)%ncent / params%natoms

C         role_string = " "

C         if ( trim(params%geom(ig)%name) .eq. trim(active_geom_name))
C     $        role_string = trim(role_string) // "<= Active "

         if (ig .eq. params%geom_active) then
            role_string = " *"
         else
            role_string = "  "
         endif

         if (me.eq.0) then

            write (luout,
C     $           "(2x,a16,2x,i6,a,i3,a,5x,i6,a,i3,a,4x,1f8.2,4x,a)")
     $           "(a,a16,i6,a,i3,a,i6,a,i3,a,2x,1f8.2,2x,3f10.2)")
     $           role_string,
     $           params%geom(ig)%name,
     $           params%geom(ig)%ncent, " (", per_cen , "% )", 
     $           params%geom(ig)%nbf,   " (", per_bf, "% )",
     $           params%geom(ig)%charge_nuc,
     $           params%geom(ig)%ndip%x,
     $           params%geom(ig)%ndip%y,
     $           params%geom(ig)%ndip%z
         endif
         
         call util_flush (luout)
      enddo


C
C     Now that we have determined the active geometry, update any
C     excitations with -1 geom (ie if exciting orbitals directly) to
C     excite the entire active geometry.
C
c$$$      do is = 1, params%nexcites
c$$$         if (params%excite(is)%geom_indx .lt. 1) then
c$$$            params%excite(is)%geom_indx = params%geom_active
c$$$         endif
c$$$      enddo



C
C     It appears that dipole integrals only work for Abelian symmetry,
C     so we enforce that for each geom.
C
C
      if (.not.sym_abelian_group(params%geom_active_handle))
     $     call errquit (pname//
     $     "active geom: non-Abelian symmetry; disable autosym",
     $     0, GEOM_ERR)

C     (debug)
C      call util_flush (luout)
C      if (.not. bas_print_all())
C     $     call errquit (pname//"basis print failed")
C      call util_flush (luout)
      
      end subroutine

c $Id$
