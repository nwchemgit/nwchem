      subroutine geom_fragment(rtdb,nfrags_return)
c
c $Id: geom_fragment.F,v 1.1 2004-02-01 22:03:50 sohirata Exp $
c
c     Subroutines and functions in this file are to help
c     perform almost automatically the fragment MO for 
c     various wave functions, if the fragments are not 
c     covalently bonded.  Fragments are separated by
c     the proximity ansatz for atoms.  See, e.g.,
c
c     K.Kitaura, T.Sawai, T.Asada, T.Nakano, & M.Uebayashi,
c     Chem.Phys.Lett. 312, 319 (1999).
c
c     Written by So Hirata in collaboration with Hideo Sekino. 
c
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "geom.fh"
#include "geom_fragment.fh"
c
      integer rtdb             ! Runtime database handle
      integer geom             ! Geometry handle
      integer ncent
      integer iatom, jatom
      logical done
      integer nfrags_return
      integer ifrag
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (.not.geom_create(geom,'geometry'))
     1  call errquit('geom_fragment: geom_create failed',0, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     1  call errquit('geom_fragment: no geometry found',0, RTDB_ERR)
      if (.not.geom_rtdb_store(rtdb,geom,'original'))
     1  call errquit('geom_fragment: failed to store geom',0, RTDB_ERR)
c     write(LuOut,9000) "Original geometry"
c     if (.not.geom_print(geom))
c    1  call errquit('geom_fragment: geom_print',0,GEOM_ERR)
      if (.not.geom_cart_get(geom,ncent,tags,coords,charges))
     1  call errquit('geom_fragment: geom_cart_get',0, RTDB_ERR)
c
c     Determine fragments
c
      do iatom = 1, ncent
        f_map(iatom) = -1
      enddo
      nfrags = 0
   10 done = .true.
      jatom = 0
      do iatom = ncent, 1, -1
        if (f_map(iatom).eq.-1) jatom = iatom
      enddo
      if (jatom.ne.0) then
        done = .false.
        nfrags = nfrags + 1
        f_map(jatom) = nfrags
        do iatom = 1, ncent
          if ((f_map(iatom).eq.-1).and.
     1       (dsqrt((coords(1,iatom)-coords(1,jatom))**2
     2             +(coords(2,iatom)-coords(2,jatom))**2
     3             +(coords(3,iatom)-coords(3,jatom))**2)
     4       .le.geom_thresh)) then
            f_map(iatom) = nfrags
          endif
        enddo
      endif
      if (.not.done) goto 10
      if (.not.geom_destroy(geom))
     1  call errquit('geom_fragment: geom_destroy failed',0,GEOM_ERR)
      if (nodezero) write(LuOut,*)
      if (nodezero) write(LuOut,9020) "@ Fragmentation"
      if (nodezero) write(LuOut,*)
      do ifrag = 1, nfrags
        if (nodezero) write(LuOut,9020) '@ Fragment',ifrag
        do iatom = 1, ncent
          if (f_map(iatom).eq.ifrag) then
            if (nodezero) write(LuOut,9020) '@ ',iatom,tags(iatom),
     1        charges(iatom),
     2        coords(1,iatom),coords(2,iatom),coords(3,iatom)
          endif
        enddo
      enddo
      if (nodezero) write(LuOut,*)
      nfrags_return = nfrags
      return
 9000 format(3x,A,i3)
 9010 format('.',i3.3)
 9020 format(1x,A,i3,3x,A5,f10.4,3f15.10)
      end



      subroutine select_fragment(rtdb,ifrag,dipole)
c
c $Id: geom_fragment.F,v 1.1 2004-02-01 22:03:50 sohirata Exp $
c
c     Overwrite the main geometry entry in RTDB by that of a fragment.
c
c     Written by So Hirata in collaboration with Hideo Sekino.
c
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "geom.fh"
#include "geom_fragment.fh"
c
      integer rtdb             ! Runtime database handle
      integer ifrag
      double precision dipole(3,*)
      integer jfrag
      integer f_geom
      integer ncent
      integer iatom
      integer jatom
      double precision cx, cy, cz, c, d
      double precision e
c
      if (.not.rtdb_get(rtdb,'fragment:dipole_length',mt_dbl,1,e))
     1  e = 0.1d0
      if (.not.geom_create(f_geom,'fragment'))
     1  call errquit('geom_fragment: geom_create failed',0, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,f_geom,'original'))
     1  call errquit('geom_fragment: no geometry found',0, RTDB_ERR)
      if (.not.geom_ncent(f_geom,ncent))
     1  call errquit('geom_fragment: no geometry found',0, RTDB_ERR)
      jatom = 0
      do jfrag = 1, nfrags
        if (jfrag.eq.ifrag) then
          do iatom = 1, ncent
            if (f_map(iatom).eq.jfrag) then
              jatom = jatom + 1
              f_coords(1,jatom) = coords(1,iatom)
              f_coords(2,jatom) = coords(2,iatom)
              f_coords(3,jatom) = coords(3,iatom)
              f_tags(jatom) = tags(iatom)
              f_charges(jatom) = charges(iatom)
            endif
          enddo
        else
          cx = 0.0d0
          cy = 0.0d0
          cz = 0.0d0
          c  = 0.0d0
          do iatom = 1, ncent
            if (f_map(iatom).eq.jfrag) then
              cx = cx + charges(iatom) * coords(1,iatom)
              cy = cy + charges(iatom) * coords(2,iatom)
              cz = cz + charges(iatom) * coords(3,iatom)
              c  = c  + charges(iatom)
            endif
          enddo
          cx = cx / c
          cy = cy / c
          cz = cz / c
          d  = dsqrt(dipole(1,jfrag)**2
     1              +dipole(2,jfrag)**2
     2              +dipole(3,jfrag)**2)
          if (d.ne.0.0d0) then
            jatom = jatom + 1
            f_coords(1,jatom) = cx + dipole(1,jfrag) * e/d
            f_coords(2,jatom) = cy + dipole(2,jfrag) * e/d
            f_coords(3,jatom) = cz + dipole(3,jfrag) * e/d
            f_tags(jatom) = 'Bq'
            f_charges(jatom) = 0.5d0 * d/e
            jatom = jatom + 1
            f_coords(1,jatom) = cx - dipole(1,jfrag) * e/d
            f_coords(2,jatom) = cy - dipole(2,jfrag) * e/d
            f_coords(3,jatom) = cz - dipole(3,jfrag) * e/d
            f_tags(jatom) = 'Bq'
            f_charges(jatom) = - 0.5d0 * d/e
          endif
        endif
      enddo
      if (.not.geom_cart_set(f_geom,jatom,
     1  f_tags,f_coords,f_charges))
     2  call errquit('geom_fragment: geom_cart_set',0, RTDB_ERR)
      if (.not.geom_print(f_geom))
     1  call errquit('geom_fragment: geom_create failed',0, GEOM_ERR)
      if (.not.geom_rtdb_store(rtdb,f_geom,'geometry'))
     1  call errquit('geom_fragment: failed to store geom',0,RTDB_ERR)
      if (.not.geom_destroy(f_geom))
     1  call errquit('geom_fragment: geom_create failed',0, GEOM_ERR)
 9010 format('.',i3.3)
      return
      end
 
 
 
      subroutine select_fragmentpair(rtdb,ifrag,jfrag,dipole)
c
c $Id: geom_fragment.F,v 1.1 2004-02-01 22:03:50 sohirata Exp $
c
c     Overwrite the main geometry entry in RTDB by that of a fragment pair.
c
c     Written by So Hirata in collaboration with Hideo Sekino.
c
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "geom.fh"
#include "geom_fragment.fh"
#include "inp.fh"
c
      integer rtdb             ! Runtime database handle
      integer ifrag
      integer jfrag
      double precision dipole(3,*)
      integer kfrag
      integer iatom
      integer jatom
      integer f_geom
      integer ncent
      double precision cx, cy, cz, c, d
      double precision e
c
      if (.not.rtdb_get(rtdb,'fragment:dipole_length',mt_dbl,1,e))
     1  e = 0.1d0
      if (.not.geom_create(f_geom,'fragmentpair'))
     1  call errquit('geom_fragment: geom_create failed',0,GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,f_geom,'original'))
     1  call errquit('geom_fragment: no geometry found',0, RTDB_ERR)
      if (.not.geom_ncent(f_geom,ncent))
     1  call errquit('geom_fragment: no geometry found',0, RTDB_ERR)
      jatom = 0
      do kfrag = 1, nfrags
        if ((kfrag.eq.ifrag).or.(kfrag.eq.jfrag)) then
          do iatom = 1, ncent
            if (f_map(iatom).eq.kfrag) then
              jatom = jatom + 1
              f_coords(1,jatom) = coords(1,iatom)
              f_coords(2,jatom) = coords(2,iatom)
              f_coords(3,jatom) = coords(3,iatom)
              f_tags(jatom) = tags(iatom)
              f_charges(jatom) = charges(iatom)
            endif
          enddo
        else
          cx = 0.0d0
          cy = 0.0d0
          cz = 0.0d0
          c  = 0.0d0
          do iatom = 1, ncent
            if (f_map(iatom).eq.kfrag) then
              cx = cx + charges(iatom) * coords(1,iatom)
              cy = cy + charges(iatom) * coords(2,iatom)
              cz = cz + charges(iatom) * coords(3,iatom)
              c  = c  + charges(iatom)
            endif
          enddo
          cx = cx / c
          cy = cy / c
          cz = cz / c
          d  = dsqrt(dipole(1,kfrag)**2
     1              +dipole(2,kfrag)**2
     2              +dipole(3,kfrag)**2)
          if (d.ne.0.0d0) then
            jatom = jatom + 1
            f_coords(1,jatom) = cx + dipole(1,kfrag) * e/d
            f_coords(2,jatom) = cy + dipole(2,kfrag) * e/d
            f_coords(3,jatom) = cz + dipole(3,kfrag) * e/d
            f_tags(jatom) = 'Bq'
            f_charges(jatom) = 0.5d0 * d/e
            jatom = jatom + 1
            f_coords(1,jatom) = cx - dipole(1,kfrag) * e/d
            f_coords(2,jatom) = cy - dipole(2,kfrag) * e/d
            f_coords(3,jatom) = cz - dipole(3,kfrag) * e/d
            f_tags(jatom) = 'Bq'
            f_charges(jatom) = - 0.5d0 * d/e
          endif
        endif
      enddo
      if (.not.geom_cart_set(f_geom,jatom,
     1  f_tags,f_coords,f_charges))
     2  call errquit('geom_fragment: geom_cart_set',0, RTDB_ERR)
      if (.not.geom_print(f_geom))
     1  call errquit('geom_fragment: geom_create failed',0,GEOM_ERR)
      if (.not.geom_rtdb_store(rtdb,f_geom,'geometry'))
     1  call errquit('geom_fragment: failed to store geom',0,RTDB_ERR)
      if (.not.geom_destroy(f_geom))
     1  call errquit('geom_fragment: geom_create failed',0,GEOM_ERR)
 9010 format('.',i3.3)
      return
      end
