      logical function hnd_property(rtdb)
*
* $Id: hnd_property.F,v 1.12 2000-10-23 15:53:25 windus Exp $
*
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "cfock.fh"
c
c     property module.
c
c     Assumes energy has been completed, MO vectors stored
c     and all information is still in the RTDB
c
c
      logical  int_normalize
      external int_normalize

      integer rtdb              ! [input] database handle
      integer geom, basis       ! handles
      logical status
      character*255 title
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
      call util_print_push
c
c     Extract high level info from the data-base setting defaults
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('property: geom_create?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('property: no geometry ', 0)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('property: bas_create?', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('property: no ao basis ', 0)
      if (.not. int_normalize(rtdb,basis))
     $     call errquit('property: normalization failed',911)
c
      if (ga_nodeid().eq.0) then
         if (util_print('information',print_low)) then
            call util_print_centered(luout,'NWChem Property Module',
     $           40,.true.)
            write(luout,*)
            call util_flush(luout)
          endif
          if (util_print('information',print_medium)) then
            write(luout,*)
            if (title .ne. ' ') then
               call util_print_centered(luout, title, 40, .false.)
               write(luout,*)
            endif
            call util_flush(luout)
         endif
         if (util_print('geometry',print_high)) then
            if (.not. geom_print(geom)) 
     $           call errquit('property: geom_print ?',0)
            call util_flush(luout)
         endif
         if (util_print('basis',print_high)) then
            if (.not. bas_print(basis))
     $           call errquit('property: bas_print ?',0)
            call util_flush(luout)
         endif
      endif
c
c     call hondo property calculation routines
c
      if(.not.bas_is_spherical(basis)) then
        if (.not.bas_any_gcorsp(basis)) then
          call hnd_propty ( rtdb, basis, geom )
        else
          write(luout,99991)
        endif
      else
         write(luout,99992)
      endif
c
c property is done. destroy basis and geometry handles  
c (e.g., preserve the memory available to other modules!!)
c
      if (.not.(
     &    (bas_destroy(basis))
     &    .and.
     &    (geom_destroy(geom))
     &    ))
     &    call errquit
     &    ('property:error destroying geom and basis handles',911)
c
      call util_print_pop
c
      hnd_property = .true.
c
      return
99991 format(' property routines lifted from -hondo- do not handle',
     &       ' sp functions or general contractions ... .',/,
     &       ' please use the "segment" keyword on the basis set ',/,
     &       ' skip and continue ')
99992 format(' property routines lifted from -hondo- do not handle',
     &       ' spherical basis sets ... .',/,
     &       ' skip and continue ')
      end
