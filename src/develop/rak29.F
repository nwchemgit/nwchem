      logical function raktask_fullsc(rtdb)
      implicit none
#include "stdio.fh"
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
      integer rtdb
c
      integer basis, geom
      integer ma_hp, ncart, nsph
      raktask_fullsc = .false.
c
      if (.not.geom_create(geom,'geometry')) call errquit
     &    ('geom create failed',911)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('geom_rtdb_load failed',911)
c
      if (.not.bas_create(basis,'ao basis')) call errquit
     &    ('bas_create failed',911)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis')) call errquit
     &    ('bas_rtdb_load failed',911)
c
      write(luout,*)' geom/basis loaded'
      if (.not. bas_print(basis))
     $    call errquit(' basis print failed', 0)
      call util_flush(luout)
c
      call spcart_get_fullcs_ma(basis,ma_hp,ncart,nsph)
      write(luout,*)' heap ma handle ',ma_hp
      write(luout,*)' cart dimension ',ncart
      write(luout,*)' spherical dim  ',nsph
c
      raktask_fullsc = ma_free_heap(ma_hp)
      raktask_fullsc = raktask_fullsc.and.bas_destroy(basis)
      raktask_fullsc = raktask_fullsc.and.geom_destroy(geom)
c
      end
