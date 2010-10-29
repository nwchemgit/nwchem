c
c $Id$
c
c
      subroutine qmmm_bq_energy_gen(rtdb,n,c,q,e)
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
c
      integer rtdb
      integer n
      double precision c(3,n)
      double precision q(n)
      double precision e

      character*32 pname
      character*3 theory
      character*255 vectors_in
c
      integer i
      logical ignore
      integer h_epot,i_epot
c
      logical property
      external property
c
      pname = "qmmm_bq_energy"
c
      if (.not. ma_push_get(mt_dbl,n,'qmmm epot',h_epot,i_epot))
     &    call errquit(pname//' ma failed',911,MA_ERR)

c
c     make bq coords available for property calcs 
c     -------------------------------------------------------- 
      if (.not. rtdb_put(rtdb,'prop:xyz',
     >                   mt_dbl,3*n,c))
     &    call errquit(pname//'prop:xyz failed',911,RTDB_ERR)
c
c     export movecs file to property to avoid doing energy again
c     ----------------------------------------------------------
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit(pname//'theory not specified',555,
     &       INPUT_ERR)

       if (theory.eq.'dft') then
          if (.not. rtdb_cget(rtdb,'dft:input vectors',1,vectors_in))
     &        call errquit(pname//'rtdb_cget failed', 100,
     &                     RTDB_ERR)
       elseif (theory.eq.'scf') then
          if (.not. rtdb_cget(rtdb,'scf:input vectors',1,vectors_in))
     &        call errquit(pname//'rtdb_cget failed', 100,
     &                     RTDB_ERR)
       endif

      if (.not. rtdb_cput(rtdb,'prop:vectors',1,vectors_in))
     &    call errquit(pname//'rtdb_cput failed', 100,
     &                 RTDB_ERR)
c
      if (.not. rtdb_put(rtdb,'prop:esp',
     >                   mt_int,1,0))
     &    call errquit(pname//'prop: efield failed',911,RTDB_ERR)
c
      if (.not. rtdb_put(rtdb,'prop:nocheck',
     >                   mt_log,1,.true.))
     &    call errquit(pname//'prop:nocheck failed',911,RTDB_ERR)
c
      if(.not.property(rtdb)) 
c
     &    call errquit(pname//'failed property call', 100,
     &                 RTDB_ERR)

c
c     get back electric field at bq positions
c     ---------------------------------------
      if (.not. rtdb_get(rtdb,'prop:epot_xyz',
     >                   mt_dbl,n,dbl_mb(i_epot)))
     &    call errquit(pname//'prop:efield failed',911,RTDB_ERR)
c
      e = 0.0d0
      do i=1,n
        e = e + q(i)*dbl_mb(i_epot+i-1)
      end do
c
      
      if(ga_nodeid().eq.0) then
        write(luout,*) " Bq interaction energy", e
      end if
c 
c     clean up rtdb
c     ------------
      ignore = rtdb_delete(rtdb,"prop:epot_xyz")
      ignore = rtdb_delete(rtdb,"prop:xyz")
      ignore = rtdb_delete(rtdb,"prop:vectors")
      ignore = rtdb_delete(rtdb,"prop:nocheck")

      if (.not. rtdb_put(rtdb,'prop:esp',
     >                   mt_int,1,1))
     &    call errquit(pname//'prop: efield failed',911,RTDB_ERR)

      if (.not.ma_pop_stack(h_epot)) call errquit
     &   (pname//'ma_pop_stack of h_epot failed',911,MA_ERR)

      end

      subroutine qmmm_bq_forces(rtdb,n,c,q,f)
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "qmmm.fh"
c
      integer rtdb
      integer n
      double precision c(3,n)
      double precision q(n)
      double precision f(3,n)

      character*32 pname
      character*32 theory
      character*255 vectors_in
      character*255 ostring
c
      integer i
      logical ignore
c
      logical property
      external property
c
      pname = "qmmm_bq_forces"
c
      if(qmmm_master()) then
        call util_print_centered(luout,
     >  "Calculating forces on classical charges",
     >  32,.true.)
        write(luout,*) "Total number of charges",n
      end if
c
c     make bq coords available for property calcs 
c     -------------------------------------------------------- 
      if (.not. rtdb_put(rtdb,'prop:xyz',
     >                   mt_dbl,3*n,c))
     &    call errquit(pname//'prop:xyz failed',911,RTDB_ERR)
c
c     export movecs file to property to avoid doing energy again
c     ----------------------------------------------------------
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit(pname//'theory not specified',555,
     &       INPUT_ERR)

       if (theory.eq.'dft') then
          if (.not. rtdb_cget(rtdb,'dft:input vectors',1,vectors_in))
     &        call errquit(pname//'rtdb_cget failed', 100,
     &                     RTDB_ERR)
       elseif (theory.eq.'scf') then
          if (.not. rtdb_cget(rtdb,'scf:input vectors',1,vectors_in))
     &        call errquit(pname//'rtdb_cget failed', 100,
     &                     RTDB_ERR)
       endif

      if (.not. rtdb_cput(rtdb,'prop:vectors',1,vectors_in))
     &    call errquit(pname//'rtdb_cput failed', 100,
     &                 RTDB_ERR)
c
      if (.not. rtdb_put(rtdb,'prop:efield',
     >                   mt_int,1,0))
     &    call errquit(pname//'prop: efield failed',911,RTDB_ERR)
c
      if (.not. rtdb_put(rtdb,'prop:nocheck',
     >                   mt_log,1,.true.))
     &    call errquit(pname//'prop:nocheck failed',911,RTDB_ERR)
c
      if(.not.property(rtdb)) 
c
     &    call errquit(pname//'failed property call', 100,
     &                 RTDB_ERR)

c
c     get back electric field at bq positions
c     ---------------------------------------
      if (.not. rtdb_get(rtdb,'prop:efield_xyz',
     >                   mt_dbl,3*n,f))
     &    call errquit(pname//'prop:efield failed',911,RTDB_ERR)
c
      do i=1,n
        f(1,i) = -q(i)*f(1,i) 
        f(2,i) = -q(i)*f(2,i) 
        f(3,i) = -q(i)*f(3,i) 
      end do
c 
c     clean up rtdb
c     ------------
      ignore = rtdb_delete(rtdb,"prop:efield_xyz")
      ignore = rtdb_delete(rtdb,"prop:xyz")
      ignore = rtdb_delete(rtdb,"prop:vectors")

      if (.not. rtdb_put(rtdb,'prop:efield',
     >                   mt_int,1,1))
     &    call errquit(pname//'prop: efield failed',911,RTDB_ERR)

      end
