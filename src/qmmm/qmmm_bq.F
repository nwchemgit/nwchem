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
#include "cfock.fh"
#include "prop.fh"
c
      integer rtdb
      integer n
      double precision c(3,n)
      double precision q(n)
      double precision f(3,n)

      character*32 pname
      character*3 theory
      character*255 vectors_in
c
      integer i
      logical ignore
c
      logical property
      external property
c
      pname = "qmmm_bq_forces"
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
        write(*,*) f(1,i),f(2,i),f(3,i)
      end do
c 
c     clean up rtdb
c     ------------
      ignore = rtdb_delete(rtdb,"prop:efield_xyz")
      ignore = rtdb_delete(rtdb,"prop:xyz")
      ignore = rtdb_delete(rtdb,"prop:vectors")

      end
