c
c $Id: qmmm_bq.F,v 1.1 2004-07-21 01:29:47 marat Exp $
c
c

      subroutine qmmm_load_bq(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
      integer irtdb
c     local variables
      integer nbqs
      integer nbqw
      integer i_cbqs,h_cbqs
      integer i_qbqs,h_qbqs
      integer i_cbqw,h_cbqw
      integer i_qbqw,h_qbqw
      character*4 bqtag
      character*16 myname

      call qmmm_unload_bq()

      call mm_get_tot_nbqs(nbqs)
      call mm_get_tot_nbqw(nbqw)
       write(*,*) "qmmm_load_bq", nbqs


      call qmmm_set_nbqw(nbqw)
      call qmmm_set_nbqs(nbqs)
      
      if(nbqs.eq.0) return

      if(.not.ma_alloc_get(mt_dbl,3*nbqs,'cbq',h_cbqs,i_cbqs))
     + call errquit('qmmm: Failed to allocate memory for cbq',3*nbqs,
     +       MA_ERR)

      if(.not.ma_alloc_get(mt_dbl,nbqs,'qbq',h_qbqs,i_qbqs))
     + call errquit('qmmm: Failed to allocate memory for qbq',nbqs,
     &       MA_ERR)

       write(*,*) "tot_nbq=",nbqs
     
      if(nbqs.ne.0) then

      call mm_get_geom_bqs(nbqs,dbl_mb(i_cbqs),dbl_mb(i_qbqs))
      write(*,*) "after mm_get_geom_bqs"
      myname = "test-bq.pdb"
      bqtag="Bq"
      call qmmm_print_bq_pdb(nbqs,myname,dbl_mb(i_cbqs),
     >      dbl_mb(i_qbqs),bqtag)

      end if

      call qmmm_set_cbqs(i_cbqs,h_cbqs)
      call qmmm_set_qbqs(i_qbqs,h_qbqs)

      end

      subroutine qmmm_unload_bq()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
c     local variables
      integer i
      integer nbqs
      integer nbqw
      integer i_cbqs,h_cbqs
      integer i_qbqs,h_qbqs
      integer i_cbqw,h_cbqw
      integer i_qbqw,h_qbqw

      call mm_get_tot_nbqs(nbqs)
      call mm_get_tot_nbqw(nbqw)

      if(nbqs.gt.0) then
        if(.not.ma_alloc_get(mt_dbl,3*nbqs,'cbq',h_cbqs,i_cbqs))
     +   call errquit('qmmm: Failed to allocate memory for cbq',3*nbqs,
     +         MA_ERR)

        if(.not.ma_alloc_get(mt_dbl,nbqs,'qbq',h_qbqs,i_qbqs))
     +   call errquit('qmmm: Failed to allocate memory for qbq',nbqs,
     &         MA_ERR)
      end if

      nbqs = 0
      nbqw = 0
      call qmmm_set_nbqw(nbqw)
      call qmmm_set_nbqs(nbqs)
      
      end
 
