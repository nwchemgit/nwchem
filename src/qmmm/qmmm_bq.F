c
c $Id: qmmm_bq.F,v 1.4 2004-08-27 22:26:59 marat Exp $
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
      integer nbq
      integer nbqs
      integer nbqw
      integer i_cbq,h_cbq
      integer i_qbq,h_qbq
      character*4 bqtag
      character*16 myname

c
c     clean up any prior bq
c     information
c     ---------------------
      call qmmm_unload_bq()
c
c     get total number of solute and
c     solvent Bq's from MM and save them
c     -----------------------------
      call mm_get_tot_nbqs(nbqs)
      call qmmm_set_nbqs(nbqs)

      call mm_get_tot_nbqw(nbqw)
      call qmmm_set_nbqw(nbqw)
c      
c     this is a total number of bq's (nbqs+nbqw)
c     ------------------------------------------
      nbq = qmmm_get_nbq()
c
c     exit if no bq's
c     ---------------
      write(*,*) "debug: nbq", nbq,nbqs,nbqw
      if(nbq.eq.0) return

c
c     allocate memory for coord and charges
c     -------------------------------------
      if(.not.ma_alloc_get(mt_dbl,3*nbq,'cbq',h_cbq,i_cbq))
     + call errquit('qmmm: Failed to allocate memory for cbq',3*nbq,
     +       MA_ERR)
      if(.not.ma_alloc_get(mt_dbl,nbq,'qbq',h_qbq,i_qbq))
     + call errquit('qmmm: Failed to allocate memory for qbq',nbq,
     &       MA_ERR)

      call qmmm_set_cbq(i_cbq,h_cbq)
      call qmmm_set_qbq(i_qbq,h_qbq)
c
c     load solute bq coord and charges from MM
c     ---------------------------------------    
      if(nbqs.ne.0) then
        call mm_get_geom_bqs(nbqs,dbl_mb(i_cbq),dbl_mb(i_qbq))
        myname = "test-bqs.pdb"
        bqtag="Bq"
        call qmmm_print_bq_pdb(nbqs,myname,dbl_mb(i_cbq),
     >      dbl_mb(i_qbq),bqtag)

      end if
c
c     append solvent bq coord and charges from MM
c     -------------------------------------------    
      if(nbqw.ne.0) then
        call mm_get_geom_bqw(nbqw,dbl_mb(i_cbq+3*nbqs),
     >                            dbl_mb(i_qbq+nbqs))
        myname = "test-bqw.pdb"
        bqtag="Bq"
        call qmmm_print_bq_pdb(nbqw,myname,dbl_mb(i_cbq+3*nbqs),
     >      dbl_mb(i_qbq+nbqs),bqtag)

      end if
      end

      subroutine qmmm_unload_bq()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
c     local variables
      integer nbq
      integer h_cbq
      integer h_qbq

      nbq   = qmmm_get_nbq()
      if(nbq.gt.0) then
        h_cbq = qmmm_get_h_cbq()
        h_qbq = qmmm_get_h_qbq()

        if(.not.ma_free_heap(h_cbq))
     +   call errquit('qmmm: memory deallocation for cbq',3*nbq,
     +         MA_ERR)

        if(.not.ma_free_heap(h_qbq))
     +   call errquit('qmmm:memory deallocation for qbq',nbq,
     &         MA_ERR)
        call qmmm_set_nbqs(0)
        call qmmm_set_nbqw(0)
      end if
     
      end
 
