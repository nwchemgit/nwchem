c
c $Id: qmmm_bq.F,v 1.7 2004-09-03 20:07:56 marat Exp $
c
c

      subroutine qmmm_load_bq(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "mm_utils.fh"
      integer irtdb
c     local variables
      integer i
      integer nbq
      integer nbqs
      integer nbqw
      integer i_cbq,h_cbq
      integer i_qbq,h_qbq
      integer i_ibq,h_ibq
      integer i_tbq
      character*4 bqtag
      character*16 myname

c
c     clean up any prior bq
c     information
c     ---------------------
      call qmmm_unload_bq()
c
c     remove bq atoms bonded to link
c     atoms per input specifications
c     -----------------------------
      call mm_prune_bqzone()
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
      if(nbq.eq.0) return

c
c     allocate memory for coord and charges
c     and indexing
c     -------------------------------------
      call qmmm_data_bq_alloc(nbq)

c
c     load solute bq coord and charges from MM
c     ---------------------------------------    
      
      i_cbq = qmmm_get_i_cbq()
      i_qbq = qmmm_get_i_qbq()
      i_ibq = qmmm_get_i_ibq()
      i_tbq = qmmm_get_i_ibq()
      if(nbqs.ne.0) then
c        call mm_get_solute_bq_ind(nbqs,int_mb(i_ibq))
        call mm_get_geom_bqs(nbqs,
     >                       int_mb(i_ibq),
     >                       dbl_mb(i_cbq),
     >                       dbl_mb(i_qbq))
        myname = "test-bqs.pdb"
        bqtag="Bq"
        call qmmm_print_bq_pdb(nbqs,myname,dbl_mb(i_cbq),
     >      dbl_mb(i_qbq),bqtag)

      end if
c
c     append solvent bq coord and charges from MM
c     -------------------------------------------    
      if(nbqw.ne.0) then
        call mm_get_geom_bqw(nbqw,
     >                       int_mb(i_ibq+nbqs),
     >                       dbl_mb(i_cbq+3*nbqs),
     >                       dbl_mb(i_qbq+nbqs))
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

      call qmmm_data_bq_dealloc()
    
      end
 
