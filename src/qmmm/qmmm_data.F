      subroutine qmmm_data_init(irtdb)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
#include "rtdb.fh"

      integer irtdb


      call mm_get_tot_nqm(nqm)
      call mm_get_tot_nlink(nlink)
c
c     allocate data arrays
c     --------------------
      call qmmm_data_alloc(nqm)

c 
c     establish persistent mapping between mm and qm atoms
c     ----------------------------------------------------
c 
c     first we do quantum non-link atoms
      call mm_get_solute_quant_ind(nqm-nlink,int_mb(i_ai))
c     and then link atoms if any
c     -------------------------
      if(nlink.gt.0) then
        call mm_get_solute_link_ind(nlink,int_mb(i_ai+nqm-nlink))
      end if


      if(.not.rtdb_get(irtdb,'charge',mt_dbl,1,qcharge)) then
        qcharge = 0.0d0    
      end if

       if (.not.rtdb_get(irtdb,'qmmm:linkatm',mt_int,1,link_atom))
     + call errquit('qmmm_data_init: qmmm:linkatm',link_atom,
     &       RTDB_ERR)

       if (.not.rtdb_get(irtdb,'qmmm:bq_exclude',mt_int,1,bq_exclude))
     + call errquit('qmmm_data_init: qmmm:bq_exclude',bq_exclude,
     &       RTDB_ERR)

      nbqs  = 0
      nbqw  = 0

      end

      subroutine qmmm_data_alloc(nat)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
      integer nat

      nqm = nat
c
c     indexing array
c     -------------
      if(.not.ma_alloc_get(MT_INT, nat, 'qmmm index array',
     &      h_ai, i_ai) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      call ifill(nqm,0,int_mb(i_ai),1)
c
c     coordinate array
c     -------------
      if(.not.ma_alloc_get(MT_DBL, 3*nat, 'qmmm coord array',
     &      h_c, i_c) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      call dfill(3*nat,0,dbl_mb(i_c),1)

c
c     charge array
c     -------------
      if(.not.ma_alloc_get(MT_DBL, nat, 'qmmm charge array',
     &      h_q, i_q) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      call dfill(nat,0,dbl_mb(i_q),1)

c
c     mass array
c     -------------
      if(.not.ma_alloc_get(MT_DBL, nat, 'qmmm mass array',
     &      h_m, i_m) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      call dfill(nat,0,dbl_mb(i_m),1)

c
c     tag array
c     -------------
      if(.not.ma_alloc_get(MT_BYTE, 16*nat, 'qmmm tag array',
     &      h_t, i_t) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)

c
c     atomic number array
c     -------------------
      if(.not.ma_alloc_get(MT_INT, nat, 'qmmm inum array',
     &      h_num, i_num) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      call ifill(nqm,0,int_mb(i_num),1)
    
       
      end

      subroutine qmmm_data_bq_alloc(nbq)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
      integer nbq

      if(.not.ma_alloc_get(mt_dbl,3*nbq,'cbq',h_cbq,i_cbq))
     + call errquit('qmmm: Failed to allocate memory for cbq',3*nbq,
     +       MA_ERR)
      call dfill(3*nbq,0.00,dbl_mb(i_cbq),1)
      if(.not.ma_alloc_get(mt_dbl,nbq,'qbq',h_qbq,i_qbq))
     + call errquit('qmmm: Failed to allocate memory for qbq',nbq,
     &       MA_ERR)
      call dfill(nbq,0.00,dbl_mb(i_qbq),1)
      if(.not.ma_alloc_get(mt_int,nbq,'ibq',h_ibq,i_ibq))
     + call errquit('qmmm: Failed to allocate memory for ibq',nbq,
     &       MA_ERR)
      call ifill(nbq,0,int_mb(i_ibq),1)
      if(.not.ma_alloc_get(MT_BYTE, 16*nbq, 'bq tag array',
     &      h_tbq, i_tbq) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nbq, MA_ERR)

      
      end

      subroutine qmmm_data_bq_dealloc()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"

      integer nbq
      nbq=nbqs+nbqw
      if(nbqs+nbqw.eq.0) return

        if(.not.ma_free_heap(h_tbq))
     +   call errquit('qmmm:memory deallocation for qbq',nbq,
     &         MA_ERR)

        if(.not.ma_free_heap(h_ibq))
     +   call errquit('qmmm:memory deallocation for qbq',nbq,
     &         MA_ERR)
      
        if(.not.ma_free_heap(h_qbq))
     +   call errquit('qmmm:memory deallocation for qbq',nbq,
     &         MA_ERR)

        if(.not.ma_free_heap(h_cbq))
     +   call errquit('qmmm: memory deallocation for cbq',3*nbq,
     +         MA_ERR)

        nbqs=0
        nbqw=0

      end
 
      subroutine qmmm_data_release()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
 
     
      if(.not.ma_free_heap(h_c))
     & call errquit('qmmm c: Failed to deallocate heap',nqm,
     &       MA_ERR)

      if(.not.ma_free_heap(h_ai))
     & call errquit('qmmm ai: Failed to deallocate heap',nqm,
     &       MA_ERR)

       return
      end

      function qmmm_get_nqm()
      implicit none
#include "qmmm_data.fh"

      integer qmmm_get_nqm

      qmmm_get_nqm = nqm

      end 

      function qmmm_get_nlink()
      implicit none
#include "qmmm_data.fh"

      integer qmmm_get_nlink

      qmmm_get_nlink = nlink

      end 

      function qmmm_get_ntot()
      implicit none
#include "qmmm_data.fh"

      integer qmmm_get_ntot

      qmmm_get_ntot = nqm+nbqs+nbqw

      end 

      function qmmm_get_nbq()
      implicit none
#include "qmmm_data.fh"

      integer qmmm_get_nbq

      qmmm_get_nbq = nbqs+nbqw

      end 

      function qmmm_get_nbqs()
      implicit none
#include "qmmm_data.fh"

      integer qmmm_get_nbqs

      qmmm_get_nbqs = nbqs

      end 

      function qmmm_get_nbqw()
      implicit none
#include "qmmm_data.fh"

      integer qmmm_get_nbqw

      qmmm_get_nbqw = nbqw

      end 

      subroutine qmmm_set_nbqw(anbqw)
      implicit none
#include "qmmm_data.fh"

      integer anbqw

      nbqw = anbqw

      end 

      subroutine qmmm_set_nbqs(anbqs)
      implicit none
#include "qmmm_data.fh"

      integer anbqs

      nbqs = anbqs

      end 

      function qmmm_get_i_ai()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_ai

      qmmm_get_i_ai = i_ai

      end 

      function qmmm_get_i_c()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_c

      qmmm_get_i_c = i_c

      end 

      function qmmm_get_h_qbq()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_h_qbq

      qmmm_get_h_qbq = h_qbq

      end 

      function qmmm_get_i_qbq()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_qbq

      qmmm_get_i_qbq = i_qbq

      end 

      function qmmm_get_h_cbq()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_h_cbq

      qmmm_get_h_cbq = h_cbq

      end 

      function qmmm_get_i_cbq()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_cbq

      qmmm_get_i_cbq = i_cbq

      end 

      function qmmm_get_i_ibq()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_ibq

      qmmm_get_i_ibq = i_ibq

      end 

      function qmmm_get_i_tbq()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_tbq

      qmmm_get_i_tbq = i_tbq

      end 

      function qmmm_get_h_cbqs()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_h_cbqs

      qmmm_get_h_cbqs = h_cbqs

      end 

      function qmmm_get_i_cbqs()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_cbqs

      qmmm_get_i_cbqs = i_cbqs

      end 

      function qmmm_get_h_cbqw()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_h_cbqw

      qmmm_get_h_cbqw = h_cbqw

      end 

      function qmmm_get_i_cbqw()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_cbqw

      qmmm_get_i_cbqw = i_cbqw

      end 

      subroutine qmmm_set_cbq(a_i_cbq,a_h_cbq)
      implicit none
#include "qmmm_data.fh"
      integer a_i_cbq,a_h_cbq

      i_cbq=a_i_cbq
      h_cbq=a_h_cbq

      end 
c
      subroutine qmmm_set_ibq(a_i_ibq,a_h_ibq)
      implicit none
#include "qmmm_data.fh"
      integer a_i_ibq,a_h_ibq

      i_ibq=a_i_ibq
      h_ibq=a_h_ibq

      end 

      subroutine qmmm_set_qbq(a_i_qbq,a_h_qbq)
      implicit none
#include "qmmm_data.fh"
      integer a_i_qbq,a_h_qbq

      i_qbq=a_i_qbq
      h_qbq=a_h_qbq

      end 

      function qmmm_get_i_t()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_t

      qmmm_get_i_t = i_t

      end 

      function qmmm_get_i_num()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_num

      qmmm_get_i_num = i_num

      end 

      function qmmm_get_i_m()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_m

      qmmm_get_i_m = i_m

      end 

      function qmmm_get_i_q()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_q

      qmmm_get_i_q = i_q

      end 

      function qmmm_get_i_qbqw()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_qbqw

      qmmm_get_i_qbqw = i_qbqw

      end 

      function qmmm_get_i_g()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_g

      qmmm_get_i_g = i_g

      end 

      function qmmm_get_i_gmm()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_gmm

      qmmm_get_i_gmm = i_gmm

      end 

      function qmmm_get_qcharge()
      implicit none
#include "qmmm_data.fh"
      double precision qmmm_get_qcharge

      qmmm_get_qcharge = qcharge

      end 

      subroutine qmmm_set_qcharge(aqcharge)
      implicit none
#include "qmmm_data.fh"
      double precision aqcharge

      qcharge = aqcharge

      end 

      function qmmm_get_bq_exclude()
      implicit none
#include "qmmm_data.fh"

      integer qmmm_get_bq_exclude

      qmmm_get_bq_exclude = bq_exclude

      end 

      function qmmm_get_link_atom()
      implicit none
#include "qmmm_data.fh"

      integer qmmm_get_link_atom

      qmmm_get_link_atom = link_atom

      end 

      subroutine qmmm_set_nlink(anlink)
      implicit none
#include "qmmm_data.fh"

      integer anlink

      nlink=anlink

      end 

