      subroutine qmmm_bq_data_init(irtdb)
      implicit none
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"
      integer irtdb
c
      nbqs=0
      nbqw=0
      nbq=nbqs+nbqw

       if (.not.rtdb_get(irtdb,'qmmm:bq_exclude',mt_int,1,bq_exclude))
     + call errquit('qmmm_data_init: qmmm:bq_exclude',bq_exclude,
     &       RTDB_ERR)


      end

      subroutine qmmm_bq_data_load()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"
#include "qmmm.fh"
      logical link

      link = qmmm_h_link()

      call qmmm_bq_data_dealloc()
c
      call mm_prune_bqzone()
c 
      call mm_get_tot_nbqs(nbqs,link)
      call mm_get_tot_nbqw(nbqw)
      nbq=nbqs+nbqw

      if(nbq.eq.0) return

      call qmmm_bq_data_alloc()

      if(nbqs.ne.0) then
        call mm_get_solute_bq_ind(nbqs,link,int_mb(i_ibq))
        call mm_get_solute_geom_bq(nbqs,link,
     >                       int_mb(i_ibq),
     >                       dbl_mb(i_cbq),
     >                       dbl_mb(i_qbq))

      end if
c
c     append solvent bq coord and charges from MM
c     -------------------------------------------    
      if(nbqw.ne.0) then
        call mm_get_solvent_ind_bq(nbqw,int_mb(i_ibq+nbqs))
        call mm_get_solvent_geom(nbqw,int_mb(i_ibq+nbqs),
     >                       dbl_mb(i_cbq+3*nbqs),
     >                       dbl_mb(i_qbq+nbqs))
      end if
     
      end

      subroutine qmmm_bq_data_alloc()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"

      if(nbq.eq.0) return

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

      subroutine qmmm_bq_data_dealloc()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"

      if(nbq.eq.0) return

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

      function qmmm_get_nbq()
      implicit none
#include "qmmm_bq_data.fh"

      integer qmmm_get_nbq

      qmmm_get_nbq = nbqs+nbqw

      end 

      function qmmm_get_nbqs()
      implicit none
#include "qmmm_bq_data.fh"

      integer qmmm_get_nbqs

      qmmm_get_nbqs = nbqs

      end 

      function qmmm_get_nbqw()
      implicit none
#include "qmmm_bq_data.fh"

      integer qmmm_get_nbqw

      qmmm_get_nbqw = nbqw

      end 

      function qmmm_get_h_qbq()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_h_qbq

      qmmm_get_h_qbq = h_qbq

      end 

      function qmmm_get_i_qbq()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_i_qbq

      qmmm_get_i_qbq = i_qbq

      end 

      function qmmm_get_h_cbq()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_h_cbq

      qmmm_get_h_cbq = h_cbq

      end 

      function qmmm_get_i_cbq()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_i_cbq

      qmmm_get_i_cbq = i_cbq

      end 

      function qmmm_get_i_ibq()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_i_ibq

      qmmm_get_i_ibq = i_ibq

      end 

      function qmmm_get_i_tbq()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_i_tbq

      qmmm_get_i_tbq = i_tbq

      end 

      function qmmm_get_h_cbqs()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_h_cbqs

      qmmm_get_h_cbqs = h_cbqs

      end 

      function qmmm_get_i_cbqs()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_i_cbqs

      qmmm_get_i_cbqs = i_cbqs

      end 

      function qmmm_get_h_cbqw()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_h_cbqw

      qmmm_get_h_cbqw = h_cbqw

      end 

      function qmmm_get_i_cbqw()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_i_cbqw

      qmmm_get_i_cbqw = i_cbqw

      end 

      function qmmm_get_i_qbqw()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_i_qbqw

      qmmm_get_i_qbqw = i_qbqw

      end 

      function qmmm_get_bq_exclude()
      implicit none
#include "qmmm_bq_data.fh"

      integer qmmm_get_bq_exclude

      qmmm_get_bq_exclude = bq_exclude

      end 


