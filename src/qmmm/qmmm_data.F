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
c     formal charge array
c     -------------
      if(.not.ma_alloc_get(MT_DBL, nat, 'qmmm formal charge array',
     &      h_qf, i_qf) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      call dfill(nat,0,dbl_mb(i_qf),1)


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

      function qmmm_get_ntot()
      implicit none
#include "qmmm_data.fh"

      integer qmmm_get_ntot

      qmmm_get_ntot = ntot

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

      subroutine qmmm_set_ntot(antot)
      implicit none
#include "qmmm_data.fh"

      integer antot

      ntot = antot

      end 

      subroutine qmmm_set_nbqs(anbqs)
      implicit none
#include "qmmm_data.fh"

      integer anbqs

      nbqs = anbqs

      end 

      function qmmm_get_i_qf()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_qf

      qmmm_get_i_qf = i_qf

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

      function qmmm_get_h_qbqs()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_h_qbqs

      qmmm_get_h_qbqs = h_qbqs

      end 

      function qmmm_get_i_qbqs()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_qbqs

      qmmm_get_i_qbqs = i_qbqs

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

      subroutine qmmm_set_cbqs(a_i_cbqs,a_h_cbqs)
      implicit none
#include "qmmm_data.fh"
      integer a_i_cbqs,a_h_cbqs

      i_cbqs=a_i_cbqs
      h_cbqs=a_h_cbqs

      end 

      subroutine qmmm_set_qbqs(a_i_qbqs,a_h_qbqs)
      implicit none
#include "qmmm_data.fh"
      integer a_i_qbqs,a_h_qbqs

      i_qbqs=a_i_qbqs
      h_qbqs=a_h_qbqs

      end 

      function qmmm_get_i_cbqw()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_cbqw

      qmmm_get_i_cbqw = i_cbqw

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


