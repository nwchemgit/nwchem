c $Id: cons_data.F,v 1.3 2004-09-21 19:12:09 marat Exp $
      subroutine cons_set_namespace(namespace)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cons_data.fh"
      character*(*) namespace
c     
      cons_name = namespace
      end
c
c      subroutine cons_get_namespace(namespace)
c      implicit none
c#include "errquit.fh"
c#include "cons_data.fh"
c      character*(*) namespace
cc     
c      namespace=cons_name
c      end
c     
      function cons_get_h_c()
      implicit none
#include "cons_data.fh"
      integer cons_get_h_c
      cons_get_h_c=h_c
      return
      end

      function cons_get_i_c()
      implicit none
#include "mafdecls.fh"
#include "cons_data.fh"
#include "errquit.fh"
      integer cons_get_i_c
c      cons_get_i_c = i_c
      if(.not. ma_get_index(h_c,cons_get_i_c))
     >     call errquit('cons_data: unable to get i_c',
     >      h_c, RTDB_ERR)
      return
      end

      subroutine cons_set_h_c(a_h_c)
      implicit none
#include "cons_data.fh"
      integer a_h_c
      h_c=a_h_c
      return
      end

      function cons_get_h_hbond_id()
      implicit none
#include "cons_data.fh"
      integer cons_get_h_hbond_id
      cons_get_h_hbond_id=h_ijhb
      return
      end

      function cons_get_h_hbond_k()
      implicit none
#include "cons_data.fh"
      integer cons_get_h_hbond_k
      cons_get_h_hbond_k=h_khb
      return
      end

      function cons_get_h_hbond_r()
      implicit none
#include "cons_data.fh"
      integer cons_get_h_hbond_r
      cons_get_h_hbond_r=h_rhb
      return
      end

      function cons_get_i_hbond_id()
      implicit none
#include "mafdecls.fh"
#include "cons_data.fh"
#include "errquit.fh"
      integer cons_get_i_hbond_id
      if(.not. ma_get_index(h_ijhb,cons_get_i_hbond_id)) 
     >    call errquit('cons_data: unable to get i_hbond_id',
     >      h_ijhb, RTDB_ERR)
      return
      end

      function cons_get_i_hbond_k()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "cons_data.fh"
      integer cons_get_i_hbond_k
      if(.not. ma_get_index(h_khb,cons_get_i_hbond_k))
     >     call errquit
     >     ('cons_data: unable to get i_hbond_k',
     >      h_khb, RTDB_ERR)

      return
      end

      function cons_get_i_hbond_r()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "cons_data.fh"
      integer cons_get_i_hbond_r
      if(.not. ma_get_index(h_rhb,cons_get_i_hbond_r)) 
     >     call errquit
     >     ('cons_data: unable to get i_hbond_r',
     >      h_rhb, RTDB_ERR)

      return
      end

      subroutine cons_set_h_hbond_r(h)
      implicit none
#include "cons_data.fh"
      integer h
      h_rhb = h
      return
      end

      subroutine cons_set_h_hbond_k(h)
      implicit none
#include "cons_data.fh"
      integer h
      h_khb = h
      return
      end

      subroutine cons_set_h_hbond_id(h)
      implicit none
#include "cons_data.fh"
      integer h
      h_ijhb = h
      return
      end

      subroutine cons_set_hbond_nhb(n)
      implicit none
#include "cons_data.fh"
      integer n
      nhb=n
      return
      end

      subroutine cons_get_hbond_nhb(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=nhb
      return
      end

      function cons_get_irtdb()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "cons_data.fh"
      integer cons_get_irtdb
      cons_get_irtdb = cons_irtdb

      return
      end

      subroutine cons_set_irtdb(artdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "cons_data.fh"
      integer artdb
      cons_irtdb = artdb

      end

