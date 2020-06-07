c $Id$
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

      subroutine cons_set_unita(unita0)
      implicit none
      real*8 unita0(3,3)

#include "cons_data.fh"

*     *** local variables ***
      integer i,j
      real*8 vol

      do j=1,3
      do i=1,3
         unita(i,j) = unita0(i,j)
      end do
      end do
      ub(1,1) = unita(2,2)*unita(3,3) - unita(3,2)*unita(2,3)
      ub(2,1) = unita(3,2)*unita(1,3) - unita(1,2)*unita(3,3)
      ub(3,1) = unita(1,2)*unita(2,3) - unita(2,2)*unita(1,3)
      ub(1,2) = unita(2,3)*unita(3,1) - unita(3,3)*unita(2,1)
      ub(2,2) = unita(3,3)*unita(1,1) - unita(1,3)*unita(3,1)
      ub(3,2) = unita(1,3)*unita(2,1) - unita(2,3)*unita(1,1)
      ub(1,3) = unita(2,1)*unita(3,2) - unita(3,1)*unita(2,2)
      ub(2,3) = unita(3,1)*unita(1,2) - unita(1,1)*unita(3,2)
      ub(3,3) = unita(1,1)*unita(2,2) - unita(2,1)*unita(1,2)
      vol = unita(1,1)*ub(1,1)+unita(2,1)*ub(2,1)+unita(3,1)*ub(3,1)
      call dscal(9,1.0d0/vol,ub,1)

      return
      end

      subroutine cons_lattice_min_diff(x,y,z)
      implicit none
      real*8 x,y,z

#include "cons_data.fh"

*     *** local variables ****
      real*8 c1,c2,c3

      if (isystype.eq.3) then
         c1 = x*ub(1,1) + y*ub(2,1) + z*ub(3,1)
         c2 = x*ub(1,2) + y*ub(2,2) + z*ub(3,2)
         c3 = x*ub(1,3) + y*ub(2,3) + z*ub(3,3)
         c1 = c1 - DNINT(c1)
         c2 = c2 - DNINT(c2)
         c3 = c3 - DNINT(c3)
         x = unita(1,1)*c1 + unita(1,2)*c2 + unita(1,3)*c3
         y = unita(2,1)*c1 + unita(2,2)*c2 + unita(2,3)*c3
         z = unita(3,1)*c1 + unita(3,2)*c2 + unita(3,3)*c3
      end if
      return
      end

      subroutine cons_set_isystype(isystype0)
      implicit none
#include "cons_data.fh"
      integer isystype0
      isystype = isystype0
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

      function cons_get_h_hdihed_id()
      implicit none
#include "cons_data.fh"
      integer cons_get_h_hdihed_id
      cons_get_h_hdihed_id=h_idhd
      return
      end

      function cons_get_h_hdihed_k()
      implicit none
#include "cons_data.fh"
      integer cons_get_h_hdihed_k
      cons_get_h_hdihed_k=h_khd
      return
      end

      function cons_get_h_hdihed_r()
      implicit none
#include "cons_data.fh"
      integer cons_get_h_hdihed_r
      cons_get_h_hdihed_r=h_rhd
      return
      end

      function cons_get_i_hdihed_id()
      implicit none
#include "mafdecls.fh"
#include "cons_data.fh"
#include "errquit.fh"
      integer cons_get_i_hdihed_id
      if(.not. ma_get_index(h_idhd,cons_get_i_hdihed_id)) 
     >    call errquit('cons_data: unable to get i_hdihed_id',
     >      h_idhd, RTDB_ERR)
      return
      end

      function cons_get_i_hdihed_k()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "cons_data.fh"
      integer cons_get_i_hdihed_k
      if(.not. ma_get_index(h_khd,cons_get_i_hdihed_k))
     >     call errquit
     >     ('cons_data: unable to get i_hdihed_k',
     >      h_khd, RTDB_ERR)

      return
      end

      function cons_get_i_hdihed_r()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "cons_data.fh"
      integer cons_get_i_hdihed_r
      if(.not. ma_get_index(h_rhd,cons_get_i_hdihed_r)) 
     >     call errquit
     >     ('cons_data: unable to get i_hdihed_r',
     >      h_rhd, RTDB_ERR)

      return
      end

      subroutine cons_set_h_hdihed_r(h)
      implicit none
#include "cons_data.fh"
      integer h
      h_rhd = h
      return
      end

      subroutine cons_set_h_hdihed_k(h)
      implicit none
#include "cons_data.fh"
      integer h
      h_khd = h
      return
      end

      subroutine cons_set_h_hdihed_id(h)
      implicit none
#include "cons_data.fh"
      integer h
      h_idhd = h
      return
      end

      subroutine cons_set_hdihed_n(n)
      implicit none
#include "cons_data.fh"
      integer n
      nhd=n
      return
      end

      subroutine cons_get_hdihed_n(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=nhd
      return
      end

      integer function cons_get_h_hbondings_n0()
      implicit none
#include "cons_data.fh"
      cons_get_h_hbondings_n0 = h_hbondings_n0
      return
      end

      integer function cons_get_i_hbondings_n0()
      implicit none
#include "mafdecls.fh"
#include "cons_data.fh"
#include "errquit.fh"
      if(.not. ma_get_index(h_hbondings_n0,cons_get_i_hbondings_n0))
     >    call errquit('cons_data: unable to get i_hbondings_n0',
     >      h_hbondings_n0, RTDB_ERR)
      return
      end


      integer function cons_get_h_hbondings_indx()
      implicit none
#include "cons_data.fh"
      cons_get_h_hbondings_indx = h_hbondings_indx
      return
      end

      integer function cons_get_i_hbondings_indx()
      implicit none
#include "mafdecls.fh"
#include "cons_data.fh"
#include "errquit.fh"
      if(.not. ma_get_index(h_hbondings_indx,cons_get_i_hbondings_indx))
     >    call errquit('cons_data: unable to get i_hbondings_indx',
     >      h_hbondings_indx, RTDB_ERR)
      return
      end

      integer function cons_get_h_hbondings_coef()
      implicit none
#include "cons_data.fh"
      cons_get_h_hbondings_coef = h_hbondings_coef
      return
      end

      integer function cons_get_i_hbondings_coef()
      implicit none
#include "mafdecls.fh"
#include "cons_data.fh"
#include "errquit.fh"
      if(.not. ma_get_index(h_hbondings_coef,cons_get_i_hbondings_coef))
     >    call errquit('cons_data: unable to get i_hbondings_coef',
     >      h_hbondings_coef, RTDB_ERR)
      return
      end

      integer function cons_get_h_hbondings_k0()
      implicit none
#include "cons_data.fh"
      cons_get_h_hbondings_k0 = h_hbondings_k0
      return
      end

      integer function cons_get_i_hbondings_k0()
      implicit none
#include "mafdecls.fh"
#include "cons_data.fh"
#include "errquit.fh"
      if(.not. ma_get_index(h_hbondings_k0,cons_get_i_hbondings_k0))
     >    call errquit('cons_data: unable to get i_hbondings_k0',
     >      h_hbondings_k0, RTDB_ERR)
      return
      end

      integer function cons_get_h_hbondings_gamma0()
      implicit none
#include "cons_data.fh"
      cons_get_h_hbondings_gamma0 = h_hbondings_gamma0
      return
      end

      integer function cons_get_i_hbondings_gamma0()
      implicit none
#include "mafdecls.fh"
#include "cons_data.fh"
#include "errquit.fh"
      if(.not. ma_get_index(h_hbondings_gamma0,
     >                      cons_get_i_hbondings_gamma0))
     >    call errquit('cons_data: unable to get i_hbondings_gamma0',
     >      h_hbondings_gamma0, RTDB_ERR)
      return
      end

      subroutine cons_set_hbondings_nhc(n)
      implicit none
#include "cons_data.fh"
      integer n
      nhc=n
      return
      end

      subroutine cons_get_hbondings_nhc(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=nhc
      return
      end

      subroutine cons_set_hbondings_n0(n)
      implicit none
#include "cons_data.fh"
      integer n
      h_hbondings_n0=n
      return
      end

      subroutine cons_get_hbondings_n0(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=h_hbondings_n0
      return
      end

      subroutine cons_set_hbondings_indx(n)
      implicit none
#include "cons_data.fh"
      integer n
      h_hbondings_indx=n
      return
      end

      subroutine cons_get_hbondings_indx(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=h_hbondings_indx
      return
      end


      subroutine cons_set_hbondings_coef(n)
      implicit none
#include "cons_data.fh"
      integer n
      h_hbondings_coef=n
      return
      end

      subroutine cons_get_hbondings_coef(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=h_hbondings_coef
      return
      end

      subroutine cons_set_hbondings_k0(n)
      implicit none
#include "cons_data.fh"
      integer n
      h_hbondings_k0=n
      return
      end

      subroutine cons_get_hbondings_k0(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=h_hbondings_k0
      return
      end


      subroutine cons_set_hbondings_gamma0(n)
      implicit none
#include "cons_data.fh"
      integer n
      h_hbondings_gamma0=n
      return
      end

      subroutine cons_get_hbondings_gamma0(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=h_hbondings_gamma0
      return
      end



      integer function cons_get_h_pbondings_n0()
      implicit none
#include "cons_data.fh"
      cons_get_h_pbondings_n0 = h_pbondings_n0
      return
      end

      integer function cons_get_i_pbondings_n0()
      implicit none
#include "mafdecls.fh"
#include "cons_data.fh"
#include "errquit.fh"
      if(.not. ma_get_index(h_pbondings_n0,cons_get_i_pbondings_n0))
     >    call errquit('cons_data: unable to get i_pbondings_n0',
     >      h_pbondings_n0, RTDB_ERR)
      return
      end


      integer function cons_get_h_pbondings_indx()
      implicit none
#include "cons_data.fh"
      cons_get_h_pbondings_indx = h_pbondings_indx
      return
      end

      integer function cons_get_i_pbondings_indx()
      implicit none
#include "mafdecls.fh"
#include "cons_data.fh"
#include "errquit.fh"
      if(.not. ma_get_index(h_pbondings_indx,cons_get_i_pbondings_indx))
     >    call errquit('cons_data: unable to get i_pbondings_indx',
     >      h_pbondings_indx, RTDB_ERR)
      return
      end

      integer function cons_get_h_pbondings_coef()
      implicit none
#include "cons_data.fh"
      cons_get_h_pbondings_coef = h_pbondings_coef
      return
      end

      integer function cons_get_i_pbondings_coef()
      implicit none
#include "mafdecls.fh"
#include "cons_data.fh"
#include "errquit.fh"
      if(.not. ma_get_index(h_pbondings_coef,cons_get_i_pbondings_coef))
     >    call errquit('cons_data: unable to get i_pbondings_coef',
     >      h_pbondings_coef, RTDB_ERR)
      return
      end

      integer function cons_get_h_pbondings_k0()
      implicit none
#include "cons_data.fh"
      cons_get_h_pbondings_k0 = h_pbondings_k0
      return
      end

      integer function cons_get_i_pbondings_k0()
      implicit none
#include "mafdecls.fh"
#include "cons_data.fh"
#include "errquit.fh"
      if(.not. ma_get_index(h_pbondings_k0,cons_get_i_pbondings_k0))
     >    call errquit('cons_data: unable to get i_pbondings_k0',
     >      h_pbondings_k0, RTDB_ERR)
      return
      end


      integer function cons_get_h_pbondings_gcut0()
      implicit none
#include "cons_data.fh"
      cons_get_h_pbondings_gcut0 = h_pbondings_gcut0
      return
      end

      integer function cons_get_i_pbondings_gcut0()
      implicit none
#include "mafdecls.fh"
#include "cons_data.fh"
#include "errquit.fh"
      if(.not. ma_get_index(h_pbondings_gcut0,
     >                      cons_get_i_pbondings_gcut0))
     >    call errquit('cons_data: unable to get i_pbondings_gcut0',
     >      h_pbondings_gcut0, RTDB_ERR)
      return
      end


      integer function cons_get_h_pbondings_gamma0()
      implicit none
#include "cons_data.fh"
      cons_get_h_pbondings_gamma0 = h_pbondings_gamma0
      return
      end

      integer function cons_get_i_pbondings_gamma0()
      implicit none
#include "mafdecls.fh"
#include "cons_data.fh"
#include "errquit.fh"
      if(.not. ma_get_index(h_pbondings_gamma0,
     >                      cons_get_i_pbondings_gamma0))
     >    call errquit('cons_data: unable to get i_pbondings_gamma0',
     >      h_pbondings_gamma0, RTDB_ERR)
      return
      end


      integer function cons_get_h_pbondings_c0()
      implicit none
#include "cons_data.fh"
      cons_get_h_pbondings_c0 = h_pbondings_c0
      return
      end

      integer function cons_get_i_pbondings_c0()
      implicit none
#include "mafdecls.fh"
#include "cons_data.fh"
#include "errquit.fh"
      if(.not. ma_get_index(h_pbondings_c0,
     >                      cons_get_i_pbondings_c0))
     >    call errquit('cons_data: unable to get i_pbondings_c0',
     >      h_pbondings_c0, RTDB_ERR)
      return
      end




      subroutine cons_set_pbondings_nhp(n)
      implicit none
#include "cons_data.fh"
      integer n
      nhp=n
      return
      end

      subroutine cons_get_pbondings_nhp(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=nhp
      return
      end

      subroutine cons_set_pbondings_n0(n)
      implicit none
#include "cons_data.fh"
      integer n
      h_pbondings_n0=n
      return
      end

      subroutine cons_get_pbondings_n0(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=h_pbondings_n0
      return
      end

      subroutine cons_set_pbondings_indx(n)
      implicit none
#include "cons_data.fh"
      integer n
      h_pbondings_indx=n
      return
      end

      subroutine cons_get_pbondings_indx(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=h_pbondings_indx
      return
      end


      subroutine cons_set_pbondings_coef(n)
      implicit none
#include "cons_data.fh"
      integer n
      h_pbondings_coef=n
      return
      end

      subroutine cons_get_pbondings_coef(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=h_pbondings_coef
      return
      end

      subroutine cons_set_pbondings_k0(n)
      implicit none
#include "cons_data.fh"
      integer n
      h_pbondings_k0=n
      return
      end

      subroutine cons_get_pbondings_k0(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=h_pbondings_k0
      return
      end

      subroutine cons_set_pbondings_gcut0(n)
      implicit none
#include "cons_data.fh"
      integer n
      h_pbondings_gcut0=n
      return
      end

      subroutine cons_get_pbondings_gcut0(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=h_pbondings_gcut0
      return
      end


      subroutine cons_set_pbondings_gamma0(n)
      implicit none
#include "cons_data.fh"
      integer n
      h_pbondings_gamma0=n
      return
      end

      subroutine cons_get_pbondings_gamma0(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=h_pbondings_gamma0
      return
      end


      subroutine cons_set_pbondings_c0(n)
      implicit none
#include "cons_data.fh"
      integer n
      h_pbondings_c0=n
      return
      end

      subroutine cons_get_pbondings_c0(n)
      implicit none
#include "cons_data.fh"
      integer n
      n=h_pbondings_c0
      return
      end





