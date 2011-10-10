      subroutine btrans1(d_a,k_a_offset,ga,atpart2,nalength2)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
      integer d_a,k_a_offset
      integer ga
      integer atpart2
      integer nalength2(*)
c
      integer d_i1,l_i1_offset,k_i1_offset,size_i1
      character*255 filename
c
c     ======================================================================
c     factorized tile structure based code use share files for intermediates
c     ======================================================================
c
      call offset_btrans1_i1(l_i1_offset,k_i1_offset,size_i1,atpart2,
     1     nalength2)
      call tce_filename('btrans1_i1',filename)
      call createfile(filename,d_i1,size_i1)
      call btrans1_i1(d_a,k_a_offset,d_i1,k_i1_offset,atpart2,
     1     nalength2)
      call reconcilefile(d_i1,size_i1)
      call btrans1_i0(d_i1,k_i1_offset,ga,atpart2,nalength2)
      call deletefile(d_i1)
      if (.not. ma_pop_stack(l_i1_offset))
     1    call errquit('btrans1: ma problem',1,ma_err)
c
      end
c
      subroutine btrans1_i1(d_a,k_a_offset,d_c,k_c_offset,atpart2,
     1           nalength2)
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
      integer d_a, k_a_offset
      integer d_c, k_c_offset
      integer atpart2,nalength2(*)
c
      integer dima,dima_sort
      integer dimb,dimb_sort
      integer dimc,dimc_sort
      integer dim_common
      integer k_a,l_a
      integer k_b,l_b
      integer k_c,l_c
      integer g1b,g2b
      integer g1,g2
      integer jb
      integer j,k,j_ind,j1
c
      integer nxtask
      integer next
      integer nprocs
      integer count
      external nxtask
c
      nprocs=ga_nnodes()
      count=0
      next=nxtask(nprocs,1)
c
      j_ind=0
      do jb=1,atpart2
      do g1b=1,noab+nvab
         if (next.eq.count) then
            dimc=int_mb(k_range+g1b-1)*nalength2(jb)
            if (.not.ma_push_get(mt_dbl,dimc,'noname',l_c,k_c))
     1          call errquit('btrans1_i1: ma problem',1,ma_err)
            call dfill(dimc,0.d0,dbl_mb(k_c),1)
            do g2b=1,noab+nvab
               if ((.not.restricted).or.(int_mb(k_spin+g2b-1)
     1              +int_mb(k_spin+g1b-1)).ne.4) then
               if (int_mb(k_spin+g2b-1).eq.int_mb(k_spin+g1b-1)) then
               if (ieor(int_mb(k_sym+g2b-1),int_mb(k_sym+g1b-1)).eq.
     1             0) then
                  dima=int_mb(k_range+g1b-1)*int_mb(k_range+g2b-1)
                  dima_sort = int_mb(k_range+g1b-1)
                  dim_common=int_mb(k_range+g2b-1)
                  if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1            call errquit('btrans1_i1: ma problem',2,ma_err)
                  call get_hash_block(d_a,dbl_mb(k_a),dima,
     1                 int_mb(k_a_offset),(g2b-1)+(g1b-1)*(noab+nvab))
                  dimb=int_mb(k_range+g2b-1)*nalength2(jb)
                  dimb_sort=nalength2(jb)
                  if(.not.ma_push_get(mt_dbl,dimb,'noname',l_b,k_b))
     1            call errquit('btrans1_i1: ma problem',3,ma_err)
                  do j=1,nalength2(jb)
                     j1=j_ind+j
                    do g2=1,int_mb(k_range+g2b-1)
                      dbl_mb(k_b+(g2-1)+(j-1)*int_mb(k_range+g2b-1))=
     1                dbl_mb(k_movecs_sorted+
     1                (j1-1)+(int_mb(k_offset+g2b-1)+g2-1)*nbf)
                    enddo
                  enddo
                  call dgemm('t','n',dima_sort,dimb_sort,dim_common,
     1                      1.d0,dbl_mb(k_a),dim_common,dbl_mb(k_b),
     1                      dim_common,1.d0,dbl_mb(k_c),dima_sort)
                  if (.not.ma_pop_stack(l_b))
     1            call errquit('btrans1_i1:ma problem',4,ma_err)
                  if (.not.ma_pop_stack(l_a))
     1            call errquit('btrans1_i1: ma problem',5,ma_err)
               endif
               endif
               endif
            enddo
            call put_hash_block(d_c,dbl_mb(k_c),dimc,
     1           int_mb(k_c_offset),
     1           (g1b-1)+(jb-1)*(noab+nvab))
            if(.not.ma_pop_stack(l_c))
     1      call errquit('btrans1_i1: ma problem',6,ma_err)
            next = nxtask(nprocs,1)
         endif
         count=count+1
      enddo
      j_ind=j_ind+nalength2(jb)
      enddo
      next = nxtask(-nprocs,1)
      call ga_sync( )
      return
      end
c
      subroutine btrans1_i0(d_a,k_a_offset,ga,atpart2,nalength2)
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
      integer d_a,k_a_offset
      integer k_1pdm_ao
      integer ga
      integer atpart2
      integer nalength2(*)
c
      integer k_b,l_b,dimb,dimb_sort
      integer k_a,l_a,dima,dima_sort
      integer k_c,l_c,dimc,dimc_sort
      integer dim_common
c
      integer ib,jb
      integer i,j,i1
      integer ihi,ilo
      integer jhi,jlo
      double precision res
      integer i_ind, j_ind
      integer g1b,g1
c
      integer nxtask
      integer next
      integer nprocs
      integer count
      external nxtask
c
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      i_ind=0
      do ib=1,atpart2
      j_ind=0
      do jb=1,atpart2
         if (next.eq.count) then
            dimc=nalength2(ib)*nalength2(jb)
            if (.not.ma_push_get(mt_dbl,dimc,'noname',l_c,k_c))
     1      call errquit('btrans1_i0: ma problem',1,ma_err)
            call dfill(dimc,0.d0,dbl_mb(k_c),1)
            do g1b=1,noab+nvab
               dima=int_mb(k_range+g1b-1)*nalength2(jb)
               dima_sort = nalength2(jb)
               if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1         call errquit('btrans1_i0: ma problem',2,ma_err)
               dimb=int_mb(k_range+g1b-1)*nalength2(ib)
               dimb_sort = nalength2(ib)
               dim_common=int_mb(k_range+g1b-1)
               if (.not.ma_push_get(mt_dbl,dimb,'noname',l_b,k_b))
     1            call errquit('btrans1_i0: ma problem',3,ma_err)
               call get_hash_block(d_a,dbl_mb(k_a),dima,
     1              int_mb(k_a_offset),
     1              (g1b-1)+(jb-1)*(noab+nvab))
               do i=1,nalength2(ib)
                  i1=i_ind+i
               do g1=1,int_mb(k_range+g1b-1)
                  dbl_mb(k_b+(g1-1)+(i-1)*int_mb(k_range+g1b-1))=
     1            dbl_mb(k_movecs_sorted+
     1            (i1-1)+(int_mb(k_offset+g1b-1)+g1-1)*nbf)
               enddo
               enddo
               call dgemm('t','n',dima_sort,dimb_sort,dim_common,1.d0,
     1              dbl_mb(k_a),dim_common,
     1              dbl_mb(k_b),dim_common,1.d0,dbl_mb(k_c),dima_sort)
               if(.not.ma_pop_stack(l_b))
     1          call errquit('btrans1_i0: ma problem',4,ma_err)
               if(.not.ma_pop_stack(l_a))
     1         call errquit('btrans1_hh_i0: ma problem',5,ma_err)
            enddo
            if (restricted) then
               do i=1,dimc
                  dbl_mb(k_c+i-1)=2.d0*dbl_mb(k_c+i-1)
               enddo
            endif
            ilo=i_ind+1
            ihi=i_ind+nalength2(ib)
            jlo=j_ind+1
            jhi=j_ind+nalength2(jb)
            call ga_put(ga,jlo,jhi,ilo,ihi,dbl_mb(k_c),nalength2(jb))
            if (.not.ma_pop_stack(l_c))
     1      call errquit('btrans1_hh_i0: ma problem',6,ma_err)
            next = nxtask(nprocs,1)
         endif
         count = count + 1
         j_ind=j_ind+nalength2(jb)
      enddo
         i_ind=i_ind+nalength2(ib)
      enddo
      next = nxtask(-nprocs,1)
      call ga_sync( )
      return
      end
c
      subroutine offset_btrans1_i1(l_a_offset,k_a_offset,size,
     1           atpart2,nalength2)
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
      integer l_a_offset
      integer k_a_offset
      integer size
      integer atpart2
      integer nalength2(*)
c
      integer length
      integer addr
      integer ib
      integer g1b
c
      length = 0
      do ib=1,atpart2
      do g1b=1,noab+nvab
         length=length+1
      enddo
      enddo
c
      if (.not. ma_push_get(mt_int,2*length+1,'noname',l_a_offset,
     1    k_a_offset))
     1  call errquit('offset_btrans1_i1',0,ma_err)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      do ib=1,atpart2
      do g1b = 1, noab+nvab
         addr = addr + 1
         int_mb(k_a_offset+addr) = (g1b-1)+(ib-1)*(noab+nvab)
         int_mb(k_a_offset+length+addr) = size
         size = size + int_mb(k_range+g1b-1)*nalength2(ib)
      enddo
      enddo
      end
c $Id$
