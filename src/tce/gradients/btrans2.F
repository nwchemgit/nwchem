      subroutine btrans2(d_2pdm,k_2pdm_offset,d_1pdm,k_1pdm_offset,
     1                   d_i0,k_i0_offset,size_i0,
     1                   atpart2,nalength2,bfglo,bfghi)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
      integer d_2pdm,k_2pdm_offset
      integer d_1pdm,k_1pdm_offset
      integer k_2pdm_ao
      integer atpart2
      integer nalength2(*),bfglo(*),bfghi(*)

      integer k_mo,k_2pdm
c
      integer ib,jb,kb,lb
      integer p,q,r,s
c
      integer nh,np
      double precision res
c
      integer k_ispq, l_ispq, dim_ispq
      integer k_ispl, l_ispl, dim_ispl
      integer k_ikpl, l_ikpl, dim_ikpl
      integer k_ikjl, l_ikjl, dim_ikjl
c
      logical one, two, three, four,five,six
      integer l_tmp,k_tmp,dim_tmp
c
      integer k_pqrs,l_pqrs,dim_pqrs
      integer h1b,h2b,h3b,h4b,h1,h2,h3,h4,h2_0,h4_0
      integer p1b,p2b,p3b,p4b
      integer g1b,g2b,g3b,g4b
      integer m
c
      integer d_i1,l_i1_offset,k_i1_offset,size_i1
      integer d_i2,l_i2_offset,k_i2_offset,size_i2
      integer d_i3,l_i3_offset,k_i3_offset,size_i3
      integer d_i0,l_i0_offset,k_i0_offset,size_i0
      integer i,j,k,l,i1,j1,k1,l1,i_ind,j_ind,k_ind,l_ind,ii
      character*255 filename
c
      integer k_a,l_a,k_b,l_b,k_c,l_c,dima,dimb,dim_common
      integer d_cc,d_cca,d_ccb,d_scf,d_scfa,d_scfb
      integer k_1pao_offset,l_1pao_offset,size_1p
      integer dima_sort,dimb_sort,dimc
      integer k_c_sort,l_c_sort,k_c0,l_c0
c
      integer nxtask
      integer nprocs
      integer count
      integer next
      external nxtask
c
c     Gamma(h1,h2,h3,h4) C(h1,i) = I(i,h2,h3,h4)
c     I(i,h2,h3,h4) C(h4,l) = I(i,h2,h3,l)
c     I(i,h2,h3,l) C(h2,j) = I(i,j,h3,l)
c     I(i,j,h3,l) C(h3,k) = I(i,j,k,l)
c
c      do i=1,atpart2
c         print *,i,nalength2(i)
c      enddo
c
c      print *,'GET INTO BTRANS2'
c      call offset_btrans2_i0(l_i0_offset,k_i0_offset,size_i0,
c     1     atpart2,nalength2)
c      call tce_filename('btrans2_i0',filename)
c      call createfile(filename,d_i0,size_i0)
c      print *,'offset_btrans2_i0 is fine!'
c     =======================================================
c     Contribution from the product of SCF 1-PDM and CC 1-PDM
c     =======================================================
c
c     1. Calculate the alpha and beta components of SCF 1-PDM
c
      call offset_btrans1_i0(l_1pao_offset,k_1pao_offset,
     1     size_1p,atpart2,nalength2)
      call tce_filename('scf_alpha',filename)
      call createfile(filename,d_scfa,size_1p)
      call tce_filename('scf_beta',filename)
      call createfile(filename,d_scfb,size_1p)
      call tce_filename('scf_total',filename)
      call createfile(filename,d_scf,size_1p)
      call tce_filename('cc_alpah',filename)
      call createfile(filename,d_cca,size_1p)
      call tce_filename('cc_beta',filename)
      call createfile(filename,d_ccb,size_1p)
      call tce_filename('cc_total',filename)
      call createfile(filename,d_cc,size_1p)
      call scf_1pao(d_scfa,d_scfb,d_scf,k_1pao_offset,size_1p,
     1     atpart2,nalength2)
      call cc_1pao(d_1pdm,k_1pdm_offset,d_cca,d_ccb,d_cc,
     1     k_1pao_offset,size_1p,atpart2,nalength2)
c      print *,'D_SCFA'
c      call print_1pdm_ao(d_scfa,k_1pao_offset,atpart2,
c     1           bfglo,bfghi)
c      print *,'D_SCFB'
c      call print_1pdm_ao(d_scfb,k_1pao_offset,atpart2,
c     1           bfglo,bfghi)
c      print *,'D_SCF'
c      call print_1pdm_ao(d_scf,k_1pao_offset,atpart2,
c     1           bfglo,bfghi)
c      print *,'D_CCA'
c      call print_1pdm_ao(d_cca,k_1pao_offset,atpart2,
c     1           bfglo,bfghi)
c      print *,'D_CCB'
c      call print_1pdm_ao(d_ccb,k_1pao_offset,atpart2,
c     1           bfglo,bfghi)
c      print *,'D_CC'
c      call print_1pdm_ao(d_cc,k_1pao_offset,atpart2,
c     1           bfglo,bfghi)
c       call ga_sync( )
c
c
c     2. calculate 2-PDM by the products of SCF and CC 1-PDM
c     2.1 (alpha+beta)x(alpha+beta)
c
      nprocs=ga_nnodes()
      count=0
      next=nxtask(nprocs,1)
c
c      i_ind=0
      do ib=1,atpart2
c      j_ind=0
      do jb=1,atpart2
         if(next.eq.count) then
           i_ind=0
           do i=1,ib-1
              i_ind=i_ind+nalength2(i)
           enddo
           j_ind=0
           do j=1,jb-1
              j_ind=j_ind+nalength2(j)
           enddo
c
         dima=nalength2(ib)*nalength2(jb)
         dima_sort=dima
         if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1   call errquit('btrans2: ma problem',3,ma_err)
         call get_hash_block(d_scf,dbl_mb(k_a),dima,
     1        int_mb(k_1pao_offset),(jb-1)+(ib-1)*atpart2)
c         do i=1,nalength2(ib)
c            i1=i_ind+i
c         do j=1,nalength2(jb)
c            j1=j_ind+j
c            write(6,'(a,i5,i5,f20.16)') 'SCF 1PDM_AO=',i1,j1,
c     1      dbl_mb(k_a+j-1+(i-1)*nalength2(jb))
c         enddo
c         enddo
c
         k_ind=0
         do kb=1,atpart2
         l_ind=0
         do lb=1,atpart2
c            if(next.eq.count) then
c
            dimb=nalength2(kb)*nalength2(lb)
            dimb_sort=dimb
            if (.not.ma_push_get(mt_dbl,dimb,'noname',l_b,k_b))
     1      call errquit('btrans2: ma problem',3,ma_err)
            call get_hash_block(d_cc,dbl_mb(k_b),dimb,
     1        int_mb(k_1pao_offset),(lb-1)+(kb-1)*atpart2)
c            do k=1,nalength2(kb)
c               k1=k_ind+k
c            do l=1,nalength2(lb)
c               l1=l_ind+l
c               write(6,'(a,i5,i5,f20.16)') 'CCSD 1PDM_AO=',k1,l1,
c     1         dbl_mb(k_b+l-1+(k-1)*nalength2(lb))
c            enddo
c            enddo
            dimc=dima*dimb
            if (.not.ma_push_get(mt_dbl,dimc,'noname',l_c,k_c))
     1      call errquit('btrans2: ma problem',4,ma_err)
            dim_common=1
            call dgemm('t','n',dima_sort,dimb_sort,dim_common,1.d0,
     1           dbl_mb(k_a),dim_common,
     1           dbl_mb(k_b),dim_common,0.d0,
     1           dbl_mb(k_c),dima_sort)
            if (.not.ma_push_get(mt_dbl,dimc,'noname',
     1           l_c_sort,k_c_sort))
     1      call errquit('btrans2: ma problem',4,ma_err)
c
            ii=0
            do i=1,nalength2(ib)
               i1=i_ind+i
            do j=1,nalength2(jb)
               j1=j_ind+j
            do k=1,nalength2(kb)
               k1=k_ind+k
            do l=1,nalength2(lb)
               l1=l_ind+l
               ii=ii+1
               dbl_mb(k_c_sort+(l-1)+(j-1)*nalength2(lb)
     1         +(k-1)*nalength2(lb)*nalength2(jb)
     1         +(i-1)*nalength2(lb)*nalength2(jb)*nalength2(kb))=
     1          dbl_mb(k_c+(j-1)+(i-1)*nalength2(jb)
     1         +(l-1)*nalength2(ib)*nalength2(jb)
     1         +(k-1)*nalength2(ib)*nalength2(jb)*nalength2(lb))
c               write(6,'(a,4i5,f20.16,f20.16)') 'I. 2PDM AO=',
c     1          i1,j1,k1,l1, 
c     1         dbl_mb(k_c_sort+(l-1)+(j-1)*nalength2(lb)
c     1         +(k-1)*nalength2(lb)*nalength2(jb)
c     1         +(i-1)*nalength2(jb)*nalength2(kb)*nalength2(lb)),
c     1         dbl_mb(k_c+ii-1)          
            enddo
            enddo
            enddo
            enddo
c               
            call put_hash_block(d_i0,dbl_mb(k_c_sort),dima*dimb,
     1           int_mb(k_i0_offset),(lb-1)+(jb-1)*atpart2+
     1           (kb-1)*atpart2*atpart2+
     1           (ib-1)*atpart2*atpart2*atpart2)
c
            if(.not.ma_pop_stack(l_c_sort))
     1      call errquit('btrans2: ma problem',101,ma_err)
c
            if (.not.ma_pop_stack(l_c))
     1      call errquit('btrans2: ma problem',101,ma_err)
            if (.not.ma_pop_stack(l_b))
     1      call errquit('btrans2: ma problem',102,ma_err)
c            next=nxtask(nprocs,1)
c            endif
c            count=count+1
         l_ind=l_ind+nalength2(lb)
         enddo
         k_ind=k_ind+nalength2(kb)
         enddo
         if (.not.ma_pop_stack(l_a))
     1   call errquit('btrans2: ma problem',103,ma_err)
         next=nxtask(nprocs,1)
         endif
         count=count+1
c      j_ind=j_ind+nalength2(jb)
      enddo
c      i_ind=i_ind+nalength2(ib)
      enddo
      next=nxtask(-nprocs,1)
      call ga_sync( )
c      print *,'(ALPHA+BETA)*(ALPHA+BETA)'
c      call print_2pdm_ao(d_i0,k_i0_offset,atpart2,
c     1     bfglo,bfghi)
      
c      print *,'2.1 is fine!'
c
c     2.2 alpha * alpha
      nprocs=ga_nnodes()
      count=0
      next=nxtask(nprocs,1)
c
c      i_ind=0
      do ib=1,atpart2
c      j_ind=0
      do jb=1,atpart2
         if(next.eq.count) then
           i_ind=0
           do i=1,ib-1
              i_ind=i_ind+nalength2(i)
           enddo
           j_ind=0
           do j=1,jb-1
              j_ind=j_ind+nalength2(j)
           enddo
c
         dima=nalength2(ib)*nalength2(jb)
         dima_sort=dima
         if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1   call errquit('btrans2: ma problem',3,ma_err)
         call get_hash_block(d_scfa,dbl_mb(k_a),dima,
     1        int_mb(k_1pao_offset),(jb-1)+(ib-1)*atpart2)
         k_ind=0
         do kb=1,atpart2
         l_ind=0
         do lb=1,atpart2
            dimb=nalength2(kb)*nalength2(lb)
            dimb_sort=dimb
            dimc=dima*dimb
            if (.not.ma_push_get(mt_dbl,dimb,'noname',l_b,k_b))
     1      call errquit('btrans2: ma problem',3,ma_err)
            call get_hash_block(d_cca,dbl_mb(k_b),dimb,
     1        int_mb(k_1pao_offset),(lb-1)+(kb-1)*atpart2)
            if (.not.ma_push_get(mt_dbl,dimc,'noname',l_c,k_c))
     1      call errquit('btrans2: ma problem',4,ma_err)
            if (.not.ma_push_get(mt_dbl,dimc,'noname',
     1           l_c_sort,k_c_sort))
     1      call errquit('btrans2: ma problem',4,ma_err)
            if (.not.ma_push_get(mt_dbl,dimc,'noname',
     1           l_c0,k_c0))
     1      call errquit('btrans2: ma problem',5,ma_err)
            call get_hash_block(d_i0,dbl_mb(k_c0),dimc,
     1           int_mb(k_i0_offset),(lb-1)+(kb-1)*atpart2+
     1           (jb-1)*atpart2*atpart2+
     1           (ib-1)*atpart2*atpart2*atpart2)
c
            dim_common=1
            call dgemm('t','n',dima_sort,dimb_sort,dim_common,1.d0,
     1           dbl_mb(k_a),dim_common,
     1           dbl_mb(k_b),dim_common,0.d0,
     1           dbl_mb(k_c),dima_sort)
c
            do i=1,dima*dimb
               dbl_mb(k_c+i-1)=-1.d0*dbl_mb(k_c+i-1)
            enddo
c           
            do i=1,nalength2(ib)
               i1=i_ind+i
            do j=1,nalength2(jb)
               j1=j_ind+j
            do k=1,nalength2(kb)
               k1=k_ind+k
            do l=1,nalength2(lb)
               l1=l_ind+l
               dbl_mb(k_c_sort+(l-1)+(k-1)*nalength2(lb)
     1         +(j-1)*nalength2(lb)*nalength2(kb)
     1         +(i-1)*nalength2(lb)*nalength2(jb)*nalength2(kb))=
     1          dbl_mb(k_c+(j-1)+(i-1)*nalength2(jb)
     1         +(l-1)*nalength2(ib)*nalength2(jb)
     1         +(k-1)*nalength2(ib)*nalength2(jb)*nalength2(lb))
            enddo
            enddo
            enddo
            enddo
c
            do i=1,dimc
               dbl_mb(k_c_sort+i-1)=dbl_mb(k_c_sort+i-1)
     1         +dbl_mb(k_c0+i-1)
            enddo
c
            call put_hash_block(d_i0,dbl_mb(k_c_sort),dima*dimb,
     1           int_mb(k_i0_offset),(lb-1)+(kb-1)*atpart2+
     1           (jb-1)*atpart2*atpart2+
     1           (ib-1)*atpart2*atpart2*atpart2)
c
            if(.not.ma_pop_stack(l_c0))
     1      call errquit('btrans2: ma problem',100,ma_err)
c
            if(.not.ma_pop_stack(l_c_sort))
     1      call errquit('btrans2: ma problem',101,ma_err)
c
            if (.not.ma_pop_stack(l_c))
     1      call errquit('btrans2: ma problem',104,ma_err)
            if (.not.ma_pop_stack(l_b))
     1      call errquit('btrans2: ma problem',105,ma_err)
c            next=nxtask(nprocs,1)
c            endif
c            count=count+1
         l_ind=l_ind+nalength2(lb)
         enddo
         k_ind=k_ind+nalength2(kb)
         enddo
         if (.not.ma_pop_stack(l_a))
     1   call errquit('btrans2: ma problem',106,ma_err)
         next=nxtask(nprocs,1)
         endif
         count=count+1
c      j_ind=j_ind+nalength2(jb)
      enddo
c      i_ind=i_ind+nalength2(ib)
      enddo
      next=nxtask(-nprocs,1)
      call ga_sync( )
c      print *,'2.2 is fine!'
c
c     2.3 beta * beta
      nprocs=ga_nnodes( )
      count=0
      next=nxtask(nprocs,1)
c
c      call tce_zero(d_i0,size_i0)
      do ib=1,atpart2
      do jb=1,atpart2
         if(next.eq.count) then
c
         dima=nalength2(ib)*nalength2(jb)
         dima_sort=dima
         if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1   call errquit('btrans2: ma problem',3,ma_err)
         call get_hash_block(d_scfb,dbl_mb(k_a),dima,
     1        int_mb(k_1pao_offset),(jb-1)+(ib-1)*atpart2)
         do kb=1,atpart2
         do lb=1,atpart2
c            if(next.eq.count) then
c
            dimb=nalength2(kb)*nalength2(lb)
            dimb_sort=dimb
            dimc=dima*dimb
            if (.not.ma_push_get(mt_dbl,dimb,'noname',l_b,k_b))
     1      call errquit('btrans2: ma problem',3,ma_err)
            call get_hash_block(d_ccb,dbl_mb(k_b),dimb,
     1        int_mb(k_1pao_offset),(lb-1)+(kb-1)*atpart2)
            if (.not.ma_push_get(mt_dbl,dimc,'noname',l_c,k_c))
     1      call errquit('btrans2: ma problem',4,ma_err)
            if (.not.ma_push_get(mt_dbl,dimc,'noname',
     1           l_c_sort,k_c_sort))
     1      call errquit('btrans2: ma problem',4,ma_err)
c
            if (.not.ma_push_get(mt_dbl,dimc,'noname',
     1           l_c0,k_c0))
     1      call errquit('btrans2: ma problem',5,ma_err)
c
            call get_hash_block(d_i0,dbl_mb(k_c0),dimc,
     1           int_mb(k_i0_offset),(lb-1)+(kb-1)*atpart2+
     1           (jb-1)*atpart2*atpart2+
     1           (ib-1)*atpart2*atpart2*atpart2)
c
            dim_common=1
            call dgemm('t','n',dima_sort,dimb_sort,dim_common,1.d0,
     1           dbl_mb(k_a),dim_common,
     1           dbl_mb(k_b),dim_common,0.d0,
     1           dbl_mb(k_c),dima_sort)
c
            do i=1,dima*dimb
               dbl_mb(k_c+i-1)=-1.d0*dbl_mb(k_c+i-1)
            enddo
c
            do i=1,nalength2(ib)
            do j=1,nalength2(jb)
            do k=1,nalength2(kb)
            do l=1,nalength2(lb)
               dbl_mb(k_c_sort+(l-1)+(k-1)*nalength2(lb)
     1         +(j-1)*nalength2(lb)*nalength2(kb)
     1         +(i-1)*nalength2(lb)*nalength2(jb)*nalength2(kb))=
     1          dbl_mb(k_c+(j-1)+(i-1)*nalength2(jb)
     1         +(l-1)*nalength2(ib)*nalength2(jb)
     1         +(k-1)*nalength2(ib)*nalength2(jb)*nalength2(lb))
            enddo
            enddo
            enddo
            enddo
c
            do i=1,dimc
               dbl_mb(k_c_sort+i-1)=dbl_mb(k_c_sort+i-1)+
     1         dbl_mb(k_c0+i-1)
            enddo
c
            call put_hash_block(d_i0,dbl_mb(k_c_sort),dima*dimb,
     1           int_mb(k_i0_offset),(lb-1)+(kb-1)*atpart2+
     1           (jb-1)*atpart2*atpart2+
     1           (ib-1)*atpart2*atpart2*atpart2)
c
            if(.not.ma_pop_stack(l_c0))
     1      call errquit('btrans2: ma problem',100,ma_err)
            if(.not.ma_pop_stack(l_c_sort))
     1      call errquit('btrans2: ma problem',101,ma_err)
            if (.not.ma_pop_stack(l_c))
     1      call errquit('btrans2: ma problem',107,ma_err)
            if (.not.ma_pop_stack(l_b))
     1      call errquit('btrans2: ma problem',108,ma_err)
c            next=nxtask(nprocs,1)
c            endif
c            count=count+1
         enddo
         enddo
         if (.not.ma_pop_stack(l_a))
     1   call errquit('btrans2: ma problem',109,ma_err)
         next=nxtask(nprocs,1)
         endif
         count=count+1
      enddo
      enddo
      next=nxtask(-nprocs,1)
      call ga_sync( )
c      print *,'2.3 is fine!'
c
      call deletefile(d_cc)
      call deletefile(d_ccb)
      call deletefile(d_cca)
      call deletefile(d_scf)
      call deletefile(d_scfb)
      call deletefile(d_scfa)
      if (.not.ma_pop_stack(l_1pao_offset))
     1   call errquit('btrans2: ma problem',1,ma_err)
c
c      print *,'1 PDM products is fine!'
c
c     =======================================================
c      call offset_btrans2_i0(l_i0_offset,k_i0_offset,size_i0,
c     1     atpart2,nalength2)
c      call tce_filename('btrans2_i0',filename)
c      call createfile(filename,d_i0,size_i0)
c
      call offset_btrans2_i1(l_i1_offset,k_i1_offset,size_i1,
     1     atpart2,nalength2)
      call tce_filename('btrans2_i1',filename)
      call createfile(filename,d_i1,size_i1)
c
      call offset_btrans2_i2(l_i2_offset,k_i2_offset,size_i2,
     1     atpart2,nalength2)
      call tce_filename('btrans2_i2',filename)
      call createfile(filename,d_i2,size_i2)
c
      call offset_btrans2_i3(l_i3_offset,k_i3_offset,size_i3,
     1     atpart2,nalength2)
      call tce_filename('btrans2_i3',filename)
      call createfile(filename,d_i3,size_i3)
      call ga_sync( )
c     =======
c     H-H-H-H
c     =======
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     1,noa,1,noa,
     1     1,noa,1,noa,
c
     1     1,noa,noa+1,noab,
     1     1,noa,noa+1,noab,
c
     1     noa+1,noab,1,noa,
     1     noa+1,noab,1,noa,
c
     1     noa+1,noab,noa+1,noab,
     1     noa+1,noab,noa+1,noab)
c      call ga_sync( )
c
c
c      print *,'H-H-H-H is fine!'
c     =======
c     H-H-H-P
c     =======
c      two=.true.
c      if(two) then
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     1,noa,1,noa,
     1     1,noa,noab+1,noab+nva,
c
     1     1,noa,noa+1,noab,
     1     1,noa,noab+nva+1,noab+nvab,
c
     1     noa+1,noab,1,noa,
     1     noa+1,noab,noab+1,noab+nva,
c
     1     noa+1,noab,noa+1,noab,
     1     noa+1,noab,noab+nva+1,noab+nvab)
c      call ga_sync( )
c
c      print *,'H-H-H-P is fine!'
c
c     =======
c     H-H-P-H
c     =======
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     1,noa,1,noa,
     1     noab+1,noab+nva,1,noa,
c
     1     1,noa,noa+1,noab,
     1     noab+1,noab+nva,noa+1,noab,
c
     1     noa+1,noab,1,noa,
     1     noab+nva+1,noab+nvab,1,noa,
c
     1     noa+1,noab,noa+1,noab,
     1     noab+nva+1,noab+nvab,noa+1,noab)
c       call ga_sync( )
c
c       print *,'H-H-P-H is fine!'
c
c     =======
c     H-H-P-P
c     =======
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     1,noa,1,noa,
     1     noab+1,noab+nva,noab+1,noab+nva,
c
     1     1,noa,noa+1,noab,
     1     noab+1,noab+nva,noab+nva+1,noab+nvab,
c
     1     noa+1,noab,1,noa,
     1     noab+nva+1,noab+nvab,noab+1,noab+nva,
c
     1     noa+1,noab,noa+1,noab,
     1     noab+nva+1,noab+nvab,noab+nva+1,noab+nvab)
c      call ga_sync( )
c
c      print *,'H-H-P-P is fine!'
c
c     =======
c     H-P-H-H
c     =======
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     1,noa,noab+1,noab+nva,
     1     1,noa,1,noa,
c
     1     1,noa,noab+nva+1,noab+nvab,
     1     1,noa,noa+1,noab,
c
     1     noa+1,noab,noab+1,noab+nva,
     1     noa+1,noab,1,noa,
c
     1     noa+1,noab,noab+nva+1,noab+nvab,
     1     noa+1,noab,noa+1,noab)
c      call ga_sync( )
c
c      print *,'H-P-H-H is fine!'
c
c     =======
c     P-H-H-H
c     =======
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     noab+1,noab+nva,1,noa,
     1     1,noa,1,noa,
c
     1     noab+1,noab+nva,noa+1,noab,
     1     1,noa,noa+1,noab,
c
     1     noab+nva+1,noab+nvab,1,noa,
     1     noa+1,noab,1,noa,
c
     1     noab+nva+1,noab+nvab,noa+1,noab,
     1     noa+1,noab,noa+1,noab)
c      call ga_sync( )
c
c      print *,'P-H-H-H is fine!'
c
c     =======
c     H-P-H-P
c     =======
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     1,noa,noab+1,noab+nva,
     1     1,noa,noab+1,noab+nva,
c
     1     1,noa,noab+nva+1,noab+nvab,
     1     1,noa,noab+nva+1,noab+nvab,
c
     1     noa+1,noab,noab+1,noab+nva,
     1     noa+1,noab,noab+1,noab+nva,
c
     1     noa+1,noab,noab+nva+1,noab+nvab,
     1     noa+1,noab,noab+nva+1,noab+nvab)
c      call ga_sync( )
c
c      print *,'H-P-H-P is fine!'
c
c     =======
c     H-P-P-H
c     =======
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     1,noa,noab+1,noab+nva,
     1     noab+1,noab+nva,1,noa,
c
     1     1,noa,noab+nva+1,noab+nvab,
     1     noab+1,noab+nva,noa+1,noab,
c
     1     noa+1,noab,noab+1,noab+nva,
     1     noab+nva+1,noab+nvab,1,noa,
c
     1     noa+1,noab,noab+nva+1,noab+nvab,
     1     noab+nva+1,noab+nvab,noa+1,noab)
c      call ga_sync( )
c
c      print *,'H-P-P-H is fine!'
c
c     =======
c     P-H-H-P
c     =======
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     noab+1,noab+nva,1,noa,
     1     1,noa,noab+1,noab+nva,
c
     1     noab+1,noab+nva,noa+1,noab,
     1     1,noa,noab+nva+1,noab+nvab,
c
     1     noab+nva+1,noab+nvab,1,noa,
     1     noa+1,noab,noab+1,noab+nva,
c
     1     noab+nva+1,noab+nvab,noa+1,noab,
     1     noa+1,noab,noab+nva+1,noab+nvab)
c      call ga_sync( )
c      print *,'P-H-H-P is fine'
c     =======
c     P-H-P-H
c     =======
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     noab+1,noab+nva,1,noa,
     1     noab+1,noab+nva,1,noa,
c
     1     noab+1,noab+nva,noa+1,noab,
     1     noab+1,noab+nva,noa+1,noab,
c
     1     noab+nva+1,noab+nvab,1,noa,
     1     noab+nva+1,noab+nvab,1,noa,
c
     1     noab+nva+1,noab+nvab,noa+1,noab,
     1     noab+nva+1,noab+nvab,noa+1,noab)
c      call ga_sync( )
c      print *,'P-H-P-H is fine'
c     =======
c     H-P-P-P
c     =======
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     1,noa,noab+1,noab+nva,
     1     noab+1,noab+nva,noab+1,noab+nva,
c
     1     1,noa,noab+nva+1,noab+nvab,
     1     noab+1,noab+nva,noab+nva+1,noab+nvab,
c
     1     noa+1,noab,noab+1,noab+nva,
     1     noab+nva+1,noab+nvab,noab+1,noab+nva,
c
     1     noa+1,noab,noab+nva+1,noab+nvab,
     1     noab+nva+1,noab+nvab,noab+nva+1,noab+nvab)
c      call ga_sync( )
c      print *,'H-P-P-P is fine'
c     =======
c     P-H-P-P
c     =======
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     noab+1,noab+nva,1,noa,
     1     noab+1,noab+nva,noab+1,noab+nva,
c
     1     noab+1,noab+nva,noa+1,noab,
     1     noab+1,noab+nva,noab+nva+1,noab+nvab,
c
     1     noab+nva+1,noab+nvab,1,noa,
     1     noab+nva+1,noab+nvab,noab+1,noab+nva,
c
     1     noab+nva+1,noab+nvab,noa+1,noab,
     1     noab+nva+1,noab+nvab,noab+nva+1,noab+nvab)
c      call ga_sync( )
c      print *,'P-H-P-P is fine!'
c     =======
c     P-P-H-H
c     =======
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     noab+1,noab+nva,noab+1,noab+nva,
     1     1,noa,1,noa,
c
     1     noab+1,noab+nva,noab+nva+1,noab+nvab,
     1     1,noa,noa+1,noab,
c
     1     noab+nva+1,noab+nvab,noab+1,noab+nva,
     1     noa+1,noab,1,noa,
c
     1     noab+nva+1,noab+nvab,noab+nva+1,noab+nvab,
     1     noa+1,noab,noa+1,noab)
c      call ga_sync( )
cc      print *,'P-P-H-H is fine!'
c     =======
c     P-P-H-P
c     =======
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     noab+1,noab+nva,noab+1,noab+nva,
     1     1,noa,noab+1,noab+nva,
c
     1     noab+1,noab+nva,noab+nva+1,noab+nvab,
     1     1,noa,noab+nva+1,noab+nvab,
c
     1     noab+nva+1,noab+nvab,noab+1,noab+nva,
     1     noa+1,noab,noab+1,noab+nva,
c
     1     noab+nva+1,noab+nvab,noab+nva+1,noab+nvab,
     1     noa+1,noab,noab+nva+1,noab+nvab)
c      call ga_sync( )
c      print *,'P-P-H-P is fine!'
c     =======
c     P-P-P-H
c     =======
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     noab+1,noab+nva,noab+1,noab+nva,
     1     noab+1,noab+nva,1,noa,
c
     1     noab+1,noab+nva,noab+nva+1,noab+nvab,
     1     noab+1,noab+nva,noa+1,noab,
c
     1     noab+nva+1,noab+nvab,noab+1,noab+nva,
     1     noab+nva+1,noab+nvab,1,noa,
c
     1     noab+nva+1,noab+nvab,noab+nva+1,noab+nvab,
     1     noab+nva+1,noab+nvab,noa+1,noab)
c      call ga_sync( )
c      print *,'P-P-P-H is fine!'
c     =======
c     P-P-P-P
c     =======
      call btrans2_all(d_2pdm,k_2pdm_offset,
     1     d_i0,k_i0_offset,size_i0,
     1     d_i1,k_i1_offset,size_i1,
     1     d_i2,k_i2_offset,size_i2,
     1     d_i3,k_i3_offset,size_i3,
     1     atpart2,nalength2,
     1     noab+1,noab+nva,noab+1,noab+nva,
     1     noab+1,noab+nva,noab+1,noab+nva,
c
     1     noab+1,noab+nva,noab+nva+1,noab+nvab,
     1     noab+1,noab+nva,noab+nva+1,noab+nvab,
c
     1     noab+nva+1,noab+nvab,noab+1,noab+nva,
     1     noab+nva+1,noab+nvab,noab+1,noab+nva,
c
     1     noab+nva+1,noab+nvab,noab+nva+1,noab+nvab,
     1     noab+nva+1,noab+nvab,noab+nva+1,noab+nvab)
c      endif !two
c
c      print *,'2PDM is fine!'
c     ======================
c     READ d_i0 to k_2pdm_ao
c     ======================
c      i_ind=0
c      do ib=1,atpart2
c      j_ind=0
c      do jb=1,atpart2
c      k_ind=0
c      do kb=1,atpart2
c      l_ind=0
c      do lb=1,atpart2
c         dima=nalength2(ib)*nalength2(jb)*nalength2(kb)*nalength2(lb)
c         if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
c     1   call errquit('btrans2: ma problem',3,ma_err)
c         call get_hash_block(d_i0,dbl_mb(k_a),dima,int_mb(k_i0_offset),
c     1        (lb-1)+(kb-1)*atpart2+(jb-1)*atpart2*atpart2+
c     1        (ib-1)*atpart2*atpart2*atpart2)
c         do i=1,nalength2(ib)
c            i1=i_ind+i
c         do j=1,nalength2(jb)
c            j1=j_ind+j
c         do k=1,nalength2(kb)
c            k1=k_ind+k
c         do l=1,nalength2(lb)
c            l1=l_ind+l
c            dbl_mb(k_2pdm_ao+(l1-1)+(j1-1)*nbf+(k1-1)*nbf*nbf
c     1      +(i1-1)*nbf*nbf*nbf) = 
c     1      dbl_mb(k_2pdm_ao+(l1-1)+(j1-1)*nbf+(k1-1)*nbf*nbf
c     1      +(i1-1)*nbf*nbf*nbf) +
c     1      dbl_mb(k_a+(l-1)+(k-1)*nalength2(lb)+
c     1      (j-1)*nalength2(lb)*nalength2(kb)+
c     1      (i-1)*nalength2(lb)*nalength2(kb)*nalength2(jb))
c         enddo
c         enddo
c         enddo
c         enddo
c         if(.not.ma_pop_stack(l_a))
c     1   call errquit('btrans2: ma problem',4,ma_err)
c         l_ind=l_ind+nalength2(lb)
c      enddo
c      k_ind=k_ind+nalength2(kb)
c      enddo
c      j_ind=j_ind+nalength2(jb)
c      enddo
c      i_ind=i_ind+nalength2(ib)
c      enddo
c
      call deletefile(d_i3)
      if (.not.ma_pop_stack(l_i3_offset))
     1   call errquit('btrans2: ma problem',1,ma_err)
      call deletefile(d_i2)
      if (.not.ma_pop_stack(l_i2_offset))
     1   call errquit('btrans2: ma problem',2,ma_err)
      call deletefile(d_i1)
      if (.not.ma_pop_stack(l_i1_offset))
     1   call errquit('btrans2: ma problem',3,ma_err)
c      call deletefile(d_i0)
c      if (.not.ma_pop_stack(l_i0_offset))
c     1   call errquit('btrans2: ma problem',4,ma_err)
c
      call ga_sync( )
c
      end
c
c
      subroutine btrans2_all(d_a,k_a_offset,
     1 d_i0,k_i0_offset,size_i0,
     1 d_i1,k_i1_offset,size_i1,
     1 d_i2,k_i2_offset,size_i2,
     1 d_i3,k_i3_offset,size_i3,
     1 atpart2,nalength2,
     1 g1b_lo_aaaa,g1b_hi_aaaa,
     1 g2b_lo_aaaa,g2b_hi_aaaa,
     1 g3b_lo_aaaa,g3b_hi_aaaa,
     1 g4b_lo_aaaa,g4b_hi_aaaa,
     1 g1b_lo_aabb,g1b_hi_aabb,
     1 g2b_lo_aabb,g2b_hi_aabb,
     1 g3b_lo_aabb,g3b_hi_aabb,
     1 g4b_lo_aabb,g4b_hi_aabb,
     1 g1b_lo_bbaa,g1b_hi_bbaa,
     1 g2b_lo_bbaa,g2b_hi_bbaa,
     1 g3b_lo_bbaa,g3b_hi_bbaa,
     1 g4b_lo_bbaa,g4b_hi_bbaa,
     1 g1b_lo_bbbb,g1b_hi_bbbb,
     1 g2b_lo_bbbb,g2b_hi_bbbb,
     1 g3b_lo_bbbb,g3b_hi_bbbb,
     1 g4b_lo_bbbb,g4b_hi_bbbb)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_a,k_a_offset
      integer d_i0,k_i0_offset,size_i0
      integer d_i1,k_i1_offset,size_i1
      integer d_i2,k_i2_offset,size_i2
      integer d_i3,k_i3_offset,size_i3
      integer atpart2,nalength2(*)
      integer g1b_lo_aaaa,g1b_hi_aaaa
      integer g2b_lo_aaaa,g2b_hi_aaaa
      integer g3b_lo_aaaa,g3b_hi_aaaa
      integer g4b_lo_aaaa,g4b_hi_aaaa
      integer g1b_lo_aabb,g1b_hi_aabb
      integer g2b_lo_aabb,g2b_hi_aabb
      integer g3b_lo_aabb,g3b_hi_aabb
      integer g4b_lo_aabb,g4b_hi_aabb
      integer g1b_lo_bbaa,g1b_hi_bbaa
      integer g2b_lo_bbaa,g2b_hi_bbaa
      integer g3b_lo_bbaa,g3b_hi_bbaa
      integer g4b_lo_bbaa,g4b_hi_bbaa
      integer g1b_lo_bbbb,g1b_hi_bbbb
      integer g2b_lo_bbbb,g2b_hi_bbbb
      integer g3b_lo_bbbb,g3b_hi_bbbb
      integer g4b_lo_bbbb,g4b_hi_bbbb
c
c     AAAA
c
      call tce_zero(d_i1,size_i1)
      call tce_zero(d_i2,size_i2)
      call tce_zero(d_i3,size_i3)
      call btrans2_i1(d_a,k_a_offset,d_i1,k_i1_offset,
     1     atpart2,nalength2,
     1     g1b_lo_aaaa,g1b_hi_aaaa,
     1     g2b_lo_aaaa,g2b_hi_aaaa,
     1     g3b_lo_aaaa,g3b_hi_aaaa,
     1     g4b_lo_aaaa,g4b_hi_aaaa)
      call btrans2_i2(d_i1,k_i1_offset,d_i2,k_i2_offset,
     1     atpart2,nalength2,
     1     g1b_lo_aaaa,g1b_hi_aaaa,
     1     g2b_lo_aaaa,g2b_hi_aaaa,
     1     g3b_lo_aaaa,g3b_hi_aaaa)
      call btrans2_i3(d_i2,k_i2_offset,d_i3,k_i3_offset,
     1     atpart2,nalength2,
     1     g1b_lo_aaaa,g1b_hi_aaaa,
     1     g2b_lo_aaaa,g2b_hi_aaaa)
      call btrans2_i0(d_i3,k_i3_offset,d_i0,k_i0_offset,
     1     atpart2,nalength2,
     1     g1b_lo_aaaa,g1b_hi_aaaa,.true.)
c
c
c     AABB
c
      call tce_zero(d_i1,size_i1)
      call tce_zero(d_i2,size_i2)
      call tce_zero(d_i3,size_i3)
      call btrans2_i1(d_a,k_a_offset,d_i1,k_i1_offset,
     1     atpart2,nalength2,
     1     g1b_lo_aabb,g1b_hi_aabb,
     1     g2b_lo_aabb,g2b_hi_aabb,
     1     g3b_lo_aabb,g3b_hi_aabb,
     1     g4b_lo_aabb,g4b_hi_aabb)
      call btrans2_i2(d_i1,k_i1_offset,d_i2,k_i2_offset,
     1     atpart2,nalength2,
     1     g1b_lo_aabb,g1b_hi_aabb,
     1     g2b_lo_aabb,g2b_hi_aabb,
     1     g3b_lo_aabb,g3b_hi_aabb)
      call btrans2_i3(d_i2,k_i2_offset,d_i3,k_i3_offset,
     1     atpart2,nalength2,
     1     g1b_lo_aabb,g1b_hi_aabb,
     1     g2b_lo_aabb,g2b_hi_aabb)
      call btrans2_i0(d_i3,k_i3_offset,d_i0,k_i0_offset,
     1     atpart2,nalength2,
     1     g1b_lo_aabb,g1b_hi_aabb,.false.)
c
c
c     BBAA
c
      call tce_zero(d_i1,size_i1)
      call tce_zero(d_i2,size_i2)
      call tce_zero(d_i3,size_i3)
      call btrans2_i1(d_a,k_a_offset,d_i1,k_i1_offset,
     1     atpart2,nalength2,
     1     g1b_lo_bbaa,g1b_hi_bbaa,
     1     g2b_lo_bbaa,g2b_hi_bbaa,
     1     g3b_lo_bbaa,g3b_hi_bbaa,
     1     g4b_lo_bbaa,g4b_hi_bbaa)
      call btrans2_i2(d_i1,k_i1_offset,d_i2,k_i2_offset,
     1     atpart2,nalength2,
     1     g1b_lo_bbaa,g1b_hi_bbaa,
     1     g2b_lo_bbaa,g2b_hi_bbaa,
     1     g3b_lo_bbaa,g3b_hi_bbaa)
      call btrans2_i3(d_i2,k_i2_offset,d_i3,k_i3_offset,
     1     atpart2,nalength2,
     1     g1b_lo_bbaa,g1b_hi_bbaa,
     1     g2b_lo_bbaa,g2b_hi_bbaa)
      call btrans2_i0(d_i3,k_i3_offset,d_i0,k_i0_offset,
     1     atpart2,nalength2,
     1     g1b_lo_bbaa,g1b_hi_bbaa,.false.)
c
c
c     BBBB
c
      if(.not.restricted) then
      call tce_zero(d_i1,size_i1)
      call tce_zero(d_i2,size_i2)
      call tce_zero(d_i3,size_i3)
      call btrans2_i1(d_a,k_a_offset,d_i1,k_i1_offset,
     1     atpart2,nalength2,
     1     g1b_lo_bbbb,g1b_hi_bbbb,
     1     g2b_lo_bbbb,g2b_hi_bbbb,
     1     g3b_lo_bbbb,g3b_hi_bbbb,
     1     g4b_lo_bbbb,g4b_hi_bbbb)
      call btrans2_i2(d_i1,k_i1_offset,d_i2,k_i2_offset,
     1     atpart2,nalength2,
     1     g1b_lo_bbbb,g1b_hi_bbbb,
     1     g2b_lo_bbbb,g2b_hi_bbbb,
     1     g3b_lo_bbbb,g3b_hi_bbbb)
      call btrans2_i3(d_i2,k_i2_offset,d_i3,k_i3_offset,
     1     atpart2,nalength2,
     1     g1b_lo_bbbb,g1b_hi_bbbb,
     1     g2b_lo_bbbb,g2b_hi_bbbb)
      call btrans2_i0(d_i3,k_i3_offset,d_i0,k_i0_offset,
     1     atpart2,nalength2,
     1     g1b_lo_bbbb,g1b_hi_bbbb,.false.)
      endif
c
      end
      
      subroutine btrans2_i0(d_a,k_a_offset,d_i0,k_i0_offset,
     1 atpart2,nalength2,
     1 g1b_lo,g1b_hi,alpha4)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
      integer d_a,k_a_offset
      integer d_i0,k_i0_offset
      integer ib,jb,kb,lb,atpart2,nalength2(*)
      integer g1b,g1b_lo,g1b_hi
      integer g1
      logical alpha4
c
      integer dimc,l_c,k_c,dimc_sort
      integer l_c0,k_c0
      integer dimb,l_b,k_b,dimb_sort
      integer dima,l_a,k_a,dima_sort
      integer i_ind,j_ind,k_ind,l_ind,i,j,k,l,i1,j1,k1,l1
      integer dim_common
      double precision res
c
      integer nxtask
      integer next
      integer count
      integer nprocs
      external nxtask
c
      nprocs=ga_nnodes( )
      count=0
      next=nxtask(nprocs,1)
c
c      l_ind=0
      do lb=1,atpart2
      do kb=1,atpart2
      do jb=1,atpart2
      do ib=1,atpart2
         if(next.eq.count) then
         l_ind=0
         do i=1,lb-1
            l_ind=l_ind+nalength2(i)
         enddo
c
         dimc=nalength2(ib)*nalength2(jb)*nalength2(kb)*nalength2(lb)
         if(.not.ma_push_get(mt_dbl,dimc,'noname',l_c,k_c))
     1   call errquit('btrans2_i0: ma problem',1,ma_err)
         do i=1,dimc
            dbl_mb(k_c+i-1)=0.d0
         enddo
c
         if(.not.ma_push_get(mt_dbl,dimc,'noname',l_c0,k_c0))
     1   call errquit('btrans2_i0: ma problem',1,ma_err)
         call get_hash_block(d_i0,dbl_mb(k_c0),dimc,
     1        int_mb(k_i0_offset),
     1        (ib-1)+(jb-1)*atpart2+(kb-1)*atpart2*atpart2
     1        +(lb-1)*atpart2*atpart2*atpart2)
c
         do g1b=g1b_lo,g1b_hi
            dima=int_mb(k_range+g1b-1)*nalength2(kb)*
     1           nalength2(jb)*nalength2(ib)
            if(.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1      call errquit('btrans2_i0: ma problem',2,ma_err)
c
            dimb=int_mb(k_range+g1b-1)*nalength2(lb)
            if(.not.ma_push_get(mt_dbl,dimb,'noname',l_b,k_b))
     1      call errquit('btrans2_i0: ma problem',3,ma_err)
c
            call get_hash_block(d_a,dbl_mb(k_a),dima,
     1           int_mb(k_a_offset),
     1           (g1b-1)+(ib-1)*(noab+nvab)+(jb-1)*(noab+nvab)*atpart2
     1          +(kb-1)*(noab+nvab)*atpart2*atpart2)
c
            do l=1,nalength2(lb)
               l1=l_ind+l
            do g1=1,int_mb(k_range+g1b-1)
               dbl_mb(k_b+(g1-1)+(l-1)*int_mb(k_range+g1b-1))=
     1         dbl_mb(k_movecs_sorted+
     1         (l1-1)+(int_mb(k_offset+g1b-1)+g1-1)*nbf)
            enddo
            enddo
c
            dima_sort=nalength2(kb)*nalength2(jb)*nalength2(ib)
            dimb_sort=nalength2(lb)
            dim_common=int_mb(k_range+g1b-1)
c
            call dgemm('t','n',dima_sort,dimb_sort,dim_common,1.d0,
     1           dbl_mb(k_a),dim_common,
     1           dbl_mb(k_b),dim_common,1.d0,dbl_mb(k_c),dima_sort)
c
            if (.not.ma_pop_stack(l_b))
     1      call errquit('btrans2_i0: ma problem',4,ma_err)
            if(.not.ma_pop_stack(l_a))
     1      call errquit('btrans2_i0: ma problem',5,ma_err)
         enddo
         if (restricted.and.alpha4) then
            do i=1,dimc
               dbl_mb(k_c+i-1)=2.d0*dbl_mb(k_c+i-1)
            enddo
         endif
         do i=1,dimc
            dbl_mb(k_c+i-1)=0.5d0*dbl_mb(k_c+i-1)
         enddo
         do i=1,dimc
            dbl_mb(k_c+i-1)= dbl_mb(k_c+i-1) + dbl_mb(k_c0+i-1)
         enddo
         call put_hash_block(d_i0,dbl_mb(k_c),dimc,
     1        int_mb(k_i0_offset),
     1        (ib-1)+(jb-1)*atpart2+(kb-1)*atpart2*atpart2
     1        +(lb-1)*atpart2*atpart2*atpart2)
         if (.not.ma_pop_stack(l_c0))
     1      call errquit('btrans2_i0: ma problem',7,ma_err)
         if (.not.ma_pop_stack(l_c))
     1      call errquit('btrans2_i0: ma problem',6,ma_err)
c
         next=nxtask(nprocs,1)
         endif
         count=count+1
      enddo
      enddo
      enddo
c      l_ind=l_ind+nalength2(lb)
      enddo
      next=nxtask(-nprocs,1)
      call ga_sync( )
      return
c
      end
c
      subroutine btrans2_i3(d_a,k_a_offset,d_c,k_c_offset,
     1   atpart2,nalength2,
     1   g1b_lo,g1b_hi,
     1   g2b_lo,g2b_hi) 
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
      integer d_a,k_a_offset,d_c,k_c_offset
      integer ib,jb,kb,atpart2,nalength2(*)
      integer g1b,g1b_lo,g1b_hi
      integer g2b,g2b_lo,g2b_hi
c
      integer dimc,l_c,k_c
      integer dimb,l_b,k_b,dimb_sort
      integer dima,l_a,k_a,dima_sort
      integer dim_common
c
      integer i_ind,j_ind,k_ind,i1,j1,k1,i,j,k
c
      integer g1,g2
c
      integer nxtask
      integer next
      integer count
      integer nprocs
      external nxtask
c
      nprocs=ga_nnodes( )
      count=0
      next=nxtask(nprocs,1)
c
c      k_ind=0
      do kb=1,atpart2
c      j_ind=0
      do jb=1,atpart2
c      i_ind=0
      do ib=1,atpart2
      do g1b=g1b_lo,g1b_hi
         if (next.eq.count) then
c
c         i_ind=0
c         do g1=1,ib-1
c            i_ind=i_ind+nalength2(g1)
c         enddo
c         j_ind=0
c         do g1=1,jb-1
c            j_ind=j_ind+nalength2(g1)
c         enddo
         k_ind=0
         do g1=1,kb-1
            k_ind=k_ind+nalength2(g1)
         enddo
c
         dimc = int_mb(k_range+g1b-1)*nalength2(ib)*
     1          nalength2(jb)*nalength2(kb)
         if(.not.ma_push_get(mt_dbl,dimc,'noname',l_c,k_c))
     1   call errquit('btrans2_i3: ma problem',1,ma_err)
c
         do g1=1,dimc
            dbl_mb(k_c+g1-1)=0.d0
         enddo
c
            do g2b=g2b_lo,g2b_hi
               dima=int_mb(k_range+g2b-1)*int_mb(k_range+g1b-1)*
     1              nalength2(jb)*nalength2(ib)
               if(.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1         call errquit('btrans2_i3: ma problem',2,ma_err)
c
               call get_hash_block(d_a,dbl_mb(k_a),dima,
     1              int_mb(k_a_offset),
     1              (g2b-1)+(g1b-1)*(noab+nvab)
     1             +(ib-1)*(noab+nvab)*(noab+nvab)
     1             +(jb-1)*(noab+nvab)*(noab+nvab)*atpart2)
c
               dimb=int_mb(k_range+g2b-1)*nalength2(kb)
               if (.not.ma_push_get(mt_dbl,dimb,'noname',l_b,k_b))
     1         call errquit('btrans2_i3: ma problem',3,ma_err)
c
               do k=1,nalength2(kb)
                  k1=k_ind+k
               do g2=1,int_mb(k_range+g2b-1)
                  dbl_mb(k_b+(g2-1)+(k-1)*int_mb(k_range+g2b-1))=
     1            dbl_mb(k_movecs_sorted+
     1            (k1-1)+(int_mb(k_offset+g2b-1)+g2-1)*nbf)
               enddo
               enddo
c
               dima_sort=int_mb(k_range+g1b-1)*nalength2(jb)*
     1         nalength2(ib)
               dimb_sort=nalength2(kb)
               dim_common=int_mb(k_range+g2b-1)
c
               call dgemm('t','n',dima_sort,dimb_sort,dim_common,1.d0,
     1              dbl_mb(k_a),dim_common,dbl_mb(k_b),dim_common,
     1              1.d0,dbl_mb(k_c),dima_sort)
c
               if (.not.ma_pop_stack(l_b))
     1         call errquit('btrans2_i3: ma problem',4,ma_err)
               if (.not.ma_pop_stack(l_a))
     1         call errquit('btrans2_i3: ma problem',5,ma_err)
            enddo
c
            call put_hash_block(d_c,dbl_mb(k_c),
     1        dimc,int_mb(k_c_offset),
     1        (g1b-1)+(ib-1)*(noab+nvab)+(jb-1)*(noab+nvab)*atpart2+
     1        (kb-1)*(noab+nvab)*atpart2*atpart2)
c
            if (.not.ma_pop_stack(l_c))
     1      call errquit('btrans2_hhpp: ma problem',6,ma_err)
c
            next=nxtask(nprocs,1)
            endif
            count=count+1
       enddo
c       i_ind=i_ind+nalength2(ib)
       enddo
c       j_ind=j_ind+nalength2(jb)
       enddo
c       k_ind=k_ind+nalength2(kb)
       enddo
       next=nxtask(-nprocs,1)
       call ga_sync( )
       return
      end
c
      subroutine btrans2_i2(d_a,k_a_offset,d_c,k_c_offset,
     1 atpart2,nalength2,
     1 g1b_lo,g1b_hi,
     1 g2b_lo,g2b_hi,
     1 g3b_lo,g3b_hi)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
      integer d_a,k_a_offset,d_c,k_c_offset
      integer jb,ib,atpart2,nalength2(*)
      integer g3b,g3b_lo,g3b_hi
      integer g2b,g2b_lo,g2b_hi
      integer g1b,g1b_lo,g1b_hi
c
      integer g3,g2,g1
      integer dima,dima_sort,l_a,k_a
      integer dimb,dimb_sort,l_b,k_b
      integer dimc,l_c,k_c
      integer dim_common
c
      integer j,i,j1,i1,j_ind,i_ind
c
      integer nxtask
      integer next
      integer count
      integer nprocs
      external nxtask
c
      nprocs=ga_nnodes( )
      count=0
      next = nxtask(nprocs,1)
c
c      j_ind=0
      do jb=1,atpart2
c      i_ind=0
      do ib=1,atpart2
      do g1b=g1b_lo,g1b_hi
      do g2b=g2b_lo,g2b_hi
         if(next.eq.count) then
c
c         i_ind=0
c         do g2=1,ib-1
c            i_ind=i_ind+nalength2(g2)
c         enddo
         j_ind=0
         do g2=1,jb-1
            j_ind=j_ind+nalength2(g2)
         enddo
c         k_ind=0
c         do g1=1,kb-1
c            k_ind=k_ind+nalength2(g1)
c         enddo
c
         dimc=int_mb(k_range+g1b-1)*int_mb(k_range+g2b-1)*
     1        nalength2(ib)*nalength2(jb)
         if(.not.ma_push_get(mt_dbl,dimc,'noname',l_c,k_c))
     1     call errquit('btrans2_i2: ma problem',1,ma_err)
         do g2=1,dimc
            dbl_mb(k_c+g2-1)=0.d0
         enddo
c
         do g3b=g3b_lo,g3b_hi
            dima=int_mb(k_range+g1b-1)*int_mb(k_range+g2b-1)*
     1           int_mb(k_range+g3b-1)*nalength2(ib)
            if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1      call errquit('btrans2_i2: ma problem',11,ma_err)
            call get_hash_block(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),
     1           (g3b-1)+(g2b-1)*(noab+nvab)
     1          +(g1b-1)*(noab+nvab)*(noab+nvab)+
     1           (ib-1)*(noab+nvab)*(noab+nvab)*(noab+nvab))
c
            dimb=int_mb(k_range+g3b-1)*nalength2(jb)
            if (.not.ma_push_get(mt_dbl,dimb,'noname',l_b,k_b))
     1         call errquit('btrans2_i2: ma problem',2,ma_err)
c
            do j=1,nalength2(jb)
               j1=j_ind+j
            do g3=1,int_mb(k_range+g3b-1)
               dbl_mb(k_b+(g3-1)+(j-1)*int_mb(k_range+g3b-1))=
     1         dbl_mb(k_movecs_sorted+
     1         (j1-1)+(int_mb(k_offset+g3b-1)+g3-1)*nbf)
            enddo
            enddo
c
            dima_sort=int_mb(k_range+g1b-1)*int_mb(k_range+g2b-1)*
     1                nalength2(ib)
            dimb_sort=nalength2(jb)
            dim_common=int_mb(k_range+g3b-1)
c
            call dgemm('t','n',dima_sort,dimb_sort,dim_common,1.d0,
     1           dbl_mb(k_a),dim_common,
     1           dbl_mb(k_b),dim_common,1.d0,
     1           dbl_mb(k_c),dima_sort)
c
            if (.not.ma_pop_stack(l_b))
     1         call errquit('btrans2_i2: ma problem',3,ma_err)
            if (.not.ma_pop_stack(l_a))
     1         call errquit('btrans2_i2: ma problem',4,ma_err)
c
         enddo
c
         call put_hash_block(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),
     1        (g2b-1)+(g1b-1)*(noab+nvab)
     1       +(ib-1)*(noab+nvab)*(noab+nvab)
     1       +(jb-1)*(noab+nvab)*(noab+nvab)*atpart2)
         if(.not.ma_pop_stack(l_c))
     1   call errquit('yes error',1,ma_err)
         next=nxtask(nprocs,1)
         endif
         count=count+1
      enddo
      enddo
c      i_ind=i_ind+nalength2(ib)
      enddo
c      j_ind=j_ind+nalength2(jb)
      enddo
      next=nxtask(-nprocs,1)
      call ga_sync( )
      return
      end
c
      subroutine btrans2_i1(d_a,k_a_offset,d_c,k_c_offset,
     1  atpart2,nalength2,
     1  g1b_lo,g1b_hi,
     1  g2b_lo,g2b_hi,
     1  g3b_lo,g3b_hi,
     1  g4b_lo,g4b_hi)
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
      integer d_a,k_a_offset,d_c,k_c_offset
      integer ib,atpart2,nalength2(*)
      integer g1b,g1b_lo,g1b_hi,g2b,g2b_lo,g2b_hi
      integer g3b,g3b_lo,g3b_hi,g4b,g4b_lo,g4b_hi
c
      integer g1,g2,g3,g4
      integer dima,k_a,l_a,dima_sort,l_a_sort,k_a_sort
      integer dimb,k_b,l_b,dimb_sort,l_b_sort,k_b_sort
      integer dimc,k_c,l_c,dimc_sort,l_c_sort,k_c_sort
      integer dim_common
      integer i,i1,i_ind
c
      integer nxtask
      integer next
      integer count
      integer nprocs
      external nxtask
c
      nprocs=ga_nnodes( )
      count=0
      next=nxtask(nprocs,1)
c
c      i_ind=0
      do ib=1,atpart2
      do g1b=g1b_lo,g1b_hi
      do g2b=g2b_lo,g2b_hi
      do g3b=g3b_lo,g3b_hi
         if (next.eq.count) then
c
         i_ind=0
         do i=1,ib-1
            i_ind=i_ind+nalength2(i)
         enddo
c
            dimc = int_mb(k_range+g1b-1)*int_mb(k_range+g2b-1)*
     1             int_mb(k_range+g3b-1)*nalength2(ib)
            if (.not.ma_push_get(mt_dbl,dimc,'noname',l_c,k_c))
     1         call errquit('btrans2_i1: ma problem',2,ma_err)
            call dfill(dimc,0.d0,dbl_mb(k_c),1)
            do g4b=g4b_lo,g4b_hi
               if (ieor(int_mb(k_sym+g1b-1),ieor(int_mb(k_sym+g2b-1),
     1             ieor(int_mb(k_sym+g3b-1),int_mb(k_sym+g4b-1))))
     1             .eq. irrep_e) then
                  dima = int_mb(k_range+g1b-1)*int_mb(k_range+g2b-1)
     1                  *int_mb(k_range+g3b-1)*int_mb(k_range+g4b-1)
                  if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1            call errquit('btrans2_i1: ma problem',10,ma_err)
                  call get_hash_block(d_a,dbl_mb(k_a),dima,
     1                int_mb(k_a_offset),
     1                g4b-1+(noab+nvab)*(g3b-1+(noab+nvab)*(g2b-1+
     1                (noab+nvab)*(g1b-1))))
                  dimb=int_mb(k_range+g4b-1)*nalength2(ib)
                  dimb_sort=nalength2(ib)
                  dim_common=int_mb(k_range+g4b-1)
                  if(.not.ma_push_get(mt_dbl,dimb,'noname',l_b,k_b))
     1            call errquit('btrans2_i1: ma problem',11,ma_err)
                  do i=1,nalength2(ib)
                     i1=i_ind+i
                  do g4=1,int_mb(k_range+g4b-1)
                     dbl_mb(k_b+(g4-1)+(i-1)*int_mb(k_range+g4b-1))=
     1               dbl_mb(k_movecs_sorted+(i1-1)
     1               +(int_mb(k_offset+g4b-1)+g4-1)*nbf)
                  enddo
                  enddo
                  dima_sort = int_mb(k_range+g1b-1)*
     1                  int_mb(k_range+g2b-1)*int_mb(k_range+g3b-1)
                  call dgemm('t','n',dima_sort,dimb_sort,dim_common,
     1                 1.d0,dbl_mb(k_a),dim_common,
     1                 dbl_mb(k_b),dim_common,1.d0,
     1                 dbl_mb(k_c),dima_sort)
                  if (.not.ma_pop_stack(l_b))
     1            call errquit('btrans2_i1: ma problem',12,ma_err)
                  if (.not.ma_pop_stack(l_a))
     1            call errquit('btrans2_i1: ma problem',14,ma_err)
               endif
            enddo
            call put_hash_block(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),
     1           (g3b-1)+(g2b-1)*(noab+nvab)
     1          +(g1b-1)*(noab+nvab)*(noab+nvab)
     1          +(ib-1)*(noab+nvab)*(noab+nvab)*(noab+nvab))
            if (.not.ma_pop_stack(l_c))
     1         call errquit('btrans2_i1',3,ma_err)
            next=nxtask(nprocs,1)
         endif
         count = count + 1
      enddo
      enddo
      enddo
c      i_ind=i_ind+nalength2(ib)
      enddo
      next = nxtask(-nprocs,1)
      call ga_sync( )
      return
      end

c     =======
c     offsets
c     =======
      subroutine offset_btrans2_i1(l_a_offset,k_a_offset,size,
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
      integer g1b,g2b,g3b
c
      length = 0
      do ib=1,atpart2
      do g1b=1,noab+nvab
      do g2b=1,noab+nvab
      do g3b=1,noab+nvab
         length=length+1
      enddo
      enddo
      enddo
      enddo
      if (.not. ma_push_get(mt_int,2*length+1,'noname',l_a_offset,
     1    k_a_offset))
     1  call errquit('offset_btrans2_i1',0,ma_err)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      do ib=1,atpart2
      do g1b=1,noab+nvab
      do g2b=1,noab+nvab
      do g3b=1,noab+nvab
         addr = addr + 1
         int_mb(k_a_offset+addr) = (g3b-1)+(g2b-1)*(noab+nvab)
     1          +(g1b-1)*(noab+nvab)*(noab+nvab)
     1          +(ib-1)*(noab+nvab)*(noab+nvab)*(noab+nvab)
         int_mb(k_a_offset+length+addr) = size
         size = size + int_mb(k_range+g1b-1)*int_mb(k_range+g2b-1)*
     1          int_mb(k_range+g3b-1)*nalength2(ib)
         enddo
      enddo
      enddo
      enddo
      end
c
      subroutine offset_btrans2_i2(l_a_offset,k_a_offset,size,
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
      integer ib,jb
      integer g1b,g2b
c
      length = 0
      do jb=1,atpart2
      do ib=1,atpart2
      do g1b=1,noab+nvab
      do g2b=1,noab+nvab
         length=length+1
      enddo
      enddo
      enddo
      enddo
      if (.not. ma_push_get(mt_int,2*length+1,'noname',l_a_offset,
     1    k_a_offset))
     1  call errquit('offset_btrans2__i1',0,ma_err)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      do jb=1,atpart2
      do ib=1,atpart2
      do g1b=1,noab+nvab
      do g2b=1,noab+nvab
            addr = addr + 1
            int_mb(k_a_offset+addr) = (g2b-1)+(g1b-1)*(noab+nvab)
     1           +(ib-1)*(noab+nvab)*(noab+nvab)
     1           +(jb-1)*(noab+nvab)*(noab+nvab)*atpart2
            int_mb(k_a_offset+length+addr) = size
            size = size + int_mb(k_range+g1b-1)*int_mb(k_range+g2b-1)*
     1             nalength2(ib)*nalength2(jb)
         enddo
         enddo
      enddo
      enddo
      end
c
      subroutine offset_btrans2_i3(l_a_offset,k_a_offset,size,
     1   atpart2,nalength2)
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
      integer ib,jb,kb
      integer g1b
c
      length = 0
      do kb=1,atpart2
      do jb=1,atpart2
      do ib=1,atpart2
      do g1b=1,noab+nvab
         length=length+1
      enddo
      enddo
      enddo
      enddo
      if (.not. ma_push_get(mt_int,2*length+1,'noname',l_a_offset,
     1    k_a_offset))
     1  call errquit('offset_btrans2_i1',0,ma_err)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      do kb=1,atpart2
      do jb=1,atpart2
      do ib=1,atpart2
      do g1b=1,noab+nvab
         addr = addr + 1
         int_mb(k_a_offset+addr) = (g1b-1)+(ib-1)*(noab+nvab)
     1       +(jb-1)*atpart2*(noab+nvab)
     1       +(kb-1)*atpart2*atpart2*(noab+nvab)
         int_mb(k_a_offset+length+addr) = size
         size = size + int_mb(k_range+g1b-1)*nalength2(ib)*
     1          nalength2(jb)*nalength2(kb)
      enddo
      enddo
      enddo
      enddo
      end
c
      subroutine offset_btrans2_i0(l_a_offset,k_a_offset,size,
     1   atpart2,nalength2)
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
      integer ib,jb,kb,lb
      integer g1b
c
      length = 0
      do lb=1,atpart2
      do kb=1,atpart2
      do jb=1,atpart2
      do ib=1,atpart2
         length=length+1
      enddo
      enddo
      enddo
      enddo
      if (.not. ma_push_get(mt_int,2*length+1,'noname',l_a_offset,
     1    k_a_offset))
     1  call errquit('offset_btrans2_i1',0,ma_err)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      do lb=1,atpart2
      do kb=1,atpart2
      do jb=1,atpart2
      do ib=1,atpart2
         addr = addr + 1
         int_mb(k_a_offset+addr) = (ib-1)+(jb-1)*atpart2
     1       +(kb-1)*atpart2*atpart2
     1       +(lb-1)*atpart2*atpart2*atpart2
         int_mb(k_a_offset+length+addr) = size
         size = size + nalength2(lb)*nalength2(ib)*
     1          nalength2(jb)*nalength2(kb)
      enddo
      enddo
      enddo
      enddo
      end
c
      subroutine offset_btrans1_i0(l_a_offset,k_a_offset,size,
     1   atpart2,nalength2)
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
      integer ib,jb
c
      length = 0
      do jb=1,atpart2
      do ib=1,atpart2
         length=length+1
      enddo
      enddo
      if (.not. ma_push_get(mt_int,2*length+1,'noname',l_a_offset,
     1    k_a_offset))
     1  call errquit('offset_btrans1_i0',0,ma_err)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      do jb=1,atpart2
      do ib=1,atpart2
         addr = addr + 1
         int_mb(k_a_offset+addr) = (ib-1)+(jb-1)*atpart2
         int_mb(k_a_offset+length+addr) = size
         size = size + nalength2(ib)*nalength2(jb)
      enddo
      enddo
      end
c
      subroutine scf_1pao(d_scfa,d_scfb,d_scf,k_1pao_offset,size_1p,
     1           atpart2,nalength2)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
      integer d_scfa,d_scfb,d_scf
      integer k_1pao_offset,size_1p
      integer atpart2,nalength2(*)
c
      integer ib,jb,i_ind,j_ind
      integer h1b,h1
      integer i,j,i1,j1
      integer dima,l_a,k_a
      integer dimb,l_b,k_b
      integer dimc,l_c,k_c
      integer dima_sort,dimb_sort,dimc_sort
      integer dim_common
c
      integer nxtask
      integer next
      integer nprocs
      integer count
      external nxtask
c
c     ==================
c     alpha component
c     =================
      nprocs = ga_nnodes( )
      count=0
      next=nxtask(nprocs,1)
c
c      i_ind = 0
      do ib=1,atpart2
c      j_ind = 0
      do jb=1,atpart2
         if (next.eq.count) then
            i_ind=0
            do i=1,ib-1
               i_ind=i_ind+nalength2(i)
            enddo
            j_ind=0
            do j=1,jb-1
               j_ind=j_ind+nalength2(j)
            enddo
c
            dimc=nalength2(ib)*nalength2(jb)
            if (.not.ma_push_get(mt_dbl,dimc,'noname',l_c,k_c))
     1      call errquit('scf_1pao: ma problem',1,ma_err)
            call dfill(dimc,0.d0,dbl_mb(k_c),1)
            do h1b=1,noa
               dima=nalength2(ib)*int_mb(k_range+h1b-1)
               if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1         call errquit('scf_1pao: ma problem',2,ma_err)
               dimb=nalength2(jb)*int_mb(k_range+h1b-1)
               if (.not.ma_push_get(mt_dbl,dimb,'noname',l_b,k_b))
     1            call errquit('scf_1pao: ma problem',3,ma_err)
               do i=1,nalength2(ib)
                  i1=i_ind+i
               do h1=1,int_mb(k_range+h1b-1)
                  dbl_mb(k_a+(h1-1)+(i-1)*int_mb(k_range+h1b-1))=
     1            dbl_mb(k_movecs_sorted+(i1-1)
     1            +(int_mb(k_offset+h1b-1)+h1-1)*nbf)
               enddo
               enddo
               do j=1,nalength2(jb)
                  j1=j_ind+j
               do h1=1,int_mb(k_range+h1b-1)
                  dbl_mb(k_b+(h1-1)+(j-1)*int_mb(k_range+h1b-1))=
     1            dbl_mb(k_movecs_sorted+(j1-1)
     1            +(int_mb(k_offset+h1b-1)+h1-1)*nbf)
               enddo
               enddo
               dima_sort=nalength2(ib)
               dimb_sort=nalength2(jb)
               dim_common=int_mb(k_range+h1b-1)
               call dgemm('t','n',dima_sort,dimb_sort,dim_common,
     1              1.d0,dbl_mb(k_a),dim_common,dbl_mb(k_b),dim_common,
     1              1.d0,dbl_mb(k_c),dima_sort)
               if (.not.ma_pop_stack(l_b))
     1            call errquit('scf_1pao',5,ma_err)
               if (.not.ma_pop_stack(l_a))
     1            call errquit('scf_1pao',6,ma_err)
            enddo
            call put_hash_block(d_scfa,dbl_mb(k_c),dimc,
     1        int_mb(k_1pao_offset),(ib-1)+(jb-1)*atpart2)
            if (.not.ma_pop_stack(l_c))
     1         call errquit('scf_1pao',4,ma_err)
            next=nxtask(nprocs,1)
         endif
         count=count+1
c         j_ind=j_ind+nalength2(jb)
      enddo
c      i_ind=i_ind+nalength2(ib)
      enddo
      next=nxtask(-nprocs,1)
      call ga_sync( )
c     ==============
c     BETA COMPONENT
c     ==============
      nprocs=ga_nnodes()
      count=0
      next=nxtask(nprocs,1)
c      i_ind = 0
      do ib=1,atpart2
c      j_ind = 0
      do jb=1,atpart2
         if (next.eq.count) then
            i_ind=0
            do i=1,ib-1
               i_ind=i_ind+nalength2(i)
            enddo
            j_ind=0
            do i=1,jb-1
               j_ind=j_ind+nalength2(i)
            enddo
c
            dimc=nalength2(ib)*nalength2(jb)
            if (.not.ma_push_get(mt_dbl,dimc,'noname',l_c,k_c))
     1         call errquit('scf_1pao: ma problem',1,ma_err)
            call dfill(dimc,0.d0,dbl_mb(k_c),1)
            do h1b=noa+1,noab
               dima=nalength2(ib)*int_mb(k_range+h1b-1)
               if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1            call errquit('scf_1pao: ma problem',2,ma_err)
               dimb=nalength2(jb)*int_mb(k_range+h1b-1)
               if (.not.ma_push_get(mt_dbl,dimb,'noname',l_b,k_b))
     1            call errquit('scf_1pao: ma problem',3,ma_err)
               do i=1,nalength2(ib)
                  i1=i_ind+i
               do h1=1,int_mb(k_range+h1b-1)
                  dbl_mb(k_a+(h1-1)+(i-1)*int_mb(k_range+h1b-1))=
     1            dbl_mb(k_movecs_sorted+(i1-1)
     1           +(int_mb(k_offset+h1b-1)+h1-1)*nbf)
               enddo
               enddo
               do j=1,nalength2(jb)
                  j1=j_ind+j
               do h1=1,int_mb(k_range+h1b-1)
                  dbl_mb(k_b+(h1-1)+(j-1)*int_mb(k_range+h1b-1))=
     1            dbl_mb(k_movecs_sorted+(j1-1)
     1            +(int_mb(k_offset+h1b-1)+h1-1)*nbf)
               enddo
               enddo
               dima_sort=nalength2(ib)
               dimb_sort=nalength2(jb)
               dim_common=int_mb(k_range+h1b-1)
               call dgemm('t','n',dima_sort,dimb_sort,dim_common,
     1              1.d0,dbl_mb(k_a),dim_common,dbl_mb(k_b),dim_common,
     1              1.d0,dbl_mb(k_c),dima_sort)
               if (.not.ma_pop_stack(l_b))
     1            call errquit('scf_1pao',5,ma_err)
               if (.not.ma_pop_stack(l_a))
     1            call errquit('scf_1pao',6,ma_err)
            enddo
            call put_hash_block(d_scfb,dbl_mb(k_c),dimc,
     1           int_mb(k_1pao_offset),(ib-1)+(jb-1)*atpart2)
            if (.not.ma_pop_stack(l_c))
     1         call errquit('scf_1pao',4,ma_err)
            next=nxtask(nprocs,1)
         endif
         count=count+1
c         j_ind=j_ind+nalength2(jb)
      enddo
c      i_ind=i_ind+nalength2(ib)
      enddo
      next=nxtask(-nprocs,1)
      call ga_sync( )
c     ====================
c     ALPHA+BETA COMPONENT
c     ====================
      nprocs=ga_nnodes()
      count=0
      next=nxtask(nprocs,1)
      do ib=1,atpart2
      do jb=1,atpart2
         if (next.eq.count) then
            dima=nalength2(ib)*nalength2(jb)
            if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1         call errquit('scf_1pao: ma problem',7,ma_err)
            if (.not.ma_push_get(mt_dbl,dima,'noname',l_b,k_b))
     1         call errquit('scf_1pao: ma problem',8,ma_err)
            if (.not.ma_push_get(mt_dbl,dima,'noname',l_c,k_c))
     1         call errquit('scf_1pao: ma problem',9,ma_err)
            call get_hash_block(d_scfa,dbl_mb(k_a),dima,
     1           int_mb(k_1pao_offset),(jb-1)+(ib-1)*atpart2)
            call get_hash_block(d_scfb,dbl_mb(k_b),dima,
     1           int_mb(k_1pao_offset),(jb-1)+(ib-1)*atpart2)
            do i=1,dima
               dbl_mb(k_c+i-1)=dbl_mb(k_a+i-1)+dbl_mb(k_b+i-1)
            enddo
            call put_hash_block(d_scf,dbl_mb(k_c),dima,
     1           int_mb(k_1pao_offset),(jb-1)+(ib-1)*atpart2)
            if (.not.ma_pop_stack(l_c))
     1         call errquit('scf_1pao',10,ma_err)
            if (.not.ma_pop_stack(l_b))
     1         call errquit('scf_1pao',11,ma_err)
            if (.not.ma_pop_stack(l_a))
     1         call errquit('scf_1pao',12,ma_err)
            next=nxtask(nprocs,1)
         endif
         count=count+1
      enddo
      enddo
      next=nxtask(-nprocs,1)
      call ga_sync( )
      return
      end
c
      subroutine cc_1pao(d_1pdm,k_1pdm_offset,d_cca,d_ccb,d_cc,
     1     k_1pao_offset,size_1p,atpart2,nalength2)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
      integer d_1pdm,k_1pdm_offset
      integer d_cca,d_ccb,d_cc,k_1pao_offset,size_1p
      integer atpart2,nalength2(*)
c
      integer d_i1,l_i1_offset,k_i1_offset,size_i1
      integer ib,jb,i
      integer dima,l_a,k_a,l_b,k_b,l_c,k_c
      character*255 filename
c
      integer nxtask
      integer next
      integer nprocs
      integer count
      external nxtask
c 
      call offset_btrans1_i1(l_i1_offset,k_i1_offset,size_i1,atpart2,
     1     nalength2)
      call tce_filename('btrans1_i1',filename)
      call createfile(filename,d_i1,size_i1)
c     =========
c     A-A (H-H)
c     =========
      call btrans1_i1_2(d_1pdm,k_1pdm_offset,d_i1,k_i1_offset,
     1     atpart2,nalength2,
     1     1,noa,1,noa)
      call btrans1_i0_2(d_i1,k_i1_offset,d_cca,k_1pao_offset,
     1     atpart2,nalength2,1,noa)
c      print *,'cc_1pao A-A hh is fine'
c     =========
c     A-A (H-P)
c     =========
      call tce_zero(d_i1,size_i1)
      call btrans1_i1_2(d_1pdm,k_1pdm_offset,d_i1,k_i1_offset,
     1     atpart2,nalength2,
     1     1,noa,noab+1,noab+nva)
      call btrans1_i0_2(d_i1,k_i1_offset,d_cca,k_1pao_offset,
     1     atpart2,nalength2,1,noa)
c      print *,'aa hp is fine'
c     =========
c     A-A (P-H)
c     =========
c      print *,'=============aa ph starts=================='
      call tce_zero(d_i1,size_i1)
      call btrans1_i1_2(d_1pdm,k_1pdm_offset,d_i1,k_i1_offset,
     1     atpart2,nalength2,
     1     noab+1,noab+nva,1,noa)
c      print *,'aa ph i1 is fine'
      call btrans1_i0_2(d_i1,k_i1_offset,d_cca,k_1pao_offset,
     1     atpart2,nalength2,noab+1,noab+nva)
c      print *,'aa ph is fine'
c     ===========
c      A-A (P-P)
c     ===========
      call tce_zero(d_i1,size_i1)
      call btrans1_i1_2(d_1pdm,k_1pdm_offset,d_i1,k_i1_offset,
     1     atpart2,nalength2,
     1     noab+1,noab+nva,noab+1,noab+nva)
      call btrans1_i0_2(d_i1,k_i1_offset,d_cca,k_1pao_offset,
     1     atpart2,nalength2,noab+1,noab+nva)
c      print *,'aa pp is fine'
c
c     =========================================================
      if (restricted) then
         call copyfile(d_cca,d_ccb,size_1p)
      endif 
c
c     =========
c     B-B (H-H)
c     =========
      if(.not.restricted) then
      call tce_zero(d_i1,size_i1)
      call btrans1_i1_2(d_1pdm,k_1pdm_offset,d_i1,k_i1_offset,
     1     atpart2,nalength2,
     1     noa+1,noab,noa+1,noab)
      call btrans1_i0_2(d_i1,k_i1_offset,d_ccb,k_1pao_offset,
     1     atpart2,nalength2,noa+1,noab)
c     =========
c     B-B (H-P)
c     =========
      call tce_zero(d_i1,size_i1)
      call btrans1_i1_2(d_1pdm,k_1pdm_offset,d_i1,k_i1_offset,
     1     atpart2,nalength2,
     1     noa+1,noab,noab+nva+1,noab+nvab)
      call btrans1_i0_2(d_i1,k_i1_offset,d_ccb,k_1pao_offset,
     1     atpart2,nalength2,noa+1,noab)
c     =========
c     B-B (P-H)
c     =========
      call tce_zero(d_i1,size_i1)
      call btrans1_i1_2(d_1pdm,k_1pdm_offset,d_i1,k_i1_offset,
     1     atpart2,nalength2,
     1     noab+nva+1,noab+nvab,noa+1,noab)
      call btrans1_i0_2(d_i1,k_i1_offset,d_ccb,k_1pao_offset,
     1     atpart2,nalength2,noab+nva+1,noab+nvab)
c     ===========
c      B-B (P-P)
c     ===========
      call tce_zero(d_i1,size_i1)
      call btrans1_i1_2(d_1pdm,k_1pdm_offset,d_i1,k_i1_offset,
     1     atpart2,nalength2,
     1     noab+nva+1,noab+nvab,noab+nva+1,noab+nvab)
      call btrans1_i0_2(d_i1,k_i1_offset,d_ccb,k_1pao_offset,
     1     atpart2,nalength2,noab+nva+1,noab+nvab)
      endif
c
c    ================================================================
c
      call deletefile(d_i1)
      if (.not. ma_pop_stack(l_i1_offset))
     1    call errquit('btrans1: ma problem',1,ma_err)
c
c     ================================================================
c
c     =============
c     ALPHA + BETA
c     =============
      nprocs=ga_nnodes()
      count=0
      next=nxtask(nprocs,1)
      do ib=1,atpart2
      do jb=1,atpart2
         if (next.eq.count) then
            dima=nalength2(ib)*nalength2(jb)
            if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1      call errquit('cc_1pao: ma problem',2,ma_err)
            if (.not.ma_push_get(mt_dbl,dima,'noname',l_b,k_b))
     1      call errquit('cc_1pao: ma problem',2,ma_err)
            if (.not.ma_push_get(mt_dbl,dima,'noname',l_c,k_c))
     1      call errquit('cc_1pao: ma problem',2,ma_err)
            call get_hash_block(d_cca,dbl_mb(k_a),dima,
     1           int_mb(k_1pao_offset),(jb-1)+(ib-1)*atpart2)
            call get_hash_block(d_ccb,dbl_mb(k_b),dima,
     1           int_mb(k_1pao_offset),(jb-1)+(ib-1)*atpart2)
            do i=1,dima
               dbl_mb(k_c+i-1)=dbl_mb(k_a+i-1)+dbl_mb(k_b+i-1)
            enddo
            call put_hash_block(d_cc,dbl_mb(k_c),dima,
     1           int_mb(k_1pao_offset),(jb-1)+(ib-1)*atpart2)
            if (.not.ma_pop_stack(l_c))
     1         call errquit('cc_1pao: ma problem',5,ma_err)
            if (.not.ma_pop_stack(l_b))
     1         call errquit('cc_1pao: ma problem',6,ma_err)
            if (.not.ma_pop_stack(l_a))
     1         call errquit('cc_1pao: ma problem',7,ma_err)
            next=nxtask(nprocs,1)
         endif
         count=count+1
      enddo
      enddo
      next=nxtask(-nprocs,1)
      call ga_sync( )
      return
      end
c
      subroutine btrans1_i1_2(d_a,k_a_offset,d_c,k_c_offset,
     1           atpart2,nalength2,
     1           g1b_lo,g1b_hi,g2b_lo,g2b_hi)
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
      integer g1b_lo,g1b_hi
      integer g2b_lo,g2b_hi
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
      double precision res
c
      integer nxtask
      integer next
      integer nprocs
      integer count
      external nxtask
c
      nprocs=ga_nnodes( )
      count=0
      next=nxtask(nprocs,1)
c
c      j_ind=0
      do jb=1,atpart2
      do g1b=g1b_lo,g1b_hi
         if (next.eq.count) then
            j_ind=0
            do j=1,jb-1
               j_ind=j_ind+nalength2(j)
            enddo
c
            dimc=int_mb(k_range+g1b-1)*nalength2(jb)
            if (.not.ma_push_get(mt_dbl,dimc,'noname',l_c,k_c))
     1      call errquit('btrans1_i1: ma problem',1,ma_err)
            call dfill(dimc,0.d0,dbl_mb(k_c),1)
            do g2b=g2b_lo,g2b_hi
               if ((.not.restricted).or.(int_mb(k_spin+g2b-1)
     1             +int_mb(k_spin+g1b-1)).ne.4) then
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
     1               dbl_mb(k_movecs_sorted+
     1               (j1-1)+(int_mb(k_offset+g2b-1)+g2-1)*nbf)
                  enddo
                  enddo
                  call dgemm('t','n',dima_sort,dimb_sort,dim_common,
     1                 1.d0,dbl_mb(k_a),dim_common,dbl_mb(k_b),
     1                 dim_common,1.d0,dbl_mb(k_c),dima_sort)
                  if (.not.ma_pop_stack(l_b))
     1               call errquit('btrans1_i1:ma problem',4,ma_err)
                  if (.not.ma_pop_stack(l_a))
     1               call errquit('btrans1_i1: ma problem',5,ma_err)
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
         count = count + 1
      enddo
c      j_ind=j_ind+nalength2(jb)
      enddo
      next = nxtask(-nprocs,1)
      call ga_sync( )
      return
      end
c
      subroutine btrans1_i0_2(d_a,k_a_offset,d_b,k_b_offset,
     1           atpart2,nalength2,
     1           g1b_lo,g1b_hi)
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
      integer d_b,k_b_offset
      integer g1b_lo,g1b_hi
      integer atpart2
      integer nalength2(*)
      logical alpha2
c
      integer k_b,l_b,dimb,dimb_sort
      integer k_a,l_a,dima,dima_sort
      integer k_c,l_c,dimc,dimc_sort
      integer k_c0,l_c0
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
      nprocs=ga_nnodes( )
      count=0
      next=nxtask(nprocs,1)
c      i_ind=0
      do ib=1,atpart2
c      j_ind=0
      do jb=1,atpart2
         if (next.eq.count) then
            i_ind=0
            do i=1,ib-1
               i_ind=i_ind+nalength2(i)
            enddo
            j_ind=0
            do i=1,jb-1
               j_ind=j_ind+nalength2(j)
            enddo
c
            dimc=nalength2(ib)*nalength2(jb)
            if (.not.ma_push_get(mt_dbl,dimc,'noname',l_c0,k_c0))
     1      call errquit('btrans1_i0_2: ma problem',11,ma_err)
            call get_hash_block(d_b,dbl_mb(k_c0),dimc,
     1           int_mb(k_b_offset),(jb-1)+(ib-1)*atpart2)
            if (.not.ma_push_get(mt_dbl,dimc,'noname',l_c,k_c))
     1      call errquit('btrans1_i0: ma problem',1,ma_err)
            call dfill(dimc,0.d0,dbl_mb(k_c),1)
            do g1b=g1b_lo,g1b_hi
               dima=int_mb(k_range+g1b-1)*nalength2(jb)
               dima_sort = nalength2(jb)
               if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1         call errquit('btrans1_i0: ma problem',2,ma_err)
               dimb=int_mb(k_range+g1b-1)*nalength2(ib)
               dimb_sort = nalength2(ib)
               dim_common=int_mb(k_range+g1b-1)
               if (.not.ma_push_get(mt_dbl,dimb,'noname',l_b,k_b))
     1         call errquit('btrans1_i0: ma problem',3,ma_err)
               call get_hash_block(d_a,dbl_mb(k_a),dima,
     1              int_mb(k_a_offset),(g1b-1)+(jb-1)*(noab+nvab))
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
     1         call errquit('btrans1_i0: ma problem',4,ma_err)
               if(.not.ma_pop_stack(l_a))
     1         call errquit('btrans1_hh_i0: ma problem',5,ma_err)
            enddo
            do i=1,dimc
               dbl_mb(k_c+i-1)=dbl_mb(k_c+i-1)+dbl_mb(k_c0+i-1)
            enddo
            call put_hash_block(d_b,dbl_mb(k_c),dimc,
     1           int_mb(k_b_offset),(jb-1)+(ib-1)*atpart2)
            if (.not.ma_pop_stack(l_c))
     1         call errquit('btrans1_hh_i0: ma problem',6,ma_err)
            if (.not.ma_pop_stack(l_c0))
     1         call errquit('btrans1_hh_i0: ma problem',7,ma_err)
            next=nxtask(nprocs,1)
         endif
         count=count+1
c         j_ind=j_ind+nalength2(jb)
      enddo
c      i_ind=i_ind+nalength2(ib)
      enddo
      next=nxtask(-nprocs,1)
      call ga_sync( )
      return
      end
c $Id$
