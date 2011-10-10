      subroutine sym_2pdm_ao(d_r2,k_r2_offset,nshblocks,
     1           bfglo,bfghi)
c
c Prints the particle-hole block of the 1-e reduced density matrix
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
      integer d_r2, k_r2_offset
      integer nshblocks, bfglo(*), bfghi(*)
      integer ib,jb,kb,lb,idim,jdim,kdim,ldim
      integer i1,j1,k1,l1, i
      integer k_r2_1,l_r2_1
      integer k_r2_2,l_r2_2
      integer k_r2_3,l_r2_3
      integer k_r2_4,l_r2_4
      integer k_r2_5,l_r2_5
      integer k_r2_6,l_r2_6
      integer k_r2_7,l_r2_7
      integer k_r2_8,l_r2_8
      integer l_a,k_a
      integer size
      integer ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8
      logical nodezero
      double precision res
c
      integer nxtask
      integer next
      integer count
      integer nprocs
      external nxtask
c

c      nodezero = (ga_nodeid().eq.0)
      nprocs=ga_nnodes( )
      count=0
      next=nxtask(nprocs,1)
c
      do ib=1,nshblocks
c         idim=bfghi(ib)-bfglo(ib)+1
c      do jb=1,nshblocks
      do jb=ib,nshblocks
c         jdim=bfghi(jb)-bfglo(jb)+1
c      do kb=1,nshblocks
      do kb=ib,nshblocks
c         kdim=bfghi(kb)-bfglo(kb)+1
c      do lb=1,nshblocks
      do lb=jb,nshblocks
         if(next.eq.count) then
         if(ib.eq.jb.and.lb.lt.kb) then
         else
c           write(6,'(a,4i4)') 'ib,jb,kb,lb=',ib,jb,kb,lb
c           write(6,'(a,4i4)') 'kb,jb,ib,lb=',kb,jb,ib,lb
c           write(6,'(a,4i4)') 'ib,lb,kb,jb=',ib,lb,kb,jb
c           write(6,'(a,4i4)') 'kb,lb,ib,jb=',kb,lb,ib,jb
c           write(6,'(a,4i4)') 'jb,ib,lb,kb=',jb,ib,lb,kb
c           write(6,'(a,4i4)') 'lb,ib,jb,kb=',lb,ib,jb,kb
c           write(6,'(a,4i4)') 'jb,kb,lb,ib=',jb,kb,lb,ib
c           write(6,'(a,4i4)') 'lb,kb,jb,ib=',lb,kb,jb,ib
c
         idim=bfghi(ib)-bfglo(ib)+1
         jdim=bfghi(jb)-bfglo(jb)+1
         kdim=bfghi(kb)-bfglo(kb)+1
         ldim=bfghi(lb)-bfglo(lb)+1
c
         size = idim*jdim*kdim*ldim
         if (.not.ma_push_get(mt_dbl,size,'r2_1',l_r2_1,k_r2_1))
     1   call errquit('sym_2pdm_ao: ma problem',0,ma_err)
         if (.not.ma_push_get(mt_dbl,size,'r2_2',l_r2_2,k_r2_2))
     1   call errquit('sym_2pdm_ao: ma problem',1,ma_err)
         if (.not.ma_push_get(mt_dbl,size,'r2_3',l_r2_3,k_r2_3))
     1   call errquit('sym_2pdm_ao: ma problem',2,ma_err)
         if (.not.ma_push_get(mt_dbl,size,'r2_4',l_r2_4,k_r2_4))
     1   call errquit('sym_2pdm_ao: ma problem',3,ma_err)
         if (.not.ma_push_get(mt_dbl,size,'r2_5',l_r2_5,k_r2_5))
     1   call errquit('sym_2pdm_ao: ma problem',4,ma_err)
         if (.not.ma_push_get(mt_dbl,size,'r2_6',l_r2_6,k_r2_6))
     1   call errquit('sym_2pdm_ao: ma problem',5,ma_err)
         if (.not.ma_push_get(mt_dbl,size,'r2_7',l_r2_7,k_r2_7))
     1   call errquit('sym_2pdm_ao: ma problem',6,ma_err)
         if (.not.ma_push_get(mt_dbl,size,'r2_8',l_r2_8,k_r2_8))
     1   call errquit('sym_2pdm_ao: ma problem',7,ma_err)
c
         call dfill(size,0.d0,dbl_mb(k_r2_1),1)
         call dfill(size,0.d0,dbl_mb(k_r2_2),1)
         call dfill(size,0.d0,dbl_mb(k_r2_3),1)
         call dfill(size,0.d0,dbl_mb(k_r2_4),1)
         call dfill(size,0.d0,dbl_mb(k_r2_5),1)
         call dfill(size,0.d0,dbl_mb(k_r2_6),1)
         call dfill(size,0.d0,dbl_mb(k_r2_7),1)
         call dfill(size,0.d0,dbl_mb(k_r2_8),1)
c
c        (i,j,k,l)
         call get_hash_block(d_r2,dbl_mb(k_r2_1),size,
     1        int_mb(k_r2_offset),(lb-1)+(kb-1)*nshblocks
     1        +(jb-1)*nshblocks*nshblocks
     1        +(ib-1)*nshblocks*nshblocks*nshblocks)
c        (k,j,i,l)
         call get_hash_block(d_r2,dbl_mb(k_r2_2),size,
     1        int_mb(k_r2_offset),(lb-1)+(ib-1)*nshblocks
     1        +(jb-1)*nshblocks*nshblocks
     1        +(kb-1)*nshblocks*nshblocks*nshblocks)
c        (i,l,k,j)
         call get_hash_block(d_r2,dbl_mb(k_r2_3),size,
     1        int_mb(k_r2_offset),(jb-1)+(kb-1)*nshblocks
     1        +(lb-1)*nshblocks*nshblocks
     1        +(ib-1)*nshblocks*nshblocks*nshblocks)
c        (k,l,i,j)
         call get_hash_block(d_r2,dbl_mb(k_r2_4),size,
     1        int_mb(k_r2_offset),(jb-1)+(ib-1)*nshblocks
     1        +(lb-1)*nshblocks*nshblocks
     1        +(kb-1)*nshblocks*nshblocks*nshblocks)
c        (j,i,l,k)
         call get_hash_block(d_r2,dbl_mb(k_r2_5),size,
     1        int_mb(k_r2_offset),(kb-1)+(lb-1)*nshblocks
     1        +(ib-1)*nshblocks*nshblocks
     1        +(jb-1)*nshblocks*nshblocks*nshblocks)
c        (l,i,j,k)
         call get_hash_block(d_r2,dbl_mb(k_r2_6),size,
     1        int_mb(k_r2_offset),(kb-1)+(jb-1)*nshblocks
     1        +(ib-1)*nshblocks*nshblocks
     1        +(lb-1)*nshblocks*nshblocks*nshblocks)
c        (j,k,l,i)
         call get_hash_block(d_r2,dbl_mb(k_r2_7),size,
     1        int_mb(k_r2_offset),(ib-1)+(lb-1)*nshblocks
     1        +(kb-1)*nshblocks*nshblocks
     1        +(jb-1)*nshblocks*nshblocks*nshblocks)
c        (l,k,j,i)
         call get_hash_block(d_r2,dbl_mb(k_r2_8),size,
     1        int_mb(k_r2_offset),(ib-1)+(jb-1)*nshblocks
     1        +(kb-1)*nshblocks*nshblocks
     1        +(lb-1)*nshblocks*nshblocks*nshblocks)
c
            do i1 = 1, idim
            do j1 = 1, jdim
            do k1 = 1, kdim
            do l1 = 1, ldim
               ind1 = (l1-1)+(k1-1)*ldim+(j1-1)*ldim*kdim
     1               +(i1-1)*ldim*kdim*jdim
               ind2 = (l1-1)+(i1-1)*ldim+(j1-1)*ldim*idim
     1               +(k1-1)*ldim*idim*jdim
               ind3 = (j1-1)+(k1-1)*jdim+(l1-1)*jdim*kdim
     1               +(i1-1)*jdim*kdim*ldim
               ind4 = (j1-1)+(i1-1)*jdim+(l1-1)*jdim*idim
     1               +(k1-1)*jdim*idim*ldim
               ind5 = (k1-1)+(l1-1)*kdim+(i1-1)*kdim*ldim
     1               +(j1-1)*kdim*ldim*idim
               ind6 = (k1-1)+(j1-1)*kdim+(i1-1)*kdim*jdim
     1               +(l1-1)*kdim*jdim*idim
               ind7 = (i1-1)+(l1-1)*idim+(k1-1)*idim*ldim
     1               +(j1-1)*idim*ldim*kdim
               ind8 = (i1-1)+(j1-1)*idim+(k1-1)*idim*jdim
     1               +(l1-1)*idim*jdim*kdim
               res = (1.d0/8.d0)*(
     1               dbl_mb(k_r2_1+ind1)+dbl_mb(k_r2_2+ind2)+
     1               dbl_mb(k_r2_3+ind3)+dbl_mb(k_r2_4+ind4)+
     1               dbl_mb(k_r2_5+ind5)+dbl_mb(k_r2_6+ind6)+
     1               dbl_mb(k_r2_7+ind7)+dbl_mb(k_r2_8+ind8) )
c
               dbl_mb(k_r2_1+ind1)=res
               dbl_mb(k_r2_2+ind2)=res
               dbl_mb(k_r2_3+ind3)=res
               dbl_mb(k_r2_4+ind4)=res
               dbl_mb(k_r2_5+ind5)=res
               dbl_mb(k_r2_6+ind6)=res
               dbl_mb(k_r2_7+ind7)=res
               dbl_mb(k_r2_8+ind8)=res
c 
            enddo
            enddo
            enddo
            enddo
c        (i,j,k,l)
         call put_hash_block(d_r2,dbl_mb(k_r2_1),size,
     1        int_mb(k_r2_offset),(lb-1)+(kb-1)*nshblocks
     1        +(jb-1)*nshblocks*nshblocks
     1        +(ib-1)*nshblocks*nshblocks*nshblocks)
c        (k,j,i,l)
         call put_hash_block(d_r2,dbl_mb(k_r2_2),size,
     1        int_mb(k_r2_offset),(lb-1)+(ib-1)*nshblocks
     1        +(jb-1)*nshblocks*nshblocks
     1        +(kb-1)*nshblocks*nshblocks*nshblocks)
c        (i,l,k,j)
         call put_hash_block(d_r2,dbl_mb(k_r2_3),size,
     1        int_mb(k_r2_offset),(jb-1)+(kb-1)*nshblocks
     1        +(lb-1)*nshblocks*nshblocks
     1        +(ib-1)*nshblocks*nshblocks*nshblocks)
c        (k,l,i,j)
         call put_hash_block(d_r2,dbl_mb(k_r2_4),size,
     1        int_mb(k_r2_offset),(jb-1)+(ib-1)*nshblocks
     1        +(lb-1)*nshblocks*nshblocks
     1        +(kb-1)*nshblocks*nshblocks*nshblocks)
c        (j,i,l,k)
         call put_hash_block(d_r2,dbl_mb(k_r2_5),size,
     1        int_mb(k_r2_offset),(kb-1)+(lb-1)*nshblocks
     1        +(ib-1)*nshblocks*nshblocks
     1        +(jb-1)*nshblocks*nshblocks*nshblocks)
c        (l,i,j,k)
         call put_hash_block(d_r2,dbl_mb(k_r2_6),size,
     1        int_mb(k_r2_offset),(kb-1)+(jb-1)*nshblocks
     1        +(ib-1)*nshblocks*nshblocks
     1        +(lb-1)*nshblocks*nshblocks*nshblocks)
c        (j,k,l,i)
         call put_hash_block(d_r2,dbl_mb(k_r2_7),size,
     1        int_mb(k_r2_offset),(ib-1)+(lb-1)*nshblocks
     1        +(kb-1)*nshblocks*nshblocks
     1        +(jb-1)*nshblocks*nshblocks*nshblocks)
c        (l,k,j,i)
         call put_hash_block(d_r2,dbl_mb(k_r2_8),size,
     1        int_mb(k_r2_offset),(ib-1)+(jb-1)*nshblocks
     1        +(kb-1)*nshblocks*nshblocks
     1        +(lb-1)*nshblocks*nshblocks*nshblocks)
c
            if (.not.ma_pop_stack(l_r2_8))
     1        call errquit('sym_2pdm_ao: MA problem',9,MA_ERR)
c
            if (.not.ma_pop_stack(l_r2_7))
     1        call errquit('sym_2pdm_ao: MA problem',10,MA_ERR)
c
            if (.not.ma_pop_stack(l_r2_6))
     1        call errquit('sym_2pdm_ao: MA problem',11,MA_ERR)
c
            if (.not.ma_pop_stack(l_r2_5))
     1        call errquit('sym_2pdm_ao: MA problem',12,MA_ERR)
c
            if (.not.ma_pop_stack(l_r2_4))
     1        call errquit('sym_2pdm_ao: MA problem',13,MA_ERR)
c
            if (.not.ma_pop_stack(l_r2_3))
     1        call errquit('sym_2pdm_ao: MA problem',14,MA_ERR)
c
            if (.not.ma_pop_stack(l_r2_2))
     1        call errquit('sym_2pdm_ao: MA problem',15,MA_ERR)
c
            if (.not.ma_pop_stack(l_r2_1))
     1        call errquit('sym_2pdm_ao: MA problem',16,MA_ERR)
            endif
            next=nxtask(nprocs,1)
            endif
            count=count+1
c
      enddo
      enddo
      enddo
      enddo
      next=nxtask(-nprocs,1)
      call ga_sync( )
c
c     check if it is symmetrized
c
c      do ib=1,nshblocks
c         idim=bfghi(ib)-bfglo(ib)+1
c      do jb=1,nshblocks
c         jdim=bfghi(jb)-bfglo(jb)+1
c      do kb=1,nshblocks
c         kdim=bfghi(kb)-bfglo(kb)+1
c      do lb=1,nshblocks
c         ldim=bfghi(lb)-bfglo(lb)+1
c         size = idim*jdim*kdim*ldim
c         if (.not.ma_push_get(mt_dbl,size,'r2_1',l_r2_1,k_r2_1))
c     1   call errquit('sym_2pdm_ao: ma problem',0,ma_err)
cc        (i,j,k,l)
c         call get_hash_block(d_r2,dbl_mb(k_r2_1),size,
c     1        int_mb(k_r2_offset),(lb-1)+(kb-1)*nshblocks
c     1        +(jb-1)*nshblocks*nshblocks
c     1        +(ib-1)*nshblocks*nshblocks*nshblocks)
cc
c         i=0
c         do i1 = 1, idim
c         do j1 = 1, jdim
c         do k1 = 1, kdim
c         do l1 = 1, ldim
c            i=i+1
c            write(6,'(a,4i5,a,f20.16)') 'sym checking 2PDM AO',
c     1               bfglo(ib)+i1-1, bfglo(jb)+j1-1,
c     1               bfglo(kb)+k1-1, bfglo(lb)+l1-1,'=',
c     1               dbl_mb(k_r2_1+i-1)
cc
c         enddo
c         enddo
c         enddo
c         enddo
c         if (.not.ma_pop_stack(l_r2_1))
c     1   call errquit('sym_2pdm_ao: MA problem',16,MA_ERR)
cc
c      enddo
c      enddo
c      enddo
c      enddo
cc
      return
      end
c $Id$
