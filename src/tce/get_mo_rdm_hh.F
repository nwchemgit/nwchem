      subroutine get_mo_rdm_hh(d_r1,k_r1_offset,k_rdm,denominator)
c
c written by Peng-Dong Fan 07/27/2007
c
c Prints the Hole-Hole block density matrix tensor
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_r1
      integer h1b
      integer h2b
      integer h1
      integer h2
      integer k_r1_offset
      integer size
      integer l_r1,k_r1
      integer i
      logical nodezero
      character*4 irrep_h1,irrep_h2
      character*5 spinnames(2),spin_h1,spin_h2
      data spinnames/'alpha','beta '/
c->pengdong.fan
      integer k_rdm
      integer ih1,ih2,ihh1,ihh2,ind
      integer nh,ind1,ind2,ind3
      double precision denominator
c<-pengdong.fan
c
      nodezero = (ga_nodeid().eq.0)
c      if (nodezero) write(LuOut,9000) 'hole-hole block'
      nh=nocc(1)+nocc(ipol)
      do h1b = 1,noab
        do h2b = 1,noab
          if (int_mb(k_spin+h1b-1) .eq. int_mb(k_spin+h2b-1)) then 
          if ((.not.restricted).or.(int_mb(k_spin+h1b-1)
     1      +int_mb(k_spin+h2b-1).ne.4)) then
          if (ieor(int_mb(k_sym+h1b-1),int_mb(k_sym+h2b-1))
     1      .eq.irrep_e) then
            size = int_mb(k_range+h1b-1) * int_mb(k_range+h2b-1)
            if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1        call errquit('tce_print_hh: MA problem',0,MA_ERR)
            call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1        int_mb(k_r1_offset),((h1b-1)*noab+h2b-1))
            i = 0
            do h1 = 1,int_mb(k_range+h1b-1)
              do h2 = 1,int_mb(k_range+h2b-1)
                i = i + 1
c                call sym_irrepname(geom,int_mb(k_sym+h1b-1)+1,
c     1               irrep_h1)
c                call sym_irrepname(geom,int_mb(k_sym+h2b-1)+1,
c     1               irrep_h2)
c                spin_h1 = spinnames(int_mb(k_spin+h1b-1))
c                spin_h2 = spinnames(int_mb(k_spin+h2b-1))
                ih1=int_mb(k_offset+h1b-1)+h1
                ih2=int_mb(k_offset+h2b-1)+h2
                ihh1=int_mb(k_mo_index+ih1-1)
                ihh2=int_mb(k_mo_index+ih2-1)
                ind=(ihh1-1)+(ihh2-1)*nh
                dbl_mb(k_rdm+ind)=dbl_mb(k_r1+i-1)/denominator
c                if (nodezero) write(LuOut,9020) 
c     1            int_mb(k_offset+h1b-1)+h1,irrep_h1,spin_h1,
c     2            int_mb(k_offset+h2b-1)+h2,irrep_h2,spin_h2,
c     3            dbl_mb(k_r1+i-1),ihh1,ihh2
              enddo
            enddo
            if (.not.ma_pop_stack(l_r1))
     1        call errquit('tce_print_hh: MA problem',1,MA_ERR)
          endif
          endif
          endif
        enddo
      enddo
      do i=1,nfc(1)+nfc(ipol)
         ind=(i-1)*nh+(i-1)
         dbl_mb(k_rdm+ind) = 1.d0
      enddo
      do i=nfc(1)+nfc(ipol)+1,nocc(1)+nocc(ipol)
         ind=(i-1)*nh+(i-1)
         dbl_mb(k_rdm+ind) = dbl_mb(k_rdm+ind) + 1.d0
      enddo
      if(restricted) then
        do ih1=1,nocc(ipol)
           ind1=(2*ih1-1)
           ind2=(2*ih1-1-1) 
           do ih2=1,nocc(ipol)
              ind = ind1 + (2*ih2-1)*nh
              ind3= ind2 + (2*ih2-1-1)*nh
              dbl_mb(k_rdm+ind) = dbl_mb(k_rdm+ind3)
           enddo
        enddo
      endif
c      do ih1=1,nh
c         ind1=(ih1-1)
c         do ih2=1,nh
c            ind=ind1+(ih2-1)*nh
c            if(nodezero) then
c              if(dabs(dbl_mb(k_rdm+ind)).gt.1.d-10) 
c     1          write(LuOut, 9010) dbl_mb(k_rdm+ind),ih1,ih2
c            endif
c         enddo
c      enddo
      call ga_sync()
 9000 format(/,1x,a)
 9010 format(f30.16,i5,i5)
 9020 format(1x,i3,a4,'(',a5,') --- ',i3,a4,'(',a5,')',f20.5,i5,i5)
      return
      end
c $Id$
