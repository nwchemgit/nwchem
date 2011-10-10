      subroutine get_mo_rdm_hp(d_r1,k_r1_offset,k_rdm,denominator)
c
c
c Prints the Hole-Particle block of 1-e reduced density matrix
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
      integer p2b
      integer h1
      integer p2
      integer k_r1_offset
      integer size
      integer l_r1,k_r1
      integer i
      logical nodezero
      character*4 irrep_h,irrep_p
      character*5 spinnames(2),spin_h,spin_p
      data spinnames/'alpha','beta '/
c->pengdong.fan
      integer k_rdm
      integer ih1,ip2,ihh1,ipp2,ind,np,nh,ind1,ind2,ind3
      double precision denominator
c<-pengdong.fan
     
c
      nodezero = (ga_nodeid().eq.0)
c      if (nodezero) write(LuOut,9000) 'hole-particle block'
      nh=nocc(1)+nocc(ipol)
      np=nmo(1)+nmo(ipol)-nocc(1)-nocc(ipol)
      do h1b = 1,noab
        do p2b = 1+noab,noab+nvab
          if (int_mb(k_spin+h1b-1) .eq. int_mb(k_spin+p2b-1)) then 
          if ((.not.restricted).or.(int_mb(k_spin+h1b-1)
     1      +int_mb(k_spin+p2b-1).ne.4)) then
          if (ieor(int_mb(k_sym+h1b-1),int_mb(k_sym+p2b-1))
     1      .eq.irrep_e) then
            size = int_mb(k_range+h1b-1) * int_mb(k_range+p2b-1)
            if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1        call errquit('tce_print_hp: MA problem',0,MA_ERR)
            call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1        int_mb(k_r1_offset),((h1b-1)*nvab+p2b-noab-1))
            i = 0
            do h1 = 1,int_mb(k_range+h1b-1)
              do p2 = 1,int_mb(k_range+p2b-1)
                i = i + 1
c                call sym_irrepname(geom,int_mb(k_sym+h1b-1)+1,irrep_h)
c                call sym_irrepname(geom,int_mb(k_sym+p2b-1)+1,irrep_p)
c                spin_h = spinnames(int_mb(k_spin+h1b-1))
c                spin_p = spinnames(int_mb(k_spin+p2b-1))
                ih1=int_mb(k_offset+h1b-1)+h1-1
                ip2=int_mb(k_offset+p2b-1)+p2-1
                ihh1=int_mb(k_mo_index+ih1)
                ipp2=int_mb(k_mo_index+ip2)
                ind=(ihh1-1) + (ipp2-nh-1)*nh
                dbl_mb(k_rdm+ind)=dbl_mb(k_r1+i-1)/denominator
c                if (nodezero) write(LuOut,9010) 
c     1            int_mb(k_offset+h1b-1)+h1,irrep_h,spin_h,
c     2            int_mb(k_offset+p2b-1)+p2,irrep_p,spin_p,
c     2            dbl_mb(k_r1+i-1),
c     3            ihh1,ipp2
              enddo
            enddo
            if (.not.ma_pop_stack(l_r1))
     1        call errquit('tce_print_hp: MA problem',1,MA_ERR)
          endif
          endif
          endif
        enddo
      enddo
      if(restricted) then
        do ih1=1,nocc(ipol)
           ind1=(2*ih1-1)
           ind2=(2*ih1-1-1)
           do ip2=1,nmo(ipol)-nocc(ipol)
              ind=ind1+(2*ip2-1)*nh
              ind3=ind2+(2*ip2-1-1)*nh
              dbl_mb(k_rdm+ind)=dbl_mb(k_rdm+ind3)
           enddo
        enddo
      endif
c      do ih1=1,nh
c         ind1=(ih1-1)
c         do ip2=1,np
c            ind=ind1+(ip2-1)*nh
c            if(nodezero) then
c               if(dabs(dbl_mb(k_rdm+ind)).gt.1.d-10) 
c     1            write(LuOut,'(f30.16,i5,i5)') 
c     2              dbl_mb(k_rdm+ind),ih1,ip2+nh
c            endif
c         enddo
c      enddo
      call ga_sync()
 9000 format(/,1x,a)
 9010 format(1x,i3,a4,'(',a5,') --- ',i3,a4,'(',a5,')',f20.5,i5,i5)
      return
      end
c $Id$
