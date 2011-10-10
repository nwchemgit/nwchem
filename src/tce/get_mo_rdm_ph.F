      subroutine get_mo_rdm_ph(d_r1,k_r1_offset,k_rdm,denominator)
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
#include "tce_main.fh"
      integer d_r1
      integer p1b
      integer h2b
      integer p1
      integer h2
      integer k_r1_offset
      integer size
      integer l_r1,k_r1
      integer i
      logical nodezero
      character*4 irrep_p,irrep_h
      character*5 spinnames(2),spin_p,spin_h
      data spinnames/'alpha','beta '/
c->pengdong.fan
      integer k_rdm
      integer ip1,ih2,ipp1,ihh2,ind,np,nh,ind1,ind2,ind3
      double precision denominator
c<-pengdong.fan
c
      nodezero = (ga_nodeid().eq.0)
c      if (nodezero) write(LuOut,9000) 'particle-hole block'
      np=nmo(1)+nmo(ipol)-nocc(1)-nocc(ipol)
      nh=nocc(1)+nocc(ipol)
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
          if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1      +int_mb(k_spin+h2b-1).ne.4)) then
          if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1      .eq.irrep_e) then
            size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
            if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1        call errquit('tce_print_ph: MA problem',0,MA_ERR)
            call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1        int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
            i = 0
            do p1 = 1,int_mb(k_range+p1b-1)
              do h2 = 1,int_mb(k_range+h2b-1)
                i = i + 1
c                call sym_irrepname(geom,int_mb(k_sym+p1b-1)+1,irrep_p)
c                call sym_irrepname(geom,int_mb(k_sym+h2b-1)+1,irrep_h)
c                spin_p = spinnames(int_mb(k_spin+p1b-1))
c                spin_h = spinnames(int_mb(k_spin+h2b-1))
                ip1=int_mb(k_offset+p1b-1)+p1-1
                ih2=int_mb(k_offset+h2b-1)+h2-1
                ipp1=int_mb(k_mo_index+ip1)
                ihh2=int_mb(k_mo_index+ih2)
                ind=(ipp1-nh-1)+(ihh2-1)*np
                dbl_mb(k_rdm+ind)=dbl_mb(k_r1+i-1)/denominator
c                if (nodezero) write(LuOut,9010) 
c     1            int_mb(k_offset+p1b-1)+p1,irrep_p,spin_p,
c     2            int_mb(k_offset+h2b-1)+h2,irrep_h,spin_h,
c     2            dbl_mb(k_r1+i-1),ipp1,ihh2
              enddo
            enddo
            if (.not.ma_pop_stack(l_r1))
     1        call errquit('tce_print_ph: MA problem',1,MA_ERR)
          endif
          endif
          endif
        enddo
      enddo
      if(restricted) then
         do ip1 = 1,nmo(ipol)-nocc(ipol)
            ind1 = (2*ip1-1)
            ind2 = (2*ip1-1-1)
            do ih2 = 1,nocc(ipol)
               ind  = ind1+(2*ih2-1)*np
               ind3 = ind2+(2*ih2-1-1)*np
               dbl_mb(k_rdm+ind) = dbl_mb(k_rdm+ind3)
            enddo
         enddo
      endif
c      do ip1=1,np
c         ind1=(ip1-1)
c         do ih2=1,nh
c            ind=ind1+(ih2-1)*np
c            if(dabs(dbl_mb(k_rdm+ind)).gt.1.d-10) then
c              if(nodezero) write(LuOut,'(f30.16,i5,i5)') 
c     1           dbl_mb(k_rdm+ind),ip1+nh,ih2
c            endif
c         enddo
c      enddo
      call ga_sync()
 9000 format(/,1x,a)
 9010 format(1x,i3,a4,'(',a5,') --- ',i3,a4,'(',a5,')',f20.5,i4,i4)
      return
      end
c $Id$
