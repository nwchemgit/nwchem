      subroutine get_mo_rdm_pp(d_r1,k_r1_offset,k_rdm,denominator)
c
c Prints the particle-particle block of 1-e reduced density matrix
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
      integer p2b
      integer p1
      integer p2
      integer k_r1_offset
      integer size
      integer l_r1,k_r1
      integer i
      logical nodezero
      double precision r1_thresh
      character*4 irrep_p1,irrep_p2
      character*5 spinnames(2),spin_p1,spin_p2
      data spinnames/'alpha','beta '/
c->pengdong.fan
      integer k_rdm
      integer ip1,ip2,ipp1,ipp2,ind,np,nh,ind1,ind2,ind3
      double precision denominator
c<-pengdong.fan
c
      nodezero = (ga_nodeid().eq.0)
c      if (nodezero) write(LuOut,9000) 'particle-particle block'
      np=nmo(1)+nmo(ipol)-nocc(1)-nocc(ipol)
      nh=nocc(1)+nocc(ipol)
      do p1b = noab+1,noab+nvab
        do p2b = noab+1,noab+nvab
          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+p2b-1)) then 
          if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1      +int_mb(k_spin+p2b-1).ne.4)) then
          if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+p2b-1))
     1      .eq.irrep_e) then
            size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
            if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1        call errquit('tce_print_pp: MA problem',0,MA_ERR)
            call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1        int_mb(k_r1_offset),((p1b-noab-1)*nvab+p2b-noab-1))
            i = 0
            do p1 = 1,int_mb(k_range+p1b-1)
              do p2 = 1,int_mb(k_range+p2b-1)
                i = i + 1
c                call sym_irrepname(geom,int_mb(k_sym+p1b-1)+1,
c     1               irrep_p1)
c                call sym_irrepname(geom,int_mb(k_sym+p2b-1)+1,
c     1               irrep_p2)
c                spin_p1 = spinnames(int_mb(k_spin+p1b-1))
c                spin_p2 = spinnames(int_mb(k_spin+p2b-1))
                ip1=int_mb(k_offset+p1b-1)+p1-1
                ip2=int_mb(k_offset+p2b-1)+p2-1
                ipp1=int_mb(k_mo_index+ip1)
                ipp2=int_mb(k_mo_index+ip2)
                ind=(ipp1-nh-1)+(ipp2-nh-1)*np
                dbl_mb(k_rdm+ind)=dbl_mb(k_r1+i-1)/denominator
c                if (nodezero) write(LuOut,9010) 
c     1            int_mb(k_offset+p1b-1)+p1,irrep_p1,spin_p1,
c     2            int_mb(k_offset+p2b-1)+p2,irrep_p2,spin_p2,
c     2            dbl_mb(k_r1+i-1),ipp1,ipp2
              enddo
            enddo
            if (.not.ma_pop_stack(l_r1))
     1        call errquit('tce_print_pp: MA problem',1,MA_ERR)
          endif
          endif
          endif
        enddo
      enddo
      if(restricted) then 
         do ip1=1,nmo(ipol)-nocc(ipol)
            ind1 = (2*ip1-1)
            ind2 = (2*ip1-2)
            do ip2=1,nmo(ipol)-nocc(ipol)
               ind  = ind1 + (2*ip2-1)*np
               ind3 = ind2 + (2*ip2-1-1)*np
               dbl_mb(k_rdm+ind)=dbl_mb(k_rdm+ind3)
            enddo
         enddo
      endif
c      do ip1=1,np
c         ind1=(ip1-1)
c         do ip2=1,np
c            ind=ind1+(ip2-1)*np
c            if(dabs(dbl_mb(k_rdm+ind)).gt.1.d-10) then
c               if(nodezero) write(LuOut,'(f30.16,i5,i5)')
c     1            dbl_mb(k_rdm+ind),ip1+nh,ip2+nh
c            endif
c         enddo
c      enddo  
      call ga_sync()
 9000 format(/,1x,a)
 9010 format(1x,i3,a4,'(',a5,') --- ',i3,a4,'(',a5,')',f20.5,i4,i4)
      return
      end
c $Id$
