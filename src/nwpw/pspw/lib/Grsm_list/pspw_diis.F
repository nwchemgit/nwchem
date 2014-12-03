*
*     $Id$                       
*

      subroutine pspw_diis_init(max_m0,g)
      implicit none
      integer max_m0
      complex*16 g(*)

#include "pspw_diis_common.fh"

*     **** local variables ****
      integer i

*     **** external functions ****
      integer  psi_neq
      external psi_neq

      max_m = max_m0
      call Pack_npack(1,npack1)
      neall = psi_neq(1)+psi_neq(2)
      nsize = neall*npack1
      call pspw_Grsm_list_init('pspwdiis',2*max_m,2*nsize)

      call dcopy(max_list*max_list,0.0d0,0,B,1)
      call dcopy(max_list,0.0d0,0,d,1)
      do i=2,max_list
        B(i,1) = 1.0d0
        B(1,i) = 1.0d0
      end do
      d(1) = 1.0d0

      m = 1
      call pspw_Grsm_list_store(m,g)
      return
      end



      subroutine pspw_diis(tmin,g,s)
      implicit none
      real*8 tmin
      complex*16 g(*)
      complex*16 s(*)


#include "bafdecls.fh"
#include "pspw_diis_common.fh"

*     **** local variables ****
      logical value
      integer i,j
      real*8 sum, A(max_list,max_list)
      integer gg(2),ipiv(max_list),info


*     **** allocate gg and ss ****
      value = BA_push_get(mt_dcpl,nsize,'gg',gg(2),gg(1))


      call pspw_Grsm_list_store(m,g)

*     **** update B matrix ****
      call Grsm_gg_trace(npack1,neall,g,g,sum)
      call D1dB_SumAll(sum)
      B(m+1,m+1) = sum
      do i=1,m-1
         call pspw_Grsm_list_load(i,dcpl_mb(gg(1)))
         call psi_1geodesic_Gtransport(tmin,dcpl_mb(gg(1)))
         call pspw_Grsm_list_store(i,dcpl_mb(gg(1)))
         call Grsm_gg_trace(npack1,neall,dcpl_mb(gg(1)),g,sum)
         call D1dB_SumAll(sum)
         B(i+1,m+1) = sum
         B(m+1,i+1) = sum
      end do

*     **** solve Bc = d ****
      call dcopy(max_list,d,1,c,1)
      call dcopy(max_list*max_list,B,1,A,1)
      call dgesv(m+1,1,A,max_list,ipiv,c,max_list,info)

*     **** s = sum(i) c(i) * gg(i) ****
      call Grsm_gg_dScale(npack1,neall,c(m+1),g,s)
      do i=1,m-1
         call pspw_Grsm_list_load(i,dcpl_mb(gg(1)))
         call Grsm_gg_daxpy(npack1,neall,c(i+1),dcpl_mb(gg(1)),s)
      end do

*     **** shift B ****
      if (m.lt.max_m) then
         m = m+1
      else
         call pspw_Grsm_list_shift()
         do j=2,m
         do i=2,m
           B(i,j) = B(i+1,j+1)
         end do
         end do
      end if

      value = BA_pop_stack(gg(2))
      return 
      end

