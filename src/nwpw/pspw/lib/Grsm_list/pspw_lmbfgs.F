*
*     $Id: pspw_lmbfgs.F,v 1.5 2001-12-27 20:54:55 bylaska Exp $                       
*

      subroutine pspw_lmbfgs_init(max_m0,g)
      implicit none
      integer max_m0
      complex*16 g(*)

#include "pspw_lmbfgs_common.fh"

*     **** local variables ****

*     **** external functions ****
      integer  psi_ne
      external psi_ne
      real*8 sum

      max_m = max_m0
      call Pack_npack(1,npack1)
      neall = psi_ne(1)+psi_ne(2)
      nsize = neall*npack1
      call pspw_Grsm_list_init('lmbfgs',max_m,2*nsize)

      m = 1
      call Grsm_gg_dScale(npack1,neall,(-1.0d0),g,g)
      call pspw_Grsm_list_store(m,g,g)
      call Grsm_gg_dScale(npack1,neall,(-1.0d0),g,g)
      return
      end





      subroutine pspw_lmbfgs(tmin,g,s)
      implicit none
      real*8 tmin
      complex*16 g(*)
      complex*16 s(*)


#include "mafdecls.fh"
#include "pspw_lmbfgs_common.fh"

*     **** local variables ****
      logical value
      integer k
      real*8 alpha(20),beta,sum,h0
      integer yy(2),ss(2)

      call Grsm_gg_dScale(npack1,neall,(-1.0d0),g,g)

*     **** allocate yy and ss ****
      value = MA_push_get(mt_dcpl,nsize,'yy',yy(2),yy(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nsize,'ss',ss(2),ss(1))


      call Grsm_gg_copy(npack1,neall,g,s)
      call pspw_Grsm_list_load(m,dcpl_mb(yy(1)),dcpl_mb(ss(1)))
      call psi_1geodesic_Gtransport(tmin,dcpl_mb(yy(1)))
      call psi_1geodesic_transport(tmin,dcpl_mb(ss(1)))
      call Grsm_ggg_Sub(npack1,neall,g,dcpl_mb(yy(1)),dcpl_mb(yy(1)))
      call pspw_Grsm_list_store(m,dcpl_mb(yy(1)),dcpl_mb(ss(1)))

      call Grsm_gg_trace(npack1,neall,
     >                   dcpl_mb(yy(1)),
     >                   dcpl_mb(ss(1)),
     >                   sum)
      rho(m) = 1.0d0/sum
c     write(*,*) "M,RHO:",m,rho(m)
c     call Grsm_gg_trace(npack1,neall,
c    >                   dcpl_mb(yy(1)),
c    >                   dcpl_mb(yy(1)),
c    >                   h0)
c     write(*,*) "h0:",h0
c     h0 = sum/h0

      call Grsm_gg_trace(npack1,neall,
     >                   dcpl_mb(ss(1)),
     >                   s,
     >                   sum)
      alpha(m) = rho(m)*sum
      call Grsm_gg_daxpy(npack1,neall,-alpha(m),dcpl_mb(yy(1)),s)

      do k=(m-1),1, -1
         call pspw_Grsm_list_load(k,dcpl_mb(yy(1)),dcpl_mb(ss(1)))
         call psi_1geodesic_Gtransport(tmin,dcpl_mb(yy(1)))
         call psi_1geodesic_Gtransport(tmin,dcpl_mb(ss(1)))
         call pspw_Grsm_list_store(k,dcpl_mb(yy(1)),dcpl_mb(ss(1)))
         call Grsm_gg_trace(npack1,neall,
     >                   dcpl_mb(ss(1)),
     >                   s,
     >                   sum)
         alpha(k) = rho(k)*sum
         call Grsm_gg_daxpy(npack1,neall,-alpha(k),dcpl_mb(yy(1)),s)

*     **** debug ****
c      call Grsm_gg_trace(npack1,neall,
c     >                   dcpl_mb(yy(1)),
c     >                   dcpl_mb(ss(1)),
c     >                   sum)
c      write(*,*) "    debug k,RHO:",k,1.0d0/sum,rho(k)

      end do

*     **** preconditioner ****
c      call Grsm_gg_dScale(npack1,neall,h0,s,s)


      do k=1,(m-1)
         call Grsm_gg_trace(npack1,neall,
     >                   dcpl_mb(yy(1)),
     >                   s,
     >                   sum)
         beta = rho(k)*sum
         sum = beta - alpha(m)
         sum = alpha(k) - beta
         call Grsm_gg_daxpy(npack1,neall,sum,dcpl_mb(ss(1)),s)
         call pspw_Grsm_list_load(k+1,dcpl_mb(yy(1)),dcpl_mb(ss(1)))
      end do
      call Grsm_gg_trace(npack1,neall,
     >                dcpl_mb(yy(1)),
     >                s,
     >                sum)
      beta = rho(m)*sum
      sum = beta - alpha(m)
      sum = alpha(m) - beta
      call Grsm_gg_daxpy(npack1,neall,sum,dcpl_mb(ss(1)),s)

      if (m.lt.max_m) then
         m = m+1
      else
         call pspw_Grsm_list_shift()
         do k=1,(m-1)
            rho(k) = rho(k+1)
         end do
      end if

      call Grsm_gg_dScale(npack1,neall,(-1.0d0),s,s)
      call pspw_Grsm_list_store(m,g,s)

      value = MA_pop_stack(ss(2))
      value = value.and.
     >        MA_pop_stack(yy(2))

      call Grsm_gg_dScale(npack1,neall,(-1.0d0),g,g)
      return 
      end

