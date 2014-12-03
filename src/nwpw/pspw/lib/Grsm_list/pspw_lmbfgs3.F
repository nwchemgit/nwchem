*
*     $Id$                       
*

      subroutine pspw_lmbfgs3_init(max_m0,y,g)
      implicit none
      integer max_m0
      complex*16 y(*)
      complex*16 g(*)

#include "pspw_lmbfgs_common.fh"


*     **** external functions ****
      integer  psi_neq
      external psi_neq

      max_m = max_m0
      call Pack_npack(1,npack1)
      neall = psi_neq(1)+psi_neq(2)
      nsize = neall*npack1
      call pspw_Grsm_list_init('lmbfgs',2*max_m,2*nsize)

      m = 1
      !call Grsm_gg_dScale(npack1,neall,(-1.0d0),g,g)
      call pspw_Grsm_list_store(2*m-1,g)
      call pspw_Grsm_list_store(2*m,y)
      !call Grsm_gg_dScale(npack1,neall,(-1.0d0),g,g)

      return
      end



      subroutine pspw_lmbfgs3(tmin,y,g,s)
      implicit none
      real*8 tmin
      complex*16 y(*)
      complex*16 g(*)
      complex*16 s(*)


#include "bafdecls.fh"
#include "pspw_lmbfgs_common.fh"

*     **** local variables ****
      logical value,ondisk
      integer k
      real*8 alpha(20),beta,sum,h0
      integer yy(2),ss(2)

*     **** external functions ****
      logical  control_lmbfgs_ondisk
      external control_lmbfgs_ondisk

      ondisk = control_lmbfgs_ondisk()
      !call Grsm_gg_dScale(npack1,neall,(-1.0d0),g,g)

*     **** allocate yy and ss ****
      if (ondisk) then
        value = BA_push_get(mt_dbl,2*nsize,'yy',yy(2),yy(1))
        value = value.and.
     >          BA_push_get(mt_dbl,2*nsize,'ss',ss(2),ss(1))
        if (.not.value) call errquit('pspw_lmbfgs: push stack',0,0)
      end if

      call Grsm_gg_copy(npack1,neall,g,s)
      if (ondisk) then
        call pspw_Grsm_list_load(2*m-1,dbl_mb(yy(1)))
        call pspw_Grsm_list_load(2*m,  dbl_mb(ss(1)))
      else
        call pspw_Grsm_list_ptr(2*m-1,yy(1))
        call pspw_Grsm_list_ptr(2*m,  ss(1))
      end if
      call psi_1geodesic_Gtransport(tmin,dbl_mb(yy(1)))
      call psi_1geodesic_Gtransport(tmin,dbl_mb(ss(1)))
      call Grsm_ggg_Sub(npack1,neall,g,dbl_mb(yy(1)),dbl_mb(yy(1)))
      call Grsm_ggg_Sub(npack1,neall,y,dbl_mb(ss(1)),dbl_mb(ss(1)))
      if (ondisk) then
        call pspw_Grsm_list_store(2*m-1,dbl_mb(yy(1)))
        call pspw_Grsm_list_store(2*m,  dbl_mb(ss(1)))
      end if

      call Grsm_gg_trace(npack1,neall,
     >                   dbl_mb(yy(1)),
     >                   dbl_mb(ss(1)),
     >                   sum)
      call D1dB_SumAll(sum)
      rho(m) = 1.0d0/sum

      call Grsm_gg_trace(npack1,neall,
     >                   dbl_mb(ss(1)),
     >                   s,
     >                   sum)
      call D1dB_SumAll(sum)
      alpha(m) = rho(m)*sum
      call Grsm_gg_daxpy(npack1,neall,-alpha(m),dbl_mb(yy(1)),s)

      do k=(m-1),1, -1
         if (ondisk) then
           call pspw_Grsm_list_load(2*k-1,dbl_mb(yy(1)))
           call pspw_Grsm_list_load(2*k,  dbl_mb(ss(1)))
         else
           call pspw_Grsm_list_ptr(2*k-1,yy(1))
           call pspw_Grsm_list_ptr(2*k,  ss(1))
         end if
         call psi_1geodesic_Gtransport(tmin,dbl_mb(yy(1)))
         call psi_1geodesic_Gtransport(tmin,dbl_mb(ss(1)))
         if (ondisk) then
           call pspw_Grsm_list_store(2*k-1,dbl_mb(yy(1)))
           call pspw_Grsm_list_store(2*k,  dbl_mb(ss(1)))
         end if
         call Grsm_gg_trace(npack1,neall,
     >                   dbl_mb(ss(1)),
     >                   s,
     >                   sum)
         call D1dB_SumAll(sum)
         alpha(k) = rho(k)*sum
         call Grsm_gg_daxpy(npack1,neall,-alpha(k),dbl_mb(yy(1)),s)

      end do

   
*     **** preconditioner ****
c      call Grsm_gg_dScale(npack1,neall,h0,s,s)


      do k=1,(m-1)
         call Grsm_gg_trace(npack1,neall,
     >                   dbl_mb(yy(1)),
     >                   s,
     >                   sum)
         call D1dB_SumAll(sum)
         beta = rho(k)*sum
         sum = alpha(k) - beta
         call Grsm_gg_daxpy(npack1,neall,sum,dbl_mb(ss(1)),s)
         if (ondisk) then
           call pspw_Grsm_list_load(2*(k+1)-1,dbl_mb(yy(1)))
           call pspw_Grsm_list_load(2*(k+1),  dbl_mb(ss(1)))
         else
           call pspw_Grsm_list_ptr(2*(k+1)-1,yy(1))
           call pspw_Grsm_list_ptr(2*(k+1),  ss(1))
         end if

      end do
      call Grsm_gg_trace(npack1,neall,
     >                dbl_mb(yy(1)),
     >                s,
     >                sum)
      call D1dB_SumAll(sum)
      beta = rho(m)*sum
      sum = alpha(m) - beta
      call Grsm_gg_daxpy(npack1,neall,sum,dbl_mb(ss(1)),s)

      if (m.lt.max_m) then
         m = m+1
      else
         call pspw_Grsm_list_shift()
         call pspw_Grsm_list_shift()
         do k=1,(m-1)
            rho(k) = rho(k+1)
         end do
      end if

      !call Grsm_gg_dScale(npack1,neall,(-1.0d0),s,s)
      call pspw_Grsm_list_store(2*m-1,g)
      call pspw_Grsm_list_store(2*m,y)

      if (ondisk) then
        value = BA_pop_stack(ss(2))
        value = value.and.
     >          BA_pop_stack(yy(2))
        if (.not.value) call errquit('pspw_lmbfgs3:pop stack',2,0)
      end if

      !call Grsm_gg_dScale(npack1,neall,(-1.0d0),g,g)
      return 
      end

