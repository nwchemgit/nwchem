*
*     $Id$                       
*

      subroutine paw_lmbfgs_init(max_m0,g)
      implicit none
      integer max_m0
      complex*16 g(*)

#include "paw_lmbfgs_common.fh"


*     **** external functions ****
      integer  paw_psi_ne
      external paw_psi_ne

      max_m = max_m0
      call Pack_npack(1,npack1)
      neall = paw_psi_ne(1)+paw_psi_ne(2)
      nsize = neall*npack1
      call paw_Grsm_list_init('lmbfgs',2*max_m,2*nsize)

      m = 1
c      call Grsm_gg_dScale(npack1,neall,(-1.0d0),g,g)
      call Grsm_gg_dScale1(npack1,neall,(-1.0d0),g)
      call paw_Grsm_list_store(2*m-1,g)
      call paw_Grsm_list_store(2*m,g)
c      call Grsm_gg_dScale(npack1,neall,(-1.0d0),g,g)
      call Grsm_gg_dScale1(npack1,neall,(-1.0d0),g)
      return
      end



      subroutine paw_lmbfgs(tmin,g,s)
      implicit none
      real*8 tmin
      complex*16 g(*)
      complex*16 s(*)


#include "bafdecls.fh"
#include "paw_lmbfgs_common.fh"

*     **** local variables ****
      logical value,ondisk
      integer k,ngras
      real*8 alpha(20),beta,sum,h0
      integer yy(2),ss(2)

*     **** external functions ****
      logical  control_lmbfgs_ondisk
      external control_lmbfgs_ondisk

      ondisk = control_lmbfgs_ondisk()
      ngras = 2*npack1*neall

      !call Grsm_gg_dScale(npack1,neall,(-1.0d0),g,g)
      call dscal(ngras,(-1.0d0),g,1)

*     **** allocate yy and ss ****
      if (ondisk) then
        value = BA_push_get(mt_dbl,2*nsize,'yy',yy(2),yy(1))
        value = value.and.
     >          BA_push_get(mt_dbl,2*nsize,'ss',ss(2),ss(1))
        if (.not.value) call errquit('paw_lmbfgs: push stack',0,0)
      end if
      call dcopy(ngras,g,1,s,1)

      if (ondisk) then
        call paw_Grsm_list_load(2*m-1,dbl_mb(yy(1)))
        call paw_Grsm_list_load(2*m,  dbl_mb(ss(1)))
      else
        call paw_Grsm_list_ptr(2*m-1,yy(1))
        call paw_Grsm_list_ptr(2*m,  ss(1))
      end if
      call paw_psi_1geodesic_Gtransport(tmin,dbl_mb(yy(1)))
      call paw_psi_1geodesic_transport(tmin,dbl_mb(ss(1)))
      call daxpy(ngras,(-1.0d0),g,1,dbl_mb(yy(1)),1)
      call dscal(ngras,(-1.0d0),dbl_mb(yy(1)),1)
      if (ondisk) then
        call paw_Grsm_list_store(2*m-1,dbl_mb(yy(1)))
        call paw_Grsm_list_store(2*m,  dbl_mb(ss(1)))
      end if

      call Grsm_gg_trace(npack1,neall,
     >                   dbl_mb(yy(1)),
     >                   dbl_mb(ss(1)),
     >                   sum)
      !*** exit if dividing by small number ***
      if (dabs(sum).lt.1.0d-15) go to 211

      rho(m) = 1.0d0/sum
      call Grsm_gg_trace(npack1,neall,
     >                   dbl_mb(ss(1)),
     >                   s,
     >                   sum)
      alpha(m) = rho(m)*sum
      call daxpy(ngras,(-alpha(m)),dbl_mb(yy(1)),1,s,1)
      do k=(m-1),1, -1
         if (ondisk) then
           call paw_Grsm_list_load(2*k-1,dbl_mb(yy(1)))
           call paw_Grsm_list_load(2*k,  dbl_mb(ss(1)))
         else
           call paw_Grsm_list_ptr(2*k-1,yy(1))
           call paw_Grsm_list_ptr(2*k,  ss(1))
         end if
         call paw_psi_1geodesic_Gtransport(tmin,dbl_mb(yy(1)))
         call paw_psi_1geodesic_Gtransport(tmin,dbl_mb(ss(1)))
         if (ondisk) then
           call paw_Grsm_list_store(2*k-1,dbl_mb(yy(1)))
           call paw_Grsm_list_store(2*k,  dbl_mb(ss(1)))
         end if

         call Grsm_gg_trace(npack1,neall,
     >                   dbl_mb(ss(1)),
     >                   s,
     >                   sum)
         alpha(k) = rho(k)*sum
         call daxpy(ngras,(-alpha(k)),dbl_mb(yy(1)),1,s,1)
      end do

   
*     **** preconditioner ****
c      call Grsm_gg_dScale(npack1,neall,h0,s,s)


      do k=1,(m-1)

         call Grsm_gg_trace(npack1,neall,
     >                   dbl_mb(yy(1)),
     >                   s,
     >                   sum)
         beta = rho(k)*sum
         sum = alpha(k) - beta
         call daxpy(ngras,sum,dbl_mb(ss(1)),1,s,1)

         if (ondisk) then
           call paw_Grsm_list_load(2*(k+1)-1,dbl_mb(yy(1)))
           call paw_Grsm_list_load(2*(k+1),  dbl_mb(ss(1)))
         else
           call paw_Grsm_list_ptr(2*(k+1)-1,yy(1))
           call paw_Grsm_list_ptr(2*(k+1),  ss(1))
         end if

      end do
      call Grsm_gg_trace(npack1,neall,
     >                dbl_mb(yy(1)),
     >                s,
     >                sum)
      beta = rho(m)*sum
      sum = alpha(m) - beta
      call daxpy(ngras,sum,dbl_mb(ss(1)),1,s,1)
      if (m.lt.max_m) then
         m = m+1
      else
         call paw_Grsm_list_shift()
         call paw_Grsm_list_shift()
         do k=1,(m-1)
            rho(k) = rho(k+1)
         end do
      end if

      call dscal(ngras,(-1.0d0),s,1)
 211  call paw_Grsm_list_store(2*m-1,g)
      call paw_Grsm_list_store(2*m,s)
      call dscal(ngras, (-1.0d0),g,1)

      if (ondisk) then
        value = BA_pop_stack(ss(2))
        value = value.and.
     >          BA_pop_stack(yy(2))
        if (.not.value) call errquit('paw_lmbfgs:pop stack',2,0)
      end if

      return 
      end

