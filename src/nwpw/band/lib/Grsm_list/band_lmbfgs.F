*
*     $Id$                       
*

      subroutine band_lmbfgs_init(max_m0,g)
      implicit none
      integer max_m0
      complex*16 g(*)

#include "band_lmbfgs_common.fh"

      integer ngras

*     **** external functions ****
      integer  cpsi_ispin,cpsi_ne,brillioun_nbrillq
      external cpsi_ispin,cpsi_ne,brillioun_nbrillq

      max_m = max_m0
      call Cram_max_npack(npack1)
      ispin = cpsi_ispin()
      ne(1) = cpsi_ne(1)
      ne(2) = cpsi_ne(2)
      neall = ne(1)+ne(2)
      nbrill = brillioun_nbrillq()
      nsize = neall*npack1*nbrill
      call band_Grsm_list_init('lmbfgs',2*max_m,2*nsize)
      m = 1
c      call BGrsm_ff_dScale(-1.0d0,g,g)
      ngras = 2*nsize
      call dscal(ngras,(-1.0d0),g,1)
      call band_Grsm_list_store(2*m-1,g)
      call band_Grsm_list_store(2*m,g)
c      call BGrsm_ff_dScale(-1.0d0,g,g)
      call dscal(ngras,(-1.0d0),g,1)
      return
      end



      subroutine band_lmbfgs(tmin,g,s)
      implicit none
      real*8 tmin
      complex*16 g(*)
      complex*16 s(*)

#include "bafdecls.fh"
#include "band_lmbfgs_common.fh"

*     **** local variables ****
      logical value,ondisk,casa1,casa2
      integer k,ngras
      real*8 alpha(20),beta,sum,h0
      integer yy(2),ss(2)

*     **** external functions ****
      logical  control_lmbfgs_ondisk
      external control_lmbfgs_ondisk

      ondisk = control_lmbfgs_ondisk()
      ngras = 2*npack1*neall*nbrill

      call dscal(ngras,(-1.0d0),g,1)

*     **** allocate yy and ss ****
      if (ondisk) then
        value = BA_push_get(mt_dbl,2*nsize,'yy',yy(2),yy(1))
        value = value.and.
     >          BA_push_get(mt_dbl,2*nsize,'ss',ss(2),ss(1))
        if (.not.value) call errquit('band_lmbfgs: push stack',0,0)
      end if

      call dcopy(ngras,g,1,s,1)

      if (ondisk) then
        call band_Grsm_list_load(2*m-1,dbl_mb(yy(1)))
        call band_Grsm_list_load(2*m,  dbl_mb(ss(1)))
      else
        call band_Grsm_list_ptr(2*m-1,yy(1))
        call band_Grsm_list_ptr(2*m,  ss(1))
      end if

      call cpsi_1geodesic_Gtransport_junk(tmin,dbl_mb(yy(1)))
      call cpsi_1geodesic_transport_junk(tmin,dbl_mb(ss(1)))
      
      call daxpy(ngras,(-1.0d0),g,1,dbl_mb(yy(1)),1)
      call dscal(ngras,(-1.0d0),dbl_mb(yy(1)),1)
      if (ondisk) then
        call band_Grsm_list_store(2*m-1,dbl_mb(yy(1)))
        call band_Grsm_list_store(2*m,  dbl_mb(ss(1)))
      end if

     
      call Pneb_ff_traceall(0,0,
     >                   dbl_mb(yy(1)),
     >                   dbl_mb(ss(1)),npack1,
     >                   sum)
      rho(m) = 1.0d0/sum
      call Pneb_ff_traceall(0,0,
     >                   dbl_mb(ss(1)),
     >                   s,npack1,
     >                   sum)
      alpha(m) = rho(m)*sum
      call daxpy(ngras,(-alpha(m)),dbl_mb(yy(1)),1,s,1)
      do k=(m-1),1, -1
         if (ondisk) then
           call band_Grsm_list_load(2*k-1,dbl_mb(yy(1)))
           call band_Grsm_list_load(2*k,  dbl_mb(ss(1)))
         else
           call band_Grsm_list_ptr(2*k-1,yy(1))
           call band_Grsm_list_ptr(2*k,  ss(1))
         end if
         call cpsi_1geodesic_Gtransport_junk(tmin,dbl_mb(yy(1)))
         call cpsi_1geodesic_Gtransport_junk(tmin,dbl_mb(ss(1)))
         if (ondisk) then
           call band_Grsm_list_store(2*k-1,dbl_mb(yy(1)))
           call band_Grsm_list_store(2*k,  dbl_mb(ss(1)))
         end if

         call Pneb_ff_traceall(0,0,
     >                   dbl_mb(ss(1)),
     >                   s,npack1,
     >                   sum)
         alpha(k) = rho(k)*sum
         call daxpy(ngras,(-alpha(k)),dbl_mb(yy(1)),1,s,1)
      end do

   
*     **** preconditioner ****


      do k=1,(m-1)

         call Pneb_ff_traceall(0,0,
     >                   dbl_mb(yy(1)),
     >                   s,npack1,
     >                   sum)
         beta = rho(k)*sum
         sum = alpha(k) - beta
      
         call daxpy(ngras,sum,dbl_mb(ss(1)),1,s,1)
         if (ondisk) then
           call band_Grsm_list_load(2*(k+1)-1,dbl_mb(yy(1)))
           call band_Grsm_list_load(2*(k+1),  dbl_mb(ss(1)))
         else
           call band_Grsm_list_ptr(2*(k+1)-1,yy(1))
           call band_Grsm_list_ptr(2*(k+1),  ss(1))
         end if

      end do
      call Pneb_ff_traceall(0,0,
     >                dbl_mb(yy(1)),
     >                s,npack1,
     >                sum)

      beta = rho(m)*sum
      sum = alpha(m) - beta
      call daxpy(ngras,sum,dbl_mb(ss(1)),1,s,1)
      if (m.lt.max_m) then
         m = m+1
      else
         call band_Grsm_list_shift()
         call band_Grsm_list_shift()
         do k=1,(m-1)
            rho(k) = rho(k+1)
         end do
      end if

      call dscal(ngras,(-1.0d0),s,1)
      call band_Grsm_list_store(2*m-1,g)
      call band_Grsm_list_store(2*m,s)
      call dscal(ngras, (-1.0d0),g,1)

      if (ondisk) then
        value = BA_pop_stack(ss(2))
        value = value.and.
     >          BA_pop_stack(yy(2))
        if (.not.value) call errquit('band_lmbfgs:pop stack',2,0)
      end if

      return 
      end

