*
*     $Id: band_lmbfgs.F,v 1.2 2006-11-29 02:33:29 d3p708 Exp $                       
*

      subroutine band_lmbfgs_init(max_m0,g)
      implicit none
      integer max_m0
      complex*16 g(*)

#include "band_lmbfgs_common.fh"


*     **** external functions ****
      integer  cpsi_ispin,cpsi_ne,brillioun_nbrillioun
      external cpsi_ispin,cpsi_ne,brillioun_nbrillioun

      max_m = max_m0
      call Cram_max_npack(npack1)
      ispin = cpsi_ispin()
      ne(1) = cpsi_ne(1)
      ne(2) = cpsi_ne(2)
      neall = ne(1)+ne(2)
      nbrill = brillioun_nbrillioun()
      nsize = neall*npack1*nbrill
      call band_Grsm_list_init('lmbfgs',2*max_m,2*nsize)
      write(0,*)"band_Grsm_init ",nsize,max_m,neall,ne(1),nbrill
      m = 1
      call Stfl_gg_dScale(ispin,ne,npack1,neall,nbrill,(-1.0d0),g,g)
      call band_Grsm_list_store(2*m-1,g)
      call band_Grsm_list_store(2*m,g)
      call Stfl_gg_dScale(ispin,ne,npack1,neall,nbrill,(-1.0d0),g,g)
      return
      end



      subroutine band_lmbfgs(tmin,g,s)
      implicit none
      real*8 tmin
      complex*16 g(*)
      complex*16 s(*)

#include "mafdecls.fh"
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

      write(0,*)"In band_lmbfgs " ,nsize
      write(0,*)ngras,ondisk
      call dscal(ngras,(-1.0d0),g,1)

*     **** allocate yy and ss ****
      if (ondisk) then
        value = MA_push_get(mt_dbl,2*nsize,'yy',yy(2),yy(1))
        value = value.and.
     >          MA_push_get(mt_dbl,2*nsize,'ss',ss(2),ss(1))
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
      call cpsi_1geodesic_Gtransport(tmin,dbl_mb(yy(1)))
      call cpsi_1geodesic_transport(tmin,dbl_mb(ss(1)))
      
      call daxpy(ngras,(-1.0d0),g,1,dbl_mb(yy(1)),1)
      call dscal(ngras,(-1.0d0),dbl_mb(yy(1)),1)
      if (ondisk) then
        call band_Grsm_list_store(2*m-1,dbl_mb(yy(1)))
        call band_Grsm_list_store(2*m,  dbl_mb(ss(1)))
      end if

     
c      casa1=yy(2)-yy(1)
c      casa2=ss(2)-ss(1)
      write(0,*)"Doing first trace in band_lmbfgs",casa2,casa1
      call Stfl_gg_traceall(ispin,ne,npack1,neall,nbrill,
     >                   dbl_mb(yy(1)),
     >                   dbl_mb(ss(1)),
     >                   sum)
      write(0,*)"did my trace"
      rho(m) = 1.0d0/sum
      write(0,*)"Doing 2nd trace in band_lmbfgs"
      call Stfl_gg_traceall(ispin,ne,npack1,neall,nbrill,
     >                   dbl_mb(ss(1)),
     >                   s,
     >                   sum)
      alpha(m) = rho(m)*sum
      call daxpy(ngras,(-alpha(m)),dbl_mb(yy(1)),1,s,1)
      write(0,*)"Doing loop in band_lmbfgs"
      do k=(m-1),1, -1
        write(0,*)"Doing load band_lmbfgs",k
         if (ondisk) then
           call band_Grsm_list_load(2*k-1,dbl_mb(yy(1)))
           call band_Grsm_list_load(2*k,  dbl_mb(ss(1)))
         else
           call band_Grsm_list_ptr(2*k-1,yy(1))
           call band_Grsm_list_ptr(2*k,  ss(1))
         end if
         write(0,*)"Doing Gtransport in band_lmbfgs",k
         call cpsi_1geodesic_Gtransport(tmin,dbl_mb(yy(1)))
         write(0,*)"Doing Gtransport in band_lmbfgs",k
         call cpsi_1geodesic_Gtransport(tmin,dbl_mb(ss(1)))
         if (ondisk) then
           call band_Grsm_list_store(2*k-1,dbl_mb(yy(1)))
           call band_Grsm_list_store(2*k,  dbl_mb(ss(1)))
         end if

         write(0,*)"Doing traceall in band_lmbfgs",k
         call Stfl_gg_traceall(ispin,ne,npack1,neall,nbrill,
     >                   dbl_mb(ss(1)),
     >                   s,
     >                   sum)
         alpha(k) = rho(k)*sum
         write(0,*)"Doing  daxpy in band_lmbfgs",k
         call daxpy(ngras,(-alpha(k)),dbl_mb(yy(1)),1,s,1)
      end do

   
*     **** preconditioner ****


      do k=1,(m-1)

         write(0,*)"Doing traceall in band_lmbfgs",k
         call Stfl_gg_traceall(ispin,ne,npack1,neall,nbrill,
     >                   dbl_mb(yy(1)),
     >                   s,
     >                   sum)
         beta = rho(k)*sum
         sum = alpha(k) - beta
      
         write(0,*)"Doing  daxpy in band_lmbfgs",k
         call daxpy(ngras,sum,dbl_mb(ss(1)),1,s,1)
         if (ondisk) then
           call band_Grsm_list_load(2*(k+1)-1,dbl_mb(yy(1)))
           call band_Grsm_list_load(2*(k+1),  dbl_mb(ss(1)))
         else
           call band_Grsm_list_ptr(2*(k+1)-1,yy(1))
           call band_Grsm_list_ptr(2*(k+1),  ss(1))
         end if

      end do
         write(0,*)"Doing  daxpy in band_lmbfgs"
      call Stfl_gg_traceall(ispin,ne,npack1,neall,nbrill,
     >                dbl_mb(yy(1)),
     >                s,
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

c      if (ondisk) then
        value = MA_pop_stack(ss(2))
        value = value.and.
     >          MA_pop_stack(yy(2))
        if (.not.value) call errquit('band_lmbfgs:pop stack',2,0)
c      end if

      write(0,*)'out of band_lmbfgs'
      return 
      end

