*
*     $Id: nwpw_scf_mixing.F,v 1.7 2007-09-14 21:15:35 bylaska Exp $                       
*

*     *****************************************
*     *                                       *
*     *        nwpw_scf_mixing_init           *
*     *                                       *
*     *****************************************

      subroutine nwpw_scf_mixing_init(algorithm0,alpha0,
     >                                max_m0,ispin0,nsize0,rho_in)
      implicit none
      integer algorithm0
      real*8  alpha0
      integer max_m0,ispin0,nsize0
      real*8 rho_in(*)

#include "mafdecls.fh"
#include "nwpw_scf_mixing.fh"

      algorithm = algorithm0
      alpha = alpha0
      max_m = max_m0
      nsize = nsize0*ispin0
      ispin = ispin0

      if (algorithm.eq.0) then
        call nwpw_list_start(1,2,nsize)
        call nwpw_scf_mixing_reset(rho_in)
      end if
      if (algorithm.eq.1) then
        call nwpw_list_start(1,8,nsize)
        call nwpw_scf_mixing_reset(rho_in)
      end if
      if (algorithm.eq.2) then
        call nwpw_list_start(1,(5+(2*max_m)),nsize)
        call nwpw_scf_mixing_reset(rho_in)
        w0 = 0.01d0
      end if
      if (algorithm.eq.3) then
        call nwpw_list_start(1,3,nsize)
        call nwpw_scf_mixing_reset(rho_in)
      end if

      return
      end


*     *****************************************
*     *                                       *
*     *        nwpw_scf_mixing_end            *
*     *                                       *
*     *****************************************

      subroutine nwpw_scf_mixing_end()
      implicit none

#include "mafdecls.fh"
#include "nwpw_scf_mixing.fh"

      call nwpw_list_end(1)
      return
      end


*     *****************************************
*     *                                       *
*     *        nwpw_scf_mixing_reset          *
*     *                                       *
*     *****************************************

      subroutine nwpw_scf_mixing_reset(rho_in)
      implicit none
      real*8 rho_in(*)

#include "mafdecls.fh"
#include "nwpw_scf_mixing.fh"

      m = 1
      call nwpw_list_store(1,1,rho_in)
      return
      end


*     *****************************************
*     *                                       *
*     *              nwpw_scf_mixing          *
*     *                                       *
*     *****************************************

      subroutine nwpw_scf_mixing(vout,vnew,deltae,scf_error)
      implicit none
      real*8 vout(*)
      real*8 vnew(*)
      real*8 deltae
      real*8 scf_error


#include "mafdecls.fh"
#include "nwpw_scf_mixing.fh"

*     **** local variables ****
      logical value
      integer i,j,info,k,ipiv
      integer rr_ptr,ff_ptr,tt_ptr
      integer V0,V1,Vout0,Vout1,F0,F1,Vbar0,Vbar1
      integer dV,U,dF,dFi
      real*8 sum0,sum1,beta,p00,p01,p11,alpha1
      real*8 BB(40,40),BBB(40,40)

*     **** external functions ****
      real*8   ddot
      external ddot

      if ((deltae.ge.0.0d0)) then
         alpha1 = alpha 
      else
         alpha1 = 1.0d0
      end if

*     **** simple mixing ****
      if (algorithm.eq.0) then
         call nwpw_list_ptr(1,1,rr_ptr)
         call nwpw_list_ptr(1,2,ff_ptr)
         call dcopy(nsize,dbl_mb(rr_ptr),1,dbl_mb(ff_ptr),1)  ! f = vm
         call daxpy(nsize,(-1.0d0),vout, 1,dbl_mb(ff_ptr),1)  ! f = f-vout = vm-vout
         call dscal(nsize,(-1.0d0),dbl_mb(ff_ptr),1)          ! f = -f = vout-vm

         scf_error = ddot(nsize,dbl_mb(ff_ptr),1,dbl_mb(ff_ptr),1)
         call D3dB_SumAll(scf_error)
         scf_error = dsqrt(scf_error)                         ! scf_error = sqrt(<f|f>)

         call dcopy(nsize,dbl_mb(rr_ptr),1,vnew,1)            !vnew = vm
         call daxpy(nsize,alpha1,dbl_mb(ff_ptr),1,vnew,1)     !vnew = vnew+alpha*f
         call dcopy(nsize,vnew,1,dbl_mb(rr_ptr),1)            !vm=vnew
      end if

*     **** Broyden potential mixing ****
cccccccccc Actually this is Broyden mixing ....
      if (algorithm.eq.1) then
         call nwpw_list_ptr(1,1,V0)
         call nwpw_list_ptr(1,2,V1)
         call nwpw_list_ptr(1,3,Vout0)
         call nwpw_list_ptr(1,4,Vout1)
         call nwpw_list_ptr(1,5,F0)
         call nwpw_list_ptr(1,6,F1)
         call nwpw_list_ptr(1,7,Vbar0)
         call nwpw_list_ptr(1,8,Vbar1)

         if (m.eq.1) then
           !**** vout0 = vout ****
           call dcopy(nsize,vout,1,dbl_mb(Vout0),1)

           !**** F0 = vout - V0 ****
           call dcopy(nsize,dbl_mb(V0),1,dbl_mb(F0),1)
           call dscal(nsize,(-1.0d0),dbl_mb(F0),1)
           call daxpy(nsize,(1.0d0),vout, 1,dbl_mb(F0),1)

           !**** V1 = V0 + alpha*F0 ****
           call dcopy(nsize,dbl_mb(V0),1,dbl_mb(V1),1)
           call daxpy(nsize,alpha1,dbl_mb(F0),1,dbl_mb(V1),1)   

           !**** vnew = V1 ****
           call dcopy(nsize,dbl_mb(V1),1,vnew,1)

         else !(m.gt.1) 
           !**** vout1 = vout ****
           call dcopy(nsize,vout,1,dbl_mb(Vout1),1)

           !**** F1 = vout - V1 ****
           call dcopy(nsize,dbl_mb(V1),1,dbl_mb(F1),1)
           call dscal(nsize,(-1.0d0),dbl_mb(F1),1)
           call daxpy(nsize,(1.0d0),vout, 1,dbl_mb(F1),1)

           !**** scf_error = sqrt(<F1|F1>) ****
           scf_error = ddot(nsize,dbl_mb(F1),1,dbl_mb(F1),1)
           call D3dB_SumAll(scf_error)
           scf_error = dsqrt(scf_error)          

       
           !**** Beta = <F1|F1-F0>/<F1-F0/F1-F0> ****
           call dcopy(nsize,dbl_mb(F1),1,dbl_mb(Vbar1),1)
           call daxpy(nsize,(-1.0d0),dbl_mb(F0),1,dbl_mb(Vbar1),1)
           sum0 = ddot(nsize,dbl_mb(F1),1,dbl_mb(Vbar1),1)
           sum1 = ddot(nsize,dbl_mb(Vbar1),1,dbl_mb(Vbar1),1)
           call D3dB_SumAll(sum0)
           call D3dB_SumAll(sum1)
           beta = sum0/sum1
           
 

           !**** Vbar1 = (1-Beta)*Vout1 + Beta*Vout0 ****
           call dcopy(nsize,dbl_mb(Vout0),1,dbl_mb(Vbar1),1)
           call dscal(nsize,(beta),dbl_mb(Vbar1),1)
           call daxpy(nsize,(1.0d0-beta),dbl_mb(Vout1),1,
     >                                   dbl_mb(Vbar1),1)

           !**** Vbar0 = (1-Beta)*V1 + Beta*V0       ****
           call dcopy(nsize,dbl_mb(V0),1,dbl_mb(Vbar0),1)
           call dscal(nsize,(beta),dbl_mb(Vbar0),1)
           call daxpy(nsize,(1.0d0-beta),dbl_mb(V1),1,dbl_mb(Vbar0),1)

           !**** F0 = F1, Vout0 = Vout1, V0 = V1 ****
           call dcopy(nsize,dbl_mb(F1),1,dbl_mb(F0),1)
           call dcopy(nsize,dbl_mb(Vout1),1,dbl_mb(Vout0),1)
           call dcopy(nsize,dbl_mb(V1),1,dbl_mb(V0),1)

           !**** V1 = (1-alpha)*Vbar0 + alpha*Vbar1 ****
           call dcopy(nsize,dbl_mb(Vbar1),1,dbl_mb(V1),1)
           call dscal(nsize,(alpha1),dbl_mb(V1),1)
           call daxpy(nsize,(1.0d0-alpha1),dbl_mb(Vbar0),1,dbl_mb(V1),1)
          

           !**** vnew = V1 ****
           call dcopy(nsize,dbl_mb(V1),1,vnew,1)

         end if

         m = m + 1

      end if

*     **** Johnson mixing ****
      if (algorithm.eq.2) then
         call nwpw_list_ptr(1,1,V0)
         call nwpw_list_ptr(1,2,V1)
         call nwpw_list_ptr(1,3,F0)
         call nwpw_list_ptr(1,4,F1)
         call nwpw_list_ptr(1,5,dV)

         if (m.eq.1) then

           !**** F0 = vout - V0 ****
           call dcopy(nsize,dbl_mb(V0),1,dbl_mb(F0),1)
           call dscal(nsize,(-1.0d0),dbl_mb(F0),1)
           call daxpy(nsize,(1.0d0),vout, 1,dbl_mb(F0),1)

           !**** scf_error = sqrt(<F1|F1>) ****
           scf_error = ddot(nsize,dbl_mb(F0),1,dbl_mb(F0),1)
           call D3dB_SumAll(scf_error)
           scf_error = dsqrt(scf_error)          

           !**** V1 = V0 + alpha*F0 ****
           call dcopy(nsize,dbl_mb(V0),1,dbl_mb(V1),1)

           call daxpy(nsize,alpha1,dbl_mb(F0),1,dbl_mb(V1),1)   

           !**** vnew = V1 ****
           call dcopy(nsize,dbl_mb(V1),1,vnew,1)

           call dcopy(max_list*max_list,0.0d0,0,A,1)

           m = m + 1

         else !(m.gt.1) 
           !**** F1 = vout - V1 ****
           call dcopy(nsize,dbl_mb(V1),1,dbl_mb(F1),1)
           call dscal(nsize,(-1.0d0),dbl_mb(F1),1)
           call daxpy(nsize,(1.0d0),vout,1,dbl_mb(F1),1)

           !**** scf_error = sqrt(<F1|F1>) ****
           scf_error = ddot(nsize,dbl_mb(F1),1,dbl_mb(F1),1)
           call D3dB_SumAll(scf_error)
           scf_error = dsqrt(scf_error)          

           !*** dF = dF(m-1), U = U(m-1) ***
           call nwpw_list_ptr(1,(5+m-1),dF)
           call nwpw_list_ptr(1,(5+max_m  +m-1),U)

           !**** dF = (F1-F0) ****
           call dcopy(nsize,dbl_mb(F1),1,dbl_mb(dF),1)
           call daxpy(nsize,(-1.0d0),dbl_mb(F0),1,dbl_mb(dF),1)

           !**** dV  = (V1-V0) ****
           call dcopy(nsize,dbl_mb(V1),1,dbl_mb(dV),1)
           call daxpy(nsize,(-1.0d0),dbl_mb(V0),1,dbl_mb(dV),1)
           
           !**** U  = alpha*dF + dV ****
           call dcopy(nsize,dbl_mb(dV),1,dbl_mb(U),1)

           call daxpy(nsize,(alpha1),dbl_mb(dF),1,dbl_mb(U),1)

           !**** Define A,c and B ****
           do i=1,m-1
              call nwpw_list_ptr(1,(5+i),dFi)
              sum0 = ddot(nsize,dbl_mb(dFi),1,dbl_mb(dF),1)
              sum1 = ddot(nsize,dbl_mb(dFi),1,dbl_mb(F1),1)
              call D3dB_SumAll(sum0)
              call D3dB_SumAll(sum1)

              A(i,m-1) = sum0
              A(m-1,i) = A(i,m-1)

              c(i) = sum1
           end do
           
           call dcopy(max_list*max_list,0.0d0,0,B,1)
           do i=1,m-1   
           do j=1,m-1   
              B(i,j) = A(i,j)
           end do
           end do

           do i=1,m-1
           do j=1,m-1
              Binv(i,j) = 0.0d0
           end do
           end do
           do i=1,m-1
              Binv(i,i) = 1.0d0
           end do

           call DGESV((m-1),(m-1),
     >                B,max_list,
     >                ipiv,
     >                Binv,max_list,
     >                info)


           !**** Define d ****
           do i=1,m-1
              d(i) = 0.0d0
              do j=1,m-1
                 d(i) = d(i) - c(j)*Binv(j,i)
              end do
           end do

           !**** F0 = F1,  V0 = V1 ****
           call dcopy(nsize,dbl_mb(F1),1,dbl_mb(F0),1)
           call dcopy(nsize,dbl_mb(V1),1,dbl_mb(V0),1)


           !**** V1 = V0 + alpha*F0 - Sum(i=1,m-1) d(i)*U(i) ****
           call daxpy(nsize,alpha1,dbl_mb(F0),1,dbl_mb(V1),1)   

           do i=1,m-1
              call nwpw_list_ptr(1,(5+max_m +i),U)
              call daxpy(nsize,(d(i)),dbl_mb(U),1,dbl_mb(V1),1)   
           end do


           if (m.lt.max_m) then
              m = m+1
           else

             !*** Shift A matrix ****
             do j=1,m-2
             do i=1,m-2
               A(i,j) = A(i+1,j+1)
             end do
             end do

              !*** Shift dF and U  ****
              call nwpw_list_shift_range(1,(5+1),(5+max_m))
              call nwpw_list_shift_range(1,(5+max_m+1),(5+2*max_m))
           end if

           call dcopy(nsize,dbl_mb(V1),1,vnew,1)

         end if

      end if


*     **** Anderson density mixing ****
      if (algorithm.eq.3) then
         
         if ((m.eq.2).or.(m.eq.1).or.(deltae.lt.0.0d0)) then
            call nwpw_list_ptr(1,1,rr_ptr)
            call nwpw_list_ptr(1,2,ff_ptr)
            call nwpw_list_ptr(1,3,tt_ptr)
ccccccccccccccccccc ff=rr
            call dcopy(nsize,dbl_mb(rr_ptr),1,dbl_mb(tt_ptr),1)  ! f = vm
ccccccccccccccccccc ff=ff-vout= vout_old - vout
            call daxpy(nsize,(-1.0d0),vout, 1,dbl_mb(tt_ptr),1)  ! f = f-vout = vm-vout
            call dscal(nsize,(-0.5d0),dbl_mb(tt_ptr),1)          ! f = -f = vout-vm

            scf_error = ddot(nsize,dbl_mb(tt_ptr),1,dbl_mb(tt_ptr),1)
            call D3dB_SumAll(scf_error)
            scf_error = 4.0*dsqrt(scf_error)                     ! scf_error = sqrt(<f|f>)
        
           call dcopy(nsize,dbl_mb(rr_ptr),1,dbl_mb(ff_ptr),1)
           call dcopy(nsize,dbl_mb(rr_ptr),1,vnew,1)            !vnew = vm
           call daxpy(nsize,alpha1,dbl_mb(tt_ptr),1,vnew,1)     !vnew = vnew+alpha*f
           call dcopy(nsize,vnew,1,dbl_mb(rr_ptr),1)            !vm=vnew       
 
         else
     
            call nwpw_list_ptr(1,1,rr_ptr)
            call nwpw_list_ptr(1,2,ff_ptr)
            call nwpw_list_ptr(1,3,tt_ptr)

            call dcopy(nsize,vout,1,vnew,1)
            call daxpy(nsize,(-1.0d0),dbl_mb(rr_ptr),1,vout,1)  ! vout = vout - rr_ptr

            call dcopy(nsize,dbl_mb(rr_ptr),1,dbl_mb(tt_ptr),1) ! tt_ptr = rr_ptr - ff_ptr
            call daxpy(nsize,(-1.0d0),dbl_mb(ff_ptr),1,dbl_mb(tt_ptr),1)

*           ***** generate beta ****
            p00 = ddot(nsize,vout,1,vout,1)                     ! p00 = <vout|vout>
            p01 = ddot(nsize,vout,1,dbl_mb(tt_ptr),1)           ! p01 = <vout|tt_ptr>
            p11 = ddot(nsize,dbl_mb(tt_ptr),1,dbl_mb(tt_ptr),1) ! p11 = <tt_ptr|tt_ptr>
            call D3dB_SumAll(p00)
            call D3dB_SumAll(p01)
            call D3dB_SumAll(p11)
            beta = (p00-p01)/(p00-2.0d0*p01+p11)
            if (beta.lt.1.0d-3) beta = 1.0d-3
            if (beta.gt.1.5d0)  beta = 1.5d0

            call dcopy(nsize,dbl_mb(rr_ptr),1,dbl_mb(ff_ptr),1) ! ff_ptr =rr_ptr
            call dcopy(nsize,vnew,1,dbl_mb(rr_ptr),1)           ! rr_ptr =vnew

            call dscal(nsize,beta,vnew,1)                       ! vnew = beta*vnew + (1-beta)*ff_ptr
            call daxpy(nsize,(1.0d0-beta),dbl_mb(ff_ptr),1,vnew,1)
            scf_error = dsqrt(p00)
         end if
         m = m + 1

      end if

      return 
      end

