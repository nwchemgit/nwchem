*
*     $Id$                       
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

#include "bafdecls.fh"
#include "nwpw_scf_mixing.fh"


!$OMP MASTER
      algorithm = algorithm0
      alpha = alpha0
      max_m = max_m0
      n2ft3d = nsize0
      nsize  = nsize0*ispin0
      ispin  = ispin0
!$OMP END MASTER
!$OMP BARRIER

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
        call nwpw_list_start(1,4,nsize)
        call nwpw_scf_mixing_reset(rho_in)
      end if

*     **** intialize kerker damping ****
      call kerker_G_init()

      return
      end


*     *****************************************
*     *                                       *
*     *        nwpw_scf_mixing_end            *
*     *                                       *
*     *****************************************

      subroutine nwpw_scf_mixing_end()
      implicit none

#include "bafdecls.fh"
#include "nwpw_scf_mixing.fh"

      call nwpw_list_end(1)
      call kerker_G_end()
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

#include "bafdecls.fh"
#include "nwpw_scf_mixing.fh"

!$OMP MASTER
      m = 1
!$OMP END MASTER
      call nwpw_list_store(1,1,rho_in)
      return
      end


*     *****************************************
*     *                                       *
*     *              nwpw_scf_mixing          *
*     *                                       *
*     *****************************************
*
*     shared output variables:
*     - scf_error

      subroutine nwpw_scf_mixing(vout,vnew,deltae,scf_error)
      implicit none
      real*8 vout(*)
      real*8 vnew(*)
      real*8 deltae
      real*8 scf_error


#include "bafdecls.fh"
#include "nwpw_scf_mixing.fh"

*     **** local variables ****
      logical value
      integer i,j,info,k,ipiv,ms,shift
      integer rr_ptr,ss_ptr,tt_ptr,ff_ptr
      integer V0,V1,Vout0,Vout1,F0,F1,Vbar0,Vbar1
      integer dV,U,dF,dFi
      real*8 sum0,sum1,beta,p00,p01,p11,alpha1,r00
      real*8 BB(40,40),BBB(40,40)

*     **** external functions ****
      real*8   ddot
      external ddot

c      if ((deltae.ge.0.0d0)) then
c         alpha1 = alpha 
c      else
c         alpha1 = 1.0d0
c      end if
      alpha1 = alpha 

*     **** simple mixing ****
      if (algorithm.eq.0) then
         call nwpw_list_ptr(1,1,rr_ptr)
         call nwpw_list_ptr(1,2,ff_ptr)
         call Parallel_shared_vector_copy(.true., nsize,
     >                dbl_mb(rr_ptr),dbl_mb(ff_ptr))  ! f = vm
         call DAXPY_OMP(nsize,(-1.0d0),vout, 1,dbl_mb(ff_ptr),1)  ! f = f-vout = vm-vout
         call DSCAL_OMP(nsize,(-1.0d0),dbl_mb(ff_ptr),1)          ! f = -f = vout-vm

!$OMP MASTER
         scf_error = ddot(nsize,dbl_mb(ff_ptr),1,dbl_mb(ff_ptr),1)
!$OMP END MASTER
         call D3dB_SumAll(scf_error)
!$OMP MASTER
         scf_error = dsqrt(scf_error)                         ! scf_error = sqrt(<f|f>)
!$OMP END MASTER

         do ms=1,ispin
            shift = (ms-1)*n2ft3d
            call kerker_G(dbl_mb(ff_ptr+shift))
         end do

         call Parallel_shared_vector_copy(.true.,nsize,
     >                dbl_mb(rr_ptr),vnew)            !vnew = vm
         call DAXPY_OMP(nsize,alpha1,dbl_mb(ff_ptr),1,vnew,1)     !vnew = vnew+alpha*f
         call Parallel_shared_vector_copy(.true.,nsize,
     >                vnew,dbl_mb(rr_ptr))            !vm=vnew
      end if

*     **** Broyden mixing ****
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
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  vout,dbl_mb(Vout0))

           !**** F0 = vout - V0 ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(V0),dbl_mb(F0))
           call DSCAL_OMP(nsize,(-1.0d0),dbl_mb(F0),1)
           call DAXPY_OMP(nsize,(1.0d0),vout, 1,dbl_mb(F0),1)

           !**** V1 = V0 + alpha*F0 ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(V0),dbl_mb(V1))
           call DAXPY_OMP(nsize,alpha1,dbl_mb(F0),1,dbl_mb(V1),1)   

           !**** vnew = V1 ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(V1),vnew)

         else !(m.gt.1) 
           !**** vout1 = vout ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  vout,dbl_mb(Vout1))

           !**** F1 = vout - V1 ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(V1),dbl_mb(F1))
           call DSCAL_OMP(nsize,(-1.0d0),dbl_mb(F1),1)
           call DAXPY_OMP(nsize,(1.0d0),vout, 1,dbl_mb(F1),1)

           !**** scf_error = sqrt(<F1|F1>) ****
!$OMP MASTER
           scf_error = ddot(nsize,dbl_mb(F1),1,dbl_mb(F1),1)
!$OMP END MASTER
           call D3dB_SumAll(scf_error)
!$OMP MASTER
           scf_error = dsqrt(scf_error)          
!$OMP END MASTER

       
           !**** Beta = <F1|F1-F0>/<F1-F0/F1-F0> ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(F1),dbl_mb(Vbar1))
           call DAXPY_OMP(nsize,(-1.0d0),dbl_mb(F0),1,dbl_mb(Vbar1),1)
           sum0 = ddot(nsize,dbl_mb(F1),1,dbl_mb(Vbar1),1)
           sum1 = ddot(nsize,dbl_mb(Vbar1),1,dbl_mb(Vbar1),1)
           call D3dB_SumAll(sum0)
           call D3dB_SumAll(sum1)
           beta = sum0/sum1
           
 

           !**** Vbar1 = (1-Beta)*Vout1 + Beta*Vout0 ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(Vout0),dbl_mb(Vbar1))
           call DSCAL_OMP(nsize,(beta),dbl_mb(Vbar1),1)
           call DAXPY_OMP(nsize,(1.0d0-beta),dbl_mb(Vout1),1,
     >                                   dbl_mb(Vbar1),1)

           !**** Vbar0 = (1-Beta)*V1 + Beta*V0       ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(V0),dbl_mb(Vbar0))
           call DSCAL_OMP(nsize,(beta),dbl_mb(Vbar0),1)
           call DAXPY_OMP(nsize,(1.0d0-beta),
     >                    dbl_mb(V1),1,dbl_mb(Vbar0),1)

           !**** F0 = F1, Vout0 = Vout1, V0 = V1 ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(F1),dbl_mb(F0))
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(Vout1),dbl_mb(Vout0))
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(V1),dbl_mb(V0))

           !**** V1 = (1-alpha)*Vbar0 + alpha*Vbar1 ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(Vbar1),dbl_mb(V1))
           call DSCAL_OMP(nsize,(alpha1),dbl_mb(V1),1)
           call DAXPY_OMP(nsize,(1.0d0-alpha1),dbl_mb(Vbar0),1,
     >                                         dbl_mb(V1),1)
          

           !**** vnew = V1 ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(V1),vnew)

         end if

!$OMP MASTER 
         m = m + 1
!$OMP END MASTER 

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
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(V0),dbl_mb(F0))
           call DSCAL_OMP(nsize,(-1.0d0),dbl_mb(F0),1)
           call DAXPY_OMP(nsize,(1.0d0),vout,1,dbl_mb(F0),1)

           !**** scf_error = sqrt(<F1|F1>) ****
!$OMP MASTER 
           scf_error = ddot(nsize,dbl_mb(F0),1,dbl_mb(F0),1)
!$OMP END MASTER 
           call D3dB_SumAll(scf_error)
!$OMP MASTER 
           scf_error = dsqrt(scf_error)          
!$OMP END MASTER 

           !**** V1 = V0 + alpha*F0 ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(V0),dbl_mb(V1))

           call DAXPY_OMP(nsize,alpha1,dbl_mb(F0),1,dbl_mb(V1),1)   

           !**** vnew = V1 ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(V1),vnew)

           call Parallel_shared_vector_zero(.true.,max_list*max_list,A)

!$OMP MASTER 
           m = m + 1
!$OMP END MASTER 

         else !(m.gt.1) 
           !**** F1 = vout - V1 ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(V1),dbl_mb(F1))
           call DSCAL_OMP(nsize,(-1.0d0),dbl_mb(F1),1)
           call DAXPY_OMP(nsize,(1.0d0),vout,1,dbl_mb(F1),1)

           !**** scf_error = sqrt(<F1|F1>) ****
!$OMP MASTER
           scf_error = ddot(nsize,dbl_mb(F1),1,dbl_mb(F1),1)
!$OMP END MASTER
           call D3dB_SumAll(scf_error)
!$OMP MASTER
           scf_error = dsqrt(scf_error)          
!$OMP END MASTER

           !*** dF = dF(m-1), U = U(m-1) ***
           call nwpw_list_ptr(1,(5+m-1),dF)
           call nwpw_list_ptr(1,(5+max_m  +m-1),U)

           !**** dF = (F1-F0) ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(F1),dbl_mb(dF))
           call DAXPY_OMP(nsize,(-1.0d0),dbl_mb(F0),1,dbl_mb(dF),1)

           !**** dV  = (V1-V0) ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(V1),dbl_mb(dV))
           call DAXPY_OMP(nsize,(-1.0d0),dbl_mb(V0),1,dbl_mb(dV),1)
           
           !**** U  = alpha*dF + dV ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(dV),dbl_mb(U))

           call DAXPY_OMP(nsize,(alpha1),dbl_mb(dF),1,dbl_mb(U),1)

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
           
           call Parallel_shared_vector_zero(.true.,max_list*max_list,B)
!$OMP MASTER
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
!$OMP END MASTER

           !**** F0 = F1,  V0 = V1 ****
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(F1),dbl_mb(F0))
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(V1),dbl_mb(V0))


           !**** V1 = V0 + alpha*F0 - Sum(i=1,m-1) d(i)*U(i) ****
           call DAXPY_OMP(nsize,alpha1,dbl_mb(F0),1,dbl_mb(V1),1)   

           do i=1,m-1
              call nwpw_list_ptr(1,(5+max_m +i),U)
              call DAXPY_OMP(nsize,(d(i)),dbl_mb(U),1,dbl_mb(V1),1)   
           end do


           if (m.lt.max_m) then
!$OMP MASTER
              m = m+1
!$OMP END MASTER
           else

             !*** Shift A matrix ****
!$OMP MASTER
             do j=1,m-2
             do i=1,m-2
               A(i,j) = A(i+1,j+1)
             end do
             end do
!$OMP END MASTER

              !*** Shift dF and U  ****
              call nwpw_list_shift_range(1,(5+1),(5+max_m))
              call nwpw_list_shift_range(1,(5+max_m+1),(5+2*max_m))
           end if

           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(V1),vnew)

         end if

      end if


*     **** Anderson density mixing ****
      if (algorithm.eq.3) then
c         if ((deltae.ge.0.0d0)) then
c            alpha1 = alpha 
c         else
c            alpha1 = 1.0d0
c         end if
         if ((m.eq.2).or.(m.eq.1).or.(deltae.gt.0.0d0)) then
            call nwpw_list_ptr(1,1,rr_ptr)
            call nwpw_list_ptr(1,2,ff_ptr)
            call nwpw_list_ptr(1,3,tt_ptr)

c           *** ff=vout-vm ***
            call Parallel_shared_vector_copy(.true.,nsize,
     >                   vout,dbl_mb(ff_ptr))  ! f = vout
            call DAXPY_OMP(nsize,(-1.0d0),dbl_mb(rr_ptr),1,
     >                                    dbl_mb(ff_ptr),1)  ! f = f-vm = vout-vm

!$OMP MASTER
            scf_error = ddot(nsize,dbl_mb(ff_ptr),1,dbl_mb(ff_ptr),1)
!$OMP END MASTER
            call D3dB_SumAll(scf_error)
!$OMP MASTER
            scf_error = dsqrt(scf_error)                         ! scf_error = sqrt(<f|f>)
!$OMP END MASTER

            do ms=1,ispin
               shift = (ms-1)*n2ft3d
               call kerker_G(dbl_mb(ff_ptr+shift))
            end do

           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(rr_ptr),vnew)            !vnew = vm
           call DAXPY_OMP(nsize,alpha1,dbl_mb(ff_ptr),1,vnew,1)     !vnew = vnew+alpha*f

           call Parallel_shared_vector_copy(.true.,nsize,
     >                  vout,dbl_mb(tt_ptr))           !vout1 = vout
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  dbl_mb(rr_ptr),dbl_mb(ff_ptr)) !vm1   = vm
           call Parallel_shared_vector_copy(.true.,nsize,
     >                  vnew,dbl_mb(rr_ptr))           !vm    = vnew       
 
         else
            call nwpw_list_ptr(1,1,rr_ptr)
     
            call nwpw_list_ptr(1,1,rr_ptr) !vm
            call nwpw_list_ptr(1,2,ss_ptr) !vm1
            call nwpw_list_ptr(1,3,tt_ptr) !vout1
            call nwpw_list_ptr(1,4,ff_ptr)

            call Parallel_shared_vector_copy(.true.,nsize,
     >                   vout,dbl_mb(ff_ptr))
            call DAXPY_OMP(nsize,(-1.0d0),dbl_mb(rr_ptr),1,
     >                                    dbl_mb(ff_ptr),1)  ! ff_ptr = vout - rr_ptr
            call DAXPY_OMP(nsize,(-1.0d0),dbl_mb(ss_ptr),1,
     >                                    dbl_mb(tt_ptr),1)  ! tt_ptr = vout1 - vm1

!$OMP MASTER
            scf_error = 0.0d0
            scf_error = ddot(nsize,dbl_mb(ff_ptr),1,dbl_mb(ff_ptr),1)
!$OMP END MASTER
            call D3dB_SumAll(scf_error)
!$OMP MASTER
            scf_error = dsqrt(scf_error)                         ! scf_error = sqrt(<f|f>)
!$OMP END MASTER

            do ms=1,ispin
               shift = (ms-1)*n2ft3d

               call kerker_G(dbl_mb(ff_ptr+shift))
               call kerker_G(dbl_mb(tt_ptr+shift))

*              ***** generate beta ****
               p00 = ddot(n2ft3d,dbl_mb(ff_ptr+shift),1,
     >                           dbl_mb(ff_ptr+shift),1)                 ! p00 = <ff_ptr|ff_ptr>

               p01 = ddot(n2ft3d,dbl_mb(ff_ptr+shift),1,
     >                           dbl_mb(tt_ptr+shift),1)                 ! p01 = <ff_ptr|tt_ptr>

               p11 = ddot(n2ft3d,dbl_mb(tt_ptr+shift),1,
     >                           dbl_mb(tt_ptr+shift),1)                 ! p11 = <tt_ptr|tt_ptr>
               call D3dB_SumAll(p00)
               call D3dB_SumAll(p01)
               call D3dB_SumAll(p11)
               r00 = p00-2.0d0*p01+p11
               beta = (p00-p01)/r00
                
               if ((r00.lt.0.0d0).or.(beta.lt.-1.5d0)) beta = 1.0d-3
               if (beta.gt.1.5d0)  beta = 1.5d0

               !*** vnew = (1-beta)*vm + beta*vm1 + alpha*((1-beta)*fm + beta*fm1) ***
               call Parallel_shared_vector_copy(.true.,n2ft3d,
     >                      dbl_mb(rr_ptr+shift),vnew(1+shift))
               call DSCAL_OMP(n2ft3d,(1.0d0-beta),vnew(1+shift),1)
               call DAXPY_OMP(n2ft3d,(beta),dbl_mb(ss_ptr+shift),1,
     >                                  vnew(1+shift),1)
               call DSCAL_OMP(n2ft3d,(1.0d0-beta),
     >                        dbl_mb(ff_ptr+shift),1)
               call DAXPY_OMP(n2ft3d,(beta),dbl_mb(tt_ptr+shift),1,
     >                                dbl_mb(ff_ptr+shift),1)
               call DAXPY_OMP(n2ft3d,alpha1,dbl_mb(ff_ptr+shift),1,
     >                                  vnew(1+shift),1)
            end do
            call Parallel_shared_vector_copy(.true.,nsize,
     >                   dbl_mb(rr_ptr),dbl_mb(ss_ptr))
            call Parallel_shared_vector_copy(.true.,nsize,
     >                   vnew,dbl_mb(rr_ptr))
            call Parallel_shared_vector_copy(.true.,nsize,
     >                   vout,dbl_mb(tt_ptr))
         end if
!$OMP MASTER
         m = m + 1
!$OMP END MASTER

      end if

      return 
      end

