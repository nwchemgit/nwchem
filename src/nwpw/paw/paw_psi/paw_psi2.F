*
* $Id$
*


************************ f_orb orbitals Part ************************

************************ virtural orbital Part ************************


************************ KS orbital Part ************************

*     ***********************************
*     *                                 *
*     *      paw_psi_KS_update          *
*     *                                 *
*     ***********************************

*    This routine (approximately) diagonalizes the KS matrix.
*
      subroutine paw_psi_KS_update(paw_psi_number,precondition,maxerror)
      implicit none
      integer paw_psi_number
      logical precondition
      real*8 maxerror
        
#include "bafdecls.fh"
#include "paw_psi.fh"
    
*     **** local variables ****
      logical done
      integer i,j,neall,maxit_orb,maxit_orbs
      real*8 error,error_out,tim1,tim2,tim,sum

*     **** external functions ****
      integer  control_ks_maxit_orb,control_ks_maxit_orbs
      external control_ks_maxit_orb,control_ks_maxit_orbs

      tim = 0.0d0
      neall = ne(1)+ne(2)
      maxit_orb  = control_ks_maxit_orb()   !*** should be read from rtdb ***
      maxit_orbs = control_ks_maxit_orbs()  !*** should be read from rtdb ***
      j = 0
 2    j = j+1
        error = 0.0d0
        !do i=neall,1,-1
        do i=1,neall
          call current_second(tim1)

         !*** orthogonalize to lower orbitals  ****
!        call paw_psi_project_out_f_orb1(
!    >           i,
!    >           dcpl_mb(psi1(1)+(i-1)*npack1))

         !*** normalize ****
         call Pack_cc_dot(1,
     >            dcpl_mb(psi1(1) +(i-1)*npack1),
     >            dcpl_mb(psi1(1) +(i-1)*npack1),
     >            sum)
         sum = 1.0d0/dsqrt(sum)
c         call Pack_c_SMul(1,sum,
c     >            dcpl_mb(psi1(1) +(i-1)*npack1),
c     >            dcpl_mb(psi1(1) +(i-1)*npack1))
         call Pack_c_SMul1(1,sum,dcpl_mb(psi1(1) +(i-1)*npack1))



          call paw_psi_KS_update_orb(paw_psi_number,precondition,
     >                         maxit_orb,maxerror,
     >                         0.1d0,i,error_out)
          call current_second(tim2)
          tim = tim + (tim2-tim1)
          error = error+error_out
        end do
        error = error/dble(neall)

        done = ((j.gt.maxit_orbs).or.(error.lt.maxerror))
      if (.not.done) go to 2
c      write(*,*) "TIME ALL=",tim

      return
      end


*     ***********************************
*     *                                 *
*     *      paw_psi_KS_update_orb      *
*     *                                 *
*     ***********************************

*    This routine performs a KS update on orbital i
*
      subroutine paw_psi_KS_update_orb(paw_psi_number,
     >                             precondition,maxiteration,
     >                             maxerror,perror,i,
     >                             error_out)
      implicit none
#include "errquit.fh"
      integer paw_psi_number
      logical precondition
      integer maxiteration
      real*8  maxerror,perror 
      integer i
      real*8 error_out
      
#include "bafdecls.fh"
#include "paw_psi.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,done,oneloop
      integer it
      real*8 e0,eold,percent_error,error0,de0,lmbda_r0,lmbda_r1
      real*8 theta,sigma
      integer r1(2),t0(2),t(2),g(2)
      integer paw_psi_ptr

      if (paw_psi_number.eq.1) then
         paw_psi_ptr=psi1(1)
      else 
         paw_psi_ptr=psi2(1)
      end if

      call Parallel_taskid(taskid)

      value = BA_push_get(mt_dcpl,npack1,'t0',t0(2),t0(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'g',g(2),g(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'t',t(2),t(1))
      if (.not. value) call errquit(
     >     'paw_psi_KS_update_orb: out of stack memory',0, MA_ERR)

      done = .false.
      error0 = 0.0d0 
      e0 = 0.0d0
      theta = -3.14159d0/600.0d0
      it = 0
 2    continue

         it = it + 1
         eold = e0

*        *** calculate residual (steepest descent) direction for a single band ***
         call paw_psi_get_gradient_orb(paw_psi_number,i,dcpl_mb(g(1)))
         call Pack_cc_dot(1,dcpl_mb(paw_psi_ptr+(i-1)*npack1),
     >                   dcpl_mb(g(1)),
     >                    e0)
         e0 = -e0
         

         done = ((it.gt.maxiteration) 
     >           .or.
     >           (dabs(e0-eold).lt.maxerror))

         if (done) go to 4

*        **** preconditioning ****
         if (precondition) then
           call ke_Precondition(npack1,1,
     >                     dcpl_mb(g(1)),
     >                     dcpl_mb(g(1)))
       
         end if

c        call Pack_c_Copy(1,dcpl_mb(g(1)),dcpl_mb(r1(1)))
c        call Pack_cc_daxpy(1,(e0),
c    >                 dcpl_mb(paw_psi_ptr+(i-1)*npack1),
c    >                 dcpl_mb(r1(1)))
         call Pack_c_Copy(1,dcpl_mb(g(1)),dcpl_mb(r1(1)))
         call paw_psi_project_out_orb(paw_psi_number,i,dcpl_mb(r1(1)))

     


*        *** determine conjuagate direction ***
         call Pack_cc_dot(1,dcpl_mb(r1(1)),
     >                   dcpl_mb(r1(1)),
     >                   lmbda_r1)
         call Pack_c_Copy(1,dcpl_mb(r1(1)),dcpl_mb(t(1)))
      
         if (it.gt.1) then
         call Pack_cc_daxpy(1,(lmbda_r1/lmbda_r0),
     >                   dcpl_mb(t0(1)),
     >                   dcpl_mb(t(1)))       
         end if
         lmbda_r0 = lmbda_r1
         oneloop = .true.
 3       call Pack_c_Copy(1,dcpl_mb(t(1)),dcpl_mb(t0(1)))



c!*        **** project out psi components from t ****
c!        call paw_psi_project_out_orb(paw_psi_number,i,dcpl_mb(t(1)))
c!        call Pack_cc_dot(1,dcpl_mb(paw_psi_ptr+(i-1)*npack1),
c!    >                   dcpl_mb(t(1)),
c!    >                    de0)
c!        de0 = -de0
c!        call Pack_cc_daxpy(1,(de0),
c!    >                 dcpl_mb(paw_psi_ptr+(i-1)*npack1),
c!    >                 dcpl_mb(t(1)))


*        *** normalize search direction, t ****
         call Pack_cc_dot(1,dcpl_mb(t(1)),
     >                   dcpl_mb(t(1)),
     >                   sigma)
         sigma = dsqrt(sigma)
         de0 = 1.0d0/sigma
c         call Pack_c_SMul(1,de0,dcpl_mb(t(1)),dcpl_mb(t(1)))
         call Pack_c_SMul1(1,de0,dcpl_mb(t(1)))


*        **** compute de0 = <t|g> ****
         call Pack_cc_dot(1,dcpl_mb(t(1)),
     >                   dcpl_mb(g(1)),
     >                   de0)

*        *** bad direction ***
         if ((de0.lt.0.0d0).and.oneloop) then
           call Pack_c_Copy(1,dcpl_mb(g(1)),dcpl_mb(t(1)))
           oneloop = .false.
           go to 3
         end if

         de0 = -2.0d0*de0
         call paw_psi_linesearch_update2(paw_psi_number,i,
     >                              theta,e0,de0,
     >                              dcpl_mb(t(1)),
     >                              sigma,
     >                              dcpl_mb(t0(1)))

      go to 2


*     **** release stack memory ****
 4    value =           BA_pop_stack(t(2)) 
      value = value.and.BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(r1(2))
      value = value.and.BA_pop_stack(t0(2))
      if (.not. value) call errquit(
     >     'paw_psi_KS_update_orb: popping stack memory',1, MA_ERR)

c      write(*,*) "iterations=",it," eig=",e0," error=",error_out,
c     >           theta
      error_out = dabs(e0-eold)
      return
      end



*     ***********************************
*     *                                 *
*     *      paw_psi_linesearch_update  *
*     *                                 *
*     ***********************************

*    This routine performs a linesearch on orbital i, in the direction t.  
* This routine is needed for a KS minimizer.
*  e0 = <orb|g>
*  de0 = 2*<t|g>
*
      subroutine paw_psi_linesearch_update(paw_psi_number,i,
     >                                     theta,e0,de0,t)
      implicit none
      integer paw_psi_number
      integer i
      real*8  theta
      real*8  e0,de0
      complex*16 t(*) !search direction
      
#include "bafdecls.fh"
#include "paw_psi.fh"
#include "errquit.fh"
   
*     **** local variables ****
      logical value
      integer orb(2),g(2),paw_psi_ptr
      real*8 x,y,pi,e1,de1
      real*8 theta2,e2,de2

      if (paw_psi_number.eq.1) then
         paw_psi_ptr=psi1(1)
      else 
         paw_psi_ptr=psi2(1)
      end if

      pi = 4.0d0*datan(1.0d0)

*     **** allocate stack memory ****
      value = BA_push_get(mt_dcpl,npack1,'orb',
     >                       orb(2),orb(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'g',
     >                       g(2),g(1))
      if (.not. value) call errquit(
     >     'paw_psi_linesearch_update: out of stack memory',0, MA_ERR)
 

      call Pack_c_Copy(1,dcpl_mb(paw_psi_ptr+(i-1)*npack1),
     >                   dcpl_mb(orb(1)))

*     **** orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
      !theta = pi/300.0d0
      x = cos(theta)
      y = sin(theta)
      call Pack_c_SMul(1,x,
     >                  dcpl_mb(orb(1)),
     >                  dcpl_mb(paw_psi_ptr+(i-1)*npack1))
      call Pack_cc_daxpy(1,y,
     >                   t,
     >                   dcpl_mb(paw_psi_ptr+(i-1)*npack1))

*     *** determine theta ***
      call paw_psi_get_gradient_orb(paw_psi_number,i,dcpl_mb(g(1)))

      call Pack_cc_dot(1,dcpl_mb(paw_psi_ptr+(i-1)*npack1),
     >                   dcpl_mb(g(1)),
     >                   e1)
      e1 = -e1
      x = (e0 - e1 + 0.5d0*de0*sin(2*theta))
     >    /(1.0d0-cos(2*theta))
      theta = 0.5d0*datan(0.5d0*de0/x) 

c     call Pack_cc_dot(1,t,
c    >                 dcpl_mb(g(1)),
c    >                 de1)
c     de1 = -2.0d0*de1
c     theta  = -de1*(pi/300.0d0)/(de1-de0)

      !write(*,*) "i,theta,e1:",i,theta,e1


*     **** orb2 = orb*cos(theta) + t*sin(theta) ****
      x = cos(theta)
      y = sin(theta)
      call Pack_c_SMul(1,x,
     >                  dcpl_mb(orb(1)),
     >                  dcpl_mb(paw_psi_ptr+(i-1)*npack1))
      call Pack_cc_daxpy(1,y,
     >                   t,
     >                   dcpl_mb(paw_psi_ptr+(i-1)*npack1))

*     **** update orb2_r and H*orb2 ****
      !call paw_electron_run_orb(i,dcpl_mb(paw_psi_ptr)) 
c     call paw_psi_get_gradient_orb(paw_psi_number,i,dcpl_mb(g(1)))
c     call Pack_cc_dot(1,dcpl_mb(paw_psi_ptr+(i-1)*npack1),
c    >                   dcpl_mb(g(1)),
c    >                   e2)
c     e2 = -e2
c     call Pack_cc_dot(1,t,
c    >                 dcpl_mb(g(1)),
c    >                 de2)
c     de2 = -2.0d0*de2

c     write(*,*) "i,theta,es:",i,theta,e0,e1,e2
c     write(*,*)

*     **** release stack memory ****
      value =           BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(orb(2))      
      if (.not. value) call errquit(
     >     'paw_psi_linesearch_update: popping stack memory',1, MA_ERR)

      return
      end

*     ***********************************
*     *                                 *
*     *   paw_psi_linesearch_update2    *
*     *                                 *
*     ***********************************

*    This routine performs a linesearch on orbital i, in the direction t.  
* This routine is needed for a KS minimizer.
*  e0 = <orb|g>
*  de0 = 2*<t|g>
*
      subroutine paw_psi_linesearch_update2(paw_psi_number,
     >                                  i,theta,e0,de0,t,
     >                                  sigma,tau_t)
      implicit none
#include "errquit.fh"
      integer paw_psi_number
      integer i
      real*8  theta
      real*8  e0,de0
      complex*16 t(*)     !search direction

      real*8     sigma
      complex*16 tau_t(*) !parallel transported search direction
      
#include "bafdecls.fh"
#include "paw_psi.fh"
   
*     **** local variables ****
      logical value
      integer orb(2),g(2),paw_psi_ptr
      real*8 x,y,pi,e1,de1
      real*8 theta2,e2,de2

      if (paw_psi_number.eq.1) then
         paw_psi_ptr=psi1(1)
      else 
         paw_psi_ptr=psi2(1)
      end if

      pi = 4.0d0*datan(1.0d0)

*     **** allocate stack memory ****
      value = BA_push_get(mt_dcpl,npack1,'orb',
     >                       orb(2),orb(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'g',
     >                       g(2),g(1))
      if (.not. value) call errquit(
     >     'paw_psi_linesearch_update: out of stack memory',0,MA_ERR)
 

      call Pack_c_Copy(1,dcpl_mb(paw_psi_ptr+(i-1)*npack1),
     >                   dcpl_mb(orb(1)))

*     **** orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
      !theta = pi/300.0d0
      x = cos(theta)
      y = sin(theta)
      call Pack_c_SMul(1,x,
     >                  dcpl_mb(orb(1)),
     >                  dcpl_mb(paw_psi_ptr+(i-1)*npack1))
      call Pack_cc_daxpy(1,y,
     >                   t,
     >                   dcpl_mb(paw_psi_ptr+(i-1)*npack1))

*     *** determine theta ***
      call paw_psi_get_gradient_orb(paw_psi_number,i,dcpl_mb(g(1)))

      call Pack_cc_dot(1,dcpl_mb(paw_psi_ptr+(i-1)*npack1),
     >                   dcpl_mb(g(1)),
     >                   e1)
      e1 = -e1
      x = (e0 - e1 + 0.5d0*de0*sin(2*theta))
     >    /(1.0d0-cos(2*theta))
      theta = 0.5d0*datan(0.5d0*de0/x) 

      x = cos(theta)
      y = sin(theta)

*     **** tau_t = (-orb*sin(theta) + t*cos(theta))*sigma ****
      call Pack_c_SMul(1,(-y),
     >                  dcpl_mb(orb(1)),
     >                  tau_t)
      call Pack_cc_daxpy(1,x,
     >                   t,
     >                   tau_t)
c      call Pack_c_SMul(1,sigma,
c     >                  tau_t,
c     >                  tau_t)
      call Pack_c_SMul1(1,sigma,tau_t)

*     **** orb2 = orb*cos(theta) + t*sin(theta) ****
      call Pack_c_SMul(1,x,
     >                  dcpl_mb(orb(1)),
     >                  dcpl_mb(paw_psi_ptr+(i-1)*npack1))
      call Pack_cc_daxpy(1,y,
     >                   t,
     >                   dcpl_mb(paw_psi_ptr+(i-1)*npack1))


*     **** release stack memory ****
      value =           BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(orb(2))      
      if (.not. value) call errquit(
     >     'paw_psi_linesearch_update: popping stack memory',1,MA_ERR)

      return
      end



*     ***************************
*     *                         *
*     *      paw_psi_set_orb    *
*     *                         *
*     ***************************

*    This routine copies an orbital, orb, into the ith psi of psi1.  
* This routine is needed for a KS minimizer.
*
      subroutine paw_psi_set_orb(paw_psi_number,i,orb)
      implicit none
      integer paw_psi_number
      integer i
      complex*16 orb(*)

#include "bafdecls.fh"
#include "paw_psi.fh"
   
*     **** local variables ****
      integer index,paw_psi_ptr

      if (paw_psi_number.eq.1) then
         paw_psi_ptr=psi1(1)
      else 
         paw_psi_ptr=psi2(1)
      end if

      index = (i-1)*npack1

      call zcopy(npack1,
     >           orb, 1,
     >           dcpl_mb(paw_psi_ptr+index),1)
      return
      end


*     ***************************
*     *                         *
*     *      paw_psi_get_orb    *
*     *                         *
*     ***************************

*    This routine copies the ith psi of psi1 into an orbital, orb.  
* This routine is needed for a KS minimizer.
*
      subroutine paw_psi_get_orb(paw_psi_number,i,orb)
      implicit none
      integer paw_psi_number
      integer i
      complex*16 orb(*)

#include "bafdecls.fh"
#include "paw_psi.fh"
   
*     **** local variables ****
      integer index,paw_psi_ptr


      if (paw_psi_number.eq.1) then
         paw_psi_ptr=psi1(1)
      else 
         paw_psi_ptr=psi2(1)
      end if

      index = (i-1)*npack1

      call zcopy(npack1,
     >           dcpl_mb(paw_psi_ptr+index), 1,
     >           orb, 1)
      return
      end

*     ***********************************
*     *                                 *
*     *      paw_psi_get_gradient_orb   *
*     *                                 *
*     ***********************************

*    This routine returns the Hpsi(i).  
* This routine is needed for a KS minimizer.
*
      subroutine paw_psi_get_gradient_orb(paw_psi_number,i,Horb)
      implicit none
      integer paw_psi_number
      integer i
      complex*16 Horb(*)

#include "bafdecls.fh"
#include "paw_psi.fh"

*     **** local variables ****
      integer paw_psi_ptr

      if (paw_psi_number.eq.1) then
         paw_psi_ptr=psi1(1)
      else 
         paw_psi_ptr=psi2(1)
      end if

      call paw_electron_run_orb(i,dcpl_mb(paw_psi_ptr))
      call paw_electron_get_gradient_orb(i,Horb)
      
      return
      end


*     *******************************************
*     *                                         *
*     *          paw_psi_project_out_orb           *
*     *                                         *
*     *******************************************
*
*    This routine projects out non-orthogonal components of Horb.  
* This routine is needed for a KS minimizer.
*
      subroutine paw_psi_project_out_orb(paw_psi_number,i,Horb)
      implicit none
#include "errquit.fh"
      integer paw_psi_number
      integer i
      complex*16 Horb(*)

#include "bafdecls.fh"
#include "paw_psi.fh"

*     **** local variables ****
      logical ok
      integer ii,n,paw_psi_ptr,np
      integer x(2)
      real*8  sum

      call Parallel_np(np)

*     **** allocate stack memory ****
      ok = BA_push_get(mt_dbl,ne(1),'x',x(2),x(1))
      if (.not.ok) 
     > call errquit('paw_psi_project_out_orb: out of stack memory',0,
     &       MA_ERR)


      if (paw_psi_number.eq.1) then
         paw_psi_ptr=psi1(1)
      else 
         paw_psi_ptr=psi2(1)
      end if

*     **** spin up orbital **** 
      if (i.le.ne(1)) then

        ii = i
!       do n=1,(ii)
!          call Pack_cc_dot(1,
!    >            dcpl_mb(paw_psi_ptr +(n-1)*npack1),
!    >            Horb,
!    >            sum)
!          call daxpy(2*npack1,
!    >               (-sum),
!    >               dcpl_mb(paw_psi_ptr+(n-1)*npack1),1,
!    >               Horb,1) 
!       end do     
        call Pack_cc_ndot(1,ii,
     >            dcpl_mb(paw_psi_ptr),
     >            Horb,
     >            dbl_mb(x(1)))
        do n=1,(ii)
           call daxpy(2*npack1,
     >               (-dbl_mb(x(1)+n-1)),
     >               dcpl_mb(paw_psi_ptr+(n-1)*npack1),1,
     >               Horb,1) 
        end do     



*     **** spin down orbital ****      
      else       


        ii = i - ne(1)
        do n=(ne(1)+1),(ne(1)+ii)
           call Pack_cc_dot(1,
     >            dcpl_mb(paw_psi_ptr +(n-1)*npack1),
     >            Horb,
     >            sum)
           call daxpy(2*npack1,
     >               (-sum),
     >               dcpl_mb(paw_psi_ptr+(n-1)*npack1),1,
     >               Horb,1) 
        end do     


      end if

*     **** release stack memory ****
      ok = BA_pop_stack(x(2))
      if (.not. ok) 
     > call errquit('paw_psi_project_out_orb: poping stack memory',0,
     &       MA_ERR)
 
      return
      end






*     ***************************
*     *                         *
*     *    paw_psi_set_density  *
*     *                         *
*     ***************************

*    This routine sets the densities and potentials in psi and electron.  
* This routine is needed for a band by band minimizer.
*
      subroutine paw_psi_set_density(paw_psi_number,rho)
      implicit none
      integer paw_psi_number
      real*8 rho(*)


#include "bafdecls.fh"
#include "paw_psi.fh"
   

*     **** local variables ****
      integer dng_ptr,rho_ptr

      if (paw_psi_number.eq.1) then
        rho_ptr     = rho1(1)
        dng_ptr     = dng1(1)
      else
        rho_ptr     = rho2(1)
        dng_ptr     = dng2(1)
      end if


      call dcopy(4*nfft3d,
     >           rho, 1,
     >           dbl_mb(rho_ptr),1)

      call paw_electron_gen_dng(dbl_mb(rho_ptr),
     >                            dcpl_mb(dng_ptr))
      call paw_electron_gen_scf_potentials(dbl_mb(rho_ptr),
     >                            dcpl_mb(dng_ptr))
      call paw_electron_gen_vall()
      return
      end


*     ***************************
*     *                         *
*     *    paw_psi_get_density  *
*     *                         *
*     ***************************

*    This routine gets the densities in psi.  
* This routine is needed for a band by band minimizer.
*
      subroutine paw_psi_get_density(paw_psi_number,rho)
      implicit none
      integer paw_psi_number
      real*8 rho(*)


#include "bafdecls.fh"
#include "paw_psi.fh"
   
*     **** local variables ****
      integer rho_ptr

      if (paw_psi_number.eq.1) then
        rho_ptr = rho1(1)
      else
        rho_ptr = rho2(1)
      end if

      call dcopy(4*nfft3d,
     >           dbl_mb(rho_ptr),1,
     >           rho,1)
      return
      end


*     **************************************
*     *                                    *
*     *   paw_psi_gen_density_potentials   *
*     *                                    *
*     **************************************

*    This routine sets the densities and potentials in psi and electron.  
* This routine is needed for a band by band minimizer.
*
      subroutine paw_psi_gen_density_potentials(paw_psi_number)
      implicit none
      integer paw_psi_number


#include "bafdecls.fh"
#include "paw_psi.fh"
   
*     **** local variables ****
      integer paw_psi_ptr,rho_ptr,dng_ptr

      if (paw_psi_number.eq.1) then
        paw_psi_ptr = psi1(1)
        rho_ptr     = rho1(1)
        dng_ptr     = dng1(1)
      else
        paw_psi_ptr = psi2(1)
        rho_ptr     = rho2(1)
        dng_ptr     = dng2(1)
      end if


      call paw_electron_gen_densities(dcpl_mb(paw_psi_ptr),
     >                            dbl_mb(rho_ptr),
     >                            dcpl_mb(dng_ptr))
      call paw_electron_gen_scf_potentials(dbl_mb(rho_ptr),
     >                            dcpl_mb(dng_ptr))
      call paw_electron_gen_vall()
      return
      end


************************ Grasmman orbitals Part ************************

*     ***************************
*     *                         *
*     *         paw_psi_1to2    *
*     *                         *
*     ***************************
      subroutine paw_psi_1to2()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"
   
      call zcopy(npack1*(ne(1)+ne(2)),
     >           dcpl_mb(psi1(1)),1,
     >           dcpl_mb(psi2(1)),1)

      return
      end


*     ***************************
*     *                         *
*     *         paw_psi_2to1    *
*     *                         *
*     ***************************
      subroutine paw_psi_2to1()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"

 
      call zcopy(npack1*(ne(1)+ne(2)),
     >           dcpl_mb(psi2(1)),1,
     >           dcpl_mb(psi1(1)),1)

c      call OrthoCheck(ispin,ne,dcpl_mb(psi1(1)))  
      return
      end


*     ***************************
*     *                         *
*     *     paw_psi_1get_psi    *
*     *                         *
*     ***************************
      subroutine paw_psi_1get_psi(rpsi)
      implicit none
      complex*16 rpsi(*)

#include "bafdecls.fh"
#include "paw_psi.fh"
   
      call zcopy(npack1*(ne(1)+ne(2)),
     >           dcpl_mb(psi1(1)),1,
     >           rpsi,1)

      return
      end


*     ***************************
*     *                         *
*     *     paw_psi_2get_psi    *
*     *                         *
*     ***************************
      subroutine paw_psi_2get_psi(rpsi)
      implicit none
      complex*16 rpsi(*)

#include "bafdecls.fh"
#include "paw_psi.fh"
   
      call zcopy(npack1*(ne(1)+ne(2)),
     >           dcpl_mb(psi2(1)),1,
     >           rpsi,1)

      return
      end

*     ***************************
*     *                         *
*     *         paw_psi_check   *
*     *                         *
*     ***************************
      subroutine paw_psi_check()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"
 

      call OrthoCheck(ispin,ne,dcpl_mb(psi1(1)))  
      return
      end



*     ***************************
*     *                         *
*     *         paw_rho_2to1    *
*     *                         *
*     ***************************
      subroutine paw_rho_2to1()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"

      call dcopy(4*nfft3d,
     >           dbl_mb(rho2(1)),1,
     >           dbl_mb(rho1(1)),1)

      return
      end

*     ***************************
*     *                         *
*     *         paw_rho_1to2    *
*     *                         *
*     ***************************
      subroutine paw_rho_1to2()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"


      call dcopy(4*nfft3d,
     >           dbl_mb(rho1(1)),1,
     >           dbl_mb(rho2(1)),1)

      return
      end

*     ***************************
*     *                         *
*     *         paw_dng_2to1    *
*     *                         *
*     ***************************
      subroutine paw_dng_2to1()
      implicit none
 
#include "bafdecls.fh"
#include "paw_psi.fh"
 
      call zcopy(npack0,
     >           dcpl_mb(dng2(1)),1,
     >           dcpl_mb(dng1(1)),1)

      return
      end




*     ***********************************
*     *                                 *
*     *         paw_psi_1toelectron             *
*     *                                 *
*     ***********************************
      subroutine paw_psi_1toelectron()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"

 
      call paw_electron_run(dcpl_mb(psi1(1)),
     >                  dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)))

      return
      end

*     ***********************************
*     *                                 *
*     *         paw_psi_1energy         *
*     *                                 *
*     ***********************************
      real*8 function paw_psi_1energy()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"

*     **** external functions ****
      real*8   paw_electron_energy
      external paw_electron_energy

      call paw_electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                   dcpl_mb(dng1(1)))
      paw_psi_1energy = paw_electron_energy(dcpl_mb(psi1(1)),
     >                               dbl_mb(rho1(1)),
     >                              dcpl_mb(dng1(1)))
      return
      end

*     ***********************************
*     *                                 *
*     *     paw_psi_1_noupdate_energy   *
*     *                                 *
*     ***********************************
      real*8 function paw_psi_1_noupdate_energy()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"

*     **** external functions ****
      real*8   paw_electron_energy
      external paw_electron_energy


      !call electron_gen_Hpaw_psi_k(dcpl_mb(psi1(1)))
      paw_psi_1_noupdate_energy = paw_electron_energy(dcpl_mb(psi1(1)),
     >                               dbl_mb(rho1(1)),
     >                              dcpl_mb(dng1(1)) )
      return
      end


*     ***********************************
*     *                                 *
*     *         paw_psi_2energy         *
*     *                                 *
*     ***********************************
      real*8 function paw_psi_2energy()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"


*     **** external functions ****
      real*8   paw_electron_energy
      external paw_electron_energy

      call paw_electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                  dcpl_mb(dng2(1)))
      paw_psi_2energy = paw_electron_energy(dcpl_mb(psi2(1)),
     >                               dbl_mb(rho2(1)),
     >                              dcpl_mb(dng2(1)))

      return
      end



*     ***********************************
*     *                                 *
*     *         paw_psi_1eorbit         *
*     *                                 *
*     ***********************************
      real*8 function paw_psi_1eorbit()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"

*     **** external functions ****
      real*8   paw_electron_eorbit
      external paw_electron_eorbit

      paw_psi_1eorbit = paw_electron_eorbit(dcpl_mb(psi1(1)))

      return
      end


*     ***********************************
*     *                                 *
*     *         paw_psi_1ke             *
*     *                                 *
*     ***********************************
      real*8 function paw_psi_1ke()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"

*     **** local variables ***
      real*8 ave,occ(1)

*     **** external functions ****
      real*8   paw_electron_eorbit
      external paw_electron_eorbit
      
      call ke_ave(ispin,ne,dcpl_mb(psi1(1)),ave,.false.,occ)
      paw_psi_1ke = ave

      return
      end


*     ***********************************
*     *                                 *
*     *         paw_psi_1vl             *
*     *                                 *
*     ***********************************
      real*8 function paw_psi_1vl()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"


*     **** external functions ****
      real*8   paw_electron_psi_vl_ave
      external paw_electron_psi_vl_ave
 
      paw_psi_1vl 
     > = paw_electron_psi_vl_ave(dcpl_mb(psi1(1)),dbl_mb(rho1(1)))

      return
      end





*     ***********************************
*     *                                 *
*     *         paw_rho_1exc            *
*     *                                 *
*     ***********************************
      real*8 function paw_rho_1exc()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"


*     **** external functions ****
      real*8   paw_electron_exc
      external paw_electron_exc

      paw_rho_1exc = paw_electron_exc(dbl_mb(rho1(1)))
      return
      end

*     ***********************************
*     *                                 *
*     *         paw_rho_1pxc            *
*     *                                 *
*     ***********************************
      real*8 function paw_rho_1pxc()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"

*     **** external functions ****
      real*8   paw_electron_pxc
      external paw_electron_pxc

      paw_rho_1pxc = paw_electron_pxc(dbl_mb(rho1(1)))
      return
      end


*     ***********************************
*     *                                 *
*     *         paw_dng_1ehartree           *
*     *                                 *
*     ***********************************
      real*8 function paw_dng_1ehartree()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"

*     **** external functions ****
      real*8   paw_electron_ehartree
      external paw_electron_ehartree

      paw_dng_1ehartree = paw_electron_ehartree(dcpl_mb(dng1(1)))
      return
      end



*     ***********************************
*     *                                 *
*     *         paw_psi_2toelectron     *
*     *                                 *
*     ***********************************
      subroutine paw_psi_2toelectron()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"


      call paw_electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                   dcpl_mb(dng2(1)))
      return
      end


*     ***********************************
*     *                                 *
*     *     paw_psi_1check_Tangent      *
*     *                                 *
*     ***********************************
*
*   This routine checks the accuracy of the tangent vector.
*   MM = Yt*S*H = Yt*S*(I-Y*Yt*S)*G = Yt*S*G - Yt*S*Y*Yt*S*G = Yt*S*G - Yt*S*G == 0

*     Updated - 5-18-2002
*
      subroutine paw_psi_1check_Tangent(H)
      implicit none
      complex*16 H(*)

#include "errquit.fh"
#include "bafdecls.fh"
#include "paw_psi.fh"

*     **** local variables ****
      logical value
      integer ms,n,indx,i,j
      integer MM(2)
      real*8 sum

      do ms=1,ispin      
         n = ne(ms)
         if (n.eq.0) go to 101  !*** ferromagnetic check ***
         value = BA_push_get(mt_dbl,n*n,'MM',MM(2),MM(1))
         if (.not. value) 
     >   call errquit('out of stack memory in paw_psi_1check_Tangent',0,
     &       MA_ERR)

         indx = (ms-1)*ne(1)*npack1

*        **** calculate MM = Yt*S*H ****
         call paw_overlap_matrix_gen(n,n,
     >                     dcpl_mb(psi1(1)+indx),
     >                     H(1+indx),
     >                     dbl_mb(MM(1)))

*        **** write out MM matrix  ****
         sum = 0.0d0
         do j=1,n
         do i=1,n
            sum = sum + dbl_mb(MM(1)+(i-1)+(j-1)*n)
         end do
         end do
         write(*,*) "paw_psi_1check_Tangent:",sum
            


         value = BA_pop_stack(MM(2))
         if (.not. value) 
     >    call errquit(
     >         'error popping stack memory in paw_psi_1check_Tangent',0,
     >        MA_ERR)

 101     continue
      end do

      return
      end



*     ***********************************
*     *                                 *
*     *     paw_psi_2check_Tangent      *
*     *                                 *
*     ***********************************
*
*   This routine checks the accuracy of the tangent vector.
*   MM = Yt*S*H = Yt*S*(I-Y*Yt*S)*G = Yt*S*G - Yt*S*Y*Yt*S*G = Yt*S*G - Yt*S*G == 0

*     Updated - 5-18-2002
*
      subroutine paw_psi_2check_Tangent(H)
      implicit none
      complex*16 H(*)

#include "bafdecls.fh"
#include "paw_psi.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer ms,n,indx,i,j
      integer MM(2)
      real*8 sum

      do ms=1,ispin
         n = ne(ms)
         if (n.eq.0) go to 101  !*** ferromagnetic check ***
         value = BA_push_get(mt_dbl,n*n,'MM',MM(2),MM(1))
         if (.not. value)
     >   call errquit('out of stack memory in paw_psi_1check_Tangent',0,
     >        MA_ERR)

         indx = (ms-1)*ne(1)*npack1

*        **** calculate MM = Yt*S*H ****
         call paw_overlap_matrix_gen(n,n,
     >                     dcpl_mb(psi2(1)+indx),
     >                     H(1+indx),
     >                     dbl_mb(MM(1)))

*        **** write out MM matrix  ****
         sum = 0.0d0
         do j=1,n
         do i=1,n
            sum = sum + dbl_mb(MM(1)+(i-1)+(j-1)*n)
         end do
         end do
         write(*,*) "paw_psi_2check_Tangent:",sum



         value = BA_pop_stack(MM(2))
         if (.not. value)
     >    call errquit(
     >         'error popping stack memory in paw_psi_2check_Tangent',0,
     &       MA_ERR)

 101     continue
      end do

      return
      end



*     ***********************************
*     *                                 *
*     *         paw_psi_1get_Tgradient  *
*     *                                 *
*     ***********************************

*     THpsi = Hpsi - Y*Y^t*S*Hpsi ! used by Grassman minimizers
*     THpsi = Hpsi - S*Y*Y^t*Hpsi ! used by Grassman minimizers
*     THpsi = S^(-1)*Hpsi - Y*Y^t*Hpsi ! used by Grassman minimizers
*    Y(t) = Y*V*Cos(Sigma*t)*Vt + U*Sin(Sigma*t)*Vt
*
      subroutine paw_psi_1get_Tgradient(THpsi,Eout)
      implicit none
      complex*16 THpsi(*)
      real*8 Eout

#include "bafdecls.fh"
#include "paw_psi.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer tmp1(2)
 
*     **** external functions ****
      real*8   paw_electron_energy
      external paw_electron_energy


      value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in paw_psi_1get_Tradient',0,
     &       MA_ERR)


      call paw_electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)))
      Eout =  paw_electron_energy(dcpl_mb(psi1(1)),
     >                        dbl_mb(rho1(1)),
     >                        dcpl_mb(dng1(1)))
      call paw_electron_gen_hml_S(dcpl_mb(psi1(1)),
     >                       dbl_mb(tmp1(1)))
      call paw_electron_get_Tgradient(dcpl_mb(psi1(1)),
     >                             dbl_mb(tmp1(1)),
     >                            THpsi)

      value = BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('paw_psi_1get_Tgradient:error popping stack',1,
     >     MA_ERR)

      return
      end

*     ***********************************
*     *                                 *
*     *         paw_psi_1get_residual   *
*     *                                 *
*     ***********************************

*
*     R = Hpsi - S*Y*Y^t*Hpsi ! used by Grassman minimizers
*
      subroutine paw_psi_1get_residual(R,Spsi,Eout)
      implicit none
      complex*16 R(*)
      complex*16 Spsi(*)
      real*8 Eout

#include "bafdecls.fh"
#include "paw_psi.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer tmp1(2)

*     **** external functions ****
      real*8   paw_electron_energy
      external paw_electron_energy


      value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
      if (.not. value)
     >   call errquit('out of stack memory in paw_psi_1gen_residual',0,
     &       MA_ERR)


      call paw_electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)))

      Eout =  paw_electron_energy(dcpl_mb(psi1(1)),
     >                        dbl_mb(rho1(1)),
     >                        dcpl_mb(dng1(1)))

      call paw_electron_gen_hml(dcpl_mb(psi1(1)),             ! tmp = <psi|Hpsi>
     >                       dbl_mb(tmp1(1)))

      call paw_ovlp_S((ne(1)+ne(2)),dcpl_mb(psi1(1)),Spsi)         ! Spsi = S*psi1

      call paw_electron_get_Tgradient(Spsi,dbl_mb(tmp1(1)),R) ! R = Hpsi - Spsi*tmp

      value = BA_pop_stack(tmp1(2))
      if (.not. value)
     > call errquit('paw_psi_1get_residual:error popping stack',1,
     >     MA_ERR)
      return
      end


*     ***********************************
*     *                                 *
*     *         paw_psi_2get_residual   *
*     *                                 *
*     ***********************************

*
*     R = Hpsi - S*Y*Y^t*Hpsi ! used by Grassman minimizers
*
      subroutine paw_psi_2get_residual(option,R,Spsi,Eout)
      implicit none
      integer option
      complex*16 R(*)
      complex*16 Spsi(*)
      real*8 Eout

#include "bafdecls.fh"
#include "paw_psi.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer tmp1(2)

*     **** external functions ****
      real*8   paw_electron_energy
      external paw_electron_energy


      value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
      if (.not. value)
     >   call errquit('out of stack memory in paw_psi_2get_residual',0,
     &       MA_ERR)


      if (option.le.1) then
        call paw_electron_run(dcpl_mb(psi2(1)),
     >                     dbl_mb(rho2(1)),
     >                     dcpl_mb(dng2(1)))
      end if
      Eout =  paw_electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)))

      call paw_electron_gen_hml(dcpl_mb(psi2(1)),             ! tmp = <psi|Hpsi>
     >                       dbl_mb(tmp1(1)))

      call paw_ovlp_S((ne(1)+ne(2)),dcpl_mb(psi2(1)),Spsi)         ! Spsi = S*psi1

      call paw_electron_get_Tgradient(Spsi,dbl_mb(tmp1(1)),R) ! R = Hpsi - Spsi*tmp

      value = BA_pop_stack(tmp1(2))
      if (.not. value)
     > call errquit('paw_psi_2get_residual:error popping stack',1,
     >     MA_ERR)
      return
      end




*     ***********************************
*     *                                 *
*     *         paw_psi_1get_Gradient   *
*     *                                 *
*     ***********************************

*     THpsi = Hpsi ! used by Projected Grassman minimizers
*
      subroutine paw_psi_1get_Gradient(THpsi,Eout)
      implicit none
      complex*16 THpsi(*)
      real*8 Eout

#include "bafdecls.fh"
#include "paw_psi.fh"


*     **** local variables ****
 
*     **** external functions ****
      real*8   paw_electron_energy
      external paw_electron_energy


      call paw_electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)))

      Eout =  paw_electron_energy(dcpl_mb(psi1(1)),
     >                        dbl_mb(rho1(1)),
     >                        dcpl_mb(dng1(1)))

      call paw_electron_get_Gradient(THpsi)

      return
      end


*     ***********************************
*     *                                 *
*     *         paw_psi_1gen_Tangent    *
*     *                                 *
*     ***********************************

*     T = T - Y*Y^t*S*T ! used by Grassman minimizers
*
      subroutine paw_psi_1gen_Tangent(T)
      implicit none
      complex*16 T(*)

#include "bafdecls.fh"
#include "paw_psi.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer tmp1(2),ii,jj


      value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >call errquit('paw_psi_1gen_Tangent:out of stack memory',0,MA_ERR)

      call paw_elecpsitang(dcpl_mb(psi1(1)),
     >                                T,
     >                                dbl_mb(tmp1(1)))
      call paw_electron_gen_Tangent(dcpl_mb(psi1(1)),
     >                          dbl_mb(tmp1(1)),
     >                          T)

      value = BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit(
     > 'error popping stack memory in paw_psi_1gen_Tangent',0,MA_ERR)

      return
      end





*     ***********************************
*     *                                 *
*     *         paw_psi_2get_Tgradient  *
*     *                                 *
*     ***********************************
      subroutine paw_psi_2get_Tgradient(option,THpsi,Eout)
      implicit none
      integer    option
      complex*16 THpsi(*)
      real*8     Eout

#include "errquit.fh"
#include "bafdecls.fh"
#include "paw_psi.fh"


*     *** local variables ****
      logical value
      integer tmp1(2)

*     **** external functions ****
      real*8   paw_electron_energy
      external paw_electron_energy


      value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in paw_psi_1get_Tradient',0,
     &       MA_ERR)

      if (option.le.1) then
        call paw_electron_run(dcpl_mb(psi2(1)),
     >                     dbl_mb(rho2(1)),
     >                    dcpl_mb(dng2(1)))
      end if

      Eout =  paw_electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)))
      call paw_electron_gen_hml_S(dcpl_mb(psi2(1)),
     >                       dbl_mb(tmp1(1)))
      call paw_electron_get_Tgradient(dcpl_mb(psi2(1)),
     >                             dbl_mb(tmp1(1)),
     >                             THpsi)
      
      value = BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit(
     > 'paw_psi_2get_Tgradient:error popping stack',1,MA_ERR)

      return
      end


*     ***********************************
*     *                                 *
*     *         paw_psi_2get_Gradient   *
*     *                                 *
*     ***********************************
      subroutine paw_psi_2get_Gradient(option,THpsi,Eout)
      implicit none
      integer    option
      complex*16 THpsi(*)
      real*8     Eout

#include "bafdecls.fh"
#include "paw_psi.fh"


*     *** local variables ****

*     **** external functions ****
      real*8   paw_electron_energy
      external paw_electron_energy


      if (option.le.1) then
        call paw_electron_run(dcpl_mb(psi2(1)),
     >                     dbl_mb(rho2(1)),
     >                    dcpl_mb(dng2(1)))
      end if

      Eout =  paw_electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)))

      call paw_electron_get_Gradient(THpsi)

      return
      end


*     ***********************************
*     *                                 *
*     *         paw_psi_2gen_Tangent    *
*     *                                 *
*     ***********************************

*     T = T - Y*Y^t*S*T ! used by Grassman minimizers
*
      subroutine paw_psi_2gen_Tangent(T)
      implicit none
      complex*16 T(*)

#include "bafdecls.fh"
#include "paw_psi.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer tmp1(2),ii,jj


      value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
      if (.not. value)
     >call errquit('paw_psi_2gen_Tangent:out of stack memory',0,MA_ERR)

      call paw_elecpsitang(dcpl_mb(psi2(1)),
     >                                T,
     >                                dbl_mb(tmp1(1)))
      call paw_electron_gen_Tangent(dcpl_mb(psi2(1)),
     >                          dbl_mb(tmp1(1)),
     >                          T)

      value = BA_pop_stack(tmp1(2))
      if (.not. value)
     > call errquit(
     > 'error popping stack memory in paw_psi_2gen_Tangent',0,MA_ERR)

      return
      end




*     ***********************************
*     *                                 *
*     *         paw_psi_1get_TSgradient *
*     *                                 *
*     ***********************************

*     THpsi = Hpsi - Y*Hpsi^t*Y ! used by Stiefel minimizers
*
      subroutine paw_psi_1get_TSgradient(THpsi,Eout)
      implicit none
      complex*16 THpsi(*)
      real*8 Eout

#include "errquit.fh"
#include "bafdecls.fh"
#include "paw_psi.fh"

*     **** local variables ****
      logical value
      integer tmp1(2)
 
*     **** external functions ****
      real*8   paw_electron_energy
      external paw_electron_energy


      value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('paw_psi_1get_TSradient:pushing stack',0, MA_ERR)


      call paw_electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)))

      Eout =  paw_electron_energy(dcpl_mb(psi1(1)),
     >                        dbl_mb(rho1(1)),
     >                        dcpl_mb(dng1(1)))


      call paw_electron_gen_hmlt(dcpl_mb(psi1(1)),
     >                       dbl_mb(tmp1(1)))
      call paw_electron_get_Tgradient(dcpl_mb(psi1(1)),
     >                             dbl_mb(tmp1(1)),
     >                            THpsi)
      

      value = BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('paw_psi_1get_TSgradient:popping stack',1, MA_ERR)

      return
      end


*     ***********************************
*     *                                 *
*     *         paw_psi_2get_TSgradient *
*     *                                 *
*     ***********************************

*     THpsi = Hpsi - Y*Hpsi^t*Y ! used by Stiefel minimizers
*
      subroutine paw_psi_2get_TSgradient(option,THpsi,Eout)
      implicit none
      integer    option
      complex*16 THpsi(*)
      real*8     Eout

#include "errquit.fh"
#include "bafdecls.fh"
#include "paw_psi.fh"

*     *** local variables ****
      logical value
      integer tmp1(2)

*     **** external functions ****
      real*8   paw_electron_energy
      external paw_electron_energy


      value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('paw_psi_2get_TSgradient:pushing stack',0, MA_ERR)

      if (option.le.1) then
        call paw_electron_run(dcpl_mb(psi2(1)),
     >                     dbl_mb(rho2(1)),
     >                    dcpl_mb(dng2(1)))
      end if

      Eout =  paw_electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)))

      call paw_electron_gen_hmlt(dcpl_mb(psi2(1)),
     >                       dbl_mb(tmp1(1)))
      call paw_electron_get_Tgradient(dcpl_mb(psi2(1)),
     >                             dbl_mb(tmp1(1)),
     >                             THpsi)
      
      value = BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('paw_psi_2get_TSgradient:popping stack',1, MA_ERR)

      return
      end




*     ***********************************
*     *                                 *
*     *         paw_psi_1get_TMgradient *
*     *                                 *
*     ***********************************
      subroutine paw_psi_1get_TMgradient(THpsi,Eout)
      implicit none
      complex*16 THpsi(*)
      real*8     Eout

#include "bafdecls.fh"
#include "paw_psi.fh"


*     **** external functions ****
      real*8   paw_electron_energy
      external paw_electron_energy


      call paw_electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)))

      Eout =  paw_electron_energy(dcpl_mb(psi1(1)),
     >                        dbl_mb(rho1(1)),
     >                        dcpl_mb(dng1(1)))

      call paw_electron_get_TMgradient(dcpl_mb(psi1(1)),
     >                            THpsi)

      return
      end



*     ***********************************
*     *                                 *
*     *         paw_psi_2get_TMgradient *
*     *                                 *
*     ***********************************
      subroutine paw_psi_2get_TMgradient(option,THpsi,Eout)
      implicit none
      integer    option
      complex*16 THpsi(*)
      real*8     Eout

#include "bafdecls.fh"
#include "paw_psi.fh"


*     **** external functions ****
      real*8   paw_electron_energy
      external paw_electron_energy

      if (option.le.1) then
        call paw_electron_run(dcpl_mb(psi2(1)),
     >                    dbl_mb(rho2(1)),
     >                    dcpl_mb(dng2(1)))
      end if

      Eout =  paw_electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)))

      call paw_electron_get_TMgradient(dcpl_mb(psi2(1)),
     >                             THpsi)
      
      return
      end

*     ***********************************
*     *                                 *
*     *    paw_psi_1ke_Precondition     *
*     *                                 *
*     ***********************************
      subroutine paw_psi_1ke_Precondition(Hpsi)
      implicit none
      complex*16 Hpsi(*)

#include "bafdecls.fh"
#include "paw_psi.fh"

*     **** local variables ****
      integer neall

      neall = ne(1)+ne(2)
      call ke_Precondition(npack1,neall,
     >                      dcpl_mb(psi1(1)),
     >                      Hpsi)
      return
      end



*     ***********************************
*     *                                 *
*     *     paw_psi_1geodesic_transport *
*     *                                 *
*     ***********************************
      subroutine paw_psi_1geodesic_transport(t,H0)
      implicit none
      real*8 t
      complex*16 H0(*)

#include "bafdecls.fh"
#include "paw_psi.fh"


      call paw_geodesic_transport(t,dcpl_mb(psi1(1)),H0)

      return
      end


*     ***********************************
*     *                                 *
*     *     paw_psi_1geodesic_Gtransport        *
*     *                                 *
*     ***********************************
      subroutine paw_psi_1geodesic_Gtransport(t,G0)
      implicit none
      real*8 t
      complex*16 G0(*)

#include "bafdecls.fh"
#include "paw_psi.fh"

      call paw_geodesic_Gtransport(t,dcpl_mb(psi1(1)),G0)

      return
      end



*     ***********************************
*     *                                 *
*     *     paw_psi_geodesic_energy     *
*     *                                 *
*     ***********************************
      real*8 function paw_psi_geodesic_energy(t)
      implicit none
      real*8 t

#include "bafdecls.fh"
#include "paw_psi.fh"


*     **** local variables ****
      real*8 e_new

*     **** external functions ****
      real*8   paw_electron_energy,paw_psi_CheckOrtho
      external paw_electron_energy,paw_psi_CheckOrtho

     
      call paw_geodesic_get(t,dcpl_mb(psi1(1)),
     >                    dcpl_mb(psi2(1)))
      call paw_electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                  dcpl_mb(dng2(1)))
      e_new =  paw_electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)))

c      write(*,*) "paw_psi_geodesic_energy:",t,e_new
      paw_psi_geodesic_energy = e_new
      return
      end

*     ***********************************
*     *                                 *
*     *    paw_psi_geodesic_denergy     *
*     *                                 *
*     ***********************************
      real*8 function paw_psi_geodesic_denergy(t)
      implicit none
      real*8 t

#include "bafdecls.fh"
#include "paw_psi.fh"

*     **** external functions ****
      real*8   paw_electron_eorbit
      external paw_electron_eorbit


      call paw_geodesic_transport(t,dcpl_mb(psi1(1)),
     >                          dcpl_mb(psi2(1)))
      paw_psi_geodesic_denergy 
     > = 2.0d0*paw_electron_eorbit(dcpl_mb(psi2(1)))

      return
      end

*     ***********************************
*     *                                 *
*     *         paw_psi_geodesic_final  *
*     *                                 *
*     ***********************************
      subroutine paw_psi_geodesic_final(t)
      implicit none
      real*8 t

#include "bafdecls.fh"
#include "paw_psi.fh"

      integer taskid,MASTER
      parameter (MASTER=0)
c     real*8 sum1,sum2
     
      call Parallel_taskid(taskid)

      call paw_geodesic_get(t,dcpl_mb(psi1(1)),
     >                    dcpl_mb(psi2(1)))          
      return
      end





*     ***********************************
*     *                                 *
*     *         paw_psito2_sd_update    *
*     *                                 *
*     ***********************************
      subroutine paw_psi1to2_sd_update(dte)
      implicit none
      real*8 dte

#include "bafdecls.fh"
#include "paw_psi.fh"
#include "frac_occ.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer nemax,ierr
      integer lmd(2),tmp_L(2)


      call paw_electron_run(dcpl_mb(psi1(1)),
     >                  dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)))

*     **** do a steepest descent step ****
      call paw_electron_sd_update(dcpl_mb(psi1(1)),
     >                        dcpl_mb(psi2(1)),
     >                        dte)

*     **** lagrange multiplier corrections ****
      nemax = ne(1)+ne(2)

*     **** allocate MA local variables ****
      value = BA_push_get(mt_dbl,(8*nemax*nemax),
     >                    'tmp_L',tmp_L(2),tmp_L(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nemax*nemax),
     >                    'lmd',lmd(2),lmd(1))

        call paw_ovlp_S(nemax,dcpl_mb(psi1(1)),dcpl_mb(psi1(1)))
        call paw_psi_lmbda(ispin,ne,nemax,npack1,
     >                 dcpl_mb(psi1(1)),dcpl_mb(psi2(1)),dte,
     >                 dbl_mb(lmd(1)),
     >                 dbl_mb(tmp_L(1)),ierr)


      value = value.and.BA_pop_stack(lmd(2))
      value = value.and.BA_pop_stack(tmp_L(2))
      if (.not. value)
     >     call errquit(
     >          'psi1to2_sd_update: stack failure', 0, MA_ERR)
      return
      end


*     ***********************************
*     *                                 *
*     *         paw_psi_1force          *
*     *                                 *
*     ***********************************
      subroutine paw_psi_1force(fion)
      implicit none
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "paw_psi.fh"

      call paw_electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)))

      call paw_electron_force(dcpl_mb(psi1(1)),dcpl_mb(dng1(1)),fion)
      return
      end



*     ***********************************
*     *                                 *
*     *         paw_psi_1Shml           *
*     *                                 *
*     ***********************************
      subroutine paw_psi_1Shml(S0,S0hml)
      implicit none
      complex*16 S0(*)
      complex*16 S0hml(*)

#include "bafdecls.fh"
#include "paw_psi.fh"

      integer ms,n,shift1,shift2

      call paw_electron_gen_hml(dcpl_mb(psi1(1)),dbl_mb(hml(1)))
      do ms=1,ispin
            n     = ne(ms)
            if (n.le.0) go to 30
            shift1 = 1 + (ms-1)*ne(1)*npack1
            shift2 =     (ms-1)*ne(1)*ne(1)
            call dgemm('N','N',2*npack1,n,n,
     >                (1.0d0),
     >                S0(shift1),            2*npack1,
     >                dbl_mb(hml(1)+shift2), n,
     >                (0.0d0),
     >                S0hml(shift1),         2*npack1)
   30       continue
      end do
      return
      end



*     ***********************************
*     *                                 *
*     *         paw_psi_1gen_hml        *
*     *                                 *
*     ***********************************
      subroutine paw_psi_1gen_hml()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"


      call paw_electron_gen_hml(dcpl_mb(psi1(1)),dbl_mb(hml(1)))

      return
      end




*     ***********************************
*     *                                 *
*     *         paw_psi_1gen_hml_g          *
*     *                                 *
*     ***********************************
      subroutine paw_psi_1gen_hml_g()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"


      call paw_electron_gen_hml_g(dcpl_mb(psi1(1)),dbl_mb(hml(1)))

      return
      end


*     ***********************************
*     *                                 *
*     *         paw_psi_2gen_hml        *
*     *                                 *
*     ***********************************
      subroutine paw_psi_2gen_hml()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"


      call paw_electron_gen_hml(dcpl_mb(psi2(1)),dbl_mb(hml(1)))

      return
      end

*     ***********************************
*     *                                 *
*     *         paw_psi_eigenvalue      *
*     *                                 *
*     ***********************************
      real*8  function paw_psi_eigenvalue(ms,i)
      implicit none
      integer ms
      integer i

#include "bafdecls.fh"
#include "paw_psi.fh"
#include "frac_occ.fh"

      real*8 sum

      sum = dbl_mb(eig(1)+(i-1)+(ms-1)*ne(1))
      paw_psi_eigenvalue = sum

      return
      end

*     ***********************************
*     *                                 *
*     *        paw_psi_virtual          *
*     *                                 *
*     ***********************************
      real*8  function paw_psi_virtual(ms,i)
      implicit none
      integer ms
      integer i

#include "bafdecls.fh"
#include "paw_psi.fh"

      paw_psi_virtual=dbl_mb(eig_excited(1)+(i-1)+(ms-1)*ne_excited(1))

      return
      end

*     ***********************************
*     *                                 *
*     *         paw_psi_hml             *
*     *                                 *
*     ***********************************
      real*8  function paw_psi_hml(ms,i,j)
      implicit none
      integer ms
      integer i,j

#include "bafdecls.fh"
#include "paw_psi.fh"

      paw_psi_hml = dbl_mb(hml(1)-1 + i 
     >                          + (j-1)*ne(ms) 
     >                          + (ms-1)*ne(1)*ne(1)) 

      return
      end


*     ***********************************
*     *                                 *
*     *         paw_psi_spin_density    *
*     *                                 *
*     ***********************************
      subroutine paw_psi_spin_density(en1,en2)
      implicit none
      real*8 en1(2),en2(2)

#include "bafdecls.fh"
#include "errquit.fh"
#include "paw_psi.fh"


*     **** local variables ****
      integer ms,nx,ny,nz,n2ft3d,tmp1(2)
      real*8  scale,sumall

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call D3dB_nfft3d(1,n2ft3d)
      n2ft3d = 2*n2ft3d
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scale = lattice_omega()/dble(nx*ny*nz)

*     **** check total number of electrons ****
      do ms =1,ispin
         call D3dB_r_dsum(1,dbl_mb(rho1(1)+(ms-1)*n2ft3d),sumall)
         en1(ms) = sumall*scale
      end do

      if (.not.BA_push_get(mt_dbl,n2ft3d,'tmp1',tmp1(2),tmp1(1)))
     >   call errquit(
     >        'paw_psi_spin_density: out of stack memory',0,MA_ERR)

      call paw_comp_charge_spin_update()

      do ms =1,ispin
         call paw_mult_dn_cmp_smooth_spin_get(ms,dbl_mb(tmp1(1)))
         call Pack_c_unpack(0,dbl_mb(tmp1(1)))
         call D3dB_cr_fft3b(1,dbl_mb(tmp1(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))
         call D3dB_r_dsum(1,dbl_mb(tmp1(1)),sumall)
         en2(ms) = sumall*scale

      end do
      if (.not.BA_pop_stack(tmp1(2)))
     >   call errquit(
     >        'paw_psi_spin_density: popping stack memory',0,MA_ERR)


      
      return
      end

*     ***********************************
*     *                                 *
*     *         paw_psi_spin2           *
*     *                                 *
*     ***********************************
      subroutine paw_psi_spin2(Sab)
      implicit none
      real*8 Sab

#include "bafdecls.fh"
#include "paw_psi.fh"

      call Calculate_psi_spin2(ispin,ne,npack1,dcpl_mb(psi1(1)),Sab)
           
      return
      end





*     ***********************************
*     *                                 *
*     *         paw_psi_1rotate2        *
*     *                                 *
*     ***********************************
      subroutine paw_psi_1rotate2()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"

*     ***** local variables *****
      integer ms,index,i,j,shift1,shift2


      do ms=1,ispin
         if (ne(ms).le.0) go to 30
         shift1 = (ms-1)*ne(1)
         shift2 = (ms-1)*ne(1)*ne(1)

         call dgemm('N','N',2*npack1,ne(ms),ne(ms),
     >              (1.0d0),
     >              dcpl_mb(psi1(1)+shift1*npack1),2*npack1,
     >              dbl_mb(hml(1)+shift2),ne(ms),
     >              (0.0d0),
     >              dcpl_mb(psi2(1)+shift1*npack1),2*npack1)

   30   continue
      end do

      return
      end

*     ***********************************
*     *                                 *
*     *         paw_psi_2rotate1        *
*     *                                 *
*     ***********************************
      subroutine paw_psi_2rotate1()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"

*     ***** local variables *****
      integer ms,index,i,j,shift1,shift2

      do ms=1,ispin
         if (ne(ms).le.0) go to 30
         shift1 = (ms-1)*ne(1)
         shift2 = (ms-1)*ne(1)*ne(1)

         call dgemm('N','N',2*npack1,ne(ms),ne(ms),
     >              (1.0d0),
     >              dcpl_mb(psi2(1)+shift1*npack1),2*npack1,
     >              dbl_mb(hml(1)+shift2),ne(ms),
     >              (0.0d0),
     >              dcpl_mb(psi1(1)+shift1*npack1),2*npack1)

   30    continue
      end do

      return
      end


*     ***********************************
*     *                                 *
*     *   paw_psi_diagonalize_hml       *
*     *                                 *
*     ***********************************
      subroutine paw_psi_diagonalize_hml()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "paw_psi.fh"


*     ***** local variables ****
      logical value
      integer ms,shift1,shift2,ierr,i,j,indx
      integer tmp1(2)


      value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     > call errquit(
     > 'out of stack memory in paw_psi_diagonalize_hml',0,MA_ERR)


*     ***** diagonalize the hamiltonian matrix *****
      call dcopy((ne(1)+ne(2)),0.0d0,0,dbl_mb(eig(1)),1)
      do ms=1,ispin
         shift1 = (ms-1)*ne(1)
         shift2 = (ms-1)*ne(1)*ne(1)
         if (ne(ms).le.0) go to 30

         call dsyev('V','U',ne(ms),
     >              dbl_mb(hml(1)+shift2),ne(ms), 
     >              dbl_mb(eig(1)+shift1),
     >              dbl_mb(tmp1(1)),2*ne(1)*ne(1),
     >              ierr)

         call eigsrt(dbl_mb(eig(1)+shift1),
     >              dbl_mb(hml(1)+shift2),
     >              ne(ms),ne(ms))

  30    continue
      end do

      
      value = BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit(
     > 'error popping stack in paw_psi_diagonalize_hml',0,MA_ERR)

      return
      end

*     *******************************************
*     *                                         *
*     *   paw_psi_diagonalize_hml_assending     *
*     *                                         *
*     *******************************************
      subroutine paw_psidiaghmlasse()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "paw_psi.fh"


*     ***** local variables ****
      logical value
      integer ms,shift1,shift2,ierr
      integer tmp1(2)


      value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit(
     >    'out of stack memory in paw_psi_diagonalize_hml_assending',0,
     >     MA_ERR)


*     ***** diagonalize the hamiltonian matrix *****
      call dcopy((ne(1)+ne(2)),0.0d0,0,dbl_mb(eig(1)),1)
      do ms=1,ispin
         shift1 = (ms-1)*ne(1)
         shift2 = (ms-1)*ne(1)*ne(1)
         if (ne(ms).le.0) go to 30

         call dsyev('V','U',ne(ms),
     >              dbl_mb(hml(1)+shift2),ne(ms), 
     >              dbl_mb(eig(1)+shift1),
     >              dbl_mb(tmp1(1)),2*ne(1)*ne(1),
     >              ierr)

   30    continue
      end do

      
      value = BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit(
     >   'error popping stack in paw_psi_diagonalize_hml_assending',0,
     >     MA_ERR)

      return
      end



*     ***************************
*     *                         *
*     *         paw_psi_error   *
*     *                         *
*     ***************************
      real*8 function paw_psi_error()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "paw_psi.fh"

*     ***** local variables ****
      logical value
      integer k,n
      real*8  error,sum,size
      integer tmp1(2)

      value = BA_push_get(mt_dcpl,(npack1),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in paw_psi_error',0, MA_ERR)


      error = 0.0d0
      size =  dble(ne(1)+ne(2))
      do n=1, (ne(1)+ne(2))
         do k=1,npack1
            dcpl_mb(tmp1(1)+k-1) = dcpl_mb(psi2(1)+k-1+(n-1)*npack1) 
     >                           - dcpl_mb(psi1(1)+k-1+(n-1)*npack1)
         end do
c         call D3dB_cc_dot(1,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)),sum)
         call Pack_cc_dot(1,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)),sum)

         error = error + sum
      end do
      error = dsqrt(error)/size
      
      value = BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit(
     > 'error popping stack memory in paw_psi_error',0,MA_ERR)


      paw_psi_error = error
      return
      end

*     ***************************
*     *                         *
*     *         paw_rho_error   *
*     *                         *
*     ***************************
      real*8 function paw_rho_error()
      implicit none

#include "errquit.fh"
#include "bafdecls.fh"
#include "paw_psi.fh"

*     ***** local variables ****
      logical value
      integer k,nx,ny,nz
      real*8  error,scale
      integer tmp1(2)

*     ***** external functions *****
      real*8   lattice_omega
      external lattice_omega

      value = BA_push_get(mt_dbl,(2*nfft3d),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in rho_error',0, MA_ERR)


      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scale = lattice_omega()

      scale = (scale)/dble(nx*ny*nz)
*     scale = (scale)/dble(nx*ny*nz)
*     scale = (scale*scale)

      do k=1,(2*nfft3d)
         dbl_mb(tmp1(1)+k-1) = (dbl_mb(rho2(1)+k-1)
     >                         -dbl_mb(rho1(1)+k-1)) 
         dbl_mb(tmp1(1)+k-1) = dbl_mb(tmp1(1)+k-1) 
     >                      + (dbl_mb(rho2(1)+k-1+(ispin-1)*(2*nfft3d))
     >                        -dbl_mb(rho1(1)+k-1+(ispin-1)*(2*nfft3d)))
      end do
      call D3dB_rr_dot(1,dbl_mb(tmp1(1)),dbl_mb(tmp1(1)),error)
      error = error*scale
*     error = dsqrt(error)

      value = BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory in rho_error',0, MA_ERR)


      paw_rho_error = error
      return
      end


*     ***************************
*     *                         *
*     *         paw_rho_dipole  *
*     *                         *
*     ***************************
*
*     Uses - Calculate_dipole (pspw/lib/psi/dipole.f)
*
      subroutine paw_rho_dipole(dipole)
      implicit none
      real*8 dipole(3)

#include "bafdecls.fh"
#include "paw_psi.fh"

      call Calculate_Dipole(ispin,ne,2*nfft3d,dbl_mb(rho1(1)),dipole)
      return
      end


*     ***************************
*     *                         *
*     *         paw_psi_ispin   *
*     *                         *
*     ***************************
      integer function paw_psi_ispin()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"

      paw_psi_ispin = ispin
      return
      end


*     ***************************
*     *                         *
*     *         paw_psi_ne              *
*     *                         *
*     ***************************
      integer function paw_psi_ne(ms)
      implicit none
      integer ms

#include "bafdecls.fh"
#include "paw_psi.fh"

      paw_psi_ne = ne(ms)
      return
      end



*     ***************************
*     *                         *
*     *    paw_psi_initialize   *
*     *                         *
*     ***************************

      logical function paw_psi_initialize()
      implicit none 


#include "bafdecls.fh"
#include "errquit.fh"
#include "btdb.fh"
#include "paw_psi.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical value,psi_nogrid
      integer nemax,ms,idum
      real*8 sum1,sum2,dum(1)
      integer hversion,hnfft(3),hispin,hne(2)
      real*8 hunita(3,3)
      integer rtdb,ind

      integer  control_rtdb,control_ngrid
      external control_rtdb,control_ngrid
      character*50 filename
      character*50 control_input_psi
      external     control_input_psi
      logical  wvfnc_expander,Dneall_m_allocate
      external wvfnc_expander,Dneall_m_allocate
      real*8   paw_psi_CheckOrtho
      external paw_psi_CheckOrtho


      call Parallel_taskid(taskid)

*     ***** get ispin, and ne, and nfft3d ****
      call psi_get_ne(ispin,ne)
      call Dneall_neq(neq)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)
      call Pack_npack(0,npack0)
      nemax = ne(1)+ne(2)

*     **** allocate memory ****
      value = BA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'psi1',psi1(2),psi1(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,4*nfft3d,
     >                     'rho1',rho1(2),rho1(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,4*nfft3d,
     >                     'rho2',rho2(2),rho2(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack0,
     >                     'dng1',dng1(2),dng1(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack0,
     >                     'dng2',dng2(2),dng2(1))
c      value = value.and.
c     >        BA_alloc_get(mt_dbl,(2*nemax*nemax),'hml',hml(2),hml(1))
      value = value.and.Dneall_m_allocate(0,hml)

      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*nemax),'eig',eig(2),eig(1))

      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

*     *****  read initial wavefunctions into psi1  ****
      rtdb = control_rtdb()
      if (.not.btdb_get(rtdb,'nwpw:psi_nogrid',
     >                  mt_log,1,psi_nogrid))
     >   psi_nogrid = .true.

      if (psi_nogrid) then
        
        call psi_get_header(hversion,hnfft,hunita,hispin,hne)
      
        if ( (hnfft(1).ne.control_ngrid(1)) .or.
     >       (hnfft(2).ne.control_ngrid(2)) .or.
     >       (hnfft(3).ne.control_ngrid(3)) ) then

        hnfft(1) = control_ngrid(1)
        hnfft(2) = control_ngrid(2)
        hnfft(3) = control_ngrid(3)

        call ga_sync()
        value = btdb_parallel(.false.)
        call ga_sync()
        if (taskid.eq.MASTER) then
       
          filename =  control_input_psi()

          ind = index(filename,' ') - 1
          if (.not. btdb_cput(rtdb,'xpndr:old_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: btdb_cput failed', 0, RTDB_ERR)

          if (.not. btdb_cput(rtdb,'xpndr:new_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: btdb_cput failed', 0, RTDB_ERR)

          if (.not. btdb_put(rtdb,'xpndr:ngrid',mt_int,3,hnfft))
     >     call errquit(
     >     'wvfnc_expander_input: btdb_put failed', 0, RTDB_ERR)

          write(*,*)
          write(*,*) "Grid is being converted:"
          write(*,*) "------------------------"
          write(*,*)
          write(*,*) "To turn off automatic grid conversion:"
          write(*,*)
          write(*,*) "set nwpw:psi_nogrid .false."
          write(*,*)
          value = wvfnc_expander(rtdb)

        end if
        call ga_sync()
        value = btdb_parallel(.true.)

      end if

      end if

      call psi_read(ispin,ne,dcpl_mb(psi1(1)),idum,dum)

c      call psi_history_read(ispin,ne,
c     >                      dcpl_mb(psi1(1)),
c     >                      dcpl_mb(psi2(1)))
     
      call paw_ovlp_coeff_set(dcpl_mb(psi1(1)))
      call paw_ovlp_weights_set()

      !**** Ortho Check ****
      do ms=1,ispin
        sum1=paw_psi_CheckOrtho(npack1,ne(ms),
     >                   dcpl_mb(psi1(1)+(ms-1)*ne(1)*npack1))

        if (sum1.gt.1.0d-10) then
          call paw_psi_MakeOrtho(npack1,ne(ms),
     >                   dcpl_mb(psi1(1)+(ms-1)*ne(1)*npack1))
          sum2=paw_psi_CheckOrtho(npack1,ne(ms),
     >                   dcpl_mb(psi1(1)+(ms-1)*ne(1)*npack1))
          if (taskid.eq.MASTER) then
            if (ms.eq.1) then
              write(*,24) sum1,sum2
            end if
            if (ms.eq.2) then
              write(*,25) sum1,sum2
            end if
  24  format('Gram-Schmidt performed on up spin  : (old error=',
     >        E10.3,' new error=',E10.3,')' )
  25  format('Gram-Schmidt performed on down spin: (old error=',
     >        E10.3,' new error=',E10.3,')' )

          end if
        end if

      end do


      paw_psi_initialize = value
      return
      end

*     ***************************
*     *                         *
*     *   paw_psi_tmp_write     *
*     *                         *
*     ***************************
      subroutine paw_psi_tmp_write()
      implicit none

#include "bafdecls.fh"
#include "paw_psi.fh"

      real*8 dum(1)

*     ***** write psi1 wavefunctions ****
      call psi_write(ispin,ne,dcpl_mb(psi1(1)),-1,dum)

      return
      end



*     ***************************
*     *                         *
*     *    paw_psi_finalize     *
*     *                         *
*     ***************************

      logical function paw_psi_finalize(wpsi)
      implicit none 
      logical wpsi

#include "bafdecls.fh"
#include "errquit.fh"
#include "paw_psi.fh"


*     **** local variables ****
      logical value
      real*8  dum

*     ***** write psi1 wavefunctions ****
      if (wpsi) then
       call psi_write(ispin,ne,dcpl_mb(psi1(1)),-1,dum)
c      call psi_history_write(ispin,ne,dcpl_mb(psi1(1)))
      end if
    
      value = BA_free_heap(eig(2))
      value = value.and.BA_free_heap(hml(2))
      value = value.and.BA_free_heap(dng2(2))
      value = value.and.BA_free_heap(dng1(2))
      value = value.and.BA_free_heap(rho2(2))
      value = value.and.BA_free_heap(rho1(2))
      value = value.and.BA_free_heap(psi2(2))
      value = value.and.BA_free_heap(psi1(2))
      if (.not. value) 
     >  call errquit('paw_psi_finalize: error freeing heap',0, MA_ERR)

      paw_psi_finalize = value
      return
      end


*     ***************************
*     *                         *
*     *    paw_psi_ne_excited   *
*     *                         *
*     ***************************
      integer function paw_psi_ne_excited(ms)
      implicit none
      integer ms

#include "bafdecls.fh"
#include "paw_psi.fh"

      paw_psi_ne_excited = ne_excited(ms)
      return
      end

*     ***************************
*     *                         *
*     *    paw_epsi_initialize  *
*     *                         *
*     ***************************

      logical function paw_epsi_initialize()
      implicit none 
#include "errquit.fh"

#include "bafdecls.fh"
#include "btdb.fh"
#include "paw_psi.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical value,psi_nogrid
      integer nemax,ispin0
      real*8 sum1,sum2

      integer hversion,hnfft(3),hispin,hne(2)
      real*8 hunita(3,3)
      integer rtdb,ind
      integer  control_rtdb,control_ngrid
      external control_rtdb,control_ngrid
      character*50 filename
      character*50 control_input_epsi
      external     control_input_epsi
      logical  wvfnc_expander
      external wvfnc_expander


*     ***** get ispin, and ne, and nfft3d ****
      call psi_get_ne_excited(ispin0,ne_excited)
      nemax = ne_excited(1)+ne_excited(2)

*     **** allocate memory ****
      value = BA_alloc_get(mt_dcpl,npack1*(nemax),
     >         'psi2_excited',psi2_excited(2),psi2_excited(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(nemax),
     >         'psi1_excited',psi1_excited(2),psi1_excited(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*nemax),'eig_excited',
     >                     eig_excited(2),eig_excited(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      call dcopy(2*npack1*nemax,0.0d0,0,dcpl_mb(psi2_excited(1)),1)
      call dcopy(2*npack1*nemax,0.0d0,0,dcpl_mb(psi1_excited(1)),1)
      call dcopy(2*nemax,0.0d0,0,dbl_mb(eig_excited(1)),1)


*     *****  read initial wavefunctions into psi1  ****
      rtdb = control_rtdb()
      if (.not.btdb_get(rtdb,'nwpw:psi_nogrid',
     >                  mt_log,1,psi_nogrid))
     >   psi_nogrid = .true.

      if (psi_nogrid) then

        
        filename =  control_input_epsi()
        call psi_get_header_filename(filename,
     >                      hversion,hnfft,hunita,hispin,hne)

        if ( (hnfft(1).ne.control_ngrid(1)) .or.
     >       (hnfft(2).ne.control_ngrid(2)) .or.
     >       (hnfft(3).ne.control_ngrid(3)) ) then

        hnfft(1) = control_ngrid(1)
        hnfft(2) = control_ngrid(2)
        hnfft(3) = control_ngrid(3)
        call Parallel_taskid(taskid)
        value = btdb_parallel(.false.)
        if (taskid.eq.MASTER) then


          ind = index(filename,' ') - 1
          if (.not. btdb_cput(rtdb,'xpndr:old_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: btdb_cput failed', 0, RTDB_ERR)

          if (.not. btdb_cput(rtdb,'xpndr:new_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: btdb_cput failed', 0, RTDB_ERR)

          if (.not. btdb_put(rtdb,'xpndr:ngrid',mt_int,3,hnfft))
     >     call errquit(
     >     'wvfnc_expander_input: btdb_put failed', 0, RTDB_ERR)

          write(*,*)
          write(*,*) "Grid is being converted:"
          write(*,*) "------------------------"
          write(*,*)
          write(*,*) "To turn off automatic grid conversion:"
          write(*,*)
          write(*,*) "set nwpw:psi_nogrid .false."
          write(*,*)
          value = wvfnc_expander(rtdb)

        end if
        value = btdb_parallel(.true.)

      end if

      end if

*     *****  read initial wavefunctions into psi1  ****
      call epsi_read(ispin0,ne_excited,dcpl_mb(psi1_excited(1)))





      paw_epsi_initialize = value
      return
      end



*     ***************************
*     *                         *
*     *    paw_epsi_finalize    *
*     *                         *
*     ***************************

      logical function paw_epsi_finalize(writepsi)
      implicit none 
#include "errquit.fh"
      logical writepsi

#include "bafdecls.fh"
#include "paw_psi.fh"

*     **** local variables ****
      logical value

*     ***** write psi1 wavefunctions ****
      if (writepsi)
     >  call epsi_write(ispin,ne_excited,dcpl_mb(psi1_excited(1)))
    
      value = BA_free_heap(eig_excited(2))
      value = value.and.BA_free_heap(psi2_excited(2))
      value = value.and.BA_free_heap(psi1_excited(2))
      if (.not. value) 
     >  call errquit('epaw_psi_finalize: error freeing heap',0, MA_ERR)

      paw_epsi_finalize = value
      return
      end



