* $Id: cpsi_KS.F,v 1.2 2004-11-13 02:22:05 bylaska Exp $
*
*
* $Log: not supported by cvs2svn $
* Revision 1.1  2003/12/02 19:16:47  bylaska
* HGH pseudpotential added.
* TM, Hamman, HGH, pspw_default, and paw_default pseudopotential libraries have been added.
* KS minimizer updates.
* ...EJB
*



************************ KS orbital Part ************************

*     ***********************************
*     *                                 *
*     *      cpsi_KS_update             *
*     *                                 *
*     ***********************************

*    This routine (approximately) diagonalizes the KS matrix.
*
      subroutine cpsi_KS_update(psi_number,precondition,maxerror)
      implicit none
      integer psi_number
      logical precondition
      real*8 maxerror

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical done
      integer nb,i,j,neall,maxit_orb,maxit_orbs,shiftb
      real*8 error,error_out,sum

*     **** external functions ****
      integer  control_ks_maxit_orb,control_ks_maxit_orbs
      external control_ks_maxit_orb,control_ks_maxit_orbs

      neall = ne(1)+ne(2)
      maxit_orb  = control_ks_maxit_orb()   !*** should be read from rtdb ***
      maxit_orbs = control_ks_maxit_orbs()  !*** should be read from rtdb ***
      j = 0
 2    j = j+1
        error = 0.0d0
        do nb=1,nbrillioun
          shiftb = (nb-1)*npack1*neall
          do i=1,neall


           !*** orthogonalize to lower orbitals  ****
           call cpsi_project_out_f_orb1(
     >             nb,i,
     >             dcpl_mb(psi1(1)+(i-1)*npack1 + shiftb))

           !*** normalize ****
           call Cram_cc_dot(nb,
     >              dcpl_mb(psi1(1) +(i-1)*npack1 + shiftb),
     >              dcpl_mb(psi1(1) +(i-1)*npack1 + shiftb),
     >              sum)
           sum = 1.0d0/dsqrt(sum)
           call Cram_c_SMul(nb,sum,
     >              dcpl_mb(psi1(1) +(i-1)*npack1 + shiftb),
     >              dcpl_mb(psi1(1) +(i-1)*npack1 + shiftb))



            call cpsi_KS_update_orb(psi_number,precondition,maxit_orb,
     >                           maxerror,
     >                           0.1d0,nb,i,error_out)

            error = error+error_out
          end do
        end do
        error = error/dble(neall)

        done = ((j.gt.maxit_orbs).or.(error.lt.maxerror))
      if (.not.done) go to 2

      return
      end


*     *******************************************
*     *                                         *
*     *          cpsi_project_out_f_orb1        *
*     *                                         *
*     *******************************************
*
*    This routine projects out non-orthogonal components of Horb.
* This routine is needed for a KS minimizer.
*
      subroutine cpsi_project_out_f_orb1(nb,ii,Horb)
      integer nb,ii
      complex*16 Horb(*)

#include "mafdecls.fh"
#include "cpsi_common.fh"

      integer ms,i,jj,kk,shifte
      complex*16  sum

*     **** spin up orbital ****
      if (ii.le.ne(1)) then

         shift  = (nb-1)*(ne(1)+ne(2))*npack1
         shifte = 0
         ms     = 1
         kk     = ii
*     **** spin down orbital ****
      else
         shifte = ne(1)*npack1+(nb-1)*(ne(1)+ne(2))*npack1
         ms     = 2
         kk     = ii-ne(1)
      end if


      !**** project out orbitals ****
      do jj=1,(kk-1)
        call Cram_cc_zdot(nb,
     >            dcpl_mb(psi1(1) +(jj-1)*npack1+shifte),
     >            Horb,
     >            sum)

        call Cram_cc_zaxpy(nb,(-sum),
     >            dcpl_mb(psi1(1) +(jj-1)*npack1+shifte),
     >            Horb)
      end do


      return
      end




*     *******************************************
*     *                                         *
*     *          cpsi_project_out_orb           *
*     *                                         *
*     *******************************************
*
*    This routine projects out non-orthogonal components of Horb.
* This routine is needed for a KS minimizer.
*
      subroutine cpsi_project_out_orb(psi_number,nb,i,Horb)
      implicit none
      integer psi_number
      integer nb,i
      complex*16 Horb(*)

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "errquit.fh"

*     **** local variables ****
      logical ok
      integer ii,n,psi_ptr
      integer x(2)


*     **** allocate stack memory ****
      ok = MA_push_get(mt_dcpl,ne(1),'x',x(2),x(1))
      if (.not.ok)
     > call errquit('cpsi_project_out_orb: out of stack memory',0,
     &       MA_ERR)

      if (psi_number.eq.1) then
         psi_ptr=psi1(1)+(nb-1)*(ne(1)+ne(2))*npack1
      else
         psi_ptr=psi2(1)+(nb-1)*(ne(1)+ne(2))*npack1
      end if


*     **** spin up orbital ****
      if (i.le.ne(1)) then

        ii = i
        call Cram_cc_nzdot(nb,ii,
     >            dcpl_mb(psi_ptr),
     >            Horb,
     >            dcpl_mb(x(1)))
        do n=1,(ii)
           call Cram_cc_zaxpy(nb,
     >               (-dcpl_mb(x(1)+n-1)),
     >               dcpl_mb(psi_ptr+(n-1)*npack1),
     >               Horb)
        end do



*     **** spin down orbital ****
      else

        psi_ptr = psi_ptr + ne(1)*npack1
        ii = i - ne(1)
        call Cram_cc_nzdot(nb,ii,
     >            dcpl_mb(psi_ptr),
     >            Horb,
     >            dcpl_mb(x(1)))
        do n=1,(ii)
           call Cram_cc_zaxpy(nb,
     >               (-dcpl_mb(x(1)+n-1)),
     >               dcpl_mb(psi_ptr+(n-1)*npack1),
     >               Horb)
        end do


      end if

*     **** release stack memory ****
      ok = MA_pop_stack(x(2))
      if (.not. ok)
     > call errquit('cpsi_project_out_orb: poping stack memory',0,
     &       MA_ERR)

      return
      end








*     ***********************************
*     *                                 *
*     *      cpsi_KS_update_orb         *
*     *                                 *
*     ***********************************

*    This routine performs a KS update on orbital i
*
      subroutine cpsi_KS_update_orb(psi_number,
     >                             precondition,maxiteration,
     >                             maxerror,perror,nb,i,
     >                             error_out)
      implicit none
      integer psi_number
      logical precondition
      integer maxiteration
      real*8  maxerror,perror
      integer nb,i
      real*8 error_out

#include "mafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,done,oneloop
      integer it
      real*8 e0,eold,percent_error,error0,de0,lmbda_r0,lmbda_r1
      real*8 theta,sigma
      integer r1(2),t0(2),t(2),g(2)
      integer psi_ptr

      if (psi_number.eq.1) then
         psi_ptr=psi1(1) + (nb-1)*npack1*(ne(1)+ne(2))
      else
         psi_ptr=psi2(1) + (nb-1)*npack1*(ne(1)+ne(2))
      end if

      call Parallel_taskid(taskid)

      value = MA_push_get(mt_dcpl,npack1,'t0',t0(2),t0(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack1,'g',g(2),g(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack1,'t',t(2),t(1))
      if (.not. value) call errquit(
     >     'cpsi_KS_update_orb: out of stack memory',0, MA_ERR)

      done = .false.
      error0 = 0.0d0
      e0 = 0.0d0
      theta = -3.14159d0/600.0d0
      it = 0
 2    continue

         it = it + 1
         eold = e0

*        *** calculate residual (steepest descent) direction for a single band ***
         call cpsi_get_gradient_orb(psi_number,nb,i,dcpl_mb(g(1)))
         call Cram_cc_dot(nb,dcpl_mb(psi_ptr+(i-1)*npack1),
     >                   dcpl_mb(g(1)),
     >                    e0)
         e0 = -e0


         done = ((it.gt.maxiteration)
     >           .or.
     >           (dabs(e0-eold).lt.maxerror))

         if (done) go to 4

c*        **** preconditioning ****
c         if (precondition) then
c           call ke_Precondition(npack1,1,
c     >                     dcpl_mb(g(1)),
c     >                     dcpl_mb(g(1)))
c
c         end if

         call Cram_c_Copy(nb,dcpl_mb(g(1)),dcpl_mb(r1(1)))
         call cpsi_project_out_orb(psi_number,nb,i,dcpl_mb(r1(1)))



*        *** determine conjuagate direction ***
         call Cram_cc_dot(nb,dcpl_mb(r1(1)),
     >                   dcpl_mb(r1(1)),
     >                   lmbda_r1)
         call Cram_c_Copy(nb,dcpl_mb(r1(1)),dcpl_mb(t(1)))

         if (it.gt.1) then
         call Cram_cc_daxpy(nb,(lmbda_r1/lmbda_r0),
     >                   dcpl_mb(t0(1)),
     >                   dcpl_mb(t(1)))
         end if
         lmbda_r0 = lmbda_r1
         oneloop = .true.
 3       call Cram_c_Copy(nb,dcpl_mb(t(1)),dcpl_mb(t0(1)))




c!*        **** project out psi components from t ****
c!        call psi_project_out_orb(psi_number,i,dcpl_mb(t(1)))
c!        call Pack_cc_dot(1,dcpl_mb(psi_ptr+(i-1)*npack1),
c!    >                   dcpl_mb(t(1)),
c!    >                    de0)
c!        de0 = -de0
c!        call Pack_cc_daxpy(1,(de0),
c!    >                 dcpl_mb(psi_ptr+(i-1)*npack1),
c!    >                 dcpl_mb(t(1)))


*        *** normalize search direction, t ****
         call Cram_cc_dot(nb,dcpl_mb(t(1)),
     >                   dcpl_mb(t(1)),
     >                   sigma)
         sigma = dsqrt(sigma)
         de0 = 1.0d0/sigma
         call Cram_c_SMul(nb,de0,dcpl_mb(t(1)),dcpl_mb(t(1)))



*        **** compute de0 = <t|g> ****
         call Cram_cc_dot(nb,dcpl_mb(t(1)),
     >                   dcpl_mb(g(1)),
     >                   de0)

*        *** bad direction ***
         if ((de0.lt.0.0d0).and.oneloop) then
           call Cram_c_Copy(nb,dcpl_mb(g(1)),dcpl_mb(t(1)))
           oneloop = .false.
           go to 3
         end if

         de0 = -2.0d0*de0
         call cpsi_linesearch_update2(psi_number,nb,i,
     >                              theta,e0,de0,
     >                              dcpl_mb(t(1)),
     >                              sigma,
     >                              dcpl_mb(t0(1)))

      go to 2


*     **** release stack memory ****
 4    value =           MA_pop_stack(t(2))
      value = value.and.MA_pop_stack(g(2))
      value = value.and.MA_pop_stack(r1(2))
      value = value.and.MA_pop_stack(t0(2))
      if (.not. value) call errquit(
     >     'psi_KS_update_orb: popping stack memory',1, MA_ERR)

c      write(*,*) "iterations=",it," eig=",e0," error=",error_out,
c     >           theta
      error_out = dabs(e0-eold)
      return
      end



*     ************************************
*     *                                  *
*     *      cpsi_get_gradient_orb       *
*     *                                  *
*     ************************************

*    This routine returns the Hpsi(i).
* This routine is needed for a KS minimizer.
*
      subroutine cpsi_get_gradient_orb(psi_number,nb,i,Horb)
      implicit none
      integer psi_number
      integer nb,i
      complex*16 Horb(*)

#include "mafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer psi_ptr

      if (psi_number.eq.1) then
         psi_ptr=psi1(1)
      else
         psi_ptr=psi2(1)
      end if

      call c_electron_run_orb(nb,i,dcpl_mb(psi_ptr))
      call c_electron_get_gradient_orb(nb,i,Horb)

      return
      end



*     ***********************************
*     *                                 *
*     *      psi_linesearch_update2     *
*     *                                 *
*     ***********************************

*    This routine performs a linesearch on orbital i, in the direction t.
* This routine is needed for a KS minimizer.
*  e0 = <orb|g>
*  de0 = 2*<t|g>
*
      subroutine cpsi_linesearch_update2(psi_number,nb,i,theta,e0,de0,t,
     >                                  sigma,tau_t)
      implicit none
      integer psi_number
      integer nb,i
      real*8  theta
      real*8  e0,de0
      complex*16 t(*)     !search direction

      real*8     sigma
      complex*16 tau_t(*) !parallel transported search direction

#include "mafdecls.fh"
#include "cpsi_common.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer orb(2),g(2),psi_ptr
      real*8 x,y,pi,e1,de1
      real*8 theta2,e2,de2

      if (psi_number.eq.1) then
         psi_ptr=psi1(1)+(i-1)*npack1+(nb-1)*(ne(1)+ne(2))*npack1
      else
         psi_ptr=psi2(1)+(i-1)*npack1+(nb-1)*(ne(1)+ne(2))*npack1
      end if

      pi = 4.0d0*datan(1.0d0)

*     **** allocate stack memory ****
      value = MA_push_get(mt_dcpl,npack1,'orb',
     >                       orb(2),orb(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack1,'g',
     >                       g(2),g(1))
      if (.not. value) call errquit(
     >     'psi_linesearch_update: out of stack memory',0, MA_ERR)


      call Cram_c_Copy(nb,dcpl_mb(psi_ptr),
     >                    dcpl_mb(orb(1)))

*     **** orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
      !theta = pi/300.0d0
      x = cos(theta)
      y = sin(theta)
      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)),
     >                  dcpl_mb(psi_ptr))
      call Cram_cc_daxpy(nb,y,
     >                   t,
     >                   dcpl_mb(psi_ptr))

*     *** determine theta ***
      call cpsi_get_gradient_orb(psi_number,nb,i,dcpl_mb(g(1)))

      call Cram_cc_dot(nb,dcpl_mb(psi_ptr),
     >                   dcpl_mb(g(1)),
     >                   e1)
      e1 = -e1
      x = (e0 - e1 + 0.5d0*de0*sin(2*theta))
     >    /(1.0d0-cos(2*theta))
      theta = 0.5d0*datan(0.5d0*de0/x)

      x = cos(theta)
      y = sin(theta)

*     **** tau_t = (-orb*sin(theta) + t*cos(theta))*sigma ****
      call Cram_c_SMul(nb,(-y),
     >                  dcpl_mb(orb(1)),
     >                  tau_t)
      call Cram_cc_daxpy(nb,x,
     >                   t,
     >                   tau_t)
      call Cram_c_SMul(nb,sigma,
     >                  tau_t,
     >                  tau_t)

*     **** orb2 = orb*cos(theta) + t*sin(theta) ****
      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)),
     >                  dcpl_mb(psi_ptr))
      call Cram_cc_daxpy(nb,y,
     >                   t,
     >                   dcpl_mb(psi_ptr))


*     **** release stack memory ****
      value =           MA_pop_stack(g(2))
      value = value.and.MA_pop_stack(orb(2))
      if (.not. value) call errquit(
     >     'psi_linesearch_update: popping stack memory',1, MA_ERR)

      return
      end
