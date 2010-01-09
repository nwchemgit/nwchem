*
* $Id: psi.F,v 1.99 2009-04-02 00:01:34 bylaska Exp $
*


* $Log: not supported by cvs2svn $
* Revision 1.98  2009/03/19 20:42:17  bylaska
* ...EJB
*
* Revision 1.97  2009/02/24 21:30:17  bert
* In psi_finalize, occ1 and occ2 were deallocated in wrong order.
*
* Revision 1.96  2009/02/07 03:50:56  bylaska
* Bassi Vectorization Fix...EJB
*
* Revision 1.95  2008/12/18 21:15:51  bylaska
* ...updates for calculating spin contaminatio....EJB
*
* Revision 1.94  2008/11/17 17:25:45  bylaska
* fractional occupation updates....EJB
*
* Revision 1.93  2008/10/22 23:56:43  bylaska
* added NWCHEM_NWPW_LIBRARY to nwchemrc. fixed bug in paw...EJB
*
* Revision 1.92  2008/09/30 19:53:35  bylaska
* Added Baden's exchange algorithm...EJB
*
* Revision 1.91  2008/09/17 00:55:36  bylaska
* ...EJB
*
* Revision 1.90  2008/09/15 20:25:33  bylaska
* ...fractional bug fix..EJB
*
* Revision 1.89  2008/09/11 21:26:51  bylaska
* ...EJB
*
* Revision 1.88  2008/06/21 19:37:16  bylaska
*
* initalization error fixed with psi_Tgradient...EJB
*
* Revision 1.87  2008/06/02 15:20:04  bylaska
* ..io fixes...EJB
*
* Revision 1.86  2008/05/13 02:10:36  bylaska
* ...EJB
*
* Revision 1.85  2008/04/21 19:34:27  bylaska
* queue fft added to cpsi_H, bug fixes in DMatrix_dgemm1_rot (MPI_Allgather routine replaced with a MPI_AllReduce for stability), np_orbital keyword replaced with np_dimensions keyword (needed for Parallel3d routines)
* ...EJB
*
* Revision 1.84  2007/11/17 00:25:26  d3p708
* ....PNJ
*
* Revision 1.83  2007/11/17 00:19:09  d3p708
* pjn
*    bugfix
*
* Revision 1.82  2007/11/16 22:32:53  d3p708
* pjn
*    stuff for berry phase dipole
*
* Revision 1.81  2007/10/01 23:02:46  bylaska
* PAW changes..EJB
*
* Revision 1.80  2007/09/29 00:33:37  bylaska
* ...EJB
*
* Revision 1.79  2007/09/24 16:58:14  bylaska
* ...preliminary PAW modifications...
*   - basis file format changed
*   - .vpp formatting routines added to pspw
*
* - zdotc routines currently modified to tzdotc.
* ...EJB
*
* Revision 1.78  2007/09/13 20:38:36  bylaska
* occupation template added to band...EJB
*
* Revision 1.77  2007/03/27 02:02:49  bylaska
* more qmmm_updates....EJB
*
* Revision 1.76  2007/03/22 20:46:22  bylaska
* New implementation of QM/MM.
* ....EJB
*
* Revision 1.75  2007/02/23 01:24:32  bylaska
* ...EJB
*
* Revision 1.74  2007/02/10 03:56:54  bylaska
* ...bug fix...
* ..EJB
*
* Revision 1.73  2007/02/10 03:40:18  bylaska
* replaced calls to Grsm_g_MakeOrtho with Dneall_f_ortho
* ...EJB
*
* Revision 1.72  2007/01/02 18:36:52  bylaska
* HGH pseudopotentials added to band.
* ...EJB
*
* Revision 1.71  2006/10/13 01:43:58  bylaska
* tcgmsg code for 2d grid distribution.  Also cleaned up Dmatrix_ calls so
* that they are Dneall_ calls instead.
* ....EJB
*
* Revision 1.70  2006/10/07 00:10:07  bylaska
* Initial implementation of 2d processor grid parallelization in pspw.  Currently works with:
*
* task pspw steepest_descent
* task pspw energy            (only minimizer 1, minimizer 2?, other minimizers not yet implemented)
*
* Currently only works with USE_MPIF option, straight tcgmsg only partially implemented.  Car-Parrinello, HFX, SIC, and various analysis codes are also not yet ported.
*
*
* The number of processors along the orbital dimension, np_orbital, is entered as follows, e.g.:
*
* nwpw
*    np_orbital  2
* end
*
* The number of processors along the grid dimension, np_grid, is currently defined using np_orbital as
*
* np_grid = np/np_orbital
*
* where np is the total number of processors.
*
* ...EJB
*
* Revision 1.69  2006/09/20 19:18:49  bylaska
* Adding Dmatrix
* ...EJB
*
* Revision 1.68  2006/08/13 01:03:28  bylaska
* Checking in code not include in 5.0 release.
* A chain algorithm was added to Nose-Hoover thermostats.
* Preliminary implementation of a processor group decomposition added to pspw, i.e. parallel decomposition is over fft grid and electrons.
* ...EJB
*
* Revision 1.67  2006/01/26 18:29:36  bylaska
* bug fix for gga checking...EJB
*
* Revision 1.66  2006/01/06 22:52:28  bylaska
* parallel io bug fix...EJB
*
* Revision 1.65  2006/01/06 21:48:43  bylaska
* io changes for inversion symmetry....EJB
*
* Revision 1.64  2005/12/29 03:06:09  marat
* qmmm interface stuff
*
* Revision 1.63  2005/12/22 01:35:07  bylaska
* revPBE added and gga logic restructured....EJB
*
* Revision 1.62  2005/10/05 21:21:30  bylaska
* psi_iptr_write added...EJB
*
* Revision 1.61  2005/05/24 17:36:27  bylaska
* Stresses added to SIC and HFX.
* ....BLYP functional updates
* ....ECCE hacks
* ...EJB
*
* Revision 1.60  2005/02/09 02:39:10  bylaska
* .............................EJB
*
* Revision 1.59  2005/01/17 20:51:33  edo
* fixed a  couple of FPEs
*
* Revision 1.58.2.1  2005/01/17 20:51:06  edo
* fixed a couple of FPEs
*
* Revision 1.58  2004/12/21 16:58:35  bylaska
* various io fixes for dos...EJB
*
* Revision 1.57  2004/12/06 20:03:25  bylaska
* RMM-DIIS diagonalizer added to pspw.
* nwpw_list updated to handle multiple lists.
* ....EJB
*
* Revision 1.56  2004/11/29 16:05:21  bylaska
* Finite difference stresses added to PSPW, BAND, and PAW modules.
*    - This is currently the default for BAND and PAW
* Fixed the analytic unrestricted gga stress term in PSPW.
* Fixed unrestricted optimization for minimizers 1 and 2 in the PSPW and PAW modules.
* Partial implementation of analytic stress in BAND module.
*    - kinetic, ewald, and coulomb stresses have been implemented
*
* ....EJB
*
* Revision 1.55  2004/11/08 01:32:49  bylaska
*
* Unrestricted and closed-shell restricted Hartree-Fock has been implemented into pspw.
*    - works with minimizers 1,2,3,4,6,and 7.
*    - band-by-band minimizer (minimizer 5) not yet implemented
*    - free-space coulomb and screened-coulomb kernels implemented.
*    - the free-space coulomb kernel has been tested.
*    - the screened-coulomb kernel needs to be debugged/tested?
*    - restricted open-shell HF has not yet been implemented
*
* ....EJB
*
* Revision 1.54  2004/03/08 22:51:27  bylaska
* Fractional occupation working in pspw with minimizer 4, steepest descent, and Car-Parrinello.
*
* Bug fix in velocity initialization in liquid and solid-state Car-Parrinello simulations...incell3 instead of incell2 was used in newton step.
*
* Added restart capabilites to thermostat masses...Qe and Qr and eke0 are now propagated to
* restart Car-Parrinello simulations.
*
* SIC input modifications.
*
* Wannier orbital output modifications.
*
* ....EJB
*
* Revision 1.53  2004/03/02 00:10:22  bylaska
* ....EJB
*
* Revision 1.52  2004/03/01 05:14:33  bylaska
* Mulliken and DOS fixes.
* Added Mulliken projections based on atomic orbitals
* Added projected density of states (based on Mulliken projections)
* ...EJB
*
* Revision 1.51  2004/02/06 01:57:22  bylaska
*
* Option added to write out temporary psi for Kiril.
* Tempory psi written if
* set nwpw:psi_tmp .true.
*
*
* ....EJB
*
* Revision 1.50  2004/01/28 00:08:47  bylaska
* Bug fixes...EJB
*
* Revision 1.49  2003/12/13 21:07:41  bylaska
* Kohn-Sham scf minimizer added to pspw. Mixing includes
* 	- simple mixing
* 	- Anderson mixing
* 	- Johnson-Pulay mixing
* 	- Kerker preconditioner
*
* ....EJB
*
* Revision 1.48  2003/12/02 23:55:34  bylaska
* density of state generation added...EJB
*
* Revision 1.47  2003/12/02 19:17:10  bylaska
* HGH pseudpotential added.
* TM, Hamman, HGH, pspw_default, and paw_default pseudopotential libraries have been added.
* KS minimizer updates.
* ...EJB
*


************************ f_orb orbitals Part ************************

*     ***********************************
*     *				        *
*     *	     psi_minimize_f_orb         *
*     *				        *
*     ***********************************
      subroutine psi_minimize_f_orb()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

      !*** local variables ***
      integer maxit_orb
      integer ii,l
      real*8  sum,maxerror,error_out,e0

      !*** external functions ***
      real*8   control_tole
      external control_tole

      !call psi_gen_density_potentials(1)
      maxit_orb=120
      maxerror = control_tole()

      do ii=1,(ne(1)+ne(2))

         !*** orthogonalize to lower orbitals  ****
         call psi_project_out_f_orb1(
     >           ii,
     >           dcpl_mb(psi1(1)+(ii-1)*npack1))

         !*** normalize ****
         call Pack_cc_dot(1,
     >            dcpl_mb(psi1(1) +(ii-1)*npack1),
     >            dcpl_mb(psi1(1) +(ii-1)*npack1),
     >            sum)
         sum = 1.0d0/dsqrt(sum)
c         call Pack_c_SMul(1,sum,
c     >            dcpl_mb(psi1(1) +(ii-1)*npack1),
c     >            dcpl_mb(psi1(1) +(ii-1)*npack1))
         call Pack_c_SMul1(1,sum,dcpl_mb(psi1(1) +(ii-1)*npack1))

         !*** minimize orbital ****
          l = 0
 2        call psi_KS_update_f_orb(maxit_orb,
     >                               maxerror,
     >                               0.001d0,ii,error_out,e0)
          !write(*,*) "e0:",ii,l,e0,error_out
          l = l+1
          if ((error_out.gt.maxerror).and.(l.le.4)) go to 2

          dbl_mb(eig(1)+ii-1) = e0

      end do
      call psi_sort_f_orb()
      
     
      return
      end

      subroutine psi_sort_f_orb()
      implicit none
#include "errquit.fh"

#include "mafdecls.fh"
#include "psi.fh"

      logical value
      integer i,j,ii,jj,ms
      integer r1(2)
      real*8  ei,ej

      value = MA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      if (.not. value) call errquit(
     >     'psi_sort_f_orb: out of stack memory',0,MA_ERR)

      do ms=1,ispin

        !*** Bubble sort ***
        do ii=1,ne(ms)
         do jj=ii+1,ne(ms)
           i = ii + (ms-1)*ne(1)
           j = jj + (ms-1)*ne(1)
           ei = dbl_mb(eig(1)+i-1)
           ej = dbl_mb(eig(1)+j-1)

           !*** swap ***
           if (ej.lt.ei) then
             dbl_mb(eig(1)+i-1) = ej
             dbl_mb(eig(1)+j-1) = ei
             call Pack_c_Copy(1,dcpl_mb(psi1(1)+(i-1)*npack1),
     >                          dcpl_mb(r1(1)))
             call Pack_c_Copy(1,dcpl_mb(psi1(1)+(j-1)*npack1),
     >                          dcpl_mb(psi1(1)+(i-1)*npack1))
             call Pack_c_Copy(1,dcpl_mb(r1(1)),
     >                          dcpl_mb(psi1(1)+(j-1)*npack1))
           end if

         end do
        end do

      end do

      value = MA_pop_stack(r1(2))
      if (.not. value) call errquit(
     >     'psi_sort_f_orb: popping stack memory',1, MA_ERR)
      return
      end

*     ***********************************
*     *				        *
*     *	     psi_KS_update_f_orb      *
*     *				        *
*     ***********************************

*    This routine performs a KS update on orbital ii
*
      subroutine psi_KS_update_f_orb(maxiteration,
     >                             maxerror,perror,ii,
     >                             error_out,e0)
      implicit none
#include "errquit.fh"
      integer maxiteration
      real*8  maxerror,perror 
      integer ii
      real*8 error_out
      real*8 e0
      
#include "mafdecls.fh"
#include "psi.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,done,oneloop
      integer it
      real*8 eold,percent_error,error0,de0,lmbda_r0,lmbda_r1
      real*8 theta
      integer r1(2),t0(2),t(2),g(2)
      integer psi_ptr

      psi_ptr=psi1(1)

      call Parallel_taskid(taskid)

      value = MA_push_get(mt_dcpl,npack1,'t0',t0(2),t0(1))
      value = value.and. 
     >        MA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      value = value.and. 
     >        MA_push_get(mt_dcpl,npack1,'g',g(2),g(1))
      value = value.and. 
     >        MA_push_get(mt_dcpl,npack1,'t',t(2),t(1))
      if (.not. value) call errquit(
     >     'psi_KS_update_f_orb: out of stack memory',0, MA_ERR)

      done = .false.
      error0 = 0.0d0 
      e0 = 0.0d0
      theta = -3.14159d0/600.0d0
      lmbda_r0 = 1.0d0
      it = 0
 2    continue

         it = it + 1
         eold = e0

*        *** calculate residual (steepest descent) direction for a single band ***
         call psi_get_gradient_f_orb(ii,dcpl_mb(g(1)))
         call Pack_cc_dot(1,dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                   dcpl_mb(g(1)),
     >                    e0)
   
         e0 = -e0
         percent_error=0d0
         if(error0.ne.0d0)
     A      percent_error = dabs(e0-eold)/error0

         done = ((it.gt.maxiteration) 
     >           .or.
     >           (dabs(e0-eold).lt.maxerror))

         if (done) go to 4


         call Pack_c_Copy(1,dcpl_mb(g(1)),dcpl_mb(r1(1)))
         call Pack_cc_daxpy(1,(e0),
     >                 dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                 dcpl_mb(r1(1)))


         !*** determine conjuagate direction ***
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


*        *** normalize search direction, t ****
         call psi_project_out_f_orb(ii,dcpl_mb(t(1)))
         call Pack_cc_dot(1,dcpl_mb(t(1)),
     >                   dcpl_mb(t(1)),
     >                   de0)
         de0 = 1.0d0/dsqrt(de0)
c         call Pack_c_SMul(1,de0,dcpl_mb(t(1)),dcpl_mb(t(1)))
         call Pack_c_SMul1(1,de0,dcpl_mb(t(1)))
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
         call psi_linesearch_f_orb(ii,
     >                               theta,e0,de0,dcpl_mb(t(1)))

      go to 2


*     **** release stack memory ****
 4    value =           MA_pop_stack(t(2)) 
      value = value.and.MA_pop_stack(g(2))
      value = value.and.MA_pop_stack(r1(2))
      value = value.and.MA_pop_stack(t0(2))
      if (.not. value) call errquit(
     >     'psi_KS_update_f_orb: popping stack memory',1, MA_ERR)

      error_out = dabs(e0-eold)
      e0        = -e0
      return
      end

*     ***********************************
*     *				        *
*     *	     psi_linesearch_f_orb
*     *				        *
*     ***********************************

*    This routine performs a linesearch on orbital ii, in the direction t.  
* This routine is needed for a KS minimizer.
*  e0 = <orb|g>
*  de0 = 2*<t|g>
*
      subroutine psi_linesearch_f_orb(ii,theta,e0,de0,t)
      implicit none
#include "errquit.fh"
      integer ii
      real*8  theta
      real*8  e0,de0
      complex*16 t(*) !search direction
      
#include "mafdecls.fh"
#include "psi.fh"
   
*     **** local variables ****
      logical value
      integer orb(2),g(2),psi_ptr
      real*8 x,y,pi,dtheta_min,e1

      psi_ptr=psi1(1)

      pi = 4.0d0*datan(1.0d0)
      !dtheta = pi/300.0d0
      dtheta_min = 0.01*theta

*     **** allocate stack memory ****
      value = MA_push_get(mt_dcpl,npack1,'orb',
     >                       orb(2),orb(1))
      value = value.and. 
     >        MA_push_get(mt_dcpl,npack1,'g',
     >                       g(2),g(1))
      if (.not. value) call errquit(
     >     'psi_linesearch_f_orb: out of stack memory',0, MA_ERR)
 

      call Pack_c_Copy(1,dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                   dcpl_mb(orb(1)))

*     **** orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
  10  x = cos(theta)
      y = sin(theta)
      call Pack_c_SMul(1,x,
     >                  dcpl_mb(orb(1)),
     >                  dcpl_mb(psi_ptr+(ii-1)*npack1))
      call Pack_cc_daxpy(1,y,
     >                   t,
     >                   dcpl_mb(psi_ptr+(ii-1)*npack1))

*     *** determine theta ***
      call psi_get_gradient_f_orb(ii,dcpl_mb(g(1)))
      call Pack_cc_dot(1,dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                   dcpl_mb(g(1)),
     >                   e1)
      e1 = -e1

     

      x = (e0 - e1 + 0.5d0*de0*sin(2*theta))
     >    /(1.0d0-cos(2*theta))
      theta = 0.5d0*datan(0.5d0*de0/x) 
    
   

*     **** orb2 = orb*cos(theta) + t*sin(theta) ****
      x = cos(theta)
      y = sin(theta)
      call Pack_c_SMul(1,x,
     >                  dcpl_mb(orb(1)),
     >                  dcpl_mb(psi_ptr+(ii-1)*npack1))
      call Pack_cc_daxpy(1,y,
     >                   t,
     >                   dcpl_mb(psi_ptr+(ii-1)*npack1))


*     **** release stack memory ****
      value =           MA_pop_stack(g(2))
      value = value.and.MA_pop_stack(orb(2))      
      if (.not. value) call errquit(
     >     'psi_linesearch_f_orb: popping stack memory',1,MA_ERR)

      return
      end


*     ***********************************
*     *				        *
*     *	     psi_get_gradient_f_orb	*
*     *				        *
*     ***********************************

*    This routine returns the Hpsi(i).  
* This routine is needed for a KS minimizer.
*
      subroutine psi_get_gradient_f_orb(ii,Horb)
      implicit none
      integer ii
      complex*16 Horb(*)

#include "mafdecls.fh"
#include "psi.fh"

*     **** local variables ****
      integer psi_ptr,ms

      psi_ptr=psi1(1)+(ii-1)*npack1

      if (ii.le.neq(1)) then
         ms = 1
      else
         ms = 2
      end if

      call electron_get_gradient_virtual(ms,dcpl_mb(psi_ptr),Horb)
      
      return
      end

*     *******************************************
*     *				                *
*     *	         psi_project_out_f_orb1        *
*     *				                *
*     *******************************************
*
*    This routine projects out non-orthogonal components of Horb.  
* This routine is needed for a KS minimizer.
*
      subroutine psi_project_out_f_orb1(ii,Horb)
      integer ii
      complex*16 Horb(*)

#include "mafdecls.fh"
#include "psi.fh"

      integer ms,jj,kk,shift,shifte
      real*8  sum

*     **** spin up orbital ****
      if (ii.le.neq(1)) then

         shift  = 0
         shifte = 0
         ms     = 1
         kk     = ii
*     **** spin down orbital ****
      else
         shift  = neq(1)*npack1
         shifte = neq(1)*npack1
         ms     = 2
         kk     = ii-neq(1)
      end if 


      !**** project out orbitals ****
      do jj=1,(kk-1)
        call Pack_cc_dot(1,
     >            dcpl_mb(psi1(1) +(jj-1)*npack1+shifte),
     >            Horb,
     >            sum)

        call Pack_cc_daxpy(1,(-sum),
     >            dcpl_mb(psi1(1) +(jj-1)*npack1+shifte),
     >            Horb)
      end do


      return
      end




*     *******************************************
*     *				                *
*     *	         psi_project_out_f_orb          *
*     *				                *
*     *******************************************
*
*    This routine projects out non-orthogonal components of Horb.  
* This routine is needed for a KS minimizer.
*
      subroutine psi_project_out_f_orb(ii,Horb)
      integer ii
      complex*16 Horb(*)

#include "mafdecls.fh"
#include "psi.fh"

      integer ms,jj,kk,shift,shifte
      real*8  sum

*     **** spin up orbital ****
      if (ii.le.neq(1)) then

         shift  = 0
         shifte = 0
         ms     = 1
         kk     = ii
*     **** spin down orbital ****
      else
         shift  = neq(1)*npack1
         shifte = neq(1)*npack1
         ms     = 2
         kk     = ii-neq(1)
      end if 


      !**** project out  orbitals ****
      do jj=1,(kk)
        call Pack_cc_dot(1,
     >            dcpl_mb(psi1(1) +(jj-1)*npack1+shifte),
     >            Horb,
     >            sum)

        call Pack_cc_daxpy(1,(-sum),
     >            dcpl_mb(psi1(1) +(jj-1)*npack1+shifte),
     >            Horb)
      end do


      return
      end




************************ virtural orbital Part ************************
*     ***********************************
*     *				        *
*     *	     psi_minimize_virtual       *
*     *				        *
*     ***********************************
      subroutine psi_minimize_virtual()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

      !*** local variables ***
      integer maxit_orb
      integer ii,l
      real*8  sum,maxerror,error_out,e0

      !*** external functions ***
      real*8   control_tole
      external control_tole

      !call psi_gen_density_potentials(1)
      maxit_orb=120
      maxerror = control_tole()

      do ii=1,(ne_excited(1)+ne_excited(2))

         !*** orthogonalize to lower orbitals  ****
         call psi_project_out_virtual1(
     >           ii,
     >           dcpl_mb(psi1_excited(1)+(ii-1)*npack1))

         !*** normalize ****
         call Pack_cc_dot(1,
     >            dcpl_mb(psi1_excited(1) +(ii-1)*npack1),
     >            dcpl_mb(psi1_excited(1) +(ii-1)*npack1),
     >            sum)
         sum = 1.0d0/dsqrt(sum)
c         call Pack_c_SMul(1,sum,
c     >            dcpl_mb(psi1_excited(1) +(ii-1)*npack1),
c     >            dcpl_mb(psi1_excited(1) +(ii-1)*npack1))
         call Pack_c_SMul1(1,sum,
     >            dcpl_mb(psi1_excited(1) +(ii-1)*npack1))

         !*** minimize orbital ****
          l = 0
 2        call psi_KS_update_virtual(maxit_orb,
     >                               maxerror,
     >                               0.001d0,ii,error_out,e0)
          l = l+1
          if ((error_out.gt.maxerror).and.(l.le.4)) go to 2

          dbl_mb(eig_excited(1)+ii-1) = e0

      end do
      call psi_sort_virtual()
      
     
      return
      end

      subroutine psi_sort_virtual()
      implicit none
#include "errquit.fh"

#include "mafdecls.fh"
#include "psi.fh"

      logical value
      integer i,j,ii,jj,ms
      integer r1(2)
      real*8  ei,ej

      value = MA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      if (.not. value) call errquit(
     >     'psi_sort_virtual: out of stack memory',0, MA_ERR)

      do ms=1,ispin

        !*** Bubble sort ***
        do ii=1,ne_excited(ms)
         do jj=ii+1,ne_excited(ms)
           i = ii + (ms-1)*ne_excited(1)
           j = jj + (ms-1)*ne_excited(1)
           ei = dbl_mb(eig_excited(1)+i-1)
           ej = dbl_mb(eig_excited(1)+j-1)

           !*** swap ***
           if (ej.lt.ei) then
             dbl_mb(eig_excited(1)+i-1) = ej
             dbl_mb(eig_excited(1)+j-1) = ei
             call Pack_c_Copy(1,dcpl_mb(psi1_excited(1)+(i-1)*npack1),
     >                          dcpl_mb(r1(1)))
             call Pack_c_Copy(1,dcpl_mb(psi1_excited(1)+(j-1)*npack1),
     >                          dcpl_mb(psi1_excited(1)+(i-1)*npack1))
             call Pack_c_Copy(1,dcpl_mb(r1(1)),
     >                          dcpl_mb(psi1_excited(1)+(j-1)*npack1))
           end if

         end do
        end do

      end do

      value = MA_pop_stack(r1(2))
      if (.not. value) call errquit(
     >     'psi_sort_virtual: popping stack memory',1, MA_ERR)
      return
      end

*     ***********************************
*     *				        *
*     *	     psi_KS_update_virtual      *
*     *				        *
*     ***********************************

*    This routine performs a KS update on virtual ii
*
      subroutine psi_KS_update_virtual(maxiteration,
     >                             maxerror,perror,ii,
     >                             error_out,e0)
      implicit none
      integer maxiteration
      real*8  maxerror,perror 
      integer ii
      real*8 error_out
      real*8 e0
      
#include "mafdecls.fh"
#include "psi.fh"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,done,oneloop
      integer it
      real*8 eold,percent_error,error0,de0,lmbda_r0,lmbda_r1
      real*8 theta
      integer r1(2),t0(2),t(2),g(2)
      integer psi_ptr

      psi_ptr=psi1_excited(1)

      call Parallel_taskid(taskid)

      value = MA_push_get(mt_dcpl,npack1,'t0',t0(2),t0(1))
      value = value.and. 
     >        MA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      value = value.and. 
     >        MA_push_get(mt_dcpl,npack1,'g',g(2),g(1))
      value = value.and. 
     >        MA_push_get(mt_dcpl,npack1,'t',t(2),t(1))
      if (.not. value) call errquit(
     >     'psi_KS_update_virtual: out of stack memory',0, MA_ERR)

      done = .false.
      error0 = 0.0d0 
      e0 = 0.0d0
      theta = -3.14159d0/600.0d0
      lmbda_r0 = 1.0d0
      it = 0
 2    continue

         it = it + 1
         eold = e0

*        *** calculate residual (steepest descent) direction for a single band ***
         call psi_get_gradient_virtual(ii,dcpl_mb(g(1)))
         call Pack_cc_dot(1,dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                   dcpl_mb(g(1)),
     >                    e0)
   
         e0 = -e0
         
         !if (it.eq.1) then
         !  percent_error = 1.0d0
         !else if (it.eq.2) then
         !  error0 = dabs(e0-eold)
         !  percent_error = 1.0d0
         !else
         percent_error=0.0d0
         if(error0.ne.0.0d0)
     A      percent_error = dabs(e0-eold)/error0
         !end if

         done = ((it.gt.maxiteration) 
     >           .or.
     >           (dabs(e0-eold).lt.maxerror))

         if (done) go to 4


         call Pack_c_Copy(1,dcpl_mb(g(1)),dcpl_mb(r1(1)))
         call Pack_cc_daxpy(1,(e0),
     >                 dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                 dcpl_mb(r1(1)))


         !*** determine conjuagate direction ***
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


*        *** normalize search direction, t ****
         call psi_project_out_virtual(ii,dcpl_mb(t(1)))
         call Pack_cc_dot(1,dcpl_mb(t(1)),
     >                   dcpl_mb(t(1)),
     >                   de0)
         de0 = 1.0d0/dsqrt(de0)
c         call Pack_c_SMul(1,de0,dcpl_mb(t(1)),dcpl_mb(t(1)))
         call Pack_c_SMul1(1,de0,dcpl_mb(t(1)))
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
         call psi_linesearch_virtual(ii,
     >                               theta,e0,de0,dcpl_mb(t(1)))

      go to 2


*     **** release stack memory ****
 4    value =           MA_pop_stack(t(2)) 
      value = value.and.MA_pop_stack(g(2))
      value = value.and.MA_pop_stack(r1(2))
      value = value.and.MA_pop_stack(t0(2))
      if (.not. value) call errquit(
     >     'psi_KS_update_virtual: popping stack memory',1, MA_ERR)

      error_out = dabs(e0-eold)
      e0        = -e0
      return
      end

*     ***********************************
*     *				        *
*     *	     psi_linesearch_virtual
*     *				        *
*     ***********************************

*    This routine performs a linesearch on orbital ii, in the direction t.  
* This routine is needed for a KS minimizer.
*  e0 = <orb|g>
*  de0 = 2*<t|g>
*
      subroutine psi_linesearch_virtual(ii,theta,e0,de0,t)
      implicit none
      integer ii
      real*8  theta
      real*8  e0,de0
      complex*16 t(*) !search direction
#include "errquit.fh"
#include "mafdecls.fh"
#include "psi.fh"
   
*     **** local variables ****
      logical value
      integer orb(2),g(2),psi_ptr
      real*8 x,y,pi,dtheta_min,e1

      psi_ptr=psi1_excited(1)

      pi = 4.0d0*datan(1.0d0)
      !dtheta = pi/300.0d0
      dtheta_min = 0.01*theta

*     **** allocate stack memory ****
      value = MA_push_get(mt_dcpl,npack1,'orb',
     >                       orb(2),orb(1))
      value = value.and. 
     >        MA_push_get(mt_dcpl,npack1,'g',
     >                       g(2),g(1))
      if (.not. value) call errquit(
     >     'psi_linesearch_virtual: out of stack memory',0, MA_ERR)
 

      call Pack_c_Copy(1,dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                   dcpl_mb(orb(1)))

*     **** orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
  10  x = cos(theta)
      y = sin(theta)
      call Pack_c_SMul(1,x,
     >                  dcpl_mb(orb(1)),
     >                  dcpl_mb(psi_ptr+(ii-1)*npack1))
      call Pack_cc_daxpy(1,y,
     >                   t,
     >                   dcpl_mb(psi_ptr+(ii-1)*npack1))

*     *** determine theta ***
      call psi_get_gradient_virtual(ii,dcpl_mb(g(1)))
      call Pack_cc_dot(1,dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                   dcpl_mb(g(1)),
     >                   e1)
      e1 = -e1

     
      !if (((-e1).gt.(-e0)).and.(theta.gt.dtheta_min)) then
      !   theta = 0.5d0*theta
      !   go to 10
      !end if

      x = (e0 - e1 + 0.5d0*de0*sin(2*theta))
     >    /(1.0d0-cos(2*theta))
      theta = 0.5d0*datan(0.5d0*de0/x) 
    
   

*     **** orb2 = orb*cos(theta) + t*sin(theta) ****
      x = cos(theta)
      y = sin(theta)
      call Pack_c_SMul(1,x,
     >                  dcpl_mb(orb(1)),
     >                  dcpl_mb(psi_ptr+(ii-1)*npack1))
      call Pack_cc_daxpy(1,y,
     >                   t,
     >                   dcpl_mb(psi_ptr+(ii-1)*npack1))


*     **** release stack memory ****
      value =           MA_pop_stack(g(2))
      value = value.and.MA_pop_stack(orb(2))      
      if (.not. value) call errquit(
     >     'psi_linesearch_virtual: popping stack memory',1, MA_ERR)

      return
      end


*     ***********************************
*     *				        *
*     *	     psi_get_gradient_virtual	*
*     *				        *
*     ***********************************

*    This routine returns the Hpsi(i).  
* This routine is needed for a KS minimizer.
*
      subroutine psi_get_gradient_virtual(ii,Horb)
      implicit none
      integer ii
      complex*16 Horb(*)

#include "mafdecls.fh"
#include "psi.fh"

*     **** local variables ****
      integer psi_ptr,ms

      psi_ptr=psi1_excited(1)+(ii-1)*npack1

      if (ii.le.ne_excited(1)) then
         ms = 1
      else
         ms = 2
      end if

      call electron_get_gradient_virtual(ms,dcpl_mb(psi_ptr),Horb)
      
      return
      end

*     *******************************************
*     *				                *
*     *	         psi_project_out_virtual1        *
*     *				                *
*     *******************************************
*
*    This routine projects out non-orthogonal components of Horb.  
* This routine is needed for a KS minimizer.
*
      subroutine psi_project_out_virtual1(ii,Horb)
      integer ii
      complex*16 Horb(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

      integer ms,i,jj,kk,shift,shifte
      integer etmp(2)
      real*8  sum

*     **** spin up orbital ****
      if (ii.le.ne_excited(1)) then

         shift  = 0
         shifte = 0
         ms     = 1
         kk     = ii
*     **** spin down orbital ****
      else
         shift  = neq(1)*npack1
         shifte = ne_excited(1)*npack1
         ms     = 2
         kk     = ii-ne_excited(1)
      end if 

      !**** project out filled orbitals ****
      if (neq(ms).eq.ne(ms)) then

         do i=1,ne(ms)
           call Pack_cc_dot(1,
     >            dcpl_mb(psi1(1) +(i-1)*npack1+shift),
     >            Horb,
     >            sum)
           call Pack_cc_daxpy(1,(-sum),
     >            dcpl_mb(psi1(1) +(i-1)*npack1+shift),
     >            Horb)
         end do      

      else

      if (.not.MA_push_get(mt_dcpl,npack1,'etmp',etmp(2),etmp(1)))
     > call errquit('psi_project_out_virtual1: out of stack',0,MA_ERR)

         call dcopy(2*npack1,0.0d0,0,dcpl_mb(etmp(1)),1)
         do i=1,neq(ms)
           call Pack_cc_dot(1,
     >            dcpl_mb(psi1(1) +(i-1)*npack1+shift),
     >            Horb,
     >            sum)
           call Pack_cc_daxpy(1,(-sum),
     >            dcpl_mb(psi1(1) +(i-1)*npack1+shift),
     >            dcpl_mb(etmp(1)))
         end do
         call D1dB_Vector_SumAll(2*npack1,dcpl_mb(etmp(1)))
         call daxpy(2*npack1,1.0d0,dcpl_mb(etmp(1)),1,Horb,1)

      if (.not.MA_pop_stack(etmp(2)))
     > call errquit('psi_project_out_virtual1:popping stack',0,MA_ERR)

      end if

      !**** project out virtual orbitals ****
      do jj=1,(kk-1)
        call Pack_cc_dot(1,
     >            dcpl_mb(psi1_excited(1) +(jj-1)*npack1+shifte),
     >            Horb,
     >            sum)

        call Pack_cc_daxpy(1,(-sum),
     >            dcpl_mb(psi1_excited(1) +(jj-1)*npack1+shifte),
     >            Horb)
      end do


      return
      end




*     *******************************************
*     *				                *
*     *	         psi_project_out_virtual        *
*     *				                *
*     *******************************************
*
*    This routine projects out non-orthogonal components of Horb.  
* This routine is needed for a KS minimizer.
*
      subroutine psi_project_out_virtual(ii,Horb)
      integer ii
      complex*16 Horb(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

      integer ms,i,jj,kk,shift,shifte
      integer etmp(2)
      real*8  sum

*     **** spin up orbital ****
      if (ii.le.ne_excited(1)) then

         shift  = 0
         shifte = 0
         ms     = 1
         kk     = ii
*     **** spin down orbital ****
      else
         shift  = neq(1)*npack1
         shifte = ne_excited(1)*npack1
         ms     = 2
         kk     = ii-ne_excited(1)
      end if 

      !**** project out filled orbitals ****
      if (neq(ms).eq.ne(ms)) then

         do i=1,ne(ms)
           call Pack_cc_dot(1,
     >            dcpl_mb(psi1(1) +(i-1)*npack1+shift),
     >            Horb,
     >            sum)

           call Pack_cc_daxpy(1,(-sum),
     >            dcpl_mb(psi1(1) +(i-1)*npack1+shift),
     >            Horb)
         end do      

      else
      if (.not.MA_push_get(mt_dcpl,npack1,'etmp',etmp(2),etmp(1)))
     > call errquit('psi_project_out_virtual1: out of stack',0,MA_ERR)

         call dcopy(2*npack1,0.0d0,0,dcpl_mb(etmp(1)),1)
         do i=1,neq(ms)
           call Pack_cc_dot(1,
     >            dcpl_mb(psi1(1) +(i-1)*npack1+shift),
     >            Horb,
     >            sum)
           call Pack_cc_daxpy(1,(-sum),
     >            dcpl_mb(psi1(1) +(i-1)*npack1+shift),
     >            dcpl_mb(etmp(1)))
         end do
         call D1dB_Vector_SumAll(2*npack1,dcpl_mb(etmp(1)))
         call daxpy(2*npack1,1.0d0,dcpl_mb(etmp(1)),1,Horb,1)

      if (.not.MA_pop_stack(etmp(2)))
     > call errquit('psi_project_out_virtual:popping stack',0,MA_ERR)

      end if

      !**** project out virtual orbitals ****
      do jj=1,(kk)
        call Pack_cc_dot(1,
     >            dcpl_mb(psi1_excited(1) +(jj-1)*npack1+shifte),
     >            Horb,
     >            sum)

        call Pack_cc_daxpy(1,(-sum),
     >            dcpl_mb(psi1_excited(1) +(jj-1)*npack1+shifte),
     >            Horb)
      end do


      return
      end



************************ KS orbital Part ************************

*     ***********************************
*     *				        *
*     *	     psi_KS_update	        *
*     *				        *
*     ***********************************

*    This routine (approximately) diagonalizes the KS matrix.
*
      subroutine psi_KS_update(psi_number,
     >                         ks_algorithm,
     >                         precondition,
     >                         maxerror)
      implicit none
      integer psi_number
      integer ks_algorithm
      logical precondition
      real*8 maxerror
        
#include "mafdecls.fh"
#include "psi.fh"
    
*     **** local variables ****
      logical done
      integer i,j,neall,maxit_orb,maxit_orbs
      real*8 error,error_out,tim1,tim2,tim,sum

*     **** external functions ****
      integer  control_ks_maxit_orb,control_ks_maxit_orbs
      external control_ks_maxit_orb,control_ks_maxit_orbs

      tim = 0.0d0
      neall = neq(1)+neq(2)
      maxit_orb  = control_ks_maxit_orb()   !*** should be read from rtdb ***
      maxit_orbs = control_ks_maxit_orbs()  !*** should be read from rtdb ***
      j = 0
 2    j = j+1
        error = 0.0d0
        !do i=neall,1,-1
        do i=1,neall

         !*** orthogonalize to lower orbitals  ****
         call psi_project_out_f_orb1(
     >           i,
     >           dcpl_mb(psi1(1)+(i-1)*npack1))

         !*** normalize ****
         call Pack_cc_dot(1,
     >            dcpl_mb(psi1(1) +(i-1)*npack1),
     >            dcpl_mb(psi1(1) +(i-1)*npack1),
     >            sum)
         sum = 1.0d0/dsqrt(sum)
c         call Pack_c_SMul(1,sum,
c     >            dcpl_mb(psi1(1) +(i-1)*npack1),
c     >            dcpl_mb(psi1(1) +(i-1)*npack1))
         call Pack_c_SMul1(1,sum,
     >            dcpl_mb(psi1(1) +(i-1)*npack1))



          if (ks_algorithm.eq.1) then
          call psi_KS_update_orb2(psi_number,precondition,maxit_orb,
     >                         maxerror,
     >                         0.1d0,i,error_out)
          else
          call psi_KS_update_orb(psi_number,precondition,maxit_orb,
     >                         maxerror,
     >                         0.1d0,i,error_out)
          end if

          error = error+error_out
        end do
        error = error/dble(neall)

        done = ((j.gt.maxit_orbs).or.(error.lt.maxerror))
      if (.not.done) go to 2

      return
      end


*     ***********************************
*     *				        *
*     *	     psi_KS_update_orb	        *
*     *				        *
*     ***********************************

*    This routine performs a KS update on orbital i
*
      subroutine psi_KS_update_orb(psi_number,
     >                             precondition,maxiteration,
     >                             maxerror,perror,i,
     >                             error_out)
      implicit none
#include "errquit.fh"
      integer psi_number
      logical precondition
      integer maxiteration
      real*8  maxerror,perror 
      integer i
      real*8 error_out
      
#include "mafdecls.fh"
#include "psi.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,done,oneloop
      integer it
      real*8 e0,eold,error0,de0,lmbda_r0,lmbda_r1
      real*8 theta,sigma
      integer r1(2),t0(2),t(2),g(2)
      integer psi_ptr

      if (psi_number.eq.1) then
         psi_ptr=psi1(1)
      else 
         psi_ptr=psi2(1)
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
     >     'psi_KS_update_orb: out of stack memory',0, MA_ERR)

      done = .false.
      error0 = 0.0d0 
      e0 = 0.0d0
      theta = -3.14159d0/600.0d0
      lmbda_r0 = 1.0d0
      it = 0
 2    continue

         it = it + 1
         eold = e0

*        *** calculate residual (steepest descent) direction for a single band ***
         call psi_get_gradient_orb(psi_number,i,dcpl_mb(g(1)))
         call Pack_cc_dot(1,dcpl_mb(psi_ptr+(i-1)*npack1),
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
c    >                 dcpl_mb(psi_ptr+(i-1)*npack1),
c    >                 dcpl_mb(r1(1)))
         call Pack_c_Copy(1,dcpl_mb(g(1)),dcpl_mb(r1(1)))
         call psi_project_out_orb(psi_number,i,dcpl_mb(r1(1)))

     


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
c!        call psi_project_out_orb(psi_number,i,dcpl_mb(t(1)))
c!        call Pack_cc_dot(1,dcpl_mb(psi_ptr+(i-1)*npack1),
c!    >                   dcpl_mb(t(1)),
c!    >                    de0)
c!        de0 = -de0
c!        call Pack_cc_daxpy(1,(de0),
c!    >                 dcpl_mb(psi_ptr+(i-1)*npack1),
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
         call psi_linesearch_update2(psi_number,i,
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




*     ***********************************
*     *				        *
*     *	     psi_KS_update_orb2	        *
*     *				        *
*     ***********************************

*    This routine performs a RMM-DIIS KS update on orbital i
*
      subroutine psi_KS_update_orb2(psi_number,
     >                             precondition,maxiteration,
     >                             maxerror,perror,i,
     >                             error_out)
      implicit none
      integer psi_number
      logical precondition
      integer maxiteration
      real*8  maxerror,perror 
      integer i
      real*8 error_out
      
#include "mafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,done
      integer it
      real*8 sigma,e0,eold,error0
      real*8 lambda
      integer r1(2),g1(2)
      integer psi_ptr

      if (psi_number.eq.1) then
         psi_ptr=psi1(1)
      else 
         psi_ptr=psi2(1)
      end if

      call Parallel_taskid(taskid)

*     **** allocate memory ****
      value =            MA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      value = value.and. MA_push_get(mt_dcpl,npack1,'g1',g1(2),g1(1))
      if (.not. value) call errquit(
     >     'psi_KS_update_orb2: out of stack memory',0, MA_ERR)

*     **** set lambda ***
      lambda = 0.1d0


*     *** calculate residual (steepest descent) direction for a single band ***
      call psi_get_gradient_orb(psi_number,i,dcpl_mb(g1(1)))
      call Pack_cc_dot(1,
     >                 dcpl_mb(psi_ptr+(i-1)*npack1),
     >                 dcpl_mb(g1(1)),
     >                 e0)
      call Pack_cc_dot(1,
     >                 dcpl_mb(psi_ptr+(i-1)*npack1),
     >                 dcpl_mb(psi_ptr+(i-1)*npack1),
     >                 sigma)
      e0 = e0/sigma
      call Pack_c_SMul(1,e0,
     >                 dcpl_mb(psi_ptr+(i-1)*npack1),
     >                 dcpl_mb(r1(1)))
      call Pack_cc_daxpy(1,(-1.0d0),
     >                 dcpl_mb(g1(1)),
     >                 dcpl_mb(r1(1)))
      
        !write(*,*) "i=",i,"it=",0, " eig=",e0,sigma

*     ***** rmmdiis start *****
      call pspw_rmmdiis_start(lambda,
     >                        dcpl_mb(r1(1)),
     >                        dcpl_mb(psi_ptr+(i-1)*npack1))


      done = .false.
      error0 = 0.0d0 
      it = 0
 2    continue

         it = it + 1
         eold = e0

*        *** calculate residual (steepest descent) direction for a single band ***
         call psi_get_gradient_orb(psi_number,i,dcpl_mb(g1(1)))
         call Pack_cc_dot(1,
     >                    dcpl_mb(psi_ptr+(i-1)*npack1),
     >                    dcpl_mb(g1(1)),
     >                    e0)
         call Pack_cc_dot(1,
     >                    dcpl_mb(psi_ptr+(i-1)*npack1),
     >                    dcpl_mb(psi_ptr+(i-1)*npack1),
     >                    sigma)
         e0 = e0/sigma
         call Pack_c_SMul(1,(e0),
     >                    dcpl_mb(psi_ptr+(i-1)*npack1),
     >                    dcpl_mb(r1(1)))
         call Pack_cc_daxpy(1,(-1.0d0),
     >                    dcpl_mb(g1(1)),
     >                    dcpl_mb(r1(1)))
         !e0 = -e0
        !write(*,*) "i=",i,"it=",it, " eig=",e0,sigma,dabs(e0-eold)

         done = ((it.gt.maxiteration) 
     >           .or.
     >           (dabs(e0-eold).lt.maxerror))

         if (done) go to 4

*        ***** rmmdiis update *****
         call pspw_rmmdiis(lambda,
     >                     dcpl_mb(r1(1)),
     >                     dcpl_mb(psi_ptr+(i-1)*npack1))

      go to 2

*     ***** normalize psi ****
 4    call Pack_cc_dot(1,
     >                 dcpl_mb(psi_ptr+(i-1)*npack1),
     >                 dcpl_mb(psi_ptr+(i-1)*npack1),
     >                 sigma)
c      call Pack_c_SMul(1,(1.0d0/dsqrt(sigma)),
c     >                 dcpl_mb(psi_ptr+(i-1)*npack1),
c     >                 dcpl_mb(psi_ptr+(i-1)*npack1))
      call Pack_c_SMul1(1,(1.0d0/dsqrt(sigma)),
     >                 dcpl_mb(psi_ptr+(i-1)*npack1))


*     **** release stack memory ****
      value =           MA_pop_stack(g1(2))
      value = value.and.MA_pop_stack(r1(2))
      if (.not. value) call errquit(
     >     'psi_KS_update_orb2: popping stack memory',1, MA_ERR)
      error_out = dabs(e0-eold)

c       write(*,*) "i=",i,"iterations=",it," eig=",e0,
c     >            " error=",error_out,
c     >            lambda
      return
      end








*     ***********************************
*     *				        *
*     *	     psi_linesearch_update	*
*     *				        *
*     ***********************************

*    This routine performs a linesearch on orbital i, in the direction t.  
* This routine is needed for a KS minimizer.
*  e0 = <orb|g>
*  de0 = 2*<t|g>
*
      subroutine psi_linesearch_update(psi_number,i,theta,e0,de0,t)
      implicit none
#include "errquit.fh"
      integer psi_number
      integer i
      real*8  theta
      real*8  e0,de0
      complex*16 t(*) !search direction
      
#include "mafdecls.fh"
#include "psi.fh"
   
*     **** local variables ****
      logical value
      integer orb(2),g(2),psi_ptr
      real*8 x,y,pi,e1

      if (psi_number.eq.1) then
         psi_ptr=psi1(1)
      else 
         psi_ptr=psi2(1)
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
 

      call Pack_c_Copy(1,dcpl_mb(psi_ptr+(i-1)*npack1),
     >                   dcpl_mb(orb(1)))

*     **** orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
      !theta = pi/300.0d0
      x = cos(theta)
      y = sin(theta)
      call Pack_c_SMul(1,x,
     >                  dcpl_mb(orb(1)),
     >                  dcpl_mb(psi_ptr+(i-1)*npack1))
      call Pack_cc_daxpy(1,y,
     >                   t,
     >                   dcpl_mb(psi_ptr+(i-1)*npack1))

*     *** determine theta ***
      call psi_get_gradient_orb(psi_number,i,dcpl_mb(g(1)))

      call Pack_cc_dot(1,dcpl_mb(psi_ptr+(i-1)*npack1),
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
     >                  dcpl_mb(psi_ptr+(i-1)*npack1))
      call Pack_cc_daxpy(1,y,
     >                   t,
     >                   dcpl_mb(psi_ptr+(i-1)*npack1))

*     **** update orb2_r and H*orb2 ****
      !call electron_run_orb(i,dcpl_mb(psi_ptr)) 
c     call psi_get_gradient_orb(psi_number,i,dcpl_mb(g(1)))
c     call Pack_cc_dot(1,dcpl_mb(psi_ptr+(i-1)*npack1),
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
      value =           MA_pop_stack(g(2))
      value = value.and.MA_pop_stack(orb(2))      
      if (.not. value) call errquit(
     >     'psi_linesearch_update: popping stack memory',1, MA_ERR)

      return
      end

*     ***********************************
*     *				        *
*     *	     psi_linesearch_update2	*
*     *				        *
*     ***********************************

*    This routine performs a linesearch on orbital i, in the direction t.  
* This routine is needed for a KS minimizer.
*  e0 = <orb|g>
*  de0 = 2*<t|g>
*
      subroutine psi_linesearch_update2(psi_number,i,theta,e0,de0,t,
     >                                  sigma,tau_t)
      implicit none
#include "errquit.fh"
      integer psi_number
      integer i
      real*8  theta
      real*8  e0,de0
      complex*16 t(*)     !search direction

      real*8     sigma
      complex*16 tau_t(*) !parallel transported search direction
      
#include "mafdecls.fh"
#include "psi.fh"
   
*     **** local variables ****
      logical value
      integer orb(2),g(2),psi_ptr
      real*8 x,y,pi,e1

      if (psi_number.eq.1) then
         psi_ptr=psi1(1)
      else 
         psi_ptr=psi2(1)
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
 

      call Pack_c_Copy(1,dcpl_mb(psi_ptr+(i-1)*npack1),
     >                   dcpl_mb(orb(1)))

*     **** orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
      !theta = pi/300.0d0
      x = cos(theta)
      y = sin(theta)
      call Pack_c_SMul(1,x,
     >                  dcpl_mb(orb(1)),
     >                  dcpl_mb(psi_ptr+(i-1)*npack1))
      call Pack_cc_daxpy(1,y,
     >                   t,
     >                   dcpl_mb(psi_ptr+(i-1)*npack1))

*     *** determine theta ***
      call psi_get_gradient_orb(psi_number,i,dcpl_mb(g(1)))

      call Pack_cc_dot(1,dcpl_mb(psi_ptr+(i-1)*npack1),
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
     >                  dcpl_mb(psi_ptr+(i-1)*npack1))
      call Pack_cc_daxpy(1,y,
     >                   t,
     >                   dcpl_mb(psi_ptr+(i-1)*npack1))


*     **** release stack memory ****
      value =           MA_pop_stack(g(2))
      value = value.and.MA_pop_stack(orb(2))      
      if (.not. value) call errquit(
     >     'psi_linesearch_update: popping stack memory',1, MA_ERR)

      return
      end



*     ***************************
*     *				*
*     *	     psi_set_orb	*
*     *				*
*     ***************************

*    This routine copies an orbital, orb, into the ith psi of psi1.  
* This routine is needed for a KS minimizer.
*
      subroutine psi_set_orb(psi_number,i,orb)
      implicit none
      integer psi_number
      integer i
      complex*16 orb(*)

#include "mafdecls.fh"
#include "psi.fh"
   
*     **** local variables ****
      integer index,psi_ptr

      if (psi_number.eq.1) then
         psi_ptr=psi1(1)
      else 
         psi_ptr=psi2(1)
      end if

      index = (i-1)*npack1

      call zcopy(npack1,
     >           orb, 1,
     >           dcpl_mb(psi_ptr+index),1)
      return
      end


*     ***************************
*     *				*
*     *	     psi_get_orb	*
*     *				*
*     ***************************

*    This routine copies the ith psi of psi1 into an orbital, orb.  
* This routine is needed for a KS minimizer.
*
      subroutine psi_get_orb(psi_number,i,orb)
      implicit none
      integer psi_number
      integer i
      complex*16 orb(*)

#include "mafdecls.fh"
#include "psi.fh"
   
*     **** local variables ****
      integer index,psi_ptr


      if (psi_number.eq.1) then
         psi_ptr=psi1(1)
      else 
         psi_ptr=psi2(1)
      end if

      index = (i-1)*npack1

      call zcopy(npack1,
     >           dcpl_mb(psi_ptr+index), 1,
     >           orb, 1)
      return
      end

*     ***********************************
*     *				        *
*     *	     psi_get_gradient_orb	*
*     *				        *
*     ***********************************

*    This routine returns the Hpsi(i).  
* This routine is needed for a KS minimizer.
*
      subroutine psi_get_gradient_orb(psi_number,i,Horb)
      implicit none
      integer psi_number
      integer i
      complex*16 Horb(*)

#include "mafdecls.fh"
#include "psi.fh"

*     **** local variables ****
      integer psi_ptr

      if (psi_number.eq.1) then
         psi_ptr=psi1(1)
      else 
         psi_ptr=psi2(1)
      end if

      call electron_run_orb(i,dcpl_mb(psi_ptr))
      call electron_get_gradient_orb(i,Horb)
      
      return
      end


*     *******************************************
*     *				                *
*     *	         psi_project_out_orb           *
*     *				                *
*     *******************************************
*
*    This routine projects out non-orthogonal components of Horb.  
* This routine is needed for a KS minimizer.
*
      subroutine psi_project_out_orb(psi_number,i,Horb)
      implicit none
#include "errquit.fh"
      integer psi_number
      integer i
      complex*16 Horb(*)

#include "mafdecls.fh"
#include "psi.fh"

*     **** local variables ****
      logical ok
      integer ii,n,psi_ptr,np
      integer x(2)
      real*8  sum

      call Parallel_np(np)

*     **** allocate stack memory ****
      ok = MA_push_get(mt_dbl,ne(1),'x',x(2),x(1))
      if (.not.ok) 
     > call errquit('psi_project_out_orb: out of stack memory',0,
     &       MA_ERR)


      if (psi_number.eq.1) then
         psi_ptr=psi1(1)
      else 
         psi_ptr=psi2(1)
      end if

*     **** spin up orbital **** 
      if (i.le.ne(1)) then

        ii = i
!       do n=1,(ii)
!          call Pack_cc_dot(1,
!    >            dcpl_mb(psi_ptr +(n-1)*npack1),
!    >            Horb,
!    >            sum)
!          call daxpy(2*npack1,
!    >               (-sum),
!    >               dcpl_mb(psi_ptr+(n-1)*npack1),1,
!    >               Horb,1) 
!       end do     
        call Pack_cc_ndot(1,ii,
     >            dcpl_mb(psi_ptr),
     >            Horb,
     >            dbl_mb(x(1)))
        do n=1,(ii)
           call daxpy(2*npack1,
     >               (-dbl_mb(x(1)+n-1)),
     >               dcpl_mb(psi_ptr+(n-1)*npack1),1,
     >               Horb,1) 
        end do     



*     **** spin down orbital ****      
      else       


        ii = i - ne(1)
        do n=(ne(1)+1),(ne(1)+ii)
           call Pack_cc_dot(1,
     >            dcpl_mb(psi_ptr +(n-1)*npack1),
     >            Horb,
     >            sum)
           call daxpy(2*npack1,
     >               (-sum),
     >               dcpl_mb(psi_ptr+(n-1)*npack1),1,
     >               Horb,1) 
        end do     


      end if

*     **** release stack memory ****
      ok = MA_pop_stack(x(2))
      if (.not. ok) 
     > call errquit('psi_project_out_orb: poping stack memory',0,
     &       MA_ERR)
 
      return
      end






*     ***************************
*     *				*
*     *	     psi_set_density	*
*     *				*
*     ***************************

*    This routine sets the densities and potentials in psi and electron.  
* This routine is needed for a band by band minimizer.
*
      subroutine psi_set_density(psi_number,rho)
      implicit none
      integer psi_number
      real*8 rho(*)


#include "mafdecls.fh"
#include "psi.fh"
   
*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** local variables ****
      integer rho_ptr,dng_ptr,rho_all_ptr

      if (psi_number.eq.1) then
        rho_ptr     = rho1(1)
        dng_ptr     = dng1(1)
        rho_all_ptr = rho1_all(1)
      else
        rho_ptr     = rho2(1)
        dng_ptr     = dng2(1)
        rho_all_ptr = rho2_all(1)
      end if


      call dcopy(4*nfft3d,
     >           rho, 1,
     >           dbl_mb(rho_ptr),1)

      call electron_gen_dng_dnall(dbl_mb(rho_ptr),
     >                            dcpl_mb(dng_ptr),
     >                            dbl_mb(rho_all_ptr))
      call electron_gen_scf_potentials(dbl_mb(rho_ptr),
     >                            dcpl_mb(dng_ptr),
     >                            dbl_mb(rho_all_ptr))
      call electron_gen_vall()
      return
      end


*     ***************************
*     *				*
*     *	     psi_get_density	*
*     *				*
*     ***************************

*    This routine gets the densities in psi.  
* This routine is needed for a band by band minimizer.
*
      subroutine psi_get_density(psi_number,rho)
      implicit none
      integer psi_number
      real*8 rho(*)


#include "mafdecls.fh"
#include "psi.fh"
   
*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** local variables ****
      integer rho_ptr

      if (psi_number.eq.1) then
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
*     *			   	           *
*     *	     psi_gen_density_potentials	   *
*     *				           *
*     **************************************

*    This routine sets the densities and potentials in psi and electron.  
* This routine is needed for a band by band minimizer.
*
      subroutine psi_gen_density_potentials(psi_number)
      implicit none
      integer psi_number


#include "mafdecls.fh"
#include "psi.fh"
   
*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** local variables ****
      integer psi_ptr,rho_ptr,dng_ptr,rho_all_ptr,occ_ptr

      if (psi_number.eq.1) then
        psi_ptr     = psi1(1)
        rho_ptr     = rho1(1)
        dng_ptr     = dng1(1)
        rho_all_ptr = rho1_all(1)
        occ_ptr     = occ1(1)
      else
        psi_ptr     = psi2(1)
        rho_ptr     = rho2(1)
        dng_ptr     = dng2(1)
        rho_all_ptr = rho2_all(1)
        occ_ptr     = occ2(1)
      end if


      call electron_gen_densities(dcpl_mb(psi_ptr),
     >                            dbl_mb(rho_ptr),
     >                            dcpl_mb(dng_ptr),
     >                            dbl_mb(rho_all_ptr),
     >                            occupation_on,dbl_mb(occ_ptr))
      call electron_gen_scf_potentials(dbl_mb(rho_ptr),
     >                            dcpl_mb(dng_ptr),
     >                            dbl_mb(rho_all_ptr))
      call electron_gen_vall()
      return
      end


************************ Grasmman orbitals Part ************************

*     ***************************
*     *				*
*     *		psi_1to2	*
*     *				*
*     ***************************
      subroutine psi_1to2()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"
   
      call zcopy(npack1*(neq(1)+neq(2)),
     >           dcpl_mb(psi1(1)),1,
     >           dcpl_mb(psi2(1)),1)

      return
      end


*     ***************************
*     *				*
*     *		psi_2to1	*
*     *				*
*     ***************************
      subroutine psi_2to1()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

 
      call zcopy(npack1*(neq(1)+neq(2)),
     >           dcpl_mb(psi2(1)),1,
     >           dcpl_mb(psi1(1)),1)

c      call OrthoCheck(ispin,ne,dcpl_mb(psi1(1)))  
      return
      end


*     ***************************
*     *                         *
*     *         epsi_2to1        *
*     *                         *
*     ***************************
      subroutine epsi_2to1()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

      call zcopy(npack1*(ne_excited(1)+ne_excited(2)),
     >           dcpl_mb(psi2_excited(1)),1,
     >           dcpl_mb(psi1_excited(1)),1)
      return
      end


*     ***************************
*     *                         *
*     *         epsi_1to2       *
*     *                         *
*     ***************************
      subroutine epsi_1to2()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

      call zcopy(npack1*(ne_excited(1)+ne_excited(2)),
     >           dcpl_mb(psi1_excited(1)),1,
     >           dcpl_mb(psi2_excited(1)),1)
      return
      end



*     ***************************
*     *				*
*     *		psi_1get_psi	*
*     *				*
*     ***************************
      subroutine psi_1get_psi(rpsi)
      implicit none
      complex*16 rpsi(*)

#include "mafdecls.fh"
#include "psi.fh"
   
      call zcopy(npack1*(neq(1)+neq(2)),
     >           dcpl_mb(psi1(1)),1,
     >           rpsi,1)

      return
      end


*     ***************************
*     *				*
*     *		psi_2get_psi	*
*     *				*
*     ***************************
      subroutine psi_2get_psi(rpsi)
      implicit none
      complex*16 rpsi(*)

#include "mafdecls.fh"
#include "psi.fh"
   
      call zcopy(npack1*(neq(1)+neq(2)),
     >           dcpl_mb(psi2(1)),1,
     >           rpsi,1)

      return
      end

*     ***************************
*     *				*
*     *		psi_check	*
*     *				*
*     ***************************
      subroutine psi_check()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"
 

      call OrthoCheck(ispin,ne,dcpl_mb(psi1(1)))  
      return
      end



*     ***************************
*     *				*
*     *		rho_2to1	*
*     *				*
*     ***************************
      subroutine rho_2to1()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

      call dcopy(4*nfft3d,
     >           dbl_mb(rho2(1)),1,
     >           dbl_mb(rho1(1)),1)

      call dcopy(4*nfft3d,
     >           dbl_mb(rho2_all(1)),1,
     >           dbl_mb(rho1_all(1)),1)

      return
      end

*     ***************************
*     *				*
*     *		rho_1to2	*
*     *				*
*     ***************************
      subroutine rho_1to2()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

      call dcopy(4*nfft3d,
     >           dbl_mb(rho1(1)),1,
     >           dbl_mb(rho2(1)),1)

      call dcopy(4*nfft3d,
     >           dbl_mb(rho1_all(1)),1,
     >           dbl_mb(rho2_all(1)),1)

      return
      end

*     ***************************
*     *				*
*     *		dng_2to1	*
*     *				*
*     ***************************
      subroutine dng_2to1()
      implicit none
 
#include "mafdecls.fh"
#include "psi.fh"
 
      call zcopy(npack0,
     >           dcpl_mb(dng2(1)),1,
     >           dcpl_mb(dng1(1)),1)

      return
      end


*     ***********************************
*     *					*
*     *		psi_1add_oep_to_vall	*
*     *					*
*     ***********************************
      subroutine psi_1add_oep_to_vall()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

 
      call electron_add_oep_to_vall(dbl_mb(rho1(1)))

      return
      end


*     ***********************************
*     *					*
*     *		psi_1toelectron		*
*     *					*
*     ***********************************
      subroutine psi_1toelectron()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all
 
      call electron_run(dcpl_mb(psi1(1)),
     >                  dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)),
     >                  dbl_mb(rho1_all(1)),
     >                  occupation_on,dbl_mb(occ1(1)))

      return
      end

*     ***********************************
*     *					*
*     *		psi_1genrho		*
*     *					*
*     ***********************************
      subroutine psi_1genrho()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all
 
      call electron_genrho(dcpl_mb(psi1(1)),
     >                     dbl_mb(rho1(1)),
     >                     occupation_on,dbl_mb(occ1(1)))

      return
      end



*     ***********************************
*     *					*
*     *		psi_1energy		*
*     *					*
*     ***********************************
      real*8 function psi_1energy()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** external functions ****
      real*8   electron_energy
      external electron_energy

      call electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                   dcpl_mb(dng1(1)),
     >                   dbl_mb(rho1_all(1)),
     >                   occupation_on,dbl_mb(occ1(1)))
      psi_1energy = electron_energy(dcpl_mb(psi1(1)),
     >                               dbl_mb(rho1(1)),
     >                              dcpl_mb(dng1(1)),
     >                              dbl_mb(rho1_all(1)),
     >                              occupation_on,dbl_mb(occ1(1)))

      return
      end

*     ***********************************
*     *					*
*     *	    psi_1_noupdate_energy	*
*     *					*
*     ***********************************
      real*8 function psi_1_noupdate_energy()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** external functions ****
      real*8   electron_energy
      external electron_energy


      !call electron_gen_Hpsi_k(dcpl_mb(psi1(1)))
      psi_1_noupdate_energy = electron_energy(dcpl_mb(psi1(1)),
     >                               dbl_mb(rho1(1)),
     >                              dcpl_mb(dng1(1)),
     >                              dbl_mb(rho1_all(1)),
     >                              occupation_on,dbl_mb(occ1(1)))

      return
      end


*     ***********************************
*     *					*
*     *		psi_2energy		*
*     *					*
*     ***********************************
      real*8 function psi_2energy()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** external functions ****
      real*8   electron_energy
      external electron_energy

      call electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                  dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)),
     >                  occupation_on,dbl_mb(occ2(1)))
      psi_2energy = electron_energy(dcpl_mb(psi2(1)),
     >                               dbl_mb(rho2(1)),
     >                              dcpl_mb(dng2(1)),
     >                              dbl_mb(rho2_all(1)),
     >                              occupation_on,dbl_mb(occ2(1)))

      return
      end



*     ***********************************
*     *					*
*     *		psi_1eorbit		*
*     *					*
*     ***********************************
      real*8 function psi_1eorbit()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

*     **** external functions ****
      real*8   electron_eorbit
      external electron_eorbit

      psi_1eorbit = electron_eorbit(dcpl_mb(psi1(1)),
     >                              occupation_on,dbl_mb(occ1(1)))

      return
      end


*     ***********************************
*     *					*
*     *		psi_1ke 		*
*     *					*
*     ***********************************
      real*8 function psi_1ke()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

*     **** local variables ***
      real*8 ave

      call ke_ave(ispin,neq,dcpl_mb(psi1(1)),ave,
     >            occupation_on,dbl_mb(occ1(1)))

      psi_1ke = ave
      return
      end




*     ***********************************
*     *                                 *
*     *         psi_1ke_atom            *
*     *                                 *
*     ***********************************
      real*8 function psi_1ke_atom()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

*     **** local variables ***
      real*8 ave

*     **** external functions ****
      real*8   psp_kinetic_atom
      external psp_kinetic_atom

      ave =  psp_kinetic_atom(ispin,neq,dcpl_mb(psi1(1)))
      call D1dB_SumAll(ave)

      psi_1ke_atom = ave
      return
      end


*     ***********************************
*     *					*
*     *		psi_1vl 		*
*     *					*
*     ***********************************
      real*8 function psi_1vl()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"


*     **** external functions ****
      real*8   electron_psi_vl_ave
      external electron_psi_vl_ave
 
      psi_1vl = electron_psi_vl_ave(dcpl_mb(psi1(1)),dbl_mb(rho1(1)))

      return
      end

*     ***********************************
*     *					*
*     *		psi_1vnl 		*
*     *					*
*     ***********************************
      real*8 function psi_1vnl()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"


*     **** external functions ****
      real*8   electron_psi_vnl_ave
      external electron_psi_vnl_ave
 
      psi_1vnl = electron_psi_vnl_ave(dcpl_mb(psi1(1)),
     >                   occupation_on,dbl_mb(occ1(1)))

      return
      end

*     *******************************
*     *				    *
*     *		psi_1v_field 	    *
*     *				    *
*     *******************************
      real*8 function psi_1v_field()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"



*     **** external functions ****
      real*8   electron_psi_v_field_ave
      external electron_psi_v_field_ave
 
      psi_1v_field = electron_psi_v_field_ave(dcpl_mb(psi1(1)),
     >                                        dbl_mb(rho1(1)))

      return
      end


*     ***********************************
*     *					*
*     *		rho_1Fcharge		*
*     *					*
*     ***********************************
      subroutine rho_1Fcharge(Fcharge)
      implicit none
      real*8 Fcharge(*)

#include "mafdecls.fh"
#include "psi.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer n2ft3d,nx,ny,nz
      integer r_grid(2),rho(2)
      real*8  dv

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

*     **** Initializationsr ****
      call D3dB_n2ft3d(1,n2ft3d) 
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      dv = lattice_omega()/dble(nx*ny*nz)


*     **** Push memory ****
      value = MA_push_get(mt_dbl,(3*n2ft3d),'r_grid',
     >                       r_grid(2),r_grid(1))
      value = value.and. 
     >        MA_push_get(mt_dbl,(3*n2ft3d),'rho',
     >                       rho(2),rho(1))
      if (.not. value) call errquit(
     >     'rho_1Fcharge: out of stack memory',0, MA_ERR)


*     **** Get r_grid and rho ****
      call lattice_r_grid(dbl_mb(r_grid(1)))
      call D3dB_rr_Sum(1,dbl_mb(rho1(1)),
     >                   dbl_mb(rho1(1)+(ispin-1)*n2ft3d),
     >                   dbl_mb(rho(1)))

*     **** Now calculate Fcharge ****
      call pspw_charge_rho_Fcharge(n2ft3d,dbl_mb(r_grid(1)),
     >                            dbl_mb(rho(1)),
     >                            dv,Fcharge)

*     **** Pop memory ****
      value =           MA_pop_stack(rho(2))
      value = value.and.MA_pop_stack(r_grid(2))
      if (.not. value) call errquit(
     >     'rho_1Fcharge: error popping stack memory',0, MA_ERR)

      return
      end



*     ***********************************
*     *					*
*     *		rho_1exc		*
*     *					*
*     ***********************************
      real*8 function rho_1exc()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** external functions ****
      real*8   electron_exc
      external electron_exc

      rho_1exc = electron_exc(dbl_mb(rho1_all(1)))
      return
      end

*     ***********************************
*     *					*
*     *		rho_1pxc		*
*     *					*
*     ***********************************
      real*8 function rho_1pxc()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

*     **** external functions ****
      real*8   electron_pxc
      external electron_pxc

      rho_1pxc = electron_pxc(dbl_mb(rho1(1)))
      return
      end


*     ***********************************
*     *					*
*     *		dng_1ehartree           *
*     *					*
*     ***********************************
      real*8 function dng_1ehartree()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

*     **** external functions ****
      integer  control_version
      real*8   electron_ehartree,electron_ehartree2
      external control_version
      external electron_ehartree,electron_ehartree2

*     **** local variables *****
      real*8 eh

      eh = 0.0d0
      if (control_version().eq.3) 
     >    eh = electron_ehartree(dcpl_mb(dng1(1)))

      if (control_version().eq.4) 
     >    eh = electron_ehartree2(dbl_mb(rho1(1)))

      dng_1ehartree = eh
      return
      end



*     ***********************************
*     *					*
*     *		psi_2toelectron		*
*     *					*
*     ***********************************
      subroutine psi_2toelectron()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

      call electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                   dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)),
     >                   occupation_on,dbl_mb(occ2(1)))

      return
      end


*     ***********************************
*     *					*
*     *		psi_1check_Tangent      *
*     *					*
*     ***********************************
*
*   This routine checks the accuracy of the tangent vector.
*   MM = Yt*H = Yt*(I-Y*Yt)*G = Yt*G - Yt*Y*Yt*G = Yt*G - Yt*G == 0

*     Updated - 5-18-2002
*
      subroutine psi_1check_Tangent(H)
      implicit none
      complex*16 H(*)

#include "errquit.fh"
#include "mafdecls.fh"
#include "psi.fh"

*     **** local variables ****
      logical value
      integer ms,n,indx,i,j
      integer MM(2)
      real*8 sum

      do ms=1,ispin      
         n = ne(ms)
         if (n.eq.0) go to 101  !*** ferromagnetic check ***
         value = MA_push_get(mt_dbl,n*n,'MM',MM(2),MM(1))
         if (.not. value) 
     >   call errquit('out of stack memory in psi_1check_Tangent',0,
     &       MA_ERR)

         indx = (ms-1)*ne(1)*npack1

*        **** calculate MM = Yt*H ****
         call Grsm_ggm_dot(npack1,n,
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
         write(*,*) "psi_1check_Tangent:",sum
            


         value = MA_pop_stack(MM(2))
         if (.not. value) 
     >    call errquit(
     >         'error popping stack memory in psi_1check_Tangent',0,
     >        MA_ERR)

 101     continue
      end do

      return
      end



*     ***********************************
*     *                                 *
*     *         psi_2check_Tangent      *
*     *                                 *
*     ***********************************
*
*   This routine checks the accuracy of the tangent vector.
*   MM = Yt*H = Yt*(I-Y*Yt)*G = Yt*G - Yt*Y*Yt*G = Yt*G - Yt*G == 0

*     Updated - 5-18-2002
*
      subroutine psi_2check_Tangent(H)
      implicit none
      complex*16 H(*)

#include "mafdecls.fh"
#include "psi.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer ms,n,indx,i,j
      integer MM(2)
      real*8 sum

      do ms=1,ispin
         n = ne(ms)
         if (n.eq.0) go to 101  !*** ferromagnetic check ***
         value = MA_push_get(mt_dbl,n*n,'MM',MM(2),MM(1))
         if (.not. value)
     >   call errquit('out of stack memory in psi_1check_Tangent',0,
     >        MA_ERR)

         indx = (ms-1)*ne(1)*npack1

*        **** calculate MM = Yt*H ****
         call Grsm_ggm_dot(npack1,n,
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
         write(*,*) "psi_2check_Tangent:",sum



         value = MA_pop_stack(MM(2))
         if (.not. value)
     >    call errquit(
     >         'error popping stack memory in psi_2check_Tangent',0,
     &       MA_ERR)

 101     continue
      end do

      return
      end



*     ***********************************
*     *					*
*     *		psi_1get_Tgradient	*
*     *					*
*     ***********************************

*     THpsi = Hpsi - Y*Y^t*Hpsi ! used by Grassman minimizers
*
      subroutine psi_1get_Tgradient(THpsi,Eout)
      implicit none
      complex*16 THpsi(*)
      real*8 Eout

#include "mafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** local variables ****
      integer tmp1(2),i,n
      logical value
 
*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack

      real*8   electron_energy
      external electron_energy

      if (.not.Dneall_m_push_get(0,tmp1))
     >   call errquit('out of stack memory in psi_1get_Tradient',0,
     >       MA_ERR)

      call electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)),
     >                   dbl_mb(rho1_all(1)),
     >                  occupation_on,dbl_mb(occ1(1)))

      Eout =  electron_energy(dcpl_mb(psi1(1)),
     >                        dbl_mb(rho1(1)),
     >                        dcpl_mb(dng1(1)),
     >                        dbl_mb(rho1_all(1)),
     >                        occupation_on,dbl_mb(occ1(1)))

      call electron_gen_hml(dcpl_mb(psi1(1)),
     >                       dbl_mb(tmp1(1)))

      call electron_get_Tgradient(dcpl_mb(psi1(1)),
     >                             dbl_mb(tmp1(1)),
     >                            THpsi)
      
      if (.not.Dneall_m_pop_stack(tmp1))
     > call errquit('psi_1get_Tgradient:error popping stack',1,
     >     MA_ERR)

      return
      end


*     ***********************************
*     *					*
*     *		psi_1get_Gradient	*
*     *					*
*     ***********************************

*     THpsi = Hpsi ! used by Projected Grassman minimizers
*
      subroutine psi_1get_Gradient(THpsi,Eout)
      implicit none
      complex*16 THpsi(*)
      real*8 Eout

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** local variables ****
 
*     **** external functions ****
      real*8   electron_energy
      external electron_energy


      call electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)),
     >                   dbl_mb(rho1_all(1)),
     >                   occupation_on,dbl_mb(occ1(1)))

      Eout =  electron_energy(dcpl_mb(psi1(1)),
     >                        dbl_mb(rho1(1)),
     >                        dcpl_mb(dng1(1)),
     >                        dbl_mb(rho1_all(1)),
     >                        occupation_on,dbl_mb(occ1(1)))

      call electron_get_Gradient(THpsi)

      return
      end


*     ***********************************
*     *					*
*     *		psi_1gen_Tangent	*
*     *					*
*     ***********************************

*     THpsi = Hpsi - Y*Y^t*Hpsi ! used by Grassman minimizers
*
      subroutine psi_1gen_Tangent(THpsi)
      implicit none
      complex*16 THpsi(*)

#include "mafdecls.fh"
#include "psi.fh"
#include "errquit.fh"

*     **** local variables ****
      integer tmp1(2)

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack

      if (.not. Dneall_m_push_get(0,tmp1)) 
     >   call errquit('psi_1gen_Tangent:out of stack memory',0, MA_ERR)

      call electron_gen_psiTangenthml(dcpl_mb(psi1(1)),
     >                                THpsi,
     >                                dbl_mb(tmp1(1)))
      call electron_gen_Tangent(dcpl_mb(psi1(1)),
     >                          dbl_mb(tmp1(1)),
     >                          THpsi)

      if (.not. Dneall_m_pop_stack(tmp1)) 
     > call errquit('error popping stack memory in psi_1get_Tradient',0,
     &       MA_ERR)

      return
      end





*     ***********************************
*     *					*
*     *		psi_2get_Tgradient	*
*     *					*
*     ***********************************
      subroutine psi_2get_Tgradient(option,THpsi,Eout)
      implicit none
      integer    option
      complex*16 THpsi(*)
      real*8     Eout

#include "errquit.fh"
#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     *** local variables ****
      integer tmp1(2)


*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack

      real*8   electron_energy
      external electron_energy

      if (.not.Dneall_m_push_get(0,tmp1))
     >   call errquit('out of stack memory in psi_2get_Tradient',0,
     >       MA_ERR)

      if (option.le.1) then
        call electron_run(dcpl_mb(psi2(1)),
     >                     dbl_mb(rho2(1)),
     >                    dcpl_mb(dng2(1)),
     >                     dbl_mb(rho2_all(1)),
     >                     occupation_on,dbl_mb(occ2(1)))
      end if

      Eout =  electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)),
     >                        dbl_mb(rho2_all(1)),
     >                        occupation_on,dbl_mb(occ2(1)))

      call electron_gen_hml(dcpl_mb(psi2(1)),
     >                       dbl_mb(tmp1(1)))
      call electron_get_Tgradient(dcpl_mb(psi2(1)),
     >                             dbl_mb(tmp1(1)),
     >                             THpsi)
      
      if (.not. Dneall_m_pop_stack(tmp1)) 
     >call errquit('psi_2get_Tgradient:error popping stack',1,MA_ERR)

      return
      end


*     ***********************************
*     *					*
*     *		psi_2get_Gradient	*
*     *					*
*     ***********************************
      subroutine psi_2get_Gradient(option,THpsi,Eout)
      implicit none
      integer    option
      complex*16 THpsi(*)
      real*8     Eout

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     *** local variables ****

*     **** external functions ****
      real*8   electron_energy
      external electron_energy


      if (option.le.1) then
        call electron_run(dcpl_mb(psi2(1)),
     >                     dbl_mb(rho2(1)),
     >                    dcpl_mb(dng2(1)),
     >                     dbl_mb(rho2_all(1)),
     >                     occupation_on,dbl_mb(occ2(1)))
      end if

      Eout =  electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)),
     >                        dbl_mb(rho2_all(1)),
     >                        occupation_on,dbl_mb(occ2(1)))

      call electron_get_Gradient(THpsi)

      return
      end

*     ***********************************
*     *					*
*     *		psi_2gen_Tangent	*
*     *					*
*     ***********************************
      subroutine psi_2gen_Tangent(THpsi)
      implicit none
      complex*16 THpsi(*)

#include "mafdecls.fh"
#include "psi.fh"
#include "errquit.fh"

*     *** local variables ****
      integer tmp1(2)

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack


      if (.not. Dneall_m_push_get(0,tmp1)) 
     >   call errquit('psi_2gen_Tangent: out of stack memory',0,MA_ERR)


      call electron_gen_psiTangenthml(dcpl_mb(psi2(1)),
     >                                THpsi,
     >                                dbl_mb(tmp1(1)))
      call electron_gen_Tangent(dcpl_mb(psi2(1)),
     >                          dbl_mb(tmp1(1)),
     >                          THpsi)
      
      if (.not. Dneall_m_pop_stack(tmp1)) 
     > call errquit('error popping stack memory in psi_1get_Tradient',0,
     &       MA_ERR)

      return
      end




*     ***********************************
*     *					*
*     *		psi_1get_TSgradient	*
*     *					*
*     ***********************************

*     THpsi = Hpsi - Y*Hpsi^t*Y ! used by Stiefel minimizers
*
      subroutine psi_1get_TSgradient(THpsi,Eout)
      implicit none
      complex*16 THpsi(*)
      real*8 Eout

#include "errquit.fh"
#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** local variables ****
      integer tmp1(2)
 
*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack

      real*8   electron_energy
      external electron_energy


      if (.not. Dneall_m_push_get(0,tmp1)) 
     >   call errquit('psi_1get_TSradient:pushing stack',0, MA_ERR)


      call electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)),
     >                   dbl_mb(rho1_all(1)),
     >                   occupation_on,dbl_mb(occ1(1)))

      Eout =  electron_energy(dcpl_mb(psi1(1)),
     >                        dbl_mb(rho1(1)),
     >                        dcpl_mb(dng1(1)),
     >                        dbl_mb(rho1_all(1)),
     >                        occupation_on,dbl_mb(occ1(1)))


      call electron_gen_hmlt(dcpl_mb(psi1(1)),
     >                       dbl_mb(tmp1(1)))
      call electron_get_Tgradient(dcpl_mb(psi1(1)),
     >                             dbl_mb(tmp1(1)),
     >                            THpsi)
      

      if (.not. Dneall_m_pop_stack(tmp1)) 
     > call errquit('psi_1get_TSgradient:popping stack',1, MA_ERR)

      return
      end


*     ***********************************
*     *					*
*     *		psi_2get_TSgradient	*
*     *					*
*     ***********************************

*     THpsi = Hpsi - Y*Hpsi^t*Y ! used by Stiefel minimizers
*
      subroutine psi_2get_TSgradient(option,THpsi,Eout)
      implicit none
      integer    option
      complex*16 THpsi(*)
      real*8     Eout

#include "errquit.fh"
#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     *** local variables ****
      integer tmp1(2)

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack

      real*8   electron_energy
      external electron_energy


      if (.not. Dneall_m_push_get(0,tmp1)) 
     >   call errquit('psi_2get_TSgradient:pushing stack',0, MA_ERR)

      if (option.le.1) then
        call electron_run(dcpl_mb(psi2(1)),
     >                     dbl_mb(rho2(1)),
     >                    dcpl_mb(dng2(1)),
     >                     dbl_mb(rho2_all(1)),
     >                    occupation_on,dbl_mb(occ2(1)))
      end if

      Eout =  electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)),
     >                        dbl_mb(rho2_all(1)),
     >                        occupation_on,dbl_mb(occ2(1)))

      call electron_gen_hmlt(dcpl_mb(psi2(1)),
     >                       dbl_mb(tmp1(1)))
      call electron_get_Tgradient(dcpl_mb(psi2(1)),
     >                             dbl_mb(tmp1(1)),
     >                             THpsi)
      
      if (.not. Dneall_m_pop_stack(tmp1)) 
     > call errquit('psi_2get_TSgradient:popping stack',1, MA_ERR)

      return
      end




*     ***********************************
*     *					*
*     *		psi_1get_TMgradient	*
*     *					*
*     ***********************************
      subroutine psi_1get_TMgradient(THpsi,Eout)
      implicit none
      complex*16 THpsi(*)
      real*8     Eout

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** external functions ****
      real*8   electron_energy
      external electron_energy


      call electron_run(dcpl_mb(psi1(1)),
     >                   dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)),
     >                   dbl_mb(rho1_all(1)),
     >                  occupation_on,dbl_mb(occ1(1)))

      Eout =  electron_energy(dcpl_mb(psi1(1)),
     >                        dbl_mb(rho1(1)),
     >                        dcpl_mb(dng1(1)),
     >                        dbl_mb(rho1_all(1)),
     >                        occupation_on,dbl_mb(occ1(1)))

      call electron_get_TMgradient(dcpl_mb(psi1(1)),
     >                            THpsi)

      return
      end



*     ***********************************
*     *					*
*     *		psi_2get_TMgradient	*
*     *					*
*     ***********************************
      subroutine psi_2get_TMgradient(option,THpsi,Eout)
      implicit none
      integer    option
      complex*16 THpsi(*)
      real*8     Eout

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** external functions ****
      real*8   electron_energy
      external electron_energy

      if (option.le.1) then
        call electron_run(dcpl_mb(psi2(1)),
     >                    dbl_mb(rho2(1)),
     >                    dcpl_mb(dng2(1)),
     >                    dbl_mb(rho2_all(1)),
     >                    occupation_on,dbl_mb(occ2(1)))
      end if

      Eout =  electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)),
     >                        dbl_mb(rho2_all(1)),
     >                        occupation_on,dbl_mb(occ2(1)))

      call electron_get_TMgradient(dcpl_mb(psi2(1)),
     >                             THpsi)
      
      return
      end

*     ***********************************
*     *					*
*     *		psi_1ke_Precondition	*
*     *					*
*     ***********************************
      subroutine psi_1ke_Precondition(Hpsi)
      implicit none
      complex*16 Hpsi(*)

#include "mafdecls.fh"
#include "psi.fh"

*     **** local variables ****
      integer neall

      neall = neq(1)+neq(2)
      call ke_Precondition(npack1,neall,
     >                      dcpl_mb(psi1(1)),
     >                      Hpsi)
      return
      end



*     ***********************************
*     *					*
*     *	    psi_1geodesic_transport	*
*     *					*
*     ***********************************
      subroutine psi_1geodesic_transport(t,H0)
      implicit none
      real*8 t
      complex*16 H0(*)

#include "mafdecls.fh"
#include "psi.fh"


      call geodesic_transport(t,dcpl_mb(psi1(1)),H0)

      return
      end


*     ***********************************
*     *					*
*     *	    psi_1geodesic_Gtransport	*
*     *					*
*     ***********************************
      subroutine psi_1geodesic_Gtransport(t,G0)
      implicit none
      real*8 t
      complex*16 G0(*)

#include "mafdecls.fh"
#include "psi.fh"

      call geodesic_Gtransport(t,dcpl_mb(psi1(1)),G0)

      return
      end



*     ***********************************
*     *					*
*     *		psi_geodesic_energy 	*
*     *					*
*     ***********************************
      real*8 function psi_geodesic_energy(t)
      implicit none
      real*8 t

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

      real*8 e_new
*     **** external functions ****
      real*8   electron_energy
      external electron_energy

     
      call geodesic_get(t,dcpl_mb(psi1(1)),
     >                    dcpl_mb(psi2(1)))
    
      call electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                  dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)),
     >                   occupation_on,dbl_mb(occ2(1)))

      e_new =  electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)),
     >                        dbl_mb(rho2_all(1)),
     >                        occupation_on,dbl_mb(occ2(1)))

      psi_geodesic_energy = e_new
      return
      end

*     ***********************************
*     *					*
*     *		psi_geodesic_denergy 	*
*     *					*
*     ***********************************
      real*8 function psi_geodesic_denergy(t)
      implicit none
      real*8 t

#include "mafdecls.fh"
#include "psi.fh"

*     **** external functions ****
      real*8   electron_eorbit_noocc
      external electron_eorbit_noocc


      call geodesic_transport(t,dcpl_mb(psi1(1)),
     >                          dcpl_mb(psi2(1)))
      psi_geodesic_denergy 
     > =  2.0d0*electron_eorbit_noocc(dcpl_mb(psi2(1)))

      return
      end

*     ***********************************
*     *					*
*     *		psi_geodesic_final 	*
*     *					*
*     ***********************************
      subroutine psi_geodesic_final(t)
      implicit none
      real*8 t

#include "mafdecls.fh"
#include "psi.fh"

      integer taskid,MASTER
      parameter (MASTER=0)
c     real*8 sum1,sum2
     
      call Parallel_taskid(taskid)

      call geodesic_get(t,dcpl_mb(psi1(1)),
     >                    dcpl_mb(psi2(1)))          
      return
      end



*     ***********************************
*     *					*
*     *	    psi_1geodesic2_start	*
*     *					*
*     ***********************************
      subroutine psi_1geodesic2_start(H0,max_sigma,dE0)
      implicit none
      complex*16 H0(*)
      real*8 max_sigma
      real*8 dE0

#include "mafdecls.fh"
#include "psi.fh"

      call geodesic2_start(dcpl_mb(psi1(1)),H0,max_sigma,dE0)

      return
      end

*     ***********************************
*     *					*
*     *	    psi_1geodesic2_transport	*
*     *					*
*     ***********************************
      subroutine psi_1geodesic2_transport(t,Hnew)
      implicit none
      real*8 t
      complex*16 Hnew(*)

#include "mafdecls.fh"
#include "psi.fh"

      call geodesic2_transport(t,dcpl_mb(psi1(1)),Hnew)

      return
      end



*     ***********************************
*     *					*
*     *		psi_geodesic2_energy 	*
*     *					*
*     ***********************************
      real*8 function psi_geodesic2_energy(t)
      implicit none
      real*8 t

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

      real*8 e_new

*     **** external functions ****
      real*8   electron_energy
      external electron_energy

     
      call geodesic2_get(t,dcpl_mb(psi1(1)),
     >                    dcpl_mb(psi2(1)))

*     **** check Orthogonality of psi2 **** !debug
*      call OrthoCheck_geo(ispin,ne,dcpl_mb(psi2(1))) !debug


      call electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                  dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)),
     >                   occupation_on,dbl_mb(occ2(1)))
      e_new =  electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)),
     >                        dbl_mb(rho2_all(1)),
     >                        occupation_on,dbl_mb(occ2(1)))

      psi_geodesic2_energy = e_new
      return
      end

*     ***********************************
*     *					*
*     *		psi_geodesic2_denergy 	*
*     *					*
*     ***********************************
      real*8 function psi_geodesic2_denergy(t)
      implicit none
      real*8 t

#include "mafdecls.fh"
#include "psi.fh"

*     **** external functions ****
      real*8   electron_eorbit
      external electron_eorbit


      call geodesic2_transport(t,dcpl_mb(psi1(1)),
     >                          dcpl_mb(psi2(1)))
      psi_geodesic2_denergy =  2.0d0*electron_eorbit(dcpl_mb(psi2(1)),
     >                                  occupation_on,dbl_mb(occ2(1)))

      return
      end

*     ***********************************
*     *					*
*     *		psi_geodesic2_final 	*
*     *					*
*     ***********************************
      subroutine psi_geodesic2_final(t)
      implicit none
      real*8 t

#include "mafdecls.fh"
#include "psi.fh"

      integer taskid,MASTER
      parameter (MASTER=0)
c     real*8 sum1,sum2
     
      call Parallel_taskid(taskid)

      call geodesic2_get(t,dcpl_mb(psi1(1)),
     >                    dcpl_mb(psi2(1)))
      return
      end



*     ***********************************
*     *					*
*     *		psito2_sd_update	*
*     *					*
*     ***********************************
      subroutine psi1to2_sd_update(dte)
      implicit none
      real*8 dte

#include "mafdecls.fh"
#include "psi.fh"
#include "errquit.fh"


*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** local variables ****
      logical value
      integer nemaxq,ierr
      integer lmd(2),tmp_L(2)

*     **** external functions ****
      logical  pspw_SIC,Dneall_m_push_get,Dneall_m_push_get_block
      logical  Dneall_m_pop_stack
      external pspw_SIC,Dneall_m_push_get,Dneall_m_push_get_block
      external Dneall_m_pop_stack

      call electron_run(dcpl_mb(psi1(1)),
     >                  dbl_mb(rho1(1)),
     >                  dcpl_mb(dng1(1)),
     >                  dbl_mb(rho1_all(1)),
     >                  occupation_on,dbl_mb(occ1(1)))

*     **** do a steepest descent step ****
      call electron_sd_update(dcpl_mb(psi1(1)),
     >                        dcpl_mb(psi2(1)),
     >			      dte)

*     **** lagrange multiplier corrections ****
      nemaxq = neq(1)+neq(2)

*     **** allocate MA local variables ****
      value =           Dneall_m_push_get_block(1,8,tmp_L)
      value = value.and.Dneall_m_push_get(0,lmd)

c        if (occupation_on) then
c        call psi_lmbda2(ispin,ne,nemaxq,npack1,
c     >                 dcpl_mb(psi1(1)),dcpl_mb(psi2(1)),
c     >                 dte,dbl_mb(occ1(1)),
c     >                 dbl_mb(lmd(1)),
c     >                 dbl_mb(tmp_L(1)),ierr)
        if (pspw_SIC().or.occupation_on) then
        call psi_lmbda_sic(ispin,ne,nemaxq,npack1,
     >                 dcpl_mb(psi1(1)),dcpl_mb(psi2(1)),dte,
     >                 dbl_mb(lmd(1)),
     >                 dbl_mb(tmp_L(1)),ierr)
        else
        call psi_lmbda(ispin,ne,nemaxq,npack1,
     >                 dcpl_mb(psi1(1)),dcpl_mb(psi2(1)),dte,
     >                 dbl_mb(lmd(1)),
     >                 dbl_mb(tmp_L(1)),ierr)
        end if


      value = value.and.Dneall_m_pop_stack(lmd)
      value = value.and.Dneall_m_pop_stack(tmp_L)
      if (.not. value)
     >     call errquit(
     >          'psi1to2_sd_update:stack failure', 0, MA_ERR)
      return
      end


*     ***********************************
*     *					*
*     *		psi_1force              *
*     *					*
*     ***********************************
      subroutine psi_1force(fion)
      implicit none
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     **** local variables ****
      logical value
      integer r_grid(2),tmp(2)

*     **** external functions ****
      integer  control_version
      external control_version

c     call electron_gen_psi_r(dcpl_mb(psi1(1)))
c     call electron_gen_densities(dcpl_mb(psi1(1)),
c    >                             dbl_mb(rho1(1)),
c    >                            dcpl_mb(dng1(1)))

      call f_vlocal(dcpl_mb(dng1(1)),fion)

      if (control_version().eq.4) then
          value = MA_push_get(mt_dbl,(2*nfft3d),'tmp',
     >                        tmp(2),tmp(1))
          value = value.and.
     >            MA_push_get(mt_dbl,(6*nfft3d),'r_grid',
     >                        r_grid(2),r_grid(1))
         if (.not. value) call errquit('out of stack memory',0, MA_ERR)

          call D3dB_rr_Sum(1,dbl_mb(rho1(1)),
     >                       dbl_mb(rho1(1)+(ispin-1)*2*nfft3d),
     >                       dbl_mb(tmp(1)))
          call lattice_r_grid(dbl_mb(r_grid(1)))
          call grad_v_lr_local(dbl_mb(r_grid(1)),
     >                         dbl_mb(tmp(1)),
     >                         fion)

          value = MA_pop_stack(r_grid(2))
          value = value.and.MA_pop_stack(tmp(2))
         if (.not. value) call errquit('error popping stack memory',0,
     &       MA_ERR)
      end if

      call f_vnonlocal(ispin,
     >                 neq,
     >                 dcpl_mb(psi1(1)),
     >                 fion,
     >                 occupation_on,dbl_mb(occ1(1)))
      return
      end

*     ***********************************
*     *                                 *
*     *         psi_1force_local        *
*     *                                 *
*     ***********************************
      subroutine psi_1force_local(fion)
      implicit none
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     **** local variables ****
      logical value
      integer r_grid(2),tmp(2)

*     **** external functions ****
      integer  control_version
      external control_version

c     call electron_gen_psi_r(dcpl_mb(psi1(1)))
c     call electron_gen_densities(dcpl_mb(psi1(1)),
c    >                             dbl_mb(rho1(1)),
c    >                            dcpl_mb(dng1(1)))

      call f_vlocal(dcpl_mb(dng1(1)),fion)

      if (control_version().eq.4) then
          value = MA_push_get(mt_dbl,(2*nfft3d),'tmp',
     >                        tmp(2),tmp(1))
          value = value.and.
     >            MA_push_get(mt_dbl,(6*nfft3d),'r_grid',
     >                        r_grid(2),r_grid(1))
         if (.not. value) call errquit('out of stack memory',0, MA_ERR)

          call D3dB_rr_Sum(1,dbl_mb(rho1(1)),
     >                       dbl_mb(rho1(1)+(ispin-1)*2*nfft3d),
     >                       dbl_mb(tmp(1)))
          call lattice_r_grid(dbl_mb(r_grid(1)))
          call grad_v_lr_local(dbl_mb(r_grid(1)),
     >                         dbl_mb(tmp(1)),
     >                         fion)

          value = MA_pop_stack(r_grid(2))
          value = value.and.MA_pop_stack(tmp(2))
         if (.not. value) call errquit('error popping stack memory',0,
     &       MA_ERR)
      end if

      return
      end

*     ***********************************
*     *                                 *
*     *         psi_1force_nonlocal     *
*     *                                 *
*     ***********************************
      subroutine psi_1force_nonlocal(fion)
      implicit none
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

      call f_vnonlocal(ispin,
     >                 neq,
     >                 dcpl_mb(psi1(1)),
     >                 fion,
     >                 occupation_on,dbl_mb(occ1(1)))
      return
      end






*     ***********************************
*     *					*
*     *		psi_1ke_stress          *
*     *					*
*     ***********************************
      subroutine psi_1ke_stress(stress)
      implicit none
      real*8 stress(3,3)

#include "mafdecls.fh"
#include "psi.fh"

      call ke_euv(ispin,neq,dcpl_mb(psi1(1)),stress)
      return
      end

*     ***********************************
*     *					*
*     *		psi_1coulomb_stress     *
*     *					*
*     ***********************************
      subroutine psi_1coulomb_stress(stress)
      implicit none
      real*8 stress(3,3)

#include "mafdecls.fh"
#include "psi.fh"

      call coulomb_euv(dcpl_mb(dng1(1)),stress)
      return
      end

*     ***********************************
*     *					*
*     *		rho_1exc_stress 	*
*     *					*
*     ***********************************
      subroutine rho_1exc_stress(stress)
      implicit none
      real*8 stress(3,3)

#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     ***** local variables ****
      integer u,v,gga
      real*8 exc,pxc
      real*8 pi,scal,hm(3,3),tstress(3,3)

*     **** external functions ****
      integer  control_gga
      real*8   rho_1exc,rho_1pxc,lattice_unitg,lattice_omega
      external control_gga
      external rho_1exc,rho_1pxc,lattice_unitg,lattice_omega

*     *** define hm ****
      pi   = 4.0d0*datan(1.0d0)
      scal = 1.0d0/(2.0d0*pi)
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
      end do
      end do

*     **** LDA part ****
      exc = rho_1exc()
      pxc = rho_1pxc()
      do v=1,3
      do u=1,3
         stress(u,v) = (exc-pxc)*hm(u,v)
      end do
      end do
      !write(*,*) "hm(1,1):",hm(1,1),1.0d0/hm(1,1)
      !write(*,*) "exc:",exc,pxc
      !write(*,*) "D:",stress(1,1)

*     **** PBE96 GGA part ****
*     **** finished? 11/24/04 - still need to test ***
      gga = control_gga()
      if ((gga.ge.10).and.(gga.lt.100)) then
       call v_bwexc_euv(gga,2*nfft3d,ispin,dbl_mb(rho1_all(1)),
     >                  1.0d0,1.0d0,tstress)
       do v=1,3
       do u=1,3
          stress(u,v) = stress(u,v) + tstress(u,v)
       end do
       end do
      end if 

      if (gga.eq.110) then
       call v_bwexc_euv(10,2*nfft3d,ispin,dbl_mb(rho1_all(1)),
     >                  0.75d0,1.0d0,tstress)
       do v=1,3
       do u=1,3
          stress(u,v) = stress(u,v) + tstress(u,v)
       end do
       end do
      end if 

      if (gga.eq.112) then
       call v_bwexc_euv(12,2*nfft3d,ispin,dbl_mb(rho1_all(1)),
     >                  0.75d0,1.0d0,tstress)
       do v=1,3
       do u=1,3
          stress(u,v) = stress(u,v) + tstress(u,v)
       end do
       end do
      end if 

      return
      end

*     ***********************************
*     *					*
*     *		rho_1semicore_stress 	*
*     *					*
*     ***********************************
      subroutine rho_1semicore_stress(stress)
      implicit none
      real*8 stress(3,3)

#include "mafdecls.fh"
#include "psi.fh"

*     **** not finished ****
      call semicore_euv(stress)

      return
      end




*     ***********************************
*     *					*
*     *		dng_1vlocal_stress      *
*     *					*
*     ***********************************

      subroutine dng_1vlocal_stress(stress)
      implicit none
      real*8 stress(3,3)


#include "mafdecls.fh"
#include "psi.fh"

      call v_local_euv(dcpl_mb(dng1(1)),stress)

      return
      end

*     ***********************************
*     *					*
*     *		psi_1vnonlocal_stress   *
*     *					*
*     ***********************************
      subroutine psi_1vnonlocal_stress(stress)
      implicit none
      real*8 stress(3,3)

#include "mafdecls.fh"
#include "psi.fh"



*     ***** local variables ****
      integer u,v
      real*8 evnl
      real*8 pi,scal,hm(3,3)

*     **** external functions ****
      real*8   psi_1vnl,lattice_unitg
      external psi_1vnl,lattice_unitg

*     *** define hm ****
      pi   = 4.0d0*datan(1.0d0)      
      scal = 1.0d0/(2.0d0*pi)
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
      end do
      end do

      call v_nonlocal_euv_2(ispin,neq,dcpl_mb(psi1(1)),stress)
      evnl = psi_1vnl()

      do v=1,3
      do u=1,3
         stress(u,v) = stress(u,v) - evnl*hm(u,v)
      end do
      end do

      return
      end




*     ***********************************
*     *					*
*     *		psi_1Orb_Analysis       *
*     *					*
*     ***********************************
      subroutine psi_1Orb_Analysis(iunit)
      implicit none
      integer iunit

#include "mafdecls.fh"
#include "psi.fh"

c      call Orb_Analysis(iunit,ispin,ne,dcpl_mb(psi1(1)))
      return
      end

*     ***********************************
*     *					*
*     *		psi_1Shml 	      	*
*     *					*
*     ***********************************
      subroutine psi_1Shml(S0,S0hml)
      implicit none
      complex*16 S0(*)
      complex*16 S0hml(*)

#include "mafdecls.fh"
#include "psi.fh"

      integer ms,n,shift1,shift2

      call electron_gen_hml(dcpl_mb(psi1(1)),dbl_mb(hml(1)))
      do ms=1,ispin
            n     = ne(ms)
            if (n.le.0) go to 30
            shift1 = 1 + (ms-1)*ne(1)*npack1
            shift2 =     (ms-1)*ne(1)*ne(1)
            call DGEMM('N','N',2*npack1,n,n,
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
*     *					*
*     *		psi_1gen_hml      	*
*     *					*
*     ***********************************
      subroutine psi_1gen_hml()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"


      call electron_gen_hml(dcpl_mb(psi1(1)),dbl_mb(hml(1)))

      return
      end




*     ***********************************
*     *                                 *
*     *         psi_1gen_hml_g          *
*     *                                 *
*     ***********************************
      subroutine psi_1gen_hml_g()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"


      call electron_gen_hml_g(dcpl_mb(psi1(1)),dbl_mb(hml(1)))

      return
      end


*     ***********************************
*     *					*
*     *		psi_2gen_hml      	*
*     *					*
*     ***********************************
      subroutine psi_2gen_hml()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"


      call electron_gen_hml(dcpl_mb(psi2(1)),dbl_mb(hml(1)))

      return
      end

*     ***********************************
*     *					*
*     *		psi_eigenvalue    	*
*     *					*
*     ***********************************
      real*8  function psi_eigenvalue(ms,i)
      implicit none
      integer ms
      integer i

#include "mafdecls.fh"
#include "psi.fh"

      real*8 sum

      sum = dbl_mb(eig(1)+(i-1)+(ms-1)*ne(1))
      psi_eigenvalue = sum

      return
      end

*     ***********************************
*     *                                 *
*     *         psi_occupation          *
*     *                                 *
*     ***********************************
      real*8  function psi_occupation(ms,i)
      implicit none
      integer ms
      integer i

#include "mafdecls.fh"
#include "psi.fh"

      if (occupation_on) then
         psi_occupation = dbl_mb(occ1(1)+(i-1)+(ms-1)*ne(1))
      else
         psi_occupation = 1.0d0
      end if
      return
      end

*     ***********************************
*     *                                 *
*     *     psi_1reverse_occupation     *
*     *                                 *
*     ***********************************
      subroutine psi_1reverse_occupation()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

      integer ms,i,indx1,indx2

      do ms=1,ispin
         indx1 = occ1(1) + ne(ms) - 1 + (ms-1)*ne(1)
         indx2 = occ2(1)              + (ms-1)*ne(1)
         do i=1,ne(ms)
            dbl_mb(indx2)=dbl_mb(indx1)
            indx1 = indx1 - 1
            indx2 = indx2 + 1
         end do
      end do
      call dcopy((ne(1)+ne(2)),dbl_mb(occ2(1)),1,dbl_mb(occ1(1)),1)
      return
      end


*     ***********************************
*     *                                 *
*     *      psi_1define_occupation    *
*     *                                 *
*     ***********************************
      subroutine psi_1define_occupation(initial)
      implicit none
      logical initial

#include "mafdecls.fh"
#include "psi.fh"

*     **** local variables ****
      integer it,itmax
      parameter (itmax=50)

      integer ms,nb,n,shift1,shift2,occ1_tag,ndiff
      real*8 e,x,kT,f,g,alpha,pi
      real*8 ZZ,Z(2),Zlower,Zmid,Zupper,elower,emid,eupper
      real*8 flower,fmid,fupper

*     **** external functions ****
      integer  control_multiplicity
      real*8   control_TotalCharge,ion_TotalCharge_qm
      real*8   psi_occ_distribution,control_ks_alpha
      external control_multiplicity
      external control_TotalCharge,ion_TotalCharge_qm
      external psi_occ_distribution,control_ks_alpha


      smearfermi(1)   = 0.0d0
      smearfermi(2)   = 0.0d0
      smearcorrection = 0.0d0

      if (occupation_on) then
         alpha = control_ks_alpha()
         kT    = smearkT
         ZZ  = ion_TotalCharge_qm() - control_TotalCharge()
         if (ispin.eq.2) then
            ndiff = control_multiplicity() - 1
            Z(1) = 0.5d0*(ZZ+ndiff)
            Z(2) = 0.5d0*(ZZ-ndiff)
         else
            Z(1) = 0.5d0*ZZ
            Z(2) = 0.0d0
         end if
       
         pi    = 4.0d0*datan(1.0d0)
         if (initial) alpha = 1.0d0

*        **** outer loop over spins ****
         smearcorrection = 0.0d0
         do ms=1,ispin

*           **** find eupper and elower ****
            elower =  9.9d12
            eupper = -9.9d12
            shift1 = eig(1) + (ms-1)*ne(1)
            do n=1,ne(ms)
              e       = dbl_mb(shift1)
              if (e.lt.elower) elower = e
              if (e.gt.eupper) eupper = e
              shift1  = shift1 + 1
            end do


*           **** find fermi level ****
            Zlower = 0.0d0
            Zupper = 0.0d0
            shift1 = eig(1) + (ms-1)*ne(1)
            do n=1,ne(ms)
              e = dbl_mb(shift1)
              Zlower = Zlower 
     >         + psi_occ_distribution(smeartype,(e-elower)/kT)
              Zupper = Zupper 
     >         + psi_occ_distribution(smeartype,(e-eupper)/kT)
              shift1  = shift1 + 1
            end do

            flower = Zlower - Z(ms)
            fupper = Zupper - Z(ms)

            if (flower*fupper.ge.0.0d0) 
     >       call errquit(
     >            'psi_1define_occupation:Fermi energy not found',ms,0)

            it = 0
  20        it = it + 1
            emid = 0.5d0*(elower + eupper)
            Zmid = 0.0d0
            shift1 = eig(1) + (ms-1)*ne(1)
            do n=1,ne(ms)
              e = dbl_mb(shift1)
              Zmid = Zmid + psi_occ_distribution(smeartype,(e-emid)/kT)
              shift1  = shift1 + 1
            end do
            fmid = Zmid - Z(ms)
            if (fmid.lt.0.0d0) then
               flower = fmid
               elower = emid
            else 
               fupper = fmid
               eupper = emid
            end if
            if ( (dabs(fmid)     .gt.1.0d-11) .and.
     >           ((eupper-elower).gt.1.0d-11) .and.
     >           (it.lt.itmax))   goto 20
   

            smearfermi(ms) = emid

*           **** determine filling and correction ****
            shift1 = eig(1)  + (ms-1)*ne(1)
            shift2 = occ1(1) + (ms-1)*ne(1)
            do n=1,ne(ms)
              e = dbl_mb(shift1)
              x = (e - smearfermi(ms))/kT
              f = psi_occ_distribution(smeartype,x)


              dbl_mb(shift2) = (1.0d0-alpha)*dbl_mb(shift2) + alpha*f

              if (smeartype.eq.1) then
                 if (  (dbl_mb(shift2)       .gt.1.0d-6) .and. 
     >               ( (1.0d0-dbl_mb(shift2)).gt.1.0d-6)) then
                smearcorrection = smearcorrection  
     >           + kT*( dbl_mb(shift2)*log(dbl_mb(shift2)) 
     >           + (1.0d0-dbl_mb(shift2))*log(1.0d0-dbl_mb(shift2)) )
                 end if
              else if (smeartype.eq.2) then
                smearcorrection 
     >              = smearcorrection 
     >              - kT*dexp(-x*x)/(4.0d0*dsqrt(pi))
              end if
              shift1  = shift1 + 1
              shift2  = shift2 + 1
            end do

         end do !** ms***
         if (ms.eq.1) smearcorrection=smearcorrection+smearcorrection

      end if

      return
      end

c  set nwpw:fractional_smeartype 1 #0-none, 1-Fermi-Dirac, 2-Gaussian, 3-Hermite

      real*8 function psi_occ_distribution(smeartype,e)
      implicit none
      integer smeartype
      real*8 e
      real*8 f

*     **** external functions ****
      real*8   util_erfc
      external util_erfc

      if (smeartype.eq.1) then
         if (e.gt.30.0d0) then
           f = 0.0d0
         else if (e.lt.(-30.0d0)) then
           f = 1.0d0
         else
           f = 1.0d0/(1.0d0+dexp(e))
         end if
      else if (smeartype.eq.2) then
         f = 0.5d0*util_erfc(e)
      else 
         if (e.gt.0.0d0) then
           f = 0.0d0
         else
           f = 1.0d0
         end if
      end if
      psi_occ_distribution = f 
      return
      end

      real*8 function psi_smearfermi(ms)
      implicit none
      integer ms
#include "psi.fh"
      psi_smearfermi = smearfermi(ms)
      return
      end
      real*8 function psi_smearcorrection()
      implicit none
#include "psi.fh"
      psi_smearcorrection = smearcorrection
      return
      end





*     ***********************************
*     *                                 *
*     *          psi_virtual            *
*     *                                 *
*     ***********************************
      real*8  function psi_virtual(ms,i)
      implicit none
      integer ms
      integer i

#include "mafdecls.fh"
#include "psi.fh"

      psi_virtual=dbl_mb(eig_excited(1)+(i-1)+(ms-1)*ne_excited(1))

      return
      end

*     ***********************************
*     *					*
*     *		psi_hml		   	*
*     *					*
*     ***********************************
      real*8  function psi_hml(ms,i,j)
      implicit none
      integer ms
      integer i,j

#include "mafdecls.fh"
#include "psi.fh"

      psi_hml = dbl_mb(hml(1)-1 + i 
     >                          + (j-1)*ne(ms) 
     >                          + (ms-1)*ne(1)*ne(1)) 

      return
      end


*     ***********************************
*     *                                 *
*     *         psi_iptr_hml            *
*     *                                 *
*     ***********************************
      integer function psi_iptr_hml(ms,i,j)
      implicit none
      integer ms
      integer i,j

#include "mafdecls.fh"
#include "psi.fh"

      psi_iptr_hml = (hml(1)-1 + i
     >                          + (j-1)*ne(ms)
     >                          + (ms-1)*ne(1)*ne(1))

      return
      end


*     ***********************************
*     *					*
*     *		psi_spin_density  	*
*     *					*
*     ***********************************
      subroutine psi_spin_density(en)
      implicit none
      real*8 en(2)

#include "mafdecls.fh"
#include "psi.fh"

*     **** local variables ****
      integer ms,nx,ny,nz,n2ft3d
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
         en(ms) = sumall*scale
      end do
      
      return
      end

*     ***********************************
*     *					*
*     *		psi_spin2     	        *
*     *					*
*     ***********************************
      subroutine psi_spin2(Sab)
      implicit none
      real*8 Sab

#include "mafdecls.fh"
#include "psi.fh"

      call Calculate_psi_spin2(ispin,ne,npack1,dcpl_mb(psi1(1)),
     >                         occupation_on,dbl_mb(occ1(1)),Sab)
      return
      end

*     ***********************************
*     *					*
*     *		psi_1rotate2       	*
*     *					*
*     ***********************************
      subroutine psi_1rotate2()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

c*     ***** local variables *****
c      integer ms,index,i,j,shift1,shift2

      call Dneall_fmf_Multiply(0,dcpl_mb(psi1(1)),npack1,
     >                          dbl_mb(hml(1)),1.0d0,
     >                          dcpl_mb(psi2(1)),0.0d0)



c      !call dcopy(2*npack1*(ne(1)+ne(2)),0.0d0,0,dcpl_mb(psi2(1)),1)
c      do ms=1,ispin
c         if (ne(ms).le.0) go to 30
c         shift1 = (ms-1)*ne(1)
c         shift2 = (ms-1)*ne(1)*ne(1)
c
c         call DGEMM('N','N',2*npack1,ne(ms),ne(ms),
c     >              (1.0d0),
c     >              dcpl_mb(psi1(1)+shift1*npack1),2*npack1,
c     >              dbl_mb(hml(1)+shift2),ne(ms),
c     >              (0.0d0),
c     >              dcpl_mb(psi2(1)+shift1*npack1),2*npack1)
cc        do j=1,ne(ms)
cc          do i=1,ne(ms)
cc             index = (i-1) + (j-1)*ne(ms) + shift2
cc            
cc              call D3dB_cc_daxpy(1,dbl_mb(hml(1)+index),
cc     >                           dcpl_mb(psi1(1)+(i-1+shift1)*nfft3d),
cc     >                           dcpl_mb(psi2(1)+(j-1+shift1)*nfft3d)) 
cc             call Pack_cc_daxpy(1,dbl_mb(hml(1)+index),
cc    >                           dcpl_mb(psi1(1)+(i-1+shift1)*npack1),
cc    >                           dcpl_mb(psi2(1)+(j-1+shift1)*npack1)) 
cc          end do
cc        end do
c
c   30   continue
c      end do

      return
      end

*     ***********************************
*     *					*
*     *		psi_2rotate1       	*
*     *					*
*     ***********************************
      subroutine psi_2rotate1()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

c*     ***** local variables *****
c      integer ms,index,i,j,shift1,shift2

      call Dneall_fmf_Multiply(0,dcpl_mb(psi2(1)),npack1,
     >                          dbl_mb(hml(1)),1.0d0,
     >                          dcpl_mb(psi1(1)),0.0d0)

c      do ms=1,ispin
c         if (ne(ms).le.0) go to 30
c         shift1 = (ms-1)*ne(1)
c         shift2 = (ms-1)*ne(1)*ne(1)
c
c         call DGEMM('N','N',2*npack1,ne(ms),ne(ms),
c     >              (1.0d0),
c     >              dcpl_mb(psi2(1)+shift1*npack1),2*npack1,
c     >              dbl_mb(hml(1)+shift2),ne(ms),
c     >              (0.0d0),
c     >              dcpl_mb(psi1(1)+shift1*npack1),2*npack1)
c
c   30    continue
c      end do

      return
      end


*     ***********************************
*     *					*
*     *		psi_diagonalize_hml	*
*     *					*
*     ***********************************
      subroutine psi_diagonalize_hml()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"


c*     ***** local variables ****
c      logical value
c      integer ms,shift1,shift2,ierr,i,j,indx
c      integer tmp1(2)


      call Dneall_m_diagonalize(0,dbl_mb(hml(1)),
     >                             dbl_mb(eig(1)),.false.)

c      value = MA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
c      if (.not. value) 
c     >   call errquit('out of stack memory in psi_diagonalize_hml',0,
c     &       MA_ERR)


c*     ***** diagonalize the hamiltonian matrix *****
c      call dcopy((ne(1)+ne(2)),0.0d0,0,dbl_mb(eig(1)),1)
c      do ms=1,ispin
c         shift1 = (ms-1)*ne(1)
c         shift2 = (ms-1)*ne(1)*ne(1)
c         if (ne(ms).le.0) go to 30

cc        call eigen(ne(ms),ne(ms),
cc    >              dbl_mb(hml(1)+shift2),
cc    >              dbl_mb(eig(1)+shift1),
cc    >              dbl_mb(tmp1(1)),ierr)
c
c         call DSYEV('V','U',ne(ms),
c     >              dbl_mb(hml(1)+shift2),ne(ms), 
c     >              dbl_mb(eig(1)+shift1),
c     >              dbl_mb(tmp1(1)),2*ne(1)*ne(1),
c     >              ierr)
c
c         call eigsrt(dbl_mb(eig(1)+shift1),
c     >              dbl_mb(hml(1)+shift2),
c     >              ne(ms),ne(ms))
c
c  30    continue
c      end do
c
c      
c      value = MA_pop_stack(tmp1(2))
c      if (.not. value) 
c     > call errquit('error popping stack in psi_diagonalize_hml',0,
c     &       MA_ERR)

      return
      end

*     ***********************************
*     *					*
*     *	  psi_diagonalize_hml_assending *
*     *					*
*     ***********************************
      subroutine psi_diagonalize_hml_assending()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"


c*     ***** local variables ****
c      logical value
c      integer ms,shift1,shift2,ierr
c      integer tmp1(2)

      call Dneall_m_diagonalize(0,dbl_mb(hml(1)),dbl_mb(eig(1)),.true.)

c      value = MA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
c      if (.not. value) 
c     >   call errquit(
c     >    'out of stack memory in psi_diagonalize_hml_assending',0,
c     >     MA_ERR)
c
c
c*     ***** diagonalize the hamiltonian matrix *****
c      call dcopy((ne(1)+ne(2)),0.0d0,0,dbl_mb(eig(1)),1)
c      do ms=1,ispin
c         shift1 = (ms-1)*ne(1)
c         shift2 = (ms-1)*ne(1)*ne(1)
c         if (ne(ms).le.0) go to 30
c
c         call DSYEV('V','U',ne(ms),
c     >              dbl_mb(hml(1)+shift2),ne(ms), 
c     >              dbl_mb(eig(1)+shift1),
c     >              dbl_mb(tmp1(1)),2*ne(1)*ne(1),
c     >              ierr)
c
c   30    continue
c      end do

      
c      value = MA_pop_stack(tmp1(2))
c      if (.not. value) 
c     > call errquit(
c     >   'error popping stack in psi_diagonalize_hml_assending',0,
c     >     MA_ERR)

      return
      end



*     ***************************
*     *				*
*     *		psi_error	*
*     *				*
*     ***************************
      real*8 function psi_error()
      implicit none
#include "errquit.fh"

#include "mafdecls.fh"
#include "psi.fh"

*     ***** local variables ****
      logical value
      integer k,n
      real*8  error,sum,size
      integer tmp1(2)

      value = MA_push_get(mt_dcpl,(npack1),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('out of stack memory in psi_error',0, MA_ERR)


      error = 0.0d0
      size =  dble(ne(1)+ne(2))
      do n=1, (neq(1)+neq(2))
         do k=1,npack1
            dcpl_mb(tmp1(1)+k-1) = dcpl_mb(psi2(1)+k-1+(n-1)*npack1) 
     >                           - dcpl_mb(psi1(1)+k-1+(n-1)*npack1)
         end do
c         call D3dB_cc_dot(1,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)),sum)
         call Pack_cc_dot(1,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)),sum)

         error = error + sum
      end do
      call D1dB_SumAll(error)
      error = dsqrt(error)/size
      
      value = MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory in psi_error',0, MA_ERR)


      psi_error = error
      return
      end

*     ***************************
*     *				*
*     *		rho_error	*
*     *				*
*     ***************************
      real*8 function rho_error()
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "psi.fh"

*     ***** local variables ****
      logical value
      integer k,nx,ny,nz
      real*8  error,scale
      integer tmp1(2)

*     ***** external functions *****
      real*8   lattice_omega
      external lattice_omega

      value = MA_push_get(mt_dbl,(2*nfft3d),'tmp1',tmp1(2),tmp1(1))
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

      value = MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory in rho_error',0, MA_ERR)


      rho_error = error
      return
      end


*     ***************************
*     *                         *
*     *         psi_a_sum       *
*     *                         *
*     ***************************
      real*8 function psi_a_sum(npack1,psi)
      implicit none
      integer npack1
      complex*16 psi(*)

      integer k
      real*8 a,tmp

      a = 0.0d0
      do k=1,npack1
         tmp = dble(psi(k))
         a = a + tmp*tmp
      end do
      call D3dB_SumAll(a)

      psi_a_sum = a
      return
      end




*     ***************************
*     *                         *
*     *         psi_b_sum       *
*     *                         *
*     ***************************
      real*8 function psi_b_sum(npack1,psi)
      implicit none
      integer npack1
      complex*16 psi(*)


      integer k
      real*8 b,tmp

      b = 0.0d0
      do k=1,npack1
         tmp = dimag(psi(k))
         b = b + tmp*tmp
      end do
      call D3dB_SumAll(b)

      psi_b_sum = b
      return
      end

*     **************************************
*     *                                    *
*     *          psi_symm_project          *
*     *                                    *
*     **************************************
      subroutine psi_a_project(npack1,psi)
      implicit none
      integer npack1
      complex*16 psi(*)
      integer k
      real*8 tmp
      do k=1,npack1
        tmp    = dble(psi(k))
        psi(k) = dcmplx(tmp,0.0d0)
      end do
      return
      end
      subroutine psi_b_project(npack1,psi)
      implicit none
      integer npack1
      complex*16 psi(*)
      integer k
      real*8 tmp
      do k=1,npack1
        tmp    = dimag(psi(k))
        psi(k) = dcmplx(0.0d0,tmp)
      end do
      return
      end

      subroutine psi_symm_project(ispin,ne,npack1,psi1)
      implicit none
      integer ispin,ne(2),npack1
      complex*16 psi1(npack1,*)

      integer i
      real*8   a,b
      real*8   psi_a_sum,psi_b_sum
      external psi_a_sum,psi_b_sum

      do i=1,(ne(1)+ne(2))
          a = psi_a_sum(npack1,psi1(1,i))
          b = psi_b_sum(npack1,psi1(1,i))
          if (a.ge.b) then
             call psi_a_project(npack1,psi1(1,i))
          else
             call psi_b_project(npack1,psi1(1,i))
          end if
      end do
      return
      end

      subroutine psi_ab_gen_irrep_names(virtual)
      implicit none
      logical virtual

#include "mafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

      integer irreps(2)
      common / ab_irrep / irreps

      integer k,n,ptr,nn
      real*8  a,b,tmpa,tmpb

      if (virtual) then
          ptr = psi1_excited(1)
          nn  = ne_excited(1)+ne_excited(2)
      else
         ptr = psi1(1)
         nn  = ne(1)+ne(2)
      end if

      if (.not.MA_alloc_get(mt_int,nn,
     >                     'irreps',irreps(2),irreps(1)))
     > call errquit('psi_ab_gen_irrep_names',0, MA_ERR)

      do n=1,nn
         a = 0.0d0
         b = 0.0d0
         do k=1,npack1
            tmpa = dble( dcpl_mb(ptr+k-1+(n-1)*npack1))
            tmpb = dimag(dcpl_mb(ptr+k-1+(n-1)*npack1))
            a = a + tmpa*tmpa
            b = b + tmpb*tmpb
         end do
         call D3dB_SumAll(a)
         call D3dB_SumAll(b)

         if      ((b .lt. 1.0d-6).and.(a .gt. 1.0d-6)) then
            int_mb(irreps(1)+n-1) = 1
         else if ((a .lt. 1.0d-6).and.(b .gt. 1.0d-6)) then
            int_mb(irreps(1)+n-1) = -1
         else
            int_mb(irreps(1)+n-1) = 0
         end if
      end do

      
      return
      end

      subroutine psi_ab_kill_irrep_names()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"

      integer irreps(2)
      common / ab_irrep / irreps
     
      if (.not.MA_free_heap(irreps(2)))
     >  call errquit('psi_ab_gen_irrep_names: error freeing heap',
     >               0, MA_ERR)

      return
      end



*     **************************************
*     *                                    *
*     *         psi_ab_irrep_name          *
*     *                                    *
*     **************************************

*     This function resturns
*        '[ag]' - if psi(n) is purely real
*        '[au]' - if psi(n) is purely imaginary
*        '    ' - if psi(n) is mixed
*
*   Not psi_ab_gen_irrep_names needs to be called before this is used.
*
      character*4 function psi_ab_irrep_name(n)
      implicit none
      integer n

#include "mafdecls.fh"

      integer irreps(2)
      common / ab_irrep / irreps

      character*4 abvalue

      if      (int_mb(irreps(1)+n-1).eq.1) then
         abvalue = '[ag]'
      else if (int_mb(irreps(1)+n-1).eq.-1) then
         abvalue = '[au]'
      else
         abvalue = '    '
      end if

      psi_ab_irrep_name = abvalue
      return
      end


       
*     ***************************
*     *				*
*     *		rho_dipole	*
*     *				*
*     ***************************
*
*     Uses - Calculate_dipole (pspw/lib/psi/dipole.f)
*
      subroutine rho_dipole(dipole)
      implicit none
      real*8 dipole(3)

#include "mafdecls.fh"
#include "psi.fh"

      call Calculate_Dipole(ispin,ne,2*nfft3d,dbl_mb(rho1(1)),dipole)
      return
      end


*     ***************************
*     *				*
*     *		psi_ispin	*
*     *				*
*     ***************************
      integer function psi_ispin()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

      psi_ispin = ispin
      return
      end


*     ***************************
*     *				*
*     *		psi_ne		*
*     *				*
*     ***************************
      integer function psi_ne(ms)
      implicit none
      integer ms

#include "mafdecls.fh"
#include "psi.fh"

      psi_ne = ne(ms)
      return
      end

*     ***************************
*     *				*
*     *		psi_neq		*
*     *				*
*     ***************************
      integer function psi_neq(ms)
      implicit none
      integer ms

#include "mafdecls.fh"
#include "psi.fh"

      psi_neq = neq(ms)
      return
      end



*     ***************************
*     *                         *
*     *    psi_cpmd_start       *
*     *                         *
*     ***************************
      subroutine psi_cpmd_start()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

      logical value

      value = MA_alloc_get(mt_dbl,2*ispin*nfft3d,'rho0',rho0(2),rho0(1))
      if (.not.value)
     >   call errquit('psi_cpmd_start',0,MA_ERR)

      call dcopy(2*ispin*nfft3d,dbl_mb(rho1(1)),1,dbl_mb(rho0(1)),1)
      return
      end

*     ***************************
*     *                         *
*     *    psi_cpmd_end         *
*     *                         *
*     ***************************
      subroutine psi_cpmd_end()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

      if (.not.MA_free_heap(rho0(2)))
     >   call errquit('psi_cpmd_end',0,MA_ERR)
      return
      end


*     ***************************
*     *                         *
*     *    psi_cpmd_step        *
*     *                         *
*     ***************************
      subroutine psi_cpmd_step(dte)
      implicit none
      real*8 dte

#include "mafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

      logical  control_precondition
      external control_precondition
      integer  control_ks_algorithm
      external control_ks_algorithm
      real*8   control_tole
      external control_tole
  

*     **** psi2 = 2*psi1 - psi0 + dt*dt/fmass*Hpsi ****
c      call electron_cpmd_update(dcpl_mb(psi0(1)),
c     >                          dcpl_mb(psi1(1)),
c     >                          dcpl_mb(psi2(1)),
c     >                          dbl_mb(hml(1)),
c     >                          dte)
c      call Dneall_f_ortho(0,dcpl_mb(psi2(1)),npack1)
c      write(*,*) "psi1 ortho:"
c      call OrthoCheck_geo(ispin,ne,dcpl_mb(psi1(1)))
c*     **** lagrange multiplier corrections ****
c      if (pspw_SIC().or.occupation_on) then
c        call psi_lmbda_sic(ispin,ne,(neq(1)+neq(2)),npack1,
c     >                 dcpl_mb(psi1(1)),dcpl_mb(psi2(1)),dte,
c     >                 dbl_mb(lmd_cpmd(1)),
c     >                 dbl_mb(tmp_L_cpmd(1)),ierr)
c      else
c        call psi_lmbda(ispin,ne,(neq(1)+neq(2)),npack1,
c     >                 dcpl_mb(psi1(1)),dcpl_mb(psi2(1)),dte,
c     >                 dbl_mb(lmd_cpmd(1)),
c     >                 dbl_mb(tmp_L_cpmd(1)),ierr)
c      end if
c      write(*,*) "psi2 ortho:"
c      call OrthoCheck_geo(ispin,ne,dcpl_mb(psi2(1)))

      call dcopy(2*ispin*nfft3d,dbl_mb(rho1(1)),1,dbl_mb(rho2(1)),1)
      call dscal(2*ispin*nfft3d, 2.0d0,dbl_mb(rho2(1)),1)
      call daxpy(2*ispin*nfft3d,-1.0d0,
     >           dbl_mb(rho0(1)),1,dbl_mb(rho2(1)),1)
      call dcopy(2*ispin*nfft3d,dbl_mb(rho1(1)),1,dbl_mb(rho0(1)),1)

      call psi_set_density(1,dbl_mb(rho2(1)))

*     **** diaganolize KS matrix ****
      call psi_KS_update(1,
     >                   control_ks_algorithm(),
     >                   control_precondition(),
     >                   control_tole())

      return
      end



*     ***************************
*     *				*
*     *	    psi_initialize 	*
*     *				*
*     ***************************

      logical function psi_initialize()
      implicit none 


#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical value,psi_nogrid
      integer nemax
      real*8 sum1,sum2,sum3
      integer hversion,hnfft(3),hispin,hne(2)
      real*8 hunita(3,3)
      integer rtdb,ind,vers
      integer  control_rtdb,control_ngrid,control_symmetry
      external control_rtdb,control_ngrid,control_symmetry
      integer  psi_get_version
      external psi_get_version
      character*30 filename
      character*30 control_input_psi
      external     control_input_psi
      logical  wvfnc_expander,Dneall_m_allocate,band_reformat_c_wvfnc
      external wvfnc_expander,Dneall_m_allocate,band_reformat_c_wvfnc
      logical  psp_pawexist,control_print
      external psp_pawexist,control_print
      integer          control_fractional_smeartype
      double precision control_fractional_kT
      external         control_fractional_smeartype
      external         control_fractional_kT


      ne_excited(1) = 0
      ne_excited(2) = 0

*     **** reformat wavefunction if it is a band wavefunction ****
      vers = psi_get_version()
      if (vers.eq.5) then
        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER) then
          value= band_reformat_c_wvfnc(1)
        end if
      end if


*     *****  get ispin,ne,neq,nfft3d,npack0,npack1 ****
      call psi_get_ne_occupation(ispin,ne,smearoccupation)
      call Dneall_neq(neq)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)
      call Pack_npack(0,npack0)
      nemax = ne(1)+ne(2)
      occupation_on = .false.
      if (smearoccupation.gt.0) occupation_on = .true.


*     **** allocate memory ****
      value = MA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'psi1',psi1(2),psi1(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,4*nfft3d,
     >                     'rho1',rho1(2),rho1(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,4*nfft3d,
     >                     'rho2',rho2(2),rho2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack0,
     >                     'dng1',dng1(2),dng1(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack0,
     >                     'dng2',dng2(2),dng2(1))
c      value = value.and.
c     >        MA_alloc_get(mt_dbl,(2*nemax*nemax),'hml',hml(2),hml(1))
      value = value.and.Dneall_m_allocate(0,hml)

      value = value.and.
     >        MA_alloc_get(mt_dbl,(2*nemax),'eig',eig(2),eig(1))

      if (occupation_on) then
        value = value.and.
     >        MA_alloc_get(mt_dbl,(2*nemax),'occ1',occ1(2),occ1(1))
        value = value.and.
     >        MA_alloc_get(mt_dbl,(2*nemax),'occ2',occ2(2),occ2(1))
        smeartype = control_fractional_smeartype()
        smearkT   = control_fractional_kT()
      end if

      value = value.and.
     >        MA_alloc_get(mt_dbl,4*nfft3d,
     >                     'rho1_all',rho1_all(2),rho1_all(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,4*nfft3d,
     >                     'rho2_all',rho2_all(2),rho2_all(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

*     *****  read initial wavefunctions into psi1  ****
      rtdb = control_rtdb()
      if (.not.rtdb_get(rtdb,'nwpw:psi_nogrid',
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
        call Parallel_taskid(taskid)

        call ga_sync()
        value = rtdb_parallel(.false.)
        call ga_sync()
        if (taskid.eq.MASTER) then
       
          filename =  control_input_psi()

          ind = index(filename,' ') - 1
          if (.not. rtdb_cput(rtdb,'xpndr:old_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: rtdb_cput failed', 0, RTDB_ERR)

          if (.not. rtdb_cput(rtdb,'xpndr:new_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: rtdb_cput failed', 0, RTDB_ERR)

          if (.not. rtdb_put(rtdb,'xpndr:ngrid',mt_int,3,hnfft))
     >     call errquit(
     >     'wvfnc_expander_input: rtdb_put failed', 0, RTDB_ERR)

          if (control_print(print_medium)) then
            write(luout,*)
            write(luout,*) "Grid is being converted:"
            write(luout,*) "------------------------"
            write(luout,*)
            write(luout,*) "To turn off automatic grid conversion:"
            write(luout,*)
            write(luout,*) "set nwpw:psi_nogrid .false."
            write(luout,*)
          endif
          value = wvfnc_expander(rtdb)

        end if
        call ga_sync()
        value = rtdb_parallel(.true.)
        value = .true.

      end if

      end if

      call psi_read(ispin,ne,dcpl_mb(psi1(1)),
     >              smearoccupation,dbl_mb(occ1(1)))


      call psi_history_read(ispin,ne,
     >                      dcpl_mb(psi1(1)),
     >                      dcpl_mb(psi2(1)))
     

*     **** force inversion symmetry ****
      if (control_symmetry().eq.1)  then
         call Parallel_taskid(taskid)
         if ((taskid.eq.MASTER).and.
     >       (control_print(print_medium))) then
         write(luout,*)
         write(luout,*) 
     >   "Projecting wavefunctions to have inversion symmetry"
         write(luout,*)
         end if
         call psi_symm_project(ispin,neq,npack1,dcpl_mb(psi1(1)))
      end if 

*     **** Ortho Check ****
      if (psp_pawexist()) then
         call psp_overlap_S(ispin,neq,dcpl_mb(psi1(1)),dcpl_mb(psi2(1)))
         call Grsm_gg_trace(npack1,(neq(1)+neq(2)),
     >                        dcpl_mb(psi1(1)),
     >                        dcpl_mb(psi2(1)),
     >                        sum2)
      else
         call Grsm_gg_trace(npack1,(neq(1)+neq(2)),
     >                        dcpl_mb(psi1(1)),
     >                        dcpl_mb(psi1(1)),
     >                        sum2)
      end if
      call D1dB_SumAll(sum2)


      sum1 = dble(ne(1) + ne(2))
      if (dabs(sum2-sum1).gt.1.0d-10) then
         call Parallel_taskid(taskid)
         if (psp_pawexist()) then
         call Dneall_f_Sortho(0,dcpl_mb(psi1(1)),
     >                          dcpl_mb(psi2(1)),npack1)
         call psp_overlap_S(ispin,neq,dcpl_mb(psi1(1)),
     >                                dcpl_mb(psi2(1)))
         call Grsm_gg_trace(npack1,(neq(1)+neq(2)),
     >                        dcpl_mb(psi1(1)),
     >                        dcpl_mb(psi2(1)),
     >                        sum3)
         else
c         call Dneall_f_ortho(0,dcpl_mb(psi1(1)),npack1)
         call Dneall_f_GramSchmidt(0,dcpl_mb(psi1(1)),npack1)
         call Grsm_gg_trace(npack1,(neq(1)+neq(2)),
     >                        dcpl_mb(psi1(1)),
     >                        dcpl_mb(psi1(1)),
     >                        sum3)
         end if
         call D1dB_SumAll(sum3)
         if ((taskid.eq.MASTER).and.(control_print(print_medium)))
     >    write(luout,*) 
     >     "Warning - Gram-Schmidt being performed on psi:",
     >               sum1,sum2,sum3,dabs(sum2-sum1)


      end if

      psi_initialize = value
      return
      end



*     ***************************
*     *				*
*     *	  psi_tmp_write  	*
*     *				*
*     ***************************
      subroutine psi_tmp_write()
      implicit none

#include "mafdecls.fh"
#include "psi.fh"

*     ***** write psi1 wavefunctions ****
      call psi_write(ispin,ne,dcpl_mb(psi1(1)),
     >               smearoccupation,dbl_mb(occ1(1)))

      return
      end



*     ***************************
*     *				*
*     *		psi_finalize	*
*     *				*
*     ***************************

      logical function psi_finalize(wpsi)
      implicit none 
      logical wpsi

#include "errquit.fh"
#include "mafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all


*     **** local variables ****
      logical value

      logical  Dneall_m_free
      external Dneall_m_free

*     ***** write psi1 wavefunctions ****
      if (wpsi) then
        call psi_write(ispin,ne,dcpl_mb(psi1(1)),
     >                 smearoccupation,dbl_mb(occ1(1)))
        call psi_history_write(ispin,ne,dcpl_mb(psi1(1)))
      end if
    
      value = MA_free_heap(eig(2))
      value = value.and.Dneall_m_free(hml)
      value = value.and.MA_free_heap(dng2(2))
      value = value.and.MA_free_heap(dng1(2))
      value = value.and.MA_free_heap(rho2(2))
      value = value.and.MA_free_heap(rho1(2))
      value = value.and.MA_free_heap(psi2(2))
      value = value.and.MA_free_heap(psi1(2))
      value = value.and.MA_free_heap(rho2_all(2))
      value = value.and.MA_free_heap(rho1_all(2))
      if (occupation_on) then
         value = value.and.MA_free_heap(occ2(2))
         value = value.and.MA_free_heap(occ1(2))
      end if

      if (.not. value) 
     >  call errquit('psi_finalize: error freeing heap',0, MA_ERR)

      psi_finalize = value
      return
      end


*     ***************************
*     *                         *
*     *      psi_ne_excited     *
*     *                         *
*     ***************************
      integer function psi_ne_excited(ms)
      implicit none
      integer ms

#include "mafdecls.fh"
#include "psi.fh"

      psi_ne_excited = ne_excited(ms)
      return
      end

*     ***************************
*     *				*
*     *     epsi_initialize 	*
*     *				*
*     ***************************

      logical function epsi_initialize()
      implicit none 
#include "errquit.fh"

#include "mafdecls.fh"
#include "rtdb.fh"
#include "psi.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical value,psi_nogrid
      integer nemax,ispin0

      integer hversion,hnfft(3),hispin,hne(2)
      real*8 hunita(3,3)
      integer rtdb,ind
      integer  control_rtdb,control_ngrid
      external control_rtdb,control_ngrid
      character*30 filename
      character*30 control_input_epsi
      external     control_input_epsi
      logical  wvfnc_expander
      external wvfnc_expander
      integer  control_symmetry
      external control_symmetry


*     ***** get ispin, and ne, and nfft3d ****
      call psi_get_ne_excited(ispin0,ne_excited)
      nemax  = ne_excited(1)  + ne_excited(2)

*     **** allocate memory ****
      value = MA_alloc_get(mt_dcpl,npack1*(nemax),
     >         'psi2_excited',psi2_excited(2),psi2_excited(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*(nemax),
     >         'psi1_excited',psi1_excited(2),psi1_excited(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(2*nemax),'eig_excited',
     >                     eig_excited(2),eig_excited(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      call dcopy(2*npack1*nemax,0.0d0,0,dcpl_mb(psi2_excited(1)),1)
      call dcopy(2*npack1*nemax,0.0d0,0,dcpl_mb(psi1_excited(1)),1)
      call dcopy(2*nemax,0.0d0,0,dbl_mb(eig_excited(1)),1)


*     *****  read initial wavefunctions into psi1  ****
      rtdb = control_rtdb()
      if (.not.rtdb_get(rtdb,'nwpw:psi_nogrid',
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
        value = rtdb_parallel(.false.)
        if (taskid.eq.MASTER) then


          ind = index(filename,' ') - 1
          if (.not. rtdb_cput(rtdb,'xpndr:old_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: rtdb_cput failed', 0, RTDB_ERR)

          if (.not. rtdb_cput(rtdb,'xpndr:new_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: rtdb_cput failed', 0, RTDB_ERR)

          if (.not. rtdb_put(rtdb,'xpndr:ngrid',mt_int,3,hnfft))
     >     call errquit(
     >     'wvfnc_expander_input: rtdb_put failed', 0, RTDB_ERR)

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
        value = rtdb_parallel(.true.)

      end if

      end if

*     *****  read initial wavefunctions into psi1  ****
      call epsi_read(ispin0,ne_excited,dcpl_mb(psi1_excited(1)))


*     **** force inversion symmetry ****
      if (control_symmetry().eq.1)  then
         call Parallel_taskid(taskid)
         if (taskid.eq.MASTER) then
         write(*,*)
         write(*,*)
     >   "Projecting virtual wavefunctions to have inversion symmetry"
         write(*,*)
         end if
         call psi_symm_project(ispin0,ne_excited,npack1,
     >                         dcpl_mb(psi1_excited(1)))
      end if



*     **** Ortho Check ****
c     call Grsm_gg_trace(npack1,(ne_excited(1)+ne_excited(2)),
c    >                        dcpl_mb(psi1_excited(1)),
c    >                        dcpl_mb(psi1_excited(1)),
c    >                        sum2)
c     
c     sum1 = dble(ne_excited(1) + ne_excited(2))
c     if (dabs(sum2-sum1).gt.1.0d-10) then
c
c        call Parallel_taskid(taskid)
c        call Grsm_g_MakeOrtho(npack1,ne_excited(1),
c    >                         dcpl_mb(psi1_excited(1)))
c        if (ispin.gt.1) then
c          call Grsm_g_MakeOrtho(npack1,ne_excited(2),
c    >                           dcpl_mb(psi1_excited(1)
c    >                                  +ne_excited(1)*npack1))
c        end if
c        call Grsm_gg_trace(npack1,(ne_excited(1)+ne_excited(2)),
c    >                        dcpl_mb(psi1_excited(1)),
c    >                        dcpl_mb(psi1_excited(1)),
c    >                        sum2)
c        if (taskid.eq.MASTER)
c    >    write(*,*) "Warning - Gram-Schmidt being performed on epsi:",
c    >               dabs(sum2-sum1)
c
c     end if



      epsi_initialize = value
      return
      end



*     ***************************
*     *				*
*     *		epsi_finalize	*
*     *				*
*     ***************************

      logical function epsi_finalize(writepsi)
      implicit none 
      logical writepsi

#include "mafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     **** local variables ****
      logical value

*     ***** write psi1 wavefunctions ****
      if (writepsi)
     >  call epsi_write(ispin,ne_excited,dcpl_mb(psi1_excited(1)))
    
      value = MA_free_heap(eig_excited(2))
      value = value.and.MA_free_heap(psi2_excited(2))
      value = value.and.MA_free_heap(psi1_excited(2))
      if (.not. value) 
     >  call errquit('epsi_finalize: error freeing heap',0, MA_ERR)

      epsi_finalize = value
      return
      end


*     ***********************
*     *			    *
*     *	    psi_Mulliken    *
*     *			    *
*     ***********************

      subroutine psi_Mulliken(rtdb)
      implicit none 
      integer rtdb

#include "mafdecls.fh"
#include "psi.fh"


*     **** Lubin Water Analysis ***
      call pspw_Lubin_water_analysis(rtdb,ispin,ne,2*nfft3d,
     >                                 dbl_mb(rho1(1)))

      call pspw_analysis(0,rtdb,ispin,ne,dcpl_mb(psi2(1)),
     >                                   dbl_mb(eig(1)))

      call pspw_gen_APC(ispin,ne,dcpl_mb(dng1(1)))
      call pspw_print_APC(6)

      call pspw_gen_Efield(rtdb,ispin,dbl_mb(rho1(1)),dcpl_mb(dng1(1)))

      return
      end

*     ***********************
*     *                     *
*     *    epsi_Mulliken    *
*     *                     *
*     ***********************

      subroutine epsi_Mulliken(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "psi.fh"

      call pspw_analysis(1,rtdb,ispin,ne_excited,
     >                   dcpl_mb(psi1_excited(1)),
     >                   dbl_mb(eig_excited(1)))
      return
      end




*     ***********************
*     *                     *
*     *     psi_DOS         *
*     *                     *
*     ***********************

      subroutine psi_DOS(rtdb)
      implicit none
      integer rtdb

#include "rtdb.fh"
#include "mafdecls.fh"
#include "psi.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer npoints,ii
      integer weight(2),nemax
      real*8 emin,emax,alpha
      character*255 filename

      nemax = ne(1)
      value = MA_push_get(mt_dbl,(nemax),'weight',weight(2),weight(1))
      if (.not. value) 
     >  call errquit('psi_dos:out of stack memory',0, MA_ERR)
      call dcopy(nemax,1.0d0,0,dbl_mb(weight(1)),1)


      if (.not.rtdb_get(rtdb,'dos:alpha',mt_dbl,1,alpha)) then
        alpha = 0.05d0/27.2116d0
      end if

      if (.not.rtdb_get(rtdb,'dos:npoints',mt_int,1,npoints)) then
        npoints = 500
      end if

      if (.not.rtdb_get(rtdb,'dos:emin',mt_dbl,1,emin)) then
         emin = 99999.0d0
         do ii=1,ne(1)+ne(2)
           if (dbl_mb(eig(1)+ii-1).lt.emin) emin = dbl_mb(eig(1)+ii-1)
         end do
         emin = emin - 0.1d0
      end if

      if (.not.rtdb_get(rtdb,'dos:emax',mt_dbl,1,emax)) then
         emax = -99999.0d0
         do ii=1,ne(1)+ne(2)
           if (dbl_mb(eig(1)+ii-1).gt.emax) emax = dbl_mb(eig(1)+ii-1)
         end do
         emax = emax + 0.1d0
      end if

*     **** generate DENSITY OF STATES *****
      if (ispin.eq.1) then
        filename = "DOS_both"
        call densityofstates(filename,
     >                     dbl_mb(eig(1)),dbl_mb(weight(1)),ne(1),
     >                     1.0d0,alpha,npoints,emin,emax)
      end if

      if (ispin.eq.2) then
        filename = "DOS_alpha"
        call densityofstates(filename,
     >                     dbl_mb(eig(1)),dbl_mb(weight(1)),ne(1),
     >                     1.0d0,alpha,npoints,emin,emax)
        filename = "DOS_beta"
        call densityofstates(filename,
     >               dbl_mb(eig(1)+ne(1)),dbl_mb(weight(1)),ne(2),
     >               -1.0d0,alpha,npoints,emin,emax)
      end if

 
      value = MA_pop_stack(weight(2))
      if (.not. value)
     >  call errquit('psi_dos: error freeing stack',0, MA_ERR)

      return
      end

*     ***********************
*     *                     *
*     *     epsi_DOS        *
*     *                     *
*     ***********************

      subroutine epsi_DOS(rtdb)
      implicit none
      integer rtdb

#include "rtdb.fh"
#include "mafdecls.fh"
#include "psi.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer weight(2),npoints,ii
      real*8 emin,emax,alpha
      character*255 filename

      value = MA_push_get(mt_dbl,(ne_excited(1)+ne_excited(2)),
     >                    'weight',weight(2),weight(1))
      if (.not. value)
     >  call errquit('epsi_dos:out of stack memory',0, MA_ERR)
      call dcopy((ne_excited(1)+ne_excited(2)),
     >          1.0d0,0,dbl_mb(weight(1)),1)

      if (.not.rtdb_get(rtdb,'dos:alpha',mt_dbl,1,alpha)) then
        alpha = 0.05d0/27.2116d0
      end if

      if (.not.rtdb_get(rtdb,'dos:npoints',mt_int,1,npoints)) then
        npoints = 500
      end if

      if (.not.rtdb_get(rtdb,'dos:emin',mt_dbl,1,emin)) then
         emin = 99999.0d0
         do ii=1,ne_excited(1)+ne_excited(2)
           if (dbl_mb(eig_excited(1)+ii-1).lt.emin) 
     >       emin = dbl_mb(eig_excited(1)+ii-1)
         end do
         emin = emin - 0.1d0
      end if

      if (.not.rtdb_get(rtdb,'dos:emax',mt_dbl,1,emax)) then
         emax = -99999.0d0
         do ii=1,ne_excited(1)+ne_excited(2)
           if (dbl_mb(eig_excited(1)+ii-1).gt.emax) 
     >       emax = dbl_mb(eig_excited(1)+ii-1)
         end do
         emax = emax + 0.1d0
      end if

*     **** generate DENSITY OF STATES *****
      if (ispin.eq.1) then
        filename = "VDOS_both"
        call densityofstates(filename,
     >                     dbl_mb(eig_excited(1)),dbl_mb(weight(1)),
     >                     ne_excited(1),
     >                     1.0d0,alpha,npoints,emin,emax)
      end if

      if (ispin.eq.2) then
        filename = "VDOS_alpha"
        call densityofstates(filename,
     >                     dbl_mb(eig_excited(1)),dbl_mb(weight(1)),
     >                     ne_excited(1),
     >                     1.0d0,alpha,npoints,emin,emax)
        filename = "VDOS_beta"
        call densityofstates(filename,
     >           dbl_mb(eig_excited(1)+ne_excited(1)),dbl_mb(weight(1)),
     >           ne_excited(2),
     >           -1.0d0,alpha,npoints,emin,emax)
      end if

      value = MA_pop_stack(weight(2))
      if (.not. value)
     >  call errquit('epsi_dos: error freeing stack',0, MA_ERR)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine psi_polariz()
      implicit none
#include "psi.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
      integer psirx(2),asize
      logical val
      asize=(ne(1)+ne(2))*nfft3d*2
      val=MA_PUSH_GET(mt_dbl,asize,"psir",psirx(2),psirx(1))
      if (.not.val) then
        call errquit("psi_polariz stack empty",0,0)
      end if
      call berry_phase_pol(psi2,psirx)
      val=MA_POP_STACK(psirx(2))
      if (.not.val) then
        call errquit("psi_polariz failed pop stack",0,0)
      end if
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine epsi_generate_kb_vnm(vnm)
      implicit none
      real*8 vnm(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     **** local variables ****
      logical value
      integer ms,i,j,vshift,ishift,neall(2)
      integer Gx(2),Gy(2),Gz(2),tmp(2),psii_ptr,psij_ptr
      real*8  vv(3)

*     **** external functions ****
      integer  Dneall_mne_size,G_indx
      external Dneall_mne_size,G_indx
       
      neall(1) = ne(1)+ne_excited(1)
      neall(2) = ne(2)+ne_excited(2)

      value = MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
      value = value.and.
     >        MA_push_get(mt_dbl, 2*npack1,'tmp',tmp(2),tmp(1))
         if (.not. value)
     >      call errquit('epsi_generate_kb_vnm:pushing stack',1,MA_ERR)


*     **** define Gx,Gy and Gz in packed space ****
      call D3dB_t_Copy(1,dbl_mb(G_indx(1)),dbl_mb(Gx(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(2)),dbl_mb(Gy(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(3)),dbl_mb(Gz(1)))
      call Pack_t_pack(1,dbl_mb(Gx(1)))
      call Pack_t_pack(1,dbl_mb(Gy(1)))
      call Pack_t_pack(1,dbl_mb(Gz(1)))


      vshift = Dneall_mne_size(0,neall)
      call dcopy(3*vshift,0.0d0,0,vnm,1)

      do ms=1,ispin
         ishift = (ms-1)*ne(1)

         do j=1,neall(ms)
            if (j.le.ne(ms)) then
               psij_ptr =(j-1+ishift)*npack1 + psi1(1)
            else
               psij_ptr =(j-ne(ms)-1+ishift)*npack1 + psi1_excited(1)
            end if

            do i=j,neall(ms)
               if (i.le.ne(ms)) then
                  psii_ptr =(i-1+ishift)*npack1 + psi1(1)
               else
                  psii_ptr =(i-ne(ms)-1+ishift)*npack1 + psi1_excited(1)
               end if

               call Pack_tc_Mul(1,
     >                          dbl_mb(Gx(1)),
     >                          dcpl_mb(psij_ptr),
     >                          dbl_mb(tmp(1)))
               call Pack_cc_dot(1,
     >                          dcpl_mb(psii_ptr),
     >                          dbl_mb(tmp(1)),
     >                          vv(1))
               call Pack_tc_Mul(1,
     >                          dbl_mb(Gy(1)),
     >                          dcpl_mb(psij_ptr),
     >                          dbl_mb(tmp(1)))
               call Pack_cc_dot(1,
     >                          dcpl_mb(psii_ptr),
     >                          dbl_mb(tmp(1)),
     >                          vv(2))
               call Pack_tc_Mul(1,
     >                          dbl_mb(Gz(1)),
     >                          dcpl_mb(psij_ptr),
     >                          dbl_mb(tmp(1)))
               call Pack_cc_dot(1,
     >                          dcpl_mb(psii_ptr),
     >                          vv(3))
   

               call Dneall_mne_set_value(vv(1),0,neall,ms,i,j,vnm)
               call Dneall_mne_set_value(vv(2),0,neall,ms,i,j,
     >                                   vnm(1+vshift))
               call Dneall_mne_set_value(vv(3),0,neall,ms,i,j,
     >                                   vnm(1+vshift+vshift))
               if (i.ne.j) then
                  call Dneall_mne_set_value(vv(1),0,neall,ms,j,i,vnm)
                  call Dneall_mne_set_value(vv(2),0,neall,ms,j,i,
     >                                      vnm(1+vshift))
                  call Dneall_mne_set_value(vv(3),0,neall,ms,j,i,
     >                                      vnm(1+vshift+vshift))
               end if
            end do
         end do
      end do

      value =           MA_pop_stack(tmp(2))
      value = value.and.MA_pop_stack(Gz(2))
      value = value.and.MA_pop_stack(Gy(2))
      value = value.and.MA_pop_stack(Gx(2))
      if (.not. value)
     >   call errquit('epsi_generate_kb_vnm:popping stack',1,MA_ERR)

      
      return
      end


