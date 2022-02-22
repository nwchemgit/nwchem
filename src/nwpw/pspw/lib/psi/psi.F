
* $Id$
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

#include "bafdecls.fh"
#include "psi.fh"

      !*** local variables ***
      integer maxit_orb
      integer ii,l,l2
      real*8  sum,maxerror,error_out,e0

      !*** external functions ***
      real*8   control_tole
      external control_tole

      !call psi_gen_density_potentials(1)
      maxit_orb=120
      maxerror = control_tole()

      do ii=1,(ne(1)+ne(2))
         l2 = 0

         !*** orthogonalize to lower orbitals  ****
 2       l2 = l2 + 1
         call psi_project_out_f_orb1(
     >           ii,
     >           dcpl_mb(psi1(1)+(ii-1)*npack1))

         !*** normalize ****
         call Pack_cc_dot(1,
     >            dcpl_mb(psi1(1) +(ii-1)*npack1),
     >            dcpl_mb(psi1(1) +(ii-1)*npack1),
     >            sum)
         sum = 1.0d0/dsqrt(sum)
         call Pack_c_SMul1(1,sum,dcpl_mb(psi1(1) +(ii-1)*npack1))

         !*** minimize orbital ****
          l = 0
 3        call psi_KS_update_f_orb(maxit_orb,
     >                               maxerror,
     >                               0.001d0,ii,error_out,e0)
          !write(*,*) "e0:",ii,l,e0,error_out
          l = l+1
          if ((error_out.gt.maxerror).and.(l.le.(1+(l2-1)*3))) go to 3
          if (((error_out.gt.maxerror).or.(e0.gt.4.0d0))
     >        .and.(l2.le.1)) then
           call Pack_c_Zero(1,dcpl_mb(psi1(1) +(ii-1)*npack1))
           call Pack_c_setzero(1,1.0d0,dcpl_mb(psi1(1) +(ii-1)*npack1))
           go to 2
          end if

          dbl_mb(eig(1)+ii-1) = e0

      end do
      call psi_sort_f_orb()
      
     
      return
      end

      subroutine psi_sort_f_orb()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "psi.fh"

      logical value
      integer i,j,ii,jj,ms
      integer r1(2)
      real*8  ei,ej

      value = BA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
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

      value = BA_pop_stack(r1(2))
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
      
#include "bafdecls.fh"
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

      value = BA_push_get(mt_dcpl,npack1,'t0',t0(2),t0(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'g',g(2),g(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'t',t(2),t(1))
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
 4    value =           BA_pop_stack(t(2)) 
      value = value.and.BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(r1(2))
      value = value.and.BA_pop_stack(t0(2))
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
      
#include "bafdecls.fh"
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
      value = BA_push_get(mt_dcpl,npack1,'orb',
     >                       orb(2),orb(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'g',
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
      value =           BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(orb(2))      
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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



************************ gradient virtural orbital Part ************************
*     ***********************************
*     *                                 *
*     *      psi_gen_hml_virtual        *
*     *                                 *
*     ***********************************
      subroutine psi_gen_hml_virtual(assending)
      implicit none
      logical assending

#include "bafdecls.fh"
#include "psi.fh"
#include "errquit.fh"

      !*** local variables ***
      integer ii,jj,jjstart,jjend,ms,mshift,indx,indxt
      integer Horb(2)
      real*8  tsum

*     **** allocate temporary memory from stack ****
      if (.not.BA_push_get(mt_dcpl,npack1,'Horb',Horb(2),Horb(1)))
     >   call errquit('psi_gen_hml_virtual:pusing stack',0,MA_ERR)

      do ii=1,(ne_excited(1)+ne_excited(2))

         !*** orthogonalize to lower orbitals  ****
         call psi_project_out_virtual1(
     >           ii,
     >           dcpl_mb(psi1_excited(1)+(ii-1)*npack1))

         !*** normalize ****
         call Pack_cc_dot(1,
     >            dcpl_mb(psi1_excited(1) +(ii-1)*npack1),
     >            dcpl_mb(psi1_excited(1) +(ii-1)*npack1),
     >            tsum)
         tsum = 1.0d0/dsqrt(tsum)
         call Pack_c_SMul1(1,tsum,
     >            dcpl_mb(psi1_excited(1) +(ii-1)*npack1))

         if (ii.le.ne_excited(1)) then
            ms = 1
            jjstart = 1
            jjend   = ii
            mshift = 0
         else
            ms = 2
            jjstart = ne_excited(1)+1
            jjend   = ii
            mshift = ne_excited(1)*ne_excited(1)
         end if
c         call Pack_cc_dot(1,
c     >               dcpl_mb(psi1_excited(1) +(ii-1)*npack1),
c     >               dcpl_mb(psi1_excited(1) +(ii-1)*npack1),
c     >               tsum)
         !write(*,*) "ii,norm=",ii,tsum

         call electron_get_gradient_virtual(ms,
     >             dcpl_mb(psi1_excited(1)+(ii-1)*npack1),
     >             dcpl_mb(Horb(1)))

         do jj=jjstart,jjend
c         call Pack_cc_dot(1,
c     >               dcpl_mb(psi1_excited(1) +(ii-1)*npack1),
c     >               dcpl_mb(psi1_excited(1) +(jj-1)*npack1),
c     >               tsum)
         !write(*,*) "ii,jj,norm=",ii,jj,tsum
            call Pack_cc_dot(1,
     >               dcpl_mb(psi1_excited(1) +(jj-1)*npack1),
     >               dcpl_mb(Horb(1)),
     >               tsum)
            !write(*,*) "ii,jj,e=",ii,jj,tsum
            indx  = (ii-jjstart)+(jj-jjstart)*ne_excited(ms)+mshift
            dbl_mb(hml_excited(1)+indx)  = tsum
            if (ii.ne.jj) then
               indxt = (jj-jjstart)+(ii-jjstart)*ne_excited(ms)+mshift
               dbl_mb(hml_excited(1)+indxt) = tsum
            end if
         end do
      end do

*     **** deallocate temporary memory from stack ****
      if (.not.BA_pop_stack(Horb(2)))
     >   call errquit('psi_gen_hml_virtual:popping stack',0,MA_ERR)

c      call Dnexall_m_diagonalize(0,ispin,ne_excited,
c     >                           dbl_mb(hml_excited(1)),
c     >                           dbl_mb(eig_excited(1)),.true.)
      call Dnexall_m_diagonalize(0,ispin,ne_excited,
     >                           dbl_mb(hml_excited(1)),
     >                           dbl_mb(eig_excited(1)),assending)
      call  Dnexall_fmf_Multiply(0,ispin,ne_excited,
     >                           dcpl_mb(psi1_excited(1)),npack1,
     >                           dbl_mb(hml_excited(1)),1.0d0,
     >                           dcpl_mb(psi2_excited(1)),0.0)
      indx = psi1_excited(1)
      psi1_excited(1) = psi2_excited(1)
      psi2_excited(1) = indx
c      call dcopy(2*npack1,
c     >           dcpl_mb(psi1_excited(1)+3*npack1),1,
c     >           dcpl_mb(psi1(1)),1)


      return
      end


      
************************ CI virtural orbital Part ************************
*     ***********************************
*     *                                 *
*     *    psi_minimize_virtual_CI      *
*     *                                 *
*     ***********************************
      subroutine psi_minimize_virtual_CI(ci_algorithm)
      implicit none
      integer ci_algorithm
      
#include "bafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"
#include "psi.fh"

      real*8 Hgg
      common / CI_Hgg_common / Hgg
      
      !*** local variables ***
      integer taskid,MASTER
      parameter (MASTER=0)

      logical oprint
      integer maxit_sweeps,maxit_orb,maxit_ls
      integer ii,l
      real*8  sum,maxerror,error_out,e0,eci
      integer g93(2),eig_eci(2)
         
      !*** external functions ***
      logical  control_print
      external control_print
      integer  control_CI_maxit_orb,control_CI_maxit_linesearch
      external control_CI_maxit_orb,control_CI_maxit_linesearch
      integer  control_CI_maxit_sweeps
      external control_CI_maxit_sweeps
      real*8   control_tole
      external control_tole

      call Parallel_taskid(taskid)
      oprint= ((taskid.eq.MASTER).and.control_print(print_medium))

      Hgg = 0.0d0

      if (.not.BA_push_get(mt_dcpl,npack1,'g93',g93(2),g93(1)))
     >   call errquit('psi_CI_minimize_virtual:out of stack',0,MA_ERR)
      if (.not.BA_push_get(mt_dbl,ne_excited(1)+ne_excited(2),
     >                     'eig_eci',eig_eci(2),eig_eci(1)))
     >   call errquit('psi_CI_minimize_virtual:out of stack',1,MA_ERR)

      !call psi_gen_density_potentials(1)
      maxit_sweeps = control_CI_maxit_sweeps()
      maxit_orb = control_CI_maxit_orb()
      maxit_ls  = control_CI_maxit_linesearch()
      maxerror  = control_tole()
      if (oprint) then
         write(luout,*) 
         write(luout,*) "COVOs Mimimization"
         write(luout,*) "------------------"
         write(luout,'(A,I9)') "CI gradient algorithm = ",ci_algorithm
         write(luout,'(A,I9)') "maxit_sweeps          = ",maxit_sweeps
         write(luout,'(A,I9)') "maxit_orb             = ",maxit_orb
         write(luout,'(A,I9)') "maxit_linesearch      = ",maxit_ls
         write(luout,'(A,E9.3)') "maxerror              = ",maxerror
         write(luout,*) 
      end if
      !maxerror = 1.0d-12
      eci = 9.99d9
      e0  = 9.99d9

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
         call Pack_c_SMul1(1,sum,
     >            dcpl_mb(psi1_excited(1) +(ii-1)*npack1))

         call psi_get_gradient_virtual(ii,dcpl_mb(g93(1)))
         call Pack_cc_dot(1,dcpl_mb(psi1_excited(1)+(ii-1)*npack1),
     >                      dcpl_mb(g93(1)),e0)
         e0 = e0

         !*** minimize orbital ****

          l = 0
 2        call psi_CI_update_virtual(ci_algorithm,
     >                               maxit_orb,maxit_ls,
     >                               maxerror,
     >                               0.001d0,ii,error_out,e0,eci)
          l = l+1
          !if ((error_out.gt.maxerror).and.(l.le.4)) go to 2
          if ((error_out.gt.maxerror).and.(l.lt.maxit_sweeps)) go to 2

          dbl_mb(eig_excited(1)+ii-1) = e0
          dbl_mb(eig_eci(1)+ii-1)     = eci

      end do
      call psi_sort_virtual_CI(dbl_mb(eig_eci(1)))

c      if (oprint) then
c         do ii=1,(ne_excited(1)+ne_excited(2))
c            write(*,*) "ii,eig=",ii,dbl_mb(eig_eci(1)+ii-1),
c     >                              dbl_mb(eig_excited(1)+ii-1)
c         end do
c      end if

      if (.not.BA_pop_stack(eig_eci(2)))
     >   call errquit('psi_CI_minimize_virtual:pop stack',1,MA_ERR)
      if (.not.BA_pop_stack(g93(2)))
     >   call errquit('psi_CI_minimize_virtual:pop stack',2,MA_ERR)


      return
      end 

      subroutine psi_sort_virtual_CI(eig_eci)
      implicit none
      real*8 eig_eci(*)

#include "bafdecls.fh"
#include "psi.fh"
#include "errquit.fh"

      logical value
      integer i,j,ii,jj,ms
      integer r1(2)
      real*8  ei,ej,eci,ecj

      value = BA_push_get(mt_dcpl,npack1,'rr1',r1(2),r1(1))
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
           eci = eig_eci(i)
           ecj = eig_eci(j)

           !*** swap ***
           if (ecj.lt.eci) then
             eig_eci(i) = ecj
             eig_eci(j) = eci
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

      value = BA_pop_stack(r1(2))
      if (.not. value) call errquit(
     >     'psi_sort_virtual: popping stack memory',1, MA_ERR)
      return
      end



*     ***********************************
*     *                                 *
*     *      psi_CI_update_virtual      *
*     *                                 *
*     ***********************************

*    This routine performs a KS update on virtual ii
*
      subroutine psi_CI_update_virtual(ci_algorithm,maxiteration,
     >                             maxit_ls,
     >                             maxerror,perror,ii,
     >                             error_out,e0,eci)
      implicit none
      integer ci_algorithm
      integer maxiteration,maxit_ls
      real*8  maxerror,perror
      integer ii
      real*8 error_out
      real*8 e0,eci

#include "bafdecls.fh"
#include "psi.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,done,oneloop,precondition,oprint
      integer it,pit,n2ft3d,nx,ny,nz
      real*8 eold,percent_error,error0,de0,lmbda_r0,lmbda_r1
      real*8 theta,ep,sp,e0old
      integer r1(2),t0(2),t(2),g(2),psig_r(2),psie_r(2)
      integer rhog(2),rhoe(2),rhoge(2),vcgg(2),vcge(2),vcee(2)
      integer vc(2),vctmp(2),te(2),vcc(2)
      integer psig_ptr,psie_ptr,psig_homo(2),hpsig(2),hpsie(2)
      integer hpsig_r(2),hpsie_r(2)
      real*8  Sgg,Sge,See
      real*8  H1gg,H1ge,H1ee
      real*8  H2gg,H2ge,H2ee
      real*8  Hgg,Hge,Hee
      real*8  scal1,scal2,dv,ehartr
      real*8  A,B,C,Elow,Ehigh,c1,c2,lmbda,sum,sigma

      logical  control_print
      external control_print
      real*8   control_Ep,control_Sp,lattice_omega
      external control_Ep,control_Sp,lattice_omega
      integer  electron_psi_r_ptr
      external electron_psi_r_ptr

      psig_ptr=psi1(1)
      psie_ptr=psi1_excited(1) + (ii-1)*npack1

      call Parallel_taskid(taskid)
      oprint= ((taskid.eq.MASTER).and.control_print(print_medium))

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv = scal1*lattice_omega()

      call D3dB_n2ft3d(1,n2ft3d)

      value = BA_push_get(mt_dcpl,npack1,'t0',t0(2),t0(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1,'g',g(2),g(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1,'t',t(2),t(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1,'te',te(2),te(1))

      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'psig_r',psig_r(2),psig_r(1))
      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'psie_r',psie_r(2),psie_r(1))
      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'rhog',rhog(2),rhog(1))
      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'rhoe',rhoe(2),rhoe(1))
      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'rhoge',rhoge(2),rhoge(1))

      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'vctmp',vctmp(2),vctmp(1))
      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'vc',vc(2),vc(1))

      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'vcc',vcc(2),vcc(1))
      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'vcgg',vcgg(2),vcgg(1))
      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'vcge',vcge(2),vcge(1))
      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'vcee',vcee(2),vcee(1))
cc      if (ci_algorithm.eq.5) then
         value = value.and.
     >           BA_push_get(mt_dcpl,npack1,'psig_homo',
     >                       psig_homo(2),psig_homo(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,(neq(1)+neq(2))*npack1,'hpsig',
     >                       hpsig(2),hpsig(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,(neq(1)+neq(2))*npack1,'hpsie',
     >                       hpsie(2),hpsie(1))
         value = value.and.
     >           BA_push_get(mt_dbl,(neq(1)+neq(2))*n2ft3d,'hpsig_r',
     >                       hpsig_r(2),hpsig_r(1))
         value = value.and.
     >           BA_push_get(mt_dbl,(neq(1)+neq(2))*n2ft3d,'hpsie_r',
     >                       hpsie_r(2),hpsie_r(1))
         call Pack_c_Copy(1,dcpl_mb(psig_ptr),dcpl_mb(psig_homo(1)))
cc      end if
      if (.not.value) call errquit(
     >     'psi_CI_update_virtual: out of stack memory',0, MA_ERR)

      call Pack_c_Copy(1,dcpl_mb(psig_ptr),dbl_mb(psig_r(1)))
      call Pack_c_unpack(1,dbl_mb(psig_r(1)))
      call D3dB_cr_fft3b(1,dbl_mb(psig_r(1)))
      call D3dB_r_Zero_Ends(1,dbl_mb(psig_r(1)))

      ep = control_Ep()
      sp = control_Sp()
      !precondition = .true.
      precondition = .false.
      done = .false.
      error0 = 0.0d0
      eci = 0.0d0
      !theta = 3.14159d0/600.0d0
      theta = 3.14159d0/10.0d0
      lmbda_r0 = 1.0d0
      it = 0
      pit = 0
 2    continue
         it = it + 1

        if (ci_algorithm.eq.2) then
*        *** calculate residual (steepest descent) direction for psie ***
         call psi_2x2_virtual_gradient(dcpl_mb(psig_ptr),
     >                        dbl_mb(psig_r(1)),
     >                        dcpl_mb(psie_ptr),
     >                        dbl_mb(psie_r(1)),
     >                        eci,dcpl_mb(g(1)),
     >                        dcpl_mb(t(1)),
     >                        dbl_mb(rhog(1)),
     >                        dbl_mb(vcc(1)),
     >                        dbl_mb(vcee(1)),dbl_mb(vcge(1)),
     >                        dcpl_mb(hpsig(1)),
     >                        dbl_mb(hpsig_r(1)),dbl_mb(hpsie_r(1)))
         else if (ci_algorithm.eq.3) then
         call psi_3x3_virtual_gradient(dcpl_mb(psig_ptr),
     >                        dbl_mb(psig_r(1)),
     >                        dcpl_mb(psie_ptr),
     >                        dbl_mb(psie_r(1)),
     >                        eci,dcpl_mb(g(1)),
     >                        dcpl_mb(t(1)),
     >                        dcpl_mb(te(1)),
     >                        dbl_mb(rhog(1)),
     >                        dbl_mb(rhoge(1)),
     >                        dbl_mb(rhoe(1)),
     >                        dbl_mb(vctmp(1)),
     >                        dbl_mb(vc(1)),
     >                        dbl_mb(vcc(1)),
     >                        dbl_mb(vcgg(1)),
     >                        dbl_mb(vcge(1)),
     >                        dbl_mb(vcee(1)),
     >                        dbl_mb(hpsig_r(1)),
     >                        dbl_mb(hpsie_r(1)))
         else if (ci_algorithm.eq.4) then
         call psi_4x4_virtual_gradient(dcpl_mb(psig_ptr),
     >                                 dbl_mb(psig_r(1)),
     >                                 dcpl_mb(psie_ptr),
     >                                 dbl_mb(psie_r(1)),
     >                                 eci,dcpl_mb(g(1)),
     >                                 dcpl_mb(t(1)),
     >                                 dcpl_mb(te(1)),
     >                                 dbl_mb(rhog(1)),
     >                                 dbl_mb(rhoge(1)),
     >                                 dbl_mb(rhoe(1)),
     >                                 dbl_mb(vctmp(1)),
     >                                 dbl_mb(vc(1)),
     >                                 dbl_mb(vcgg(1)),
     >                                 dbl_mb(vcge(1)),
     >                                 dbl_mb(vcee(1)))
         else if (ci_algorithm.eq.5) then
         call psi_2x2ne_virtual_gradient(dcpl_mb(psig_ptr),
     >                          dbl_mb(electron_psi_r_ptr()),
     >                          dcpl_mb(psie_ptr),
     >                          dbl_mb(psie_r(1)),
     >                          dcpl_mb(psig_homo(1)),
     >                          dbl_mb(psig_r(1)),
     >                          dcpl_mb(psie_ptr),
     >                          dbl_mb(psie_r(1)),
     >                          eci,dcpl_mb(g(1)),
     >                          dcpl_mb(t(1)),
     >                          dbl_mb(rhog(1)),
     >                          dbl_mb(vcc(1)),
     >                          dbl_mb(vcee(1)),dbl_mb(vcge(1)),
     >                          dcpl_mb(hpsig(1)),
     >                          dbl_mb(hpsig_r(1)),dbl_mb(hpsie_r(1)))
         else if (ci_algorithm.eq.6) then
         call psi_3x3ne_virtual_gradient(dcpl_mb(psig_ptr),
     >                          dbl_mb(electron_psi_r_ptr()),
     >                          dcpl_mb(psie_ptr),
     >                          dbl_mb(psie_r(1)),
     >                          dcpl_mb(psig_homo(1)),
     >                          dbl_mb(psig_r(1)),
     >                          dcpl_mb(psie_ptr),
     >                          dbl_mb(psie_r(1)),
     >                          eci,dcpl_mb(g(1)),
     >                          dcpl_mb(t(1)),
     >                          dcpl_mb(te(1)),
     >                          dbl_mb(rhog(1)),
     >                          dbl_mb(rhoge(1)),
     >                          dbl_mb(rhoe(1)),
     >                          dbl_mb(vctmp(1)),
     >                          dbl_mb(vc(1)),
     >                          dbl_mb(vcc(1)),
     >                          dbl_mb(vcgg(1)),
     >                          dbl_mb(vcge(1)),
     >                          dbl_mb(vcee(1)),
     >                          dcpl_mb(hpsig(1)),
     >                          dcpl_mb(hpsie(1)),
     >                          dbl_mb(hpsig_r(1)),
     >                          dbl_mb(hpsie_r(1)))
         end if


         eold = eci
c         write(*,*) "it,eold=",it,eold

c         if (oprint) 
c     >   write(*,*) "orb,it,eci,eold,eci-eold=",ii,it,eci,eold,eci-eold

         !*** project out to lower orbitals  ****
         call Pack_c_Copy(1,dcpl_mb(g(1)),dcpl_mb(r1(1)))

*        **** preconditioning ****
         if (precondition) then
            pit = pit + 1
            call ke_Precondition(npack1,1,
     >                           dcpl_mb(psie_ptr),
     >                           dcpl_mb(r1(1)))
         end if

         !call psi_project_out_virtual1(ii,dcpl_mb(r1(1)))
         call psi_project_out_virtual(ii,dcpl_mb(r1(1)))

*        *** determine conjuagate direction ***
         call Pack_cc_dot(1,dcpl_mb(r1(1)),
     >                   dcpl_mb(r1(1)),
     >                   lmbda_r1)

         call Pack_c_Copy(1,dcpl_mb(r1(1)),dcpl_mb(t(1)))

         if (it.gt.1) then
            call Pack_cc_daxpy(1,(lmbda_r1/lmbda_r0),
     >                      dcpl_mb(t0(1)),
     >                      dcpl_mb(t(1)))
         end if
         lmbda_r0 = lmbda_r1
         oneloop = .true.
 3       call Pack_c_Copy(1,dcpl_mb(t(1)),dcpl_mb(t0(1)))

*        *** normalize search direction, t ****
         call psi_project_out_virtual(ii,dcpl_mb(t(1)))
         call Pack_cc_dot(1,dcpl_mb(t(1)),dcpl_mb(t(1)),sigma)
         sigma = dsqrt(sigma)
         de0 = 1.0d0/sigma
         call Pack_c_SMul1(1,de0,dcpl_mb(t(1)))

*        **** compute de0 = <t|g> ****
         call Pack_cc_dot(1,dcpl_mb(t(1)),
     >                      dcpl_mb(g(1)),
     >                      de0)

*        *** bad direction ***
         if ((de0.lt.0.0d0).and.oneloop) then
           call Pack_c_Copy(1,dcpl_mb(g(1)),dcpl_mb(t(1)))
           oneloop = .false.
           go to 3
         end if

         call psi_linesearch_CI_virtual(ci_algorithm,maxit_ls,
     >                                   ii,theta,eci,de0,dcpl_mb(t(1)),
     >                                  dcpl_mb(psig_ptr),
     >                                  dbl_mb(psig_r(1)),
     >                                  dbl_mb(psie_r(1)),
     >                                  dcpl_mb(r1(1)),dcpl_mb(te(1)),
     >                                  dbl_mb(rhog(1)),
     >                                  dbl_mb(rhoe(1)),
     >                                  dbl_mb(rhoge(1)),
     >                                  dbl_mb(vctmp(1)),
     >                                  dbl_mb(vc(1)),
     >                                  dbl_mb(vcc(1)),
     >                                  dbl_mb(vcgg(1)),
     >                                  dbl_mb(vcee(1)),
     >                                  dbl_mb(vcge(1)),
     >                                  psig_homo(1),
     >                                  electron_psi_r_ptr(),
     >                                  hpsig(1),
     >                                  hpsie(1),
     >                                  hpsig_r(1),hpsie_r(1),
     >                                  maxerror)

         !precondition = (dabs(eci-eold).gt.(sp*maxerror))

         if (oprint) 
     >   write(*,*) "    orb,it,eci,eold,eci-eold=",
     >              ii,it,eci,eold,eci-eold,precondition

         done = ((it.gt.maxiteration) 
     >           .or.
     >           (dabs(eci-eold).lt.maxerror))

         !write(*,*) "done=",done,it,maxiteration,dabs(eci-eold),maxerror

         if (done) go to 4
      go to 2

 4    e0old = e0
      call psi_get_gradient_virtual(ii,dcpl_mb(g(1)))
      call Pack_cc_dot(1,dcpl_mb(psie_ptr),dcpl_mb(g(1)),e0)
      e0 = e0
      error_out = dabs(e0-e0old)

*     **** release stack memory ****
      value = .true.
cc      if (ci_algorithm.eq.5) then
         value = BA_pop_stack(hpsie_r(2))
         value = value.and.BA_pop_stack(hpsig_r(2))
         value = value.and.BA_pop_stack(hpsie(2))
         value = value.and.BA_pop_stack(hpsig(2))
         value = value.and.BA_pop_stack(psig_homo(2))
cc      end if
      value = value.and.BA_pop_stack(vcee(2))
      value = value.and.BA_pop_stack(vcge(2))
      value = value.and.BA_pop_stack(vcgg(2))
      value = value.and.BA_pop_stack(vcc(2))
      value = value.and.BA_pop_stack(vc(2))
      value = value.and.BA_pop_stack(vctmp(2))
      value = value.and.BA_pop_stack(rhoge(2))
      value = value.and.BA_pop_stack(rhoe(2))
      value = value.and.BA_pop_stack(rhog(2))
      value = value.and.BA_pop_stack(psie_r(2))
      value = value.and.BA_pop_stack(psig_r(2))
      value = value.and.BA_pop_stack(te(2))
      value = value.and.BA_pop_stack(t(2))
      value = value.and.BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(r1(2))
      value = value.and.BA_pop_stack(t0(2))
      if (.not. value) call errquit(
     >     'psi_2x2_update_virtual: popping stack memory',1, MA_ERR)

      if (oprint) then
       write(luout,921) ii,e0,error_out,eci,dabs(eci-eold),it,pit,ep,sp
  921  format(5x,"orbital",I4," current e=",E14.6," error=",E9.3,
     >       " (eci=",E17.9,
     >       " error=",E9.3,")",
     >       " iterations",I4,"(",I4,
     >       " preconditioned, Ep,Sp=",F5.1,F7.1,")")
      end if


      return
      end


*     ***********************************
*     *				        *
*     *	     psi_linesearch_CI_virtual  *
*     *				        *
*     ***********************************

*    This routine performs a linesearch on orbital ii, in the direction t.  
* This routine is needed for a KS minimizer.
*  e0 = <orb|g>
*  de0 = 2*<t|g>
*

      subroutine psi_linesearch_CI_virtual(ci_algorithm,maxit,
     >                                     ii,theta,e0,de0,t,
     >                                     psig,psig_r,psie_r,
     >                                     h1psig,h1psie,
     >                                     rho,rhoge,rhoee,
     >                                     vctmp,vc,vcc,
     >                                     vcgg,vcee,vcge,
     >                                     psig_homo_ptr,psiall_r_ptr,
     >                                     hpsig_ptr,
     >                                     hpsie_ptr,
     >                                     hpsig_r_ptr,hpsie_r_ptr,
     >                                     maxerror)
      implicit none
      integer ci_algorithm,maxit
      integer ii
      real*8  theta
      real*8  e0,de0
      complex*16  t(*) !search direction
      complex*16  psig(*)
      real*8      psig_r(*)
      real*8      psie_r(*)
      complex*16  h1psig(*),h1psie(*)
      real*8      rho(*),rhoge(*),rhoee(*)
      real*8      vctmp(*),vc(*)
      complex*16  vcc(*)
      real*8      vcgg(*),vcee(*),vcge(*)
      integer     psig_homo_ptr,psiall_r_ptr
      integer     hpsig_ptr,hpsie_ptr
      integer     hpsig_r_ptr,hpsie_r_ptr
      real*8      maxerror

#include "errquit.fh"
#include "bafdecls.fh"
#include "psi.fh"
   
*     **** local variables ****
      !integer maxit
      !parameter (maxit=15)

      logical value,done,first
      integer orb(2),g(2),psi_ptr,it,iteration
      real*8 x,y,pi,dtheta_min,e1,de1,de1a,de1b
      real*8 theta0,theta1,theta2,dtheta,emin,theta_min,demin,tnorm
      real*8 ea,eb,ec,e2,de2
      real*8 emid,demid,thetatmp

      psi_ptr=psi1_excited(1)

      pi = 4.0d0*datan(1.0d0)
      dtheta_min = 0.01*theta

*     **** allocate stack memory ****
      value = BA_push_get(mt_dcpl,npack1,'orb',
     >                       orb(2),orb(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'gg2',
     >                       g(2),g(1))
      if (.not. value) call errquit(
     >     'psi_linesearch_CI_virtual: out of stack memory',0, MA_ERR)
 

      call Pack_c_Copy(1,dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                   dcpl_mb(orb(1)))

      de0 = -de0
      call Pack_cc_dot(1,dcpl_mb(orb(1)),t,tnorm)

      theta0 = 0.0d0
      iteration = 0
c      write(*,*) 
c      write(*,*) "e0,de0, <psi|t> =",e0,de0,tnorm
c      write(*,*) "start line search:"
      emin = e0
      demin = de0
      dtheta = (2.0d0*pi/100.0d0)/99.0
      first = .true.

c      if (ii.eq.2) then
c         write(*,*) 
c         write(*,*) "e0,de0, <psi|t> =",e0,de0,tnorm
c         write(*,*) "start line search:"
c         dtheta = 0.05d0
c         thetatmp = theta 
c         theta = 0.0d0
c      end if

 10   iteration = iteration + 1
      x = cos(theta)
      y = sin(theta)
      call Pack_c_SMul(1,x,
     >                  dcpl_mb(orb(1)),
     >                  dcpl_mb(psi_ptr+(ii-1)*npack1))
       call Pack_cc_daxpy(1,y,
     >                   t,
     >                   dcpl_mb(psi_ptr+(ii-1)*npack1))

*     *** determine e1,g1, and de1/dtheta ***
*     **** de/dtheta = -sin(theta)*<orb|g> + cos(theta)*<t|g> ****

      if (ci_algorithm.eq.2) then
      call psi_2x2_virtual_gradient(psig,psig_r,
     >                              dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                              psie_r,
     >                              e1,dcpl_mb(g(1)),
     >                              h1psie,rho,vcc,
     >                              vcee,vcge,
     >                              dcpl_mb(hpsig_ptr),
     >                              dbl_mb(hpsig_r_ptr),
     >                              dbl_mb(hpsie_r_ptr))
      else if (ci_algorithm.eq.3) then
      call psi_3x3_virtual_gradient(psig,psig_r,
     >                              dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                              psie_r,
     >                              e1,dcpl_mb(g(1)),
     >                              h1psig,h1psie,
     >                              rho,rhoge,rhoee,
     >                              vctmp,vc,vcc,
     >                              vcgg,vcge,vcee,
     >                              dbl_mb(hpsig_r_ptr),
     >                              dbl_mb(hpsie_r_ptr))
      else if (ci_algorithm.eq.4) then
      call psi_4x4_virtual_gradient(psig,psig_r,
     >                              dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                              psie_r,
     >                              e1,dcpl_mb(g(1)),
     >                              h1psig,h1psie,
     >                              rho,rhoge,rhoee,
     >                              vctmp,vc,
     >                              vcgg,vcge,vcee)

      else if (ci_algorithm.eq.5) then
      call psi_2x2ne_virtual_gradient(psig,dbl_mb(psiall_r_ptr),
     >                              dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                              psie_r,
     >                              dcpl_mb(psig_homo_ptr),psig_r,
     >                              dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                              psie_r,
     >                              e1,dcpl_mb(g(1)),
     >                              h1psie,rho,vcc,
     >                              vcee,vcge,
     >                              dcpl_mb(hpsig_ptr),
     >                              dbl_mb(hpsig_r_ptr),
     >                              dbl_mb(hpsie_r_ptr))
      else if (ci_algorithm.eq.6) then
      call psi_3x3ne_virtual_gradient(psig,dbl_mb(psiall_r_ptr),
     >                              dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                              psie_r,
     >                              dcpl_mb(psig_homo_ptr),psig_r,
     >                              dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                              psie_r,
     >                              e1,dcpl_mb(g(1)),
     >                              h1psig,h1psie,
     >                              rho,rhoge,rhoee,
     >                              vctmp,vc,vcc,
     >                              vcgg,vcge,vcee,
     >                              dcpl_mb(hpsig_ptr),
     >                              dcpl_mb(hpsie_ptr),
     >                              dbl_mb(hpsig_r_ptr),
     >                              dbl_mb(hpsie_r_ptr))
      end if
      call Pack_cc_dot(1,dcpl_mb(orb(1)),dcpl_mb(g(1)),de1a)
      call Pack_cc_dot(1,t,              dcpl_mb(g(1)),de1b)
      de1 = y*de1a - x*de1b
c      write(*,*) "theta,e0,e1,de0,de1 = ",
c     >           ii,iteration,theta,e0,e1,de0,de1

c      if (e1.lt.emin) then
c         emin  = e1
c         demin = de1
c         theta_min = theta
c      end if

c      if ((first).and.(ii.eq.2)) then
c         write(*,*) "theta,e0,e1,de0,de1 = ",
c     >           ii,iteration,theta,e0,e1,de0,de1
c         theta = theta + dtheta
c         if (theta.gt.2.0d0*pi) then
c            first = .false.
c           dtheta = (2.0d0*pi/100.0d0)/99.0
c           theta = thetatmp
c        
c         end if
c
c         go to 10
c      end if

      if (e1.lt.emin) then
         emin  = e1
         demin = de1
         theta_min = theta
      end if

      done = (de0.le.0.0d0).and.(de1.ge.0.0d0).or.(iteration.gt.maxit)

      if (.not.done) then
         theta = theta + theta     
         go to 10
      end if

      theta0 = 0.0d0
      theta1 = theta

      theta = theta1 - de0*(theta1-theta0)/(de1-de0)

cccccc 20   theta = 0.5d0*(theta0+theta1)
 20   iteration = iteration + 1
      theta1 = theta
      x = cos(theta1)
      y = sin(theta1)
      call Pack_c_SMul(1,x,
     >                  dcpl_mb(orb(1)),
     >                  dcpl_mb(psi_ptr+(ii-1)*npack1))
       call Pack_cc_daxpy(1,y,
     >                   t,
     >                   dcpl_mb(psi_ptr+(ii-1)*npack1))

*     *** determine e1,g1, and de1/dtheta ***
*     **** de/dtheta = -sin(theta)*<orb|g> + cos(theta)*<t|g> ****
      if (ci_algorithm.eq.2) then
      call psi_2x2_virtual_gradient(psig,psig_r,
     >                              dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                              psie_r,
     >                              e1,dcpl_mb(g(1)),
     >                              h1psie,rho,vcc,
     >                              vcee,vcge,
     >                              dcpl_mb(hpsig_ptr),
     >                              dbl_mb(hpsig_r_ptr),
     >                              dbl_mb(hpsie_r_ptr))
      else if (ci_algorithm.eq.3) then
      call psi_3x3_virtual_gradient(psig,psig_r,
     >                              dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                              psie_r,
     >                              e1,dcpl_mb(g(1)),
     >                              h1psig,h1psie,
     >                              rho,rhoge,rhoee,
     >                              vctmp,vc,vcc,
     >                              vcgg,vcge,vcee,
     >                              dbl_mb(hpsig_r_ptr),
     >                              dbl_mb(hpsie_r_ptr))
      else if (ci_algorithm.eq.4) then
      call psi_4x4_virtual_gradient(psig,psig_r,
     >                              dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                              psie_r,
     >                              e1,dcpl_mb(g(1)),
     >                              h1psig,h1psie,
     >                              rho,rhoge,rhoee,
     >                              vctmp,vc,
     >                              vcgg,vcge,vcee)
      else if (ci_algorithm.eq.5) then
      call psi_2x2ne_virtual_gradient(psig,dbl_mb(psiall_r_ptr),
     >                              dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                              psie_r,
     >                              dcpl_mb(psig_homo_ptr),psig_r,
     >                              dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                              psie_r,
     >                              e1,dcpl_mb(g(1)),
     >                              h1psie,rho,vcc,
     >                              vcee,vcge,
     >                              dcpl_mb(hpsig_ptr),
     >                              dbl_mb(hpsig_r_ptr),
     >                              dbl_mb(hpsie_r_ptr))
      else if (ci_algorithm.eq.6) then
      call psi_3x3ne_virtual_gradient(psig,dbl_mb(psiall_r_ptr),
     >                              dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                              psie_r,
     >                              dcpl_mb(psig_homo_ptr),psig_r,
     >                              dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                              psie_r,
     >                              e1,dcpl_mb(g(1)),
     >                              h1psig,h1psie,
     >                              rho,rhoge,rhoee,
     >                              vctmp,vc,vcc,
     >                              vcgg,vcge,vcee,
     >                              dcpl_mb(hpsig_ptr),
     >                              dcpl_mb(hpsie_ptr),
     >                              dbl_mb(hpsig_r_ptr),
     >                              dbl_mb(hpsie_r_ptr))
      end if
      call Pack_cc_dot(1,dcpl_mb(orb(1)),dcpl_mb(g(1)),de1a)
      call Pack_cc_dot(1,t,              dcpl_mb(g(1)),de1b)
      de1 = y*de1a - x*de1b
      !dtheta =  -de0*(theta-theta0)/(demid-de0)

      !write(*,*) "theta,e0,e1,(e1-e0)=",iteration,theta,e0,e1,e1-e0
      !write(*,*) "            de0,de1=",de0,de1
      if (e1.lt.emin) then
         emin  = e1
         demin = de1
         theta_min = theta1
      end if

      if (e1.le.e0) then
         theta  = theta1 - de0*(theta1-theta0)/(de1-de0)
         e0     = e1
         de0    = de0
         theta0 = theta1
         done = .false.
      else if (e1.gt.e0) then
         theta = theta0 + 0.5d0*(theta1-theta0)
         done = .false.
      end if

      done = (dabs(e1-e0).lt.maxerror).or.(iteration.gt.maxit)

      if (.not.done) go to 20




c         theta = -pi/100.0 + (it-1)*(2.0d0*pi/100.0d0)/99.0
c
c*        **** orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
c  10     x = cos(theta)
c         y = sin(theta)
c         call Pack_c_SMul(1,x,
c     >                  dcpl_mb(orb(1)),
c     >                  dcpl_mb(psi_ptr+(ii-1)*npack1))
c         call Pack_cc_daxpy(1,y,
c     >                   t,
c     >                   dcpl_mb(psi_ptr+(ii-1)*npack1))
c
c*        *** determine theta ***
c         call psi_2x2_virtual_gradient(psig,psig_r,
c     >                              dcpl_mb(psi_ptr+(ii-1)*npack1),
c     >                              psie_r,
c     >                              e1,dcpl_mb(g(1)),
c     >                              h1psie,rho,vcee,vcge)
c
c*        **** compute de/dtheta = -sin(theta)*<orb|g> + cos(theta)*<t|g> ****
c         call Pack_cc_dot(1,dcpl_mb(orb(1)),dcpl_mb(g(1)),de1a)
c         call Pack_cc_dot(1,t,              dcpl_mb(g(1)),de1b)
c         de1 = y*de1a - x*de1b
c
c         if (e1<emin) then
c            emin  = e1
c            demin = de1
c            theta_min = theta
c         end if
c
c
c         write(*,*) theta,e1,de1,(e1-ea)/(2.0d0*dtheta)
c         ea = eb
c         eb = e1
c      end do
c      theta = theta_min
      !write(*,*) ":end line search"
      !write(*,*)
      !write(*,*) "theta_min,emin,demin = ",theta_min,emin,demin
      !write(*,*)

c*     **** compute de/dtheta = -sin(theta)*<orb|g> + cos(theta)*<t|g> ****
c      call Pack_cc_dot(1,dcpl_mb(orb(1)),dcpl_mb(g(1)),de1a)
c      call Pack_cc_dot(1,t,              dcpl_mb(g(1)),de1b)
c      de1 = y*de1a - x*de1b
c
c      write(*,*) "GERDD,ii,it,theta,e0,e1=",ii,it,theta,e0,e1,de0,de1
c     >                      
c
c      it = it + 1
c      if ((e1>e0).and.(it<5)) then
c         theta = 0.5d0*theta
c         goto 10
c      end if
c      dtheta = de0*(theta-0.0d0)/(de1-de0)
c      theta = theta + dtheta
c      write(*,*) "NEW theta = ",theta
       


c      x = (e0 - e1 + 0.5d0*de0*sin(2*theta))
c     >    /(1.0d0-cos(2*theta))
c      theta = 0.5d0*datan(0.5d0*de0/x) 
c      !write(*,*) "GERE,x,theta,e0,e1=",x,theta,e0,e1

c*     **** orb2 = orb*cos(theta) + t*sin(theta) ****
c      x = cos(theta)
c      y = sin(theta)
c      call Pack_c_SMul(1,x,
c     >                  dcpl_mb(orb(1)),
c     >                  dcpl_mb(psi_ptr+(ii-1)*npack1))
c      call Pack_cc_daxpy(1,y,
c     >                   t,
c     >                   dcpl_mb(psi_ptr+(ii-1)*npack1))


*     **** release stack memory ****
      value =           BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(orb(2))      
      if (.not. value) call errquit(
     >     'psi_linesearch_CI_virtual: popping stack memory',1, MA_ERR)

      return
      end

c      subroutine psi_follow_geovirt(orb1,t,orb2)
c      implicit none
c      complex*16 orb(*),t(*),orb2(*)
c
c      integer i
c      real*8 theta x,y
c
c      do i=1,50
c         theta = (i-1)/50.0 * 2.0d0*datan(1.0d0)
c         call psi_linesearch_genpsi(theta,orb,t,orb2)
c
c         call psi_2x2_virtual_gradient(psig,psig_r,
c     >                              orb2,
c     >                              psie_r,
c     >                              e1,dcpl_mb(g(1)),
c     >                              h1psie,rho,vcee,vcge)
c      end do

c      return
c      end


      subroutine psi_linesearch_genpsi(theta,orb,t,orb2)
      implicit none
      real*8 theta
      complex*16 orb(*),t(*),orb2(*)
      real*8 x,y

      x = cos(theta)
      y = sin(theta)
      call Pack_c_SMul(1,x,orb,orb2)
      call Pack_cc_daxpy(1,y,t,orb2)
      return
      end


*     ***********************************
*     *                                 *
*     *      psi_2x2_virtual_gradient   *
*     *                                 *
*     ***********************************
*
*    This routine calculates the 2x2 CI Energy and its gradient wrt to psie.
*
*    Entry - psig,psig_r
*            psie 
*
*    Exit - E,dEpsie
*
*    Use - psie_r
*          h1psie,rho,vcee,vcge,hpsig,hpsig_r,hpsie_r
* 
      subroutine psi_2x2_virtual_gradient(psig,psig_r,psie,psie_r,
     >                                    E,dEdpsie,
     >                                    h1psie,rho,vcc,vcee,vcge,
     >                                    hpsig,
     >                                    hpsig_r,hpsie_r)

      implicit none
      complex*16 psig(*)
      real*8     psig_r(*)
      complex*16 psie(*)
      real*8     psie_r(*)
      real*8     E
      complex*16 dEdpsie(*)

      complex*16 h1psie(*)
      complex*16 vcc(*)
      real*8     rho(*),vcee(*),vcge(*)
      complex*16 hpsig(*)
      real*8     hpsig_r(*)
      real*8     hpsie_r(*)

#include "bafdecls.fh"
#include "psi.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"

*     **** local variables ****
      logical oprint
      integer neq1(2),i,q,ms,n1q(2),n2q(2)
      integer nx,ny,nz,n2ft3d
      real*8 Sgg,Sge,See,A,B,C,c1,c2
      real*8 H1gg,H1ge,H1ee
      real*8 H2gg,H2ge,H2ee
      real*8 Hgg,Hge,Hee
      real*8 ehartr,scal1,scal2,dv,Elow,Ehigh,lmbda,sum1
      real*8 ehfx,phfx,Eion
      integer  control_version
      external control_version
      real*8   coulomb_e,ewald_e,ion_ion_e
      external coulomb_e,ewald_e,ion_ion_e

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      !write(*,*) "ENTER 2x2 GRADIENT"
      oprint = .false.

      !*** generate psie_r ****
      call Pack_c_Copy(1,psie,psie_r)
      call Pack_c_unpack(1,psie_r)
      call D3dB_cr_fft3b(1,psie_r)
      call D3dB_r_Zero_Ends(1,psie_r)

      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv = scal1*lattice_omega()

      !*** calculate Sgg, Sge, See ***
      call Pack_cc_dot(1,psig,psig,Sgg)
      call Pack_cc_dot(1,psig,psie,Sge)
      call Pack_cc_dot(1,psie,psie,See)
      !if (oprint) write(*,*) "Sgg=",Sgg," Sge=",Sge," See=",See

*     **** ion-ion part ****
      if (control_version().eq.3) Eion = ewald_e()   !**** get ewald energy ****
      if (control_version().eq.4) Eion = ion_ion_e() !**** get free-space ion-ion energy ****
c     if (oprint) write(*,*) "Eion=",Eion
      write(*,*) "Eion=",Eion

*     ***************************
*     **** one-electron part ****
*     ***************************
*     **** apply H1 operator ****
      call  psi_H1psi(ispin,neq,npack1,n2ft3d,psig,psig_r,dEdpsie)
      call  psi_H1psi(ispin,neq,npack1,n2ft3d,psie,psie_r,h1psie)
 
      call Pack_cc_dot(1,psig,dEdpsie,H1gg)
      call Pack_cc_dot(1,psig,h1psie,H1ge)
      call Pack_cc_dot(1,psie,h1psie,H1ee)
      H1gg = -H1gg
      H1ge = -H1ge
      H1ee = -H1ee
      if (ispin.eq.1) then
         H1gg = H1gg + H1gg
         H1ge = H1ge + H1ge
         H1ee = H1ee + H1ee
      end if
      H1ge = H1ge*Sge

      if (oprint) write(*,*) "H1gg=",H1gg," H1ge=",H1ge," H1ee=",H1ee
      
*     **** apply H2 operator  - note rho = 2*psi*psi ****
      call pspw_et_gen_rho(ispin,neq,n2ft3d,psig_r,rho)
      if (control_version().eq.4) then
         call coulomb2_v(rho,vcee)
         call D3dB_rr_dot(1,rho,vcee,ehartr)
         H2gg = 0.5d0*(0.5d0*ehartr*dv) 
         if (oprint) write(*,*) "Coulombgg=",2.0d0*H2gg,
     >                          "Exchangegg=",-H2gg
      else
         call D3dB_r_SMul1(1,scal1,rho)
         call D3dB_rc_fft3f(1,rho)
         call Pack_c_pack(0,rho)
         H2gg = coulomb_e(rho)
         call D3dB_r_nZero(1,neq(1)+neq(2),hpsig_r)
         call pspw_potential_HFX(ispin,psig_r,hpsig_r)
         call pspw_energy_HFX(ispin,psig_r,ehfx,phfx)
         if (oprint) write(*,*) "Coulombgg=",H2gg,
     >                          "Exchangegg=",ehfx
         H2gg = H2gg + ehfx
      end if


      call pspw_et_gen_rho(ispin,neq,n2ft3d,psie_r,rho)
      if (control_version().eq.4) then
         call coulomb2_v(rho,vcee)
         call D3dB_rr_dot(1,rho,vcee,ehartr)
         H2ee = 0.5d0*(0.5d0*ehartr*dv)
         if (oprint) write(*,*) "Coulombee=",2.0d0*H2ee,
     >                          "Exchangeee=",-H2ee
      else
         call D3dB_r_SMul1(1,scal1,rho)
         call D3dB_rc_fft3f(1,rho)
         call Pack_c_pack(0,rho)
         H2ee = coulomb_e(rho)
         call coulomb_v(rho,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_c_Copy(1,vcc,vcee)
         call D3dB_r_nZero(1,neq(1)+neq(2),hpsig_r)
         call pspw_potential_HFX(ispin,psie_r,hpsig_r)
         call pspw_energy_HFX(ispin,psie_r,ehfx,phfx)
         if (oprint) write(*,*) "Coulombee=",H2ee,
     >                          "Exchangeee=",ehfx
         H2ee = H2ee + ehfx
         call D3dB_rr_Mul2(1,psie_r,vcee)
         call D3dB_rr_Sum2(1,hpsig_r,vcee) 
      end if 

      call pspw_et_gen_rho12(ispin,neq,n2ft3d,psig_r,psie_r,rho)
      if (control_version().eq.4) then
         call coulomb2_v(rho,vcge)
         call D3dB_rr_dot(1,rho,vcge,ehartr)
         H2ge = 0.5d0*(0.5d0*ehartr*dv)
         if (oprint) write(*,*) "Coulombge=",2.0d0*H2ge,
     >                          "Exchangege=",-H2ge
      else
         call D3dB_r_SMul1(1,scal1,rho)
         call D3dB_rc_fft3f(1,rho)
         call Pack_c_pack(0,rho)
c        H2ge = coulomb_e(rho)
c        call coulomb_v(rho,vcc)
c        call Pack_c_unpack(0,vcc)
c        call D3dB_cr_fft3b(1,vcc)
c        call D3dB_c_Copy(1,vcc,vcge)
         call D3dB_r_nZero(1,neq(1)+neq(2),hpsig_r)
         call D3dB_r_nZero(1,neq(1)+neq(2),hpsie_r)
         call pspw_potential_HFX2_dnc(ispin,psig_r,psie_r,
     >                                   hpsig_r,hpsie_r,
     >                                   ehfx,phfx)
c        if (oprint) write(*,*) "HFX2_dnc:ehfx,phfx",ehfx,phfx
c        call Pack_cc_dot(1,hpsig_r,hpsig_r,c1)
c        call Pack_cc_dot(1,hpsie_r,hpsie_r,c2)
c        if (oprint) write(*,*) "HFX2_dnc:hspi",c1,c2
c
c        call D3dB_r_nZero(1,neq(1)+neq(2),hpsig_r)
c        call D3dB_r_nZero(1,neq(1)+neq(2),hpsie_r)
c        call pspw_potential_HFX2(ispin,psig_r,psie_r,hpsig_r,hpsie_r)
c        call pspw_energy_HFX2(ispin,psig_r,psie_r,ehfx,phfx)
c        if (oprint) write(*,*) "HFX2:ehfx,phfx",ehfx,phfx
c        call Pack_cc_dot(1,hpsig_r,hpsig_r,c1)
c        call Pack_cc_dot(1,hpsie_r,hpsie_r,c2)
c        if (oprint) write(*,*) "HFX2:hspi",c1,c2
         if (oprint) write(*,*) "Coulombge=",H2ge,
     >                          "Exchangege=",ehfx
c        H2ge = H2ge + ehfx 
         H2ge = -ehfx 
      end if 

      if (oprint) write(*,*) "H2gg=",H2gg," H2ge=",H2ge," H2ee=",H2ee
*     **** generate and diagonalize 2x2 CI matrix ****
      Hgg = H1gg + H2gg + Eion
      Hge = H1ge + H2ge 
      Hee = H1ee + H2ee + Eion
c     if (oprint) write(*,*) "Hgg=",Hgg," Hge=",Hge," Hee=",Hee
      write(*,*) "Hgg=",Hgg," Hge=",Hge," Hee=",Hee
      A = 1.0d0
      B = -Hgg - Hee
      C = Hgg*Hee - Hge*Hge
      Elow  = (-B-dsqrt(B*B-4.0d0*A*C))/(2*A)
      Ehigh = (-B+dsqrt(B*B-4.0d0*A*C))/(2*A)
      lmbda = Hge/(Elow-Hee)
      c1 =  1.0d0/dsqrt(1.0d0+lmbda**2)
      !c2 =  dsqrt(1.0d0-c1*c1)
      c2 =  lmbda/dsqrt(1.0d0+lmbda**2)
      if (oprint) write(*,*) "Elow=",Elow," Ehigh=",Ehigh
      if (oprint) write(*,*) "c1,c2=",c1,c2

*     **** generate dElow/dpsie ****
      if (control_version().eq.4) then
         call D3dB_rr_Mul2(1,psie_r,vcee)
         call D3dB_rr_Mul2(1,psig_r,vcge)
         call D3dB_r_SMul1(1,(-1.0d0*c2*c2),vcee)
         call D3dB_rr_daxpy(1,(-1.0d0*c1*c2),vcge,vcee)
      else
c        call D3dB_rr_Mul2(1,psig_r,vcge)
c        call D3dB_rr_Sum2(1,hpsie_r,vcge) 
         call D3dB_r_SMul1(1,(-2.0d0*c2*c2),vcee)
c        call D3dB_rr_daxpy(1,(-2.0d0*c1*c2),vcge,vcee)
         call D3dB_rr_daxpy(1,(2.0d0*c1*c2),hpsie_r,vcee)
      end if

      call D3dB_rc_fft3f(1,vcee)
      call Pack_c_pack(1,vcee)
      call Pack_c_SMul(1,scal1,vcee,dEdpsie)
     
      call Pack_cc_daxpy(1,2.0d0*c2*c2,h1psie,dEdpsie)
      E = Elow

      return
      end

*     ***********************************
*     *                                 *
*     *    psi_2x2ne_virtual_gradient   *
*     *                                 *
*     ***********************************
*
*    This routine calculates the 2x2 CI Energy and its gradient wrt to psie.
*
*    Entry - psig,psig_r
*            psie 
*
*    Exit - E,dEpsie
*
*    Use - psie_r
*          h1psie,rho,vcee,vcge
* 
      subroutine psi_2x2ne_virtual_gradient(psig,psig_r,
     >                                    psie,psie_r,
     >                                    psihomo,psihomo_r,
     >                                    psilumo,psilumo_r,
     >                                    E,dEdpsie,
     >                                    h1psie,rho,vcc,
     >                                    vcee,vcge,
     >                                    hpsig,
     >                                    hpsig_r,hpsie_r)
      implicit none
      complex*16 psig(*)
      real*8     psig_r(*)

      complex*16 psie(*)
      real*8     psie_r(*)

      complex*16 psihomo(*)
      real*8     psihomo_r(*)

      complex*16 psilumo(*)
      real*8     psilumo_r(*)

      real*8     E
      complex*16 dEdpsie(*)

      complex*16 h1psie(*)
      complex*16 vcc(*)

      real*8     rho(*),vcee(*),vcge(*)

      complex*16 hpsig(*)
      real*8     hpsig_r(*)
      real*8     hpsie_r(*)


#include "bafdecls.fh"
#include "psi.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"

      real*8 Hgg
      common / CI_Hgg_common / Hgg

*     **** local variables ****
      logical oprint
      integer neq1(2),i,q,ms,n1q(2),n2q(2)
      integer nx,ny,nz,n2ft3d
      real*8 Sgg,Sge,See,A,B,C,c1,c2
      real*8 H1gg,H1ge,H1ee
      real*8 H2gg,H2ge,H2ee
      real*8 Hge,Hee
      real*8 ehartr,scal1,scal2,dv,Elow,Ehigh,lmbda,sum1
      real*8 ehfx,phfx,Eion

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega
      integer  control_version
      external control_version
      real*8   coulomb_e,ewald_e,ion_ion_e
      external coulomb_e,ewald_e,ion_ion_e

      !write(*,*) "ENTER 2x2 ne GRADIENT"
      oprint = .false.

      !*** generate psilumo_r ****
      call Pack_c_Copy(1,psilumo,psilumo_r)
      call Pack_c_unpack(1,psilumo_r)
      call D3dB_cr_fft3b(1,psilumo_r)
      call D3dB_r_Zero_Ends(1,psilumo_r)

      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv = scal1*lattice_omega()

      call Pack_c_Zero(1,dEdpsie)

      !do i=1,5
      !   write(*,*) "i,psi_r,psi_homo=",psig_r(i),psihomo_r(i)
      !end do

      !*** calculate Sgg, Sge, See ***
      Sgg = 1.0d0
      See = 1.0d0
      Sge = 0.0d0
      !if (oprint) write(*,*) "Sgg=",Sgg," Sge=",Sge," See=",See
*     **** ion-ion part ****
      if (control_version().eq.3) Eion = ewald_e()   !**** get ewald energy ****
      if (control_version().eq.4) Eion = ion_ion_e() !**** get free-space ion-ion energy ****
      if (oprint) write(*,*) "Eion=",Eion
      write(*,*) "Eion=",Eion


*     ***************************
*     **** one-electron part ****
*     ***************************
      n1q(1) = 1
      n2q(1) = neq(1)
      n1q(2) = neq(1)+1
      n2q(2) = neq(1)+neq(2)

      !*** generate HF energy and gradient of state psig ****
      !*** generate Hgg if Hgg==0 ***
*     *** generate <g|H|g> ***
      if (dabs(Hgg).lt.1.0d-9) then
*        **** apply H1 operator ****
         call psi_H1psi(ispin,neq,npack1,n2ft3d,psig,psig_r,hpsig)
         H1gg=0.0d0
         do ms=1,ispin
         if (neq(ms).gt.0) then
            do q=n1q(ms),n2q(ms)
               call Pack_cc_idot(1,psig(1+npack1*(q-1)),
     >                            hpsig(1+npack1*(q-1)),sum1)
               H1gg = H1gg - sum1
            end do
         end if
         end do
         call Parallel_SumAll(H1gg)
         if (ispin.eq.1) H1gg = H1gg + H1gg

         !write(*,*) "H1gg=",H1gg

*        **** apply H2 operator  - note rho = 2*psi*psi ****
         H2gg=0.0d0
*        **** coulomb part ****
         call pspw_et_gen_rho(ispin,neq,n2ft3d,psig_r,rho)
         if (control_version().eq.4) then
            call coulomb2_v(rho,vcee)
            call D3dB_rr_dot(1,rho,vcee,ehartr)
            H2gg = (0.5d0*ehartr*dv)
*        **** exchange part ****
            call D3dB_r_nZero(1,neq(1)+neq(2),hpsig_r)
            call pspw_potential_HFX(ispin,psig_r,hpsig_r)
            call pspw_energy_HFX(ispin,psig_r,ehfx,phfx)
            if (oprint) write(*,*) "Coulombgg=",H2gg,
     >                             "Exchangegg=",ehfx
            H2gg = H2gg + ehfx
            Hgg = H1gg + H2gg + Eion
            !write(*,*) "H2gg+ehfx=",H2gg
            !write(*,*) "Hgg=",Hgg
         else
            call D3dB_r_SMul1(1,scal1,rho)
            call D3dB_rc_fft3f(1,rho)
            call Pack_c_pack(0,rho)
            H2gg = coulomb_e(rho)
            call D3dB_r_nZero(1,neq(1)+neq(2),hpsig_r)
            call pspw_potential_HFX(ispin,psig_r,hpsig_r)
            call pspw_energy_HFX(ispin,psig_r,ehfx,phfx)
            if (oprint) write(*,*) "Coulombgg=",H2gg,
     >                             "Exchangegg=",ehfx
            H2gg = H2gg + ehfx 
            Hgg = H1gg + H2gg + Eion
          end if

         !write(*,*) "H2gg=",H2gg

      end if
      !write(*,*) "Hgg=",Hgg
      !write(*,*)


      !*** generate HF energy and gradient of state psie ****
*     *** generate <e|H|e> ***
      call Pack_c_Copy(1,psilumo,psig)
      call D3dB_c_Copy(1,psilumo_r,psig_r)

*     **** apply H1 operator ****
      call psi_H1psi(ispin,neq,npack1,n2ft3d,psig,psig_r,hpsig)
      call Pack_c_Copy(1,hpsig,h1psie)

      H1ee=0.0d0
      do ms=1,ispin
         if (neq(ms).gt.0) then
            do q=n1q(ms),n2q(ms)
               call Pack_cc_idot(1,psig(1+npack1*(q-1)),
     >                            hpsig(1+npack1*(q-1)),sum1)
               H1ee = H1ee - sum1
            end do
         end if
      end do
      call Parallel_SumAll(H1ee)
      if (ispin.eq.1) H1ee = H1ee + H1ee
      !write(*,*) "H1ee=",H1ee

*     **** apply H2 operator  - note rho = 2*psi*psi ****
      H2ee=0.0d0
*     **** coulomb part ****
      call pspw_et_gen_rho(ispin,neq,n2ft3d,psig_r,rho)
      if (control_version().eq.4) then
         call coulomb2_v(rho,vcee)
         call D3dB_rr_dot(1,rho,vcee,ehartr)
         H2ee = (0.5d0*ehartr*dv)
         !write(*,*) "H2ee=",H2ee

*        **** exchange part ****
         call D3dB_r_nZero(1,neq(1)+neq(2),hpsig_r)
         call pspw_potential_HFX(ispin,psig_r,hpsig_r)
         call pspw_energy_HFX(ispin,psig_r,ehfx,phfx)
         if (oprint) write(*,*) "Coulombee=",H2ee,
     >                          "Exchangeee=",ehfx
         H2ee = H2ee + ehfx
         Hee = H1ee + H2ee + Eion
         !write(*,*) "H2ee+ehfx=",H2ee,ehfx
         !write(*,*) "Hee=",Hee
         !write(*,*)
      else
         call D3dB_r_SMul1(1,scal1,rho)
         call D3dB_rc_fft3f(1,rho)
         call Pack_c_pack(0,rho)
         H2ee = coulomb_e(rho)
         call coulomb_v(rho,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_c_Copy(1,vcc,vcee)
         call D3dB_r_nZero(1,neq(1)+neq(2),hpsig_r)
         call pspw_potential_HFX(ispin,psig_r,hpsig_r)
         call pspw_energy_HFX(ispin,psig_r,ehfx,phfx)
         if (oprint) write(*,*) "Coulombee=",H2ee,
     >                          "Exchangeee=",ehfx
         H2ee = H2ee + ehfx
         Hee = H1ee + H2ee + Eion
         call D3dB_rr_Mul2(1,psilumo_r,vcee)
         call D3dB_rr_Sum2(1,hpsig_r,vcee)
      end if

      call Pack_c_Copy(1,psihomo,psig)
      call D3dB_c_Copy(1,psihomo_r,psig_r)
      call Pack_c_Copy(1,psilumo,psie)
      call D3dB_c_Copy(1,psilumo_r,psie_r)


*     *** generate <e|H|g> ***
      H1ge = 0.0d0

      neq1(1) = 1
      neq1(2) = 0

      H2ge = 0.0d0
      call pspw_et_gen_rho12(ispin,neq1,n2ft3d,psihomo_r,psilumo_r,rho)
      if (control_version().eq.4) then
         call coulomb2_v(rho,vcge)
         call D3dB_rr_dot(1,rho,vcge,ehartr)
         H2ge = 0.5d0*(0.5d0*ehartr*dv)
         if (oprint) write(*,*) "Coulombge=",2.0d0*H2ge,
     >                          "Exchangege=",-H2ge
         Hge = H1ge + H2ge
      else
c        call D3dB_r_SMul1(1,scal1,rho)
c        call D3dB_rc_fft3f(1,rho)
c        call Pack_c_pack(0,rho)
c        H2ge = coulomb_e(rho)
c        call coulomb_v(rho,vcc)
c        call Pack_c_unpack(0,vcc)
c        call D3dB_cr_fft3b(1,vcc)
c        call D3dB_c_Copy(1,vcc,vcge)
         call D3dB_r_nZero(1,neq(1)+neq(2),hpsie_r)
         call pspw_potential_HFX2_dnc(ispin,psig_r,psie_r,
     >                                   hpsig_r,hpsie_r,
     >                                   ehfx,phfx)
c        call pspw_potential_HFX2(ispin,psig_r,psie_r,
c    >                            hpsig_r,hpsie_r)
c        call pspw_energy_HFX2(ispin,psig_r,psie_r,ehfx,phfx)
         if (oprint) write(*,*) !"Coulombge=",H2ge, 
     >                          "Exchangege=",ehfx
c        H2ge = H2ge + ehfx
         H2ge = -ehfx
         Hge = H1ge + H2ge
c        call D3dB_rr_Mul2(1,psihomo_r,vcge)
c        call D3dB_rr_Sum2(1,hpsie_r,vcge)
      end if
      
      !write(*,*) "H1ge=",H1ge
      !write(*,*) "H2ge=",H2ge
      !write(*,*) "Hge=",Hge
      !write(*,*)

      write(*,*) "H1gg,H1ge,H1ee",H1gg,H1ge,H1ee
      write(*,*) "H2gg,H2ge,H2ee",H2gg,H2ge,H2ee
      write(*,*) "Hgg,Hge,Hee",Hgg,Hge,Hee
*     **** generate and diagonalize 2x2 CI matrix ****
      A = 1.0d0
      B = -Hgg - Hee
      C = Hgg*Hee - Hge*Hge
      Elow  = (-B-dsqrt(B*B-4.0d0*A*C))/(2*A)
      Ehigh = (-B+dsqrt(B*B-4.0d0*A*C))/(2*A)
      lmbda = Hge/(Elow-Hee)
      c1 =  1.0d0/dsqrt(1.0d0+lmbda**2)
      !c2 =  dsqrt(1.0d0-c1*c1)
      c2 =  lmbda/dsqrt(1.0d0+lmbda**2)
      !write(*,*) "Elow=",Elow," Ehigh=",Ehigh
      !write(*,*)

*     **** generate dElow/dpsie ****
      if (control_version().eq.4) then

         call D3dB_rr_Mul2(1,psilumo_r,vcee)
         call D3dB_rr_Sum2(1,hpsig_r,vcee)

         call D3dB_rr_Mul2(1,psihomo_r,vcge)

         call D3dB_r_SMul1(1,(-2.0d0*c2*c2),vcee)
         call D3dB_rr_daxpy(1,(-1.0d0*c1*c2),vcge,vcee)
         !original call D3dB_rr_daxpy(1,(-1.0d0*c1*c2),vcge,vcee)
      else
         call D3dB_r_SMul1(1,(-2.0d0*c2*c2),vcee)
c        call D3dB_rr_daxpy(1,(-2.0d0*c1*c2),vcge,vcee)
         call D3dB_rr_daxpy(1,(2.0d0*c1*c2),hpsie_r,vcee)
      end if

        
      call D3dB_rc_fft3f(1,vcee)
      call Pack_c_pack(1,vcee)
      call Pack_c_SMul(1,scal1,vcee,dEdpsie)

      call Pack_cc_daxpy(1,2.0d0*c2*c2,h1psie,dEdpsie)
      E = Elow
      !if (oprint) write(*,*) "2x2ne c1,c2=",c1,c2
      write(*,*) "2x2ne c1,c2=",c1,c2,Elow

c     call Pack_cc_dot(1,psie,dEdpsie,c1)
c     call Pack_cc_dot(1,psie,h1psie,c1)
c     write(*,*) "check last gradient=",c1
c     call Pack_cc_dot(1,psig,dEdpsie,c1)
c     call Pack_cc_dot(1,psig,h1psie,c1)
c     write(*,*) "check last gradient=",c1
c     call Pack_cc_dot(1,dEdpsie,dEdpsie,c1)
c     call Pack_cc_dot(1,h1psie,h1psie,c1)
c     write(*,*) "check last gradient=",c1

c

      return
      end 





*     ***********************************
*     *                                 *
*     *      psi_3x3_virtual_gradient   *
*     *                                 *
*     ***********************************
*
*    This routine calculates the 3x3 CI Energy and its gradient wrt to psie.
*
*    Entry - psig,psig_r
*            psie 
*
*    Exit - E,dEpsie
*
*    Use - psie_r
*          h1psig,h1psie,rhogg,rhoge,rhoee,
*          vctmp,vc,vcc,vcgg,vcge,vcee,
*          hpsig_r,hpsie_r
*
*    The matrix is
*        g    e    m    
*    g  Hgg  Hge  Hgm  
*
*    e  Heg  Hee  Hem  
*
*    m  Hmg  Hme  Hmm  

      subroutine psi_3x3_virtual_gradient(psig,psig_r,psie,psie_r,
     >                                    E,dEdpsie,
     >                                    h1psig,h1psie,
     >                                    rhogg,rhoge,rhoee,
     >                                    vctmp,vc,
     >                                    vcc,
     >                                    vcgg,vcge,vcee,
     >                                    hpsig_r,hpsie_r)
      implicit none
      complex*16 psig(*)
      real*8     psig_r(*)
      complex*16 psie(*)
      real*8     psie_r(*)
      real*8     E
      complex*16 dEdpsie(*)

      complex*16 h1psie(*),h1psig(*)
      real*8     rhogg(*),rhoee(*),rhoge(*)
      real*8     vcgg(*),vcee(*),vcge(*)
      real*8     vctmp(*),vc(*)
      complex*16 vcc(*)
      real*8     hpsig_r(*)
      real*8     hpsie_r(*)

#include "bafdecls.fh"
#include "psi.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical oprint
      integer nx,ny,nz,n2ft3d
      real*8 Sgg,Sge,See,c1,c2,c3,A,B,C
      real*8 H1gg,H1ge,H1gm,H1eg,H1ee,H1em,H1mg,H1me,H1mm
      real*8 H2gg,H2ge,H2gm,H2eg,H2ee,H2em,H2mg,H2me,H2mm
      real*8 Hgg,Hge,Hgm,Heg,Hee,Hem,Hmg,Hme,Hmm
      real*8 ehartr,scal1,scal2,dv,Elow,Ehigh,lmbda
      real*8 ehfx,phfx,Eion
      real*8 Hci(3,3),Eci(3),wkopt(1),work(999)
      !real*8 Hci(2,2),Eci(2),wkopt(1),work(999)
      integer lwork,INFO

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega
      integer  control_version
      external control_version
      real*8   coulomb_e,ewald_e,ion_ion_e
      external coulomb_e,ewald_e,ion_ion_e

c     call Parallel_taskid(taskid)
      !write(*,*) "ENTER 3x3 GRADIENT"
c     oprint = .true.
c     oprint = (taskid.eq.MASTER)
      oprint = .false.
c     if (oprint) then
c        write(*,*) "here here taskid=",taskid
c     end if

      !*** generate psie_r ****
      call Pack_c_Copy(1,psie,psie_r)
      call Pack_c_unpack(1,psie_r)
      call D3dB_cr_fft3b(1,psie_r)
      call D3dB_r_Zero_Ends(1,psie_r)

      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv = scal1*lattice_omega()

*     **** ion-ion part ****
      if (control_version().eq.3) Eion = ewald_e()   !**** get ewald energy ****
      if (control_version().eq.4) Eion = ion_ion_e() !**** get free-space ion-ion energy ****
      if (oprint) write(*,*) "Eion=",Eion

      !*** calculate Sgg, Sge, See ***
      call Pack_cc_dot(1,psig,psig,Sgg)
      call Pack_cc_dot(1,psig,psie,Sge)
      call Pack_cc_dot(1,psie,psie,See)
      if (oprint) write(*,*) "Sgg=",Sgg," Sge=",Sge," See=",See


*     **** apply H1 operator ****
      call  psi_H1psi(ispin,neq,npack1,n2ft3d,psig,psig_r,h1psig)
      call  psi_H1psi(ispin,neq,npack1,n2ft3d,psie,psie_r,h1psie)

      call Pack_cc_dot(1,psig,h1psig,H1gg)
      call Pack_cc_dot(1,psig,h1psie,H1ge)
      call Pack_cc_dot(1,psie,h1psie,H1ee)
      H1gg = -H1gg
      H1ge = -H1ge
      H1ee = -H1ee
      H1eg =  H1ge
      H1gm = dsqrt(2.0d0)*H1ge
      H1em = H1gm
      H1mg = H1gm
      H1me = H1em
      H1mm = H1gg + H1ee
      if (ispin.eq.1) then
         H1gg = H1gg + H1gg
         H1ge = H1ge + H1ge
         H1ee = H1ee + H1ee
         H1eg = H1ge
      end if
      H1ge = H1ge*Sge
      H1eg = H1eg*Sge

      if (oprint) write(*,*) "H1gg=",H1gg," H1ge=",H1ge," H1gm=",H1gm
      if (oprint) write(*,*) "H1eg=",H1eg," H1ee=",H1ee," H1em=",H1em
      if (oprint) write(*,*) "H1mg=",H1mg," H1me=",H1me," H1mm=",H1mm

*     **** apply H2 operator ****

      call pspw_et_gen_rho(ispin,neq,n2ft3d,psig_r,rhogg)
      call pspw_et_gen_rho12(ispin,neq,n2ft3d,psig_r,psie_r,rhoge)
      call pspw_et_gen_rho(ispin,neq,n2ft3d,psie_r,rhoee)

      if (control_version().eq.4) then
         call coulomb2_v(rhogg,vcgg)
         call coulomb2_v(rhoge,vcge)
         call coulomb2_v(rhoee,vcee)

         call D3dB_rr_dot(1,rhogg,vcgg,ehartr)
         H2gg = 0.5d0*(0.5d0*ehartr*dv)

         call D3dB_rr_dot(1,rhoge,vcge,ehartr)
         H2ge = 0.5d0*(0.5d0*ehartr*dv)

         call D3dB_rr_dot(1,rhogg,vcge,ehartr)
         H2gm = dsqrt(2.0d0)*0.5d0*(0.5d0*ehartr*dv)

         call D3dB_rr_dot(1,rhoge,vcge,ehartr)
         H2eg = 0.5d0*(0.5d0*ehartr*dv)

         call D3dB_rr_dot(1,rhoee,vcee,ehartr)
         H2ee = 0.5d0*(0.5d0*ehartr*dv)
      
         call D3dB_rr_dot(1,rhoee,vcge,ehartr)
         H2em = dsqrt(2.0d0)*0.5d0*(0.5d0*ehartr*dv)
 
         call D3dB_rr_dot(1,rhoge,vcgg,ehartr)
         H2mg = dsqrt(2.0d0)*0.5d0*(0.5d0*ehartr*dv)

         call D3dB_rr_dot(1,rhoge,vcee,ehartr)
         H2me = dsqrt(2.0d0)*0.5d0*(0.5d0*ehartr*dv)

         call D3dB_rr_dot(1,rhogg,vcee,ehartr)
         H2mm = 0.5d0*(0.5d0*ehartr*dv)
         H2mm = H2mm + H2ge
      else
*        *** generate <g|H2|g> ***
         call D3dB_r_SMul1(1,scal1,rhogg)
         call D3dB_rc_fft3f(1,rhogg)
         call Pack_c_pack(0,rhogg)
         H2gg = coulomb_e(rhogg)
         call coulomb_v(rhogg,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_r_Copy(1,vcc,vcgg)
         call D3dB_r_nZero(1,neq(1)+neq(2),hpsig_r)
         call pspw_potential_HFX(ispin,psig_r,hpsig_r)
         call pspw_energy_HFX(ispin,psig_r,ehfx,phfx)
         H2gg = H2gg + ehfx

*        *** generate <g|H2|e> ***
         call D3dB_r_nZero(1,neq(1)+neq(2),hpsig_r)
         call D3dB_r_nZero(1,neq(1)+neq(2),hpsie_r)
         call pspw_potential_HFX2_dnc(ispin,psig_r,psie_r,
     >                                   hpsig_r,hpsie_r,
     >                                   ehfx,phfx)
         H2ge = -ehfx 

*        *** generate <g|H2|m> ***
         call D3dB_r_SMul1(1,scal1,rhoge)
         call D3dB_rc_fft3f(1,rhoge)
         call Pack_c_pack(0,rhoge)
         call coulomb_v(rhoge,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_r_Copy(1,vcc,vcge)
         call coulomb_screened_v(rhoge,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_rr_daxpy(1,-0.5d0,vcc,vcge)
         call pspw_et_gen_rho(ispin,neq,n2ft3d,psig_r,rhogg)
         call D3dB_rr_dot(1,rhogg,vcge,H2gm)
         H2gm = dsqrt(2.0d0)*(0.50d0*H2gm*dv)

*        *** generate <e|H2|e> ***
         call D3dB_r_SMul1(1,scal1,rhoee)
         call D3dB_rc_fft3f(1,rhoee)
         call Pack_c_pack(0,rhoee)
         H2ee = coulomb_e(rhoee)
         call coulomb_v(rhoee,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_r_Copy(1,vcc,vcee)
         call D3dB_r_nZero(1,neq(1)+neq(2),hpsig_r)
         call pspw_potential_HFX(ispin,psie_r,hpsig_r)
         call pspw_energy_HFX(ispin,psie_r,ehfx,phfx)
         H2ee = H2ee + ehfx

*        *** generate <e|H2|m> ***
         call pspw_et_gen_rho(ispin,neq,n2ft3d,psie_r,rhoee)
         call D3dB_rr_dot(1,rhoee,vcge,H2em)
         H2em = dsqrt(2.0d0)*(0.50d0*H2em*dv)

*        *** generate <m|H2|m> ***
         call D3dB_r_SMul1(1,scal1,rhoee)
         call D3dB_rc_fft3f(1,rhoee)
         call Pack_c_pack(0,rhoee)
         call coulomb_v(rhoee,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_r_Copy(1,vcc,vcee)
         call coulomb_screened_v(rhoee,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_rr_daxpy(1,-0.5d0,vcc,vcee)
         call pspw_et_gen_rho(ispin,neq,n2ft3d,psig_r,rhogg)
         call D3dB_rr_dot(1,rhogg,vcee,H2mm)
         H2mm = 0.50d0*H2mm*dv
         H2mm = H2mm + H2ge
      end if
         !if (oprint) write(*,*) "H2gg=",H2gg," H2ge=",H2ge," H2ee=",H2ee
      if (oprint) write(*,*) "H2gg=",H2gg," H2ge=",H2ge," H2gm=",H2gm
      if (oprint) write(*,*) "H2eg=",H2eg," H2ee=",H2ee," H2em=",H2em
      if (oprint) write(*,*) "H2mg=",H2mg," H2me=",H2me," H2mm=",H2mm


*     **** generate and diagonalize 3x3 CI matrix ****
      Hgg = H1gg + H2gg + Eion
      Hge = H1ge + H2ge 
      Hgm = H1gm + H2gm 
c     Heg = H1eg + H2eg 
      Heg = Hge
      Hee = H1ee + H2ee + Eion
      Hem = H1em + H2em 
c     Hmg = H1mg + H2mg 
      Hmg = Hgm
c     Hme = H1me + H2me 
      Hme = Hem
      Hmm = H1mm + H2mm + Eion

      Hci(1,1) = Hgg
      Hci(1,2) = Hge
      Hci(1,3) = Hgm
      Hci(2,1) = Heg
      Hci(2,2) = Hee
      Hci(2,3) = Hem
      Hci(3,1) = Hmg
      Hci(3,2) = Hme
      Hci(3,3) = Hmm

      if (oprint) then 
         write(*,*)"Hgg,Hge,Hgm",Hci(1,1),Hci(1,2),Hci(1,3)
         write(*,*)"Heg,Hee,Hem",Hci(2,1),Hci(2,2),Hci(2,3)
         write(*,*)"Hmg,Hme,Hmm",Hci(3,1),Hci(3,2),Hci(3,3)
      end if

      call DSYEV('V','U',3,Hci,3,Eci,wkopt,-1,INFO )
      lwork = wkopt(1)
      call DSYEV('V','U',3,Hci,3,Eci,work,lwork,INFO )

      c1 = Hci(1,1)
      c2 = Hci(2,1)
      c3 = Hci(3,1)
      if (oprint) write(*,*) "c1,c2,c3",c1,c2,c3

*     **** generate dElow/dpsie ****
      if (control_version().eq.4) then
         call D3dB_r_Copy(1,vcee,vc)
         call D3dB_rr_Mul2(1,psie_r,vc)
         call D3dB_r_SMul1(1,(-2.0d0*c2*c2),vc)

         call D3dB_r_Copy(1,vcge,vctmp)
         call D3dB_rr_Mul2(1,psig_r,vctmp)
         call D3dB_rr_daxpy(1,(-2.0d0*c1*c2),vctmp,vc)

         call D3dB_r_Copy(1,vcee,vctmp)
         call D3dB_rr_Mul2(1,psig_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*c2*c3),vctmp,vc)
 
         call D3dB_r_Copy(1,vcge,vctmp)
         call D3dB_rr_Mul2(1,psie_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*c2*c3),vctmp,vc)
 
         call D3dB_r_Copy(1,vcgg,vctmp)
         call D3dB_rr_Mul2(1,psig_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*c3*c1),vctmp,vc)
 
         call D3dB_r_Copy(1,vcge,vctmp)
         call D3dB_rr_Mul2(1,psie_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*c3*c2),vctmp,vc)
 
         call D3dB_r_Copy(1,vcgg,vctmp)
         call D3dB_rr_Mul2(1,psie_r,vctmp)
         call D3dB_rr_daxpy(1,(-c3*c3),vctmp,vc)
 
         call D3dB_r_Copy(1,vcge,vctmp)
         call D3dB_rr_Mul2(1,psig_r,vctmp)
         call D3dB_rr_daxpy(1,(-c3*c3),vctmp,vc)
      else
         call coulomb_v(rhoee,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_r_Copy(1,vcc,vcee)
 
         call D3dB_r_Copy(1,vcee,vc)
         call D3dB_rr_Mul2(1,psie_r,vc)
         call D3dB_rr_Sum2(1,hpsig_r,vc) 
         call D3dB_r_SMul1(1,(-4.0d0*c2*c2),vc)
 
         call D3dB_rr_daxpy(1,(4.0d0*c1*c2),hpsie_r,vc)
 
         call pspw_et_gen_rho(ispin,neq,n2ft3d,psig_r,rhogg)
         call D3dB_r_SMul1(1,scal1,rhogg)
         call D3dB_rc_fft3f(1,rhogg)
         call Pack_c_pack(0,rhogg)
         call coulomb_v(rhogg,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_r_Copy(1,vcc,vcgg)
         call coulomb_screened_v(rhogg,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_rr_daxpy(1,-0.5d0,vcc,vcgg)

         call pspw_et_gen_rho(ispin,neq,n2ft3d,psie_r,rhoee)
         call D3dB_r_SMul1(1,scal1,rhoee)
         call D3dB_rc_fft3f(1,rhoee)
         call Pack_c_pack(0,rhoee)
         call coulomb_v(rhoee,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_r_Copy(1,vcc,vcee)
         call coulomb_screened_v(rhoee,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_rr_daxpy(1,-0.5d0,vcc,vcee)

         call pspw_et_gen_rho12(ispin,neq,n2ft3d,psig_r,psie_r,rhoge)
         call D3dB_r_SMul1(1,scal1,rhoge)
         call D3dB_rc_fft3f(1,rhoge)
         call Pack_c_pack(0,rhoge)
         call coulomb_v(rhoge,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_r_Copy(1,vcc,vcge)
         call coulomb_screened_v(rhoge,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_rr_daxpy(1,-0.5d0,vcc,vcge)

         call D3dB_r_Copy(1,vcee,vctmp)
         call D3dB_rr_Mul2(1,psig_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*2.0d0*c2*c3),vctmp,vc)
 
         call D3dB_r_Copy(1,vcge,vctmp)
         call D3dB_rr_Mul2(1,psie_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*2.0d0*c2*c3),vctmp,vc)
c
         call D3dB_r_Copy(1,vcgg,vctmp)
         call D3dB_rr_Mul2(1,psig_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*2.0d0*c3*c1),vctmp,vc)
c
         call D3dB_r_Copy(1,vcge,vctmp)
         call D3dB_rr_Mul2(1,psie_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*2.0d0*c3*c2),vctmp,vc)

         call D3dB_r_Copy(1,vcgg,vctmp)
         call D3dB_rr_Mul2(1,psie_r,vctmp)
         call D3dB_rr_daxpy(1,(-2.0d0*c3*c3),vctmp,vc)
 
         call D3dB_rr_daxpy(1,(2.0d0*c3*c3),hpsie_r,vc)
      end if

      call D3dB_rc_fft3f(1,vc)
      call Pack_c_pack(1,vc)
      call Pack_c_SMul(1,scal1,vc,dEdpsie)

      call Pack_cc_daxpy(1,4.0d0*c2*c2,h1psie,dEdpsie)
      call Pack_cc_daxpy(1,2.0d0*dsqrt(2.0d0)*c2*c3,h1psig,dEdpsie)
      call Pack_cc_daxpy(1,2.0d0*dsqrt(2.0d0)*c3*c1,h1psig,dEdpsie)
      call Pack_cc_daxpy(1,2.0d0*c3*c3,h1psie,dEdpsie)

      E = Eci(1)
      !if (oprint) write(*,*) "c1,c2=",c1,c2

      return
      end


*     ***********************************
*     *                                 *
*     *    psi_3x3ne_virtual_gradient   *
*     *                                 *
*     ***********************************
*
*    This routine calculates the 3x3 CI Energy and its gradient wrt to psie.
*
*    Entry - psig,psig_r
*            psie 
*
*    Exit - E,dEpsie
*
*    Use - psie_r
*          h1psig,h1psie,rhogg,rhoge,rhoee,
*          vctmp,vc,vcc,vcgg,vcge,vcee,
*          hpsig_r,hpsie_r
*
*    The matrix is
*        g    e    m    
*    g  Hgg  Hge  Hgm  
*
*    e  Heg  Hee  Hem  
*
*    m  Hmg  Hme  Hmm  

      subroutine psi_3x3ne_virtual_gradient(psig,psig_r,
     >                                    psie,psie_r,
     >                                    psihomo,psihomo_r,
     >                                    psilumo,psilumo_r,
     >                                    E,dEdpsie,
     >                                    h1psig,h1psie,
     >                                    rhogg,rhoge,rhoee,
     >                                    vctmp,vc,
     >                                    vcc,
     >                                    vcgg,vcge,vcee,
     >                                    hpsig,hpsie,
     >                                    hpsig_r,hpsie_r)
      implicit none
      complex*16 psig(*)
      real*8     psig_r(*)
      complex*16 psie(*)
      real*8     psie_r(*)

      complex*16 psihomo(*)
      real*8     psihomo_r(*)

      complex*16 psilumo(*)
      real*8     psilumo_r(*)

      real*8     E
      complex*16 dEdpsie(*)

      complex*16 h1psie(*),h1psig(*)
      real*8     rhogg(*),rhoee(*),rhoge(*)
      real*8     vcgg(*),vcee(*),vcge(*)
      real*8     vctmp(*),vc(*)
      complex*16 vcc(*)
      complex*16 hpsig(*)
      complex*16 hpsie(*)
      real*8     hpsig_r(*)
      real*8     hpsie_r(*)

#include "bafdecls.fh"
#include "psi.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"

c     real*8 Hgg
c     common / CI_Hgg_common / Hgg

*     **** local variables ****
      logical oprint
      integer neq1(2),i,q,ms,n1q(2),n2q(2)
      integer nx,ny,nz,n2ft3d
      real*8 Sgg,Sge,See,c1,c2,c3,A,B,C
      real*8 H1gg,H1ge,H1gm,H1eg,H1ee,H1em,H1mg,H1me,H1mm
      real*8 H2gg,H2ge,H2gm,H2eg,H2ee,H2em,H2mg,H2me,H2mm
      real*8 Hgg,Hge,Hgm,Heg,Hee,Hem,Hmg,Hme,Hmm
      real*8 ehartr,scal1,scal2,dv,Elow,Ehigh,lmbda,sum1
      real*8 ehfx,phfx,Eion
      real*8 Hci(3,3),Eci(3),wkopt(1),work(999)
      integer lwork,INFO

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega
      integer  control_version
      external control_version
      real*8   coulomb_e,ewald_e,ion_ion_e
      external coulomb_e,ewald_e,ion_ion_e

      !write(*,*) "ENTER 3x3 GRADIENT"
c     oprint = .true.
      oprint = .false.

      !*** generate psilumo_r ****
      call Pack_c_Copy(1,psilumo,psilumo_r)
      call Pack_c_unpack(1,psilumo_r)
      call D3dB_cr_fft3b(1,psilumo_r)
      call D3dB_r_Zero_Ends(1,psilumo_r)

      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv = scal1*lattice_omega()

      call Pack_c_Zero(1,dEdpsie)

*     **** ion-ion part ****
      if (control_version().eq.3) Eion = ewald_e()   !**** get ewald energy ****
      if (control_version().eq.4) Eion = ion_ion_e() !**** get free-space ion-ion energy ****
      if (oprint) write(*,*) "Eion=",Eion

      !*** calculate Sgg, Sge, See ***
      Sgg = 1.0d0
      See = 1.0d0
      Sge = 0.0d0
      if (oprint) write(*,*) "Sgg=",Sgg," Sge=",Sge," See=",See

*     ***************************
*     **** one-electron part ****
*     ***************************
      n1q(1) = 1
      n2q(1) = neq(1)
      n1q(2) = neq(1)+1
      n2q(2) = neq(1)+neq(2)

      !*** generate HF energy and gradient of state psig ****
      !*** generate Hgg if Hgg==0 ***
*     *** generate <g|H1|g> ***
c     if (dabs(Hgg).lt.1.0d-9) then
*        **** apply H1 operator ****
         call psi_H1psi(ispin,neq,npack1,n2ft3d,psig,psig_r,hpsig)
         call Pack_c_Copy(1,hpsig,h1psig)
         H1gg=0.0d0
         do ms=1,ispin
         if (neq(ms).gt.0) then
            do q=n1q(ms),n2q(ms)
               call Pack_cc_idot(1,psig(1+npack1*(q-1)),
     >                            hpsig(1+npack1*(q-1)),sum1)
               H1gg = H1gg - sum1
            end do
         end if
         end do
         call Parallel_SumAll(H1gg)

         !write(*,*) "H1gg=",H1gg
c     end if

      call Pack_c_Copy(1,psihomo,psig)
      call D3dB_c_Copy(1,psihomo_r,psig_r)
      call Pack_c_Copy(1,psilumo,psie)
      call D3dB_c_Copy(1,psilumo_r,psie_r)

*     *** generate <e|H1|e> ***
      call psi_H1psi(ispin,neq,npack1,n2ft3d,psie,psie_r,hpsie)
      call Pack_c_Copy(1,hpsie,h1psie)

      H1ee=0.0d0
      do ms=1,ispin
         if (neq(ms).gt.0) then
            do q=n1q(ms),n2q(ms)
               call Pack_cc_idot(1,psie(1+npack1*(q-1)),
     >                            hpsie(1+npack1*(q-1)),sum1)
               H1ee = H1ee - sum1
            end do
         end if
      end do
      call Parallel_SumAll(H1ee)
      !write(*,*) "H1ee=",H1ee

      H1ge=0.0d0
      do ms=1,ispin
         if (neq(ms).gt.0) then
            do q=n1q(ms),n2q(ms)
               call Pack_cc_idot(1,psig(1+npack1*(q-1)),
     >                            hpsie(1+npack1*(q-1)),sum1)
               H1ge = H1ge - sum1
            end do
         end if
      end do
      call Parallel_SumAll(H1ge)

      H1eg =  H1ge
      H1gm = dsqrt(2.0d0)*H1ge
      H1em = H1gm
      H1mg = H1gm
      H1me = H1em
      H1mm = H1gg + H1ee
      if (ispin.eq.1) then
         H1gg = H1gg + H1gg
         H1ge = H1ge + H1ge
         H1ee = H1ee + H1ee
      end if
      H1ge = H1ge*Sge
      H1eg = H1eg*Sge

      if (oprint) write(*,*) "H1gg=",H1gg," H1ge=",H1ge," H1gm=",H1gm
      if (oprint) write(*,*) "H1eg=",H1eg," H1ee=",H1ee," H1em=",H1em
      if (oprint) write(*,*) "H1mg=",H1mg," H1me=",H1me," H1mm=",H1mm

*     **** apply H2 operator ****
      neq1(1) = 1
      neq1(2) = 0

      call pspw_et_gen_rho(ispin,neq,n2ft3d,psig_r,rhogg)
      call pspw_et_gen_rho12(ispin,neq1,n2ft3d,psihomo_r,psilumo_r,
     >                       rhoge)
      call pspw_et_gen_rho(ispin,neq,n2ft3d,psie_r,rhoee)

      if (control_version().eq.4) then
         call coulomb2_v(rhogg,vcgg)
         call coulomb2_v(rhoge,vcge)
         call coulomb2_v(rhoee,vcee)

c        if (dabs(Hgg).lt.1.0d-9) then
            call D3dB_rr_dot(1,rhogg,vcgg,ehartr)
            H2gg = (0.5d0*ehartr*dv)
            call D3dB_r_nZero(1,neq(1)+neq(2),hpsig_r)
            call pspw_potential_HFX(ispin,psig_r,hpsig_r)
            call pspw_energy_HFX(ispin,psig_r,ehfx,phfx)
            H2gg = H2gg + ehfx
c        end if

         call D3dB_rr_dot(1,rhoge,vcge,ehartr)
         H2ge = 0.5d0*(0.5d0*ehartr*dv)

         call D3dB_rr_dot(1,rhogg,vcge,ehartr)
         H2gm = dsqrt(2.0d0)*0.5d0*(0.5d0*ehartr*dv)

         call D3dB_rr_dot(1,rhoge,vcge,ehartr)
         H2eg = 0.5d0*(0.5d0*ehartr*dv)

         call D3dB_rr_dot(1,rhoee,vcee,ehartr)
         H2ee = (0.5d0*ehartr*dv)
         call D3dB_r_nZero(1,neq(1)+neq(2),hpsie_r)
         call pspw_potential_HFX(ispin,psie_r,hpsie_r)
         call pspw_energy_HFX(ispin,psie_r,ehfx,phfx)
         H2ee = H2ee + ehfx

         call D3dB_rr_dot(1,rhoee,vcge,ehartr)
         H2em = dsqrt(2.0d0)*0.5d0*(0.5d0*ehartr*dv)
 
         call D3dB_rr_dot(1,rhoge,vcgg,ehartr)
         H2mg = dsqrt(2.0d0)*0.5d0*(0.5d0*ehartr*dv)

         call D3dB_rr_dot(1,rhoge,vcee,ehartr)
         H2me = dsqrt(2.0d0)*0.5d0*(0.5d0*ehartr*dv)

         call D3dB_rr_dot(1,rhogg,vcee,ehartr)
         H2mm = 0.5d0*(0.5d0*ehartr*dv)
         H2mm = H2mm + H2ge
      else
*        *** generate <g|H2|g> ***
c        if (dabs(Hgg).lt.1.0d-9) then
            call D3dB_r_SMul1(1,scal1,rhogg)
            call D3dB_rc_fft3f(1,rhogg)
            call Pack_c_pack(0,rhogg)
            H2gg = coulomb_e(rhogg)
            call coulomb_v(rhogg,vcc)
            call Pack_c_unpack(0,vcc)
            call D3dB_cr_fft3b(1,vcc)
            call D3dB_r_Copy(1,vcc,vcgg)
            call D3dB_r_nZero(1,neq(1)+neq(2),hpsig_r)
            call pspw_potential_HFX(ispin,psig_r,hpsig_r)
            call pspw_energy_HFX(ispin,psig_r,ehfx,phfx)
            H2gg = H2gg + ehfx
c        end if
*        *** generate <g|H2|e> ***
         call D3dB_r_nZero(1,neq(1)+neq(2),hpsig_r)
         call D3dB_r_nZero(1,neq(1)+neq(2),hpsie_r)
         call pspw_potential_HFX2_dnc(ispin,psig_r,psie_r,
     >                                   hpsig_r,hpsie_r,
     >                                   ehfx,phfx)
         H2ge = -ehfx 

*        *** generate <g|H2|m> ***
         call D3dB_r_SMul1(1,scal1,rhoge)
         call D3dB_rc_fft3f(1,rhoge)
         call Pack_c_pack(0,rhoge)
         call coulomb_v(rhoge,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_r_Copy(1,vcc,vcge)
         call coulomb_screened_v(rhoge,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_rr_daxpy(1,-0.5d0,vcc,vcge)
         call pspw_et_gen_rho(ispin,neq,n2ft3d,psig_r,rhogg)
         call D3dB_rr_dot(1,rhogg,vcge,H2gm)
         H2gm = dsqrt(2.0d0)*(0.50d0*H2gm*dv)

*        *** generate <e|H2|e> ***
         call D3dB_r_SMul1(1,scal1,rhoee)
         call D3dB_rc_fft3f(1,rhoee)
         call Pack_c_pack(0,rhoee)
         H2ee = coulomb_e(rhoee)
         call coulomb_v(rhoee,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_r_Copy(1,vcc,vcee)
         call D3dB_r_nZero(1,neq(1)+neq(2),hpsig_r)
         call pspw_potential_HFX(ispin,psie_r,hpsig_r)
         call pspw_energy_HFX(ispin,psie_r,ehfx,phfx)
         H2ee = H2ee + ehfx

*        *** generate <e|H2|m> ***
         call pspw_et_gen_rho(ispin,neq,n2ft3d,psie_r,rhoee)
         call D3dB_rr_dot(1,rhoee,vcge,H2em)
         H2em = dsqrt(2.0d0)*(0.50d0*H2em*dv)

*        *** generate <m|H2|m> ***
         call D3dB_r_SMul1(1,scal1,rhoee)
         call D3dB_rc_fft3f(1,rhoee)
         call Pack_c_pack(0,rhoee)
         call coulomb_v(rhoee,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_r_Copy(1,vcc,vcee)
         call coulomb_screened_v(rhoee,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_rr_daxpy(1,-0.5d0,vcc,vcee)
         call pspw_et_gen_rho(ispin,neq,n2ft3d,psig_r,rhogg)
         call D3dB_rr_dot(1,rhogg,vcee,H2mm)
         H2mm = 0.50d0*H2mm*dv
         H2mm = H2mm + H2ge
      end if
         !if (oprint) write(*,*) "H2gg=",H2gg," H2ge=",H2ge," H2ee=",H2ee
      if (oprint) write(*,*) "H2gg=",H2gg," H2ge=",H2ge," H2gm=",H2gm
      if (oprint) write(*,*) "H2eg=",H2eg," H2ee=",H2ee," H2em=",H2em
      if (oprint) write(*,*) "H2mg=",H2mg," H2me=",H2me," H2mm=",H2mm


*     **** generate and diagonalize 3x3 CI matrix ****
      Hgg = H1gg + H2gg + Eion
      Hge = H1ge + H2ge 
      Hgm = H1gm + H2gm 
c     Heg = H1eg + H2eg 
      Heg = Hge
      Hee = H1ee + H2ee + Eion
      Hem = H1em + H2em 
c     Hmg = H1mg + H2mg 
      Hmg = Hgm
c     Hme = H1me + H2me 
      Hme = Hem
      Hmm = H1mm + H2mm + Eion

      Hci(1,1) = Hgg
      Hci(1,2) = Hge
      Hci(1,3) = Hgm
      Hci(2,1) = Heg
      Hci(2,2) = Hee
      Hci(2,3) = Hem
      Hci(3,1) = Hmg
      Hci(3,2) = Hme
      Hci(3,3) = Hmm

      if (oprint) then 
         write(*,*)"Hgg,Hge,Hgm",Hci(1,1),Hci(1,2),Hci(1,3)
         write(*,*)"Heg,Hee,Hem",Hci(2,1),Hci(2,2),Hci(2,3)
         write(*,*)"Hmg,Hme,Hmm",Hci(3,1),Hci(3,2),Hci(3,3)
      end if

      call DSYEV('V','U',3,Hci,3,Eci,wkopt,-1,INFO )
      lwork = wkopt(1)
      call DSYEV('V','U',3,Hci,3,Eci,work,lwork,INFO )

      c1 = Hci(1,1)
      c2 = Hci(2,1)
      c3 = Hci(3,1)
      if (oprint) write(*,*) "c1,c2,c3",c1,c2,c3

*     **** generate dElow/dpsie ****
      if (control_version().eq.4) then
         call D3dB_r_Copy(1,vcee,vc)
         call D3dB_rr_Mul2(1,psilumo_r,vc)
         call D3dB_rr_Sum2(1,hpsie_r,vc)
         call D3dB_r_SMul1(1,(-4.0d0*c2*c2),vc)

         call D3dB_r_Copy(1,vcge,vctmp)
         call D3dB_rr_Mul2(1,psihomo_r,vctmp)
         call D3dB_rr_daxpy(1,(-2.0d0*c1*c2),vctmp,vc)

         call D3dB_r_Copy(1,vcee,vctmp)
         call D3dB_rr_Mul2(1,psihomo_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*c2*c3),vctmp,vc)
 
         call D3dB_r_Copy(1,vcge,vctmp)
         call D3dB_rr_Mul2(1,psilumo_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*c2*c3),vctmp,vc)
 
         call D3dB_r_Copy(1,vcgg,vctmp)
         call D3dB_rr_Mul2(1,psihomo_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*c3*c1),vctmp,vc)
 
         call D3dB_r_Copy(1,vcge,vctmp)
         call D3dB_rr_Mul2(1,psilumo_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*c3*c2),vctmp,vc)
 
         call D3dB_r_Copy(1,vcgg,vctmp)
         call D3dB_rr_Mul2(1,psilumo_r,vctmp)
         call D3dB_rr_daxpy(1,(-c3*c3),vctmp,vc)
 
         call D3dB_r_Copy(1,vcge,vctmp)
         call D3dB_rr_Mul2(1,psihomo_r,vctmp)
         call D3dB_rr_daxpy(1,(-c3*c3),vctmp,vc)
      else
         call coulomb_v(rhoee,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_r_Copy(1,vcc,vcee)

         call D3dB_r_Copy(1,vcee,vc)
         call D3dB_rr_Mul2(1,psilumo_r,vc)
         call D3dB_rr_Sum2(1,hpsig_r,vc) 
         call D3dB_r_SMul1(1,(-4.0d0*c2*c2),vc)

         call D3dB_rr_daxpy(1,(4.0d0*c1*c2),hpsie_r,vc)

         call pspw_et_gen_rho(ispin,neq,n2ft3d,psig_r,rhogg)
         call D3dB_r_SMul1(1,scal1,rhogg)
         call D3dB_rc_fft3f(1,rhogg)
         call Pack_c_pack(0,rhogg)
         call coulomb_v(rhogg,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_r_Copy(1,vcc,vcgg)
         call coulomb_screened_v(rhogg,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_rr_daxpy(1,-0.5d0,vcc,vcgg)

         call pspw_et_gen_rho(ispin,neq,n2ft3d,psie_r,rhoee)
         call D3dB_r_SMul1(1,scal1,rhoee)
         call D3dB_rc_fft3f(1,rhoee)
         call Pack_c_pack(0,rhoee)
         call coulomb_v(rhoee,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_r_Copy(1,vcc,vcee)
         call coulomb_screened_v(rhoee,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_rr_daxpy(1,-0.5d0,vcc,vcee)

         call pspw_et_gen_rho12(ispin,neq1,n2ft3d,psihomo_r,psilumo_r,
     >                          rhoge)
         call D3dB_r_SMul1(1,scal1,rhoge)
         call D3dB_rc_fft3f(1,rhoge)
         call Pack_c_pack(0,rhoge)
         call coulomb_v(rhoge,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_r_Copy(1,vcc,vcge)
         call coulomb_screened_v(rhoge,vcc)
         call Pack_c_unpack(0,vcc)
         call D3dB_cr_fft3b(1,vcc)
         call D3dB_rr_daxpy(1,-0.5d0,vcc,vcge)

         call D3dB_r_Copy(1,vcee,vctmp)
         call D3dB_rr_Mul2(1,psihomo_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*2.0d0*c2*c3),vctmp,vc)

         call D3dB_r_Copy(1,vcge,vctmp)
         call D3dB_rr_Mul2(1,psilumo_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*2.0d0*c2*c3),vctmp,vc)

         call D3dB_r_Copy(1,vcgg,vctmp)
         call D3dB_rr_Mul2(1,psihomo_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*2.0d0*c3*c1),vctmp,vc)

         call D3dB_r_Copy(1,vcge,vctmp)
         call D3dB_rr_Mul2(1,psilumo_r,vctmp)
         call D3dB_rr_daxpy(1,(-dsqrt(2.0d0)*2.0d0*c3*c2),vctmp,vc)

         call D3dB_r_Copy(1,vcgg,vctmp)
         call D3dB_rr_Mul2(1,psilumo_r,vctmp)
         call D3dB_rr_daxpy(1,(-2.0d0*c3*c3),vctmp,vc)

         call D3dB_rr_daxpy(1,(2.0d0*c3*c3),hpsie_r,vc)
      end if

      call D3dB_rc_fft3f(1,vc)
      call Pack_c_pack(1,vc)
      call Pack_c_SMul(1,scal1,vc,dEdpsie)

      call Pack_cc_daxpy(1,4.0d0*c2*c2,h1psie,dEdpsie)
      call Pack_cc_daxpy(1,2.0d0*dsqrt(2.0d0)*c2*c3,h1psig,dEdpsie)
      call Pack_cc_daxpy(1,2.0d0*dsqrt(2.0d0)*c3*c1,h1psig,dEdpsie)
      call Pack_cc_daxpy(1,2.0d0*c3*c3,h1psie,dEdpsie)

      E = Eci(1)
      !if (oprint) write(*,*) "c1,c2=",c1,c2

      return
      end

*     ***********************************
*     *                                 *
*     *      psi_4x4_virtual_gradient   *
*     *                                 *
*     ***********************************

*     This routine calculates the 4x4 CI Energy and its gradient wrt to psie.
*
*       g    e    a    b
*   g  Hgg  Hge  Hga  Hgb
*
*   e  Heg  Hee  Hea  Heb
*
*   a  Hag  Hae  Haa  Hab
*
*   b  Hbg  Hbe  Hba  Hbb

      subroutine psi_4x4_virtual_gradient(psig,psig_r,psie,psie_r,
     >                                    E,dEdpsie,
     >                                    h1psig,h1psie,
     >                                    rhogg,rhoge,rhoee,
     >                                    vctmp,vc,
     >                                    vcgg,vcge,vcee)
      implicit none
      complex*16 psig(*)
      real*8     psig_r(*)
      complex*16 psie(*)
      real*8     psie_r(*)
      real*8     E
      complex*16 dEdpsie(*)

      complex*16 h1psie(*),h1psig(*)
      real*8     rhogg(*),rhoee(*),rhoge(*)
      real*8     vcgg(*),vcee(*),vcge(*)
      real*8     vctmp(*),vc(*)

#include "bafdecls.fh"
#include "psi.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"

*     **** local variables ****
      integer nx,ny,nz,n2ft3d
      real*8 Sgg,Sge,See,c1,c2,c3,c4
      real*8 H1gg,H1ge,H1ga,H1gb,H1ee,H1ea,H1eb,H1aa,H1ab,H1bb
      real*8 H2gg,H2ge,H2ga,H2gb,H2ee,H2ea,H2eb,H2aa,H2ab,H2bb
      real*8 Hgg,Hge,Hga,Hgb,Hee,Hea,Heb,Haa,Hab,Hbb
      real*8 ehartr,scal1,scal2,dv,Elow,Ehigh,lmbda
      real*8 Hci(4,4),ciE(4),wkopt(1),work(999)
      integer lwork,INFO

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      !write(*,*) "ENTER 4x4 GRADIENT"

      !*** generate psie_r ****
      call Pack_c_Copy(1,psie,psie_r)
      call Pack_c_unpack(1,psie_r)
      call D3dB_cr_fft3b(1,psie_r)
      call D3dB_r_Zero_Ends(1,psie_r)


      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv = scal1*lattice_omega()

      !*** calculate Sgg, Sge, See ***
      call Pack_cc_dot(1,psig,psig,Sgg)
      call Pack_cc_dot(1,psig,psie,Sge)
      call Pack_cc_dot(1,psie,psie,See)
      !if (oprint) write(*,*) "Sgg=",Sgg," Sge=",Sge," See=",See

*     **** apply H1 operator ****
      call  psi_H1psi(ispin,neq,npack1,n2ft3d,psig,psig_r,h1psig)
      call  psi_H1psi(ispin,neq,npack1,n2ft3d,psie,psie_r,h1psie)

      call Pack_cc_dot(1,psig,h1psig,H1gg)
      call Pack_cc_dot(1,psig,h1psie,H1ge)
      call Pack_cc_dot(1,psie,h1psie,H1ee)
      H1gg = -H1gg
      H1ge = -H1ge
      H1ee = -H1ee
      H1ga =  H1ge
      H1gb = -H1ge
      H1ea =  H1ge
      H1eb = -H1ge
      H1aa =  H1gg + H1ee
      H1ab =  H1gg + H1ee
      H1bb =  H1gg + H1ee
      if (ispin.eq.1) then
         H1gg =  H1gg + H1gg
         H1ge =  H1ge + H1ge
         H1ee =  H1ee + H1ee
      end if
      H1ge = H1ge*Sge
      H1ab = H1ab*Sge

*     **** apply H2 operator ****
      call pspw_et_gen_rho(ispin,neq,n2ft3d,psig_r,rhogg)
      call coulomb2_v(rhogg,vcgg)
      call pspw_et_gen_rho(ispin,neq,n2ft3d,psie_r,rhoee)
      call coulomb2_v(rhoee,vcee)
      call pspw_et_gen_rho12(ispin,neq,n2ft3d,psig_r,psie_r,rhoge)
      call coulomb2_v(rhoge,vcge)

      call D3dB_rr_dot(1,rhogg,vcgg,ehartr)
      H2gg = 0.5d0*(0.5d0*ehartr*dv)

      call D3dB_rr_dot(1,rhoge,vcge,ehartr)
      H2ge = 0.5d0*(0.5d0*ehartr*dv)

      call D3dB_rr_dot(1,rhoge,vcgg,ehartr)
      H2ga = 0.5d0*(0.5d0*ehartr*dv)
      
      H2gb = -H2ga

      call D3dB_rr_dot(1,rhoee,vcee,ehartr)
      H2ee = 0.5d0*(0.5d0*ehartr*dv)

      call D3dB_rr_dot(1,rhoge,vcee,ehartr)
      H2ea = 0.5d0*(0.5d0*ehartr*dv)

      H2eb = -H2ea

      call D3dB_rr_dot(1,rhogg,vcee,ehartr)
      H2aa = 0.5d0*(0.5d0*ehartr*dv)

      H2ab = -H2ge

      H2bb = H2aa

      !if (oprint) write(*,*) "H2gg=",H2gg," H2ge=",H2ge," H2ee=",H2ee

*     **** generate and diagonalize 4x4 CI matrix ****
      Hgg = H1gg + H2gg
      Hge = H1ge + H2ge
      Hga = H1ga + H2ga
      Hgb = H1gb + H2gb
      Hee = H1ee + H2ee
      Hea = H1ea + H2ea
      Heb = H1eb + H2eb
      Haa = H1aa + H2aa
      Hab = H1ab + H2ab
      Hbb = H1bb + H2bb

      Hci(1,1) = Hgg
      Hci(1,2) = Hge
      Hci(1,3) = Hga
      Hci(1,4) = Hgb
      Hci(2,1) = Hci(1,2)
      Hci(2,2) = Hee
      Hci(2,3) = Hea
      Hci(2,4) = Heb
      Hci(3,1) = Hci(1,3)
      Hci(3,2) = Hci(2,3)
      Hci(3,3) = Haa
      Hci(3,4) = Hab
      Hci(4,1) = Hci(1,4)
      Hci(4,2) = Hci(2,4)
      Hci(4,3) = Hci(3,4)
      Hci(4,4) = Hbb
      
c     write(*,*) "Matrix"
c     write(*,*) Hci(1,1),Hci(1,2),Hci(1,3),Hci(1,4) 
c     write(*,*) Hci(2,1),Hci(2,2),Hci(2,3),Hci(2,4) 
c     write(*,*) Hci(3,1),Hci(3,2),Hci(3,3),Hci(3,4) 
c     write(*,*) Hci(4,1),Hci(4,2),Hci(4,3),Hci(4,4) 
c
      call DSYEV('V','U',4,Hci,4,ciE,wkopt,-1,INFO )
      lwork = wkopt(1)
      call DSYEV('V','U',4,Hci,4,ciE,work,lwork,INFO )

c     write(*,*) "Eigen"
c     write(*,*) Hci(1,1),Hci(1,2),Hci(1,3),Hci(1,4) 
c     write(*,*) Hci(2,1),Hci(2,2),Hci(2,3),Hci(2,4) 
c     write(*,*) Hci(3,1),Hci(3,2),Hci(3,3),Hci(3,4) 
c     write(*,*) Hci(4,1),Hci(4,2),Hci(4,3),Hci(4,4) 
c     write(*,*) ciE(1),ciE(2),ciE(3),ciE(4)

      c1 = Hci(1,1)
      c2 = Hci(2,1)
      c3 = Hci(3,1)
      c4 = Hci(4,1)
      
*     **** generate dElow/dpsie ****
      call D3dB_r_Copy(1,vcee,vc)
      call D3dB_rr_Mul2(1,psie_r,vc)
      call D3dB_r_SMul1(1,(-2.0d0*c2*c2),vc)

      call D3dB_r_Copy(1,vcge,vctmp)
      call D3dB_rr_Mul2(1,psig_r,vctmp)
      call D3dB_rr_daxpy(1,(-2.0d0*c2*c1),vctmp,vc)

      call D3dB_r_Copy(1,vcge,vctmp)
      call D3dB_rr_Mul2(1,psie_r,vctmp)
      call D3dB_rr_daxpy(1,(-c2*c3),vctmp,vc)

      call D3dB_r_Copy(1,vcee,vctmp)
      call D3dB_rr_Mul2(1,psig_r,vctmp)
      call D3dB_rr_daxpy(1,(-c2*c3),vctmp,vc)

      call D3dB_r_Copy(1,vcge,vctmp)
      call D3dB_rr_Mul2(1,psie_r,vctmp)
      call D3dB_rr_daxpy(1,(c2*c4),vctmp,vc)

      call D3dB_r_Copy(1,vcee,vctmp)
      call D3dB_rr_Mul2(1,psig_r,vctmp)
      call D3dB_rr_daxpy(1,(c2*c4),vctmp,vc)

      call D3dB_r_Copy(1,vcgg,vctmp)
      call D3dB_rr_Mul2(1,psig_r,vctmp)
      call D3dB_rr_daxpy(1,(-c3*c1),vctmp,vc)

      call D3dB_r_Copy(1,vcge,vctmp)
      call D3dB_rr_Mul2(1,psie_r,vctmp)
      call D3dB_rr_daxpy(1,(-c3*c2),vctmp,vc)

      call D3dB_r_Copy(1,vcgg,vctmp)
      call D3dB_rr_Mul2(1,psie_r,vctmp)
      call D3dB_rr_daxpy(1,(-c3*c3),vctmp,vc)

      call D3dB_r_Copy(1,vcge,vctmp)
      call D3dB_rr_Mul2(1,psig_r,vctmp)
      call D3dB_rr_daxpy(1,(c3*c4),vctmp,vc)

      call D3dB_r_Copy(1,vcgg,vctmp)
      call D3dB_rr_Mul2(1,psig_r,vctmp)
      call D3dB_rr_daxpy(1,(c4*c1),vctmp,vc)

      call D3dB_r_Copy(1,vcge,vctmp)
      call D3dB_rr_Mul2(1,psie_r,vctmp)
      call D3dB_rr_daxpy(1,(c4*c2),vctmp,vc)

      call D3dB_r_Copy(1,vcge,vctmp)
      call D3dB_rr_Mul2(1,psig_r,vctmp)
      call D3dB_rr_daxpy(1,(c4*c3),vctmp,vc)

      call D3dB_r_Copy(1,vcgg,vctmp)
      call D3dB_rr_Mul2(1,psie_r,vctmp)
      call D3dB_rr_daxpy(1,(-c4*c4),vctmp,vc)

      call D3dB_rc_fft3f(1,vc)
      call Pack_c_pack(1,vc)
      call Pack_c_SMul(1,scal1,vc,dEdpsie)

      call Pack_cc_daxpy(1, 4.0d0*c2*c2,h1psie,dEdpsie)
      call Pack_cc_daxpy(1, 2.0d0*c2*c3,h1psig,dEdpsie)
      call Pack_cc_daxpy(1,-2.0d0*c2*c4,h1psig,dEdpsie)
      call Pack_cc_daxpy(1, 2.0d0*c3*c1,h1psig,dEdpsie)
      call Pack_cc_daxpy(1, 2.0d0*c3*c3,h1psie,dEdpsie)
      call Pack_cc_daxpy(1,-2.0d0*c4*c1,h1psig,dEdpsie)
      call Pack_cc_daxpy(1, 2.0d0*c4*c4,h1psie,dEdpsie)

      E = ciE(1)

      return
      end





*     ****************************************************
*     *                                                  *
*     *                psi_H1psi                         *
*     *                                                  *
*     ****************************************************
      subroutine psi_H1psi(ispin,neq,npack1,n2ft3d,
     >                        psi,psi_r,
     >                        H1psi)
      implicit none
      integer ispin,neq(2),npack1,n2ft3d
      complex*16 psi(npack1,*)
      real*8     psi_r(n2ft3d,*)
      real*8     H1psi(npack1,*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value,aperiodic,field_exist
      integer q,n,m,ms,nfft3d,nemaxq,npack0,n1q(2),n2q(2)
      integer Hpsi2(2),vlr_l(2),r_grid(2),v_field(2),vl(2)
      real*8 tmp1(2),tmp2(2),sum1,sum2

*     **** external functions ****
      logical  pspw_charge_found,pspw_Efield_found
      external pspw_charge_found,pspw_Efield_found
      integer  control_version
      external control_version

      aperiodic   = (control_version().eq.4)
      field_exist = pspw_charge_found().or.pspw_Efield_found()
      nemaxq = neq(1)+neq(2)

      call Pack_npack(0,npack0)
      call D3dB_nfft3d(1,nfft3d)

*     **** allocate memory ****
      value = .true.
      if (aperiodic) then
       value = value.and.
     >        BA_push_get(mt_dbl,(n2ft3d),'vlr_l',vlr_l(2),vlr_l(1))
      end if
      if (field_exist.or.aperiodic) then
         value = value.and.
     >            BA_push_get(mt_dbl,(3*n2ft3d),'r_grid',
     >                        r_grid(2),r_grid(1))
         value = value.and.
     >           BA_push_get(mt_dbl,(n2ft3d),'v_field',
     >                       v_field(2),v_field(1))
      end if
      value = value.and.
     >        BA_push_get(mt_dcpl,(npack0),'vloc',vl(2),vl(1))
      if (.not. value)
     >   call errquit('pspw_et_sub1:out of stack',0,MA_ERR)


*     **** generate r_grid ****
      if (aperiodic.or.field_exist)
     >   call lattice_r_grid(dbl_mb(r_grid(1)))

*     **** generate local pseudopotential  ****
      call v_local(dcpl_mb(vl(1)),.false.,tmp1,tmp2)

*     *** long-range psp for charge systems ***
      if (control_version().eq.4) then
         call v_lr_local(dbl_mb(r_grid(1)),dbl_mb(vlr_l(1)))
      end if

*     ***** generate other real-space fields ****
      if (field_exist) then
         !call dcopy(n2ft3d,0.0d0,0,dbl_mb(v_field(1)),1)
         call Parallel_shared_vector_zero(.true.,n2ft3d,
     >                                    dbl_mb(v_field(1)))
         call pspw_charge_Generate_V(n2ft3d,
     >                               dbl_mb(r_grid(1)),
     >                               dbl_mb(v_field(1)))
         call pspw_Efield_Generate_V(n2ft3d,
     >                               dbl_mb(r_grid(1)),
     >                               dbl_mb(v_field(1)))
      end if


*     **** get H1psi  ****
      !call dcopy(2*nemaxq*npack1,0.0d0,0,H1psi,1)
      call Parallel_shared_vector_zero(.true.,2*nemaxq*npack1,H1psi)
      if (aperiodic) then
         call psi_H1v4(ispin,neq,psi,psi_r,
     >             dcpl_mb(vl(1)),dbl_mb(vlr_l(1)),
     >             dbl_mb(v_field(1)),field_exist,
     >             H1psi)
      else
         call psi_H1(ispin,neq,psi,psi_r,
     >             dcpl_mb(vl(1)),
     >             dbl_mb(v_field(1)),field_exist,
     >             H1psi)
      end if


*     **** deallocate memory ****
      value = BA_pop_stack(vl(2))
      if (field_exist.or.aperiodic) then
         value = value.and.BA_pop_stack(v_field(2))
         value = value.and.BA_pop_stack(r_grid(2))
      end if
      if (aperiodic) then
         value = value.and.BA_pop_stack(vlr_l(2))
      end if
      if (.not.value)
     >   call errquit('psi_H1psi:pop stack',0,MA_ERR)

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

#include "bafdecls.fh"
#include "psi.fh"

      !*** local variables ***
      integer maxit_orb
      integer ii,l,l2
      real*8  sum,maxerror,error_out,e0

      !*** external functions ***
      real*8   control_tole
      external control_tole

      !call psi_gen_density_potentials(1)
      maxit_orb=120
      maxerror = control_tole()

      do ii=1,(ne_excited(1)+ne_excited(2))
         l2= 0

         !*** orthogonalize to lower orbitals  ****
 2       l2 = l2 + 1
         call psi_project_out_virtual1(
     >           ii,
     >           dcpl_mb(psi1_excited(1)+(ii-1)*npack1))

         !*** normalize ****
         call Pack_cc_dot(1,
     >            dcpl_mb(psi1_excited(1) +(ii-1)*npack1),
     >            dcpl_mb(psi1_excited(1) +(ii-1)*npack1),
     >            sum)
         sum = 1.0d0/dsqrt(sum)
         call Pack_c_SMul1(1,sum,
     >            dcpl_mb(psi1_excited(1) +(ii-1)*npack1))

         !*** minimize orbital ****
          l = 0
 3        call psi_KS_update_virtual(maxit_orb,
     >                               maxerror,
     >                               0.001d0,ii,error_out,e0)
          l  = l+1
          if ((error_out.gt.maxerror).and.(l.le.(1+(l2-1)*3))) go to 3
          if (((error_out.gt.maxerror).or.(e0.gt.4.0d0))
     >        .and.(l2.le.1)) then
            call Pack_c_Zero(1,
     >               dcpl_mb(psi1_excited(1) +(ii-1)*npack1))
            call Pack_c_setzero(1,1.0d0,
     >               dcpl_mb(psi1_excited(1) +(ii-1)*npack1))
             go to 2
          end if

          dbl_mb(eig_excited(1)+ii-1) = e0

      end do
      call psi_sort_virtual()
      
     
      return
      end

      subroutine psi_sort_virtual()
      implicit none

#include "bafdecls.fh"
#include "psi.fh"
#include "errquit.fh"

      logical value
      integer i,j,ii,jj,ms
      integer r1(2)
      real*8  ei,ej

      value = BA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
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

      value = BA_pop_stack(r1(2))
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
      
#include "bafdecls.fh"
#include "psi.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,done,oneloop,precondition,oprint
      integer it,pit
      real*8 eold,percent_error,error0,de0,lmbda_r0,lmbda_r1
      real*8 theta,ep,sp
      integer r1(2),t0(2),t(2),g(2)
      integer psi_ptr

      logical  control_print
      external control_print
      real*8   control_Ep,control_Sp
      external control_Ep,control_Sp

      psi_ptr=psi1_excited(1)

      call Parallel_taskid(taskid)
      oprint= ((taskid.eq.MASTER).and.control_print(print_medium))


      value = BA_push_get(mt_dcpl,npack1,'t0',t0(2),t0(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'g',g(2),g(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'t',t(2),t(1))
      if (.not. value) call errquit(
     >     'psi_KS_update_virtual: out of stack memory',0, MA_ERR)

      ep = control_Ep()
      sp = control_Sp()
      precondition = .true.
      done = .false.
      error0 = 0.0d0 
      e0 = 0.0d0
      theta = -3.14159d0/600.0d0
      lmbda_r0 = 1.0d0
      it = 0
      pit = 0
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

         precondition = (dabs(e0-eold).gt.(sp*maxerror))

         done = ((it.gt.maxiteration) 
     >           .or.
     >           (dabs(e0-eold).lt.maxerror))

         if (done) go to 4


         call Pack_c_Copy(1,dcpl_mb(g(1)),dcpl_mb(r1(1)))
         call Pack_cc_daxpy(1,(e0),
     >                 dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                 dcpl_mb(r1(1)))

*        **** preconditioning ****
         if (precondition) then
            pit = pit + 1
            call ke_Precondition(npack1,1,
     >                           dcpl_mb(psi_ptr+(ii-1)*npack1),
     >                           dcpl_mb(r1(1)))
         end if 

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
 4    value =           BA_pop_stack(t(2)) 
      value = value.and.BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(r1(2))
      value = value.and.BA_pop_stack(t0(2))
      if (.not. value) call errquit(
     >     'psi_KS_update_virtual: popping stack memory',1, MA_ERR)

      if (oprint) then
         write(luout,921) ii,-e0,dabs(e0-eold),it,pit,ep,sp
  921 format(5x,"orbital",I4," current e=",E10.3,
     >       " (error=",E9.3,")",
     >       " iterations",I4,"(",I4,
     >       " preconditioned, Ep,Sp=",F5.1,F7.1,")")
      end if

      error_out = dabs(e0-eold)
      e0        = -e0
      return
      end

*     ***********************************
*     *				        *
*     *	     psi_linesearch_virtual     *
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
#include "bafdecls.fh"
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
      value = BA_push_get(mt_dcpl,npack1,'orb',
     >                       orb(2),orb(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'g',
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
      value =           BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(orb(2))      
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

      if (.not.BA_push_get(mt_dcpl,npack1,'etmp',etmp(2),etmp(1)))
     > call errquit('psi_project_out_virtual1: out of stack',0,MA_ERR)

         !call dcopy(2*npack1,0.0d0,0,dcpl_mb(etmp(1)),1)
         call Parallel_shared_vector_zero(.true.,
     >                                    2*npack1,dcpl_mb(etmp(1)))
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
         call daxpy_omp(2*npack1,1.0d0,dcpl_mb(etmp(1)),1,Horb,1)

      if (.not.BA_pop_stack(etmp(2)))
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

#include "bafdecls.fh"
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
      if (.not.BA_push_get(mt_dcpl,npack1,'etmp',etmp(2),etmp(1)))
     > call errquit('psi_project_out_virtual1: out of stack',0,MA_ERR)

         !call dcopy(2*npack1,0.0d0,0,dcpl_mb(etmp(1)),1)
         call Parallel_shared_vector_zero(.true.,
     >                       2*npack1,dcpl_mb(etmp(1)))
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

      if (.not.BA_pop_stack(etmp(2)))
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

*
      subroutine psi_KS_update00(psi_number,
     >                         maxit_orbs,maxit_orb,
     >                         ks_algorithm,
     >                         precondition,
     >                         maxerror)
      implicit none
      integer psi_number
      integer maxit_orbs,maxit_orb
      integer ks_algorithm
      logical precondition
      real*8  maxerror

#include "bafdecls.fh"
#include "psi.fh"

*     **** local variables ****
      logical done
      integer i,j,neall
      real*8 error,error_out,tim1,tim2,tim,sum

*     **** external functions ****

      tim = 0.0d0
      neall = neq(1)+neq(2)
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
        
#include "bafdecls.fh"
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
*     *                                 *
*     *      psi_KS_block_update        *
*     *                                 *
*     ***********************************
 
*    This routine (approximately) diagonalizes the KS matrix.
*
*    input variables
*    - current_iteration
*    - maxerror
*    global shared variables output 
*    - Enew
*    - deltae
*
      subroutine psi_KS_block_update(Enew,deltae,current_iteration,
     >                               maxit_orbs,maxerror)
      implicit none
      real*8  Enew,deltae,Ermdr
      integer current_iteration
      integer maxit_orbs
      real*8  maxerror
        
#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      real*8  deltat_min
      parameter (deltat_min=1.0d-2)

      integer G0(2),S0(2),G1(2),G2(2)
      real*8  E0,dE0,Enew0

      logical precondition
      common / cgsd_block2 / precondition

      real*8 tmin,dte,sum0,sum1
      common / bfgs_block / tmin,dte,sum0,sum1

      integer it,it_in
      real*8 deltat,tmin0,deltae0,deltac0
      real*8 max_sigma,dt,kappa

      logical value
      integer neall

*     **** define pointers ****
      integer Y,U,HY,HU

*     **** external functions ****
      logical  control_precondition,psp_pawexist
      external control_precondition,psp_pawexist
      integer  control_lmbfgs_size
      external control_lmbfgs_size
      !integer  control_ks_maxit_orbs,psi_neq
      !external control_ks_maxit_orbs,psi_neq
      integer  psi_neq
      external psi_neq
      real*8   control_time_step
      external control_time_step
      real*8   psi_geodesic_energy0,psi_geodesic_denergy
      external psi_geodesic_energy0,psi_geodesic_denergy
      real*8   linesearch
      external linesearch
      integer  geodesic_get_U_ptr
      external geodesic_get_U_ptr
      logical  Dneall_m_pop_stack,Dneall_m_push_get_block
      external Dneall_m_pop_stack,Dneall_m_push_get_block
      real*8   geodesic_ABC_energy,geodesic_ABC_denergy
      external geodesic_ABC_energy,geodesic_ABC_denergy



      !maxit_orb  = control_ks_maxit_orb()   !*** should be read from rtdb ***
      !maxit_orbs = control_ks_maxit_orbs()  !*** should be read from rtdb ***
      it_in = maxit_orbs

      dt = control_time_step()

      call Parallel_taskid(taskid) 

!$OMP MASTER
      if (current_iteration.eq.1) then
         precondition = control_precondition()
         tmin  = 10*deltat_min
      end if
!$OMP END MASTER

      neall = neq(1)+neq(2)

      call geodesic_ABC_init()

*     **** allocate G0, S0 ****
      value = BA_push_get(mt_dcpl,npack1*neall,
     >                     'S0',S0(2),S0(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1*neall,
     >                     'G0',G0(2),G0(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1*neall,
     >                     'G1',G1(2),G1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1*neall,
     >                     'G2',G2(2),G2(1))
      if (.not.value) 
     >   call errquit('psi_KS_block_update:out of heap memory',0,MA_ERR)

*     **** set ptrs ****
      Y = psi1(1)
      U = geodesic_get_U_ptr()
      HY = G0(1)
      HU = S0(1)

      call Parallel_shared_vector_zero(.true.,2*npack1*neall,
     >                                 dcpl_mb(G0(1)))

*     ***** get the initial direction ****
      if (pawexist) then
         call psi_1get_STgradient(dcpl_mb(S0(1)),dcpl_mb(G0(1)),Enew0)
      else
         if (precondition) then
            call psi_1get_TMgradient0(dcpl_mb(G0(1)),Enew0)
         else
            call psi_1get_Tgradient0(dcpl_mb(G0(1)),Enew0)
         end if
      end if
      call psi_1get_remainder0(Ermdr)
      E0    = Enew0

*     ***** use the initial gradient for the direction ****
      call pspw_lmbfgs_init(control_lmbfgs_size(),dcpl_mb(G0(1)))
      call Grsm_gg_Copy(npack1,neall,
     >                  dcpl_mb(G0(1)),
     >                  dcpl_mb(S0(1)))

      do it=2,it_in
*        **** initialize the geoedesic line data structure ****
         call geodesic_start(dcpl_mb(S0(1)),max_sigma,dE0)

*        **** Copy Hpsi_k to HY then ****
*        **** generate Hpsi using psi_k=U and copy Hpsi_k to HU ****
*        **** compute A,B,C, <S0|H|S0>, <U|H|S0> and <U|H|U> matrices ****
         call geodesic_ABC_start(ispin,ne,npack1,Ermdr,
     >                        dcpl_mb(Y), dcpl_mb(U),
     >                        dcpl_mb(HY),dcpl_mb(HU))

*        ******* line search *********
         if ((tmin.gt.deltat_min).and.(tmin.lt.1.0d4)) then
            deltat = tmin
         else
            deltat = deltat_min
         end if
 20      continue
c         write(*,*) "     --- tmin0=",tmin
         tmin0 = tmin
         deltae0 = deltae
         Enew0 = linesearch(0.0d0,E0,dE0,deltat,
     >                        psi_geodesic_energy0,
     >                        psi_geodesic_denergy,
     >                        0.50d0,tmin0,deltae0,2)

!$OMP MASTER
         tmin = tmin0
         deltae = deltae0
         Enew = Enew0 + Ermdr
!$OMP END MASTER
!$OMP BARRIER
c         write(*,*) "     --- tmin=",tmin
c         write(*,*) "     --- it,Enew,deltae=",it,Enew,deltae

         call psi_geodesic_final(tmin)

*        **** exit loop early ****
         if (dabs(deltae).lt.maxerror) then
            if (.not.precondition) go to 30
            precondition = .false.
         end if

*        **** get the new gradient ****
         if (pawexist) then
              call psi_2get_STgradient(2,dcpl_mb(S0(1)),
     >                                   dcpl_mb(G0(1)),Enew0)
         else
            if (precondition) then
               call psi_2get_TMgradient0(2,dcpl_mb(G0(1)),Enew0)
            else
               call psi_2get_Tgradient0(2,dcpl_mb(G0(1)),Enew0)
            end if
         end if
!$OMP MASTER
      E0 = Enew0
!$OMP END MASTER
!$OMP BARRIER

         call pspw_lmbfgs(tmin,dcpl_mb(G0(1)),dcpl_mb(S0(1)))
         call psi_2to1()

      end do

*     **** initialize the geoedesic line data structure ****
      call geodesic_start(dcpl_mb(S0(1)),max_sigma,dE0)


*     **** Copy Hpsi_k to HY then ****
*     **** generate Hpsi using psi_k=U and copy Hpsi_k to HU ****
*     **** compute A,B,C, <S0|H|S0>, <U|H|S0> and <U|H|U> matrices ****
      call geodesic_ABC_start(ispin,ne,npack1,Ermdr,
     >                        dcpl_mb(Y), dcpl_mb(U),
     >                        dcpl_mb(HY),dcpl_mb(HU))

*     ******* line search *********
      if ((tmin.gt.deltat_min).and.(tmin.lt.1.0d4)) then
         deltat = tmin
      else
         deltat = deltat_min
      end if

 25   continue
      tmin0 = tmin
      deltae0 = deltae
c      Enew0 = linesearch(0.0d0,E0,dE0,deltat,
c     >                        psi_geodesic_energy0,
c     >                        psi_geodesic_denergy,
c     >                        0.50d0,tmin0,deltae0,2)
c      Enew0 = linesearch(0.0d0,E0,dE0,deltat,
c     >                        geodesic_ABC_energy,
c     >                        geodesic_ABC_denergy,
c     >                        0.50d0,tmin0,deltae0,2)
      Enew0 = linesearch(0.0d0,E0+Ermdr,dE0,deltat,
     >                        geodesic_ABC_energy,
     >                        geodesic_ABC_denergy,
     >                        maxerror,tmin0,deltae0,1)

c       if (taskid.eq.MASTER) then
c          write(*,*) "ABC energy E0+Ermdr,Enew0,deltae0 =", 
c     >            E0+Ermdr,Enew0,deltae0,
c     >            tmin0,maxerror,it
c       end if


!$OMP MASTER
         tmin = tmin0
         deltae = deltae0
         Enew = Enew0 + Ermdr
!$OMP END MASTER
!$OMP BARRIER

c      write(*,*) "psi_geodesic_energy0=",tmin,
c     >            psi_geodesic_energy0(tmin),Enew0,
c     >            geodesic_ABC_energy(tmin),
c     >            geodesic_ABC_energy(tmin)-psi_geodesic_energy0(tmin)
c
c      write(*,*) "psi_geodesic_denergy=",tmin,
c     >            psi_geodesic_denergy(tmin),
c     >            geodesic_ABC_denergy(tmin),
c     >            geodesic_ABC_denergy(tmin)-psi_geodesic_denergy(tmin)
 
      call psi_geodesic_final(tmin)

 30   CONTINUE


*     **** free memory ****
      value =           BA_pop_stack(G2(2))
      value = value.and.BA_pop_stack(G1(2))
      value = value.and.BA_pop_stack(G0(2))
      value = value.and.BA_pop_stack(S0(2))
      if (.not. value)
     >  call errquit('psi_KS_block_update:freeing heap',0, MA_ERR)

      call geodesic_ABC_finalize()

      call psi_2to1()

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
      integer psi_number
      logical precondition
      integer maxiteration
      real*8  maxerror,perror 
      integer i
      real*8 error_out
      
#include "bafdecls.fh"
#include "psi.fh"
#include "errquit.fh"

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

      value = BA_push_get(mt_dcpl,npack1,'t0',t0(2),t0(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'g',g(2),g(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'t',t(2),t(1))
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
         !dbl_mb(eig(1)+i-1) = e0
         e0 = -e0

         done = ((it.gt.maxiteration) 
     >           .or.
     >           (dabs(e0-eold).lt.maxerror))

         if (done) go to 4

*        **** preconditioning ****
         if (precondition) then
           call ke_Precondition(npack1,1,
     >                     dcpl_mb(psi_ptr+(i-1)*npack1),
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



*        **** project out psi components from t - may not be needed! ****
         !call psi_project_out_orb(psi_number,i,dcpl_mb(t(1)))
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
 4    value =           BA_pop_stack(t(2)) 
      value = value.and.BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(r1(2))
      value = value.and.BA_pop_stack(t0(2))
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
      
#include "bafdecls.fh"
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
      value =            BA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      value = value.and. BA_push_get(mt_dcpl,npack1,'g1',g1(2),g1(1))
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
      value =           BA_pop_stack(g1(2))
      value = value.and.BA_pop_stack(r1(2))
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
      
#include "bafdecls.fh"
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
      value = BA_push_get(mt_dcpl,npack1,'orb',
     >                       orb(2),orb(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'g',
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
      value =           BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(orb(2))      
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
      
#include "bafdecls.fh"
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
      value = BA_push_get(mt_dcpl,npack1,'orb',
     >                       orb(2),orb(1))
      value = value.and. 
     >        BA_push_get(mt_dcpl,npack1,'g',
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
      value =           BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(orb(2))      
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "psi.fh"

*     **** local variables ****
      logical ok
      integer ii,n,psi_ptr,np
      integer x(2)
      real*8  sum

      call Parallel_np(np)

*     **** allocate stack memory ****
      ok = BA_push_get(mt_dbl,ne(1),'x',x(2),x(1))
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
      ok = BA_pop_stack(x(2))
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


#include "bafdecls.fh"
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

c      call dcopy(4*nfft3d,
c     >           rho, 1,
c     >           dbl_mb(rho_ptr),1)

      call Parallel_shared_vector_copy(.true.,4*nfft3d,
     >                                 rho,dbl_mb(rho_ptr))

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


#include "bafdecls.fh"
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

c      call dcopy(4*nfft3d,
c     >           dbl_mb(rho_ptr),1,
c     >           rho,1)
      call Parallel_shared_vector_copy(.true.,4*nfft3d,
     >                                 dbl_mb(rho_ptr),rho)
      return
      end

*     ***************************
*     *                         *
*     *      psi_write_density  *
*     *                         *
*     ***************************

*    This routine writes the densities in psi to disk.  
* This routine is needed for a band by band minimizer.
*
      subroutine psi_write_density(psi_number)
      implicit none
      integer psi_number

#include "bafdecls.fh"
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

      call rho_write(ispin,dbl_mb(rho_ptr))
      return
      end

*     ***************************
*     *                         *
*     *   psi_try_read_density  *
*     *                         *
*     ***************************

*    This routine reads the densities from disk to psi.
* This routine is needed for a band by band minimizer.
*
      logical function psi_try_read_density(psi_number)
      implicit none
      integer psi_number

#include "bafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** local variables ****
      logical ok
      integer rho_ptr

*     **** external functions ****
      logical  rho_check_header
      external rho_check_header

      if (psi_number.eq.1) then
        rho_ptr = rho1(1)
      else
        rho_ptr = rho2(1)
      end if

      if (rho_check_header(ispin,.false.)) then
         call rho_read(ispin,dbl_mb(rho_ptr))
         ok = .true.
      else
         ok = .false.
      end if

      psi_try_read_density = ok
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
      subroutine psi_gen_density_potentials(psi_number0)
      implicit none
      integer psi_number0


#include "bafdecls.fh"
#include "psi.fh"
   
*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** local variables ****
      logical zeroscf
      integer psi_ptr,rho_ptr,dng_ptr,rho_all_ptr,occ_ptr
      integer psi_number

*     *** hacky for now ***
      psi_number = psi_number0
      zeroscf = .false.
      if (psi_number.gt.2) then
         psi_number = psi_number - 2
         zeroscf = .true.
      end if
      !write(*,*) "psi_number,zeroscf=",psi_number,zeroscf

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


      if (zeroscf) then
          call electron_gen_vall0()
      else
         call electron_gen_densities(dcpl_mb(psi_ptr),
     >                            dbl_mb(rho_ptr),
     >                            dcpl_mb(dng_ptr),
     >                            dbl_mb(rho_all_ptr),
     >                            occupation_on,dbl_mb(occ_ptr))
         call electron_gen_scf_potentials(dbl_mb(rho_ptr),
     >                            dcpl_mb(dng_ptr),
     >                            dbl_mb(rho_all_ptr))
         !if (zeroscf) then
         !   call electron_zero_scf_potentials()
         !end if
         call electron_gen_vall()
      end if
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

#include "bafdecls.fh"
#include "psi.fh"
   
c      call zcopy(npack1*(neq(1)+neq(2)),
c     >           dcpl_mb(psi1(1)),1,
c     >           dcpl_mb(psi2(1)),1)
      call Parallel_shared_vector_copy(.true.,2*npack1*(neq(1)+neq(2)),
     >                                 dcpl_mb(psi1(1)),
     >                                 dcpl_mb(psi2(1)))

      return
      end


*     ***************************
*     *				*
*     *		psi_2to1	*
*     *				*
*     ***************************
      subroutine psi_2to1()
      implicit none

#include "bafdecls.fh"
#include "psi.fh"

 
c      call zcopy(npack1*(neq(1)+neq(2)),
c     >           dcpl_mb(psi2(1)),1,
c     >           dcpl_mb(psi1(1)),1)
      call Parallel_shared_vector_copy(.true.,2*npack1*(neq(1)+neq(2)),
     >                                  dcpl_mb(psi2(1)),
     >                                  dcpl_mb(psi1(1)))

c      call OrthoCheck(ispin,ne,dcpl_mb(psi1(1)))  
      return
      end

*     ***************************
*     *                         *
*     *       psiocc_1to2       *
*     *                         *
*     ***************************
      subroutine psiocc_1to2()
      implicit none

#include "bafdecls.fh"
#include "psi.fh"
  
      call Parallel_shared_vector_copy(.false.,2*npack1*(neq(1)+neq(2)),
     >                                 dcpl_mb(psi1(1)),
     >                                 dcpl_mb(psi2(1)))
      call Parallel_shared_vector_copy(.true.,(ne(1)+ne(2)),
     >                                 dbl_mb(occ1(1)),
     >                                 dbl_mb(occ2(1)))
      return
      end

*     ***************************
*     *                         *
*     *       psiocc_2to1       *
*     *                         *
*     ***************************
      subroutine psiocc_2to1()
      implicit none

#include "bafdecls.fh"
#include "psi.fh"
  
      call Parallel_shared_vector_copy(.false.,2*npack1*(neq(1)+neq(2)),
     >                                 dcpl_mb(psi2(1)),
     >                                 dcpl_mb(psi1(1)))
      call Parallel_shared_vector_copy(.true.,(ne(1)+ne(2)),
     >                                 dbl_mb(occ2(1)),
     >                                 dbl_mb(occ1(1)))
      return
      end


*     ***************************
*     *                         *
*     *         epsi_2to1        *
*     *                         *
*     ***************************
      subroutine epsi_2to1()
      implicit none

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

c      call dcopy(4*nfft3d,
c     >           dbl_mb(rho2(1)),1,
c     >           dbl_mb(rho1(1)),1)
      call Parallel_shared_vector_copy(.true.,4*nfft3d,
     >                    dbl_mb(rho2(1)),
     >                    dbl_mb(rho1(1)))

c      call dcopy(4*nfft3d,
c     >           dbl_mb(rho2_all(1)),1,
c     >           dbl_mb(rho1_all(1)),1)
      call Parallel_shared_vector_copy(.true.,4*nfft3d,
     >           dbl_mb(rho2_all(1)),
     >           dbl_mb(rho1_all(1)))

      return
      end

*     ***************************
*     *				*
*     *		rho_1to2	*
*     *				*
*     ***************************
      subroutine rho_1to2()
      implicit none

#include "bafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

c      call dcopy(4*nfft3d,
c     >           dbl_mb(rho1(1)),1,
c     >           dbl_mb(rho2(1)),1)
c
c      call dcopy(4*nfft3d,
c     >           dbl_mb(rho1_all(1)),1,
c     >           dbl_mb(rho2_all(1)),1)

      call Parallel_shared_vector_copy(.false.,4*nfft3d,
     >           dbl_mb(rho1(1)),
     >           dbl_mb(rho2(1)))
      call Parallel_shared_vector_copy(.true.,4*nfft3d,
     >           dbl_mb(rho1_all(1)),
     >           dbl_mb(rho2_all(1)))

      return
      end

*     ***************************
*     *				*
*     *		dng_2to1	*
*     *				*
*     ***************************
      subroutine dng_2to1()
      implicit none
 
#include "bafdecls.fh"
#include "psi.fh"
 
      call zcopy(npack0,
     >           dcpl_mb(dng2(1)),1,
     >           dcpl_mb(dng1(1)),1)

      return
      end

*     ***************************
*     *                         *
*     *         dng_1to2        *
*     *                         *
*     ***************************
      subroutine dng_1to2()
      implicit none

#include "bafdecls.fh"
#include "psi.fh"

      call zcopy(npack0,
     >           dcpl_mb(dng1(1)),1,
     >           dcpl_mb(dng2(1)),1)

      return
      end



*     ***********************************
*     *					*
*     *		psi_1add_oep_to_vall	*
*     *					*
*     ***********************************
      subroutine psi_1add_oep_to_vall()
      implicit none

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "psi.fh"

*     **** local variables ***
      real*8 ave

*     **** external functions ****
      real*8   psp_kinetic_atom
      external psp_kinetic_atom

      ave = psp_kinetic_atom(ispin,neq,dcpl_mb(psi1(1)))

      psi_1ke_atom = ave
      return
      end

*     ***********************************
*     *                                 *
*     *     psi_1valence_core_atom      *
*     *                                 *
*     ***********************************
      real*8 function psi_1valence_core_atom()
      implicit none

#include "bafdecls.fh"
#include "psi.fh"

*     **** local variables ***
      real*8 ave

*     **** external functions ****
      real*8   psp_valence_core_atom
      external psp_valence_core_atom

      ave =  psp_valence_core_atom(ispin,neq,dcpl_mb(psi1(1)))

      psi_1valence_core_atom = ave
      return
      end



*     ***********************************
*     *                                 *
*     *         psi_1vloc_atom          *
*     *                                 *
*     ***********************************
      real*8 function psi_1vloc_atom()
      implicit none

#include "bafdecls.fh"
#include "psi.fh"

*     **** local variables ***
      real*8 ave

*     **** external functions ****
      real*8   psp_vloc_atom
      external psp_vloc_atom

      ave =  psp_vloc_atom(ispin,neq,dcpl_mb(psi1(1)))

      psi_1vloc_atom = ave
      return
      end

        

*     ***********************************  
*     *                                 *  
*     *       psi_1ncmp_vloc            *
*     *                                 *
*     ***********************************
      real*8 function psi_1ncmp_vloc()
      implicit none
         
#include "bafdecls.fh"
#include "psi.fh"
         
*     **** local variables ***    
      real*8 ave 

*     **** external functions ****
      real*8   psp_ncmp_vloc
      external psp_ncmp_vloc
     
      ave =  psp_ncmp_vloc(ispin)

      psi_1ncmp_vloc = ave
      return
      end



*     ***********************************
*     *                                 *
*     *         psi_1hartree_atom      *  
*     *                                 *  
*     ***********************************
      real*8 function psi_1hartree_atom()     
      implicit none
         
#include "bafdecls.fh"
#include "psi.fh"
         
*     **** local variables ***
      real*8 ave
         
*     **** external functions ****
      real*8   psp_hartree_atom
      external psp_hartree_atom

      ave =  psp_hartree_atom(ispin,neq,dcpl_mb(psi1(1)))
         
      psi_1hartree_atom = ave
      return
      end


*     ***********************************
*     *                                 *
*     *      psi_1hartree_cmp_cmp       *
*     *                                 *
*     ***********************************
      real*8 function psi_1hartree_cmp_cmp()
      implicit none
         
#include "psi.fh"
         
*     **** external functions ****
      real*8   psp_hartree_cmp_cmp
      external psp_hartree_cmp_cmp
         
      psi_1hartree_cmp_cmp = psp_hartree_cmp_cmp(ispin)
      return
      end

*     ***********************************
*     *                                 *
*     *         psi_1hartree_cmp_pw     *
*     *                                 *
*     ***********************************
      real*8 function psi_1hartree_cmp_pw()
      implicit none
         
#include "bafdecls.fh"
#include "psi.fh"
         
*     **** external functions ****
      real*8   psp_hartree_cmp_pw
      external psp_hartree_cmp_pw

      psi_1hartree_cmp_pw = psp_hartree_cmp_pw(ispin,dcpl_mb(dng1(1)),
     >                                               dbl_mb(rho1(1)))
      return
      end


c*     ***********************************
c*     *                                 *
c*     *         dng_1vlpaw_pw           *
c*     *                                 *
c*     ***********************************
c      real*8 function dng_1vlpaw_pw()
c      implicit none
c
c#include "bafdecls.fh"
c#include "psi.fh"
c
c*     **** external functions ****
c      real*8   electron_dng_vlpaw_ave
c      external electron_dng_vlpaw_ave
c
c      dng_1vlpaw_pw = electron_dng_vlpaw_ave(dcpl_mb(dng1(1)))
c
c      return
c      end



*     ***********************************
*     *                                 *
*     *         psi_1xc_atom            *
*     *                                 *
*     ***********************************
      subroutine psi_1xc_atom(exc,pxc)
      implicit none
      real*8 exc,pxc

#include "bafdecls.fh"
#include "psi.fh"

      call psp_xc_atom(ispin,neq,dcpl_mb(psi1(1)),exc,pxc)
      !call D1dB_SumAll(exc)
      !call D1dB_SumAll(pxc)
      return
      end


*     ***********************************
*     *                                 *
*     *         psi_1qlm_atom           *
*     *                                 *
*     ***********************************
      subroutine psi_1qlm_atom()
      implicit none

#include "bafdecls.fh"
#include "psi.fh"

      call  psp_qlm_atom(ispin,neq,dcpl_mb(psi1(1)))
      return
      end




*     ***********************************
*     *					*
*     *		psi_1vl 		*
*     *					*
*     ***********************************
      real*8 function psi_1vl()
      implicit none

#include "bafdecls.fh"
#include "psi.fh"


*     **** external functions ****
      real*8   electron_psi_vl_ave
      external electron_psi_vl_ave
 
      psi_1vl = electron_psi_vl_ave(dcpl_mb(psi1(1)),dbl_mb(rho1(1)))

      return
      end



*     ***********************************
*     *                                 *
*     *         dng_1apc                *
*     *                                 *
*     ***********************************
      subroutine dng_1apc(Eapc,Papc)
      implicit none
      real*8 Eapc,Papc

#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     **** local variables ****
      integer vtmp(2)
      real*8  ftmp(3)

      if (.not.BA_push_get(mt_dcpl,npack0,'vtmp',vtmp(2),vtmp(1)))
     >   call errquit('dng_1apc: out of stack memory',0,MA_ERR)

      call pspw_V_APC(ispin,ne,dcpl_mb(dng1(1)),
     >                         dcpl_mb(vtmp(1)),
     >                         Eapc,Papc,.false.,ftmp)

      if (.not.BA_pop_stack(vtmp(2)))
     >  call errquit('dng_1apc popping stack',1,MA_ERR)

      return
      end


*     ***********************************
*     *                                 *
*     *         dng_1vl_mm              *
*     *                                 *
*     ***********************************
      real*8 function dng_1vl_mm()
      implicit none
      
#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     **** local variables ****
      integer nx,ny,nz,n2ft3d,vltmp(2),r_grid(2)
      real*8  elocal,esum,dv
         
*     **** external functions ****
      integer  control_version
      external control_version
      real*8   lattice_omega
      external lattice_omega

      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      dv = lattice_omega()/dble(nx*ny*nz)

      if (.not.BA_push_get(mt_dbl,(n2ft3d),'vltmp',vltmp(2),vltmp(1)))
     >   call errquit('dng_1_vl_mm: out of stack memory',0,MA_ERR)

*     **** average Kohn-Sham v_local energy ****
      call v_local_mm(dbl_mb(vltmp(1)))
      call Pack_cc_dot(0,dcpl_mb(dng1(1)),dbl_mb(vltmp(1)),elocal)

*     *** add in long range part ****
      if (control_version().eq.4) then
         if (.not.BA_push_get(mt_dbl,(3*n2ft3d),'r_grid',
     >                       r_grid(2),r_grid(1)))
     >      call errquit('dng_1_vl_mm: out of stack memory',1,MA_ERR)
         call lattice_r_grid(dbl_mb(r_grid(1)))

         call v_lr_local_mm(dbl_mb(r_grid(1)),dbl_mb(vltmp(1)))
         call D3dB_rr_dot(1,dbl_mb(rho1(1)),dbl_mb(vltmp(1)),esum)
         elocal = elocal + esum*dv
         if (.not.BA_pop_stack(r_grid(2)))
     >      call errquit('dng_1_vl_mm: popping stack',0,MA_ERR)
      end if

      if (.not.BA_pop_stack(vltmp(2)))
     >  call errquit('dng_1_vl_mm: popping stack',1,MA_ERR)

         
      dng_1vl_mm = elocal
      return
      end





*     ***********************************
*     *					*
*     *		psi_1vnl 		*
*     *					*
*     ***********************************
      real*8 function psi_1vnl()
      implicit none

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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
      value = BA_push_get(mt_dbl,(3*n2ft3d),'r_grid',
     >                       r_grid(2),r_grid(1))
      value = value.and. 
     >        BA_push_get(mt_dbl,(3*n2ft3d),'rho',
     >                       rho(2),rho(1))
      if (.not. value) call errquit(
     >     'rho_1Fcharge: out of stack memory',0, MA_ERR)
      call Parallel_shared_vector_zero(.true.,3*n2ft3d,dbl_mb(rho(1)))


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
      value =           BA_pop_stack(rho(2))
      value = value.and.BA_pop_stack(r_grid(2))
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "psi.fh"

*     **** external functions ****
      real*8   electron_pxc
      external electron_pxc

      rho_1pxc = electron_pxc(dbl_mb(rho1(1)))
      return
      end


*     ***********************************
*     *                                 *
*     *         psi_1meta_gga_pxc       *
*     *                                 *
*     ***********************************
      real*8 function psi_1meta_gga_pxc()
      implicit none

#include "bafdecls.fh"
#include "psi.fh"

*     **** external functions ****
      real*8   nwpw_meta_gga_pxc
      external nwpw_meta_gga_pxc

      psi_1meta_gga_pxc = nwpw_meta_gga_pxc(ispin,neq,dcpl_mb(psi1(1)))
      return
      end



*     ***********************************
*     *					*
*     *		dng_1ehartree           *
*     *					*
*     ***********************************
      real*8 function dng_1ehartree()
      implicit none

#include "bafdecls.fh"
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
*     * 	 psi_1vl_cosmo          *
*     *					*
*     ***********************************
      real*8 function psi_1vl_cosmo()
      implicit none

#include "bafdecls.fh"
#include "psi.fh"
#include "errquit.fh"

*     **** local variables ****
      logical ok
      integer r_grid(2),tmp1(2)
      integer nx,ny,nz,n2ft3d
      real*8  elocal,e1,e2,dv

*     **** external functions ****
      logical  nwpw_cosmo1_on,nwpw_cosmo2_on
      external nwpw_cosmo1_on,nwpw_cosmo2_on
      real*8   lattice_omega,nwpw_cosmo_EQelcq
      external lattice_omega,nwpw_cosmo_EQelcq

      if (nwpw_cosmo1_on()) then
         call D3dB_n2ft3d(1,n2ft3d)
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         dv = lattice_omega()/dble(nx*ny*nz)

         ok=BA_push_get(mt_dbl,(3*n2ft3d),'r_grid',r_grid(2),r_grid(1))
         ok=ok.and.BA_push_get(mt_dbl,(n2ft3d),'tmp1',tmp1(2),tmp1(1))
         if (.not.ok) call errquit("psi_1vl_cosmo: push stack",0,MA_ERR)

         call lattice_r_grid(dbl_mb(r_grid(1)))

         call v_local_cosmo(dbl_mb(tmp1(1)))
         call Pack_cc_dot(0,dcpl_mb(dng1(1)),dbl_mb(tmp1(1)),elocal)

         call v_lr_local_cosmo(dbl_mb(r_grid(1)),dbl_mb(tmp1(1)))
         call D3dB_rr_dot(1,dbl_mb(rho1(1)),
     >                    dbl_mb(tmp1(1)),e1)
         call D3dB_rr_dot(1,dbl_mb(rho1(1)+(ispin-1)*n2ft3d),
     >                    dbl_mb(tmp1(1)),e2)
         elocal = elocal + (e1+e2)*dv

         ok =        BA_pop_stack(tmp1(2))
         ok = ok.and.BA_pop_stack(r_grid(2))
         if (.not.ok) call errquit("psi_1vl_cosmo: pop stack",0,MA_ERR)
      else if (nwpw_cosmo2_on()) then
         elocal = nwpw_cosmo_EQelcq()
      else
         elocal = 0.0d0
      end if

      psi_1vl_cosmo = elocal
      return
      end 


*     ***********************************
*     *					*
*     *		rho_1dngen_APC   	*
*     *					*
*     ***********************************
      subroutine rho_1dngen_APC()
      implicit none

#include "bafdecls.fh"
#include "psi.fh"

      call pspw_dngen_APC(ispin,ne,dbl_mb(rho1(1)),.false.)
      return
      end 


*     ***********************************
*     *					*
*     *		psi_2toelectron		*
*     *					*
*     ***********************************
      subroutine psi_2toelectron()
      implicit none

#include "bafdecls.fh"
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
#include "bafdecls.fh"
#include "psi.fh"

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
            


         value = BA_pop_stack(MM(2))
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

#include "bafdecls.fh"
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
         value = BA_push_get(mt_dbl,n*n,'MM',MM(2),MM(1))
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



         value = BA_pop_stack(MM(2))
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

#include "bafdecls.fh"
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

*     **** get the tangent THpsi = Hpsi - psi1*tmp1 , TY = HY - Y*Y'HY ****
      call electron_get_Tgradient(dcpl_mb(psi1(1)),
     >                             dbl_mb(tmp1(1)),
     >                            THpsi)
      
      if (.not.Dneall_m_pop_stack(tmp1))
     > call errquit('psi_1get_Tgradient:error popping stack',1,
     >     MA_ERR)

      return
      end


*     ***********************************
*     *                                 *
*     *         psi_1get_remainder0     *
*     *                                 *
*     ***********************************
*
*     Must be called after call to psi_1get_Tgradient0 or psi_1get_Tgradient
*
      subroutine psi_1get_remainder0(Eremainder)
      implicit none
      real*8 Eremainder

#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

*     **** external functions ****
      real*8   electron_energy_remainder
      external electron_energy_remainder

      Eremainder =  electron_energy_remainder(dcpl_mb(psi1(1)),
     >                        dbl_mb(rho1(1)),
     >                        dcpl_mb(dng1(1)),
     >                        dbl_mb(rho1_all(1)))

      return
      end

*     ***********************************
*     *                                 *
*     *         psi_2get_remainder0     *
*     *                                 *
*     ***********************************
*     
*     Must be called after call to psi_2get_Tgradient0 or psi_2get_Tgradient
*     
      subroutine psi_2get_remainder0(Eremainder)
      implicit none
      real*8 Eremainder
      
#include "bafdecls.fh"                 
#include "errquit.fh"
#include "psi.fh"
      
*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all
            
*     **** external functions ****
      real*8   electron_energy_remainder
      external electron_energy_remainder
      
      Eremainder =  electron_energy_remainder(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)),
     >                        dbl_mb(rho2_all(1)))
      return
      end


*     ***********************************
*     *                                 *
*     *         psi_1get_Tgradient0     *
*     *                                 *
*     ***********************************

*     THpsi = Hpsi - Y*Y^t*Hpsi ! used by Grassman minimizers
*
      subroutine psi_1get_Tgradient0(THpsi,Eout)
      implicit none
      complex*16 THpsi(*)
      real*8 Eout

#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"


*     **** local variables ****
      integer tmp1(2),i,n
      logical value

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack

      real*8   electron_energy0
      external electron_energy0

      if (.not.Dneall_m_push_get(0,tmp1))
     >   call errquit('out of stack memory in psi_1get_Tradient0',0,
     >       MA_ERR)

      call electron_run0(dcpl_mb(psi1(1)))
      Eout =  electron_energy0(dcpl_mb(psi1(1)))

      call electron_gen_hml(dcpl_mb(psi1(1)),
     >                       dbl_mb(tmp1(1)))

*     **** get the tangent THpsi = Hpsi - psi1*tmp1 , TY = HY - Y*Y'HY ****
      call electron_get_Tgradient(dcpl_mb(psi1(1)),
     >                             dbl_mb(tmp1(1)),
     >                            THpsi)

      if (.not.Dneall_m_pop_stack(tmp1))
     > call errquit('psi_1get_Tgradient0:error popping stack',1,
     >     MA_ERR)

      return
      end


*     ***********************************
*     *                                 *
*     *         psi_1get_STgradient     *
*     *                                 *
*     ***********************************

*     THpsi = Hpsi - Y*Y^t*Hpsi ! used by Grassman minimizers
*
      subroutine psi_1get_STgradient(Rpsi,THpsi,Eout)
      implicit none
      complex*16 Rpsi(*)
      complex*16 THpsi(*)
      real*8 Eout

#include "bafdecls.fh"
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
     >   call errquit('out of stack memory in psi_1get_STradient',0,
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

*     **** get the tangent Rpsi = Hpsi - Spsi1*tmp1 , TY = HY - SY*Y'HY ****
      call psp_overlap_S(ispin,neq,dcpl_mb(psi1(1)),dcpl_mb(spsi1(1)))
      call electron_get_Tgradient(dcpl_mb(spsi1(1)),
     >                             dbl_mb(tmp1(1)),
     >                            Rpsi)

*     **** THpsi = Rpsi - psi1 * <spsi1|Rpsi> ****
      call Grsm_gg_Copy(npack1,neq(1)+neq(2),Rpsi,THpsi)
      call Dneall_ffm_sym_Multiply(0,dcpl_mb(spsi1(1)),Rpsi,npack1,
     >                            dbl_mb(tmp1(1)))
      call Dneall_fmf_Multiply(0,dcpl_mb(psi1(1)),npack1,
     >                         dbl_mb(tmp1(1)),
     >                         -1.0d0,THpsi,1.0d0)

      if (.not.Dneall_m_pop_stack(tmp1))
     > call errquit('psi_1get_STgradient:error popping stack',1,
     >     MA_ERR)

      return
      end


*     ***********************************
*     *                                 *
*     *         psi_2get_STgradient     *
*     *                                 *
*     ***********************************

*     THpsi = Hpsi - Y*Y^t*Hpsi ! used by Grassman minimizers
*
      subroutine psi_2get_STgradient(option,Rpsi,THpsi,Eout)
      implicit none
      integer    option
      complex*16 Rpsi(*)
      complex*16 THpsi(*)
      real*8 Eout

#include "bafdecls.fh"
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
     >   call errquit('out of stack memory in psi_2get_STradient',0,
     >       MA_ERR)

      if (option.le.1) then
         call electron_run(dcpl_mb(psi2(1)),
     >                   dbl_mb(rho2(1)),
     >                  dcpl_mb(dng2(1)),
     >                   dbl_mb(rho2_all(1)),
     >                  occupation_on,dbl_mb(occ2(1)))
      end if
      Eout =  electron_energy(dcpl_mb(psi2(1)),
     >                        dbl_mb(rho2(1)),
     >                        dcpl_mb(dng2(1)),
     >                        dbl_mb(rho2_all(1)),
     >                        occupation_on,dbl_mb(occ2(1)))
      call electron_gen_hml(dcpl_mb(psi2(1)),
     >                       dbl_mb(tmp1(1)))

*     **** get the tangent Rpsi = Hpsi - Spsi2*tmp1 , TY = HY - SY*Y'HY ****
      call psp_overlap_S(ispin,neq,dcpl_mb(psi2(1)),dcpl_mb(spsi1(1)))
      call electron_get_Tgradient(dcpl_mb(spsi1(1)),
     >                             dbl_mb(tmp1(1)),
     >                            Rpsi)

*     **** THpsi = Rpsi - psi2 * <spsi2|Rpsi> ****
      call Grsm_gg_Copy(npack1,neq(1)+neq(2),Rpsi,THpsi)
      call Dneall_ffm_sym_Multiply(0,dcpl_mb(spsi1(1)),Rpsi,npack1,
     >                            dbl_mb(tmp1(1)))
      call Dneall_fmf_Multiply(0,dcpl_mb(psi2(1)),npack1,
     >                         dbl_mb(tmp1(1)),
     >                         -1.0d0,THpsi,1.0d0)

      if (.not.Dneall_m_pop_stack(tmp1))
     > call errquit('psi_2get_STgradient:error popping stack',1,
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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
#include "bafdecls.fh"
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


!$OMP BARRIER

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
!$OMP BARRIER

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
*     *                                 *
*     *         psi_2get_Tgradient0     *
*     *                                 *
*     ***********************************
      subroutine psi_2get_Tgradient0(option,THpsi,Eout)
      implicit none
      integer    option
      complex*16 THpsi(*)
      real*8     Eout

#include "errquit.fh"
#include "bafdecls.fh"
#include "psi.fh"


*     **** local variables ****
      integer tmp1(2)

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack

      real*8   electron_energy0
      external electron_energy0


!$OMP BARRIER

      if (.not.Dneall_m_push_get(0,tmp1))
     >   call errquit('out of stack memory in psi_2get_Tradient0',0,
     >       MA_ERR)

      if (option.le.1) then
        call electron_run0(dcpl_mb(psi2(1)))
      end if

      Eout =  electron_energy0(dcpl_mb(psi2(1)))
!$OMP BARRIER

      call electron_gen_hml(dcpl_mb(psi2(1)),
     >                       dbl_mb(tmp1(1)))
      call electron_get_Tgradient(dcpl_mb(psi2(1)),
     >                             dbl_mb(tmp1(1)),
     >                             THpsi)

      if (.not. Dneall_m_pop_stack(tmp1))
     >call errquit('psi_2get_Tgradient0:error popping stack',1,MA_ERR)

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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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
#include "bafdecls.fh"
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
#include "bafdecls.fh"
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

#include "bafdecls.fh"
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
*     *                                 *
*     *        psi_1get_TMgradient0     *
*     *                                 *
*     ***********************************
      subroutine psi_1get_TMgradient0(THpsi,Eout)
      implicit none
      complex*16 THpsi(*)
      real*8     Eout

#include "bafdecls.fh"
#include "psi.fh"

*     **** external functions ****
      real*8   electron_energy0
      external electron_energy0


      call electron_run0(dcpl_mb(psi1(1)))
      Eout =  electron_energy0(dcpl_mb(psi1(1)))

      call electron_get_TMgradient(dcpl_mb(psi1(1)),
     >                             THpsi)

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

#include "bafdecls.fh"
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
*     *                                 *
*     *         psi_2get_TMgradient0    *
*     *                                 *
*     ***********************************
      subroutine psi_2get_TMgradient0(option,THpsi,Eout)
      implicit none
      integer    option
      complex*16 THpsi(*)
      real*8     Eout

#include "bafdecls.fh"
#include "psi.fh"

*     **** external functions ****
      real*8   electron_energy0
      external electron_energy0

      if (option.le.1) then
        call electron_run0(dcpl_mb(psi2(1)))
      end if

      Eout =  electron_energy0(dcpl_mb(psi2(1)))

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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "psi.fh"

      call geodesic_Gtransport(t,dcpl_mb(psi1(1)),G0)

      return
      end



*     ***********************************
*     *                                 *
*     *         psi_geodesic_energy0    *
*     *                                 *
*     ***********************************
*
*    This function follows a geodesic but without updating
* the densities.

      real*8 function psi_geodesic_energy0(t)
      implicit none
      real*8 t

#include "bafdecls.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all

      real*8 e_new
*     **** external functions ****
      real*8   electron_energy0
      external electron_energy0

      call geodesic_get(t,dcpl_mb(psi1(1)),
     >                    dcpl_mb(psi2(1)))

      call electron_run0(dcpl_mb(psi2(1)))

      e_new =  electron_energy0(dcpl_mb(psi2(1)))
        
      psi_geodesic_energy0 = e_new
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "psi.fh"

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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

      if (occupation_on) 
     >   call dcopy((ne(1)+ne(2)),dbl_mb(occ1(1)),1,dbl_mb(occ2(1)),1)

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

#include "bafdecls.fh"
#include "psi.fh"

*     **** external functions ****
      real*8   electron_eorbit
      external electron_eorbit


      call geodesic2_transport(t,dcpl_mb(psi1(1)),
     >                           dcpl_mb(psi2(1)))
      if (occupation_on)
     >   call dcopy((ne(1)+ne(2)),dbl_mb(occ1(1)),1,dbl_mb(occ2(1)),1)

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

#include "bafdecls.fh"
#include "psi.fh"

      integer taskid,MASTER
      parameter (MASTER=0)
c     real*8 sum1,sum2
     
      call Parallel_taskid(taskid)

      call geodesic2_get(t,dcpl_mb(psi1(1)),
     >                    dcpl_mb(psi2(1)))
      if (occupation_on)
     >   call dcopy((ne(1)+ne(2)),dbl_mb(occ1(1)),1,dbl_mb(occ2(1)),1)
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

#include "bafdecls.fh"
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
     >                        dte)

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

        if (pawexist) then
           call psp_overlap_S(ispin,neq,
     >                        dcpl_mb(psi1(1)),
     >                        dcpl_mb(spsi1(1)))
           call psi_lmbda_paw(ispin,neq,nemaxq,npack1,
     >                        dcpl_mb(spsi1(1)),
     >                        dcpl_mb(psi2(1)),
     >                        dte,
     >                        dbl_mb(lmd(1)),
     >                        dbl_mb(tmp_L(1)),ierr)
        else if (occupation_on) then
           call psi_lmbda2(ispin,ne,nemaxq,npack1,
     >                     dcpl_mb(psi1(1)),
     >                     dcpl_mb(psi2(1)),
     >                     dte,dbl_mb(occ1(1)),
     >                     dbl_mb(lmd(1)),
     >                     dbl_mb(tmp_L(1)),ierr)

        else if (pspw_SIC()) then
           call psi_lmbda_sic(ispin,ne,nemaxq,npack1,
     >                 dcpl_mb(psi1(1)),dcpl_mb(psi2(1)),dte,
     >                 dbl_mb(lmd(1)),
     >                 dbl_mb(tmp_L(1)),ierr)
        else
           call psi_lmbda(ispin,neq,nemaxq,npack1,
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     **** local variables ****
      logical value
      integer r_grid(2),tmp(2)

*     **** external functions ****
      integer  control_version
      external control_version
      logical  psp_U_psputerm,pspw_V_APC_on
      external psp_U_psputerm,pspw_V_APC_on
      logical  Dneall_m_push_get, Dneall_m_pop_stack
      external Dneall_m_push_get, Dneall_m_pop_stack

c     call electron_gen_psi_r(dcpl_mb(psi1(1)))
c     call electron_gen_densities(dcpl_mb(psi1(1)),
c    >                             dbl_mb(rho1(1)),
c    >                            dcpl_mb(dng1(1)))

      call f_vlocal(dcpl_mb(dng1(1)),fion)

      if (control_version().eq.4) then
          value = BA_push_get(mt_dbl,(2*nfft3d),'tmp',
     >                        tmp(2),tmp(1))
          value = value.and.
     >            BA_push_get(mt_dbl,(6*nfft3d),'r_grid',
     >                        r_grid(2),r_grid(1))
         if (.not. value) 
     >      call errquit('psi_1force:out of stack memory',0, MA_ERR)
          !call dcopy(2*nfft3d,0.0d0,0,dbl_mb(tmp(1)),1)
          call Parallel_shared_vector_zero(.true.,
     >                                     2*nfft3d,dbl_mb(tmp(1)))

          call D3dB_rr_Sum(1,dbl_mb(rho1(1)),
     >                       dbl_mb(rho1(1)+(ispin-1)*2*nfft3d),
     >                       dbl_mb(tmp(1)))
          call lattice_r_grid(dbl_mb(r_grid(1)))
          call grad_v_lr_local(dbl_mb(r_grid(1)),
     >                         dbl_mb(tmp(1)),
     >                         fion)

          value = BA_pop_stack(r_grid(2))
          value = value.and.BA_pop_stack(tmp(2))
         if (.not.value)
     >      call errquit('psi_1force:error popping stack memory',0,
     &                   MA_ERR)
      end if

*     *** APC force here ***
      if (pspw_V_APC_on()) then
         call pspw_force_APC(ispin,ne,dcpl_mb(dng1(1)),fion)
      end if

      call f_vnonlocal(ispin,
     >                 neq,
     >                 dcpl_mb(psi1(1)),
     >                 fion,
     >                 occupation_on,dbl_mb(occ1(1)))

      if (pawexist) then
         value = Dneall_m_push_get(0,tmp)
         if (.not.value)
     >      call errquit('psi_1force:out of stack memory',0,MA_ERR)

         call psi_1toelectron()
         call electron_gen_hml(dcpl_mb(psi1(1)),dbl_mb(tmp(1)))

         call psp_paw_overlap_fion(ispin,
     >                             dbl_mb(tmp(1)),
     >                             dcpl_mb(psi1(1)),
     >                             fion)

         value = Dneall_m_pop_stack(tmp)
         if (.not.value)
     >      call errquit('psi_1force:error popping stack',0,MA_ERR)
      end if

      if (psp_U_psputerm()) then
         call f_psp_U_v_nonlocal(ispin,
     >                 neq,
     >                 dcpl_mb(psi1(1)),
     >                 fion,
     >                 occupation_on,dbl_mb(occ1(1)),.true.)
      end if

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

#include "bafdecls.fh"
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
          value = BA_push_get(mt_dbl,(2*nfft3d),'tmp',
     >                        tmp(2),tmp(1))
          value = value.and.
     >            BA_push_get(mt_dbl,(6*nfft3d),'r_grid',
     >                        r_grid(2),r_grid(1))
         if (.not. value) call errquit('out of stack memory',0, MA_ERR)
          call dcopy(2*nfft3d,0.0d0,0,dbl_mb(tmp(1)),1)

          call D3dB_rr_Sum(1,dbl_mb(rho1(1)),
     >                       dbl_mb(rho1(1)+(ispin-1)*2*nfft3d),
     >                       dbl_mb(tmp(1)))
          call lattice_r_grid(dbl_mb(r_grid(1)))
          call grad_v_lr_local(dbl_mb(r_grid(1)),
     >                         dbl_mb(tmp(1)),
     >                         fion)

          value = BA_pop_stack(r_grid(2))
          value = value.and.BA_pop_stack(tmp(2))
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

#include "bafdecls.fh"
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
*     *                                 *
*     *   psi_1force_psp_U_v_nonlocal   *
*     *                                 *
*     ***********************************
      subroutine psi_1force_psp_U_v_nonlocal(fion)
      implicit none
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

c     **** external functions ****
      logical  psp_U_psputerm
      external psp_U_psputerm

      if (psp_U_psputerm()) then
         call f_psp_U_v_nonlocal(ispin,
     >                 neq,
     >                 dcpl_mb(psi1(1)),
     >                 fion,
     >                 occupation_on,dbl_mb(occ1(1)),.true.)
      end if

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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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


#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "psi.fh"


      call electron_gen_hml(dcpl_mb(psi2(1)),dbl_mb(hml(1)))

      return
      end

*     ***********************************
*     *					*
*     *		psi_hml_value     	*
*     *					*
*     ***********************************
      real*8  function psi_hml_value(ms,i,j)
      implicit none
      integer ms
      integer i,j

#include "bafdecls.fh"
#include "psi.fh"

      psi_hml_value = dbl_mb(hml(1) + (i-1)+(j-1)*ne(ms)
     >                              + (ms-1)*ne(1)*ne(1))
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "psi.fh"

      integer ms,i,indx1,indx2

!$OMP MASTER
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
!$OMP END MASTER
!$OMP BARRIER
      return
      end

*     ***********************************
*     *                                 *
*     *     psi_1assending_occupation   *
*     *                                 *
*     ***********************************
*
*     makes sure the occupation follows an assending order 
*
      subroutine psi_1assending_occupation()
      implicit none
         
#include "bafdecls.fh"
#include "psi.fh"
     
      if (dbl_mb(occ1(1)).lt.dbl_mb(occ1(1)+ne(1)-1)) then
         call psi_1reverse_occupation()
      end if
      return
      end

*     ***********************************
*     *                                 *
*     *     psi_1desending_occupation   *
*     *                                 *
*     ***********************************
*
*     makes sure the occupation follows a desending order 
*
      subroutine psi_1desending_occupation()
      implicit none
         
#include "bafdecls.fh"
#include "psi.fh"
     
      if (dbl_mb(occ1(1)+ne(1)-1).lt.dbl_mb(occ1(1))) then
         call psi_1reverse_occupation()
      end if
      return
      end


*     ***********************************
*     *                                 *
*     *      psi_0define_occupation     *
*     *                                 *
*     ***********************************
      subroutine psi_0define_occupation(initial_alpha,use_hml,
     >                                  ispin,ne,eig,hml,occ,
     >                                  smeartype,smearkT,
     >                                  smearfermi,smearcorrection)
      implicit none
      real*8  initial_alpha
      logical use_hml
      integer ispin,ne(2)
      real*8  eig(*),hml(*),occ(*)
      integer smeartype
      real*8 smearkT,smearfermi(2),smearcorrection

*     **** local variables ****
      integer it,itmax
      parameter (itmax=50)

      integer ms,nb,n,shift1,shift2,occ1_tag,ndiff
      real*8 e,x,kT,f,g,alpha,pi,f0
      real*8 ZZ,Z(2),Zlower,Zmid,Zupper,elower,emid,eupper
      real*8 flower,fmid,fupper,lmbda

*     **** external functions ****
      integer  control_multiplicity
      real*8   control_TotalCharge,ion_TotalCharge_qm
      real*8   psi_occ_distribution,control_fractional_alpha
      external control_multiplicity
      external control_TotalCharge,ion_TotalCharge_qm
      external psi_occ_distribution,control_fractional_alpha

      ZZ  = ion_TotalCharge_qm() - control_TotalCharge()

!$OMP MASTER
      if (use_hml) then
         do ms=1,ispin
            shift1 = (ms-1)*ne(1)
            shift2 = (ms-1)*ne(1)*ne(1)
            do n=1,ne(ms)
               eig(n+shift1) = hml(n+(n-1)*ne(ms)+shift2)
            end do
         end do
      end if
      smearfermi(1)   = 0.0d0
      smearfermi(2)   = 0.0d0
      smearcorrection = 0.0d0

      
         if (initial_alpha.lt.0.0d0) then
            alpha = control_fractional_alpha()
         else
            alpha = initial_alpha
         end if
         kT    = smearkT
         !ZZ  = ion_TotalCharge_qm() - control_TotalCharge()

         if (dabs(ZZ).lt.1.0d-9) go to 98

         if (ispin.eq.2) then
            ndiff = control_multiplicity() - 1
            Z(1) = 0.5d0*(ZZ+ndiff)
            Z(2) = 0.5d0*(ZZ-ndiff)
         else
            Z(1) = 0.5d0*ZZ
            Z(2) = 0.0d0
         end if

         pi    = 4.0d0*datan(1.0d0)
         !if (initial) alpha = 1.0d0
         if (smeartype.le.0) alpha = 0.0d0

*        **** outer loop over spins ****
         smearcorrection = 0.0d0
         do ms=1,ispin

*           **** find eupper and elower ****
            elower =  9.9d12
            eupper = -9.9d12
            shift1 = (ms-1)*ne(1) + 1
            do n=1,ne(ms)
              e       = eig(shift1)
              if (e.lt.elower) elower = e
              if (e.gt.eupper) eupper = e
              shift1  = shift1 + 1
            end do

*           **** find fermi level ****
            Zlower = 0.0d0
            Zupper = 0.0d0
            shift1 = (ms-1)*ne(1) + 1
            do n=1,ne(ms)
              e = eig(shift1)
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
     >            'psi_0define_occupation:Fermi energy not found',ms,0)

            it = 0
  20        it = it + 1
            emid = 0.5d0*(elower + eupper)
            Zmid = 0.0d0
            shift1 = (ms-1)*ne(1) + 1
            do n=1,ne(ms)
              e = eig(shift1)
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
            shift1 = (ms-1)*ne(1) + 1
            shift2 = (ms-1)*ne(1) + 1
            do n=1,ne(ms)
              e = eig(shift1)
              x = (e - smearfermi(ms))/kT
              f = psi_occ_distribution(smeartype,x)
              f0 = occ(shift2)

              occ(shift2) = (1.0d0-alpha)*f0 + alpha*f

              if (smeartype.eq.1) then
                 if (  (occ(shift2)       .gt.1.0d-6) .and.
     >               ( (1.0d0-occ(shift2)).gt.1.0d-6)) then
                smearcorrection = smearcorrection
     >           + kT*( occ(shift2)*log(occ(shift2))
     >           + (1.0d0-occ(shift2))*log(1.0d0-occ(shift2)) )
                 end if
              else if (smeartype.eq.2) then
                smearcorrection
     >              = smearcorrection
     >              - kT*dexp(-x*x)/(4.0d0*dsqrt(pi))
              else if (smeartype.eq.4) then
                smearcorrection
     >              = smearcorrection
     >              - kT*dexp(-(x+dsqrt(0.5d0))*(x+dsqrt(0.5d0)))
     >              * (1.0d0 + dsqrt(2.0d0) * x)
     >              / (2.0d0 * dsqrt(pi))
              end if

              shift1  = shift1 + 1
              shift2  = shift2 + 1
            end do

         end do !** ms***
         if (ms.eq.1) smearcorrection=smearcorrection+smearcorrection

         go to 99


  98     continue
         smearcorrection = 0.0d0
         call dcopy(ne(1)+ne(2),0.0d0,0,occ,1)

  99     continue

!$OMP END MASTER
!$OMP BARRIER

      return
      end






*     ***********************************
*     *                                 *
*     *      psi_1define_occupation     *
*     *                                 *
*     ***********************************
      subroutine psi_1define_occupation(initial_alpha,use_hml)
      implicit none
      real*8  initial_alpha
      logical use_hml

#include "bafdecls.fh"
#include "psi.fh"

*     **** local variables ****
      integer it,itmax
      parameter (itmax=50)

      integer ms,nb,n,shift1,shift2,occ1_tag,ndiff
      real*8 e,x,kT,f,g,alpha,pi,f0
      real*8 ZZ,Z(2),Zlower,Zmid,Zupper,elower,emid,eupper
      real*8 flower,fmid,fupper,lmbda

*     **** external functions ****
      integer  control_multiplicity
      real*8   control_TotalCharge,ion_TotalCharge_qm
      real*8   psi_occ_distribution,control_fractional_alpha
      external control_multiplicity
      external control_TotalCharge,ion_TotalCharge_qm
      external psi_occ_distribution,control_fractional_alpha

      ZZ  = ion_TotalCharge_qm() - control_TotalCharge()

!$OMP MASTER
      if (use_hml) then
         do ms=1,ispin
            shift1 = eig(1) + (ms-1)*ne(1)
            shift2 = hml(1) + (ms-1)*ne(1)*ne(1)
            do n=1,ne(ms)
               dbl_mb(shift1+n-1) = dbl_mb(shift2+(n-1)+(n-1)*ne(ms))
            end do
         end do
      end if
      smearfermi(1)   = 0.0d0
      smearfermi(2)   = 0.0d0
      smearcorrection = 0.0d0

      if (occupation_on) then
         if (initial_alpha.lt.0.0d0) then
            alpha = control_fractional_alpha()
         else
            alpha = initial_alpha
         end if
         kT    = smearkT
         !ZZ  = ion_TotalCharge_qm() - control_TotalCharge()

         if (dabs(ZZ).lt.1.0d-9) go to 98

         if (ispin.eq.2) then
            ndiff = control_multiplicity() - 1
            Z(1) = 0.5d0*(ZZ+ndiff)
            Z(2) = 0.5d0*(ZZ-ndiff)
         else
            Z(1) = 0.5d0*ZZ
            Z(2) = 0.0d0
         end if
       
         pi    = 4.0d0*datan(1.0d0)
         !if (initial) alpha = 1.0d0
         if (smeartype.le.0) alpha = 0.0d0

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
              f0 = dbl_mb(shift2)

              dbl_mb(shift2) = (1.0d0-alpha)*f0 + alpha*f

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
              else if (smeartype.eq.4) then
                smearcorrection
     >              = smearcorrection
     >              - kT*dexp(-(x+dsqrt(0.5d0))*(x+dsqrt(0.5d0)))
     >              * (1.0d0 + dsqrt(2.0d0) * x)
     >              / (2.0d0 * dsqrt(pi))
              end if

              shift1  = shift1 + 1
              shift2  = shift2 + 1
            end do

         end do !** ms***
         if (ms.eq.1) smearcorrection=smearcorrection+smearcorrection

         go to 99

  98     continue
         smearcorrection = 0.0d0
         do ms=1,ispin
            shift2 = occ1(1) + (ms-1)*ne(1)
            call dcopy(ne(ms),0.0d0,0,dbl_mb(shift2),1)
         end do

  99     continue

      end if
!$OMP END MASTER
!$OMP BARRIER

      return
      end


c  set nwpw:fractional_smeartype 1 #0-none, 1-Fermi-Dirac, 2-Gaussian, 3-Hermite
c                                   4-Marzari-Vanderbilt
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
      else if (smeartype.eq.4) then
         f = dexp(-(e + dsqrt(0.5d0)) * (e + dsqrt(0.5d0)))
     >     * dsqrt(0.125d0 / datan(1.0d0))
     >     + 0.5d0 * util_erfc(e + dsqrt(0.5d0))
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "psi.fh"

*     **** local variables ****
      integer ms,nx,ny,nz,n2ft3d
      real*8  scale,sumall

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call D3dB_n2ft3d(1,n2ft3d)
      !n2ft3d = 2*n2ft3d
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scale = lattice_omega()/dble(nx*ny*nz)

*     **** check total number of electrons ****
      en(2) = 0.0d0
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "psi.fh"


c*     ***** local variables ****
c      logical value
c      integer ms,shift1,shift2,ierr,i,j,indx
c      integer tmp1(2)


      call Dneall_m_diagonalize(0,dbl_mb(hml(1)),
     >                             dbl_mb(eig(1)),.false.)

c      value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
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
c      value = BA_pop_stack(tmp1(2))
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

#include "bafdecls.fh"
#include "psi.fh"


c*     ***** local variables ****
c      logical value
c      integer ms,shift1,shift2,ierr
c      integer tmp1(2)

      call Dneall_m_diagonalize(0,dbl_mb(hml(1)),dbl_mb(eig(1)),.true.)

c      value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
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

      
c      value = BA_pop_stack(tmp1(2))
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

#include "bafdecls.fh"
#include "psi.fh"

*     ***** local variables ****
      logical value
      integer k,n
      real*8  error,sum,size
      integer tmp1(2)

      value = BA_push_get(mt_dcpl,(npack1),'tmp1',tmp1(2),tmp1(1))
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
      
      value = BA_pop_stack(tmp1(2))
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
#include "bafdecls.fh"
#include "psi.fh"

*     ***** local variables ****
      logical value
      integer k,nx,ny,nz
      real*8  error,scale
      integer tmp1(2)

      real*8 e1,e2
      common /eenergy_tmp_common/ e1,e2

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

!$OMP DO 
      do k=1,(2*nfft3d)
         dbl_mb(tmp1(1)+k-1) = (dbl_mb(rho2(1)+k-1)
     >                         -dbl_mb(rho1(1)+k-1)) 
         dbl_mb(tmp1(1)+k-1) = dbl_mb(tmp1(1)+k-1) 
     >                      + (dbl_mb(rho2(1)+k-1+(ispin-1)*(2*nfft3d))
     >                        -dbl_mb(rho1(1)+k-1+(ispin-1)*(2*nfft3d)))
      end do
!$OMP END DO
      call D3dB_rr_dot(1,dbl_mb(tmp1(1)),dbl_mb(tmp1(1)),e1)
      error = e1*scale
*     error = dsqrt(error)


      value = BA_pop_stack(tmp1(2))
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

#include "bafdecls.fh"
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

      if (.not.BA_alloc_get(mt_int,nn,
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

#include "bafdecls.fh"
#include "errquit.fh"

      integer irreps(2)
      common / ab_irrep / irreps
     
      if (.not.BA_free_heap(irreps(2)))
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

#include "bafdecls.fh"

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
*     *                         *
*     *   psi1_crystal_dipole   *
*     *                         *
*     ***************************
*     
*     Uses - electron_crystal_dipole
*
      subroutine psi1_crystal_dipole(dipole)
      implicit none
      real*8 dipole(3)

#include "bafdecls.fh"
#include "psi.fh" 
      
      call Calculate_Resta_Dipole(.true.,ispin,ne,neq,npack1,nfft3d,
     >                            dcpl_mb(psi1(1)),dipole)

      return
      end


*     **********************************
*     *                                *
*     *   psi1_crystal_polarizability  *
*     *                                *
*     **********************************
*     
*     Uses - electron_crystal_polarizability
*       
      subroutine psi1_crystal_polarizability(alpha)
      implicit none
      real*8 alpha(3,3)
      
#include "bafdecls.fh"
#include "psi.fh" 
      
      call Calculate_Resta_Polarizability(.true.,
     >                   ispin,ne,neq,npack1,nfft3d,
     >                   dcpl_mb(psi1(1)),alpha)
     
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

      logical value

      value = BA_alloc_get(mt_dbl,2*ispin*nfft3d,'rho0',rho0(2),rho0(1))
      if (.not.value)
     >   call errquit('psi_cpmd_start',0,MA_ERR)

      call Parallel_shared_vector_copy(.true.,2*ispin*nfft3d,
     >                                 dbl_mb(rho1(1)),dbl_mb(rho0(1)))
      return
      end

*     ***************************
*     *                         *
*     *    psi_cpmd_end         *
*     *                         *
*     ***************************
      subroutine psi_cpmd_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

      if (.not.BA_free_heap(rho0(2)))
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

#include "bafdecls.fh"
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

      call Parallel_shared_vector_copy(.true.,2*ispin*nfft3d,
     >                          dbl_mb(rho1(1)),dbl_mb(rho2(1)))
      call DSCAL_OMP(2*ispin*nfft3d, 2.0d0,dbl_mb(rho2(1)),1)
      call DAXPY_OMP(2*ispin*nfft3d,-1.0d0,
     >               dbl_mb(rho0(1)),1,dbl_mb(rho2(1)),1)
      call Parallel_shared_vector_copy(.true.,2*ispin*nfft3d,
     >                          dbl_mb(rho1(1)),dbl_mb(rho0(1)))

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


#include "bafdecls.fh"
#include "btdb.fh"
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
      character*50 filename
      character*50 control_input_psi
      external     control_input_psi
      logical  wvfnc_expander,Dneall_m_allocate,band_reformat_c_wvfnc
      external wvfnc_expander,Dneall_m_allocate,band_reformat_c_wvfnc
      logical  psp_pawexist,control_print
      external psp_pawexist,control_print
      integer          control_fractional_smeartype
      double precision control_fractional_kT
      external         control_fractional_smeartype
      external         control_fractional_kT
      logical          control_ortho
      external         control_ortho


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
      pawexist = psp_pawexist()


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
      value = BA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'psi1',psi1(2),psi1(1))
      if (pawexist) then
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'spsi1',spsi1(2),spsi1(1))
      end if
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

      if (occupation_on) then
        value = value.and.
     >        BA_alloc_get(mt_dbl,(2*nemax),'occ1',occ1(2),occ1(1))
        value = value.and.
     >        BA_alloc_get(mt_dbl,(2*nemax),'occ2',occ2(2),occ2(1))
        smeartype = control_fractional_smeartype()
        smearkT   = control_fractional_kT()
      end if

      value = value.and.
     >        BA_alloc_get(mt_dbl,4*nfft3d,
     >                     'rho1_all',rho1_all(2),rho1_all(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,4*nfft3d,
     >                     'rho2_all',rho2_all(2),rho2_all(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)
c      call dcopy(4*nfft3d,0.0d0,0,dbl_mb(rho1_all(1)),1)
c      call dcopy(4*nfft3d,0.0d0,0,dbl_mb(rho2_all(1)),1)
      call Parallel_shared_vector_zero(.true.,4*nfft3d,
     >                    dbl_mb(rho1_all(1)))
      call Parallel_shared_vector_zero(.true.,4*nfft3d,
     >                    dbl_mb(rho2_all(1)))

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
        call Parallel_taskid(taskid)

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
        value = btdb_parallel(.true.)
        value = .true.

      end if

      end if

      call psi_read(ispin,ne,dcpl_mb(psi1(1)),
     >              smearoccupation,dbl_mb(occ1(1)))

      if (occupation_on) then
         call frac_occ_set(rtdb,ispin,ne,dbl_mb(occ1(1)))
      end if

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
      call psi_1ortho_check_fix()
c      if (pawexist) then
c         call psp_overlap_S(ispin,neq,dcpl_mb(psi1(1)),dcpl_mb(psi2(1)))
c         call Grsm_gg_trace(npack1,(neq(1)+neq(2)),
c     >                        dcpl_mb(psi1(1)),
c     >                        dcpl_mb(psi2(1)),
c     >                        sum2)
c      else
c         call Grsm_gg_trace(npack1,(neq(1)+neq(2)),
c     >                        dcpl_mb(psi1(1)),
c     >                        dcpl_mb(psi1(1)),
c     >                        sum2)
c      end if
c      call D1dB_SumAll(sum2)
c
c
c      sum1 = dble(ne(1) + ne(2))
c      if ((control_ortho()).and.(dabs(sum2-sum1).gt.1.0d-10)) then
c         call Parallel_taskid(taskid)
c         if (pawexist) then
c         call Dneall_f_Sortho(0,dcpl_mb(psi1(1)),
c     >                          dcpl_mb(psi2(1)),npack1)
c         call psp_overlap_S(ispin,neq,dcpl_mb(psi1(1)),
c     >                                dcpl_mb(psi2(1)))
c         call Grsm_gg_trace(npack1,(neq(1)+neq(2)),
c     >                        dcpl_mb(psi1(1)),
c     >                        dcpl_mb(psi2(1)),
c     >                        sum3)
c         else
cc         call Dneall_f_ortho(0,dcpl_mb(psi1(1)),npack1)
c         call Dneall_f_GramSchmidt(0,dcpl_mb(psi1(1)),npack1)
c         call Grsm_gg_trace(npack1,(neq(1)+neq(2)),
c     >                        dcpl_mb(psi1(1)),
c     >                        dcpl_mb(psi1(1)),
c     >                        sum3)
c         end if
c         call D1dB_SumAll(sum3)
c         if ((taskid.eq.MASTER).and.(control_print(print_medium)))
c     >    write(luout,*) 
c     >     "Warning - Gram-Schmidt being performed on psi:",
c     >               sum1,sum2,sum3,dabs(sum2-sum1)
c
c
c      end if

      psi_initialize = value
      return
      end

*     ***************************
*     *				*
*     *	 psi_1ortho_check_fix   *
*     *				*
*     ***************************
      subroutine psi_1ortho_check_fix()
      implicit none

#include "bafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "psi.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      real*8 sum2,sum1,sum3

*     **** external functions ****
      logical  control_ortho,control_print
      external control_ortho,control_print

*     **** Ortho Check ****
      if (pawexist) then
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
      if ((control_ortho()).and.(dabs(sum2-sum1).gt.1.0d-10)) then
         call Parallel_taskid(taskid)
         if (pawexist) then
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
         if ((taskid.eq.MASTER).and.(control_print(print_medium))) then
             write(luout,321) sum1,sum2,sum3,dabs(sum2-sum1)
         end if
 321     format(/" Warning - K.S. orbitals are not orthonormal. ",
     >      "Applying Gram-Schmidt orthonormalization."/,
     >       "         - exact norm=",E12.6," norm=",E12.6,
     >       " corrected norm=",E12.6,
     >       " (error=",E12.6,")"/)

      end if
      return
      end



*     ***************************
*     *				*
*     *	  psi_tmp_write  	*
*     *				*
*     ***************************
      subroutine psi_tmp_write()
      implicit none

#include "bafdecls.fh"
#include "psi.fh"

*     ***** write psi1 wavefunctions ****
      call psi_write(ispin,ne,dcpl_mb(psi1(1)),
     >               smearoccupation,dbl_mb(occ1(1)))

      return
      end


*     ************************************
*     *				         *
*     *	    psi_tmp_write_full_filename  *
*     *				         *
*     ************************************
      subroutine psi_tmp_write_full_filename(full_filename)
      implicit none
      character*(*) full_filename

#include "bafdecls.fh"
#include "psi.fh"

*     ***** write psi1 wavefunctions ****
      call psi_write_full_filename(full_filename,
     >                             ispin,ne,dcpl_mb(psi1(1)),
     >                             smearoccupation,dbl_mb(occ1(1)))

      return
      end

*     ************************************
*     *				         *
*     *	    psi_tmp_read_full_filename   *
*     *				         *
*     ************************************
      subroutine psi_tmp_read_full_filename(full_filename)
      implicit none
      character*(*) full_filename

#include "bafdecls.fh"
#include "psi.fh"

      call psi_read_full_filename(full_filename,
     >                            ispin,ne,dcpl_mb(psi1(1)),
     >                            smearoccupation,dbl_mb(occ1(1)))
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
#include "bafdecls.fh"
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
    
      value = BA_free_heap(eig(2))
      value = value.and.Dneall_m_free(hml)
      value = value.and.BA_free_heap(dng2(2))
      value = value.and.BA_free_heap(dng1(2))
      value = value.and.BA_free_heap(rho2(2))
      value = value.and.BA_free_heap(rho1(2))
      value = value.and.BA_free_heap(psi2(2))
      value = value.and.BA_free_heap(psi1(2))
      if (pawexist) then
          value = value.and.BA_free_heap(spsi1(2))
      end if
      value = value.and.BA_free_heap(rho2_all(2))
      value = value.and.BA_free_heap(rho1_all(2))
      if (occupation_on) then
         value = value.and.BA_free_heap(occ2(2))
         value = value.and.BA_free_heap(occ1(2))
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

#include "bafdecls.fh"
#include "psi.fh"

      psi_ne_excited = ne_excited(ms)
      return
      end

*     ***************************
*     *				*
*     *   psi_2epsi_gradients	*
*     *				*
*     ***************************
      logical function psi_2epsi_gradients(gn)
      implicit none
      integer gn

#include "bafdecls.fh"
#include "btdb.fh"
#include "psi.fh"
#include "errquit.fh"
      
*     **** local variables ****
      logical value
      integer nfac,nex(2),nemax,nexmax,ii,n
      real*8 tsum,tsum0
      
      nfac = 3
      if (gn.eq.2) then
         nfac = 9
      else if (gn.eq.3) then
         nfac = 19
      else 
         nfac = gn
      end if
      nex(1) = ne(1)*nfac
      nex(2) = ne(2)*nfac
      nemax   = ne(1) + ne(2)
      nexmax  = nex(1) + nex(2)

*     **** allocate memory ****
      value = BA_alloc_get(mt_dcpl,npack1*(nexmax),
     >         'psi1_excited',psi1_excited(2),psi1_excited(1))
      if (.not.value) 
     >  call errquit('psi_2epsi_gradients: out of heap memory',0,MA_ERR)

      do ii=1,nfac
         do n=1,nemax
            call Pack_cc_multiplegradients(1,ii,
     >             dcpl_mb(psi1(1)+(n-1)*npack1),
     >             dcpl_mb(psi1_excited(1)+((n-1)+(ii-1)*nemax)*npack1))
            call Pack_cc_dot(1,
     >             dcpl_mb(psi1_excited(1)+((n-1)+(ii-1)*nemax)*npack1),
     >             dcpl_mb(psi1_excited(1)+((n-1)+(ii-1)*nemax)*npack1),
     >             tsum)
            call Pack_cc_dot(1,
     >             dcpl_mb(psi1(1)+(n-1)*npack1),
     >             dcpl_mb(psi1(1)+(n-1)*npack1),
     >             tsum0)
            write(*,*) "ii,n,tsum=",ii,n,tsum,tsum0
         end do
      end do

      call epsi_write(ispin,nex,dcpl_mb(psi1_excited(1)))
      value = BA_free_heap(psi1_excited(2))
      if (.not.value) 
     >  call errquit('psi_2epsi_gradients:freeing heap',0,MA_ERR)

      return
      end


*     ***************************
*     *				*
*     *     epsi_initialize 	*
*     *				*
*     ***************************

      logical function epsi_initialize()
      implicit none 

#include "bafdecls.fh"
#include "btdb.fh"
#include "psi.fh"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical value,psi_nogrid
      integer nemax,ispin0

      integer hversion,hnfft(3),hispin,hne(2)
      real*8 hunita(3,3),sum1,sum2
      integer rtdb,ind,nn
      integer  control_rtdb,control_ngrid
      external control_rtdb,control_ngrid
      character*50 filename
      character*50 control_input_epsi
      external     control_input_epsi
      logical  wvfnc_expander
      external wvfnc_expander
      integer  control_symmetry
      external control_symmetry


*     ***** get ispin, and ne, and nfft3d ****
      call psi_get_ne_excited(ispin0,ne_excited)
      nemax  = ne_excited(1)  + ne_excited(2)
      nn = ne_excited(1)*ne_excited(1) + ne_excited(2)*ne_excited(2)

*     **** allocate memory ****
      value = BA_alloc_get(mt_dcpl,npack1*(nemax),
     >         'psi2_excited',psi2_excited(2),psi2_excited(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(nemax),
     >         'psi1_excited',psi1_excited(2),psi1_excited(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*nemax),'eig_excited',
     >                     eig_excited(2),eig_excited(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,nn,'hml_excited',
     >                     hml_excited(2),hml_excited(1))
      if (.not.value) 
     >   call errquit('epsi_initialize: out of heap memory',0,MA_ERR)

      !call dcopy(2*npack1*nemax,0.0d0,0,dcpl_mb(psi2_excited(1)),1)
      !call dcopy(2*npack1*nemax,0.0d0,0,dcpl_mb(psi1_excited(1)),1)
      !call dcopy(2*nemax,0.0d0,0,dbl_mb(eig_excited(1)),1)
      !call dcopy(nn,0.0d0,0,dbl_mb(hml_excited(1)),1)
      call Parallel_shared_vector_zero(.true.,2*npack1*nemax,
     >                    dcpl_mb(psi2_excited(1)))
      call Parallel_shared_vector_zero(.true.,2*npack1*nemax,
     >                    dcpl_mb(psi1_excited(1)))
      call Parallel_shared_vector_zero(.true.,2*nemax,
     >                    dbl_mb(eig_excited(1)))
      call Parallel_shared_vector_zero(.true.,nn,
     >                    dbl_mb(hml_excited(1)))


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


c*     **** Ortho Check ****
c      call Grsm_gg_trace(npack1,(ne_excited(1)+ne_excited(2)),
c     >                        dcpl_mb(psi1_excited(1)),
c     >                        dcpl_mb(psi1_excited(1)),
c     >                        sum2)
c      
c      sum1 = dble(ne_excited(1) + ne_excited(2))
c      if (dabs(sum2-sum1).gt.1.0d-10) then
c 
c         call Parallel_taskid(taskid)
c         call Grsm_g_MakeOrtho(npack1,ne_excited(1),
c     >                         dcpl_mb(psi1_excited(1)))
c         if (ispin.gt.1) then
c           call Grsm_g_MakeOrtho(npack1,ne_excited(2),
c     >                           dcpl_mb(psi1_excited(1)
c     >                                  +ne_excited(1)*npack1))
c         end if
c         call Grsm_gg_trace(npack1,(ne_excited(1)+ne_excited(2)),
c     >                        dcpl_mb(psi1_excited(1)),
c     >                        dcpl_mb(psi1_excited(1)),
c     >                        sum2)
c         if (taskid.eq.MASTER)
c     >    write(*,*) "Warning - Gram-Schmidt being performed on epsi:",
c     >               dabs(sum2-sum1)
c 
c      end if


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

#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     **** local variables ****
      logical value

*     ***** write psi1 wavefunctions ****
      if (writepsi)
     >  call epsi_write(ispin,ne_excited,dcpl_mb(psi1_excited(1)))
    
      value = BA_free_heap(hml_excited(2))
      value = value.and.BA_free_heap(eig_excited(2))
      value = value.and.BA_free_heap(psi2_excited(2))
      value = value.and.BA_free_heap(psi1_excited(2))
      if (.not. value) 
     >  call errquit('epsi_finalize: error freeing heap',0, MA_ERR)

      epsi_finalize = value
      return
      end




*     ***************************
*     *                         *
*     *     qpsi_initialize     *
*     *                         *
*     ***************************
      logical function qpsi_initialize()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"
#include "psi.fh"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical value,psi_nogrid
      integer nemax,ispin0

      integer hversion,hnfft(3),hispin,hne(2)
      real*8 hunita(3,3),sum1,sum2
      integer rtdb,ind,nn
      integer  control_rtdb,control_ngrid
      external control_rtdb,control_ngrid
      character*50 filename
      character*50 control_input_qpsi
      external     control_input_qpsi
      logical  wvfnc_expander
      external wvfnc_expander
      integer  control_symmetry
      external control_symmetry


      nemax  = neq(1)  + neq(2)
      nn = ne(1)*ne(1) + ne(2)*ne(2)

*     **** allocate memory ****
      value = BA_alloc_get(mt_dcpl,npack1*(nemax),
     >         'psi2_pertab',psi2_pertab(2),psi2_pertab(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(nemax),
     >         'psi1_pertab',psi1_pertab(2),psi1_pertab(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*nemax),'eig_excited',
     >                     eig_pertab(2),eig_pertab(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,nn,'hml_pertab',
     >                     hml_pertab(2),hml_pertab(1))
      if (.not.value)
     >   call errquit('qpsi_initialize: out of heap memory',0,MA_ERR)

      call Parallel_shared_vector_zero(.true.,2*npack1*nemax,
     >                    dcpl_mb(psi2_pertab(1)))
      call Parallel_shared_vector_zero(.true.,2*npack1*nemax,
     >                    dcpl_mb(psi1_pertab(1)))
      call Parallel_shared_vector_zero(.true.,2*nemax,
     >                    dbl_mb(eig_pertab(1)))
      call Parallel_shared_vector_zero(.true.,nn,
     >                    dbl_mb(hml_pertab(1)))


*     *****  read initial wavefunctions into psi1  ****
      rtdb = control_rtdb()
      if (.not.btdb_get(rtdb,'nwpw:psi_nogrid',
     >                  mt_log,1,psi_nogrid))
     >   psi_nogrid = .true.

      if (psi_nogrid) then


        filename =  control_input_qpsi()
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
      call qpsi_read(ispin,ne,dcpl_mb(psi1_pertab(1)))


*     **** force inversion symmetry ****
      if (control_symmetry().eq.1)  then
         call Parallel_taskid(taskid)
         if (taskid.eq.MASTER) then
         write(*,*)
         write(*,*)
     >   "Projecting reponse wavefunctions to have inversion symmetry"
         write(*,*)
         end if
         call psi_symm_project(ispin,neq,npack1,
     >                         dcpl_mb(psi1_pertab(1)))
      end if

      qpsi_initialize = value
      return
      end


*     ***************************
*     *                         *
*     *      qpsi_finalize      *
*     *                         *
*     ***************************

      logical function qpsi_finalize(writeqpsi)
      implicit none
      logical writeqpsi

#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     **** local variables ****
      logical value

*     ***** write psi1 wavefunctions ****
      if (writeqpsi)
     >  call qpsi_write(ispin,ne_excited,dcpl_mb(psi1_excited(1)))

      value = BA_free_heap(hml_pertab(2))
      value = value.and.BA_free_heap(eig_pertab(2))
      value = value.and.BA_free_heap(psi2_pertab(2))
      value = value.and.BA_free_heap(psi1_pertab(2))
      if (.not. value)
     >  call errquit('qpsi_finalize: error freeing heap',0, MA_ERR)

      qpsi_finalize = value
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

#include "bafdecls.fh"
#include "psi.fh"
#include "stdio.fh"


*     **** Lubin Water Analysis ****
      call pspw_Lubin_water_analysis(rtdb,ispin,ne,2*nfft3d,
     >                                 dbl_mb(rho1(1)))

*     **** Atom Analysis ****
      call pspw_atom_analysis(rtdb,ispin,2*nfft3d,dbl_mb(rho1(1)))

*     **** Mulliken Analysis ****
      call pspw_analysis(0,rtdb,ispin,ne,dcpl_mb(psi2(1)),
     >                                   dbl_mb(eig(1)))

      call pspw_gen_APC(ispin,ne,dcpl_mb(dng1(1)),.false.)
      call pspw_print_APC(luout)

      call pspw_gen_atom_Efield(rtdb,ispin,dbl_mb(rho1(1)),
     >                          dcpl_mb(dng1(1)))

      call pspw_gen_atom_Efield_grad(rtdb,ispin,ne,
     >                          dcpl_mb(psi1(1)),
     >                          dcpl_mb(dng1(1)))

      return
      end


*     ***********************
*     *                     *
*     *     psi_Born        *
*     *                     *
*     ***********************

      subroutine psi_Born()
      implicit none

#include "bafdecls.fh"
#include "psi.fh"
#include "stdio.fh"

      call pspw_gen_APC(ispin,ne,dcpl_mb(dng1(1)),.false.)
      call pspw_print_APC(luout)

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

#include "bafdecls.fh"
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

#include "btdb.fh"
#include "bafdecls.fh"
#include "psi.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer npoints,ii
      integer weight(2),nemax
      real*8 emin,emax,alpha
      character*255 filename

      nemax = ne(1)
      value = BA_push_get(mt_dbl,(nemax),'weight',weight(2),weight(1))
      if (.not. value) 
     >  call errquit('psi_dos:out of stack memory',0, MA_ERR)
      call dcopy(nemax,1.0d0,0,dbl_mb(weight(1)),1)


      if (.not.btdb_get(rtdb,'dos:alpha',mt_dbl,1,alpha)) then
        alpha = 0.05d0/27.2116d0
      end if

      if (.not.btdb_get(rtdb,'dos:npoints',mt_int,1,npoints)) then
        npoints = 500
      end if

      if (.not.btdb_get(rtdb,'dos:emin',mt_dbl,1,emin)) then
         emin = 99999.0d0
         do ii=1,ne(1)+ne(2)
           if (dbl_mb(eig(1)+ii-1).lt.emin) emin = dbl_mb(eig(1)+ii-1)
         end do
         emin = emin - 0.1d0
      end if

      if (.not.btdb_get(rtdb,'dos:emax',mt_dbl,1,emax)) then
         emax = -99999.0d0
         do ii=1,ne(1)+ne(2)
           if (dbl_mb(eig(1)+ii-1).gt.emax) emax = dbl_mb(eig(1)+ii-1)
         end do
         emax = emax + 0.1d0
      end if

*     **** generate DENSITY OF STATES *****
      if (ispin.eq.1) then
        filename = "smear_dos_both"
        call densityofstates(filename,.false.,
     >                     dbl_mb(eig(1)),dbl_mb(weight(1)),ne(1),
     >                     1.0d0,alpha,npoints,emin,emax)
        filename = "smear_fdos_both"
        call densityofstates(filename,.false.,
     >                     dbl_mb(eig(1)),dbl_mb(weight(1)),ne(1),
     >                     1.0d0,alpha,npoints,emin,emax)
      end if

      if (ispin.eq.2) then
        filename = "smear_dos_alpha"
        call densityofstates(filename,.false.,
     >                     dbl_mb(eig(1)),dbl_mb(weight(1)),ne(1),
     >                     1.0d0,alpha,npoints,emin,emax)
        filename = "smear_dos_beta"
        call densityofstates(filename,.false.,
     >               dbl_mb(eig(1)+ne(1)),dbl_mb(weight(1)),ne(2),
     >               -1.0d0,alpha,npoints,emin,emax)
        filename = "smear_fdos_alpha"
        call densityofstates(filename,.false.,
     >                     dbl_mb(eig(1)),dbl_mb(weight(1)),ne(1),
     >                     1.0d0,alpha,npoints,emin,emax)
        filename = "smear_fdos_beta"
        call densityofstates(filename,.false.,
     >               dbl_mb(eig(1)+ne(1)),dbl_mb(weight(1)),ne(2),
     >               -1.0d0,alpha,npoints,emin,emax)
      end if

 
      value = BA_pop_stack(weight(2))
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

#include "btdb.fh"
#include "bafdecls.fh"
#include "psi.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer weight(2),npoints,ii
      real*8 emin,emax,alpha
      character*255 filename

      value = BA_push_get(mt_dbl,(ne_excited(1)+ne_excited(2)),
     >                    'weight',weight(2),weight(1))
      if (.not. value)
     >  call errquit('epsi_dos:out of stack memory',0, MA_ERR)
      call dcopy((ne_excited(1)+ne_excited(2)),
     >          1.0d0,0,dbl_mb(weight(1)),1)

      if (.not.btdb_get(rtdb,'dos:alpha',mt_dbl,1,alpha)) then
        alpha = 0.05d0/27.2116d0
      end if

      if (.not.btdb_get(rtdb,'dos:npoints',mt_int,1,npoints)) then
        npoints = 500
      end if

      if (.not.btdb_get(rtdb,'dos:emin',mt_dbl,1,emin)) then
         emin = 99999.0d0
         do ii=1,ne_excited(1)+ne_excited(2)
           if (dbl_mb(eig_excited(1)+ii-1).lt.emin) 
     >       emin = dbl_mb(eig_excited(1)+ii-1)
         end do
         emin = emin - 0.1d0
      end if

      if (.not.btdb_get(rtdb,'dos:emax',mt_dbl,1,emax)) then
         emax = -99999.0d0
         do ii=1,ne_excited(1)+ne_excited(2)
           if (dbl_mb(eig_excited(1)+ii-1).gt.emax) 
     >       emax = dbl_mb(eig_excited(1)+ii-1)
         end do
         emax = emax + 0.1d0
      end if

*     **** generate DENSITY OF STATES *****
      if (ispin.eq.1) then
        filename = "smear_vdos_both"
        call densityofstates(filename,.false.,
     >                     dbl_mb(eig_excited(1)),dbl_mb(weight(1)),
     >                     ne_excited(1),
     >                     1.0d0,alpha,npoints,emin,emax)
        filename = "smear_dos_both"
        call densityofstates(filename,.true.,
     >                     dbl_mb(eig_excited(1)),dbl_mb(weight(1)),
     >                     ne_excited(1),
     >                     1.0d0,alpha,npoints,emin,emax)
      end if

      if (ispin.eq.2) then
        filename = "smear_vdos_alpha"
        call densityofstates(filename,.false.,
     >                     dbl_mb(eig_excited(1)),dbl_mb(weight(1)),
     >                     ne_excited(1),
     >                     1.0d0,alpha,npoints,emin,emax)
        filename = "smear_vdos_beta"
        call densityofstates(filename,.false.,
     >           dbl_mb(eig_excited(1)+ne_excited(1)),dbl_mb(weight(1)),
     >           ne_excited(2),
     >           -1.0d0,alpha,npoints,emin,emax)
        filename = "smear_dos_alpha"
        call densityofstates(filename,.true.,
     >                     dbl_mb(eig_excited(1)),dbl_mb(weight(1)),
     >                     ne_excited(1),
     >                     1.0d0,alpha,npoints,emin,emax)
        filename = "smear_dos_beta"
        call densityofstates(filename,.true.,
     >           dbl_mb(eig_excited(1)+ne_excited(1)),dbl_mb(weight(1)),
     >           ne_excited(2),
     >           -1.0d0,alpha,npoints,emin,emax)
      end if

      value = BA_pop_stack(weight(2))
      if (.not. value)
     >  call errquit('epsi_dos: error freeing stack',0, MA_ERR)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccc

*     ***********************
*     *                     *
*     *     psi_polariz     *
*     *                     *
*     ***********************

*** Early version of periodic dipole written by Patrick Nichols                                               ***
*** This routine has been replace by Calculate_Resta_Dipole, but keeping around as a reference implementation ***
*** this routine only works for cubic unit cells                                                              ***

      subroutine psi_polariz()
      implicit none

#include "psi.fh"
#include "bafdecls.fh"
#include "errquit.fh"
#include "util.fh"

*     **** local variables ****
      logical value
      integer psirx(2),tsirx(2),asize

      asize=(neq(1)+neq(2))*2*nfft3d
      value=          BA_push_get(mt_dbl,asize,"psir",psirx(2),psirx(1))
      value=value.and.BA_push_get(mt_dbl,asize,"tsir",tsirx(2),tsirx(1))
      if (.not.value) call errquit("psi_polariz stack empty",0,MA_ERR)

      call berry_phase_pol(ispin,ne,neq,npack1,nfft3d,
     >                     dcpl_mb(psi2(1)),
     >                     dbl_mb(psirx(1)),
     >                     dbl_mb(tsirx(1)))

      value=          BA_pop_stack(tsirx(2))
      value=value.and.BA_pop_stack(psirx(2))
      if (.not.value) call errquit("psi_polariz: pop stack",0,MA_ERR)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine epsi_generate_kb_vnm(vnm)
      implicit none
      real*8 vnm(*)

#include "bafdecls.fh"
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

      value = BA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
      value = value.and.
     >        BA_push_get(mt_dbl, 2*npack1,'tmp',tmp(2),tmp(1))
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
      !call dcopy(3*vshift,0.0d0,0,vnm,1)
      call Parallel_shared_vector_zero(.true.,3*vshift,vnm)

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
     >                          dbl_mb(tmp(1)),
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

      value =           BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(Gz(2))
      value = value.and.BA_pop_stack(Gy(2))
      value = value.and.BA_pop_stack(Gx(2))
      if (.not. value)
     >   call errquit('epsi_generate_kb_vnm:popping stack',1,MA_ERR)

      
      return
      end



*     *********************************
*     *                               *
*     *     psi_1pressure_stress      *
*     *                               *
*     *********************************

      subroutine psi_1pressure_stress(pressure,p1,p2,stress)
      implicit none
      real*8 pressure,p1,p2,stress(3,3)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     ***** rhoall common block ****
      integer rho1_all(2)
      integer rho2_all(2)
      common / rhoall_block / rho1_all,rho2_all
 
*     ***** external functions *****
      integer  electron_xcp_ptr
      external electron_xcp_ptr
      real*8   psi_1vnl,rho_1exc,rho_1pxc
      external psi_1vnl,rho_1exc,rho_1pxc

      call cgsd_pressure_stress(ispin,neq,
     >                          dcpl_mb(psi1(1)),
     >                          dbl_mb(rho1_all(1)),
     >                          dcpl_mb(dng1(1)),
     >                          dbl_mb(electron_xcp_ptr()),
     >                          psi_1vnl(),rho_1exc(),rho_1pxc(),
     >                          pressure,p1,p2,stress)


      return
      end



*     ***********************
*     *                     *
*     *    psi_MP2_energy   *
*     *                     *
*     ***********************
      subroutine psi_MP2_energy(rtdb)
      implicit none
      integer rtdb

#include "stdio.fh"
#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical oprint,value
      integer ms,a,b,r,s,n2ft3d,icount
      real*8 ea,eb,er,es
      real*8 e2,d2,tmp2
      integer vpsi_r(2),v1h(2),v2h(2)

*     **** allocate memory from heap ****
      n2ft3d = 2*nfft3d
      value = BA_alloc_get(mt_dcpl,nfft3d,'v1h',v1h(2),v1h(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d,'v2h',v2h(2),v2h(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,n2ft3d*(ne_excited(1)+ne_excited(2)),
     >                     'vpsi_r',vpsi_r(2),vpsi_r(1))
      if (.not. value)
     >  call errquit('psi_MP2_energy: error allocating heap memory',0,
     &       MA_ERR)


      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

      icount = 0
      do ms=1,ispin
         do a=1,ne(ms)-1
            ea = dbl_mb(eig(1)+a-1+(ms-1)*ne(1))
            do r=1,ne_excited(ms)-1
               er = dbl_mb(eig_excited(1)+r-1+(ms-1)*ne(1))

               do b=a+1,ne(ms)
                  eb = dbl_mb(eig(1)+b-1+(ms-1)*ne(1))
                  do s=r+1,ne_excited(ms)
                     es = dbl_mb(eig_excited(1)+s-1+(ms-1)*ne(1))
                     d2 = ea+eb-er-es
                     tmp2 = 1.0d0
                     e2 = e2 + tmp2*tmp2/d2
                     icount = icount + 1
                  if (oprint) 
     >            write(luout,*) "a,b,r,s, Esub=",a,b,r,s,(tmp2*tmp2/d2)
                  end do
               end do
            end do
         end do
      end do
      if (ispin.eq.1) e2=e2+e2

*     **** deallocate memory from heap ****
      value =     BA_free_heap(v1h(2))
     >       .and.BA_free_heap(v2h(2))
     >       .and.BA_free_heap(vpsi_r(2))
      if (.not. value)
     >  call errquit('psi_MP2_energy: error freeing heap memory',1,
     &       MA_ERR)



      if (oprint) then
         write(luout,*) "EMP2 = ", e2, " icount=",icount
      end if

      return
      end


*     ***********************
*     *                     *
*     *  psi_2q_integrals   *
*     *                     *
*     ***********************
      subroutine psi_2q_integrals(rtdb)
      implicit none
      integer rtdb
      
#include "stdio.fh"
#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"
      
*     **** local variables ****
      integer taskid,MASTER,taskid_j,pj
      parameter (MASTER=0)
      
      logical oprint,value
      integer icount,iicount,version
      integer i,j,k,l,ij,kl
      integer nall,nnall,n4all,nx,ny,nz
      integer psiij
      integer ipackm(2),jpackm(2)
      integer h1_integrals(2),h2_integrals(2)
      integer vij(2),dnij(2),orbi(2),orbj(2),vhij(2)
      real*8 e1,e2,e2h,scal1,scal2,dv

*     **** external functions ****
      integer  control_version
      external control_version
      real*8   lattice_omega
      external lattice_omega

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)
      oprint = (taskid.eq.MASTER)
      version = control_version()

      if (oprint) then
         write(luout,*) 
     >   "== Generating One-Electron and Two-Electron Integrals =="
         write(luout,*) 
      end if
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = scal1*lattice_omega()

      nall = ne(1) + ne_excited(1)
      nnall = (nall*(nall+1))/2
      n4all = (nnall*(nnall+1))/2

      value = BA_alloc_get(mt_dcpl,nfft3d,'vij',vij(2),vij(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d,'dnij',dnij(2),dnij(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d,'vhij',vhij(2),vhij(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d,'orbi',orbi(2),orbi(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d,'orbj',orbj(2),orbj(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,nnall,'h1_integrals',
     >                     h1_integrals(2),h1_integrals(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,n4all,'h2_integrals',
     >                     h2_integrals(2),h2_integrals(1))
      value = value.and.
     >        BA_alloc_get(mt_int,nnall,'ipackm',
     >                     ipackm(2),ipackm(1))
      value = value.and.
     >        BA_alloc_get(mt_int,nnall,'jpackm',
     >                     jpackm(2),jpackm(1))
      if (.not.value)
     >   call errquit('psi_2q_integrals:error allocating heap',0,MA_ERR)

      call Parallel_shared_vector_zero(.true.,nnall,
     >                                 dbl_mb(h1_integrals(1)))
      call Parallel_shared_vector_zero(.true.,n4all,
     >                                 dbl_mb(h2_integrals(1)))

      !**** generate 1e integrals - H1 only contains kinetic + ion-electron potentials ****
      call electron_gen_vall0()
      icount = 0
      do i=1,nall
         call psi_fetch_orbi_replicated(i,dcpl_mb(orbi(1)))
         call electron_get_H0psi_k_orb(dcpl_mb(orbi(1)),
     >                                 dcpl_mb(vij(1)))

         do j=i,nall
            call psi_fetch_orbj_indx(j,psiij,pj)
            if (pj.eq.taskid_j) then
               call Pack_cc_idot(1,dcpl_mb(psiij),dcpl_mb(vij(1)),e1)
               dbl_mb(h1_integrals(1)+icount) = e1
            end if
            int_mb(ipackm(1) + icount) = i
            int_mb(jpackm(1) + icount) = j
            icount = icount + 1
         end do
         call ga_sync()
      end do
      call Parallel_Vector_SumAll(nnall,dbl_mb(h1_integrals(1)))

      if (oprint) then
         write(luout,*) "begin_one_electron_integrals"
         icount = 0
         do i=1,nall
            do j=i,nall
               write(luout,'(I5,I5,F21.10)') 
     >         j,i,dbl_mb(h1_integrals(1)+icount)
               icount = icount + 1
            end do
         end do
         write(luout,*) "end_one_electron_integrals"
      end if

      !*** generate 2e integrals ***
      iicount = 0
      do ij=1,nnall
         i = int_mb(ipackm(1)+ij-1)
         j = int_mb(jpackm(1)+ij-1)
         call psi_fetch_orbi_replicated(i,dcpl_mb(orbi(1)))
         call Pack_c_unpack(1,   dcpl_mb(orbi(1)))
         call D3dB_cr_pfft3b(1,1,dcpl_mb(orbi(1)))
         call psi_fetch_orbi_replicated(j,dcpl_mb(orbj(1)))
         call Pack_c_unpack(1,   dcpl_mb(orbj(1)))
         call D3dB_cr_pfft3b(1,1,dcpl_mb(orbj(1)))

*        **** generate dnij for Vij  ****
         if (version.eq.4) then
            call D3dB_rr_Mul(1,dcpl_mb(orbi(1)),
     >                       dcpl_mb(orbj(1)), 
     >                       dcpl_mb(dnij(1)))
            call D3dB_r_SMul1(1,scal2,dcpl_mb(dnij(1)))
            call coulomb2_v(dcpl_mb(dnij(1)),dcpl_mb(vij(1)))
         else
            call D3dB_rr_Mul(1,dcpl_mb(orbi(1)),
     >                         dcpl_mb(orbj(1)), 
     >                         dcpl_mb(dnij(1)))
            call D3dB_r_SMul1(1,scal2,dcpl_mb(dnij(1)))
            call D3dB_r_SMul1(1,scal1,dcpl_mb(dnij(1)))
            call D3dB_rc_fft3f(1,dcpl_mb(dnij(1)))
            call Pack_c_pack(0,dcpl_mb(dnij(1)))

            call coulomb_v(dcpl_mb(dnij(1)),dcpl_mb(vhij(1)))
            call Pack_c_unpack(0,dcpl_mb(vhij(1)))
            call D3dB_cr_fft3b(1,dcpl_mb(vhij(1)))

            call coulomb_screened_v(dcpl_mb(dnij(1)),dcpl_mb(vij(1)))
            call Pack_c_unpack(0,dcpl_mb(vij(1)))
            call D3dB_cr_fft3b(1,dcpl_mb(vij(1)))
         end if

         do kl=ij,nnall
            k = int_mb(ipackm(1)+kl-1)
            l = int_mb(jpackm(1)+kl-1)
            call psi_fetch_orbi_replicated(k,dcpl_mb(orbi(1)))
            call Pack_c_unpack(1,   dcpl_mb(orbi(1)))
            call D3dB_cr_pfft3b(1,1,dcpl_mb(orbi(1)))
            call psi_fetch_orbi_replicated(l,dcpl_mb(orbj(1)))
            call Pack_c_unpack(1,   dcpl_mb(orbj(1)))
            call D3dB_cr_pfft3b(1,1,dcpl_mb(orbj(1)))

            call D3dB_rr_Mul(1,dcpl_mb(orbi(1)),
     >                         dcpl_mb(orbj(1)), 
     >                         dcpl_mb(dnij(1)))
            call D3dB_r_SMul1(1,scal2,dcpl_mb(dnij(1)))

            if (version .eq. 4) then
               call D3dB_rr_dot(1,dcpl_mb(dnij(1)),dcpl_mb(vij(1)),e2)
               !e2 = 0.5d0*e2*dv
               e2 = e2*dv
               dbl_mb(h2_integrals(1)+iicount) = e2

            else 
               if (i.eq.j .and. k.eq.l) then
                  call D3dB_rr_dot(1,dcpl_mb(dnij(1)),dcpl_mb(vhij(1)),
     >                             e2h)
                  e2h = e2h*dv
                  call D3dB_rr_dot(1,dcpl_mb(dnij(1)),dcpl_mb(vij(1)),
     >                             e2)
                  e2 = e2*dv
                  e2 = 2.0d0*e2h - e2
                  dbl_mb(h2_integrals(1)+iicount) = e2
               else if (i.eq.j .and. k.ne.l) then
                  call D3dB_rr_dot(1,dcpl_mb(dnij(1)),dcpl_mb(vhij(1)),
     >                             e2h)
                  e2h = e2h*dv
                  call D3dB_rr_dot(1,dcpl_mb(dnij(1)),dcpl_mb(vij(1)),
     >                             e2)
                  e2 = e2*dv
                  e2 = 2.0d0*e2h - e2
                  dbl_mb(h2_integrals(1)+iicount) = e2
               else if (i.ne.j .and. k.eq.l) then
                  call D3dB_rr_dot(1,dcpl_mb(dnij(1)),dcpl_mb(vhij(1)),
     >                             e2h)
                  e2h = e2h*dv
                  call D3dB_rr_dot(1,dcpl_mb(dnij(1)),dcpl_mb(vij(1)),
     >                             e2)
                  e2 = e2*dv
                  e2 = 2.0d0*e2h - e2
                  dbl_mb(h2_integrals(1)+iicount) = e2
               else
                  call D3dB_rr_dot(1,dcpl_mb(dnij(1)),dcpl_mb(vij(1)),
     >                             e2)
                  e2 = e2*dv
                  dbl_mb(h2_integrals(1)+iicount) = e2
               end if
            end if

            !write(*,*) "i,j,k,l,e2=",i,j,k,l,e2
            iicount = iicount + 1
         end do
      end do


      if (oprint) then
         write(luout,*)
         write(luout,*) "begin_two_electron_integrals"
         iicount = 0
         do ij=1,nnall
            i = int_mb(ipackm(1)+ij-1)
            j = int_mb(jpackm(1)+ij-1)
            do kl=ij,nnall
               k = int_mb(ipackm(1)+kl-1)
               l = int_mb(jpackm(1)+kl-1)
               write(luout,'(I5,I5,I5,I5,F20.10)') 
     >         j,i,l,k,dbl_mb(h2_integrals(1)+iicount)
               iicount = iicount + 1
            end do
         end do
         write(luout,*) "end_two_electron_integrals"
         write(luout,*)

      end if

*     **** deallocate memory from heap ****
      value = BA_free_heap(vij(2))
      value = value.and.BA_free_heap(dnij(2))
      value = value.and.BA_free_heap(vhij(2))
      value = value.and.BA_free_heap(orbi(2))
      value = value.and.BA_free_heap(orbj(2))
      value = value.and.BA_free_heap(h1_integrals(2))
      value = value.and.BA_free_heap(h2_integrals(2))
      value = value.and.BA_free_heap(ipackm(2))
      value = value.and.BA_free_heap(jpackm(2))
      if (.not. value)
     >   call errquit('psi_2q_integrals:error freeing heap',1,MA_ERR)

      return
      end



      subroutine psi_fetch_orbi_replicated(i,orbi)
      implicit none
      integer i
      complex*16 orbi(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     **** local variables ****
      integer pi,q,taskid_j,psii_ptr

      call Parallel2d_taskid_j(taskid_j)

      if (i.le.ne(1)) then
         call Dneall_ntoqp(i,q,pi)
         call Pack_c_Zero(1,orbi)
         if (pi.eq.taskid_j) then
            psii_ptr = psi1(1) + (q-1)*npack1
            call Pack_c_Copy(1,dcpl_mb(psii_ptr),orbi)
         end if
         call D1dB_Vector_SumAll(2*npack1,orbi)
      else
         psii_ptr = psi1_excited(1) + (i-1-ne(1))*npack1
         call Pack_c_Copy(1,dcpl_mb(psii_ptr),orbi)
      end if

      return
      end

      subroutine psi_fetch_orbj_indx(j,orbj_indx,pj)
      implicit none
      integer j
      integer orbj_indx
      integer pj

#include "bafdecls.fh"
#include "errquit.fh"
#include "psi.fh"

*     **** local variables ****
      integer q,np_j

      call Parallel2d_np_j(np_j)
      if (j.le.ne(1)) then
         call Dneall_ntoqp(j,q,pj)
         orbj_indx = psi1(1) + (q-1)*npack1
      else
         orbj_indx = psi1_excited(1) + (j-1-ne(1))*npack1
         pj = mod(pj+1,np_j)
      end if

      return
      end

*     ***********************************
*     *                                 *
*     *         psi_get_psi_ptr         *
*     *                                 *
*     ***********************************
      integer function psi_get_psi_ptr(i)
      implicit none
      integer i

#include "psi.fh"

      if (i.eq.1) then
         psi_get_psi_ptr = psi1(1)
      else
         psi_get_psi_ptr = psi2(1)
      end if

      return
      end




