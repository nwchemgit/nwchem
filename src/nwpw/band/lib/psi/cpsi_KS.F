* $Id$
*
*
* $Log: not supported by cvs2svn $
* Revision 1.16  2007/09/12 17:28:19  bylaska
* cpsi_data_ routines have been added...EJB
*
* Revision 1.15  2007/08/24 02:25:01  bylaska
* added density mixing to band...EJB
*
* Revision 1.14  2007/03/14 20:35:52  d3p708
* pjn
*    in cpsi_spin2 we should have cram_cc_zdot not cram_cc_dot. The overlap
* matrix is _complex_ for different complex wavefunctions.
*
* Revision 1.13  2007/02/23 01:18:11  bylaska
* ...EJB
*
* Revision 1.12  2007/01/24 22:32:51  d3p708
* pjn
*
* Revision 1.11  2007/01/23 23:03:58  d3p708
* PJN
*
* Revision 1.10  2007/01/22 21:28:59  d3p708
* pjn oops forget to check in proper file
*
* Revision 1.9  2007/01/22 02:23:19  d3p708
* PJN CHANGES TO SORT ROUTINE
*
* Revision 1.8  2007/01/22 02:17:38  d3p708
* pjn added missing shift
*
* Revision 1.7  2007/01/20 23:52:46  d3p708
* pjn removed extraneous output
*
* Revision 1.6  2007/01/19 22:34:53  d3p708
* pjn added spin-orbit
*
* Revision 1.5  2005/07/09 22:44:21  bylaska
* adding Louie FFT....EJB
* flag added for PAW xc and comp angular integration.
* ....EJB
*
* Revision 1.4  2005/02/01 02:50:06  bylaska
* Various updates...EJB
*
* Revision 1.3  2004/11/14 23:31:50  bylaska
* Band structure plotting code has been added to the BAND module.  The code the following task is used to run the code,
*
* task band structure
*
* The line traversed thru the Brillouin zone is defined by setting structure_zone_name to a previously
* defined zone_name,i.e.
*
* nwpw
*    ...
*    #define path in Brillouin zone
*    brillouin_zone
*      zone_name path1
*      kvector 0.00 0.00 0.0
*      kvector 0.10 0.00 0.0
*      kvector 0.20 0.00 0.0
*      kvector 0.30 0.00 0.0
*      kvector 0.40 0.00 0.0
*      kvector 0.50 0.00 0.0
*      kvector 0.40 0.00 0.0
*      kvector 0.30 0.00 0.0
*      kvector 0.20 0.00 0.0
*      kvector 0.10 0.00 0.0
*      kvector 0.00 0.00 0.0
*      ...
*    end
*    ...
*    structure_zone_name path1
*    ...
* end
*
*
* This code uses results from a task band energy calculation, i.e. the following sequence needs to be run for a reasonable band structure plot to be generated.
*
*      task band energy
*      task band structure
*
*
*   - This code has not been fully tested, but for silicon-carbide.
*   - The zone weights are ignored in a band structure plotting calculation
*   - kp extensions not yet added.
*   - automatic brillouin zone generation for band structure plotting  not yet added.
*
* ....EJB
*
* Revision 1.2  2004/11/13 02:22:05  bylaska
* Band-by-band minimizer implemented into BAND.
*   - works for Gamma point
*   - needs to be debugged for multiple k-points.
* ...EJB
*
* Revision 1.1  2003/12/02 19:16:47  bylaska
* HGH pseudpotential added.
* TM, Hamman, HGH, pspw_default, and paw_default pseudopotential libraries have been added.
* KS minimizer updates.
* ...EJB
*



*     ***************************
*     *                         *
*     *      cpsi_set_density   *
*     *                         *
*     ***************************

*    This routine sets the densities and potentials in psi and electron.
* This routine is needed for a band by band minimizer.
*
      subroutine cpsi_set_density(psi_number,rho)
      implicit none
      integer psi_number
      real*8 rho(*)


#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"


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

c      call dcopy(2*nfft3d,
c     >           rho, 1,
c     >           dbl_mb(rho_ptr),1)

      call Parallel_shared_vector_copy(.true.,2*nfft3d,
     >                                 rho,dbl_mb(rho_ptr))

      call c_electron_gen_dng_dnall(dbl_mb(rho_ptr),
     >                            dcpl_mb(dng_ptr),
     >                            dbl_mb(rho_all_ptr))
      call c_electron_gen_scf_potentials(dbl_mb(rho_ptr),
     >                            dcpl_mb(dng_ptr),
     >                            dbl_mb(rho_all_ptr))
      call c_electron_gen_vall()
      return
      end



*     ***************************
*     *                         *
*     *      cpsi_get_density   *
*     *                         *
*     ***************************

*    This routine gets the densities in psi.
* This routine is needed for a band by band minimizer.
*
      subroutine cpsi_get_density(psi_number,rho)
      implicit none
      integer psi_number
      real*8 rho(*)

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer rho_ptr

      if (psi_number.eq.1) then
        rho_ptr = rho1(1)
      else
        rho_ptr = rho2(1)
      end if

c      call dcopy(2*nfft3d,
c     >           dbl_mb(rho_ptr),1,
c     >           rho,1)
      call Parallel_shared_vector_copy(.true.,2*nfft3d,
     >                                 dbl_mb(rho_ptr),rho)
      return
      end


*     ***************************
*     *                         *
*     *    cpsi_write_density   *
*     *                         *
*     ***************************

*    This routine gets the densities in psi.
* This routine is needed for a band by band minimizer.
*
      subroutine cpsi_write_density(psi_number)
      implicit none
      integer psi_number

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer rho_ptr

      if (psi_number.eq.1) then
        rho_ptr = rho1(1)
      else
        rho_ptr = rho2(1)
      end if

      call crho_write(ispin,dbl_mb(rho_ptr))
      return
      end


*     ***************************
*     *                         *
*     *  cpsi_try_read_density  *
*     *                         *
*     ***************************

*    This routine reads the densities from disk to psi.
* This routine is needed for a band by band minimizer.
*
      logical function cpsi_try_read_density(psi_number)
      implicit none
      integer psi_number

#include "bafdecls.fh"
#include "cpsi_common.fh"

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

      if (rho_check_header(ispin,.true.)) then
         call crho_read(ispin,dbl_mb(rho_ptr))
         ok = .true.
      else
         ok = .false.
      end if

      cpsi_try_read_density = ok
      return
      end




*     **************************************
*     *                                    *
*     *     cpsi_gen_density_potentials    *
*     *                                    *
*     **************************************

*    This routine sets the densities and potentials in psi and electron.  
* This routine is needed for a band by band minimizer.
*
      subroutine cpsi_gen_density_potentials(psi_number)
      implicit none
      integer psi_number

#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "rhoall_common.fh"

*     **** local variables ****
      integer rho_ptr,dng_ptr,rho_all_ptr,psi_ptr

      if (psi_number.eq.1) then
        psi_ptr     = psi1_tag
        rho_ptr     = rho1(1)
        dng_ptr     = dng1(1)
        rho_all_ptr = rho1_all(1)
      else
        psi_ptr     = psi2_tag
        rho_ptr     = rho2(1)
        dng_ptr     = dng2(1)
        rho_all_ptr = rho2_all(1)
      end if

c      call c_electron_gen_psi_r(psi_ptr)
      call c_electron_gen_densities(dbl_mb(rho_ptr),
     >                              dcpl_mb(dng_ptr),
     >                              dbl_mb(rho_all_ptr))
      call c_electron_gen_scf_potentials(dbl_mb(rho_ptr),
     >                                 dcpl_mb(dng_ptr),
     >                                 dbl_mb(rho_all_ptr))
      call c_electron_gen_vall()
      return
      end




************************ KS orbital Part ************************

*     ***********************************
*     *                                 *
*     *        cpsi_KS_update00         *
*     *                                 *
*     ***********************************

*    This routine (approximately) diagonalizes the KS matrix.
*
      subroutine cpsi_KS_update00(psi_number,
     >                            maxit_orbs,maxit_orb,
     >                            precondition,
     >                            maxerror)
      implicit none
      integer psi_number
      integer maxit_orbs,maxit_orb
      logical precondition
      real*8  maxerror

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical done
      integer nb,i,j,neall,psi1_shift
      real*8 error,error_out,sum

*     **** external functions ****
      integer  control_ks_maxit_orb,control_ks_maxit_orbs
      integer  cpsi_data_get_ptr
      external control_ks_maxit_orb,control_ks_maxit_orbs
      external cpsi_data_get_ptr

      neall = neq(1)+neq(2)
      j = 0
 2    j = j+1
        error = 0.0d0
        do nb=1,nbrillq
        do i=1,neall
          psi1_shift = cpsi_data_get_ptr(psi1_tag,nb,i)

           !*** orthogonalize to lower orbitals  ****
           if (i.gt.1)
     >        call cpsi_project_out_f_orb1(nb,i,dbl_mb(psi1_shift))

           !*** normalize ****
           call Cram_cc_dot(nb,
     >                      dbl_mb(psi1_shift),
     >                      dbl_mb(psi1_shift),
     >                      sum)
           sum = 1.0d0/dsqrt(sum)
c           call Cram_c_SMul(nb,sum,
c     >                      dbl_mb(psi1_shift),
c     >                      dbl_mb(psi1_shift))
           call Cram_c_SMul1(nb,sum,dbl_mb(psi1_shift))


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

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical done
      integer nb,i,j,neall,maxit_orb,maxit_orbs,psi1_shift
      real*8 error,error_out,sum

*     **** external functions ****
      integer  control_ks_maxit_orb,control_ks_maxit_orbs
      integer  cpsi_data_get_ptr
      external control_ks_maxit_orb,control_ks_maxit_orbs
      external cpsi_data_get_ptr

      if (spin_orbit) then
        call cpsi_KS_update2com(psi_number,precondition,maxerror)
        return 
      end if

      neall = neq(1)+neq(2)
      maxit_orb  = control_ks_maxit_orb()   !*** should be read from rtdb ***
      maxit_orbs = control_ks_maxit_orbs()  !*** should be read from rtdb ***
      j = 0
 2    j = j+1
        error = 0.0d0
        do nb=1,nbrillq
        do i=1,neall
          psi1_shift = cpsi_data_get_ptr(psi1_tag,nb,i)

           !*** orthogonalize to lower orbitals  ****
           if (i.gt.1)
     >        call cpsi_project_out_f_orb1(nb,i,dbl_mb(psi1_shift))

           !*** normalize ****
           call Cram_cc_dot(nb,
     >                      dbl_mb(psi1_shift),
     >                      dbl_mb(psi1_shift),
     >                      sum)
           sum = 1.0d0/dsqrt(sum)
c           call Cram_c_SMul(nb,sum,
c     >                      dbl_mb(psi1_shift),
c     >                      dbl_mb(psi1_shift))
           call Cram_c_SMul1(nb,sum,dbl_mb(psi1_shift))


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
*     ***********************************
*     *                                 *
*     *      cpsi_KS_minimize           *
*     *                                 *
*     ***********************************

*    This routine (approximately) diagonalizes the KS matrix.
*
      subroutine cpsi_KS_minimize(psi_number,precondition,
     >                            maxerror,maxpsi_error)
      implicit none
      integer psi_number
      logical precondition
      real*8 maxerror
      real*8 maxpsi_error

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical stalled
      integer nb,i,j,neall,maxit_orb,maxit_orbs
      integer psi1_shift,eig_shift
      real*8 error_out,psi_error_out,e0,sum

*     **** external functions ****
      integer  control_ks_maxit_orb,control_ks_maxit_orbs
      integer  cpsi_data_get_ptr
      external control_ks_maxit_orb,control_ks_maxit_orbs
      external cpsi_data_get_ptr


      if (spin_orbit) then
        call  cpsi_KS_minimize2com(psi_number,precondition,
     >                            maxerror,maxpsi_error)
        return 
      end  if
       
      neall = neq(1)+neq(2)
      maxit_orb  = 120
      maxit_orbs = control_ks_maxit_orbs()  !*** should be read from rtdb ***
      stalled = .false.
  
        do nb=1,nbrillq
        do i=1,neall
          psi1_shift = cpsi_data_get_ptr(psi1_tag,nb,i)
          eig_shift  = cpsi_data_get_ptr(eig_tag,nb,i)

           if (stalled) then
              !*** intialize with steepest descent ****
              call c_sdminimize_noscf(0)

*             **** diagonalize current result ***
              call cpsi_1gen_hml()
              call cpsi_diagonalize_hml()
              call cpsi_1rotate2()
              call cpsi_2to1()
              stalled = .false.
           end if

           !*** orthogonalize to lower orbitals  ****
           call cpsi_project_out_f_orb1(nb,i,dbl_mb(psi1_shift))

           !*** normalize ****
           call Cram_cc_dot(nb,
     >              dbl_mb(psi1_shift),
     >              dbl_mb(psi1_shift),
     >              sum)
           sum = 1.0d0/dsqrt(sum)
c           call Cram_c_SMul(nb,sum,
c     >              dbl_mb(psi1_shift),
c     >              dbl_mb(psi1_shift))
           call Cram_c_SMul1(nb,sum,dbl_mb(psi1_shift))


            !**** minimize orbital ****
            j = 0
 2          call cpsi_KS_minimize_orb(psi_number,precondition,maxit_orb,
     >                       maxerror,maxpsi_error,0.001d0,
     >                       nb,i,error_out,psi_error_out,e0)
            j = j+1

            if (((error_out.gt.maxerror).or.
     >           (psi_error_out.gt.maxpsi_error))
     >          .and.(j.le.24)) go to 2
            if (j.gt.1) stalled = .true.

            dbl_mb(eig_shift) = e0

        end do
        end do

      call cpsi_sort_minimize()

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
      subroutine cpsi_project_out_f_orb1(nb,i,Horb)
      implicit none
      integer nb,i
      complex*16 Horb(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical ok
      integer n,ii,psi_ptr,x(2),nshift

*     **** external functions ***
      integer  cpsi_data_get_chnk,cpsi_data_nsize
      external cpsi_data_get_chnk,cpsi_data_nsize

      if (spin_orbit) then
        call cpsi_project_out_f_orb1_2com(nb,i,Horb)
        return 
      end  if


*     **** allocate stack memory ****
      ok = BA_push_get(mt_dcpl,ne(1),'x',x(2),x(1))
      if (.not.ok)
     > call errquit('cpsi_project_out_orb: out of stack memory',0,
     >       MA_ERR)

c      psi_ptr=psi1(1)+(nb-1)*(ne(1)+ne(2))*npack1
      psi_ptr = cpsi_data_get_chnk(psi1_tag,nb)
      nshift  = cpsi_data_nsize(psi1_tag)
      if (i.le.ne(1)) then
        ii = i-1
      else
        ii = i-ne(1)-1
        psi_ptr = psi_ptr + 2*ne(1)*npack1
      end if

      call Cram_cc_nzdot(nb,ii,
     >            dbl_mb(psi_ptr),
     >            Horb,
     >            dcpl_mb(x(1)))
      do n=1,(ii)
           call Cram_cc_zaxpy(nb,
     >               (-dcpl_mb(x(1)+n-1)),
     >               dbl_mb(psi_ptr),
     >               Horb)
           psi_ptr = psi_ptr + nshift
      end do

*     **** release stack memory ****
      ok = BA_pop_stack(x(2))
      if (.not. ok)
     > call errquit('cpsi_project_out_orb: poping stack memory',0,
     &       MA_ERR)

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

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical ok
      integer ii,n,psi_ptr,nsize
      integer x(2)
*     ***** external functions ****
      integer  cpsi_data_get_chnk,cpsi_data_nsize
      external cpsi_data_get_chnk,cpsi_data_nsize


      if (spin_orbit) then
        call cpsi_project_out_orb_2com(psi_number,nb,i,Horb)
        return 
      end if

*     **** allocate stack memory ****
      ok = BA_push_get(mt_dcpl,ne(1),'x',x(2),x(1))
      if (.not.ok)
     > call errquit('cpsi_project_out_orb: out of stack memory',0,
     &       MA_ERR)

      if (psi_number.eq.1) then
         psi_ptr = cpsi_data_get_chnk(psi1_tag,nb)
         nsize   = cpsi_data_nsize(psi1_tag)
      else
         psi_ptr = cpsi_data_get_chnk(psi2_tag,nb)
         nsize   = cpsi_data_nsize(psi2_tag)
      end if


      if (i.le.ne(1)) then
        ii = i
      else
        psi_ptr = psi_ptr + 2*ne(1)*npack1
        ii = i - ne(1)
      end if
      
      call Cram_cc_nzdot(nb,ii,
     >            dbl_mb(psi_ptr),
     >            Horb,
     >            dcpl_mb(x(1)))
      do n=1,(ii)
           call Cram_cc_zaxpy(nb,
     >               (-dcpl_mb(x(1)+n-1)),
     >               dbl_mb(psi_ptr),
     >               Horb)
         psi_ptr = psi_ptr + nsize
      end do


*     **** release stack memory ****
      ok = BA_pop_stack(x(2))
      if (.not. ok)
     > call errquit('cpsi_project_out_orb: poping stack memory',0,
     &       MA_ERR)

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

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer psi_ptr

      if (psi_number.eq.1) then
         psi_ptr=psi1_tag
      else
         psi_ptr=psi2_tag
      end if

      call c_electron_run_orb(nb,i,psi_ptr)
      call c_electron_get_gradient_orb(nb,i,Horb)

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

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical value,done,oneloop
      integer it
      real*8 e0,eold,error0,de0,lmbda_r0,lmbda_r1
      real*8 theta,sigma
      integer r1(2),t0(2),t(2),g(2)
      integer psi_ptr

*     **** external functions ****
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr

      if (spin_orbit) then
        call cpsi_KS_update_orb_2com(psi_number,
     >                             precondition,maxiteration,
     >                             maxerror,perror,nb,i,
     >                             error_out)
        return
      end if

      lmbda_r0 = 1.0d0
      if (psi_number.eq.1) then
         psi_ptr=cpsi_data_get_ptr(psi1_tag,nb,i)
      else
         psi_ptr=cpsi_data_get_ptr(psi2_tag,nb,i)
      end if


      value = BA_push_get(mt_dcpl,npack1,'t0',t0(2),t0(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1,'g',g(2),g(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1,'t',t(2),t(1))
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
ccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Check 
ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccc this cc_dot is okay! cccccccccccccccccc
*        *** calculate residual (steepest descent) direction for a single band ***
         call cpsi_get_gradient_orb(psi_number,nb,i,dcpl_mb(g(1)))
         call Cram_cc_dot(nb,dbl_mb(psi_ptr),
     >                   dcpl_mb(g(1)),
     >                    e0)
         e0 = -e0

         done = ((it.gt.maxiteration)
     >           .or.
     >           (dabs(e0-eold).lt.maxerror))

         if (done) go to 4

*        **** preconditioning ****
         if (precondition) then
           call cke_Precondition(nb,npack1,1,
     >                     dbl_mb(psi_ptr),
     >                     dcpl_mb(g(1)))

         end if

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
c!        call cpsi_project_out_orb(psi_number,i,dcpl_mb(t(1)))
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
c         call Cram_c_SMul(nb,de0,dcpl_mb(t(1)),dcpl_mb(t(1)))
         call Cram_c_SMul1(nb,de0,dcpl_mb(t(1)))



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
 4    value =           BA_pop_stack(t(2))
      value = value.and.BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(r1(2))
      value = value.and.BA_pop_stack(t0(2))
      if (.not. value) call errquit(
     >     'psi_KS_update_orb: popping stack memory',1, MA_ERR)

      error_out = dabs(e0-eold)
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

#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer orb(2),g(2),psi_ptr
      real*8 x,y,pi,e1

*     **** external functions ****
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr

      if (spin_orbit) then
        call cpsi_linesearch_update2_2com(psi_number,nb,i,theta,e0,
     >                                  de0,t,
     >                                  sigma,tau_t)
        return
      end if 

      if (psi_number.eq.1) then
         psi_ptr=cpsi_data_get_ptr(psi1_tag,nb,i)
      else
         psi_ptr=cpsi_data_get_ptr(psi2_tag,nb,i)
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


      call Cram_c_Copy(nb,dbl_mb(psi_ptr),dcpl_mb(orb(1)))

*     **** orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
      !theta = pi/300.0d0
      x = cos(theta)
      y = sin(theta)
      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)),
     >                  dbl_mb(psi_ptr))
      call Cram_cc_daxpy(nb,y,
     >                   t,
     >                   dbl_mb(psi_ptr))

*     *** determine theta ***
      call cpsi_get_gradient_orb(psi_number,nb,i,dcpl_mb(g(1)))

      call Cram_cc_dot(nb,dbl_mb(psi_ptr),
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
c      call Cram_c_SMul(nb,sigma,
c     >                  tau_t,
c     >                  tau_t)
      call Cram_c_SMul1(nb,sigma,tau_t)

*     **** orb2 = orb*cos(theta) + t*sin(theta) ****
      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)),
     >                  dbl_mb(psi_ptr))
      call Cram_cc_daxpy(nb,y,
     >                   t,
     >                   dbl_mb(psi_ptr))


*     **** release stack memory ****
      value =           BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(orb(2))
      if (.not. value) call errquit(
     >     'psi_linesearch_update: popping stack memory',1, MA_ERR)

      return
      end





*     ***********************************
*     *                                 *
*     *      cpsi_KS_minimize_orb      *
*     *                                 *
*     ***********************************

*    This routine performs a KS update on orbital i
*
      subroutine cpsi_KS_minimize_orb(psi_number,
     >                             precondition,maxiteration,
     >                             maxerror,maxpsi_error,perror,nb,i,
     >                             error_out,psi_error,e0)
      implicit none
      integer psi_number
      logical precondition
      integer maxiteration
      real*8  maxerror,maxpsi_error,perror
      integer nb,i
      real*8 error_out,psi_error
      real*8 e0

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical value,done,oneloop
      integer it
      real*8 eold,error0,de0,lmbda_r0,lmbda_r1
      real*8 theta,sigma
      integer r1(2),t0(2),t(2),g(2)
      integer psi_ptr

*     **** external functions ****
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr

     
      if (spin_orbit) then
        call cpsi_KS_minimize_orb_2com(psi_number,
     >                             precondition,maxiteration,
     >                             maxerror,maxpsi_error,perror,nb,i,
     >                             error_out,psi_error,e0)
        return
      end if 


      lmbda_r0 = 1.0d0
      if (psi_number.eq.1) then
         psi_ptr=cpsi_data_get_ptr(psi1_tag,nb,i)
      else
         psi_ptr=cpsi_data_get_ptr(psi2_tag,nb,i)
      end if


      value = BA_push_get(mt_dcpl,npack1,'t0',t0(2),t0(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1,'g',g(2),g(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1,'t',t(2),t(1))
      if (.not. value) call errquit(
     >     'cpsi_KS_minimize_orb: out of stack memory',0, MA_ERR)

      done = .false.
      error0 = 0.0d0
      e0 = 0.0d0
      theta = -3.14159d0/600.0d0
      it = 0
      psi_error = 10.0d0
 2    continue

         it = it + 1
         eold = e0

*        *** calculate residual (steepest descent) direction for a single band ***
         call cpsi_get_gradient_orb(psi_number,nb,i,dcpl_mb(g(1)))
         call Cram_cc_dot(nb,dbl_mb(psi_ptr),
     >                   dcpl_mb(g(1)),
     >                    e0)
         e0 = -e0

         done = ((it.gt.maxiteration)
     >           .or.
     >           ((dabs(e0-eold).lt.maxerror).and.
     >           (psi_error.lt.maxpsi_error)))

         if (done) go to 4

*        **** preconditioning ****
         if (precondition) then
           call cke_Precondition(nb,npack1,1,
     >                     dbl_mb(psi_ptr),
     >                     dcpl_mb(g(1)))

         end if

         call Cram_c_Copy(nb,dcpl_mb(g(1)),dcpl_mb(r1(1)))
         call Cram_cc_daxpy(nb,(e0),
     >                 dbl_mb(psi_ptr),
     >                 dcpl_mb(r1(1)))


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


*        *** normalize search direction, t ****
         call cpsi_project_out_orb(psi_number,nb,i,dcpl_mb(t(1)))
         call Cram_cc_dot(nb,dcpl_mb(t(1)),
     >                   dcpl_mb(t(1)),
     >                   sigma)
         sigma = dsqrt(sigma)
         de0 = 1.0d0/sigma
c         call Cram_c_SMul(nb,de0,dcpl_mb(t(1)),dcpl_mb(t(1)))
         call Cram_c_SMul1(nb,de0,dcpl_mb(t(1)))



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
         call cpsi_linesearch_minimize(psi_number,nb,i,
     >                              theta,e0,de0,
     >                              dcpl_mb(t(1)),psi_error)

      go to 2


*     **** release stack memory ****
 4    value =           BA_pop_stack(t(2))
      value = value.and.BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(r1(2))
      value = value.and.BA_pop_stack(t0(2))
      if (.not. value) call errquit(
     >     'psi_KS_minimize_orb: popping stack memory',1, MA_ERR)

c      write(*,*) "iterations=",it," eig=",e0," error=",error_out,
c     >           theta
      error_out = dabs(e0-eold)
      e0 = -e0
      return
      end



*     ***********************************
*     *                                 *
*     *      psi_linesearch_minimize    *
*     *                                 *
*     ***********************************

*    This routine performs a linesearch on orbital i, in the direction t.
* This routine is needed for a KS minimizer.
*  e0 = <orb|g>
*  de0 = 2*<t|g>
*
      subroutine cpsi_linesearch_minimize(psi_number,nb,i,theta,e0,de0,
     >                                    t,psi_error)
      implicit none
      integer psi_number
      integer nb,i
      real*8  theta
      real*8  e0,de0
      complex*16 t(*)     !search direction
      real*8     psi_error


#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer orb(2),g(2),psi_ptr
      real*8 x,y,pi,e1

*     **** external functions ****
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr

      if (spin_orbit) then
        call cpsi_linesearch_minimize2com(psi_number,nb,i,theta,e0,de0,
     >                                    t,psi_error)
        return
      end if 

      if (psi_number.eq.1) then
         psi_ptr=cpsi_data_get_ptr(psi1_tag,nb,i)
      else
         psi_ptr=cpsi_data_get_ptr(psi2_tag,nb,i)
      end if

      pi = 4.0d0*datan(1.0d0)

*     **** allocate stack memory ****
      value = BA_push_get(mt_dcpl,npack1,'orb',
     >                       orb(2),orb(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1,'g',
     >                       g(2),g(1))
      if (.not. value) call errquit(
     >     'psi_linesearch_minimize: out of stack memory',0, MA_ERR)


      call Cram_c_Copy(nb,dbl_mb(psi_ptr),
     >                    dcpl_mb(orb(1)))

*     **** orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
      !theta = pi/300.0d0
      x = cos(theta)
      y = sin(theta)
      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)),
     >                  dbl_mb(psi_ptr))
      call Cram_cc_daxpy(nb,y,
     >                   t,
     >                   dbl_mb(psi_ptr))

*     *** determine theta ***
      call cpsi_get_gradient_orb(psi_number,nb,i,dcpl_mb(g(1)))

      call Cram_cc_dot(nb,dbl_mb(psi_ptr),
     >                   dcpl_mb(g(1)),
     >                   e1)
      e1 = -e1
      x = (e0 - e1 + 0.5d0*de0*sin(2*theta))
     >    /(1.0d0-cos(2*theta))
      theta = 0.5d0*datan(0.5d0*de0/x)


*     **** orb2 = orb*cos(theta) + t*sin(theta) ****
      x = cos(theta)
      y = sin(theta)
      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)),
     >                  dbl_mb(psi_ptr))
      call Cram_cc_daxpy(nb,y,
     >                   t,
     >                   dbl_mb(psi_ptr))

*     **** calculated psi_error = <(orb-psi)|(orb-psi)>  ****
      call Cram_cc_daxpy(nb,(-1.0d0),dbl_mb(psi_ptr),dcpl_mb(orb(1)))
      call Cram_cc_dot(nb,dcpl_mb(orb(1)),dcpl_mb(orb(1)),psi_error)


*     **** release stack memory ****
      value =           BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(orb(2))
      if (.not. value) call errquit(
     >     'psi_linesearch_minimize: popping stack memory',1, MA_ERR)

      return
      end


      subroutine cpsi_sort_minimize()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"


      logical value
      integer i,j,ii,jj,ms,nb,nshift
      integer r1(2),psi1_shift,eig_shift
      real*8  ei,ej

*      ***** external functions ****
       integer  cpsi_data_get_chnk
       external cpsi_data_get_chnk

      if (spin_orbit) then
        call cpsi_sort_minimize2com()
        return
      end if

      value = BA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      if (.not. value) call errquit(
     >     'psi_sort_sort_minimize: out of stack memory',0,MA_ERR)

      nshift = 2*npack1
      do nb=1,nbrillq
        psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
        eig_shift  = cpsi_data_get_chnk(eig_tag,nb)
        do ms=1,ispin

        !*** Bubble sort ***
        do ii=1,ne(ms)
         do jj=ii+1,ne(ms)
           i = ii + (ms-1)*ne(1)
           j = jj + (ms-1)*ne(1)
           ei = dbl_mb(eig_shift+i-1)
           ej = dbl_mb(eig_shift+j-1)

           !*** swap ***
           if (ej.lt.ei) then
             dbl_mb(eig_shift+i-1) = ej
             dbl_mb(eig_shift+j-1) = ei
             call Cram_c_Copy(nb,dbl_mb(psi1_shift+(i-1)*nshift),
     >                           dcpl_mb(r1(1)))
             call Cram_c_Copy(nb,dbl_mb(psi1_shift+(j-1)*nshift),
     >                           dbl_mb(psi1_shift+(i-1)*nshift))
             call Cram_c_Copy(nb,dcpl_mb(r1(1)),
     >                           dbl_mb(psi1_shift+(j-1)*nshift))
           end if

         end do
        end do

      end do
      end do

      value = BA_pop_stack(r1(2))
      if (.not. value) call errquit(
     >     'cpsi_sort_minimize: popping stack memory',1, MA_ERR)
      return
      end
************************ KS orbital Part ************************
c  routines for KS MINIMIZATION WITH 2 COMPONENT WAVEFUNCTIONS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

*     ***********************************
*     *                                 *
*     *      cpsi_KS_update2com         *
*     *                                 *
*     ***********************************

*    This routine (approximately) diagonalizes the KS matrix.
*
      subroutine cpsi_KS_update2com(psi_number,precondition,maxerror)
      implicit none
      integer psi_number
      logical precondition
      real*8 maxerror

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical done
      integer nb,i,j,neall,maxit_orb,maxit_orbs,shifts,nshift,psi1_shift
      real*8 error,error_out,sum,sum1,sum2

*     **** external functions ****
      integer  control_ks_maxit_orb,control_ks_maxit_orbs
      integer  cpsi_data_get_chnk
      external control_ks_maxit_orb,control_ks_maxit_orbs
      external cpsi_data_get_chnk

      nshift=2*npack1
      shifts=nshift*ne(1)
      neall = neq(1)
      maxit_orb  = control_ks_maxit_orb()   !*** should be read from rtdb ***
      maxit_orbs = control_ks_maxit_orbs()  !*** should be read from rtdb ***
      j = 0
 2    j = j+1
        error = 0.0d0
        do nb=1,nbrillq
          psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
          do i=1,neall

           !*** orthogonalize to lower orbitals  ****
           call cpsi_project_out_f_orb1_2com(nb,i,
     >             dbl_mb(psi1_shift+(i-1)*nshift))

           !*** normalize ****
           call Cram_cc_dot(nb,
     >              dbl_mb(psi1_shift+(i-1)*nshift),
     >              dbl_mb(psi1_shift+(i-1)*nshift),
     >              sum1)
           call Cram_cc_dot(nb,
     >              dbl_mb(psi1_shift+(i-1)*nshift + shifts),
     >              dbl_mb(psi1_shift+(i-1)*nshift + shifts),
     >              sum2)
           sum = sum1 + sum2
           sum = 1.0d0/dsqrt(sum)
c           call Cram_c_SMul(nb,sum,
c     >              dbl_mb(psi1_shift+(i-1)*nshift),
c     >              dbl_mb(psi1_shift+(i-1)*nshift))
           call Cram_c_SMul1(nb,sum,
     >              dbl_mb(psi1_shift+(i-1)*nshift))
c           call Cram_c_SMul(nb,sum,
c     >              dbl_mb(psi1_shift+(i-1)*nshift + shifts),
c     >              dbl_mb(psi1_shift+(i-1)*nshift + shifts))
           call Cram_c_SMul1(nb,sum,
     >              dbl_mb(psi1_shift+(i-1)*nshift + shifts))


            call cpsi_KS_update_orb_2com(psi_number,precondition,
     >                           maxit_orb,
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





*     ***********************************
*     *                                 *
*     *      cpsi_KS_minimize2com       *
*     *                                 *
*     ***********************************

*    This routine (approximately) diagonalizes the KS matrix.
*
      subroutine cpsi_KS_minimize2com(psi_number,precondition,
     >                            maxerror,maxpsi_error)
      implicit none
      integer psi_number
      logical precondition
      real*8 maxerror
      real*8 maxpsi_error

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer nb,i,j,ne1,maxit_orb,maxit_orbs,shifts
      integer psi_ptr,eig_ptr
      real*8 error_out,psi_error_out,e0,sum,sum1,sum2

*     **** external functions ****
      integer  control_ks_maxit_orb,control_ks_maxit_orbs
      integer  cpsi_data_get_ptr
      external control_ks_maxit_orb,control_ks_maxit_orbs
      external cpsi_data_get_ptr

      ne1 = ne(1)
      maxit_orb  = 120
      maxit_orbs = control_ks_maxit_orbs()  !*** should be read from rtdb ***

      shifts=2*npack1*ne(1)   

c      do nb=1,nbrillioun
c           do i=1,ne(1)
c           do korb=1,ne(1)
c              psi_ptr=psi1(1)+(i-1)*npack1+(nb-1)*shifts*2
c              kpsi=psi1(1)+(korb-1)*npack1+(nb-1)*shifts*2
c              call Cram_cc_zdot(nb,
c     >              dcpl_mb(psi_ptr),
c     >              dcpl_mb(kpsi),zsum1)
c              call Cram_cc_zdot(nb,
c     >              dcpl_mb(psi_ptr+shifts),
c     >              dcpl_mb(kpsi+shifts),zsum2)
c              zsum1=zsum1+zsum2
c              write(*,*)"AAOvlp(",i,",",korb,")= ",zsum1
c           end do
c           end do
c      end do
      do nb=1,nbrillq
      do i=1,ne1
         psi_ptr=cpsi_data_get_ptr(psi1_tag,nb,i)
         eig_ptr=cpsi_data_get_ptr(eig_tag,nb,i)


            !*** orthogonalize to lower orbitals  ****
            call cpsi_project_out_f_orb1_2com(
     >             nb,i,
     >             dbl_mb(psi_ptr))

           !*** normalize ****
           call Cram_cc_dot(nb,
     >              dbl_mb(psi_ptr),
     >              dbl_mb(psi_ptr),
     >              sum1)
           call Cram_cc_dot(nb,
     >              dbl_mb(psi_ptr+shifts),
     >              dbl_mb(psi_ptr+shifts),
     >              sum2)
           sum = sum1 + sum2
           sum = 1.0d0/dsqrt(sum)

c           call Cram_c_SMul(nb,sum,
c     >              dbl_mb(psi_ptr),
c     >              dbl_mb(psi_ptr))
c           call Cram_c_SMul(nb,sum,
c     >              dbl_mb(psi_ptr+shifts),
c     >              dbl_mb(psi_ptr+shifts))
           call Cram_c_SMul1(nb,sum,dbl_mb(psi_ptr))
           call Cram_c_SMul1(nb,sum,dbl_mb(psi_ptr+shifts))


            !**** minimize orbital ****
            j = 0
 2          call cpsi_KS_minimize_orb_2com(psi_number,precondition,
     >                       maxit_orb,
     >                       maxerror,maxpsi_error,0.001d0,
     >                       nb,i,error_out,psi_error_out,e0)

            j = j+1
            call flush(6)
            if (((error_out.gt.maxerror).or.
     >           (psi_error_out.gt.maxpsi_error))
     >          .and.(j.le.24)) go to 2

            dbl_mb(eig_ptr) = e0


c           write(*,*)"THIS IS I IN KS_MINIM",i," ",e0
           call flush(6)
          end do
        end do
        call cpsi_sort_minimize2com()

      return
      end


*     *******************************************
*     *                                         *
*     *          cpsi_project_out_f_orb1_2com   *
*     *                                         *
*     *******************************************
*
*    This routine projects out non-orthogonal components of Horb.
* This routine is needed for a KS minimizer.
*
      subroutine cpsi_project_out_f_orb1_2com(nb,i,Horb)
      implicit none
      integer nb,i
      complex*16 Horb(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical ok
      integer n,ii,psi_ptr,x(2),shifts,shift
      complex*16 sum1,sum2

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

*     **** allocate stack memory ****
      ok = BA_push_get(mt_dcpl,2*ne(1),'x',x(2),x(1))
      if (.not.ok)
     > call errquit('cpsi_project_out_orb: out of stack memory',0,
     &       MA_ERR)

      shifts=npack1*ne(1)
      psi_ptr=cpsi_data_get_chnk(psi1_tag,nb)

      ii=i-1
      shift = psi_ptr
      do n=1,ii
        call Cram_cc_zdot(nb,dbl_mb(shift),Horb,sum1)
        call Cram_cc_zdot(nb,dbl_mb(shift+2*shifts),Horb(1+shifts),sum2)
        dcpl_mb(x(1)+n-1)=(sum1+sum2)
        shift = shift + 2*npack1
      end do

      shift = psi_ptr
      do n=1,ii
         sum1= (- dcpl_mb(x(1)+n-1))
         call Cram_cc_zaxpy(nb,
     >             sum1,
     >             dbl_mb(shift),
     >             Horb)
         call Cram_cc_zaxpy(nb,
     >             sum1,
     >             dbl_mb(shift+2*shifts),
     >             Horb(shifts+1))
         shift = shift + 2*npack1
      end do

*     **** release stack memory ****
      ok = BA_pop_stack(x(2))
      if (.not. ok)
     > call errquit('cpsi_project_out_orb: poping stack memory',0,
     &       MA_ERR)

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
      subroutine cpsi_project_out_orb_2com(psi_number,nb,i,Horb)
      implicit none
      integer psi_number
      integer nb,i
      complex*16 Horb(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical ok
      integer ii,n,psi_ptr,shifts,shift
      integer x(2)
      complex*16 zsum1,zsum2
*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

*     **** allocate stack memory ****
      ok = BA_push_get(mt_dcpl,ne(1),'x',x(2),x(1))
      if (.not.ok)
     > call errquit('cpsi_project_out_orb: out of stack memory',0,
     &       MA_ERR)

      shifts=2*npack1*ne(1)
      if (psi_number.eq.1) then
        psi_ptr=cpsi_data_get_chnk(psi1_tag,nb)
      else
        psi_ptr=cpsi_data_get_chnk(psi2_tag,nb)
      end if
      ii=i
      shift = psi_ptr
      do n=1,ii
         call Cram_cc_zdot(nb,dbl_mb(shift),Horb,zsum1)
         call Cram_cc_zdot(nb,dbl_mb(shift+shifts),Horb(1+npack1),zsum2)
         dcpl_mb(x(1)+n-1)=zsum1+zsum2
         shift = shift + 2*npack1
      end do

      shift = psi_ptr
      do n=1,ii
         zsum1 = -dcpl_mb(x(1)+n-1)
         call Cram_cc_zaxpy(nb,
     >               zsum1,
     >               dbl_mb(shift),
     >               Horb)
         call Cram_cc_zaxpy(nb,
     >               zsum1,
     >               dbl_mb(shift+shifts),
     >               Horb(1+npack1))
         shift = shift + 2*npack1
      end do


*     **** release stack memory ****
      ok = BA_pop_stack(x(2))
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
      subroutine cpsi_KS_update_orb_2com(psi_number,
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical value,done,oneloop
      integer it,shifts
      real*8 e0,eold,error0,de0,lmbda_r0,lmbda_r1
      real*8 theta,sigma,tmp
      integer r1(2),t0(2),t(2),g(2)
      integer psi_ptr

*     **** external functions ****
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr


      lmbda_r0 = 1.0d0
      shifts=2*npack1*ne(1)
      if (psi_number.eq.1) then
         psi_ptr=cpsi_data_get_ptr(psi1_tag,nb,i)
      else
         psi_ptr=cpsi_data_get_ptr(psi2_tag,nb,i)
      end if


      value = BA_push_get(mt_dcpl,2*npack1,'t0',t0(2),t0(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,2*npack1,'r1',r1(2),r1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,2*npack1,'g',g(2),g(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,2*npack1,'t',t(2),t(1))
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
         call Cram_cc_dot(nb,dbl_mb(psi_ptr),
     >                   dcpl_mb(g(1)),
     >                    e0)
         call Cram_cc_dot(nb,dbl_mb(psi_ptr+shifts),
     >                   dcpl_mb(g(1)+npack1),
     >                    tmp)
         e0 = -(e0+tmp)

         done = ((it.gt.maxiteration)
     >           .or.
     >           (dabs(e0-eold).lt.maxerror))

         if (done) go to 4

*        **** preconditioning ****
         if (precondition) then
           call cke_Precondition(nb,npack1,1,
     >                     dbl_mb(psi_ptr),
     >                     dcpl_mb(g(1)))

         end if

         call Cram_c_Copy(nb,dcpl_mb(g(1)),dcpl_mb(r1(1)))
         call Cram_c_Copy(nb,dcpl_mb(g(1)+npack1),
     >       dcpl_mb(r1(1)+npack1))
         call cpsi_project_out_orb_2com(psi_number,nb,i,
     >       dcpl_mb(r1(1)))



*        *** determine conjuagate direction ***
         call Cram_cc_dot(nb,dcpl_mb(r1(1)),
     >                   dcpl_mb(r1(1)),
     >                   lmbda_r1)
         call Cram_cc_dot(nb,dcpl_mb(r1(1)+npack1),
     >                   dcpl_mb(r1(1)+npack1),
     >                   tmp)
         lmbda_r1=lmbda_r1+tmp
         call Cram_c_Copy(nb,dcpl_mb(r1(1)),dcpl_mb(t(1)))
         call Cram_c_Copy(nb,dcpl_mb(r1(1)+npack1),
     >       dcpl_mb(t(1)+npack1))

         if (it.gt.1) then
         call Cram_cc_daxpy(nb,(lmbda_r1/lmbda_r0),
     >                   dcpl_mb(t0(1)),
     >                   dcpl_mb(t(1)))
         call Cram_cc_daxpy(nb,(lmbda_r1/lmbda_r0),
     >                   dcpl_mb(t0(1)+npack1),
     >                   dcpl_mb(t(1)+npack1))
         end if
         lmbda_r0 = lmbda_r1
         oneloop = .true.
 3       call Cram_c_Copy(nb,dcpl_mb(t(1)),dcpl_mb(t0(1)))
         call Cram_c_Copy(nb,dcpl_mb(t(1)+npack1),
     >      dcpl_mb(t0(1)+npack1))

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
         call Cram_cc_dot(nb,dcpl_mb(t(1)+npack1),
     >                   dcpl_mb(t(1)+npack1),
     >                   tmp)
         sigma=sigma+tmp
         sigma = dsqrt(sigma)
         de0 = 1.0d0/sigma
c         call Cram_c_SMul(nb,de0,dcpl_mb(t(1)),dcpl_mb(t(1)))
c         call Cram_c_SMul(nb,de0,dcpl_mb(t(1)+npack1),
c     >        dcpl_mb(t(1)+npack1))
         call Cram_c_SMul1(nb,de0,dcpl_mb(t(1)))
         call Cram_c_SMul1(nb,de0,dcpl_mb(t(1)+npack1))



*        **** compute de0 = <t|g> ****
         call Cram_cc_dot(nb,dcpl_mb(t(1)),
     >                   dcpl_mb(g(1)),
     >                   de0)
         call Cram_cc_dot(nb,dcpl_mb(t(1)+npack1),
     >                   dcpl_mb(g(1)+npack1),
     >                   tmp)
         de0=de0+tmp

*        *** bad direction ***
         if ((de0.lt.0.0d0).and.oneloop) then
           call Cram_c_Copy(nb,dcpl_mb(g(1)),dcpl_mb(t(1)))
           call Cram_c_Copy(nb,dcpl_mb(g(1)+npack1),
     >         dcpl_mb(t(1)+npack1))
           oneloop = .false.
           go to 3
         end if

         de0 = -2.0d0*de0
         call cpsi_linesearch_update2_2com(psi_number,nb,i,
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
*     *                                 *
*     *      psi_linesearch_update2     *
*     *                                 *
*     ***********************************

*    This routine performs a linesearch on orbital i, in the direction t.
* This routine is needed for a KS minimizer.
*  e0 = <orb|g>
*  de0 = 2*<t|g>
*
      subroutine cpsi_linesearch_update2_2com(psi_number,nb,i,theta,
     >                                  e0,de0,t,
     >                                  sigma,tau_t)
      implicit none
      integer psi_number
      integer nb,i
      real*8  theta
      real*8  e0,de0
      complex*16 t(*)     !search direction

      real*8     sigma
      complex*16 tau_t(*) !parallel transported search direction

#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer orb(2),g(2),psi_ptr,shifts
      real*8 x,y,pi,e1
      real*8 tmp
*     **** external functions ****
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr

      shifts=2*ne(1)*npack1
      if (psi_number.eq.1) then
         psi_ptr=cpsi_data_get_ptr(psi1_tag,nb,i)
      else
         psi_ptr=cpsi_data_get_ptr(psi2_tag,nb,i)
      end if

      pi = 4.0d0*datan(1.0d0)

*     **** allocate stack memory ****
      value = BA_push_get(mt_dcpl,npack1*2,'orb',
     >                       orb(2),orb(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1*2,'g',
     >                       g(2),g(1))
      if (.not. value) call errquit(
     >     'psi_linesearch_update2com: out of stack memory',0, MA_ERR)

      call Cram_c_Copy(nb,dbl_mb(psi_ptr),dcpl_mb(orb(1)))
      call Cram_c_Copy(nb,dbl_mb(psi_ptr+shifts),dcpl_mb(orb(1)+npack1))

*     **** orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
      !theta = pi/300.0d0
      x = dcos(theta)
      y = dsin(theta)

      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)),
     >                  dbl_mb(psi_ptr))
      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)+npack1),
     >                  dbl_mb(psi_ptr+shifts))
      call Cram_cc_daxpy(nb,y,
     >                   t,
     >                   dcpl_mb(psi_ptr))
      call Cram_cc_daxpy(nb,y,
     >                   t(1+npack1),
     >                   dbl_mb(psi_ptr+shifts))

*     *** determine theta ***
      call cpsi_get_gradient_orb(psi_number,nb,i,dcpl_mb(g(1)))

      call Cram_cc_dot(nb,dbl_mb(psi_ptr),
     >                   dcpl_mb(g(1)),
     >                   e1)
      call Cram_cc_dot(nb,dbl_mb(psi_ptr+shifts),
     >                   dcpl_mb(g(1)+npack1),
     >                   tmp)
      e1 = -(e1+tmp)
      x = (e0 - e1 + 0.5d0*de0*dsin(2*theta))
     >    /(1.0d0-dcos(2*theta))
      theta = 0.5d0*datan(0.5d0*de0/x)

      x = dcos(theta)
      y = dsin(theta)

*     **** tau_t = (-orb*sin(theta) + t*cos(theta))*sigma ****
      call Cram_c_SMul(nb,(-y),
     >                  dcpl_mb(orb(1)),
     >                  tau_t)
      call Cram_c_SMul(nb,(-y),
     >                  dcpl_mb(orb(1)+npack1),
     >                  tau_t(npack1+1))
      call Cram_cc_daxpy(nb,x,
     >                   t,
     >                   tau_t)
      call Cram_cc_daxpy(nb,x,
     >                   t(npack1+1),
     >                   tau_t(npack1+1))
c      call Cram_c_SMul(nb,sigma,
c     >                  tau_t,
c     >                  tau_t)
c      call Cram_c_SMul(nb,sigma,
c     >                  tau_t(npack1+1),
c     >                  tau_t(npack1+1))
      call Cram_c_SMul1(nb,sigma,tau_t)
      call Cram_c_SMul1(nb,sigma,tau_t(npack1+1))

*     **** orb2 = orb*cos(theta) + t*sin(theta) ****
      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)),
     >                  dbl_mb(psi_ptr))
      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)+npack1),
     >                  dbl_mb(psi_ptr+shifts))
      call Cram_cc_daxpy(nb,y,
     >                   t,
     >                   dbl_mb(psi_ptr))
      call Cram_cc_daxpy(nb,y,
     >                   t(npack1+1),
     >                   dbl_mb(psi_ptr+shifts))


*     **** release stack memory ****
      value =           BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(orb(2))
      if (.not. value) call errquit(
     >     'psi_linesearch_update2com: popping stack memory',
     >     1, MA_ERR)
      return
      end

*     ***********************************
*     *                                 *
*     *      cpsi_KS_minimize_orb2_com  *
*     *                                 *
*     ***********************************

*    This routine performs a KS update on orbital i
*
      subroutine cpsi_KS_minimize_orb_2com(psi_number,
     >                             precondition,maxiteration,
     >                             maxerror,maxpsi_error,perror,nb,i,
     >                             error_out,psi_error,e0)
      implicit none
      integer psi_number
      logical precondition
      integer maxiteration
      real*8  maxerror,maxpsi_error,perror
      integer nb,i
      real*8 error_out,psi_error
      real*8 e0

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical value,done,oneloop
      integer it,shifts,ne1
      real*8 eold,error0,de0,lmbda_r0,lmbda_r1
      real*8 theta,sigma,tmp
      integer r1(2),t0(2),t(2),g(2)
      integer psi_ptr

*     **** external functions ****
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr

      lmbda_r0 = 1.0d0
      ne1=ne(1)
      shifts=2*ne1*npack1
      if (psi_number.eq.1) then
         psi_ptr=cpsi_data_get_ptr(psi1_tag,nb,i)
      else
         psi_ptr=cpsi_data_get_ptr(psi2_tag,nb,i)
      end if


      value = BA_push_get(mt_dcpl,2*npack1,'t0',t0(2),t0(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,2*npack1,'r1',r1(2),r1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,2*npack1,'g',g(2),g(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,2*npack1,'t',t(2),t(1))
      if (.not. value) call errquit(
     >     'cpsi_KS_minimize_orb: out of stack memory',0, MA_ERR)

      done = .false.
      error0 = 0.0d0
      e0 = 0.0d0
      theta = -3.14159d0/600.0d0
      it = 0
      psi_error = 10.0d0
 2    continue

c         write(*,*)"it=",it," I= ",i
         it = it + 1
         eold = e0

*        *** calculate residual (steepest descent) direction for a single band ***
         call cpsi_get_gradient_orb(psi_number,nb,i,dcpl_mb(g(1)))
         call Cram_cc_dot(nb,dbl_mb(psi_ptr),
     >                   dcpl_mb(g(1)),
     >                    e0)
         call Cram_cc_dot(nb,dbl_mb(psi_ptr+shifts),
     >                   dcpl_mb(g(1)+npack1),
     >                    tmp)
         e0 = -(e0+tmp)
c         write(*,*) "it,e0,eold,psi_error=",it,(-e0),(-eold),psi_error

         done = ((it.gt.maxiteration)
     >           .or.
     >           ((dabs(e0-eold).lt.maxerror).and.
     >           (psi_error.lt.maxpsi_error)))

         if (done) go to 4

         call Cram_c_Copy(nb,dcpl_mb(g(1)),dcpl_mb(r1(1)))
         call Cram_c_Copy(nb,dcpl_mb(g(1)+npack1),
     >          dcpl_mb(r1(1)+npack1))
         call Cram_cc_daxpy(nb,(e0),
     >                 dbl_mb(psi_ptr),
     >                 dcpl_mb(r1(1)))
         call Cram_cc_daxpy(nb,(e0),
     >                 dbl_mb(psi_ptr+shifts),
     >                 dcpl_mb(r1(1)+npack1))


*        *** determine conjuagate direction ***
         call Cram_cc_dot(nb,dcpl_mb(r1(1)),
     >                   dcpl_mb(r1(1)),
     >                   lmbda_r1)
         call Cram_cc_dot(nb,dcpl_mb(r1(1)+npack1),
     >                   dcpl_mb(r1(1)+npack1),
     >                   tmp)
         lmbda_r1=lmbda_r1+tmp
         call Cram_c_Copy(nb,dcpl_mb(r1(1)),dcpl_mb(t(1)))
         call Cram_c_Copy(nb,dcpl_mb(r1(1)+npack1),
     >       dcpl_mb(t(1)+npack1))

         if (it.gt.1) then
         call Cram_cc_daxpy(nb,(lmbda_r1/lmbda_r0),
     >                   dcpl_mb(t0(1)),
     >                   dcpl_mb(t(1)))
         call Cram_cc_daxpy(nb,(lmbda_r1/lmbda_r0),
     >                   dcpl_mb(t0(1)+npack1),
     >                   dcpl_mb(t(1)+npack1))
         end if
         lmbda_r0 = lmbda_r1
         oneloop = .true.
 3       call Cram_c_Copy(nb,dcpl_mb(t(1)),dcpl_mb(t0(1)))
         call Cram_c_Copy(nb,dcpl_mb(t(1)+npack1),
     >         dcpl_mb(t0(1)+npack1))



*        *** normalize search direction, t ****
         call cpsi_project_out_orb_2com(psi_number,nb,i,dcpl_mb(t(1)))
         call Cram_cc_dot(nb,dcpl_mb(t(1)),
     >                   dcpl_mb(t(1)),
     >                   sigma)
         call Cram_cc_dot(nb,dcpl_mb(t(1)+npack1),
     >                   dcpl_mb(t(1)+npack1),
     >                   tmp)
         sigma = sigma + tmp
         sigma = dsqrt(sigma)
         de0 = 1.0d0/sigma
c         call Cram_c_SMul(nb,de0,dcpl_mb(t(1)),dcpl_mb(t(1)))
c         call Cram_c_SMul(nb,de0,dcpl_mb(t(1)+npack1),
c     >       dcpl_mb(t(1)+npack1))
         call Cram_c_SMul1(nb,de0,dcpl_mb(t(1)))
         call Cram_c_SMul1(nb,de0,dcpl_mb(t(1)+npack1))


*        **** compute de0 = <t|g> ****
         call Cram_cc_dot(nb,dcpl_mb(t(1)),
     >                   dcpl_mb(g(1)),
     >                   de0)
         call Cram_cc_dot(nb,dcpl_mb(t(1)+npack1),
     >                   dcpl_mb(g(1)+npack1),
     >                   tmp)
         de0=de0+tmp
c         write(*,*)"de0=",de0,(de0-tmp),tmp
*        *** bad direction ***
         if ((de0.lt.0.0d0).and.oneloop) then
           call Cram_c_Copy(nb,dcpl_mb(g(1)),dcpl_mb(t(1)))
           call Cram_c_Copy(nb,dcpl_mb(g(1)+npack1),
     >           dcpl_mb(t(1)+npack1))
           oneloop = .false.
           go to 3
         end if

         de0 = -2.0d0*de0
         call cpsi_linesearch_minimize2com(psi_number,nb,i,
     >                              theta,e0,de0,
     >                              dcpl_mb(t(1)),psi_error)

      go to 2


*     **** release stack memory ****
 4    value =           BA_pop_stack(t(2))
      value = value.and.BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(r1(2))
      value = value.and.BA_pop_stack(t0(2))
      if (.not. value) call errquit(
     >     'psi_KS_minimize_orb: popping stack memory',1, MA_ERR)

c      write(*,*) "iterations=",it," eig=",e0," error=",error_out,
c     >           theta
      error_out = dabs(e0-eold)
      e0 = -e0
      return
      end



*     ***********************************
*     *                                 *
*     *      psi_linesearch_minimize    *
*     *                                 *
*     ***********************************

*    This routine performs a linesearch on orbital i, in the direction t.
* This routine is needed for a KS minimizer.
*  e0 = <orb|g>
*  de0 = 2*<t|g>
*
      subroutine cpsi_linesearch_minimize2com(psi_number,nb,i,
     >                                   theta,e0,de0,
     >                                    t,psi_error)
      implicit none
      integer psi_number
      integer nb,i
      real*8  theta
      real*8  e0,de0
      complex*16 t(*)     !search direction
      real*8     psi_error


#include "bafdecls.fh"
#include "cpsi_common.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer orb(2),g(2),psi_ptr,shifts
      real*8 x,y,pi,e1
      real*8 tmp
*     **** external functions ****
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr

      shifts=2*npack1*ne(1)
      if (psi_number.eq.1) then
         psi_ptr=cpsi_data_get_ptr(psi1_tag,nb,i)
      else
         psi_ptr=cpsi_data_get_ptr(psi2_tag,nb,i)
      end if

      pi = 4.0d0*datan(1.0d0)
 
*     **** allocate stack memory ****
      value = BA_push_get(mt_dcpl,2*npack1,'orb',orb(2),orb(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,2*npack1,'g',g(2),g(1))
      if (.not. value) call errquit(
     >     'psi_linesearch_minimize: out of stack memory',0, MA_ERR)

  
      call Cram_c_Copy(nb,dbl_mb(psi_ptr),dcpl_mb(orb(1)))
      call Cram_c_Copy(nb,dbl_mb(psi_ptr+shifts),dcpl_mb(orb(1)+npack1))

*     **** orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
      !theta = pi/300.0d0
      x = cos(theta)
      y = sin(theta)
      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)),
     >                  dbl_mb(psi_ptr))
      call Cram_cc_daxpy(nb,y,
     >                   t,
     >                   dbl_mb(psi_ptr))
      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)+npack1),
     >                  dbl_mb(psi_ptr+shifts))
      call Cram_cc_daxpy(nb,y,
     >                   t(npack1+1),
     >                   dbl_mb(psi_ptr+shifts))

*     *** determine theta ***
      call cpsi_get_gradient_orb(psi_number,nb,i,dcpl_mb(g(1)))

      call Cram_cc_dot(nb,dbl_mb(psi_ptr),
     >                   dcpl_mb(g(1)),
     >                   e1)
      call Cram_cc_dot(nb,dbl_mb(psi_ptr+shifts),
     >                   dcpl_mb(g(1)+npack1),
     >                   tmp)
      e1 = e1 + tmp
      e1 = -e1
c      write(*,*)"line nb,i,e1= ",nb,i,e1
      call flush(6)
      x = (e0 - e1 + 0.5d0*de0*dsin(2*theta))
     >    /(1.0d0-dcos(2*theta))
      theta = 0.5d0*datan(0.5d0*de0/x)


*     **** orb2 = orb*cos(theta) + t*sin(theta) ****
      x = cos(theta)
      y = sin(theta)
      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)),
     >                  dbl_mb(psi_ptr))
      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)+npack1),
     >                  dbl_mb(psi_ptr+shifts))
      call Cram_cc_daxpy(nb,y,
     >                   t,
     >                   dbl_mb(psi_ptr))
      call Cram_cc_daxpy(nb,y,
     >                   t(npack1+1),
     >                   dbl_mb(psi_ptr+shifts))

*     **** calculated psi_error = <(orb-psi)|(orb-psi)>  ****
      call Cram_cc_daxpy(nb,(-1.0d0),dbl_mb(psi_ptr),dcpl_mb(orb(1)))
      call Cram_cc_dot(nb,dcpl_mb(orb(1)),dcpl_mb(orb(1)),psi_error)
      call Cram_cc_daxpy(nb,(-1.0d0),dbl_mb(psi_ptr+shifts),
     >                               dcpl_mb(orb(1)+npack1))
      call Cram_cc_dot(nb,dcpl_mb(orb(1)+npack1),
     >                    dcpl_mb(orb(1)+npack1),
     >                    tmp)
      psi_error=psi_error+tmp

*     **** release stack memory ****
      value =           BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(orb(2))
      if (.not. value) call errquit(
     >     'psi_linesearch_minimize: popping stack memory',1, MA_ERR)

      return
      end

      subroutine cpsi_sort_minimize2com()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical value
      integer i,j,nb,shifts,ne1,indx,jndx
      integer r1(2),psi1,eig1
      real*8  ei,ej

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      value = BA_push_get(mt_dcpl,npack1*2,'r1',r1(2),r1(1))
      if (.not. value) call errquit(
     >     'psi_sort_sort_minimize2com: out of stack memory',0,MA_ERR)

      ne1=ne(1)
      shifts=2*ne1*npack1
      do nb=1,nbrillq

        psi1 = cpsi_data_get_chnk(psi1_tag,nb)
        eig1 = cpsi_data_get_chnk( eig_tag,nb)

        !*** Bubble sort ***
        do i=1,ne1
         do j=i+1,ne1
           ei = dbl_mb(eig1+i-1)
           ej = dbl_mb(eig1+j-1)

           !*** swap ***
           if (ej.lt.ei) then
             dbl_mb(eig1+i-1) = ej
             dbl_mb(eig1+j-1) = ei
             indx=(i-1)*2*npack1
             jndx=(j-1)*2*npack1
             call Cram_c_Copy(nb,dbl_mb(psi1+indx),
     >                           dcpl_mb(r1(1)))
             call Cram_c_Copy(nb,dbl_mb(psi1+indx+shifts),
     >                           dcpl_mb(r1(1)+npack1))
             call Cram_c_Copy(nb,dbl_mb(psi1+jndx),
     >                           dbl_mb(psi1+indx))
             call Cram_c_Copy(nb,dbl_mb(psi1+jndx+shifts),
     >                           dbl_mb(psi1+indx+shifts))
             call Cram_c_Copy(nb,dcpl_mb(r1(1)),
     >                           dbl_mb(psi1+jndx))
             call Cram_c_Copy(nb,dcpl_mb(r1(1)+npack1),
     >                           dbl_mb(psi1+jndx+shifts))
           end if

         end do
        end do

      end do

      value = BA_pop_stack(r1(2))
      if (.not. value) call errquit(
     >     'cpsi_sort_minimize: popping stack memory',1, MA_ERR)
      return
      end



************************ virtural orbital Part ************************


*     ***********************************
*     *                                 *
*     *     cpsi_minimize_virtual       *
*     *                                 *
*     ***********************************

      subroutine cpsi_minimize_virtual()
      implicit none

#include "bafdecls.fh"
#include "cpsi_common.fh"
      
      !*** local variables ***
      integer maxit_orb,taskid_k
      integer ii,l,nb,epsi_ptr,eig_ptr
      real*8  sum,maxerror,error_out,e0

      !*** external functions ***
      integer  cpsi_data_get_ptr
      real*8   control_tole
      external cpsi_data_get_ptr
      external control_tole

      maxit_orb=120
      maxerror = control_tole()
      call Parallel3d_taskid_k(taskid_k)

      do nb=1,nbrillq
      do ii=1,(ne_excited(1)+ne_excited(2))
         epsi_ptr = cpsi_data_get_ptr(psi1_excited_tag,nb,ii)
         eig_ptr  = cpsi_data_get_ptr(eig_excited_tag,nb,ii)
     
         !*** orthogonalize to lower orbitals  ****
         call cpsi_project_out_virtual1(nb,ii,dbl_mb(epsi_ptr))

         !*** normalize ****
         call Cram_cc_dot(nb,
     >            dbl_mb(epsi_ptr),
     >            dbl_mb(epsi_ptr),
     >            sum)
         sum = 1.0d0/dsqrt(sum)
         call Cram_c_SMul1(nb,sum,dbl_mb(epsi_ptr))


         !*** minimize orbital ****
          l = 0
 2        call cpsi_KS_update_virtual(maxit_orb,
     >                               maxerror,
     >                               0.001d0,nb,ii,error_out,e0)
          l = l+1
          if ((error_out.gt.maxerror).and.(l.le.4)) go to 2

          dbl_mb(eig_ptr) = e0

      end do
      end do
      call cpsi_sort_virtual()

      return
      end
*
c      subroutine cpsi_check_orthodebug(nb,i,Horb)
c      implicit none
c      integer nb,i
c      complex*16 Horb(*)
c
c#include "bafdecls.fh"
c#include "errquit.fh"
c#include "cpsi_common.fh"
c
c      integer x(2),ii,ms,n,psi_ptr,epsi_ptr,nshift
c
c*     **** external functions ****
c      integer  cpsi_data_get_chnk,cpsi_data_nsize
c      external cpsi_data_get_chnk,cpsi_data_nsize
c
c*     **** allocate stack memory ****
c      if(.not.BA_push_get(mt_dcpl,(neq(1)+ne_excited(1)),'x',x(2),x(1)))
c     > call errquit('cpsi_project_out_virtual1: out of stack memory',0,
c     >       MA_ERR)
c
c      psi_ptr  = cpsi_data_get_chnk(psi1_tag,nb)
c      epsi_ptr = cpsi_data_get_chnk(psi1_excited_tag,nb)
c      nshift   = cpsi_data_nsize(psi1_tag)
c      if (i.le.ne_excited(1)) then
c        ii = i-1
c        ms = 1
c      else
c        ms = 2
c        ii = i-ne_excited(1)-1
c        psi_ptr  = psi_ptr  + 2*ne(1)*npack1
c        epsi_ptr = epsi_ptr + 2*ne_excited(1)*npack1
c      end if
c      call Cram_cc_nzdot(nb,ne(ms),
c     >            dbl_mb(psi_ptr),
c     >            Horb,
c     >            dcpl_mb(x(1)))
c      do n=1,ne(ms)
c         write(*,*) "occ,n,x=",n,dcpl_mb(x(1)+n-1)
c      end do
c
c      !**** project out virtual orbitals ****
c      if (ii.gt.0) then
c         call Cram_cc_nzdot(nb,ii,
c     >               dbl_mb(epsi_ptr),
c     >               Horb,
c     >               dcpl_mb(x(1)))
c         do n=1,ii
c            write(*,*) "virt,n,x=",n,dcpl_mb(x(1)+n-1)
c         end do
c      end if
c
c*     **** release stack memory ****
c      if (.not.BA_pop_stack(x(2)))
c     > call errquit('cpsi_project_out_virtual1:popping stack memory',0,
c     &       MA_ERR)
c
c      return
c      end



*     *******************************************
*     *                                         *
*     *          cpsi_project_out_virtual1      *
*     *                                         *
*     *******************************************
*
*    This routine projects out non-orthogonal components of Horb.
* This routine is needed for a KS minimizer.
*
      subroutine cpsi_project_out_virtual1(nb,i,Horb)
      implicit none
      integer nb,i
      complex*16 Horb(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

      integer x(2),ii,ms,n,psi_ptr,epsi_ptr,nshift
      complex*16  xxx

*     **** external functions ****
      integer  cpsi_data_get_chnk,cpsi_data_nsize
      external cpsi_data_get_chnk,cpsi_data_nsize

*     **** allocate stack memory ****
      if(.not.BA_push_get(mt_dcpl,(neq(1)+ne_excited(1)),'x',x(2),x(1)))
     > call errquit('cpsi_project_out_virtual1: out of stack memory',0,
     >       MA_ERR)

      psi_ptr  = cpsi_data_get_chnk(psi1_tag,nb)
      epsi_ptr = cpsi_data_get_chnk(psi1_excited_tag,nb)
      nshift   = cpsi_data_nsize(psi1_tag)
      if (i.le.ne_excited(1)) then
        ms = 1
        ii = i-1
      else
        ms = 2
        ii = i-ne_excited(1)-1
        psi_ptr  = psi_ptr  + 2*ne(1)*npack1
        epsi_ptr = epsi_ptr + 2*ne_excited(1)*npack1
      end if

      !**** project out occupied orbitals ****
c      call Cram_cc_nzdot(nb,ne(ms),
c     >            dbl_mb(psi_ptr),
c     >            Horb,
c     >            dcpl_mb(x(1)))
      do n=1,ne(ms)
         call Cram_cc_zdot(nb,dbl_mb(psi_ptr),Horb,xxx)
         call Cram_cc_zaxpy(nb,(-xxx),dbl_mb(psi_ptr),Horb)
         psi_ptr = psi_ptr + nshift
      end do

      !**** project out virtual orbitals ****
      if (ii.gt.0) then
c         call Cram_cc_nzdot(nb,ii,
c     >               dbl_mb(epsi_ptr),
c     >               Horb,
c     >               dcpl_mb(x(1)))
         do n=1,ii
            call Cram_cc_zdot(nb,dbl_mb(epsi_ptr),Horb,xxx)
            call Cram_cc_zaxpy(nb,(-xxx),dbl_mb(epsi_ptr),Horb)
            epsi_ptr = epsi_ptr + nshift
         end do
      end if

*     **** release stack memory ****
      if (.not.BA_pop_stack(x(2)))
     > call errquit('cpsi_project_out_virtual1:popping stack memory',0,
     &       MA_ERR)

      return
      end 


*     *******************************************
*     *                                         *
*     *          cpsi_project_out_virtual      *
*     *                                         *
*     *******************************************
*     
*    This routine projects out non-orthogonal components of Horb.
* This routine is needed for a KS minimizer.
* 
      subroutine cpsi_project_out_virtual(nb,i,Horb)
      implicit none
      integer nb,i
      complex*16 Horb(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer x(2),ii,n,psi_ptr,epsi_ptr,nshift,ms
      complex*16 xxx

*     **** external functions ****
      integer  cpsi_data_get_chnk,cpsi_data_nsize
      external cpsi_data_get_chnk,cpsi_data_nsize

*     **** allocate stack memory ****
      if (.not.BA_push_get(mt_dcpl,neq(1)+ne_excited(1),'x',x(2),x(1)))
     > call errquit('cpsi_project_out_virtual: out of stack memory',0,
     >       MA_ERR)

      psi_ptr  = cpsi_data_get_chnk(psi1_tag,nb)
      epsi_ptr = cpsi_data_get_chnk(psi1_excited_tag,nb)
      nshift   = cpsi_data_nsize(psi1_tag)
      if (i.le.ne_excited(1)) then
        ii = i
        ms = 1
      else
        ii = i-ne_excited(1)
        ms = 2
        psi_ptr  = psi_ptr  + 2*ne(1)*npack1
        epsi_ptr = epsi_ptr + 2*ne_excited(1)*npack1
      end if

      !**** project out occupied orbitals ****
c      call Cram_cc_nzdot(nb,ne(ms),
c     >            dbl_mb(psi_ptr),
c     >            Horb,
c     >            dcpl_mb(x(1)))
      do n=1,ne(ms)
         call Cram_cc_zdot(nb,dbl_mb(psi_ptr),Horb,xxx)
         call Cram_cc_zaxpy(nb,(-xxx),dbl_mb(psi_ptr),Horb)
         psi_ptr = psi_ptr + nshift
      end do

      !**** project out virtual orbitals ****
c      call Cram_cc_nzdot(nb,ii, 
c     >            dbl_mb(epsi_ptr),
c     >            Horb,
c     >            dcpl_mb(x(1)))
      do n=1,ii
         call Cram_cc_zdot(nb,dbl_mb(epsi_ptr),Horb,xxx)
         call Cram_cc_zaxpy(nb,(-xxx),dbl_mb(epsi_ptr),Horb)
         epsi_ptr = epsi_ptr + nshift
      end do

*     **** release stack memory ****
      if (.not.BA_pop_stack(x(2)))
     > call errquit('cpsi_project_out_virtual:popping stack memory',0,
     &       MA_ERR)

      return
      end



      

*     ***********************************
*     *                                 *
*     *      cpsi_KS_update_virtual     *
*     *                                 *
*     ***********************************

*    This routine performs a KS update on orbital i
*
      subroutine cpsi_KS_update_virtual(maxiteration,
     >                             maxerror,perror,nb,i,
     >                             error_out,e0)
      implicit none
      integer maxiteration
      real*8  maxerror,perror
      integer nb,i
      real*8 error_out,eo

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"


*     **** local variables ****
      logical value,done,oneloop
      integer it
      real*8 e0,eold,percent_error,error0,de0,lmbda_r0,lmbda_r1
      real*8 theta
      integer r1(2),t0(2),t(2),g(2)
      integer psi_ptr

*     **** external functions ****
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr

      psi_ptr=cpsi_data_get_ptr(psi1_excited_tag,nb,i)

      lmbda_r0 = 1.0d0

      value = BA_push_get(mt_dcpl,npack1,'t0',t0(2),t0(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1,'g',g(2),g(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1,'t',t(2),t(1))
      if (.not. value) call errquit(
     >     'cpsi_KS_update_virtual: out of stack memory',0, MA_ERR)

      done = .false.
      error0 = 0.0d0
      e0 = 0.0d0
      theta = -3.14159d0/600.0d0
      it = 0
 2    continue
 
         it = it + 1
         eold = e0 

*        *** calculate residual (steepest descent) direction for a single band ***
         call cpsi_get_gradient_virtual(nb,i,dcpl_mb(g(1)))
         call Cram_cc_dot(nb,dbl_mb(psi_ptr),
     >                   dcpl_mb(g(1)),
     >                    e0)
         e0 = -e0

         percent_error=0.0d0
         if(error0.ne.0.0d0)
     A      percent_error = dabs(e0-eold)/error0

         done = ((it.gt.maxiteration)
     >           .or.
     >           (dabs(e0-eold).lt.maxerror))

         if (done) go to 4

         call Cram_c_Copy(nb,dcpl_mb(g(1)),dcpl_mb(r1(1)))
         call Cram_cc_daxpy(nb,e0,dbl_mb(psi_ptr),dcpl_mb(r1(1)))


*        *** determine conjuagate direction ***
         call Cram_cc_dot(nb,dcpl_mb(r1(1)),
     >                       dcpl_mb(r1(1)),
     >                       lmbda_r1)
         call Cram_c_Copy(nb,dcpl_mb(r1(1)),dcpl_mb(t(1)))

         if (it.gt.1) then
         call Cram_cc_daxpy(nb,(lmbda_r1/lmbda_r0),
     >                         dcpl_mb(t0(1)),
     >                         dcpl_mb(t(1)))
         end if
         lmbda_r0 = lmbda_r1
         oneloop = .true.
 3       call Cram_c_Copy(nb,dcpl_mb(t(1)),dcpl_mb(t0(1)))

*        *** normalize search direction, t ****
         call cpsi_project_out_virtual(nb,i,dcpl_mb(t(1)))
         call Cram_cc_dot(nb,dcpl_mb(t(1)),
     >                       dcpl_mb(t(1)),
     >                       de0)
         de0 = 1.0d0/dsqrt(de0)
         call Cram_c_SMul1(nb,de0,dcpl_mb(t(1)))
         call Cram_cc_dot(nb,dcpl_mb(t(1)),
     >                       dcpl_mb(g(1)),
     >                       de0)

*        *** bad direction ***
         if ((de0.lt.0.0d0).and.oneloop) then
           call Cram_c_Copy(nb,dcpl_mb(g(1)),dcpl_mb(t(1)))
           oneloop = .false.
           go to 3
         end if

         de0 = -2.0d0*de0
         call cpsi_linesearch_virtual(nb,i,
     >                              theta,e0,de0,
     >                              dcpl_mb(t(1)))

      go to 2


*     **** release stack memory ****
 4    value =           BA_pop_stack(t(2))
      value = value.and.BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(r1(2))
      value = value.and.BA_pop_stack(t0(2))
      if (.not.value) call errquit(
     >     'cpsi_KS_update_virtual: popping stack memory',1,MA_ERR)

      error_out = dabs(e0-eold)
      e0 = -e0
      return
      end


*     ***********************************
*     *                                 *
*     *      cpsi_get_gradient_virtual   *
*     *                                 *
*     ***********************************

*    This routine returns the Hpsi(i).  
* This routine is needed for a KS minimizer.
*
      subroutine cpsi_get_gradient_virtual(nb,ii,Horb)
      implicit none
      integer nb,ii
      complex*16 Horb(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer psi_ptr,ms

*     **** external functions ****
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr

      psi_ptr=cpsi_data_get_ptr(psi1_excited_tag,nb,ii)
      if (ii.le.ne_excited(1)) then
         ms = 1
      else
         ms = 2
      end if

      call c_electron_get_gradient_virtual(nb,ms,dbl_mb(psi_ptr),Horb)
      return
      end






*     ***********************************
*     *                                 *
*     *      cpsi_linesearch_virtual    *
*     *                                 *
*     ***********************************

*    This routine performs a linesearch on orbital ii, in the direction t.  
* This routine is needed for a KS minimizer.
*  e0 = <orb|g>
*  de0 = 2*<t|g>
*
      subroutine cpsi_linesearch_virtual(nb,ii,theta,e0,de0,t)
      implicit none
      integer nb,ii
      real*8  theta
      real*8  e0,de0
      complex*16 t(*) !search direction

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      logical value
      integer orb(2),g(2),psi_ptr
      real*8 x,y,pi,dtheta_min,e1

*     **** external functions ****
      integer  cpsi_data_get_ptr
      external cpsi_data_get_ptr

      psi_ptr=cpsi_data_get_ptr(psi1_excited_tag,nb,ii)

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
     >  'cpsi_linesearch_virtual: out of stack memory',0,MA_ERR)


      call Cram_c_Copy(nb,dbl_mb(psi_ptr),
     >                    dcpl_mb(orb(1)))

*     **** orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
  10  x = cos(theta)
      y = sin(theta)
      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)),
     >                  dbl_mb(psi_ptr))
      call Cram_cc_daxpy(nb,y,
     >                   t,
     >                   dbl_mb(psi_ptr))

*     *** determine theta ***
      call cpsi_get_gradient_virtual(nb,ii,dcpl_mb(g(1)))
      call Cram_cc_dot(nb,dbl_mb(psi_ptr),
     >                   dcpl_mb(g(1)),
     >                   e1)
      e1 = -e1

      x = (e0 - e1 + 0.5d0*de0*sin(2*theta))
     >    /(1.0d0-cos(2*theta))
      theta = 0.5d0*datan(0.5d0*de0/x)

*     **** orb2 = orb*cos(theta) + t*sin(theta) ****
      x = cos(theta)
      y = sin(theta)
      call Cram_c_SMul(nb,x,
     >                  dcpl_mb(orb(1)),
     >                  dbl_mb(psi_ptr))
      call Cram_cc_daxpy(nb,y,
     >                   t,
     >                   dbl_mb(psi_ptr))

*     **** release stack memory ****
      value =           BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(orb(2))
      if (.not. value) call errquit(
     >     'cpsi_linesearch_virtual: popping stack memory',1, MA_ERR)

      return
      end


*     ***********************************
*     *                                 *
*     *        cpsi_sort_virtual        *
*     *                                 *
*     ***********************************

      subroutine cpsi_sort_virtual()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsi_common.fh"

      logical value
      integer i,j,ii,jj,ms,nb,nshift
      integer r1(2),psi1_shift,eig_shift
      real*8  ei,ej

*      ***** external functions ****
       integer  cpsi_data_get_chnk
       external cpsi_data_get_chnk

c      if (spin_orbit) then
c        call cpsi_sort_minimize2com()
c        return
c      end if

      if (.not.BA_push_get(mt_dcpl,npack1,'r1',r1(2),r1(1)))
     >  call errquit(
     >     'cpsi_sort_sort_virtual: out of stack memory',0,MA_ERR)

      nshift = 2*npack1
      do nb=1,nbrillq
        psi1_shift = cpsi_data_get_chnk(psi1_excited_tag,nb)
        eig_shift  = cpsi_data_get_chnk(eig_excited_tag,nb)
        do ms=1,ispin

        !*** Bubble sort ***
        do ii=1,ne_excited(ms)
         do jj=ii+1,ne_excited(ms)
           i = ii + (ms-1)*ne_excited(1)
           j = jj + (ms-1)*ne_excited(1)
           ei = dbl_mb(eig_shift+i-1)
           ej = dbl_mb(eig_shift+j-1)

           !*** swap ***
           if (ej.lt.ei) then
             dbl_mb(eig_shift+i-1) = ej
             dbl_mb(eig_shift+j-1) = ei
             call Cram_c_Copy(nb,dbl_mb(psi1_shift+(i-1)*nshift),
     >                           dcpl_mb(r1(1)))
             call Cram_c_Copy(nb,dbl_mb(psi1_shift+(j-1)*nshift),
     >                           dbl_mb(psi1_shift+(i-1)*nshift))
             call Cram_c_Copy(nb,dcpl_mb(r1(1)),
     >                           dbl_mb(psi1_shift+(j-1)*nshift))
           end if

         end do
        end do
      end do
      end do

      if (.not.BA_pop_stack(r1(2)))
     >  call errquit(
     >     'cpsi_sort_virtual: popping stack memory',1, MA_ERR)
      return
      end





