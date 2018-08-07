*
* $Id$
*

* $Log: not supported by cvs2svn $
* Revision 1.44  2008/09/15 20:52:32  bylaska
* ..paw fixes...EJB
*
* Revision 1.43  2007/09/24 16:58:12  bylaska
* ...preliminary PAW modifications...
*   - basis file format changed
*   - .vpp formatting routines added to pspw
*
* - zdotc routines currently modified to tzdotc.
* ...EJB
*
* Revision 1.42  2006/02/11 02:50:46  bylaska
* GGA's using 1st derivative formulas have been added in core part of PAW....EJB
*
* Revision 1.41  2005/02/09 02:38:57  bylaska
* ..............EJB
*
* Revision 1.40  2004/11/08 23:37:41  bylaska
* Bug fix in pspw_hfx found by M. Hackler.
* PBE0 has been implemented.
*
*  ........EJB
*
* Revision 1.39  2004/08/12 18:39:41  bylaska
* A prototype of a Grassmann CG paw minimizer (i.e. nwpw:minimizer 1) has been added.
* The code is similar to the CG minimizer in pspw, but differences exist
* because the residual |R> = (1 - S|psi><psi|)|Hpsi> is not the same as the
* tangent vector |T> = (1 - |psi><psi|S)|R>.
*
* Forces still need to be implemented.
*
* ...EJB
*
* Revision 1.38  2003/10/24 18:45:24  bylaska
*
* Aperiodic convolution capability has been added to PAW .... EJB
*
* Revision 1.37  2003/10/21 02:05:15  marat
* switched to new errquit by running global replace operation
* see the script below (note it will not work on multiline errquit calls)
* *********************************************************
* #!/bin/sh
*
* e=`find . -name "*F" -print`
*
* for f in $e
* do
* cp $f $f.bak
* sed  's|\(^[ ].*call[ ]*errquit([^,]*\)\(,[^,]*\)\()\)|\1,0\2\3|' $f.bak > $f
* #rm $f.bak
* done
* **********************************************************
*
* Revision 1.36  2003/03/25 19:43:31  bylaska
* bug fix...EJB
*
* Revision 1.35  2003/03/22 02:30:01  bylaska
* paw cpsd program finished....
* The nwpw directory structure is ready to be checked into 4.5 release tree.
*
* ....EJB
*
* Revision 1.34  2003/03/21 23:41:13  bylaska
*
* paw updates ...EJB
*
* Revision 1.33  2003/03/15 02:14:44  bylaska
* orthonormalization checking fixed to work with forces...EJB
*
* Revision 1.32  2003/03/15 01:47:43  bylaska
* steepest descent loop has been modified for the inclusion of forces....
* Lagrange Multipliers require a recalculation of the phase factors
* after call to paw_overlap_S and before paw_psi_lagrange.....EJB
*
* Revision 1.31  2003/03/14 01:20:59  marat
* moved call to paw_force_solve after nonlocal
* matrices have been calculated
* MV
*
* Revision 1.30  2003/03/11 17:57:10  bylaska
* updates...EJB
*
* Revision 1.29  2003/03/07 20:51:10  bylaska
* Code cleanup...0.0 changed to 0.0d0 in paw_xc.F
* Tangent vector now used for SD with Gram-schmidt.
* ....EJB
*
* Revision 1.28  2003/03/06 01:46:39  bylaska
* bug fix in paw_energy_core_atom...and ma_chop_stack changed to ma_pop_stack
* ...EJB
*
* Revision 1.27  2003/03/05 23:16:31  bylaska
* Commented out write statements and other minor fixes.....
* self-consistent loop looks like it is working.....
* ....EJB
*
* Revision 1.26  2003/03/04 00:04:03  marat
* added printouts for atomic potenitials
* for debug purposes
* MV
*
* Revision 1.25  2003/02/24 22:38:51  bylaska
* Fixed bugs in ehartr_pw calculation....EJB
*
* Revision 1.24  2003/02/24 21:58:59  marat
* ...
* MV
*
* Revision 1.23  2003/02/24 21:05:00  bylaska
*  $Log: added to cvs output....EJB
*



      subroutine paw_inner_loop(ispin,ne,
     >                      npack1,nfft3d,nemax,
     >                      psi1,psi2,dn,
     >                      dn_cmp_smooth,
     >                      it_in,E,deltae,deltac,deltar,
     >                      hml,lmd,lmd1,first_iteration,
     >                      psi_r,Hpsi)
      implicit none
      integer    ispin,ne(2)
      integer    npack1,nfft3d,nemax
      complex*16 psi1(npack1,nemax)
      complex*16 psi2(npack1,nemax)
      real*8     dn(2*nfft3d,2)
      real*8     dn_cmp_smooth(2*nfft3d)
      integer    it_in
      real*8     E(*)
      real*8     deltae,deltac,deltar
      real*8     hml(2*nemax*nemax)
      real*8     lmd(2*nemax*nemax),lmd1(2*nemax*nemax)
      logical    first_iteration

*     **** very big workspace variables ****
      real*8     psi_r(2*nfft3d,nemax)
      complex*16 Hpsi(npack1,nemax)


#include "bafdecls.fh"
#include "paw_energy_kin_atom.fh"
#include "paw_energy_vloc_atom.fh"
#include "paw_energy_ion_atom.fh"
#include "paw_energy_core_atom.fh"
#include "paw_energy_hartree_atom.fh"
#include "paw_xc.fh"
#include "nwpwxc.fh"
#include "util.fh"


*     **** local variables ****
      logical move
      logical use_lda
      integer n2ft3d,np
      integer i,j,ii,jj,n,n1(2),n2(2),it,ms,nn,ierr
      integer nx,ny,nz
      integer index,indext
      double precision evloc_pw,evloc_atom,occ(1)
      real*8  sum,Eold,eorbit,ehartr_pw,eke,enlocal
      real*8  exc,exc2,pxc,pxc2,dte,scal1,scal2,dv,dt
      real*8  deltamm,vzero
      double precision ekin_atom
      double precision eion_atom
      double precision ecore_atom
      double precision ecore_ion_atom
      double precision ecore_self_atom
      double precision ehartree_atom
      double precision exc_atom


*     **** MA local variables ****
      logical value,gram_schmidt
      integer tmp_L(2)
      integer tmp1(2),tmp2(2)
      integer vl(2),vh(2),vc(2),vcomp(2),dng(2)
      integer rho(2)
      integer xcp(2),xce(2),dnall(2)
      integer natmx,fion(2),ftest(2)
      integer sumi(2)
      integer npack0,gga

*     ***** external functions ****
      logical  control_move,control_gram_schmidt
      integer  ion_nion,control_gga
      real*8   control_time_step,control_fake_mass,ion_dti
      real*8   lattice_omega,coulomb_e,ewald_e
      external control_move,control_gram_schmidt
      external ion_nion,control_gga
      external control_time_step,control_fake_mass,ion_dti
      external lattice_omega,coulomb_e,ewald_e
      integer  control_version
      external control_version
      real*8   ion_ion_e
      external ion_ion_e
      real*8   paw_mult_energy_atom_comp !**no header file for paw_mult**
      real*8   paw_mult_energy_atom_self
      real*8   paw_mult_energy_atom_mult
      external paw_mult_energy_atom_comp
      external paw_mult_energy_atom_self
      external paw_mult_energy_atom_mult


      call Parallel_np(np)
      call Pack_npack(0,npack0)
      n2ft3d = 2*nfft3d
      deltamm = 0.0d0
      gga = control_gga()


      call nwpw_timing_start(12)
*     **** allocate MA local variables ****
      value = BA_push_get(mt_dbl,(8*nemax*nemax),
     >                     'tmp_L',tmp_L(2),tmp_L(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'tmp2',tmp2(2),tmp2(1))

      if (control_version().eq.3) then
       value = value.and.
     >        BA_push_get(mt_dcpl,(npack0),'vcomp',vcomp(2),vcomp(1))
       value = value.and.
     >        BA_push_get(mt_dcpl,(npack0),'vh',vh(2),vh(1))
       value = value.and.
     >        BA_push_get(mt_dcpl,(npack0),'vc',vc(2),vc(1))
      end if

      if (control_version().eq.4) then
       value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'vcomp',vcomp(2),vcomp(1))
       value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'vh',vh(2),vh(1))
       value = value.and.
     >        BA_push_get(mt_dcpl,(npack0),'vc',vc(2),vc(1))
      end if

      value = value.and.
     >        BA_push_get(mt_dcpl,(npack0),'vloc', vl(2), vl(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(n2ft3d),'rho',rho(2),rho(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(npack0),'dng',dng(2), dng(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'xcp',xcp(2), xcp(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'xce',xce(2), xce(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'dnall',dnall(2),dnall(1))
      natmx = ion_nion()
      value = value.and.
     >        BA_push_get(mt_dbl,(3*natmx),'fion',fion(2),fion(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(3*natmx),'ftest',ftest(2),ftest(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nemax),'sumi',sumi(2),sumi(1))

      if (.not. value) call errquit('out of stack memory',0,0)

      call nwpw_timing_end(12)
     
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      move         = control_move()
      gram_schmidt = control_gram_schmidt()

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

      dt = control_time_step()
      dte = dt/dsqrt(control_fake_mass())
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = scal1*lattice_omega()


*     ******************************************
*     ****                                  ****
*     ****   Start of steepest descent loop ****
*     ****                                  ****
*     ******************************************
      do it=1,it_in

*       **** shift wavefunction and atoms ****
        call dcopy(2*npack1*nemax,psi2,1,psi1,1)
        if (move) call ion_shift()
        !if (move) call phafac()
        if (move) call paw_set_mult_energy_coeff()

        call nwpw_timing_start(11)
*       *******************
*       **** get psi_r ****
*       *******************
        do n=n1(1),n2(ispin)
           call Pack_c_Copy(1,psi1(1,n),psi_r(1,n))
           call Pack_c_unpack(1,psi_r(1,n))
           call D3dB_cr_fft3b(1,psi_r(1,n))
           call D3dB_r_Zero_Ends(1,psi_r(1,n))
        end do

*       *******************
*       **** set overlaps *
*       *******************
        call paw_ovlp_coeff_set(psi1)
        call paw_ovlp_weights_set()
        !call paw_ovlp_weights_write(89)

*       *************************************
*       ****generate comp charge potential***
*       *************************************
        call paw_comp_charge_update()
        call paw_pot_comp_solve()
        !call paw_pot_comp_print()

       
*       *********************
*       **** generate dn ****
*       *********************
        call dcopy(ispin*n2ft3d,0.0d0,0,dn,1)
        do ms=1,ispin
           do n=n1(ms),n2(ms)
              do i=1,n2ft3d
                 dn(i,ms) = dn(i,ms) + scal2*(psi_r(i,n)**2)
              end do
           end do
           call D3dB_r_Zero_Ends(1,dn(1,ms))
        end do


*       **********************
*       **** generate dng ****
*       **********************
        call D3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dbl_mb(rho(1)))
        call D3dB_r_SMul(1,scal1,dbl_mb(rho(1)),dcpl_mb(tmp1(1)))
        call D3dB_rc_fft3f(1,dcpl_mb(tmp1(1)))
        call Pack_c_pack(0,dcpl_mb(tmp1(1)))
        call Pack_c_Copy(0,dcpl_mb(tmp1(1)),dcpl_mb(dng(1)))



*       *****************************************
*       **** generate local pseudopotential  ****
*       **** and also get force if move true ****
*       *****************************************
       call paw_vloc(dcpl_mb(vl(1)),
     >               move,
     >               dcpl_mb(dng(1)),
     >               dbl_mb(fion(1)))
       call Pack_cc_dot(0,dcpl_mb(dng(1)),dcpl_mb(vl(1)),evloc_pw)


*      ************************************
*      **** generate coulomb potential ****
*      ************************************
        
*      *** atomic portion ***
       call paw_pot_hartree_solve()
           
        call paw_mult_dn_cmp_get(dcpl_mb(tmp1(1)),
     >                           dn_cmp_smooth)
        if (control_version().eq.3)  then
           call Pack_cc_Sub(0,
     >                      dcpl_mb(tmp1(1)),
     >                      dn_cmp_smooth,
     >                      dcpl_mb(tmp2(1))) !** tmp2 = dn_cmp - dn_cmp_smooth
           call Pack_cc_Sum(0,
     >                      dn_cmp_smooth,
     >                      dcpl_mb(dng(1)),
     >                      dcpl_mb(tmp1(1))) !** tmp1 = dng+dn_cmp_smooth **

           !**** vh *****
           call coulomb_v(dcpl_mb(tmp1(1)), 
     >                    dcpl_mb(vh(1))) 

           !**** vcmp *****
           call coulomb_v(dcpl_mb(tmp2(1)),
     >                    dcpl_mb(vcomp(1)))
           call paw_mult_vzero(vzero)
           call Pack_c_setzero(0,vzero,dcpl_mb(vcomp(1)))

           call paw_mult_coeff_set(dcpl_mb(vh(1)),dcpl_mb(vcomp(1)))
           
           call Pack_cc_Sum(0,
     >                      dcpl_mb(vh(1)),
     >                      dcpl_mb(vcomp(1)),
     >                      dcpl_mb(vc(1)))
                   
        end if

        if (control_version().eq.4) then
        
           call Pack_cc_Sub(0,
     >                      dcpl_mb(tmp1(1)),
     >                      dn_cmp_smooth,
     >                      dcpl_mb(tmp2(1))) !** tmp2 = dn_cmp - dn_cmp_smooth
           call Pack_cc_Sum(0,
     >                      dn_cmp_smooth,
     >                      dcpl_mb(dng(1)),
     >                      dcpl_mb(tmp1(1))) !** tmp1 = dng+dn_cmp_smooth **

           call Pack_c_unpack(0,dcpl_mb(tmp1(1)))
           call D3dB_cr_fft3b(1,dcpl_mb(tmp1(1)))
           call D3dB_r_Zero_Ends(1,dcpl_mb(tmp1(1)))

           call coulomb2_v(dcpl_mb(tmp1(1)),dcpl_mb(vh(1)))

           call D3dB_rc_fft3f(1,dcpl_mb(vh(1)))
c           call D3dB_r_SMul(1,scal1,dcpl_mb(vh(1)),dcpl_mb(vh(1)))
           call D3dB_r_SMul1(1,scal1,dcpl_mb(vh(1)))
           call Pack_c_pack(0,dcpl_mb(vh(1)))


           call Pack_c_unpack(0,dcpl_mb(tmp2(1)))
           call D3dB_cr_fft3b(1,dcpl_mb(tmp2(1)))
           call D3dB_r_Zero_Ends(1,dcpl_mb(tmp2(1)))

           call coulomb2_v(dcpl_mb(tmp2(1)),dcpl_mb(vcomp(1)))

           call D3dB_rc_fft3f(1,dcpl_mb(vcomp(1)))
c           call D3dB_r_SMul(1,scal1,dcpl_mb(vcomp(1)),dcpl_mb(vcomp(1)))
           call D3dB_r_SMul1(1,scal1,dcpl_mb(vcomp(1)))
           call Pack_c_pack(0,dcpl_mb(vcomp(1)))

           call paw_mult_vzero(vzero)
           call Pack_c_setzero(0,vzero,dcpl_mb(vcomp(1)))

           call paw_mult_coeff_set(dcpl_mb(vh(1)),dcpl_mb(vcomp(1)))

           call Pack_cc_Sum(0,
     >                      dcpl_mb(vh(1)),
     >                      dcpl_mb(vcomp(1)),
     >                      dcpl_mb(vc(1)))

        end if


*       *************************************************
*       **** generate exchange-correlation potential ****
*       *************************************************

*       *** local portion ***
c        call paw_density_solve()
        call paw_xc_solve()
        !call paw_xc_print()
        
*       *** plane wave ***
       
        use_lda = ((.not.nwpwxc_is_on().and.gga.eq.0).or.
     &             (nwpwxc_is_on().and.nwpwxc_is_lda()))
        
        if (use_lda) then
          call vxc(n2ft3d,ispin,dn,
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)),
     >                      dcpl_mb(tmp1(1)))
        else
          call v_bwexc(gga,n2ft3d,ispin,dn,
     >                      1.0d0,1.0d0,
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)))
        end if
     

*       ******************
*       **** get Hpsi ****
*       ******************
        call nwpw_timing_start(13)
        call paw_psi_H(ispin,ne,psi1,psi_r,
     >             dcpl_mb(vl(1)),
     >             dcpl_mb(vc(1)),dbl_mb(xcp(1)),Hpsi,
     >             move,dbl_mb(fion(1)))
 
        
        !call paw_Gop_print()

*       ************************************
*       **** do a steepest descent step ****
*       ************************************
*   
*       **** if gram-schmidt make Hpsi a tangent vecotor ****
        if (gram_schmidt) then

          call paw_ovlp_S(n2(ispin),psi1,psi2)   ! psi2 = S*psi1
          do ms=1,ispin
            call Grsm_ggm_sym_dot(npack1,ne(ms),
     >                            psi1(1,n1(ms)),
     >                            Hpsi(1,n1(ms)),
     >                            dbl_mb(tmp_L(1)))
            call dscal(ne(ms)*ne(ms),(-1.0d0),dbl_mb(tmp_L(1)),1)
            call Grsm_gmg_Mul(npack1,ne(ms),
     >                          psi2(1,n1(ms)),
     >                          dbl_mb(tmp_L(1)),
     >                          psi_r)            ! psi_r =  -S*psi1*<psi1|Hpsi2>
            call Grsm_ggg_Sum(npack1,ne(ms),
     >                        Hpsi(1,n1(ms)),
     >                        psi_r,
     >                        psi2(1,n1(ms)))    ! psi2 = Hpsi - S*psi1*<psi1|Hpsi1>

c            call Grsm_gg_dscale(npack1,ne(ms),dte,
c     >                         psi2(1,n1(ms)),
c     >                         psi2(1,n1(ms)))
            call Grsm_gg_dScale1(npack1,ne(ms),dte,psi2(1,n1(ms)))

c            call Grsm_ggg_Sum(npack1,ne(ms),
c     >                        psi2(1,n1(ms)), 
c     >                        psi1(1,n1(ms)),
c     >                        psi2(1,n1(ms)))
            call Grsm_ggg_Sum2(npack1,ne(ms),
     >                        psi1(1,n1(ms)),
     >                        psi2(1,n1(ms)))
          end do !*ms*
       
*       **** else don't change Hpsi ***
        else 

          do n=1,n2(ispin)
            call Pack_c_SMul(1,dte,Hpsi(1,n),psi2(1,n))
c            call Pack_cc_Sum(1,psi2(1,n),psi1(1,n),psi2(1,n))
            call Pack_cc_Sum2(1,psi1(1,n),psi2(1,n))
          end do
        end if

        call nwpw_timing_end(13)

*       *******************************************
*       **** get ion forces and do steepest    ****
*       **** descent on ions                   ****
*       *******************************************

*       *********************
*       **** generate force *
*       *********************
         if (move) then

           call paw_mult_pw_force(dcpl_mb(vh(1)),
     >                            dcpl_mb(vcomp(1)),
     >                            dbl_mb(fion(1)))

           call paw_force_solve(psi1,dbl_mb(fion(1))) 


*           *** compute hamiltonian matrix if first iteration ****
           if (first_iteration) then
             n = ne(1)
             nn = n*n
             do ms=1,ispin
                do ii=n1(ms),n2(ms)
                  i = ii-n1(ms)
                  index = (i+1) + i*n + (ms-1)*nn
                  call Pack_cc_idot(1,psi1(1,ii),Hpsi(1,ii),sum)
 
                  hml(index) =  -sum
                  do jj=ii+1,n2(ms)
                     j = jj-n1(ms)
                     index  = (i+1) + j*n + (ms-1)*nn
                     indext = (j+1) + i*n + (ms-1)*nn
                     call Pack_cc_idot(1,psi1(1,ii),Hpsi(1,jj),sum)
 
                     hml(index)  =  -sum
                     hml(indext) =  -sum
                  end do 
                end do
             end do
             if (np.gt.1)  call D3dB_Vector_SumAll((ispin*nn),hml)
             call dcopy(2*nemax*nemax,hml,1,lmd1,1)
             call dcopy(2*nemax*nemax,hml,1,lmd,1)
             first_iteration = .false.
           end if

           call dcopy(2*nemax*nemax,lmd,1,dbl_mb(tmp_L(1)),1)
           call dscal(2*nemax*nemax,2.0d0,dbl_mb(tmp_L(1)),1)
           call daxpy(2*nemax*nemax,-1.0d0,lmd1,1,dbl_mb(tmp_L(1)),1)
           call paw_force_constraint(dbl_mb(tmp_L(1)),dbl_mb(fion(1))) 




*          **** remove ion forces using ion_FixIon ****
           call ion_FixIon(dbl_mb(fion(1)))

           call ion_optimize_step(dbl_mb(fion(1)))
        end if


*       *****************************************        
*       **** lagrange multiplier corrections ****
*       *****************************************        
        if (gram_schmidt) then
          if (move) call phafac2()
          do ms=1,ispin
            call paw_psi_MakeOrtho(npack1,ne(ms),psi2(1,n1(ms)))
          end do
        else

          call paw_ovlp_S(n2(ispin),psi1,psi_r)

          if (move) call phafac2()
          call dcopy(2*nemax*nemax,lmd,1,lmd1,1)
          call paw_psi_lmbda(ispin,ne,nemax,npack1,psi_r,psi2,dte,
     >                 lmd,dbl_mb(tmp_L(1)),ierr)
        end if

!*        ***** debug ***
!         do ms=1,ispin
!           write(23,*)
!           write(23,*)
!           write(23,*) "DEBUG: iteration=",it
!           write(23,*) "DEBUG: overlap matrix"
!           call paw_overlap_matrix_gen(ne(ms),ne(ms),
!     >                                 psi2(1,n1(ms)),
!     >                                 psi2(1,n1(ms)),
!     >                                 dbl_mb(tmp_L(1)))
!           write(23,*) "Overlap matrix, spin:",ms
!           do i=1,ne(ms) 
!             write(23,*) (dbl_mb(tmp_L(1)+(i-1) +(j-1)*ne(ms)),
!     >                   j=1,ne(ms))
!           end do
!         end do
!*        ***** debug ***


      end do

*     *************************************
*     ***** total energy calculation ******
*     *************************************
      call nwpw_timing_start(10)

      !if (move) call phafac() !*** reset phase factors to r1 ***

*     *** get orbital energies ****
      n = ne(1)
      nn = n*n
      do ms=1,ispin
         do ii=n1(ms),n2(ms)
           i = ii-n1(ms)
           index = (i+1) + i*n + (ms-1)*nn
           call Pack_cc_idot(1,psi1(1,ii),Hpsi(1,ii),sum)

           hml(index) =  -sum
           do jj=ii+1,n2(ms)
              j = jj-n1(ms)
              index  = (i+1) + j*n + (ms-1)*nn
              indext = (j+1) + i*n + (ms-1)*nn
              call Pack_cc_idot(1,psi1(1,ii),Hpsi(1,jj),sum)

              hml(index)  =  -sum
              hml(indext) =  -sum
           end do 
         end do
      end do
      if (np.gt.1)  call D3dB_Vector_SumAll((ispin*nn),hml)
      eorbit = 0.0d0
      do ms=1,ispin
         do ii=1,ne(ms)
            index = (ii) + (ii-1)*n + (ms-1)*nn
            eorbit = eorbit + hml(index)
         end do
      end do
      if (ispin.eq.1) eorbit = eorbit+eorbit
        

    
*     **** get coulomb energy ****
      call Pack_cc_Sum(0,
     >                   dcpl_mb(dng(1)),
     >                   dn_cmp_smooth,
     >                   dcpl_mb(tmp1(1)))
      call Pack_c_Copy(0,
     >                   dcpl_mb(vcomp(1)),
     >                   dcpl_mb(tmp2(1)))
      call Pack_cc_daxpy(0,0.5d0,
     >                   dcpl_mb(vh(1)),
     >                   dcpl_mb(tmp2(1)))
      call Pack_cc_dot(0,
     >                  dcpl_mb(tmp1(1)),
     >                  dcpl_mb(tmp2(1)),
     >                  ehartr_pw)
      ehartr_pw = ehartr_pw*lattice_omega()
     

        

*     **** get exchange-correlation energy ****
      call D3dB_rr_dot(1,dn(1,1),dbl_mb(xce(1)),exc)
      call D3dB_rr_dot(1,dn(1,1),dbl_mb(xcp(1)),pxc)
      if (ispin.eq.1) then
         exc= exc + exc 
         pxc= pxc + pxc 
      else
         call D3dB_rr_dot(1,dn(1,2),dbl_mb(xce(1)),exc2)
         call D3dB_rr_dot(1,dn(1,2),dbl_mb(xcp(1)+n2ft3d),pxc2)
         exc= exc + exc2
         pxc= pxc + pxc2
      end if
      exc = exc*dv
      pxc = pxc*dv
        
  
*     ***** average Kohn-Sham kinetic energy ****
      call ke_ave(ispin,ne,psi1,eke,.false.,occ)


*     **** average Kohn-Sham v_local energy ****
      call Pack_cc_dot(0,dcpl_mb(dng(1)),dcpl_mb(vl(1)),evloc_pw)



 
*     ***** average Kohn-Sham v_nonlocal energy ****
c     call dcopy(2*npack1*nemax,0.0d0,0,Hpsi,1)
c     call v_nonlocal(ispin,ne,psi1,Hpsi,
c    >                move,dbl_mb(ftest(1)))
      enlocal = 0.0d0
c     do ms=1,ispin
c     do n=n1(ms),n2(ms)
c        call Pack_cc_idot(1,psi1(1,n),Hpsi(1,n),sum)
c        enlocal = enlocal - sum
c     end do
c     end do
c     if (np.gt.1) call D3dB_SumAll(enlocal)
c     if (ispin.eq.1) enlocal = 2.0d0*enlocal

*     **** atomic energies ***
      ehartree_atom = paw_energy_hartree_atom()
      ekin_atom = paw_energy_kin_atom()
      evloc_atom = paw_energy_vloc_atom()
      eion_atom = paw_energy_ion_atom()
      ecore_atom = paw_energy_core_atom()
      ecore_ion_atom = paw_energy_core_ion_atom()
      ecore_self_atom = paw_energy_core_self_atom()
      exc_atom = paw_energy_xc_atom()

*?????????????????????? what is this ??????????????      
      call Pack_c_unpack(0,dn_cmp_smooth)
      call D3dB_cr_fft3b(1,dn_cmp_smooth)
      call D3dB_r_Zero_Ends(1,dn_cmp_smooth)


*     *** fill in total energy array ***

*     *** kinetic energy 
      E(2) = eke
      E(3) = ekin_atom 
      E(4) = ehartr_pw

*     *** coulomb contributions 
      E(5) = eion_atom + ecore_atom + ehartree_atom +
     >       ecore_ion_atom + ecore_self_atom +
     >       paw_mult_energy_atom_self() +
     >       paw_mult_energy_atom_comp() 

      E(6)=paw_mult_energy_atom_mult()
      
*     *** exch-correlation
      E(7) = exc
      E(8) = exc_atom
      
*     *** local pseudopot ***
      E(9)  = evloc_pw 
      E(10) = evloc_atom


*     *** total energy ***
      Eold=E(1)
      E(1) = 0.0d0
      do i=2,10
       E(1) = E(1) + E(i)
      end do

      E(11) = eorbit      

*     **** set convergence variables ****
      deltae = (E(1)-Eold)/(dt*dble(it_in))

*     *** deltac ***
      do n=n1(1),n2(ispin)
         do i=1,npack1
            Hpsi(i,n) = psi2(i,n) - psi1(i,n)
         end do
      end do

      do n=n1(1),n2(ispin)
         call Pack_cc_idot(1,Hpsi(1,n),Hpsi(1,n),dbl_mb(sumi(1)+n-1))
      end do
      if (np.gt.1) 
     >     call D3dB_Vector_SumAll((n2(ispin)-n1(1)),
     >                             dbl_mb(sumi(1)))
      deltac = 0.0d0
      do n=n1(1),n2(ispin)
         if (dbl_mb(sumi(1)+n-1).gt.deltac) deltac=dbl_mb(sumi(1)+n-1)
      end do
      deltac = deltac/dte


*     *** deltar ***
      deltar = deltamm
      if (move) then
        do i=1,ion_nion()
           sum = dsqrt( dbl_mb(fion(1)+(i-1)*3  )**2 
     >                + dbl_mb(fion(1)+(i-1)*3+1)**2 
     >                + dbl_mb(fion(1)+(i-1)*3+2)**2)
           if (sum.gt.deltar) deltar = sum
        end do
      end if
       
      call nwpw_timing_end(10)

*     **** dealocate MA local variables ****
      call nwpw_timing_start(12)
      value = BA_pop_stack(sumi(2))
      value = BA_pop_stack(ftest(2))
      value = BA_pop_stack(fion(2))
      value = BA_pop_stack(dnall(2))
      value = BA_pop_stack(xce(2))
      value = BA_pop_stack(xcp(2))
      value = BA_pop_stack(dng(2))
      value = BA_pop_stack(rho(2))
      value = BA_pop_stack(vl(2))


      value = BA_pop_stack(vc(2))
      value = BA_pop_stack(vh(2))
      value = BA_pop_stack(vcomp(2))
      value = BA_pop_stack(tmp2(2))
      value = BA_pop_stack(tmp1(2))
      value = BA_pop_stack(tmp_L(2))

      call nwpw_timing_end(12)



      return
      end
 
