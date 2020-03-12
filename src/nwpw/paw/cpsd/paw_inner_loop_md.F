*
* $Id$
*

* $Log: not supported by cvs2svn $
* Revision 1.14  2008/10/22 23:56:40  bylaska
* added NWCHEM_NWPW_LIBRARY to nwchemrc. fixed bug in paw...EJB
*
* Revision 1.13  2008/09/15 20:52:32  bylaska
* ..paw fixes...EJB
*
* Revision 1.12  2007/03/22 20:46:20  bylaska
* New implementation of QM/MM.
* ....EJB
*
* Revision 1.11  2006/02/11 02:50:46  bylaska
* GGA's using 1st derivative formulas have been added in core part of PAW....EJB
*
* Revision 1.10  2005/02/09 02:38:57  bylaska
* ..............EJB
*
* Revision 1.9  2005/01/31 16:38:45  bylaska
* ...EJB
*
* Revision 1.8  2005/01/31 03:58:32  bylaska
* aperiodic bug fixes in PAW-CPMD ... EJB
*
* Revision 1.7  2005/01/31 03:28:37  bylaska
* Thermostat bug fix in PAW-CPMD...EJB
*
* Revision 1.6  2005/01/25 22:32:04  bylaska
* Various bug fixes for release, containing band-DOS and input fixes.
* Also, center of rotation constrained added to pspw-cpmd and paw-cpmd (not in release).
* ...EJB
*
* Revision 1.5  2004/11/08 23:37:41  bylaska
* Bug fix in pspw_hfx found by M. Hackler.
* PBE0 has been implemented.
*
*  ........EJB
*
* Revision 1.4  2004/03/08 22:51:25  bylaska
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
* Revision 1.3  2003/10/21 02:05:15  marat
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
* Revision 1.2  2003/09/26 16:37:32  bylaska
* bug fix...EJB
*
* Revision 1.1  2003/09/25 23:46:24  bylaska
* PAW Car-Parrinello added...EJB
*
*

      subroutine paw_inner_loop_md(verlet,sa_alpha,
     >                      ispin,ne,
     >                      npack1,nfft3d,nemax,
     >                      psi0,psi1,psi2,
     >                      dn,dn_cmp_smooth,
     >                      it_in,it_sum,E,
     >                      hml,lmd,lmd1,
     >                      Hpsi,psi_r)
      implicit none
      logical    verlet
      real*8     sa_alpha(2)
      integer    ispin,ne(2)
      integer    npack1,nfft3d,nemax
      complex*16 psi0(npack1,nemax)
      complex*16 psi1(npack1,nemax)
      complex*16 psi2(npack1,nemax)
      real*8     dn(2*nfft3d,2)
      real*8     dn_cmp_smooth(2*nfft3d)
      integer    it_in,it_sum
      real*8     E(*)
      real*8     hml(2*nemax*nemax)
      real*8     lmd(2*nemax*nemax),lmd1(2*nemax*nemax)

*     **** very big workspace variables ****
      complex*16 Hpsi(npack1,nemax)
      real*8     psi_r(2*nfft3d,nemax)


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
      logical move,allow_translation
      logical use_lda
      integer n2ft3d,np
      integer i,j,ii,jj,n,n1(2),n2(2),it,ms,nn,ierr
      integer nx,ny,nz
      integer index,indext
      double precision evloc_pw,evloc_atom,occ(1)
      real*8  sum,Eold,eorbit,ehartr_pw,enlocal
      real*8  eke,eki,sse,ssr,sa1,sa2
      real*8  exc,exc2,pxc,pxc2,dte,dte0,scal1,scal2,dv,dt,h
      real*8  deltamm,vzero,fmass
      double precision ekin_atom,ekin_pw
      double precision eion_atom
      double precision ecore_atom
      double precision ecore_ion_atom
      double precision ecore_self_atom
      double precision ehartree_atom
      double precision exc_atom
      double precision cpu1,cpu2


*     **** MA local variables ****
      logical value,gram_schmidt,nose,field_exist
      integer tmp_L(2)
      integer tmp1(2),tmp2(2)
      integer vl(2),vh(2),vc(2),vcomp(2),dng(2)
      integer rho(2),r_grid(2)
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
      real*8   ion_ion_e,ion_ke
      external ion_ion_e,ion_ke
      real*8   paw_mult_energy_atom_comp !**no header file for paw_mult**
      real*8   paw_mult_energy_atom_self
      real*8   paw_mult_energy_atom_mult
      external paw_mult_energy_atom_comp
      external paw_mult_energy_atom_self
      external paw_mult_energy_atom_mult

      logical  control_Nose,control_allow_translation
      external control_Nose,control_allow_translation
      real*8   Nose_e_energy,Nose_r_energy,Nose_sse,Nose_ssr
      external Nose_e_energy,Nose_r_energy,Nose_sse,Nose_ssr

*     ***** QM/MM external functions ****
      logical  dplot_iteration_check
      external dplot_iteration_check




      call Parallel_np(np)
      call Pack_npack(0,npack0)
      n2ft3d = 2*nfft3d
      deltamm = 0.0d0
      gga = control_gga()

      !field_exist = pspw_charge_found()
c      allow_translation = control_allow_translation()



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

      if (.not. value) 
     > call errquit('paw_inner_loop_md:out of stack memory',0,0)

      call nwpw_timing_end(12)
     
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      move = .true.

      nose = control_Nose()
      sse = 1.0d0
      ssr = 1.0d0

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

      dt    = control_time_step()
      fmass = control_fake_mass()
      dte   = dt*dt/fmass
      if (.not. verlet) dte=0.5d0*dte
      if (.not.nose) then
        sa1 =    1.0d0/(2.0d0-sa_alpha(1))
        sa2 = sa_alpha(1)/(2.0d0-sa_alpha(1))
      end if

      h = 1.0d0/(2.0d0*dt)

      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = scal1*lattice_omega()



*     ******************************************
*     ****                                  ****
*     **** Start of molecular dynamics loop ****
*     ****                                  ****
*     ******************************************
      do it=1,it_in

*       **** shift wavefunction and atoms ****

        call dcopy(2*npack1*nemax,psi1,1,psi0,1)
        call dcopy(2*npack1*nemax,psi2,1,psi1,1)

*       *** skip ion_shift if newton step ***
        if (verlet) call ion_shift()
        if (nose)   call Nose_Shift()


*       ********************************
*       **** generate phaze factors ****
*       ********************************
        call phafac()
        call paw_set_mult_energy_coeff()

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

*       ********************
*       **** call dplot ****
*       ********************
        if (dplot_iteration_check(it+it_sum)) then
         call dplot_iteration((it+it_sum),ispin,ne,psi1,dn,psi_r)
        end if


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
       call paw_pot_hartree_solve()
        
*      *** atomic portion ***
           
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

!*       **********************************************
!*       **** generate other real-space potentials ****
!*       **********************************************
!        if (field_exist) then
!
!           call dcopy(n2ft3d,0.0d0,0,dbl_mb(v_field(1)),1)
!
!*          **** generate QM/MM potential ****
!           if (pspw_qmmm_found()) then
!            call pspw_qmmm_Generate_V(n2ft3d,
!     >                               dbl_mb(r_grid(1)),
!     >                               dbl_mb(v_field(1)))
!           end if
!
!*          **** generate charge potential ****
!           if (pspw_charge_found()) then
!            call pspw_charge_Generate_V(n2ft3d,
!     >                               dbl_mb(r_grid(1)),
!     >                               dbl_mb(v_field(1)))
!           end if
!
!        end if



     

*       ******************
*       **** get Hpsi ****
*       ******************
        call nwpw_timing_start(13)
        call paw_psi_H(ispin,ne,psi1,psi_r,
     >             dcpl_mb(vl(1)),
     >             dcpl_mb(vc(1)),dbl_mb(xcp(1)),Hpsi,
     >             move,dbl_mb(fion(1)))
        call nwpw_timing_end(13)

        

*       *******************************************
*       **** get ion forces and do update      ****
*       ****  on ions  and psi                 ****
*       *******************************************

*       *********************
*       **** generate force *
*       *********************

        call paw_mult_pw_force(dcpl_mb(vh(1)),
     >                         dcpl_mb(vcomp(1)),
     >                         dbl_mb(fion(1)))

        call paw_force_solve(psi1,dbl_mb(fion(1))) 


*       *** compute hamiltonian matrix if first iteration ****
        if (.not.verlet) then
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
        end if

        call dcopy(2*nemax*nemax,lmd,1,dbl_mb(tmp_L(1)),1)
        call dscal(2*nemax*nemax,2.0d0,dbl_mb(tmp_L(1)),1)
        call daxpy(2*nemax*nemax,-1.0d0,lmd1,1,dbl_mb(tmp_L(1)),1)
        call paw_force_constraint(dbl_mb(tmp_L(1)),dbl_mb(fion(1))) 


*       **** remove ion forces using ion_FixIon ****
        call ion_FixIon(dbl_mb(fion(1)))

c        !**** center of mass constraint ****
c        if (.not.allow_translation) then
c          call remove_center_F_mass(dbl_mb(fion(1)))
c        end if



*       *********************************
*       **** do a verlet step        ****
*       *********************************
        if (verlet) then

*          **** constant temperature ****
           if (nose) then
             sse = Nose_sse()
             ssr = Nose_ssr()
             do n=1,n2(ispin)
              call Pack_c_SMul(1,0.5d0*dte,Hpsi(1,n),psi2(1,n))
              call Pack_cc_daxpy(1,-1.0d0,psi0(1,n),psi2(1,n))
              call Pack_cc_daxpy(1,1.0d0,psi1(1,n),psi2(1,n))
c              call Pack_c_SMul(1,2.0d0*sse,psi2(1,n),psi2(1,n))
              call Pack_c_SMul1(1,2.0d0*sse,psi2(1,n))
              call Pack_cc_daxpy(1,1.0d0,psi0(1,n),psi2(1,n))
             end do
             call ion_nose_step(ssr,dbl_mb(fion(1)))

*          **** constant energy ****
           else
             do n=1,n2(ispin)
              call Pack_c_SMul(1,dte*sa1,Hpsi(1,n),psi2(1,n))
              call Pack_cc_daxpy(1,-1.0d0*sa2,psi0(1,n),psi2(1,n))
              call Pack_cc_daxpy(1,2.0d0*sa1,psi1(1,n),psi2(1,n))
             end do
             call ion_verlet_step(dbl_mb(fion(1)),sa_alpha(2))

           end if

*       *********************************
*       **** do a newton step        ****
*       *********************************
        else
           do n=1,n2(ispin)
              call Pack_c_SMul(1,dte,Hpsi(1,n),psi2(1,n))
              call Pack_cc_daxpy(1,dt*sa_alpha(1),psi0(1,n),psi2(1,n))
c              call Pack_cc_Sum(1,psi2(1,n),psi1(1,n),psi2(1,n))
              call Pack_cc_Sum2(1,psi1(1,n),psi2(1,n))
           end do
           call ion_newton_step(dbl_mb(fion(1)),sa_alpha(2))

        end if


*       *****************************************        
*       **** lagrange multiplier corrections ****
*       *****************************************        
        call paw_ovlp_S(n2(ispin),psi1,psi_r)

        call phafac2()
        call dcopy(2*nemax*nemax,lmd,1,lmd1,1)
        dte0 = dte
        if (nose) dte0 = dte*sse
        call paw_psi_lmbda(ispin,ne,nemax,npack1,psi_r,psi2,dte0,
     >                 lmd,dbl_mb(tmp_L(1)),ierr)

*       **************************
*       *** update thermostats ***
*       **************************
        if (nose) then
          if (verlet) then
             eke = 0.0d0
             do i=1,n2(ispin)
c                call Pack_c_SMul(1,-h,psi0(1,i),psi0(1,i))
                call Pack_c_SMul1(1,-h,psi0(1,i))
                call Pack_cc_daxpy(1,h,psi2(1,i),psi0(1,i))
                call Pack_cc_idot(1,psi0(1,i),psi0(1,i),sum)
                eke = eke+sum
             end do
             if (np.gt.1) call D3dB_SumAll(eke)
             eke = eke*fmass
             if (ispin.eq.1) eke = 2.0d0*eke
             eki = ion_ke()
             call Nose_Verlet_Step(eke,eki)
          else
              eke = 0.0d0
              do i=1,n2(ispin)
                call Pack_cc_idot(1,psi0(1,i),psi0(1,i),sum)
                eke = eke+sum
              end do
              if (np.gt.1) call D3dB_SumAll(eke)
              if (ispin.eq.1) eke = 2.0d0*eke
              eke = eke*fmass
              eki = ion_ke()
              call Nose_Newton_Step(eke,eki)
          end if
        end if


      end do
*     ******************************************************
*     ***** end main loop **********************************
*     ******************************************************

*     **** if newton then skip energy calculations ****
      if (.not. verlet) goto 333


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
     >                 dcpl_mb(dng(1)),
     >                 dn_cmp_smooth,
     >                 dcpl_mb(tmp1(1)))
      call Pack_c_Copy(0,
     >                 dcpl_mb(vcomp(1)),
     >                 dcpl_mb(tmp2(1)))
      call Pack_cc_daxpy(0,0.5d0,
     >                 dcpl_mb(vh(1)),
     >                 dcpl_mb(tmp2(1)))
      call Pack_cc_dot(0,
     >                dcpl_mb(tmp1(1)),
     >                dcpl_mb(tmp2(1)),
     >                ehartr_pw)
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
      call ke_ave(ispin,ne,psi1,ekin_pw,.false.,occ)


*     **** average Kohn-Sham v_local energy ****
      call Pack_cc_dot(0,dcpl_mb(dng(1)),dcpl_mb(vl(1)),evloc_pw)

 
*     ***** average Kohn-Sham v_nonlocal energy ****
      enlocal = 0.0d0

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



*     **** velocity an kinetic energy of psi ****
      if (.not.nose) then
c        h = 1.0d0/(2.0d0*dt)
         eke = 0.0d0
         do i=1,n2(ispin)
c            call Pack_c_SMul(1,-h,psi0(1,i),psi0(1,i))
            call Pack_c_SMul1(1,-h,psi0(1,i))
            call Pack_cc_daxpy(1,h,psi2(1,i),psi0(1,i))
            call Pack_cc_idot(1,psi0(1,i),psi0(1,i),sum)
            eke = eke+sum
         end do
         if (np.gt.1) call D3dB_SumAll(eke)
         eke = eke*fmass
         if (ispin.eq.1) eke = 2.0d0*eke
      end if


*     *** fill in total energy array ***

*     **** total energy ****
      Eold=E(1)
      E(2) = ekin_pw + ekin_atom
     >     + ehartr_pw 
     >     + eion_atom + ecore_atom + ehartree_atom 
     >     + ecore_ion_atom + ecore_self_atom 
     >     + paw_mult_energy_atom_self() 
     >     + paw_mult_energy_atom_comp()
     >     + paw_mult_energy_atom_mult()
     >     + exc + exc_atom
     >     + evloc_pw + evloc_atom
      E(3) = eke
      E(4) = ion_ke()

      E(5) = eorbit
      E(6) = ehartr_pw 
     >     + eion_atom + ecore_atom + ehartree_atom
     >     + ecore_ion_atom + ecore_self_atom
     >     + paw_mult_energy_atom_self()
     >     + paw_mult_energy_atom_comp()
      E(7) = exc+exc_atom
      E(8) = 0.0d0



*     **** Energy and Energy**2 sum ***
      E(25) = E(25) + E(2)
      E(26) = E(26) + E(2)*E(2)

*     **** output Forces for Fei ***
      call fei_output(E(2),dbl_mb(fion(1)))


      if (nose) then
        E(9) = Nose_e_energy()
        E(10) = Nose_r_energy()
        E(1) = E(2)+E(3)+E(4)+E(9)+E(10)
      else
        E(1) = E(2)+E(3)+E(4)
      end if

      call nwpw_timing_end(10)

*     **** dealocate MA local variables ****
 333  continue
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
 
