*
* $Id$
* $Log: not supported by cvs2svn $
* Revision 1.26  2008/09/15 20:52:31  bylaska
* ..paw fixes...EJB
*
* Revision 1.25  2007/11/17 22:45:48  bylaska
* ...EJB
*
* Revision 1.24  2007/04/02 15:38:39  bylaska
* Symbol names increased to length 4 instead of 2.
* standalone  rdf calculation added.
* lattice vectors added to ion_motion files
* changed the qmmm psp to have short range repulsion for positive mm charges
* ...EJB
*
* Revision 1.23  2007/03/22 20:46:20  bylaska
* New implementation of QM/MM.
* ....EJB
*
* Revision 1.22  2006/08/13 01:03:26  bylaska
* Checking in code not include in 5.0 release.
* A chain algorithm was added to Nose-Hoover thermostats.
* Preliminary implementation of a processor group decomposition added to pspw, i.e. parallel decomposition is over fft grid and electrons.
* ...EJB
*
* Revision 1.20  2006/02/11 02:50:46  bylaska
* GGA's using 1st derivative formulas have been added in core part of PAW....EJB
*
* Revision 1.19  2006/01/07 00:55:33  marat
* removing ma_auto_verify
* added some comments to paw
*
* Revision 1.18  2005/12/22 01:35:05  bylaska
* revPBE added and gga logic restructured....EJB
*
* Revision 1.17  2005/07/09 22:44:22  bylaska
* adding Louie FFT....EJB
* flag added for PAW xc and comp angular integration.
* ....EJB
*
* Revision 1.16  2005/02/09 02:38:57  bylaska
* ..............EJB
*
* Revision 1.15  2005/01/31 23:15:16  bylaska
* ...EJB
*
* Revision 1.14  2005/01/31 16:38:45  bylaska
* ...EJB
*
* Revision 1.13  2005/01/31 03:58:32  bylaska
* aperiodic bug fixes in PAW-CPMD ... EJB
*
* Revision 1.12  2005/01/31 02:09:12  bylaska
* Large nested loops in PAW  parallelized using an openMP-like strategy.
* ...EJB
*
* Revision 1.11  2004/09/04 17:56:21  bylaska
* Added local potential to the projector file (.jpp).
* More updates to constraint force.
* ...EJB
*
* Revision 1.10  2004/08/01 02:02:45  bylaska
* updates...EJB
*
* Revision 1.9  2004/07/29 15:55:31  bylaska
* Temporary array for Gaunt coefficients added.  Speeds up the program considerably, but it is extrememely memory intensive.  Also added timing routines to multipole calculations and fixed the initial total charge calculation.
*
*  ...EJB
*
* Revision 1.8  2004/05/16 05:38:07  edo
* *** empty log message ***
*
* Revision 1.7  2004/05/05 19:33:50  bylaska
* hilbert mapping added....EJB
*
* Revision 1.6  2004/03/15 15:17:41  bylaska
* Restructuring for 2d Hilbert mapping.
*
* D3dB_Init(nb,n1,n2,n2) changed to D3dB_Init(nb,n1,n2,n3,map)
*
* Calls to D3dB_ktoqp changed to either D3dB_ijktoindexp or D3dB_ijktoindex2p
* depending the context of the call.
*
* Calls to D3dB_qtok have been eliminated from coulomb2.F and the code
* has been restructured to used D3dB_ijktoindexp and D3dB_ijktoindex2p
*
*
* ...EJB
*
* Revision 1.5  2004/02/10 04:27:25  edo
* spellinggggggggggggg
*
* Revision 1.4  2003/10/21 02:05:15  marat
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
* Revision 1.3  2003/09/26 16:31:16  bylaska
* bug fix....EJB
*
* Revision 1.2  2003/09/26 01:16:48  bylaska
* io bug fix...EJB
*
* Revision 1.1  2003/09/25 23:46:23  bylaska
* PAW Car-Parrinello added...EJB
*


***********************************************************************
*                      paw_cpmd                                       *
*                                                                     *
*     This is a developing PAW Car-Parrinello MD code for NWChem.     *
*                                                                     *
*                                                                     * 
*  Authors: Marat Valiev and Eric J. Bylaska                          *
*                                                                     *
***********************************************************************

      logical function paw_cpmd(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "stdio.fh"
      
      logical value

      real*8 kb
      parameter (kb=3.16679d-6)


*     **** parallel variables ****
      integer  taskid,np,np_i,np_j
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,nfft3d,n2ft3d
      integer npack1

*     **** electronic variables ****
      logical first_iteration,psi_nogrid
      integer ispin
      integer ne(2),n1(2),n2(2),nemax,neall,neq(2),nemaxq
      real*8  r_charge,icharge
      real*8  dipole(3),occ(1)

      integer psi0(2),psi1(2),psi2(2)
      integer dn(2),dn_cmp_smooth(2)
      integer Hpsi(2),psir(2)
    

*     ***** energy variables ****
      real*8  E(30),eke,eave,evar,cv

*     real*8  eig(2*nemax)
*     real*8  hml(2*nemax*nemax)
*     real*8  lmd(2*nemax*nemax)
      integer eig(2),hml(2),lmd(2),lmd1(2)
      real*8 Te_init,Tr_init,Te_new,Tr_new,sa_decay(2),sa_alpha(2)




*     **** error variables ****
      integer ierr

*     **** local variables ****
      logical verlet,mulliken,SA,found
      integer ms,lmax,idum
      real*8  deltae,deltac,deltar,dum(1)
      real*8  gx,gy,gz,cx,cy,cz,sum1,sum2
      real*8 vcx,vcy,vcz,vgx,vgy,vgz
      real*8  ekg,eki0,eki1,sum
      real*8  eke0,eke1
      real*8  EV,pi,dt
      real*8  emotion_time_shift
      integer i,j,k,ia,n,nn
      integer ii,jj,indx
      integer icount,it_in,it_out,icount_shift
      real*8 w,sumall,virial
      integer nfft3,mapping,mapping1d
      parameter (nfft3=32)
      character*255 full_filename
      character*50 filename

      integer hversion,hnfft(3),hispin,hne(2)
      real*8 hunita(3,3)
      integer ind

      character*50 control_input_psi
      external     control_input_psi
      logical  wvfnc_expander
      external wvfnc_expander

  


*     **** external functions ****
      real*8      paw_mult_rcut,ion_amass
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg
      integer     paw_mult_ncut
      character   spdf_name
      character*4 ion_aname,ion_atom
      external    paw_mult_rcut,ion_amass
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg
      external    paw_mult_ncut
      external    spdf_name
      external    ion_aname,ion_atom



      real*8   control_rti,control_rte
      real*8   ion_vion,ion_com_ke,ion_ke
      real*8   ion_Temperature,ion_com_Temperature
      external control_rti,control_rte
      external ion_vion,ion_com_ke,ion_ke
      external ion_Temperature,ion_com_Temperature

      real*8   control_tole,control_tolc,control_tolr,ion_rion
      external control_tole,control_tolc,control_tolr,ion_rion
      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      logical  control_read,ion_init,ion_q_FixIon
      external control_read,ion_init,ion_q_FixIon

      integer  pack_nwave_all
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,control_lmax_multipole
      integer  ion_nion,ion_natm,ion_katm,ion_nkatm
      external pack_nwave_all
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,control_lmax_multipole
      external ion_nion,ion_natm,ion_katm,ion_nkatm

      character*12 control_boundry
      external     control_boundry

      logical      control_Mulliken
      external     control_Mulliken
      logical      psi_filefind
      external     psi_filefind
      real*8       nwpw_timing,paw_psi_CheckOrtho
      external     nwpw_timing,paw_psi_CheckOrtho

      logical  control_Nose,Nose_restart
      external control_Nose,Nose_restart

      integer  Nose_Mchain,Nose_Nchain
      external Nose_Mchain,Nose_Nchain
 
      real*8   control_Nose_Te,Nose_Qe,Nose_Pe,Nose_Ee0
      external control_Nose_Te,Nose_Qe,Nose_Pe,Nose_Ee0

      real*8   control_Nose_Tr,Nose_Qr,Nose_Pr,Nose_Er0
      external control_Nose_Tr,Nose_Qr,Nose_Pr,Nose_Er0
      logical      v_psi_filefind
      external     v_psi_filefind

      logical  control_out_of_time,control_new_vpsi
      external control_out_of_time,control_new_vpsi


      logical  control_SA,control_Fei
      real*8   control_SA_decay
      external control_SA,control_Fei
      external control_SA_decay
      integer  control_np_orbital,control_mapping,control_mapping1d
      external control_np_orbital,control_mapping,control_mapping1d

      logical  control_translation,control_rotation
      external control_translation,control_rotation



*                            |************|
*****************************|  PROLOGUE  |****************************
*                            |************|

      value = .true.
      pi = 4.0d0*datan(1.0d0)

      call nwpw_timing_init()
      call dcopy(30,0.0d0,0,E,1)


*     **** get parallel variables ****
      call Parallel_Init()
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) call current_second(cpu1)

*     ***** print out header ****
      if (taskid.eq.MASTER) then
         write(luout,1000)
         write(luout,1010)
         write(luout,1020)
         write(luout,1010)
         write(luout,1030)
         write(luout,1010)
         write(luout,1035)
         write(luout,1036)
         write(luout,1010)
         write(luout,1040)
         write(luout,1010)
         write(luout,1041)
         write(luout,1010)
         write(luout,1042)
         write(luout,1043)
         write(luout,1010)
         write(luout,1000)
         call nwpw_message(1)
         write(luout,1110)
      end if
      
      value = control_read(7,rtdb)
      call Parallel2d_Init(control_np_orbital())
      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)

      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0
      mapping = control_mapping()


*     **** initialize D3dB data structure ****
      call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d

*     ***** Initialize double D3dB data structure ****
      if (control_version().eq.4) 
     >   call D3dB_Init(2,2*ngrid(1),2*ngrid(2),2*ngrid(3),mapping)


*     **** initialize lattice and packing data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_init()
      call Pack_npack(1,npack1)      

      call D3dB_pfft_init()


*     **** initialize Gaunt array ****
      call paw_gaunt_init()

 
*     **** read ions ****
      value = ion_init(rtdb)
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      first_iteration = .true.


*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)


*     **** allocate paw data structure and read in paw basis into it ****
c      value = BA_set_auto_verify(.true.)
      call init_paw_basis()

*     *** initialize paw matrices ***      
      call init_paw_kin_matrix()
      call init_paw_ion_matrix()
      call init_paw_vloc_matrix()
      call init_paw_core_matrix()
      call init_paw_hartree_matrix()
      call init_paw_overlap_matrix()
      call init_paw_comp_charge_matrix()
      call init_paw_comp_pot_matrix()

*     *** initialize paw projectors ***      
      call paw_proj_init()
c      call paw_vloc_init()
      
*     *** initialize compensation charge ***
      call paw_comp_charge_init()
      
*     *** initialize paw atomc potentials ***      
      call init_paw_pot_hartree()
      call init_paw_pot_comp()
      call paw_mult_init()


*     **** initialize G,mask,ke,and coulomb data structures ****
      call ke_init()
      if (control_version().eq.3) call coulomb_init()
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()



*     ***** allocate psi2,and psi1 wavefunctions ****
      call psi_get_ne(ispin,ne)
      mapping1d = control_mapping1d()
      call Dne_init(ispin,ne,mapping1d)
      call Dneall_neq(neq)
      nemaxq = neq(1)+neq(2)

      value = BA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'psi1',psi1(2),psi1(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'psi0',psi0(2),psi0(1))
      if (.not. value) call errquit('paw_cpmd:out of heap memory',0,0)


*     *****  read psi2 wavefunctions ****
      call psi_read(ispin,ne,dcpl_mb(psi2(1)),idum,dum)


*     **** move  wavefunction velocities ****
      if (control_new_vpsi()) then
        call v_psi_delete()
      end if

*    **** generate initial wavefunction velocities if it does not exist ****
      if (.not.v_psi_filefind()) then
        call v_psi_new(ispin,ne)
      end if

*     *****  read psi0 wavefunctions ****
      call v_psi_read(ispin,ne,dcpl_mb(psi1(1)))



      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)
      nemax = ne(1)+ne(2)


*     **** allocate other variables *****
      value = BA_alloc_get(mt_dbl,(2*nemax),'eig',eig(2),eig(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*nemax*nemax),'hml',hml(2),hml(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*nemax*nemax),'lmd',lmd(2),lmd(1))
      value = value.and.
     >     BA_alloc_get(mt_dbl,(2*nemax*nemax),'lmd1',lmd1(2),lmd1(1))
      call dcopy(2*nemax*nemax,0.0d0,0,dbl_mb(lmd(1)), 1)
      call dcopy(2*nemax*nemax,0.0d0,0,dbl_mb(lmd1(1)),1)

      value = value.and.
     >        BA_alloc_get(mt_dbl,(4*nfft3d),
     >                     'dn',dn(2),dn(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*nfft3d),
     >             'dn_cmp_smooth',dn_cmp_smooth(2),dn_cmp_smooth(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'Hpsi',Hpsi(2),Hpsi(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d*(ne(1)+ne(2)),
     >                     'psir',psir(2),psir(1))
      if (.not. value) call errquit('paw_cpmd: out of heap memory',0,1)

*     *** intialize overlap coefficient data structure ***
      call phafac()
      call paw_ovlp_init(ispin,ne)
      call paw_nonlocal_init(ispin,ne)
      lmax = control_lmax_multipole()
      if (lmax.lt.0) lmax =  paw_basis_max_mult_l()
c      call init_paw_density(ispin,lmax)
      call init_paw_xc(ispin,lmax)

*     *** intialize paw force ***      
      call paw_force_init()

      call paw_ovlp_coeff_set(dcpl_mb(psi2(1)))
      call paw_ovlp_weights_set()


      !**** Ortho Check ****
      do ms=1,ispin
        deltae=paw_psi_CheckOrtho(npack1,ne(ms),
     >                   dcpl_mb(psi2(1)+(n1(ms)-1)*npack1))

        if (deltae.gt.1.0d-10) then
          call paw_psi_MakeOrtho(npack1,ne(ms),
     >                   dcpl_mb(psi2(1)+(n1(ms)-1)*npack1))
          deltac=paw_psi_CheckOrtho(npack1,ne(ms),
     >                   dcpl_mb(psi2(1)+(n1(ms)-1)*npack1))
          if (taskid.eq.MASTER) then
            if (ms.eq.1) then
              write(luout,*) "Warning: ",
     >                   "Gram-Schmidt performed on up spin of psi2 "
              write(luout,*) "       : (old error=",deltae,
     >                   " new error=",deltac,")"
            end if
            if (ms.eq.2) then
              write(luout,*) "Warning: ",
     >                   "Gram-Schmidt performed on down spin of psi2 "
              write(luout,*) "       : (old error=",deltae,
     >                   " new error=",deltac,")"
            end if

          end if
        end if
     
      end do


*     ******************************
*     **** scaling psi velocity ****
*     ******************************
      call dcopy(2*(ne(1)+ne(2))*npack1,dcpl_mb(psi1(1)),1,
     >                                  dcpl_mb(psi0(1)),1)
      call dscal(2*(ne(1)+ne(2))*npack1,control_rte(),
     >           dcpl_mb(psi1(1)),1)
      eke0 = 0.0d0
      eke1 = 0.0d0
      do i=1,(ne(1)+ne(2))
         call Pack_cc_dot(1,dcpl_mb(psi0(1)+(i-1)*npack1),
     >                      dcpl_mb(psi0(1)+(i-1)*npack1),
     >                     sum)
         eke0 = eke0 + sum
         call Pack_cc_dot(1,dcpl_mb(psi1(1)+(i-1)*npack1),
     >                      dcpl_mb(psi1(1)+(i-1)*npack1),
     >                    sum)
         eke1 = eke1 + sum
      end do
      eke0 = control_fake_mass()*eke0
      eke1 = control_fake_mass()*eke1
      call ion_init_ke(ekg,eki0,eki1)

*     **** Initialize thermostats ****
      if (control_Nose()) then
         call ke_ave(ispin,ne,dcpl_mb(psi2(1)),w,.false.,occ)
         call Nose_Init((ne(1)+ne(2)),w)
      end if


*     **** Initialize simulated annealing ****
      SA=.false.
      Te_init     = 0.0d0
      Tr_init     = 0.0d0
      sa_alpha(1) = 1.0d0
      sa_alpha(2) = 1.0d0
      if (control_SA()) then
         if (control_Nose()) then
            SA          = .true.
            sa_decay(1) = control_SA_decay(1)
            sa_decay(2) = control_SA_decay(2)
            Te_init     = control_Nose_Te()
            Tr_init     = control_Nose_Tr()
         else
            dt = control_time_step()
            SA          = .false.
            sa_decay(1) = control_SA_decay(1)
            sa_decay(2) = control_SA_decay(2)
            sa_alpha(1) = dexp( -(dt/control_SA_decay(1)) )
            sa_alpha(2) = dexp( -(dt/control_SA_decay(2)) )
         end if
      end if


*     **** initialize QM/MM ****
c      call pspw_qmmm_init(rtdb)

*     **** initialize dplot ****
      call dplot_iteration_init()




*                |**************************|
******************   summary of input data  **********************
*                |**************************|

      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      call center_v_geom(vcx,vcy,vcz)
      call center_v_mass(vgx,vgy,vgz)
      mulliken = control_Mulliken()

      if (taskid.eq.MASTER) then
         write(luout,1111) np
         write(luout,1117) np_i,np_j
         if (mapping.eq.1) write(luout,1112)
         if (mapping.eq.2) write(luout,1113)
         if (mapping.eq.3) write(luout,1118)

         write(luout,1115)
         write(luout,1121) control_boundry(),control_version()
         if (ispin.eq.1) write(luout,1130) 'restricted'
         if (ispin.eq.2) write(luout,1130) 'unrestricted'
         IF (control_gga().eq.0) THEN
            write(luout,1131) 'Vosko et al parameterization'
         ELSE IF (control_gga().eq.10) THEN
            write(luout,1131) 
     >      'PBE96 (White and Bird) parameterization'
         ELSE IF (control_gga().eq.11) THEN
            write(luout,1131) 
     >      'BLYP (White and Bird) parameterization'
         ELSE IF (control_gga().eq.12) THEN
            write(luout,1131) 
     >      'revPBE (White and Bird) parameterization'
         ELSE
            write(luout,1131) 'unknown parameterization'
            call errquit('bad exchange_correlation',0,0)
         END IF

         write(luout,1140)
         do ia = 1,ion_nkatm()
            write(luout,1141) ia,ion_atom(ia),
     >                    paw_basis_ion_charge(ia),
     >                    paw_basis_core_charge(ia)
            write(luout,1143) paw_basis_sphere_radius(ia)
            !write(luout,1144) paw_basis_sigma(ia),paw_basis_sigma(ia)
            write(luout,1144) paw_basis_sigma(ia)
            write(luout,1150) paw_proj_nbasis(ia)
            write(luout,1151)
            do i=1,paw_basis_nbasis(ia)
              write(luout,1152) paw_basis_n_ps(i,ia),
     >                      paw_basis_n(i,ia),
     >                      spdf_name(paw_basis_orb_l(i,ia)),
     >                      paw_basis_eig(i,ia),
     >                      2*paw_basis_orb_l(i,ia)+1
            end do
         end do

         icharge = -(ne(1)+ne(ispin))
         do ia=1,ion_nkatm()
           icharge = icharge + ion_natm(ia)*
     >                        (paw_basis_ion_charge(ia)
     >                        -paw_basis_core_charge(ia))
         end do
         write(luout,1159) icharge

         write(luout,1160)
         write(luout,1170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
         write(luout,1180)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0
           else
           write(luout,1190) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz



         write(luout,1181)
         write(luout,1192) (I,ion_aname(I),
     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
         write(luout,1200) vcx,vcy,vcz
         write(luout,1210) vgx,vgy,vgz


         write(luout,1220) ne(1),ne(ispin),' ( fourier space)'
         write(luout,1230)
         write(luout,1241) lattice_unita(1,1),
     >                 lattice_unita(2,1),
     >                 lattice_unita(3,1)
         write(luout,1242) lattice_unita(1,2),
     >                 lattice_unita(2,2),
     >                 lattice_unita(3,2)
         write(luout,1243) lattice_unita(1,3),
     >                 lattice_unita(2,3),
     >                 lattice_unita(3,3)
         write(luout,1244) lattice_unitg(1,1),
     >                 lattice_unitg(2,1),
     >                 lattice_unitg(3,1)
         write(luout,1245) lattice_unitg(1,2),
     >                 lattice_unitg(2,2),
     >                 lattice_unitg(3,2)
         write(luout,1246) lattice_unitg(1,3),
     >                 lattice_unitg(2,3),
     >                 lattice_unitg(3,3)
         write(luout,1231) lattice_omega()
         write(luout,1250) lattice_ecut(),ngrid(1),ngrid(2),ngrid(3),
     >                 pack_nwave_all(0),pack_nwave(0)
         write(luout,1251) lattice_wcut(),ngrid(1),ngrid(2),ngrid(3),
     >                 pack_nwave_all(1),pack_nwave(1)
         write(luout,1260) paw_mult_rcut(),paw_mult_ncut()
         write(luout,1270)
         write(luout,1262) lmax
         if (.not.control_translation()) write(luout,1271)
         if (.not.control_rotation())    write(luout,1272)
         write(luout,1280) control_time_step(),control_fake_mass()
         write(luout,1290) control_rte(),control_rti()
         call ion_scaling_atoms_print(luout)
         write(luout,1222) eke0,eki0,ekg
         write(luout,1223) eke1,eki1
         write(luout,1224) (eke1-eke0),(eki1-eki0)
         if (control_Nose()) then
           write(luout,1295)
           if (Nose_restart()) then
              write(luout,*) "    thermostats resused"
           else
              write(luout,*) "    thermostats initialized"
           end if
           do i=1,Nose_Mchain()
             write(luout,1297) i,control_Nose_Te(),Nose_Qe(i),
     >                     Nose_Pe(i),Nose_Ee0(i)
           end do
           do i=1,Nose_Nchain()
             write(luout,1298) i,control_Nose_Tr(),Nose_Qr(i),
     >                     Nose_Pr(i),Nose_Er0(i)
           end do
         else
           write(luout,1294)
         end if
        if (control_SA()) then
           write(luout,1296) sa_decay(1),sa_decay(2)
         end if

         if (mulliken) write(luout,1299)
         write(luout,1300)
         write(luout,1305)
         call util_flush(luout)
      end if


*                |***************************|
******************     start iterations      **********************
*                |***************************|
*     **** open xyz and MOTION file ****
      call xyz_init()          ! unit=18
      call MOTION_init(rtdb)   ! unit=19

*     *** fei io ****
      call fei_init()

*     ************************************
*     **** open up other MOTION files ****
*     ************************************


*     **** open EMOTION file ****
      if (.not.btdb_cget(rtdb,'cpmd:emotion_filename',1,filename))
     >  call util_file_prefix('emotion',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)
      if (taskid.eq.MASTER) then

         emotion_time_shift = 0.0d0
         icount_shift       = 0
         inquire(file=full_filename,exist=found)
         if (found) then
           open(unit=31,file=full_filename,form='formatted',
     >          status='old')
           do while (found)
           read(31,*,end=100) emotion_time_shift,w,sum
           E(25) = E(25) + sum                          !*** take care of running sums ***
           E(26) = E(26) + sum*sum
           icount_shift = icount_shift + 1
           end do
  100      continue
#if defined(FUJITSU) || defined(__crayx1) || defined(GCC46)
           backspace 31
#endif
         else
           open(unit=31,file=full_filename,form='formatted',
     >          status='new')
         end if
      end if

*     **** open EIGMOTION file ****
      if (mulliken) then
      value = btdb_cget(rtdb,'cpmd:eigmotion_filename',1,filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)
      if (taskid.eq.MASTER)
     >   open(unit=32,file=full_filename,form='formatted')
      end if

*     **** open HMOTION file ****
      if (mulliken) then
      value = btdb_cget(rtdb,'cpmd:hmotion_filename',1,filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)
      if (taskid.eq.MASTER)
     >   open(unit=34,file=full_filename,form='formatted')
      end if

*     **** open OMOTION file ****
      if (mulliken) call Orb_Init(rtdb,ispin,ne) !unit=33

*     **** write initial position to xyz data ****
      call xyz_write()


*     ***** first step using velocity ****
      verlet = .false.
      call paw_inner_loop_md(verlet,sa_alpha,ispin,ne,
     >             npack1,nfft3d,nemax,
     >             dcpl_mb(psi0(1)),
     >             dcpl_mb(psi1(1)),
     >             dcpl_mb(psi2(1)),
     >             dbl_mb(dn(1)),dbl_mb(dn_cmp_smooth(1)),
     >             1,0,E,
     >             dbl_mb(hml(1)),dbl_mb(lmd(1)),dbl_mb(lmd1(1)),
     >             dcpl_mb(Hpsi(1)),dcpl_mb(psir(1)))


      if (taskid.eq.MASTER) call current_second(cpu2)
      if (taskid.eq.MASTER) CALL nwpw_message(6)
      it_in  = control_it_in()
      it_out = control_it_out()
      icount = 0
      verlet = .true.
      eke    = 0.0d0
      if (it_out.lt.1) goto 102



      dt = control_time_step()

      Te_new = Te_init
      Tr_new = Tr_init
  101 continue
         icount = icount + 1
         call paw_inner_loop_md(verlet,sa_alpha,ispin,ne,
     >             npack1,nfft3d,nemax,
     >             dcpl_mb(psi0(1)),
     >             dcpl_mb(psi1(1)),
     >             dcpl_mb(psi2(1)),
     >             dbl_mb(dn(1)),dbl_mb(dn_cmp_smooth(1)),
     >             it_in,((icount-1)*it_in),
     >             E,
     >             dbl_mb(hml(1)),dbl_mb(lmd(1)),dbl_mb(lmd1(1)),
     >             dcpl_mb(Hpsi(1)), dcpl_mb(psir(1)))

         eke = eke + E(3)

         if (taskid.eq.MASTER) then

           if (SA) then
             write(6,1309) icount*it_in,E(1),E(2),E(3),E(4),
     >                     Te_new,Tr_new
           else
             write(6,1310) icount*it_in,E(1),E(2),E(3),E(4),
     >                     ion_Temperature()
           end if
           call util_flush(6)

*          **** write out EMOTION data ****
           eave = E(25)/dble(icount+icount_shift)
           evar = E(26)/dble(icount+icount_shift)
           evar = evar - eave*eave
           if (control_Nose()) then
             write(31,1311) icount*it_in*dt + emotion_time_shift,
     >                    e(1),e(2),e(3),e(4),e(5),e(6),
     >                    e(7),e(8),e(9),e(10),
     >                    eave,evar,ion_Temperature()
           else
             write(31,1311) icount*it_in*dt + emotion_time_shift,
     >                    e(1),e(2),e(3),e(4),e(5),e(6),
     >                    e(7),e(8),
     >                    eave,evar,ion_Temperature()
           end if
           call util_flush(31)


*          **** write out EIGMOTION data -diagonal hml matrix ****
           if (mulliken) then
           write(32,1311) icount*it_in*dt,
     >       (( dbl_mb(hml(1)+ii-1+(ii-1)*ne(1)+(ms-1)*ne(1)*ne(1)),
     >         ii=1,ne(ms)),ms=1,ispin)
           call util_flush(32)
           end if

*          **** write out HMOTION data - hml matrix ****
           if (mulliken) then
           write(34,1312) icount*it_in*dt,ispin
           do ms=1,ispin
             write(34,1313) ms,ne(ms),ne(ms)
             do ii=1,ne(ms)
               write(34,1311)
     >         (dbl_mb(hml(1)+ii-1+(jj-1)*ne(1)+(ms-1)*ne(1)*ne(1)),
     >          jj=1,ne(ms))
             end do
           end do
           call util_flush(34)
           end if

         end if


*        **** write xyz data ****
         call xyz_write()
         call MOTION_write((icount*it_in*control_time_step()))

*        **** write OMOTION data ****
         if (mulliken) call Orb_Write(dcpl_mb(psi1(1)))

*        **** update thermostats using SA decay ****
         if (SA) then
           t1 = icount*it_in*dt/sa_decay(1)
           t2 = icount*it_in*dt/sa_decay(2)
           Te_new = Te_init*dexp(-t1)
           Tr_new = Tr_init*dexp(-t2)
           call Nose_reset_T(Te_new,Tr_new)
         end if


*        **** exit early ****
         if (control_out_of_time()) then
            if (taskid.eq.MASTER)
     >       write(6,*) ' *** out of time. iteration terminated'
            go to 102
         end if
      if (icount.lt.it_out) go to 101
      if (taskid.eq.MASTER)
     > write(6,*) '*** Arrived at the maximum iteration. Terminated.'

*::::::::::::::::::::  end of iteration loop  :::::::::::::::::::::::::

  102 continue


*     **** close xyz and MOTION files ****
      call xyz_end()
      call MOTION_end()
      if (taskid.eq.MASTER) then
        close(unit=31)
        close(unit=32)
        close(unit=34)
      end if

*     *** close fei io ****
      call fei_end()

*     **** close OMOTION file ****
      if (mulliken) call Orb_End()


      if (taskid.eq.MASTER) CALL nwpw_message(3)
      if (taskid.eq.MASTER) call current_second(cpu3)


*



*         |****************************************|
*********** produce CHECK file and diagonalize hml *****************
*         |****************************************|

*     **** produce CHECK FILE ****
      if (taskid.eq.MASTER) then
         call util_file_name('CHECK',.true.,
     >                               .false.,
     >                        full_filename)
         open(unit=17,file=full_filename,form='formatted')
      end if


*     **** check total number of electrons ****

      !*** psi charge ****
      r_charge = 0.0d0
      do ms =1,ispin
         call D3dB_r_dsum(1,dbl_mb(dn(1)+(ms-1)*n2ft3d),sumall)
         r_charge = r_charge
     >          + sumall*lattice_omega()
     >             /dble(ngrid(1)*ngrid(2)*ngrid(3))
      end do
      if (ispin.eq.1) r_charge = 2.0d0*r_charge

      !*** add comp charge ****
      call D3dB_r_dsum(1,dbl_mb(dn_cmp_smooth(1)),sum1)
      sum1 = sum1*lattice_omega()
     >      /dble(ngrid(1)*ngrid(2)*ngrid(3))
      r_charge = r_charge + sum1


      if (taskid.eq.MASTER) then
         write(17,1321) r_charge
      end if

*     **** comparison between hamiltonian an lambda matrix ****
      n = ne(1)
      nn = n*n
      if (taskid.eq.MASTER) then
         write(17,1330)
         do ms=1,ispin
            do i=n1(ms),n2(ms)
               ii = i-n1(ms)
               do j=n1(ms),n2(ms)
                  jj = j-n1(ms)
                  indx = (ii+1) + jj*n +(ms-1)*nn
                  write(17,1340) ms,ii+1,jj+1,
     >                           dbl_mb(hml(1)+indx-1),
     >                           dbl_mb(lmd(1)+indx-1),
     >             dbl_mb(hml(1)+indx-1)-dbl_mb(lmd(1)+indx-1)
               end do
            end do
         end do
      end if



*     **** check orthonormality ****
      if (taskid.eq.MASTER) then
         write(17,1350)
      end if

      call phafac()  !*** reset phase factors to r1 ***
      do ms=1,ispin
         do i=n1(ms),n2(ms)
            ii = i-n1(ms)+1
            do j=i,n2(ms)
               jj = j-n1(ms)+1
               call paw_overlap_matrix_gen(1,1,
     >                          dcpl_mb(psi1(1)+(i-1)*npack1),
     >                          dcpl_mb(psi1(1)+(j-1)*npack1),
     >                          w)
               if (taskid.eq.MASTER) then
                  write(17,1360) ms,ii,jj,w
               end if
            end do
         end do
      end do

*     **** close check file ****
      if (taskid.eq.MASTER) then
         close(17)
      end if


*     ***** diagonalize the hamiltonian matrix ****
      n = ne(1)
      nn = n*n
      call dcopy(2*nemax,0.0d0,0,dbl_mb(eig(1)),1)
      do ms=1,ispin
c        call eigen(n,ne(ms),
c    >              dbl_mb(hml(1)+(ms-1)*nn),
c    >              dbl_mb(eig(1)+(ms-1)*n),
c    >              dbl_mb(lmd(1)),ierr)
         call DSYEV('V','U',ne(ms),
     >              dbl_mb(hml(1)+(ms-1)*nn),n,
     >              dbl_mb(eig(1)+(ms-1)*n),
     >              dbl_mb(lmd(1)),(2*nemax*nemax),
     >              ierr)
        call eigsrt(dbl_mb(eig(1)+(ms-1)*n),
     >              dbl_mb(hml(1)+(ms-1)*nn),
     >              ne(ms),n)

      end do

*     *** rotate current psi ***
      call dcopy(2*npack1*nemax,0.0d0,0,dcpl_mb(psi2(1)),1)
      do ms=1,ispin
         do j=n1(ms),n2(ms)
            jj = j-n1(ms)
            do i=n1(ms),n2(ms)
               ii = i-n1(ms)
               indx = (ii+1) + jj*n + (ms-1)*nn

               call Pack_cc_daxpy(1,dbl_mb(hml(1)+indx-1),
     >                            dcpl_mb(psi1(1)+(i-1)*npack1),
     >                            dcpl_mb(psi2(1)+(j-1)*npack1))
            end do
         end do
      end do


*     *** rotate current v_psi ***
      call dcopy(2*npack1*nemax,dcpl_mb(psi0(1)),1,
     >                          dcpl_mb(psi1(1)),1)
      call dcopy(2*npack1*nemax,0.0d0,0,dcpl_mb(psi0(1)),1)
      do ms=1,ispin
         do j=n1(ms),n2(ms)
            jj = j-n1(ms)
            do i=n1(ms),n2(ms)
               ii = i-n1(ms)
               indx = (ii+1) + jj*n + (ms-1)*nn

               call Pack_cc_daxpy(1,dbl_mb(hml(1)+indx-1),
     >                            dcpl_mb(psi1(1)+(i-1)*npack1),
     >                            dcpl_mb(psi0(1)+(j-1)*npack1))
            end do
         end do
      end do






*                |***************************|
****************** report summary of results **********************
*                |***************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      call center_v_geom(vcx,vcy,vcz)
      call center_v_mass(vgx,vgy,vgz)


      if (taskid.eq.MASTER) then
         write(6,1300)
         write(6,1410)
         write(6,1420)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(6,1191) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           else
           write(6,1190) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           end if
         end do
         write(6,1200) cx,cy,cz
         write(6,1210) gx,gy,gz


         write(6,1421)
         write(6,1192) (I,ion_aname(I),
     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
         write(6,1200) vcx,vcy,vcz
         write(6,1210) vgx,vgy,vgz



         write(6,*)
         write(6,1321) r_charge,' (real space)'


         write(6,*)
         write(6,1430) E(2),E(2)/ion_nion()
         write(6,1440) E(5),E(5)/n2(ispin)
         write(6,1450) E(6),E(6)/n2(ispin)
         write(6,1460) E(7),E(7)/n2(ispin)
         write(6,1470) E(8),E(8)/ion_nion()
         write(6,1471) E(3),E(3)/n2(ispin)
         write(6,1472) ion_ke(),ion_ke()/ion_nion()


c         if (pspw_qmmm_found()) then
c            write(6,1700)
c            write(6,1701)
c            write(6,1702) E(11)
c            write(6,1703) E(12)
c            write(6,1704) E(13)
c            write(6,1705) E(14)
c            write(6,1706) E(15)
c         end if


         if (control_Nose()) then
           write(6,1473) E(9),E(9)/n2(ispin)
           write(6,1474) E(10),E(10)/ion_nion()
         end if
         write(6,1226) E(3),ion_ke(),ion_com_ke()
         eke = eke/dble(it_out)
         eke = 2.0d0*eke/kb/(ne(1)+ne(ispin))/pack_nwave_all(1)
         !eke = 2.0d0*eke/kb/(ne(1)+ne(ispin))

*       **** write out Temperatures ****
         write(6,1491) eke
         eki0 = ion_Temperature()
         write(6,1480) eki0
         write(6,1490) ion_com_Temperature()

         eave = E(25)/dble(icount)
         evar = E(26)/dble(icount)
         evar = evar - eave*eave
         cv = (evar)/(kb*ion_Temperature()**2)
         cv = cv/dble(ion_nion())
         write(6,1492) eave
         write(6,1493) evar
         write(6,1494) cv

         write(6,1500)
         NN=NE(1)-NE(2)
         EV=27.2116d0
         do i=1,NN
           write(6,1510) dbl_mb(EIG(1)+i-1),dbl_mb(EIG(1)+i-1)*EV
         end do
         do i=1,ne(2)
           write(6,1510) dbl_mb(EIG(1)+i-1+NN),
     >                    dbl_mb(EIG(1)+i-1+NN)*EV,
     >                    dbl_mb(EIG(1)+i-1+n1(2)-1),
     >                    dbl_mb(EIG(1)+i-1+n1(2)-1)*EV
         end do

*        *** Extra energy output added for QA test ****
         write(6,1600) E(2)
      end if

*                |***************************|
******************         Prologue          **********************
*                |***************************|

!*     **** calculate spin contamination ****
!      call Calculate_psi_spin2(ispin,ne,npack1,dcpl_mb(psi2(1)),w)
!
!*     **** calculate the Dipole ***
!      call Calculate_Dipole(ispin,ne,n2ft3d,dbl_mb(dn(1)),dipole)
!      
!
!*     ***** write psi2 wavefunctions ****
       call psi_write(ispin,ne,dcpl_mb(psi2(1)),-1,dum)
      call v_psi_write(ispin,ne,dcpl_mb(psi0(1)))

       !call psi_marat_write("marat3.elc",ispin,ne,dcpl_mb(psi2(1)))


!*     **** write geometry to rtdb ****
       call ion_write(rtdb)


*     **** deallocate heap memory ****
      call strfac_end()
      if (control_version().eq.3) call coulomb_end()
      if (control_version().eq.4) call coulomb2_end()
      call ke_end()
      call mask_end()
      call Pack_end()
      call G_end()
      call ion_end()
     
      call dealloc_paw_basis_data()
      call paw_proj_end()
      call paw_ovlp_end()
      call paw_nonlocal_end()
      call paw_comp_charge_end()
      call paw_mult_end()
c      call paw_vloc_end()
      call end_paw_kin_matrix()
      call end_paw_vloc_matrix()
      call end_paw_ion_matrix()
      call end_paw_overlap_matrix()
      call end_paw_hartree_matrix()
      call end_paw_core_matrix()
      call end_paw_comp_pot_matrix()
      call end_paw_comp_charge_matrix()
      call end_paw_pot_comp()
      call end_paw_pot_hartree()
c      call paw_density_end()
      call paw_xc_end()
      call paw_force_end()
      call paw_gaunt_end()
      if (control_Nose()) call Nose_end()


      value =           BA_free_heap(psir(2))
      value = value.and.BA_free_heap(Hpsi(2))
      value = value.and.BA_free_heap(dn(2))
      value = value.and.BA_free_heap(dn_cmp_smooth(2))
      value = value.and.BA_free_heap(eig(2))
      value = value.and.BA_free_heap(hml(2))
      value = value.and.BA_free_heap(lmd(2))
      value = value.and.BA_free_heap(lmd1(2))
      value = value.and.BA_free_heap(psi0(2))
      value = value.and.BA_free_heap(psi1(2))
      value = value.and.BA_free_heap(psi2(2))
      if (.not. value) call errquit('paw_cpmd:error freeing heap',0,2)


      call D3dB_pfft_end()

      call D3dB_end(1)
      if (control_version().eq.4) call D3dB_end(2)
      call Dne_end()


*     **** do anaylysis on MOTION files ****
      call cpmd_properties(rtdb)

*                |***************************|
****************** report consumed cputime   **********************
*                |***************************|
      if (taskid.eq.MASTER) then
         CALL current_second(cpu4)

         T1=CPU2-CPU1
         T2=CPU3-CPU2
         T3=CPU4-CPU3
         T4=CPU4-CPU1
         AV=T2/dble(icount*it_in)
         write(6,*)
         write(6,*) '-----------------'
         write(6,*) 'cputime in seconds'
         write(6,*) 'prologue    : ',T1
         write(6,*) 'main loop   : ',T2
         write(6,*) 'epilogue    : ',T3
         write(6,*) 'total       : ',T4
         write(6,*) 'cputime/step: ',AV
         write(6,*)
         write(6,*) '-------------------------------'
         write(6,*) 'Time spent doing:'
         write(6,*) '  FFTs                       : ', nwpw_timing(1)
         write(6,*) '  dot products               : ', nwpw_timing(2)
         write(6,*) '  orthonormalization         : ', nwpw_timing(3)
         write(6,*) '  exchange correlation       : ', nwpw_timing(4)
         write(6,*) '  local pseudopotentials     : ', nwpw_timing(5)
         write(6,*) '  non-local pseudopotentials : ', nwpw_timing(6)
         write(6,*) '  hartree potentials         : ', nwpw_timing(7)
         write(6,*) '  structure factors          : ', nwpw_timing(8)
         write(6,*) '  masking and packing        : ', nwpw_timing(9)
         write(6,*) '  total energy evaluation    : ',nwpw_timing(10)
         write(6,*) '  density                    : ',nwpw_timing(11)
         write(6,*) '  allocate and deallocate    : ',nwpw_timing(12)
         write(6,*) '  Hpsi and update            : ',nwpw_timing(13)
         write(6,*)
         CALL nwpw_MESSAGE(4)
      end if 


      call Parallel_Finalize()
      paw_cpmd = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*  PAW Car-Parrinello microcluster calculation     *')
 1030 FORMAT(10X,'*     [    extended Lagrangian molecular    ]      *')
 1035 FORMAT(10x,'*     [         dynamics simulation         ]      *')
 1036 FORMAT(10x,'*     [  Northwest Chemistry implementation ]      *')
 1040 FORMAT(10X,'*            version #1.00   10/01/03              *')
 1041 FORMAT(10X,'*    Authors: Marat Valiev and Eric J. Bylaska     *')
 1042 FORMAT(10X,'*    This code is based upon algorithms and code   *')
 1043 FORMAT(10X,'*    developed by the group of Prof. John H. Weare *')
 1044 FORMAT(10X,'*                                                  *')
 1045 FORMAT(10X,'*    References:                                   *')
 1046 FORMAT(10X,'*                                                  *')
 1047 FORMAT(10X,'*    M. Valiev, E. J. Bylaska, A. Gramada,         *')
 1048 FORMAT(10X,'*    and J. H. Weare,                              *')
 1049 FORMAT(10X,'*    Reviews in Modern  Quantum Chemistry,         *')
 1050 FORMAT(10X,'*    1684 (World Scientific, Singapore, 2002)      *')
 1051 FORMAT(10X,'*                                                  *')
 1052 FORMAT(10X,'*    E. J. Bylaska, M. Valiev, R. Kawai,           *')
 1053 FORMAT(10X,'*    and J. H. Weare,                              *')
 1054 FORMAT(10X,'*    Computer Physics  Communications, 143 (2002)  *')
 1055 FORMAT(10X,'*                                                  *')
 1056 FORMAT(10X,'*    M. Valiev and J. H. Weare,                    *')
 1057 FORMAT(10X,'*    J. Phys. Chem. A 103, 10588 (1999).           *')
 1058 FORMAT(10X,'*                                                  *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ PAW input data ===================')
 1111 FORMAT(/' number of processors used:',I3)
 1117 FORMAT( ' processor grid           :',I4,' x',I4)
 1112 FORMAT( ' parallel mapping         : slab')
 1113 FORMAT( ' parallel mapping         : hilbert')
 1115 FORMAT(/' options:')
 1118 FORMAT( ' parallel mapping         : hcurve')
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1140 FORMAT(/' elements involved in the cluster:')
 1141 FORMAT(5X,I2,': ',A4,'  ion charge:',F4.1,'  core charge:',F4.1)

 1143 FORMAT(12x,' augmentation sphere radius  :',F6.3)
c1144 FORMAT(12x,' compensation sigma          :',F6.3,
c    .   ' (',F6.3,' smooth)')
 1144 FORMAT(12x,' compensation sigma          :',F6.3)

 1150 FORMAT(12x,' total number of projectors  :',I3)

 1151 FORMAT(12x,' n_ps (n) l          eig    #projector')
 1152 FORMAT(14X,I3,' (',I1,') ',A,F13.6,I14)

 1153 FORMAT(12X,' local potential used           : ',i2)
 1154 FORMAT(12X,' number of non-local projections: ',i2)
 1155 FORMAT(12X,' semicore corrections included  : ',
     >       F6.3,' (radius) ',F6.3,' (charge)')
 1156 FORMAT(12X,' aperiodic cutoff radius        : ',F6.3)
 1159 FORMAT(/' total charge:',F8.3)
 1160 FORMAT(/' atomic composition:')
 1170 FORMAT(7(5X,A4,':',I3))
 1180 FORMAT(/' initial position of ions:')
 1181 FORMAT(/' initial velocity of ions:')
 1190 FORMAT(5X, I4, A5, ' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1191 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F6.3,' - fixed')
 1192 FORMAT(5X, I4, A5  ,' (',3F11.5,' )')
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1220 FORMAT(/' number of electrons: spin up=',I3,'  spin down=',I3,A)
 1222 format(5x,' initial kinetic energy: ',e12.5,' (psi)', 2x,
     >                                      e12.5,' (ion)',/50x,
     >                                      e12.5,' (c.o.m.)')
 1223 format(5x,' after scaling:          ',e12.5,' (psi)', 2x,
     >                                      e12.5,' (ion)')
 1224 format(5x,' increased energy:       ',e12.5,' (psi)', 2x,
     >                                      e12.5,' (ion)')
 1226 format(/' final kinetic energy:  ',e12.5,' (psi)', 2x,
     >                                      e12.5,' (ion)',/44x,
     >                                      e12.5,' (c.o.m.)')
 1230 FORMAT(/' supercell:')
 1231 FORMAT(5x,' volume : ',F10.1)
 1241 FORMAT(5x,' lattice:    a1=<',3f8.3,' >')
 1242 FORMAT(5x,'             a2=<',3f8.3,' >')
 1243 FORMAT(5x,'             a3=<',3f8.3,' >')
 1244 FORMAT(5x,' reciprocal: b1=<',3f8.3,' >')
 1245 FORMAT(5x,'             b2=<',3f8.3,' >')
 1246 FORMAT(5x,'             b3=<',3f8.3,' >')

 1250 FORMAT(5X,' density cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
 1251 FORMAT(5X,' wavefnc cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
 1260 FORMAT(5X,' smooth compensation (ewald) summation: cut radius=',
     &       F8.2,'  and',I3)
 1261 FORMAT(5X,'                   madelung=',f11.8)
 1262 FORMAT(5X,' core integration lmax=',I2)
 1270 FORMAT(/' technical parameters:')
 1271 FORMAT(5x, ' translation constrained')
 1272 FORMAT(5x, ' rotation constrained')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'fictitious mass=',F10.1)
 1290 FORMAT(5X, ' cooling/heatting rates: ',e12.5,' (psi)',2x,
     >                                       e12.5,' (ion)')
 1294 format(/' Constant Energy Simulation                     ')
 1295 format(/' Nose-Hoover Simulation - thermostat parameters:')
 1296 format(5x, 'SA decay rates  =',e8.3,' (elc)',e8.3,' (ion)')
 1297 format(5x, 'link = ',I3,
     >       ' Te =',f8.2,' Qe =',e8.3,' 2*pi/we=',e8.3,' Ee0=',e8.3)
 1298 format(5x, 'link = ',I3,
     >       ' Tr =',f8.2,' Qr =',e8.3,' 2*pi/wr=',e8.3,' Er0=',e8.3)
 1299 format(//' Mulliken analysis output generated            ')
 1300 FORMAT(//)
 1305 FORMAT(10X,'============ Car-Parrinello iteration ==============')
 1309 FORMAT(I8,2E19.10,2E14.5,2F9.1)
 1310 FORMAT(I8,2E19.10,2E14.5,F14.2)
 1311 format(100e19.10)
 1312 format(e14.6,i3)
 1313 format(3i4)
 1320 FORMAT(' number of electrons: spin up=',F11.5,'  down=',F11.5,A)
 1321 FORMAT(' total charge of system:',F11.5,A)
 1330 FORMAT(/' comparison between hamiltonian and lambda matrix')
 1331 FORMAT(/' Hamiltonian matrix')
 1340 FORMAT(I3,2I3,' H=',E16.7,', L=',E16.7,', H-L=',E16.7)
 1341 FORMAT(I3,2I3,' H=',E16.6)
 1350 FORMAT(/' orthonormality')
 1360 FORMAT(I3,2I3,E18.7)
 1370 FORMAT(I3)
 1380 FORMAT(' ''',a,'''',I4)
 1390 FORMAT(I3)
 1400 FORMAT(I3,3E18.8/3X,3E18.8)
 1410 FORMAT(10X,'=============  summary of results  =================')
 1420 FORMAT(/' final position of ions:')
 1421 FORMAT(/' final velocity of ions:')
 1430 FORMAT(/' total     energy    :',E19.10,' (',E15.5,'/ion)')
 1440 FORMAT( ' total orbital energy:',E19.10,' (',E15.5,'/electron)')
 1450 FORMAT( ' hartree   energy    :',E19.10,' (',E15.5,'/electron)')
 1460 FORMAT( ' exc-corr  energy    :',E19.10,' (',E15.5,'/electron)')
 1470 FORMAT( ' ion-ion   energy    :',E19.10,' (',E15.5,'/ion)')
 1471 FORMAT(/' Kinetic energy (elc)    :',E19.10,' (',E15.5,'/elc)')
 1472 FORMAT( ' Kinetic energy (ion)    :',E19.10,' (',E15.5,'/ion)')
 1473 FORMAT( ' thermostat energy (elc) :',E19.10,' (',E15.5,'/elc)')
 1474 FORMAT( ' thermostat energy (ion) :',E19.10,' (',E15.5,'/ion)')
 1480 FORMAT(' Temperature :    ',F10.1,' K (ion)')
 1490 FORMAT('             :    ',F10.1,' K (c.o.m.)')
 1491 FORMAT(' Temperature :    ',F10.1,' K (elc)')
 1492 FORMAT(/' Eaverage           :    ',E19.10)
 1493 FORMAT( ' Evariance          :    ',E19.10)
 1494 FORMAT( ' Cv - f*kb/(2*nion) :    ',E19.10)

 1500 FORMAT(/' orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))

 1600 FORMAT(/' Total PAW energy   :',E19.10)

 1700 FORMAT(/' QM/MM-pol-vib/CAV Energies')
 1701 FORMAT( ' --------------------------')
 1702 FORMAT( ' QM/MM energy           :',E19.10)
 1703 FORMAT( ' MM/MM energy           :',E19.10)
 1704 FORMAT( ' MM Polarization energy :',E19.10)
 1705 FORMAT( ' MM Vibration energy    :',E19.10)
 1706 FORMAT( ' (QM+MM)/Cavity energy  :',E19.10)

 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(6,9010) ierr
      call Parallel_Finalize()

      paw_cpmd = value
      return
      end

