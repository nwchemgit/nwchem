*
* $Id$
*

***********************************************************************
*                      cpmdv5-mpi  (MPI code)                         *
*                                                                     *
*     This is a developing cpsdv5 parallel code for NWChem            *
*       + mpi message passing library used                            *
*       + ngp is used instead of nfft in this proceudure              *
*       + error checking is based on aimd.h parameters                *
*         then control file                                           *
*       + my own slap-decomposed parallel 3d-FFT(real->complex) used  *
*                                                                     *
*                                                                     *
***********************************************************************

      logical function cpmdv5(rtdb)
CDEC$ OPTIMIZE:0
      implicit none
      integer rtdb

#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "errquit.fh"
cccc#include "frac_occ.fh"
      
      logical value

      real*8 kb
      parameter (kb=3.16679d-6)
      real*8 autoatm
      parameter (autoatm =290.360032539d6)


      
*     **** parallel variables ****
      integer  taskid,np,np_i,np_j
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,nfft3d,n2ft3d,ngrid_small(3)
      integer npack1
      real*8  a,b,c,alpha,beta,gamma

*     **** electronic variables ****
      real*8 icharge
      integer ispin
      integer ne(2),n1(2),n2(2),nemax,neall,neq(2),nemaxq
      real*8  en(2),en1(2),en2(2)
      real*8 dipole(3)

*     complex*16 psi1(nfft3d,nemax)
*     complex*16 psi2(nfft3d,nemax)
*     real*8     dn(n2ft3d,2)
*     complex*16 Hpsi(nfft3d,nemax)
*     complex*16 psir(nfft3d,nemax)
      integer psi0(2),psi1(2),psi2(2)
      integer occ0(2),occ1(2),occ2(2)
      integer dn(2)
      integer Hpsi(2),psir(2)
    

*     ***** energy variables ****
      real*8  E(60),eke,eave,evar,cv,have,hvar,qave,qvar,Egas

*     real*8  eig(2*nemax)
*     real*8  hml(2*nemax*nemax)
*     real*8  lmd(2*nemax*nemax)
      integer eig(2),hml(2),lmd(2),lmd1(2)

*     **** psi smearing block ****
      logical fractional
      integer smearoccupation,smeartype
      real*8 smearfermi(2),smearcorrection,smearkT

*     **** error variables ****
      integer ierr

*     **** local variables ****
      logical verlet,mulliken,SA,found,calc_pressure,found_bak
      logical field_exist
      integer ms
      real*8  gx,gy,gz,cx,cy,cz
      real*8  vgx,vgy,vgz,vcx,vcy,vcz
      real*8  ekg,eki0,eki1,sum
      real*8  eke0,eke1
      real*8  EV,pi,dt
      real*8  emotion_time_shift
      integer i,j,k,ia,n,nn
      integer ii,jj,index,indx
      integer icount,it_in,it_out,icount_shift
      real*8 w,sumall,pressure,stress(3,3),p1,p2
      real*8 Te_init,Tr_init,Te_new,Tr_new,sa_decay(2),sa_alpha(2)
      double precision tollz
      parameter(tollz=1d-16)
      integer mapping,mapping1d
      character*50 filename
      character*255 full_filename,full_bak
      integer tmp1(2)
  


*     **** external functions ****
      real*8      psp_zv,psp_rc,ewald_rcut,ion_amass
      real*8      ewald_mandelung,lattice_omega_small
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg,lattice_unitg_small,lattice_unita_small
      integer     ewald_ncut,ewald_nshl3d
      integer     psp_lmmax,psp_lmax,psp_locp,ion_nkatm
      character*4 ion_atom,ion_aname
      external    psp_zv,psp_rc,ewald_rcut,ion_amass
      external    ewald_mandelung,lattice_omega_small
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg,lattice_unitg_small,lattice_unita_small
      external    ewald_ncut,ewald_nshl3d
      external    psp_lmmax,psp_lmax,psp_locp,ion_nkatm
      external    ion_atom,ion_aname

      real*8   control_rti,control_rte,ion_rion
      real*8   ion_vion,ion_com_ke,ion_ke
      real*8   ion_Temperature,ion_com_Temperature
      external control_rti,control_rte,ion_rion
      external ion_vion,ion_com_ke,ion_ke
      external ion_Temperature,ion_com_Temperature
      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      logical  control_read,control_move,ion_init,ion_q_FixIon
      external control_read,control_move,ion_init,ion_q_FixIon
      logical  ion_q_xyzFixIon
      external ion_q_xyzFixIon
      character*14 ion_q_xyzFixIon_label
      external     ion_q_xyzFixIon_label

      integer  pack_nwave_all
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      external pack_nwave_all
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm

      character*12 control_boundry
      external     control_boundry

      logical  psp_semicore,pspw_qmmm_found,pspw_charge_found
      external psp_semicore,pspw_qmmm_found,pspw_charge_found
      real*8   psp_rcore,psp_ncore,psp_rlocal
      external psp_rcore,psp_ncore,psp_rlocal
      logical  pspw_Efield_found
      external pspw_Efield_found

      logical  control_Nose,control_Mulliken,Nose_restart
      external control_Nose,control_Mulliken,Nose_restart

      integer  Nose_Mchain,Nose_Nchain
      external Nose_Mchain,Nose_Nchain
  
      real*8   control_Nose_Te,Nose_Qe,Nose_Pe,Nose_Ee0
      external control_Nose_Te,Nose_Qe,Nose_Pe,Nose_Ee0

      real*8   control_Nose_Tr,Nose_Qr,Nose_Pr,Nose_Er0
      external control_Nose_Tr,Nose_Qr,Nose_Pr,Nose_Er0
      logical      v_psi_filefind
      external     v_psi_filefind
      real*8   nwpw_timing
      external nwpw_timing

      logical  control_out_of_time,control_new_vpsi
      external control_out_of_time,control_new_vpsi

      logical  pspw_HFX_localize2
      logical  control_SA,control_Fei,pspw_SIC,pspw_HFX,control_pressure
      real*8   control_SA_decay
      external pspw_HFX_localize2
      external control_SA,control_Fei,pspw_SIC,pspw_HFX,control_pressure
      external control_SA_decay

      integer  control_np_orbital,control_mapping,control_mapping1d
      external control_np_orbital,control_mapping,control_mapping1d


      logical  control_translation,control_rotation,control_balance
      external control_translation,control_rotation,control_balance
     
      logical  Dneall_m_allocate,Dneall_m_free,control_parallel_io
      external Dneall_m_allocate,Dneall_m_free,control_parallel_io

      real*8   Dneall_m_value,pspw_qmmm_lambda
      external Dneall_m_value,pspw_qmmm_lambda

      logical  meta_found,tamd_found,psp_U_psputerm,ion_disp_on
      external meta_found,tamd_found,psp_U_psputerm,ion_disp_on

      integer  ion_nconstraints,ion_ndof,control_ngrid_small
      external ion_nconstraints,ion_ndof,control_ngrid_small

      logical  psp_pawexist,ion_makehmass2,control_has_ngrid_small
      external psp_pawexist,ion_makehmass2,control_has_ngrid_small

      logical  nwpw_born_on
      external nwpw_born_on
      real*8   nwpw_born_screen
      external nwpw_born_screen
      logical  pspw_V_APC_on
      external pspw_V_APC_on
      real*8   control_gas_energy
      external control_gas_energy

      character*50 control_cell_name
      external     control_cell_name
      integer  Parallel_maxthreads
      external Parallel_maxthreads


*                            |************|
*****************************|  PROLOGUE  |****************************
*                            |************|

      value = .true.
      pi = 4.0d0*datan(1.0d0)

      call nwpw_timing_init()
      call dcopy(60,0.0d0,0,E,1)


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
         write(luout,1031)
         write(luout,1010)
         write(luout,1035)
         write(luout,1010)
         write(luout,1040)
         write(luout,1010)
         write(luout,1041)
         write(luout,1042)
         write(luout,1043)
         write(luout,1010)
         write(luout,1000)
         call nwpw_message(1)
         write(luout,1110)
         call util_flush(luout)
      end if
      call ga_sync()
      
      value = control_read(2,rtdb)
      call Parallel2d_Init(control_np_orbital())
      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)

      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0
      mapping = control_mapping()

*     **** initialize psi_data ****
      call psi_data_init(100)


*     **** initialize D3dB data structure ****
      call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      if (control_version().eq.4) 
     >   call D3dB_Init(2,2*ngrid(1),2*ngrid(2),2*ngrid(3),mapping)

      if (control_has_ngrid_small()) then
         ngrid_small(1) = control_ngrid_small(1)
         ngrid_small(2) = control_ngrid_small(2)
         ngrid_small(3) = control_ngrid_small(3)
         call D3dB_Init(3,ngrid_small(1),ngrid_small(2),ngrid_small(3),
     >                  mapping)
      end if


*     **** initialize lattice data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_init()
      call Pack_npack(1,npack1)

      call D3dB_pfft_init()
      call ga_sync()

*     ***** allocate psi2, psi1, and psi0 wavefunctions ****
      call psi_get_ne_occupation(ispin,ne,smearoccupation)
      if (smearoccupation.gt.0) then
         fractional = .true.
      else
         fractional = .false.
      end if
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
      if (fractional) then
      value = value.and.
     >        BA_alloc_get(mt_dbl,(ne(1)+ne(2)),'occ0',occ0(2),occ0(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(ne(1)+ne(2)),'occ1',occ1(2),occ1(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(ne(1)+ne(2)),'occ2',occ2(2),occ2(1))
      end if
      if (.not.value) call errquit('out of heap memory',0, MA_ERR)



*     *****  read psi2 wavefunctions ****
      call psi_read(ispin,ne,dcpl_mb(psi2(1)),
     >              smearoccupation,dbl_mb(occ2(1)))


*     **** move  wavefunction velocities ****
      if (control_new_vpsi()) then
        call v_psi_delete()
      end if


*     **** generate initial wavefunction velocities if it does not exist ****
      if (.not.v_psi_filefind()) then
        call v_psi_new(ispin,ne)
      end if


*     *****  read psi0 wavefunctions ****
      call dcopy(2*npack1*(neq(1)+neq(2)),0.0d0,0,dcpl_mb(psi1(1)),1)
      call v_psi_read(ispin,ne,dcpl_mb(psi1(1)))
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)
      nemax = ne(1)+ne(2)

*     **** allocate other variables *****
      value = BA_alloc_get(mt_dbl,(2*nemax),'eig',eig(2),eig(1))
      value = value.and.Dneall_m_allocate(0,hml)
      value = value.and.Dneall_m_allocate(0,lmd)
      value = value.and.Dneall_m_allocate(0,lmd1)

      value = value.and.
     >        BA_alloc_get(mt_dbl,(4*nfft3d),
     >                     'dn',dn(2),dn(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'Hpsi',Hpsi(2),Hpsi(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d*(neq(1)+neq(2)),
     >                     'psir',psir(2),psir(1))
      if (.not.value) call errquit('out of heap memory',0, MA_ERR)



*     **** read ions ****
      value = ion_init(rtdb)


*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)


*     **** allocate psp data structure and read in psedupotentials into it ****
      call psp_init()
      call psp_readall()
      if (psp_semicore(0)) call semicore_check()


*     **** initialize G,mask,ke,and coulomb data structures ****
      call ke_init()
      if (control_version().eq.3) call coulomb_init()
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()
      if (control_version().eq.3) call ewald_init()
      
*     **** initialize QM/MM ****
      call pspw_init_APC(rtdb)
      call pspw_qmmm_init(rtdb)
      call pspw_charge_init(rtdb)
      call pspw_Efield_init(rtdb,ispin,ne)
      field_exist = pspw_charge_found().or.pspw_Efield_found()

*     ******************************
*     **** scaling psi velocity ****
*     ******************************
      call dcopy(2*(neq(1)+neq(2))*npack1,dcpl_mb(psi1(1)),1,
     >                                  dcpl_mb(psi0(1)),1)
      call dscal(2*(neq(1)+neq(2))*npack1,control_rte(),
     >           dcpl_mb(psi1(1)),1)
      eke0 = 0.0d0
      eke1 = 0.0d0
      do i=1,(neq(1)+neq(2))
         call Pack_cc_dot(1,dcpl_mb(psi0(1)+(i-1)*npack1),
     >                      dcpl_mb(psi0(1)+(i-1)*npack1),
     >                     sum)
         eke0 = eke0 + sum
         call Pack_cc_dot(1,dcpl_mb(psi1(1)+(i-1)*npack1),
     >                      dcpl_mb(psi1(1)+(i-1)*npack1),
     >                    sum)
         eke1 = eke1 + sum
      end do

      call D1dB_SumAll(eke0)
      call D1dB_SumAll(eke1)
      eke0 = control_fake_mass()*eke0
      eke1 = control_fake_mass()*eke1
      call ion_init_ke(ekg,eki0,eki1)


*     **** Initialize thermostats ****
      if (control_Nose()) then
         call ke_ave(ispin,neq,dcpl_mb(psi2(1)),w,
     >               fractional,dbl_mb(occ2(1)))
         call Nose_Init((ne(1)+ne(2)),w)
      end if


*     **** Initialize simulated annealing ****
      SA       = .false.
      Te_init  = 0.0d0
      Tr_init  = 0.0d0
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


*     **** initialize two-electron Gaussian integrals ****
*     **** initialize paw ncmp*Vloc ****
      if (psp_pawexist()) then
         call nwpw_gintegrals_init()
         call nwpw_gintegrals_set(.true.)
         call psp_dE_ncmp_vloc_Qlm(ispin,.false.,dipole)
      end if

*     **** initialize metadynamics and tamd ****
      call meta_initialize(rtdb)
      call tamd_initialize(rtdb)

*     **** initialize dplot ****
      call dplot_iteration_init()

c*     **** initialize frac_occ data structure ****
c      call frac_occ_init(rtdb,ispin,ne)

*     **** initialize SIC and HFX ****
      call pspw_init_SIC(rtdb,ne)
      call pspw_init_HFX(rtdb,ispin,ne)


*     **** initialize DFT+U ****
      call psp_U_init()

*     **** initialize META GGA ****
      call nwpw_meta_gga_init(control_gga())

*     **** initialize vdw ****
      call vdw_DF_init()


*     **** initialize pressure ****
      calc_pressure = control_pressure().and.(control_version().eq.3)
      pressure      = 0.0d0
      p1            = 0.0d0
      p2            = 0.0d0
      if (calc_pressure) then
         call psp_stress_init()
         call psp_stress_readall()
      end if



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
         if (control_balance()) then
           write(luout,1114)
         else
           write(luout,1116)
         end if
         if (control_parallel_io()) then
           write(luout,1119)
         else
           write(luout,1120)
         end if
         write(luout,1123) Parallel_maxthreads()

         write(luout,1115)
         write(luout,1121) control_boundry(),control_version()
         if (ispin.eq.1) write(luout,1130) 'restricted'
         if (ispin.eq.2) write(luout,1130) 'unrestricted'

         call v_bwexc_print(luout,control_gga())

         if (fractional) write(luout,1132)
         call pspw_print_SIC(luout)
         call pspw_print_HFX(luout)
         if (ion_makehmass2()) write(luout,1135)
         write(luout,1140)
         do ia = 1,ion_nkatm()
            call psp_print(ia)
c           write(luout,1150) ia,ion_atom(ia),
c     >                    psp_zv(ia),psp_lmax(ia)
c           write(luout,1152) psp_lmax(ia)
c           write(luout,1153) psp_locp(ia)
c           write(luout,1154) psp_lmmax(ia)
c           if (control_version().eq.4) write(luout,1156) psp_rlocal(ia)
c           if (psp_semicore(ia))
c     >         write(luout,1155) psp_rcore(ia),psp_ncore(ia)
c           write(luout,1151) (psp_rc(i,ia),i=0,psp_lmax(ia))
         end do

         icharge = -(ne(1)+ne(ispin))
         en(1)     = ne(1)
         en(ispin) = ne(ispin)
         if (fractional) then
            icharge = 0.0d0
            do ms=1,ispin
            en(ms) =0.0
            do i=n1(ms),n2(ms)
              icharge = icharge - (3-ispin)*dbl_mb(occ2(1)+i-1)
              en(ms) = en(ms) + dbl_mb(occ2(1)+i-1)
            end do
            end do
         end if

         do ia=1,ion_nkatm()
           icharge = icharge + ion_natm(ia)*psp_zv(ia)
         end do
         write(luout,1159) icharge

         write(luout,1160)
         write(luout,1170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
         write(luout,1180)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),
     >                    (ion_rion(K,I),K=1,3),ion_amass(i)/1822.89d0
           else if (ion_q_xyzFixIon(I)) then
           write(luout,1194) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_q_xyzFixIon_label(I)
           else
           write(luout,1190) I,ion_aname(I),
     >                    (ion_rion(K,I),K=1,3),ion_amass(i)/1822.89d0
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz


         write(luout,1181)
         write(luout,1192) (I,ion_aname(I),
     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
         write(luout,1200) vcx,vcy,vcz
         write(luout,1210) vgx,vgy,vgz
         write(luout,1211) ion_nconstraints(),ion_ndof()

         call pspw_charge_Print(luout)
         call pspw_Efield_Print(luout)

         if (fractional) then
           write(luout,1219) en(1),en(ispin),' (   fractional)'
           write(luout,1221) ne(1),neq(1),
     >                   ne(ispin),neq(ispin),' (Fourier space)'
         else
           write(luout,1220) ne(1),neq(1),
     >                   ne(ispin),neq(ispin),' (Fourier space)'
           write(luout,1221) ne(1),neq(1),
     >                   ne(ispin),neq(ispin),' (Fourier space)'
         end if
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
         if (control_version().eq.3) then
         write(luout,1260) ewald_rcut(),ewald_ncut()
         write(luout,1261) ewald_mandelung()
         end if

         if (control_has_ngrid_small()) then
            write(luout,1229)
            write(luout,1233) control_cell_name()
            write(luout,1241) lattice_unita_small(1,1),
     >                    lattice_unita_small(2,1),
     >                    lattice_unita_small(3,1)
            write(luout,1242) lattice_unita_small(1,2),
     >                    lattice_unita_small(2,2),
     >                    lattice_unita_small(3,2)
            write(luout,1243) lattice_unita_small(1,3),
     >                    lattice_unita_small(2,3),
     >                    lattice_unita_small(3,3)
            write(luout,1244) lattice_unitg_small(1,1),
     >                    lattice_unitg_small(2,1),
     >                    lattice_unitg_small(3,1)
            write(luout,1245) lattice_unitg_small(1,2),
     >                    lattice_unitg_small(2,2),
     >                    lattice_unitg_small(3,2)
            write(luout,1246) lattice_unitg_small(1,3),
     >                    lattice_unitg_small(2,3),
     >                    lattice_unitg_small(3,3)
            call lattice_small_abc_abg(a,b,c,alpha,beta,gamma)
            write(luout,1232) a,b,c,alpha,beta,gamma
            write(luout,1231) lattice_omega_small()
            write(luout,1250) lattice_ecut(),
     >                 ngrid_small(1),ngrid_small(2),ngrid_small(3),
     >                 pack_nwave_all(2),pack_nwave(2)
            write(luout,1251) lattice_wcut(),
     >                 ngrid_small(1),ngrid_small(2),ngrid_small(3),
     >                 pack_nwave_all(3),pack_nwave(3)
         end if

         write(luout,1270)
         if (.not.control_translation()) write(luout,1271)
         if (.not.control_rotation())    write(luout,1272)
         write(luout,1280) control_time_step(),control_fake_mass()
         write(luout,1290) control_rte(),control_rti()
         call ion_scaling_atoms_print(luout)
         write(luout,1281) control_it_in()*control_it_out(),
     >                 control_it_in(),control_it_out()
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
     >                   Nose_Pr(i),Nose_Er0(i)
           end do
         else
           write(luout,1294)
         end if
        if (calc_pressure) write(luout,1293)
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
*     **** open xyz and MOTION and dipole file ****
      call xyz_init()          ! unit=18
      call MOTION_init(rtdb)   ! unit=19
      call dipole_motion_init(rtdb)   ! unit=36

*     *** fei io ****
      call fei_init(rtdb)


*     **** ecce print ****
      call ecce_print_module_entry('task Car-Parrinello')
      !call ecce_print_module_entry('driver')
      call movecs_ecce_print_off()




*     ************************************
*     **** open up other MOTION files ****
*     ************************************
      icount_shift = 0

*     **** open EMOTION file ****
      if (.not.btdb_cget(rtdb,'cpmd:emotion_filename',1,filename))
     >  call util_file_prefix('emotion',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)
      if (taskid.eq.MASTER) then

*        **** check for backup file ****
         call util_file_name_noprefix('EMOTION99-bak',.false.,
     >                                .false.,
     >                                full_bak)
         inquire(file=full_bak,exist=found_bak)
         if (found_bak) then
            write(*,*) 
            write(*,*) "EMOTION99-bak exists:"
            i=index(full_bak,' ')
            j=index(full_filename,' ')
            write(*,*) "   Copying ",full_bak(1:i),
     >                 " to ",full_filename(1:j)
            write(*,*)
            call util_file_copy(full_bak,full_filename)
         end if

         emotion_time_shift = 0.0d0
         icount_shift       = 0
         inquire(file=full_filename,exist=found)
         if (found) then

*          **** make a new backup file ***
           call util_file_copy(full_filename,full_bak)

           open(unit=31,file=full_filename,form='formatted',
     >          status='old')
           do while (found)
           read(31,*,end=100) emotion_time_shift,w,sum,gx,gy,gz
           E(25) = E(25) + sum                          !*** take care of running sums ***
           E(26) = E(26) + sum*sum
           E(27) = E(27) + (sum+gx+gy)
           E(28) = E(28) + (sum+gx+gy)**2
           E(23) = E(23) + gz
           E(24) = E(24) + gz*gz
           icount_shift = icount_shift + 1
           end do
  100      continue
#if defined(FUJITSU) || defined(PSCALE) || defined(__crayx1) || defined(GCC46)
           backspace 31
#endif
         else
           open(unit=31,file=full_filename,form='formatted',
     >          status='new')
         end if
      end if


*     **** open EIGMOTION file ****
      if (mulliken) then
        if (.not.btdb_cget(rtdb,'cpmd:eigmotion_filename',1,filename))
     >    call util_file_prefix('eigmotion',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)
      if (taskid.eq.MASTER) 
     >   open(unit=32,file=full_filename,form='formatted')
      end if

*     **** open HMOTION file ****
      if (mulliken) then
       if (.not.btdb_cget(rtdb,'cpmd:hmotion_filename',1,filename))
     >  call util_file_prefix('hmotion',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)
      if (taskid.eq.MASTER) 
     >   open(unit=34,file=full_filename,form='formatted')
      end if

*     **** open OMOTION file ****
      if (mulliken) call Orb_Init(rtdb,ispin,ne) !unit=33

c*     **** write initial position to xyz data ****
c      call xyz_write()

*     ***** first step using velocity ****
      verlet = .false.
      call inner_loop_md(verlet,sa_alpha,ispin,ne,neq,
     >             npack1,nfft3d,nemaxq,
     >             dcpl_mb(psi0(1)),
     >             dcpl_mb(psi1(1)),
     >             dcpl_mb(psi2(1)),
     >             dbl_mb(dn(1)),
     >             1,0,E,
     >             dbl_mb(hml(1)),dbl_mb(lmd(1)),dbl_mb(lmd1(1)),
     >             dcpl_mb(Hpsi(1)),dcpl_mb(psir(1)),
     >             calc_pressure,pressure,p1,p2,
     >             fractional,
     >             dbl_mb(occ0(1)),dbl_mb(occ1(1)),dbl_mb(occ2(1)))


      if (taskid.eq.MASTER) call current_second(cpu2)
      if ((taskid.eq.MASTER).and.(.not.calc_pressure)) 
     >   CALL nwpw_message(6)
      if ((taskid.eq.MASTER).and.(calc_pressure)) 
     >   CALL nwpw_message(9)
      
      it_in  = control_it_in()
      it_out = control_it_out()
      icount = 0
      verlet = .true.
      eke    = 0.0d0
      dt = control_time_step()

      if (it_out.lt.1) goto 102

      Te_new = Te_init
      Tr_new = Tr_init
  101 continue
         icount = icount + 1
         call inner_loop_md(verlet,sa_alpha,ispin,ne,neq,
     >             npack1,nfft3d,nemaxq,
     >             dcpl_mb(psi0(1)),
     >             dcpl_mb(psi1(1)),
     >             dcpl_mb(psi2(1)),
     >             dbl_mb(dn(1)),
     >             it_in,((icount-1)*it_in),
     >             E,
     >             dbl_mb(hml(1)),dbl_mb(lmd(1)),dbl_mb(lmd1(1)),
     >             dcpl_mb(Hpsi(1)), dcpl_mb(psir(1)),
     >             calc_pressure,pressure,p1,p2,
     >             fractional,
     >             dbl_mb(occ0(1)),dbl_mb(occ1(1)),dbl_mb(occ2(1)))

         eke = eke + E(3)


         !**** metadynamics and tamd update ****
         call meta_update(ispin,neq,dcpl_mb(psi1(1)),E)
         call tamd_update(ispin,neq,dcpl_mb(psi1(1)),E)

         if (taskid.eq.MASTER) then 

           if (calc_pressure) then
             if (SA) then
             write(luout,1309) icount*it_in,E(1),E(2),E(3),E(4),
     >                     Te_new,Tr_new,pressure*autoatm
             else
             write(luout,1310) icount*it_in,E(1),E(2),E(3),E(4),
     >                     ion_Temperature(),pressure*autoatm
             end if
           else
             if (SA) then
             write(luout,1309) icount*it_in,E(1),E(2),E(3),E(4),
     >                     Te_new,Tr_new
             else
             write(luout,1310) icount*it_in,E(1),E(2),E(3),E(4),
     >                     ion_Temperature()
             end if
           end if
           call util_flush(luout)

*          **** write out EMOTION data ****
           qave = E(23)/dble(icount+icount_shift)
           qvar = E(24)/dble(icount+icount_shift)
           qvar = qvar - qave*qave
           eave = E(25)/dble(icount+icount_shift)
           evar = E(26)/dble(icount+icount_shift)
           evar = evar - eave*eave
           have = E(27)/dble(icount+icount_shift)
           hvar = E(28)/dble(icount+icount_shift)
           hvar = hvar - have*have
           if (control_Nose()) then
             write(31,1311) icount*it_in*dt + emotion_time_shift,
     >                    e(1),e(2),e(3),e(4),e(14),e(5),e(6),
     >                    e(7),e(8),e(9),e(10),
     >                    eave,evar,have,hvar,qave,qvar,
     >                    ion_Temperature(),pressure
           else
             write(31,1311) icount*it_in*dt + emotion_time_shift,
     >                    e(1),e(2),e(3),e(4),e(14),e(5),e(6),
     >                    e(7),e(8),
     >                    eave,evar,have,hvar,qave,qvar,
     >                    ion_Temperature(),pressure
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


*        **** write xyz, MOTION, and dipole data ****
         call xyz_write()
         call MOTION_write(icount*it_in*control_time_step())

         call dipole_motion_write((control_version().eq.3),
     >                     (icount*it_in*control_time_step()),
     >                     ispin,ne,neq,npack1,nfft3d,
     >                     dbl_mb(dn(1)),
     >                     dcpl_mb(psi1(1)))


*        **** write OMOTION data ****
         if (mulliken) call Orb_Write(dcpl_mb(psi1(1)))

*        **** update thermostats using SA decay ****
         if (SA) then
           if(abs(sa_decay(1)).gt.tollz) 
     *     t1 = icount*it_in*dt/sa_decay(1)
           if(abs(sa_decay(2)).gt.tollz) 
     *     t2 = icount*it_in*dt/sa_decay(2)
           Te_new = Te_init*dexp(-t1)
           Tr_new = Tr_init*dexp(-t2)
           call Nose_reset_T(Te_new,Tr_new)
         end if


*        **** exit early ****
         if (control_out_of_time()) then
            if (taskid.eq.MASTER) 
     >       write(luout,*) ' *** out of time. iteration terminated'
            go to 102
         end if
      if (icount.lt.it_out) go to 101
      if (taskid.eq.MASTER) 
     > write(luout,*)
     > '*** arrived at the Maximum iteration.   terminated.'
 
*::::::::::::::::::::  end of iteration loop  :::::::::::::::::::::::::

  102 continue

*     **** close xyz,MOTION and dipole files ****
      call xyz_end()
      call MOTION_end()
      call dipole_motion_end()
      if (taskid.eq.MASTER) then
        close(unit=31)
        close(unit=32)
        close(unit=34)

*        **** remove EMOTION backup file ***
         call util_file_name_noprefix('EMOTION99-bak',.false.,
     >                                .false.,
     >                                full_bak)
         call util_file_unlink(full_bak)
      end if

*     *** close fei io ****
      call fei_end()

*     **** close OMOTION file ****
      if (mulliken) call Orb_End()

*     **** ecce print ****
      !call ecce_print_module_exit('driver', 'ok')
      call ecce_print_module_exit('task Car-Parrinello', 'ok')


*     **** finalize pressure ****
      if (calc_pressure) then
         call psp_stress_end()
      end if


      if (taskid.eq.MASTER) CALL nwpw_message(3)
      if (taskid.eq.MASTER) call current_second(cpu3)


*         |****************************************|
*********** produce CHECK file and diagonalize hml *****************
*         |****************************************|

*     **** produce CHECK FILE ****
      if (taskid.eq.MASTER) then
         call util_file_name('CHECK',.true.,
     >                               .false.,
     >                       full_filename)
         open(unit=17,file=full_filename,form='formatted')
      end if

*     **** check total number of electrons ****
      do ms =1,ispin
         call D3dB_r_dsum(1,dbl_mb(dn(1)+(ms-1)*n2ft3d),sumall)
         en1(ms) = sumall*lattice_omega()
     >             /dble(ngrid(1)*ngrid(2)*ngrid(3))
      end do

      if (psp_pawexist()) then
         if (.not.BA_push_get(mt_dbl,n2ft3d,'tmp1',tmp1(2),tmp1(1)))
     >   call errquit(
     >        'cpmdv5: out of stack memory',0,MA_ERR)

         call psp_qlm_atom(ispin,neq,dcpl_mb(psi1(1)))
         do ms=1,ispin
           call nwpw_compcharge_gen_dn_cmp_smooth_ms(ms,dbl_mb(tmp1(1)))
           call Pack_c_unpack(0,dbl_mb(tmp1(1)))
           call D3dB_cr_fft3b(1,dbl_mb(tmp1(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))
           call D3dB_r_dsum(1,dbl_mb(tmp1(1)),sumall)
           en2(ms) = sumall*lattice_omega()
     >              /dble(ngrid(1)*ngrid(2)*ngrid(3))
         end do
         if (.not.BA_pop_stack(tmp1(2)))
     >   call errquit(
     >        'cpmdv5: popping stack memory',0,MA_ERR)
      else
         en2(1) = 0.0d0
         en2(2) = 0.0d0
      end if
      en(1) = en1(1)+en2(1)
      en(2) = en1(2)+en2(2)

      if (taskid.eq.MASTER) then
         write(17,1320) (en(ms),ms=1,ispin)
         if (psp_pawexist()) then
            write(17,1322) (en1(ms),ms=1,ispin)
            write(17,1323) (en2(ms),ms=1,ispin)
         end if
      end if

*     **** comparison between hamiltonian an lambda matrix ****
      if (taskid.eq.MASTER) write(17,1330)
      do ms=1,ispin
         do i=1,ne(ms)
         do j=1,ne(ms)
            w   = Dneall_m_value(0,ms,i,j,dbl_mb(hml(1)))
            sum = Dneall_m_value(0,ms,i,j,dbl_mb(lmd(1)))

            if (taskid.eq.MASTER)
     >      write(17,1340) ms,i,j,w,sum,w-sum

         end do
         end do
      end do

*     **** check orthonormality ****
      if (taskid.eq.MASTER) then
         write(17,1350)
      end if

      call Dneall_ffm_Multiply(0,dcpl_mb(psi1(1)),
     >                           dcpl_mb(psi1(1)),npack1,
     >                           dbl_mb(lmd(1)))
      do ms=1,ispin
         do j=1,ne(ms)
         do i=j,ne(ms)
            w  = Dneall_m_value(0,ms,i,j,dbl_mb(lmd(1)))
            if (taskid.eq.MASTER) write(17,1360) ms,i,j,w
         end do
         end do
      end do

*     **** close check file ****
      if (taskid.eq.MASTER) then
         close(17)
      end if



*     ***** do not diagonalize the hamiltonian matrix *****
      if (pspw_SIC()) then
        call dcopy(2*npack1*nemaxq,
     >             dcpl_mb(psi1(1)),1,
     >             dcpl_mb(psi2(1)),1)
      
*     ***** diagonalize the hamiltonian matrix *****
      else

c         if (fractional) then
c           call Dneall_m_HmltimesSA(0,dbl_mb(hml(1)),dbl_mb(fweight(1)))
c         end if

         call Dneall_m_diagonalize(0,dbl_mb(hml(1)),
     >                               dbl_mb(eig(1)),.false.)


c         if (fractional) then
c         do ii=1,ne(ms)
c           dbl_mb(eig(1)+(ii-1)+(ms-1)*n)
c     >       =dbl_mb(eig(1)+(ii-1)+(ms-1)*n)
c     >       /dbl_mb(fweight(1)+(ii-1)+(ms-1)*n)
c         end do
c         end if


*        **** do not rotate for wannier localization algorithm ****
         if (.not.pspw_HFX_localize2()) then
*           *** rotate current psi ***
            call Dneall_fmf_Multiply(0,dcpl_mb(psi1(1)),npack1,
     >                                 dbl_mb(hml(1)),  1.0d0,
     >                                 dcpl_mb(psi2(1)),0.0d0)


*           *** rotate current v_psi ***
            call dcopy(2*npack1*nemaxq,dcpl_mb(psi0(1)),1,
     >                                 dcpl_mb(psi1(1)),1)

            call Dneall_fmf_Multiply(0,dcpl_mb(psi1(1)),npack1,
     >                                 dbl_mb(hml(1)),  1.0d0,
     >                                 dcpl_mb(psi0(1)),0.0d0)
         end if

      end if



*                |***************************|
****************** report summary of results **********************
*                |***************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      call center_v_geom(vcx,vcy,vcz)
      call center_v_mass(vgx,vgy,vgz)

      if (taskid.eq.MASTER) then
         call print_elapsed_time(icount*it_in*dt)
         write(luout,1300)
         write(luout,1410)
         write(luout,1420)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           else if (ion_q_xyzFixIon(I)) then
           write(6,1194) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_q_xyzFixIon_label(I)
           else
           write(luout,1190) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz


         write(luout,1421)
         write(luout,1192) (I,ion_aname(I),
     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
         write(luout,1200) vcx,vcy,vcz
         write(luout,1210) vgx,vgy,vgz
         write(luout,1211) ion_nconstraints(),ion_ndof()

         call pspw_charge_Print(luout)
         call pspw_Efield_Print(luout)

         write(luout,*)
         write(luout,1320) en(1),en(ispin),' (real space)'
         if (psp_pawexist()) then
            write(luout,1322) en1(1),en1(ispin),' (real space)'
            write(luout,1323) en2(1),en2(ispin),' (real space)'
         end if

         if (psp_pawexist()) then
            write(luout,1434) (E(1)+E(36)+E(45)),
     >                        (E(1)+E(36)+E(45))/ion_nion()
         end if

*       **** write APC potential and charges ***
        if (pspw_V_APC_on()) call pspw_shortprint_APC(luout)

         write(luout,1430) E(2),E(2)/ion_nion()
         write(luout,1440) E(5),E(5)/n2(ispin)
         write(luout,1450) E(6),E(6)/n2(ispin)
         write(luout,1460) E(7),E(7)/n2(ispin)
         if (pspw_SIC()) then
           write(luout,1455) E(16),E(16)/n2(ispin)
           write(luout,1456) E(17),E(17)/n2(ispin)
         end if
         if (pspw_HFX()) then
           write(luout,1457) E(20),E(20)/n2(ispin)
         end if
         if (psp_U_psputerm()) then
           write(luout,1458) E(29),E(29)/n2(ispin)
         end if
         if (meta_found()) then
           write(luout,1459) E(31),E(31)/ion_nion()
         end if
         if (tamd_found()) then
           write(luout,1461) E(34),E(34)/ion_nion()
         end if
         if (pspw_V_APC_on()) then
           write(luout,1505) E(52),E(52)/ion_nion()
         end if
         write(luout,1470) E(8),E(8)/ion_nion()
         write(luout,1471) E(3),E(3)/n2(ispin)
         write(luout,1472) ion_ke(),ion_ke()/ion_nion()


         if (pspw_qmmm_found()) then
            write(luout,1700)
            write(luout,1701)
            write(luout,1702) E(11)
            write(luout,1703) E(12)
            write(luout,1704) E(13)
            qave = E(23)/dble(icount+icount_shift)
            qvar = E(24)/dble(icount+icount_shift)
            qvar = qvar - qave*qave
            write(luout,1707) pspw_qmmm_lambda()
            write(luout,1705) E(14),qave,qvar
            !write(luout,1706) qave,qvar
         end if
        if (ion_disp_on()) then
            write(luout,1720) E(33)
        end if

        if (field_exist) then
           write(luout,1800)
           write(luout,1801)
           write(luout,1805) E(19)+E(20)+E(21)
           write(luout,1802) E(19)
           write(luout,1803) E(20)
           write(luout,1804) E(21)
        end if

         if (control_Nose()) then
           write(luout,1473) E(9),E(9)/n2(ispin)
           write(luout,1474) E(10),E(10)/ion_nion()
         end if
         write(luout,1226) E(3),ion_ke(),ion_com_ke()
         eke = eke/dble(it_out)
         eke = 2.0d0*eke/kb/(ne(1)+ne(ispin))/pack_nwave_all(1)
         !eke = 2.0d0*eke/kb/(ne(1)+ne(ispin))

*       **** write out Temperatures ****
         write(luout,1491) eke
         eki0 = ion_Temperature()
c         if (pspw_qmmm_found()) then
c            eki1 =pspw_qmmm_Temperature()
c            sum = ion_nion() + pspw_qmmm_nion() - 2.0d0
c            eki0 = eki0*((ion_nion()-2.0d0)/sum)
c     >           + eki1*((pspw_qmmm_nion()-2.0d0)/sum)
c         end if
         write(luout,1480) eki0
         write(luout,1490) ion_com_Temperature()

         eave = E(25)/dble(icount+icount_shift)
         evar = E(26)/dble(icount+icount_shift)
         evar = evar - eave*eave
         have = E(27)/dble(icount+icount_shift)
         hvar = E(28)/dble(icount+icount_shift)
         hvar = hvar - have*have
         cv = (evar)/(kb*ion_Temperature()**2)
         cv = cv/dble(ion_nion())
         write(luout,1492) eave,have
         write(luout,1493) evar,hvar
         write(luout,1494) cv

*        **** write out diagonal <psi|H|psi> matrix ****
         if (pspw_SIC()) then

           n = ne(1)
           nn = n*n
           do ms=1,ispin
             if (ms.eq.1) write(luout,1331)
             if (ms.eq.2) write(luout,1332)
             !*** call Gainsville matrix output ***
             call output(dbl_mb(hml(1)+(ms-1)*nn),
     >                    1,ne(ms),1,ne(ms),
     >                    n,n,1)
           end do


*        **** write out KS eigenvalues ****
         else
         write(luout,1500)
         NN=NE(1)-NE(2)
         EV=27.2116d0

         if (fractional) then
           do i=1,NN
             write(luout,1511) dbl_mb(EIG(1)+i-1),
     >                     dbl_mb(EIG(1)+i-1)*EV,
     >                     dbl_mb(occ2(1)+i-1)
           end do
           do i=1,ne(2)
             write(luout,1511) dbl_mb(EIG(1)+i-1+NN),
     >                     dbl_mb(EIG(1)+i-1+NN)*EV,
     >                     dbl_mb(occ2(1)+i-1+NN),
     >                     dbl_mb(EIG(1)+i-1+n1(2)-1),
     >                     dbl_mb(EIG(1)+i-1+n1(2)-1)*EV,
     >                     dbl_mb(occ2(1)+i-1+n1(2)-1)
           end do
         else
           do i=1,NN
             write(luout,1510) dbl_mb(EIG(1)+i-1),dbl_mb(EIG(1)+i-1)*EV
           end do
           do i=1,ne(2)
             write(luout,1510) dbl_mb(EIG(1)+i-1+NN),
     >                     dbl_mb(EIG(1)+i-1+NN)*EV,
     >                     dbl_mb(EIG(1)+i-1+n1(2)-1),
     >                     dbl_mb(EIG(1)+i-1+n1(2)-1)*EV
           end do
         end if

         end if
      end if

*     **** write out extended Born solvation energies ****
      if (nwpw_born_on()) then
         Egas = control_gas_energy()
         if (taskid.eq.MASTER) then
            write(luout,1740)
            write(luout,1741) nwpw_born_screen()
            write(luout,1745) E(52),E(52)*27.2116d0*23.06d0
            if (dabs(Egas).gt.1.0d-6)
     >            write(luout,1746) E(2)-Egas,
     >                             (E(2)-Egas)*27.2116d0*23.06d0
            call nwpw_born_print(luout,Egas,E(2))
         end if
      end if

      if (taskid.eq.MASTER) then
*        *** Extra energy output added for QA test ****
         write(luout,1600) E(2)
      end if

*                |***************************|
******************         Prologue          **********************
*                |***************************|

*     **** calculate spin contamination ****
      call Calculate_psi_spin2(ispin,ne,npack1,dcpl_mb(psi2(1)),
     >                         fractional,dbl_mb(occ2(1)),
     >                         dipole)

*     **** calculate the Dipole ***
      call Calculate_Dipole(ispin,ne,n2ft3d,dbl_mb(dn(1)),dipole)

      
*     ***** write wavefunctions and v_wavefunctions ****
      call psi_write(ispin,ne,dcpl_mb(psi2(1)),
     >               smearoccupation,dbl_mb(occ2(1)))
      call v_psi_write(ispin,ne,dcpl_mb(psi0(1)))

*     **** write geometry to rtdb ****
      call pspw_charge_write(rtdb)
      call ion_write(rtdb)

*     **** deallocate heap memory ****
      if (control_version().eq.3) call ewald_end()
      call strfac_end()
      if (controL_version().eq.3) call coulomb_end()
      if (controL_version().eq.4) call coulomb2_end()
      call ke_end()
      call mask_end()
      call Pack_end()
      call G_end()
      call psp_U_end()
      call vdw_DF_end()
      call nwpw_meta_gga_end()
      call pspw_end_SIC()
      call pspw_end_HFX()
      call pspw_end_APC()
      call pspw_qmmm_end()
      call meta_finalize(rtdb)
      call tamd_finalize(rtdb)
      call dplot_iteration_end()
      call pspw_charge_end()
      call pspw_Efield_end()
c      call frac_occ_end()
      if (control_Nose()) call Nose_end()
      if (psp_pawexist()) call nwpw_gintegrals_end()
      !if (psp_pawexist()) call nwpw_cgintegrals_end()

      call ion_end()
      call psp_end()
      call ion_end_FixIon()

      value = BA_free_heap(psir(2))
      value = BA_free_heap(Hpsi(2))
      value = BA_free_heap(dn(2))
      value = BA_free_heap(eig(2))
      value = Dneall_m_free(hml)
      value = Dneall_m_free(lmd)
      value = Dneall_m_free(lmd1)
      value = BA_free_heap(psi0(2))
      value = BA_free_heap(psi1(2))
      value = BA_free_heap(psi2(2))
      if (fractional) then
      value = BA_free_heap(occ0(2))
      value = BA_free_heap(occ1(2))
      value = BA_free_heap(occ2(2))
      end if
      call D3dB_pfft_end()
      call D3dB_end(1)
      if (control_version().eq.4) call D3dB_end(2)
      if (control_has_ngrid_small()) call D3dB_end(3)
      call Dne_end()
      call psi_data_end()

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
         write(luout,*)
         write(luout,*) '-----------------'
         write(luout,*) 'cputime in seconds'
         write(luout,*) 'prologue    : ',T1
         write(luout,*) 'main loop   : ',T2
         write(luout,*) 'epilogue    : ',T3
         write(luout,*) 'total       : ',T4
         write(luout,*) 'cputime/step: ',AV
         write(luout,*)

         call nwpw_timing_print_final(.true.,(icount*it_in))
         CALL nwpw_message(4)
      end if 


      call Parallel2d_Finalize()
      call Parallel_Finalize()
      cpmdv5 = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*     Car-Parrinello microcluster calculation      *')
 1030 FORMAT(10X,'*      [   extended Lagrangian molecular   ]       *')
 1031 FORMAT(10X,'*      [        dynamics simulation        ]       *')
 1035 FORMAT(10x,'*      [ NorthWest Chemistry implementation ]      *')
 1040 FORMAT(10X,'*            version #5.00   06/01/99              *')
 1041 FORMAT(10X,'*    This code was developed by Eric J. Bylaska,   *')
 1042 FORMAT(10X,'*    and was based upon algorithms and code        *')
 1043 FORMAT(10X,'*    developed by the group of Prof. John H. Weare *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1111 FORMAT(/' number of processors used:',I10)
 1112 FORMAT( ' parallel mapping         :      1d-slab')
 1113 FORMAT( ' parallel mapping         :   2d-hilbert')
 1114 FORMAT( ' parallel mapping         :     balanced')
 1115 FORMAT(/' options:')
 1116 FORMAT( ' parallel mapping         : not balanced')
 1117 FORMAT( ' processor grid           :',I4,' x',I4)
 1118 FORMAT( ' parallel mapping         :    2d-hcurve')
 1119 FORMAT( ' parallel io              :        on')
 1120 FORMAT( ' parallel io              :       off')
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1123 FORMAT( ' number of threads        :',I10)
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1132 FORMAT(5X,' using fractional occupation')
c 1135 FORMAT(/' The masses of QM H atoms converted to 2.0 amu. ',
c     >       /' To turn off this default',
c     >        ' set nwpw:makehmass2 .false.')
 1135 FORMAT(/' The masses of QM H atoms converted to 2.0 amu. ',
     >       /' To turn off this default',
     >       /' nwpw',
     >       /'    makehmass2 off',
     >       /' end')
 1140 FORMAT(/' elements involved in the cluster:')
 1150 FORMAT(5X,I2,': ',A4,'  core charge:',F4.1,'  lmax=',I1)
 1151 FORMAT(5X,'        cutoff =',4F8.3)
 1152 FORMAT(12X,' highest angular component      : ',i2)
 1153 FORMAT(12X,' local potential used           : ',i2)
 1154 FORMAT(12X,' number of non-local projections: ',i2)
 1155 FORMAT(12X,' semicore corrections included  : ',
     >       F6.3,' (radius) ',F6.3,' (charge)')
 1156 FORMAT(12X,' aperiodic cutoff radius        : ',F6.3)
 1159 FORMAT(/' total charge=',F8.3)
 1160 FORMAT(/' atomic composition:')
 1170 FORMAT(7(5X,A4,':',I5))
 1180 FORMAT(/' initial position of ions:')
 1181 FORMAT(/' initial velocity of ions:')
 1190 FORMAT(5X, I4, A5  ,' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1191 FORMAT(5X, I4, A5  ,' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - fixed')
 1192 FORMAT(5X, I4, A5  ,' (',3F11.5,' )')
 1193 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - z fixed')
 1194 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,A)

 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1211 FORMAT(5X,'   number of constraints = ', I6,' ( DOF = ',I6,' )' )
 1219 FORMAT(/' number of electrons: spin up=',F6.2,'  down=',F6.2,A)
 1220 FORMAT(/' number of electrons: spin up=',I6,
     >        ' (',I4,' per task)',
     >        '  down=',I6,
     >        ' (',I4,' per task)',
     >        A)
 1221 FORMAT( ' number of orbitals : spin up=',I6,
     >        ' (',I4,' per task)',
     >        '  down=',I6,
     >        ' (',I4,' per task)',
     >        A)
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
 1229 FORMAT(/' small supercell:')
 1230 FORMAT(/' supercell:')
 1231 FORMAT(5x,' volume : ',F12.1)
 1232 FORMAT(5x,' lattice:    a=    ',f8.3,' b=   ',f8.3,' c=    ',f8.3,
     >      /5x,'             alpha=',f8.3,' beta=',f8.3,' gamma=',f8.3)
 1233 FORMAT(5x,' cell_name:  ',A)
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
  
 1260 FORMAT(5X,' Ewald summation: cut radius=',F8.2,'  and',I3)
 1261 FORMAT(5X,'                   madelung=',f14.8)
 1270 FORMAT(/' technical parameters:')
 1271 FORMAT(5x, ' translation constrained')
 1272 FORMAT(5x, ' rotation constrained')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'fictitious mass=',F10.1)
 1281 FORMAT(5X, ' maximum iterations =',I10,
     >           ' ( ',I4,' inner ',I6,' outer )')
 1290 FORMAT(5X, ' cooling/heatting rates: ',e12.5,' (psi)',2x,
     >                                       e12.5,' (ion)')
 
 1293 format(/' Pressure Output Generated         ')
 1294 format(/' Constant Energy Simulation                     ')
 1295 format(/' Nose-Hoover Simulation - Thermostat Parameters:')
 1296 format(5x, 'SA decay rates  =',e10.3,' (elc)',e10.3,' (ion)')
 1297 format(5x, 'link = ',I3,
     > ' Te =',f8.2,' Qe =',e10.3,' 2*pi/we=',e10.3,' Ee0=',e10.3)
 1298 format(5x, 'link = ',I3,
     > ' Tr =',f8.2,' Qr =',e10.3,' 2*pi/wr=',e10.3,' Er0=',e10.3)
 1299 format(//' Mulliken Analysis Output Generated            ')
 1300 FORMAT(//)
 1305 FORMAT(10X,'============ Car-Parrinello iteration ==============')
 1309 FORMAT(I8,2E19.10,2E14.5,2F9.1,3E11.3)
 1310 FORMAT(I8,2E19.10,2E14.5,F14.2,3E11.3)
 1311 format(100e19.10)
 1312 format(e14.6,i3)
 1313 format(3i4)
 1320 FORMAT(' number of electrons: spin up=',F11.5,'  down=',F11.5,A)
 1321 FORMAT(' total charge of system:',F11.5,A)
 1322 FORMAT('     plane-wave part:         ',F11.5,'       ',F11.5,A)
 1323 FORMAT('      augmented part:         ',F11.5,'       ',F11.5,A)
 1330 FORMAT(/' comparison between hamiltonian and lambda matrix')
 1331 FORMAT(/' Elements of Hamiltonian matrix (up/restricted)')
 1332 FORMAT(/' Elements of Hamiltonian matrix (down)')
 1340 FORMAT(I5,2I5,' H=',E16.7,', L=',E16.7,', H-L=',E16.7)
 1341 FORMAT(I5,2I5,' H=',E16.6)
 1350 FORMAT(/' orthonormality')
 1360 FORMAT(I5,2I5,E18.7)
 1370 FORMAT(I3)
 1380 FORMAT(' ''',a,'''',I4)
 1390 FORMAT(I3)
 1400 FORMAT(I3,3E18.8/3X,3E18.8)
 1410 FORMAT(10X,'=============  summary of results  =================')
 1420 FORMAT(/' final position of ions:')
 1421 FORMAT(/' final velocity of ions:')
 1430 FORMAT(/' total     energy    :',E19.10,' (',E15.5,'/ion)')
 1431 FORMAT(/' QM Energies')
 1432 FORMAT( '------------')
 1434 FORMAT(//' total paw energy    :',E19.10,' (',E15.5,'/ion)')
 1440 FORMAT( ' total orbital energy:',E19.10,' (',E15.5,'/electron)')
 1450 FORMAT( ' hartree   energy    :',E19.10,' (',E15.5,'/electron)')
 1455 FORMAT( ' SIC-hartree energy  :',E19.10,' (',E15.5,'/electron)')
 1456 FORMAT( ' SIC-exc-corr energy :',E19.10,' (',E15.5,'/electron)')
 1457 FORMAT( ' HF exchange energy  :',E19.10,' (',E15.5,'/electron)')
 1458 FORMAT( ' DFT+U     energy    :',E19.10,' (',E15.5,'/ion)')
 1459 FORMAT( ' Metadynamics energy :',E19.10,' (',E15.5,'/ion)')
 1460 FORMAT( ' exc-corr  energy    :',E19.10,' (',E15.5,'/electron)')
 1461 FORMAT( ' TAMD energy         :',E19.10,' (',E15.5,'/ion)')
 1470 FORMAT( ' ion-ion   energy    :',E19.10,' (',E15.5,'/ion)')
 1471 FORMAT(/' Kinetic energy (elc)    :',E19.10,' (',E15.5,'/elc)')
 1472 FORMAT( ' Kinetic energy (ion)    :',E19.10,' (',E15.5,'/ion)')
 1473 FORMAT( ' thermostat energy (elc) :',E19.10,' (',E15.5,'/elc)')
 1474 FORMAT( ' thermostat energy (ion) :',E19.10,' (',E15.5,'/ion)')
 1480 FORMAT(' Temperature :    ',F10.1,' K (ion)')
 1490 FORMAT('             :    ',F10.1,' K (c.o.m.)')
 1491 FORMAT(' Temperature :    ',F10.1,' K (elc)')
 1492 FORMAT(/' Vaverage  Eaverage :    ',E19.10, E19.10)
 1493 FORMAT( ' Vvariance Evariance:    ',E19.10, E19.10)
 1494 FORMAT( ' Cv - f*kb/(2*nion) :    ',E19.10)
 1499 FORMAT( ' K.S. SIC-hartree energy  :',E19.10,
     >        ' (',E15.5,'/electron)')
 1501 FORMAT( ' K.S. SIC-exc-corr energy :',E19.10,
     >        ' (',E15.5,'/electron)')
 
 1500 FORMAT(/' orbital energies:')
 1505 FORMAT( ' APC energy          :',E19.10,' (',E15.5,'/ion)')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1511 FORMAT(2(E18.7,' (',F8.3,'eV) occ=',F6.3))
 1600 FORMAT(/' Total PSPW energy   :',E19.10)

 1700 FORMAT(/' QM/MM-pol-vib/CAV Energies')
 1701 FORMAT( ' --------------------------')
 1702 FORMAT( ' LJ energy              :',E19.10)
 1703 FORMAT( ' Residual Coulomb energy:',E19.10)
 1704 FORMAT( ' MM Vibration energy    :',E19.10)
 1705 FORMAT( ' QM/MM coupling energy  :',E19.10, 
     >        ' (average=',E19.10,' variance=',E19.10,')'/)
 1706 FORMAT( ' Average and Variance of QM/MM coupling energy :',
     >       E19.10,E19.10)
 1707 FORMAT( ' QM/MM coupling param.  :',E19.10)

 1720 FORMAT(/' Dispersion energy   :',E19.10)


 1740 FORMAT(/' extended Born solvation energies:')
 1741 FORMAT(5x,' screen=(epsilon-1)/(epsilon):',F11.6)
 1745 FORMAT(5x,' solvation energy (w/o QM polarization) :',E19.10,
     >   ' (',F8.3,' kcal/mol)')
 1746 FORMAT(5x,' solvation energy (w/  QM polarization) :',E19.10,
     >   ' (',F8.3,' kcal/mol)')

 1800 FORMAT(/' Charge Field Energies')
 1801 FORMAT( ' ---------------------')
 1802 FORMAT( ' - Charge Field/Electron    :',E19.10)
 1803 FORMAT( ' - Charge Field/Ion         :',E19.10)
 1804 FORMAT( ' - Charge Field/Charge Field:',E19.10)
 1805 FORMAT( ' Charge Field Energy        :',E19.10)

 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(luout,9010) ierr
      call Parallel2d_Finalize()
      call Parallel_Finalize()

      cpmdv5 = value
      return
      END

c
c Now in nwpw/utilities
c
c      subroutine print_elapsed_time(autime)
c      implicit none
c      real*8 autime
c
c#include "stdio.fh"
c
c      real*8 sectime
c
c      sectime = autime*2.41889d-17
c
c      if (sectime.lt.1.0d-12) then
c         write(luout,1800) (sectime/1.0d-15)," fs"
c      else if (sectime.lt.1.0d-9) then
c         write(luout,1800) (sectime/1.0d-12)," ps"
c      else 
c         write(luout,1800) (sectime/1.0d-9 )," ns"
c      end if
c
c      return
c 1800 format(//' Elapsed time of simulation was',F8.3,A)
c      end


