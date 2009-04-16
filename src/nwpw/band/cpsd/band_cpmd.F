*
* $Id: cpmdv5.F,v 1.76 2009-01-31 01:18:00 bylaska Exp $
*

***********************************************************************
*                                                                     *
*                         band_cpmd                                   *
*                                                                     *
***********************************************************************

      logical function band_cpmd(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
      
      logical value

      real*8 kb
      parameter (kb=3.16679d-6)
c      real*8 autoatm
c      parameter (autoatm =290.360032539d6)


      
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
      real*8 icharge
      integer ispin
      integer ne(2),n1(2),n2(2),nemax,neall,neq(2),nemaxq
      real*8  en(2)
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
      real*8  E(30),eke,eave,evar,cv

*     real*8  eig(2*nemax)
*     real*8  hml(2*nemax*nemax)
*     real*8  lmd(2*nemax*nemax)
      integer eig(2),hml(2),lmd(2)

*     **** psi smearing block ****
      logical fractional
      integer smearoccupation,smeartype
      real*8 smearfermi(2),smearcorrection,smearkT

*     **** error variables ****
      integer ierr

*     **** local variables ****
      logical verlet,mulliken,SA,found,calc_pressure
      logical lprint,mprint,hprint
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
      integer nfft3,mapping,mapping1d
      parameter (nfft3=32)
      character*30 filename
      character*255 full_filename
  


*     **** external functions ****
      real*8      psp_zv,psp_rc,ewald_rcut,ion_amass
      real*8      ewald_mandelung
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg
      integer     ewald_ncut,ewald_nshl3d
      integer     psp_lmmax,psp_lmax,psp_locp,ion_nkatm
      character*4 ion_atom,ion_aname
      external    psp_zv,psp_rc,ewald_rcut,ion_amass
      external    ewald_mandelung
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg
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

      integer  pack_nwave_all
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      external pack_nwave_all
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm

      character*12 control_boundry
      external     control_boundry

      logical      psp_semicore,pspw_qmmm_found
      real*8       psp_rcore,psp_ncore,psp_rlocal
      external     psp_semicore,pspw_qmmm_found
      external     psp_rcore,psp_ncore,psp_rlocal

      logical  control_Nose,control_Mulliken
      external control_Nose,control_Mulliken

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

      logical  control_SA,control_Fei,pspw_SIC,pspw_HFX,control_pressure
      real*8   control_SA_decay
      external control_SA,control_Fei,pspw_SIC,pspw_HFX,control_pressure
      external control_SA_decay

      integer  control_np_orbital,control_mapping,control_mapping1d
      external control_np_orbital,control_mapping,control_mapping1d


      logical  control_translation,control_rotation,control_balance
      external control_translation,control_rotation,control_balance
     
      logical  Dneall_m_allocate,Dneall_m_free
      external Dneall_m_allocate,Dneall_m_free

      real*8   Dneall_m_value
      external Dneall_m_value




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

      if (.not.control_read(5,rtdb))
     >   call errquit('band_sd:error reading control',0,DISK_ERR)

      lprint = ((taskid.eq.MASTER).and.(control_print(print_low)))
      mprint = ((taskid.eq.MASTER).and.(control_print(print_medium)))
      hprint = ((taskid.eq.MASTER).and.(control_print(print_high)))

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
      end if
      
      call Parallel3d_Init(control_np_orbital())
      call Parallel3d_np_i(np_i)
      call Parallel3d_np_j(np_j)
      call Parallel3d_np_j(np_k)

      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0
      mapping = control_mapping()

*     **** initialize D3dB data structure ****
      call C3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call C3dB_nfft3d(1,nfft3d)

*     **** initialize psi_data ****
      call cpsi_data_init(20)

*     **** initialize lattice and packing data structure ****
      call lattice_init()
      call c_G_init()
      call brillioun_init()
      call Cram_Init()
      call C3dB_pfft_init()

*     ***** Initialize double D3dB data structure ****
      if ((control_gga().ge.10).and.(control_gga().lt.100)) then
         call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
         call G_init()
         call mask_init()
      end if

*     **** read ions ****
      value = ion_init(rtdb)
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

*     **** allocate psp data structure and read in psedupotentials into it ****
      call cpsp_init()
      call cpsp_readall()
      if (cpsp_semicore(0)) call c_semicore_check()


*     **** initialize G,mask,ke,and coulomb data structures ****
      call cstrfac_init()
      call cke_init()
      call c_coulomb_init()
      call ewald_init()

c*     ***** allocate psi2, psi1, and psi0 wavefunctions ****
c      call psi_get_ne_occupation(ispin,ne,smearoccupation)
c      if (smearoccupation.gt.0) then
c         fractional = .true.
c      else
c         fractional = .false.
c      end if
c      mapping1d = control_mapping1d()
c      call Dne_init(ispin,ne,mapping1d)
c      call Dneall_neq(neq)
c      nemaxq = neq(1)+neq(2)
c
c      value = MA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
c     >                     'psi2',psi2(2),psi2(1))
c      value = value.and.
c     >        MA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
c     >                     'psi1',psi1(2),psi1(1))
c      value = value.and.
c     >        MA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
c     >                     'psi0',psi0(2),psi0(1))
c      if (fractional) then
c      value = value.and.
c     >        MA_alloc_get(mt_dbl,(ne(1)+ne(2)),'occ0',occ0(2),occ0(1))
c      value = value.and.
c     >        MA_alloc_get(mt_dbl,(ne(1)+ne(2)),'occ1',occ1(2),occ1(1))
c      value = value.and.
c     >        MA_alloc_get(mt_dbl,(ne(1)+ne(2)),'occ2',occ2(2),occ2(1))
c      end if
c      if (.not.value) call errquit('out of heap memory',0, MA_ERR)
c
c
c
c*     *****  read psi2 wavefunctions ****
c      call psi_read(ispin,ne,dcpl_mb(psi2(1)),
c     >              smearoccupation,dbl_mb(occ2(1)))
c
c
c*     **** move  wavefunction velocities ****
c      if (control_new_vpsi()) then
c        call v_psi_delete()
c      end if
c
c
c*     **** generate initial wavefunction velocities if it does not exist ****
c      if (.not.v_psi_filefind()) then
c        call v_psi_new(ispin,ne)
c      end if
c
c
c*     *****  read psi0 wavefunctions ****
c      call dcopy(2*npack1*(neq(1)+neq(2)),0.0d0,0,dcpl_mb(psi1(1)),1)
c      call v_psi_read(ispin,ne,dcpl_mb(psi1(1)))
c      n1(1) = 1
c      n2(1) = ne(1)
c      n1(2) = ne(1)+1
c      n2(2) = ne(1)+ne(2)
c      nemax = ne(1)+ne(2)
c
c*     **** allocate other variables *****
c      value = MA_alloc_get(mt_dbl,(2*nemax),'eig',eig(2),eig(1))
c      value = value.and.Dneall_m_allocate(0,hml)
c      value = value.and.Dneall_m_allocate(0,lmd)
c
c      value = value.and.
c     >        MA_alloc_get(mt_dbl,(4*nfft3d),
c     >                     'dn',dn(2),dn(1))
c      value = value.and.
c     >        MA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
c     >                     'Hpsi',Hpsi(2),Hpsi(1))
c      value = value.and.
c     >        MA_alloc_get(mt_dcpl,nfft3d*(neq(1)+neq(2)),
c     >                     'psir',psir(2),psir(1))
c      if (.not.value) call errquit('out of heap memory',0, MA_ERR)
c
c
c
c*     **** read ions ****
c      value = ion_init(rtdb)
c
c*     **** initialize FixIon constraint ****
c      call ion_init_FixIon(rtdb)
c
c*     **** allocate psp data structure and read in psedupotentials into it ****
c      call psp_init()
c      call psp_readall()
c      if (psp_semicore(0)) call semicore_check()
c
c
c*     **** initialize G,mask,ke,and coulomb data structures ****
c      call ke_init()
c      if (control_version().eq.3) call coulomb_init()
c      if (control_version().eq.4) call coulomb2_init()
c      call strfac_init()
c      if (control_version().eq.3) call ewald_init()
c      
c
c*     ******************************
c*     **** scaling psi velocity ****
c*     ******************************
c      call dcopy(2*(neq(1)+neq(2))*npack1,dcpl_mb(psi1(1)),1,
c     >                                  dcpl_mb(psi0(1)),1)
c      call dscal(2*(neq(1)+neq(2))*npack1,control_rte(),
c     >           dcpl_mb(psi1(1)),1)
c      eke0 = 0.0d0
c      eke1 = 0.0d0
c      do i=1,(neq(1)+neq(2))
c         call Pack_cc_dot(1,dcpl_mb(psi0(1)+(i-1)*npack1),
c     >                      dcpl_mb(psi0(1)+(i-1)*npack1),
c     >                     sum)
c         eke0 = eke0 + sum
c         call Pack_cc_dot(1,dcpl_mb(psi1(1)+(i-1)*npack1),
c     >                      dcpl_mb(psi1(1)+(i-1)*npack1),
c     >                    sum)
c         eke1 = eke1 + sum
c      end do
c
c      call D1dB_SumAll(eke0)
c      call D1dB_SumAll(eke1)
c      eke0 = control_fake_mass()*eke0
c      eke1 = control_fake_mass()*eke1
c      call ion_init_ke(ekg,eki0,eki1)
c
c
c*     **** Initialize thermostats ****
c      if (control_Nose()) then
c         call ke_ave(ispin,neq,dcpl_mb(psi2(1)),w,
c     >               fractional,dbl_mb(occ2(1)))
c         call Nose_Init((ne(1)+ne(2)),w)
c      end if
c
c
c*     **** Initialize simulated annealing ****
c      SA       = .false.
c      Te_init  = 0.0d0
c      Tr_init  = 0.0d0
c      sa_alpha(1) = 1.0d0
c      sa_alpha(2) = 1.0d0
c      if (control_SA()) then
c         if (control_Nose()) then
c            SA          = .true.
c            sa_decay(1) = control_SA_decay(1)
c            sa_decay(2) = control_SA_decay(2)
c            Te_init     = control_Nose_Te()
c            Tr_init     = control_Nose_Tr()
c         else
c            dt = control_time_step()
c            SA          = .false.
c            sa_decay(1) = control_SA_decay(1)
c            sa_decay(2) = control_SA_decay(2)
c            sa_alpha(1) = dexp( -(dt/control_SA_decay(1)) ) 
c            sa_alpha(2) = dexp( -(dt/control_SA_decay(2)) ) 
c         end if
c      end if
c
c
c*     **** initialize QM/MM ****
c      call pspw_qmmm_init(rtdb)
c
c*     **** initialize dplot ****
c      call dplot_iteration_init()
c
cc*     **** initialize frac_occ data structure ****
cc      call frac_occ_init(rtdb,ispin,ne)
c
c*     **** initialize SIC and HFX ****
c      call pspw_init_SIC(rtdb,ne)
c      call pspw_init_HFX(rtdb,ispin,ne)
c
c
c
c*     **** initialize pressure ****
c      calc_pressure = control_pressure().and.(control_version().eq.3)
c      pressure      = 0.0d0
c      p1            = 0.0d0
c      p2            = 0.0d0
c      if (calc_pressure) then
c         call psp_stress_init()
c         call psp_stress_readall()
c      end if
c
c
c
c
c
c
c
c
c*                |**************************|
c******************   summary of input data  **********************
c*                |**************************|
c      call center_geom(cx,cy,cz)
c      call center_mass(gx,gy,gz)
c      call center_v_geom(vcx,vcy,vcz)
c      call center_v_mass(vgx,vgy,vgz)
c      mulliken = control_Mulliken()
c
c      if (taskid.eq.MASTER) then
c         write(6,1111) np
c         write(6,1117) np_i,np_j
c         if (mapping.eq.1) write(6,1112)
c         if (mapping.eq.2) write(6,1113)
c         if (control_balance()) then
c           write(6,1114)
c         else
c           write(6,1116)
c         end if
c
c         write(6,1115)
c         write(6,1121) control_boundry(),control_version()
c         if (ispin.eq.1) write(6,1130) 'restricted'
c         if (ispin.eq.2) write(6,1130) 'unrestricted'
c         IF (control_gga().eq.-1) THEN
c            write(6,1131) 'off'
c         ELSE IF (control_gga().eq.0) THEN
c            write(6,1131) 'LDA (Vosko et al) parameterization'
c         ELSE IF (control_gga().eq.10) THEN
c            write(6,1131)
c     >      'PBE96 (White and Bird) parameterization'
c         ELSE IF (control_gga().eq.11) THEN
c            write(6,1131)
c     >      'BLYP (White and Bird) parameterization'
c         ELSE IF (control_gga().eq.12) THEN
c            write(6,1131)
c     >      'revPBE (White and Bird) parameterization'
c
c         ELSE IF (control_gga().eq.110) THEN
c            write(6,1131)
c     >      'PBE0 (White and Bird) parameterization'
c         ELSE IF (control_gga().eq.112) THEN
c            write(6,1131)
c     >      'revPBE0 (White and Bird) parameterization'
c         ELSE IF (control_gga().eq.113) THEN
c            write(6,1131)
c     >      'BNL (White and Bird) parameterization'
c
c         ELSE IF (control_gga().eq.200) THEN
c            write(6,1131) 'Hartree-Fock'
c         ELSE
c            write(6,1131) 'unknown parameterization'
c            call errquit('bad exchange_correlation',0, INPUT_ERR)
c         END IF
c         if (fractional) write(6,1132)
c         call pspw_print_SIC(6)
c         call pspw_print_HFX(6)
c         write(6,1140)
c         do ia = 1,ion_nkatm()
c           write(6,1150) ia,ion_atom(ia),
c     >                    psp_zv(ia),psp_lmax(ia)
c           write(6,1152) psp_lmax(ia)
c           write(6,1153) psp_locp(ia)
c           write(6,1154) psp_lmmax(ia)
c           if (control_version().eq.4) write(6,1156) psp_rlocal(ia)
c           if (psp_semicore(ia))
c     >         write(6,1155) psp_rcore(ia),psp_ncore(ia)
c           write(6,1151) (psp_rc(i,ia),i=0,psp_lmax(ia))
c         end do
c
c         icharge = -(ne(1)+ne(ispin))
c         en(1)     = ne(1)
c         en(ispin) = ne(ispin)
c         if (fractional) then
c            icharge = 0.0d0
c            do ms=1,ispin
c            en(ms) =0.0
c            do i=n1(ms),n2(ms)
c              icharge = icharge - (3-ispin)*dbl_mb(occ2(1)+i-1)
c              en(ms) = en(ms) + dbl_mb(occ2(1)+i-1)
c            end do
c            end do
c         end if
c
c         do ia=1,ion_nkatm()
c           icharge = icharge + ion_natm(ia)*psp_zv(ia)
c         end do
c         write(6,1159) icharge
c
c         write(6,1160)
c         write(6,1170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
c         write(6,1180)
c         do I=1,ion_nion()
c           if (ion_q_FixIon(I)) then
c           write(6,1191) I,ion_aname(I),
c     >                    (ion_rion(K,I),K=1,3),ion_amass(i)/1822.89d0
c           else
c           write(6,1190) I,ion_aname(I),
c     >                    (ion_rion(K,I),K=1,3),ion_amass(i)/1822.89d0
c           end if
c         end do
c         write(6,1200) cx,cy,cz
c         write(6,1210) gx,gy,gz
c
c
c         write(6,1181)
c         write(6,1192) (I,ion_aname(I),
c     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
c         write(6,1200) vcx,vcy,vcz
c         write(6,1210) vgx,vgy,vgz
c
c
c         if (fractional) then
c           write(6,1219) en(1),en(ispin),' (   fractional)'
c           write(6,1221) ne(1),neq(1),
c     >                   ne(ispin),neq(ispin),' (fourier space)'
c         else
c           write(6,1220) ne(1),neq(1),
c     >                   ne(ispin),neq(ispin),' (fourier space)'
c           write(6,1221) ne(1),neq(1),
c     >                   ne(ispin),neq(ispin),' (fourier space)'
c         end if
c         write(6,1230)
c         write(6,1241) lattice_unita(1,1),
c     >                 lattice_unita(2,1),
c     >                 lattice_unita(3,1)
c         write(6,1242) lattice_unita(1,2),
c     >                 lattice_unita(2,2),
c     >                 lattice_unita(3,2)
c         write(6,1243) lattice_unita(1,3),
c     >                 lattice_unita(2,3),
c     >                 lattice_unita(3,3)
c         write(6,1244) lattice_unitg(1,1),
c     >                 lattice_unitg(2,1),
c     >                 lattice_unitg(3,1)
c         write(6,1245) lattice_unitg(1,2),
c     >                 lattice_unitg(2,2),
c     >                 lattice_unitg(3,2)
c         write(6,1246) lattice_unitg(1,3),
c     >                 lattice_unitg(2,3),
c     >                 lattice_unitg(3,3)
c         write(6,1231) lattice_omega()
c         write(6,1250) lattice_ecut(),ngrid(1),ngrid(2),ngrid(3),
c     >                 pack_nwave_all(0),pack_nwave(0)
c         write(6,1251) lattice_wcut(),ngrid(1),ngrid(2),ngrid(3),
c     >                 pack_nwave_all(1),pack_nwave(1)
c         if (control_version().eq.3) then
c         write(6,1260) ewald_rcut(),ewald_ncut()
c         write(6,1261) ewald_mandelung()
c         end if
c         write(6,1270)
c         if (.not.control_translation()) write(6,1271)
c         if (.not.control_rotation())    write(6,1272)
c         write(6,1280) control_time_step(),control_fake_mass()
c         write(6,1290) control_rte(),control_rti()
c         write(6,1222) eke0,eki0,ekg
c         write(6,1223) eke1,eki1
c         write(6,1224) (eke1-eke0),(eki1-eki0)
c         if (control_Nose()) then
c           write(6,1295)
c           do i=1,Nose_Mchain()
c             write(6,1297) i,control_Nose_Te(),Nose_Qe(i),
c     >                     Nose_Pe(i),Nose_Ee0(i)
c           end do
c           do i=1,Nose_Nchain()
c             write(6,1298) i,control_Nose_Tr(),Nose_Qr(i),
c     >                   Nose_Pr(i),Nose_Er0(i)
c           end do
c         else
c           write(6,1294)
c         end if
c        if (calc_pressure) write(6,1293)
c        if (control_SA()) then
c           write(6,1296) sa_decay(1),sa_decay(2)
c         end if
c
c
c         if (mulliken) write(6,1299)
c         write(6,1300)
c         write(6,1305)
c         call util_flush(6)
c      end if
c
c*                |***************************|
c******************     start iterations      **********************
c*                |***************************|
c*     **** open xyz and MOTION file ****
c      call xyz_init()          ! unit=18
c      call MOTION_init(rtdb)   ! unit=19
c
c*     *** fei io ****
c      call fei_init(rtdb)
c
c*     **** ecce print ****
c      call ecce_print_module_entry('task Car-Parrinello')
c      !call ecce_print_module_entry('driver')
c      call movecs_ecce_print_off()
c
c
c
c
c*     ************************************
c*     **** open up other MOTION files ****
c*     ************************************
c
c
c*     **** open EMOTION file ****
c      if (.not.rtdb_cget(rtdb,'cpmd:emotion_filename',1,filename))
c     >  call util_file_prefix('emotion',filename)
c      call util_file_name_noprefix(filename,.false.,
c     >                             .false.,
c     >                    full_filename)
c      if (taskid.eq.MASTER) then
c
c         emotion_time_shift = 0.0d0
c         icount_shift       = 0
c         inquire(file=full_filename,exist=found)
c         if (found) then
c           open(unit=31,file=full_filename,form='formatted',
c     >          status='old')
c           do while (found)
c           read(31,*,end=100) emotion_time_shift,w,sum
c           E(25) = E(25) + sum                          !*** take care of running sums ***
c           E(26) = E(26) + sum*sum
c           icount_shift = icount_shift + 1
c           end do
c  100      continue
c#if defined(FUJITSU_SOLARIS) || defined(PSCALE) || defined(__crayx1)
c           backspace 31
c#endif
c         else
c           open(unit=31,file=full_filename,form='formatted',
c     >          status='new')
c         end if
c      end if
c
c
c*     **** open EIGMOTION file ****
c      if (mulliken) then
c        if (.not.rtdb_cget(rtdb,'cpmd:eigmotion_filename',1,filename))
c     >    call util_file_prefix('eigmotion',filename)
c      call util_file_name_noprefix(filename,.false.,
c     >                             .false.,
c     >                    full_filename)
c      if (taskid.eq.MASTER) 
c     >   open(unit=32,file=full_filename,form='formatted')
c      end if
c
c*     **** open HMOTION file ****
c      if (mulliken) then
c       if (.not.rtdb_cget(rtdb,'cpmd:hmotion_filename',1,filename))
c     >  call util_file_prefix('hmotion',filename)
c      call util_file_name_noprefix(filename,.false.,
c     >                             .false.,
c     >                    full_filename)
c      if (taskid.eq.MASTER) 
c     >   open(unit=34,file=full_filename,form='formatted')
c      end if
c
c*     **** open OMOTION file ****
c      if (mulliken) call Orb_Init(rtdb,ispin,ne) !unit=33
c
c*     **** write initial position to xyz data ****
c      call xyz_write()
c
c*     ***** first step using velocity ****
c      verlet = .false.
c      call inner_loop_md(verlet,sa_alpha,ispin,ne,neq,
c     >             npack1,nfft3d,nemaxq,
c     >             dcpl_mb(psi0(1)),
c     >             dcpl_mb(psi1(1)),
c     >             dcpl_mb(psi2(1)),
c     >             dbl_mb(dn(1)),
c     >             1,0,E,
c     >             dbl_mb(hml(1)),dbl_mb(lmd(1)),
c     >             dcpl_mb(Hpsi(1)),dcpl_mb(psir(1)),
c     >             calc_pressure,pressure,p1,p2,
c     >             fractional,
c     >             dbl_mb(occ0(1)),dbl_mb(occ1(1)),dbl_mb(occ2(1)))
c
c
c      if (taskid.eq.MASTER) call current_second(cpu2)
c      if ((taskid.eq.MASTER).and.(.not.calc_pressure)) 
c     >   CALL nwpw_message(6)
c      if ((taskid.eq.MASTER).and.(calc_pressure)) 
c     >   CALL nwpw_message(9)
c      
c      it_in  = control_it_in()
c      it_out = control_it_out()
c      icount = 0
c      verlet = .true.
c      eke    = 0.0d0
c      if (it_out.lt.1) goto 102
c
c
c
c      dt = control_time_step()
c
c      Te_new = Te_init
c      Tr_new = Tr_init
c  101 continue
c         icount = icount + 1
c         call inner_loop_md(verlet,sa_alpha,ispin,ne,neq,
c     >             npack1,nfft3d,nemaxq,
c     >             dcpl_mb(psi0(1)),
c     >             dcpl_mb(psi1(1)),
c     >             dcpl_mb(psi2(1)),
c     >             dbl_mb(dn(1)),
c     >             it_in,((icount-1)*it_in),
c     >             E,
c     >             dbl_mb(hml(1)),dbl_mb(lmd(1)),
c     >             dcpl_mb(Hpsi(1)), dcpl_mb(psir(1)),
c     >             calc_pressure,pressure,p1,p2,
c     >             fractional,
c     >             dbl_mb(occ0(1)),dbl_mb(occ1(1)),dbl_mb(occ2(1)))
c
c         eke = eke + E(3)
c
c         !**** calculate pressure ****
c
c         if (taskid.eq.MASTER) then 
c
c           if (calc_pressure) then
c             if (SA) then
c             write(6,1309) icount*it_in,E(1),E(2),E(3),E(4),
c     >                     Te_new,Tr_new,pressure
c             else
c             write(6,1310) icount*it_in,E(1),E(2),E(3),E(4),
c     >                     ion_Temperature(),pressure
c             end if
c           else
c             if (SA) then
c             write(6,1309) icount*it_in,E(1),E(2),E(3),E(4),
c     >                     Te_new,Tr_new
c             else
c             write(6,1310) icount*it_in,E(1),E(2),E(3),E(4),
c     >                     ion_Temperature()
c             end if
c           end if
c           call util_flush(6)
c
c*          **** write out EMOTION data ****
c           eave = E(25)/dble(icount+icount_shift)
c           evar = E(26)/dble(icount+icount_shift)
c           evar = evar - eave*eave
c           if (control_Nose()) then
c             write(31,1311) icount*it_in*dt + emotion_time_shift,
c     >                    e(1),e(2),e(3),e(4),e(5),e(6),
c     >                    e(7),e(8),e(9),e(10),
c     >                    eave,evar,ion_Temperature(),pressure
c           else
c             write(31,1311) icount*it_in*dt + emotion_time_shift,
c     >                    e(1),e(2),e(3),e(4),e(5),e(6),
c     >                    e(7),e(8),
c     >                    eave,evar,ion_Temperature(),pressure
c           end if
c           call util_flush(31)
c
c*          **** write out EIGMOTION data -diagonal hml matrix ****
c           if (mulliken) then
c           write(32,1311) icount*it_in*dt,
c     >       (( dbl_mb(hml(1)+ii-1+(ii-1)*ne(1)+(ms-1)*ne(1)*ne(1)), 
c     >         ii=1,ne(ms)),ms=1,ispin)
c           call util_flush(32)
c           end if
c
c*          **** write out HMOTION data - hml matrix ****
c           if (mulliken) then
c           write(34,1312) icount*it_in*dt,ispin
c           do ms=1,ispin
c             write(34,1313) ms,ne(ms),ne(ms)
c             do ii=1,ne(ms)
c               write(34,1311) 
c     >         (dbl_mb(hml(1)+ii-1+(jj-1)*ne(1)+(ms-1)*ne(1)*ne(1)), 
c     >          jj=1,ne(ms))
c             end do
c           end do
c           call util_flush(34)
c           end if
c
c         end if
c
c
c*        **** write xyz data ****
c         call xyz_write()
c         call MOTION_write((icount*it_in*control_time_step()))
c
c*        **** write OMOTION data ****
c         if (mulliken) call Orb_Write(dcpl_mb(psi1(1)))
c
c*        **** update thermostats using SA decay ****
c         if (SA) then
c           t1 = icount*it_in*dt/sa_decay(1)
c           t2 = icount*it_in*dt/sa_decay(2)
c           Te_new = Te_init*dexp(-t1)
c           Tr_new = Tr_init*dexp(-t2)
c           call Nose_reset_T(Te_new,Tr_new)
c         end if
c
c
c*        **** exit early ****
c         if (control_out_of_time()) then
c            if (taskid.eq.MASTER) 
c     >       write(6,*) ' *** out of time. iteration terminated'
c            go to 102
c         end if
c      if (icount.lt.it_out) go to 101
c      if (taskid.eq.MASTER) 
c     > write(6,*) '*** arived at the Maximum iteration.   terminated.'
c 
c*::::::::::::::::::::  end of iteration loop  :::::::::::::::::::::::::
c
c  102 continue
c
c*     **** close xyz and MOTION files ****
c      call xyz_end()
c      call MOTION_end()
c      if (taskid.eq.MASTER) then
c        close(unit=31)
c        close(unit=32)
c        close(unit=34)
c      end if
c
c*     *** close fei io ****
c      call fei_end()
c
c*     **** close OMOTION file ****
c      if (mulliken) call Orb_End()
c
c*     **** ecce print ****
c      !call ecce_print_module_exit('driver', 'ok')
c      call ecce_print_module_exit('task Car-Parrinello', 'ok')
c
c
c*     **** finalize pressure ****
c      if (calc_pressure) then
c         call psp_stress_end()
c      end if
c
c
c      if (taskid.eq.MASTER) CALL nwpw_message(3)
c      if (taskid.eq.MASTER) call current_second(cpu3)
c
c
c*         |****************************************|
c*********** produce CHECK file and diagonalize hml *****************
c*         |****************************************|
c
c*     **** produce CHECK FILE ****
c      if (taskid.eq.MASTER) then
c         call util_file_name('CHECK',.true.,
c     >                               .false.,
c     >                       full_filename)
c         open(unit=17,file=full_filename,form='formatted')
c      end if
c
c*     **** check total number of electrons ****
c      do ms =1,ispin
c         call D3dB_r_dsum(1,dbl_mb(dn(1)+(ms-1)*n2ft3d),sumall)
c         en(ms) = sumall*lattice_omega()
c     >             /dble(ngrid(1)*ngrid(2)*ngrid(3))
c      end do
c      if (taskid.eq.MASTER) then
c         write(17,1320) (en(ms),ms=1,ispin)
c      end if
c
c*     **** comparison between hamiltonian an lambda matrix ****
c      if (taskid.eq.MASTER) write(17,1330)
c      do ms=1,ispin
c         do i=1,ne(ms)
c         do j=1,ne(ms)
c            w   = Dneall_m_value(0,ms,i,j,dbl_mb(hml(1)))
c            sum = Dneall_m_value(0,ms,i,j,dbl_mb(lmd(1)))
c
c            if (taskid.eq.MASTER)
c     >      write(17,1340) ms,i,j,w,sum,w-sum
c
c         end do
c         end do
c      end do
c
c*     **** check orthonormality ****
c      if (taskid.eq.MASTER) then
c         write(17,1350)
c      end if
c
c      call Dneall_ffm_Multiply(0,dcpl_mb(psi1(1)),
c     >                           dcpl_mb(psi1(1)),npack1,
c     >                           dbl_mb(lmd(1)))
c      do ms=1,ispin
c         do j=1,ne(ms)
c         do i=j,ne(ms)
c            w  = Dneall_m_value(0,ms,i,j,dbl_mb(lmd(1)))
c            if (taskid.eq.MASTER) write(17,1360) ms,i,j,w
c         end do
c         end do
c      end do
c
c*     **** close check file ****
c      if (taskid.eq.MASTER) then
c         close(17)
c      end if
c
c
c
c*     ***** do not diagonalize the hamiltonian matrix *****
c      if (pspw_SIC()) then
c        call dcopy(2*npack1*nemaxq,
c     >             dcpl_mb(psi1(1)),1,
c     >             dcpl_mb(psi2(1)),1)
c      
c*     ***** diagonalize the hamiltonian matrix *****
c      else
c
cc         if (fractional) then
cc           call Dneall_m_HmltimesSA(0,dbl_mb(hml(1)),dbl_mb(fweight(1)))
cc         end if
c
c         call Dneall_m_diagonalize(0,dbl_mb(hml(1)),
c     >                               dbl_mb(eig(1)),.false.)
c
c
cc         if (fractional) then
cc         do ii=1,ne(ms)
cc           dbl_mb(eig(1)+(ii-1)+(ms-1)*n)
cc     >       =dbl_mb(eig(1)+(ii-1)+(ms-1)*n)
cc     >       /dbl_mb(fweight(1)+(ii-1)+(ms-1)*n)
cc         end do
cc         end if
c
c
c*        *** rotate current psi ***
c         call Dneall_fmf_Multiply(0,dcpl_mb(psi1(1)),npack1,
c     >                              dbl_mb(hml(1)),  1.0d0,
c     >                              dcpl_mb(psi2(1)),0.0d0)
c
c
c*        *** rotate current v_psi ***
c         call dcopy(2*npack1*nemaxq,dcpl_mb(psi0(1)),1,
c     >                              dcpl_mb(psi1(1)),1)
c
c         call Dneall_fmf_Multiply(0,dcpl_mb(psi1(1)),npack1,
c     >                              dbl_mb(hml(1)),  1.0d0,
c     >                              dcpl_mb(psi0(1)),0.0d0)
c
c      end if
c
c
c
c*                |***************************|
c****************** report summary of results **********************
c*                |***************************|
c      call center_geom(cx,cy,cz)
c      call center_mass(gx,gy,gz)
c      call center_v_geom(vcx,vcy,vcz)
c      call center_v_mass(vgx,vgy,vgz)
c
c      if (taskid.eq.MASTER) then
c         call print_elapsed_time(icount*it_in*dt)
c         write(6,1300)
c         write(6,1410)
c         write(6,1420)
c         do I=1,ion_nion()
c           if (ion_q_FixIon(I)) then
c           write(6,1191) I,ion_aname(I),(ion_rion(k,i),K=1,3),
c     >                   ion_amass(I)/1822.89d0
c           else
c           write(6,1190) I,ion_aname(I),(ion_rion(k,i),K=1,3),
c     >                   ion_amass(I)/1822.89d0
c           end if
c         end do
c         write(6,1200) cx,cy,cz
c         write(6,1210) gx,gy,gz
c
c
c         write(6,1421)
c         write(6,1192) (I,ion_aname(I),
c     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
c         write(6,1200) vcx,vcy,vcz
c         write(6,1210) vgx,vgy,vgz
c
c
c         write(6,*)
c         write(6,1320) en(1),en(ispin),' (real space)'
c         write(6,1430) E(2),E(2)/ion_nion()
c         write(6,1440) E(5),E(5)/n2(ispin)
c         write(6,1450) E(6),E(6)/n2(ispin)
c         write(6,1460) E(7),E(7)/n2(ispin)
c         if (pspw_SIC()) then
c           write(6,1455) E(16),E(16)/n2(ispin)
c           write(6,1456) E(17),E(17)/n2(ispin)
c         end if
c         if (pspw_HFX()) then
c           write(6,1457) E(20),E(20)/n2(ispin)
c         end if
c         write(6,1470) E(8),E(8)/ion_nion()
c         write(6,1471) E(3),E(3)/n2(ispin)
c         write(6,1472) ion_ke(),ion_ke()/ion_nion()
c
c
c         if (pspw_qmmm_found()) then
c            write(6,1700)
c            write(6,1701)
c            write(6,1702) E(11)
c            write(6,1703) E(12)
c            write(6,1704) E(13)
ccc            write(6,1705) E(14)
c            write(6,1706) E(15)
c         end if
c
c
c         if (control_Nose()) then
c           write(6,1473) E(9),E(9)/n2(ispin)
c           write(6,1474) E(10),E(10)/ion_nion()
c         end if
c         write(6,1226) E(3),ion_ke(),ion_com_ke()
c         eke = eke/dble(it_out)
c         eke = 2.0d0*eke/kb/(ne(1)+ne(ispin))/pack_nwave_all(1)
c         !eke = 2.0d0*eke/kb/(ne(1)+ne(ispin))
c
c*       **** write out Temperatures ****
c         write(6,1491) eke
c         eki0 = ion_Temperature()
cc         if (pspw_qmmm_found()) then
cc            eki1 =pspw_qmmm_Temperature()
cc            sum = ion_nion() + pspw_qmmm_nion() - 2.0d0
cc            eki0 = eki0*((ion_nion()-2.0d0)/sum)
cc     >           + eki1*((pspw_qmmm_nion()-2.0d0)/sum)
cc         end if
c         write(6,1480) eki0
c         write(6,1490) ion_com_Temperature()
c
c         eave = E(25)/dble(icount+icount_shift)
c         evar = E(26)/dble(icount+icount_shift)
c         evar = evar - eave*eave
c         cv = (evar)/(kb*ion_Temperature()**2)
c         cv = cv/dble(ion_nion())
c         write(6,1492) eave
c         write(6,1493) evar
c         write(6,1494) cv
c
c*        **** write out diagonal <psi|H|psi> matrix ****
c         if (pspw_SIC()) then
c
c           n = ne(1)
c           nn = n*n
c           do ms=1,ispin
c             if (ms.eq.1) write(6,1331)
c             if (ms.eq.2) write(6,1332)
c             !*** call Gainsville matrix output ***
c             call output(dbl_mb(hml(1)+(ms-1)*nn),
c     >                    1,ne(ms),1,ne(ms),
c     >                    n,n,1)
c           end do
c
c
c*        **** write out KS eigenvalues ****
c         else
c         write(6,1500)
c         NN=NE(1)-NE(2)
c         EV=27.2116d0
c
c         if (fractional) then
c           do i=1,NN
c             write(6,1511) dbl_mb(EIG(1)+i-1),
c     >                     dbl_mb(EIG(1)+i-1)*EV,
c     >                     dbl_mb(occ2(1)+i-1)
c           end do
c           do i=1,ne(2)
c             write(6,1511) dbl_mb(EIG(1)+i-1+NN),
c     >                     dbl_mb(EIG(1)+i-1+NN)*EV,
c     >                     dbl_mb(occ2(1)+i-1+NN),
c     >                     dbl_mb(EIG(1)+i-1+n1(2)-1),
c     >                     dbl_mb(EIG(1)+i-1+n1(2)-1)*EV,
c     >                     dbl_mb(occ2(1)+i-1+n1(2)-1)
c           end do
c         else
c           do i=1,NN
c             write(6,1510) dbl_mb(EIG(1)+i-1),dbl_mb(EIG(1)+i-1)*EV
c           end do
c           do i=1,ne(2)
c             write(6,1510) dbl_mb(EIG(1)+i-1+NN),
c     >                     dbl_mb(EIG(1)+i-1+NN)*EV,
c     >                     dbl_mb(EIG(1)+i-1+n1(2)-1),
c     >                     dbl_mb(EIG(1)+i-1+n1(2)-1)*EV
c           end do
c         end if
c
c         end if
c
c*        *** Extra energy output added for QA test ****
c         write(6,1600) E(2)
c      end if
c
c*                |***************************|
c******************         Prologue          **********************
c*                |***************************|
c
c*     **** calculate spin contamination ****
c      call Calculate_psi_spin2(ispin,ne,npack1,dcpl_mb(psi2(1)),
c     >                         fractional,dbl_mb(occ2(1)),
c     >                         dipole)
c
c*     **** calculate the Dipole ***
c      call Calculate_Dipole(ispin,ne,n2ft3d,dbl_mb(dn(1)),dipole)
c
c      
c*     ***** write wavefunctions and v_wavefunctions ****
c      call psi_write(ispin,ne,dcpl_mb(psi2(1)),
c     >               smearoccupation,dbl_mb(occ2(1)))
c      call v_psi_write(ispin,ne,dcpl_mb(psi0(1)))
c
c*     **** write geometry to rtdb ****
c      call ion_write(rtdb)
c
c*     **** deallocate heap memory ****
c      if (control_version().eq.3) call ewald_end()
c      call strfac_end()
c      if (controL_version().eq.3) call coulomb_end()
c      if (controL_version().eq.4) call coulomb2_end()
c      call ke_end()
c      call mask_end()
c      call Pack_end()
c      call G_end()
c      call ion_end()
c      call psp_end()
c      call ion_end_FixIon()
c      call pspw_end_SIC()
c      call pspw_end_HFX()
c      call pspw_qmmm_end()
c      call dplot_iteration_end()
cc      call frac_occ_end()
c      if (control_Nose()) call Nose_end()
c
c      value = MA_free_heap(psir(2))
c      value = MA_free_heap(Hpsi(2))
c      value = MA_free_heap(dn(2))
c      value = MA_free_heap(eig(2))
c      value = Dneall_m_free(hml)
c      value = Dneall_m_free(lmd)
c      value = MA_free_heap(psi0(2))
c      value = MA_free_heap(psi1(2))
c      value = MA_free_heap(psi2(2))
c      if (fractional) then
c      value = MA_free_heap(occ0(2))
c      value = MA_free_heap(occ1(2))
c      value = MA_free_heap(occ2(2))
c      end if
c      call D3dB_pfft_end()
c      call D3dB_end(1)
c      if (control_version().eq.4) call D3dB_end(2)
c      call Dne_end()
c      call psi_data_end()
c
c*     **** do anaylysis on MOTION files ****
c      call cpmd_properties(rtdb)
c
c
c*                |***************************|
c****************** report consumed cputime   **********************
c*                |***************************|
c      if (taskid.eq.MASTER) then
c         CALL current_second(cpu4)
c
c         T1=CPU2-CPU1
c         T2=CPU3-CPU2
c         T3=CPU4-CPU3
c         T4=CPU4-CPU1
c         AV=T2/dble(icount*it_in)
c         write(6,*)
c         write(6,*) '-----------------'
c         write(6,*) 'cputime in seconds'
c         write(6,*) 'prologue    : ',T1
c         write(6,*) 'main loop   : ',T2
c         write(6,*) 'epilogue    : ',T3
c         write(6,*) 'total       : ',T4
c         write(6,*) 'cputime/step: ',AV
c         write(6,*)
c
c         call nwpw_timing_print_final(.true.,(icount*it_in))
c         CALL nwpw_message(4)
c      end if 
c
c
      call Parallel3d_Finalize()
      call Parallel_Finalize()
      band_cpmd = value
      return
cc
ccc
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
 1112 FORMAT( ' parallel mapping         : slab')
 1113 FORMAT( ' parallel mapping         : hilbert')
 1114 FORMAT( ' parallel mapping         : balanced')
 1115 FORMAT(/' options:')
 1116 FORMAT( ' parallel mapping         : not balanced')
 1117 FORMAT( ' processor grid           :',I4,' x',I4)
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1132 FORMAT(5X,' using fractional occupation')
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
 1170 FORMAT(7(5X,A4,':',I3))
 1180 FORMAT(/' initial position of ions:')
 1181 FORMAT(/' initial velocity of ions:')
 1190 FORMAT(5X, I4, A5  ,' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1191 FORMAT(5X, I4, A5  ,' (',3F11.5,
     >       ' ) - atomic mass= ',F6.3,' - fixed')
 1192 FORMAT(5X, I4, A5  ,' (',3F11.5,' )')
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
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
  
 1260 FORMAT(5X,' ewald summation: cut radius=',F8.2,'  and',I3)
 1261 FORMAT(5X,'                   madelung=',f11.8)
 1270 FORMAT(/' technical parameters:')
 1271 FORMAT(5x, ' translation contrained')
 1272 FORMAT(5x, ' rotation constrained')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'ficticious mass=',F10.1)
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
 1330 FORMAT(/' comparison between hamiltonian and lambda matrix')
 1331 FORMAT(/' Elements of Hamiltonian matrix (up/restricted)')
 1332 FORMAT(/' Elements of Hamiltonian matrix (down)')
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
 1455 FORMAT( ' SIC-hartree energy  :',E19.10,' (',E15.5,'/electron)')
 1456 FORMAT( ' SIC-exc-corr energy :',E19.10,' (',E15.5,'/electron)')
 1457 FORMAT( ' HF exchange energy  :',E19.10,' (',E15.5,'/electron)')
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
 1499 FORMAT( ' K.S. SIC-hartree energy  :',E19.10,
     >        ' (',E15.5,'/electron)')
 1501 FORMAT( ' K.S. SIC-exc-corr energy :',E19.10,
     >        ' (',E15.5,'/electron)')
 
 1500 FORMAT(/' orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1511 FORMAT(2(E18.7,' (',F8.3,'eV) occ=',F5.3))
 1600 FORMAT(/' Total PSPW energy   :',E19.10)

 1700 FORMAT(/' QM/MM-pol-vib/CAV Energies')
 1701 FORMAT( ' --------------------------')
 1702 FORMAT( ' LJ energy              :',E19.10)
 1703 FORMAT( ' Residual Coulomb energy:',E19.10)
 1704 FORMAT( ' MM Vibration energy    :',E19.10)
 1705 FORMAT( ' MM Vibration energy    :',E19.10)
 1706 FORMAT( ' (QM+MM)/Cavity energy  :',E19.10)

 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(6,9010) ierr
      call Parallel3d_Finalize()
      call Parallel_Finalize()

      band_cpmd = value
      return
      END

      subroutine print_elapsed_time(autime)
      implicit none
      real*8 autime

      real*8 sectime

      sectime = autime*2.41889d-17

      if (sectime.lt.1.0d-12) then
         write(6,1800) (sectime/1.0d-15)," fs"
      else if (sectime.lt.1.0d-9) then
         write(6,1800) (sectime/1.0d-12)," ps"
      else 
         write(6,1800) (sectime/1.0d-9 )," ns"
      end if

      return
 1800 format(//' Elapsed time of simulation was',F8.3,A)
      end


