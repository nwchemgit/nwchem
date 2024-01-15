
***********************************************************************
*                      cpsdv5-mpi  (MPI code)                         *
*                                                                     *
*     This is a developing cpsdv5 parallel code for the SP2           *
*       + mpi message passing library used                            *
*       + ngp is used instead of nfft in this proceudure              *
*       + error checking is based on aimd.h parameters                *
*         then control file                                           *
*       + my own slap-decomposed parallel 3d-FFT(real->complex) used  *
*                                                                     *
*                                                                     *
***********************************************************************

      logical function cpsdv5(rtdb)
*
* $Id$
*
      implicit none
      integer rtdb

#include "global.fh"
#include "bafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "btdb.fh"
#include "errquit.fh"
cccc#include "frac_occ.fh"
      
*     **** parallel variables ****
      logical  omaster
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
      logical first_iteration
      real*8  icharge
      integer ispin
      integer ne(2),n1(2),n2(2),nemax,neq(2),nemaxq
      real*8  en(2),en1(2),en2(2)
      real*8  dipole(3)

*     complex*16 psi1(nfft3d,nemax)
*     complex*16 psi2(nfft3d,nemax)
*     real*8     dn(n2ft3d,2)
*     complex*16 Hpsi(nfft3d,nemax)
*     complex*16 psir(nfft3d,nemax)
      integer psi1(2),psi2(2)
      integer occ1(2),occ2(2)
      integer dn(2)
      integer Hpsi(2),psir(2)
    

*     ***** energy variables ****
      real*8  E(60),Egas,dt,emotion_time_shift

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
      logical oprint,lprint,hprint,mparallelized,found,found_bak
      integer ms,mapping,mapping1d
      real*8  deltae,deltac,deltar
      real*8  gx,gy,gz,cx,cy,cz,sum1,sum2
      real*8  EV,pi
      integer i,j,k,ia,n,nn
      integer ii,jj,indx,indx1
      integer icount,it_in,it_out
      real*8 w,sumall,virial,dv
      integer nfft3
      parameter (nfft3=32)
      character*255 full_filename,full_bak
      integer tmp1(2)

      logical value,psi_nogrid,field_exist
      integer hversion,hnfft(3),hispin,hne(2)
      real*8 hunita(3,3)
      integer ind
      character*50 filename
      character*50 control_input_psi
      external     control_input_psi
      logical  wvfnc_expander
      external wvfnc_expander


*     **** external functions ****
      real*8      psp_zv,psp_rc,ewald_rcut,ion_amass,nwpw_cosmo_qc0
      real*8      ewald_mandelung,lattice_omega_small
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg,lattice_unitg_small,lattice_unita_small
      integer     ewald_ncut,ewald_nshl3d
      integer     psp_lmmax,psp_lmax,psp_locp
      character*4 ion_aname,ion_atom
      external    psp_zv,psp_rc,ewald_rcut,ion_amass,nwpw_cosmo_qc0
      external    ewald_mandelung,lattice_omega_small
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg,lattice_unitg_small,lattice_unita_small
      external    ewald_ncut,ewald_nshl3d
      external    psp_lmmax,psp_lmax,psp_locp
      external    ion_aname,ion_atom

      character*50 control_cell_name
      external     control_cell_name

      real*8   control_tole,control_tolc,control_tolr,ion_rion
      external control_tole,control_tolc,control_tolr,ion_rion
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
      integer  control_ngrid,pack_nwave
      integer  ion_nion,ion_natm,ion_katm,ion_nkatm,ion_nkatm0
      integer  ion_nion_qm,ion_nion_mm,ion_nion_q
      external pack_nwave_all
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave
      external ion_nion,ion_natm,ion_katm,ion_nkatm,ion_nkatm0
      external ion_nion_qm,ion_nion_mm,ion_nion_q

      character*12 control_boundry
      external     control_boundry

      logical  pspw_HFX_localize2,pspw_charge_found,pspw_Efield_found
      external pspw_HFX_localize2,pspw_charge_found,pspw_Efield_found
      logical  pspw_SIC,pspw_SIC_relaxed,pspw_qmmm_found
      logical  pspw_HFX,pspw_HFX_relaxed,psp_pawexist
      logical  psp_semicore,control_Mulliken
      real*8   psp_rcore,psp_ncore,psp_rlocal
      external pspw_SIC,pspw_SIC_relaxed,pspw_qmmm_found
      external pspw_HFX,pspw_HFX_relaxed,psp_pawexist
      external psp_semicore,control_Mulliken
      external psp_rcore,psp_ncore,psp_rlocal
      logical  control_check_charge_multiplicity
      external control_check_charge_multiplicity
      real*8   nwpw_timing
      external nwpw_timing
      integer  control_np_orbital,control_mapping,control_mapping1d
      external control_np_orbital,control_mapping,control_mapping1d

      logical  control_translation,control_rotation,control_balance
      external control_translation,control_rotation,control_balance

      logical  Dneall_m_allocate,Dneall_m_free,control_parallel_io
      external Dneall_m_allocate,Dneall_m_free,control_parallel_io

      real*8   Dneall_m_value
      external Dneall_m_value
      character*9 ion_amm
      external    ion_amm
      logical  psp_U_psputerm,meta_found,ion_disp_on
      external psp_U_psputerm,meta_found,ion_disp_on
      integer  ion_nconstraints,ion_ndof
      external ion_nconstraints,ion_ndof
      logical  control_print,ion_makehmass2,control_has_ngrid_small
      external control_print,ion_makehmass2,control_has_ngrid_small
      integer  control_ngrid_small
      external control_ngrid_small
c#include "Parallel.fh"
      integer  Parallel_threadid,Parallel_nthreads,Parallel_maxthreads
      external Parallel_threadid,Parallel_nthreads,Parallel_maxthreads
      integer tid

      logical  control_mparallelized
      external control_mparallelized
      logical  nwpw_cosmo_on,nwpw_born_on
      external nwpw_cosmo_on,nwpw_born_on
      real*8   nwpw_cosmo_screen,nwpw_born_screen
      external nwpw_cosmo_screen,nwpw_born_screen
      logical  pspw_V_APC_on
      external pspw_V_APC_on
      real*8   control_gas_energy
      external control_gas_energy

      logical  control_fractional
      external control_fractional
      integer  control_fractional_smeartype
      external control_fractional_smeartype
      real*8   control_fractional_temperature
      external control_fractional_temperature
      real*8   control_fractional_kT
      external control_fractional_kT
      real*8   control_fractional_alpha
      external control_fractional_alpha
      logical  control_gram_schmidt
      external control_gram_schmidt

      integer  pspw_Efield_type
      external pspw_Efield_type



*                            |************|
*****************************|  PROLOGUE  |****************************
*                            |************|

      value = .true.
      pi = 4.0d0*datan(1.0d0)
      dt = 1.0d0

      call nwpw_timing_init()
      call dcopy(60,0.0d0,0,E,1)


*     **** get parallel variables ****
      call Parallel_Init()
      call Parallel_np(np)
      call Parallel_taskid(taskid)

      value = MA_set_hard_fail(.true.)

      tid = Parallel_threadid()
      omaster =  (taskid.eq.MASTER).and.(Parallel_threadid().eq.MASTER)

      if (omaster) call current_second(cpu1)
      if (omaster) call current_second(cpu1)

*     ***** print out header ****
      if (omaster) then
         write(luout,1000)
         write(luout,1010)
         write(luout,1020)
         write(luout,1010)
         write(luout,1030)
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
      
      value = control_read(1,rtdb)


      mparallelized = control_mparallelized()
      call Parallel2d_Init(control_np_orbital())
      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)

      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0
      mapping = control_mapping()

      oprint = omaster.and.control_print(print_medium)
      lprint = omaster.and.control_print(print_low)
      hprint = omaster.and.control_print(print_high)


*     **** initialize psi_data ****
      call psi_data_init(100)

*     **** initialize D3dB data structure ****
      call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d

*     ***** Initialize double D3dB data structure ****
      if (control_version().eq.4) 
     >   call D3dB_Init(2,2*ngrid(1),2*ngrid(2),2*ngrid(3),mapping)

      if (control_has_ngrid_small()) then
         ngrid_small(1) = control_ngrid_small(1)
         ngrid_small(2) = control_ngrid_small(2)
         ngrid_small(3) = control_ngrid_small(3)
         call D3dB_Init(3,ngrid_small(1),ngrid_small(2),ngrid_small(3),
     >                  mapping)
      end if


*     **** initialize lattice and packing data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_Init()
      call Pack_npack(1,npack1)      
      call D3dB_pfft_init()
      !call ga_sync()

*     **** read ions ****
      value = ion_init(rtdb)
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      first_iteration = .true.

*     **** allocate psp data structure and read in psedupotentials into it ****
      call psp_init()
      call psp_readall()
      if (psp_semicore(0)) call semicore_check()

*     **** initialize G,mask,ke,and coulomb data structures ****
      call ke_init()
      if (control_version().eq.3) call coulomb_init()
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()
      call phafac()
      if (control_version().eq.3) call ewald_init()

*     **** generate initial wavefunction if it does not exist ****
      if (.not.control_check_charge_multiplicity()) then
        call psi_new()
      end if

*     ***** allocate psi2,and psi1 wavefunctions ****
      call psi_get_ne_occupation(ispin,ne,smearoccupation)
      if (smearoccupation.gt.0) then
         fractional = .true.
      else
         fractional = .false.
      end if
      smearcorrection = 0.0d0
      smeartype = control_fractional_smeartype()
      smearkT   = control_fractional_kT()
      
      mapping1d = control_mapping1d()
      call Dne_init(ispin,ne,mapping1d)
      call Dneall_neq(neq)
      nemaxq = neq(1)+neq(2)
      
      value = BA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'psi1',psi1(2),psi1(1))
      if (fractional) then
         value = value.and.
     >        BA_alloc_get(mt_dbl,(ne(1)+ne(2)),'occ1',occ1(2),occ1(1))
         value = value.and.
     >        BA_alloc_get(mt_dbl,(ne(1)+ne(2)),'occ2',occ2(2),occ2(1))
      end if
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)


*     *****  read initial wavefunctions into psi2  ****
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
        if (oprint) then

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

*     *****  read psi2 wavefunctions ****
      call psi_read(ispin,ne,dcpl_mb(psi2(1)),
     >              smearoccupation,dbl_mb(occ2(1)))
      if (smearoccupation.gt.0) 
     >   call Parallel_shared_vector_copy(.true.,ne(1)+ne(2),
     >                                    dbl_mb(occ2(1)),
     >                                    dbl_mb(occ1(1)))


      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)
      nemax = ne(1)+ne(2)

*     **** Ortho Check ****
      if (psp_pawexist()) then
      call psp_overlap_S(ispin,neq,dcpl_mb(psi2(1)),dcpl_mb(psi1(1)))
      call Grsm_gg_trace(npack1,(neq(1)+neq(2)),
     >                        dcpl_mb(psi1(1)),
     >                        dcpl_mb(psi2(1)),
     >                        sum2)
      else
      call Grsm_gg_trace(npack1,(neq(1)+neq(2)),
     >                        dcpl_mb(psi2(1)),
     >                        dcpl_mb(psi2(1)),
     >                        sum2)
      end if
      call D1dB_SumAll(sum2)

      sum1 = dble(ne(1) + ne(2))
      if (dabs(sum2-sum1).gt.1.0d-10) then
         if (oprint)
     >      write(*,*) "Warning: Gram-Schmidt being performed on psi2"
         if (psp_pawexist()) then
            call Dneall_f_Sortho(0,dcpl_mb(psi2(1)),
     >                             dcpl_mb(psi1(1)),npack1)
         else
            call Dneall_f_GramSchmidt(0,dcpl_mb(psi2(1)),npack1)
         end if
c         call Dneall_f_ortho(0,dcpl_mb(psi2(1)),npack1)
c         call Grsm_g_MakeOrtho(npack1,ne(1),dcpl_mb(psi2(1)))
c         if (ispin.gt.1) then
c           call Grsm_g_MakeOrtho(npack1,ne(2),
c     >                           dcpl_mb(psi2(1) + ne(1)*npack1))
c         end if
      end if

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

      if (.not. value) 
     >   call errquit('cpsdv5:out of heap memory',0,MA_ERR)



*     **** initialize two-electron Gaussian integrals ****
*     **** initialize paw ncmp*Vloc ****
      if (psp_pawexist()) then
         call nwpw_gintegrals_init()
         call nwpw_gintegrals_set(control_move())
         call psp_dE_ncmp_vloc_Qlm(ispin,.false.,hunita)
      end if


*     **** initialize SIC and HFX  ****
      call pspw_init_SIC(rtdb,ne)
      call pspw_init_HFX(rtdb,ispin,ne)

*     **** initialize DFT+U ****
      call psp_U_init()

*     **** initialize META GGA ****
      call nwpw_meta_gga_init(control_gga())

*     **** initialize vdw ****
      call vdw_DF_init()

*     **** initialize rho_symmetry  ****
      call rho_symmetrizer_init()

*     **** initialize metadynamics ****
      call meta_initialize(rtdb)

*     **** initialize QM/MM ****
      call pspw_init_APC(rtdb)
      call pspw_qmmm_init(rtdb)
      call pspw_charge_init(rtdb)
      call pspw_Efield_init(rtdb,ispin,ne)
      !call pspw_cdft_init(rtdb)
      field_exist = pspw_charge_found().or.pspw_Efield_found()

*     **** initialize frac_occ data structure ****
c      call frac_occ_init(rtdb,ispin,ne)

*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)



*                |**************************|
******************   summary of input data  **********************
*                |**************************|

      if (oprint) then
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
           write(luout,1122)
         end if
         write(luout,1123) Parallel_maxthreads()


         write(luout,1115)
         IF(control_move()) THEN
           write(luout,1120) 'yes'
         ELSE
           write(luout,1120) 'no'
         ENDIF
         write(luout,1121) control_boundry(),control_version()
         if (ispin.eq.1) write(luout,1130) 'restricted'
         if (ispin.eq.2) write(luout,1130) 'unrestricted'

         call v_bwexc_print(luout,control_gga())

         if (fractional) write(luout,1132)
         call pspw_print_SIC(luout)
         call pspw_print_HFX(luout)
         call nwpw_meta_gga_print(luout)
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
              icharge = icharge - (3-ispin)*dbl_mb(occ1(1)+i-1)
              en(ms) = en(ms) + dbl_mb(occ1(1)+i-1)
            end do
            end do
         end if
            
         do ia=1,ion_nkatm0()
           icharge = icharge + ion_natm(ia)*psp_zv(ia)
         end do
         write(luout,1159) icharge

         write(luout,1160)
         write(luout,1170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
         write(luout,1180)
         do I=1,ion_nion()+ion_nion_q()
           if (I.gt.ion_nion()) then
           write(luout,1195) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                       nwpw_cosmo_qc0(I-ion_nion()),ion_amm(i)
           else if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
           else if (ion_q_xyzFixIon(I)) then
           write(luout,1194) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_q_xyzFixIon_label(I)
           else
           write(luout,1190) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz
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

         call rho_symmetrizer_print(luout)

         write(luout,1270)
         if (.not.control_translation()) write(luout,1271)
         if (.not.control_rotation())    write(luout,1272)
         write(luout,1280) control_time_step(),control_fake_mass()
         write(luout,1290) control_tole(),control_tolc(),control_tolr()
         write(luout,1281) control_it_in()*control_it_out(),
     >                 control_it_in(),control_it_out()
         if (control_gram_schmidt()) write(luout,1282)

        if (control_fractional()) then
           write(6,1297)
           if (control_fractional_smeartype().eq.-1)
     >       write(6,1298) "fixed occupation"
           if (control_fractional_smeartype().eq.0)
     >       write(6,1298) "step function"
           if (control_fractional_smeartype().eq.1)
     >       write(6,1298) "Fermi-Dirac"
           if (control_fractional_smeartype().eq.2)
     >       write(6,1298) "Gaussian"
           if (control_fractional_smeartype().eq.4)
     >       write(6,1298) "Marzari-Vanderbilt"
           if (control_fractional_smeartype().ge.0)
     >       write(6,1299) control_fractional_kT(),
     >                     control_fractional_temperature(),
     >                     control_fractional_alpha()
         end if

         write(luout,1300)
         write(luout,1305)
         call util_flush(luout)
      end if

*                |***************************|
******************     start iterations      **********************
*                |***************************|

*     **** open xyz and MOTION and fei file ****
      call xyz_init()          ! unit=18
      call MOTION_init(rtdb)   ! unit=19
      call fei_init(rtdb)

*     **** open EMOTION file ****
      if (.not.btdb_cget(rtdb,'cpsd:emotion_filename',1,filename))
     >  call util_file_prefix('emotion',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)
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
         inquire(file=full_filename,exist=found)
         if (found) then

*          **** make a new backup file ***
           call util_file_copy(full_filename,full_bak)

           open(unit=31,file=full_filename,form='formatted',
     >          status='old')
           do while (found)
           read(31,*,end=100) emotion_time_shift,w,sumall,gx,gy,gz
           end do
  100      continue
#if defined(FUJITSU) || defined(PSCALE) || defined(__crayx1) || defined(GCC46)
#warning backspace           
           backspace 31
#endif
         else
           open(unit=31,file=full_filename,form='formatted',
     >          status='new')
         end if
      end if


      if (oprint) call current_second(cpu2)
      if (oprint) call nwpw_message(14)
     
      it_in  = control_it_in()
      it_out = control_it_out()
      icount = 0
   1  continue
         icount = icount + 1

         call inner_loop(ispin,neq,
     >             npack1,nfft3d,nemaxq,
     >             dcpl_mb(psi1(1)),dcpl_mb(psi2(1)),dbl_mb(dn(1)),
     >             it_in,E,deltae,deltac,deltar,
     >             dbl_mb(hml(1)),dbl_mb(lmd(1)),dbl_mb(lmd1(1)),
     >             first_iteration,
     >             dcpl_mb(psir(1)),dcpl_mb(Hpsi(1)),
     >             fractional,dbl_mb(occ1(1)),dbl_mb(occ2(1)))
         if (fractional) E(1) = E(1) + smearcorrection

         if (oprint) then 
           write(luout,1310) icount*it_in,E(1),deltae,deltac,deltar
           call util_flush(luout)
         end if
         if ((deltae.gt.0.0d0).and.(icount.gt.1)) then
            if ((icount.ge.it_out).or.(.not.fractional)) then
               if (oprint) 
     >          write(luout,*) 
     >          ' *** Energy going up.  iteration terminated.'
               go to 2
            end if
         end if
         deltae = dabs(deltae)
         if ((deltae.lt.control_tole()).and.
     >       (deltac.lt.control_tolc()).and.
     >       (deltar.lt.control_tolr())) then
            if (oprint) 
     >       write(luout,*) 
     >       ' *** tolerance ok.     iteration terminated.'
            go to 2
         end if


         !**** emotion write ****
         if (taskid.eq.MASTER) then
            write(31,1311) icount*it_in*dt + emotion_time_shift,
     >                     e(1),deltae,deltac,deltar,e(2),e(3),e(4),
     >                     e(5),e(6),e(7),e(8),e(9),e(10)
            call util_flush(31)
         end if

         !**** xyz write ****
         call xyz_novelocity_write()

*        ***** define fractional occupation ****
         if ((icount.lt.it_out).and.fractional.and.
     >       (smeartype.ge.0)) then

            call Dneall_m_diagonalize(0,dbl_mb(hml(1)),
     >                                  dbl_mb(eig(1)),.true.)
c            call Dneall_fmf_Multiply(0,dcpl_mb(psi2(1)),npack1,
c     >                                 dbl_mb(hml(1)),  1.0d0,
c     >                                 dcpl_mb(psi1(1)),0.0d0)

            call psi_0define_occupation(-1.0d0,.false.,
     >                               ispin,ne,
     >                               dbl_mb(eig(1)),dbl_mb(hml(1)),
     >                               dbl_mb(occ2(1)),
     >                               smeartype,smearkT,
     >                               smearfermi,smearcorrection)
c            write(*,*) "definefrac smearfermi=",smearfermi
c            do i=1,ne(1)
c               write(*,*) "i,occ=",i,dbl_mb(eig(1)+i-1),
c     >                               dbl_mb(occ1(1)+i-1)
c            end do
         end if



      if (icount.lt.it_out) go to 1
      if (oprint) 
     > write(luout,*) 
     > '*** arrived at the Maximum iteration.   terminated.'

*::::::::::::::::::::  end of iteration loop  :::::::::::::::::::::::::

   2  continue
      if (oprint) CALL nwpw_message(3)
      if (oprint) call current_second(cpu3)

      if (taskid.eq.MASTER) then
        close(unit=31)
*        **** remove EMOTION backup file ***
         call util_file_name_noprefix('EMOTION99-bak',.false.,
     >                                .false.,
     >                                full_bak)
         call util_file_unlink(full_bak)
      end if

*     *** close xyz and fei io ****
      call xyz_end()
      call MOTION_end()
      call fei_end()

*         |****************************************|
*********** produce CHECK file and diagonalize hml *****************
*         |****************************************|

*     **** produce CHECK FILE ****
      if (oprint) then
         call util_file_name('CHECK',.true.,
     >                               .false.,
     >                        full_filename)
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
     >        'cgsdv5: out of stack memory',0,MA_ERR)

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
     >        'cgsdv5: popping stack memory',0,MA_ERR)
      else
         en2(1) = 0.0d0
         en2(2) = 0.0d0
      end if
      en(1) = en1(1)+en2(1)
      en(2) = en1(2)+en2(2)
      if (omaster) then
         write(17,1320) (en(ms),ms=1,ispin)
         if (psp_pawexist()) then
            write(17,1322) (en1(ms),ms=1,ispin)
            write(17,1323) (en2(ms),ms=1,ispin)
         end if
      end if



*     **** comparison between hamiltonian an lambda matrix ****
      if (omaster) write(17,1330)
      do ms=1,ispin
         do j=1,ne(ms)
         do i=1,ne(ms)
            w      = Dneall_m_value(0,ms,i,j,dbl_mb(hml(1)))
            virial = Dneall_m_value(0,ms,i,j,dbl_mb(lmd(1)))

            if (omaster)
     >      write(17,1340) ms,i,j,w,virial,w-virial
         end do
         end do
      end do


*     **** check orthonormality ****
      if (omaster) then
         write(17,1350)
      end if
      call Dneall_ffm_Multiply(0,dcpl_mb(psi1(1)),
     >                           dcpl_mb(psi1(1)),npack1,
     >                           dbl_mb(lmd(1)))
      do ms=1,ispin
         do j=1,ne(ms)
         do i=j,ne(ms)
            w  = Dneall_m_value(0,ms,i,j,dbl_mb(lmd(1)))
            if (omaster) write(17,1360) ms,i,j,w
         end do
         end do
      end do

*     **** close check file ****
      if (omaster) then
         close(17)
      end if

*     ***** do not diagonalize the hamiltonian matrix *****
      if (pspw_SIC().or.mparallelized) then
        call dcopy(2*npack1*nemaxq,
     >             dcpl_mb(psi1(1)),1,
     >             dcpl_mb(psi2(1)),1)

*     ***** diagonalize the hamiltonian matrix but don't rotate ****
      else if (fractional) then
         n = ne(1)
c         call Dneall_m_HmltimesSA(0,dbl_mb(hml(1)),dbl_mb(fweight(1)))
         call Dneall_m_diagonalize(0,dbl_mb(hml(1)),
     >                              dbl_mb(eig(1)),.false.)
c         do ms=1,ispin
c         do ii=1,ne(ms)
c            dbl_mb(eig(1)+(ii-1)+(ms-1)*n)
c     >      =dbl_mb(eig(1)+(ii-1)+(ms-1)*n)
c     >      /dbl_mb(fweight(1)+(ii-1)+(ms-1)*n)
c         end do
c         end do
         !*** reverse the occupation ***
         do ms=1,ispin
            do i=1,ne(ms)
                 dbl_mb(occ1(1)+(ms-1)*ne(1) + i-1) = 
     >           dbl_mb(occ2(1)+(ms-1)*ne(1) + ne(ms)-i)
            end do
         end do
         !call dcopy(ne(1)+ne(2),dbl_mb(occ2(1)),1,dbl_mb(occ1(1)),1)
         call dcopy(2*npack1*nemaxq,
     >             dcpl_mb(psi1(1)),1,
     >             dcpl_mb(psi2(1)),1)

*     ***** diagonalize and rotate the hamiltonian matrix ****
      else
        call Dneall_m_diagonalize(0,dbl_mb(hml(1)),
     >                              dbl_mb(eig(1)),.false.)

*       **** do not rotate if wannier localization ***
        if (.not.pspw_HFX_localize2())
     >      call Dneall_fmf_Multiply(0,dcpl_mb(psi1(1)),npack1,
     >                                 dbl_mb(hml(1)),  1.0d0,
     >                                 dcpl_mb(psi2(1)),0.0d0)
      end if

*                |***************************|
****************** report summary of results **********************
*                |***************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

      if (oprint) then
         write(luout,1300)
         write(luout,1410)
         write(luout,1420)
         do I=1,ion_nion()+ion_nion_q()
           if (I.gt.ion_nion()) then
           write(luout,1195) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   nwpw_cosmo_qc0(I-ion_nion()),ion_amm(i)
           else if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
           else if (ion_q_xyzFixIon(I)) then
           write(6,1194) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_q_xyzFixIon_label(I)
           else
           write(luout,1190) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz
         write(luout,1211) ion_nconstraints(),ion_ndof()

         call pspw_charge_Print(luout)

         write(luout,*)
         write(luout,1320) en(1),en(ispin),' (real space)'
         if (psp_pawexist()) then
            write(luout,1322) en1(1),en1(ispin),' (real space)'
            write(luout,1323) en2(1),en2(ispin),' (real space)'
         end if

*       **** write APC potential and charges ***
        if (pspw_V_APC_on()) call pspw_shortprint_APC(luout)


         write(luout,1430) E(1),E(1)/ion_nion()
         if (pspw_qmmm_found()) then
            write(luout,1431) 
            write(luout,1432) 
            write(luout,1433) (E(1)-E(11)),
     >                    (E(1)-E(11))/ion_nion()
         end if

         if (field_exist) then
            write(luout,1431)
            write(luout,1432)
            if (pspw_Efield_type().eq.0) then
               write(luout,1433) (E(1)-E(49)-E(51)),
     >         (E(1)-E(49)-E(51))/ion_nion()
            else
               write(luout,1433) (E(1)-E(49)-E(50)-E(51)),
     >         (E(1)-E(49)-E(50)-E(51))/ion_nion()
            end if
         end if

         !*** print out total paw energy including core energies ****
         if (psp_pawexist()) 
     >      write(luout,1434) (E(1)+E(36)+E(45)),
     >                        (E(1)+E(36)+E(45))/ion_nion()

         write(luout,1440) E(2),E(2)/n2(ispin)
         write(luout,1450) E(3),E(3)/n2(ispin)
         write(luout,1460) E(4),E(4)/n2(ispin)
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
         if (pspw_V_APC_on()) then
           write(luout,1505) E(52),E(52)/ion_nion()
         end if


         write(luout,1470) E(5),E(5)/ion_nion()
         write(luout,1480) E(6),E(6)/n2(ispin)
         write(luout,1490) E(7),E(7)/n2(ispin)
         write(luout,1495) E(8),E(8)/n2(ispin)
         write(luout,1496) E(9),E(9)/n2(ispin)
         write(luout,1497) E(10),E(10)/n2(ispin)
         if (pspw_SIC().and.pspw_SIC_relaxed())  then
           write(luout,1499) E(18),E(18)/n2(ispin)
           write(luout,1501) E(19),E(19)/n2(ispin)
         end if
         if (pspw_HFX().and.pspw_HFX_relaxed())  then
           write(luout,1502) E(21),E(21)/n2(ispin)
         end if
         if (psp_U_psputerm()) then
           write(luout,1503) E(30),E(30)/n2(ispin)
         end if
         if (meta_found()) then
           write(luout,1504) E(32),E(32)/n2(ispin)
         end if
         if (pspw_V_APC_on()) then
           write(luout,1506) E(53),E(53)/n2(ispin)
         end if

         virial = (E(10)+E(9)+E(8)+E(7))/E(6)
         write(luout,1498) virial

c        **** paw energies ****
         if (psp_pawexist()) then
         write(luout,1680) E(35),E(35)/n2(ispin)
         write(luout,1681) E(36),E(36)/ion_nion()
         write(luout,1682) E(37),E(37)/n2(ispin)

         write(luout,1683) E(38),E(38)/n2(ispin)
         write(luout,1684) E(39),E(39)/n2(ispin)

         write(luout,1685) E(40),E(40)/n2(ispin)
         write(luout,1686) E(41),E(41)/n2(ispin)
         write(luout,1687) E(42),E(42)/n2(ispin)
         write(luout,1688) E(43),E(43)/n2(ispin)
         write(luout,1689) E(44),E(44)/ion_nion()
         write(luout,1690) E(45),E(45)/ion_nion()
         end if

      if (pspw_qmmm_found()) then
            write(luout,1700)
            write(luout,1701)
            write(luout,1702) E(11)
            write(luout,1703) E(12)
            write(luout,1704) E(13)
            write(luout,1705) E(14)
c
c            write(luout,1703) E(14)+E(15)
c            write(luout,1710) E(14)
c            write (luout,1711) E(15)
c
c            write(luout,1704) E(16)
c            write(luout,1705) E(17)
c            write(luout,1706) E(18)
        end if
        if (ion_disp_on()) then
            write(luout,1720) E(33)
        end if

        if (field_exist) then
           write(luout,1800)
           write(luout,1801)
           if (pspw_Efield_type().eq.0) then
              write(luout,1805) E(49)+E(51)
              write(luout,1806) E(49)
              write(luout,1807) E(50)
              !write(luout,1804) E(51)
              call dipole_Efield_print(luout)
           else
              write(luout,1805) E(49)+E(50)+E(51)
              write(luout,1802) E(49)
              write(luout,1803) E(50)
              write(luout,1804) E(51)
           end if
        end if



         if (.not.mparallelized) then
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
         EV=27.2116d0
         if (control_fractional()) then
            if (ispin.eq.1) then
               write(luout,1507) smearfermi(1),smearfermi(1)*EV
            else
               write(luout,1507) smearfermi(1),smearfermi(1)*EV,
     >                           smearfermi(2),smearfermi(2)*EV
            end if
         end if
         write(luout,1500)
         NN=NE(1)-NE(2)
         if (fractional) then
           do i=1,NN
             write(luout,1511) dbl_mb(EIG(1)+i-1),
     >                     dbl_mb(EIG(1)+i-1)*EV,
     >                     dbl_mb(occ1(1)+i-1)
           end do
           do i=1,ne(2)
             write(luout,1511) dbl_mb(EIG(1)+i-1+NN),
     >                     dbl_mb(EIG(1)+i-1+NN)*EV,
     >                     dbl_mb(occ1(1)+i-1+NN),
     >                     dbl_mb(EIG(1)+i-1+n1(2)-1),
     >                     dbl_mb(EIG(1)+i-1+n1(2)-1)*EV,
     >                     dbl_mb(occ1(1)+i-1+n1(2)-1)
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
      end if

*     **** write out COSMO energies ****
      if (nwpw_cosmo_on()) then
         Egas = control_gas_energy()
         if (oprint) then
            write(luout,1730) 
            write(luout,1731) nwpw_cosmo_screen()
            write(luout,1732) E(46)
            write(luout,1733) E(47)
            write(luout,1734) E(48)
            write(luout,1735) E(46)+E(47)+E(48),
     >                        (E(46)+E(47)+E(48))*27.2116d0*23.06d0
            if (dabs(Egas).gt.1.0d-6)
     >         write(luout,1736)  E(1)-Egas,
     >                           (E(1)-Egas)*27.2116d0*23.06d0
            call nwpw_cosmo_print(luout,Egas,E(1))
         end if
      end if

*     **** write out extended Born solvation energies ****
      if (nwpw_born_on()) then
         if (pspw_V_APC_on()) then
            Egas = control_gas_energy()
            if (oprint) then
               write(luout,1740)
               write(luout,1741) nwpw_born_screen()
               write(luout,1745) E(52),E(52)*27.2116d0*23.06d0
               if (dabs(Egas).gt.1.0d-6)
     >            write(luout,1746)  E(1)-Egas,
     >                           (E(1)-Egas)*27.2116d0*23.06d0
               call nwpw_born_print(luout,Egas,E(1))
            end if
         else
            call pspw_dngen_APC(ispin,ne,dbl_mb(dn(1)),.false.)
            call pspw_print_APC(luout)
         end if
      end if


      if (mparallelized) then
          do ms=1,ispin
             if ((ms.eq.1).and.(oprint)) write(luout,1333)
             if ((ms.eq.2).and.(oprint)) write(luout,1334)
             do i=1,ne(ms)
                w = Dneall_m_value(0,ms,i,i,dbl_mb(hml(1)))
                if (oprint) write(luout,1510) w,w*EV
             end do
          end do
      end if

      if (oprint) then
*        ***** extra energy output for QA test ****
         write(luout,1600) E(1)
      end if

   

*                |***************************|
******************         Prologue          **********************
*                |***************************|

*     **** calculate spin contamination ****
      call Calculate_psi_spin2(ispin,ne,npack1,dcpl_mb(psi2(1)),
     >                         fractional,dbl_mb(occ2(1)),w)

*     **** calculate the Dipole ***
      call Calculate_Dipole(ispin,ne,n2ft3d,dbl_mb(dn(1)),dipole)
      
*     **** perfom Lubin and Mulliken analysis ***
      if (control_Mulliken()) then

*       **** Lubin Water Analysis ***
        call pspw_Lubin_water_analysis(rtdb,ispin,ne,n2ft3d,
     >                                 dbl_mb(dn(1)))

*       **** Analysis ***
        call pspw_analysis(0,rtdb,ispin,ne,dcpl_mb(psi2(1)),
     >                                   dbl_mb(eig(1)))

*       **** generate APC *****
        call pspw_dngen_APC(ispin,ne,dbl_mb(dn(1)),.false.)
        call pspw_print_APC(luout)

      end if

*     ***** write psi2 wavefunctions ****
      call psi_write(ispin,ne,dcpl_mb(psi2(1)),
     >               smearoccupation,dbl_mb(occ2(1)))
 
*     **** write geometry to rtdb ****
      call pspw_charge_write(rtdb)
      call ion_write(rtdb)


*     **** deallocate heap memory ****
      if (control_version().eq.3) call ewald_end()
      call strfac_end()
      if (control_version().eq.3) call coulomb_end()
      if (control_version().eq.4) call coulomb2_end()
      call ke_end()
      call mask_end()
      call Pack_end()
      call G_end()
      call psp_U_end()
      call vdw_DF_end()
      call nwpw_meta_gga_end()
      call rho_symmetrizer_end()
      call pspw_end_SIC()
      call pspw_end_HFX()
      call pspw_end_APC()
      call pspw_qmmm_end()
      call pspw_charge_end()
      call pspw_Efield_end()
      call meta_finalize(rtdb)
      if (psp_pawexist()) call nwpw_gintegrals_end()

      call ion_end()
      call psp_end()
      call ion_end_FixIon()

c      call frac_occ_end()

      value = BA_free_heap(psir(2))
      value = BA_free_heap(Hpsi(2))
      value = BA_free_heap(dn(2))
      value = BA_free_heap(eig(2))
      value = Dneall_m_free(hml)
      value = Dneall_m_free(lmd)
      value = Dneall_m_free(lmd1)
      value = BA_free_heap(psi1(2))
      value = BA_free_heap(psi2(2))
      if (fractional) then
         value = BA_free_heap(occ1(2))
         value = BA_free_heap(occ2(2))
      end if
      call D3dB_pfft_end()
      call D3dB_end(1)
      if (control_version().eq.4) call D3dB_end(2)
      if (control_has_ngrid_small()) call D3dB_end(3)
      call Dne_end()
      call psi_data_end()

*                |***************************|
****************** report consumed cputime   **********************
*                |***************************|
      if (oprint) then
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
         CALL nwpw_MESSAGE(4)
      end if 


      call Parallel2d_Finalize()
      call Parallel_Finalize()
      cpsdv5 = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*     Car-Parrinello microcluster calculation      *')
 1030 FORMAT(10X,'*     [     steepest descent minimization   ]      *')
 1035 FORMAT(10x,'*     [ NorthWest Chemistry implementation ]       *')
 1040 FORMAT(10X,'*            version #5.00   06/01/99              *')
 1041 FORMAT(10X,'*    This code was developed by Eric J. Bylaska,   *')
 1042 FORMAT(10X,'*    and was based upon algorithms and code        *')
 1043 FORMAT(10X,'*    developed by the group of Prof. John H. Weare *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ PSPW input data ===================')
 1111 FORMAT(/' number of processors used:',I10)
 1112 FORMAT( ' parallel mapping         :      1d-slab')
 1113 FORMAT( ' parallel mapping         :   2d-hilbert')
 1114 FORMAT( ' parallel mapping         :     balanced')
 1115 FORMAT(/' options:')
 1116 FORMAT( ' parallel mapping         : not balanced')
 1117 FORMAT( ' processor grid           :',I4,' x',I4)
 1118 FORMAT( ' parallel mapping         :    2d-hcurve')
 1119 FORMAT( ' parallel io              :        on')
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1122 FORMAT( ' parallel io              :       off')
 1123 FORMAT( ' number of threads        :',I10)
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1132 FORMAT(5X,' using fractional occupation')
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
 1170 FORMAT(7(5X,A2,':',I5))
 1180 FORMAT(/' initial position of ions:')
 1190 FORMAT(5X, I4, A5, ' (',3F11.5,' ) - atomic mass= ',F7.3,' ',A)
 1191 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - fixed ',A)
 1193 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - z fixed')
 1194 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,A)
 1195 FORMAT(5X, I4, A5, ' (',3F11.5,' ) -      charge= ',E13.6,' ',A)
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1211 FORMAT(5X,'   number of constraints = ', I6,' ( DOF = ',I6,' )' )
 1219 FORMAT(/' number of electrons: spin up=',F6.2, 16x,
     >                               '  down=',F6.2,A)
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
 1282 FORMAT(5X, ' Gram Schmidt orthonormality used')
 1290 FORMAT(5X, ' tolerance=',E9.3,' (energy)',E12.3,
     &        ' (electron)',E12.3,' (ion)')
 1297 FORMAT(/' fractional smearing parameters:')
 1298 FORMAT(5X, ' smearing algorithm   = ',A)
 1299 FORMAT(5X, ' smearing parameter   = ',E9.3,' (',F7.1,' K)'/,
     >       5X, ' mixing parameter     =',F7.4)
 1300 FORMAT(//)
 1305 FORMAT(10X,
     >       '==================== iteration =========================')
 1310 FORMAT(I8,E20.10,3E15.5)
 1311 format(100e19.10)
 1320 FORMAT(' number of electrons: spin up=',F11.5,'  down=',F11.5,A)
 1321 FORMAT(' total charge of system:',F11.5,A)
 1322 FORMAT('     plane-wave part:         ',F11.5,'       ',F11.5,A)
 1323 FORMAT('      augmented part:         ',F11.5,'       ',F11.5,A)

 1330 FORMAT(/' comparison between hamiltonian and lambda matrix')
 1331 FORMAT(/' Elements of Hamiltonian matrix (up/restricted)')
 1332 FORMAT(/' Elements of Hamiltonian matrix (down)')
 1333 FORMAT(/' Diagonal Elements of Hamiltonian matrix(up/restricted)')
 1334 FORMAT(/' Diagonal Elements of Hamiltonian matrix(down)')
 1340 FORMAT(I5,2I5,' H=',E16.7,', L=',E16.7,', H-L=',E16.7)
 1341 FORMAT(I5,2I5,' H=',E16.6)
 1350 FORMAT(/' orthonormality')
 1360 FORMAT(I5,2I5,E18.7)
 1370 FORMAT(I3)
 1380 FORMAT(' ''',a,'''',I4)
 1390 FORMAT(I3)
 1400 FORMAT(I3,3E18.8/3X,3E18.8)
 1410 FORMAT(10X,'=============  summary of results  =================')
 1420 FORMAT( ' final position of ions:')
 1430 FORMAT(//' total     energy    :',E19.10,' (',E15.5,'/ion)')
 1431 FORMAT(/' QM Energies')
 1432 FORMAT( '------------')
 1433 FORMAT( ' total  QM energy    :',E19.10,' (',E15.5,'/ion)')
 1434 FORMAT(//' total paw energy    :',E19.10,' (',E15.5,'/ion)')
 1440 FORMAT( ' total orbital energy:',E19.10,' (',E15.5,'/electron)')
 1450 FORMAT( ' hartree   energy    :',E19.10,' (',E15.5,'/electron)')
 1455 FORMAT( ' SIC-hartree energy  :',E19.10,' (',E15.5,'/electron)')
 1456 FORMAT( ' SIC-exc-corr energy :',E19.10,' (',E15.5,'/electron)')
 1457 FORMAT( ' HF exchange energy  :',E19.10,' (',E15.5,'/electron)')
 1458 FORMAT( ' DFT+U     energy    :',E19.10,' (',E15.5,'/electron)')
 1459 FORMAT( ' Metadynamics energy :',E19.10,' (',E15.5,'/ion)')
 1460 FORMAT( ' exc-corr  energy    :',E19.10,' (',E15.5,'/electron)')
 1470 FORMAT( ' ion-ion   energy    :',E19.10,' (',E15.5,'/ion)')
 1480 FORMAT(/' K.S. kinetic energy :',E19.10,' (',E15.5,'/electron)')
 1490 FORMAT( ' K.S. V_l  energy    :',E19.10,' (',E15.5,'/electron)')
 1495 FORMAT( ' K.S. V_nl energy    :',E19.10,' (',E15.5,'/electron)')
 1496 FORMAT( ' K.S. V_Hart energy  :',E19.10,' (',E15.5,'/electron)')
 1497 FORMAT( ' K.S. V_xc energy    :',E19.10,' (',E15.5,'/electron)')
 1498 FORMAT( ' Virial Coefficient  :',E19.10)
 1499 FORMAT( ' K.S. SIC-hartree energy  :',E19.10,
     >        ' (',E15.5,'/electron)')
 1500 FORMAT(/' orbital energies:')
 1501 FORMAT( ' K.S. SIC-exc-corr energy :',E19.10,
     >        ' (',E15.5,'/electron)')
 1502 FORMAT( ' K.S. HFX energy     :',E19.10,
     >        ' (',E15.5,'/electron)')
 1503 FORMAT( ' K.S. DFT+U energy   :',E19.10,' (',E15.5,'/electron)')
 1504 FORMAT( ' K.S. Metadynamics energy :',E19.10,' (',E15.5,'/ion)')
 1505 FORMAT( ' APC energy          :',E19.10,' (',E15.5,'/ion)')
 1506 FORMAT( ' K.S. V_APC energy   :',E19.10,' (',E15.5,'/ion)')
 1507 FORMAT(/' Fermi energy =',2(E18.7,' (',F8.3,'eV)'))


 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1511 FORMAT(2(E18.7,' (',F8.3,'eV) occ=',F5.3))
 1600 FORMAT(/' Total PSPW energy   :',E19.10)

 1680 FORMAT(/' kinetic (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1681 FORMAT( ' kinetic       (core):',E19.10,' (',E15.5,'/ion)')
 1682 FORMAT( ' valence-core        :',E19.10,' (',E15.5,'/electron)')
 1683 FORMAT( ' ncmp*V_local        :',E19.10,' (',E15.5,'/electron)')
 1684 FORMAT( ' V_lpaw  (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1685 FORMAT( ' coulomb (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1686 FORMAT( ' coulomb  (ncmp-ncmp):',E19.10,' (',E15.5,'/electron)')
 1687 FORMAT( ' coulomb    (ncmp-pw):',E19.10,' (',E15.5,'/electron)')
 1688 FORMAT( ' exc-cor (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1689 FORMAT( ' pxc-cor (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1690 FORMAT( ' ion-core            :',E19.10,' (',E15.5,'/ion)')
 1691 FORMAT( ' V_local PAW residual:',E19.10,' (',E15.5,'/electron)')

 1700 FORMAT(/' QM/MM-pol-vib/CAV Energies')
 1701 FORMAT( ' --------------------------')
 1702 FORMAT( ' LJ energy              :',E19.10)
 1703 FORMAT( ' Residual Coulomb energy:',E19.10)
 1704 FORMAT( ' MM Vibration energy    :',E19.10)
 1705 FORMAT( ' QM/MM coupling energy  :',E19.10)
 1706 FORMAT( ' (QM+MM)/Cavity energy  :',E19.10)

 1720 FORMAT(/' Dispersion energy   :',E19.10)

 1730 FORMAT(/' COSMO energies:')
 1731 FORMAT(5x,' screen=(epsilon-1)/(epsilon+1/2):',F11.6)
 1732 FORMAT(5x,' screen*Qelc*B*q energy  :',E19.10)
 1733 FORMAT(5x,' screen*Qion*B*q energy  :',E19.10)
 1734 FORMAT(5x,' screen*0.5*q*A*q energy :',E19.10)
 1735 FORMAT(5x,' solvation energy (w/o QM polarization) :',E19.10,
     >   ' (',F8.3,' kcal/mol)')
 1736 FORMAT(5x,' solvation energy (w/  QM polarization) :',E19.10,
     >   ' (',F8.3,' kcal/mol)')

 1740 FORMAT(/' extended Born solvation energies:')
 1741 FORMAT(5x,' screen=(epsilon-1)/(epsilon):',F11.6)
 1745 FORMAT(5x,' solvation energy (w/o QM polarization) :',E19.10,
     >   ' (',F8.3,' kcal/mol)')
 1746 FORMAT(5x,' solvation energy (w/  QM polarization) :',E19.10,
     >   ' (',F8.3,' kcal/mol)')


 1800 FORMAT(/' Charge+Electric Field Energies')
 1801 FORMAT( ' ------------------------------')
 1802 FORMAT( ' - Charge+Electric Field/Electron    :',E19.10)
 1803 FORMAT( ' - Charge+Electric Field/Ion         :',E19.10)
 1804 FORMAT( ' - Charge+Electric Field/Charge Field:',E19.10)
 1805 FORMAT( ' Charge+Electric Field Energy        :',E19.10)
 1806 FORMAT( ' - Electric Field/Resta_Dipole       :',E19.10)
 1807 FORMAT( ' - K.S. Electric Field/Resta_Dipole  :',E19.10)

 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (omaster) write(luout,9010) ierr
      call Parallel2d_Finalize()
      call Parallel_Finalize()

      cpsdv5 = value
      return
      END
