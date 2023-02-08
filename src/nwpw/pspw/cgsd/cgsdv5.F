*
* $Id$
*

***********************************************************************
*                      cgsdv5-mpi  (MPI code)                         *
*                                                                     *
*     This is a developing cgsdv5 parallel code for NWCHEM            *
*       + tcgmsg message passing library used                         *
*       + my own slap-decomposed parallel 3d-FFT(real->complex) used  *
*                                                                     *
*                                                                     *
***********************************************************************

      logical function cgsdv5(rtdb,flag)
      implicit none
      integer rtdb
      integer flag

#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "inp.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"

#include "nwpw_timing.fh"
      
*     **** parallel variables ****
      integer  taskid,np,np_i,np_j
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,nfft3d,n2ft3d
      integer ngrid_small(3)
      real*8  a,b,c,alpha,beta,gamma
    

*     ***** energy variables ****
      integer ispin,ne(2),ms
      real*8  E(60),en(2)
      real*8  dipole(3),cdipole(3)
      real*8  stress(3,3),lstress(6)

*     **** gradient variables ****
      integer fion(2)

*     **** error variables ****
      logical value,newpsi
      integer ierr

*     **** local variables ****
      logical oprint,lprint,hprint,qmmm
      real*8  gx,gy,gz,cx,cy,cz
      real*8  EV,EVQ,pi,rs
      real*8  icharge
      integer i,k,ia,nion,minimizer,mapping,mapping1d


*     **** external functions ****
      real*8      psp_zv,psp_rc,ewald_rcut
      real*8      ewald_mandelung
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg,ion_amass
      integer     ewald_ncut,ewald_nshl3d
      integer     psp_lmmax,psp_lmax,psp_locp,ion_nkatm
      integer     psp_nprj,psp_psp_type
      character*4 ion_atom,ion_aname
      external    psp_zv,psp_rc,ewald_rcut
      external    ewald_mandelung
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg,ion_amass
      external    ewald_ncut,ewald_nshl3d
      external    psp_lmmax,psp_lmax,psp_locp,ion_nkatm
      external    psp_nprj,psp_psp_type
      external    ion_atom,ion_aname
    
      external      psp_comment
      character*255 psp_comment,comment

      real*8   control_tole,control_tolc,control_tolr,ion_rion
      external control_tole,control_tolc,control_tolr,ion_rion
      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      logical  control_read,control_move,ion_init
      external control_read,control_move,ion_init

      integer  pack_nwave_all
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      integer  pspw_charge_nion
      external pack_nwave_all
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      external pspw_charge_nion
      integer  ion_nion_q
      external ion_nion_q
      real*8   nwpw_cosmo_qc0
      external nwpw_cosmo_qc0
      character*9 ion_amm
      external    ion_amm

      character*12 control_boundry
      external     control_boundry
      character*50 control_cell_name
      external     control_cell_name

      logical      psp_semicore
      real*8       psp_rcore,psp_ncore,psp_rlocal
      external     psp_semicore
      external     psp_rcore,psp_ncore,psp_rlocal

      logical  psi_initialize, psi_finalize
      integer  psi_ispin,psi_ne,psi_neq,electron_count,linesearch_count
      external psi_initialize, psi_finalize
      external psi_ispin,psi_ne,psi_neq,electron_count,linesearch_count
      real*8   cgsd_energy,cgsd_noit_energy
      external cgsd_energy,cgsd_noit_energy
      logical  control_Mulliken,control_DOS,nwpw_born_on
      external control_Mulliken,control_DOS,nwpw_born_on
      real*8   ion_TotalCharge,psi_occupation
      external ion_TotalCharge,psi_occupation
      logical   control_check_charge_multiplicity
      external  control_check_charge_multiplicity
      logical   pspw_charge_found,ion_q_FixIon
      external  pspw_charge_found,ion_q_FixIon
      integer  control_minimizer,control_scf_algorithm
      external control_minimizer,control_scf_algorithm
      integer  control_ks_algorithm,control_diis_histories
      external control_ks_algorithm,control_diis_histories
      real*8   control_ks_alpha,control_kerker_g0
      external control_ks_alpha,control_kerker_g0
      logical  control_print,control_balance
      external control_print,control_balance
      integer  control_mapping,control_mapping1d,control_np_orbital
      external control_mapping,control_mapping1d,control_np_orbital
      integer  control_ks_maxit_orb,control_ks_maxit_orbs
      external control_ks_maxit_orb,control_ks_maxit_orbs
      logical  control_fractional,control_parallel_io
      external control_fractional,control_parallel_io
      integer  control_fractional_smeartype,control_excited_ne
      external control_fractional_smeartype,control_excited_ne
      real*8   control_fractional_kT,control_fractional_temperature
      external control_fractional_kT,control_fractional_temperature
      real*8   control_fractional_alpha
      external control_fractional_alpha
      logical  psp_pawexist,psp_paw_use_grid_cmp
      external psp_pawexist,psp_paw_use_grid_cmp
      integer  psp_paw_mult_l_max
      external psp_paw_mult_l_max
      logical  control_gram_schmidt
      external control_gram_schmidt

      integer  Parallel_maxthreads
      external Parallel_maxthreads


      logical pspw_bqext,control_fast_erf
      external pspw_bqext,control_fast_erf
      integer  ion_nconstraints,ion_ndof
      external ion_nconstraints,ion_ndof
      logical  control_smooth_cutoff,control_periodic_dipole
      external control_smooth_cutoff,control_periodic_dipole
      real*8   control_smooth_cutoff_values
      external control_smooth_cutoff_values
      logical  control_has_ngrid_small
      external control_has_ngrid_small
      integer  control_ngrid_small
      external control_ngrid_small
      real*8   lattice_unita_small,lattice_unitg_small
      real*8   lattice_omega_small
      external lattice_unita_small,lattice_unitg_small
      external lattice_omega_small
      integer  ion_rion_ptr
      external ion_rion_ptr
      logical  control_runsocket,control_precondition
      external control_runsocket,control_precondition
      logical  control_gas_energy_set
      external control_gas_energy_set


*                            |************|
*****************************|  PROLOGUE  |****************************
*                            |************|

      call ecce_print_module_entry('pspw')

      value = .true.
      pi = 4.0d0*datan(1.0d0)

      call nwpw_timing_init()
      call dcopy(10,0.0d0,0,E,1)


*     **** get parallel variables ****
      call Parallel_Init()
      call Parallel_np(np)
      call Parallel_taskid(taskid)

      value = control_read(3,rtdb)
      if (.not. value) 
     > call errquit('error reading control',0, INPUT_ERR)

      call Parallel2d_Init(control_np_orbital())
      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)



      oprint= ((taskid.eq.MASTER).and.control_print(print_medium))
      lprint= ((taskid.eq.MASTER).and.control_print(print_low))
      hprint= ((taskid.eq.MASTER).and.control_print(print_high))


      if (taskid.eq.MASTER) call current_second(cpu1)

*     ***** print out header ****
      if (oprint) then
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
         call flush(luout)
      end if
      
      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0
      minimizer = control_minimizer()
      mapping   = control_mapping()


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


 
*     **** initialize lattice data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_Init()
      call D3dB_pfft_init()
      call ga_sync()


*     **** read ions ****
      value = ion_init(rtdb)
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

*     **** allocate psp data structure and read in psedupotentials into it ****
      call psp_init()
      call psp_readall()
      if (psp_semicore(0)) call semicore_check()

*     **** initialize G,mask,ke,coulomb,and ewald data structures ****
      call ke_init()
      if (control_version().eq.3) call coulomb_init()
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()
      call phafac()
      if (control_version().eq.3) then
         call ewald_init()
         call ewald_phafac()
      end if


*     **** read in wavefunctions and initialize psi ****
      if (.not.control_check_charge_multiplicity()) then
        call psi_new()
        newpsi = .true.
      else
        newpsi = .false.
      end if

*     **** Initialize 1d mappings for ne spaces ****
      call psi_get_ne(ispin,ne)
      mapping1d   = control_mapping1d()
      call Dne_init(ispin,ne,mapping1d)


*     **** read in wavefunctions and initialize psi ****
      value = psi_initialize()

      !call D3dB_n_fft_init(1,psi_ne(1)+psi_ne(2))

*     **** electron and geodesic data structures ****
      call electron_init()
      if ((minimizer.eq.4).or.(minimizer.eq.7)) then  !debug
        call geodesic2_init()
      else
         call geodesic_init()
      end if

*     **** initialize two-electron Gaussian integrals ****
      if (psp_pawexist()) then
         call nwpw_gintegrals_init()
         call nwpw_gintegrals_set((flag.gt.0))
      end if


*     **** initialize QM/MM ****
      call pspw_init_APC(rtdb)
      call pspw_qmmm_init(rtdb)
      call pspw_charge_init(rtdb)
      call pspw_Efield_init(rtdb,psi_ispin(),ne)

*     **** initialize SIC and HFX ****
      ne(1) = psi_ne(1)
      ne(2) = psi_ne(2)
      call pspw_init_SIC(rtdb,ne)
      call pspw_init_HFX(rtdb,psi_ispin(),ne)

*     **** initialize DFT+U ****
      call psp_U_init()

*     **** initialize  META GGA ****
      call nwpw_meta_gga_init(control_gga())

*     **** initialize vdw ****
      call vdw_DF_init()

*     **** initialize metadynamics ****
      call meta_initialize(rtdb)

*     **** initialize  META GGA ****

*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)

*     **** initialize linesearching ****
      call linesearch_init()





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
         write(luout,1119)  Parallel_maxthreads()

         if (control_parallel_io()) then
           write(luout,1123)
         else
           write(luout,1124)
         end if
         if (control_fast_erf()) write(luout,1125)

         write(luout,1115)
         write(luout,1121) control_boundry(),control_version()
         if (psi_ispin().eq.1) write(luout,1130) "restricted"
         if (psi_ispin().eq.2) write(luout,1130) "unrestricted"
         !if (qmmm) write(luout,1122)

         call v_bwexc_print(luout,control_gga())

         if (control_fractional()) write(luout,1132)
         call pspw_print_SIC(luout)
         call pspw_print_HFX(luout)
         call nwpw_meta_gga_print(luout)
         write(luout,1140)
         do ia = 1,ion_nkatm()
           call psp_print(ia)
c           call psp_check_print(ia)
         end do

         icharge = -(psi_ne(1)+psi_ne(psi_ispin()))
         en(1)     = psi_ne(1)
         en(psi_ispin()) = psi_ne(psi_ispin())
         if (control_fractional()) then
            icharge = 0
            do ms=1,psi_ispin()
            en(ms) =0.0d0
            do i=1,ne(ms)
              !!!!k = fweight(1)+(i-1)+(ms-1)*ne(1)
              icharge = icharge - (3-psi_ispin())*psi_occupation(ms,i)
              en(ms) = en(ms) + psi_occupation(ms,i)
            end do
            end do
         end if
         icharge = icharge + ion_TotalCharge()
         write(luout,1159) icharge

         if (control_version().eq.3) then
            rs  = (3.0d0*lattice_omega()/(4.0d0*pi))**(1.0d0/3.0d0)
            EVQ = 0.5d0*ewald_mandelung()*(icharge**2)/rs
         else
            EVQ = 0.0d0
         end if

         write(luout,1160)
         write(luout,1170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
         if (hprint) then
         write(luout,1180)
         do I=1,ion_nion()+ion_nion_q()
           if (I.gt.ion_nion()) then
           write(luout,1195) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                       nwpw_cosmo_qc0(I-ion_nion()),ion_amm(i)
           else if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
           else
           write(luout,1190) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz
         write(luout,1211) ion_nconstraints(),ion_ndof()

         call pspw_charge_Print(luout)
         end if
         call pspw_Efield_Print(luout)

         !write(luout,1220) psi_ne(1),psi_ne(psi_ispin()),' (Fourier space)'
         if (control_fractional()) then
          write(luout,1219) en(1),en(psi_ispin()),' (   fractional)'
          write(luout,1221) psi_ne(1),psi_neq(1),
     >                  psi_ne(psi_ispin()),psi_neq(psi_ispin()),
     >                  ' (Fourier space)'
         else
          write(luout,1220) psi_ne(1),psi_neq(1),
     >                  psi_ne(psi_ispin()),psi_neq(psi_ispin()),
     >                 ' (Fourier space)'
          write(luout,1221) psi_ne(1),psi_neq(1),
     >                  psi_ne(psi_ispin()),psi_neq(psi_ispin()),
     >                  ' (Fourier space)'
         end if
         write(luout,1230)
         write(luout,1233) control_cell_name()
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
         call lattice_abc_abg(a,b,c,alpha,beta,gamma)
         write(luout,1232) a,b,c,alpha,beta,gamma
         write(luout,1231) lattice_omega()
         write(luout,1250) lattice_ecut(),ngrid(1),ngrid(2),ngrid(3),
     >                 pack_nwave_all(0),pack_nwave(0)
         write(luout,1251) lattice_wcut(),ngrid(1),ngrid(2),ngrid(3),
     >                 pack_nwave_all(1),pack_nwave(1)
         if (control_smooth_cutoff()) 
     >      write(luout,1252) 
     >      control_smooth_cutoff_values(1)*lattice_wcut(),
     >      control_smooth_cutoff_values(2) 
         if (psp_pawexist()) then
            if (psp_paw_use_grid_cmp()) then
               write(luout,1262) psp_paw_mult_l_max(),
     >                           "grid"
            else
               write(luout,1262) psp_paw_mult_l_max(),
     >                           "two-electron Gaussian"
            end if
         end if
         if (control_version().eq.3) then
         write(luout,1260) ewald_rcut(),ewald_ncut()
         write(luout,1261) ewald_mandelung(),
     >                     ewald_mandelung()
     >                    *(4.0d0*pi/3.0)**(1.0d0/3.0d0),rs
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
         write(luout,1280) control_time_step(),control_fake_mass()
         write(luout,1290) control_tole(),control_tolc()
         write(luout,1281) control_it_in()*control_it_out(),
     >                 control_it_in(),control_it_out()
         if (control_gram_schmidt()) write(luout,1282)
         if ((minimizer.eq.5).or.(minimizer.eq.8)) then
           write(luout,1291)
           if (control_ks_algorithm().eq.-1)  then
             if (control_precondition()) then
                write(luout,1292) 
     >             "block preconditioned conjugate gradient"
             else
                write(luout,1292) "block conjugate gradient"
             end if
           end if
           if (control_ks_algorithm().eq.0)  then
             if (control_precondition()) then
                write(luout,1292) "preconditioned conjugate gradient"
             else
                write(luout,1292) "conjugate gradient"
             end if
           end if
           if (control_ks_algorithm().eq.1) 
     >       write(luout,1292) "rmm-diis"

           write(luout,1295) control_ks_maxit_orb(),
     >                       control_ks_maxit_orbs()
           if (control_scf_algorithm().eq.0) 
     >       write(luout,1293) "simple mixing"
           if (control_scf_algorithm().eq.1) 
     >       write(luout,1293) "Broyden mixing"
           if (control_scf_algorithm().eq.2) 
     >       write(luout,1289) "Johnson-Pulay mixing",
     >                         control_diis_histories()
           if (control_scf_algorithm().eq.3) 
     >       write(luout,1293) "Anderson density mixing"
           if (minimizer.eq.5) write(luout,1296) "potential"
           if (minimizer.eq.8) write(luout,1296) "density"
           write(luout,1294) control_ks_alpha()
           write(luout,1301) control_kerker_g0()
         end if
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
c         call ion_Print_neighborlist()
         call util_flush(luout)
      end if

*                |***************************|
******************     call CG minimizer     **********************
*                |***************************|
      if (taskid.eq.MASTER) call current_second(cpu2)

c     if (flag.eq.1) then
c       call control_gradient_iterations()
c     end if

*     **** calculate energy ****
      if (flag.eq.-1) then
        EV= cgsd_noit_energy()
      else
        EV= cgsd_energy(newpsi)
      end if

*     **** calculate excited state orbitals ****
      call ga_sync()
      call cgsd_excited()

*     **** calculate oep eigenvalues ****
      call cgsd_oep_eigenvalues()

*     **** extra energy output for QA test ****
      if (lprint) write(LuOut,1600) EV
      if (lprint) then
         if (EVQ.gt.0.0d0) write(LuOut,1601) EV+EVQ
      end if

*     **** calculate the spin contamination ****
      if (flag.gt.-1) call psi_spin2(dipole(1))

      
*     **** calculate the dipole ***
      dipole(1) = 0.0d0
      dipole(2) = 0.0d0
      dipole(3) = 0.0d0
      cdipole(1) = 0.0d0
      cdipole(2) = 0.0d0
      cdipole(3) = 0.0d0
      if (flag.gt.-1) then
        call rho_dipole(dipole)
        if (control_periodic_dipole()) call psi1_crystal_dipole(cdipole)

c        !call psi_polariz() !early cubic implementation
      end if

*     ****  calculate gradient ***
      if (flag.gt.0) then
      nion = ion_nion()
      if (pspw_charge_found().and.
     >    (.not.pspw_bqext())) nion = nion + pspw_charge_nion()

      value = BA_push_get(mt_dbl,(3*nion),
     >                       'fion',fion(2),fion(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
     
      call cgsd_energy_gradient(dbl_mb(fion(1)))
      call dscal(3*nion,(-1.0d0),dbl_mb(fion(1)),1)
      end if

*     **** calculate the stress tensor ****
      call dcopy(9,0.0d0,0,stress,1)
      call dcopy(6,0.0d0,0,lstress,1)
      if (flag.eq.3) then
         call psp_stress_init()
         call psp_stress_readall()
         call cgsd_energy_stress(stress,lstress)
         call psp_stress_end()
      end if

*     **** run socket client ****
      if (control_runsocket().and.(flag.gt.0)) then
         call runsocket_ipi(nion,EV,
     >                      dbl_mb(ion_rion_ptr()),dbl_mb(fion(1)))
      end if


*     *************************************************************
*     **** output energy, dipole, and gradient to rtdb for use ****
*     **** by task_energy and task_gradient                    ****
*     *************************************************************
      value = .true.
      if (flag.gt.-1) then
      value = btdb_put(rtdb,'pspw:energy',mt_dbl,1,EV)
      value = value.and.
     >        btdb_put(rtdb,'pspw:dipole',mt_dbl,
     >                 3,dipole)
      value = value.and.
     >        btdb_put(rtdb,'pspw:crystal_dipole',mt_dbl,
     >                 3,cdipole)
      value = value.and.control_gas_energy_set(EV)
      end if
      if (flag.gt.0) then
        value = value.and.
     >        btdb_put(rtdb,'pspw:gradient',mt_dbl,
     >                 3*nion,dbl_mb(fion(1)))
        value = value.and.BA_pop_stack(fion(2))
      end if
      if (flag.eq.3) then
        value = value.and.
     >        btdb_put(rtdb,'pspw:stress',mt_dbl,
     >                 9,stress)
        value = value.and.
     >        btdb_put(rtdb,'pspw:lstress',mt_dbl,
     >                 6,lstress)
      end if
      if (.not. value) call errquit('cgsdv5: error writing rtdb',0,
     &       RTDB_ERR)

      if (taskid.eq.MASTER) call current_second(cpu3)

*                |***************************|
******************         Epilogue          **********************
*                |***************************|

*     **** calculate Mulliken Populations ***
      if (control_Mulliken()) then
         call psi_Mulliken(rtdb)
         call psi_Born()
         i = control_excited_ne(1)+control_excited_ne(2)
         if (i.eq.0) call kp_structure(rtdb)
      end if


*     **** write wavefunctions to file and finalize psi ****
      if (flag.eq.-1) then
        value = psi_finalize(.false.)
      else
        value = psi_finalize(.true.)
      end if


*     **** write geometry to rtdb ****
      call pspw_charge_write(rtdb)
      call ion_write(rtdb)


*     **** deallocate heap memory ****
      call electron_finalize()
      if ((minimizer.eq.4).or.(minimizer.eq.7)) then
        call geodesic2_finalize()
      else
        call geodesic_finalize()
      end if
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
      call pspw_end_SIC()
      call pspw_end_HFX()
      call pspw_end_APC()
      call pspw_charge_end()
      call pspw_Efield_end()
      call pspw_qmmm_end()
c      call frac_occ_end()
      call meta_finalize(rtdb)
      if (psp_pawexist()) call nwpw_gintegrals_end()

      call ion_end()
      call psp_end()
      call ion_end_FixIon()
      !call D3dB_n_fft_end(1)
      call D3dB_pfft_end()
      call D3dB_end(1)
      if (control_version().eq.4)    call D3dB_end(2)
      if (control_has_ngrid_small()) call D3dB_end(3)
      call Dne_end()
      call psi_data_end()

*                |***************************|
****************** report consumed cputime   **********************
*                |***************************|
      if (lprint) then
         CALL current_second(cpu4)

         T1=CPU2-CPU1
         T2=CPU3-CPU2
         T3=CPU4-CPU3
         T4=CPU4-CPU1
         AV=T2/dble(electron_count())
         write(luout,1801)
         write(luout,1802)
         write(luout,1803) T1
         write(luout,1804) T2
         write(luout,1805) T3
         write(luout,1806) T4
         write(luout,1807) AV,electron_count(),linesearch_count()

         call nwpw_timing_print_final(oprint,electron_count())
         write(luout,*)
         call pspw_print_HFX_counters(luout)
         CALL nwpw_MESSAGE(4)
      end if 



      call Parallel2d_Finalize()
      call Parallel_Finalize()
      call ecce_print_module_exit('pspw', 'ok')
      cgsdv5 = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*               NWPW PSPW Calculation              *')
 1030 FORMAT(10X,'* [ (Grassmann/Stiefel manifold implementation)  ] *')
 1035 FORMAT(10x,'*      [ NorthWest Chemistry implementation ]      *')
 1040 FORMAT(10X,'*            version #5.10   06/12/02              *')
 1041 FORMAT(10X,'*    This code was developed by Eric J. Bylaska,   *')
 1042 FORMAT(10X,'*    and was based upon algorithms and code        *')
 1043 FORMAT(10X,'*    developed by the group of Prof. John H. Weare *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1111 FORMAT(/' number of processors used:',I10)
 1112 FORMAT( ' parallel mapping         :   1d-slab')
 1113 FORMAT( ' parallel mapping         :2d-hilbert')
 1114 FORMAT( ' parallel mapping         :  balanced')
 1115 FORMAT(/' options:')
 1116 FORMAT( ' parallel mapping         : not balanced')
 1117 FORMAT( ' processor grid           :',I4,' x',I4)
 1118 FORMAT( ' parallel mapping         : 2d-hcurve')
 1119 FORMAT( ' number of threads        :',I10)
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1122 FORMAT(5X,' qmmm simulation')
 1123 FORMAT( ' parallel io              :        on')
 1124 FORMAT( ' parallel io              :       off')
 1125 FORMAT( ' fast erf                 :        on')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1132 FORMAT(5X,' using fractional occupation')
 1140 FORMAT(/' elements involved in the cluster:')
 1150 FORMAT(5X,I2,': ',A4,'  core charge:',F8.4,'  lmax=',I3)
 1151 FORMAT(5X,'        cutoff =',4F8.3)
 1152 FORMAT(12X,' highest angular component      : ',i3)
 1153 FORMAT(12X,' local potential used           : ',i3)
 1154 FORMAT(12X,' number of non-local projections: ',i3)
 1155 FORMAT(12X,' semicore corrections included  : ',
     >       F6.3,' (radius) ',F6.3,' (charge)')
 1156 FORMAT(12X,' aperiodic cutoff radius        : ',F6.3)
 1157 FORMAT(12X,' comment    : ',A)
 1158 FORMAT(12X,' pseudpotential type            : ',i3)

 1159 FORMAT(/' total charge:',F8.3)
 1160 FORMAT(/' atomic composition:')
 1170 FORMAT(7(5X,A4,':',I5))
 1180 FORMAT(/' position of ions (au):')
 1190 FORMAT(5X, I4, A5  ,' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1191 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - fixed',A)
 1193 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - z fixed')
 1194 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,A)
 1195 FORMAT(5X, I4, A5, ' (',3F11.5,' ) -      charge= ',E13.6,' ',A)
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1211 FORMAT(5X,'   number of constraints = ', I6,' ( DOF = ',I6,' )' )
 1219 FORMAT(/' number of electrons: spin up=',F6.2,16x,
     >        '  down=',F6.2,A)
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
 1231 FORMAT(5x,'             omega=',F12.1)
 1232 FORMAT(5x,' lattice:    a=    ',f8.3,' b=   ',f8.3,' c=    ',f8.3,
     >      /5x,'             alpha=',f8.3,' beta=',f8.3,' gamma=',f8.3)
 1233 FORMAT(5x,' cell_name:  ',A)
 1241 FORMAT(5x,' lattice:    a1=<',3f8.3,' >')
 1242 FORMAT(5x,'             a2=<',3f8.3,' >')
 1243 FORMAT(5x,'             a3=<',3f8.3,' >')
 1244 FORMAT(5x,' reciprocal: b1=<',3f8.3,' >')
 1245 FORMAT(5x,'             b2=<',3f8.3,' >')
 1246 FORMAT(5x,'             b3=<',3f8.3,' >')

 1250 FORMAT(/5X,' density cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
 1251 FORMAT(5X,' wavefnc cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
 1252 FORMAT(5X,' smooth wavefnc cutoff:  A=',F7.3,' sigma = ',F7.3)
 1260 FORMAT(5X,' Ewald summation: cut radius=',F8.2,'  and',I3)
 1261 FORMAT(5X,'                  Madelung Wigner-Seitz=',f14.8,
     >                              ' (alpha=',f14.8,' rs=',f11.8,')')
 1262 FORMAT(5X,' paw integration: xc lmax=',I3,
     >          ' compensation charge= ',A)
 1270 FORMAT(/' technical parameters:')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'fictitious mass=',F10.1)
 1281 FORMAT(5X, ' maximum iterations =',I10,
     >           ' ( ',I4,' inner ',I6,' outer )')
 1282 FORMAT(5X, ' Gram Schmidt orthonormality used')
 1289 FORMAT(5X, ' SCF algorithm        = ',A,' (',I2,' histories)')
 1290 FORMAT(5X, ' tolerance=',E9.3,' (energy)',E12.3,
     &        ' (density)')
 1291 FORMAT(/' Kohn-Sham scf parameters:')
 1292 FORMAT(5X, ' Kohn-Sham algorithm  = ',A)
 1293 FORMAT(5X, ' SCF algorithm        = ',A)
 1294 FORMAT(5X, ' SCF mixing parameter =',F7.4)
 1295 FORMAT(5X, ' Kohn-Sham iterations = ',I3,
     >           ' (',I3,' outer)')
 1296 FORMAT(5X, ' SCF mixing type      = ',A)
 1297 FORMAT(/' fractional smearing parameters:')
 1298 FORMAT(5X, ' smearing algorithm   = ',A)
 1299 FORMAT(5X, ' smearing parameter   = ',E9.3,' (',F7.1,' K)'/,
     >       5X, ' mixing parameter     =',F7.4) 
 1300 FORMAT(//)
 1301 FORMAT(5X, ' Kerker damping       =',F7.4)
 1305 FORMAT(10X,'================ iteration =========================')
 1310 FORMAT(I8,E20.10,3E15.5)
 1320 FORMAT(' number of electrons: spin up=',F11.5,'  down=',F11.5,A)
 1330 FORMAT(/' comparison between hamiltonian and lambda matrix')
 1340 FORMAT(I3,2I3,' H=',E16.7,', L=',E16.7,', H-L=',E16.7)
 1350 FORMAT(/' orthonormality')
 1360 FORMAT(I3,2I3,E18.7)
 1370 FORMAT(I3)
 1380 FORMAT(' ''',a,'''',I4)
 1390 FORMAT(I3)
 1400 FORMAT(I3,3E18.8/3X,3E18.8)
 1410 FORMAT(10X,'=============  summary of results  =================')
 1420 FORMAT( ' final position of ions:')
 1430 FORMAT(/' total     energy    :',E19.10,' (',E15.5,'/ion)')
 1431 FORMAT(/' total paw energy    :',E19.10,' (',E15.5,'/ion)')
 1440 FORMAT( ' total orbital energy:',E19.10,' (',E15.5,'/electron)')
 1450 FORMAT( ' hartree   energy    :',E19.10,' (',E15.5,'/electron)')
 1460 FORMAT( ' exc-corr  energy    :',E19.10,' (',E15.5,'/electron)')
 1470 FORMAT( ' ion-ion   energy    :',E19.10,' (',E15.5,'/ion)')
 1480 FORMAT(/' K.S. kinetic energy :',E19.10,' (',E15.5,'/electron)')
 1490 FORMAT( ' K.S. V_l  energy    :',E19.10,' (',E15.5,'/electron)')
 1495 FORMAT( ' K.S. V_nl energy    :',E19.10,' (',E15.5,'/electron)')
 1496 FORMAT( ' K.S. V_Hart energy  :',E19.10,' (',E15.5,'/electron)')
 1497 FORMAT( ' K.S. V_xc energy    :',E19.10,' (',E15.5,'/electron)')
 1498 FORMAT( ' Virial Coefficient  :',E19.10)
 1500 FORMAT(/' orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1600 FORMAT(/' Total PSPW energy     :',E19.10)
 1601 FORMAT(/' Total PSPW energy ',
     >        '(charge correction - E + 0.5*q**2*mandelung/rs) :',
     >        E19.10)
 1801 FORMAT(//'== Timing ==')
 1802 FORMAT(/'cputime in seconds')
 1803 FORMAT( '  prologue    : ',E14.6)
 1804 FORMAT( '  main loop   : ',E14.6)
 1805 FORMAT( '  epilogue    : ',E14.6)
 1806 FORMAT( '  total       : ',E14.6)
 1807 FORMAT( '  cputime/step: ',E14.6,
     >        '       (',I8,' evalulations,', I8,' linesearches)')
 1808 FORMAT(A,E14.6,E14.6)
 1809 FORMAT(//A,2A14)

 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(luout,9010) ierr
      call Parallel_Finalize()
      call ecce_print_module_exit('pspw', 'failed')

      cgsdv5 = value
      return
      END
