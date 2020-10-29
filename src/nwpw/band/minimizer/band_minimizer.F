*
* $Id$
*
***********************************************************************
*                     band_minimizer								  *
*                                                                     *
*     This is a developing band structure parallel code for NWCHEM    *
*       + tcgmsg message passing library used                         *
*       + my own slap-decomposed parallel 3d-FFT(real->complex) used  *
*                                                                     *
*                                                                     *
***********************************************************************

      logical function band_minimizer(rtdb,flag)
      implicit none
      integer rtdb
      integer flag

#include "global.fh"
#include "bafdecls.fh"
#include "inp.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"      

*     **** parallel variables ****
      integer  taskid,np,np_i,np_j,np_k
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,nfft3d
      real*8  a,b,c,alpha,beta,gamma

*     ***** energy variables ****
      logical spin_orbit
      integer ispin,ne(2)
      real*8  E(50)
      real*8  dipole(3)
      real*8  stress(3,3),lstress(6)

*     **** gradient variables ****
      integer fion(2)

*     **** error variables ****
      logical value
      integer ierr

*     **** local variables ****
      logical newpsi,lprint,mprint,hprint
      real*8  gx,gy,gz,cx,cy,cz
      real*8  EV,pi,weight
      real*8  icharge,en(2)
      real*8  f0,f1,f2,f3,f4,f5,f6
      integer if1,if2
      integer i,k,ia,nion,vers,nbrillioun,nb,ms
      integer mapping,minimizer,mapping1d

*     **** external functions ****
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg,ion_amass,ion_TotalCharge
      character*4 ion_aname
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg,ion_amass,ion_TotalCharge
      external    ion_aname


      real*8   control_tole,control_tolc,control_tolr,ion_rion
      external control_tole,control_tolc,control_tolr,ion_rion
      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      logical  control_read,control_move,ion_init,ion_q_FixIon
      external control_read,control_move,ion_init,ion_q_FixIon
      logical control_spin_orbit
      external control_spin_orbit
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      integer  ion_nkatm
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      external ion_nkatm

      character*12 control_boundry
      external     control_boundry

      logical  brillioun_print,control_Mulliken
      integer  brillioun_nbrillioun,brillioun_nbrillq
      integer  brillioun_nbrillioun0
      real*8   brillioun_weight_brdcst,brillioun_ks_brdcst
      real*8   brillioun_k_brdcst,brillioun_weight
      external brillioun_print,control_Mulliken
      external brillioun_nbrillioun,brillioun_nbrillq
      external brillioun_nbrillioun0
      external brillioun_weight_brdcst,brillioun_ks_brdcst 
      external brillioun_k_brdcst,brillioun_weight
      integer  c_electron_count,linesearch_count
      external c_electron_count,linesearch_count

      real*8   nwpw_timing
      external nwpw_timing
      integer  Cram_nwave_all_brdcst,Cram_nwave_brdcst
      external Cram_nwave_all_brdcst,Cram_nwave_brdcst

      integer  ewald_ncut
      real*8   ewald_rcut,ewald_mandelung,ewald_e
      external ewald_ncut
      external ewald_rcut,ewald_mandelung,ewald_e
      logical  cpsp_semicore,psi_filefind,cpsi_initialize,cpsi_finalize
      external cpsp_semicore,psi_filefind,cpsi_initialize,cpsi_finalize
      real*8   c_cgsd_energy,c_cgsd_noit_energy
      external c_cgsd_energy,c_cgsd_noit_energy
      integer  cpsp_lmax,cpsp_locp,cpsp_lmmax
      external cpsp_lmax,cpsp_locp,cpsp_lmmax
      integer  cpsp_nprj,cpsp_psp_type
      external cpsp_nprj,cpsp_psp_type
      real*8   cpsp_rcore,cpsp_rc,cpsp_ncore,cpsp_zv
      external cpsp_rcore,cpsp_rc,cpsp_ncore,cpsp_zv
      character*4 ion_atom
      external    ion_atom

      character*255 cpsp_comment,comment
      external      cpsp_comment

      integer  cpsi_ispin,cpsi_ne,psi_get_version
      external cpsi_ispin,cpsi_ne,psi_get_version
      logical  pspw_reformat_c_wvfnc,control_check_charge_multiplicity
      external pspw_reformat_c_wvfnc,control_check_charge_multiplicity
      integer  control_mapping,control_minimizer,control_np_dimensions
      external control_mapping,control_minimizer,control_np_dimensions
      integer  control_mapping1d
      external control_mapping1d
      integer  control_scf_algorithm,control_ks_algorithm
      external control_scf_algorithm,control_ks_algorithm
      real*8   control_ks_alpha,control_kerker_g0
      external control_ks_alpha,control_kerker_g0
      integer  control_ks_maxit_orb,control_ks_maxit_orbs
      external control_ks_maxit_orb,control_ks_maxit_orbs
      real*8   cpsi_occupation,control_fractional_temperature
      external cpsi_occupation,control_fractional_temperature
      logical  control_fractional,control_print
      external control_fractional,control_print
      integer  control_fractional_smeartype
      external control_fractional_smeartype
      real*8   control_fractional_kT,control_fractional_alpha
      external control_fractional_kT,control_fractional_alpha

c      character*255 cpsp_comment,comment
c      external      cpsp_comment

      integer  ion_nconstraints,ion_ndof
      external ion_nconstraints,ion_ndof
      logical  control_smooth_cutoff
      external control_smooth_cutoff
      real*8   control_smooth_cutoff_values
      external control_smooth_cutoff_values
      logical  nwpwxc_is_on,control_only_lda,control_precondition
      external nwpwxc_is_on,control_only_lda,control_precondition


*****************************|  PROLOGUE  |****************************

      value = .true.
      pi = 4.0d0*datan(1.0d0)

      call nwpw_timing_init()
      call dcopy(50,0.0d0,0,E,1)


*     **** get parallel variables ****
      call Parallel_Init()
      call Parallel_np(np)
      call Parallel_taskid(taskid)

      value = control_read(5,rtdb)
      if (.not. value) 
     > call errquit('error reading control',0, DISK_ERR)

      lprint = ((taskid.eq.MASTER).and.(control_print(print_low)))
      mprint = ((taskid.eq.MASTER).and.(control_print(print_medium)))
      hprint = ((taskid.eq.MASTER).and.(control_print(print_high)))

      if (taskid.eq.MASTER) call current_second(cpu1)
*     ***** print out header ****
      if (mprint) then
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
         write(luout,1044)
         write(luout,1010)
         write(luout,1000)
         call nwpw_message(1)
         write(luout,1110)
         call flush(luout)
      end if
      

      call Parallel3d_Init(control_np_dimensions(2),
     >                     control_np_dimensions(3))
      call Parallel3d_np_i(np_i)
      call Parallel3d_np_j(np_j)
      call Parallel3d_np_k(np_k)

      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0
      minimizer = control_minimizer()
      mapping   = control_mapping()
      ierr = 0 

*     **** initialize C3dB data structure ****
      call C3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call C3dB_nfft3d(1,nfft3d)
       
      call cpsi_data_init(20)

*     **** read ions ****
      value = ion_init(rtdb)
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

*     **** initialize lattice data structure ****
      call lattice_init()
      call c_G_init()

*     **** nitialize brillion ***
      call brillioun_init()
      call Cram_Init()
      call C3dB_pfft_init()

*     **** initialize D3dB data structure and mask for GGA ****
c      if (((control_gga().ge.10).and.(control_gga().le.200))
c     >     .or.nwpwxc_is_on()) THEN
      if (.not.control_only_lda()) then
         call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
         call G_init()
         call mask_init()
      end if

c*     **** read ions ****
c      value = ion_init(rtdb)
c      call center_geom(cx,cy,cz)
c      call center_mass(gx,gy,gz)

*     **** allocate psp data structure and read in psedupotentials into it ****
      call cpsp_init()
      call cpsp_readall()
      if (cpsp_semicore(0)) call c_semicore_check()


*     **** initialize ke,and coulomb data structures ****
      call cstrfac_init()
      call cke_init()
      call c_coulomb_init()
      call ewald_init()

*     **** set up phase factors at the current geometry  ****
      call cphafac()
      call cphafac_k()
      call ewald_phafac()

*     **** read in wavefunctions and initialize psi ****
      if (.not.control_check_charge_multiplicity()) then
         call cpsi_new()
         newpsi = .true.
      else
         newpsi = .false.
*        **** convert from pspw format to band format ****
         vers = psi_get_version()
         if ((vers.eq.3).or.(vers.eq.4)) then
            newpsi = .true.
            value = btdb_parallel(.false.)
            if (taskid.eq.MASTER) then
              value= pspw_reformat_c_wvfnc(1)  
            end if
            value = btdb_parallel(.true.)
         end if
      end if

      call psi_get_ne(ispin,ne)
      if (ispin.eq.3) then
         spin_orbit = .true.
         ispin=2
      else
         spin_orbit = .false.
      end if
      nbrillioun = brillioun_nbrillioun()
      mapping1d  = control_mapping1d()
      call Pneb_init(ispin,ne,nbrillioun,mapping1d,spin_orbit)
      value = cpsi_initialize(.true.)
      !newpsi = newpsi.or.value

*     **** electron and geodesic data structures ****
      call c_electron_init()
      call c_geodesic_init()


*     **** initialize QM/MM ****
      call band_init_APC(rtdb)


*     **** initialize HFX ****
      call band_init_HFX(rtdb,nbrillioun,ispin,ne)


*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)
 
*     **** initialize linesearching ****
      call linesearch_init()


*                |**************************|
******************   summary of input data  **********************
*                |**************************|


*     **** determine en ****
      if (.not.control_spin_orbit()) then
        icharge = 0.0d0
        en(1)   = 0.0d0
        en(2)   = 0.0d0
        b = dble(3-cpsi_ispin())
        do nb=1,brillioun_nbrillq()
        weight = brillioun_weight(nb)
        do ms=1,cpsi_ispin()
          do i=1,ne(ms)
            a = cpsi_occupation(nb,ms,i)
            icharge = icharge - b*a*weight
            en(ms)  = en(ms) + a*weight
          end do
        end do
        end do
        call K1dB_Vector_SumAll(2,en)
        call K1dB_SumAll(icharge)
      else
        icharge          = -cpsi_ne(1)
        en(1)            =  cpsi_ne(1)
        en(cpsi_ispin()) =  cpsi_ne(cpsi_ispin())
      end if


      if (mprint) then
         write(luout,1111) np
         write(luout,1117) np_i,np_j,np_k
         if (mapping.eq.1) write(luout,1112)
         if (mapping.eq.2) write(luout,1113)
         if (mapping.eq.3) write(luout,1118)

         write(luout,1115)
         write(luout,1121) control_boundry(),control_version()
         if (.not.spin_orbit) then
           if (cpsi_ispin().eq.1) write(luout,1130) "restricted"
           if (cpsi_ispin().eq.2) write(luout,1130) "unrestricted"
         else
           write(luout,1130) "spin orbit"
         end if

         call v_bwexc_print(luout,control_gga())
         
         call band_print_HFX(luout)
         write(luout,1140)
         do ia = 1,ion_nkatm()
           write(luout,1150) ia,ion_atom(ia),
     >                    cpsp_zv(ia),cpsp_lmax(ia)
           comment = cpsp_comment(ia)
           i = inp_strlen(comment)
           write(luout,1157) comment(1:i)
           write(luout,1158) cpsp_psp_type(ia)
           write(luout,1152) cpsp_lmax(ia)
           write(luout,1153) cpsp_locp(ia)
c           write(luout,1154) cpsp_lmmax(ia)
           write(luout,1154) cpsp_nprj(ia)
           if (cpsp_semicore(ia))
     >         write(luout,1155) cpsp_rcore(ia),cpsp_ncore(ia)
           write(luout,1151) (cpsp_rc(i,ia),i=0,cpsp_lmax(ia))
         end do

         icharge = icharge + ion_TotalCharge()
         write(luout,1159) icharge


         write(luout,1160)
         write(luout,1170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
         if (hprint) then
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
         write(luout,1211) ion_nconstraints(),ion_ndof()
         endif
         write(luout,1220) en(1),en(cpsi_ispin()),
     >                 ' (Fourier space)'
         write(luout,1221) cpsi_ne(1),cpsi_ne(cpsi_ispin()),
     >                 ' (Fourier space)'

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
         call lattice_abc_abg(a,b,c,alpha,beta,gamma)
         write(luout,1232) a,b,c,alpha,beta,gamma
         write(luout,1231) lattice_omega()
         write(luout,1260) ewald_rcut(),ewald_ncut()
         write(luout,1261) ewald_mandelung()

         write(luout,1255)
         write(luout,1256) brillioun_nbrillioun(),
     >                     brillioun_nbrillioun0()
      end if

c     **** print brillioun zone - extra logic for distributed kpoints ****
      if (brillioun_print()) then
         do i=1,brillioun_nbrillioun()
            f0 = brillioun_weight_brdcst(i) 
            f1 = brillioun_ks_brdcst(1,i) 
            f2 = brillioun_ks_brdcst(2,i) 
            f3 = brillioun_ks_brdcst(3,i) 
            f4 = brillioun_k_brdcst(1,i)  
            f5 = brillioun_k_brdcst(2,i)  
            f6 = brillioun_k_brdcst(3,i)
            if (mprint) write(luout,1257) f0,f1,f2,f3,f4,f5,f6
         end do
      else
        if (mprint) write(luout,1258) 
      end if

      if1 = Cram_nwave_all_brdcst(0)
      if2 = Cram_nwave_brdcst(0)
      if (mprint) then
         write(luout,1249)
         write(luout,1250) lattice_ecut(),ngrid(1),ngrid(2),ngrid(3),
     >                     if1,if2
      end if

      if (brillioun_print()) then
        do i=1,brillioun_nbrillioun()
          if1 = Cram_nwave_all_brdcst(i)
          if2 = Cram_nwave_brdcst(i)
          if (mprint) then
          write(luout,1251) i,lattice_wcut(),ngrid(1),ngrid(2),ngrid(3),
     >                      if1,if2
          end if
        end do
      else
        if (mprint) write(luout,1252) lattice_wcut()
      end if

      if (mprint) then
         if (control_smooth_cutoff())
     >      write(luout,1253)
     >      control_smooth_cutoff_values(1)*lattice_wcut(),
     >      control_smooth_cutoff_values(2)

         call c_rho_symmetrizer_print(luout)

         write(luout,1270)
         write(luout,1280) control_time_step(),control_fake_mass()
         write(luout,1290) control_tole(),control_tolc()
         write(luout,1281) control_it_in()*control_it_out(),
     >                 control_it_in(),control_it_out()

         if ((minimizer.eq.5).or.(minimizer.eq.8)) then
           write(luout,1291)
           if (control_ks_algorithm().eq.0)  then
              if (control_precondition()) then
                 write(luout,1292) "preconditioned conjugate gradient"
              else
                 write(luout,1292) "conjugate gradient"
              end if
           end if
           if (control_ks_algorithm().eq.1) write(luout,1292) "rmm-diis"
           write(luout,1295) control_ks_maxit_orb(),
     >                       control_ks_maxit_orbs()
           if (control_scf_algorithm().eq.0)
     >       write(luout,1293) "simple mixing"
           if (control_scf_algorithm().eq.1)
     >       write(luout,1293) "Anderson potential mixing"
           if (control_scf_algorithm().eq.2)
     >       write(luout,1293) "Johnson-Pulay mixing"
           if (control_scf_algorithm().eq.3)
     >       write(luout,1293) "Anderson density mixing"
           if (minimizer.eq.5) write(luout,1296) "potential"
           if (minimizer.eq.8) write(luout,1296) "density"
           write(luout,1294) control_ks_alpha()
           write(luout,1301) control_kerker_g0()
         end if
         if (control_fractional()) then
           write(luout,1297)
           if (control_fractional_smeartype().eq.0)
     >       write(luout,1298) "step function"
           if (control_fractional_smeartype().eq.1)
     >       write(luout,1298) "Fermi-Dirac"
           if (control_fractional_smeartype().eq.2)
     >       write(luout,1298) "Gaussian"
           if (control_fractional_smeartype().eq.4)
     >       write(luout,1298) "Marzari-Vanderbilt"
           write(luout,1299) control_fractional_kT(),
     >                   control_fractional_temperature(),
     >                   control_fractional_alpha()
         end if
         write(luout,1300)
         call util_flush(luout)
      end if



*                |***************************|
******************     call CG minimizer     **********************
*                |***************************|
      if (taskid.eq.MASTER) call current_second(cpu2)


*     **** calculate energy ****
      if (flag.eq.-1) then

        EV= c_cgsd_noit_energy()
      else

        EV= c_cgsd_energy(newpsi)
      end if

*     **** calculate excited state orbitals ****
      call ga_sync()
      call c_cgsd_excited()

*     **** extra energy output for QA test ****
      if (lprint) write(luout,1600) EV

*     **** calculate the spin contamination ****
      if (flag.gt.-1) call cpsi_spin2(dipole(1))
      
*     **** calculate the dipole ***
      dipole(1) = 0.0d0
      dipole(2) = 0.0d0
      dipole(3) = 0.0d0

*     ****  calculate gradient ***
      if (flag.gt.0) then
      nion = ion_nion()
      value = BA_push_get(mt_dbl,(3*nion),
     >                       'fion',fion(2),fion(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


      call c_cgsd_energy_gradient(dbl_mb(fion(1)))
      call dscal(3*nion,(-1.0d0),dbl_mb(fion(1)),1)
      end if
     
*     **** calculate the stress tensor ****
      call dcopy(9,0.0d0,0,stress,1)
      call dcopy(6,0.0d0,0,lstress,1)
      if (flag.eq.3) then

         call cpsp_stress_init()
         call cpsp_stress_readall()
         call c_cgsd_energy_stress(stress,lstress)
         call cpsp_stress_end()
      end if


*     *************************************************************
*     **** output energy, dipole, and gradient to rtdb for use ****
*     **** by task_energy and task_gradient                    ****
*     *************************************************************
      value = .true.
      if (flag.gt.-1) then
      value = btdb_put(rtdb,'band:energy',mt_dbl,1,EV)
      value = value.and.
     >        btdb_put(rtdb,'band:dipole',mt_dbl,
     >                 3,dipole)
      end if
      if (flag.gt.0) then
        value = value.and.btdb_put(rtdb,'band:gradient',mt_dbl,
     >                             3*nion,dbl_mb(fion(1)))
        value = value.and.BA_pop_stack(fion(2))
      end if
      if (flag.eq.3) then
        value = value.and.
     >        btdb_put(rtdb,'band:stress',mt_dbl,
     >                 9,stress)
        value = value.and.
     >        btdb_put(rtdb,'band:lstress',mt_dbl,
     >                 6,lstress)
      end if
      if (.not. value) 
     > call errquit('band_minimizer: error writing rtdb',0, RTDB_ERR)
     


      if (taskid.eq.MASTER) call current_second(cpu3)

*                |***************************|
******************         Epilogue          **********************
*                |***************************|

*     **** generate eig grid data file ****
      call cpsi_printeigall(rtdb)

*     **** calculate interpolated band structure plot ***
      call band_interpolate_structure(rtdb)
      call band_fft_structure(rtdb)
      !call band_kp_structure(rtdb)


*     **** calculate Mulliken Populations ***
      if (control_Mulliken()) then
         call cpsi_Mulliken(rtdb)
      end if


*     **** write geometry to rtdb ****
      call ion_write(rtdb)

*     **** write wavefunctions to file and finalize psi ****
      if (flag.eq.-1) then
        value = cpsi_finalize(.false.)
      else
        value = cpsi_finalize(.true.)
      end if

*     **** deallocate heap memory ****
      call c_electron_finalize()
      call c_geodesic_finalize()
      call ewald_end()
      call cstrfac_end()
      call c_coulomb_end()
      call cke_end()
      call cpsp_end()
      call Cram_end()
      call c_G_end()
      call band_end_HFX()
      call brillioun_end()
      call band_end_APC()
      call ion_end()
      call ion_end_FixIon()
      call cpsi_data_end()
      call C3dB_pfft_end()
      call C3dB_end(1)
      call Pneb_end()
c      IF (((control_gga().ge.10).and.(control_gga().le.200)) 
c     >     .or.nwpwxc_is_on()) THEN
      if (.not.control_only_lda()) then
         call mask_end()
         call G_end()
         call D3dB_end(1)
      end if

*                |***************************|
****************** report consumed cputime   **********************
*                |***************************|
      if (lprint) then
         CALL current_second(cpu4)

         T1=CPU2-CPU1
         T2=CPU3-CPU2
         T3=CPU4-CPU3
         T4=CPU4-CPU1
         AV=T2/dble(c_electron_count())
         write(luout,1801)
         write(luout,1802)
         write(luout,1803) T1
         write(luout,1804) T2
         write(luout,1805) T3
         write(luout,1806) T4
         write(luout,1807) AV,c_electron_count(),linesearch_count()

         call nwpw_timing_print_final(mprint,c_electron_count())
         write(luout,*)
         CALL nwpw_MESSAGE(4)
      end if 


      call Parallel3d_Finalize()
      call Parallel_Finalize()
      band_minimizer = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,
     > '**********************************************************')
 1010 FORMAT(10X,
     > '*                                                        *')
 1020 FORMAT(10X,
     > '*                   NWPW BAND Calculation                *')
 1030 FORMAT(10X,
     > '*  [(bundled Grassmann/Stiefel manifold implementation)] *')
 1035 FORMAT(10x,
     > '*         [ NorthWest Chemistry implementation ]         *')
 1040 FORMAT(10X,
     > '*                version #1.10   01/31/03                *')
 1041 FORMAT(10X,
     > '*  A pseudopotential plane-wave band structure program   *')
 1042 FORMAT(10X,
     > '*  with Brillouin zone sampling for optimizing crystals, *')
 1043 FORMAT(10X,
     > '*  slabs, and polymers.  Developed by Eric J. Bylaska,   *')
 1044 FORMAT(10X,
     > '*  Edoardo Apra, and Patrick Nichols.                    *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1111 FORMAT(/' number of processors used:',I16)
 1112 FORMAT( ' parallel mapping         :         1d slab')
 1113 FORMAT( ' parallel mapping         :     2d  hilbert')
 1115 FORMAT(/' options:')
 1117 FORMAT( ' processor grid           :',I4,' x',I4,' x',I4)
 1118 FORMAT( ' parallel mapping         :       2d hcurve')
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1140 FORMAT(/' elements involved in the cluster:')
 1150 FORMAT(5X,I2,': ',A4,'  core charge:',F4.1,'  lmax=',I1)
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
 1170 FORMAT(7(5X,A4,':',I3))
 1180 FORMAT(/' initial position of ions:')
 1190 FORMAT(5X, I4, A5  ,' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1191 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - fixed')
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1211 FORMAT(5X,'   number of constraints = ', I6,' ( DOF = ',I6,' )' )
 1220 FORMAT(/' number of electrons: spin up=',F8.2,
     >                          '  spin down=',F8.2,A)
 1221 FORMAT( ' number of orbitals:  spin up=',I8,
     >                          '  spin down=',I8,A)
 1230 FORMAT(/' supercell:')
 1231 FORMAT(5x,' volume : ',F10.1)
 1232 FORMAT(/5x,' lattice:    a=',f8.3,'    b=',f8.3,'     c=',f8.3,
     >       /5x,'         alpha=',f8.3,' beta=',f8.3,' gamma=',f8.3)
 1241 FORMAT(5x,' lattice:    a1=<',3f8.3,' >')
 1242 FORMAT(5x,'             a2=<',3f8.3,' >')
 1243 FORMAT(5x,'             a3=<',3f8.3,' >')
 1244 FORMAT(5x,' reciprocal: b1=<',3f8.3,' >')
 1245 FORMAT(5x,'             b2=<',3f8.3,' >')
 1246 FORMAT(5x,'             b3=<',3f8.3,' >')

 1249 FORMAT(/' computational grids:')
 1250 FORMAT(5X,' density     cutoff=',F7.3,'  fft=',I4,'x',I4,'x',I4,
     &       '( ',I8,' waves ',I8,' per task)')
 1251 FORMAT(5X,' wavefnc ',I3,' cutoff=',F7.3,
     &        '  fft=',I4,'x',I4,'x',I4,
     &       '( ',I8,' waves ',I8,' per task)')
 1252 FORMAT(5x,' wavefnc     cutoff=',F7.3,
     >       ' wavefunction grids not printed - ',
     >       'number of k-points is very large')
 1253 FORMAT(5X,' smooth wavefnc cutoff:  A=',F7.3,' sigma = ',F7.3)

  
 1255 FORMAT(/' brillouin zone:')
 1256 FORMAT(5x,' number of zone points:',I6,
     >          ' (reduced from ',I0,' zone points without symmetry)')
 1257 FORMAT(5x,' weight=',f8.3,'  ks=<',3f8.3,' >, k=<',3f8.3,'>')
 1258 FORMAT(5x,' number of k-points is very large')

 1260 FORMAT(5X,' Ewald summation: cut radius=',F8.2,'  and',I3)
 1261 FORMAT(5X,'                   madelung=',f14.8)
 
 1270 FORMAT(/' technical parameters:')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'fictitious mass=',F10.1)
 1281 FORMAT(5X, ' maximum iterations =',I10,
     >           ' ( ',I4,' inner ',I6,' outer )')
 1290 FORMAT(5X, ' tolerance=',E8.3,' (energy)',E12.3,
     &        ' (density)')
 1291 FORMAT(/' Kohn-Sham scf parameters:')
 1292 FORMAT(5X, ' Kohn-Sham algorithm  = ',A)
 1293 FORMAT(5X, ' scf algorithm        = ',A)
 1294 FORMAT(5X, ' scf mixing parameter =',F7.4)
 1295 FORMAT(5X, ' Kohn-Sham iterations = ',I3,
     >           ' (',I3,' outer)')
 1296 FORMAT(5X, ' scf mixing type      = ',A)
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
 1600 FORMAT(/' Total BAND energy   :',E19.10)

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

      band_minimizer = value
      return
      END

