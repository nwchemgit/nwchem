*
* $Id$
*

***********************************************************************
*                      cgmontecarlov1                                 *
*                                                                     *
*     This is a developing NVT and NPT Monte-Carlo code for NWCHEM    *
*                                                                     *
***********************************************************************

      logical function cgmontecarlov1(rtdb,flag)
      implicit none
      integer rtdb
      integer flag

#include "global.fh"
#include "errquit.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "inp.fh"
#include "util.fh"
#include "geom.fh"
#include "stdio.fh"

#include "nwpw_timing.fh"

*     *** local variables and parameters ****
      double precision kb
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
      integer ngrid(3)
      real*8  a,b,c,alpha,beta,gamma,mc_beta
    

*     ***** energy variables ****
      integer ms,ispin
      real*8  E(40),en(2)
      real*8  dipole(3),cdipole(3)
      real*8  stress(3,3),lstress(6)

*     **** gradient variables ****
      integer fion(2)

*     **** error variables ****
      logical value,newpsi
      integer ierr

*     **** local variables ****
      logical verlet,mulliken,SA,calc_pressure,change_volume
      logical oprint,lprint,hprint,qmmm,found,notfirststep,found_bak
      real*8  gx,gy,gz,cx,cy,cz,ddx,c0,c1,c2,dA,ddv,v0,v1,v2
      real*8  vgx,vgy,vgz,vcx,vcy,vcz,ekg,eki0,eki1,dt
      real*8  sum,w,emotion_time_shift,wa,wr,wa_r,wr_r,wa_v,wr_v
      real*8  EV,EV0,pi,esum1,esum2,eave,evar,E25,E26
      real*8  mc_unita1(3,3),mc_unita2(3,3)
      real*8  mc_dr,mc_dV,mc_atom_direction(3)
      real*8  ratio,ratio_v,ratio_r,aratio,omega1,omega2
      real*8  omega_all,omega_all2,omega_ave,omega_var
      real*8  icharge,cv
      real*8  sa_alpha(2),sa_decay(2),tt1,tt2,ssr,r,p1,p2,pressure
      real*8 Te_new,Tr_new,Te_init,Tr_init,mc_temperature,mc_pressure
      integer i,k,ia,nion,minimizer,mapping,icount
      integer nbq,mapping1d,it_in,it_out,seed
      integer mc_algorithm
      integer mc_napply,mc_ngroups,mc_group_size
      integer mc_group_start(2),mc_group_end(2),mc_group(2)
      integer vcount,geom1,frac1(2),nion1

      character*50 filename
      character*255 full_filename,full_bak


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

      real*8   control_tole,control_tolc,control_tolr,ion_rion,ion_vion
      external control_tole,control_tolc,control_tolr,ion_rion,ion_vion
      real*8   control_time_step,control_fake_mass,control_bo_fake_mass
      external control_time_step,control_fake_mass,control_bo_fake_mass
      logical  control_read,control_move,ion_init,control_out_of_time
      external control_read,control_move,ion_init,control_out_of_time
      logical  control_translation,control_rotation,control_parallel_io
      external control_translation,control_rotation,control_parallel_io
      integer  pack_nwave_all
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      integer  pspw_charge_nion
      external pack_nwave_all
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      external pspw_charge_nion

      character*12 control_boundry
      external     control_boundry
      character*50 control_cell_name
      external     control_cell_name

      logical      psp_semicore
      real*8       psp_rcore,psp_ncore,psp_rlocal
      external     psp_semicore
      external     psp_rcore,psp_ncore,psp_rlocal

      logical  psi_initialize, psi_finalize
      integer  psi_ispin,psi_ne,electron_count,linesearch_count
      external psi_initialize, psi_finalize
      external psi_ispin,psi_ne,electron_count,linesearch_count
      real*8   cgsd_energy,cgsd_noit_energy
      external cgsd_energy,cgsd_noit_energy
      logical  control_Mulliken,control_DOS
      external control_Mulliken,control_DOS
      real*8   ion_TotalCharge
      external ion_TotalCharge
      logical  control_check_charge_multiplicity
      external control_check_charge_multiplicity
      logical   pspw_charge_found,ion_q_FixIon,ion_q_xyzFixIon
      external  pspw_charge_found,ion_q_FixIon,ion_q_xyzFixIon
      integer  control_minimizer,control_scf_algorithm
      external control_minimizer,control_scf_algorithm
      integer  control_ks_algorithm
      external control_ks_algorithm
      real*8   control_ks_alpha,control_kerker_g0
      external control_ks_alpha,control_kerker_g0
      logical  control_print,control_balance
      external control_print,control_balance
      integer  control_mapping,control_mapping1d,control_np_orbital
      external control_mapping,control_mapping1d,control_np_orbital
      integer  control_ks_maxit_orb,control_ks_maxit_orbs
      external control_ks_maxit_orb,control_ks_maxit_orbs
      character*14 ion_q_xyzFixIon_label
      external     ion_q_xyzFixIon_label

      logical  pspw_bqext,control_Fei
      external pspw_bqext,control_Fei

      integer  control_bo_steps_in,control_bo_steps_out
      integer  control_bo_algorithm,control_mc_seed,control_mc_algorithm
      real*8   control_mc_step_size,control_rti,ion_ke,ion_Temperature
      external control_bo_steps_in,control_bo_steps_out 
      external control_bo_algorithm,control_mc_seed,control_mc_algorithm
      external control_mc_step_size,control_rti,ion_ke,ion_Temperature
      real*8   control_mc_aratio,control_mc_Temperature,control_mc_ddx
      external control_mc_aratio,control_mc_Temperature,control_mc_ddx
      real*8   control_mc_ddv,control_mc_volume_step
      external control_mc_ddv,control_mc_volume_step
      real*8   control_mc_pressure
      external control_mc_pressure
      logical  control_SA,control_pressure,ion_makehmass2
      external control_SA,control_pressure,ion_makehmass2
      real*8   control_SA_decay
      external control_SA_decay
      integer  ion_nconstraints,ion_ndof,ion_rion_ptr
      external ion_nconstraints,ion_ndof,ion_rion_ptr
      logical  control_periodic_dipole
      external control_periodic_dipole


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

      value = control_read(11,rtdb)
      if (.not. value) 
     > call errquit('error reading control',0, INPUT_ERR)

      call Parallel2d_Init(control_np_orbital())
      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)

      oprint= ((taskid.eq.MASTER).and.control_print(print_medium))
      lprint= ((taskid.eq.MASTER).and.control_print(print_low))
      hprint= ((taskid.eq.MASTER).and.control_print(print_high))


      if (oprint) call current_second(cpu1)

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
         write(luout,1010)
         write(luout,1000)
         call nwpw_message(1)
         write(luout,1110)
      end if
      
      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      minimizer = control_minimizer()
      mapping   = control_mapping()

*     **** initialize psi_data ****
      call psi_data_init(100)


*     **** initialize a bunch of stuff ***
      call cgmontecarlo_start(rtdb,newpsi)

*     **** initialize DFT+U ****
      call psp_U_init()

*     **** initialize metadynamics ****
      call meta_initialize(rtdb)


*     **** initialize pressure ****
      calc_pressure = control_pressure().and.(control_version().eq.3)
      pressure      = 0.0d0
      p1            = 0.0d0
      p2            = 0.0d0

*     **** initialize linesearching ****
      call linesearch_init()


*     ******************************
*     **** scaling ion velocity ****
*     ******************************
      call ion_init_ke(ekg,eki0,eki1)


*     **** Initialize simulated annealing ****
      SA       = .false.
      Tr_init  = 0.0d0
      sa_alpha(2) = 1.0d0
      if (control_SA()) then
         dt = 1.0d0
         SA          = .false.
         sa_decay(2) = control_SA_decay(2)
         sa_alpha(2) = dexp( -(dt/control_SA_decay(2)) )
      end if



*     **** Initialize Monte-Carlo Parameters ****
      mc_dr  = control_mc_step_size()
      mc_dV  = control_mc_volume_step()
      mc_algorithm  = control_mc_algorithm()
      call control_mc_atom_direction(mc_atom_direction)
      call control_mc_ngroups(mc_napply,mc_ngroups,mc_group_size)
      if (mc_ngroups.lt.1) then
      
         value = BA_alloc_get(mt_int,1,
     >                    'mc_group_start',
     >                     mc_group_start(2),
     >                     mc_group_start(1))
         value = value.and.BA_alloc_get(mt_int,1,
     >                    'mc_group_end',
     >                     mc_group_end(2),
     >                     mc_group_end(1))
         value = value.and.BA_alloc_get(mt_int,1,
     >                    'mc_group',
     >                     mc_group(2),
     >                     mc_group(1))
         if (.not.value) 
     >     call errquit('cgmontecarlov1:out of heap memory',0,MA_ERR)
      else
         value = BA_alloc_get(mt_int,mc_ngroups,
     >                    'mc_group_start',
     >                     mc_group_start(2),
     >                     mc_group_start(1))
         value = value.and.BA_alloc_get(mt_int,mc_ngroups,
     >                    'mc_group_end',
     >                     mc_group_end(2),
     >                     mc_group_end(1))
         value = value.and.BA_alloc_get(mt_int,mc_group_size,
     >                    'mc_group',
     >                     mc_group(2),
     >                     mc_group(1))
         if (.not.value) 
     >      call errquit('cgmontecarlov1:out of heap memory',1,MA_ERR)
         call control_mc_groups(int_mb(mc_group_start(1)),
     >                          int_mb(mc_group_end(1)),
     >                          int_mb(mc_group(1)))
      end if


      seed           = control_mc_seed()
      mc_temperature = control_mc_Temperature()
      mc_pressure    = control_mc_pressure()
      mc_beta = util_random(control_mc_seed()) !**seed the random number generator
      mc_beta = 1.0d0/(kb*mc_temperature)

      aratio      = control_mc_aratio()

      ddx = control_mc_ddx()
      dA  = 1.0d0/(aratio*aratio-aratio)
      c0  = 1.0d0-ddx
      c1  = (aratio*aratio*(2*ddx) - (ddx))*dA
      c2  = (-aratio      *(2*ddx) + (ddx))*dA

      if (mc_algorithm.eq.2) then
         ddv = control_mc_ddv()
         v0  = 1.0d0-ddv
         v1  = (aratio*aratio*(2*ddv) - (ddv))*dA
         v2  = (-aratio      *(2*ddv) + (ddv))*dA
      end if




*                |**************************|
******************   summary of input data  **********************
*                |**************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      call center_v_geom(vcx,vcy,vcz)
      call center_v_mass(vgx,vgy,vgz)
      mulliken = control_Mulliken()


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
         write(luout,1115)
         if (control_parallel_io()) then
           write(luout,1123)
         else
           write(luout,1124)
         end if

         write(luout,1121) control_boundry(),control_version()
         if (psi_ispin().eq.1) write(luout,1130) "restricted"
         if (psi_ispin().eq.2) write(luout,1130) "unrestricted"
         !if (qmmm) write(luout,1122)
         IF (control_gga().eq.-1) THEN
            write(luout,1131) 'off'
         ELSE IF (control_gga().eq.0) THEN
            write(luout,1131) 'LDA (Vosko et al) parameterization'
         ELSE IF (control_gga().eq.10) THEN
            write(luout,1131) 
     >      'PBE96 (White and Bird) parameterization'
         ELSE IF (control_gga().eq.11) THEN
            write(luout,1131) 
     >      'BLYP (White and Bird) parameterization'
         ELSE IF (control_gga().eq.12) THEN
            write(luout,1131) 
     >      'revPBE (White and Bird) parameterization'
         ELSE IF (control_gga().eq.13) THEN
            write(luout,1131) 
     >      'PBEsol (White and Bird) parameterization'

         ELSE IF (control_gga().eq.110) THEN
            write(luout,1131) 
     >      'PBE0 (White and Bird) parameterization'
         ELSE IF (control_gga().eq.112) THEN
            write(luout,1131) 
     >      'revPBE0 (White and Bird) parameterization'
         ELSE IF (control_gga().eq.113) THEN
            write(luout,1131) 
     >      'BNL (White and Bird) parameterization'
         ELSE IF (control_gga().eq.114) THEN
            write(luout,1131) 
     >      'HSE (White and Bird) parameterization'
         ELSE IF (control_gga().eq.115) THEN
            write(luout,1131) 
     >      'B3LYP (White and Bird) parameterization'

         ELSE IF (control_gga().eq.200) THEN
            write(luout,1131) 'Hartree-Fock'
         ELSE
            write(luout,1131) 'unknown parameterization'
            call errquit('bad exchange_correlation',0, INPUT_ERR)
         END IF
         call pspw_print_SIC(luout)
         call pspw_print_HFX(luout)
         if (ion_makehmass2()) write(luout,1135)
         write(luout,1140)
         do ia = 1,ion_nkatm()
            call psp_print(ia)
         end do

         icharge = -(psi_ne(1)+psi_ne(psi_ispin()))
         en(1)     = psi_ne(1)
         en(psi_ispin()) = psi_ne(psi_ispin())

         icharge = icharge + ion_TotalCharge()
         write(luout,1159) icharge
         write(luout,1160)
         write(luout,1170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
         !if (hprint) then
         write(luout,1180)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0
           else if (ion_q_xyzFixIon(I)) then
           write(luout,1194) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_q_xyzFixIon_label(I)
           else
           write(luout,1190) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz

         call pspw_charge_Print(luout)

         write(luout,1181)
         write(luout,1192) (I,ion_aname(I),
     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
         write(luout,1200) vcx,vcy,vcz
         write(luout,1210) vgx,vgy,vgz
         write(luout,1211) ion_nconstraints(),ion_ndof()

         !end if

         !write(6,1220) psi_ne(1),psi_ne(psi_ispin()),' (Fourier space)'
c         if (fractional) then
c          write(6,1219) en(1),en(psi_ispin()),' (   fractional)'
c          write(6,1221) psi_ne(1),psi_ne(psi_ispin()),' (Fourier space)'
c         else
          write(luout,1220) psi_ne(1),psi_ne(psi_ispin()),
     >                      ' (Fourier space)'
          write(luout,1221) psi_ne(1),psi_ne(psi_ispin()),
     >                      ' (Fourier space)'
c         end if
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
         if (control_version().eq.3) then
         write(luout,1260) ewald_rcut(),ewald_ncut()
         write(luout,1261) ewald_mandelung()
         end if

         write(luout,1270)
         write(luout,1280) control_time_step(),control_fake_mass()
         write(luout,1290) control_tole(),control_tolc()
         write(luout,1281) control_it_in()*control_it_out(),
     >                 control_it_in(),control_it_out()
         if ((minimizer.eq.5).or.(minimizer.eq.8)) then
           write(luout,1291)
           if (control_ks_algorithm().eq.0) 
     >       write(luout,1292) "conjugate gradient"
           if (control_ks_algorithm().eq.1) 
     >       write(luout,1292) "rmm-diis"
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
           write(luout,1294) control_ks_alpha()
           write(luout,1301) control_kerker_g0()
         end if
         write(luout,1310)
         if (.not.control_translation()) write(luout,1271)
         if (.not.control_rotation())    write(luout,1272)
         if (mc_algorithm.eq.0) then
            write(luout,1311) 'NVE Microcanonical Ensemble - ',
     >                        'HyperVolume Monte Carlo Algorithm'
         else if (mc_algorithm.eq.1) then
            write(luout,1311) 'NVT Canonical Ensemble - ',
     >                        'Metropolis Algorithm'
         else if (mc_algorithm.eq.2) then
            write(luout,1311) 'NPT Canonical Ensemble - ',
     >                        'Metropolis Algorithm'
         end if
         
         write(luout,1320) mc_dr,
     >    control_bo_steps_in()*control_bo_steps_out(),
     >    control_bo_steps_in(),control_bo_steps_out() 

         write(luout,1322) mc_temperature,mc_beta,seed
         if (mc_algorithm.eq.2) 
     >      write(luout,1326) mc_pressure,mc_pressure*autoatm,
     >                        ddv,v0,v1,v2
         write(luout,1323) aratio,ddx,c0,c1,c2
         write(luout,1324) mc_napply,mc_atom_direction
         write(luout,1325) mc_ngroups

         write(luout,1340) control_rti()
         call ion_scaling_atoms_print(luout)
         write(luout,1222) eki0,ekg
         write(luout,1223) eki1
         write(luout,1224) (eki1-eki0)

c         write(luout,1394)
       
        if (calc_pressure) write(luout,1393)
         if (control_SA()) then
           write(luout,1396) sa_decay(2)
         end if
         if (mulliken) write(luout,1399)
         write(luout,1300)
         call util_flush(luout)
      end if

*                |***************************|
******************  simple Monte-Carlo loop  **********************
*                |***************************|
      if (taskid.eq.MASTER) call current_second(cpu2)
      if ((taskid.eq.MASTER).and.(.not.calc_pressure))  then
         if (mc_algorithm.lt.2) then
            call nwpw_message(12)
         else
            call nwpw_message(13)
         end if
      end if

      if ((taskid.eq.MASTER).and.(calc_pressure)) 
     >   call nwpw_message(11)


*     **** write initial position to xyz data ****
      call xyz_init()          ! unit=18
      call MOTION_init(rtdb)   ! unit=19


*     *** fei io ****
      call fei_init(rtdb)

*     **** ecce print ****
      call ecce_print_module_entry('task pspw Monte-Carlo')
      !call ecce_print_module_entry('driver')
      call movecs_ecce_print_off()



*     ************************************
*     **** open up other MOTION files ****
*     ************************************
*     **** open EMOTION file ****
      E25 = 0.0d0
      E26 = 0.0d0
      omega_all  = 0.0d0
      omega_all2 = 0.0d0
      wa  = 0.0d0
      wr  = 0.0d0
      wa_r  = 0.0d0
      wr_r  = 0.0d0
      wa_v  = 0.0d0
      wr_v  = 0.0d0
      if (.not.btdb_cget(rtdb,'nwpw:emotion_filename',1,filename))
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
            k=index(full_filename,' ')
            write(*,*) "   Copying ",full_bak(1:i),
     >                 " to ",full_filename(1:k)
            write(*,*)
            call util_file_copy(full_bak,full_filename)
         end if

         emotion_time_shift = 0.0d0
         inquire(file=full_filename,exist=found)
         if (found) then
           open(unit=31,file=full_filename,form='formatted',
     >          status='old')
           do while (found)
           if (mc_algorithm.lt.2) then
              read(31,*,end=100) emotion_time_shift,
     >                        EV,EV0,eave,evar,ratio,wa,wr,mc_dr
           else
              read(31,*,end=100) emotion_time_shift,
     >                EV,EV0,eave,evar,ratio,wa,wr,
     >                mc_dr,ratio_r,wa_r,wr_r,
     >                mc_dV,ratio_v,wa_v,wr_v,
     >                omega1,omega_ave,omega_var
              omega_all  = omega_all  + omega1
              omega_all2 = omega_all2 + omega1*omega1
           end if

           E25 = E25 + EV   !*** take care of running sums ***
           E26 = E26 + EV*EV
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
      call Parallel_Brdcst_value(MASTER,mc_dr)
      call Parallel_Brdcst_value(MASTER,E25)
      call Parallel_Brdcst_value(MASTER,E26)
      call Parallel_Brdcst_value(MASTER,wa)
      call Parallel_Brdcst_value(MASTER,wr)
      call Parallel_Brdcst_value(MASTER,wa_r)
      call Parallel_Brdcst_value(MASTER,wr_r)

      if (mc_algorithm.eq.2) then
         call Parallel_Brdcst_value(MASTER,mc_dV)
         call Parallel_Brdcst_value(MASTER,omega_all)
         call Parallel_Brdcst_value(MASTER,omega_all2)
         call Parallel_Brdcst_value(MASTER,wa_v)
         call Parallel_Brdcst_value(MASTER,wr_v)
      end if

*     **** open EIGMOTION file ****
      if (mulliken) then
        if (.not.btdb_cget(rtdb,'nwpw:eigmotion_filename',1,filename))
     >    call util_file_prefix('eigmotion',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)
      if (taskid.eq.MASTER)
     >   open(unit=32,file=full_filename,form='formatted')
      end if

      call xyz_write()

*     ****  allocate fion ***
      nion = ion_nion()
      if (pspw_charge_found().and.
     >    (.not.pspw_bqext())) nion = nion + pspw_charge_nion()

      value = BA_push_get(mt_dbl,(3*nion),
     >                       'fion',fion(2),fion(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


*     *******************************
*     ***** do a simple mc loop *****
*     *******************************
      notfirststep = .false.



      call control_reduce_print()
      EV = cgsd_energy(newpsi)
      if (control_Fei())
     >   call cgsd_energy_gradient_md(dbl_mb(fion(1)))

      call xyz_write()

      it_out = control_bo_steps_out()
      it_in  = control_bo_steps_in()
      icount = 0
      if (it_out.lt.1) goto 102


      omega1 = lattice_omega()
      vcount = 0
      change_volume = .false.
  101 continue
         icount = icount + 1

c       **** inner loop ****
        do i=1,it_in

           !*** compute the trial energy ****
           if (mc_algorithm.eq.2) then
              if (taskid.eq.MASTER) alpha = util_random(0)
              call Parallel_Brdcst_value(MASTER,alpha)
              change_volume = alpha.le.(1.0d0/dble(ion_nion()))
           end if

*          **** change volume ****
           if (change_volume) then
              vcount = vcount + 1
              omega1 = lattice_omega()

              !*** rion --> frac1  ***
              call lattice_r1_to_frac(nion,
     >                                dbl_mb(ion_rion_ptr()),
     >                                dbl_mb(frac1(1)))
              call control_get_unita(mc_unita1)
              call cgmontecarlo_end(rtdb,flag)

*             **** volume change to unita1 --> unita2 ****
              if (taskid.eq.MASTER) alpha = (1.0d0*util_random(0)-0.5d0)
              call Parallel_Brdcst_value(MASTER,alpha)
              alpha = dexp(alpha*mc_dV/3.0d0)
              call dcopy(9,mc_unita1,1,mc_unita2,1)
              call dscal(9,alpha,mc_unita2,1)
              call Parallel_Brdcst_values(MASTER,9,mc_unita2)


*             **** change unitcell and reinitialize ****
              call control_set_unita(mc_unita2)
              call control_set_unita_frozen(mc_unita2)
              call cgmontecarlo_start(rtdb,newpsi)

              !*** frac1 -- > rion ***
              call lattice_frac_to_r1(nion,
     >                                dbl_mb(frac1(1)),
     >                                dbl_mb(ion_rion_ptr()))
              omega2 = lattice_omega()


*          **** change atom positions ****
           else
              call ion_MC_step(mc_dr,
     >                         mc_atom_direction,
     >                         mc_napply,
     >                         mc_ngroups,
     >                         int_mb(mc_group_start(1)),
     >                         int_mb(mc_group_end(1)),
     >                         int_mb(mc_group(1)))
           end if
           EV0 = cgsd_energy(.false.)

           if (change_volume) then
              w = dexp(-mc_beta*(EV0-EV) 
     >                 -mc_beta*mc_pressure*(omega2-omega1) 
     >                 +(ion_nion()+1)*dlog(omega2/omega1))
           else
              w = dexp(-mc_beta*(EV0-EV))
           end if

           if (taskid.eq.MASTER) alpha = util_random(0)
           call Parallel_Brdcst_value(MASTER,alpha)
        
           !*** accept the sequence ***
           if (alpha.lt.w) then
              EV = EV0
              if (change_volume) then
                 wa_v = wa_v + 1.0d0
                 omega1 = omega2
              else
                 wa_r = wa_r + 1.0d0
              end if
              wa = wa + 1.0d0
              if (control_Fei())
     >           call cgsd_energy_gradient_md(dbl_mb(fion(1)))

           !*** reject the sequence ***
           else
              if (change_volume) then
                 call cgmontecarlo_end(rtdb,flag)
                 call Parallel_Brdcst_values(MASTER,9,mc_unita1)
                 call control_set_unita(mc_unita1)
                 call control_set_unita_frozen(mc_unita1)
                 call cgmontecarlo_start(rtdb,newpsi)
                 call lattice_frac_to_r1(nion,
     >                                   dbl_mb(frac1(1)),
     >                                   dbl_mb(ion_rion_ptr()))
                 wr_v = wr_v + 1.0d0
              else
                 call ion_MC_reject_step()
                 wr_r = wr_r + 1.0d0
              end if
              wr = wr + 1.0d0
           end if

           E25 = E25 + EV
           E26 = E26 + EV*EV
           eave = E25/(wa+wr)
           evar = E26/(wa+wr)
           evar = evar - eave*eave
           ratio_r = wa_r/(wa_r+wr_r)
           ratio   = wa/(wa+wr)

           if (mc_algorithm.eq.2) then
              omega_all  = omega_all  + omega1
              omega_all2 = omega_all2 + omega1*omega1
              omega_ave  = omega_all/(wa+wr)
              omega_var  = omega_all2/(wa+wr)
              omega_var  = omega_var - omega_ave*omega_ave
              if ((wa_v+wr_v).gt.1.0d-6) then
                 ratio_v = wa_v/(wa_v+wr_v)
              else
                 ratio_v = 0.0d0
              end if
           end if

           if (oprint) then
             if (mc_algorithm.lt.2) then
                write(luout,1350) (icount-1)*it_in+i,EV,EV0,eave,
     >                              ratio,mc_dr
                write(31,1360) 
     >                dble((icount-1)*it_in+i)+emotion_time_shift,
     >                EV,EV0,eave,evar,ratio,wa,wr,mc_dr
             else
                write(luout,1351) (icount-1)*it_in+i,
     >                        EV,EV0,eave,
     >                        omega1,omega_ave,
     >                        ratio,mc_dr,mc_dv
c     >                        ratio_v,wa,wr,wa_v,wr_v
                write(31,1360) 
     >                dble((icount-1)*it_in+i)+emotion_time_shift,
     >                EV,EV0,eave,evar,ratio,wa,wr,
     >                mc_dr,ratio_r,wa_r,wr_r,
     >                mc_dv,ratio_v,wa_v,wr_v,
     >                omega1,omega_ave,omega_var
             end if
             call util_flush(luout)
             call util_flush(31)
           end if
           call fei_output(EV,dbl_mb(fion(1)))
           call xyz_write()
           call MOTION_write(dble((icount-1)*it_in+i)
     >                       +emotion_time_shift)
           call CIF_write(rtdb) !*** write out CIF file ***

        end do

*       **** adjust MC steplength ****
        mc_dr = mc_dr*(c0 + c1*ratio_r + c2*ratio_r*ratio_r)
        !if ((ratio/aratio).gt.1.05d0) then
        !   dr = dr*(1.05d0)
        !else if ((ratio/aratio).lt.0.95d0) then
        !   dr = dr*(0.95d0)
        !else
        !   dr = dr*(ratio/aratio)
        !end if
        !if (ratio.lt.aratio) dr = dr*0.95d0
        !if (ratio.gt.aratio) dr = dr*1.05d0
        call Parallel_Brdcst_value(MASTER,mc_dr)

        if (mc_algorithm.eq.2) then
           if (vcount.ge.it_in) then
              mc_dV = mc_dV*(v0 + v1*ratio_v + v2*ratio_v*ratio_v)
              call Parallel_Brdcst_value(MASTER,mc_dV)
              vcount = 0
           end if
        end if

*        **** update thermostats using SA decay ****
         if (SA) then
           tt2 = icount*it_in*dt/sa_decay(2)
           Tr_new = Tr_init*dexp(-tt2)
         end if

*        **** exit early ****
         if (control_out_of_time()) then
            if (oprint)
     >       write(luout,*) ' *** out of time. iteration terminated'
            go to 102
         end if
      if (icount.lt.it_out) go to 101
      if (oprint) write(luout,*) 
     > '*** arrived at the Maximum iteration.   terminated.'

*::::::::::::::::::::  end of iteration loop  :::::::::::::::::::::::::
      call control_up_print()

  102 continue
  

*     **** close xyz and MOTION files ****
      call xyz_end()
      call MOTION_end()
      if (taskid.eq.MASTER) then
         close(unit=31)
         if (mulliken) close(unit=32)

*        **** remove EMOTION backup file ***
         call util_file_name_noprefix('EMOTION99-bak',.false.,
     >                                .false.,
     >                                full_bak)
         call util_file_unlink(full_bak)

      end if

*     *** close fei io ****
      call fei_end()




*                |***************************|
******************     simple MD loop        **********************
*                |***************************|


*                |***************************|
****************** report summary of results **********************
*                |***************************|

      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      call center_v_geom(vcx,vcy,vcz)
      call center_v_mass(vgx,vgy,vgz)

      if (oprint) then
         call print_elapsed_time(icount*it_in*dt)
         write(luout,1300)
         write(luout,1410)
         write(luout,1420)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           else if (ion_q_xyzFixIon(I)) then
           write(luout,1194) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_q_xyzFixIon_label(I)
           else
           write(luout,1190) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz
         write(luout,1211) ion_nconstraints(),ion_ndof()
      end if

      call control_reduce_print()
      EV0 = cgsd_energy(.false.)
      call control_up_print()
      EV  = cgsd_noit_energy()

      if (oprint) then


*       **** write out Temperatures ****
         !cv = (evar)/(kb*ion_Temperature()**2)
         cv = (evar)/(kb*mc_temperature**2)
         cv = cv/dble(ion_nion())
         write(luout,1492) eave
         write(luout,1493) evar
         write(luout,1494) cv

      end if
      call cgsd_energy_gradient(dbl_mb(fion(1)))

*     **** calculate excited state orbitals ****
      call cgsd_excited()

*     **** calculate oep eigenvalues ****
      call cgsd_oep_eigenvalues()

*     **** extra energy output for QA test ****
      if (lprint) write(6,1600) EV

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


*     *************************************************************
*     **** output energy, dipole, and gradient to rtdb for use ****
*     **** by task_energy and task_gradient                    ****
*     *************************************************************
      if (flag.gt.-1) then
      value = btdb_put(rtdb,'pspw:energy',mt_dbl,1,EV)
      value = value.and.
     >        btdb_put(rtdb,'pspw:dipole',mt_dbl,
     >                 3,dipole)
      value = value.and.
     >        btdb_put(rtdb,'pspw:crystal_dipole',mt_dbl,
     >                 3,cdipole)
      end if
      if (flag.gt.0) then
        value = value.and.
     >        btdb_put(rtdb,'pspw:gradient',mt_dbl,
     >                 3*nion,dbl_mb(fion(1)))
      end if
      if (flag.eq.3) then
        value = value.and.
     >        btdb_put(rtdb,'pspw:stress',mt_dbl,
     >                 9,stress)
        value = value.and.
     >        btdb_put(rtdb,'pspw:lstress',mt_dbl,
     >                 6,lstress)
      end if
      value = value.and.BA_pop_stack(fion(2))
      if (.not. value) 
     > call errquit('cgmontecarlov1: error writing rtdb',0,
     >       RTDB_ERR)

      if (taskid.eq.MASTER) call current_second(cpu3)

*                |***************************|
******************         Epilogue          **********************
*                |***************************|

*     **** calculate Mulliken Populations ***
      if (control_Mulliken()) call psi_Mulliken(rtdb)


*     **** deallocate memory and data structures ****
      call cgmontecarlo_end(rtdb,flag)


*     **** deallocate heap memory ****
      value =           BA_free_heap(mc_group_start(2))
      value = value.and.BA_free_heap(mc_group_end(2))
      value = value.and.BA_free_heap(mc_group(2))
      if (.not. value) 
     >  call errquit('cgmontecarlov1:freeing heap',3,MA_ERR)

      call meta_finalize(rtdb)
      call psp_U_end()
      call psi_data_end()


*     **** do anaylysis on MOTION files ****
      call cpmd_properties(rtdb)



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
         write(6,1801)
         write(6,1802)
         write(6,1803) T1
         write(6,1804) T2
         write(6,1805) T3
         write(6,1806) T4
         write(6,1807) AV,electron_count(),linesearch_count()

         call nwpw_timing_print_final(oprint,electron_count())

         write(6,*)
         CALL nwpw_MESSAGE(4)
      end if 



      call Parallel_Finalize()
      cgmontecarlov1 = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*               NWPW PSPW Calculation              *')
 1030 FORMAT(10X,'*     [  NVT and NPT Monte-Carlo simulation   ]    *')
 1035 FORMAT(10x,'*     [  NorthWest Chemistry implementation   ]    *')
 1040 FORMAT(10X,'*            version #1.00   07/18/06              *')
 1041 FORMAT(10X,'*    This code was developed by Eric J. Bylaska.   *')
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
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1122 FORMAT(5X,' qmmm simulation')
 1123 FORMAT( ' parallel io              :        on')
 1124 FORMAT( ' parallel io              :       off')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1132 FORMAT(5X,' using fractional occupation')
 1135 FORMAT(/' The masses of QM H atoms converted to 2.0 amu.',
     >       /' To turn off this default',
     >        ' set nwpw:makehmass2 .false.')
 1140 FORMAT(/' elements involved in the cluster:')
 1150 FORMAT(5X,I2,': ',A4,'  core charge:',F4.1,'  lmax=',I3)
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
 1180 FORMAT(/' initial position of ions (au):')
 1181 FORMAT(/' initial velocity of ions after scaling (au):')

 1190 FORMAT(5X, I4, A5  ,' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1191 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - fixed')
 1192 FORMAT(5X, I4, A5  ,' (',3F11.5,' )')
 1193 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - z fixed')
 1194 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,A)

 1200 FORMAT(5X,'    G.C. ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1211 FORMAT(5X,'   number of constraints = ', I6,' ( DOF = ',I6,' )' )
 1219 FORMAT(/' number of active electrons: spin up=',F6.2,
     >        '  down=',F6.2,A)
 1220 FORMAT(/' number of active electrons: spin up=',I6,  
     >        '  down=',I6,A)
 1221 FORMAT( ' number of active orbitals : spin up=',I6,  
     >        '  down=',I6,A)

 1222 format(5x,' initial kinetic energy=',e12.5,' (ion)',2x,
     >                                     e12.5,' (c.o.m.)')
 1223 format(5x,' after scaling=         ',e12.5,' (ion)')
 1224 format(5x,' increased energy=      ',e12.5,' (ion)')
 1226 format(/' final kinetic energy= ',   e12.5,' (ion)',2x,
     >                                     e12.5,' (c.o.m.)')

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
  
 1260 FORMAT(5X,' Ewald summation: cut radius=',F8.2,'  and',I3)
 1261 FORMAT(5X,'                   madelung=',f14.8)
 1270 FORMAT(/' technical parameters for minimizer:')
 1271 FORMAT(5x, ' translation constrained')
 1272 FORMAT(5x, ' rotation constrained')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'fictitious mass=',F10.1)
 1281 FORMAT(5X, ' maximum iterations =',I10, 
     >           ' ( ',I4,' inner ',I6,' outer )')
 1290 FORMAT(5X, ' tolerance=',E8.3,' (energy)',E12.3,
     &        ' (density)')
 1291 FORMAT(/' Kohn-Sham scf parameters:')
 1292 FORMAT(5X, ' Kohn-Sham algorithm  = ',A)
 1293 FORMAT(5X, ' SCF algorithm        = ',A)
 1294 FORMAT(5X, ' SCF mixing parameter =',F7.4)
 1295 FORMAT(5X, ' Kohn-Sham iterations = ',I4,
     >           ' (',I4,' outer)')
 1300 FORMAT(//)
 1301 FORMAT(5X, ' Kerker damping       =',F7.4)
 1310 FORMAT(/' Monte-Carlo parameters:')
 1311 FORMAT(5X, ' MC Algorithm= ',A,A)
 1320 FORMAT(5X, ' MC displacement=',F10.2,5X,' iterations=',I10,
     >           ' ( ',I4,' inner ',I6,' outer )')
 1321 FORMAT(5X, ' time step=',F10.2,5X,'fictitious mass=',F10.1,
     >       /5X,' iterations=',I10,
     >           ' ( ',I4,' inner ',I6,' outer )')
 1322 FORMAT(5X, ' temperature='F10.3,5x,' beta=',E10.3,' seed=',I10)
 1323 FORMAT(5X, ' targeted accept.=',F10.3,9X,
     >           ' ddx=',F10.3,' c0=',F10.3,' c1=',F10.3,' c2=',F10.3)
 1324 FORMAT(5X, ' number atom changes = 'I10,
     >           ' change vector = <',3F10.3'>')
 1325 FORMAT(5X, ' number of groups = 'I10)
 1326 FORMAT(5X, ' pressure='E10.3,' (',F10.3,' atm)',
     >           ' ddv=',F10.3,' v0=',F10.3,' v1=',F10.3,' v2=',F10.3)

 1330 FORMAT(5X, ' integration algorithm=   ',A)
 1340 FORMAT(/5X, ' cooling/heatting rate= ',e12.5,' (ion)')
 1350 FORMAT(I8,3E19.10,F16.3,F10.3,3E19.10,2F8.4)
 1351 FORMAT(I8,3E19.10,F10.2,3x,F10.2,4x,F8.3,2x,F8.3,3x,F8.3,9E19.10)
 1360 format(100e19.10)

 1393 format(/' Pressure Output Generated         ')
 1394 format(/' Constant Energy Simulation                     ')
 1395 format(/' Nose-Hoover Simulation - Thermostat Parameters:')
 1396 format(5x, 'SA decay rates  =',e10.3,' (ion)')
 1397 format(5x, 'link = ',I3,
     > ' Te =',f8.2,' Qe =',e10.3,' 2*pi/we=',e10.3,' Ee0=',e10.3)
 1398 format(5x, 'link = ',I3,
     > ' Tr =',f8.2,' Qr =',e10.3,' 2*pi/wr=',e10.3,' Er0=',e10.3)

 1399 format(//' Mulliken Analysis Output Generated            ')
 1400 FORMAT(I3,3E18.8/3X,3E18.8)
 1410 FORMAT(10X,'=============  summary of results  =================')
 1420 FORMAT(/' final position of ions:')
 1421 FORMAT(/' final velocity of ions:')
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
 1600 FORMAT(/' Total PSPW energy   :',E19.10)
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

 9000 if (taskid.eq.MASTER) write(6,9010) ierr
      call Parallel_Finalize()

      cgmontecarlov1 = value
      return
      END



*     ***********************************************
*     *                                             *
*     *             cgmontecarlo_start              *
*     *                                             *
*     ***********************************************

      subroutine cgmontecarlo_start(rtdb,newpsi)
      implicit none
      integer rtdb
      logical newpsi

*     **** local variables ****
      logical value,calc_pressure
      integer minimizer,mapping,mapping1d,ngrid(3),ne(2),ispin

*     **** external functions ****
      logical  control_pressure,control_check_charge_multiplicity
      logical  psp_semicore,ion_init,psi_initialize
      integer  control_version,control_ngrid,control_minimizer
      integer  control_mapping,control_mapping1d,psi_ne,psi_ispin
      external control_pressure,control_check_charge_multiplicity
      external psp_semicore,ion_init,psi_initialize
      external control_version,control_ngrid,control_minimizer
      external control_mapping,control_mapping1d,psi_ne,psi_ispin
      integer  control_gga
      external control_gga



      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      minimizer = control_minimizer()
      mapping   = control_mapping()
      calc_pressure = control_pressure().and.(control_version().eq.3)

*     **** initialize D3dB data structure ****
      call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)

*     ***** Initialize double D3dB data structure ****
      if (control_version().eq.4)
     >   call D3dB_Init(2,2*ngrid(1),2*ngrid(2),2*ngrid(3),mapping)

*     **** initialize lattice data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_init()
      call D3dB_pfft_init()
      call ga_sync()

*     **** read ions ****
      value = ion_init(rtdb)

*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)

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

*     **** initialize QM/MM ****
      call pspw_init_APC(rtdb)
      call pspw_qmmm_init(rtdb)
      call pspw_charge_init(rtdb)

*     **** initialize SIC and HFX ****
      ne(1) = psi_ne(1)
      ne(2) = psi_ne(2)
      call pspw_init_SIC(rtdb,ne)
      call pspw_init_HFX(rtdb,psi_ispin(),ne)

*     **** initialize Meta GGA ****
      call nwpw_meta_gga_init(control_gga())


      if (calc_pressure) then
         call psp_stress_init()
         call psp_stress_readall()
      end if

      return
      end


*     ***********************************************
*     *                                             *
*     *             cgmontecarlo_end                *
*     *                                             *
*     ***********************************************

      subroutine cgmontecarlo_end(rtdb,flag)
      implicit none
      integer rtdb
      integer flag

*     **** local variables ****
      logical value,calc_pressure
      integer minimizer


*     **** external functions ****
      logical  control_pressure,psi_finalize
      integer  control_minimizer,control_version
      external control_pressure,psi_finalize
      external control_minimizer,control_version


      minimizer     = control_minimizer()
      calc_pressure = control_pressure().and.(control_version().eq.3)

*     **** finalize pressure ****
      if (calc_pressure) then
         call psp_stress_end()
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
      call nwpw_meta_gga_end()
      call pspw_end_SIC()
      call pspw_end_HFX()
      call pspw_end_APC()
      call pspw_charge_end()
      call pspw_qmmm_end()

      call ion_end()
      call psp_end()
      call ion_end_FixIon()
      !call D3dB_n_fft_end(1)
      call D3dB_pfft_end()
      call D3dB_end(1)
      if (control_version().eq.4) call D3dB_end(2)
      call Dne_end()

      return
      end
