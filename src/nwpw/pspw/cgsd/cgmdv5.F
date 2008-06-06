*
* $Id: cgmdv5.F,v 1.9 2008-06-06 18:06:03 d3p708 Exp $
*

***********************************************************************
*                      cgmdv5                                         *

*     This is a developing Born-Oppenheimer MD code for NWCHEM        *
*                                                                     *
***********************************************************************

      logical function cgmdv5(rtdb,flag)
      implicit none
#include "errquit.fh"
      integer rtdb
      integer flag

#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "frac_occ.fh"

#include "nwpw_timing.fh"

*     *** local variables and parameters ****
      double precision kb
      parameter (kb=3.16679d-6)

      
*     **** parallel variables ****
      integer  taskid,np
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,nfft3d,n2ft3d
      real*8  a,b,c,alpha,beta,gamma
    

*     ***** energy variables ****
      integer ne(2),ms
      real*8  E(30),en(2)
      real*8  dipole(3)
      real*8  stress(3,3),lstress(6)

*     **** gradient variables ****
      integer fion(2)

*     **** error variables ****
      logical value,newpsi
      integer ierr

*     **** local variables ****
      logical oprint,lprint,hprint,qmmm,found
      real*8  gx,gy,gz,cx,cy,cz
      real*8  vgx,vgy,vgz,vcx,vcy,vcz,ekg,eki0,eki1,dt
      real*8  sum,w,emotion_time_shift,wa,wr
      real*8  EV,EV0,pi,esum1,esum2,eave,evar,dr
      real*8  ratio,aratio
      real*8  icharge
      integer i,k,ia,nion,minimizer,mapping,it,icount_shift
      integer nbq,egrid(3)

      character*30 filename
      character*255 full_filename


  


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
      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      logical  control_read,control_move,ion_init
      external control_read,control_move,ion_init
      integer  control_ewald_ngrid
      external control_ewald_ngrid
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
      character*30 control_cell_name
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
      logical   pspw_charge_found,ion_q_FixIon
      external  pspw_charge_found,ion_q_FixIon
      integer  control_minimizer,control_scf_algorithm
      external control_minimizer,control_scf_algorithm
      integer  control_ks_algorithm
      external control_ks_algorithm
      real*8   control_ks_alpha
      external control_ks_alpha
      logical  control_print,control_balance
      external control_print,control_balance
      integer  control_mapping
      external control_mapping
      integer  control_ks_maxit_orb,control_ks_maxit_orbs
      external control_ks_maxit_orb,control_ks_maxit_orbs

      logical pspw_bqext
      external pspw_bqext

      integer  control_bo_steps,control_bo_algorithm
      real*8   control_bo_time_step,control_rti,ion_ke,ion_Temperature
      external control_bo_steps,control_bo_algorithm
      external control_bo_time_step,control_rti,ion_ke,ion_Temperature

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

      oprint= ((taskid.eq.MASTER).and.control_print(print_medium))
      lprint= ((taskid.eq.MASTER).and.control_print(print_low))
      hprint= ((taskid.eq.MASTER).and.control_print(print_high))


      if (oprint) call current_second(cpu1)

*     ***** print out header ****
      if (oprint) then
         write(6,1000)
         write(6,1010)
         write(6,1020)
         write(6,1010)
         write(6,1030)
         write(6,1031)
         write(6,1010)
         write(6,1035)
         write(6,1010)
         write(6,1040)
         write(6,1010)
         write(6,1041)
         write(6,1010)
         write(6,1000)
         call nwpw_message(1)
         write(6,1110)
      end if
      
      egrid(1) = control_ewald_ngrid(1)
      egrid(2) = control_ewald_ngrid(2)
      egrid(3) = control_ewald_ngrid(3)
 
      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0
      minimizer = control_minimizer()
      mapping   = control_mapping()
      !qmmm      = control_qmmm()


*     **** initialize D3dB data structure ****
      call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d

*     ***** Initialize double D3dB data structure ****
      if (control_version().eq.4)
     >   call D3dB_Init(2,2*ngrid(1),2*ngrid(2),2*ngrid(3),mapping)


 
*     **** initialize lattice data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_init()

      call D3dB_pfft_init()


*     **** read in wavefunctions and initialize psi ****
c     value = psi_initialize()

*     **** read ions ****
      value = ion_init(rtdb)

*     **** allocate psp data structure and read in psedupotentials into it ****
      call psp_init()
      call psp_readall()
      if (psp_semicore(0)) call semicore_check()

*     **** initialize G,mask,ke,coulomb,and ewald data structures ****
      call ke_init()
      if (control_version().eq.3) call coulomb_init()
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()
      if (control_version().eq.3) call ewald_init(egrid)

*     **** read in wavefunctions and initialize psi ****
      if (.not.control_check_charge_multiplicity()) then
        call psi_new()
        newpsi = .true.
      else
        newpsi = .false.
      end if
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

*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)

*     **** initialize frac_occ data structure ****
      call frac_occ_init(rtdb,psi_ispin(),ne)

*     **** initialize linesearching ****
      call linesearch_init()


*     ******************************
*     **** scaling ion velocity ****
*     ******************************
      call ion_init_ke(ekg,eki0,eki1)

c     **** Initialize thermostats ****
c      if (control_Nose()) then
c         !call Nose_Init(w)
c      end if
c
c*     **** Initialize simulated annealing ****
c      SA       = .false.
c      Tr_init  = 0.0d0
c      sa_alpha(1) = 1.0d0
c      if (control_SA()) then
c         if (control_Nose()) then
c            SA          = .true.
c            sa_decay(1) = control_SA_decay(1)
c            Tr_init     = control_Nose_Tr()
c         else
c            dt = control_time_step()
c            SA          = .false.
c            sa_decay(1) = control_SA_decay(1)
c            sa_alpha(1) = dexp( -(dt/control_SA_decay(1)) )
c         end if
c      end if






*                |**************************|
******************   summary of input data  **********************
*                |**************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      call center_v_geom(vcx,vcy,vcz)
      call center_v_mass(vgx,vgy,vgz)

      if (oprint) then
         write(6,1111) np
         if (mapping.eq.1) write(6,1112) 
         if (mapping.eq.2) write(6,1113) 
         if (control_balance()) then
           write(6,1114)
         else
           write(6,1116)
         end if
         write(6,1115)

         write(6,1121) control_boundry(),control_version()
         if (psi_ispin().eq.1) write(6,1130) "restricted"
         if (psi_ispin().eq.2) write(6,1130) "unrestricted"
         !if (qmmm) write(6,1122)
         IF (control_gga().eq.-1) THEN
            write(6,1131) 'off'
         ELSE IF (control_gga().eq.0) THEN
            write(6,1131) 'LDA (Vosko et al) parameterization'
         ELSE IF (control_gga().eq.10) THEN
            write(6,1131) 
     >      'PBE96 (White and Bird) parameterization'
         ELSE IF (control_gga().eq.11) THEN
            write(6,1131) 
     >      'BLYP (White and Bird) parameterization'
         ELSE IF (control_gga().eq.12) THEN
            write(6,1131) 
     >      'revPBE (White and Bird) parameterization'

         ELSE IF (control_gga().eq.110) THEN
            write(6,1131) 
     >      'PBE0 (White and Bird) parameterization'
         ELSE IF (control_gga().eq.112) THEN
            write(6,1131) 
     >      'revPBE0 (White and Bird) parameterization'

         ELSE IF (control_gga().eq.200) THEN
            write(6,1131) 'Hartree-Fock'
         ELSE
            write(6,1131) 'unknown parameterization'
            call errquit('bad exchange_correlation',0, INPUT_ERR)
         END IF
         if (fractional) write(6,1132)
         call pspw_print_SIC(6)
         call pspw_print_HFX(6)
         write(6,1140)
         do ia = 1,ion_nkatm()
           write(6,1150) ia,ion_atom(ia),
     >                    psp_zv(ia),psp_lmax(ia)
           comment = psp_comment(ia)
           i = inp_strlen(comment)
           write(6,1157) comment(1:i)
           write(6,1158) psp_psp_type(ia)
           write(6,1152) psp_lmax(ia)
           write(6,1153) psp_locp(ia)
           write(6,1154) psp_nprj(ia)
           if (control_version().eq.4) write(6,1156) psp_rlocal(ia)
           if (psp_semicore(ia))
     >         write(6,1155) psp_rcore(ia),psp_ncore(ia)
           write(6,1151) (psp_rc(i,ia),i=0,psp_lmax(ia))
           write(6,*)
         end do

         icharge = -(psi_ne(1)+psi_ne(psi_ispin()))
         en(1)     = psi_ne(1)
         en(psi_ispin()) = psi_ne(psi_ispin())
         if (fractional) then
            icharge = 0
            do ms=1,psi_ispin()
            en(ms) =0.0
            do i=1,ne(ms)
              k = fweight(1)+(i-1)+(ms-1)*ne(1)
              icharge = icharge - (3-psi_ispin())*dbl_mb(k)
              en(ms) = en(ms) + dbl_mb(k)
            end do
            end do
         end if
         icharge = icharge + ion_TotalCharge()
         write(6,1159) icharge

         write(6,1160)
         write(6,1170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
         !if (hprint) then
         write(6,1180)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(6,1191) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0
           else
           write(6,1190) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0
           end if
         end do
         write(6,1200) cx,cy,cz
         write(6,1210) gx,gy,gz

         call pspw_charge_Print(6)

         write(6,1181)
         write(6,1192) (I,ion_aname(I),
     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
         write(6,1200) vcx,vcy,vcz
         write(6,1210) vgx,vgy,vgz

         !end if

         !write(6,1220) psi_ne(1),psi_ne(psi_ispin()),' (fourier space)'
         if (fractional) then
          write(6,1219) en(1),en(psi_ispin()),' (   fractional)'
          write(6,1221) psi_ne(1),psi_ne(psi_ispin()),' (fourier space)'
         else
          write(6,1220) psi_ne(1),psi_ne(psi_ispin()),' (fourier space)'
          write(6,1221) psi_ne(1),psi_ne(psi_ispin()),' (fourier space)'
         end if
         write(6,1230)
         write(6,1233) control_cell_name()
         write(6,1241) lattice_unita(1,1),
     >                 lattice_unita(2,1),
     >                 lattice_unita(3,1)
         write(6,1242) lattice_unita(1,2),
     >                 lattice_unita(2,2),
     >                 lattice_unita(3,2)
         write(6,1243) lattice_unita(1,3),
     >                 lattice_unita(2,3),
     >                 lattice_unita(3,3)
         write(6,1244) lattice_unitg(1,1),
     >                 lattice_unitg(2,1),
     >                 lattice_unitg(3,1)
         write(6,1245) lattice_unitg(1,2),
     >                 lattice_unitg(2,2),
     >                 lattice_unitg(3,2)
         write(6,1246) lattice_unitg(1,3),
     >                 lattice_unitg(2,3),
     >                 lattice_unitg(3,3)
         call lattice_abc_abg(a,b,c,alpha,beta,gamma)
         write(6,1232) a,b,c,alpha,beta,gamma
         write(6,1231) lattice_omega()
         write(6,1250) lattice_ecut(),ngrid(1),ngrid(2),ngrid(3),
     >                 pack_nwave_all(0),pack_nwave(0)
         write(6,1251) lattice_wcut(),ngrid(1),ngrid(2),ngrid(3),
     >                 pack_nwave_all(1),pack_nwave(1)
         if (control_version().eq.3) then
         write(6,1260) ewald_rcut(),ewald_ncut()
         write(6,1261) ewald_mandelung()
         end if

         write(6,1270)
         write(6,1280) control_time_step(),control_fake_mass()
         write(6,1290) control_tole(),control_tolc()
         write(6,1281) control_it_in()*control_it_out(),
     >                 control_it_in(),control_it_out()
         if (minimizer.eq.5) then
           write(6,1291)
           if (control_ks_algorithm().eq.0) 
     >       write(6,1292) "conjugate gradient"
           if (control_ks_algorithm().eq.1) 
     >       write(6,1292) "rmm-diis"
           write(6,1295) control_ks_maxit_orb(),control_ks_maxit_orbs()
           if (control_scf_algorithm().eq.0) 
     >       write(6,1293) "simple mixing"
           if (control_scf_algorithm().eq.1) 
     >       write(6,1293) "Anderson potential mixing"
           if (control_scf_algorithm().eq.2) 
     >       write(6,1293) "Johnson-Pulay mixing"
           if (control_scf_algorithm().eq.3) 
     >       write(6,1293) "Anderson density mixing"
           write(6,1294) control_ks_alpha()
         end if
         write(6,1310)
         write(6,1320) control_bo_time_step(), control_bo_steps()
         if (control_bo_algorithm().eq.0) write(6,1330) "Verlet"

         write(6,1340) control_rti()
         write(6,1222) eki0,ekg
         write(6,1223) eki1
         write(6,1224) (eki1-eki0)
         write(6,1300)
         call util_flush(6)
      end if

*                |***************************|
******************     simple MD loop        **********************
*                |***************************|
      if (taskid.eq.MASTER) call current_second(cpu2)
      if (taskid.eq.MASTER) call nwpw_message(10)


*     **** write initial position to xyz data ****
      call xyz_init()          ! unit=18
      call MOTION_init(rtdb)   ! unit=19



*     ************************************
*     **** open up other MOTION files ****
*     ************************************
*     **** open EMOTION file ****
      if (.not.rtdb_cget(rtdb,'nwpw:emotion_filename',1,filename))
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
           !E(25) = E(25) + sum                          !*** take care of running sums ***
           !E(26) = E(26) + sum*sum
           icount_shift = icount_shift + 1
           end do
  100      continue
#if defined(FUJITSU_SOLARIS) || defined(PSCALE) || defined(__crayx1)
           backspace 31
#endif
         else
           open(unit=31,file=full_filename,form='formatted',
     >          status='new')
         end if
      end if


      if (.not.rtdb_cget(rtdb,'nwpw:emotion_filename',1,filename))
     >  call util_file_prefix('emotion',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)
      if (taskid.eq.MASTER) then
         open(unit=31,file=full_filename,form='formatted',
     >        status='new')
      end if

      call xyz_write()


*     ****  allocate fion ***
      nion = ion_nion()
      if (pspw_charge_found().and.
     >    (.not.pspw_bqext())) nion = nion + pspw_charge_nion()

      value = MA_push_get(mt_dbl,(3*nion),
     >                       'fion',fion(2),fion(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


      dt = control_bo_time_step()


c*     ***** do a simple md loop ***
c      call control_reduce_print()
c      EV = cgsd_energy(newpsi)
c      call cgsd_energy_gradient(dbl_mb(fion(1)))
c      call ion_newton_step(dbl_mb(fion(1)),1.0d0)
c      call xyz_write()
c
c      do it=1,control_bo_steps()
c        call ion_shift()
c        EV = cgsd_energy(.false.)
c        call cgsd_energy_gradient(dbl_mb(fion(1)))
c        call ion_verlet_step(dbl_mb(fion(1)),1.0d0)
c        eki1 = ion_ke()
c        write(6,1350) it,EV+eki1,EV,eki1,ion_Temperature()
c
c        write(31,1360) it*dt+emotion_time_shift,
c     >                 EV+eki1,EV,eki1,ion_Temperature()
c
c        call util_flush(6)
c        call util_flush(31)
c
c        call xyz_write()
c        call MOTION_write(it*dt)
c
c
c      end do
c*     ***** do a simple md loop ***
c      call control_up_print()


*     ***** do a simple mc loop ***
      call control_reduce_print()
      dr = 0.08d0
      aratio = 0.5d0
      wa = 0.0d0
      wr = 0.0d0
      esum1 = 0.0d0
      esum2 = 0.0d0
      eave = 0.0d0
      evar = 0.0d0
      beta = util_random(39) !**seed the random number generator
      beta = 1.0d0/(kb*298.15)
      call xyz_write()
      EV0 = cgsd_energy(newpsi)
      do it=1,control_bo_steps()

        call ion_MC_step(dr)
        EV = cgsd_energy(.false.)
        w = dexp(-beta*(EV-EV0))
        alpha = util_random(0)

        !*** accept the sequence ***
        if (alpha.lt.w) then
          EV0 = EV
          !write(6,*) "   ----ACCEPTED----"
          if (taskid.eq.MASTER)  then
             write(6,1350) it,EV,EV0,w,alpha,wa,wr,eave,evar,dr
             call util_flush(6)
          end if
          call xyz_write()
          call MOTION_write(it*dt)
          wa = wa + 1.0d0
          esum1 = esum1 + EV
          esum2 = esum2 + EV*EV
          eave = esum1/wa
          evar = esum2/wa
          evar = evar - eave*eave

        !*** reject the sequence ***
        else
           call ion_MC_reject_step()
          !write(6,*) "   ----REJECTED----"
          wr = wr + 1.0d0
        end if

        !*** adjust ***
        if (mod(it,10).eq.0) then
            ratio = wa/(wa+wr)
            dr = dr*(ratio/aratio)
            !if (ratio.lt.aratio) dr = dr*0.95d0
            !if (ratio.gt.aratio) dr = dr*1.05d0
        end if

      end do

      call control_up_print()
*     ***** do a simple mc loop ***

      call xyz_end()
      call MOTION_end()
      if (taskid.eq.MASTER) then
        close(unit=31)
      end if

      EV = cgsd_energy(.false.)
      call cgsd_energy_gradient(dbl_mb(fion(1)))


*                |***************************|
******************     simple MD loop        **********************
*                |***************************|




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
      if (flag.gt.-1) call rho_dipole(dipole)


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
      value = rtdb_put(rtdb,'pspw:energy',mt_dbl,1,EV)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:dipole',mt_dbl,
     >                 3,dipole)
      end if
      if (flag.gt.0) then
        value = value.and.
     >        rtdb_put(rtdb,'pspw:gradient',mt_dbl,
     >                 3*nion,dbl_mb(fion(1)))
      end if
      if (flag.eq.3) then
        value = value.and.
     >        rtdb_put(rtdb,'pspw:stress',mt_dbl,
     >                 9,stress)
        value = value.and.
     >        rtdb_put(rtdb,'pspw:lstress',mt_dbl,
     >                 6,lstress)
      end if
      value = value.and.MA_pop_stack(fion(2))
      if (.not. value) call errquit('cgmdv5: error writing rtdb',0,
     &       RTDB_ERR)

      if (taskid.eq.MASTER) call current_second(cpu3)

*                |***************************|
******************         Epilogue          **********************
*                |***************************|

*     **** calculate Mulliken Populations ***
      if (control_Mulliken()) call psi_Mulliken(rtdb)




*     **** calculate Mulliken Populations ***

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
      call pspw_end_SIC()
      call pspw_end_HFX()
      call pspw_end_APC()
      call pspw_charge_end()
      call pspw_qmmm_end()
      call frac_occ_end()
      call ion_end()
      call psp_end()
      call ion_end_FixIon()
      !call D3dB_n_fft_end(1)
      call D3dB_pfft_end()
      call D3dB_end(1)
      if (control_version().eq.4) call D3dB_end(2)


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
      cgmdv5 = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*               NWPW PSPW Calculation              *')
 1030 FORMAT(10X,'*     [      Born-Oppenheimer molecular       ]    *')
 1031 FORMAT(10X,'*     [          dynamics simulation          ]    *')
 1035 FORMAT(10x,'*     [  NorthWest Chemistry implementation   ]    *')
 1040 FORMAT(10X,'*            version #1.00   07/18/06              *')
 1041 FORMAT(10X,'*    This code was developed by Eric J. Bylaska.   *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1111 FORMAT(/' number of processors used:',I3)
 1112 FORMAT( ' parallel mapping         : slab')
 1113 FORMAT( ' parallel mapping         : hilbert')
 1114 FORMAT( ' parallel mapping         : balanced')
 1115 FORMAT(/' options:')
 1116 FORMAT( ' parallel mapping         : not balanced')
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1122 FORMAT(5X,' qmmm simulation')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1132 FORMAT(5X,' using fractional occupation')
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
 1170 FORMAT(7(5X,A4,':',I3))
 1180 FORMAT(/' initial position of ions (au):')
 1181 FORMAT(/' initial velocity of ions after scaling (au):')

 1190 FORMAT(5X, I4, A4  ,' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1191 FORMAT(5X, I4, A4, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - fixed')
 1192 FORMAT(5X, I4, A4  ,' (',3F11.5,' )')

 1200 FORMAT(5X,'  G.C. ',' (',3F11.5,' )')
 1210 FORMAT(5X,' C.O.M.',' (',3F11.5,' )')
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
 1231 FORMAT(5x,'             omega=',F8.1)
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
  
 1260 FORMAT(5X,' ewald summation: cut radius=',F8.2,'  and',I3)
 1261 FORMAT(5X,'                   madelung=',f11.8)
 1270 FORMAT(/' technical parameters for minimizer:')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'ficticious mass=',F10.1)
 1281 FORMAT(5X, ' maximum iterations =',I8, 
     >           ' ( ',I4,' inner ',I4,' outer )')
 1290 FORMAT(5X, ' tolerance=',E8.3,' (energy)',E12.3,
     &        ' (density)')
 1291 FORMAT(/' Kohn-Sham scf parameters:')
 1292 FORMAT(5X, ' Kohn-Sham algorithm  = ',A)
 1293 FORMAT(5X, ' SCF algorithm        = ',A)
 1294 FORMAT(5X, ' SCF mixing parameter =',F7.4)
 1295 FORMAT(5X, ' Kohn-Sham iterations = ',I3,
     >           ' (',I3,' outer)')
 1300 FORMAT(//)
 1310 FORMAT(/' molecular dynamics parameters:')
 1320 FORMAT(5X, ' time step=',F10.2,5X,'iterations=',I8)
 1330 FORMAT(5X, ' integration agorithm=   ',A)
 1340 FORMAT(/5X, ' cooling/heatting rate= ',e12.5,' (ion)')
 1350 FORMAT('@',I8,2E19.10,4F8.2,2E19.10,F8.4)
 1360 format(100e19.10)
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

      cgmdv5 = value
      return
      END
