*
* $Id$
*

***********************************************************************
*                      mmmdv1                                         *

*     This is a developing Born-Oppenheimer MD code for NWCHEM        *
*                                                                     *
***********************************************************************

      logical function mmmdv1(rtdb,flag)
      implicit none
      integer rtdb
      integer flag

#include "global.fh"
#include "bafdecls.fh"
#include "errquit.fh"
#include "btdb.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"

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
      real*8  en(2),E(60)
      real*8  dipole(3)
      real*8  stress(3,3),lstress(6)

*     **** gradient variables ****
      integer fion(2),fion1(2)

*     **** error variables ****
      logical value,newpsi
      integer ierr

*     **** local variables ****
      logical oprint,lprint,hprint,qmmm,found,nose,SA,vverlet
      logical calc_pressure,field_exist
      real*8  gx,gy,gz,cx,cy,cz
      real*8  vgx,vgy,vgz,vcx,vcy,vcz,ekg,eki0,eki1,dt
      real*8  sum0,sum1,sum2,sum3,w,emotion_time_shift,wa,wr
      real*8  EV,EV0,pi,esum1,esum2,eave,evar,dr,qave,qvar,have,hvar
      real*8  ratio,aratio
      real*8  icharge,r,cv,E25,E26
      real*8  sa_alpha(2),sa_decay(2),ssr,tt2
      real*8  Tr_new,Tr_init,p1,p2,pressure
      integer i,k,ia,nion,minimizer,mapping,it,icount_shift,icount
      integer nbq,it_in,it_out

      character*50 filename
      character*255 full_filename


*     **** external functions ****
      real*8      mmq_zv,ewald_rcut
      real*8      ewald_mandelung
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg,ion_amass
      integer     ewald_ncut,ewald_nshl3d,ewald_npack
      integer     ewald_grid_nx,ewald_grid_ny,ewald_grid_nz
      integer     psp_lmmax,psp_lmax,psp_locp,ion_nkatm
      integer     psp_nprj,psp_psp_type
      character*4 ion_atom,ion_aname
      external    mmq_zv,ewald_rcut
      external    ewald_mandelung
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg,ion_amass
      external    ewald_ncut,ewald_nshl3d,ewald_npack
      external    ewald_grid_nx,ewald_grid_ny,ewald_grid_nz
      external    psp_lmmax,psp_lmax,psp_locp,ion_nkatm
      external    psp_nprj,psp_psp_type
      external    ion_atom,ion_aname
      
      external      psp_comment
      character*255 psp_comment,comment

      real*8   control_tole,control_tolc,control_tolr,ion_rion,ion_vion
      external control_tole,control_tolc,control_tolr,ion_rion,ion_vion
      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      logical  control_read,control_move,ion_init,control_out_of_time
      external control_read,control_move,ion_init,control_out_of_time
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      integer  pspw_charge_nion
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      external pspw_charge_nion

      character*12 control_boundry
      external     control_boundry
      character*50 control_cell_name
      external     control_cell_name

      real*8   md_energy
      external md_energy
      real*8   ion_TotalCharge
      external ion_TotalCharge
      logical  control_check_charge_multiplicity,control_rotation
      external control_check_charge_multiplicity,control_rotation
      logical   pspw_charge_found,ion_q_FixIon,pspw_Efield_found
      external  pspw_charge_found,ion_q_FixIon,pspw_Efield_found
      integer  control_minimizer,control_scf_algorithm
      external control_minimizer,control_scf_algorithm
      integer  control_ks_algorithm
      external control_ks_algorithm
      real*8   control_ks_alpha
      external control_ks_alpha
      logical  control_print,control_balance,control_translation
      external control_print,control_balance,control_translation
      integer  control_mapping
      external control_mapping

      logical pspw_bqext,pspw_qmmm_found
      external pspw_bqext,pspw_qmmm_found

      integer  control_bo_steps_in,control_bo_steps_out
      integer  control_bo_algorithm
      real*8   control_bo_time_step,control_rti,ion_ke,ion_Temperature
      external control_bo_steps_in,control_bo_steps_out
      external  control_bo_algorithm
      external control_bo_time_step,control_rti,ion_ke,ion_Temperature
      real*8   ion_com_Temperature
      external ion_com_Temperature
      logical  control_Nose,control_SA,control_pressure,Nose_restart
      real*8   control_Nose_Te,control_Nose_Tr,control_SA_decay,Nose_ssr
      real*8   Nose_Pr,Nose_Pe,Nose_Ee0,Nose_Er0,Nose_dXr
      real*8   Nose_Tr,Nose_Te,Nose_Qr,Nose_Qe,Nose_r_energy
      integer  Nose_Mchain,Nose_Nchain
      external control_Nose,control_SA,control_pressure,Nose_restart
      external control_Nose_Te,control_Nose_Tr,control_SA_decay,Nose_ssr
      external Nose_Pr,Nose_Pe,Nose_Ee0,Nose_Er0,Nose_dXr
      external Nose_Tr,Nose_Te,Nose_Qr,Nose_Qe,Nose_r_energy
      external Nose_Mchain,Nose_Nchain
      integer  ion_nconstraints,ion_ndof
      external ion_nconstraints,ion_ndof
      real*8   pspw_qmmm_lambda
      external pspw_qmmm_lambda
      integer  pspw_Efield_type
      external pspw_Efield_type


*                            |************|
*****************************|  PROLOGUE  |****************************
*                            |************|

      value = .true.
      pi = 4.0d0*datan(1.0d0)


      call nwpw_timing_init()


*     **** get parallel variables ****
      call Parallel_Init()
      call Parallel_np(np)
      call Parallel_taskid(taskid)

      value = control_read(12,rtdb)
      if (.not. value) 
     > call errquit('error reading control',0, INPUT_ERR)

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
         write(luout,1031)
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
      

*     **** initialize lattice data structure ****
      call lattice_init()

*     **** read ions ****
      value = ion_init(rtdb)

*     **** read in ion charges ****
      call mmq_init(rtdb)

*     **** initialize ewald data structures ****
      if (control_version().eq.3) call ewald_init()

*     **** initialize QM/MM ****
      call pspw_qmmm_init(rtdb)
      call pspw_charge_init(rtdb)


      field_exist = pspw_charge_found()

*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)


*     **** initialize pressure ****
      calc_pressure = control_pressure().and.(control_version().eq.3)
      pressure      = 0.0d0
      p1            = 0.0d0
      p2            = 0.0d0


*     ******************************
*     **** scaling ion velocity ****
*     ******************************
      call ion_init_ke(ekg,eki0,eki1)


*     **** Initialize thermostats ****
      nose = control_Nose()
      if (nose) call Nose_Init(2,0.01d0)

*     **** Initialize simulated annealing ****
      SA       = .false.
      Tr_init  = 0.0d0
      sa_alpha(2) = 1.0d0
      if (control_SA()) then
         if (nose) then
            SA          = .true.
            sa_decay(2) = control_SA_decay(2)
            Tr_init     = control_Nose_Tr()
         else
            dt = control_bo_time_step()
            SA          = .false.
            sa_decay(2) = control_SA_decay(2)
            sa_alpha(2) = dexp( -(dt/control_SA_decay(2)) )
         end if
      end if

      vverlet = (control_bo_algorithm().eq.1)



*                |**************************|
******************   summary of input data  **********************
*                |**************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      call center_v_geom(vcx,vcy,vcz)
      call center_v_mass(vgx,vgy,vgz)

      if (oprint) then
         write(luout,1111) np
         write(luout,1115)

         write(luout,1121) control_boundry(),control_version()

         write(luout,1140)
         do ia = 1,ion_nkatm()
           write(luout,1150) ia,ion_atom(ia),mmq_zv(ia)
         end do

         icharge = 0.0d0
         icharge = icharge + ion_TotalCharge()
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
         write(luout,1211) ion_nconstraints(),ion_ndof()

         call pspw_charge_Print(luout)

         write(luout,1181)
         write(luout,1192) (I,ion_aname(I),
     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
         write(luout,1200) vcx,vcy,vcz
         write(luout,1210) vgx,vgy,vgz

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
         if (control_version().eq.3) then
         write(luout,1260) ewald_rcut(),ewald_ncut()
         write(luout,1261) ewald_mandelung()
         write(luout,1252) 
     >        ewald_grid_nx(),ewald_grid_ny(),ewald_grid_nz(),
     >        ewald_npack(),ewald_npack()
         end if

c         write(luout,1270)
c         write(luout,1280) control_time_step(),control_fake_mass()
c         write(luout,1290) control_tole(),control_tolc()
c         write(luout,1281) control_it_in()*control_it_out(),
c     >                 control_it_in(),control_it_out()

         write(luout,1310)
         if (.not.control_translation()) write(luout,1271)
         if (.not.control_rotation())    write(luout,1272)
         write(luout,1320) control_bo_time_step(),
     >    control_bo_steps_in()*control_bo_steps_out(),
     >    control_bo_steps_in(),control_bo_steps_out()
         if(control_bo_algorithm().eq.0) write(6,1330) "Position Verlet"
         if(control_bo_algorithm().eq.1) write(6,1330) "Velocity Verlet"
         if(control_bo_algorithm().eq.2) write(6,1330) "Leap Frog"


         write(luout,1340) control_rti()
         call ion_scaling_atoms_print(luout)
         write(luout,1222) eki0,ekg
         write(luout,1223) eki1
         write(luout,1224) (eki1-eki0)
         if (nose) then
           write(luout,1395)
           if (Nose_restart()) then
              write(luout,*) "    thermostats resused"
           else
              write(luout,*) "    thermostats initialized"
           end if
           do i=1,Nose_Nchain()
             write(luout,1398) i,control_Nose_Tr(),Nose_Qr(i),
     >                   Nose_Pr(i),Nose_Er0(i)
           end do
         else
           write(luout,1394)
         end if
        if (calc_pressure) write(6,1293)
         if (control_SA()) then
           write(luout,1396) sa_decay(2)
         end if
         write(luout,1300)
         call util_flush(luout)
      end if

*                |***************************|
******************     simple MD loop        **********************
*                |***************************|
      call dcopy(50,0.0d0,0,E,1)
      if (taskid.eq.MASTER) call current_second(cpu2)
      if (taskid.eq.MASTER) call nwpw_message(10)


*     **** write initial position to xyz data ****
      call xyz_init()          ! unit=18
      call MOTION_init(rtdb)   ! unit=19

*     ************************************
*     **** open up other MOTION files ****
*     ************************************
*     **** open EMOTION file ****
      E(23)=0.0d0
      E(24)=0.0d0
      E(25)=0.0d0
      E(26)=0.0d0
      icount_shift = 0
      if (.not.btdb_cget(rtdb,'nwpw:emotion_filename',1,filename))
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
           read(31,*,end=100) emotion_time_shift,w,sum0,sum1,sum2,sum3
           E(25) = E(25) + sum0   !*** take care of running sums ***
           E(26) = E(26) + sum0*sum0
           E(27) = E(27) + (sum0+sum1+sum2)
           E(28) = E(28) + (sum0+sum1+sum2)**2
           E(23) = E(23) + sum3
           E(24) = E(24) + sum3*sum3
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

      call xyz_write()


*     ****  allocate fion ***
      nion = ion_nion()
      if (pspw_charge_found().and.
     >    (.not.pspw_bqext())) nion = nion + pspw_charge_nion()

      value = BA_push_get(mt_dbl,(3*nion),
     >                       'fion',fion(2),fion(1))
      if (vverlet) then
         value = value.and.BA_push_get(mt_dbl,(3*nion),
     >                                'fion1',fion1(2),fion1(1))
      end if
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)



*     ***** do a simple md loop ***
      dt = control_bo_time_step()
      call control_reduce_print()
      !EV = md_energy(E)
      call md_energy_gradient(dbl_mb(fion(1)))

      r = 1.0d0 
      if (nose) r =  (1.0d0-0.5d0*dt*Nose_dXr())
      call ion_newton_step(dbl_mb(fion(1)),sa_alpha(2)*r)
      eki1 = ion_ke()
      if (nose) call Nose_Newton_Step(0.01d0,eki1)

      call xyz_write()

      it_out = control_bo_steps_out()
      it_in  = control_bo_steps_in()
      icount = 0
      if (it_out.lt.1) goto 102


  101 continue
         icount = icount + 1

c       **** inner loop ****
        do i=1,it_in
           if (vverlet) then
              call dcopy(3*nion,dbl_mb(fion(1)),1,dbl_mb(fion1(1)),1)
              call ion_shift21()
           else
              call ion_shift()
              if (nose) call Nose_shift()
           end if

           !EV = md_energy(E)
           call md_energy_gradient(dbl_mb(fion(1)))
           if (nose) then
              ssr = Nose_ssr()
              call ion_nose_step(ssr,dbl_mb(fion(1)))
              eki1 = ion_ke()
              call Nose_Verlet_Step(0.01d0,eki1)
           else if (vverlet) then
              call ion_vverlet_step(dbl_mb(fion(1)),
     >                              dbl_mb(fion1(1)))
              call ion_vshift()
              call ion_newton_step(dbl_mb(fion(1)),sa_alpha(2))
              eki1 = ion_ke()
           else
              call ion_verlet_step(dbl_mb(fion(1)),sa_alpha(2))
              eki1 = ion_ke()
           end if
        end do

        !**** calculate energies ****
        EV = md_energy(E)
        E(1) = E(2) + eki1
        E(4) = eki1
        if (nose) then
           E(10) = Nose_r_energy()
           E(1) = E(1) + E(10)
        end if

        !*** take care of running sums ***
        E(23) = E(23) + E(14)
        E(24) = E(24) + E(14)*E(14)
        E(25) = E(25) + E(2)
        E(26) = E(26) + E(2)*E(2)
        E(27) = E(27) + E(1)
        E(28) = E(28) + E(1)*E(1)

        if (oprint) then
           if (nose) then
              write(luout,1350) icount*it_in,E(1),
     >                          E(2),E(4),ion_Temperature()
           else
              write(luout,1350) icount*it_in,E(1),E(2),E(4),
     >                          ion_Temperature()
           end if
           call util_flush(luout)
        end if

*       **** write out EMOTION data ****
        qave = E(23)/dble(icount+icount_shift)
        qvar = E(24)/dble(icount+icount_shift)
        qvar = qvar - qave*qave
        eave = E(25)/dble(icount+icount_shift)
        evar = E(26)/dble(icount+icount_shift)
        evar = evar - eave*eave
        have = E(27)/dble(icount+icount_shift)
        hvar = E(28)/dble(icount+icount_shift)
        hvar = hvar - have*have
        if (nose) then
           write(31,1360) icount*it_in*dt+emotion_time_shift,
     >                    E(1),E(2),E(3),E(4),E(14),E(5),E(6),
     >                    E(7),E(8),E(9),E(10),
     >                    eave,evar,have,hvar,qave,qvar,
     >                    ion_Temperature()
        else
           write(31,1360) icount*it_in*dt+emotion_time_shift,
     >                    E(1),E(2),E(3),E(4),E(14),E(5),E(6),
     >                    E(7),E(8),
     >                    eave,evar,have,hvar,qave,qvar,
     >                    ion_Temperature()

        end if
        call util_flush(31)

*       **** write xyz, MOTION data ****
        call xyz_write()
        call MOTION_write(icount*it_in*dt)


*        **** update thermostats using SA decay ****
         if (SA) then
           tt2 = icount*it_in*dt/sa_decay(2)
           Tr_new = Tr_init*dexp(-tt2)
           call Nose_reset_T(Tr_new,Tr_new)
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


*::::::::::::::::::::  end of iteration loop  :::::::::::::::::::::::::

*     **** close xyz and MOTION files ****
      call xyz_end()
      call MOTION_end()
      if (taskid.eq.MASTER) then
        close(unit=31)
      end if

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
           else
           write(luout,1190) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz
         write(luout,1211) ion_nconstraints(),ion_ndof()

         write(luout,1421)
         write(luout,1192) (I,ion_aname(I),
     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
         write(luout,1200) vcx,vcy,vcz
         write(luout,1210) vgx,vgy,vgz

         write(luout,1304)
         write(luout,1410)
         write(luout,*)
         write(luout,1430) E(1),E(1)/ion_nion()


         if (field_exist) then
            write(luout,1431)
            write(luout,1432)
            write(luout,1433) (E(1)-E(49)-E(50)-E(51)),
     >         (E(1)-E(49)-E(50)-E(51))/ion_nion()
         end if


         write(luout,1470) E(5),E(5)/ion_nion()
         if (pspw_qmmm_found()) then
            write(luout,1700)
            write(luout,1701)
            write(luout,1702) E(11)
            write(luout,1703) E(12)
            write(luout,1704) E(13)
            qave = E(23)/dble(icount+icount_shift)
            qvar = E(24)/dble(icount+icount_shift)
            qvar = qvar - qave*qave
            write(luout,1706) pspw_qmmm_lambda()
            write(luout,1705) E(14),qave,qvar
         end if
      end if
      call ecce_print1 ('total energy', mt_dbl, E(1), 1)
      call ecce_print1 ('nuclear repulsion energy', mt_dbl, E(5), 1)


      if (oprint) then

         write(luout,1300)
         write(luout,1472) ion_ke(),ion_ke()/ion_nion()

*       **** write out Temperatures ****
         eki0 = ion_Temperature()
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
      end if


*     **** extra energy output for QA test ****
      if (lprint) write(luout,1600) EV

*     **** calculate the spin contamination ****

*     **** calculate the dipole ***

*     **** calculate the stress tensor ****


*     *************************************************************
*     **** output energy, dipole, and gradient to rtdb for use ****
*     **** by task_energy and task_gradient                    ****
*     *************************************************************
      if (flag.gt.-1) then
      value = btdb_put(rtdb,'pspw:energy',mt_dbl,1,EV)
      value = value.and.
     >        btdb_put(rtdb,'pspw:dipole',mt_dbl,
     >                 3,dipole)
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
      if (vverlet) value = value.and.BA_pop_stack(fion1(2))
      value = value.and.BA_pop_stack(fion(2))
      if (.not. value) call errquit('cgmdv5: error writing rtdb',0,
     &       RTDB_ERR)

      if (taskid.eq.MASTER) call current_second(cpu3)

*                |***************************|
******************         Epilogue          **********************
*                |***************************|

      call ion_write(rtdb)

*     **** deallocate heap memory ****
      if (control_version().eq.3) call ewald_end()
   
      call mmq_end()
      call pspw_charge_end()
      call pspw_qmmm_end()
      if (nose) call Nose_end()
      call ion_end()
      call ion_end_FixIon()

*     **** do anaylysis on MOTION files ****
c      call cpmd_properties(rtdb)



*                |***************************|
****************** report consumed cputime   **********************
*                |***************************|
      if (lprint) then
         CALL current_second(cpu4)

         T1=CPU2-CPU1
         T2=CPU3-CPU2
         T3=CPU4-CPU3
         T4=CPU4-CPU1
         AV=T2/dble(icount*it_in)
         write(luout,1801)
         write(luout,1802)
         write(luout,1803) T1
         write(luout,1804) T2
         write(luout,1805) T3
         write(luout,1806) T4
         write(luout,1807) AV,icount*it_in

         call nwpw_timing_print_final(oprint,icount*it_in)

         write(luout,*)
         CALL nwpw_MESSAGE(4)
      end if 


      call Parallel_Finalize()
      mmmdv1 = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*               NWPW MD   Calculation              *')
 1030 FORMAT(10X,'*     [      Born-Oppenheimer molecular       ]    *')
 1031 FORMAT(10X,'*     [          dynamics simulation          ]    *')
 1035 FORMAT(10x,'*     [  NorthWest Chemistry implementation   ]    *')
 1040 FORMAT(10X,'*            version #1.00   07/18/06              *')
 1041 FORMAT(10X,'*    This code was developed by Eric J. Bylaska.   *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1111 FORMAT(/' number of processors used:',I3)
 1115 FORMAT(/' options:')
 1116 FORMAT( ' parallel mapping         : not balanced')
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1122 FORMAT(5X,' qmmm simulation')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1132 FORMAT(5X,' using fractional occupation')
 1140 FORMAT(/' elements involved in the cluster:')
 1150 FORMAT(5X,I2,': ',A4,'  core charge:',F9.6)
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

 1200 FORMAT(5X,'  G.C. ',' (',3F11.5,' )')
 1210 FORMAT(5X,' C.O.M.',' (',3F11.5,' )')
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
 1252 FORMAT(5X,' ewald fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
  
 1260 FORMAT(5X,' Ewald summation: cut radius=',F8.2,'  and',I3)
 1261 FORMAT(5X,'                   madelung=',f14.8)
 1270 FORMAT(/' technical parameters for minimizer:')
 1271 FORMAT(5x, ' translation constrained')
 1272 FORMAT(5x, ' rotation constrained')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'fictitious mass=',F10.1)
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
 1304 FORMAT(/)
 1310 FORMAT(/' molecular dynamics parameters:')
 1320 FORMAT(5X, ' time step=',F10.2,5X,'iterations=',I8,
     >           ' ( ',I4,' inner ',I4,' outer )')
 1330 FORMAT(5X, ' integration agorithm=   ',A)
 1340 FORMAT(/5X, ' cooling/heatting rate= ',e12.5,' (ion)')
 1350 FORMAT(I8,2E19.10,E14.5,4F14.2,2E19.10,F8.4)
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
 1430 FORMAT(/' total     energy    :',E19.10,' (',E15.5,'/ion)')
 1431 FORMAT(/' MD Energies')
 1432 FORMAT( '------------')
 1433 FORMAT( ' total  MD energy    :',E19.10,' (',E15.5,'/ion)')
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
 1600 FORMAT(/' Total MD   energy   :',E19.10)
 1700 FORMAT(/' QM/MM-pol-vib/CAV Energies')
 1701 FORMAT( ' --------------------------')
 1702 FORMAT( ' LJ energy              :',E19.10)
 1703 FORMAT( ' Residual Coulomb energy:',E19.10)
 1704 FORMAT( ' MM Vibration energy    :',E19.10)
 1705 FORMAT( ' QM/MM coupling energy  :',E19.10,
     >        ' (average=',E19.10,' variance=',E19.10,')'/)
 1706 FORMAT(/' QM/MM coupling lambda  :',E19.10)
 1707 FORMAT( ' - MM Charge Field/QM Electron    :',E19.10)
 1708 FORMAT( ' - MM Charge Field/QM Ion         :',E19.10)
 1709 FORMAT( ' - MM LJ/QM LJ                    :',E19.10)
 1710 FORMAT( ' - MM Charge Field/MM Charge Field:',E19.10)
 1711 FORMAT( ' - MM LJ/MM LJ                    :',E19.10)
 1801 FORMAT(//'== Timing ==')
 1802 FORMAT(/'cputime in seconds')
 1803 FORMAT( '  prologue    : ',E14.6)
 1804 FORMAT( '  main loop   : ',E14.6)
 1805 FORMAT( '  epilogue    : ',E14.6)
 1806 FORMAT( '  total       : ',E14.6)
 1807 FORMAT( '  cputime/step: ',E14.6,
     >        '       (',I8,' evalulations)')
 1808 FORMAT(A,E14.6,E14.6)
 1809 FORMAT(//A,2A14)

 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(6,9010) ierr
      call Parallel_Finalize()

      mmmdv1 = value
      return
      END
