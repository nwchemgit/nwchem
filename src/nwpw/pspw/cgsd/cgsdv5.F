
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
*
* $Id: cgsdv5.F,v 1.1 2001-08-30 00:25:11 edo Exp $
*
      implicit none
      integer rtdb
      integer flag

#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"


*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - strfac_time
*     9 - maskC_time
*    10 - geodesic_time
      real*8 times(50)
      common / timing / times


      
*     **** parallel variables ****
      integer  taskid,np
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,nfft3d,n2ft3d
    

*     ***** energy variables ****
      real*8  E(10)
      real*8  dipole(3)
      real*8  stress(3,3)

*     **** gradient variables ****
      integer fion(2)

*     **** error variables ****
      logical value
      integer ierr

*     **** local variables ****
      real*8  gx,gy,gz,cx,cy,cz
      real*8  EV,pi
      integer i,k,ia,nion
  


*     **** external functions ****
      real*8      psp_zv,psp_amass,psp_rc,ewald_rcut
      real*8      ewald_mandelung
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg
      integer     ewald_ncut,ewald_nshl3d,psp_npsp
      integer     psp_lmmax,psp_lmax,psp_locp
      character*2 psp_atom
      external    psp_zv,psp_amass,psp_rc,ewald_rcut
      external    ewald_mandelung
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg
      external    ewald_ncut,ewald_nshl3d,psp_npsp
      external    psp_lmmax,psp_lmax,psp_locp
      external    psp_atom


      real*8   control_tole,control_tolc,control_tolr,ion_rion
      external control_tole,control_tolc,control_tolr,ion_rion
      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      logical  control_read,control_move,ion_init
      external control_read,control_move,ion_init

      integer  pack_nwave_all
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      integer  Waterpsp_nwater
      external pack_nwave_all
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      external Waterpsp_nwater

      character*12 control_boundry
      external     control_boundry

      logical      psp_semicore
      real*8       psp_rcore,psp_ncore,psp_rlocal
      external     psp_semicore
      external     psp_rcore,psp_ncore,psp_rlocal

      logical  psi_initialize, psi_finalize, Waterpsp_found
      integer  psi_ispin,psi_ne,electron_count
      external psi_initialize, psi_finalize, Waterpsp_found
      external psi_ispin,psi_ne,electron_count
      real*8   cgsd_energy
      external cgsd_energy
      logical  control_Mulliken
      external control_Mulliken
*                            |************|
*****************************|  PROLOGUE  |****************************
*                            |************|

      value = .true.
      pi = 4.0d0*datan(1.0d0)

      call dcopy(50,0.0d0,0,times,1)
      call dcopy(10,0.0d0,0,E,1)


*     **** get parallel variables ****
      call Parallel_Init()
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) call current_second(cpu1)

*     ***** print out header ****
      if (taskid.eq.MASTER) then
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
         write(6,1042)
         write(6,1043)
         write(6,1010)
         write(6,1000)
         call message(1)
         write(6,1110)
      end if
      
      value = control_read(3,rtdb)
      if (.not. value) 
     > call errquit('error reading control',0)

      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0


      ierr = 0 
      if ((np.gt.(control_ngrid(3)))) then
        if (taskid.eq.MASTER)
     >     write(6,*) 'Error: the number of processors must ',
     >                 ' be in the range ( 1 ...',
     >                   control_ngrid(3),')'
        ierr=100
        value = .false.
        go to 9000
      end if


*     **** initialize D3dB data structure ****
      call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3))
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d

*     ***** Initialize double D3dB data structure ****
      if (control_version().eq.4)
     >   call D3dB_Init(2,2*ngrid(1),2*ngrid(2),2*ngrid(3))


 
*     **** initialize lattice data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_init()


*     **** read in wavefunctions and initialize psi ****
      value = psi_initialize()

*     **** allocate psp data structure and read in psedupotentials into it ****
      call psp_init()
      call psp_readall()

*     **** read ions ****
      value = ion_init(rtdb)
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

*     **** initialize G,mask,ke,and coulomb data structures ****
      call ke_init()
      if (control_version().eq.3) call coulomb_init()
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()
      if (control_version().eq.3) call ewald_init()
      call electron_init()
      call geodesic_init()
      if (psp_semicore(0)) call semicore_check()

*     **** initialize QM/MM ****
      call pspw_init_APC(rtdb)
      call Waterpsp_init(rtdb)


*                |**************************|
******************   summary of input data  **********************
*                |**************************|

      if (taskid.eq.MASTER) then
         write(6,1111) np
         write(6,1115)
c         IF(control_move()) THEN
c           write(6,1120) 'yes'
c         ELSE
c           write(6,1120) 'no'
c         ENDIF
         write(6,1121) control_boundry(),control_version()
         if (psi_ispin().eq.1) write(6,1130) "restricted"
         if (psi_ispin().eq.2) write(6,1130) "unrestricted"
         IF (control_gga().eq.0) THEN
            write(6,1131) 'Vosko et al parameterization'
         ELSE IF (control_gga().eq.1) THEN
            write(6,1131) 
     >      'PBE96 (White and Bird) parameterization'
         ELSE
            write(6,1131) 'unknown parameterization'
            call errquit('bad exchange_correlation',0)
         END IF
         write(6,1140)
         do ia = 1,psp_npsp()
           write(6,1150) ia,psp_atom(ia),psp_amass(ia)/1822.89d0,
     >                    psp_zv(ia),psp_lmax(ia)
           write(6,1152) psp_lmax(ia)
           write(6,1153) psp_locp(ia)
           write(6,1154) psp_lmmax(ia)
           if (control_version().eq.4) write(6,1156) psp_rlocal(ia)
           if (psp_semicore(ia))
     >         write(6,1155) psp_rcore(ia),psp_ncore(ia)
           write(6,1151) (psp_rc(i,ia),i=0,psp_lmax(ia))
         end do
         write(6,1160)
         write(6,1170) (psp_atom(K),ion_natm(K),K=1,psp_npsp())
         write(6,1180)
         write(6,1190) (I,psp_atom(ion_katm(I)),
     >                  (ion_rion(K,I),K=1,3),I=1,ion_nion())
         write(6,1200) cx,cy,cz
         write(6,1210) gx,gy,gz

         call Waterpsp_Print(6)

         write(6,1220) psi_ne(1),psi_ne(psi_ispin()),' (fourier space)'
         write(6,1230)
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
         write(6,1290) control_tole(),control_tolc(),control_tolr()
         write(6,1300)
         call util_flush(6)
      end if

*                |***************************|
******************     call CG minimizer     **********************
*                |***************************|
      if (taskid.eq.MASTER) call current_second(cpu2)

c     if (flag.eq.1) then
c       call control_gradient_iterations()
c     end if

*     **** calculate energy ****
      EV= cgsd_energy()

*     **** extra energy output for QA test ****
      if (taskid.eq.MASTER) write(6,1600) EV

*     **** calculate the spin contamination ****
      call psi_spin2(dipole(1))
      
*     **** calculate the dipole ***
      dipole(1) = 0.0d0
      dipole(2) = 0.0d0
      dipole(3) = 0.0d0
      call rho_dipole(dipole)

*     ****  calculate gradient ***
      nion = ion_nion()
      if (Waterpsp_found()) nion = nion + 3*Waterpsp_nwater()

      value = MA_push_get(mt_dbl,(3*nion),
     >                       'fion',fion(2),fion(1))
      if (.not. value) call errquit('out of stack memory',0)
     
      call cgsd_energy_gradient(dbl_mb(fion(1)))
      call dscal(3*nion,(-1.0d0),dbl_mb(fion(1)),1)

*     **** calculate the stress tensor ****
      call dcopy(9,0.0d0,0,stress,1)
      if (flag.eq.3) then
         call psp_stress_init()
         call psp_stress_readall()
         call cgsd_energy_stress(stress)
         call psp_stress_end()
      end if


*     *************************************************************
*     **** output energy, dipole, and gradient to rtdb for use ****
*     **** by task_energy and task_gradient                    ****
*     *************************************************************
      value = rtdb_put(rtdb,'cg_pspw:energy',mt_dbl,1,EV)
      value = value.and.
     >        rtdb_put(rtdb,'cg_pspw:dipole',mt_dbl,
     >                 3,dipole)
      value = value.and.
     >        rtdb_put(rtdb,'cg_pspw:stress',mt_dbl,
     >                 9,stress)
      value = value.and.
     >        rtdb_put(rtdb,'cg_pspw:gradient',mt_dbl,
     >                 3*nion,dbl_mb(fion(1)))
      if (.not. value) call errquit('cgsdv5: error writing rtdb',0)
      value = MA_pop_stack(fion(2))
      if (.not. value) call errquit('erroring popping stack memory',0)    



      if (taskid.eq.MASTER) call current_second(cpu3)

*                |***************************|
******************         Epilogue          **********************
*                |***************************|

*     **** calculate Mulliken Populations ***
      if (control_Mulliken()) call psi_Mulliken(rtdb)

*     **** calculate Mulliken Populations ***

*     **** write geometry to rtdb ****
      call Waterpsp_write(rtdb)
      call ion_write(rtdb)

*     **** write wavefunctions to file and finalize psi ****
      value = psi_finalize()

*     **** deallocate heap memory ****
      call electron_finalize()
      call geodesic_finalize()
      if (control_version().eq.3) call ewald_end()
      call strfac_end()
      if (control_version().eq.3) call coulomb_end()
      if (control_version().eq.4) call coulomb2_end()
      call ke_end()
      call mask_end()
      call Pack_end()
      call G_end()
      call pspw_end_APC()
      call Waterpsp_end()
      call ion_end()
      call psp_end()
      call D3dB_end(1)
      if (control_version().eq.4) call D3dB_end(2)

*                |***************************|
****************** report consumed cputime   **********************
*                |***************************|
      if (taskid.eq.MASTER) then
         CALL current_second(cpu4)

         T1=CPU2-CPU1
         T2=CPU3-CPU2
         T3=CPU4-CPU3
         T4=CPU4-CPU1
         AV=T2/dble(electron_count())
         write(6,*)
         write(6,*) '-----------------'
         write(6,*) 'cputime in seconds'
         write(6,*) 'prologue    : ',T1
         write(6,*) 'main loop   : ',T2
         write(6,*) 'epilogue    : ',T3
         write(6,*) 'total       : ',T4
         write(6,*) 'cputime/step: ',AV,
     >              '   (',electron_count(),' evalulations)'
         write(6,*)
         write(6,*) '-------------------------------'
         write(6,*) 'Time spent doing:'
         write(6,*) '  FFTs                       : ', times(1)
         write(6,*) '  dot products               : ', times(2)
         write(6,*) '  geodesic                   : ', times(10)
         write(6,*) '  exchange correlation       : ', times(4)
         write(6,*) '  local pseudopotentials     : ', times(5)
         write(6,*) '  non-local pseudopotentials : ', times(6)
         write(6,*) '  hartree potentials         : ', times(7)
         write(6,*) '  structure factors          : ', times(8)
         write(6,*) '  masking and packing        : ', times(9)
         write(6,*)
         CALL MESSAGE(4)
      end if 



      call Parallel_Finalize()
      cgsdv5 = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*     Car-Parrinello microcluster calculation      *')
 1030 FORMAT(10X,'*     [     conjugate gradient minimization   ]    *')
 1031 FORMAT(10X,'*     [   (Grassman manifold implementation)  ]    *')
 1035 FORMAT(10x,'*      [ NorthWest Chemistry implementation ]      *')
 1040 FORMAT(10X,'*            version #5.00   06/01/00              *')
 1041 FORMAT(10X,'*    This code was developed by Eric J. Bylaska,   *')
 1042 FORMAT(10X,'*    and was based upon algorithms and code        *')
 1043 FORMAT(10X,'*    developed by the group of Prof. John H. Weare *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1111 FORMAT(/' number of processors used:',I3)
 1115 FORMAT(/' options:')
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundry conditions   = ',A,'(version', I1,')')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1140 FORMAT(/' elements involved in the cluster:')
 1150 FORMAT(5X,I2,': ',A2,'  mass no.:',F6.1,'   core charge:',F4.0,
     &       '  lmax=',I1)
 1151 FORMAT(5X,'        cutoff =',4F8.3)
 1152 FORMAT(12X,' highest angular component      : ',i2)
 1153 FORMAT(12X,' local potential used           : ',i2)
 1154 FORMAT(12X,' number of non-local projections: ',i2)
 1155 FORMAT(12X,' semicore corrections included  : ',
     >       F6.3,' (radius) ',F6.3,' (charge)')
 1156 FORMAT(12X,' aperiodic cutoff radius        : ',F6.3)
 1160 FORMAT(/' atomic composition:')
 1170 FORMAT(7(5X,A2,':',I3))
 1180 FORMAT(/' initial position of ions:')
 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' )')
 1200 FORMAT(5X,'  G.C. ',' (',3F11.5,' )')
 1210 FORMAT(5X,' C.O.M.',' (',3F11.5,' )')
 1220 FORMAT(/' number of electrons: spin up=',I3,'  spin down=',I3,A)
 1230 FORMAT(/' supercell:')
 1231 FORMAT(5x,' volume : ',F10.1)
 1241 FORMAT(5x,' lattice: a1=<',3f8.3,' >')
 1242 FORMAT(5x,'          a2=<',3f8.3,' >')
 1243 FORMAT(5x,'          a3=<',3f8.3,' >')
 1244 FORMAT(5x,'          b1=<',3f8.3,' >')
 1245 FORMAT(5x,'          b2=<',3f8.3,' >')
 1246 FORMAT(5x,'          b3=<',3f8.3,' >')

 1250 FORMAT(5X,' density cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
 1251 FORMAT(5X,' wavefnc cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
  
 1260 FORMAT(5X,' ewald summation: cut radius=',F8.2,'  and',I3)
 1261 FORMAT(5X,'                   mandelung=',f11.8)
 1270 FORMAT(/' technical parameters:')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'fictacious mass=',F10.1)
 1290 FORMAT(5X, ' tolerance=',E8.3,' (energy)',E12.3,
     &        ' (density)',E12.3,' (ion)')
 1300 FORMAT(//)
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
 1600 FORMAT(/' Total PSPW energy   :',E19.10)
 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(6,9010) ierr
      call Parallel_Finalize()

      cgsdv5 = value
      return
      END
