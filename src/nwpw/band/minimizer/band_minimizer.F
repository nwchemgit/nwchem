*
* $Id: band_minimizer.F,v 1.3 2001-12-12 03:17:02 bylaska Exp $
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
      real*8 ks(3),kv(3)

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
*     **** external functions ****
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg,ion_amass
      character*2 ion_aname
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg,ion_amass
      external    ion_aname


      real*8   control_tole,control_tolc,control_tolr,ion_rion
      external control_tole,control_tolc,control_tolr,ion_rion
      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      logical  control_read,control_move,ion_init
      external control_read,control_move,ion_init

      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm

      character*12 control_boundry
      external     control_boundry

      integer  control_num_kvectors
      external control_num_kvectors

      integer  Cram_nwave_all,Cram_nwave
      external Cram_nwave_all,Cram_nwave

      integer  c_ewald_ncut
      real*8   c_ewald_rcut,c_ewald_mandelung,c_ewald_e
      external c_ewald_ncut
      external c_ewald_rcut,c_ewald_mandelung,c_ewald_e



*****************************|  PROLOGUE  |****************************

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
         write(6,1031)
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
      
      value = control_read(5,rtdb)
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
      call C3dB_Init(1,ngrid(1),ngrid(2),ngrid(3))
      call C3dB_nfft3d(1,nfft3d)

 
*     **** initialize lattice data structure ****
      call lattice_init()
      call G_init()
      call Cram_init()



*     **** read ions ****
      value = ion_init(rtdb)
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

*     **** allocate psp data structure and read in psedupotentials into it ****
c     call cpsp_init()
c     call cpsp_readall()
c     if (cpsp_semicore(0)) call semicore_check()

*     **** initialize ke,and coulomb data structures ****
      call c_ewald_init()



*                |**************************|
******************   summary of input data  **********************
*                |**************************|

      value = rtdb_parallel(.false.)
      if (taskid.eq.MASTER) then
         write(6,1111) np
         write(6,1115)
         write(6,1121) control_boundry(),control_version()
         IF (control_gga().eq.0) THEN
            write(6,1131) 'Vosko et al parameterization'
         ELSE IF (control_gga().eq.1) THEN
            write(6,1131)
     >      'PBE96 (White and Bird) parameterization'
         ELSE
            write(6,1131) 'unknown parameterization'
            call errquit('bad exchange_correlation',0)
         END IF
         write(6,1180)
         write(6,1190) (I,ion_aname(I),
     >                  (ion_rion(K,I),K=1,3),
     >                  ion_amass(I)/1822.89d0,
     >                 I=1,ion_nion())
         write(6,1200) cx,cy,cz
         write(6,1210) gx,gy,gz


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
         write(6,1260) c_ewald_rcut(),c_ewald_ncut()
         write(6,1261) c_ewald_mandelung()

         ia = control_num_kvectors()
         write(6,1255)
         write(6,1256) ia
         do i=1,ia
           call control_ksvector(i,ks)
           call control_kvector(i,kv)
           write(6,1257) ks,kv
         end do

         write(6,1249)
         write(6,1250) lattice_ecut(),ngrid(1),ngrid(2),ngrid(3),
     >                 Cram_nwave_all(0),Cram_nwave(0)
         do i=1,ia
           write(6,1251) i,lattice_wcut(),ngrid(1),ngrid(2),ngrid(3),
     >                 Cram_nwave_all(i),Cram_nwave(i)
         end do


         write(6,1270)
         write(6,1280) control_time_step(),control_fake_mass()
         write(6,1290) control_tole(),control_tolc(),control_tolr()
         write(6,1300)
         call util_flush(6)
      end if
      value = rtdb_parallel(.true.)






*                |***************************|
******************     call CG minimizer     **********************
*                |***************************|
      if (taskid.eq.MASTER) call current_second(cpu2)


*     **** calculate energy ****
      EV= 0.0d0

*     **** extra energy output for QA test ****
      if (taskid.eq.MASTER) write(6,1600) EV

*     **** calculate the spin contamination ****
      
*     **** calculate the dipole ***
      dipole(1) = 0.0d0
      dipole(2) = 0.0d0
      dipole(3) = 0.0d0

*     ****  calculate gradient ***
      nion = 3

      value = MA_push_get(mt_dbl,(3*nion),
     >                       'fion',fion(2),fion(1))
      if (.not. value) call errquit('out of stack memory',0)
     

*     **** calculate the stress tensor ****
      call dcopy(9,0.0d0,0,stress,1)


*     *************************************************************
*     **** output energy, dipole, and gradient to rtdb for use ****
*     **** by task_energy and task_gradient                    ****
*     *************************************************************
      value = rtdb_put(rtdb,'band:energy',mt_dbl,1,EV)
      value = value.and.
     >        rtdb_put(rtdb,'band:dipole',mt_dbl,
     >                 3,dipole)
      value = value.and.
     >        rtdb_put(rtdb,'band:stress',mt_dbl,
     >                 9,stress)
      value = value.and.
     >        rtdb_put(rtdb,'band:gradient',mt_dbl,
     >                 3*nion,dbl_mb(fion(1)))
      if (.not. value) 
     > call errquit('band_minimizer: error writing rtdb',0)
      value = MA_pop_stack(fion(2))
      if (.not. value) call errquit('erroring popping stack memory',0)    



      if (taskid.eq.MASTER) call current_second(cpu3)

*                |***************************|
******************         Epilogue          **********************
*                |***************************|

      
      call cphafac()
      E(5) = c_ewald_e()
      if (taskid.eq.MASTER) then
      write(6,1470) E(5),E(5)/ion_nion()
      end if

*     **** calculate Mulliken Populations ***

*     **** write geometry to rtdb ****
      call ion_write(rtdb)

*     **** write wavefunctions to file and finalize psi ****

*     **** deallocate heap memory ****
      call c_ewald_end()
c     call cpsp_end()
      call Cram_end()
      call G_end()
      call ion_end()
      call C3dB_end(1)

*                |***************************|
****************** report consumed cputime   **********************
*                |***************************|
      if (taskid.eq.MASTER) then
         CALL current_second(cpu4)

         T1=CPU2-CPU1
         T2=CPU3-CPU2
         T3=CPU4-CPU3
         T4=CPU4-CPU1
c        AV=T2/dble(electron_count())
         write(6,*)
         write(6,*) '-----------------'
         write(6,*) 'cputime in seconds'
         write(6,*) 'prologue    : ',T1
         write(6,*) 'main loop   : ',T2
         write(6,*) 'epilogue    : ',T3
         write(6,*) 'total       : ',T4
c        write(6,*) 'cputime/step: ',AV,
c    >              '   (',electron_count(),' evalulations, ', 
c    >                     linesearch_count(),' line searches)'
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
      band_minimizer = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*           NWPW Band Structure Calculation        *')
 1031 FORMAT(10X,'*     [   (Stiefel manifold implementation)   ]    *')
 1040 FORMAT(10X,'*            version #1.00   11/16/01              *')
 1041 FORMAT(10X,'*    This code was developed by Eric J. Bylaska    *')
 1042 FORMAT(10X,'*    and Edoardo Apra                              *')
 1043 FORMAT(10X,'*                                                  *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1111 FORMAT(/' number of processors used:',I3)
 1115 FORMAT(/' options:')
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundry conditions   = ',A,'(version', I1,')')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1140 FORMAT(/' elements involved in the cluster:')
 1150 FORMAT(5X,I2,': ',A2,'  core charge:',F4.1,'  lmax=',I1)
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
 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' ) - atomic mass= ',F6.3,' ')
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

 1249 FORMAT(/' computation grids:')
 1250 FORMAT(5X,' density     cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
 1251 FORMAT(5X,' wavefnc ',I3,' cutoff=',F7.3,
     &        '  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
  
 1255 FORMAT(/' brillouin zone:')
 1256 FORMAT(5x,'number of zone points:',I3)
 1257 FORMAT(5x,'   ks=<',3f8.3,' >, k=<',3f8.3,'>')

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

      band_minimizer = value
      return
      END
