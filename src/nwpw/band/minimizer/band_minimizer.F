*
* $Id: band_minimizer.F,v 1.16 2002-01-16 20:51:13 bylaska Exp $
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


      
*     **** parallel variables ****
      integer  taskid,np
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,nfft3d

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
      logical newpsi
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

      integer  brillioun_nbrillioun
      real*8   brillioun_weight,brillioun_ks,brillioun_k
      external brillioun_nbrillioun
      external brillioun_weight,brillioun_ks,brillioun_k
      integer  c_electron_count,linesearch_count
      external c_electron_count,linesearch_count

      real*8   nwpw_timing
      external nwpw_timing
      integer  Cram_nwave_all,Cram_nwave
      external Cram_nwave_all,Cram_nwave

      integer  c_ewald_ncut
      real*8   c_ewald_rcut,c_ewald_mandelung,c_ewald_e
      external c_ewald_ncut
      external c_ewald_rcut,c_ewald_mandelung,c_ewald_e
      logical  cpsp_semicore,psi_filefind,cpsi_initialize
      external cpsp_semicore,psi_filefind,cpsi_initialize
      real*8   c_cgsd_energy
      external c_cgsd_energy
      integer  cpsp_lmax,cpsp_locp,cpsp_lmmax,cpsp_npsp
      external cpsp_lmax,cpsp_locp,cpsp_lmmax,cpsp_npsp
      real*8   cpsp_rcore,cpsp_rc,cpsp_ncore,cpsp_zv
      external cpsp_rcore,cpsp_rc,cpsp_ncore,cpsp_zv
      character*2 cpsp_atom
      external    cpsp_atom
      integer  cpsi_ispin,cpsi_ne
      external cpsi_ispin,cpsi_ne


*****************************|  PROLOGUE  |****************************

      value = .true.
      pi = 4.0d0*datan(1.0d0)

      call nwpw_timing_init()
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
      call brillioun_init()
      call Cram_init()



*     **** read ions ****
      value = ion_init(rtdb)
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

*     **** allocate psp data structure and read in psedupotentials into it ****
      call cpsp_init()
      call cpsp_readall()
      if (cpsp_semicore(0)) call c_semicore_check()


*     **** initialize ke,and coulomb data structures ****
      call cstrfac_init()
      call cke_init()
      call c_coulomb_init()
      call c_ewald_init()

*     **** set up phase factors at the current geometry  ****
      call cphafac()
      call cphafac_k()


*     **** read in wavefunctions and initialize psi ****
      if (.not.psi_filefind()) then
        call cpsi_new()
        newpsi = .true.
      else
        newpsi = .false.
      end if
      value = cpsi_initialize()
      newpsi = .true.


*     **** electron and geodesic data structures ****
      call c_electron_init()
      call c_geodesic_init()
      call linesearch_init()



*                |**************************|
******************   summary of input data  **********************
*                |**************************|

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
         write(6,1140)
         do ia = 1,cpsp_npsp()
           write(6,1150) ia,cpsp_atom(ia),
     >                    cpsp_zv(ia),cpsp_lmax(ia)
           write(6,1152) cpsp_lmax(ia)
           write(6,1153) cpsp_locp(ia)
           write(6,1154) cpsp_lmmax(ia)
           if (cpsp_semicore(ia))
     >         write(6,1155) cpsp_rcore(ia),cpsp_ncore(ia)
           write(6,1151) (cpsp_rc(i,ia),i=0,cpsp_lmax(ia))
         end do

         write(6,1180)
         write(6,1190) (I,ion_aname(I),
     >                  (ion_rion(K,I),K=1,3),
     >                  ion_amass(I)/1822.89d0,
     >                 I=1,ion_nion())
         write(6,1200) cx,cy,cz
         write(6,1210) gx,gy,gz
         write(6,1220) cpsi_ne(1),cpsi_ne(cpsi_ispin()),
     >                 ' (fourier space)'


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

         ia = brillioun_nbrillioun()
         write(6,1255)
         write(6,1256) ia
         do i=1,ia
           write(6,1257) brillioun_weight(i),
     >        brillioun_ks(1,i),brillioun_ks(2,i),brillioun_ks(3,i),
     >        brillioun_k(1,i), brillioun_k(2,i), brillioun_k(3,i)
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






*                |***************************|
******************     call CG minimizer     **********************
*                |***************************|
      if (taskid.eq.MASTER) call current_second(cpu2)


*     **** calculate energy ****
      EV= c_cgsd_energy(newpsi)

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

      
*     **** calculate Mulliken Populations ***

*     **** write geometry to rtdb ****
      call ion_write(rtdb)

*     **** write wavefunctions to file and finalize psi ****
      call cpsi_finalize()

*     **** deallocate heap memory ****
      call c_electron_finalize()
      call c_geodesic_finalize()
      call c_ewald_end()
      call cstrfac_end()
      call c_coulomb_end()
      call cke_end()
      call cpsp_end()
      call Cram_end()
      call G_end()
      call brillioun_end()
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
         AV=T2/dble(c_electron_count())
         write(6,*)
         write(6,*) '-----------------'
         write(6,*) 'cputime in seconds'
         write(6,*) 'prologue    : ',T1
         write(6,*) 'main loop   : ',T2
         write(6,*) 'epilogue    : ',T3
         write(6,*) 'total       : ',T4
         write(6,*) 'cputime/step: ',AV,
     >              '   (',c_electron_count(),' evalulations, ', 
     >                     linesearch_count(),' line searches)'

         write(6,*)
         write(6,*) '-------------------------------'
         write(6,*) 'Time spent doing:'
         write(6,*) '  FFTs                       : ', nwpw_timing(1),
     >                          nwpw_timing(1)/dble(c_electron_count())
         write(6,*) '  dot products               : ', nwpw_timing(2),
     >                          nwpw_timing(2)/dble(c_electron_count())
         write(6,*) '  geodesic                   : ', nwpw_timing(10),
     >                          nwpw_timing(10)/dble(c_electron_count())
         write(6,*) '  exchange correlation       : ', nwpw_timing(4),
     >                          nwpw_timing(4)/dble(c_electron_count())
         write(6,*) '  local pseudopotentials     : ', nwpw_timing(5),
     >                          nwpw_timing(5)/dble(c_electron_count())
         write(6,*) '  non-local pseudopotentials : ', nwpw_timing(6),
     >                          nwpw_timing(6)/dble(c_electron_count())
         write(6,*) '  hartree potentials         : ', nwpw_timing(7),
     >                          nwpw_timing(7)/dble(c_electron_count())
         write(6,*) '  structure factors          : ', nwpw_timing(8),
     >                          nwpw_timing(8)/dble(c_electron_count())
         write(6,*) '  masking and packing        : ', nwpw_timing(9),
     >                          nwpw_timing(9)/dble(c_electron_count())
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
 1040 FORMAT(10X,'*            version #1.00   01/01/02              *')
 1041 FORMAT(10X,'*  Developed by Eric J. Bylaska and Edoardo Apra   *')
 1043 FORMAT(10X,'*                                                  *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1111 FORMAT(/' number of processors used:',I3)
 1115 FORMAT(/' options:')
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
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
 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1200 FORMAT(5X,'  G.C. ',' (',3F11.5,' )')
 1210 FORMAT(5X,' C.O.M.',' (',3F11.5,' )')
 1220 FORMAT(/' number of electrons: spin up=',I3,'  spin down=',I3,A)
 1230 FORMAT(/' supercell:')
 1231 FORMAT(5x,' volume : ',F10.1)
 1241 FORMAT(5x,' lattice:    a1=<',3f8.3,' >')
 1242 FORMAT(5x,'             a2=<',3f8.3,' >')
 1243 FORMAT(5x,'             a3=<',3f8.3,' >')
 1244 FORMAT(5x,' reciprocal: b1=<',3f8.3,' >')
 1245 FORMAT(5x,'             b2=<',3f8.3,' >')
 1246 FORMAT(5x,'             b3=<',3f8.3,' >')

 1249 FORMAT(/' computational grids:')
 1250 FORMAT(5X,' density     cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
 1251 FORMAT(5X,' wavefnc ',I3,' cutoff=',F7.3,
     &        '  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
  
 1255 FORMAT(/' brillouin zone:')
 1256 FORMAT(5x,'number of zone points:',I3)
 1257 FORMAT(5x,' weight=',f8.3,'  ks=<',3f8.3,' >, k=<',3f8.3,'>')

 1260 FORMAT(5X,' ewald summation: cut radius=',F8.2,'  and',I3)
 1261 FORMAT(5X,'                   mandelung=',f11.8)
 
 1270 FORMAT(/' technical parameters:')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'ficticious mass=',F10.1)
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
 1600 FORMAT(/' Total BAND energy   :',E19.10)
 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(6,9010) ierr
      call Parallel_Finalize()

      band_minimizer = value
      return
      END

