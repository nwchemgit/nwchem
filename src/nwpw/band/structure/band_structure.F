*
* $Id: band_structure.F,v 1.7 2004-12-21 16:58:35 bylaska Exp $
*
***********************************************************************
*                     band_structure								  *
*                                                                     *
*     This is a developing band structure parallel code for NWCHEM    *
*       + tcgmsg message passing library used                         *
*       + my own slap-decomposed parallel 3d-FFT(real->complex) used  *
*                                                                     *
*                                                                     *
***********************************************************************

      logical function band_structure(rtdb,flag)
      implicit none
      integer rtdb
      integer flag


#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"


      
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
      integer vall(2),nn,nb,ispin,ne(2)
      real*8  E(10),en(2)
      real*8  dipole(3)
      real*8  stress(3,3)

      integer eigs_dos(2),dosgrid(3)

*     **** gradient variables ****
      integer fion(2)

*     **** error variables ****
      logical value,ortho
      integer ierr

*     **** local variables ****
      logical newpsi
      real*8  gx,gy,gz,cx,cy,cz
      real*8  EV,pi
      real*8  pathlength,dist,kold(3)
      integer i,k,ia,nion,vers,nbrillioun,icharge,isize,indx
      integer mapping
      character*255 full_filename

  


*     **** external functions ****
*     **** external functions ****
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg,ion_amass,ion_TotalCharge
      character*2 ion_aname
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg,ion_amass,ion_TotalCharge
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
      logical  cpsp_semicore,psi_filefind,cpsi_initialize,cpsi_finalize
      external cpsp_semicore,psi_filefind,cpsi_initialize,cpsi_finalize
      real*8   c_cgsd_noit_energy,cpsi_1energy,cpsi_eigenvalue
      external c_cgsd_noit_energy,cpsi_1energy,cpsi_eigenvalue
      integer  cpsp_lmax,cpsp_locp,cpsp_lmmax,cpsp_npsp
      external cpsp_lmax,cpsp_locp,cpsp_lmmax,cpsp_npsp
      real*8   cpsp_rcore,cpsp_rc,cpsp_ncore,cpsp_zv
      external cpsp_rcore,cpsp_rc,cpsp_ncore,cpsp_zv
      character*2 cpsp_atom
      external    cpsp_atom
      integer  cpsi_ispin,cpsi_ne,psi_get_version
      external cpsi_ispin,cpsi_ne,psi_get_version
      logical  pspw_reformat_c_wvfnc
      external pspw_reformat_c_wvfnc
      integer  control_mapping
      external control_mapping
      integer  control_num_kvectors_structure,control_excited_ne
      external control_num_kvectors_structure,control_excited_ne
    


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
         write(luout,1000)
         write(luout,1010)
         write(luout,1020)
         write(luout,1010)
         write(luout,1040)
         write(luout,1010)
         write(luout,1041)
         write(luout,1043)
         write(luout,1010)
         write(luout,1000)
         call nwpw_message(1)
         write(luout,1110)
      end if
      
      value = control_read(5,rtdb)
      if (.not. value) 
     > call errquit('error reading control',0, DISK_ERR)

      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0
      mapping = control_mapping()




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


*     **** initialize C3dB data structure ****
      call C3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call C3dB_nfft3d(1,nfft3d)


 
*     **** initialize lattice data structure ****
      call lattice_init()
      call c_G_init()
      call brillioun_init()
      call Cram_init()


*     **** initialize D3dB data structure and mask for GGA ****
      IF (control_gga().eq.1) THEN
      call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),1)
      call G_init()
      call mask_init()
      END IF


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

*       **** convert from pspw format to band format ****
        vers = psi_get_version()
        if ((vers.eq.3).or.(vers.eq.4)) then
           nbrillioun = brillioun_nbrillioun()
           newpsi = .true.
           if (taskid.eq.MASTER) then
             value= pspw_reformat_c_wvfnc(1)
           end if
        end if
      end if
      value = cpsi_initialize(.true.)



*     **** electron and geodesic data structures ****
      call c_electron_init()
      call c_geodesic_init()
      call linesearch_init()





*                |**************************|
******************   summary of input data  **********************
*                |**************************|

      if (taskid.eq.MASTER) then
         write(luout,1111) np
         write(luout,1115)
         write(luout,1121) control_boundry(),control_version()
         IF (control_gga().eq.0) THEN
            write(luout,1131) 'Vosko et al parameterization'
         ELSE IF (control_gga().eq.1) THEN
            write(luout,1131)
     >      'PBE96 (White and Bird) parameterization'
         ELSE
            write(luout,1131) 'unknown parameterization'
            call errquit('bad exchange_correlation',0, INPUT_ERR)
         END IF
         write(luout,1140)
         do ia = 1,cpsp_npsp()
           write(luout,1150) ia,cpsp_atom(ia),
     >                    cpsp_zv(ia),cpsp_lmax(ia)
           write(luout,1152) cpsp_lmax(ia)
           write(luout,1153) cpsp_locp(ia)
           write(luout,1154) cpsp_lmmax(ia)
           if (cpsp_semicore(ia))
     >         write(luout,1155) cpsp_rcore(ia),cpsp_ncore(ia)
           write(luout,1151) (cpsp_rc(i,ia),i=0,cpsp_lmax(ia))
         end do


         icharge = -(cpsi_ne(1)+cpsi_ne(cpsi_ispin()))
     >           + ion_TotalCharge()
         write(luout,1159) icharge


         write(luout,1180)
         write(luout,1190) (I,ion_aname(I),
     >                  (ion_rion(K,I),K=1,3),
     >                  ion_amass(I)/1822.89d0,
     >                 I=1,ion_nion())
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz
         write(luout,1220) cpsi_ne(1),cpsi_ne(cpsi_ispin()),
     >                 ' (fourier space)'


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
         write(luout,1260) c_ewald_rcut(),c_ewald_ncut()
         write(luout,1261) c_ewald_mandelung()

         ia = brillioun_nbrillioun()
         write(luout,1255)
         write(luout,1256) ia
         do i=1,ia
           write(luout,1257) brillioun_weight(i),
     >        brillioun_ks(1,i),brillioun_ks(2,i),brillioun_ks(3,i),
     >        brillioun_k(1,i), brillioun_k(2,i), brillioun_k(3,i)
         end do

         write(luout,1249)
         write(luout,1250) lattice_ecut(),ngrid(1),ngrid(2),ngrid(3),
     >                 Cram_nwave_all(0),Cram_nwave(0)
         do i=1,ia
          write(luout,1251) i,lattice_wcut(),ngrid(1),ngrid(2),ngrid(3),
     >                 Cram_nwave_all(i),Cram_nwave(i)
         end do


         write(luout,1270)
         write(luout,1280) control_time_step(),control_fake_mass()
         write(luout,1290) control_tole(),control_tolc()
         write(luout,1300)
         call util_flush(luout)
      end if



      !**** set the size of the band - include virtual orbitals ****
      ispin = cpsi_ispin()
      ne(1) = cpsi_ne(1)+control_excited_ne(1)
      ne(2) = 0
      if (ispin.gt.1) ne(2) = cpsi_ne(2)+control_excited_ne(2)


      if (taskid.eq.MASTER) call current_second(cpu2)

*     **** allocate vall ****
      value = MA_alloc_get(mt_dcpl,2*nfft3d,'vall',vall(2),vall(1))
      if (.not. value)
     >  call errquit('band_structure:out of heap memory',0, MA_ERR)



      EV = c_cgsd_noit_energy()
      call c_electron_gen_vall()
      call c_electron_get_vall(dcpl_mb(vall(1)))

      if (taskid.eq.MASTER) then
         write(luout,1600) EV
         write(luout,*)
         write(luout,*) "Self-Consistent Potential Generated"
      end if



*     **** deallocate heap memory ****
      value=cpsi_finalize(.false.)

      call c_electron_finalize()
      call c_geodesic_finalize()
      call c_ewald_end()
      call cstrfac_end()
      call c_coulomb_end()
      call cke_end()
      call cpsp_end()
      call Cram_end()
      call c_G_end()
      call brillioun_end()


*     **** produce eigenvalue band file(s) ****
      if (ispin.eq.1) then
        call util_file_name('restricted_band',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=58,file=full_filename,form='formatted')
        end if
      else
        call util_file_name('alpha_band',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=58,file=full_filename,form='formatted')
        end if
        call util_file_name('beta_band',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=59,file=full_filename,form='formatted')
        end if
      end if


*     **** DOS calculation ****
      if (flag.eq.1) then
        if (taskid.eq.MASTER) write(luout,*) "DOS of states calculation"
        call control_dos_grid_structure(dosgrid)


*     **** band structure calculation ****
      else
       if (taskid.eq.MASTER) write(luout,*) "band structure calculation"
      end if


      call control_reset_band_structure()
      do k=1,control_num_kvectors_structure()
        ortho = .true.

        if (taskid.eq.MASTER) then
          write(luout,1301) k
        end if


*     **** initialize lattice data structure ****
      call lattice_init()
      call c_G_init()
      call brillioun_structure_init(k)
      call Cram_init()

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
        call cpsi_new_ne(ispin,ne)
        newpsi = .true.
        ortho  = .false.

      else
        newpsi = .false.

*       **** convert from pspw format to band format ****
        vers = psi_get_version()
        if ((vers.eq.3).or.(vers.eq.4)) then
           nbrillioun = brillioun_nbrillioun()
           newpsi = .true.
           if (taskid.eq.MASTER) then
             value= pspw_reformat_c_wvfnc(1)
           end if
        end if
      end if
      value = cpsi_initialize(ortho)


*     **** allocate eigs_dos if first iteration ****
      if ((flag.eq.1).and.(k.eq.1).and.(taskid.eq.MASTER)) then
        isize = dosgrid(1)*dosgrid(2)*dosgrid(3)*(cpsi_ne(1)+cpsi_ne(2))
        value = MA_alloc_get(mt_dbl,isize,
     >                       'eigs_dos',eigs_dos(2),eigs_dos(1))
        if (.not. value)
     >   call errquit('band_structure:out of heap memory',0, MA_ERR)
      end if



*     **** electron and geodesic data structures ****
      call c_electron_init()
      call c_geodesic_init()
      call linesearch_init()




*     **** diagonalize hamiltonian and rotate psi ****
      call c_electron_set_vall(dcpl_mb(vall(1)))
      call cpsi_KS_minimize(1,.false.,control_tole())
      call cpsi_check_indx(k)


      if (taskid.eq.MASTER) then

         !*** determine pathlength ****
         if (k.eq.1) then
           pathlength = 0.0d0
         else
            dist = (brillioun_k(1,1) - kold(1))**2
     >           + (brillioun_k(2,1) - kold(2))**2
     >           + (brillioun_k(3,1) - kold(3))**2
            dist = dsqrt(dist)     
            pathlength = pathlength + dist
         end if
         kold(1) = brillioun_k(1,1)
         kold(2) = brillioun_k(2,1)
         kold(3) = brillioun_k(3,1)

         NN=cpsi_ne(1)-cpsi_ne(2)
         EV=27.2116
         nb=1
           write(luout,1508) k,
     >                   pathlength,
     >                   brillioun_ks(1,nb),
     >                   brillioun_ks(2,nb),
     >                   brillioun_ks(3,nb),
     >                   brillioun_k(1,nb),
     >                   brillioun_k(2,nb),
     >                   brillioun_k(3,nb)
           write(luout,1500)
           do i=0,NN-1
             write(luout,1510) cpsi_eigenvalue(nb,1,cpsi_ne(1)-i),
     >                     cpsi_eigenvalue(nb,1,cpsi_ne(1)-i)*EV
           end do
           do i=0,cpsi_ne(2)-1
             write(luout,1510)  cpsi_eigenvalue(nb,1,cpsi_ne(1)-i-NN),
     >                      cpsi_eigenvalue(nb,1,cpsi_ne(1)-i-NN)*EV,
     >                      cpsi_eigenvalue(nb,2,cpsi_ne(2)-i),
     >                      cpsi_eigenvalue(nb,2,cpsi_ne(2)-i)*EV
           end do
           write(58,'(1000E14.6)') pathlength,
     >          (cpsi_eigenvalue(nb,1,i),i=1,cpsi_ne(1))
           if (ispin.eq.2)
     >      write(59,'(1000E14.6)') pathlength,
     >          (cpsi_eigenvalue(nb,2,i),i=1,cpsi_ne(2))
         
           if (flag.eq.1) then
            do i=1,cpsi_ne(1)
             indx = eigs_dos(1) 
     >            + (k-1) 
     >            + (i-1)*dosgrid(1)*dosgrid(2)*dosgrid(3)
             dbl_mb(indx) = cpsi_eigenvalue(nb,1,i)
            end do
            do i=1,cpsi_ne(2)
             indx = eigs_dos(1) 
     >            + (k-1) 
     >            + (i-1+cpsi_ne(1))*dosgrid(1)*dosgrid(2)*dosgrid(3)
             dbl_mb(indx) = cpsi_eigenvalue(nb,2,i)
            end do
           end if
                  
      end if



*     **** deallocate heap memory ****
      value = cpsi_finalize(.true.)

      call c_electron_finalize()
      call c_geodesic_finalize()
      call c_ewald_end()
      call cstrfac_end()
      call c_coulomb_end()
      call cke_end()
      call cpsp_end()
      call Cram_end()
      call c_G_end()
      call brillioun_end()

      end do
      if (taskid.eq.MASTER) then
        close(58)
        if (ispin.eq.2) close(59)
      end if


      if (taskid.eq.MASTER) call current_second(cpu3)

*                |***************************|
******************         Epilogue          **********************
*                |***************************|

      

*     **** deallocate heap memory ****

      value = MA_free_heap(vall(2)) 
      if ((flag.eq.1).and.(taskid.eq.MASTER))      ! deallocate eigs_dos
     >  value = value.and.MA_free_heap(eigs_dos(2)) 
      
      call ion_end()
      call C3dB_end(1)
      IF (control_gga().eq.1) THEN
      call mask_end()
      call G_end()
      call D3dB_end(1)
      end if

*                |***************************|
****************** report consumed cputime   **********************
*                |***************************|
      if (taskid.eq.MASTER) then
         CALL current_second(cpu4)

         T1=CPU2-CPU1
         T2=CPU3-CPU2
         T3=CPU4-CPU3
         T4=CPU4-CPU1
         write(luout,*)
         write(luout,*) '-----------------'
         write(luout,*) 'cputime in seconds'
         write(luout,*) 'prologue    : ',T1
         write(luout,*) 'main loop   : ',T2
         write(luout,*) 'epilogue    : ',T3
         write(luout,*) 'total       : ',T4

         write(luout,*)
         write(luout,*) '-------------------------------'
         write(luout,*) 'Time spent doing:'
         write(luout,*) '  FFTs                       : ', 
     >                          nwpw_timing(1)
         write(luout,*) '  dot products               : ', 
     >                          nwpw_timing(2)
         write(luout,*) '  geodesic                   : ', 
     >                          nwpw_timing(10)
         write(luout,*) '  exchange correlation       : ', 
     >                          nwpw_timing(4)
         write(luout,*) '  local pseudopotentials     : ', 
     >                          nwpw_timing(5)
         write(luout,*) '  non-local pseudopotentials : ', 
     >                          nwpw_timing(6)
         write(luout,*) '  hartree potentials         : ', 
     >                          nwpw_timing(7)
         write(luout,*) '  structure factors          : ', 
     >                          nwpw_timing(8)
         write(luout,*) '  masking and packing        : ', 
     >                          nwpw_timing(9)
         write(luout,*)
         CALL nwpw_MESSAGE(4)
      end if 

      call Parallel_Finalize()
      band_structure = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*           NWPW Band Structure Calculation        *')
 1040 FORMAT(10X,'*            version #1.00   11/02/03              *')
 1041 FORMAT(10X,'*          Developed by Eric J. Bylaska            *')
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
 1159 FORMAT(/' total charge:',I2)
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
 1261 FORMAT(5X,'                   madelung=',f11.8)
 
 1270 FORMAT(/' technical parameters:')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'ficticious mass=',F10.1)
 1290 FORMAT(5X, ' tolerance=',E8.3,' (energy)',E12.3,
     &        ' (density)')
 1300 FORMAT(//)
 1301 FORMAT(//'== Optimizing Brillouin Zone Point:',I6,' =='/)
 1304 FORMAT(/)
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
 1508 FORMAT(/' Brillouin zone point: ',i2,
     >       /'pathlength=',f10.6,
     >       /'    k     =<',3f8.3,'> . <b1,b2,b3> ',
     >       /'          =<',3f8.3,'>')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1600 FORMAT(/' Total BAND energy   :',E19.10)
 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(6,9010) ierr
      call Parallel_Finalize()

      band_structure = value
      return
      END

