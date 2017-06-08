*
* $Id$
*

***********************************************************************
*                      paw_minimizer                                  *
*                                                                     *
*     This is a developing paw  parallel code for NWCHEM              *
*       + tcgmsg message passing library used                         *
*       + my own slap-decomposed and hilbert-decompled parallel       *
*         3d-FFT(real->complex) are used                              *
*                                                                     *
*                                                                     *
***********************************************************************

      logical function paw_minimizer(rtdb,flag)
      implicit none
      integer rtdb
      integer flag

#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "inp.fh"
#include "util.fh"
#include "errquit.fh"

#include "nwpw_timing.fh"
      
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
      integer ne(2),ms,ispin
      real*8  E(10),en(2)
      real*8  dipole(3)
      real*8  stress(3,3)

*     **** gradient variables ****
      integer fion(2)

*     **** error variables ****
      logical value,newpsi
      integer ierr

*     **** local variables ****
      logical oprint,lprint,hprint
      real*8  gx,gy,gz,cx,cy,cz
      real*8  EV,pi
      real*8  icharge
      integer i,k,ia,nion,minimizer,mapping,lmax,mapping1d
  

*     **** external functions ****
      real*8      paw_mult_rcut,ion_amass
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg
      integer     paw_mult_ncut
      character   spdf_name
      character*4 ion_aname,ion_atom
      external    paw_mult_rcut,ion_amass
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg
      external    paw_mult_ncut
      external    spdf_name
      external    ion_aname,ion_atom


      real*8   control_tole,control_tolc,control_tolr,ion_rion
      external control_tole,control_tolc,control_tolr,ion_rion
      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      logical  control_read,control_move,ion_init,ion_q_FixIon
      external control_read,control_move,ion_init,ion_q_FixIon

      integer  pack_nwave_all
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,control_lmax_multipole
      integer  ion_nion,ion_natm,ion_katm,ion_nkatm
      external pack_nwave_all
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,control_lmax_multipole
      external ion_nion,ion_natm,ion_katm,ion_nkatm

      character*12 control_boundry
      external     control_boundry
      character*50 control_cell_name
      external     control_cell_name

      logical  paw_psi_initialize, paw_psi_finalize
      external paw_psi_initialize, paw_psi_finalize


      logical      control_Mulliken
      external     control_Mulliken
      logical      psi_filefind
      external     psi_filefind
      real*8       paw_psi_CheckOrtho
      external     paw_psi_CheckOrtho
      integer      control_mapping,control_mapping1d
      external     control_mapping,control_mapping1d

      integer  paw_electron_count,linesearch_count
      external paw_electron_count,linesearch_count
      
      real*8   ion_TotalCharge,paw_cgsd_energy,paw_cgsd_noit_energy
      external ion_TotalCharge,paw_cgsd_energy,paw_cgsd_noit_energy
      integer  control_minimizer,control_scf_algorithm
      external control_minimizer,control_scf_algorithm
      real*8   control_ks_alpha
      external control_ks_alpha
      logical  control_print
      external control_print

*



*                            |************|
*****************************|  PROLOGUE  |****************************
*                            |************|

      value = .true.
      pi = 4.0d0*datan(1.0d0)

      call nwpw_timing_init()
c
c     energy array
c     ------------
      call dcopy(10,0.0d0,0,E,1)


*     **** get parallel variables ****
      call Parallel_Init()
      call Parallel_np(np)
      call Parallel_taskid(taskid)

c
c     get stuff from rtdb    
c     -------------------
      value = control_read(8,rtdb)
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
         write(6,1010)
         write(6,1035)
         write(6,1010)
         write(6,1040)
         write(6,1010)
         write(6,1041)
         write(6,1042)
         write(6,1043)
         write(6,1044)
         write(6,1010)
         write(6,1000)
         call nwpw_message(1)
         write(6,1110)
      end if
      

      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0
      minimizer = control_minimizer()
      mapping   = control_mapping()


*     **** initialize D3dB data structure ****
      call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d

*     ***** Initialize double D3dB data structure ****
c     for aperiodic
      if (control_version().eq.4)
     >   call D3dB_Init(2,2*ngrid(1),2*ngrid(2),2*ngrid(3),mapping)


*     **** initialize lattice data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_init()

c
c     FFT for fermi sphere
c     --------------------
      call D3dB_pfft_init()


*     **** initialize Gaunt array ****
      call paw_gaunt_init()


*     **** read ions ****
      value = ion_init(rtdb)
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)


*     **** allocate paw data structure and read in paw basis into it ****
      call init_paw_basis()
      call init_paw_small_basis()

*     *** initialize paw matrices ***
      call init_paw_kin_matrix()
      call init_paw_ion_matrix()
      call init_paw_vloc_matrix()
      call init_paw_core_matrix()
      call init_paw_hartree_matrix()
      call init_paw_overlap_matrix()
      call init_paw_comp_charge_matrix()
      call init_paw_comp_pot_matrix()

*     *** initialize paw projectors ***
      call paw_proj_init()

*     *** initialize compensation charge ***
      call paw_comp_charge_init()

c*     *** initialize paw atomic potentials ***
c      call current_second(t1)
c      call paw_vloc_init()
c      call current_second(t2)
c      write(*,*) "time6:",t2-t1

c     call current_second(t1)
      call init_paw_pot_hartree()
c     call current_second(t2)
c     write(*,*) "time7:",t2-t1

c     call current_second(t1)
      call init_paw_pot_comp()
c     call current_second(t2)
c     write(*,*) "time8:",t2-t1


c     call current_second(t1)
      call paw_mult_init()
c     call current_second(t2)
c     write(*,*) "time9:",t2-t1


*     **** initialize G,mask,ke,and coulomb data structures ****
c     call current_second(t1)
      call ke_init()
      if (control_version().eq.3) call coulomb_init()
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()


*     ***** create wavefunctions if they don't exist ****
      if (.not.psi_filefind()) then
        call paw_psi_new()
        newpsi = .true.
      else
        newpsi = .false.
      end if

*     *** intialize overlap coefficient data structure ***
      call psi_get_ne(ispin,ne)
      mapping1d   = control_mapping1d()
      call Dne_init(ispin,ne,mapping1d)

      call phafac()
      call paw_ovlp_init(ispin,ne)
      call paw_nonlocal_init(ispin,ne)
      lmax = control_lmax_multipole()
      if (lmax.lt.0) lmax =  paw_basis_max_mult_l()
c      call init_paw_density(ispin,lmax)
      call init_paw_xc(ispin,lmax)

      !call paw_ovlp_coeff_set(dcpl_mb(psi2(1)))
      !call paw_ovlp_weights_set()


*     **** read in wavefunctions and initialize psi ****
      value = paw_psi_initialize()

*     **** electron and geodesic data structures ****
      call paw_electron_init()
      call paw_geodesic_init()

*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)

*     **** initialize linesearching ****
      call linesearch_init()

c     call current_second(t2)
c     write(*,*) "time10:",t2-t1



*                |**************************|
******************   summary of input data  **********************
*                |**************************|

      if (oprint) then
         write(6,1111) np
         if (mapping.eq.1) write(6,1112) 
         if (mapping.eq.2) write(6,1113) 
         if (mapping.eq.3) write(6,1118) 
         write(6,1115)

         write(6,1121) control_boundry(),control_version()
         if (ispin.eq.1) write(6,1130) "restricted"
         if (ispin.eq.2) write(6,1130) "unrestricted"
         IF (control_gga().eq.0) THEN
            write(6,1131) 'Vosko et al parameterization'
         ELSE IF (control_gga().eq.10) THEN
            write(6,1131) 
     >      'PBE96 (White and Bird) parameterization'
         ELSE IF (control_gga().eq.11) THEN
            write(6,1131) 
     >      'BLYP (White and Bird) parameterization'
         ELSE IF (control_gga().eq.12) THEN
            write(6,1131) 
     >      'revPBE (White and Bird) parameterization'
         ELSE
            write(6,1131) 'unknown parameterization'
            call errquit('bad exchange_correlation',0, INPUT_ERR)
         END IF


         write(6,1140)
         do ia = 1,ion_nkatm()
            write(6,1141) ia,ion_atom(ia),
     >                    paw_basis_ion_charge(ia),
     >                    paw_basis_core_charge(ia)
            
            write(6,1142) dbl_mb(paw_basis_i_r1(ia))
            write(6,1143) dbl_mb(paw_basis_i_rmax(ia))
            write(6,1144) paw_basis_ngrid(ia)
            write(6,1145) paw_basis_sphere_radius(ia), 
     >                    paw_basis_i_cut(ia),
     >                    paw_basis_npoints(ia)
            write(6,1146) paw_basis_sigma(ia)

            write(6,1150) paw_proj_nbasis(ia)
            write(6,1151)
            do i=1,paw_basis_nbasis(ia)
              write(6,1152) paw_basis_n_ps(i,ia),
     >                      paw_basis_n(i,ia),
     >                      spdf_name(paw_basis_orb_l(i,ia)),
     >                      paw_basis_eig(i,ia),
     >                      2*paw_basis_orb_l(i,ia)+1
            end do
         end do

         icharge = -(ne(1)+ne(ispin))
         do ia=1,ion_nkatm()
           icharge = icharge + ion_natm(ia)*
     >                        (paw_basis_ion_charge(ia)
     >                        -paw_basis_core_charge(ia))
         end do
         write(6,1159) icharge


         write(6,1160)
         write(6,1170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
         if (hprint) then
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

         end if

         
         write(6,1220) ne(1),ne(ispin),' (fourier space)'
         write(6,1221) ne(1),ne(ispin),' (fourier space)'
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
         write(6,1260) paw_mult_rcut(),paw_mult_ncut()
         write(6,1270)
         write(6,1262) lmax
         write(6,1280) control_time_step(),control_fake_mass()
         write(6,1290) control_tole(),control_tolc()
         if (minimizer.eq.5) then
           write(6,1291)
           write(6,1292)
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
      if (flag.eq.-1) then
         EV= paw_cgsd_noit_energy()
      else
         EV= paw_cgsd_energy(newpsi)
      end if


*     **** extra energy output for QA test ****
      if (lprint) write(6,1600) EV

      
*     ****  calculate gradient ***
      if (flag.gt.0) then
        nion = ion_nion()
        value = BA_push_get(mt_dbl,(3*nion),
     >                       'fion',fion(2),fion(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
        call paw_force_init()
        call paw_cgsd_energy_gradient(dbl_mb(fion(1)))
        call dscal(3*nion,(-1.0d0),dbl_mb(fion(1)),1)
        call paw_force_end()
      end if


*     *************************************************************
*     **** output energy, dipole, and gradient to rtdb for use ****
*     **** by task_energy and task_gradient                    ****
*     *************************************************************
      if (flag.gt.-1) then
      value = btdb_put(rtdb,'paw:energy',mt_dbl,1,EV)
      value = value.and.
     >        btdb_put(rtdb,'paw:dipole',mt_dbl,
     >                 3,dipole)
      end if
      if (flag.gt.0) then
        value = value.and.
     >        btdb_put(rtdb,'paw:gradient',mt_dbl,
     >                 3*nion,dbl_mb(fion(1)))
        value = value.and.BA_pop_stack(fion(2))
      end if
      if (flag.eq.3) then
        value = value.and.
     >        btdb_put(rtdb,'paw:stress',mt_dbl,
     >                 9,stress)
      end if
      if (.not. value) call errquit('cgsdv5: error writing rtdb',0,
     &       RTDB_ERR)

      if (taskid.eq.MASTER) call current_second(cpu3)

*                |***************************|
******************         Epilogue          **********************
*                |***************************|


*     **** calculate Mulliken Populations ***


*     **** write wavefunctions to file and finalize psi ****
      if (flag.eq.-1) then
        value = paw_psi_finalize(.false.)
      else
        value = paw_psi_finalize(.true.)
      end if




*     **** write geometry to rtdb ****
      call ion_write(rtdb)


*     **** deallocate heap memory ****
      call paw_electron_finalize()
      call paw_geodesic_finalize()

      call strfac_end()
      if (control_version().eq.3) call coulomb_end()
      if (control_version().eq.4) call coulomb2_end()
      call ke_end()
      call mask_end()
      call Pack_end()
      call G_end()
      call ion_end()
      call ion_end_FixIon()

      call dealloc_paw_basis_data()
      call paw_proj_end()
      call paw_ovlp_end()
      call paw_nonlocal_end()
      call paw_comp_charge_end()
      call paw_mult_end()
c      call paw_vloc_end()
      call end_paw_kin_matrix()
      call end_paw_vloc_matrix()
      call end_paw_ion_matrix()
      call end_paw_overlap_matrix()
      call end_paw_hartree_matrix()
      call end_paw_core_matrix()
      call end_paw_comp_pot_matrix()
      call end_paw_comp_charge_matrix()
      call end_paw_pot_comp()
      call end_paw_pot_hartree()
c      call paw_density_end()
      call paw_xc_end()
      call paw_gaunt_end()

      call D3dB_pfft_end()
      call D3dB_end(1)
      if (control_version().eq.4) call D3dB_end(2)
      call Dne_end()

*                |***************************|
****************** report consumed cputime   **********************
*                |***************************|
      if (lprint) then
         CALL current_second(cpu4)

         T1=CPU2-CPU1
         T2=CPU3-CPU2
         T3=CPU4-CPU3
         T4=CPU4-CPU1
         AV=T2/dble(paw_electron_count())
         write(6,1801)
         write(6,1802)
         write(6,1803) T1
         write(6,1804) T2
         write(6,1805) T3
         write(6,1806) T4
         write(6,1807) AV,paw_electron_count(),linesearch_count()
         if (oprint) then
         write(6,1809) 'Time spent doing               ',
     >              'total','step'
         write(6,1808) '  FFTs                       : ', 
     >                 nwpw_timing(1),
     >                 nwpw_timing(1)/dble(paw_electron_count())
         write(6,1808) '  dot products               : ', 
     >                 nwpw_timing(2),
     >                 nwpw_timing(2)/dble(paw_electron_count())
         write(6,1808) '  geodesic                   : ', 
     >                 nwpw_timing(10),
     >                 nwpw_timing(10)/dble(paw_electron_count())
         write(6,1808) '  exchange correlation       : ', 
     >                 nwpw_timing(4),
     >                 nwpw_timing(4)/dble(paw_electron_count())
         write(6,1808) '  local pseudopotentials     : ', 
     >                 nwpw_timing(5),
     >                 nwpw_timing(5)/dble(paw_electron_count())
         write(6,1808) '  non-local pseudopotentials : ', 
     >                 nwpw_timing(6),
     >                 nwpw_timing(6)/dble(paw_electron_count())
         write(6,1808) '  hartree potentials         : ', 
     >                 nwpw_timing(7),
     >                 nwpw_timing(7)/dble(paw_electron_count())
         write(6,1808) '  structure factors          : ', 
     >                 nwpw_timing(8),
     >                 nwpw_timing(8)/dble(paw_electron_count())
         write(6,1808) '  masking and packing        : ',
     >                 nwpw_timing(9),
     >                 nwpw_timing(9)/dble(paw_electron_count())
         write(6,1808) '  atomic density generation  : ',
     >                 nwpw_timing(21),
     >                 nwpw_timing(21)/dble(paw_electron_count())
         write(6,1808) '  atomic xc matrix elements  : ',
     >                 nwpw_timing(22),
     >                 nwpw_timing(22)/dble(paw_electron_count())
         end if
         write(6,*)
         CALL nwpw_message(4)
      end if 



      call Parallel_Finalize()
      paw_minimizer = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*               NWPW PAW Calculation               *')
 1030 FORMAT(10X,'*  [ (Grassman/Stiefel manifold implementation) ]  *')
 1035 FORMAT(10x,'*      [ NorthWest Chemistry implementation ]      *')
 1040 FORMAT(10X,'*            version #5.10   06/12/02              *')
 1041 FORMAT(10X,'*      This code was developed by Marat Valiev     *')
 1042 FORMAT(10X,'*      and Eric J. Bylaska and was based upon      *')
 1043 FORMAT(10X,'*      algorithms and code developed by the        *')
 1044 FORMAT(10X,'*      group of Prof. John H. Weare                *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1111 FORMAT(/' number of processors used:',I3)
 1112 FORMAT( ' parallel mapping         : slab')
 1113 FORMAT( ' parallel mapping         : hilbert')
 1115 FORMAT(/' options:')
 1118 FORMAT( ' parallel mapping         : hcurve')
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1132 FORMAT(5X,' using fractional occupation')
 1140 FORMAT(/' elements involved in the cluster:')
 1141 FORMAT(5X,I2,': ',A4,'  ion charge:',F4.1,'  core charge:',F4.1)
 1142 FORMAT(12x,' loggrid parameter r0        :',E10.3)
 1143 FORMAT(12x,' loggrid parameter rmax      :',E10.3)
 1144 FORMAT(12x,' loggrid parameter npoints   :',I10)
 1145 FORMAT(12x,' augmentation sphere radius  :',F10.3,
     >           ' (',I5,' npoints',I5,' per task)')
 1146 FORMAT(12x,' compensation sigma          :',F10.3)
 1150 FORMAT(12x,' total number of projectors  :',I10)
 1151 FORMAT(12x,' n_ps (n) l          eig      #projector')
 1152 FORMAT(14X,I3,' (',I1,') ',A,F13.6,I16)
 1153 FORMAT(12X,' local potential used           : ',i2)
 1154 FORMAT(12X,' number of non-local projections: ',i2)
 1155 FORMAT(12X,' semicore corrections included  : ',
     >       F6.3,' (radius) ',F6.3,' (charge)')
 1156 FORMAT(12X,' aperiodic cutoff radius        : ',F6.3)
 1159 FORMAT(/' total charge:',F8.3)
 1160 FORMAT(/' atomic composition:')
 1170 FORMAT(7(5X,A4,':',I3))
 1180 FORMAT(/' position of ions (au):')
 1190 FORMAT(5X, I4, A5  ,' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1191 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - fixed')
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1219 FORMAT(/' number of active electrons: spin up=',F6.2,
     >        '  down=',F6.2,A)
 1220 FORMAT(/' number of active electrons: spin up=',I6,  
     >        '  down=',I6,A)
 1221 FORMAT( ' number of active orbitals : spin up=',I6,  
     >        '  down=',I6,A)
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
 1262 FORMAT(5X,' core integration lmax=',I2)
 1270 FORMAT(/' technical parameters:')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'fictitious mass=',F10.1)
 1290 FORMAT(5X, ' tolerance=',E8.3,' (energy)',E12.3,
     &        ' (density)')
 1291 FORMAT(/' Kohn-Sham scf parameters:')
 1292 FORMAT(5X, ' Kohn-Sham algorithm  = conjugate gradient')
 1293 FORMAT(5X, ' scf algorithm        = ',A)
 1294 FORMAT(5X, ' scf mixing parameter =',F7.4)
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
 1600 FORMAT(/' Total PAW energy   :',E19.10)

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

      paw_minimizer = value
      return
      END
