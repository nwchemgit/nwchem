*
* $Id$
*

***********************************************************************
*                      ofpw  (MPI code)                               *
*                                                                     *
*     This is a developing orbital free dft parallel code for NWCHEM  *
*                                                                     *
***********************************************************************

      logical function ofpw_minimizer(rtdb,flag)
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
      real*8  a,b,c,alpha,beta,gamma
    

*     ***** energy variables ****
      integer ispin,ne(2),ms
      real*8     sum0,sum1,sum3,scale,tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion

      real*8  E(10),en(2)
      real*8  dipole(3)
      real*8  stress(3,3),lstress(6)

*     **** gradient variables ****
      integer fion(2)

*     **** error variables ****
      logical value,newpsi
      integer ierr

*     **** local variables ****
      logical oprint,lprint,hprint,qmmm
      real*8  gx,gy,gz,cx,cy,cz
      real*8  EV,pi
      real*8  icharge
      integer i,k,ia,nion,minimizer,mapping,mapping1d
  
      integer neall,npack1
      real*8 e_lj,e_q,e_spring


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
      real*8   ofpw_cgsd_energy
      external ofpw_cgsd_energy
      logical  control_Mulliken
      external control_Mulliken
      real*8   ion_TotalCharge
      external ion_TotalCharge
      logical   control_check_charge_multiplicity
      external  control_check_charge_multiplicity
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
      integer  control_mapping,control_mapping1d,control_np_orbital
      external control_mapping,control_mapping1d,control_np_orbital
      integer  control_ks_maxit_orb,control_ks_maxit_orbs
      external control_ks_maxit_orb,control_ks_maxit_orbs

      logical pspw_bqext
      external pspw_bqext

*                            |************|
*****************************|  PROLOGUE  |****************************
*                            |************|

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
         write(luout,1042)
         write(luout,1043)
         write(luout,1044)
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
      call Pack_Init()
      call D3dB_pfft_init()


*     **** read ions ****
      value = ion_init(rtdb)
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

*     **** allocate psp data structure and read in psedupotentials into it ****
      call psp_init()
      call psp_readall()
      if (psp_semicore(0)) call semicore_check()


*     **** initialize G,mask,ke,coulomb,and ewald data structures ****
      !call ke_init()
      if (control_version().eq.3) call coulomb_init()
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()
      if (control_version().eq.3) call ewald_init()

*     **** read in wavefunctions and initialize psi ****
c      chi(r) = ifft(chi(G)
c      loop
c         chi(r) = ifft(chi(G))
c         rho(r) = chi(r)*chi(r)
c         rho(G) = fft(rho(r))
c
c         Hchi(r) = 0
c
c         *** local psp term ****
c         Vloc(G) = Struct(I,G)*vlocal(G)
c         Vloc(r) = ifft(Vloc(G)
c         Hchi(r) += Vloc(r)*chi(r)
c
c         *** hartree term ***
c         Vh(G) = 4pi/G2 * rho(G)
c         Vh(r) = ifft(Vh(G))
c         Hchi(r) += Vh(r)*chi(r)
c
c         *** xc term ***
c         vxc(r) = VXC(rho(r))
c         Hchi(r) += vxc(r)*chi(r)
c
c         **** KE term ****
c         vtfvw(r) = VTFVW(rho(r))
c         Hchi(r) += vtfvw(r)*chi(r)
c        
c         vwgc(G) = Kernal(G)*rho(G)
c         Vwgc(r) = ifft(Vwgc(G))
c         Hchi(r) += Vwgc(r)*chi(r)
c
c
c         chi(r) = chi(r) +alpha*Hchi(r)
c       end loop
         


*     **** Initialize 1d mappings for ne spaces ****



*     **** electron and geodesic data structures ****

*     **** initialize QM/MM ****
      call pspw_init_APC(rtdb)
      call pspw_qmmm_init(rtdb)
      call pspw_charge_init(rtdb)

*     **** initialize SIC and HFX ****

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
         write(luout,1115)

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

         ELSE IF (control_gga().eq.110) THEN
            write(luout,1131) 
     >      'PBE0 (White and Bird) parameterization'
         ELSE IF (control_gga().eq.112) THEN
            write(luout,1131) 
     >      'revPBE0 (White and Bird) parameterization'

         ELSE IF (control_gga().eq.200) THEN
            write(luout,1131) 'Hartree-Fock'
         ELSE
            write(luout,1131) 'unknown parameterization'
            call errquit('bad exchange_correlation',0, INPUT_ERR)
         END IF

         write(luout,1140)
         do ia = 1,ion_nkatm()
           write(luout,1150) ia,ion_atom(ia),
     >                    psp_zv(ia),psp_lmax(ia)
           comment = psp_comment(ia)
           i = inp_strlen(comment)
           write(luout,1157) comment(1:i)
           write(luout,1158) psp_psp_type(ia)
           write(luout,1152) psp_lmax(ia)
           write(luout,1153) psp_locp(ia)
           write(luout,1154) psp_nprj(ia)
           if (control_version().eq.4) write(luout,1156) psp_rlocal(ia)
           if (psp_semicore(ia))
     >         write(luout,1155) psp_rcore(ia),psp_ncore(ia)
           write(luout,1151) (psp_rc(i,ia),i=0,psp_lmax(ia))
           write(luout,*)
         end do

         !icharge = -(psi_ne(1)+psi_ne(psi_ispin()))
         icharge = 0.0d0
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

         call pspw_charge_Print(luout)
         end if

         !write(luout,1220) psi_ne(1),psi_ne(psi_ispin()),' (fourier space)'
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

         write(luout,1300)
c         call ion_Print_neighborlist()
         call util_flush(luout)
      end if

*                |***************************|
******************     call CG minimizer     **********************
*                |***************************|
      if (taskid.eq.MASTER) call current_second(cpu2)


*     **** calculate energy ****
      if (flag.eq.-1) then
        !EV= ofpw_cgsd_noit_energy()
        EV= 0.0d0
      else
        EV= ofpw_cgsd_energy()
      end if


*     **** extra energy output for QA test ****
      if (lprint) write(6,1600) EV

      
*     **** calculate the dipole ***
      dipole(1) = 0.0d0
      dipole(2) = 0.0d0
      dipole(3) = 0.0d0
      !if (flag.gt.-1) call rho_dipole(dipole)

*     ****  calculate gradient ***
      if (flag.gt.0) then
      nion = ion_nion()
      if (pspw_charge_found().and.
     >    (.not.pspw_bqext())) nion = nion + pspw_charge_nion()

      value = BA_push_get(mt_dbl,(3*nion),
     >                       'fion',fion(2),fion(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
     
      !call cgsd_energy_gradient(dbl_mb(fion(1)))
      call dcopy(3*nion,0.0d0,0, dbl_mb(fion(1)),1)
      call dscal(3*nion,(-1.0d0),dbl_mb(fion(1)),1)
      end if

*     **** calculate the stress tensor ****
      call dcopy(9,0.0d0,0,stress,1)
      call dcopy(6,0.0d0,0,lstress,1)
      if (flag.eq.3) then
         !call psp_stress_init()
         !call psp_stress_readall()
         !call cgsd_energy_stress(stress,lstress)
         !call psp_stress_end()
      end if


*     *************************************************************
*     **** output energy, dipole, and gradient to rtdb for use ****
*     **** by task_energy and task_gradient                    ****
*     *************************************************************
      if (flag.gt.-1) then
      value = btdb_put(rtdb,'ofpw:energy',mt_dbl,1,EV)
      value = value.and.
     >        btdb_put(rtdb,'ofpw:dipole',mt_dbl,
     >                 3,dipole)
      end if
      if (flag.gt.0) then
        value = value.and.
     >        btdb_put(rtdb,'ofpw:gradient',mt_dbl,
     >                 3*nion,dbl_mb(fion(1)))
        value = value.and.BA_pop_stack(fion(2))
      end if
      if (flag.eq.3) then
        value = value.and.
     >        btdb_put(rtdb,'ofpw:stress',mt_dbl,
     >                 9,stress)
        value = value.and.
     >        btdb_put(rtdb,'ofpw:lstress',mt_dbl,
     >                 6,lstress)
      end if
      if (.not. value) call errquit('ofpw: error writing rtdb',0,
     &       RTDB_ERR)

      if (taskid.eq.MASTER) call current_second(cpu3)

*                |***************************|
******************         Epilogue          **********************
*                |***************************|


*     **** calculate Mulliken Populations ***

*     **** write wavefunctions to file and finalize psi ****

*     **** write geometry to rtdb ****
      call pspw_charge_write(rtdb)
      call ion_write(rtdb)


*     **** deallocate heap memory ****
      if (control_version().eq.3) call ewald_end()
      call strfac_end()
      if (control_version().eq.3) call coulomb_end()
      if (control_version().eq.4) call coulomb2_end()
      !call ke_end()
      call mask_end()
      call Pack_end()
      call G_end()
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
      

*                |***************************|
****************** report consumed cputime   **********************
*                |***************************|
      if (lprint) then
         CALL current_second(cpu4)

         T1=CPU2-CPU1
         T2=CPU3-CPU2
         T3=CPU4-CPU3
         T4=CPU4-CPU1
         AV=T2/dble(1)
         AV= 0.0
         write(luout,1801)
         write(luout,1802)
         write(luout,1803) T1
         write(luout,1804) T2
         write(luout,1805) T3
         write(luout,1806) T4
         write(luout,1807) AV,1,linesearch_count()

         !call nwpw_timing_print_final(oprint,1)
         write(luout,*)
         CALL nwpw_MESSAGE(4)
      end if 



      call Parallel_Finalize()
      ofpw_minimizer = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*        NWPW Orbital Free DFT Calculation         *')
 1030 FORMAT(10X,'*         [ (unknown implementation) ]             *')
 1035 FORMAT(10x,'*      [ NorthWest Chemistry implementation ]      *')
 1040 FORMAT(10X,'*            version #1.00   06/30/07              *')
 1041 FORMAT(10X,'*    This code was developed by Niri Govind, and   *')
 1042 FORMAT(10X,'*    Eric J. Bylaska and is based upon the         *')
 1043 FORMAT(10X,'*    algorithms and code developed by the group of *')
 1044 FORMAT(10X,'*    Prof. Emily A. Carter.                        *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1111 FORMAT(/' number of processors used:',I10)
 1112 FORMAT( ' parallel mapping         :      slab')
 1113 FORMAT( ' parallel mapping         :   hilbert')
 1114 FORMAT( ' parallel mapping         :  balanced')
 1115 FORMAT(/' options:')
 1116 FORMAT( ' parallel mapping         : not balanced')
 1117 FORMAT( ' processor grid           :',I4,' x',I4)
 1118 FORMAT( ' parallel mapping         :    hcurve')
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
 1180 FORMAT(/' position of ions (au):')
 1190 FORMAT(5X, I4, A5  ,' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1191 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - fixed')
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1219 FORMAT(/' number of electrons: spin up=',F6.2,
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
 1270 FORMAT(/' technical parameters:')
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

      ofpw_minimizer = value
      return
      END
