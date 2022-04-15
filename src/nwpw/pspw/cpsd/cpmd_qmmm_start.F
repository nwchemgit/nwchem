*
* $Id$
*

      logical function cpmd_qmmm_start(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
ccc#include "frac_occ.fh"
#include "errquit.fh"
      
      logical value

      real*8 kb
      parameter (kb=3.16679d-6)

      
*     **** parallel variables ****
      integer  taskid,np
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,n2ft3d

*     **** electronic variables ****
      real*8 icharge
      integer n1(2),n2(2),neall
      real*8  en(2)
      real*8 dipole(3)


*     ***** energy variables ****
      real*8  E(30),cv

      integer eig(2)



*     **** cpmd_qmmm_block ****
      integer nfft3d
      integer npack1
      integer ispin
      integer ne(2),nemax
      integer psi0(2),psi1(2),psi2(2)
      integer dn(2)
      integer Hpsi(2),psir(2)
      integer hml(2),lmd(2)
      integer icount,icount_shift
      real*8  emotion_time_shift,dt,eave,evar,eke
      logical mulliken
      common /cpmd_qmmm_block/ psi0,psi1,psi2,dn,hml,lmd,Hpsi,psir,
     >                         ispin,ne,npack1,nfft3d,nemax,
     >                         icount,icount_shift,
     >                         emotion_time_shift,dt,eave,evar,eke,
     >                         mulliken



*     **** error variables ****
      integer ierr

*     **** local variables ****
      logical verlet,SA,found
      integer ms
      real*8  gx,gy,gz,cx,cy,cz
      real*8  vgx,vgy,vgz,vcx,vcy,vcz
      real*8  ekg,eki0,eki1,sum
      real*8  eke0,eke1
      real*8  EV,pi
      integer i,j,k,ia,n,nn
      integer ii,jj,index,indx
      integer it_in,it_out
      real*8 w,sumall
      real*8 Te_init,Tr_init,Te_new,Tr_new,sa_decay(2),sa_alpha(2)
      integer nfft3,mapping
      parameter (nfft3=32)
      character*50 filename
      character*255 full_filename
  


*     **** external functions ****
      real*8      psp_zv,psp_rc,ewald_rcut,ion_amass
      real*8      ewald_mandelung
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg
      integer     ewald_ncut,ewald_nshl3d
      integer     psp_lmmax,psp_lmax,psp_locp,ion_nkatm
      character*4 ion_atom,ion_aname
      external    psp_zv,psp_rc,ewald_rcut,ion_amass
      external    ewald_mandelung
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg
      external    ewald_ncut,ewald_nshl3d
      external    psp_lmmax,psp_lmax,psp_locp,ion_nkatm
      external    ion_atom,ion_aname


      real*8   control_rti,control_rte,ion_rion
      real*8   ion_vion,ion_com_ke,ion_ke
      real*8   ion_Temperature,ion_com_Temperature
      real*8   pspw_qmmm_Temperature
      external control_rti,control_rte,ion_rion
      external ion_vion,ion_com_ke,ion_ke
      external ion_Temperature,ion_com_Temperature
      external pspw_qmmm_Temperature
      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      logical  control_read,control_move,ion_init,ion_q_FixIon
      external control_read,control_move,ion_init,ion_q_FixIon

      integer  pack_nwave_all
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      integer  pspw_qmmm_nion
      external pack_nwave_all
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      external pspw_qmmm_nion

      character*12 control_boundry
      external     control_boundry

      logical      psp_semicore
      real*8       psp_rcore,psp_ncore,psp_rlocal
      external     psp_semicore
      external     psp_rcore,psp_ncore,psp_rlocal

      logical  control_Nose,control_Mulliken,pspw_qmmm_found
      external control_Nose,control_Mulliken,pspw_qmmm_found

      integer  Nose_Mchain,Nose_Nchain
      external Nose_Mchain,Nose_Nchain

      real*8   control_Nose_Te,Nose_Qe,Nose_Pe,Nose_Ee0
      external control_Nose_Te,Nose_Qe,Nose_Pe,Nose_Ee0

      real*8   control_Nose_Tr,Nose_Qr,Nose_Pr,Nose_Er0
      external control_Nose_Tr,Nose_Qr,Nose_Pr,Nose_Er0
      logical      v_psi_filefind
      external     v_psi_filefind
      real*8   nwpw_timing
      external nwpw_timing

      logical  control_out_of_time,control_new_vpsi
      external control_out_of_time,control_new_vpsi

      logical  control_SA,control_Fei,pspw_SIC,pspw_HFX
      real*8   control_SA_decay
      external control_SA,control_Fei,pspw_SIC,pspw_HFX
      external control_SA_decay
      integer  control_mapping
      external control_mapping

      logical  control_translation,control_rotation
      external control_translation,control_rotation
     



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
         call nwpw_message(1)
         write(6,1110)
      end if
      
      value = control_read(2,rtdb)
      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0
      mapping = control_mapping()


*     **** initialize D3dB data structure ****
      call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      if (control_version().eq.4) 
     >   call D3dB_Init(2,2*ngrid(1),2*ngrid(2),2*ngrid(3),mapping)


*     **** initialize lattice data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_init()
      call Pack_npack(1,npack1)

      call D3dB_pfft_init()


*     ***** allocate psi2, psi1, and psi0 wavefunctions ****
      call psi_get_ne(ispin,ne)
      value = BA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'psi1',psi1(2),psi1(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'psi0',psi0(2),psi0(1))
      if (.not.value) call errquit('out of heap memory',0, MA_ERR)



*     *****  read psi2 wavefunctions ****
      call psi_read(ispin,ne,dcpl_mb(psi2(1)))


*     **** move  wavefunction velocities ****
      if (control_new_vpsi()) then
        call v_psi_delete()
      end if

*     **** generate initial wavefunction velocities if it does not exist ****
      if (.not.v_psi_filefind()) then
        call v_psi_new(ispin,ne)
      end if

*     *****  read psi0 wavefunctions ****
      call v_psi_read(ispin,ne,dcpl_mb(psi1(1)))
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)
      nemax = ne(1)+ne(2)

*     **** allocate other variables *****
      value = BA_alloc_get(mt_dbl,(2*nemax),'eig',eig(2),eig(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*nemax*nemax),'hml',hml(2),hml(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*nemax*nemax),'lmd',lmd(2),lmd(1))

      value = value.and.
     >        BA_alloc_get(mt_dbl,(4*nfft3d),
     >                     'dn',dn(2),dn(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'Hpsi',Hpsi(2),Hpsi(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d*(ne(1)+ne(2)),
     >                     'psir',psir(2),psir(1))
      if (.not.value) call errquit('out of heap memory',0, MA_ERR)



*     **** read ions ****
      value = ion_init(rtdb)

*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)

*     **** allocate psp data structure and read in psedupotentials into it ****
      call psp_init()
      call psp_readall()
      if (psp_semicore(0)) call semicore_check()


*     **** initialize G,mask,ke,and coulomb data structures ****
      call ke_init()
      if (control_version().eq.3) call coulomb_init()
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()
      if (control_version().eq.3) call ewald_init()


*     ******************************
*     **** scaling psi velocity ****
*     ******************************
      call dcopy(2*(ne(1)+ne(2))*npack1,dcpl_mb(psi1(1)),1,
     >                                  dcpl_mb(psi0(1)),1)
      call dscal(2*(ne(1)+ne(2))*npack1,control_rte(),
     >           dcpl_mb(psi1(1)),1)
      eke0 = 0.0d0
      eke1 = 0.0d0
      do i=1,(ne(1)+ne(2))
         call Pack_cc_dot(1,dcpl_mb(psi0(1)+(i-1)*npack1),
     >                      dcpl_mb(psi0(1)+(i-1)*npack1),
     >                     sum)
         eke0 = eke0 + sum
         call Pack_cc_dot(1,dcpl_mb(psi1(1)+(i-1)*npack1),
     >                      dcpl_mb(psi1(1)+(i-1)*npack1),
     >                    sum)
         eke1 = eke1 + sum
      end do
      eke0 = control_fake_mass()*eke0
      eke1 = control_fake_mass()*eke1
      call ion_init_ke(ekg,eki0,eki1)

*     **** Initialize thermostats ****
      if (control_Nose()) then
         call ke_ave(ispin,ne,dcpl_mb(psi2(1)),w)
         call Nose_Init((ne(1)+ne(2)),w)
      end if

*     **** Initialize simulated annealing ****
      SA       = .false.
      Te_init  = 0.0d0
      Tr_init  = 0.0d0
      sa_alpha(1) = 1.0d0
      sa_alpha(2) = 1.0d0
      if (control_SA()) then
         if (control_Nose()) then
            SA          = .true.
            sa_decay(1) = control_SA_decay(1)
            sa_decay(2) = control_SA_decay(2)
            Te_init     = control_Nose_Te()
            Tr_init     = control_Nose_Tr()
         else
            dt = control_time_step()
            SA          = .false.
            sa_decay(1) = control_SA_decay(1)
            sa_decay(2) = control_SA_decay(2)
            sa_alpha(1) = dexp( -(dt/control_SA_decay(1)) ) 
            sa_alpha(2) = dexp( -(dt/control_SA_decay(2)) ) 
         end if
      end if


*     **** initialize QM/MM ****
      call pspw_qmmm_init(rtdb)

*     **** initialize dplot ****
      call dplot_iteration_init()

*     **** initialize frac_occ data structure ****
ccc      call frac_occ_init(rtdb,ispin,ne)

*     **** initialize SIC and HFX ****
      call pspw_init_SIC(rtdb,ne)
      call pspw_init_HFX(rtdb,ispin,ne)







*                |**************************|
******************   summary of input data  **********************
*                |**************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      call center_v_geom(vcx,vcy,vcz)
      call center_v_mass(vgx,vgy,vgz)
      mulliken = control_Mulliken()

      if (taskid.eq.MASTER) then
         write(6,1111) np
         if (mapping.eq.1) write(6,1112)
         if (mapping.eq.2) write(6,1113)
         if (mapping.eq.3) write(6,1118)

         write(6,1115)
         write(6,1121) control_boundry(),control_version()
         if (ispin.eq.1) write(6,1130) 'restricted'
         if (ispin.eq.2) write(6,1130) 'unrestricted'
         IF (control_gga().eq.-1) THEN
            write(6,1131) 'off'
         ELSE IF (control_gga().eq.0) THEN
            write(6,1131) 'LDA (Vosko et al) parameterization'
         ELSE IF (control_gga().eq.1) THEN
            write(6,1131)
     >      'PBE96 (White and Bird) parameterization'
         ELSE IF (control_gga().eq.2) THEN
            write(6,1131)
     >      'BLYP (White and Bird) parameterization'
         ELSE IF (control_gga().eq.3) THEN
            write(6,1131)
     >      'PBE0 (White and Bird) parameterization'
         ELSE IF (control_gga().eq.4) THEN
            write(6,1131) 'Hartree-Fock'
         ELSE
            write(6,1131) 'unknown parameterization'
            call errquit('bad exchange_correlation',0, INPUT_ERR)
         END IF
c         if (fractional) write(6,1132)
         call pspw_print_SIC(6)
         call pspw_print_HFX(6)
         write(6,1140)
         do ia = 1,ion_nkatm()
           write(6,1150) ia,ion_atom(ia),
     >                    psp_zv(ia),psp_lmax(ia)
           write(6,1152) psp_lmax(ia)
           write(6,1153) psp_locp(ia)
           write(6,1154) psp_lmmax(ia)
           if (control_version().eq.4) write(6,1156) psp_rlocal(ia)
           if (psp_semicore(ia))
     >         write(6,1155) psp_rcore(ia),psp_ncore(ia)
           write(6,1151) (psp_rc(i,ia),i=0,psp_lmax(ia))
         end do

         icharge = -(ne(1)+ne(ispin))
         en(1)     = ne(1)
         en(ispin) = ne(ispin)
c         if (fractional) then
c            icharge = 0.0d0
c            do ms=1,ispin
c            en(ms) =0.0
c            do i=n1(ms),n2(ms)
c              icharge = icharge - (3-ispin)*dbl_mb(fweight(1)+i-1)
c              en(ms) = en(ms) + dbl_mb(fweight(1)+i-1)
c            end do
c            end do
c         end if

         do ia=1,ion_nkatm()
           icharge = icharge + ion_natm(ia)*psp_zv(ia)
         end do
         write(6,1159) icharge

         write(6,1160)
         write(6,1170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
         write(6,1180)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(6,1191) I,ion_aname(I),
     >                    (ion_rion(K,I),K=1,3),ion_amass(i)/1822.89d0
           else
           write(6,1190) I,ion_aname(I),
     >                    (ion_rion(K,I),K=1,3),ion_amass(i)/1822.89d0
           end if
         end do
         write(6,1200) cx,cy,cz
         write(6,1210) gx,gy,gz

         call pspw_qmmm_Print(6)

         write(6,1181)
         write(6,1192) (I,ion_aname(I),
     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
         write(6,1200) vcx,vcy,vcz
         write(6,1210) vgx,vgy,vgz

         call pspw_qmmm_Print2(6)

c         if (fractional) then
c           write(6,1219) en(1),en(ispin),' (   fractional)'
c           write(6,1221) ne(1),ne(ispin),' (fourier space)'
c         else
           write(6,1220) ne(1),ne(ispin),' (fourier space)'
           write(6,1221) ne(1),ne(ispin),' (fourier space)'
c         end if
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
         if (.not.control_translation()) write(6,1271)
         if (.not.control_rotation())    write(6,1272)
         write(6,1280) control_time_step(),control_fake_mass()
         write(6,1290) control_rte(),control_rti()
         call ion_scaling_atoms_print(6)
         write(6,1222) eke0,eki0,ekg
         write(6,1223) eke1,eki1
         write(6,1224) (eke1-eke0),(eki1-eki0)
         if (control_Nose()) then
           write(6,1295)
           do i=1,Nose_Mchain()
             write(6,1297) i,control_Nose_Te(),Nose_Qe(i),
     >                     Nose_Pe(i),Nose_Ee0(i)
           end do
           do i=1,Nose_Nchain()
             write(6,1298) i,control_Nose_Tr(),Nose_Qr(i),
     >                     Nose_Pr(i),Nose_Er0(i)
           end do
         else
           write(6,1294)
         end if
        if (control_SA()) then
           write(6,1296) sa_decay(1),sa_decay(2)
         end if


         if (mulliken) write(6,1299)
         write(6,1300)
         write(6,1305)
         call util_flush(6)
      end if


*                |***************************|
******************     start iterations      **********************
*                |***************************|
*     **** open xyz and MOTION file ****
      call xyz_init()          ! unit=18
      call MOTION_init(rtdb)   ! unit=19

*     *** fei io ****
      call fei_init()




*     ************************************
*     **** open up other MOTION files ****
*     ************************************


*     **** open EMOTION file ****
      value = btdb_cget(rtdb,'cpmd:emotion_filename',1,filename)
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
           E(25) = E(25) + sum                          !*** take care of running sums ***
           E(26) = E(26) + sum*sum
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


*     **** open EIGMOTION file ****
      if (mulliken) then
      value = btdb_cget(rtdb,'cpmd:eigmotion_filename',1,filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)
      if (taskid.eq.MASTER) 
     >   open(unit=32,file=full_filename,form='formatted')
      end if

*     **** open HMOTION file ****
      if (mulliken) then
      value = btdb_cget(rtdb,'cpmd:hmotion_filename',1,filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)
      if (taskid.eq.MASTER) 
     >   open(unit=34,file=full_filename,form='formatted')
      end if

*     **** open OMOTION file ****
      if (mulliken) call Orb_Init(rtdb,ispin,ne) !unit=33

*     **** write initial position to xyz data ****
      call xyz_write()

      eke = 0.0d0
      dt  = control_time_step()




      cpmd_qmmm_start = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*     Car-Parrinello microcluster calculation      *')
 1030 FORMAT(10X,'*      [   extended Lagrangian molecular   ]       *')
 1031 FORMAT(10X,'*      [        dynamics simulation        ]       *')
 1035 FORMAT(10x,'*      [ NorthWest Chemistry implementation ]      *')
 1040 FORMAT(10X,'*            version #5.00   06/01/99              *')
 1041 FORMAT(10X,'*    This code was developed by Eric J. Bylaska,   *')
 1042 FORMAT(10X,'*    and was based upon algorithms and code        *')
 1043 FORMAT(10X,'*    developed by the group of Prof. John H. Weare *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1111 FORMAT(/' number of processors used:',I3)
 1112 FORMAT( ' parallel mapping         :     1d-slab')
 1113 FORMAT( ' parallel mapping         :  2d-hilbert')
 1115 FORMAT(/' options:')
 1118 FORMAT( ' parallel mapping         :   2d-hcurve')
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1132 FORMAT(5X,' using fractional occupation')
 1140 FORMAT(/' elements involved in the cluster:')
 1150 FORMAT(5X,I2,': ',A4,'  core charge:',F4.1,'  lmax=',I1)
 1151 FORMAT(5X,'        cutoff =',4F8.3)
 1152 FORMAT(12X,' highest angular component      : ',i2)
 1153 FORMAT(12X,' local potential used           : ',i2)
 1154 FORMAT(12X,' number of non-local projections: ',i2)
 1155 FORMAT(12X,' semicore corrections included  : ',
     >       F6.3,' (radius) ',F6.3,' (charge)')
 1156 FORMAT(12X,' aperiodic cutoff radius        : ',F6.3)
 1159 FORMAT(/' total charge=',F8.3)
 1160 FORMAT(/' atomic composition:')
 1170 FORMAT(7(5X,A4,':',I3))
 1180 FORMAT(/' initial position of ions:')
 1181 FORMAT(/' initial velocity of ions:')
 1190 FORMAT(5X, I4, A4  ,' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1191 FORMAT(5X, I4, A4  ,' (',3F11.5,
     >       ' ) - atomic mass= ',F6.3,' - fixed')
 1192 FORMAT(5X, I4, A4  ,' (',3F11.5,' )')
 1200 FORMAT(5X,'  G.C. ',' (',3F11.5,' )')
 1210 FORMAT(5X,' C.O.M.',' (',3F11.5,' )')
 1219 FORMAT(/' number of electrons: spin up=',F6.2,'  down=',F6.2,A)
 1220 FORMAT(/' number of electrons: spin up=',I6,  '  down=',I6,A)
 1221 FORMAT( ' number of orbitals : spin up=',I6,  '  down=',I6,A)
 1222 format(5x,' initial kinetic energy: ',e12.5,' (psi)', 2x,
     >                                      e12.5,' (ion)',/50x,
     >                                      e12.5,' (c.o.m.)')
 1223 format(5x,' after scaling:          ',e12.5,' (psi)', 2x,
     >                                      e12.5,' (ion)')
 1224 format(5x,' increased energy:       ',e12.5,' (psi)', 2x,
     >                                      e12.5,' (ion)')
 1226 format(/' final kinetic energy:  ',e12.5,' (psi)', 2x,
     >                                      e12.5,' (ion)',/44x,
     >                                      e12.5,' (c.o.m.)')
 1230 FORMAT(/' supercell:')
 1231 FORMAT(5x,' volume : ',F10.1)
 1241 FORMAT(5x,' lattice:    a1=<',3f8.3,' >')
 1242 FORMAT(5x,'             a2=<',3f8.3,' >')
 1243 FORMAT(5x,'             a3=<',3f8.3,' >')
 1244 FORMAT(5x,' reciprocal: b1=<',3f8.3,' >')
 1245 FORMAT(5x,'             b2=<',3f8.3,' >')
 1246 FORMAT(5x,'             b3=<',3f8.3,' >')

 1250 FORMAT(5X,' density cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
 1251 FORMAT(5X,' wavefnc cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
  
 1260 FORMAT(5X,' Ewald summation: cut radius=',F8.2,'  and',I3)
 1261 FORMAT(5X,'                   madelung=',f14.8)
 1270 FORMAT(/' technical parameters:')
 1271 FORMAT(5x, ' translation constrained')
 1272 FORMAT(5x, ' rotation constrained')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'fictitious mass=',F10.1)
 1290 FORMAT(5X, ' cooling/heatting rates: ',e12.5,' (psi)',2x,
     >                                       e12.5,' (ion)')
 
 1294 format(/' Constant Energy Simulation                     ')
 1295 format(/' Nose-Hoover Simulation - thermostat parameters:')
 1296 format(5x, 'SA decay rates  =',e10.3,' (elc)',e10.3,' (ion)')
 1297 format(5x, 'link = ',I3,
     > ' Te =',f8.2,' Qe =',e10.3,' 2*pi/we=',e10.3,' Ee0=',e10.3)
 1298 format(5x, 'link = ',I3,
     > ' Tr =',f8.2,' Qr =',e10.3,' 2*pi/wr=',e10.3,' Er0=',e10.3)
 1299 format(//' Mulliken analysis output generated            ')
 1300 FORMAT(//)
 1305 FORMAT(10X,'============ Car-Parrinello iteration ==============')
 1309 FORMAT(I8,2E19.10,2E14.5,2F9.1)
 1310 FORMAT(I8,2E19.10,2E14.5,F14.2)
 1311 format(100e19.10)
 1312 format(e14.6,i3)
 1313 format(3i4)
 1320 FORMAT(' number of electrons: spin up=',F11.5,'  down=',F11.5,A)
 1330 FORMAT(/' comparison between hamiltonian and lambda matrix')
 1331 FORMAT(/' Diagonal Elements of Hamiltonian matrix')
 1340 FORMAT(I3,2I3,' H=',E16.7,', L=',E16.7,', H-L=',E16.7)
 1341 FORMAT(I3,2I3,' H=',E16.6)
 1350 FORMAT(/' orthonormality')
 1360 FORMAT(I3,2I3,E18.7)
 1370 FORMAT(I3)
 1380 FORMAT(' ''',a,'''',I4)
 1390 FORMAT(I3)
 1400 FORMAT(I3,3E18.8/3X,3E18.8)
 1410 FORMAT(10X,'=============  summary of results  =================')
 1420 FORMAT(/' final position of ions:')
 1421 FORMAT(/' final velocity of ions:')
 1430 FORMAT(/' total     energy    :',E19.10,' (',E15.5,'/ion)')
 1440 FORMAT( ' total orbital energy:',E19.10,' (',E15.5,'/electron)')
 1450 FORMAT( ' hartree   energy    :',E19.10,' (',E15.5,'/electron)')
 1455 FORMAT( ' SIC-hartree energy  :',E19.10,' (',E15.5,'/electron)')
 1456 FORMAT( ' SIC-exc-corr energy :',E19.10,' (',E15.5,'/electron)')
 1457 FORMAT( ' HF exchange energy  :',E19.10,' (',E15.5,'/electron)')
 1460 FORMAT( ' exc-corr  energy    :',E19.10,' (',E15.5,'/electron)')
 1470 FORMAT( ' ion-ion   energy    :',E19.10,' (',E15.5,'/ion)')
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
 1499 FORMAT( ' K.S. SIC-hartree energy  :',E19.10,
     >        ' (',E15.5,'/electron)')
 1501 FORMAT( ' K.S. SIC-exc-corr energy :',E19.10,
     >        ' (',E15.5,'/electron)')
 
 1500 FORMAT(/' orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1511 FORMAT(2(E18.7,' (',F8.3,'eV) occ=',F5.3))
 1600 FORMAT(/' Total PSPW energy   :',E19.10)

 1700 FORMAT(/' QM/MM-pol-vib/CAV Energies')
 1701 FORMAT( ' --------------------------')
 1702 FORMAT( ' QM/MM energy           :',E19.10)
 1703 FORMAT( ' MM/MM energy           :',E19.10)
 1704 FORMAT( ' MM Polarization energy :',E19.10)
 1705 FORMAT( ' MM Vibration energy    :',E19.10)
 1706 FORMAT( ' (QM+MM)/Cavity energy  :',E19.10)

 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(6,9010) ierr
      call Parallel_Finalize()

      cpmd_qmmm_start = value
      return
      END



