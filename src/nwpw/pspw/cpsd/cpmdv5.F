***********************************************************************
*                      cpmdv5-mpi  (MPI code)                         *
*                                                                     *
*     This is a developing cpsdv5 parallel code for NWChem            *
*       + mpi message passing library used                            *
*       + ngp is used instead of nfft in this proceudure              *
*       + error checking is based on aimd.h parameters                *
*         then control file                                           *
*       + my own slap-decomposed parallel 3d-FFT(real->complex) used  *
*                                                                     *
*                                                                     *
***********************************************************************

      logical function cpmdv5(rtdb)
*
* $Id: cpmdv5.F,v 1.6 2001-11-26 21:34:22 bylaska Exp $
*
      implicit none
      integer rtdb

#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      
      logical value

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
*    10 - energy_time
*    11 - density_time
*    12 - allocate/deallocate_time
      real*8 times(50)
      common / timing / times



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
      integer ngrid(3),nwave,nfft3d,n2ft3d
      integer npack1

*     **** electronic variables ****
      integer ispin,icharge
      integer ne(2),n1(2),n2(2),nemax
      real*8  en(2)
      real*8 dipole(3)

*     complex*16 psi1(nfft3d,nemax)
*     complex*16 psi2(nfft3d,nemax)
*     real*8     dn(n2ft3d,2)
*     complex*16 Hpsi(nfft3d,nemax)
*     complex*16 psir(nfft3d,nemax)
      integer psi0(2),psi1(2),psi2(2)
      integer dn(2)
      integer Hpsi(2),psir(2)
    

*     ***** energy variables ****
      real*8  E(20),eke

*     real*8  eig(2*nemax)
*     real*8  hml(2*nemax*nemax)
*     real*8  lmd(2*nemax*nemax)
      integer eig(2),hml(2),lmd(2)




*     **** error variables ****
      integer ierr

*     **** local variables ****
      logical verlet,mulliken
      integer ms
      real*8  gx,gy,gz,cx,cy,cz
      real*8  vgx,vgy,vgz,vcx,vcy,vcz
      real*8  ekg,eki0,eki1,sum
      real*8  eke0,eke1
      real*8  EV,pi,dt
      integer i,j,k,ia,n,nn
      integer ii,jj,index
      integer icount,it_in,it_out
      real*8 w,sumall
      integer nfft3
      parameter (nfft3=32)
      character*30 filename
      character*255 full_filename
  


*     **** external functions ****
      real*8      psp_zv,psp_rc,ewald_rcut,ion_amass
      real*8      ewald_mandelung
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg
      integer     ewald_ncut,ewald_nshl3d,psp_npsp
      integer     psp_lmmax,psp_lmax,psp_locp,ion_nkatm
      character*2 psp_atom,ion_aname
      external    psp_zv,psp_rc,ewald_rcut,ion_amass
      external    ewald_mandelung
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg
      external    ewald_ncut,ewald_nshl3d,psp_npsp
      external    psp_lmmax,psp_lmax,psp_locp,ion_nkatm
      external    psp_atom,ion_aname


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

      real*8   control_Nose_Te,Nose_Qe,control_Nose_Pe,Nose_Ee0
      external control_Nose_Te,Nose_Qe,control_Nose_Pe,Nose_Ee0

      real*8   control_Nose_Tr,Nose_Qr,control_Nose_Pr,Nose_Er0
      external control_Nose_Tr,Nose_Qr,control_Nose_Pr,Nose_Er0


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
      
      value = control_read(2,rtdb)
      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0


      ierr = 0 
      if ((np.gt.(control_ngrid(3)))) then
        if (taskid.eq.MASTER)
     >     write(*,*) 'Error: the number of processors must ',
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
      if (control_version().eq.4) 
     >   call D3dB_Init(2,2*ngrid(1),2*ngrid(2),2*ngrid(3))


*     **** initialize lattice data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_init()
      call Pack_npack(1,npack1)
*     ***** allocate psi2, psi1, and psi0 wavefunctions ****
      call psi_get_ne(ispin,ne)
      value = MA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'psi1',psi1(2),psi1(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'psi0',psi0(2),psi0(1))
      if (.not.value) call errquit('out of heap memory',0)



*     *****  read psi2 and psi0 wavefunctions ****
      call psi_read(ispin,ne,dcpl_mb(psi2(1)))
      call v_psi_read(ispin,ne,dcpl_mb(psi1(1)))
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)
      nemax = ne(1)+ne(2)

*     **** allocate other variables *****
      value = MA_alloc_get(mt_dbl,(2*nemax),'eig',eig(2),eig(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(2*nemax*nemax),'hml',hml(2),hml(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(2*nemax*nemax),'lmd',lmd(2),lmd(1))

      value = value.and.
     >        MA_alloc_get(mt_dbl,(4*nfft3d),
     >                     'dn',dn(2),dn(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'Hpsi',Hpsi(2),Hpsi(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,nfft3d*(ne(1)+ne(2)),
     >                     'psir',psir(2),psir(1))
      if (.not.value) call errquit('out of heap memory',0)



*     **** read ions ****
      value = ion_init(rtdb)

*     **** allocate psp data structure and read in psedupotentials into it ****
      call psp_init()
      call psp_readall()


*     **** initialize G,mask,ke,and coulomb data structures ****
      call ke_init()
      if (control_version().eq.3) call coulomb_init()
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()

      if (control_version().eq.3) call ewald_init()
      if (psp_semicore(0)) call semicore_check()

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
         call Nose_Init(w)
      end if

*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)

*     **** initialize QM/MM ****
      call pspw_qmmm_init(rtdb)




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
         write(6,1115)
         write(6,1121) control_boundry(),control_version()
         if (ispin.eq.1) write(6,1130) 'restricted'
         if (ispin.eq.2) write(6,1130) 'unrestricted'
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
           write(6,1150) ia,psp_atom(ia),
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
         do ia=1,ion_nkatm()
           icharge = icharge + ion_natm(ia)*psp_zv(ia)
         end do
         write(6,1159) icharge

         write(6,1160)
         write(6,1170) (psp_atom(K),ion_natm(K),K=1,psp_npsp())
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


         write(6,1220) ne(1),ne(ispin),' (fourier space)'
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
         write(6,1290) control_rte(),control_rti()
         write(6,1222) eke0,eki0,ekg
         write(6,1223) eke1,eki1
         write(6,1224) (eke1-eke0),(eki1-eki0)
         if (control_Nose()) then
           write(6,1295)
           write(6,1297) control_Nose_Te(),Nose_Qe(),
     >                   control_Nose_Pe(),Nose_Ee0()
           write(6,1298) control_Nose_Tr(),Nose_Qr(),
     >                   control_Nose_Pr(),Nose_Er0()
         else
           write(6,1294)
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
      call xyz_init()
      call MOTION_init(rtdb)


*     **** open up other MOTION files ****
      value = rtdb_cget(rtdb,'cpmd:emotion_filename',1,filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)
      if (taskid.eq.MASTER) 
     >   open(unit=31,file=full_filename,form='formatted')

*     **** open EIGMOTION file ****
      if (mulliken) then
      value = rtdb_cget(rtdb,'cpmd:eigmotion_filename',1,filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)
      if (taskid.eq.MASTER) 
     >   open(unit=32,file=full_filename,form='formatted')
      end if

*     **** open HMOTION file ****
      if (mulliken) then
      value = rtdb_cget(rtdb,'cpmd:hmotion_filename',1,filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)
      if (taskid.eq.MASTER) 
     >   open(unit=34,file=full_filename,form='formatted')
      end if

*     **** open OMOTION file ****
      if (mulliken) call Orb_Init(rtdb,ispin,ne)

*     **** write initial position to xyz data ****
      call xyz_write()

*     ***** first step using velocity ****
      verlet = .false.
      call inner_loop_md(verlet,ispin,ne,
     >             npack1,nfft3d,nemax,
     >             dcpl_mb(psi0(1)),
     >             dcpl_mb(psi1(1)),
     >             dcpl_mb(psi2(1)),
     >             dbl_mb(dn(1)),
     >             1,E,
     >             dbl_mb(hml(1)),dbl_mb(lmd(1)),
     >             dcpl_mb(Hpsi(1)),dcpl_mb(psir(1)))


      if (taskid.eq.MASTER) call current_second(cpu2)
      if (taskid.eq.MASTER) CALL MESSAGE(2)
      it_in  = control_it_in()
      it_out = control_it_out()
      icount = 0
      verlet = .true.
      eke    = 0.0d0
      if (it_out.lt.1) goto 2



      dt = control_time_step()

   1  continue
         icount = icount + 1
         call inner_loop_md(verlet,ispin,ne,
     >             npack1,nfft3d,nemax,
     >             dcpl_mb(psi0(1)),
     >             dcpl_mb(psi1(1)),
     >             dcpl_mb(psi2(1)),
     >             dbl_mb(dn(1)),
     >             it_in,E,
     >             dbl_mb(hml(1)),dbl_mb(lmd(1)),
     >             dcpl_mb(Hpsi(1)), dcpl_mb(psir(1)))

         eke = eke + E(3)

         if (taskid.eq.MASTER) then 

           write(6,1310) icount*it_in,E(1),E(2),E(3),E(4)
           call util_flush(6)

*          **** write out EMOTION data ****
           if (control_Nose()) then
             write(31,1311) icount*it_in*dt,
     >                    e(1),e(2),e(3),e(4),e(5),e(6),
     >                    e(7),e(8),e(9),e(10)
           else
             write(31,1311) icount*it_in*dt,
     >                    e(1),e(2),e(3),e(4),e(5),e(6),
     >                    e(7),e(8)
           end if
           call util_flush(31)

*          **** write out EIGMOTION data -diagonal hml matrix ****
           if (mulliken) then
           write(32,1311) icount*it_in*dt,
     >       (( dbl_mb(hml(1)+ii-1+(ii-1)*ne(1)+(ms-1)*ne(1)*ne(1)), 
     >         ii=1,ne(ms)),ms=1,ispin)
           call util_flush(32)
           end if

*          **** write out HMOTION data - hml matrix ****
           if (mulliken) then
           write(34,1312) icount*it_in*dt,ispin
           do ms=1,ispin
             write(34,1313) ms,ne(ms),ne(ms)
             do ii=1,ne(ms)
               write(34,1311) 
     >         (dbl_mb(hml(1)+ii-1+(jj-1)*ne(1)+(ms-1)*ne(1)*ne(1)), 
     >          jj=1,ne(ms))
             end do
           end do
           call util_flush(34)
           end if

         end if


*        **** write xyz data ****
         call xyz_write()
         call MOTION_write((icount*it_in*control_time_step()))

*        **** write OMOTION data ****
         if (mulliken) call Orb_Write(dcpl_mb(psi1(1)))

      if (icount.lt.it_out) go to 1
      if (taskid.eq.MASTER) 
     > write(6,*) '*** arived at the Maximum iiteration.   terminated.'
 
*::::::::::::::::::::  end of iteration loop  :::::::::::::::::::::::::

   2  continue

*     **** close xyz and MOTION files ****
      call xyz_end()
      call MOTION_end()
      if (taskid.eq.MASTER) then
        close(unit=31)
        close(unit=32)
        close(unit=34)
      end if
*     **** clsoe OMOTION file ****
      if (mulliken) call Orb_End()


      if (taskid.eq.MASTER) CALL MESSAGE(3)
      if (taskid.eq.MASTER) call current_second(cpu3)


*         |****************************************|
*********** produce CHECK file and diagonalize hml *****************
*         |****************************************|

*     **** produce CHECK FILE ****
      if (taskid.eq.MASTER) then
         call util_file_name('CHECK',.true.,
     >                               .false.,
     >                       full_filename)
         open(unit=17,file=full_filename,form='formatted')
      end if

*     **** check total number of electrons ****
      do ms =1,ispin
         call D3dB_r_dsum(1,dbl_mb(dn(1)+(ms-1)*n2ft3d),sumall)
         en(ms) = sumall*lattice_omega()
     >             /dble(ngrid(1)*ngrid(2)*ngrid(3))
      end do
      if (taskid.eq.MASTER) then
         write(17,1320) (en(ms),ms=1,ispin)
      end if

*     **** comparison between hamiltonian an lambda matrix ****
      n = ne(1)
      nn = n*n
      if (taskid.eq.MASTER) then
         write(17,1330)
         do ms=1,ispin
            do i=n1(ms),n2(ms)
               ii = i-n1(ms)
               do j=n1(ms),n2(ms)
                  jj = j-n1(ms)
                  index = (ii+1) + jj*n +(ms-1)*nn
                  write(17,1340) ms,ii+1,jj+1,
     >                           dbl_mb(hml(1)+index-1),
     >                           dbl_mb(lmd(1)+index-1),
     >             dbl_mb(hml(1)+index-1)-dbl_mb(lmd(1)+index-1)
               end do
            end do
         end do
      end if

*     **** check orthonormality ****
      if (taskid.eq.MASTER) then
         write(17,1350)
      end if

      do ms=1,ispin
         do i=n1(ms),n2(ms)
            ii = i-n1(ms)+1
            do j=i,n2(ms)
               jj = j-n1(ms)+1
               call Pack_cc_dot(1,dcpl_mb(psi1(1)+(i-1)*npack1),
     >                          dcpl_mb(psi1(1)+(j-1)*npack1),
     >                          w)
               if (taskid.eq.MASTER) then
                  write(17,1360) ms,ii,jj,w
               end if
            end do
         end do
      end do

*     ***** diagonalize the hamiltonian matrix *****
      n = ne(1)
      nn = n*n
      call dcopy(2*nemax,0.0d0,0,dbl_mb(eig(1)),1)
      do ms=1,ispin
         call eigen(n,ne(ms),
     >              dbl_mb(hml(1)+(ms-1)*nn),
     >              dbl_mb(eig(1)+(ms-1)*n),
     >              dbl_mb(lmd(1)),ierr)
      end do

*     *** rotate current psi ***
      call dcopy(2*npack1*nemax,0.0d0,0,dcpl_mb(psi2(1)),1)
      do ms=1,ispin
         do j=n1(ms),n2(ms)
            jj = j-n1(ms)
            do i=n1(ms),n2(ms)
               ii = i-n1(ms)
               index = (ii+1) + jj*n + (ms-1)*nn
             
               call Pack_cc_daxpy(1,dbl_mb(hml(1)+index-1),
     >                            dcpl_mb(psi1(1)+(i-1)*npack1),
     >                            dcpl_mb(psi2(1)+(j-1)*npack1))
            end do
         end do
      end do

*     *** rotate current v_psi ***
      call dcopy(2*npack1*nemax,dcpl_mb(psi0(1)),1,
     >                          dcpl_mb(psi1(1)),1)
      call dcopy(2*npack1*nemax,0.0d0,0,dcpl_mb(psi0(1)),1)
      do ms=1,ispin
         do j=n1(ms),n2(ms)
            jj = j-n1(ms)
            do i=n1(ms),n2(ms)
               ii = i-n1(ms)
               index = (ii+1) + jj*n + (ms-1)*nn
             
               call Pack_cc_daxpy(1,dbl_mb(hml(1)+index-1),
     >                            dcpl_mb(psi1(1)+(i-1)*npack1),
     >                            dcpl_mb(psi0(1)+(j-1)*npack1))
            end do
         end do
      end do




*     **** close check file ****
      if (taskid.eq.MASTER) then
         close(17)
      end if





*                |***************************|
****************** report summary of results **********************
*                |***************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      call center_v_geom(vcx,vcy,vcz)
      call center_v_mass(vgx,vgy,vgz)

      if (taskid.eq.MASTER) then
         write(6,1300)
         write(6,1410)
         write(6,1420)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(6,1191) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           else
           write(6,1190) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           end if
         end do
         write(6,1200) cx,cy,cz
         write(6,1210) gx,gy,gz

         call pspw_qmmm_Print(6)

         write(6,1421)
         write(6,1192) (I,ion_aname(I),
     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
         write(6,1200) vcx,vcy,vcz
         write(6,1210) vgx,vgy,vgz

         call pspw_qmmm_Print2(6)

         write(6,*)
         write(6,1320) en(1),en(ispin),' (real space)'
         write(6,1430) E(2),E(2)/ion_nion()
         write(6,1440) E(5),E(5)/n2(ispin)
         write(6,1450) E(6),E(6)/n2(ispin)
         write(6,1460) E(7),E(7)/n2(ispin)
         write(6,1470) E(8),E(8)/ion_nion()
         write(6,1471) E(3),E(3)/n2(ispin)
         write(6,1472) ion_ke(),ion_ke()/ion_nion()


         if (pspw_qmmm_found()) then
            write(6,1700)
            write(6,1701)
            write(6,1702) E(11)
            write(6,1703) E(12)
            write(6,1704) E(13)
            write(6,1705) E(14)
            write(6,1706) E(15)
         end if


         if (control_Nose()) then
           write(6,1473) E(9),E(9)/n2(ispin)
           write(6,1474) E(10),E(10)/ion_nion()
         end if
         write(6,1226) E(3),ion_ke(),ion_com_ke()
         eke = eke/dble(it_out)
         eke = 2.0d0*eke/kb/n2(ispin)/pack_nwave_all(1)

*       **** write out Temperatures ****
         write(6,1491) eke
         eki0 = ion_Temperature()
         if (pspw_qmmm_found()) then
            eki1 =pspw_qmmm_Temperature()
            sum = ion_nion() + pspw_qmmm_nion() - 2.0d0
            eki0 = eki0*((ion_nion()-2.0d0)/sum)
     >           + eki1*((pspw_qmmm_nion()-2.0d0)/sum)
         end if
         write(6,1480) eki0
         write(6,1490) ion_com_Temperature()

         write(6,1500)
         NN=NE(1)-NE(2)
         EV=27.2116d0
         do i=1,NN
           write(6,1510) dbl_mb(EIG(1)+i-1),dbl_mb(EIG(1)+i-1)*EV
         end do
         do i=1,ne(2)
           write(6,1510) dbl_mb(EIG(1)+i-1+NN),
     >                    dbl_mb(EIG(1)+i-1+NN)*EV,
     >                    dbl_mb(EIG(1)+i-1+n1(2)-1),
     >                    dbl_mb(EIG(1)+i-1+n1(2)-1)*EV
         end do

*        *** Extra energy output added for QA test ****
         write(6,1600) E(2)
      end if

*                |***************************|
******************         Prologue          **********************
*                |***************************|

*     **** calculate spin contamination ****
      call Calculate_psi_spin2(ispin,ne,nfft3d,dcpl_mb(psi2(1)),
     >                         dipole)

*     **** calculate the Dipole ***
      call Calculate_Dipole(ispin,ne,n2ft3d,dbl_mb(dn(1)),dipole)

      
*     ***** write wavefunctions and v_wavefunctions ****
      call psi_write(ispin,ne,dcpl_mb(psi2(1)))
      call v_psi_write(ispin,ne,dcpl_mb(psi0(1)))

*     **** write geometry to rtdb ****
      call pspw_qmmm_write(rtdb)
      call ion_write(rtdb)

*     **** deallocate heap memory ****
      if (control_version().eq.3) call ewald_end()
      call strfac_end()
      if (controL_version().eq.3) call coulomb_end()
      if (controL_version().eq.4) call coulomb2_end()
      call ke_end()
      call mask_end()
      call Pack_end()
      call G_end()
      call ion_end()
      call psp_end()
      call ion_end_FixIon()
      call pspw_qmmm_end()

      value = MA_free_heap(psir(2))
      value = MA_free_heap(Hpsi(2))
      value = MA_free_heap(dn(2))
      value = MA_free_heap(eig(2))
      value = MA_free_heap(hml(2))
      value = MA_free_heap(lmd(2))
      value = MA_free_heap(psi0(2))
      value = MA_free_heap(psi1(2))
      value = MA_free_heap(psi2(2))
      call D3dB_end(1)
      if (control_version().eq.4) call D3dB_end(2)

*     **** do anaylysis on MOTION files ****
      call cpmd_properties(rtdb)


*                |***************************|
****************** report consumed cputime   **********************
*                |***************************|
      if (taskid.eq.MASTER) then
         CALL current_second(cpu4)

         T1=CPU2-CPU1
         T2=CPU3-CPU2
         T3=CPU4-CPU3
         T4=CPU4-CPU1
         AV=T2/dble(icount*it_in)
         write(6,*)
         write(6,*) '-----------------'
         write(6,*) 'cputime in seconds'
         write(6,*) 'prologue    : ',T1
         write(6,*) 'main loop   : ',T2
         write(6,*) 'epilogue    : ',T3
         write(6,*) 'total       : ',T4
         write(6,*) 'cputime/step: ',AV
         write(6,*)
         write(6,*) '-------------------------------'
         write(6,*) 'Time spent doing:'
         write(6,*) '  FFTs                       : ', times(1)
         write(6,*) '  dot products               : ', times(2)
         write(6,*) '  orthonormalization         : ', times(3)
         write(6,*) '  exchange correlation       : ', times(4)
         write(6,*) '  local pseudopotentials     : ', times(5)
         write(6,*) '  non-local pseudopotentials : ', times(6)
         write(6,*) '  hartree potentials         : ', times(7)
         write(6,*) '  structure factors          : ', times(8)
         write(6,*) '  masking and packing        : ', times(9)
         write(6,*) '  total energy evalulation   : ', times(10)
         write(6,*) '  density                    : ', times(11)
         write(6,*) '  allocate and deallocate    : ', times(12)
         write(6,*)
         CALL MESSAGE(4)
      end if 


      call Parallel_Finalize()
      cpmdv5 = value
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
 1115 FORMAT(/' options:')
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
 1159 FORMAT(/' total charge:',I2)
 1160 FORMAT(/' atomic composition:')
 1170 FORMAT(7(5X,A2,':',I3))
 1180 FORMAT(/' initial position of ions:')
 1181 FORMAT(/' initial velocity of ions:')
 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' ) - atomic mass= ',F6.3,' ')
 1191 FORMAT(5X, I4, A3  ,' (',3F11.5,
     >       ' ) - atomic mass= ',F6.3,' - fixed')
 1192 FORMAT(5X, I4, A3  ,' (',3F11.5,' )')
 1200 FORMAT(5X,'  G.C. ',' (',3F11.5,' )')
 1210 FORMAT(5X,' C.O.M.',' (',3F11.5,' )')
 1220 FORMAT(/' number of electrons: spin up=',I3,'  spin down=',I3,A)
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
 1290 FORMAT(5X, ' cooling/heatting rates: ',e12.5,' (psi)',2x,
     >                                       e12.5,' (ion)')
 1294 format(/' Constant Energy Simulation                     ')
 1295 format(/' Nose-Hoover Simulation - thermostat parameters:')
 1297 format(5x, 'Te =',f8.2,' Qe =',e8.3,' 2*pi/we=',e8.3,' Ee0=',e8.3)
 1298 format(5x, 'Tr =',f8.2,' Qr =',e8.3,' 2*pi/wr=',e8.3,' Er0=',e8.3)
 1299 format(//' Mulliken analysis output generated            ')
 1300 FORMAT(//)
 1305 FORMAT(10X,'================ iteration =========================')
 1310 FORMAT(I8,2E19.10,2E14.5)
 1311 format(100e19.10)
 1312 format(e14.6,i3)
 1313 format(3i4)
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
 1420 FORMAT(/' final position of ions:')
 1421 FORMAT(/' final velocity of ions:')
 1430 FORMAT(/' total     energy    :',E19.10,' (',E15.5,'/ion)')
 1440 FORMAT( ' total orbital energy:',E19.10,' (',E15.5,'/electron)')
 1450 FORMAT( ' hartree   energy    :',E19.10,' (',E15.5,'/electron)')
 1460 FORMAT( ' exc-corr  energy    :',E19.10,' (',E15.5,'/electron)')
 1470 FORMAT( ' ion-ion   energy    :',E19.10,' (',E15.5,'/ion)')
 1471 FORMAT(/' Kinetic energy (elc)    :',E19.10,' (',E15.5,'/elc)')
 1472 FORMAT( ' Kinetic energy (ion)    :',E19.10,' (',E15.5,'/ion)')
 1473 FORMAT( ' thermostat energy (elc) :',E19.10,' (',E15.5,'/elc)')
 1474 FORMAT( ' thermostat energy (ion) :',E19.10,' (',E15.5,'/ion)')
 1480 FORMAT(' Temperature :    ',F10.1,' K (ion)')
 1490 FORMAT('             :    ',F10.1,' K (c.o.m.)')
 1491 FORMAT(' Temperature :    ',F10.1,' K (elc)')
 1500 FORMAT(/' orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
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

      cpmdv5 = value
      return
      END
