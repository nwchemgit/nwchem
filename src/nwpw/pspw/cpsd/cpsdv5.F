
***********************************************************************
*                      cpsdv5-mpi  (MPI code)                         *
*                                                                     *
*     This is a developing cpsdv5 parallel code for the SP2           *
*       + mpi message passing library used                            *
*       + ngp is used instead of nfft in this proceudure              *
*       + error checking is based on aimd.h parameters                *
*         then control file                                           *
*       + my own slap-decomposed parallel 3d-FFT(real->complex) used  *
*                                                                     *
*                                                                     *
***********************************************************************

      logical function cpsdv5(rtdb)
*
* $Id: cpsdv5.F,v 1.2 2001-09-29 00:28:06 bylaska Exp $
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
*    13 - hpsi and update
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
      integer npack1

*     **** electronic variables ****
      integer ispin
      integer ne(2),n1(2),n2(2),nemax
      real*8  en(2)
      real*8  dipole(3)

*     complex*16 psi1(nfft3d,nemax)
*     complex*16 psi2(nfft3d,nemax)
*     real*8     dn(n2ft3d,2)
*     complex*16 Hpsi(nfft3d,nemax)
*     complex*16 psir(nfft3d,nemax)
      integer psi1(2),psi2(2)
      integer dn(2)
      integer Hpsi(2),psir(2)
    

*     ***** energy variables ****
      real*8  E(20)

*     real*8  eig(2*nemax)
*     real*8  hml(2*nemax*nemax)
*     real*8  lmd(2*nemax*nemax)
      integer eig(2),hml(2),lmd(2)




*     **** error variables ****
      integer ierr

*     **** local variables ****
      integer ms
      real*8  deltae,deltac,deltar
      real*8  gx,gy,gz,cx,cy,cz,sum1,sum2
      real*8  EV,pi
      integer i,j,k,ia,n,nn
      integer ii,jj,index
      integer icount,it_in,it_out
      real*8 w,sumall,virial
      integer nfft3
      parameter (nfft3=32)
      character*255 full_filename
  


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
      logical  control_read,control_move,ion_init,ion_q_FixIon
      external control_read,control_move,ion_init,ion_q_FixIon
 
      integer  pack_nwave_all
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      external pack_nwave_all
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm

      character*12 control_boundry
      external     control_boundry

      logical      pspw_SIC,pspw_SIC_relaxed
      logical      psp_semicore,control_Mulliken,Waterpsp_found
      real*8       psp_rcore,psp_ncore,psp_rlocal
      external     pspw_SIC,pspw_SIC_relaxed
      external     psp_semicore,control_Mulliken,Waterpsp_found
      external     psp_rcore,psp_ncore,psp_rlocal

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
      
      value = control_read(1,rtdb)
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


*     **** initialize lattice and packing data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_init()
      call Pack_npack(1,npack1)      

*     ***** allocate psi2,and psi1 wavefunctions ****
      call psi_get_ne(ispin,ne)
      value = MA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'psi1',psi1(2),psi1(1))
      if (.not. value) call errquit('out of heap memory',0)


*     *****  read psi2 wavefunctions ****
      call psi_read(ispin,ne,dcpl_mb(psi2(1)))

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)
      nemax = ne(1)+ne(2)

*     **** Ortho Check ****
      call Grsm_gg_trace(npack1,(ne(1)+ne(2)),
     >                        dcpl_mb(psi2(1)),
     >                        dcpl_mb(psi2(1)),
     >                        sum2)

      sum1 = dble(ne(1) + ne(2))
      if (dabs(sum2-sum1).gt.1.0d-10) then
         if (taskid.eq.MASTER)
     >      write(*,*) "Waring: Gram-Schmidt being performed on psi2"
         call Grsm_g_MakeOrtho(npack1,ne(1),dcpl_mb(psi2(1)))
         if (ispin.gt.1) then
           call Grsm_g_MakeOrtho(npack1,ne(2),
     >                           dcpl_mb(psi2(1) + ne(1)*npack1))
         end if
      end if




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

      if (.not. value) call errquit('out of heap memory',0)



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
      if (psp_semicore(0)) call semicore_check()

*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)

*     **** initialize SIC ****
      call pspw_init_SIC(ne,rtdb)

*     **** initialize QM/MM ****
      call pspw_init_APC(rtdb)
      call Waterpsp_init(rtdb)



*                |**************************|
******************   summary of input data  **********************
*                |**************************|

      if (taskid.eq.MASTER) then
         write(6,1111) np
         write(6,1115)
         IF(control_move()) THEN
           write(6,1120) 'yes'
         ELSE
           write(6,1120) 'no'
         ENDIF
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
         call pspw_print_SIC(6)
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
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(6,1191) I,psp_atom(ion_katm(I)),(ion_rion(K,I),K=1,3)
           else
           write(6,1190) I,psp_atom(ion_katm(I)),(ion_rion(K,I),K=1,3)
           end if
         end do
         write(6,1200) cx,cy,cz
         write(6,1210) gx,gy,gz

         call Waterpsp_Print(6)

         write(6,1220) ne(1),ne(ispin),' ( fourier space)'
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
         write(6,1305)
         call util_flush(6)
      end if

*                |***************************|
******************     start iterations      **********************
*                |***************************|

      if (taskid.eq.MASTER) call current_second(cpu2)
      if (taskid.eq.MASTER) CALL MESSAGE(2)
      it_in  = control_it_in()
      it_out = control_it_out()
      icount = 0
   1  continue
         icount = icount + 1
         call inner_loop(ispin,ne,
     >             npack1,nfft3d,nemax,
     >             dcpl_mb(psi1(1)),dcpl_mb(psi2(1)),dbl_mb(dn(1)),
     >             it_in,E,deltae,deltac,deltar,
     >             dbl_mb(hml(1)),dbl_mb(lmd(1)),
     >             dcpl_mb(psir(1)),dcpl_mb(Hpsi(1)))
         if (taskid.eq.MASTER) then 
           write(6,1310) icount*it_in,E(1),deltae,deltac,deltar
           call util_flush(6)
         end if
         if ((deltae.gt.0.0d0).and.(icount.gt.1)) then
            if (taskid.eq.MASTER) 
     >       write(6,*) ' *** Energy going up.  iteration terminated.'
            go to 2
         end if
         deltae = dabs(deltae)
         if ((deltae.lt.control_tole()).and.
     >       (deltac.lt.control_tolc()).and.
     >       (deltar.lt.control_tolr())) then
            if (taskid.eq.MASTER) 
     >       write(6,*) ' *** tolerance ok.     iteration terminated.'
            go to 2
         end if
      if (icount.lt.it_out) go to 1
      if (taskid.eq.MASTER) 
     > write(6,*) '*** arived at the Maximum iiteration.   terminated.'

*::::::::::::::::::::  end of iteration loop  :::::::::::::::::::::::::

   2  continue
      if (taskid.eq.MASTER) CALL MESSAGE(3)
      if (taskid.eq.MASTER) call current_second(cpu3)



*         |****************************************|
*********** produce CHECK file and diagonalize hml *****************
*         |****************************************|

*     **** produce CHECK FILE ****
      if (taskid.eq.MASTER) then
         call util_file_name('CHECK',.true.,
     >                               .false.,
     >                        full_filename)
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

*     **** close check file ****
      if (taskid.eq.MASTER) then
         close(17)
      end if





*                |***************************|
****************** report summary of results **********************
*                |***************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

      if (taskid.eq.MASTER) then
         write(6,1300)
         write(6,1410)
         write(6,1420)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(6,1191) I,psp_atom(ion_katm(I)),(ion_rion(k,i),K=1,3)
           else
           write(6,1190) I,psp_atom(ion_katm(I)),(ion_rion(k,i),K=1,3)
           end if
         end do
         write(6,1200) cx,cy,cz
         write(6,1210) gx,gy,gz

         call Waterpsp_Print(6)

         write(6,*)
         write(6,1320) en(1),en(ispin),' (real space)'

         write(6,1430) E(1),E(1)/ion_nion()
         write(6,1440) E(2),E(2)/n2(ispin)
         write(6,1450) E(3),E(3)/n2(ispin)
         write(6,1460) E(4),E(4)/n2(ispin)
         if (pspw_SIC())
     >     write(6,1455) E(14),E(14)/n2(ispin)
         write(6,1470) E(5),E(5)/ion_nion()
         write(6,1480) E(6),E(6)/n2(ispin)
         write(6,1490) E(7),E(7)/n2(ispin)
         write(6,1495) E(8),E(8)/n2(ispin)
         write(6,1496) E(9),E(9)/n2(ispin)
         write(6,1497) E(10),E(10)/n2(ispin)
         if (pspw_SIC().and.pspw_SIC_relaxed())
     >     write(6,1499) E(15),E(15)/n2(ispin)
         virial = (E(10)+E(9)+E(8)+E(7))/E(6)
         write(6,1498) virial
 
         if (Waterpsp_found()) then
            write(6,1700) E(11),E(11)/n2(ispin)
            write(6,1701) E(12),E(12)/ion_nion()
            write(6,1702) E(13)
         end if



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

*        ***** extra energy output for QA test ****
         write(6,1600) E(1)
      end if

*                |***************************|
******************         Prologue          **********************
*                |***************************|

*     **** calculate spin contamination ****
      call Calculate_psi_spin2(ispin,ne,npack1,dcpl_mb(psi2(1)),w)

*     **** calculate the Dipole ***
      call Calculate_Dipole(ispin,ne,n2ft3d,dbl_mb(dn(1)),dipole)
      
*     **** perfom Lubin and Mulliken analysis ***
      if (control_Mulliken()) then

*       **** Lubin Water Analysis ***
        call pspw_Lubin_water_analysis(rtdb,ispin,ne,n2ft3d,
     >                                 dbl_mb(dn(1)))

*       **** Analysis ***
        call pspw_analysis(rtdb,ispin,ne,dcpl_mb(psi2(1)),
     >                                   dbl_mb(eig(1)))

*       **** generate APC *****
        call pspw_dngen_APC(ispin,ne,dbl_mb(dn(1)))
        call pspw_print_APC(6)

      end if

*     ***** write psi2 wavefunctions ****
      call psi_write(ispin,ne,dcpl_mb(psi2(1)))
 
*     **** write geometry to rtdb ****
      call Waterpsp_write(rtdb)
      call ion_write(rtdb)


*     **** deallocate heap memory ****
      if (control_version().eq.3) call ewald_end()
      call strfac_end()
      if (control_version().eq.3) call coulomb_end()
      if (control_version().eq.4) call coulomb2_end()
      call ke_end()
      call mask_end()
      call Pack_end()
      call G_end()
      call ion_end()
      call psp_end()
      call ion_end_FixIon()
      call pspw_end_SIC()
      call pspw_end_APC()
      call Waterpsp_end()

      value = MA_free_heap(psir(2))
      value = MA_free_heap(Hpsi(2))
      value = MA_free_heap(dn(2))
      value = MA_free_heap(eig(2))
      value = MA_free_heap(hml(2))
      value = MA_free_heap(lmd(2))
      value = MA_free_heap(psi1(2))
      value = MA_free_heap(psi2(2))
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
         write(6,*) '  total energy evaluation    : ',times(10)
         write(6,*) '  density                    : ',times(11)
         write(6,*) '  allocate and deallocate    : ',times(12)
         write(6,*) '  Hpsi and update            : ',times(13)
         write(6,*)
         CALL MESSAGE(4)
      end if 


      call Parallel_Finalize()
      cpsdv5 = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*     Car-Parrinello microcluster calculation      *')
 1030 FORMAT(10X,'*     [     steepest descent minimization   ]      *')
 1035 FORMAT(10x,'*     [ NorthWest Chemistry implementation ]       *')
 1040 FORMAT(10X,'*            version #5.00   06/01/99              *')
 1041 FORMAT(10X,'*    This code was developed by Eric J. Bylaska,   *')
 1042 FORMAT(10X,'*    and was based upon algorithms and code        *')
 1043 FORMAT(10X,'*    developed by the group of Prof. John H. Weare *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ PSPW input data ===================')
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
 1191 FORMAT(5X, I4, A3  ,' (',3F11.5,' ) - fixed')
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
     &        ' (electron)',E12.3,' (ion)')
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
 1455 FORMAT( ' SIC       energy    :',E19.10,' (',E15.5,'/electron)')
 1460 FORMAT( ' exc-corr  energy    :',E19.10,' (',E15.5,'/electron)')
 1470 FORMAT( ' ion-ion   energy    :',E19.10,' (',E15.5,'/ion)')
 1480 FORMAT(/' K.S. kinetic energy :',E19.10,' (',E15.5,'/electron)')
 1490 FORMAT( ' K.S. V_l  energy    :',E19.10,' (',E15.5,'/electron)')
 1495 FORMAT( ' K.S. V_nl energy    :',E19.10,' (',E15.5,'/electron)')
 1496 FORMAT( ' K.S. V_Hart energy  :',E19.10,' (',E15.5,'/electron)')
 1497 FORMAT( ' K.S. V_xc energy    :',E19.10,' (',E15.5,'/electron)')
 1498 FORMAT( ' Virial Coefficient  :',E19.10)
 1499 FORMAT( ' K.S. SIC  energy    :',E19.10,' (',E15.5,'/electron)')
 1500 FORMAT(/' orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1600 FORMAT(/' Total PSPW energy   :',E19.10)
 1700 FORMAT(/' water-elc energy    :',E19.10,' (',E15.5,'/electron)')
 1701 FORMAT( ' water-ion energy    :',E19.10,' (',E15.5,'/ion)')
 1702 FORMAT( ' water-water energy  :',E19.10)

 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(6,9010) ierr
      call Parallel_Finalize()

      cpsdv5 = value
      return
      END
