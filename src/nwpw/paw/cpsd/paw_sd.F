*
* $Id: paw_sd.F,v 1.43 2003-02-23 23:32:05 bylaska Exp $
*
***********************************************************************
*                      paw_sd                                         *
*                                                                     *
*     This is a developing PAW steepest descent code for NWChem.      *
*                                                                     *
*                                                                     * 
*  Authors: Marat Valiev and Eric J. Bylaska                          *
*                                                                     *
***********************************************************************

      logical function paw_sd(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
      
      logical value


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
      real*8  dipole(3)

*     complex*16 psi1(nfft3d,nemax)
*     complex*16 psi2(nfft3d,nemax)
*     real*8     dn(n2ft3d,2)
*     complex*16 Hpsi(nfft3d,nemax)
*     complex*16 psir(nfft3d,nemax)
      integer psi1(2),psi2(2)
      integer dn(2),dn_cmp_smooth(2)
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
      real*8      ewald_rcut,ion_amass
      real*8      ewald_mandelung
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg
      integer     ewald_ncut,ewald_nshl3d
      character   spdf_name
      character*2 ion_aname,ion_atom
      external    ewald_rcut,ion_amass
      external    ewald_mandelung
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg
      external    ewald_ncut,ewald_nshl3d
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
      integer  control_ngrid,pack_nwave
      integer  ion_nion,ion_natm,ion_katm,ion_nkatm
      external pack_nwave_all
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave
      external ion_nion,ion_natm,ion_katm,ion_nkatm

      character*12 control_boundry
      external     control_boundry

      logical      control_Mulliken
      external     control_Mulliken
      logical      psi_filefind
      external     psi_filefind
      real*8       nwpw_timing,paw_psi_CheckOrtho
      external     nwpw_timing,paw_psi_CheckOrtho

*                            |************|
*****************************|  PROLOGUE  |****************************
*                            |************|

      value = .true.
      pi = 4.0d0*datan(1.0d0)

      call nwpw_timing_init()
      call dcopy(20,0.0d0,0,E,1)


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
         write(6,1010)
         write(6,1042)
         write(6,1043)
         write(6,1044)
         write(6,1045)
         write(6,1046)
         write(6,1047)
         write(6,1048)
         write(6,1049)
         write(6,1050)
         write(6,1051)
         write(6,1052)
         write(6,1053)
         write(6,1054)
         write(6,1055)
         write(6,1056)
         write(6,1057)
         write(6,1058)
         write(6,1010)
         write(6,1000)
         call nwpw_message(1)
         write(6,1110)
      end if
      
      value = control_read(6,rtdb)
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

 
*     **** read ions ****
      value = ion_init(rtdb)
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)


*     **** allocate paw data structure and read in paw basis into it ****
      write(*,*) "generating PAW basis - debug"
      value = MA_set_auto_verify(.true.)
      call init_paw_basis()

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
      
*     *** initialize paw atomc potentials ***      
      call paw_vloc_init()
      call init_paw_pot_hartree()
      call init_paw_pot_comp()
      call paw_mult_init()


*     **** initialize G,mask,ke,and coulomb data structures ****
      call ke_init()
      if (control_version().eq.3) call coulomb_init()
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()
      if (control_version().eq.3) call ewald_init()


*     **** generate initial wavefunction if it does not exist ****
!      if (.not.psi_filefind()) then
!        call psi_new()
!      end if

*     ***** allocate psi2,and psi1 wavefunctions ****
      call psi_get_ne(ispin,ne)
      value = MA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'psi1',psi1(2),psi1(1))
      if (.not. value) call errquit('paw_sd:out of heap memory',0)


*     *****  read psi2 wavefunctions ****
      call psi_read(ispin,ne,dcpl_mb(psi2(1)))

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
     >        MA_alloc_get(mt_dbl,(2*nfft3d),
     >             'dn_cmp_smooth',dn_cmp_smooth(2),dn_cmp_smooth(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'Hpsi',Hpsi(2),Hpsi(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,nfft3d*(ne(1)+ne(2)),
     >                     'psir',psir(2),psir(1))
      if (.not. value) call errquit('paw_sd: out of heap memory',1)

*     *** intialize overlap coefficient data structure ***
      call phafac()
      call paw_ovlp_init(ispin,ne)
      call init_paw_density(0)
      call init_paw_xc(ispin,0)
      call paw_ovlp_coeff_set(dcpl_mb(psi2(1)))
      call paw_ovlp_weights_set()


      !**** Ortho Check ****
      do ms=1,ispin
        deltae=paw_psi_CheckOrtho(npack1,ne(ms),
     >                   dcpl_mb(psi2(1)+(n1(ms)-1)*npack1))

        if (deltae.gt.1.0d-10) then
          call paw_psi_MakeOrtho(npack1,ne(ms),
     >                   dcpl_mb(psi2(1)+(n1(ms)-1)*npack1))
          deltac=paw_psi_CheckOrtho(npack1,ne(ms),
     >                   dcpl_mb(psi2(1)+(n1(ms)-1)*npack1))
          if (taskid.eq.MASTER) then
            if (ms.eq.1) then
              write(*,*) "Warning: ",
     >                   "Gram-Schmidt performed on up spin of psi2 "
              write(*,*) "       : (old error=",deltae,
     >                   " new error=",deltac,")"
            end if
            if (ms.eq.2) then
              write(*,*) "Warning: ",
     >                   "Gram-Schmidt performed on down spin of psi2 "
              write(*,*) "       : (old error=",deltae,
     >                   " new error=",deltac,")"
            end if

          end if
        end if
     
      end do


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

         write(6,1140)
         do ia = 1,ion_nkatm()
            write(6,1141) ia,ion_atom(ia),
     >                    paw_basis_ion_charge(ia),
     >                    paw_basis_core_charge(ia)
            write(6,1143) paw_basis_sphere_radius(ia)
            write(6,1144) paw_basis_sigma(ia),paw_basis_sigma(ia)
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
      if (taskid.eq.MASTER) CALL nwpw_MESSAGE(2)
      it_in  = control_it_in()
      it_out = control_it_out()
      icount = 0
   1  continue
         icount = icount + 1
         
         call paw_inner_loop(ispin,ne,
     >             npack1,nfft3d,nemax,
     >             dcpl_mb(psi1(1)),dcpl_mb(psi2(1)),
     >             dbl_mb(dn(1)),dbl_mb(dn_cmp_smooth(1)),
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
      if (taskid.eq.MASTER) CALL nwpw_MESSAGE(3)
      if (taskid.eq.MASTER) call current_second(cpu3)


      call psi_marat_write("marat.elc",ispin,ne,dcpl_mb(psi2(1)))


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

      !*** comp charge ****
      call D3dB_r_dsum(1,dbl_mb(dn_cmp_smooth(1)),sum1)
      sum1 = sum1*lattice_omega()
     >      /dble(ngrid(1)*ngrid(2)*ngrid(3))

      !*** psi charge ****
      do ms =1,ispin
         call D3dB_r_dsum(1,dbl_mb(dn(1)+(ms-1)*n2ft3d),sumall)
         en(ms) = ne(ms) + 0.5d0*sum1
     >          + sumall*lattice_omega()
     >             /dble(ngrid(1)*ngrid(2)*ngrid(3))
      end do

      if (taskid.eq.MASTER) then
         write(17,1321) en(1)+en(ispin)
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
               call paw_overlap_matrix_gen(1,1,
     >                          dcpl_mb(psi1(1)+(i-1)*npack1),
     >                          dcpl_mb(psi1(1)+(j-1)*npack1),
     >                          w)
               if (taskid.eq.MASTER) then
                  write(17,1360) ms,ii,jj,w
               end if
            end do
         end do
      end do

*     **** close check file ****
      if (taskid.eq.MASTER) then
         close(17)
      end if


*     ***** diagonalize the hamiltonian matrix ****
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
           write(6,1191) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           else
           write(6,1190) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           end if
         end do
         write(6,1200) cx,cy,cz
         write(6,1210) gx,gy,gz


         write(6,*)
         write(6,1321) en(1)+en(ispin),' (real space)'

         call paw_energy_output(6,ion_nion(),n2(ispin),E)
c         write(6,1430) E(1),E(1)/ion_nion()
c         write(6,1440) E(2),E(2)/n2(ispin)
c         write(6,1450) E(3),E(3)/n2(ispin)
c         write(6,1460) E(4),E(4)/n2(ispin)
c         write(6,1470) E(5),E(5)/ion_nion()
c         write(6,1480) E(6),E(6)/n2(ispin)
c         write(6,1490) E(7),E(7)/n2(ispin)
c         write(6,1495) E(8),E(8)/n2(ispin)
c         write(6,1496) E(9),E(9)/n2(ispin)
c         write(6,1497) E(10),E(10)/n2(ispin)
c         virial = (E(10)+E(9)+E(8)+E(7))/E(6)
c         write(6,1498) virial


*        **** write out KS eigenvalues ****
         write(6,1500)
         NN=NE(1)-NE(2)
         EV=27.2116d0
         do i=1,NN
           write(6,1510) dbl_mb(EIG(1)+i-1),dbl_mb(EIG(1)+i-1)*EV
         end do
         do i=1,ne(2)
           write(6,1510) dbl_mb(EIG(1)+i-1+NN),
     >                   dbl_mb(EIG(1)+i-1+NN)*EV,
     >                   dbl_mb(EIG(1)+i-1+n1(2)-1),
     >                   dbl_mb(EIG(1)+i-1+n1(2)-1)*EV
         end do

*        ***** extra energy output for QA test ****
         write(6,1600) E(1)
      end if

*                |***************************|
******************         Prologue          **********************
*                |***************************|

!*     **** calculate spin contamination ****
!      call Calculate_psi_spin2(ispin,ne,npack1,dcpl_mb(psi2(1)),w)
!
!*     **** calculate the Dipole ***
!      call Calculate_Dipole(ispin,ne,n2ft3d,dbl_mb(dn(1)),dipole)
!      
!
!*     ***** write psi2 wavefunctions ****
!      call psi_write(ispin,ne,dcpl_mb(psi2(1)))

*       call psi_marat_write(ispin,ne,dcpl_mb(psi2(1)))


!*     **** write geometry to rtdb ****
!      call ion_write(rtdb)


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
     
      call dealloc_paw_basis_data()
      call paw_proj_end()
      call paw_ovlp_end()
      call paw_comp_charge_end()
      call paw_mult_end()
      call paw_vloc_end()
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
      call paw_density_end()
      call paw_xc_end()

      value =           MA_free_heap(psir(2))
      value = value.and.MA_free_heap(Hpsi(2))
      value = value.and.MA_free_heap(dn(2))
      value = value.and.MA_free_heap(dn_cmp_smooth(2))
      value = value.and.MA_free_heap(eig(2))
      value = value.and.MA_free_heap(hml(2))
      value = value.and.MA_free_heap(lmd(2))
      value = value.and.MA_free_heap(psi1(2))
      value = value.and.MA_free_heap(psi2(2))
      if (.not. value) call errquit('paw_sd:error freeing heap',2)

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
         write(6,*) '  FFTs                       : ', nwpw_timing(1)
         write(6,*) '  dot products               : ', nwpw_timing(2)
         write(6,*) '  orthonormalization         : ', nwpw_timing(3)
         write(6,*) '  exchange correlation       : ', nwpw_timing(4)
         write(6,*) '  local pseudopotentials     : ', nwpw_timing(5)
         write(6,*) '  non-local pseudopotentials : ', nwpw_timing(6)
         write(6,*) '  hartree potentials         : ', nwpw_timing(7)
         write(6,*) '  structure factors          : ', nwpw_timing(8)
         write(6,*) '  masking and packing        : ', nwpw_timing(9)
         write(6,*) '  total energy evaluation    : ',nwpw_timing(10)
         write(6,*) '  density                    : ',nwpw_timing(11)
         write(6,*) '  allocate and deallocate    : ',nwpw_timing(12)
         write(6,*) '  Hpsi and update            : ',nwpw_timing(13)
         write(6,*)
         CALL nwpw_MESSAGE(4)
      end if 


      call Parallel_Finalize()
      paw_sd = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*           PAW microcluster calculation           *')
 1030 FORMAT(10X,'*     [     steepest descent minimization   ]      *')
 1035 FORMAT(10x,'*     [ NorthWest Chemistry implementation ]       *')
 1040 FORMAT(10X,'*            version #1.00   08/01/02              *')
 1041 FORMAT(10X,'*    Authors: Marat Valiev and Eric J. Bylaska     *')
 1042 FORMAT(10X,'*    This code is based upon algorithms and code   *')
 1043 FORMAT(10X,'*    developed by the group of Prof. John H. Weare *')
 1044 FORMAT(10X,'*                                                  *')
 1045 FORMAT(10X,'*    References:                                   *')
 1046 FORMAT(10X,'*                                                  *')
 1047 FORMAT(10X,'*    M. Valiev, E. J. Bylaska, A. Gramada,         *')
 1048 FORMAT(10X,'*    and J. H. Weare,                              *')
 1049 FORMAT(10X,'*    Reviews in Modern  Quantum Chemistry,         *')
 1050 FORMAT(10X,'*    1684 (World Scientific, Singapore, 2002)      *')
 1051 FORMAT(10X,'*                                                  *')
 1052 FORMAT(10X,'*    E. J. Bylaska, M. Valiev, R. Kawai,           *')
 1053 FORMAT(10X,'*    and J. H. Weare,                              *')
 1054 FORMAT(10X,'*    Computer Physics  Communications, 143 (2002)  *')
 1055 FORMAT(10X,'*                                                  *')
 1056 FORMAT(10X,'*    M. Valiev and J. H. Weare,                    *')
 1057 FORMAT(10X,'*    J. Phys. Chem. A 103, 10588 (1999).           *')
 1058 FORMAT(10X,'*                                                  *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ PAW input data ===================')
 1111 FORMAT(/' number of processors used:',I3)
 1115 FORMAT(/' options:')
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundary conditions   = ',A,'(version', I1,')')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1140 FORMAT(/' elements involved in the cluster:')
 1141 FORMAT(5X,I2,': ',A2,'  ion charge:',F4.1,'  core charge:',F4.1)

 1143 FORMAT(12x,' augmentation sphere radius  :',F6.3)
 1144 FORMAT(12x,' compensation sigma          :',F6.3,
     .   ' (',F6.3,' smooth)')

 1150 FORMAT(12x,' total number of projectors  :',I3)

 1151 FORMAT(12x,' n_ps (n) l          eig    #projector')
 1152 FORMAT(14X,I3,' (',I1,') ',A,F13.6,I14)

 1153 FORMAT(12X,' local potential used           : ',i2)
 1154 FORMAT(12X,' number of non-local projections: ',i2)
 1155 FORMAT(12X,' semicore corrections included  : ',
     >       F6.3,' (radius) ',F6.3,' (charge)')
 1156 FORMAT(12X,' aperiodic cutoff radius        : ',F6.3)
 1159 FORMAT(/' total charge:',I2)
 1160 FORMAT(/' atomic composition:')
 1170 FORMAT(7(5X,A2,':',I3))
 1180 FORMAT(/' initial position of ions:')
 1190 FORMAT(5X, I4, A3, ' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1191 FORMAT(5X, I4, A3, ' (',3F11.5,
     >       ' ) - atomic mass= ',F6.3,' - fixed')
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

 1250 FORMAT(5X,' density cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
 1251 FORMAT(5X,' wavefnc cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
 1260 FORMAT(5X,' ewald summation: cut radius=',F8.2,'  and',I3)
 1261 FORMAT(5X,'                   mandelung=',f11.8)
 1270 FORMAT(/' technical parameters:')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'ficticious mass=',F10.1)
 1290 FORMAT(5X, ' tolerance=',E8.3,' (energy)',E12.3,
     &        ' (electron)',E12.3,' (ion)')
 1300 FORMAT(//)
 1305 FORMAT(10X,'================ iteration =========================')
 1310 FORMAT(I8,E20.10,3E15.5)
 1320 FORMAT(' number of electrons: spin up=',F11.5,'  down=',F11.5,A)
 1321 FORMAT(' total charge of system:',F11.5,A)
 1330 FORMAT(/' comparison between hamiltonian and lambda matrix')
 1331 FORMAT(/' Hamiltonian matrix')
 1340 FORMAT(I3,2I3,' H=',E16.7,', L=',E16.7,', H-L=',E16.7)
 1341 FORMAT(I3,2I3,' H=',E16.6)
 1350 FORMAT(/' orthonormality')
 1360 FORMAT(I3,2I3,E18.7)
 1370 FORMAT(I3)
 1380 FORMAT(' ''',a,'''',I4)
 1390 FORMAT(I3)
 1400 FORMAT(I3,3E18.8/3X,3E18.8)
 1410 FORMAT(10X,'=============  summary of results  =================')
 1420 FORMAT( ' final position of ions:')
 1430 FORMAT(//' total     energy    :',E19.10,' (',E15.5,'/ion)')
 1431 FORMAT(/' QM Energies')
 1432 FORMAT( '------------')
 1433 FORMAT( ' total  QM energy    :',E19.10,' (',E15.5,'/ion)')
 1440 FORMAT( ' total orbital energy:',E19.10,' (',E15.5,'/electron)')
 1450 FORMAT( ' hartree   energy    :',E19.10,' (',E15.5,'/electron)')
 1455 FORMAT( ' SIC-hartree energy  :',E19.10,' (',E15.5,'/electron)')
 1456 FORMAT( ' SIC-exc-corr energy :',E19.10,' (',E15.5,'/electron)')
 1460 FORMAT( ' exc-corr  energy    :',E19.10,' (',E15.5,'/electron)')
 1470 FORMAT( ' ion-ion   energy    :',E19.10,' (',E15.5,'/ion)')
 1480 FORMAT(/' K.S. kinetic energy :',E19.10,' (',E15.5,'/electron)')
 1490 FORMAT( ' K.S. V_l  energy    :',E19.10,' (',E15.5,'/electron)')
 1495 FORMAT( ' K.S. V_nl energy    :',E19.10,' (',E15.5,'/electron)')
 1496 FORMAT( ' K.S. V_Hart energy  :',E19.10,' (',E15.5,'/electron)')
 1497 FORMAT( ' K.S. V_xc energy    :',E19.10,' (',E15.5,'/electron)')
 1498 FORMAT( ' Virial Coefficient  :',E19.10)
 1499 FORMAT( ' K.S. SIC-hartree energy  :',E19.10,
     >        ' (',E15.5,'/electron)')
 1501 FORMAT( ' K.S. SIC-exc-corr energy :',E19.10,
     >        ' (',E15.5,'/electron)')
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

      paw_sd = value
      return
      end

      character function spdf_name(l)
      implicit none
      integer l
      character name
      name = '?'
      if (l.eq.0) name = 's'
      if (l.eq.1) name = 'p'
      if (l.eq.2) name = 'd'
      if (l.eq.3) name = 'f'
      if (l.eq.4) name = 'g'
      if (l.eq.5) name = 'h'
      if (l.eq.6) name = 'i'
      if (l.eq.7) name = 'j'
      if (l.eq.8) name = 'k'
      if (l.eq.9) name = 'l'
      spdf_name = name
      return
      end
