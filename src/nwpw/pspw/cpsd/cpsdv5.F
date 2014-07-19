
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
* $Id$
*
      implicit none
      integer rtdb

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
cccc#include "frac_occ.fh"
      
*     **** parallel variables ****
      integer  taskid,np,np_i,np_j
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,nfft3d,n2ft3d
      integer npack1

*     **** electronic variables ****
      real*8 icharge
      integer ispin
      integer ne(2),n1(2),n2(2),nemax,neq(2),nemaxq
      real*8  en(2),en1(2),en2(2)
      real*8  dipole(3)

*     complex*16 psi1(nfft3d,nemax)
*     complex*16 psi2(nfft3d,nemax)
*     real*8     dn(n2ft3d,2)
*     complex*16 Hpsi(nfft3d,nemax)
*     complex*16 psir(nfft3d,nemax)
      integer psi1(2),psi2(2)
      integer occ1(2),occ2(2)
      integer dn(2)
      integer Hpsi(2),psir(2)
    

*     ***** energy variables ****
      real*8  E(50)

*     real*8  eig(2*nemax)
*     real*8  hml(2*nemax*nemax)
*     real*8  lmd(2*nemax*nemax)
      integer eig(2),hml(2),lmd(2)

*     **** psi smearing block ****
      logical fractional
      integer smearoccupation,smeartype
      real*8 smearfermi(2),smearcorrection,smearkT




*     **** error variables ****
      integer ierr

*     **** local variables ****
      integer ms,mapping,mapping1d
      real*8  deltae,deltac,deltar
      real*8  gx,gy,gz,cx,cy,cz,sum1,sum2
      real*8  EV,pi
      integer i,j,k,ia,n,nn
      integer ii,jj,indx,indx1
      integer icount,it_in,it_out
      real*8 w,sumall,virial,dv
      integer nfft3
      parameter (nfft3=32)
      character*255 full_filename
      integer tmp1(2)

      logical value,psi_nogrid
      integer hversion,hnfft(3),hispin,hne(2)
      real*8 hunita(3,3)
      integer ind
      character*50 filename
      character*50 control_input_psi
      external     control_input_psi
      logical  wvfnc_expander
      external wvfnc_expander

  


*     **** external functions ****
      real*8      psp_zv,psp_rc,ewald_rcut,ion_amass
      real*8      ewald_mandelung
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg
      integer     ewald_ncut,ewald_nshl3d
      integer     psp_lmmax,psp_lmax,psp_locp
      character*4 ion_aname,ion_atom
      external    psp_zv,psp_rc,ewald_rcut,ion_amass
      external    ewald_mandelung
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg
      external    ewald_ncut,ewald_nshl3d
      external    psp_lmmax,psp_lmax,psp_locp
      external    ion_aname,ion_atom


      real*8   control_tole,control_tolc,control_tolr,ion_rion
      external control_tole,control_tolc,control_tolr,ion_rion
      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      logical  control_read,control_move,ion_init,ion_q_FixIon
      external control_read,control_move,ion_init,ion_q_FixIon
      logical  ion_q_zFixIon
      external ion_q_zFixIon
 
      integer  pack_nwave_all
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave
      integer  ion_nion,ion_natm,ion_katm,ion_nkatm
      external pack_nwave_all
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave
      external  ion_nion,ion_natm,ion_katm,ion_nkatm

      character*12 control_boundry
      external     control_boundry

      logical      pspw_SIC,pspw_SIC_relaxed,pspw_qmmm_found
      logical      pspw_HFX,pspw_HFX_relaxed,psp_pawexist
      logical      psp_semicore,control_Mulliken
      real*8       psp_rcore,psp_ncore,psp_rlocal
      external     pspw_SIC,pspw_SIC_relaxed,pspw_qmmm_found
      external     pspw_HFX,pspw_HFX_relaxed,psp_pawexist
      external     psp_semicore,control_Mulliken
      external     psp_rcore,psp_ncore,psp_rlocal
      logical      control_check_charge_multiplicity
      external     control_check_charge_multiplicity
      real*8       nwpw_timing
      external     nwpw_timing
      integer      control_np_orbital,control_mapping,control_mapping1d
      external     control_np_orbital,control_mapping,control_mapping1d

      logical  control_translation,control_rotation,control_balance
      external control_translation,control_rotation,control_balance

      logical  Dneall_m_allocate,Dneall_m_free,control_parallel_io
      external Dneall_m_allocate,Dneall_m_free,control_parallel_io

      real*8   Dneall_m_value
      external Dneall_m_value
      character*9 ion_amm
      external    ion_amm
      logical  psp_U_psputerm,meta_found,ion_disp_on
      external psp_U_psputerm,meta_found,ion_disp_on
      integer  ion_nconstraints,ion_ndof
      external ion_nconstraints,ion_ndof



*                            |************|
*****************************|  PROLOGUE  |****************************
*                            |************|

      value = .true.
      pi = 4.0d0*datan(1.0d0)

      call nwpw_timing_init()
      call dcopy(50,0.0d0,0,E,1)


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
         write(luout,1030)
         write(luout,1010)
         write(luout,1035)
         write(luout,1010)
         write(luout,1040)
         write(luout,1010)
         write(luout,1041)
         write(luout,1042)
         write(luout,1043)
         write(luout,1010)
         write(luout,1000)
         call nwpw_message(1)
         write(luout,1110)
      end if
      
      value = control_read(1,rtdb)
      call Parallel2d_Init(control_np_orbital())
      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)

      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0
      mapping = control_mapping()

*     **** initialize psi_data ****
      call psi_data_init(20)

*     **** initialize D3dB data structure ****
      call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d

*     ***** Initialize double D3dB data structure ****
      if (control_version().eq.4) 
     >   call D3dB_Init(2,2*ngrid(1),2*ngrid(2),2*ngrid(3),mapping)


*     **** initialize lattice and packing data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_Init()
      call Pack_npack(1,npack1)      

      call D3dB_pfft_init()
      call ga_sync()


 
*     **** read ions ****
      value = ion_init(rtdb)
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)


*     **** allocate psp data structure and read in psedupotentials into it ****
      call psp_init()
      call psp_readall()
      if (psp_semicore(0)) call semicore_check()

*     **** initialize G,mask,ke,and coulomb data structures ****
      call ke_init()
      if (control_version().eq.3) call coulomb_init()
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()
      call phafac()
      if (control_version().eq.3) call ewald_init()

*     **** generate initial wavefunction if it does not exist ****
      if (.not.control_check_charge_multiplicity()) then
        call psi_new()
      end if

*     ***** allocate psi2,and psi1 wavefunctions ****
      call psi_get_ne_occupation(ispin,ne,smearoccupation)
      if (smearoccupation.gt.0) then
         fractional = .true.
      else
         fractional = .false.
      end if
      mapping1d = control_mapping1d()
      call Dne_init(ispin,ne,mapping1d)
      call Dneall_neq(neq)
      nemaxq = neq(1)+neq(2)
      
      value = MA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'psi2',psi2(2),psi2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'psi1',psi1(2),psi1(1))
      if (fractional) then
         value = value.and.
     >        MA_alloc_get(mt_dbl,(ne(1)+ne(2)),'occ1',occ1(2),occ1(1))
         value = value.and.
     >        MA_alloc_get(mt_dbl,(ne(1)+ne(2)),'occ2',occ2(2),occ2(1))
      end if
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)


*     *****  read initial wavefunctions into psi2  ****
      if (.not.rtdb_get(rtdb,'nwpw:psi_nogrid',
     >                  mt_log,1,psi_nogrid))
     >   psi_nogrid = .true.

      if (psi_nogrid) then

        call psi_get_header(hversion,hnfft,hunita,hispin,hne)

        if ( (hnfft(1).ne.control_ngrid(1)) .or.
     >       (hnfft(2).ne.control_ngrid(2)) .or.
     >       (hnfft(3).ne.control_ngrid(3)) ) then

        hnfft(1) = control_ngrid(1)
        hnfft(2) = control_ngrid(2)
        hnfft(3) = control_ngrid(3)
        call Parallel_taskid(taskid)

        call ga_sync()
        value = rtdb_parallel(.false.)
        call ga_sync()
        if (taskid.eq.MASTER) then

          filename =  control_input_psi()

          ind = index(filename,' ') - 1
          if (.not. rtdb_cput(rtdb,'xpndr:old_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: rtdb_cput failed', 0, RTDB_ERR)

          if (.not. rtdb_cput(rtdb,'xpndr:new_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: rtdb_cput failed', 0, RTDB_ERR)

          if (.not. rtdb_put(rtdb,'xpndr:ngrid',mt_int,3,hnfft))
     >     call errquit(
     >     'wvfnc_expander_input: rtdb_put failed', 0, RTDB_ERR)

          write(*,*)
          write(*,*) "Grid is being converted:"
          write(*,*) "------------------------"
          write(*,*)
          write(*,*) "To turn off automatic grid conversion:"
          write(*,*)
          write(*,*) "set nwpw:psi_nogrid .false."
          write(*,*)
          value = wvfnc_expander(rtdb)

        end if
        call ga_sync()
        value = rtdb_parallel(.true.)

      end if

      end if

*     *****  read psi2 wavefunctions ****
      call psi_read(ispin,ne,dcpl_mb(psi2(1)),
     >              smearoccupation,dbl_mb(occ2(1)))
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)
      nemax = ne(1)+ne(2)

*     **** Ortho Check ****
      if (psp_pawexist()) then
      call psp_overlap_S(ispin,neq,dcpl_mb(psi2(1)),dcpl_mb(psi1(1)))
      call Grsm_gg_trace(npack1,(neq(1)+neq(2)),
     >                        dcpl_mb(psi1(1)),
     >                        dcpl_mb(psi2(1)),
     >                        sum2)
      else
      call Grsm_gg_trace(npack1,(neq(1)+neq(2)),
     >                        dcpl_mb(psi2(1)),
     >                        dcpl_mb(psi2(1)),
     >                        sum2)
      end if
      call D1dB_SumAll(sum2)

      sum1 = dble(ne(1) + ne(2))
      if (dabs(sum2-sum1).gt.1.0d-10) then
         if (taskid.eq.MASTER)
     >      write(*,*) "Warning: Gram-Schmidt being performed on psi2"
         if (psp_pawexist()) then
            call Dneall_f_Sortho(0,dcpl_mb(psi2(1)),
     >                             dcpl_mb(psi1(1)),npack1)
         else
            call Dneall_f_GramSchmidt(0,dcpl_mb(psi2(1)),npack1)
         end if
c         call Dneall_f_ortho(0,dcpl_mb(psi2(1)),npack1)
c         call Grsm_g_MakeOrtho(npack1,ne(1),dcpl_mb(psi2(1)))
c         if (ispin.gt.1) then
c           call Grsm_g_MakeOrtho(npack1,ne(2),
c     >                           dcpl_mb(psi2(1) + ne(1)*npack1))
c         end if
      end if

*     **** allocate other variables *****
      value = MA_alloc_get(mt_dbl,(2*nemax),'eig',eig(2),eig(1))
      value = value.and.Dneall_m_allocate(0,hml)
      value = value.and.Dneall_m_allocate(0,lmd)

      value = value.and.
     >        MA_alloc_get(mt_dbl,(4*nfft3d),
     >                     'dn',dn(2),dn(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'Hpsi',Hpsi(2),Hpsi(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,nfft3d*(neq(1)+neq(2)),
     >                     'psir',psir(2),psir(1))

      if (.not. value) 
     >   call errquit('cpsdv5:out of heap memory',0,MA_ERR)


*     **** initialize two-electron Gaussian integrals ****
      if (psp_pawexist()) then
         call nwpw_gintegrals_init()
         call nwpw_gintegrals_set(control_move())
      end if

*     **** initialize SIC and HFX  ****
      call pspw_init_SIC(rtdb,ne)
      call pspw_init_HFX(rtdb,ispin,ne)

*     **** initialize DFT+U ****
      call psp_U_init()

*     **** initialize META GGA ****
      call nwpw_meta_gga_init()

*     **** initialize metadynamics ****
      call meta_initialize(rtdb)

*     **** initialize QM/MM ****
      call pspw_init_APC(rtdb)
      call pspw_qmmm_init(rtdb)

*     **** initialize frac_occ data structure ****
c      call frac_occ_init(rtdb,ispin,ne)

*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)


*                |**************************|
******************   summary of input data  **********************
*                |**************************|

      if (taskid.eq.MASTER) then
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
         if (control_parallel_io()) then
           write(luout,1119)
         else
           write(luout,1122)
         end if


         write(luout,1115)
         IF(control_move()) THEN
           write(luout,1120) 'yes'
         ELSE
           write(luout,1120) 'no'
         ENDIF
         write(luout,1121) control_boundry(),control_version()
         if (ispin.eq.1) write(luout,1130) 'restricted'
         if (ispin.eq.2) write(luout,1130) 'unrestricted'

         call v_bwexc_print(luout,control_gga())

         if (fractional) write(luout,1132)
         call pspw_print_SIC(luout)
         call pspw_print_HFX(luout)
         write(luout,1140)
         do ia = 1,ion_nkatm()
           call psp_print(ia)
c           write(luout,1150) ia,ion_atom(ia),
c     >                    psp_zv(ia),psp_lmax(ia)
c           write(luout,1152) psp_lmax(ia)
c           write(luout,1153) psp_locp(ia)
c           write(luout,1154) psp_lmmax(ia)
c           if (control_version().eq.4) write(luout,1156) psp_rlocal(ia)
c           if (psp_semicore(ia)) 
c     >         write(luout,1155) psp_rcore(ia),psp_ncore(ia)
c           write(luout,1151) (psp_rc(i,ia),i=0,psp_lmax(ia))
         end do

         icharge = -(ne(1)+ne(ispin))
         en(1)     = ne(1)
         en(ispin) = ne(ispin)
         if (fractional) then
            icharge = 0.0d0
            do ms=1,ispin
            en(ms) =0.0
            do i=n1(ms),n2(ms)
              icharge = icharge - (3-ispin)*dbl_mb(occ1(1)+i-1)
              en(ms) = en(ms) + dbl_mb(occ1(1)+i-1)
            end do
            end do
         end if
            
         do ia=1,ion_nkatm()
           icharge = icharge + ion_natm(ia)*psp_zv(ia)
         end do
         write(luout,1159) icharge

         write(luout,1160)
         write(luout,1170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
         write(luout,1180)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
           else if (ion_q_zFixIon(I)) then
           write(luout,1193) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0
           else
           write(luout,1190) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz
         write(luout,1211) ion_nconstraints(),ion_ndof()


         if (fractional) then
           write(luout,1219) en(1),en(ispin),' (   fractional)'
           write(luout,1221) ne(1),neq(1),
     >                   ne(ispin),neq(ispin),' (Fourier space)'
         else
           write(luout,1220) ne(1),neq(1),
     >                   ne(ispin),neq(ispin),' (Fourier space)'
           write(luout,1221) ne(1),neq(1),
     >                   ne(ispin),neq(ispin),' (Fourier space)'
         end if
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
         write(luout,1250) lattice_ecut(),ngrid(1),ngrid(2),ngrid(3),
     >                 pack_nwave_all(0),pack_nwave(0)
         write(luout,1251) lattice_wcut(),ngrid(1),ngrid(2),ngrid(3),
     >                 pack_nwave_all(1),pack_nwave(1)
         if (control_version().eq.3) then
         write(luout,1260) ewald_rcut(),ewald_ncut()
         write(luout,1261) ewald_mandelung()
         end if
         write(luout,1270)
         if (.not.control_translation()) write(luout,1271)
         if (.not.control_rotation())    write(luout,1272)
         write(luout,1280) control_time_step(),control_fake_mass()
         write(luout,1290) control_tole(),control_tolc(),control_tolr()
        write(luout,1281) control_it_in()*control_it_out(),
     >                 control_it_in(),control_it_out()
         write(luout,1300)
         write(luout,1305)
         call util_flush(luout)
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
         call inner_loop(ispin,neq,
     >             npack1,nfft3d,nemaxq,
     >             dcpl_mb(psi1(1)),dcpl_mb(psi2(1)),dbl_mb(dn(1)),
     >             it_in,E,deltae,deltac,deltar,
     >             dbl_mb(hml(1)),dbl_mb(lmd(1)),
     >             dcpl_mb(psir(1)),dcpl_mb(Hpsi(1)),
     >             fractional,dbl_mb(occ1(1)),dbl_mb(occ2(1)))
         if (taskid.eq.MASTER) then 
           write(luout,1310) icount*it_in,E(1),deltae,deltac,deltar
           call util_flush(luout)
         end if
         if ((deltae.gt.0.0d0).and.(icount.gt.1)) then
            if (taskid.eq.MASTER) 
     >       write(luout,*) 
     >       ' *** Energy going up.  iteration terminated.'
            go to 2
         end if
         deltae = dabs(deltae)
         if ((deltae.lt.control_tole()).and.
     >       (deltac.lt.control_tolc()).and.
     >       (deltar.lt.control_tolr())) then
            if (taskid.eq.MASTER) 
     >       write(luout,*) 
     >       ' *** tolerance ok.     iteration terminated.'
            go to 2
         end if
      if (icount.lt.it_out) go to 1
      if (taskid.eq.MASTER) 
     > write(luout,*) 
     > '*** arived at the Maximum iiteration.   terminated.'

*::::::::::::::::::::  end of iteration loop  :::::::::::::::::::::::::

   2  continue
      if (taskid.eq.MASTER) CALL nwpw_MESSAGE(3)
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
         en1(ms) = sumall*lattice_omega()
     >             /dble(ngrid(1)*ngrid(2)*ngrid(3))
      end do

      if (psp_pawexist()) then
         if (.not.MA_push_get(mt_dbl,n2ft3d,'tmp1',tmp1(2),tmp1(1)))
     >   call errquit(
     >        'cgsdv5: out of stack memory',0,MA_ERR)

         call psp_qlm_atom(ispin,neq,dcpl_mb(psi1(1)))
         do ms=1,ispin
           call nwpw_compcharge_gen_dn_cmp_smooth_ms(ms,dbl_mb(tmp1(1)))
           call Pack_c_unpack(0,dbl_mb(tmp1(1)))
           call D3dB_cr_fft3b(1,dbl_mb(tmp1(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))
           call D3dB_r_dsum(1,dbl_mb(tmp1(1)),sumall)
           en2(ms) = sumall*lattice_omega()
     >              /dble(ngrid(1)*ngrid(2)*ngrid(3))
         end do
         if (.not.MA_pop_stack(tmp1(2)))
     >   call errquit(
     >        'cgsdv5: popping stack memory',0,MA_ERR)
      else
         en2(1) = 0.0d0
         en2(2) = 0.0d0
      end if
      en(1) = en1(1)+en2(1)
      en(2) = en1(2)+en2(2)
      if (taskid.eq.MASTER) then
         write(17,1320) (en(ms),ms=1,ispin)
         if (psp_pawexist()) then
            write(17,1322) (en1(ms),ms=1,ispin)
            write(17,1323) (en2(ms),ms=1,ispin)
         end if
      end if



*     **** comparison between hamiltonian an lambda matrix ****
      if (taskid.eq.MASTER) write(17,1330)
      do ms=1,ispin
         do j=1,ne(ms)
         do i=1,ne(ms)
            w      = Dneall_m_value(0,ms,i,j,dbl_mb(hml(1)))
            virial = Dneall_m_value(0,ms,i,j,dbl_mb(lmd(1)))

            if (taskid.eq.MASTER)
     >      write(17,1340) ms,i,j,w,virial,w-virial
         end do
         end do
      end do


*     **** check orthonormality ****
      if (taskid.eq.MASTER) then
         write(17,1350)
      end if
      call Dneall_ffm_Multiply(0,dcpl_mb(psi1(1)),
     >                           dcpl_mb(psi1(1)),npack1,
     >                           dbl_mb(lmd(1)))
      do ms=1,ispin
         do j=1,ne(ms)
         do i=j,ne(ms)
            w  = Dneall_m_value(0,ms,i,j,dbl_mb(lmd(1)))
            if (taskid.eq.MASTER) write(17,1360) ms,i,j,w
         end do
         end do
      end do

*     **** close check file ****
      if (taskid.eq.MASTER) then
         close(17)
      end if

*     ***** do not diagonalize the hamiltonian matrix *****
      if (pspw_SIC()) then
        call dcopy(2*npack1*nemaxq,
     >             dcpl_mb(psi1(1)),1,
     >             dcpl_mb(psi2(1)),1)

*     ***** diagonalize the hamiltonian matrix but don't rotate ****
      else if (fractional) then
         n = ne(1)
c         call Dneall_m_HmltimesSA(0,dbl_mb(hml(1)),dbl_mb(fweight(1)))
         call Dneall_m_diagonalize(0,dbl_mb(hml(1)),
     >                              dbl_mb(eig(1)),.false.)
c         do ms=1,ispin
c         do ii=1,ne(ms)
c            dbl_mb(eig(1)+(ii-1)+(ms-1)*n)
c     >      =dbl_mb(eig(1)+(ii-1)+(ms-1)*n)
c     >      /dbl_mb(fweight(1)+(ii-1)+(ms-1)*n)
c         end do
c         end do
         call dcopy(2*npack1*nemaxq,
     >             dcpl_mb(psi1(1)),1,
     >             dcpl_mb(psi2(1)),1)

*     ***** diagonalize and rotate the hamiltonian matrix ****
      else
        call Dneall_m_diagonalize(0,dbl_mb(hml(1)),
     >                              dbl_mb(eig(1)),.false.)

        call Dneall_fmf_Multiply(0,dcpl_mb(psi1(1)),npack1,
     >                            dbl_mb(hml(1)),  1.0d0,
     >                            dcpl_mb(psi2(1)),0.0d0)
      end if

*                |***************************|
****************** report summary of results **********************
*                |***************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

      if (taskid.eq.MASTER) then
         write(luout,1300)
         write(luout,1410)
         write(luout,1420)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
           else if (ion_q_zFixIon(I)) then
           write(6,1193) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           else
           write(luout,1190) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz
         write(luout,1211) ion_nconstraints(),ion_ndof()


         write(luout,*)
         write(luout,1320) en(1),en(ispin),' (real space)'
         if (psp_pawexist()) then
            write(luout,1322) en1(1),en1(ispin),' (real space)'
            write(luout,1323) en2(1),en2(ispin),' (real space)'
         end if


         write(luout,1430) E(1),E(1)/ion_nion()
         if (pspw_qmmm_found()) then
            write(luout,1431) 
            write(luout,1432) 
            write(luout,1433) (E(1)-E(11)),
     >                    (E(1)-E(11))/ion_nion()
         end if
         write(luout,1440) E(2),E(2)/n2(ispin)
         write(luout,1450) E(3),E(3)/n2(ispin)
         write(luout,1460) E(4),E(4)/n2(ispin)
         if (pspw_SIC()) then
           write(luout,1455) E(16),E(16)/n2(ispin)
           write(luout,1456) E(17),E(17)/n2(ispin)
         end if
         if (pspw_HFX()) then
           write(luout,1457) E(20),E(20)/n2(ispin)
         end if
         if (psp_U_psputerm()) then
           write(luout,1458) E(29),E(29)/n2(ispin)
         end if
         if (meta_found()) then
           write(luout,1459) E(31),E(31)/ion_nion()
         end if


         write(luout,1470) E(5),E(5)/ion_nion()
         write(luout,1480) E(6),E(6)/n2(ispin)
         write(luout,1490) E(7),E(7)/n2(ispin)
         write(luout,1495) E(8),E(8)/n2(ispin)
         write(luout,1496) E(9),E(9)/n2(ispin)
         write(luout,1497) E(10),E(10)/n2(ispin)
         if (pspw_SIC().and.pspw_SIC_relaxed())  then
           write(luout,1499) E(18),E(18)/n2(ispin)
           write(luout,1501) E(19),E(19)/n2(ispin)
         end if
         if (pspw_HFX().and.pspw_HFX_relaxed())  then
           write(luout,1502) E(21),E(21)/n2(ispin)
         end if
         if (psp_U_psputerm()) then
           write(luout,1503) E(30),E(30)/n2(ispin)
         end if
         if (meta_found()) then
           write(luout,1504) E(32),E(32)/n2(ispin)
         end if
         virial = (E(10)+E(9)+E(8)+E(7))/E(6)
         write(luout,1498) virial

c        **** paw energies ****
         if (psp_pawexist()) then
         write(luout,1680) E(35),E(35)/n2(ispin)
         write(luout,1681) E(36),E(36)/ion_nion()
         write(luout,1682) E(37),E(37)/n2(ispin)

         write(luout,1683) E(38),E(38)/n2(ispin)
         write(luout,1684) E(39),E(39)/n2(ispin)

         write(luout,1685) E(40),E(40)/n2(ispin)
         write(luout,1686) E(41),E(41)/n2(ispin)
         write(luout,1687) E(42),E(42)/n2(ispin)
         write(luout,1688) E(43),E(43)/n2(ispin)
         write(luout,1689) E(44),E(44)/ion_nion()
         write(luout,1690) E(45),E(45)/ion_nion()
         end if

 
      if (pspw_qmmm_found()) then
            write(luout,1700)
            write(luout,1701)
            write(luout,1702) E(11)
            write(luout,1703) E(12)
            write(luout,1704) E(13)
c
c            write(luout,1703) E(14)+E(15)
c            write(luout,1710) E(14)
c            write (luout,1711) E(15)
c
c            write(luout,1704) E(16)
c            write(luout,1705) E(17)
c            write(luout,1706) E(18)
        end if
        if (ion_disp_on()) then
            write(luout,1720) E(33)
        end if


*        **** write out diagonal <psi|H|psi> matrix ****
         if (pspw_SIC()) then

          n = ne(1)
          nn = n*n
          do ms=1,ispin
             if (ms.eq.1) write(luout,1331)
             if (ms.eq.2) write(luout,1332)
             !*** call Gainsville matrix output ***
             call output(dbl_mb(hml(1)+(ms-1)*nn),
     >                    1,ne(ms),1,ne(ms),
     >                    n,n,1)

           end do

*        **** write out KS eigenvalues ****
         else
         write(luout,1500)
         NN=NE(1)-NE(2)
         EV=27.2116d0
         if (fractional) then
           do i=1,NN
             write(luout,1511) dbl_mb(EIG(1)+i-1),
     >                     dbl_mb(EIG(1)+i-1)*EV,
     >                     dbl_mb(occ1(1)+i-1)
           end do
           do i=1,ne(2)
             write(luout,1511) dbl_mb(EIG(1)+i-1+NN),
     >                     dbl_mb(EIG(1)+i-1+NN)*EV,
     >                     dbl_mb(occ1(1)+i-1+NN),
     >                     dbl_mb(EIG(1)+i-1+n1(2)-1),
     >                     dbl_mb(EIG(1)+i-1+n1(2)-1)*EV,
     >                     dbl_mb(occ1(1)+i-1+n1(2)-1)
           end do
         else
           do i=1,NN
             write(luout,1510) dbl_mb(EIG(1)+i-1),dbl_mb(EIG(1)+i-1)*EV
           end do
           do i=1,ne(2)
             write(luout,1510) dbl_mb(EIG(1)+i-1+NN),
     >                     dbl_mb(EIG(1)+i-1+NN)*EV,
     >                     dbl_mb(EIG(1)+i-1+n1(2)-1),
     >                     dbl_mb(EIG(1)+i-1+n1(2)-1)*EV
           end do
         end if

         end if

*        ***** extra energy output for QA test ****
         write(luout,1600) E(1)
      end if

*                |***************************|
******************         Prologue          **********************
*                |***************************|

*     **** calculate spin contamination ****
      call Calculate_psi_spin2(ispin,ne,npack1,dcpl_mb(psi2(1)),
     >                         fractional,dbl_mb(occ2(1)),w)

*     **** calculate the Dipole ***
      call Calculate_Dipole(ispin,ne,n2ft3d,dbl_mb(dn(1)),dipole)
      
*     **** perfom Lubin and Mulliken analysis ***
      if (control_Mulliken()) then

*       **** Lubin Water Analysis ***
        call pspw_Lubin_water_analysis(rtdb,ispin,ne,n2ft3d,
     >                                 dbl_mb(dn(1)))

*       **** Analysis ***
        call pspw_analysis(0,rtdb,ispin,ne,dcpl_mb(psi2(1)),
     >                                   dbl_mb(eig(1)))

*       **** generate APC *****
        call pspw_dngen_APC(ispin,ne,dbl_mb(dn(1)))
        call pspw_print_APC(6)

      end if

*     ***** write psi2 wavefunctions ****
      call psi_write(ispin,ne,dcpl_mb(psi2(1)),
     >               smearoccupation,dbl_mb(occ2(1)))
 
*     **** write geometry to rtdb ****
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
      call psp_U_end()
      call nwpw_meta_gga_end()
      call pspw_end_SIC()
      call pspw_end_HFX()
      call pspw_end_APC()
      call pspw_qmmm_end()
      call meta_finalize(rtdb)
      if (psp_pawexist()) call nwpw_gintegrals_end()

c      call frac_occ_end()

      value = MA_free_heap(psir(2))
      value = MA_free_heap(Hpsi(2))
      value = MA_free_heap(dn(2))
      value = MA_free_heap(eig(2))
      value = Dneall_m_free(hml)
      value = Dneall_m_free(lmd)
      value = MA_free_heap(psi1(2))
      value = MA_free_heap(psi2(2))
      if (fractional) then
         value = MA_free_heap(occ1(2))
         value = MA_free_heap(occ2(2))
      end if
      call D3dB_pfft_end()
      call D3dB_end(1)
      if (control_version().eq.4) call D3dB_end(2)
      call Dne_end()
      call psi_data_end()

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
         write(luout,*)
         write(luout,*) '-----------------'
         write(luout,*) 'cputime in seconds'
         write(luout,*) 'prologue    : ',T1
         write(luout,*) 'main loop   : ',T2
         write(luout,*) 'epilogue    : ',T3
         write(luout,*) 'total       : ',T4
         write(luout,*) 'cputime/step: ',AV
         write(luout,*)
         call nwpw_timing_print_final(.true.,(icount*it_in))
         CALL nwpw_MESSAGE(4)
      end if 


      call Parallel2d_Finalize()
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
 1111 FORMAT(/' number of processors used:',I10)
 1112 FORMAT( ' parallel mapping         :      1d slab')
 1113 FORMAT( ' parallel mapping         :   2d hilbert')
 1114 FORMAT( ' parallel mapping         :     balanced')
 1115 FORMAT(/' options:')
 1116 FORMAT( ' parallel mapping         : not balanced')
 1117 FORMAT( ' processor grid           :',I4,' x',I4)
 1118 FORMAT( ' parallel mapping         :    2d hcurve')
 1119 FORMAT( ' parallel io              :        on')
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1122 FORMAT( ' parallel io              :       off')
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
 1170 FORMAT(7(5X,A2,':',I5))
 1180 FORMAT(/' initial position of ions:')
 1190 FORMAT(5X, I4, A5, ' (',3F11.5,' ) - atomic mass= ',F7.3,' ',A)
 1191 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - fixed ',A)
 1193 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - z fixed')
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1211 FORMAT(5X,'   number of constraints = ', I6,' ( DOF = ',I6,' )' )
 1219 FORMAT(/' number of electrons: spin up=',F6.2, 16x,
     >                               '  down=',F6.2,A)
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
 1231 FORMAT(5x,' volume : ',F12.1)
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
 1261 FORMAT(5X,'                   madelung=',f11.8)
 1270 FORMAT(/' technical parameters:')
 1271 FORMAT(5x, ' translation constrained')
 1272 FORMAT(5x, ' rotation constrained')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'ficticious mass=',F10.1)
 1281 FORMAT(5X, ' maximum iterations =',I10,
     >           ' ( ',I4,' inner ',I6,' outer )')
 1290 FORMAT(5X, ' tolerance=',E8.3,' (energy)',E12.3,
     &        ' (electron)',E12.3,' (ion)')
 1300 FORMAT(//)
 1305 FORMAT(10X,'================ iteration =========================')
 1310 FORMAT(I8,E20.10,3E15.5)
 1320 FORMAT(' number of electrons: spin up=',F11.5,'  down=',F11.5,A)
 1321 FORMAT(' total charge of system:',F11.5,A)
 1322 FORMAT('     plane-wave part:         ',F11.5,'       ',F11.5,A)
 1323 FORMAT('      augmented part:         ',F11.5,'       ',F11.5,A)

 1330 FORMAT(/' comparison between hamiltonian and lambda matrix')
 1331 FORMAT(/' Elements of Hamiltonian matrix (up/restricted)')
 1332 FORMAT(/' Elements of Hamiltonian matrix (down)')
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
 1457 FORMAT( ' HF exchange energy  :',E19.10,' (',E15.5,'/electron)')
 1458 FORMAT( ' DFT+U     energy    :',E19.10,' (',E15.5,'/electron)')
 1459 FORMAT( ' Metadynamics energy :',E19.10,' (',E15.5,'/ion)')
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
 1500 FORMAT(/' orbital energies:')
 1501 FORMAT( ' K.S. SIC-exc-corr energy :',E19.10,
     >        ' (',E15.5,'/electron)')
 1502 FORMAT( ' K.S. HFX energy     :',E19.10,
     >        ' (',E15.5,'/electron)')
 1503 FORMAT( ' K.S. DFT+U energy   :',E19.10,' (',E15.5,'/electron)')
 1504 FORMAT( ' K.S. Metadynamics energy :',E19.10,' (',E15.5,'/ion)')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1511 FORMAT(2(E18.7,' (',F8.3,'eV) occ=',F5.3))
 1600 FORMAT(/' Total PSPW energy   :',E19.10)

 1680 FORMAT(/' kinetic (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1681 FORMAT( ' kinetic       (core):',E19.10,' (',E15.5,'/ion)')
 1682 FORMAT( ' valence-core        :',E19.10,' (',E15.5,'/electron)')
 1683 FORMAT( ' ncmp*V_local        :',E19.10,' (',E15.5,'/electron)')
 1684 FORMAT( ' V_lpaw  (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1685 FORMAT( ' coulomb (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1686 FORMAT( ' coulomb  (ncmp-ncmp):',E19.10,' (',E15.5,'/electron)')
 1687 FORMAT( ' coulomb    (ncmp-pw):',E19.10,' (',E15.5,'/electron)')
 1688 FORMAT( ' exc-cor (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1689 FORMAT( ' pxc-cor (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1690 FORMAT( ' ion-core            :',E19.10,' (',E15.5,'/ion)')
 1691 FORMAT( ' V_local PAW residual:',E19.10,' (',E15.5,'/electron)')

 1700 FORMAT(/' QM/MM-pol-vib/CAV Energies')
 1701 FORMAT( ' --------------------------')
 1702 FORMAT( ' LJ energy              :',E19.10)
 1703 FORMAT( ' Residual Coulomb energy:',E19.10)
 1704 FORMAT( ' MM Vibration energy    :',E19.10)
 1705 FORMAT( ' MM Vibration energy    :',E19.10)
 1706 FORMAT( ' (QM+MM)/Cavity energy  :',E19.10)

 1720 FORMAT(/' Dispersion energy   :',E19.10)

 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(luout,9010) ierr
      call Parallel2d_Finalize()
      call Parallel_Finalize()

      cpsdv5 = value
      return
      END
