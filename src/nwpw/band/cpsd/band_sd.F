*
* $Id$
*

***********************************************************************
*                      band_sd                                        *
*                                                                     *
*     This is a developing steepest descent code for BAND             *
*                                                                     *
*                                                                     *
***********************************************************************

      logical function band_sd(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "bafdecls.fh"
#include "inp.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
      

*     **** parallel variables ****
      integer  taskid,np,np_i,np_j,np_k
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,nfft3d
      integer npack1,npack0

*     **** electronic variables ****
      logical spin_orbit,newpsi
      real*8 icharge
      integer ispin,ispinq,nbrillioun,nbrillq
      integer ne(2),n1(2),n2(2),nemax,neq(2),nemaxq,neall
      real*8  en(2)
      real*8  dipole(3)

      integer psi1_tag,psi2_tag,psir_tag,Hpsi_tag,next
      integer psi1_shift,psi2_shift
      integer dn(2)
    

*     ***** energy variables ****
      real*8  E(50)

      integer eig_tag,hml_tag,lmd_tag,svec_tag
      integer eig_shift,hml_shift,lmd_shift,svec_shift

*     **** psi smearing block ****
      logical fractional
      integer smearoccupation,smeartype
      real*8 smearfermi(2),smearcorrection,smearkT


*     **** error variables ****
      integer ierr

*     **** local variables ****
      complex*16 zero,one
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))

      logical lprint,mprint,hprint
      integer ms,mapping,mapping1d
      real*8  deltae,deltac,deltar
      real*8  gx,gy,gz,cx,cy,cz,sum1,sum2
      real*8  EV,pi,f0,f1,f2,f3,f4,f5,f6
      integer i,j,k,ia,n,nn,nb
      integer ii,jj,indx,indx1,if1,if2
      integer icount,it_in,it_out
      real*8 w,sumall,virial,a,b,c,alpha,beta,gamma
      integer nfft3
      parameter (nfft3=32)
      character*255 full_filename

      logical value,psi_nogrid
      integer hversion,hnfft(3),hispin,hne(2)
      real*8 hunita(3,3)
      integer ind,vers

      character*255 filename
      character*50 control_input_psi
      external     control_input_psi
      logical  c_wvfnc_expander
      external c_wvfnc_expander


*     **** external functions ****
      real*8      cpsp_zv,cpsp_rc,ewald_rcut,ion_amass,ion_TotalCharge
      real*8      ewald_mandelung
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg,brillioun_weight
      integer     ewald_ncut,ewald_nshl3d
      integer     cpsp_nprj,cpsp_lmax,cpsp_locp,cpsp_psp_type
      character*4 ion_aname,ion_atom
      external    cpsp_zv,cpsp_rc,ewald_rcut,ion_amass,ion_TotalCharge
      external    ewald_mandelung
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg,brillioun_weight
      external    ewald_ncut,ewald_nshl3d
      external    cpsp_nprj,cpsp_lmax,cpsp_locp,cpsp_psp_type
      external    ion_aname,ion_atom

      logical  control_fractional
      external control_fractional
      real*8   control_fractional_temperature
      external control_fractional_temperature
      real*8   control_fractional_kT,control_fractional_smeartype
      external control_fractional_kT,control_fractional_smeartype
      real*8   control_tole,control_tolc,control_tolr,ion_rion
      external control_tole,control_tolc,control_tolr,ion_rion
      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      logical  control_read,control_move,ion_init,ion_q_FixIon
      external control_read,control_move,ion_init,ion_q_FixIon
      logical  ion_q_xyzFixIon,band_HFX,band_HFX_relaxed
      external ion_q_xyzFixIon,band_HFX,band_HFX_relaxed
      character*14 ion_q_xyzFixIon_label
      external     ion_q_xyzFixIon_label

      logical  brillioun_print
      external brillioun_print
      real*8   brillioun_k_brdcst,brillioun_ks_brdcst
      external brillioun_k_brdcst,brillioun_ks_brdcst
      real*8   brillioun_weight_brdcst
      external brillioun_weight_brdcst
      real*8   cpsi_eig_brdcst_tag,cpsi_occ_brdcst_tag
      external cpsi_eig_brdcst_tag,cpsi_occ_brdcst_tag
      real*8   cpsi_sv_brdcst_tag
      external cpsi_sv_brdcst_tag

      integer  Cram_nwave_brdcst,Cram_nwave_all_brdcst
      external Cram_nwave_brdcst,Cram_nwave_all_brdcst
 
      integer  psi_get_version,brillioun_nbrillioun
      integer  Pneb_ispinq,Pneb_nbrillq,Pneb_w_size
      integer  cpsi_data_alloc,cpsi_data_get_next,cpsi_data_get_chnk
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave
      integer  ion_nion,ion_natm,ion_katm,ion_nkatm
      external psi_get_version,brillioun_nbrillioun
      external Pneb_ispinq,Pneb_nbrillq,Pneb_w_size
      external cpsi_data_alloc,cpsi_data_get_next,cpsi_data_get_chnk
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave
      external  ion_nion,ion_natm,ion_katm,ion_nkatm

      character*12 control_boundry
      external     control_boundry

      logical      pspw_reformat_c_wvfnc
      logical      pspw_HFX,pspw_HFX_relaxed
      logical      cpsp_semicore,control_Mulliken
      real*8       cpsp_rcore,cpsp_ncore
      external     pspw_reformat_c_wvfnc
      external     pspw_HFX,pspw_HFX_relaxed
      external     cpsp_semicore,control_Mulliken
      external     cpsp_rcore,cpsp_ncore
      logical      control_check_charge_multiplicity
      external     control_check_charge_multiplicity
      real*8       nwpw_timing
      external     nwpw_timing
      integer      control_np_dimensions
      integer      control_mapping,control_mapping1d
      external     control_np_dimensions
      external     control_mapping,control_mapping1d

      logical  control_print
      logical  control_translation,control_rotation,control_balance
      external control_print
      external control_translation,control_rotation,control_balance

      logical  Dneall_m_allocate,Dneall_m_free,ion_disp_on
      external Dneall_m_allocate,Dneall_m_free,ion_disp_on

      complex*16 Pneb_w_value
      external   Pneb_w_value
      character*9 ion_amm
      external    ion_amm

      character*255 cpsp_comment,comment
      external      cpsp_comment
      integer  ion_nconstraints,ion_ndof
      external ion_nconstraints,ion_ndof
      logical  ion_makehmass2,control_only_lda
      external ion_makehmass2,control_only_lda



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
      
      if (.not.control_read(13,rtdb))
     >   call errquit('band_sd:error reading control',0,DISK_ERR)


      lprint = ((taskid.eq.MASTER).and.(control_print(print_low)))
      mprint = ((taskid.eq.MASTER).and.(control_print(print_medium)))
      hprint = ((taskid.eq.MASTER).and.(control_print(print_high)))

      if (taskid.eq.MASTER) call current_second(cpu1)

*     ***** print out header ****
      if (mprint) then
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

      call Parallel3d_Init(control_np_dimensions(2),
     >                     control_np_dimensions(3))
      call Parallel3d_np_i(np_i)
      call Parallel3d_np_j(np_j)
      call Parallel3d_np_k(np_k)

      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0
      mapping = control_mapping()

*     **** initialize D3dB data structure ****
      call C3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call C3dB_nfft3d(1,nfft3d)

*     **** initialize psi_data ****
      call cpsi_data_init(20)

*     **** read ions ****
      value = ion_init(rtdb)
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

*     **** initialize lattice and packing data structure ****
      call lattice_init()
      call c_G_init()
      call brillioun_init()
      call Cram_Init()
      call C3dB_pfft_init()

*     ***** Initialize double D3dB data structure ****
      !if ((control_gga().ge.10).and.(control_gga().le.200)) then
      if (.not.control_only_lda()) then
         call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
         call G_init()
         call mask_init()
      end if

c*     **** read ions ****
c      value = ion_init(rtdb)
c      call center_geom(cx,cy,cz)
c      call center_mass(gx,gy,gz)

*     **** allocate psp data structure and read in psedupotentials into it ****
      call cpsp_init()
      call cpsp_readall()
      if (cpsp_semicore(0)) call c_semicore_check()


*     **** initialize G,mask,ke,and coulomb data structures ****
      call cstrfac_init()
      call cke_init()
      call c_coulomb_init()
      call ewald_init()


*     **** generate initial wavefunction if it does not exist ****
      if (.not.control_check_charge_multiplicity()) then
        call cpsi_new()
        newpsi = .true.
      else
         newpsi = .false.

*        **** convert from pspw format to band format ****
         vers = psi_get_version()
         if ((vers.eq.3).or.(vers.eq.4)) then
           newpsi = .true.
           value = btdb_parallel(.false.)
           if (taskid.eq.MASTER) then
             value= pspw_reformat_c_wvfnc(1)
           end if
           value = btdb_parallel(.true.)
         end if
      end if

      call psi_get_ne(ispin,ne)
      if (ispin.eq.3) then
         spin_orbit = .true.
         ispin=2
      else
         spin_orbit = .false.
      end if
      nbrillioun = brillioun_nbrillioun()
      mapping1d = control_mapping1d()
      call Pneb_init(ispin,ne,nbrillioun,mapping1d,spin_orbit)
      call Pneb_neq(neq)


*     ***** allocate psi2,and psi1 wavefunctions ****
      call psi_get_ne_occupation(ispin,ne,smearoccupation)
      if (smearoccupation.gt.0) then
         fractional = .true.
      else
         fractional = .false.
      end if
      ispinq  = Pneb_ispinq()
      nbrillq = Pneb_nbrillq()

      call Cram_npack(0,npack0)
      call Cram_max_npack(npack1)
      nemaxq = neq(1)+neq(2)
      neall  = ne(1) +ne(2)

      psi2_tag = cpsi_data_alloc(nbrillq,nemaxq,2*npack1)
      psi1_tag = cpsi_data_alloc(nbrillq,nemaxq,2*npack1)

*     *** fractional orbitals ***
      if (smearoccupation.gt.0) then
        call cpsi_data_set_next(psi2_tag,
     >                 cpsi_data_alloc(nbrillq,nemaxq,1))
        call cpsi_data_set_next(psi1_tag,
     >                 cpsi_data_alloc(nbrillq,nemaxq,1))
        smeartype = control_fractional_smeartype()
        smearkT   = control_fractional_kT()
      end if

c     **** allocate other variables ****
      Hpsi_tag = cpsi_data_alloc(nbrillq,nemaxq,2*npack1)
      psir_tag = cpsi_data_alloc(nbrillq,nemaxq,2*nfft3d)
      hml_tag  = cpsi_data_alloc(nbrillq,1,2*Pneb_w_size(0,1))
      eig_tag  = cpsi_data_alloc(nbrillq,ne(1)+ne(2),1)
      svec_tag = cpsi_data_alloc(nbrillq,neq(1),3)
      value = BA_alloc_get(mt_dbl,2*nfft3d,'dn',dn(2),dn(1))
      if (.not. value)
     >  call errquit('band_sd:out of heap memory',0,MA_ERR)


*     *****  read initial wavefunctions into psi2  ****
      if (.not.btdb_get(rtdb,'nwpw:psi_nogrid',
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
        value = btdb_parallel(.false.)
        call ga_sync()
        if (hprint) then

          filename =  control_input_psi()

          ind = index(filename,' ') - 1
          if (.not. btdb_cput(rtdb,'c_xpndr:old_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'c_wvfnc_expander_input: btdb_cput failed', 0, RTDB_ERR)

          if (.not. btdb_cput(rtdb,'c_xpndr:new_wavefunction_filename',
     >                    1,filename(1:ind)))
     >     call errquit(
     >     'c_wvfnc_expander_input: btdb_cput failed', 0, RTDB_ERR)

          if (.not. btdb_put(rtdb,'c_xpndr:ngrid',mt_int,3,hnfft))
     >     call errquit(
     >     'c_wvfnc_expander_input: btdb_put failed', 0, RTDB_ERR)

          write(*,*)
          write(*,*) "Grid is being converted:"
          write(*,*) "------------------------"
          write(*,*)
          write(*,*) "To turn off automatic grid conversion:"
          write(*,*)
          write(*,*) "set nwpw:psi_nogrid .false."
          write(*,*)
          value = c_wvfnc_expander(rtdb)

        end if
        call ga_sync()
        value = btdb_parallel(.true.)

        end if
      end if

*     *****  read psi2 wavefunctions ****
      call cpsi_read(spin_orbit,ispin,ne,nbrillioun,psi2_tag)

*     **** Ortho Check ****
      call Pneb_orthoCheckMake_tag(.true.,0,0,npack1,psi2_tag)

*     **** initialize SIC and HFX  ****
c      call pspw_init_SIC(rtdb,ne)
       call band_init_HFX(rtdb,nbrillioun,ispin,ne)

*     **** initialize rho_symmetry  ****
      call c_rho_symmetrizer_init()

*     **** initialize QM/MM ****
c      call pspw_init_APC(rtdb)
c      call pspw_qmmm_init(rtdb)


*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)


*                |**************************|
******************   summary of input data  **********************
*                |**************************|


*     **** determine en ****
      if (.not.spin_orbit) then
        next    = cpsi_data_get_next(psi2_tag)
        icharge = 0.0d0
        en(1)   = 0.0d0
        en(2)   = 0.0d0
        b = dble(3-ispin)
        do nb=1,nbrillq
        w = brillioun_weight(nb)
        do ms=1,ispin
          do i=1,ne(ms)
            if (next.lt.0) then
               a = 1.0d0
            else
               indx = cpsi_data_get_chnk(next,nb)+i-1
               if (.not.spin_orbit) indx = indx + (ms-1)*ne(1)
               a = dbl_mb(indx)
            end if
            icharge = icharge - b*a*w
            en(ms)  = en(ms) + a*w
          end do
        end do
        end do
        call K1dB_Vector_SumAll(2,en)
        call K1dB_SumAll(icharge)
      else
        icharge   = -ne(1)
        en(1)     =  ne(1)
        en(ispin) =  ne(ispin)
      end if


      if (mprint) then
         write(luout,1111) np
         write(luout,1117) np_i,np_j,np_k
         if (mapping.eq.1) write(luout,1112)
         if (mapping.eq.2) write(luout,1113)
         if (mapping.eq.3) write(luout,1118)

         write(luout,1115)
         IF(control_move()) THEN
           write(luout,1120) 'yes'
         ELSE
           write(luout,1120) 'no'
         ENDIF
         write(luout,1121) control_boundry(),control_version()
         if (.not.spin_orbit) then
            if (ispin.eq.1) write(luout,1130) 'restricted'
            if (ispin.eq.2) write(luout,1130) 'unrestricted'
         else
            write(luout,1130) 'spin orbit'
         end if

         call v_bwexc_print(luout,control_gga())

c         if (fractional) write(6,1132)
c         call pspw_print_SIC(6)
c         call pspw_print_HFX(6)
         if (ion_makehmass2()) write(luout,1135)
         write(luout,1140)
         do ia = 1,ion_nkatm()
           write(luout,1150) ia,ion_atom(ia),
     >                    cpsp_zv(ia),cpsp_lmax(ia)

           comment = cpsp_comment(ia)
           i = inp_strlen(comment)
           write(luout,1157) comment(1:i)
           write(luout,1158) cpsp_psp_type(ia)
           write(luout,1152) cpsp_lmax(ia)
           write(luout,1153) cpsp_locp(ia)
           write(luout,1154) cpsp_nprj(ia)
           if (cpsp_semicore(ia)) 
     >         write(luout,1155) cpsp_rcore(ia),cpsp_ncore(ia)
           write(luout,1151) (cpsp_rc(i,ia),i=0,cpsp_lmax(ia))
         end do


         icharge = icharge + ion_TotalCharge()
         write(6,1159) icharge

         write(luout,1160)
         write(luout,1170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
         write(luout,1180)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
           else if (ion_q_xyzFixIon(I)) then
           write(luout,1194) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_q_xyzFixIon_label(I)
           else
           write(luout,1190) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz
         write(luout,1211) ion_nconstraints(),ion_ndof()

         write(luout,1220) en(1),en(ispin),' (Fourier space)'
         write(luout,1221) ne(1),neq(1),
     >                     ne(ispin),neq(ispin),' (Fourier space)'

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
         call lattice_abc_abg(a,b,c,alpha,beta,gamma)
         write(luout,1232) a,b,c,alpha,beta,gamma
         write(luout,1260) ewald_rcut(),ewald_ncut()
         write(luout,1261) ewald_mandelung()

         write(luout,1255)
         write(luout,1256) brillioun_nbrillioun()
      end if


c         write(6,1250) lattice_ecut(),ngrid(1),ngrid(2),ngrid(3),
c     >                 pack_nwave_all(0),pack_nwave(0)
c         write(6,1251) lattice_wcut(),ngrid(1),ngrid(2),ngrid(3),
c     >                 pack_nwave_all(1),pack_nwave(1)

c     **** print brillioun zone - extra logic for distributed kpoints ****
      if (brillioun_print()) then
         do i=1,brillioun_nbrillioun()
            f0 = brillioun_weight_brdcst(i)
            f1 = brillioun_ks_brdcst(1,i)
            f2 = brillioun_ks_brdcst(2,i)
            f3 = brillioun_ks_brdcst(3,i)
            f4 = brillioun_k_brdcst(1,i)
            f5 = brillioun_k_brdcst(2,i)
            f6 = brillioun_k_brdcst(3,i)
            if (mprint) write(luout,1257) f0,f1,f2,f3,f4,f5,f6
         end do
      else
        if (mprint) write(luout,1258)
      end if

      if1 = Cram_nwave_all_brdcst(0)
      if2 = Cram_nwave_brdcst(0)
      if (mprint) then
         write(luout,1249)
         write(luout,1250) lattice_ecut(),ngrid(1),ngrid(2),ngrid(3),
     >                     if1,if2
      end if

      if (brillioun_print()) then
        do i=1,brillioun_nbrillioun()
          if1 = Cram_nwave_all_brdcst(i)
          if2 = Cram_nwave_brdcst(i)
          if (mprint) then
          write(luout,1251) i,lattice_wcut(),ngrid(1),ngrid(2),ngrid(3),
     >                      if1,if2
          end if
        end do
      else
        if (mprint) write(luout,1252)
      end if


      if (mprint) then


         call c_rho_symmetrizer_print(luout)
         
         write(luout,1270)
         if (.not.control_translation()) write(luout,1271)
         if (.not.control_rotation())    write(luout,1272)
         write(luout,1280) control_time_step(),control_fake_mass()
         write(luout,1290) control_tole(),control_tolc(),control_tolr()
         write(luout,1281) control_it_in()*control_it_out(),
     >                 control_it_in(),control_it_out()

         if (control_fractional()) then
           write(6,1297)
           if (control_fractional_smeartype().eq.0)
     >       write(6,1298) "step function"
           if (control_fractional_smeartype().eq.1)
     >       write(6,1298) "Fermi-Dirac"
           if (control_fractional_smeartype().eq.2)
     >       write(6,1298) "Gaussian"
           write(6,1299) control_fractional_kT(),
     >                   control_fractional_temperature()
         end if
         write(luout,1300)
         write(luout,1305)
         call util_flush(luout)
      end if
c
c*                |***************************|
c******************     start iterations      **********************
c*                |***************************|
c
      if (taskid.eq.MASTER) call current_second(cpu2)
      if (taskid.eq.MASTER) CALL nwpw_MESSAGE(2)
      it_in  = control_it_in()
      it_out = control_it_out()
      icount = 0
   1  continue
         icount = icount + 1
         call band_inner_loop(ispin,ispinq,ne,neq,nbrillioun,nbrillq,
     >                      nfft3d,
     >                      psi1_tag,psi2_tag,dbl_mb(dn(1)),
     >                      it_in,E,deltae,deltac,deltar,
     >                      hml_tag,
     >                      psir_tag,Hpsi_tag)

         if (mprint) then 
           write(luout,1310) icount*it_in,E(1),deltae,deltac,deltar
           call util_flush(luout)
         end if
         if ((deltae.gt.0.0d0).and.(icount.gt.1)) then
            if (mprint) 
     >       write(luout,*) 
     >       ' *** Energy going up.  iteration terminated.'
            go to 2
         end if
         deltae = dabs(deltae)
         if ((deltae.lt.control_tole()).and.
     >       (deltac.lt.control_tolc()).and.
     >       (deltar.lt.control_tolr())) then
            if (mprint) 
     >       write(luout,*) 
     >       ' *** tolerance ok.     iteration terminated.'
            go to 2
         end if
      if (icount.lt.it_out) go to 1
      if (mprint) 
     > write(luout,*) 
     > '*** arrived at the Maximum iteration.   terminated.'

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
         call C3dB_r_dsum(1,dbl_mb(dn(1)+(ms-1)*nfft3d),sumall)
         en(ms) = sumall*lattice_omega()
     >             /dble(ngrid(1)*ngrid(2)*ngrid(3))
      end do
      if (taskid.eq.MASTER) then
         write(17,1320) (en(ms),ms=1,ispin)
      end if

*     **** comparison between hamiltonian an lambda matrix ****
*     **** not done - because this can generate a very large data file ****

*     **** check orthonormality ****
      if (taskid.eq.MASTER) then
         write(17,1350)
      end if
      
c      lmd_tag   = cpsi_data_alloc(1,1,2*Pneb_w_size(0,1))
c      lmd_shift = cpsi_data_get_chnk(lmd_tag,1)
c      do nb=1,nbrillq
c        psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
c        call Pneb_ffw_Multiply(nb,0,dbl_mb(psi1_shift),
c     >                              dbl_mb(psi1_shift),
c     >                              dbl_mb(lmd_shift))
c         do ms=1,ispin
c            do j=1,ne(ms)
c            do i=j,ne(ms)
c               if (taskid.eq.MASTER) 
c     >         write(17,1360) ms,i,j,
c     >                        Pneb_w_value(nb,ms,i,j,dbl_mb(lmd_shift))
c            end do
c            end do
c         end do
c      end do
c      call cpsi_data_dealloc(lmd_tag)

*     **** close check file ****
      if (taskid.eq.MASTER) then
         close(17)
      end if



*     ***** diagonalize the hamiltonian matrix but don't rotate ****
      call cpsi_data_update(hml_tag)
      call cpsi_data_update(eig_tag)
      do nb=1,nbrillq
        hml_shift = cpsi_data_get_chnk(hml_tag,nb)
        eig_shift = cpsi_data_get_chnk(eig_tag,nb)
        call Pneb_w_diag(0,nb,dbl_mb(eig_shift),dbl_mb(hml_shift))
      end do
      call cpsi_data_noupdate(hml_tag)
      call cpsi_data_noupdate(eig_tag)


*     ***** diagonalize and rotate the hamiltonian matrix ****
      call cpsi_data_update(psi2_tag)
      do nb=1,nbrillq
        psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
        psi2_shift = cpsi_data_get_chnk(psi2_tag,nb)
        hml_shift  = cpsi_data_get_chnk(hml_tag,nb)
        call Pneb_fwf_Multiply(0,nb,
     >                       one,
     >                       dbl_mb(psi1_shift),npack1,
     >                       dbl_mb(hml_shift),
     >                       zero,
     >                       dbl_mb(psi2_shift))
      end do
      call cpsi_data_noupdate(psi2_tag)



*                |***************************|
****************** report summary of results **********************
*                |***************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

      if (mprint) then
         write(luout,1300)
         write(luout,1410)
         write(luout,1420)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_amm(i)
           else if (ion_q_xyzFixIon(I)) then
           write(6,1194) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_q_xyzFixIon_label(I)
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

         write(luout,1430) E(1),E(1)/ion_nion()
         write(luout,1440) E(2),E(2)/neall
         write(luout,1450) E(3),E(3)/neall
         write(luout,1460) E(4),E(4)/neall
         if (band_HFX()) then
           write(luout,1457) E(20),E(20)/n2(ispin)
         end if

         write(luout,1470) E(5),E(5)/ion_nion()
         write(luout,1480) E(6),E(6)/neall
         write(luout,1490) E(7),E(7)/neall
         write(luout,1495) E(8),E(8)/neall
         write(luout,1496) E(9),E(9)/neall
         write(luout,1497) E(10),E(10)/neall
         if (band_HFX().and.band_HFX_relaxed())  then
           write(luout,1502) E(21),E(21)/n2(ispin)
         end if
         virial = (E(10)+E(9)+E(8)+E(7))/E(6)
         write(luout,1498) virial

        if (ion_disp_on()) then
            write(luout,1720) E(33)
        end if

      end if

      NN=ne(1)-ne(2)
      EV=27.2116d0
      if (mprint) then
        if (control_fractional()) then
          if (ispin.eq.1) then
            write(luout,1507) smearfermi(1),smearfermi(1)*EV
          else
            write(luout,1507) smearfermi(1),smearfermi(1)*EV,
     >                        smearfermi(2),smearfermi(2)*EV
          end if
        end if
      end if

*     *** generate spinorbit vector ***
      if (spin_orbit) call Pneb_f_SOSpins_tag(psi2_tag,svec_tag)

*     *** printout the eigenvalue spetra ****
      if (brillioun_print()) then
      do nb=1,brillioun_nbrillioun()
        f0 = brillioun_weight_brdcst(nb)
        f1 = brillioun_ks_brdcst(1,nb)
        f2 = brillioun_ks_brdcst(2,nb)
        f3 = brillioun_ks_brdcst(3,nb)
        f4 = brillioun_k_brdcst(1,nb)
        f5 = brillioun_k_brdcst(2,nb)
        f6 = brillioun_k_brdcst(3,nb)
        if (mprint) then
          write(luout,1508) nb,f0,f1,f2,f3,f4,f5,f6
          write(luout,1500)
        end if
        next = cpsi_data_get_next(psi1_tag)
        if (spin_orbit) then
c          if (mprint) write(luout,1511)
          do i=0,ne(1)-1
            f1=cpsi_eig_brdcst_tag(ne,spin_orbit,eig_tag, nb,1,ne(1)-i)
            f2=cpsi_sv_brdcst_tag( spin_orbit,svec_tag,nb,ne(1)-i,1)
            f3=cpsi_sv_brdcst_tag( spin_orbit,svec_tag,nb,ne(1)-i,2)
            f4=cpsi_sv_brdcst_tag( spin_orbit,svec_tag,nb,ne(1)-i,3)
            f0 = dsqrt(f2*f2 + f3*f3 + f4*f4)
            f5=cpsi_occ_brdcst_tag(ne,spin_orbit,next,    nb,1,ne(1)-i)
            if (mprint) write(luout,1512) f1,f1*EV,f0,f2,f3,f4,f5
          end do
        else
          do i=0,NN-1
            f1=cpsi_eig_brdcst_tag(ne,spin_orbit,eig_tag,nb,1,ne(1)-i)
            f2=cpsi_occ_brdcst_tag(ne,spin_orbit,next,   nb,1,ne(1)-i)
            if (mprint) write(luout,1510) f1,f1*EV,f2
          end do
          do i=0,ne(2)-1
           f1=cpsi_eig_brdcst_tag(ne,spin_orbit,eig_tag,nb,1,ne(1)-i-NN)
           f2=cpsi_occ_brdcst_tag(ne,spin_orbit,next,   nb,1,ne(1)-i-NN)
           f3=cpsi_eig_brdcst_tag(ne,spin_orbit,eig_tag,nb,2,ne(2)-i)
           f4=cpsi_occ_brdcst_tag(ne,spin_orbit,next,   nb,2,ne(2)-i)
           if (mprint) write(luout,1510) f1,f1*EV,f2,f3,f3*EV,f4
          end do
        end if
      end do
      endif



*     ***** extra energy output for QA test ****
      if (mprint) then
         write(luout,1600) E(1)
      end if

*                |***************************|
******************         Prologue          **********************
*                |***************************|
c
c*     **** calculate spin contamination ****
c      call Calculate_psi_spin2(ispin,ne,npack1,dcpl_mb(psi2(1)),
c     >                         fractional,dbl_mb(occ2(1)),w)
c
c*     **** calculate the Dipole ***
c      call Calculate_Dipole(ispin,ne,n2ft3d,dbl_mb(dn(1)),dipole)
c      
c*     **** perfom Lubin and Mulliken analysis ***
c      if (control_Mulliken()) then
c
c*       **** Lubin Water Analysis ***
c        call pspw_Lubin_water_analysis(rtdb,ispin,ne,n2ft3d,
c     >                                 dbl_mb(dn(1)))
c
c*       **** Analysis ***
c        call pspw_analysis(0,rtdb,ispin,ne,dcpl_mb(psi2(1)),
c     >                                   dbl_mb(eig(1)))
c
c*       **** generate APC *****
c        call pspw_dngen_APC(ispin,ne,dbl_mb(dn(1)))
c        call pspw_print_APC(6)
c
c      end if

*     *****  write psi2 wavefunctions ****
      call cpsi_write(spin_orbit,ispin,ne,nbrillioun,psi2_tag)

*     **** write geometry to rtdb ****
      call ion_write(rtdb)


*     **** deallocate heap memory ****
      call ewald_end()
      call cstrfac_end()
      call c_coulomb_end()
      call cke_end()
      call cpsp_end()
      call Cram_end()
      call c_G_end()
      call brillioun_end()

      call ion_end()
      call ion_end_FixIon()
      call c_rho_symmetrizer_end()
c      call pspw_end_SIC()
      call band_end_HFX()
c      call pspw_end_APC()
c      call pspw_qmmm_end()
c
      call cpsi_data_dealloc(Hpsi_tag)
      call cpsi_data_dealloc(psir_tag)
      next = cpsi_data_get_next(psi1_tag)
      if (next.ge.0) call cpsi_data_dealloc(next)
      call cpsi_data_dealloc(psi1_tag)

      next = cpsi_data_get_next(psi2_tag)
      if (next.ge.0) call cpsi_data_dealloc(next)
      call cpsi_data_dealloc(psi2_tag)
      call cpsi_data_dealloc(hml_tag)
      call cpsi_data_dealloc(eig_tag)
      call cpsi_data_dealloc(svec_tag)
      value = BA_free_heap(dn(2))
      if (.not. value)
     >  call errquit('band_sd:freeing heap memory',0,MA_ERR)


      call C3dB_pfft_end()
      call cpsi_data_end()
      call C3dB_end(1)
      call Pneb_end()
      !if ((control_gga().ge.10).and.(control_gga().le.200)) then
      if (.not.control_only_lda()) then
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
         call nwpw_timing_print_final(.true.,(icount*it_in))
         CALL nwpw_MESSAGE(4)
      end if 


      call Parallel3d_Finalize()
      call Parallel_Finalize()
      band_sd = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*     Car-Parrinello solid-state calculation       *')
 1030 FORMAT(10X,'*     [     steepest descent minimization   ]      *')
 1035 FORMAT(10x,'*     [ NorthWest Chemistry implementation ]       *')
 1040 FORMAT(10X,'*            version #1.00   03/14/09              *')
 1041 FORMAT(10X,'*    This code was developed by Eric J. Bylaska    *')
 1042 FORMAT(10X,'*                                                  *')
 1043 FORMAT(10X,'*                                                  *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ BAND input data ===================')
 1111 FORMAT(/' number of processors used:',I16)
 1112 FORMAT( ' parallel mapping         :         1d slab')
 1113 FORMAT( ' parallel mapping         :      2d hilbert')
 1114 FORMAT( ' parallel mapping         :        balanced')
 1115 FORMAT(/' options:')
 1116 FORMAT( ' parallel mapping         : not balanced')
 1117 FORMAT( ' processor grid           :',I4,' x',I4,' x',I4)
 1118 FORMAT( ' parallel mapping         :       2d hcurve')
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1132 FORMAT(5X,' using fractional occupation')
 1135 FORMAT(/' The masses of QM H atoms converted to 2.0 amu.',
     >       /' To turn off this default',
     >        ' set nwpw:makehmass2 .false.')
 1140 FORMAT(/' elements involved in the cluster:')
 1150 FORMAT(5X,I2,': ',A4,'  core charge:',F4.1,'  lmax=',I1)
 1151 FORMAT(5X,'        cutoff =',4F8.3)
 1152 FORMAT(12X,' highest angular component      : ',i3)
 1153 FORMAT(12X,' local potential used           : ',i3)
 1154 FORMAT(12X,' number of non-local projections: ',i3)
 1155 FORMAT(12X,' semicore corrections included  : ',
     >       F6.3,' (radius) ',F6.3,' (charge)')
 1156 FORMAT(12X,' aperiodic cutoff radius        : ',F6.3)
 1157 FORMAT(12X,' comment    : ',A)
 1158 FORMAT(12X,' pseudpotential type            : ',i3)

 1159 FORMAT(/' total charge=',F8.3)
 1160 FORMAT(/' atomic composition:')
 1170 FORMAT(7(5X,A2,':',I3))
 1180 FORMAT(/' initial position of ions:')
 1190 FORMAT(5X, I4, A5, ' (',3F11.5,' ) - atomic mass= ',F7.3,' ',A)
 1191 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F6.3,' - fixed ',A)
 1193 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - z fixed')
 1194 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,A)
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1211 FORMAT(5X,'   number of constraints = ', I6,' ( DOF = ',I6,' )' )
 1220 FORMAT(/' number of electrons: spin up=',F6.2, 16x,
     >                               '  down=',F6.2,A)
c 1220 FORMAT(/' number of electrons: spin up=',I6,
c     >        ' (',I4,' per task)',
c     >        '  down=',I6,
c     >        ' (',I4,' per task)',
c     >        A)
 1221 FORMAT( ' number of orbitals : spin up=',I6,  
     >        ' (',I4,' per task)',
     >        '  down=',I6,
     >        ' (',I4,' per task)',
     >        A)
 1230 FORMAT(/' supercell:')
 1231 FORMAT(5x,' volume : ',F10.1)
 1232 FORMAT(/5x,' lattice:    a=',f8.3,'    b=',f8.3,'     c=',f8.3,
     >       /5x,'         alpha=',f8.3,' beta=',f8.3,' gamma=',f8.3)
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
 1252 FORMAT(5x, ' wavefunction grids not printed - ',
     >           'number of k-point is very large')

 1255 FORMAT(/' brillouin zone:')
 1256 FORMAT(5x,' number of zone points:',I3)
 1257 FORMAT(5x,' weight=',f8.3,'  ks=<',3f8.3,' >, k=<',3f8.3,'>')
 1258 FORMAT(5x,' number of k-point is very large or distributed')

 1260 FORMAT(5X,' Ewald summation: cut radius=',F8.2,'  and',I3)
 1261 FORMAT(5X,'                   madelung=',f14.8)
 1270 FORMAT(/' technical parameters:')
 1271 FORMAT(5x, ' translation constrained')
 1272 FORMAT(5x, ' rotation constrained')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'fictitious mass=',F10.1)
 1281 FORMAT(5X, ' maximum iterations =',I10,
     >           ' ( ',I4,' inner ',I6,' outer )')
 1290 FORMAT(5X, ' tolerance=',E8.3,' (energy)',E12.3,
     &        ' (electron)',E12.3,' (ion)')
 1297 FORMAT(/' fractional smearing parameters:')
 1298 FORMAT(5X, ' smearing algorithm   = ',A)
 1299 FORMAT(5X, ' smearing parameter   = ',E9.3,' (',F7.1,' K)')
 1300 FORMAT(//)
 1305 FORMAT(10X,'================ iteration =========================')
 1310 FORMAT(I8,E20.10,3E15.5)
 1320 FORMAT(' number of electrons: spin up=',F11.5,'  down=',F11.5,A)
 1330 FORMAT(/' comparison between hamiltonian and lambda matrix')
 1331 FORMAT(/' Elements of Hamiltonian matrix (up/restricted)')
 1332 FORMAT(/' Elements of Hamiltonian matrix (down)')
 1340 FORMAT(I3,2I3,' H=',E16.7,', L=',E16.7,', H-L=',E16.7)
 1341 FORMAT(I3,2I3,' H=',E16.6)
 1350 FORMAT(/' orthonormality')
 1360 FORMAT(I3,2I3,'(',2E18.7,')')
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
 1502 FORMAT( ' K.S. HFX energy     :',E19.10,
     >        ' (',E15.5,'/electron)')
 1500 FORMAT(/' orbital energies:')
 1507 FORMAT(/' Fermi energy =',2(E18.7,' (',F8.3,'eV)'))
 1508 FORMAT(/' Brillouin zone point: ',i3,
     >       /'    weight=',f10.6,
     >       /'    k     =<',3f8.3,'> . <b1,b2,b3> ',
     >       /'          =<',3f8.3,'>')
 1510 FORMAT(4(E18.7,' (',F8.3,'eV) occ=',F5.3))
 1511 FORMAT(33x,"Spin(Sz,Sy,Sz)")
 1512 FORMAT(E18.7,' (',F8.3,' eV) (|s| =',F6.3, 
     >       ', s = <',F7.3,',',F7.3,',',F7.3,'> ) occ=',F5.3)

 1600 FORMAT(/' Total BAND energy   :',E19.10)

 1700 FORMAT(/' QM/MM-pol-vib/CAV Energies')
 1701 FORMAT( ' --------------------------')
 1702 FORMAT( ' LJ energy              :',E19.10)
 1703 FORMAT( ' Residual Coulomb energy:',E19.10)
 1704 FORMAT( ' MM Vibration energy    :',E19.10)
 1705 FORMAT( ' MM Vibration energy    :',E19.10)
 1706 FORMAT( ' (QM+MM)/Cavity energy  :',E19.10)

 1720 FORMAT(/' Dispersion energy   :',E19.10)

 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(6,9010) ierr
      call Parallel_Finalize()

      band_sd = value
      return
      END
