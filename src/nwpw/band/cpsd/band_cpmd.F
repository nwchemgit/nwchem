*
* $Id$
*

***********************************************************************
*                                                                     *
*                         band_cpmd                                   *
*                                                                     *
***********************************************************************

      logical function band_cpmd(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "inp.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
      
      logical value

      real*8 kb
      parameter (kb=3.16679d-6)
c      real*8 autoatm
c      parameter (autoatm =290.360032539d6)


      
*     **** parallel variables ****
      integer  taskid,np,np_i,np_j,np_k
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,nfft3d,n2ft3d
      integer npack0,npack1

*     **** electronic variables ****
      logical spin_orbit
      real*8 icharge
      integer ispin,ispinq
      integer ne(2),nemax,neall,neq(2),nemaxq
      integer nbrillioun,nbrillq
      real*8  en(2)
      real*8 dipole(3)

      integer psi0_tag,psi1_tag,psi2_tag,next
      integer psi0_shift,psi1_shift,psi2_shift
      integer dn(2)
      integer Hpsi_tag,psir_tag
    

*     ***** energy variables ****
      real*8  E(50),eke,eave,evar,cv

      integer eig_tag,hml_tag,svec_tag
      integer eig_shift,hml_shift,svec_shift

*     **** psi smearing block ****
      logical fractional
      integer smearoccupation,smeartype
      real*8 smearfermi(2),smearcorrection,smearkT

*     **** error variables ****
      integer ierr

*     **** local variables ****
      complex*16 zero,one
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))

      logical verlet,mulliken,SA,found,calc_pressure
      logical lprint,mprint,hprint
      integer ms,if1,if2,vers
      real*8  gx,gy,gz,cx,cy,cz
      real*8  vgx,vgy,vgz,vcx,vcy,vcz
      real*8  ekg,eki0,eki1,sum
      real*8  eke0,eke1,f0,f1,f2,f3,f4,f5,f6
      real*8  EV,pi,dt,a,b,c,alpha,beta,gamma
      real*8  emotion_time_shift
      integer i,j,k,ia,n,nn,nb
      integer ii,jj,index,indx
      integer icount,it_in,it_out,icount_shift
      real*8 w,sumall,pressure,stress(3,3),p1,p2
      real*8 Te_init,Tr_init,Te_new,Tr_new,sa_decay(2),sa_alpha(2)
      integer nfft3,mapping,mapping1d
      parameter (nfft3=32)
      character*50 filename
      character*255 full_filename
  


*     **** external functions ****
      real*8      cpsp_zv,cpsp_rc,ewald_rcut,ion_amass
      real*8      ewald_mandelung
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg
      integer     ewald_ncut,ewald_nshl3d
      integer     cpsp_nprj,cpsp_lmax,cpsp_locp,ion_nkatm,cpsp_psp_type
      character*4 ion_atom,ion_aname
      external    cpsp_zv,cpsp_rc,ewald_rcut,ion_amass
      external    ewald_mandelung
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg
      external    ewald_ncut,ewald_nshl3d
      external    cpsp_nprj,cpsp_lmax,cpsp_locp,ion_nkatm,cpsp_psp_type
      external    ion_atom,ion_aname

      real*8   control_rti,control_rte,ion_rion
      real*8   ion_vion,ion_com_ke,ion_ke,control_fractional_kT
      real*8   ion_Temperature,ion_com_Temperature
      external control_rti,control_rte,ion_rion
      external ion_vion,ion_com_ke,ion_ke,control_fractional_kT
      external ion_Temperature,ion_com_Temperature
      real*8   control_time_step,control_fake_mass,brillioun_weight
      external control_time_step,control_fake_mass,brillioun_weight
      logical  control_read,control_move,ion_init,ion_q_FixIon
      external control_read,control_move,ion_init,ion_q_FixIon
      logical  ion_q_xyzFixIon
      external ion_q_xyzFixIon
      character*14 ion_q_xyzFixIon_label
      external     ion_q_xyzFixIon_label

      logical  brillioun_print
      external brillioun_print
      real*8   brillioun_weight_brdcst
      external brillioun_weight_brdcst
      real*8   brillioun_ks_brdcst,brillioun_k_brdcst
      external brillioun_ks_brdcst,brillioun_k_brdcst
      integer  Cram_nwave_all_brdcst,Cram_nwave_brdcst
      external Cram_nwave_all_brdcst,Cram_nwave_brdcst
      real*8   cpsi_eig_brdcst_tag,cpsi_occ_brdcst_tag
      external cpsi_eig_brdcst_tag,cpsi_occ_brdcst_tag
      real*8   cpsi_sv_brdcst_tag
      external cpsi_sv_brdcst_tag

      integer  brillioun_nbrillioun,control_fractional_smeartype
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      external brillioun_nbrillioun,control_fractional_smeartype
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      integer  cpsi_data_alloc,cpsi_data_get_next,cpsi_data_get_chnk
      external cpsi_data_alloc,cpsi_data_get_next,cpsi_data_get_chnk
      integer  Pneb_nbrillq,Pneb_ispinq,Pneb_w_size,Cram_nwave_all
      external Pneb_nbrillq,Pneb_ispinq,Pneb_w_size,Cram_nwave_all

      character*12 control_boundry
      external     control_boundry

      logical  pspw_reformat_c_wvfnc,pspw_reformat_c_v_wvfnc
      external pspw_reformat_c_wvfnc,pspw_reformat_c_v_wvfnc
      logical      cpsp_semicore,pspw_qmmm_found
      real*8       cpsp_rcore,cpsp_ncore,cpsp_rlocal,ion_TotalCharge
      external     cpsp_semicore,pspw_qmmm_found
      external     cpsp_rcore,cpsp_ncore,cpsp_rlocal,ion_TotalCharge

      logical  control_Nose,control_Mulliken,control_print,Nose_restart
      external control_Nose,control_Mulliken,control_print,Nose_restart

      integer  Nose_Mchain,Nose_Nchain,psi_get_version,v_psi_get_version
      external Nose_Mchain,Nose_Nchain,psi_get_version,v_psi_get_version
  
      real*8   control_Nose_Te,Nose_Qe,Nose_Pe,Nose_Ee0
      external control_Nose_Te,Nose_Qe,Nose_Pe,Nose_Ee0

      real*8   control_Nose_Tr,Nose_Qr,Nose_Pr,Nose_Er0
      external control_Nose_Tr,Nose_Qr,Nose_Pr,Nose_Er0
      logical      v_psi_filefind
      external     v_psi_filefind
      real*8   nwpw_timing
      external nwpw_timing

      logical  control_out_of_time,control_new_vpsi,control_fractional
      external control_out_of_time,control_new_vpsi

      logical  control_SA,control_Fei,pspw_SIC,pspw_HFX,control_pressure
      real*8   control_SA_decay
      external control_SA,control_Fei,pspw_SIC,pspw_HFX,control_pressure
      external control_SA_decay

      integer  control_np_dimensions,control_mapping,control_mapping1d
      external control_np_dimensions,control_mapping,control_mapping1d


      logical  control_translation,control_rotation,control_balance
      external control_translation,control_rotation,control_balance
     
      character*255 cpsp_comment,comment
      external      cpsp_comment
      integer  ion_nconstraints,ion_ndof
      external ion_nconstraints,ion_ndof
      logical  ion_disp_on,ion_makehmass2,control_only_lda
      external ion_disp_on,ion_makehmass2,control_only_lda




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

      if (.not.control_read(14,rtdb))
     >   call errquit('band_cpmd:error reading control',0,DISK_ERR)

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
         write(luout,1031)
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

*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)

*     **** initialize lattice and packing data structure ****
      call lattice_init()
      call c_G_init()
      call brillioun_init()
      call Cram_Init()
      call C3dB_pfft_init()

*     ***** Initialize double D3dB data structure ****
      !if ((control_gga().ge.10).and.(control_gga().lt.100)) then
      if (.not.control_only_lda()) then
         call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
         call G_init()
         call mask_init()
      end if

c*     **** read ions ****
c      value = ion_init(rtdb)
c      call center_geom(cx,cy,cz)
c      call center_mass(gx,gy,gz)

c*     **** initialize FixIon constraint ****
c      call ion_init_FixIon(rtdb)

*     **** allocate psp data structure and read in psedupotentials into it ****
      call cpsp_init()
      call cpsp_readall()
      if (cpsp_semicore(0)) call c_semicore_check()


*     **** initialize G,mask,ke,and coulomb data structures ****
      call cstrfac_init()
      call cke_init()
      call c_coulomb_init()
      call ewald_init()


      call psi_get_ne(ispin,ne)
      if (ispin.eq.3) then
         spin_orbit = .true.
         ispin=2
      else
         spin_orbit = .false.
      end if
      nbrillioun = brillioun_nbrillioun()
      mapping1d  = control_mapping1d()
      call Pneb_init(ispin,ne,nbrillioun,mapping1d,spin_orbit)
      call Pneb_neq(neq)


*     ***** allocate psi2, psi1, and psi0 wavefunctions ****
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
      psi0_tag = cpsi_data_alloc(nbrillq,nemaxq,2*npack1)

*     *** fractional orbitals ***
      if (smearoccupation.gt.0) then
        call cpsi_data_set_next(psi2_tag,
     >                 cpsi_data_alloc(nbrillq,nemaxq,1))
        call cpsi_data_set_next(psi1_tag,
     >                 cpsi_data_alloc(nbrillq,nemaxq,1))
        call cpsi_data_set_next(psi0_tag,
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
     >  call errquit('band_cpmd:out of heap memory',0,MA_ERR)

*     **** convert from pspw format to band format ****
      vers = psi_get_version()
      if ((vers.eq.3).or.(vers.eq.4)) then
        value = btdb_parallel(.false.)
        if (taskid.eq.MASTER) then
          value= pspw_reformat_c_wvfnc(1)
        end if
        value = btdb_parallel(.true.)
      end if

*     *****  read psi2 wavefunctions ****
      call cpsi_read(spin_orbit,ispin,ne,nbrillioun,psi2_tag)

*     **** move  wavefunction velocities ****
      if (control_new_vpsi()) then
        call v_cpsi_delete()
      end if

*     **** generate initial wavefunction velocities if it does not exist ****
      if (.not.v_psi_filefind()) then
        call v_cpsi_new(spin_orbit,ispin,ne,nbrillioun)
      end if

*     **** convert from pspw format to band format ****
      vers = v_psi_get_version()
      if ((vers.eq.3).or.(vers.eq.4)) then
        value = btdb_parallel(.false.)
        if (taskid.eq.MASTER) then
          value= pspw_reformat_c_v_wvfnc(1)
        end if
        value = btdb_parallel(.true.)
      end if


*     *****  read psi0 wavefunctions ****
      call v_cpsi_read(spin_orbit,ispin,ne,nbrillioun,psi1_tag)


c*     **** initialize QM/MM ****
c      call pspw_qmmm_init(rtdb)
c


*     ******************************
*     **** scaling psi velocity ****
*     ******************************
      call cpsi_data_copyall(psi1_tag,psi0_tag)
      call cpsi_data_scalall(control_rte(),psi1_tag)
      call Pneb_ff_traceall_tag(0,0,psi0_tag,psi0_tag,eke0)
      call Pneb_ff_traceall_tag(0,0,psi1_tag,psi1_tag,eke1)
      eke0 = control_fake_mass()*eke0
      eke1 = control_fake_mass()*eke1
      call ion_init_ke(ekg,eki0,eki1)

*     **** Initialize thermostats ****
      if (control_Nose()) then
         call cke_ave(ispin,neq,psi2_tag,w)
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


c*     **** initialize dplot ****
c      call dplot_iteration_init()


c*     **** initialize SIC and HFX ****
c      call pspw_init_SIC(rtdb,ne)
c      call pspw_init_HFX(rtdb,ispin,ne)


*     **** initialize pressure ****
      calc_pressure = control_pressure()
      pressure      = 0.0d0
      p1            = 0.0d0
      p2            = 0.0d0
      if (calc_pressure) then
         call cpsp_stress_init()
         call cpsp_stress_readall()
      end if



*                |**************************|
******************   summary of input data  **********************
*                |**************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      call center_v_geom(vcx,vcy,vcz)
      call center_v_mass(vgx,vgy,vgz)
      mulliken = control_Mulliken()

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
               indx = indx + (ms-1)*ne(1)
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
         if (control_balance()) then
           write(luout,1114)
         else
           write(luout,1116)
         end if

         write(luout,1115)
         write(luout,1121) control_boundry(),control_version()
         if (ispin.eq.1) write(luout,1130) 'restricted'
         if (ispin.eq.2) write(luout,1130) 'unrestricted'

         call v_bwexc_print(luout,control_gga())

c         if (fractional) write(6,1132)
c         call pspw_print_SIC(6)
c         call pspw_print_HFX(6)

         if (ion_makehmass2()) write(luout,1135)
         write(luout,1140)
         do ia = 1,ion_nkatm()
           write(luout,1150) ia,ion_atom(ia),cpsp_zv(ia),cpsp_lmax(ia)
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
         write(luout,1159) icharge

         write(luout,1160)
         write(luout,1170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
         write(luout,1180)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),
     >                    (ion_rion(K,I),K=1,3),ion_amass(i)/1822.89d0
           else if (ion_q_xyzFixIon(I)) then
           write(luout,1194) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_q_xyzFixIon_label(I)
           else
           write(luout,1190) I,ion_aname(I),
     >                    (ion_rion(K,I),K=1,3),ion_amass(i)/1822.89d0
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz

         write(luout,1181)
         write(luout,1192) (I,ion_aname(I),
     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
         write(luout,1200) vcx,vcy,vcz
         write(luout,1210) vgx,vgy,vgz
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
         write(luout,1270)
         if (.not.control_translation()) write(luout,1271)
         if (.not.control_rotation())    write(luout,1272)
         write(luout,1280) control_time_step(),control_fake_mass()
         write(luout,1290) control_rte(),control_rti()
         call ion_scaling_atoms_print(luout)
         write(luout,1281) control_it_in()*control_it_out(),
     >                 control_it_in(),control_it_out()

         write(luout,1222) eke0,eki0,ekg
         write(luout,1223) eke1,eki1
         write(luout,1224) (eke1-eke0),(eki1-eki0)
         if (control_Nose()) then
           write(luout,1295)
           if (Nose_restart()) then
              write(luout,*) "    thermostats resused"
           else
              write(luout,*) "    thermostats initialized"
           end if
           do i=1,Nose_Mchain()
             write(luout,1297) i,control_Nose_Te(),Nose_Qe(i),
     >                     Nose_Pe(i),Nose_Ee0(i)
           end do
           do i=1,Nose_Nchain()
             write(luout,1298) i,control_Nose_Tr(),Nose_Qr(i),
     >                   Nose_Pr(i),Nose_Er0(i)
           end do
         else
           write(luout,1294)
         end if
        if (calc_pressure) write(luout,1293)
        if (control_SA()) then
           write(luout,1296) sa_decay(1),sa_decay(2)
         end if


         if (mulliken) write(luout,1299)
         write(luout,1300)
         write(luout,1305)
         call util_flush(luout)
      end if

*                |***************************|
******************     start iterations      **********************
*                |***************************|
*     **** open xyz and MOTION file ****
      call xyz_init()          ! unit=18
      call MOTION_init(rtdb)   ! unit=19

*     *** fei io ****
      call fei_init(rtdb)

*     **** ecce print ****
      call ecce_print_module_entry('task Car-Parrinello')
      !call ecce_print_module_entry('driver')
      call movecs_ecce_print_off()




*     ************************************
*     **** open up other MOTION files ****
*     ************************************


*     **** open EMOTION file ****
      if (.not.btdb_cget(rtdb,'cpmd:emotion_filename',1,filename))
     >  call util_file_prefix('emotion',filename)
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
c      if (mulliken) then
c        if (.not.btdb_cget(rtdb,'cpmd:eigmotion_filename',1,filename))
c     >    call util_file_prefix('eigmotion',filename)
c      call util_file_name_noprefix(filename,.false.,
c     >                             .false.,
c     >                    full_filename)
c      if (taskid.eq.MASTER) 
c     >   open(unit=32,file=full_filename,form='formatted')
c      end if

c*     **** open HMOTION file ****
c      if (mulliken) then
c       if (.not.btdb_cget(rtdb,'cpmd:hmotion_filename',1,filename))
c     >  call util_file_prefix('hmotion',filename)
c      call util_file_name_noprefix(filename,.false.,
c     >                             .false.,
c     >                    full_filename)
c      if (taskid.eq.MASTER) 
c     >   open(unit=34,file=full_filename,form='formatted')
c      end if

c*     **** open OMOTION file ****
c      if (mulliken) call Orb_Init(rtdb,ispin,ne) !unit=33

*     **** write initial position to xyz data ****
      call xyz_write()

*     ***** first step using velocity ****
      verlet = .false.
      call band_inner_loop_md(verlet,sa_alpha,
     >                      ispin,ispinq,ne,neq,nbrillioun,nbrillq,
     >                      nfft3d,
     >                      psi0_tag,psi1_tag,psi2_tag,dbl_mb(dn(1)),
     >                      1,0,
     >                      E,
     >                      hml_tag,
     >                      psir_tag,Hpsi_tag,
     >                      calc_pressure,pressure,p1,p2)

      if (taskid.eq.MASTER) call current_second(cpu2)
      if ((mprint).and.(.not.calc_pressure))  call nwpw_message(6)
      if ((mprint).and.     (calc_pressure))  call nwpw_message(9)
      
       it_in  = control_it_in()
       it_out = control_it_out()
       icount = 0
       verlet = .true.
       eke    = 0.0d0
       if (it_out.lt.1) goto 102
 
 
       dt = control_time_step()
 
       Te_new = Te_init
       Tr_new = Tr_init
  101 continue
         icount = icount + 1
         call  band_inner_loop_md(verlet,sa_alpha,
     >                      ispin,ispinq,ne,neq,nbrillioun,nbrillq,
     >                      nfft3d,
     >                      psi0_tag,psi1_tag,psi2_tag,dbl_mb(dn(1)),
     >                      it_in,((icount-1)*it_in),
     >                      E,
     >                      hml_tag,
     >                      psir_tag,Hpsi_tag,
     >                      calc_pressure,pressure,p1,p2)
         eke = eke + E(3)

         !**** calculate pressure ****

         if (taskid.eq.MASTER) then 

           if (mprint) then
              if (calc_pressure) then
                if (SA) then
                write(luout,1309) icount*it_in,E(1),E(2),E(3),E(4),
     >                        Te_new,Tr_new,pressure
                else
                write(luout,1310) icount*it_in,E(1),E(2),E(3),E(4),
     >                        ion_Temperature(),pressure
                end if
              else
                if (SA) then
                write(luout,1309) icount*it_in,E(1),E(2),E(3),E(4),
     >                        Te_new,Tr_new
                else
                write(luout,1310) icount*it_in,E(1),E(2),E(3),E(4),
     >                        ion_Temperature()
                end if
              end if
              call util_flush(luout)
           end if

*          **** write out EMOTION data ****
           eave = E(25)/dble(icount+icount_shift)
           evar = E(26)/dble(icount+icount_shift)
           evar = evar - eave*eave
           if (control_Nose()) then
             write(31,1311) icount*it_in*dt + emotion_time_shift,
     >                    e(1),e(2),e(3),e(4),e(5),e(6),
     >                    e(7),e(8),e(9),e(10),
     >                    eave,evar,ion_Temperature(),pressure
           else
             write(31,1311) icount*it_in*dt + emotion_time_shift,
     >                    e(1),e(2),e(3),e(4),e(5),e(6),
     >                    e(7),e(8),
     >                    eave,evar,ion_Temperature(),pressure
           end if
           call util_flush(31)

c*          **** write out EIGMOTION data -diagonal hml matrix ****
c           if (mulliken) then
c           write(32,1311) icount*it_in*dt,
c     >       (( dbl_mb(hml(1)+ii-1+(ii-1)*ne(1)+(ms-1)*ne(1)*ne(1)), 
c     >         ii=1,ne(ms)),ms=1,ispin)
c           call util_flush(32)
c           end if

c*          **** write out HMOTION data - hml matrix ****
c           if (mulliken) then
c           write(34,1312) icount*it_in*dt,ispin
c           do ms=1,ispin
c             write(34,1313) ms,ne(ms),ne(ms)
c             do ii=1,ne(ms)
c               write(34,1311) 
c     >         (dbl_mb(hml(1)+ii-1+(jj-1)*ne(1)+(ms-1)*ne(1)*ne(1)), 
c     >          jj=1,ne(ms))
c             end do
c           end do
c           call util_flush(34)
c           end if

         end if

*        **** write xyz data ****
         call xyz_write()
         call MOTION_write((icount*it_in*control_time_step()))

c*        **** write OMOTION data ****
c         if (mulliken) call Orb_Write(dcpl_mb(psi1(1)))

*        **** update thermostats using SA decay ****
         if (SA) then
           t1 = icount*it_in*dt/sa_decay(1)
           t2 = icount*it_in*dt/sa_decay(2)
           Te_new = Te_init*dexp(-t1)
           Tr_new = Tr_init*dexp(-t2)
           call Nose_reset_T(Te_new,Tr_new)
         end if

*        **** exit early ****
         if (control_out_of_time()) then
            if (mprint)
     >       write(luout,*) 
     >       ' *** out of time. iteration terminated'
            go to 102
         end if
      if (icount.lt.it_out) go to 101
      if (mprint)
     > write(luout,*) 
     > '*** arrived at the Maximum iteration.   terminated.'
 
*::::::::::::::::::::  end of iteration loop  :::::::::::::::::::::::::

  102 continue

*     **** close xyz and MOTION files ****
      call xyz_end()
      call MOTION_end()
      if (taskid.eq.MASTER) then
        close(unit=31)
c        close(unit=32)
c        close(unit=34)
      end if

*     *** close fei io ****
      call fei_end()

c*     **** close OMOTION file ****
c      if (mulliken) call Orb_End()

*     **** ecce print ****
      !call ecce_print_module_exit('driver', 'ok')
      call ecce_print_module_exit('task Car-Parrinello', 'ok')


*     **** finalize pressure ****
      if (calc_pressure) then
         call cpsp_stress_end()
      end if


      if (mprint) CALL nwpw_message(3)
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
         call C3dB_r_dsum(1,dbl_mb(dn(1)+(ms-1)*nfft3d),sumall)
         en(ms) = sumall*lattice_omega()
     >             /dble(ngrid(1)*ngrid(2)*ngrid(3))
      end do
      if (taskid.eq.MASTER) then
         write(17,1320) (en(ms),ms=1,ispin)
      end if

c*     **** comparison between hamiltonian an lambda matrix ****
c      if (taskid.eq.MASTER) write(17,1330)
c      do ms=1,ispin
c         do i=1,ne(ms)
c         do j=1,ne(ms)
c            w   = Dneall_m_value(0,ms,i,j,dbl_mb(hml(1)))
c            sum = Dneall_m_value(0,ms,i,j,dbl_mb(lmd(1)))
c
c            if (taskid.eq.MASTER)
c     >      write(17,1340) ms,i,j,w,sum,w-sum
c
c         end do
c         end do
c      end do
c
c*     **** check orthonormality ****
c      if (taskid.eq.MASTER) then
c         write(17,1350)
c      end if
c
c      call Dneall_ffm_Multiply(0,dcpl_mb(psi1(1)),
c     >                           dcpl_mb(psi1(1)),npack1,
c     >                           dbl_mb(lmd(1)))
c      do ms=1,ispin
c         do j=1,ne(ms)
c         do i=j,ne(ms)
c            w  = Dneall_m_value(0,ms,i,j,dbl_mb(lmd(1)))
c            if (taskid.eq.MASTER) write(17,1360) ms,i,j,w
c         end do
c         end do
c      end do
c
*     **** close check file ****
      if (taskid.eq.MASTER) then
         close(17)
      end if


*     ***** diagonalize the hamiltonian matrix ****
      call cpsi_data_update(hml_tag)
      call cpsi_data_update(eig_tag)
      do nb=1,nbrillq
        hml_shift = cpsi_data_get_chnk(hml_tag,nb)
        eig_shift = cpsi_data_get_chnk(eig_tag,nb)
        call Pneb_w_diag(0,nb,dbl_mb(eig_shift),dbl_mb(hml_shift))
      end do
      call cpsi_data_noupdate(hml_tag)
      call cpsi_data_noupdate(eig_tag)


*     ***** rotate the psi and v_psi ****
      call cpsi_data_update(psi2_tag)
      call cpsi_data_update(psi0_tag)
      do nb=1,nbrillq
        psi0_shift = cpsi_data_get_chnk(psi0_tag,nb)
        psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
        psi2_shift = cpsi_data_get_chnk(psi2_tag,nb)
        hml_shift  = cpsi_data_get_chnk(hml_tag,nb)

*       **** rotate psi ****
        call Pneb_fwf_Multiply(0,nb,
     >                       one,
     >                       dbl_mb(psi1_shift),npack1,
     >                       dbl_mb(hml_shift),
     >                       zero,
     >                       dbl_mb(psi2_shift))

*       **** rotate v_psi ****
        call dcopy(2*npack1*nemaxq,dbl_mb(psi0_shift),1,
     >                             dbl_mb(psi1_shift),1)
        call Pneb_fwf_Multiply(0,nb,
     >                       one,
     >                       dbl_mb(psi1_shift),npack1,
     >                       dbl_mb(hml_shift),
     >                       zero,
     >                       dbl_mb(psi0_shift))

      end do
      call cpsi_data_noupdate(psi2_tag)
      call cpsi_data_noupdate(psi0_tag)




*                |***************************|
****************** report summary of results **********************
*                |***************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      call center_v_geom(vcx,vcy,vcz)
      call center_v_mass(vgx,vgy,vgz)

      if (mprint) then
         call print_elapsed_time(icount*it_in*dt)
         write(luout,1300)
         write(luout,1410)
         write(luout,1420)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           else if (ion_q_xyzFixIon(I)) then
           write(6,1194) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_q_xyzFixIon_label(I)
           else
           write(luout,1190) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz


         write(luout,1421)
         write(luout,1192) (I,ion_aname(I),
     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
         write(luout,1200) vcx,vcy,vcz
         write(luout,1210) vgx,vgy,vgz
         write(luout,1211) ion_nconstraints(),ion_ndof()

         write(luout,*)
         write(luout,1320) en(1),en(ispin),' (real space)'
         write(luout,1430) E(2),E(2)/ion_nion()
         write(luout,1440) E(5),E(5)/neall
         write(luout,1450) E(6),E(6)/neall
         write(luout,1460) E(7),E(7)/neall
         write(luout,1470) E(8),E(8)/neall
         write(luout,1471) E(3),E(3)/neall
         write(luout,1472) ion_ke(),ion_ke()/ion_nion()

c         if (pspw_qmmm_found()) then
c            write(6,1700)
c            write(6,1701)
c            write(6,1702) E(11)
c            write(6,1703) E(12)
c            write(6,1704) E(13)
ccc            write(6,1705) E(14)
c            write(6,1706) E(15)
c         end if
         if (ion_disp_on()) then
           write(luout,1720) E(33)
         end if

         if (control_Nose()) then
           write(luout,1473) E(9),E(9)/neall
           write(luout,1474) E(10),E(10)/ion_nion()
         end if
         write(luout,1226) E(3),ion_ke(),ion_com_ke()
         eke = eke/dble(it_out)

*        **** !!!!Cram_nwave_all probably not correct!!!!
         eke = 2.0d0*eke/kb/(ne(1)+ne(ispin))/Cram_nwave_all(1)
         !eke = 2.0d0*eke/kb/(ne(1)+ne(ispin))

*       **** write out Temperatures ****
         write(luout,1491) eke
         eki0 = ion_Temperature()
         write(luout,1480) eki0
         write(luout,1490) ion_com_Temperature()

         eave = E(25)/dble(icount+icount_shift)
         evar = E(26)/dble(icount+icount_shift)
         evar = evar - eave*eave
         cv = (evar)/(kb*ion_Temperature()**2)
         cv = cv/dble(ion_nion())
         write(luout,1492) eave
         write(luout,1493) evar
         write(luout,1494) cv

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

*     *** Extra energy output added for QA test ****
      if (mprint) then
         write(luout,1600) E(2)
      end if

*                |***************************|
******************         Prologue          **********************
*                |***************************|
c
c*     **** calculate spin contamination ****
c      call Calculate_psi_spin2(ispin,ne,npack1,dcpl_mb(psi2(1)),
c     >                         fractional,dbl_mb(occ2(1)),
c     >                         dipole)
c
c*     **** calculate the Dipole ***
c      call Calculate_Dipole(ispin,ne,n2ft3d,dbl_mb(dn(1)),dipole)

*     ***** write wavefunctions and v_wavefunctions ****
      call cpsi_write(  spin_orbit,ispin,ne,nbrillioun,psi2_tag)
      call v_cpsi_write(spin_orbit,ispin,ne,nbrillioun,psi0_tag)

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
       if (control_Nose()) call Nose_end()

c      call pspw_end_SIC()
c      call pspw_end_HFX()
c      call pspw_qmmm_end()

      call cpsi_data_dealloc(Hpsi_tag)
      call cpsi_data_dealloc(psir_tag)
      next = cpsi_data_get_next(psi0_tag)
      if (next.ge.0) call cpsi_data_dealloc(next)
      call cpsi_data_dealloc(psi0_tag)
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
     >  call errquit('band_cpmd:error freeing heap',0,MA_ERR)



      call C3dB_pfft_end()
      call cpsi_data_end()
      call C3dB_end(1)
      call Pneb_end()
      !if ((control_gga().ge.10).and.(control_gga().lt.100)) then
      if (.not.control_only_lda()) then
         call mask_end()
         call G_end()
         call D3dB_end(1)
      end if

*     **** do anaylysis on MOTION files ****
      call cpmd_properties(rtdb)


*                |***************************|
****************** report consumed cputime   **********************
*                |***************************|
      if (mprint) then
         CALL current_second(cpu4)

         t1=cpu2-cpu1
         t2=cpu3-cpu2
         t3=cpu4-cpu3
         t4=cpu4-cpu1
         av=t2/dble(icount*it_in)
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
         CALL nwpw_message(4)
      end if 


      call Parallel3d_Finalize()
      call Parallel_Finalize()
      band_cpmd = value
      return
cc
ccc
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
 1111 FORMAT(/' number of processors used:',I16)
 1112 FORMAT( ' parallel mapping         :         1d-slab')
 1113 FORMAT( ' parallel mapping         :      2d-hilbert')
 1114 FORMAT( ' parallel mapping         :        balanced')
 1115 FORMAT(/' options:')
 1116 FORMAT( ' parallel mapping         :    not balanced')
 1118 FORMAT( ' parallel mapping         :       2d-hcurve')
 1117 FORMAT( ' processor grid           :',I4,' x',I4,' x',I4)
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
 1170 FORMAT(7(5X,A4,':',I3))
 1180 FORMAT(/' initial position of ions:')
 1181 FORMAT(/' initial velocity of ions:')
 1190 FORMAT(5X, I4, A5  ,' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1191 FORMAT(5X, I4, A5  ,' (',3F11.5,
     >       ' ) - atomic mass= ',F6.3,' - fixed')
 1192 FORMAT(5X, I4, A5  ,' (',3F11.5,' )')
 1193 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - z fixed')
 1194 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,A)
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1211 FORMAT(5X,'   number of constraints = ', I6,' ( DOF = ',I6,' )' )
 1220 FORMAT(/' number of electrons: spin up=',F6.2,16x,
     >                               '  down=',F6.2,A)
 1221 FORMAT( ' number of orbitals : spin up=',I6,
     >        ' (',I4,' per task)',
     >        '  down=',I6,
     >        ' (',I4,' per task)',
     >        A)
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
 1290 FORMAT(5X, ' cooling/heatting rates: ',e12.5,' (psi)',2x,
     >                                       e12.5,' (ion)')
 1293 format(/' Pressure Output Generated         ')
 1294 format(/' Constant Energy Simulation                     ')
 1295 format(/' Nose-Hoover Simulation - Thermostat Parameters:')
 1296 format(5x, 'SA decay rates  =',e10.3,' (elc)',e10.3,' (ion)')
 1297 format(5x, 'link = ',I3,
     > ' Te =',f8.2,' Qe =',e10.3,' 2*pi/we=',e10.3,' Ee0=',e10.3)
 1298 format(5x, 'link = ',I3,
     > ' Tr =',f8.2,' Qr =',e10.3,' 2*pi/wr=',e10.3,' Er0=',e10.3)
 1299 format(//' Mulliken Analysis Output Generated            ')
 1300 FORMAT(//)
 1305 FORMAT(10X,'============ Car-Parrinello iteration ==============')
 1309 FORMAT(I8,2E19.10,2E14.5,2F9.1,3E11.3)
 1310 FORMAT(I8,2E19.10,2E14.5,F14.2,3E11.3)
 1311 format(100e19.10)
 1312 format(e14.6,i3)
 1313 format(3i4)
 1320 FORMAT(' number of electrons: spin up=',F11.5,'  down=',F11.5,A)
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
      call Parallel3d_Finalize()
      call Parallel_Finalize()

      band_cpmd = value
      return
      END
