*
* $Id: pspw_et.F 22726 2012-08-16 02:15:40Z bylaska $
*

***********************************************************************
*                      pspw_et  (MPI code)                            *
*                                                                     *
***********************************************************************

      logical function pspw_et(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
      
*     **** parallel variables ****
      integer  taskid,np,np_i,np_j
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,nfft3d,n2ft3d
      integer npack1,npack0

*     **** electronic variables ****
      real*8 icharge
      integer ispin
      integer ne(2),n1(2),n2(2),nemax,neq(2),nemaxq
      real*8  en(2)
      real*8  dipole(3)

*     complex*16 psi1(nfft3d,nemax)
*     complex*16 psi2(nfft3d,nemax)
*     real*8     dn(n2ft3d,2)
*     complex*16 Hpsi(nfft3d,nemax)
*     complex*16 psir(nfft3d,nemax)
      integer psi1(2),psi2(2)
      integer occ1(2),occ2(2)
      integer rhoii(2),rhojj(2),rhoij(2)
      integer Hpsi1(2),psi1r(2)
      integer Hpsi2(2),psi2r(2)
      integer vl(2),vlr_l(2),v_field(2),r_grid(2),vc(2)
    

*     ***** energy variables ****
      real*8  E(40)
      integer eig(2),hml1(2),hml2(2),overlap(2),itmp(2)
      integer U_ab(2),V_ab(2),overlap_ab(2)
      integer hml_aa(2),hml_bb(2),hml_ab(2)

*     **** psi smearing block ****
      logical fractional,field_exist
      integer smearoccupation,smeartype
      real*8 smearfermi(2),smearcorrection,smearkT


*     **** error variables ****
      integer ierr

*     **** local variables ****
      real*8 Saa,Sbb,Sab,H1aa,H1bb,H1ab,H2aa,H2bb,H2ab,Vab,Haa,Hbb,Hab
      integer ms,mapping,mapping1d
      real*8  deltae,deltac,deltar
      real*8  gx,gy,gz,cx,cy,cz,sum1,sum2,Eion
      real*8  EV,pi
      integer i,j,k,ia,n,nn
      integer ii,jj,indx,indx1
      real*8 w,sumall,virial,scal2
      integer nfft3
      parameter (nfft3=32)
      character*255 full_filename

      logical value,psi_nogrid
      integer hversion,hnfft(3),hispin,hne(2)
      real*8 hunita(3,3)
      integer ind
      character*50 filename1,filename2

  


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
      logical      pspw_HFX,pspw_HFX_relaxed
      logical      psp_semicore,control_Mulliken
      real*8       psp_rcore,psp_ncore,psp_rlocal
      external     pspw_SIC,pspw_SIC_relaxed,pspw_qmmm_found
      external     pspw_HFX,pspw_HFX_relaxed
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

      real*8   Dneall_m_value,ewald_e,ion_ion_e
      external Dneall_m_value,ewald_e,ion_ion_e
      character*9 ion_amm
      external    ion_amm
      logical  pspw_charge_found
      external pspw_charge_found


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
      call Pack_npack(0,npack0)

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
      if (control_version().eq.3) then
         call coulomb_init()
         call coulomb_screened_init(0,8.0d0,8.0d0)
      end if
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()
      if (control_version().eq.3) call ewald_init()

      if (.not.rtdb_cget(rtdb,'pspw:et:movecs_a',1,filename1))
     >   call util_file_prefix('movecs',filename1)
      if (.not.rtdb_cget(rtdb,'pspw:et:movecs_b',1,filename2))
     >   call util_file_prefix('movecs',filename2)


      call psi_get_header_filename(filename1,
     >                             hversion,hnfft,hunita,ispin,ne)

      if (taskid.eq.MASTER) then
         write(*,*) "filename1=",filename1
         write(*,*) "hversion =",hversion
         write(*,*) "hnfft    =",hnfft
         write(*,*) "hunita   =",hunita
         write(*,*) "ispin    =",ispin
         write(*,*) "ne       =",ne
         write(*,*) 
      end if
      call psi_get_header_filename(filename2,
     >                             hversion,hnfft,hunita,hispin,hne)
      if (taskid.eq.MASTER) then
         write(*,*) "filename2=",filename2
         write(*,*) "hversion =",hversion
         write(*,*) "hnfft    =",hnfft
         write(*,*) "hunita   =",hunita
         write(*,*) "ispin    =",hispin
         write(*,*) "ne       =",hne
         write(*,*) 
      end if
      

   

*     ***** allocate psi2,and psi1 wavefunctions ****
      smearoccupation=0
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



*     *****  read psi1 and psi2 wavefunctions ****
      call psi_read_filename(filename1,ispin,ne,dcpl_mb(psi1(1)))
      call psi_read_filename(filename2,ispin,ne,dcpl_mb(psi2(1)))
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)
      nemax = ne(1)+ne(2)


*     **** allocate other variables *****
      field_exist = pspw_charge_found()
      value = MA_alloc_get(mt_dbl,(2*nemax),'eig',eig(2),eig(1))
      value = value.and.Dneall_m_allocate(0,hml1)
      value = value.and.Dneall_m_allocate(0,hml2)
      value = value.and.Dneall_m_allocate(0,itmp)
      value = value.and.Dneall_m_allocate(0,overlap)
      value = value.and.Dneall_m_allocate(0,U_ab)
      value = value.and.Dneall_m_allocate(0,V_ab)
      value = value.and.Dneall_m_allocate(0,hml_aa)
      value = value.and.Dneall_m_allocate(0,hml_bb)
      value = value.and.Dneall_m_allocate(0,hml_ab)

      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'Hpsi1',Hpsi1(2),Hpsi1(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'Hpsi2',Hpsi2(2),Hpsi2(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,n2ft3d*(neq(1)+neq(2)),
     >                     'psi1r',psi1r(2),psi1r(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,n2ft3d*(neq(1)+neq(2)),
     >                     'psi2r',psi2r(2),psi2r(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,n2ft3d,'rhoii',rhoii(2),rhoii(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,n2ft3d,'rhojj',rhojj(2),rhojj(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,n2ft3d,'rhoij',rhoij(2),rhoij(1))


      if (control_version().eq.4) then
       value = value.and.
     >        MA_alloc_get(mt_dbl,(n2ft3d),'vlr_l',vlr_l(2),vlr_l(1))
       value = value.and.
     >        MA_alloc_get(mt_dbl,(n2ft3d),'vc',vc(2),vc(1))
      end if

      if ((field_exist).or.(control_version().eq.4)) then
       value = value.and.
     >    MA_alloc_get(mt_dbl,(3*n2ft3d),'r_grid',r_grid(2),r_grid(1))
      end if
      value = value.and.
     >  MA_alloc_get(mt_dbl,(n2ft3d),'v_field',v_field(2),v_field(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,(npack0),'vloc',vl(2),vl(1))

      if (.not. value) 
     >   call errquit('pspw_et:out of heap memory',0,MA_ERR)



*     **** initialize QM/MM ****
      call pspw_qmmm_init(rtdb)


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
         write(luout,1121) control_boundry(),control_version()
         if (ispin.eq.1) write(luout,1130) 'restricted'
         if (ispin.eq.2) write(luout,1130) 'unrestricted'

         if (fractional) write(luout,1132)

         write(luout,1140)
         do ia = 1,ion_nkatm()
           write(luout,1150) ia,ion_atom(ia),
     >                    psp_zv(ia),psp_lmax(ia)
           write(luout,1152) psp_lmax(ia)
           write(luout,1153) psp_locp(ia)
           write(luout,1154) psp_lmmax(ia)
           if (control_version().eq.4) write(luout,1156) psp_rlocal(ia)
           if (psp_semicore(ia)) 
     >         write(luout,1155) psp_rcore(ia),psp_ncore(ia)
           write(luout,1151) (psp_rc(i,ia),i=0,psp_lmax(ia))
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


         if (fractional) then
           write(luout,1219) en(1),en(ispin),' (   fractional)'
           write(luout,1221) ne(1),neq(1),
     >                   ne(ispin),neq(ispin),' (fourier space)'
         else
           write(luout,1220) ne(1),neq(1),
     >                   ne(ispin),neq(ispin),' (fourier space)'
           write(luout,1221) ne(1),neq(1),
     >                   ne(ispin),neq(ispin),' (fourier space)'
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

         write(luout,1300)
         write(luout,1800) filename1
         write(luout,1801) filename2
         call util_flush(luout)
      end if

*                |***************************|
******************     start iterations      **********************
*                |***************************|

      if (taskid.eq.MASTER) call current_second(cpu2)

      Saa = 1.0d0
      Sbb = 1.0d0
      call Dneall_ffm_Multiply(0,dcpl_mb(psi1(1)),
     >                           dcpl_mb(psi2(1)),npack1,
     >                           dbl_mb(itmp(1)))
      call Dneall_m_SVD(0,dbl_mb(itmp(1)),dbl_mb(U_ab(1)),
     >                    dbl_mb(eig(1)), dbl_mb(V_ab(1)))

      write(*,*) "V_ab=",dbl_mb(V_ab(1))
      write(*,*) "U_ab=",dbl_mb(U_ab(1))

c      do n=1,nemaxq
c         call Pack_c_Copy(1,dcpl_mb(psi1(1) +(n-1)*npack1),
c     >                      dcpl_mb(Hpsi1(1)+(n-1)*npack1))
c         call Pack_c_Copy(1,dcpl_mb(psi2(1) +(n-1)*npack1),
c     >                      dcpl_mb(Hpsi2(1)+(n-1)*npack1))
c      end do
c      call Dneall_fmf_Multiply(0,dcpl_mb(Hpsi1(1)),npack1,
c     >                            dbl_mb(V_ab(1)),  1.0d0,
c     >                            dcpl_mb(psi1(1)),0.0d0)
c      call Dneall_fmf_Multiply(0,dcpl_mb(Hpsi2(1)),npack1,
c     >                            dbl_mb(U_ab(1)),  1.0d0,
c     >                            dcpl_mb(psi2(1)),0.0d0)

      Sab = 1.0d0
      do n=1,nemax
         Sab = Sab*dbl_mb(eig(1)+n-1)
      end do
      call dcopy(2*nemaxq*npack1,0.0d0,0,dcpl_mb(Hpsi1(1)),1)
      call dcopy(2*nemaxq*npack1,0.0d0,0,dcpl_mb(Hpsi2(1)),1)


      if ((control_version().eq.4).or.(field_exist))
     >   call lattice_r_grid(dbl_mb(r_grid(1)))

*     **** generate local pseudopotential  ****
*     **** and also get force if move true ****
      call v_local(dcpl_mb(vl(1)),.false.,en,dipole)

*     *** long-range psp for charge systems ***
      if (control_version().eq.4) then
         call v_lr_local(dbl_mb(r_grid(1)),dbl_mb(vlr_l(1)))
      end if

*     ***** generate other real-space fields ****
      if (field_exist) then
         call dcopy(n2ft3d,0.0d0,0,dbl_mb(v_field(1)),1)

*        **** generate charge potential ****
         if (pspw_charge_found()) then
            call pspw_charge_Generate_V(n2ft3d,
     >                                  dbl_mb(r_grid(1)),
     >                                  dbl_mb(v_field(1)))
           end if
      end if


      do n=1,nemaxq
         call Pack_c_Copy(1,dcpl_mb(psi1(1)+(n-1)*npack1),
     >                      dbl_mb(psi1r(1)+(n-1)*n2ft3d))
      end do
      call Grsm_gh_fftb(nfft3d,nemaxq,dbl_mb(psi1r(1)))
      call Grsm_h_Zero_Ends(nfft3d,nemaxq,dbl_mb(psi1r(1)))  !*** probably not neeeded!

      do n=1,nemaxq
         call Pack_c_Copy(1,dcpl_mb(psi2(1)+(n-1)*npack1),
     >                      dbl_mb(psi2r(1)+(n-1)*n2ft3d))
      end do
      call Grsm_gh_fftb(nfft3d,nemaxq,dbl_mb(psi2r(1)))
      call Grsm_h_Zero_Ends(nfft3d,nemaxq,dbl_mb(psi2r(1)))  !*** probably not neeeded!

*       **** get Hpsi1,Hpsi2 ****
      if (control_version().eq.3) then
         call psi_H1(ispin,neq,dcpl_mb(psi1(1)),dbl_mb(psi1r(1)),
     >             dcpl_mb(vl(1)),
     >             dbl_mb(v_field(1)),field_exist,
     >             dcpl_mb(Hpsi1(1)))
         call psi_H1(ispin,neq,dcpl_mb(psi2(1)),dbl_mb(psi2r(1)),
     >             dcpl_mb(vl(1)),
     >             dbl_mb(v_field(1)),field_exist,
     >             dcpl_mb(Hpsi2(1)))
      end if

      if (control_version().eq.4) then
         call psi_H1v4(ispin,neq,dcpl_mb(psi1(1)),dbl_mb(psi1r(1)),
     >             dcpl_mb(vl(1)),dbl_mb(vlr_l(1)),
     >             dbl_mb(v_field(1)),field_exist,
     >             dcpl_mb(Hpsi1(1)))
         call psi_H1v4(ispin,neq,dcpl_mb(psi2(1)),dbl_mb(psi2r(1)),
     >             dcpl_mb(vl(1)),dbl_mb(vlr_l(1)),
     >             dbl_mb(v_field(1)),field_exist,
     >             dcpl_mb(Hpsi2(1)))
      end if

*     **** one-electron part ****
      H1aa = 0.0d0
      H1bb = 0.0d0
      H1ab = 0.0d0
      do ms=1,ispin
         do i=n1(ms),n2(ms)
            call Pack_cc_dot(1,dcpl_mb(psi1(1) +(i-1)*npack1),
     >                         dcpl_mb(Hpsi1(1)+(i-1)*npack1),sum1)
            H1aa = H1aa - sum1
           
            call Pack_cc_dot(1,dcpl_mb(psi2(1) +(i-1)*npack1),
     >                         dcpl_mb(Hpsi2(1)+(i-1)*npack1),sum1)
            H1bb = H1bb - sum1

            call Pack_cc_dot(1,dcpl_mb(psi2(1) +(i-1)*npack1),
     >                         dcpl_mb(Hpsi1(1)+(i-1)*npack1),sum1)
            sum2 = 1.0d0
            do j=n1(ms),n2(ms)
               if (j.ne.i) sum2 = sum2*dbl_mb(eig(1)+j-1)
            end do
            H1ab = H1ab - sum1*sum2
         end do
      end do


*     **** ion-ion part ****
      if (control_version().eq.3) Eion = ewald_e()   !**** get ewald energy ****
      if (control_version().eq.4) Eion = ion_ion_e() !**** get free-space ion-ion energy ****



*     **** two-electron part ****
      scal2 = 1.0d0/lattice_omega()
      H2aa = 0.0d0
      H2bb = 0.0d0
      H2ab = 0.0d0
      do ms=1,ispin
         do i=n1(ms),n2(ms)-1
            do j=i+1,n2(ms)
               
               call dcopy(n2ft3d,0.0d0,0,dbl_mb(rhoii(1)),1)
               call dcopy(n2ft3d,0.0d0,0,dbl_mb(rhojj(1)),1)
               call dcopy(n2ft3d,0.0d0,0,dbl_mb(rhoij(1)),1)
               do k=1,n2ft3d
                  dbl_mb(rhoii(1)+k-1) = dbl_mb(rhoii(1)+k-1)
     >                         + scal2*dbl_mb(psi1r(1)+(i-1)*n2ft3d+k-1)
     >                                *dbl_mb(psi1r(1)+(i-1)*n2ft3d+k-1)
                  dbl_mb(rhojj(1)+k-1) = dbl_mb(rhojj(1)+k-1)
     >                         + scal2*dbl_mb(psi1r(1)+(j-1)*n2ft3d+k-1)
     >                                *dbl_mb(psi1r(1)+(j-1)*n2ft3d+k-1)
                  dbl_mb(rhoij(1)+k-1) = dbl_mb(rhojj(1)+k-1)
     >                         + scal2*dbl_mb(psi1r(1)+(i-1)*n2ft3d+k-1)
     >                                *dbl_mb(psi1r(1)+(j-1)*n2ft3d+k-1)
               end do
               if (control_version().eq.4) then
                 call coulomb2_v(dbl_mb(rhojj(1)),dbl_mb(vc(1)))
                 call D3dB_rr_dot(1,dbl_mb(rhoii(1)),dbl_mb(vc(1)),sum1)
                 call coulomb2_v(dbl_mb(rhoij(1)),dbl_mb(vc(1)))
                 call D3dB_rr_dot(1,dbl_mb(rhoij(1)),dbl_mb(vc(1)),sum2)
                 H2aa = H2aa + sum1 - sum2
               end if

               call dcopy(n2ft3d,0.0d0,0,dbl_mb(rhoii(1)),1)
               call dcopy(n2ft3d,0.0d0,0,dbl_mb(rhojj(1)),1)
               call dcopy(n2ft3d,0.0d0,0,dbl_mb(rhoij(1)),1)
               do k=1,n2ft3d
                  dbl_mb(rhoii(1)+k-1) = dbl_mb(rhoii(1)+k-1)
     >                         + scal2*dbl_mb(psi2r(1)+(i-1)*n2ft3d+k-1)
     >                                *dbl_mb(psi2r(1)+(i-1)*n2ft3d+k-1)
                  dbl_mb(rhojj(1)+k-1) = dbl_mb(rhojj(1)+k-1)
     >                         + scal2*dbl_mb(psi2r(1)+(j-1)*n2ft3d+k-1)
     >                                *dbl_mb(psi2r(1)+(j-1)*n2ft3d+k-1)
                  dbl_mb(rhoij(1)+k-1) = dbl_mb(rhojj(1)+k-1)
     >                         + scal2*dbl_mb(psi2r(1)+(i-1)*n2ft3d+k-1)
     >                                *dbl_mb(psi2r(1)+(j-1)*n2ft3d+k-1)
               end do
               if (control_version().eq.4) then
                 call coulomb2_v(dbl_mb(rhojj(1)),dbl_mb(vc(1)))
                 call D3dB_rr_dot(1,dbl_mb(rhoii(1)),dbl_mb(vc(1)),sum1)
                 call coulomb2_v(dbl_mb(rhoij(1)),dbl_mb(vc(1)))
                 call D3dB_rr_dot(1,dbl_mb(rhoij(1)),dbl_mb(vc(1)),sum2)
                 H2bb = H2bb + sum1 - sum2
               end if

               call dcopy(n2ft3d,0.0d0,0,dbl_mb(rhoii(1)),1)
               call dcopy(n2ft3d,0.0d0,0,dbl_mb(rhojj(1)),1)
               call dcopy(n2ft3d,0.0d0,0,dbl_mb(rhoij(1)),1)
               do k=1,n2ft3d
                  dbl_mb(rhoii(1)+k-1) = dbl_mb(rhoii(1)+k-1)
     >                         + scal2*dbl_mb(psi1r(1)+(i-1)*n2ft3d+k-1)
     >                                *dbl_mb(psi2r(1)+(i-1)*n2ft3d+k-1)
                  dbl_mb(rhojj(1)+k-1) = dbl_mb(rhojj(1)+k-1)
     >                         + scal2*dbl_mb(psi1r(1)+(j-1)*n2ft3d+k-1)
     >                                *dbl_mb(psi2r(1)+(j-1)*n2ft3d+k-1)
                  dbl_mb(rhoij(1)+k-1) = dbl_mb(rhojj(1)+k-1)
     >                         + scal2*dbl_mb(psi1r(1)+(i-1)*n2ft3d+k-1)
     >                                *dbl_mb(psi2r(1)+(j-1)*n2ft3d+k-1)
               end do
               if (control_version().eq.4) then
                 call coulomb2_v(dbl_mb(rhojj(1)),dbl_mb(vc(1)))
                 call D3dB_rr_dot(1,dbl_mb(rhoii(1)),dbl_mb(vc(1)),sum1)
                 call coulomb2_v(dbl_mb(rhoij(1)),dbl_mb(vc(1)))
                 call D3dB_rr_dot(1,dbl_mb(rhoij(1)),dbl_mb(vc(1)),sum2)
                 H2ab = H2ab + sum1 - sum2
               end if

            end do
         end do
      end do
 

      Haa = H1aa+H2aa
      Hbb = H1bb+H2bb
      Hab = H1ab+H2ab

      Vab = (Hab-Sab*(Haa+Hbb)/2.0d0)/(1.0d0-Sab*Sab);

      
*::::::::::::::::::::  end of iteration loop  :::::::::::::::::::::::::

      if (taskid.eq.MASTER) call current_second(cpu3)




*                |***************************|
****************** report summary of results **********************
*                |***************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

      if (taskid.eq.MASTER) then
         write(luout,1300)
         write(luout,1410)

         write(luout,*)

         write(luout,1810)
         write(luout,1811) Saa
         write(luout,1812) Sbb
         write(luout,1813) Sab

         write(luout,1820)
         write(luout,1821) H1aa
         write(luout,1822) H1bb
         write(luout,1823) H1ab
         write(luout,1824) Eion

         write(luout,1830)
         write(luout,1831) H2aa
         write(luout,1832) H2bb
         write(luout,1833) H2ab

         write(luout,1840)
         write(luout,1841) Haa
         write(luout,1842) Hbb
         write(luout,1843) Hab

         write(luout,1850) 
         write(luout,1851) Vab


*        **** write out diagonal <psi1|H|psi2> matrix ****
          n = ne(1)
          nn = n*n
          do ms=1,ispin
             if (ms.eq.1) write(luout,1331)
             if (ms.eq.2) write(luout,1332)
             !*** call Gainsville matrix output ***
             call output(dbl_mb(hml1(1)+(ms-1)*nn),
     >                    1,ne(ms),1,ne(ms),
     >                    n,n,1)

           end do

      end if

*                |***************************|
******************         Prologue          **********************
*                |***************************|



*     **** deallocate heap memory ****
      if (control_version().eq.3) call ewald_end()
      call strfac_end()
      if (control_version().eq.3) then
         call coulomb_end()
         call coulomb_screened_end()
      end if
      if (control_version().eq.4) call coulomb2_end()
      call ke_end()
      call mask_end()
      call Pack_end()
      call G_end()
      call ion_end()
      call psp_end()
      call ion_end_FixIon()
      call pspw_qmmm_end()

      value = .true.
      if ((field_exist).or.(control_version().eq.4))
     >   value = value.and.MA_free_heap(r_grid(2))
      value = value.and.MA_free_heap(v_field(2))
      if (control_version().eq.4)
     >   value = value.and.MA_free_heap(vlr_l(2))
      value = value.and.MA_free_heap(vl(2))
      value = value.and.MA_free_heap(psi1r(2))
      value = value.and.MA_free_heap(psi2r(2))
      value = value.and.MA_free_heap(Hpsi1(2))
      value = value.and.MA_free_heap(Hpsi2(2))
      value = value.and.MA_free_heap(eig(2))
      value = value.and.Dneall_m_free(hml1)
      value = value.and.Dneall_m_free(hml2)
      value = value.and.Dneall_m_free(overlap)
      value = value.and.Dneall_m_free(itmp)
      value = value.and.Dneall_m_free(U_ab)
      value = value.and.Dneall_m_free(V_ab)
      value = value.and.Dneall_m_free(hml_aa)
      value = value.and.Dneall_m_free(hml_bb)
      value = value.and.Dneall_m_free(hml_ab)
      value = value.and.MA_free_heap(psi1(2))
      value = value.and.MA_free_heap(psi2(2))
      if (fractional) then
         value = value.and.MA_free_heap(occ1(2))
         value = value.and.MA_free_heap(occ2(2))
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
         write(luout,*)
         write(luout,*) '-----------------'
         write(luout,*) 'cputime in seconds'
         write(luout,*) 'prologue    : ',T1
         write(luout,*) 'main loop   : ',T2
         write(luout,*) 'epilogue    : ',T3
         write(luout,*) 'total       : ',T4
         write(luout,*)
         call nwpw_timing_print_final(.true.,1)
         CALL nwpw_MESSAGE(4)
      end if 


      call Parallel2d_Finalize()
      call Parallel_Finalize()
      pspw_et = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*        Electron-Tranfer Calculation (ET)         *')
 1030 FORMAT(10X,'*     [     Algorithm of Farazdel et al.   ]       *')
 1035 FORMAT(10x,'*     [ NorthWest Chemistry implementation ]       *')
 1040 FORMAT(10X,'*            version #1.00   09/23/2012            *')
 1041 FORMAT(10X,'*    This code was developed by Eric J. Bylaska,   *')
 1042 FORMAT(10X,'*    and Duo Song                                  *')
 1043 FORMAT(10X,'*                                                  *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'=============== PSPW ET input data =================')
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
 1180 FORMAT(/' position of ions:')
 1190 FORMAT(5X, I4, A5, ' (',3F11.5,' ) - atomic mass= ',F7.3,' ',A)
 1191 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F6.3,' - fixed ',A)
 1193 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - z fixed')
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
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
 1271 FORMAT(5x, ' translation contrained')
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

 1700 FORMAT(/' QM/MM-pol-vib/CAV Energies')
 1701 FORMAT( ' --------------------------')
 1702 FORMAT( ' LJ energy              :',E19.10)
 1703 FORMAT( ' Residual Coulomb energy:',E19.10)
 1704 FORMAT( ' MM Vibration energy    :',E19.10)
 1705 FORMAT( ' MM Vibration energy    :',E19.10)
 1706 FORMAT( ' (QM+MM)/Cavity energy  :',E19.10)
                        
 1800 format( ' Reactant  psi filename: ',A)
 1801 format( ' Product   psi filename: ',A)
 1810 format(/' Electron Overlaps:')
 1811 format( '    Reactants/Reactant overlap          S(RR) : ',E19.10)
 1812 format( '    Reactants/Reactant overlap          S(PP) : ',E19.10)
 1813 format( '    Reactants/Products overlap          S(RP) : ',E19.10)
 1820 format(/' One-electron reactants/product interation energies:')
 1821 format( '    Electron energy of reactant         H1(RR): ',E19.10)
 1822 format( '    Electron energy of product          H1(PP): ',E19.10)
 1823 format( '    Electron energy of reactant/product H1(RP): ',E19.10)
 1824 format( '    Ion-Ion energy                      Eion  : ',E19.10)
 1830 format(/' Two-electron reactants/product interation energies:')
 1831 format( '    Electron energy of reactant         H2(RR): ',E19.10)
 1832 format( '    Electron energy of product          H2(PP): ',E19.10)
 1833 format( '    Electron energy of reactant/product H2(RP): ',E19.10)
 1840 format(/' Total reactants/product interation energies:')
 1841 format( '    Electron energy of reactant         H(RR) : ',E19.10)
 1842 format( '    Electron energy of product          H(PP) : ',E19.10)
 1843 format( '    Electron energy of reactant/product H(RP) : ',E19.10)
 1850 format(/' Electron Transfer Coupling:')
 1851 format( '    reactant/product coupling energy    V(RP) : ',E19.10)


!Electronic energy of reactants     H(RR)      -5.3402392824
!Electronic energy of products      H(PP)      -5.3402392824
!
!Reactants/Products overlap         S(RP)      -0.0006033839
!
!Reactants/Products interaction energy:
!-------------------------------------
!One-electron contribution         H1(RP)       0.0040314092
!
!Beginning calculation of 2e contribution
!Two-electron integral screening (tol2e) : 6.03E-11
!
!Two-electron contribution         H2(RP)      -0.0007837138
!Total interaction energy           H(RP)       0.0032476955
!
!Electron Transfer Coupling Energy |V(RP)|      0.0000254810
!                                                      5.592 cm-1
!                                                   0.000693 eV
!                                                      0.016 kcal/mol

 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(luout,9010) ierr
      call Parallel2d_Finalize()
      call Parallel_Finalize()

      pspw_et = value
      return
      END
