*
* $Id$
*
***********************************************************************
*                     band_structure				      *
*                                                                     *
*     This is a developing band structure parallel code for NWCHEM    *
*       + tcgmsg message passing library used                         *
*       + my own slap-decomposed parallel 3d-FFT(real->complex) used  *
*                                                                     *
*                                                                     *
***********************************************************************

      logical function band_structure(rtdb,flag)
      implicit none
      integer rtdb
      integer flag


#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"

      
*     **** parallel variables ****
      integer  taskid,taskid_k,np,np_i,np_j,np_k
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,nfft3d

*     ***** energy variables ****
      integer vall(2),nn,nb,ispin,ne(2),rho(2),neall
      real*8  E(20),en(2),ein,eke
      real*8  dipole(3)
      real*8  stress(3,3)

      integer eigs_dos(2),dosgrid(3),weight_dos(2),pweight_dos(2)
      integer pweight_lmax

*     **** gradient variables ****
      integer fion(2)

*     **** error variables ****
      logical value,ortho,mulliken
      integer ierr

*     **** local variables ****
      logical newpsi,grid3d,spin_orbit,lprint,mprint,hprint
      real*8  gx,gy,gz,cx,cy,cz
      real*8  EV,pi,e1,e2,f0,f1,f2,f3,f4,f5,f6,ttl1
      real*8  pathlength,dist,kold(3),emin,emax,lmbda,rcut
      integer ii,jj,kk,ll,i,k,ia,nion,vers,nbrillioun,icharge,isize,indx
      integer npoints,nbrillall,if1,if2,nbrillq,l3
      integer mapping,tmp(2),norbs_dos,mapping1d
      character*255 full_filename
      character*50 filename
      character*72 cube_comment


*     **** external functions ****
*     **** external functions ****
      integer     cpsp_nprj
      external    cpsp_nprj
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg,ion_amass,ion_TotalCharge
      logical     cpsi_spin_orbit,control_spin_orbit,control_print
      character*4 ion_aname
      integer     control_ispin
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg,ion_amass,ion_TotalCharge
      external    ion_aname,control_print
      external    cpsi_spin_orbit,control_spin_orbit,control_ispin


      real*8   control_tole,control_tolc,control_tolr,ion_rion
      external control_tole,control_tolc,control_tolr,ion_rion
      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      logical  control_read,control_move,ion_init
      external control_read,control_move,ion_init
      integer  cpsp_psp_type
      external cpsp_psp_type
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      integer  ion_nkatm
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      external ion_nkatm
      integer  control_mapping1d
      external control_mapping1d

      character*12 control_boundry
      external     control_boundry

      logical  brillioun_print
      integer  brillioun_nbrillioun,Pneb_nbrillq
      real*8   brillioun_weight_brdcst
      real*8   brillioun_ks_brdcst,brillioun_k_brdcst
      external brillioun_print
      external brillioun_nbrillioun,Pneb_nbrillq
      external brillioun_weight_brdcst
      external brillioun_ks_brdcst,brillioun_k_brdcst
      integer  c_electron_count,linesearch_count
      external c_electron_count,linesearch_count

      real*8   nwpw_timing
      external nwpw_timing
      integer  Cram_nwave_all_brdcst,Cram_nwave_brdcst
      external Cram_nwave_all_brdcst,Cram_nwave_brdcst

      integer  ewald_ncut
      real*8   ewald_rcut,ewald_mandelung,ewald_e
      external ewald_ncut
      external ewald_rcut,ewald_mandelung,ewald_e
      logical  cpsp_semicore,psi_filefind,cpsi_initialize,cpsi_finalize
      external cpsp_semicore,psi_filefind,cpsi_initialize,cpsi_finalize
      logical  cpsi_band_finalize
      external cpsi_band_finalize
      real*8   c_cgsd_noit_energy,cpsi_1energy,cpsi_eigenvalue_brdcst
      external c_cgsd_noit_energy,cpsi_1energy,cpsi_eigenvalue_brdcst
      integer  cpsp_lmax,cpsp_locp,cpsp_lmmax
      external cpsp_lmax,cpsp_locp,cpsp_lmmax
      real*8   cpsp_rcore,cpsp_rc,cpsp_ncore,cpsp_zv
      external cpsp_rcore,cpsp_rc,cpsp_ncore,cpsp_zv
      character*4 ion_atom
      external    ion_atom
      integer  cpsi_ispin,cpsi_ne,psi_get_version
      external cpsi_ispin,cpsi_ne,psi_get_version
      logical  pspw_reformat_c_wvfnc
      external pspw_reformat_c_wvfnc
      integer  control_mapping,control_np_dimensions
      external control_mapping,control_np_dimensions
      integer  control_num_kvectors_structure,control_excited_ne
      external control_num_kvectors_structure,control_excited_ne
      logical  band_dplot_iteration_check,control_Mulliken
      external band_dplot_iteration_check,control_Mulliken
      integer  cpsi_iptr_psi, cpsi_iptr_dn, c_electron_iptr_psir
      external cpsi_iptr_psi, cpsi_iptr_dn, c_electron_iptr_psir
      character   spdf_name
      external    spdf_name
      character*7 c_index_name
      external    c_index_name
      logical     control_only_lda
      external    control_only_lda
    
*****************************|  PROLOGUE  |****************************

      value = .true.
      pi = 4.0d0*datan(1.0d0)

      call nwpw_timing_init()
      call dcopy(20,0.0d0,0,E,1)


*     **** get parallel variables ****
      call Parallel_Init()
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) call current_second(cpu1)

      
      value = control_read(5,rtdb)
      if (.not. value) 
     > call errquit('error reading control',0, DISK_ERR)

      lprint = ((taskid.eq.MASTER).and.(control_print(print_low)))
      mprint = ((taskid.eq.MASTER).and.(control_print(print_medium)))
      hprint = ((taskid.eq.MASTER).and.(control_print(print_high)))
      mulliken = control_Mulliken()

*     ***** print out header ****
      if (mprint) then
         write(luout,1000)
         write(luout,1010)
         write(luout,1020)
         write(luout,1010)
         write(luout,1040)
         write(luout,1010)
         write(luout,1041)
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
      call Parallel3d_taskid_k(taskid_k)


      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0
      mapping = control_mapping()


      ierr = 0 

*     **** initialize C3dB data structure ****
      call C3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call C3dB_nfft3d(1,nfft3d)

      call cpsi_data_init(20)

*     **** read ions ****
      value = ion_init(rtdb)
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
 
*     **** initialize lattice data structure ****
      call lattice_init()
      call c_G_init()

*     **** initalize brillioun ****
      call brillioun_init()
      call Cram_Init()
      call C3dB_pfft_init()


*     **** initialize D3dB data structure and mask for GGA ****
      !if ((control_gga().ge.10).and.(control_gga().lt.100)) THEN
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


*     **** initialize ke,and coulomb data structures ****
      call cstrfac_init()
      call cke_init()
      call c_coulomb_init()
      
      call ewald_init()

*     **** set up phase factors at the current geometry  ****
      call cphafac()
      call cphafac_k()
      call ewald_phafac()

*     **** read in wavefunctions and initialize psi ****
      if (.not.psi_filefind()) then
        call cpsi_new()
        newpsi = .true.

      else
        newpsi = .false.

*       **** convert from pspw format to band format ****
        vers = psi_get_version()
        if ((vers.eq.3).or.(vers.eq.4)) then
           nbrillioun = brillioun_nbrillioun()
           newpsi = .true.
           if (taskid.eq.MASTER) then
             value= pspw_reformat_c_wvfnc(1)
           end if
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
      mapping1d  = control_mapping1d()
      call Pneb_init(ispin,ne,nbrillioun,mapping1d,spin_orbit)
      value = cpsi_initialize(.true.) 

*     **** electron and geodesic data structures ****
      call c_electron_init()
      call c_geodesic_init()
      call linesearch_init()
      call band_dplot_iteration_init()





*                |**************************|
******************   summary of input data  **********************
*                |**************************|

      if (mprint) then
         write(luout,1111) np
         write(luout,1117) np_i,np_j,np_k
         if (mapping.eq.1) write(luout,1112)
         if (mapping.eq.2) write(luout,1113)
         if (mapping.eq.3) write(luout,1118)

         write(luout,1115)
         write(luout,1121) control_boundry(),control_version()

         call v_bwexc_print(luout,control_gga())

         write(luout,1140)
         do ia = 1,ion_nkatm()
           write(luout,1150) ia,ion_atom(ia),
     >                    cpsp_zv(ia),cpsp_lmax(ia)
           write(luout,2000) cpsp_psp_type(ia)
           write(luout,1152) cpsp_lmax(ia)
           write(luout,1153) cpsp_locp(ia)
           write(luout,1154) cpsp_nprj(ia)
           if (cpsp_semicore(ia))
     >         write(luout,1155) cpsp_rcore(ia),cpsp_ncore(ia)
           write(luout,1151) (cpsp_rc(i,ia),i=0,cpsp_lmax(ia))
         end do
         if (control_spin_orbit()) then
           icharge=-cpsi_ne(1)+ion_TotalCharge()
         else
           icharge = -(cpsi_ne(1)+cpsi_ne(cpsi_ispin()))
     >           + ion_TotalCharge()
         end if
         write(luout,1159) icharge

         write(luout,1180)
         write(luout,1190) (I,ion_aname(I),
     >                  (ion_rion(K,I),K=1,3),
     >                  ion_amass(I)/1822.89d0,
     >                 I=1,ion_nion())
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz
         write(luout,1220) cpsi_ne(1),cpsi_ne(cpsi_ispin()),
     >                 ' (Fourier space)'


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
         write(luout,1260) ewald_rcut(),ewald_ncut()
         write(luout,1261) ewald_mandelung()

         ia = brillioun_nbrillioun()
         write(luout,1255)
         write(luout,1256) ia
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
        if (mprint) write(luout,1252) lattice_wcut()
      end if


      if (mprint) then
         write(luout,1270)
         write(luout,1280) control_time_step(),control_fake_mass()
         write(luout,1290) control_tole(),control_tolc()
         write(luout,1300)
         call util_flush(luout)
         call flush(luout)
      end if



      !**** set the size of the band - include virtual orbitals ****
      ispin = cpsi_ispin()
      ne(1) = cpsi_ne(1)+control_excited_ne(1)
      ne(2) = 0
      if (ispin.gt.1) ne(2) = cpsi_ne(2)+control_excited_ne(2)
      

      if (taskid.eq.MASTER) call current_second(cpu2)

*     **** allocate vall ****
      value = BA_alloc_get(mt_dcpl,2*nfft3d,'vall',vall(2),vall(1))
      if (.not. value)
     >  call errquit('band_structure:out of heap memory',0, MA_ERR)



      EV = c_cgsd_noit_energy()
      call c_electron_gen_vall()
      call c_electron_get_vall(dcpl_mb(vall(1)))

      if (taskid.eq.MASTER) then
         write(luout,1600) EV
         write(luout,*)
         write(luout,*) "Self-Consistent Potential Generated"
      end if

      value=cpsi_finalize(.true.)
      call c_electron_finalize()
      call c_geodesic_finalize()
      call ewald_end()
      call cstrfac_end()
      call c_coulomb_end()
      call cke_end()
      call cpsp_end()
      call C3dB_pfft_end()
      call Cram_end()
      call c_G_end()
      call brillioun_end()
      call Pneb_end()


*     **** produce eigenvalue band file(s) ****
      if (ispin.eq.1) then
        call util_file_name('restricted_band',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=58,file=full_filename,form='formatted')
        end if
        call util_file_name('restricted_band_with_k',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=68,file=full_filename,form='formatted')
        end if

      else
        if (cpsi_spin_orbit()) then
        call util_file_name('spinor_band',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=58,file=full_filename,form='formatted')
        end if
        call util_file_name('spinor_band_with_k',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=68,file=full_filename,form='formatted')
        end if

        else
        call util_file_name('alpha_band',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=58,file=full_filename,form='formatted')
        end if
        call util_file_name('beta_band',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=59,file=full_filename,form='formatted')
        end if

        call util_file_name('alpha_band_with_k',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=68,file=full_filename,form='formatted')
        end if
        call util_file_name('beta_band_with_k',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=69,file=full_filename,form='formatted')
        end if

        end if
      end if



*     **** DOS calculation ****
      if (flag.eq.1) then
        if (taskid.eq.MASTER) write(luout,*) "DOS of states calculation"
        call control_dos_grid_structure(dosgrid)



*     **** DOS_dplot calculation ****
      else if (flag.eq.2) then
        if (taskid.eq.MASTER) 
     >   write(luout,*) "DOS_dplot of states calculation"

        call control_dos_grid_structure(dosgrid)
        isize = dosgrid(1)*dosgrid(2)*dosgrid(3)*(ne(1)+ne(2))
        value =           BA_alloc_get(mt_dbl,isize,
     >                       'weight_dos',weight_dos(2),weight_dos(1))
        value = value.and.BA_alloc_get(mt_dbl,isize,
     >                       'eigs_dos',eigs_dos(2),eigs_dos(1))
        value = value.and.BA_alloc_get(mt_dbl,ispin*nfft3d,
     >                       'rho',rho(2),rho(1))
        if (.not. value)
     >   call errquit('band_structure:out of heap memory',0, MA_ERR)

       !**** get eigs_dos from rtdb ****
        if (.not.btdb_get(rtdb,'dos:eigs_dos',mt_dbl,
     >                    isize,dbl_mb(eigs_dos(1)))) then
         call errquit('band_structure:cannot read eigs_dos from rtdb',
     >                0,RTDB_ERR)
        end if

       !**** get dos:ein from rtdb ****
        if (.not.btdb_get(rtdb,'dos:ein',mt_dbl,1,ein)) then
         call errquit('band_structure:cannot read dos:ein from rtdb',
     >                0,RTDB_ERR)
        end if

        call band_dos_weights_generate(dosgrid(1),dosgrid(2),dosgrid(3),
     >                                 dbl_mb(eigs_dos(1)),ne(1),
     >                                 ein,
     >                                 dbl_mb(weight_dos(1)))
        ii = dosgrid(1)*dosgrid(2)*dosgrid(3)*ne(1)
        if (ispin.eq.2)
     >  call band_dos_weights_generate(dosgrid(1),dosgrid(2),dosgrid(3),
     >                                dbl_mb(eigs_dos(1)+ii),ne(2),
     >                                ein,
     >                                dbl_mb(weight_dos(1)+ii))
        call dcopy(ispin*nfft3d,0.0d0,0,dbl_mb(rho(1)),1)

        if (taskid.eq.MASTER)  then
        do k=1,dosgrid(1)*dosgrid(2)*dosgrid(3)
           write(*,*)
           write(*,*) "brillioun k=",k
           write(*,*) "---------------------"
           do ii=1,ne(1)
              write(*,*) "weight=",k,ii,
     >         dbl_mb(weight_dos(1)
     >                +(ii-1)
     >                +(k-1)*ne(1))
              call flush(6)
           end do
        end do
        end if



*     **** band structure calculation ****
      else
         if (taskid.eq.MASTER) 
     >     write(luout,*) "band structure calculation"
      end if

      nbrillall = control_num_kvectors_structure()
      call control_reset_band_structure()
      kk = 0
      do k=1,nbrillall,np_k
        ortho = .true.
        kk = kk + np_k
        if (kk.gt.nbrillall) kk = nbrillall

        if (taskid.eq.MASTER) then
          do ii=k,kk
             write(luout,1301) ii
          end do
        end if


*     **** initialize lattice data structure ****
      call lattice_init()
      call c_G_init()
      call brillioun_structure_init(k,kk-k+1)
      call Cram_Init()
      call C3dB_pfft_init()

*     **** allocate psp data structure and read in psedupotentials into it ****
      call cpsp_init()
      call cpsp_readall()
      if (cpsp_semicore(0)) call c_semicore_check()

*     **** initialize ke,and coulomb data structures ****
      call cstrfac_init()
      call cke_init()
      call c_coulomb_init()
      call ewald_init()

*     **** set up phase factors at the current geometry  ****
      call cphafac()
      call cphafac_k()
      call ewald_phafac()

*     **** read in wavefunctions and initialize psi ****
      if (.not.psi_filefind()) then
        call cpsi_new_ne(ispin,ne)
        newpsi = .true.
        ortho  = .false.

      else
        newpsi = .false.

*       **** convert from pspw format to band format ****
        vers = psi_get_version()
        if ((vers.eq.3).or.(vers.eq.4)) then
           nbrillioun = brillioun_nbrillioun()
           newpsi = .true.
           if (taskid.eq.MASTER) then
             value= pspw_reformat_c_wvfnc(1)
           end if
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
      mapping1d  = control_mapping1d()
      call Pneb_init(ispin,ne,nbrillioun,mapping1d,spin_orbit)
      value = cpsi_initialize(ortho)

c      if (flag.eq.2) call cpsi_dospsi_read(k)


*     **** allocate eigs_dos if first iteration and band structure calculation ****
      if ((flag.eq.1).and.(k.eq.1)) then
        isize = dosgrid(1)*dosgrid(2)*dosgrid(3)*(cpsi_ne(1)+cpsi_ne(2))
        value = BA_alloc_get(mt_dbl,isize,
     >                       'eigs_dos',eigs_dos(2),eigs_dos(1))
        if (mulliken) then
           if (.not.btdb_get(rtdb,'nwpw:dos:orb:norb',
     >                       mt_int,1,norbs_dos)) 
     >        norbs_dos = 0
           value = value.and.BA_alloc_get(mt_dbl,isize*(4+norbs_dos),
     >                      'pweight_dos',pweight_dos(2),pweight_dos(1))
        end if
        if (.not. value)
     >   call errquit('band_structure:out of heap memory',0, MA_ERR)
      end if


*     **** electron and geodesic data structures ****
      call c_electron_init()
      call c_geodesic_init()
      call linesearch_init()

*     **** diagonalize hamiltonian and rotate psi ****
      call c_electron_set_vall(dcpl_mb(vall(1)))

*     **** initialize with steepest descent ***
      call c_sdminimize_noscf(0)

*     **** diagonalize current result ***
      call cpsi_1gen_hml()
      call cpsi_diagonalize_hml()
      call cpsi_1rotate2()
      call cpsi_2to1()

      call cpsi_KS_minimize(1,.false.,control_tole(),control_tolc())

      if (control_print(print_high)) call cpsi_check_indx(k)

      if (band_dplot_iteration_check(k)) then
        call band_dplot_iteration(k,ispin,ne,1,
     >                            dcpl_mb(cpsi_iptr_psi(1)),
     >                            dbl_mb(cpsi_iptr_dn(1)),
     >                            dcpl_mb(c_electron_iptr_psir()))
      end if

     
      NN=cpsi_ne(1)-cpsi_ne(2)
      EV=27.2116d0
      do nb=1,brillioun_nbrillioun()
         
         f1 = brillioun_ks_brdcst(1,nb)
         f2 = brillioun_ks_brdcst(2,nb)
         f3 = brillioun_ks_brdcst(3,nb)
         f4 = brillioun_k_brdcst(1,nb)
         f5 = brillioun_k_brdcst(2,nb)
         f6 = brillioun_k_brdcst(3,nb)
         if ((k+nb).eq.2) then
            pathlength = 0.0d0
         else
            dist=dsqrt((f4-kold(1))**2+(f5-kold(2))**2+(f6-kold(3))**2)
            pathlength = pathlength + dist
         end if
         kold(1) = f4
         kold(2) = f5
         kold(3) = f6

         if (taskid.eq.MASTER) then
            write(luout,1508) k+nb-1,pathlength,f1,f2,f3,f4,f5,f6
            write(luout,1500)
         end if

         !*** flag==2 ***
         if (flag.eq.2) then

            !*** spin-orbit ****
            if (cpsi_spin_orbit()) then
               do i=0,cpsi_ne(1)-1
                  ii=cpsi_ne(1)-i
                  e1 = cpsi_eigenvalue_brdcst(nb,1,cpsi_ne(1)-i)
                  if (taskid.eq.MASTER)
     >            write(luout,1511)  e1,e1*EV,
     >                     dbl_mb(weight_dos(1)+ii+(k+nb-1)*cpsi_ne(1))
               end do

            !*** not spin-orbit ****
            else
               do i=0,NN-1
                  ii = cpsi_ne(1)-i
                  e1 = cpsi_eigenvalue_brdcst(nb,1,cpsi_ne(1)-i)
                  if (taskid.eq.MASTER) 
     >            write(luout,1511) e1,e1*EV,
     >                     dbl_mb(weight_dos(1)+(ii-1)+(k+nb-1)*ne(1))
               end do
               do i=0,cpsi_ne(2)-1
                  ii = cpsi_ne(1)-i
                  jj = cpsi_ne(2)-i +
     >               dosgrid(1)*dosgrid(2)*dosgrid(3)*ne(1)
                  e1 = cpsi_eigenvalue_brdcst(nb,1,cpsi_ne(1)-i-NN)
                  e2 = cpsi_eigenvalue_brdcst(nb,2,cpsi_ne(2)-i)
                  if (taskid.eq.MASTER)
     >            write(luout,1511)  e1,e1*EV,
     >                      dbl_mb(weight_dos(1)+(ii-1)+(k+nb-2)*ne(1)),
     >                      e2,e2*EV,
     >                      dbl_mb(weight_dos(1)+(jj-1)+(k+nb-2)*ne(2))
               end do
            end if

         !*** flag!=2 ***
         else

            !*** spin-orbit ****
            if (cpsi_spin_orbit()) then
               do i=0,cpsi_ne(1)-1
                  e1 = cpsi_eigenvalue_brdcst(nb,1,cpsi_ne(1)-i)
                  if (taskid.eq.MASTER)
     >            write(luout,1510)  e1,e1*EV
               end do

            !*** not spin-orbit ****
            else
               do i=0,NN-1
                  e1 = cpsi_eigenvalue_brdcst(nb,1,cpsi_ne(1)-i)
                  if (taskid.eq.MASTER)
     >            write(luout,1510) e1,e1*EV
               end do
               do i=0,cpsi_ne(2)-1
                  e1 = cpsi_eigenvalue_brdcst(nb,1,cpsi_ne(1)-i-NN)
                  e2 = cpsi_eigenvalue_brdcst(nb,2,cpsi_ne(2)-i)
                  if (taskid.eq.MASTER)
     >            write(luout,1510)  e1,e1*EV,e2,e2*EV
               end do
            end if

         end if

         !*** set eigs_dos ***
         if (flag.eq.1) then
            do i=1,cpsi_ne(1)
               indx = eigs_dos(1) 
     >              + (k+nb-2) 
     >              + (i-1)*dosgrid(1)*dosgrid(2)*dosgrid(3)
               dbl_mb(indx) = cpsi_eigenvalue_brdcst(nb,1,i)
            end do
            do i=1,cpsi_ne(2)
               indx = eigs_dos(1) 
     >              + (k+nb-2) 
     >              + (i-1+cpsi_ne(1))*dosgrid(1)*dosgrid(2)*dosgrid(3)
               dbl_mb(indx) = cpsi_eigenvalue_brdcst(nb,2,i)
            end do
         end if
         
         if (.not.BA_push_get(mt_dbl,cpsi_ne(1),'tmp',tmp(2),tmp(1)))
     >      call errquit('band_structure:push stack',99,MA_ERR)

         do i=1,cpsi_ne(1)
            dbl_mb(tmp(1)+i-1) = cpsi_eigenvalue_brdcst(nb,1,i)
         end do
         if (taskid.eq.MASTER)
     >   write(58,'(1000E14.6)') pathlength,
     >          (dbl_mb(tmp(1)+i-1),i=1,cpsi_ne(1))
         if (taskid.eq.MASTER)
     >   write(68,'(1000E14.6)') pathlength,f1,f2,f3,f4,f5,f6,
     >          (dbl_mb(tmp(1)+i-1),i=1,cpsi_ne(1))

         if ((.not.cpsi_spin_orbit()).and.(ispin.eq.2)) then
           do i=1,cpsi_ne(2)
              dbl_mb(tmp(1)+i-1)=cpsi_eigenvalue_brdcst(nb,2,i)
           end do
           if (taskid.eq.MASTER)
     >     write(59,'(1000E14.6)') pathlength,
     >          (dbl_mb(tmp(1)+i-1),i=1,cpsi_ne(2))
           if (taskid.eq.MASTER)
     >     write(69,'(1000E14.6)') pathlength,f1,f2,f3,f4,f5,f6,
     >          (dbl_mb(tmp(1)+i-1),i=1,cpsi_ne(2))
         end if         


         if (.not.BA_pop_stack(tmp(2)))
     >      call errquit('band_structure:pop stack',99,MA_ERR)

      end do !*** nb ***

      !*** set rho ***
      if ((flag.eq.2).and.((k+taskid_k).le.nbrillall)) then
        ii = (k+taskid_k-1)*ne(1)
        call c_electron_gen_weighted_density(1,
     >                               dbl_mb(weight_dos(1)+ii),
     >                               dbl_mb(rho(1)))
        ii = (k+taskid_k-1)*ne(2)+dosgrid(1)*dosgrid(2)*dosgrid(3)*ne(1)
        if (ispin.eq.2)
     >  call c_electron_gen_weighted_density(2,
     >                               dbl_mb(weight_dos(1)+ii),
     >                               dbl_mb(rho(1)+nfft3d))
      end if

*     **** set pweight_dos ****
      if ((flag.eq.1).and.mulliken) then
         call cpsi_projected_dos_weights(rtdb,dosgrid,k,
     >                                   dbl_mb(pweight_dos(1)),
     >                                   pweight_lmax,norbs_dos)
      end if


*     **** writeout dospsi -- needed for task band dos_dplot ****
cccc      if (flag.eq.1) call cpsi_dospsi_write(k)


*     **** finalize and deallocate cpsi ****
      value = cpsi_finalize(.true.)

*     **** deallocate heap memory ****
      call c_electron_finalize()
      call c_geodesic_finalize()
      call ewald_end()
      call cstrfac_end()
      call c_coulomb_end()
      call cke_end()
      call cpsp_end()
      call C3dB_pfft_end()
      call Cram_end()
      call c_G_end()
      call brillioun_end()
      call Pneb_end()

      end do
      if (taskid.eq.MASTER) then
        close(58)
        close(68)
        if (ispin.eq.2) then
           close(59)
           close(69)
        end if
      end if

      if (taskid.eq.MASTER) call current_second(cpu3)

*                |***************************|
******************       DOS plotting        **********************
*                |***************************|

      value = btdb_parallel(.false.)
      if ((flag.eq.1).and.(taskid.eq.MASTER)) then


        if (.not.btdb_get(rtdb,'dos:npoints',mt_int,1,npoints)) then
          npoints = 500
        end if
  
        if (.not.btdb_get(rtdb,'dos:emin',mt_dbl,1,emin)) then
           emin = 99999.0d0
           do ii=1,(ne(1)+ne(2))*dosgrid(1)*dosgrid(2)*dosgrid(3)
             if (dbl_mb(eigs_dos(1)+ii-1).lt.emin) 
     >         emin = dbl_mb(eigs_dos(1)+ii-1)
           end do
           emin = emin - 0.1d0
        end if
  
        if (.not.btdb_get(rtdb,'dos:emax',mt_dbl,1,emax)) then
           emax = -99999.0d0
           do ii=1,(ne(1)+ne(2))*dosgrid(1)*dosgrid(2)*dosgrid(3)
             if (dbl_mb(eigs_dos(1)+ii-1).gt.emax) 
     >         emax = dbl_mb(eigs_dos(1)+ii-1)
           end do
           emax = emax + 0.1d0
        end if

        call util_file_name('dos',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        open(unit=58,file=full_filename,form='formatted')
        call band_dos_generate(58,dosgrid(1),dosgrid(2),dosgrid(3),
     >                           dbl_mb(eigs_dos(1)),ne(1)+ne(2),
     >                           (2.0d0*(3-ispin)),
     >                           npoints,emin,emax)
        close(58)

        if (ispin.eq.2) then
           call util_file_name('dos_alpha',
     >                       .false.,
     >                       .false.,
     >                       full_filename)
           open(unit=58,file=full_filename,form='formatted')
           call band_dos_generate(58,dosgrid(1),dosgrid(2),dosgrid(3),
     >                              dbl_mb(eigs_dos(1)),ne(1),
     >                              (1.0d0),
     >                              npoints,emin,emax)
           close(58)
           call util_file_name('dos_beta',
     >                       .false.,
     >                       .false.,
     >                       full_filename)
           open(unit=58,file=full_filename,form='formatted')
           ii = dosgrid(1)*dosgrid(2)*dosgrid(3)*ne(1)
           call band_dos_generate(58,dosgrid(1),dosgrid(2),dosgrid(3),
     >                              dbl_mb(eigs_dos(1)+ii),ne(2),
     >                              (-1.0d0),
     >                              npoints,emin,emax)
           close(58)
        end if

        if (mulliken) then
        do ll=0,pweight_lmax
           call util_file_name('dos_both_'//spdf_name(ll),
     >                       .false.,
     >                       .false.,
     >                       full_filename)
           open(unit=58,file=full_filename,form='formatted')
           jj=dosgrid(1)*dosgrid(2)*dosgrid(3)*(ne(1)+ne(2))*ll
           call band_projected_dos_generate(58,
     >                           dosgrid(1),dosgrid(2),dosgrid(3),
     >                           dbl_mb(eigs_dos(1)),
     >                           dbl_mb(pweight_dos(1)+jj),ne(1)+ne(2),
     >                           (1.0d0*(3-ispin)),
     >                           npoints,emin,emax)
           close(58)

           if (ispin.eq.2) then
              call util_file_name('dos_alpha_'//spdf_name(ll),
     >                       .false.,
     >                       .false.,
     >                       full_filename)
              open(unit=58,file=full_filename,form='formatted')
              jj=dosgrid(1)*dosgrid(2)*dosgrid(3)*(ne(1)+ne(2))*ll
              call band_projected_dos_generate(58,
     >                              dosgrid(1),dosgrid(2),dosgrid(3),
     >                              dbl_mb(eigs_dos(1)),
     >                              dbl_mb(pweight_dos(1)+jj),ne(1),
     >                              (1.0d0),
     >                              npoints,emin,emax)
              close(58)
              call util_file_name('dos_beta_'//spdf_name(ll),
     >                       .false.,
     >                       .false.,
     >                       full_filename)
              open(unit=58,file=full_filename,form='formatted')
              ii = dosgrid(1)*dosgrid(2)*dosgrid(3)*ne(1)
              jj=dosgrid(1)*dosgrid(2)*dosgrid(3)*(ne(1)+ne(2))*ll + ii
              call band_projected_dos_generate(58,
     >                              dosgrid(1),dosgrid(2),dosgrid(3),
     >                              dbl_mb(eigs_dos(1)+ii),
     >                              dbl_mb(pweight_dos(1)+jj),ne(2),
     >                              (-1.0d0),
     >                              npoints,emin,emax)
              close(58)
           end if
        end do

        !*** ORBITAL DOS ***
         do ll=1,norbs_dos
            l3 = ll+pweight_lmax+1
            call util_file_name('dos_both_orb'//c_index_name(ll),
     >                       .false.,
     >                       .false.,
     >                       full_filename)
             open(unit=58,file=full_filename,form='formatted')
             jj=dosgrid(1)*dosgrid(2)*dosgrid(3)*(ne(1)+ne(2))*l3
             call band_projected_dos_generate(58,
     >                           dosgrid(1),dosgrid(2),dosgrid(3),
     >                           dbl_mb(eigs_dos(1)),
     >                           dbl_mb(pweight_dos(1)+jj),ne(1)+ne(2),
     >                           (1.0d0*(3-ispin)),
     >                           npoints,emin,emax)
             close(58)

             if (ispin.eq.2) then
                call util_file_name('dos_alpha_orb'//c_index_name(ll),
     >                       .false.,
     >                       .false.,
     >                       full_filename)
                open(unit=58,file=full_filename,form='formatted')
                jj=dosgrid(1)*dosgrid(2)*dosgrid(3)*(ne(1)+ne(2))*l3
                call band_projected_dos_generate(58,
     >                              dosgrid(1),dosgrid(2),dosgrid(3),
     >                              dbl_mb(eigs_dos(1)),
     >                              dbl_mb(pweight_dos(1)+jj),ne(1),
     >                              (1.0d0),
     >                              npoints,emin,emax)
                close(58)
                call util_file_name('dos_beta_orb'//c_index_name(ll),
     >                       .false.,
     >                       .false.,
     >                       full_filename)
                open(unit=58,file=full_filename,form='formatted')
                ii = dosgrid(1)*dosgrid(2)*dosgrid(3)*ne(1)
                jj=dosgrid(1)*dosgrid(2)*dosgrid(3)*(ne(1)+ne(2))*l3+ii
                call band_projected_dos_generate(58,
     >                              dosgrid(1),dosgrid(2),dosgrid(3),
     >                              dbl_mb(eigs_dos(1)+ii),
     >                              dbl_mb(pweight_dos(1)+jj),ne(2),
     >                              (-1.0d0),
     >                              npoints,emin,emax)
                close(58)
             end if

         end do


        end if


        !**** put eigs_dos on rtdb for use by task band dos_dplot ***
        if (.not.btdb_put(rtdb,'dos:eigs_dos',mt_dbl,
     >                    isize,dbl_mb(eigs_dos(1)))) then
         call errquit('band_structure:cannot write eigs_dos to rtdb',
     >                0,RTDB_ERR)
        end if
      end if
      value = btdb_parallel(.true.)


*                |***************************|
******************     DOS_dplot plotting    **********************
*                |***************************|
      if (flag.eq.2) then


        grid3d = .false.
        if (btdb_get(rtdb,'band_dplot:3d_grid:nx',mt_int,1,i))
     >    grid3d = .true.

        if (.not.btdb_cget(rtdb,'dos:dplot_up',1,filename)) then
           filename     = 'dos_up.cube '
        end if
        indx = index(filename,' ') - 1
        write(cube_comment,'(A,F8.3)') "dos up density, e=",ein
        write(*,*) '   writing dos up density E=',ein,
     >                  ' to filename: ',filename(1:11)
        if (grid3d) then
             call band_dplot_gcube_write3d(rtdb,filename,
     >                         -1,cube_comment,dbl_mb(rho(1)))
        else
             call band_dplot_gcube_write(rtdb,filename,
     >                         -1,cube_comment,dbl_mb(rho(1)))
        endif

        if (ispin.eq.2) then
          if (.not.btdb_cget(rtdb,'dos:dplot_dn',1,filename)) then
             filename     = 'dos_dn.cube '
          end if
          indx = index(filename,' ') - 1
          write(cube_comment,'(A,F8.3)') "dos down density, e=",ein
          write(*,*) '   writing dos down density E=',ein,
     >                  ' to filename: ',filename(1:11)
          if (grid3d) then
             call band_dplot_gcube_write3d(rtdb,filename,
     >                        -2,cube_comment,dbl_mb(rho(1)+nfft3d))
          else
             call band_dplot_gcube_write(rtdb,filename,
     >                        -2,cube_comment,dbl_mb(rho(1)+nfft3d))
          endif
        end if

      end if
    

*                |***************************|
******************         Epilogue          **********************
*                |***************************|

      

*     **** deallocate heap memory ****

      value = BA_free_heap(vall(2)) 

      !*** deallocate eigs_dos and pweight_dos ****
      if (flag.eq.1)  then 
         value = value.and.BA_free_heap(eigs_dos(2)) 
         if (mulliken) then
            value = value.and.BA_free_heap(pweight_dos(2)) 
         end if
      end if

      if (flag.eq.2) then
        value = value.and.BA_free_heap(weight_dos(2)) 
        value = value.and.BA_free_heap(eigs_dos(2)) 
        value = value.and.BA_free_heap(rho(2)) 
      end if
      
      call ion_write(rtdb)
      call ion_end()
      call cpsi_data_end()
      call C3dB_end(1)
      !IF ((control_gga().ge.10).and.(control_gga().lt.100)) THEN
      if (.not.control_only_lda()) then
         call mask_end()
         call G_end()
         call D3dB_end(1)
      end if


*                |***************************|
****************** report consumed cputime   **********************
*                |***************************|
      if (lprint) then
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
         write(luout,*) '-------------------------------'
         write(luout,*) 'Time spent doing:'
         write(luout,*) '  FFTs                       : ', 
     >                          nwpw_timing(1)
         write(luout,*) '  dot products               : ', 
     >                          nwpw_timing(2)
         write(luout,*) '  geodesic                   : ', 
     >                          nwpw_timing(10)
         write(luout,*) '  exchange correlation       : ', 
     >                          nwpw_timing(4)
         write(luout,*) '  local pseudopotentials     : ', 
     >                          nwpw_timing(5)
         write(luout,*) '  non-local pseudopotentials : ', 
     >                          nwpw_timing(6)
         write(luout,*) '  hartree potentials         : ', 
     >                          nwpw_timing(7)
         write(luout,*) '  structure factors          : ', 
     >                          nwpw_timing(8)
         write(luout,*) '  masking and packing        : ', 
     >                          nwpw_timing(9)
         write(luout,*)
         CALL nwpw_MESSAGE(4)
      end if 

      call Parallel3d_Finalize()
      call Parallel_Finalize()
      band_structure = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*           NWPW Band Structure Calculation        *')
 1040 FORMAT(10X,'*            version #2.00   1/20/07               *')
 1041 FORMAT(10X,'*          Developed by Eric J. Bylaska            *')
 1043 FORMAT(10X,'*          and Patrick Nichols                     *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1111 FORMAT(/' number of processors used:',I16)
 1112 FORMAT( ' parallel mapping         :         1d slab')
 1113 FORMAT( ' parallel mapping         :      2d hilbert')
 1115 FORMAT(/' options:')
 1117 FORMAT( ' processor grid           :',I4,' x',I4,' x',I4)
 1118 FORMAT( ' parallel mapping         :       2d hcurve')
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1140 FORMAT(/' elements involved in the cluster:')
 1150 FORMAT(5X,I2,': ',A4,'  core charge:',F4.1,'  lmax=',I1)
 1151 FORMAT(5X,'        cutoff =',4F8.3)
 1152 FORMAT(12X,' highest angular component      : ',i2)
 1153 FORMAT(12X,' local potential used           : ',i2)
 1154 FORMAT(12X,' number of non-local projections: ',i2)
 1155 FORMAT(12X,' semicore corrections included  : ',
     >       F6.3,' (radius) ',F6.3,' (charge)')
 1156 FORMAT(12X,' aperiodic cutoff radius        : ',F6.3)
 1159 FORMAT(/' total charge:',I2)
 1160 FORMAT(/' atomic composition:')
 1170 FORMAT(7(5X,A4,':',I3))
 1180 FORMAT(/' initial position of ions:')
 1190 FORMAT(5X, I4, A5  ,' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1220 FORMAT(/' number of electrons: spin up=',I4,'  spin down=',I4,A)
 1230 FORMAT(/' supercell:')
 1231 FORMAT(5x,' volume : ',F10.1)
 1241 FORMAT(5x,' lattice:    a1=<',3f8.3,' >')
 1242 FORMAT(5x,'             a2=<',3f8.3,' >')
 1243 FORMAT(5x,'             a3=<',3f8.3,' >')
 1244 FORMAT(5x,' reciprocal: b1=<',3f8.3,' >')
 1245 FORMAT(5x,'             b2=<',3f8.3,' >')
 1246 FORMAT(5x,'             b3=<',3f8.3,' >')

 1249 FORMAT(/' computational grids:')
 1250 FORMAT(5X,' density      cutoff=',F7.3,'  fft=',I4,'x',I4,'x',I4,
     &       '( ',I8,' waves ',I8,' per task)')
 1251 FORMAT(5X,' wavefnc ',I4,' cutoff=',F7.3,
     &        '  fft=',I4,'x',I4,'x',I4,
     &       '( ',I8,' waves ',I8,' per task)')
 1252 FORMAT(5x,' wavefnc     cutoff=',F7.3,
     >       ' wavefunction grids not printed - ',
     >       'number of k-points is very large')
 1255 FORMAT(/' brillouin zone:')
 1256 FORMAT(5x,'number of zone points:',I6)
 1257 FORMAT(5x,' weight=',f8.3,'  ks=<',3f8.3,' >, k=<',3f8.3,'>')
 1258 FORMAT(5x,' number of k-points is very large')

 1260 FORMAT(5X,' Ewald summation: cut radius=',F8.2,'  and',I3)
 1261 FORMAT(5X,'                   madelung=',f14.8)
 
 1270 FORMAT(/' technical parameters:')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'fictitious mass=',F10.1)
 1290 FORMAT(5X, ' tolerance=',E8.3,' (energy)',E12.3,
     &        ' (density)')
 1300 FORMAT(//)
 1301 FORMAT(//'== Optimizing Brillouin Zone Point:',I6,' =='/)
 1304 FORMAT(/)
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
 1508 FORMAT(/' Brillouin zone point: ',i6,
     >       /'pathlength=',f10.6,
     >       /'    k     =<',3f8.3,'> . <b1,b2,b3> ',
     >       /'          =<',3f8.3,'>')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1511 FORMAT(2(E18.7,' (',F8.3,'eV)  dplot weight=',F8.3))
 1600 FORMAT(/' Total BAND energy   :',E19.10)
 2000 FORMAT(12X,' pseudpotential type            : ',i2)
 2005 FORMAT(12x,' number of non local projectors : ',i3)
 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(6,9010) ierr
      call Parallel_Finalize()

      band_structure = value
      return
      END

**************** Definitions of Cubes and Tetrahedrons *****************
*                                                                      *
*                                                                      *
*     (011)------------(111)                    ( 3 )------------( 7 ) *
*       +                +                        +                +   *
*      /.               /|                       /.               /|   *
*     / .              / |                      / .              / |   *
*    /  .             /  |                     /  .             /  |   *
*   /   .            /   |                    /   .            /   |   *
* (001)------------(101) |      <====>      ( 1 )------------( 5 ) |   *
*   |   .            |   |                    |   .            |   |   *
*   | (010)..........|.(110)                  | ( 2 )..........|.( 6 ) *
*   |   .            |   /                    |   .            |   /   *
*   |  .             |  /                     |  .             |  /    *
*   | .              | /                      | .              | /     *
*   |.               |/                       |.               |/      *
*   +                +                        +                +       *
* (000)------------(100)                    ( 0 )------------( 4 )     *
*                                                                      *
*                                                                      *
* Algorithm to find diagaonals                                         *
*                                                                      *
*  Given a cube vertice d1                                             *
*  then d2 = d1^(111) = d1^7                                           *
*                                                                      *
*   Where the cOR bit operator "^" is defined as follows:              *
*      0^0 = 0                                                         *
*      1^1 = 0                                                         *
*      1^0 = 1                                                         *
*      0^1 = 1                                                         *
*                                                                      *
* The four possible cube diagonals are                                 *
*     (000) --- (111)                              (0, 7)              *
*     (001) --- (110)           <====>  2-tuple    (1, 6)              *
*     (010) --- (101)                   rep.       (2, 5)              *
*     (011) --- (100)                              (3, 4)              *
*                                                                      *
* Given a 2-tuple (d1,d2) that defines the diagonal of the cube,       *
* six tetrahedrons are defined, e.g.                                   *
*                                                                      *
*                      (111)                                           *
*                     .  / .                                           *
*                   .   /  .                                           *
*                 .    /   .                                           *
*                .    /   .                                            *
*              .     /    .                                            *
*             .    (101)  .                                            *
*           .     /  |   .        <====> 4-tuple (0, 7, 4, 5)          *
*          .    /    |   .               rep.                          *
*        .    /      |   .                                             *
*       .   /        |  .                                              *
*     .   /          |  .                                              *
*    .  /            |  .                                              *
*  .  /              | .                                               *
* (000)------------(100)                                               *
*                                                                      *
*                                                                      *
* Algorithm to find the six tetradedrons                               *
*                                                                      *
*  Given the diagonals vertices d1 and d2 such that d2=d1^7, the six   * 
*  tetradedrons (six 4-tuples) can be found using the following        *
*  algorithm:                                                          *
*                                                                      *
*   shift(0) = (001) = 1                                               *
*   shift(1) = (010) = 2                                               *
*   shift(2) = (100) = 4                                               *
*   tcount = 0                                                         *
*   For i=0,2                                                          *
*   For j=0,2                                                          *
*     c1 = d1^shift(i)                                                 *
*     c2 = c1^shift(j)                                                 *
*     If (c1 != d1) and (c1 != d2) and (c2!=d1) and (c2!=d2) Then      *
*       tetra(tcount) = (d1,d2,c1,c2)                                  *
*       tcount = tcount + 1                                            *
*     End If                                                           *
*   End For                                                            *
*   End For                                                            *
*                                                                      *
**************** Definitions of Cubes and Tetrahedrons *****************


*     *********************************************
*     *                                           *
*     *            band_dos_generate              *
*     *                                           *
*     *********************************************

      subroutine band_dos_generate(unit,idx,idy,idz,eigs,neigs,
     >                             sign,npoints,emin,emax)
      implicit none
      integer unit
      integer idx,idy,idz
      real*8 eigs(idx,idy,idz,*)
      integer neigs
      real*8  sign
      integer npoints
      real*8 emin,emax

*     **** local variables ****
      integer dosgrid(3)
      integer i,j,k,ii,jj,kk,ncubes,ntetra,count
      integer ishft,jshft,kshft
      integer k1_d(4),k2_d(4),k3_d(4),k1_dd(4),k2_dd(4),k3_dd(4)
      integer id,d1(4),d2(4)
      integer itetra(4,6)
      real*8  VT,VG
      real*8  B(3,3),unitg(3,3),e,ecube(8),f,g,de
      real*8  k1,k2,k3,kx,ky,kz,kkx,kky,kkz,r,rmax

*     **** external functions ****
      real*8   lattice_unitg,Dstates_Cube,Nstates_Cube
      external lattice_unitg,Dstates_Cube,Nstates_Cube

       dosgrid(1) = idx
       dosgrid(2) = idy
       dosgrid(3) = idz

c      write(unit,*) "dosgrid:",dosgrid
c      write(unit,*) "neigs:     ",neigs
c      write(unit,*) "sign:      ",sign
c      write(unit,*) "npoints:   ",npoints
c      write(unit,*) "emin:      ", emin
c      write(unit,*) "emax:      ", emax

      do j=1,3
      do i=1,3
        B(i,j) = lattice_unitg(i,j)
      end do  
      end do  

*     **** volume of reciprocal unit cell, VG ****
      unitg(1,1) = B(2,2)*B(3,3) - B(3,2)*B(2,3)
      unitg(2,1) = B(3,2)*B(1,3) - B(1,2)*B(3,3)
      unitg(3,1) = B(1,2)*B(2,3) - B(2,2)*B(1,3)

      unitg(1,2) = B(2,3)*B(3,1) - B(3,3)*B(2,1)
      unitg(2,2) = B(3,3)*B(1,1) - B(1,3)*B(3,1)
      unitg(3,2) = B(1,3)*B(2,1) - B(2,3)*B(1,1)

      unitg(1,3) = B(2,1)*B(3,2) - B(3,1)*B(2,2)
      unitg(2,3) = B(3,1)*B(1,2) - B(1,1)*B(3,2)
      unitg(3,3) = B(1,1)*B(2,2) - B(2,1)*B(1,2)
      VG = B(1,1)*unitg(1,1)
     >   + B(2,1)*unitg(2,1)
     >   + B(3,1)*unitg(3,1)
      
      ncubes = dosgrid(1)*dosgrid(2)*dosgrid(3)
      ntetra = ncubes*6
      VT = VG/dble(ntetra)

c      write(unit,*) "VG:     ",VG
c      write(unit,*) "number of cubes:",ncubes
c      write(unit,*) "number of tetra:",ntetra
c      write(unit,*) "VT:     ",VT
c
c      count = 0
c      do k=0,dosgrid(3)-1
c      do j=0,dosgrid(2)-1
c      do i=0,dosgrid(1)-1
c         count = count + 1
c         k1 = (dble(i)/dble(dosgrid(1)))
c         k2 = (dble(j)/dble(dosgrid(2)))
c         k3 = (dble(k)/dble(dosgrid(3)))
c         kx = k1*B(1,1) + k2*B(1,2) + k3*B(1,3)
c         ky = k1*B(2,1) + k2*B(2,2) + k3*B(2,3)
c         kz = k1*B(3,1) + k2*B(3,2) + k3*B(3,3)
c         write(unit,*) i,j,k
c         write(unit,3508) count,k1,k2,k3,kx,ky,kz
c         write(unit,*)
c      end do
c      end do
c      end do

*     ********************************
*     **** find shortest diagonal ****
*     ********************************

*     **** (000) ---- (111) ****
      k1_d(1) = 0
      k2_d(1) = 0
      k3_d(1) = 0
      k1_dd(1) = 1
      k2_dd(1) = 1
      k3_dd(1) = 1
      d1(1) = 0
      d2(1) = 7

*     **** (001) ---- (110) ****
      k1_d(2) = 1
      k2_d(2) = 0
      k3_d(2) = 0
      k1_dd(2) = 0
      k2_dd(2) = 1
      k3_dd(2) = 1
      d1(2) = 1
      d2(2) = 6

*     **** (010) ---- (101) ****
      k1_d(3) = 0
      k2_d(3) = 1
      k3_d(3) = 0
      k1_dd(3) = 1
      k2_dd(3) = 0
      k3_dd(3) = 1
      d1(3) = 2
      d2(3) = 5

*     **** (011) ---- (100) ****
      k1_d(4) = 1
      k2_d(4) = 1
      k3_d(4) = 0
      k1_dd(4) = 0
      k2_dd(4) = 0
      k3_dd(4) = 1
      d1(4) = 3
      d2(4) = 4

      id = 1
      rmax = 9.99d9
      do i=1,4
         kx = k1_d(i)*B(1,1) + k2_d(i)*B(1,2) + k3_d(i)*B(1,3)
         ky = k1_d(i)*B(2,1) + k2_d(i)*B(2,2) + k3_d(i)*B(2,3)
         kz = k1_d(i)*B(3,1) + k2_d(i)*B(3,2) + k3_d(i)*B(3,3)

         kkx = k1_dd(i)*B(1,1) + k2_dd(i)*B(1,2) + k3_dd(i)*B(1,3)
         kky = k1_dd(i)*B(2,1) + k2_dd(i)*B(2,2) + k3_dd(i)*B(2,3)
         kkz = k1_dd(i)*B(3,1) + k2_dd(i)*B(3,2) + k3_dd(i)*B(3,3)
         r = (kx-kkx)**2 + (ky-kky)**2 + (kz-kkz)**2
         !write(unit,*) "diagonal distance:",i,r,rmax
         if (r.lt.rmax) then
           rmax = r
           id = i
         end if
      end do

      !write(unit,*) "diagonal d1,d2 =",d1(id),d2(id)

*     **** define six tetradrons - clunky but don't know defn of cOR in fortran ****
      if (id.eq.1) then
        itetra(1,1) = 0   +1
        itetra(2,1) = 7   +1
        itetra(3,1) = 1   +1
        itetra(4,1) = 3   +1

        itetra(1,2) = 0   +1
        itetra(2,2) = 7   +1
        itetra(3,2) = 1   +1
        itetra(4,2) = 5   +1

        itetra(1,3) = 0   +1
        itetra(2,3) = 7   +1
        itetra(3,3) = 2   +1
        itetra(4,3) = 3   +1

        itetra(1,4) = 0   +1
        itetra(2,4) = 7   +1
        itetra(3,4) = 2   +1
        itetra(4,4) = 6   +1

        itetra(1,5) = 0   +1
        itetra(2,5) = 7   +1
        itetra(3,5) = 4   +1
        itetra(4,5) = 5   +1

        itetra(1,6) = 0   +1
        itetra(2,6) = 7   +1
        itetra(3,6) = 4   +1
        itetra(4,6) = 6   +1
      else if (id.eq.2) then
        itetra(1,1) = 1   +1
        itetra(2,1) = 6   +1
        itetra(3,1) = 0   +1
        itetra(4,1) = 2   +1

        itetra(1,2) = 1   +1
        itetra(2,2) = 6   +1
        itetra(3,2) = 0   +1
        itetra(4,2) = 4   +1

        itetra(1,3) = 1   +1
        itetra(2,3) = 6   +1
        itetra(3,3) = 3   +1
        itetra(4,3) = 2   +1

        itetra(1,4) = 1   +1
        itetra(2,4) = 6   +1
        itetra(3,4) = 3   +1
        itetra(4,4) = 7   +1

        itetra(1,5) = 1   +1
        itetra(2,5) = 6   +1
        itetra(3,5) = 5   +1
        itetra(4,5) = 4   +1

        itetra(1,6) = 1   +1
        itetra(2,6) = 6   +1
        itetra(3,6) = 5   +1
        itetra(4,6) = 7   +1
      else if (id.eq.3) then
        itetra(1,1) = 2   +1
        itetra(2,1) = 5   +1
        itetra(3,1) = 3   +1
        itetra(4,1) = 1   +1

        itetra(1,2) = 2   +1
        itetra(2,2) = 5   +1
        itetra(3,2) = 3   +1
        itetra(4,2) = 7   +1

        itetra(1,3) = 2   +1
        itetra(2,3) = 5   +1
        itetra(3,3) = 0   +1
        itetra(4,3) = 1   +1

        itetra(1,4) = 2   +1
        itetra(2,4) = 5   +1
        itetra(3,4) = 0   +1
        itetra(4,4) = 4   +1

        itetra(1,5) = 2   +1
        itetra(2,5) = 5   +1
        itetra(3,5) = 6   +1
        itetra(4,5) = 7   +1

        itetra(1,6) = 2   +1
        itetra(2,6) = 5   +1
        itetra(3,6) = 6   +1
        itetra(4,6) = 4   +1
      else if (id.eq.4) then
        itetra(1,1) = 3   +1
        itetra(2,1) = 4   +1
        itetra(3,1) = 2    +1
        itetra(4,1) = 0   +1

        itetra(1,2) = 3   +1
        itetra(2,2) = 4   +1
        itetra(3,2) = 2   +1
        itetra(4,2) = 6   +1

        itetra(1,3) = 3   +1
        itetra(2,3) = 4   +1
        itetra(3,3) = 1   +1
        itetra(4,3) = 0   +1

        itetra(1,4) = 3   +1
        itetra(2,4) = 4   +1
        itetra(3,4) = 1   +1
        itetra(4,4) = 5   +1

        itetra(1,5) = 3   +1
        itetra(2,5) = 4   +1
        itetra(3,5) = 7   +1
        itetra(4,5) = 6   +1

        itetra(1,6) = 3   +1
        itetra(2,6) = 4   +1
        itetra(3,6) = 7   +1
        itetra(4,6) = 5   +1
      end if


c      do i=1,6
c        write(unit,*) id,"tetra :",i,"(",(itetra(j,i),j=1,4),")"
c      end do

      de = (emax-emin)/dble(npoints-1)
      do k=1,npoints
        e = emin + (k-1)*de

        f = 0.0d0
        g = 0.0d0
        do kk=1,dosgrid(3)
        do jj=1,dosgrid(2)
        do ii=1,dosgrid(1)
          ishft = ii+1
          jshft = jj+1
          kshft = kk+1
          if (ishft.gt.dosgrid(1)) ishft=1
          if (jshft.gt.dosgrid(2)) jshft=1
          if (kshft.gt.dosgrid(3)) kshft=1
          do i=1,neigs
            ecube(1) = eigs(ii,       jj,    kk, i)  ! (000)
            ecube(2) = eigs(ishft,    jj,    kk, i)  ! (001)
            ecube(3) = eigs(ii,    jshft,    kk, i)  ! (010)
            ecube(4) = eigs(ishft, jshft,    kk, i)  ! (011)
            ecube(5) = eigs(ii,       jj, kshft, i)  ! (100)
            ecube(6) = eigs(ishft,    jj, kshft, i)  ! (101)
            ecube(7) = eigs(   ii, jshft, kshft, i)  ! (110)
            ecube(8) = eigs(ishft, jshft, kshft, i)  ! (111)
           
            f = f + Dstates_Cube(e,itetra,ecube)
            g = g + Nstates_Cube(e,itetra,ecube)
          end do
        end do
        end do
        end do
        f = f*(VT/VG)
        g = g*(VT/VG)

        write(unit,1310) e,f*sign,g*sign
      end do


      return
 1310 FORMAT(3E15.5)
 3508 FORMAT(/' Brillouin zone point: ',i5,
     >       /'    k     =<',3f8.3,'> . <b1,b2,b3> ',
     >       /'          =<',3f8.3,'>')
      end


      real*8 function Dstates_Cube(e,itetra,ecube)
      implicit none
      real*8  e
      integer itetra(4,6)
      real*8  ecube(8)

*     **** local variables ****
      integer i,j,k
      real*8 ds,etetra(4),swap

      real*8   Dstates_Tetra
      external Dstates_Tetra

*     **** sum over 6 tetrahedrons ****
      ds = 0.0d0
      do k=1,6
        etetra(1) = ecube(itetra(1,k))
        etetra(2) = ecube(itetra(2,k))
        etetra(3) = ecube(itetra(3,k))
        etetra(4) = ecube(itetra(4,k))

*       **** bubble sort ****
        do j=1,3
        do i=j+1,4
          if (etetra(j).gt.etetra(i)) then
            swap      = etetra(i)
            etetra(i) = etetra(j)
            etetra(j) = swap
          end if
        end do
        end do
        ds = ds + Dstates_Tetra(e,etetra)
      end do

      Dstates_cube = ds
      return
      end


      real*8 function Dstates_Tetra(e,ee)
      implicit none
      real*8 e
      real*8 ee(4)

*     **** local variables ****
      real*8 ds
      real*8 e1,e2,e4
      real*8 e21,e31,e41,e32,e42,e43

      if ((ee(1).le.e).and.(e.lt.ee(2))) then
        e1 = e-ee(1)
        e21 = ee(2) - ee(1)
        e31 = ee(3) - ee(1)
        e41 = ee(4) - ee(1)
        ds = 3.0d0*e1*e1/(e21*e31*e41)
      else if ((ee(2).le.e).and.(e.lt.ee(3))) then
        e2 = e-ee(2) 
        e21 = ee(2) - ee(1) 
        e31 = ee(3) - ee(1) 
        e41 = ee(4) - ee(1) 
        e32 = ee(3) - ee(2) 
        e42 = ee(4) - ee(2) 
        ds = (3.0d0*e21+6.0d0*e2-3.0d0*(e31+e42)*e2*e2/(e32*e42))
     >       /(e31*e41)
      else if ((ee(3).le.e).and.(e.lt.ee(4))) then
        e4 = ee(4)-e 
        e41 = ee(4) - ee(1) 
        e42 = ee(4) - ee(2) 
        e43 = ee(4) - ee(3) 
        ds = 3.0d0*e4*e4/(e41*e42*e43)
      else
        ds = 0.0d0
      end if


      Dstates_Tetra = ds
      return
      end
      

      real*8 function Nstates_Cube(e,itetra,ecube)
      implicit none
      real*8  e
      integer itetra(4,6)
      real*8  ecube(8)

*     **** local variables ****
      integer i,j,k
      real*8 ds,etetra(4),swap

      real*8   Nstates_Tetra
      external Nstates_Tetra

*     **** sum over 6 tetrahedrons ****
      ds = 0.0d0
      do k=1,6
        etetra(1) = ecube(itetra(1,k))
        etetra(2) = ecube(itetra(2,k))
        etetra(3) = ecube(itetra(3,k))
        etetra(4) = ecube(itetra(4,k))

*       **** bubble sort ****
        do j=1,3
        do i=j+1,4
          if (etetra(j).gt.etetra(i)) then
            swap      = etetra(i)
            etetra(i) = etetra(j)
            etetra(j) = swap
          end if
        end do
        end do
        ds = ds + Nstates_Tetra(e,etetra)
      end do

      Nstates_cube = ds
      return
      end


      real*8 function Nstates_Tetra(e,ee)
      implicit none
      real*8 e
      real*8 ee(4)

*     **** local variables ****
      real*8 ds
      real*8 e1,e2,e4
      real*8 e21,e31,e41,e32,e42,e43

      if ((ee(1).le.e).and.(e.lt.ee(2))) then
        e1 = e-ee(1)
        e21 = ee(2) - ee(1)
        e31 = ee(3) - ee(1)
        e41 = ee(4) - ee(1)
        ds = e1*e1*e1/(e21*e31*e41)
      else if ((ee(2).le.e).and.(e.lt.ee(3))) then
        e2 = e-ee(2) 
        e21 = ee(2) - ee(1) 
        e31 = ee(3) - ee(1) 
        e41 = ee(4) - ee(1) 
        e32 = ee(3) - ee(2) 
        e42 = ee(4) - ee(2) 
        ds = (e21*e21 
     >        + 3.0d0*e21*e2
     >        + 3.0d0*e2*e2
     >        - (e31+e42)*e2*e2*e2/(e32*e42))
     >       /(e31*e41)
      else if ((ee(3).le.e).and.(e.lt.ee(4))) then
        e4 = ee(4)-e 
        e41 = ee(4) - ee(1) 
        e42 = ee(4) - ee(2) 
        e43 = ee(4) - ee(3) 
        ds = 1.0d0 - e4*e4*e4/(e41*e42*e43)
      else if (e.ge.ee(4)) then
        ds = 1.0d0
      else
        ds = 0.0d0
      end if


      Nstates_Tetra = ds
      return
      end



*     *********************************************
*     *                                           *
*     *            band_dos_weights_generate      *
*     *                                           *
*     *********************************************

      subroutine band_dos_weights_generate(idx,idy,idz,
     >                                     eigs,neigs,
     >                                     ein,weight)
      implicit none
      integer idx,idy,idz
      real*8 eigs(idx,idy,idz,*)
      integer neigs
      real*8 ein
      real*8 weight(neigs,idx,idy,idz)

*     **** local variables ****
      integer dosgrid(3)
      integer i,j,k,ii,jj,kk,ncubes,ntetra,count
      integer ishft,jshft,kshft
      integer k1_d(4),k2_d(4),k3_d(4),k1_dd(4),k2_dd(4),k3_dd(4)
      integer id,d1(4),d2(4)
      integer itetra(4,6)
      real*8  VT,VG
      real*8  B(3,3),unitg(3,3),e,ecube(8),wcube(8)
      real*8  k1,k2,k3,kx,ky,kz,kkx,kky,kkz,r,rmax

*     **** external functions ****
      real*8   lattice_unitg
      external lattice_unitg

      dosgrid(1) = idx
      dosgrid(2) = idy
      dosgrid(3) = idz


      do j=1,3
      do i=1,3
        B(i,j) = lattice_unitg(i,j)
      end do  
      end do  

*     **** volume of reciprocal unit cell, VG ****
      unitg(1,1) = B(2,2)*B(3,3) - B(3,2)*B(2,3)
      unitg(2,1) = B(3,2)*B(1,3) - B(1,2)*B(3,3)
      unitg(3,1) = B(1,2)*B(2,3) - B(2,2)*B(1,3)

      unitg(1,2) = B(2,3)*B(3,1) - B(3,3)*B(2,1)
      unitg(2,2) = B(3,3)*B(1,1) - B(1,3)*B(3,1)
      unitg(3,2) = B(1,3)*B(2,1) - B(2,3)*B(1,1)

      unitg(1,3) = B(2,1)*B(3,2) - B(3,1)*B(2,2)
      unitg(2,3) = B(3,1)*B(1,2) - B(1,1)*B(3,2)
      unitg(3,3) = B(1,1)*B(2,2) - B(2,1)*B(1,2)
      VG = B(1,1)*unitg(1,1)
     >   + B(2,1)*unitg(2,1)
     >   + B(3,1)*unitg(3,1)
      
      ncubes = dosgrid(1)*dosgrid(2)*dosgrid(3)
      ntetra = ncubes*6
      VT = VG/dble(ntetra)


*     ********************************
*     **** find shortest diagonal ****
*     ********************************

*     **** (000) ---- (111) ****
      k1_d(1) = 0
      k2_d(1) = 0
      k3_d(1) = 0
      k1_dd(1) = 1
      k2_dd(1) = 1
      k3_dd(1) = 1
      d1(1) = 0
      d2(1) = 7

*     **** (001) ---- (110) ****
      k1_d(2) = 1
      k2_d(2) = 0
      k3_d(2) = 0
      k1_dd(2) = 0
      k2_dd(2) = 1
      k3_dd(2) = 1
      d1(2) = 1
      d2(2) = 6

*     **** (010) ---- (101) ****
      k1_d(3) = 0
      k2_d(3) = 1
      k3_d(3) = 0
      k1_dd(3) = 1
      k2_dd(3) = 0
      k3_dd(3) = 1
      d1(3) = 2
      d2(3) = 5

*     **** (011) ---- (100) ****
      k1_d(4) = 1
      k2_d(4) = 1
      k3_d(4) = 0
      k1_dd(4) = 0
      k2_dd(4) = 0
      k3_dd(4) = 1
      d1(4) = 3
      d2(4) = 4

      id = 1
      rmax = 9.99d9
      do i=1,4
         kx = k1_d(i)*B(1,1) + k2_d(i)*B(1,2) + k3_d(i)*B(1,3)
         ky = k1_d(i)*B(2,1) + k2_d(i)*B(2,2) + k3_d(i)*B(2,3)
         kz = k1_d(i)*B(3,1) + k2_d(i)*B(3,2) + k3_d(i)*B(3,3)

         kkx = k1_dd(i)*B(1,1) + k2_dd(i)*B(1,2) + k3_dd(i)*B(1,3)
         kky = k1_dd(i)*B(2,1) + k2_dd(i)*B(2,2) + k3_dd(i)*B(2,3)
         kkz = k1_dd(i)*B(3,1) + k2_dd(i)*B(3,2) + k3_dd(i)*B(3,3)
         r = (kx-kkx)**2 + (ky-kky)**2 + (kz-kkz)**2
         !write(unit,*) "diagonal distance:",i,r,rmax
         if (r.lt.rmax) then
           rmax = r
           id = i
         end if
      end do


*     **** define six tetradrons - clunky but don't know defn of cOR in fortran ****
      if (id.eq.1) then
        itetra(1,1) = 0   +1
        itetra(2,1) = 7   +1
        itetra(3,1) = 1   +1
        itetra(4,1) = 3   +1

        itetra(1,2) = 0   +1
        itetra(2,2) = 7   +1
        itetra(3,2) = 1   +1
        itetra(4,2) = 5   +1

        itetra(1,3) = 0   +1
        itetra(2,3) = 7   +1
        itetra(3,3) = 2   +1
        itetra(4,3) = 3   +1

        itetra(1,4) = 0   +1
        itetra(2,4) = 7   +1
        itetra(3,4) = 2   +1
        itetra(4,4) = 6   +1

        itetra(1,5) = 0   +1
        itetra(2,5) = 7   +1
        itetra(3,5) = 4   +1
        itetra(4,5) = 5   +1

        itetra(1,6) = 0   +1
        itetra(2,6) = 7   +1
        itetra(3,6) = 4   +1
        itetra(4,6) = 6   +1
      else if (id.eq.2) then
        itetra(1,1) = 1   +1
        itetra(2,1) = 6   +1
        itetra(3,1) = 0   +1
        itetra(4,1) = 2   +1

        itetra(1,2) = 1   +1
        itetra(2,2) = 6   +1
        itetra(3,2) = 0   +1
        itetra(4,2) = 4   +1

        itetra(1,3) = 1   +1
        itetra(2,3) = 6   +1
        itetra(3,3) = 3   +1
        itetra(4,3) = 2   +1

        itetra(1,4) = 1   +1
        itetra(2,4) = 6   +1
        itetra(3,4) = 3   +1
        itetra(4,4) = 7   +1

        itetra(1,5) = 1   +1
        itetra(2,5) = 6   +1
        itetra(3,5) = 5   +1
        itetra(4,5) = 4   +1

        itetra(1,6) = 1   +1
        itetra(2,6) = 6   +1
        itetra(3,6) = 5   +1
        itetra(4,6) = 7   +1
      else if (id.eq.3) then
        itetra(1,1) = 2   +1
        itetra(2,1) = 5   +1
        itetra(3,1) = 3   +1
        itetra(4,1) = 1   +1

        itetra(1,2) = 2   +1
        itetra(2,2) = 5   +1
        itetra(3,2) = 3   +1
        itetra(4,2) = 7   +1

        itetra(1,3) = 2   +1
        itetra(2,3) = 5   +1
        itetra(3,3) = 0   +1
        itetra(4,3) = 1   +1

        itetra(1,4) = 2   +1
        itetra(2,4) = 5   +1
        itetra(3,4) = 0   +1
        itetra(4,4) = 4   +1

        itetra(1,5) = 2   +1
        itetra(2,5) = 5   +1
        itetra(3,5) = 6   +1
        itetra(4,5) = 7   +1

        itetra(1,6) = 2   +1
        itetra(2,6) = 5   +1
        itetra(3,6) = 6   +1
        itetra(4,6) = 4   +1
      else if (id.eq.4) then
        itetra(1,1) = 3   +1
        itetra(2,1) = 4   +1
        itetra(3,1) = 2    +1
        itetra(4,1) = 0   +1

        itetra(1,2) = 3   +1
        itetra(2,2) = 4   +1
        itetra(3,2) = 2   +1
        itetra(4,2) = 6   +1

        itetra(1,3) = 3   +1
        itetra(2,3) = 4   +1
        itetra(3,3) = 1   +1
        itetra(4,3) = 0   +1

        itetra(1,4) = 3   +1
        itetra(2,4) = 4   +1
        itetra(3,4) = 1   +1
        itetra(4,4) = 5   +1

        itetra(1,5) = 3   +1
        itetra(2,5) = 4   +1
        itetra(3,5) = 7   +1
        itetra(4,5) = 6   +1

        itetra(1,6) = 3   +1
        itetra(2,6) = 4   +1
        itetra(3,6) = 7   +1
        itetra(4,6) = 5   +1
      end if



      
      e = ein

      call dcopy(dosgrid(1)*dosgrid(2)*dosgrid(3)*neigs,
     >           0.0d0,0,weight,1)
      do kk=1,dosgrid(3)
      do jj=1,dosgrid(2)
      do ii=1,dosgrid(1)
        ishft = ii+1
        jshft = jj+1
        kshft = kk+1
        if (ishft.gt.dosgrid(1)) ishft=1
        if (jshft.gt.dosgrid(2)) jshft=1
        if (kshft.gt.dosgrid(3)) kshft=1
        do i=1,neigs
          ecube(1) = eigs(ii,       jj,    kk, i)  ! (000)
          ecube(2) = eigs(ishft,    jj,    kk, i)  ! (001)
          ecube(3) = eigs(ii,    jshft,    kk, i)  ! (010)
          ecube(4) = eigs(ishft, jshft,    kk, i)  ! (011)
          ecube(5) = eigs(ii,       jj, kshft, i)  ! (100)
          ecube(6) = eigs(ishft,    jj, kshft, i)  ! (101)
          ecube(7) = eigs(   ii, jshft, kshft, i)  ! (110)
          ecube(8) = eigs(ishft, jshft, kshft, i)  ! (111)
         
          call Dstates_weight_Cube(e,itetra,ecube,wcube)
          weight(i,ii,      jj,   kk) = weight(i,ii,      jj,   kk) 
     >                                + wcube(1)
          weight(i,ishft,   jj,   kk) = weight(i,ishft,   jj,   kk)
     >                                + wcube(2)
          weight(i,ii,   jshft,   kk) = weight(i,ii,   jshft,   kk) 
     >                                + wcube(3)
          weight(i,ishft,jshft,   kk) = weight(i,ishft,jshft,   kk)
     >                                + wcube(4)
          weight(i,ii,      jj,kshft) = weight(i,ii,      jj,kshft)
     >                                + wcube(5)
          weight(i,ishft,   jj,kshft) = weight(i,ishft,   jj,kshft)
     >                                + wcube(6)
          weight(i,   ii,jshft,kshft) = weight(i,   ii,jshft,kshft)
     >                                + wcube(7)
          weight(i,ishft,jshft,kshft) = weight(i,ishft,jshft,kshft)
     >                                + wcube(8)
        end do
      end do
      end do
      end do
      call dscal(dosgrid(1)*dosgrid(2)*dosgrid(3)*neigs,
     >          (VT/VG),weight,1)


      return
      end


      subroutine Dstates_weight_Cube(e,itetra,ecube,wcube)
      implicit none
      real*8  e
      integer itetra(4,6)
      real*8  ecube(8)
      real*8  wcube(8)

*     **** local variables ****
      integer i,j,k
      real*8 ds,etetra(4),swap

*     **** external functions ****
      real*8   Dstates_Tetra
      external Dstates_Tetra


*     **** sum over 6 tetrahedrons ****
      call dcopy(8,0.0d0,0,wcube,1)
      do k=1,6
        etetra(1) = ecube(itetra(1,k))
        etetra(2) = ecube(itetra(2,k))
        etetra(3) = ecube(itetra(3,k))
        etetra(4) = ecube(itetra(4,k))

*       **** bubble sort ****
        do j=1,3
        do i=j+1,4
          if (etetra(j).gt.etetra(i)) then
            swap      = etetra(i)
            etetra(i) = etetra(j)
            etetra(j) = swap
          end if
        end do
        end do

        ds = Dstates_Tetra(e,etetra)

        wcube(itetra(1,k)) = wcube(itetra(1,k)) + 0.25d0*ds
        wcube(itetra(2,k)) = wcube(itetra(2,k)) + 0.25d0*ds
        wcube(itetra(3,k)) = wcube(itetra(3,k)) + 0.25d0*ds
        wcube(itetra(4,k)) = wcube(itetra(4,k)) + 0.25d0*ds
      end do

      return
      end



*     *********************************************
*     *                                           *
*     *            band_projected_dos_generate    *
*     *                                           *
*     *********************************************

      subroutine band_projected_dos_generate(
     >                             unit,idx,idy,idz,eigs,pweight,neigs,
     >                             sign,npoints,emin,emax)
      implicit none
      integer unit
      integer idx,idy,idz
      real*8 eigs(idx,idy,idz,*)
      real*8 pweight(idx,idy,idz,*)
      integer neigs
      real*8  sign
      integer npoints
      real*8 emin,emax

*     **** local variables ****
      integer dosgrid(3)
      integer i,j,k,ii,jj,kk,ncubes,ntetra,count
      integer ishft,jshft,kshft
      integer k1_d(4),k2_d(4),k3_d(4),k1_dd(4),k2_dd(4),k3_dd(4)
      integer id,d1(4),d2(4)
      integer itetra(4,6)
      real*8  VT,VG
      real*8  B(3,3),unitg(3,3),e,ecube(8),pcube(8),f,g,de,ff
      real*8  k1,k2,k3,kx,ky,kz,kkx,kky,kkz,r,rmax

*     **** external functions ****
      real*8   lattice_unitg
      real*8   Dstates_Cube_projected,Nstates_Cube_projected
      external lattice_unitg
      external Dstates_Cube_projected,Nstates_Cube_projected

      dosgrid(1) = idx
      dosgrid(2) = idy
      dosgrid(3) = idz

      do j=1,3
      do i=1,3
        B(i,j) = lattice_unitg(i,j)
      end do  
      end do  

*     **** volume of reciprocal unit cell, VG ****
      unitg(1,1) = B(2,2)*B(3,3) - B(3,2)*B(2,3)
      unitg(2,1) = B(3,2)*B(1,3) - B(1,2)*B(3,3)
      unitg(3,1) = B(1,2)*B(2,3) - B(2,2)*B(1,3)

      unitg(1,2) = B(2,3)*B(3,1) - B(3,3)*B(2,1)
      unitg(2,2) = B(3,3)*B(1,1) - B(1,3)*B(3,1)
      unitg(3,2) = B(1,3)*B(2,1) - B(2,3)*B(1,1)

      unitg(1,3) = B(2,1)*B(3,2) - B(3,1)*B(2,2)
      unitg(2,3) = B(3,1)*B(1,2) - B(1,1)*B(3,2)
      unitg(3,3) = B(1,1)*B(2,2) - B(2,1)*B(1,2)
      VG = B(1,1)*unitg(1,1)
     >   + B(2,1)*unitg(2,1)
     >   + B(3,1)*unitg(3,1)
      
      ncubes = dosgrid(1)*dosgrid(2)*dosgrid(3)
      ntetra = ncubes*6
      VT = VG/dble(ntetra)


*     ********************************
*     **** find shortest diagonal ****
*     ********************************

*     **** (000) ---- (111) ****
      k1_d(1) = 0
      k2_d(1) = 0
      k3_d(1) = 0
      k1_dd(1) = 1
      k2_dd(1) = 1
      k3_dd(1) = 1
      d1(1) = 0
      d2(1) = 7

*     **** (001) ---- (110) ****
      k1_d(2) = 1
      k2_d(2) = 0
      k3_d(2) = 0
      k1_dd(2) = 0
      k2_dd(2) = 1
      k3_dd(2) = 1
      d1(2) = 1
      d2(2) = 6

*     **** (010) ---- (101) ****
      k1_d(3) = 0
      k2_d(3) = 1
      k3_d(3) = 0
      k1_dd(3) = 1
      k2_dd(3) = 0
      k3_dd(3) = 1
      d1(3) = 2
      d2(3) = 5

*     **** (011) ---- (100) ****
      k1_d(4) = 1
      k2_d(4) = 1
      k3_d(4) = 0
      k1_dd(4) = 0
      k2_dd(4) = 0
      k3_dd(4) = 1
      d1(4) = 3
      d2(4) = 4

      id = 1
      rmax = 9.99d9
      do i=1,4
         kx = k1_d(i)*B(1,1) + k2_d(i)*B(1,2) + k3_d(i)*B(1,3)
         ky = k1_d(i)*B(2,1) + k2_d(i)*B(2,2) + k3_d(i)*B(2,3)
         kz = k1_d(i)*B(3,1) + k2_d(i)*B(3,2) + k3_d(i)*B(3,3)

         kkx = k1_dd(i)*B(1,1) + k2_dd(i)*B(1,2) + k3_dd(i)*B(1,3)
         kky = k1_dd(i)*B(2,1) + k2_dd(i)*B(2,2) + k3_dd(i)*B(2,3)
         kkz = k1_dd(i)*B(3,1) + k2_dd(i)*B(3,2) + k3_dd(i)*B(3,3)
         r = (kx-kkx)**2 + (ky-kky)**2 + (kz-kkz)**2
         !write(unit,*) "diagonal distance:",i,r,rmax
         if (r.lt.rmax) then
           rmax = r
           id = i
         end if
      end do

      !write(unit,*) "diagonal d1,d2 =",d1(id),d2(id)

*     **** define six tetradrons - clunky but don't know defn of cOR in fortran ****
      if (id.eq.1) then
        itetra(1,1) = 0   +1
        itetra(2,1) = 7   +1
        itetra(3,1) = 1   +1
        itetra(4,1) = 3   +1

        itetra(1,2) = 0   +1
        itetra(2,2) = 7   +1
        itetra(3,2) = 1   +1
        itetra(4,2) = 5   +1

        itetra(1,3) = 0   +1
        itetra(2,3) = 7   +1
        itetra(3,3) = 2   +1
        itetra(4,3) = 3   +1

        itetra(1,4) = 0   +1
        itetra(2,4) = 7   +1
        itetra(3,4) = 2   +1
        itetra(4,4) = 6   +1

        itetra(1,5) = 0   +1
        itetra(2,5) = 7   +1
        itetra(3,5) = 4   +1
        itetra(4,5) = 5   +1

        itetra(1,6) = 0   +1
        itetra(2,6) = 7   +1
        itetra(3,6) = 4   +1
        itetra(4,6) = 6   +1
      else if (id.eq.2) then
        itetra(1,1) = 1   +1
        itetra(2,1) = 6   +1
        itetra(3,1) = 0   +1
        itetra(4,1) = 2   +1

        itetra(1,2) = 1   +1
        itetra(2,2) = 6   +1
        itetra(3,2) = 0   +1
        itetra(4,2) = 4   +1

        itetra(1,3) = 1   +1
        itetra(2,3) = 6   +1
        itetra(3,3) = 3   +1
        itetra(4,3) = 2   +1

        itetra(1,4) = 1   +1
        itetra(2,4) = 6   +1
        itetra(3,4) = 3   +1
        itetra(4,4) = 7   +1

        itetra(1,5) = 1   +1
        itetra(2,5) = 6   +1
        itetra(3,5) = 5   +1
        itetra(4,5) = 4   +1

        itetra(1,6) = 1   +1
        itetra(2,6) = 6   +1
        itetra(3,6) = 5   +1
        itetra(4,6) = 7   +1
      else if (id.eq.3) then
        itetra(1,1) = 2   +1
        itetra(2,1) = 5   +1
        itetra(3,1) = 3   +1
        itetra(4,1) = 1   +1

        itetra(1,2) = 2   +1
        itetra(2,2) = 5   +1
        itetra(3,2) = 3   +1
        itetra(4,2) = 7   +1

        itetra(1,3) = 2   +1
        itetra(2,3) = 5   +1
        itetra(3,3) = 0   +1
        itetra(4,3) = 1   +1

        itetra(1,4) = 2   +1
        itetra(2,4) = 5   +1
        itetra(3,4) = 0   +1
        itetra(4,4) = 4   +1

        itetra(1,5) = 2   +1
        itetra(2,5) = 5   +1
        itetra(3,5) = 6   +1
        itetra(4,5) = 7   +1

        itetra(1,6) = 2   +1
        itetra(2,6) = 5   +1
        itetra(3,6) = 6   +1
        itetra(4,6) = 4   +1
      else if (id.eq.4) then
        itetra(1,1) = 3   +1
        itetra(2,1) = 4   +1
        itetra(3,1) = 2    +1
        itetra(4,1) = 0   +1

        itetra(1,2) = 3   +1
        itetra(2,2) = 4   +1
        itetra(3,2) = 2   +1
        itetra(4,2) = 6   +1

        itetra(1,3) = 3   +1
        itetra(2,3) = 4   +1
        itetra(3,3) = 1   +1
        itetra(4,3) = 0   +1

        itetra(1,4) = 3   +1
        itetra(2,4) = 4   +1
        itetra(3,4) = 1   +1
        itetra(4,4) = 5   +1

        itetra(1,5) = 3   +1
        itetra(2,5) = 4   +1
        itetra(3,5) = 7   +1
        itetra(4,5) = 6   +1

        itetra(1,6) = 3   +1
        itetra(2,6) = 4   +1
        itetra(3,6) = 7   +1
        itetra(4,6) = 5   +1
      end if


c      do i=1,6
c        write(unit,*) id,"tetra :",i,"(",(itetra(j,i),j=1,4),")"
c      end do

      ff = 0.0d0
      de = (emax-emin)/dble(npoints-1)
      do k=1,npoints
        e = emin + (k-1)*de

        f = 0.0d0
        g = 0.0d0
        do kk=1,dosgrid(3)
        do jj=1,dosgrid(2)
        do ii=1,dosgrid(1)
          ishft = ii+1
          jshft = jj+1
          kshft = kk+1
          if (ishft.gt.dosgrid(1)) ishft=1
          if (jshft.gt.dosgrid(2)) jshft=1
          if (kshft.gt.dosgrid(3)) kshft=1
          do i=1,neigs
            ecube(1) = eigs(ii,       jj,    kk, i)  ! (000)
            ecube(2) = eigs(ishft,    jj,    kk, i)  ! (001)
            ecube(3) = eigs(ii,    jshft,    kk, i)  ! (010)
            ecube(4) = eigs(ishft, jshft,    kk, i)  ! (011)
            ecube(5) = eigs(ii,       jj, kshft, i)  ! (100)
            ecube(6) = eigs(ishft,    jj, kshft, i)  ! (101)
            ecube(7) = eigs(   ii, jshft, kshft, i)  ! (110)
            ecube(8) = eigs(ishft, jshft, kshft, i)  ! (111)

            pcube(1) = pweight(ii,       jj,    kk, i)  ! (000)
            pcube(2) = pweight(ishft,    jj,    kk, i)  ! (001)
            pcube(3) = pweight(ii,    jshft,    kk, i)  ! (010)
            pcube(4) = pweight(ishft, jshft,    kk, i)  ! (011)
            pcube(5) = pweight(ii,       jj, kshft, i)  ! (100)
            pcube(6) = pweight(ishft,    jj, kshft, i)  ! (101)
            pcube(7) = pweight(   ii, jshft, kshft, i)  ! (110)
            pcube(8) = pweight(ishft, jshft, kshft, i)  ! (111)
     
           
            f = f + Dstates_Cube_projected(e,itetra,ecube,pcube)
            !g = g + Nstates_Cube_projected(e,itetra,ecube,pcube)
          end do
        end do
        end do
        end do
        f = f*(VT/VG)
        !g = g*(VT/VG)
        if ((k.eq.1).or.(k.eq.npoints)) then
           ff = ff + 0.5d0*f*de
        else
           ff = ff + f*de
        end if

        !write(unit,1310) e,f*sign,g
        write(unit,1310) e,f*sign,ff*sign
      end do


      return
 1310 FORMAT(3E15.5)
 3508 FORMAT(/' Brillouin zone point: ',i5,
     >       /'    k     =<',3f8.3,'> . <b1,b2,b3> ',
     >       /'          =<',3f8.3,'>')
      end


      real*8 function Dstates_Cube_projected(e,itetra,ecube,pcube)
      implicit none
      real*8  e
      integer itetra(4,6)
      real*8  ecube(8),pcube(8)

*     **** local variables ****
      integer i,j,k
      real*8 ds,etetra(4),ptetra(4),swap

      real*8   Dstates_Tetra_projected
      external Dstates_Tetra_projected

*     **** sum over 6 tetrahedrons ****
      ds = 0.0d0
      do k=1,6
        etetra(1) = ecube(itetra(1,k))
        etetra(2) = ecube(itetra(2,k))
        etetra(3) = ecube(itetra(3,k))
        etetra(4) = ecube(itetra(4,k))

        ptetra(1) = pcube(itetra(1,k))
        ptetra(2) = pcube(itetra(2,k))
        ptetra(3) = pcube(itetra(3,k))
        ptetra(4) = pcube(itetra(4,k))

*       **** bubble sort ****
        do j=1,3
        do i=j+1,4
          if (etetra(j).gt.etetra(i)) then
            swap      = etetra(i)
            etetra(i) = etetra(j)
            etetra(j) = swap
            swap      = ptetra(i)
            ptetra(i) = ptetra(j)
            ptetra(j) = swap
          end if
        end do
        end do
        ds = ds + Dstates_Tetra_projected(e,etetra,ptetra)
      end do

      Dstates_cube_projected = ds
      return
      end


      real*8 function Dstates_Tetra_projected(e,ee,pp)
      implicit none
      real*8 e
      real*8 ee(4),pp(4)

*     **** local variables ****
      real*8 ds
      real*8 e1,e2,e4
      real*8 e21,e31,e41,e32,e42,e43
      real*8 points(3,3)

*     **** external functions ****
      real*8   Dstate_triangle
      external Dstate_triangle

      if ((ee(1).le.e).and.(e.lt.ee(2))) then
c        e1 = e-ee(1)
c        e21 = ee(2) - ee(1)
c        e31 = ee(3) - ee(1)
c        e41 = ee(4) - ee(1)
c        ds = 3.0d0*e1*e1/(e21*e31*e41)
        points(1,1) = (e-ee(1))/(ee(2)-ee(1))
        points(2,1) = 0.0d0
        points(3,1) = 0.0d0
        points(1,2) = 0.0d0
        points(2,2) = (e-ee(1))/(ee(3)-ee(1))
        points(3,2) = 0.0d0
        points(1,3) = 0.0d0 
        points(2,3) = 0.0d0
        points(3,3) = (e-ee(1))/(ee(4)-ee(1))
        ds = Dstate_triangle(points,ee,pp)

      else if ((ee(2).le.e).and.(e.lt.ee(3))) then
c        e2 = e-ee(2) 
c        e21 = ee(2) - ee(1) 
c        e31 = ee(3) - ee(1) 
c        e41 = ee(4) - ee(1) 
c        e32 = ee(3) - ee(2) 
c        e42 = ee(4) - ee(2) 
c        ds = (3.0d0*e21+6.0d0*e2-3.0d0*(e31+e42)*e2*e2/(e32*e42))
c     >       /(e31*e41)
        points(1,1) = 0.0d0
        points(2,1) = (e-ee(1))/(ee(3)-ee(1))
        points(3,1) = 0.0d0
        points(1,2) = 1.0d0 - (e-ee(2))/(ee(3)-ee(2))
        points(2,2) =         (e-ee(2))/(ee(3)-ee(2))
        points(3,2) = 0.0d0
        points(1,3) = 0.0d0
        points(2,3) = 0.0d0
        points(3,3) = (e-ee(1))/(ee(4)-ee(1))
        ds = Dstate_triangle(points,ee,pp)
        points(1,1) = 1.0d0 - (e-ee(2))/(ee(4)-ee(2))
        points(2,1) = 0.0d0
        points(3,1) =         (e-ee(2))/(ee(4)-ee(2))
        points(1,2) = 1.0d0 - (e-ee(2))/(ee(3)-ee(2))
        points(2,2) = (e-ee(2))/(ee(3)-ee(2))
        points(3,2) = 0.0d0
        points(1,3) = 0.0d0
        points(2,3) = 0.0d0
        points(3,3) = (e-ee(1))/(ee(4)-ee(1))
        ds = ds + Dstate_triangle(points,ee,pp)

      else if ((ee(3).le.e).and.(e.lt.ee(4))) then
c        e4 = ee(4)-e 
c        e41 = ee(4) - ee(1) 
c        e42 = ee(4) - ee(2) 
c        e43 = ee(4) - ee(3) 
c        ds = 3.0d0*e4*e4/(e41*e42*e43)
        points(1,1) = 1.0d0 - (e-ee(2))/(ee(4)-ee(2))
        points(2,1) = 0.0d0
        points(3,1) = (e-ee(2))/(ee(4)-ee(2))
        points(1,2) = 0.0d0
        points(2,2) = 1.0d0 - (e-ee(3))/(ee(4)-ee(3)) 
        points(3,2) = (e-ee(3))/(ee(4)-ee(3))
        points(1,3) = 0.0d0
        points(2,3) = 0.0d0
        points(3,3) = (e-ee(1))/(ee(4)-ee(1))
        ds = Dstate_triangle(points,ee,pp)

      else
        ds = 0.0d0
      end if


      Dstates_Tetra_projected = ds
      return
      end



*     ******************************************
*     *                                        *
*     *          Dstate_triangle               *
*     *                                        *
*     ******************************************
      real*8 function Dstate_triangle(points,ee,ff)
      implicit none
      real*8 points(3,3)
      real*8 ee(4),ff(4)

      real*8 p13x,p13y,p13z,p23x,p23y,p23z
      real*8 f0,f1,f2,e10,e20,e30,nde,n2
      real*8 n(3)

      p13x = points(1,1) - points(1,3)
      p13y = points(2,1) - points(2,3)
      p13z = points(3,1) - points(3,3)

      p23x = points(1,2) - points(1,3)
      p23y = points(2,2) - points(2,3)
      p23z = points(3,2) - points(3,3)

      n(1) = p13y*p23z - p13z*p23y
      n(2) = p13z*p23x - p13x*p23z
      n(3) = p13x*p23y - p13y*p23x
      n2 = n(1)*n(1) + n(2)*n(2) + n(3)*n(3)

      e10 = ee(2)-ee(1)
      e20 = ee(3)-ee(1)
      e30 = ee(4)-ee(1)
      nde = dabs(e10*n(1) + e20*n(2) + e30*n(3))

      f0  = ff(1) + (ff(2)-ff(1))*points(1,3)
     >            + (ff(3)-ff(1))*points(2,3)
     >            + (ff(4)-ff(1))*points(3,3)

      f1 = (ff(2)-ff(1))*(points(1,1)-points(1,3))
     >   + (ff(3)-ff(1))*(points(2,1)-points(2,3))
     >   + (ff(4)-ff(1))*(points(3,1)-points(3,3))

      f2 = (ff(2)-ff(1))*(points(1,2)-points(1,3))
     >   + (ff(3)-ff(1))*(points(2,2)-points(2,3))
     >   + (ff(4)-ff(1))*(points(3,2)-points(3,3))

      Dstate_triangle = 6.0d0*(n2/nde)*(f0/2.0d0 + f1/6.0d0 + f2/6.0d0)
      return
      end



      

      real*8 function Nstates_Cube_projected(e,itetra,ecube,pcube)
      implicit none
      real*8  e
      integer itetra(4,6)
      real*8  ecube(8),pcube(8)

*     **** local variables ****
      integer i,j,k
      real*8 ds,etetra(4),ptetra(4),swap

      real*8   Nstates_Tetra_projected
      external Nstates_Tetra_projected

*     **** sum over 6 tetrahedrons ****
      ds = 0.0d0
      do k=1,6
        etetra(1) = ecube(itetra(1,k))
        etetra(2) = ecube(itetra(2,k))
        etetra(3) = ecube(itetra(3,k))
        etetra(4) = ecube(itetra(4,k))
        ptetra(1) = pcube(itetra(1,k))
        ptetra(2) = pcube(itetra(2,k))
        ptetra(3) = pcube(itetra(3,k))
        ptetra(4) = pcube(itetra(4,k))

*       **** bubble sort ****
        do j=1,3
        do i=j+1,4
          if (etetra(j).gt.etetra(i)) then
            swap      = etetra(i)
            etetra(i) = etetra(j)
            etetra(j) = swap
            swap      = ptetra(i)
            ptetra(i) = ptetra(j)
            ptetra(j) = swap
          end if
        end do
        end do
        ds = ds + Nstates_Tetra_projected(e,etetra,ptetra)
      end do

      Nstates_cube_projected = ds
      return
      end


      real*8 function Nstates_Tetra_projected(e,ee,pp)
      implicit none
      real*8 e
      real*8 ee(4),pp(4)

*     **** local variables ****
      real*8 ds
      real*8 e1,e2,e4
      real*8 e21,e31,e41,e32,e42,e43

      if ((ee(1).le.e).and.(e.lt.ee(2))) then
        e1 = e-ee(1)
        e21 = ee(2) - ee(1)
        e31 = ee(3) - ee(1)
        e41 = ee(4) - ee(1)
        ds = e1*e1*e1/(e21*e31*e41)
      else if ((ee(2).le.e).and.(e.lt.ee(3))) then
        e2 = e-ee(2) 
        e21 = ee(2) - ee(1) 
        e31 = ee(3) - ee(1) 
        e41 = ee(4) - ee(1) 
        e32 = ee(3) - ee(2) 
        e42 = ee(4) - ee(2) 
        ds = (e21*e21 
     >        + 3.0d0*e21*e2
     >        + 3.0d0*e2*e2
     >        - (e31+e42)*e2*e2*e2/(e32*e42))
     >       /(e31*e41)
      else if ((ee(3).le.e).and.(e.lt.ee(4))) then
        e4 = ee(4)-e 
        e41 = ee(4) - ee(1) 
        e42 = ee(4) - ee(2) 
        e43 = ee(4) - ee(3) 
        ds = 1.0d0 - e4*e4*e4/(e41*e42*e43)
      else if (e.ge.ee(4)) then
        ds = 1.0d0
      else
        ds = 0.0d0
      end if


      Nstates_Tetra_projected = ds
      return
      end


