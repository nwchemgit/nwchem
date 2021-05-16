*
* $Id$
*

*     ****************************************************************
*     *                                                              *
*     *                  wvfnc_new_lcao_density                      *
*     *                                                              *
*     ****************************************************************
      logical function wvfnc_new_lcao_density(oprint_in,
     >                                        wavefunction_filename,
     >                                        version,ngrid,unita,
     >                                        ispin,ne,oddelcfill)

      implicit none
      logical     oprint_in
      character*50 wavefunction_filename
      integer      version
      integer      ngrid(3)
      real*8       unita(3,3)
      integer      ispin,ne(2)
      logical      oddelcfill

#include "bafdecls.fh"
#include "errquit.fh"

 
*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      real*8 ALPHA,talpha
      parameter (ALPHA=0.30d0)

      logical value,value2,oprint,field_exist
      integer ms,i,ia,icharge,ne_excited(2),nn,l,m,n,ii
      integer dn(2),dnall(2),phi1(2),rho(2),vc(2),xcp(2),xce(2)
      integer psi1(2),occ1(2),smearoccupation
      integer ee(2),eigs(2),neq(2),mapping1d
      integer dng(2),tmp1(2),vall(2),v_field(2),vl_lr(2),vl(2)
      integer nx,ny,nz
      integer nbasis,nbasis2,n2ft3d,npack0,npack1,it
      real*8  rho_error,rho_error_old,total_dn,sum
      real*8  Etotal,Eorbs,Ehart,Eexc,Eion,EV,de,olde,pxc
      real*8  cpu1,cpu2,cpu3,cpu4,cpu5,t1,t2,t3,t4,av
      real*8  scal1,scal2,dv
      character*50 filename

*     **** external functions ****
      logical     psp_semicore,aorbs_init,aorbs_readall
      character*4 ion_atom_qm,ion_atom
      integer     aorbs_norbs,aorbs_nbasis,ion_nkatm,psp_lmax,aorbs_katm
      integer     psp_locp,psp_lmmax,ion_natm,ion_nion,ion_nkatm_qm
      integer     control_version,aorbs_l,aorbs_m,aorbs_ii
      integer     ewald_ncut,ewald_nshl3d
      real*8      psp_ncore,psp_rcore,psp_rlocal,psp_rc,psp_zv
      real*8      ewald_rcut,ewald_e,ion_ion_e,coulomb_e
      real*8      ewald_mandelung
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg,control_tole,control_tolc
      external    psp_semicore,aorbs_init,aorbs_readall
      external    ion_atom_qm,ion_atom
      external    aorbs_norbs,aorbs_nbasis,ion_nkatm,psp_lmax,aorbs_katm
      external    psp_locp,psp_lmmax,ion_natm,ion_nion,ion_nkatm_qm
      external    control_version,aorbs_l,aorbs_m,aorbs_ii
      external    ewald_ncut,ewald_nshl3d
      external    psp_ncore,psp_rcore,psp_rlocal,psp_rc,psp_zv
      external    ewald_rcut,ewald_e,ion_ion_e,coulomb_e
      external    ewald_mandelung
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg,control_tole,control_tolc
      integer  pack_nwave_all
      integer  control_gga
      integer  control_ngrid,pack_nwave
      external pack_nwave_all
      external control_gga
      external control_ngrid,pack_nwave
      character*12 control_boundry
      external     control_boundry

      integer      control_excited_ne,control_mapping1d
      character*50 control_input_epsi
      external     control_excited_ne,control_mapping1d
      external     control_input_epsi
      logical      ion_chargeexist,ion_mmexist
      external     ion_chargeexist,ion_mmexist




      Etotal = 0.0d0
      oprint = oprint_in

      call Parallel_taskid(taskid)

      if ((taskid.eq.MASTER).and.oprint) call current_second(cpu1)

      call D3dB_n2ft3d(1,n2ft3d)
      call Pack_npack(0,npack0)
      call Pack_npack(1,npack1)

      field_exist = ion_chargeexist().or.ion_mmexist()

      mapping1d = control_mapping1d()
      call Dne_init(ispin,ne,mapping1d)
      call Dneall_neq(neq)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)

*     **** initialize atomic orbitals ****
      value = aorbs_init()
      value = value.and.aorbs_readall()
      if (.not.value) go to 101
      nbasis = aorbs_nbasis()
      nbasis2 = nbasis**2

*     **** basis set not big enough to meeting electronic filling ****
      if ((nbasis.lt.ne(1)).or.(nbasis.lt.ne(2))) then
        value = .false.
        go to 101
      end if


*     ***** allocate memory from heap memory ****
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(neq(1)+neq(2)),
     >                     'psi1',psi1(2),psi1(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*n2ft3d),'dn',dn(2),dn(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*n2ft3d),'dnall',dnall(2),dnall(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,(npack0),'dng',dng(2),dng(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,(npack1),'phi1',phi1(2),phi1(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*nbasis),'ee',ee(2),ee(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(n2ft3d),'rho',rho(2),rho(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(n2ft3d),'vc',vc(2),vc(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*n2ft3d),'xce',xce(2),xce(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*n2ft3d),'xcp',xcp(2),xcp(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*n2ft3d),'vall',vall(2),vall(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack0,
     >                     'vl2',vl(2),vl(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,n2ft3d,
     >                     'vl_lr',vl_lr(2),vl_lr(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,n2ft3d,
     >                     'v_field',v_field(2),v_field(1))
      eigs(1) = ee(1)
      eigs(2) = ee(1) + nbasis

*     ***** allocate matrices using GA ****
      if (.not. value) then
        go to 100 
      end if

      if ((taskid.eq.MASTER).and.(oprint)) then
        write(*,1000)
        write(*,1010)
        write(*,1020)
        write(*,1010)
        write(*,1030)
        write(*,1010)
        write(*,1040)
        write(*,1010)
        write(*,1000)
        write(*,1110) 
        call nwpw_message(1)


         write(6,1121) control_boundry(),control_version()
         if (ispin.eq.1) write(6,1130) 'restricted'
         if (ispin.eq.2) write(6,1130) 'unrestricted'
      
         IF (control_gga().eq.0) THEN
            write(6,1131) 'Vosko et al parameterization'
         ELSE IF (control_gga().eq.10) THEN
            write(6,1131)
     >      'PBE96 (White and Bird) parameterization'
         ELSE IF (control_gga().eq.11) THEN
            write(6,1131)
     >      'BLYP (White and Bird) parameterization'
         ELSE IF (control_gga().eq.12) THEN
            write(6,1131)
     >      'revPBE (White and Bird) parameterization'
         ELSE
            write(6,1131) 'unknown parameterization'
            call errquit('bad exchange_correlation',0, INPUT_ERR)
         END IF

        write(*,1117) 
        do ia=1,ion_nkatm_qm()
           write(*,1118) ion_atom_qm(ia),aorbs_norbs(ia)
        end do
        write(*,1119) nbasis

        do n = 1,nbasis
          ii = aorbs_ii(n)
          ia = aorbs_katm(n)
          l = aorbs_l(ia,n)
          m = aorbs_m(ia,n)
          write(*,*) "basis function=",n,"  ii,ia,l,m=",ii,ia,l,m
        end do


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
         do ia=1,ion_nkatm()
           icharge = icharge + ion_natm(ia)*psp_zv(ia)
         end do
         write(6,1159) icharge

        write(*,1220) ispin, ne(1), ne(ispin)
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
         write(6,1300)
         write(6,1305)
         call util_flush(6)

        write(6,*)
      end if
      if ((taskid.eq.MASTER).and.oprint) call current_second(cpu2)



*     ********************************
*     **** generate phaze factors ****
*     ********************************
      call phafac()
      if (psp_semicore(0)) call semicore_density_update()


*     ********************************
*     **** intialize the density  ****
*     ********************************
      call lcao_init_dn(ispin,ne,n2ft3d,dbl_mb(dn(1)),dcpl_mb(phi1(1)))


*     **********************
*     **** generate dng ****
*     **********************
      value = BA_push_get(mt_dbl,n2ft3d,'tmp1',tmp1(2),tmp1(1))
         if (.not. value) call errquit(
     >     'electron_gen_dng_dnall: out of stack memory',0, MA_ERR)

       
      call D3dB_rr_Sum(1,dbl_mb(dn(1)),
     >                   dbl_mb(dn(1)+(ispin-1)*n2ft3d),
     >                   dbl_mb(tmp1(1)))
      call D3dB_r_SMul1(1,scal1,dbl_mb(tmp1(1)))
      call D3dB_rc_fft3f(1,dbl_mb(tmp1(1)))
      call Pack_c_pack(0,dbl_mb(tmp1(1)))
      call Pack_c_Copy(0,dbl_mb(tmp1(1)),dcpl_mb(dng(1)))

*     ********************************************************
*     **** generate dnall - used for semicore corrections ****
*     ********************************************************
      if (psp_semicore(0)) then
         call semicore_density(dbl_mb(tmp1(1)))
         call D3dB_r_SMul1(1,0.5d0,dbl_mb(tmp1(1)))
      else
         call D3dB_r_Zero(1,dbl_mb(tmp1(1)))
      end if
      do ms=1,ispin
        call D3dB_rr_Sum(1,dbl_mb(dn(1)+(ms-1)*n2ft3d),
     >                     dbl_mb(tmp1(1)),
     >                     dbl_mb(dnall(1)+(ms-1)*n2ft3d))
      end do

*     *****************************
*     **** generate potentials ****
*     *****************************
      if (control_version().eq.3) then
            call coulomb_v(dcpl_mb(dng(1)),dbl_mb(vc(1)))
      end if

      if (control_version().eq.4)  then
         call D3dB_rr_Sum(1,dbl_mb(dn(1)),
     >                      dbl_mb(dn(1)+(ispin-1)*n2ft3d),
     >                      dbl_mb(tmp1(1)))

         call coulomb2_v(dbl_mb(tmp1(1)),dbl_mb(vc(1)))
      end if

*    **** xc potential ****
      call v_bwexc_all(control_gga(),n2ft3d,ispin,dbl_mb(dnall(1)),
     >                 dbl_mb(xcp(1)),dbl_mb(xce(1)))
      value = BA_pop_stack(tmp1(2))
      if (.not.value) call errquit(
     >     'electron_gen_dng_dnall: poping stack',1, MA_ERR)



      scal2 = 1.0d0/lattice_omega()

      if (control_version().eq.3) then

*       **** add up k-space potentials, vall = scal2*vl + vc  ****
        call Pack_c_SMul(0,scal2,dcpl_mb(vl(1)),
     >                           dbl_mb(vall(1)))

        call Pack_cc_Sum2(0,dcpl_mb(vc(1)),dbl_mb(vall(1)))

*       **** fourier transform k-space potentials ****
        call Pack_c_unpack(0,dbl_mb(vall(1)))
        call D3dB_cr_fft3b(1,dbl_mb(vall(1)))
        if (field_exist)
     >    call D3dB_rr_Sum2(1,dbl_mb(v_field(1)),dbl_mb(vall(1)))

      else

*       **** add up k-space potentials, vall = scal2*vsr_l    ****
        call Pack_c_SMul(0,scal2,dcpl_mb(vl(1)),
     >                           dbl_mb(vall(1)))

*        **** fourier transform k-space potentials ****
         call Pack_c_unpack(0,dbl_mb(vall(1)))
         call D3dB_cr_fft3b(1,dbl_mb(vall(1)))

         call D3dB_rr_Sum2(1,dbl_mb(vl_lr(1)),dbl_mb(vall(1)))
         call D3dB_rr_Sum2(1,dcpl_mb(vc(1)),dbl_mb(vall(1)))
         if (field_exist)
     >     call D3dB_rr_Sum2(1,dbl_mb(v_field(1)),dbl_mb(vall(1)))

      end if

      if (ispin.eq.2) then
        call D3dB_rr_Sum(1,dbl_mb(vall(1)),
     >                   dbl_mb(xcp(1) +n2ft3d),
     >                   dbl_mb(vall(1)+n2ft3d))
        call D3dB_r_Zero_Ends(1,dbl_mb(vall(1)+n2ft3d))
      end if
      call D3dB_rr_Sum2(1,dbl_mb(xcp(1)),dbl_mb(vall(1)))
      call D3dB_r_Zero_Ends(1,dbl_mb(vall(1)))


*     *******************************
*     *** read psi1 wavefunctions ***
*     *******************************
      write(*,*) "HERA"
      occ1(1) = rho(1)
      call psi_get_ne_occupation(ispin,ne,smearoccupation)
      if (smearoccupation.gt.0) then
      write(*,*) "HERAaa"
         value = value.and.
     >        BA_alloc_get(mt_dbl,(ne(1)+ne(2)),'occ1',occ1(2),occ1(1))
      end if
      write(*,*) "HERB,ispin,ne,smearocc=",ispin,ne,smearoccupation
      call psi_read(ispin,ne,dcpl_mb(psi1(1)),
     >              smearoccupation,dbl_mb(occ1(1)))

      write(*,*) "HERC,ispin,ne,smearocc=",ispin,ne,smearoccupation
*     ********************************
*     *** write psi1 wavefunctions ***
*     ********************************
      call psi_write(ispin,ne,dcpl_mb(psi1(1)),
     >               smearoccupation,dbl_mb(occ1(1)))






      rho_error = 10000.0d0
      if ((taskid.eq.MASTER).and.oprint) call current_second(cpu3)
      if ((taskid.eq.MASTER).and.oprint) CALL nwpw_MESSAGE(8)


 50   if ((taskid.eq.MASTER).and.oprint) call current_second(cpu4)
      if ((taskid.eq.MASTER).and.oprint) CALL nwpw_MESSAGE(3)
      if ((taskid.eq.MASTER).and.oprint) write(*,*)
      if ((taskid.eq.MASTER).and.oprint) write(*,*)


c 50   call lcao_diis_end()
c      call lcao_write_psi(wavefunction_filename,
c     >                    version,
c     >                    ngrid,
c     >                    unita,
c     >                    ispin,ne,
c     >                    psimatrix,dcpl_mb(phi2(1)))
c

*     **** write out excited orbitals ****
      ne_excited(1) = 0
      ne_excited(2) = 0
      ne_excited(1) = control_excited_ne(1)
      if (ispin.eq.2) ne_excited(2) = control_excited_ne(2)
      filename = control_input_epsi()

c      if ((ne_excited(1)+ne_excited(2)).gt.0) 
c     >  call lcao_write_epsi(filename,
c     >                    version,
c     >                    ngrid,
c     >                    unita,
c     >                    ispin,ne,ne_excited,
c     >                    psimatrix,dcpl_mb(phi2(1)))
c

*:::::::::::::::::   report summary of results  :::::::::::::::::::::::

      if ((taskid.eq.MASTER).and.oprint) then
      write(*,*)
      write(*,1410)
      end if

C     **** caluclate number of Density ****
      do ms=1,ispin
         call D3dB_r_dsum(1,dbl_mb(dn(1)+(ms-1)*n2ft3d),total_dn)
         total_dn = total_dn*lattice_omega()
         total_dn = total_dn/dble(ngrid(1)*ngrid(2)*ngrid(3))
         if ((taskid.eq.MASTER).and.oprint) then
         write(*,*) 'Spin=',MS,' Total Orbital Density:', total_dn
         end if
      end do





*:::::::::::::::::::   report consumed cputime  :::::::::::::::::::::::
      if ((taskid.eq.MASTER).and.oprint) call current_second(cpu5)


      if ((taskid.eq.MASTER).and.oprint) then
      t1 = cpu2 - cpu1
      t2 = cpu3 - cpu2
      t3 = cpu4 - cpu3
      t4 = cpu5 - cpu1
      WRITE(*,*)
      WRITE(*,*) '----------------------'
      WRITE(*,*) 'Proglogue        : ',t1
      WRITE(*,*) 'Setup Matrices   : ',t2
      WRITE(*,*) 'Convergence Loop : ',t3
      WRITE(*,*) 'Total            : ',t4
      WRITE(*,*) 'CpuTime/Step     : ',av
      end if
      if ((taskid.eq.MASTER).and.oprint) CALL nwpw_MESSAGE(4)
      




*     ***** free heap memory ****
  100 continue
      value2 = BA_free_heap(vc(2))
      value2 = value2.and.BA_free_heap(xce(2))
      value2 = value2.and.BA_free_heap(xcp(2))
      value2 = value2.and.BA_free_heap(vall(2))
      value2 = value2.and.BA_free_heap(vl_lr(2))
      value2 = value2.and.BA_free_heap(vl(2))
      value2 = value2.and.BA_free_heap(v_field(2))
      value2 = value2.and.BA_free_heap(rho(2))
      value2 = value2.and.BA_free_heap(dn(2))
      value2 = value2.and.BA_free_heap(psi1(2))
      value2 = value2.and.BA_free_heap(dng(2))
      value2 = value2.and.BA_free_heap(phi1(2))
      value2 = value2.and.BA_free_heap(ee(2))
      call Dne_end()
  101 call aorbs_end()

      value = .false.
      wvfnc_new_lcao_density = value
      return

 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*                 LCAO Calculations                *')
 1030 FORMAT(10X,'*       [ NorthWest Chemistry implementation ]     *')
 1040 FORMAT(10X,'*                  Version #2.00                   *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1111 FORMAT(' Reading in ELCIN file(1-yes, 0-no): ',I1)
 1115 FORMAT(/' Pseudopotentials used:')
 1116 FORMAT(5X,A,'  # of L: ',I1)
 1117 FORMAT(/' Atomic orbitals used:')
 1118 FORMAT(5X,A,'  # of orbitals: ',I6)
 1119 FORMAT(/' Number of basis functions:',I4)
 1121 FORMAT(/5X,' boundry conditions   = ',A,'(version', I1,')')
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
 1159 FORMAT(/' total charge=',I2)

 1160 FORMAT(/' ATOMIC COMPOSITION:')
 1170 FORMAT(7(5X,A2,':',I3))
 1180 FORMAT(/' INITIAL POSITION OF IONS:')
 1190 FORMAT(5X, I4, A4  ,' (',3F11.5,' )')
 1200 FORMAT(5X,'  G.C. ',' (',3F11.5,' )')
 1210 FORMAT(5X,' C.O.M.',' (',3F11.5,' )')
 1220 FORMAT(' # of Electrons(ispin=',I1,')  :',I4,' up',I4,' down')
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
 1261 FORMAT(5X,'                   madelung=',f14.8)
 1270 FORMAT(/' technical parameters:')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'fictacious mass=',F10.1)
 1290 FORMAT(5X, ' tolerance=',E8.3,' (energy)',E12.3,
     &        ' (electron)',E12.3,' (ion)')
 1300 FORMAT(//)

 1305 FORMAT(10X,'================ ITERATION =========================')
 1310 FORMAT(I8,E20.10,2E15.5)
 1320 FORMAT(' NUMBER OF ELECTRONS: SPIN UP=',F11.5,'  DOWN=',F11.5)
 1330 FORMAT(/' COMPARISON BETWEEN HAMILTONIAN AND LAMBDA MATRIX')
 1340 FORMAT(I3,2I3,' H=',E16.7,', L=',E16.7,', H-L=',E16.7)
 1350 FORMAT(/' ORTHONORNALITY')
 1360 FORMAT(I3,2I3,E18.7)
 1370 FORMAT(I3)
 1380 FORMAT(' ''',a,'''',I4)
 1390 FORMAT(I3)
 1400 FORMAT(I3,3E18.8/3X,3E18.8)
 1410 FORMAT(10X,'=============  SUMMARY OF RESULTS  =================')
 1420 FORMAT( ' FINAL POSITION OF IONS:')
 1430 FORMAT(/' TOTAL     ENERGY    :',E19.10,' (',E15.5,'/ION)')
 1440 FORMAT( ' TOTAL ORBITAL ENERGY:',E19.10,' (',E15.5,'/ION)')
 1450 FORMAT( ' HARTREE   ENERGY    :',E19.10,' (',E15.5,'/ELECTRON)')
 1460 FORMAT( ' EXC-CORR  ENERGY    :',E19.10,' (',E15.5,'/ELECTRON)')
 1470 FORMAT( ' ION-ION   ENERGY    :',E19.10,' (',E15.5,'/ION)')
 1480 FORMAT( ' PXC       ENERGY    :',E19.10,' (',E15.5,'/ION)')
 1500 FORMAT(/' ORBITAL ENERGIES:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

      end



