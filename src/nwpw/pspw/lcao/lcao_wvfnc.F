*
* $Id$
*
      logical function lcao_wvfnc(oprint_in,it_max,
     >                      wavefunction_filename,
     >             			version,
     >						ngrid,
     >						unita,
     >						ispin,
     >						ne)

      implicit none
#include "errquit.fh"
      logical oprint_in
      integer it_max
      character*50 wavefunction_filename
      integer      version
      integer      ngrid(3)
      real*8       unita(3,3)
      integer      ispin,ne(2)

#include "bafdecls.fh"
#include "global.fh"
 
*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      real*8 ALPHA,talpha
      parameter (ALPHA=0.30d0)

      logical value,value2,oprint
      integer ms,i,ia,icharge,ne_excited(2),nn,l,m,n,ii
      integer dn(2),phi1(2),phi2(2),rho(2),vc(2)
      integer ee(2),eigs(2)
      integer smatrix,kematrix
      integer vlmatrix,vnlmatrix,vhxcmatrix(2),vhxc_new(2)
      integer psimatrix(2),hmatrix(2)
      integer nbasis,nbasis2,n2ft3d,npack1,it
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

      integer      control_excited_ne
      character*50 control_input_epsi
      external     control_excited_ne
      external     control_input_epsi


      Etotal = 0.0d0
      oprint = oprint_in

      call Parallel_taskid(taskid)

      if ((taskid.eq.MASTER).and.oprint) call current_second(cpu1)

      call D3dB_n2ft3d(1,n2ft3d)
      call Pack_npack(1,npack1)

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
     >        BA_alloc_get(mt_dbl,(2*n2ft3d),'dn',dn(2),dn(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,(npack1),'phi1',phi1(2),phi1(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,(npack1),'phi2',phi2(2),phi2(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(2*nbasis),'ee',ee(2),ee(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(n2ft3d),'rho',rho(2),rho(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(n2ft3d),'vc',vc(2),vc(1))
      eigs(1) = ee(1)
      eigs(2) = ee(1) + nbasis

*     ***** allocate matrices using GA ****
      value = value.and.
     >        ga_create(mt_dbl,nbasis,nbasis,'smatrix',
     >                  nbasis,0,smatrix)
      value = value.and.
     >        ga_create(mt_dbl,nbasis,nbasis,'kematrix',
     >                  nbasis,0,kematrix)
      value = value.and.
     >        ga_create(mt_dbl,nbasis,nbasis,'vlmatrix',
     >                  nbasis,0,vlmatrix)
      value = value.and.
     >        ga_create(mt_dbl,nbasis,nbasis,'vnlmatrix',
     >                  nbasis,0,vnlmatrix)
      value = value.and.
     >        ga_create(mt_dbl,nbasis,nbasis,'vhxcmatrix1',
     >                  nbasis,0,vhxcmatrix(1))
      value = value.and.
     >        ga_create(mt_dbl,nbasis,nbasis,'vhxc_new',
     >                  nbasis,0,vhxc_new(1))
      value = value.and.
     >        ga_create(mt_dbl,nbasis,nbasis,'hmatrix1',
     >                  nbasis,0,hmatrix(1))
      value = value.and.
     >        ga_create(mt_dbl,nbasis,nbasis,'psimatrix1',
     >                  nbasis,0,psimatrix(1))

      if (ispin.gt.1) then
      value = value.and.
     >        ga_create(mt_dbl,nbasis,nbasis,'vhxcmatrix2',
     >                  nbasis,0,vhxcmatrix(2))
      value = value.and.
     >        ga_create(mt_dbl,nbasis,nbasis,'vhxc_new2',
     >                  nbasis,0,vhxc_new(2))
      value = value.and.
     >        ga_create(mt_dbl,nbasis,nbasis,'hmatrix2',
     >                  nbasis,0,hmatrix(2))
      value = value.and.
     >        ga_create(mt_dbl,nbasis,nbasis,'psimatrix2',
     >                  nbasis,0,psimatrix(2))
      end if
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

*     **** get ewald energy ****
      Eion = 0.0d0
      if (control_version().eq.3) Eion = ewald_e()

*     **** get free-space ion-ion energy ****
      if (control_version().eq.4) Eion = ion_ion_e()


*     **********************************
*     **** initialize lcao matrices ****
*     **********************************
      call lcao_generate_smatrix(smatrix,
     >                           dcpl_mb(phi1(1)),
     >                           dcpl_mb(phi2(1)))
      call lcao_generate_kematrix(kematrix,
     >                           dcpl_mb(phi1(1)),
     >                           dcpl_mb(phi2(1)))
      call lcao_generate_vlmatrix(vlmatrix,
     >                           dcpl_mb(phi1(1)),
     >                           dcpl_mb(phi2(1)))
      call lcao_generate_vnlmatrix(vnlmatrix,
     >                           dcpl_mb(phi1(1)),
     >                           dcpl_mb(phi2(1)))
      call lcao_init_dn(ispin,ne,n2ft3d,
     >                  dbl_mb(dn(1)),
     >                  dcpl_mb(phi2(1)))
      call lcao_generate_vhxcmatrix(ispin,n2ft3d,
     >                           dbl_mb(dn(1)),
     >                           vhxcmatrix,
     >                           dcpl_mb(phi1(1)),
     >                           dcpl_mb(phi2(1)),
     >                           Ehart,Eexc,pxc)
       sum = ga_ddot(vhxcmatrix(1),vhxcmatrix(1))

c     call lcao_diis_init(5,ispin,nbasis)
      rho_error = 10000.0d0
      if ((taskid.eq.MASTER).and.oprint) call current_second(cpu3)
      if ((taskid.eq.MASTER).and.oprint) CALL nwpw_MESSAGE(8)

      do it=0,it_max-1

         do ms=1,ispin

*          **** generate spin hamiltonian ****
           call ga_add(1.0d0,kematrix,
     >                 1.0d0,vlmatrix,
     >                 hmatrix(ms))
           call ga_add(1.0d0,vnlmatrix,
     >                 1.0d0,hmatrix(ms),
     >                 hmatrix(ms))
           call ga_add(1.0d0,vhxcmatrix(ms),
     >                 1.0d0,hmatrix(ms),
     >                 hmatrix(ms))
   
*          **** diagonalize spin hamiltonian ****
#if defined(PARALLEL_DIAG) & !defined(SCALAPACK)
           call ga_diag(
     >                  hmatrix(ms),
     >                  smatrix,
     >                  psimatrix(ms),
     >                  dbl_mb(eigs(ms)))
#else
           call ga_diag_seq(
     >                  hmatrix(ms),
     >                  smatrix,
     >                  psimatrix(ms),
     >                  dbl_mb(eigs(ms)))
#endif

         end do

*        **** calculate new density ****
         call lcao_new_dn(ispin,ne,n2ft3d,
     >                    dbl_mb(dn(1)),
     >                    psimatrix,
     >                    dcpl_mb(phi2(1)))
         call lcao_generate_vhxcmatrix(ispin,n2ft3d,
     >                           dbl_mb(dn(1)),
     >                           vhxc_new,
     >                           dcpl_mb(phi1(1)),
     >                           dcpl_mb(phi2(1)),
     >                           Ehart,Eexc,pxc)
       sum = ga_ddot(vhxc_new(1),vhxc_new(1))
c       write(*,*) "SUM new:",sum

         call ga_add((1.0d0),vhxc_new(1),
     >              (-1.0d0),vhxcmatrix(1),
     >               hmatrix(1))
         if (ispin.gt.1)
     >   call ga_add((1.0d0),vhxc_new(2),
     >              (-1.0d0),vhxcmatrix(2),
     >               hmatrix(2))

        
         rho_error_old = rho_error
         rho_error = ga_ddot(hmatrix(1),hmatrix(1))
         if (ispin.gt.1)
     >    rho_error = rho_error + ga_ddot(hmatrix(2),hmatrix(2))
         rho_error = rho_error/dble(nbasis*nbasis)

         if (it.eq.0) talpha=ALPHA
         if (rho_error.gt.rho_error_old) talpha=0.5d0*talpha

c        *** orbital energies ***
         Eorbs = 0.0d0
         do ms=1,ispin
         do i=1,ne(ms)
            Eorbs = Eorbs + dbl_mb(eigs(ms)+i-1)
         end do
         end do

C        *** Calculate total energy ***
         olde   = Etotal
         Etotal = 0.0d0
         Etotal = Eorbs + Eion + Eexc - pxc - Ehart 
         de     = dabs(Etotal-olde)



         if ((oprint).and.(taskid.eq.MASTER))
     >     write(*,510) it,Etotal,Eorbs,de,rho_error,talpha
  510 FORMAT(I8,2E20.10,3E15.5)

         if ((de.lt.control_tole()).and.
     >       (rho_error.lt.control_tolc())) then
           if ((oprint).and.(taskid.eq.MASTER))
     >        write(*,*) ' *** tolerance ok.     iteration terminated.'
           go to 50
         end if
c         if ((oprint)) then
c            do ms=1,ispin
c              call D3dB_r_dsum(1,dbl_mb(dn(1)+(ms-1)*n2ft3d),total_dn)
c              total_dn = total_dn*lattice_omega()
c              total_dn = total_dn/dble(ngrid(1)*ngrid(2)*ngrid(3))
c              if (taskid.eq.MASTER) 
c     >         write(*,*) "ms,dn:",ms,total_dn
c            end do
c         end if

         call ga_add(talpha,vhxc_new(1),
     >               (1.0d0-talpha),vhxcmatrix(1),
     >                vhxcmatrix(1))
         if (ispin.gt.1)
     >    call ga_add(talpha,vhxc_new(2),
     >                (1.0d0-talpha),vhxcmatrix(2),
     >                 vhxcmatrix(2))


c        reset_diis = (rho_error.gt.rho_error_old).or.
c    >                (rho_error.gt.1.0d0)
c        call lcao_diis_generate(talpha,vhxc_new,vhxcmatrix,reset_diis)
c        if (reset_diis) then
c           call lcao_diis_end()
c           call lcao_diis_init(5,ispin,nbasis)
c        end if



      end do
      if ((oprint).and.(taskid.eq.MASTER))
     > write(*,*) '*** arrived at the Maximum iteration.   terminated.'

 50   if ((taskid.eq.MASTER).and.oprint) call current_second(cpu4)
      if ((taskid.eq.MASTER).and.oprint) CALL nwpw_MESSAGE(3)
      if ((taskid.eq.MASTER).and.oprint) write(*,*)
      if ((taskid.eq.MASTER).and.oprint) write(*,*)


c 50   call lcao_diis_end()
      call lcao_write_psi(wavefunction_filename,
     >                    version,
     >                    ngrid,
     >                    unita,
     >                    ispin,ne,
     >                    psimatrix,dcpl_mb(phi2(1)))


*     **** write out excited orbitals ****
      ne_excited(1) = 0
      ne_excited(2) = 0
      ne_excited(1) = control_excited_ne(1)
      if (ispin.eq.2) ne_excited(2) = control_excited_ne(2)
      filename = control_input_epsi()

      if ((ne_excited(1)+ne_excited(2)).gt.0) 
     >  call lcao_write_epsi(filename,
     >                    version,
     >                    ngrid,
     >                    unita,
     >                    ispin,ne,ne_excited,
     >                    psimatrix,dcpl_mb(phi2(1)))


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


      if ((taskid.eq.MASTER).and.oprint) then
      write(*,1430) Etotal,Etotal/ion_nion()
      write(*,1440) Eorbs, Eorbs/(ne(1)+ne(ispin))
      write(*,1450) Ehart, Ehart/(ne(1)+ne(ispin))
      write(*,1460) Eexc,  Eexc/(ne(1)+ne(ispin))
      write(*,1470) Eion,  Eion/ion_nion()
      write(*,1480) pxc,   pxc/(ne(1)+ne(ispin))

C     * write out the eigenvalues *
      nn = ne(1) - ne(2)
      EV = 27.2116d0
      WRITE(*,*) ' '
      WRITE(*,*) 'Filled Eigenvalues ...'
      do i=0,nn-1
         WRITE(*,1510) dbl_mb(eigs(1)+ne(1)-i),
     >                 dbl_mb(eigs(1)+ne(1)-i)*EV
      end do
      do i=0,ne(2)-1
         WRITE(*,1510) dbl_mb(eigs(1)+ne(1)-i-nn),
     >                 dbl_mb(eigs(1)+ne(1)-i-nn)*EV,
     >                 dbl_mb(eigs(2)+ne(2)-i),   
     >                 dbl_mb(eigs(2)+ne(2)-i)*EV
      end do
      if (ispin.eq.1) then
       WRITE(*,*)
       WRITE(*,*) 'Virtual Eigenvalues ...'
       do i=0,(nbasis-ne(1))-1
         WRITE(*,1510) dbl_mb(eigs(1)+nbasis-i),
     >                 dbl_mb(eigs(1)+nbasis-i)*EV
       end do
      else
       write(*,*)
       write(*,*) "Unrestricted Virtuals Eigenvalues ..."
       do i=1,(nbasis-ne(1))-1
         WRITE(*,1510) dbl_mb(eigs(1)+ne(1)+i),
     >                 dbl_mb(eigs(1)+ne(1)+i)*EV,
     >                 dbl_mb(eigs(2)+ne(2)+i),   
     >                 dbl_mb(eigs(2)+ne(2)+i)*EV
       end do
      end if

      end if




*:::::::::::::::::::   report consumed cputime  :::::::::::::::::::::::
      if ((taskid.eq.MASTER).and.oprint) call current_second(cpu5)


      if ((taskid.eq.MASTER).and.oprint) then
      t1 = cpu2 - cpu1
      t2 = cpu3 - cpu2
      t3 = cpu4 - cpu3
      t4 = cpu5 - cpu1
      av = t3/it_max
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
      value2 = BA_free_heap(rho(2))
      value2 = BA_free_heap(dn(2))
      value2 = BA_free_heap(phi1(2))
      value2 = BA_free_heap(phi2(2))
      value2 = BA_free_heap(ee(2))
      value2 = ga_destroy(smatrix)
      value2 = ga_destroy(kematrix)
      value2 = ga_destroy(vlmatrix)
      value2 = ga_destroy(vnlmatrix)
      value2 = ga_destroy(vhxcmatrix(1))
      value2 = ga_destroy(vhxc_new(1))
      value2 = ga_destroy(hmatrix(1))
      value2 = ga_destroy(psimatrix(1))
      if (ispin.gt.1) then
      value2 = ga_destroy(vhxcmatrix(2))
      value2 = ga_destroy(vhxc_new(2))
      value2 = ga_destroy(hmatrix(2))
      value2 = ga_destroy(psimatrix(2))
      end if
  101 call aorbs_end()

c     value = .false.
      lcao_wvfnc = value
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



