      subroutine aor_get_giaoints(rtdb,basis,geom, nspin, lifetime,
     &  g_vecF1, g_vecF1_im, g_vecB1, g_smat1, g_sket1, giaowsre,
     &  giaowsim) 


c     This routine calculates the GIAO fock matrix derivative wrt the
c     components of the magnetic field, and the GIAO overlap matrix
c     derivatives needed for optical rotation calculations per Eq. (57)
c     of JCP 123 (2005), 114103.  It also calculates the third term in
c     W(beta) from the W(beta)S(alpha) contribution in that equation,
c     which can be neglected for non-hybrid functionals.
c
c     the output GAs should have dimensions of (nbf, nbf, 3).  all
c     output arrays must be initialized with zeros by the calling
c     routine.
c
c     Written by J. Autschbach, SUNY Buffalo

      implicit none
      
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"
#include "case.fh"

c     subroutine arguments:
      
      integer rtdb              ! [in] rtdb handle
      integer basis             ! [in] basis handle
      integer geom              ! [in] geometry handle
      integer nspin             ! [in] number of spin components
      logical lifetime          ! [in] response with damping?
      integer g_vecF1(2,2)      ! [in] E-field perturbed MOs Re
      integer g_vecF1_im(2,2)   ! [in] E-field perturbed MOs Im
      
      integer g_vecB1(2)        ! [out] GIAO Fock matrix derivative
      integer g_smat1           ! [out] GIAO overlap derivative
      integer g_sket1           ! [out] GIAO overlap ket derivative

      double precision
     &  giaowsre(3,3),    ! [out] 
     &  giaowsim(3,3)     ! [out]

c     local variables:
      
      integer nclosed(2), nopen(2), nvirt(2), ndens, nbf, nmo
      
      integer ixy, ix, iy, iatom, iocc, ifld, ioff, idir, jdir
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer dlo(3), dhi(3), elo(3), ehi(3)
      
      integer l_occ, k_occ, l_eval, k_eval
      
      integer g_dens(3), g_s10, g_d1, g_d1_im, g_fock, g_fock_im

      integer vectors(2), geomnew, i, j, ij, g_xc(3)
      integer g_s1mo, g_f1mo, g_f1mo_im
      integer g_temp
      
      integer ndata,m1,m2
      double precision atn, tol2e, val, rtemp
      double precision jfac(3),kfac(3),a(6),xfac
      character*3 scftyp
      character*16 tag
      character*32 element
      character*256 cphf_rhs, cphf_sol
      character*2 symbol

      integer ndir,    ! nr directions (x,y,z)
     &        ntot,    ! sum_{i=1,npol} nocc(i)*nvirt(i)
     &        ispin,disp,shift,nind_jk,
     &        nocc(2), ! store nr occupations 
     &        npol     ! nr of polarizations =1 (  restricted) 
                       !                     =2 (unrestricted) calc    

      integer nbq, nextbq, ncosbq

      integer nat
      parameter (nat=1)
      double precision origin(3)
      data origin/0d0,0d0,0d0/
   
      logical  oskel, status, me
      integer dbg
      logical  xc_gotxc
      external xc_gotxc
      double precision ppm
      data tol2e   /1.0d-10/
      
c     ==================================================================

      ndir=3                    ! nr directions (x,y,z)
      
      if (ga_nodeid().eq.0) write(luout,'(/1x,a/)')
     &  'Calculating some arrays needed for GIAO optical rotation'

      me = (ga_nodeid().eq.0)
      dbg = 0                   ! debug level

      if (dbg>0 .and. me) then
        write (luout,*) 'aor_giaoints: xc_gotxc =',xc_gotxc()
        write (luout,*) 'aor_giaoints: use_theory =',use_theory
        write (luout,*) 'nspin, lifetime=',nspin,lifetime
      end if

c     there is a possibility that this routine gets called from tddft
c     in which case use_theory is not set. We set it to 'dft',
c     assuming that we are indeed calling from some DFT response code

      if ((use_theory.ne.'scf').and.(use_theory.ne.'dft')) then
        if (ga_nodeid().eq.0) write (luout,*)
     &     'aor_giaoints: assuming DFT/TDDFT'
        use_theory = 'dft'
      end if
c
c     Current CPHF does not handle symmetry 
c     Making C1 geometry and store it on rtdb
c     
      oskel = .false.
c
c     If DFT get part of the exact exchange defined
      
      xfac = 1.0d0
      if (use_theory.eq.'dft') xfac = bgj_kfac()

      do idir = 1,ndir
        jfac(idir) =  0.0d0       
        kfac(idir) = -1.0d0*xfac  
      enddo
c     
c     Integral initialization
      
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
      call scf_get_fock_param(rtdb,tol2e)
      
c     status = rtdb_parallel(.true.)

c     ---------------------------------------------------
c     calculate the overlap matrix derivatives needed for
c     optical rotation calculations with GIAOs
c     ---------------------------------------------------
      
c     smat1 from d/dB <mu | nu> :
      
      call int_giao_1ega(basis,basis,g_smat1,'s10',origin,nat,oskel)
      
c     sket1 from <mu | d/dB nu> :
      call int_giao_1ega(basis,basis,g_sket1,'srxRb',origin,nat,oskel)
      
      if (dbg>0 .and. me) write (luout,*)
     &  'aor_giaoints: GIAO S(B), S-ket done' 

c     ----------------------------------------
c     calculate the GIAO B-field Fock operator
c     derivative and store in g_vecB1:
c     -----------------------------------------
      
c     Get Unperturbed MOs and energies
      
      if (.not. bas_numbf(basis,nbf)) call
     &  errquit('giaoints: could not get nbf',0, BASIS_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'occ num',l_occ,k_occ)) call
     &  errquit('giaoints: ma_push_get failed k_occ',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'eigenval',l_eval,k_eval)) call
     &  errquit('giaoints: ma_push_get failed k_eval',0,MA_ERR)
      
      call hnd_prp_vec_read(rtdb,geom,basis, ! in : handles
     &  nbf,                    ! out: nr basis fcts.
     &  nclosed,nopen,nvirt,    ! out: occ/unocc
     &  scftyp,                 ! out: type calc
     &  vectors,                ! out: MO vectors
     &  dbl_mb(k_occ),          ! out: occupations
     &  dbl_mb(k_eval),         ! out: DFT energies
     &  nmo)                    ! out: nr MOs
      
      call get_nocc(rtdb,       ! in : rtdb handle
     &  nocc,                   ! out: nr occupations
     &  npol,                   ! out: nr of polarization
     &  nclosed,                ! in : nr closed shells
     &  nopen,                  ! in : nr open shells
     &  nvirt,                  ! in : nr virtual MOs
     &  scftyp,                 ! in : string = UHF or RHF
     &  ntot)                   ! out: sum_{i,npol} nocc(i)*nvirt(i)
      
c     get the unperturbed density matrix
      
      call hnd_prp_get_dens(rtdb,geom,basis,     
     &  g_dens,ndens,        
     &  scftyp,              
     &  nclosed,nopen,nvirt) 
      
      if (dbg>0 .and. me) write (luout,*)
     &  'aor_giaoints: unpertubed Pmat done'

c     sanity checks:
      
      if (npol.lt.1 .or. npol.gt.2) call errquit(
     &  'aor_giaoints: npol out of range',0,CALC_ERR)
      
      if (npol.ne.nspin) call errquit(
     &  'aor_giaoints: npol not in sysnc with calling routine',
     &  npol,CALC_ERR)
         
c     calculate the one-electron integrals for the GIAO
c     Fock operator derivative and store in g_vecB1.
      
c     This includes potential energy terms from point charges and
c     a solvent model, if necessary.
      
      call int_giao_1ega(basis,basis,g_vecB1(1),'l10',origin,nat,oskel)
      
      call int_giao_1ega(basis,basis,g_vecB1(1),'tv10',origin,nat,oskel)

c     Get external and cosmo bq contribution
      
      nbq    = 0
      nextbq = 0
      ncosbq = 0
      if(geom_extbq_on()) nextbq = geom_extbq_ncenter()
      nbq = nextbq ! external bq's
      if (rtdb_get(rtdb,'cosmo:nefc',mt_int,1,ncosbq))
     &    nbq = ncosbq ! cosmo bq's
      if (nextbq.gt.0.and.ncosbq.gt.0)
     &  nbq = nextbq + ncosbq   ! tally up cosmo and external bqs
      
      if (nbq.gt.0) then
        call int_giao_1ega(basis,basis,g_vecB1(1),'bq10',
     &    origin,nat,oskel)
        if (dbg>0 .and. me) write (luout,*)
     &    'aor_giaoints: bq terms added'
      end if
      
      if (npol.eq.2) call ga_copy(g_vecB1(1), g_vecB1(2))
      
      if (dbg>0 .and. me) write (luout,*)
     &  'aor_giaoints: F(B) GIAO 1e done'
      
c     add the two-electron part of the GIAO Fock matrix
c     derivative to g_vecB1:
      
      alo(1) = -1 
      alo(2) = -1
      alo(3) =  1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = ndir*npol
      
      if (.not.nga_create(MT_DBL,ndir,ahi,'f_giao',
     &  alo,g_fock)) call 
     &  errquit('aor_giaoints: nga_create failed g_fock',1,GA_ERR)
      call ga_zero(g_fock)
      
      if(use_theory.eq.'dft') then
        if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     &    call errquit('aor_giaoints: rtdb_put of xc_active failed',0,
     &    RTDB_ERR)
        
        ifld = 4
        if (.not. rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,ifld))
     &    call errquit('aor_giaoints: rtdb_put failed',1,RTDB_ERR)
      endif
      
      call new_giao_2e(geom,basis,nbf,tol2e,
     &  g_dens, g_fock, xfac, npol)      
      
      do ispin = 1,npol
        alo(1) = 1
        ahi(1) = nbf
        alo(2) = 1
        ahi(2) = nbf
        alo(3) = 1
        ahi(3) = 3
        disp=ndir*(ispin-1)
        blo(1) = 1
        bhi(1) = nbf
        blo(2) = 1
        bhi(2) = nbf
        blo(3) = disp+1
        bhi(3) = disp+ndir
        
       call nga_add_patch(
     &   1.0d0, g_vecB1(ispin),alo,ahi,
     &   1.0d0,g_fock,blo,bhi,
     &   g_vecB1(ispin),alo,ahi)
       
      end do                    ! ispin

      if (dbg>0 .and. me) write (luout,*)
     &  'aor_giaoints: F(B) GIAO 2e done'

      do ispin = 1,npol
        if (.not.ga_destroy(g_dens(ispin))) call 
     &    errquit('giaoints: ga_destroy failed g_dens',0,GA_ERR)
      end do
      
      if (.not.ga_destroy(g_fock)) call 
     &  errquit('giaoints: ga_destroy failed g_fock',1,GA_ERR)      
      
c     done with Fock matrix B-field derivative. 
      
c     ----------------------------------------------------
c     For HF and hybrid functionals, we need an additional
c     term that depends on the occ-occ block of the
c     electric-field perturbed Fock matrix:
c     -----------------------------------------------------


      clo(1) = 3
      clo(2) = nbf
      clo(3) = nbf
      chi(1) = 1  
      chi(2) = -1 
      chi(3) = -1      
      
      if (.not.nga_create(MT_DBL,ndir,clo,'d1', chi,g_d1)) call 
     &  errquit('aor_giaoints: nga_create failed g_d1',1,GA_ERR)
      call ga_zero(g_d1)

      if (.not.nga_create(MT_DBL,ndir,clo,'f_fock', chi,g_fock)) call 
     &  errquit('aor_giaoints: nga_create failed g_fock',1,GA_ERR)
      call ga_zero(g_fock)

      if (lifetime) then
        if (.not.nga_create(MT_DBL,ndir,clo,'d1_im', chi,g_d1_im)) call 
     &    errquit('aor_giaoints: nga_create failed g_d1_im',1,GA_ERR)
        call ga_zero(g_d1_im)
        
        if (.not.nga_create(MT_DBL,ndir,clo,'f_fock_im',
     &    chi,g_fock_im)) call 
     &    errquit('aor_giaoints: nga_create failed g_fock_im',1,GA_ERR)
        call ga_zero(g_fock_im)
      end if                    ! lifetime

      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 2))
     &  call errquit('hess_cphf: rtdb_put of calc_type failed',0,
     &  RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     &  call errquit('hess_cphf: rtdb_put of j_derfit failed',0,
     &  RTDB_ERR)

      if (dbg>1) then
        clo(1) = nbf
        clo(2) = nbf
        clo(3) = 1
        chi(1) = 1  
        chi(2) = -1 
        chi(3) = -1

        if (.not.nga_create(MT_DBL,ndir,clo,'temp', chi,g_temp)) call 
     &    errquit('aor_giaoints: nga_create failed g_temp',1,GA_ERR)
        call ga_zero(g_temp)
      end if                    ! dbg>1

      do ispin = 1,npol
        
        if (dbg>0 .and. me) write(luout,*)
     &    'goaoints: nocc =',nocc(ispin)
        
        alo(1) = -1 
        alo(2) = -1
        alo(3) =  1
        ahi(1) = nocc(ispin)
        ahi(2) = nocc(ispin)
        ahi(3) = ndir
        
        if (.not.nga_create(MT_DBL,ndir,ahi,'s1mo',
     &    alo,g_s1mo)) call 
     &    errquit('aor_giaoints: nga_create failed s1mo',1,GA_ERR)
        call ga_zero(g_s1mo)
        
        if (.not.nga_create(MT_DBL,ndir,ahi,'f1mo', alo,g_f1mo)) call 
     &    errquit('aor_giaoints: nga_create failed f1mo',1,GA_ERR)
        call ga_zero(g_f1mo)

        if (lifetime) then
          if (.not.nga_create(MT_DBL,ndir,ahi,'f1mo_im', alo,
     &      g_f1mo_im)) call 
     &    errquit('aor_giaoints: nga_create failed f1mo_im',1,GA_ERR)
        call ga_zero(g_f1mo_im)
        end if                  ! lifetime

c       calculate S(B) derivative in MO basis, occ-occ block only.
c       use g_dens and g_fock as scratch space
        
        alo(1) = 1
        ahi(1) = nbf
        alo(2) = 1
        ahi(2) = nbf
        
        blo(1) = 1
        bhi(1) = nbf
        blo(2) = 1
        bhi(2) = nocc(ispin)
        blo(3) = 1
        bhi(3) = 1

c       this is for the transpose of the MO matrix
        dlo(1) = 1
        dhi(1) = nocc(ispin)
        dlo(2) = 1
        dhi(2) = nbf
        dlo(3) = 1
        dhi(3) = 1

        clo(1) = 1
        chi(1) = nocc(ispin)
        clo(2) = 1
        chi(2) = nocc(ispin)

        elo(2) = 1
        ehi(2) = nbf
        elo(3) = 1
        ehi(3) = nocc(ispin)
        elo(1) = 1
        ehi(1) = 1        
        
        do idir = 1,ndir
          alo(3) = idir
          ahi(3) = idir
          clo(3) = idir
          chi(3) = idir

          call ga_zero(g_fock)

          call nga_matmul_patch ('n','n', 1.d0, 0.d0,
     &      g_smat1       , alo, ahi,
     &      vectors(ispin), blo, bhi,
     &      g_fock        , elo, ehi)

          if (dbg>0 .and. me) write(luout,*) 'goaoints: s1 trafo step 1'
          
          call nga_matmul_patch ('t','n', 1.d0, 0.d0,
     &      vectors(ispin), dlo, dhi,
     &      g_fock        , elo, ehi,          
     &      g_s1mo        , clo, chi)
          
          if (dbg>0 .and. me) write(luout,*) 'goaoints: s1 trafo step 2'
          
          if (dbg>1) call nga_print_patch(g_s1mo,clo,chi,1)

        end do                  ! idir

c       calculate the asymmetric part of the density
c       matrix perturbation due to the dynamic electric field

        call ga_zero(g_d1)
        if (lifetime)  call ga_zero(g_d1_im)

        alo(1) = 1
        ahi(1) = nbf
        alo(2) = 1
        ahi(2) = nocc(ispin)
        alo(3) = 1
        ahi(3) = 1

        blo(1) = 1
        bhi(1) = nocc(ispin)
        blo(2) = 1
        bhi(2) = nbf
        blo(3) = 1
        bhi(3) = 1

        dlo(1) = 1
        dhi(1) = nbf
        dlo(2) = 1
        dhi(2) = nocc(ispin)

        elo(1) = 1
        ehi(1) = nocc(ispin)
        elo(2) = 1
        ehi(2) = nbf
        
        clo(2) = 1
        chi(2) = nbf
        clo(3) = 1
        chi(3) = nbf        
                
        do idir = 1,ndir
          
          dlo(3) = idir
          dhi(3) = idir

          elo(3) = idir
          ehi(3) = idir

c         note: in shell_fock_build the index for the set of
c         density and fock matrices is the fist, not last
          
          clo(1) = idir
          chi(1) = idir

          call nga_matmul_patch ('n','t', -0.5d0, 0.d0,
     &      vectors(ispin) , alo, ahi,
     &      g_vecF1(ispin,2)   , elo, ehi,
     &      g_d1           , clo, chi)
          
          if (dbg>0.and.me) write(luout,*) 'goaoints: dens trafo step 1'

          call nga_matmul_patch ('n','t', -0.5d0, 1.d0,
     &      g_vecF1(ispin,1)   , dlo, dhi,
     &      vectors(ispin) , blo, bhi,
     &      g_d1           , clo, chi)

          if (dbg>0.and.me) write(luout,*) 'goaoints: dens trafo step 2'

          call nga_matmul_patch ('n','t', 0.5d0, 1.d0,
     &      vectors(ispin) , alo, ahi,
     &      g_vecF1(ispin,1)   , elo, ehi,
     &      g_d1           , clo, chi)
          
          if (dbg>0.and.me) write(luout,*) 'goaoints: dens trafo step 3'

          call nga_matmul_patch ('n','t', 0.5d0, 1.d0,
     &      g_vecF1(ispin,2)   , dlo, dhi,
     &      vectors(ispin) , blo, bhi,
     &      g_d1           , clo, chi)
          
          if (dbg>0.and.me) write(luout,*) 'goaoints: dens trafo step 4'

         if (lifetime) then
c          if (.false.) then
            
            
            call nga_matmul_patch ('n','t', 0.5d0, 0.d0,
     &        vectors(ispin) , alo, ahi,
     &        g_vecF1_im(ispin,2)   , elo, ehi,
     &        g_d1_im           , clo, chi)
            
            if (dbg>0.and.me) write(luout,*) 'densim trafo step 1'
            
            call nga_matmul_patch ('n','t', -0.5d0, 1.d0,
     &        g_vecF1_im(ispin,1)   , dlo, dhi,
     &        vectors(ispin) , blo, bhi,
     &        g_d1_im           , clo, chi)
            
            if (dbg>0.and.me) write(luout,*) 'densim trafo step 2'
            
            call nga_matmul_patch ('n','t', -0.5d0, 1.d0,
     &        vectors(ispin) , alo, ahi,
     &        g_vecF1_im(ispin,1)   , elo, ehi,
     &        g_d1_im           , clo, chi)
            
            if (dbg>0.and.me) write(luout,*) 'densim trafo step 3'
            
            call nga_matmul_patch ('n','t', 0.5d0, 1.d0,
     &        g_vecF1_im(ispin,2)   , dlo, dhi,
     &        vectors(ispin) , blo, bhi,
     &        g_d1_im           , clo, chi)
            
            if (dbg>0.and.me) write(luout,*) 'densim trafo step 4'

          end if                ! lifetime
          
        end do                  ! idir

        if (dbg>1) then
          alo(2) = 1
          ahi(2) = nbf
          alo(3) = 1
          ahi(3) = nbf
          blo(1) = 1
          bhi(1) = nbf
          blo(2) = 1
          bhi(2) = nbf
          blo(3) = 1
          bhi(3) = 1
          
          do idir = 1,3
            alo(1) = idir
            ahi(1) = idir
            call nga_copy_patch('n',g_d1,alo,ahi, g_temp,blo,bhi)
            if (me) write(luout,*) 'D1 Re AO, idir = ',idir
            call nga_print_patch(g_temp,blo,bhi,1)
           if (lifetime) then
c            if (.false.) then
              call nga_copy_patch('n',g_d1_im,alo,ahi, g_temp,blo,bhi)
              if (me) write(luout,*) 'D1 Im AO, idir = ',idir
              call nga_print_patch(g_temp,blo,bhi,1)
            end if
          end do
        end if
        
c       calculate the Fock matrix
c       Note: Just the exchange: jfac = 0.d0 (see above)

        call ga_zero(g_fock)
        if (lifetime) call ga_zero(g_fock_im)

        ifld = 3
        
        if (.not.cam_exch) then
          call shell_fock_build(geom, basis, 0, ifld,
     &      jfac, kfac, tol2e, g_d1, g_fock, .false.)
        else
          call shell_fock_build_cam(geom, basis, 0, ifld,
     %      jfac, kfac, tol2e, g_d1, g_fock, .false.)
        end if

        if (lifetime) then
          if (.not.cam_exch) then
            call shell_fock_build(geom, basis, 0, ifld,
     &        jfac, kfac, tol2e, g_d1_im, g_fock_im, .false.)
          else
            call shell_fock_build_cam(geom, basis, 0, ifld,
     %        jfac, kfac, tol2e, g_d1_im, g_fock_im, .false.)
          end if
        end if                  ! lifetime

        if (dbg>0.and.me) write(luout,*) 'goaoints: fock_build done'

        if (dbg>1) then
          alo(2) = 1
          ahi(2) = nbf
          alo(3) = 1
          ahi(3) = nbf
          blo(1) = 1
          bhi(1) = nbf
          blo(2) = 1
          bhi(2) = nbf
          blo(3) = 1
          bhi(3) = 1
          
          do idir = 1,3
            alo(1) = idir
            ahi(1) = idir
            call nga_copy_patch('n',g_fock,alo,ahi, g_temp,blo,bhi)
            if (me) write(luout,*) 'F1 Re AO, idir = ',idir
            call nga_print_patch(g_temp,blo,bhi,1)
           if (lifetime) then
c            if (.false.) then
              call nga_copy_patch('n',g_fock,alo,ahi, g_temp,blo,bhi)
              if (me) write(luout,*) 'F1 Im AO, idir = ',idir
              call nga_print_patch(g_temp,blo,bhi,1)
            end if
          end do
        end if
          

c       transform to the MO basis, occ-occ block. use g_d1
c       for scratch space

        alo(2) = 1
        ahi(2) = nbf
        alo(3) = 1
        ahi(3) = nbf
        
        blo(1) = 1
        bhi(1) = nbf
        blo(2) = 1
        bhi(2) = nocc(ispin)
        blo(3) = 1
        bhi(3) = 1
        
c       this is for the transpose of the MO matrix
        dlo(1) = 1
        dhi(1) = nocc(ispin)
        dlo(2) = 1
        dhi(2) = nbf
        dlo(3) = 1
        dhi(3) = 1

        elo(2) = 1
        ehi(2) = nbf
        elo(3) = 1
        ehi(3) = nocc(ispin)
        elo(1) = 1
        ehi(1) = 1        

        clo(1) = 1
        chi(1) = nocc(ispin)
        clo(2) = 1
        chi(2) = nocc(ispin)
        
        do idir = 1,ndir
          
          alo(1) = idir
          ahi(1) = idir
          
          clo(3) = idir
          chi(3) = idir
          
          call ga_zero(g_d1)

          call nga_matmul_patch ('n','n', 1.d0, 0.d0,
     &      g_fock       , alo, ahi,
     &      vectors(ispin), blo, bhi,
     &      g_d1        , elo, ehi)

          if (dbg>0.and.me) write(luout,*) 'goaoints: f1 trafo step 1'
          
          call nga_matmul_patch ('t','n', 1.d0, 0.d0,
     &      vectors(ispin), dlo, dhi,
     &      g_d1        , elo, ehi,          
     &      g_f1mo        , clo, chi)
          
          if (dbg>0.and.me) write(luout,*) 'goaoints: f1 trafo step 2'
          
          if (dbg>1.and.me)
     &      call nga_print_patch(g_f1mo,clo,chi,1)

          if (lifetime) then
c          if (.false.) then            
            
            call ga_zero(g_d1_im)
            
            call nga_matmul_patch ('n','n', 1.d0, 0.d0,
     &        g_fock_im        , alo, ahi,
     &        vectors(ispin), blo, bhi,
     &        g_d1_im          , elo, ehi)

            if (dbg>0.and.me) write(luout,*) 'f1im trafo step 1'
            
            call nga_matmul_patch ('t','n', 1.d0, 0.d0,
     &        vectors(ispin), dlo, dhi,
     &        g_d1_im          , elo, ehi,          
     &        g_f1mo_im     , clo, chi)
          
          if (dbg>0.and.me) write(luout,*) 'f1im trafo step 2'
          
          if (dbg>1.and.me)
     &      call nga_print_patch(g_f1mo_im,clo,chi,1)
          end if                ! lifetime

        end do                  ! idir


        alo(1) = 1
        ahi(1) = nocc(ispin)
        alo(2) = 1
        ahi(2) = nocc(ispin)
        
        blo(1) = 1
        bhi(1) = nocc(ispin)
        blo(2) = 1
        bhi(2) = nocc(ispin)

        do idir = 1,ndir
         do jdir = 1,ndir

            alo(3) = idir
            ahi(3) = idir     
            blo(3) = jdir
            bhi(3) = jdir

            rtemp = nga_ddot_patch(
     &        g_f1mo, 't', alo, ahi,
     &        g_s1mo, 'n', blo, bhi)

            giaowsre(idir,jdir) = giaowsre(idir,jdir) + rtemp

           if (lifetime) then
c            if (.false.) then
              
              rtemp = nga_ddot_patch(
     &          g_f1mo_im, 't', alo, ahi,
     &          g_s1mo, 'n', blo, bhi)
              
              giaowsim(idir,jdir) = giaowsim(idir,jdir) + rtemp 
              
            end if              ! lifetime

          end do                ! jdir
          
        end do                  ! idir

        
        if (.not.ga_destroy(g_f1mo)) call 
     &    errquit('giaoints: ga_destroy failed g_f1mo',1,GA_ERR)
        if (.not.ga_destroy(g_s1mo)) call 
     &    errquit('giaoints: ga_destroy failed g_s1mo',1,GA_ERR)

        if (lifetime) then
          if (.not.ga_destroy(g_f1mo_im)) call 
     &      errquit('giaoints: ga_destroy failed g_f1moim',1,GA_ERR)
        end if
        
      end do                    ! ispin

      if (dbg>0.and.me) then
        write(luout,*) 'goaoints: giaowsre: '
        call output (giaowsre,1,3,1,3,3,3,1)
        write(luout,*) 'goaoints: giaowsim: '
        call output (giaowsim,1,3,1,3,3,3,1)
      end if
        

c     done with Fock matrices. deallocate GAs

      if (dbg>1) then
        if (.not.ga_destroy(g_temp)) call 
     &    errquit('giaoints: ga_destroy failed g_temp',1,GA_ERR)
      end if                    ! dbg>1

      if (.not.ga_destroy(g_fock)) call 
     &  errquit('giaoints: ga_destroy failed g_fock',1,GA_ERR)

      if (.not.ga_destroy(g_d1)) call 
     &  errquit('giaoints: ga_destroy failed g_d1',1,GA_ERR)

      if (lifetime) then
        if (.not.ga_destroy(g_fock_im)) call 
     &    errquit('giaoints: ga_destroy failed g_fock_im',1,GA_ERR)
        
        if (.not.ga_destroy(g_d1_im)) call 
     &    errquit('giaoints: ga_destroy failed g_d1_im',1,GA_ERR)
      end if                    ! lifetime
        
      do ispin=1,npol
        if (.not.ga_destroy(vectors(ispin))) call 
     &    errquit('giaoints: ga_destroy failed vectors',0,GA_ERR)
      enddo
      
      if(use_theory.eq.'dft') then
        ifld = 0
        if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, ifld))
     $    call errquit('aor_giaoints: rtdb_put failed',1,RTDB_ERR)
        if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $    call errquit('aor_giaoints: rtdb_put of xc_active failed',1,
     &    RTDB_ERR)
      endif
     
      

      
      call schwarz_tidy()
      call int_terminate()
      
      call ga_sync()

c     all done.
      
      return

c     =================================================================

 7000 format(/,10x,'B-field perturbed MOs cannot be calculated for',
     1      ' UHF or ROHF wave functions')
 9983 format(6x,i1,3x,3f12.4)
 9985 format(10x,3f12.4,/)
      
      end
