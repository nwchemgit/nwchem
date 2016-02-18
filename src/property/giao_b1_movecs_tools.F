      subroutine update_rhs_fock2e(
     &                      g_rhs,  ! in/out: 
     &                      vectors,! in : MO vectors
     &                      rtdb,   ! in : rtdb  handle
     &                      basis,  ! in : basis handle
     &                      geom,   ! in : geom  handle
     &                      g_dens, ! in : e-density
     &                      nocc,   ! in : nr. occ  shells
     &                      nvirt,  ! in : nr. virt shells
     &                      npol,   ! in : nr. polarizations
     &                      nbf,    ! in : nr. basis functions
     &                      nmo,    ! in : nr. MOs   
     &                      xfac,   ! in : exchange factor
     &                      tol2e,  ! in : tolerance coeff.
     &                      debug)  ! in : logical for debugging
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)

       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "prop.fh"
      integer g_fock,g_rhs
      integer rtdb,basis,geom
      integer npol,ispin
      integer vectors(npol),g_dens(3)
      integer ifld,ndir,nbf,nmo,disp,shift,      
     &        nocc(npol),nvirt(npol)
      double precision tol2e,xfac
      integer alo(3), ahi(3),
     &        blo(3), bhi(3) 
      logical debug
      external new_giao_2e,giao_aotomo
      ndir=3 ! = nr directions (x,y,z)
c     Remaining term is Perturbed (GIAO) two-electron term times
c     Unperturbed density Calculate Sum(r,s) D0(r,s) * G10(m,n,r,s) in
c     AO basis
      alo(1) = -1 
      alo(2) = -1
      alo(3) =  1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = ndir*npol
      if (.not.nga_create(MT_DBL,ndir,ahi,'Fock matrix',
     &    alo,g_fock)) call 
     &    errquit('giao_b1: nga_create failed g_fock',0,GA_ERR)
      call ga_zero(g_fock)
      if(use_theory.eq.'dft') then
         ifld = 4
         if (.not. rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,ifld))
     $      call errquit('giao_b1: rtdb_put failed',0,RTDB_ERR)
      endif

       call new_giao_2e(geom,basis,nbf,tol2e,
     &                  g_dens, !  in: e-denstiy 
     &                  g_fock, ! out: fock matrix
     &                  xfac,
     &                  npol)

      if(use_theory.eq.'dft') then
         ifld = 0
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, ifld))
     $      call errquit('giao_b1: rtdb_put failed',0,RTDB_ERR)
         if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $       call errquit('giao_b1: rtdb_put of xc_active failed',0,
     &       RTDB_ERR)
      endif
c
c     Transform to MO basis and add to right-hand-side
      call giao_aotomo(g_fock,vectors,nocc,nvirt,npol,ndir,nbf)
      do ispin=1,npol
       disp=ndir*(ispin-1)
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = disp+1
       ahi(3) = disp+ndir
       shift=nocc(1)*nvirt(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = ndir
       call nga_add_patch(1.0d0, g_rhs,blo,bhi,
     &                    1.0d0,g_fock,alo,ahi,
     &                           g_rhs,blo,bhi)
      enddo ! end-loop-ispin
      if (debug) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '------- g_fock2e-nw ---- START'
        call ga_print(g_fock)
       if (ga_nodeid().eq.0)
     &  write(*,*) '------- g_fock2e-nw ---- END'
      endif ! end-if-debug
      if (.not.ga_destroy(g_fock)) call 
     &    errquit('giao_b1: ga_destroy failed g_fock',0,GA_ERR)
      return
      end

      subroutine update_rhs_shfock(g_rhs,  ! in/out: RHS used for cphf2/3
     &                             g_d1,   ! in    :
     &                             vectors,! in    : MO vectors
     &                             rtdb,   ! in    : rtdb  handle
     &                             geom,   ! in    : geom  handle
     &                             basis,  ! in    : basis handle 
     &                             jfac,   ! in    : exch factors
     &                             kfac,   ! in    : exch factors
     &                             tol2e,  ! in    : tolerance coeff
     &                             nocc,   ! in    : nr occ shells
     &                             nvirt,  ! in    : nr vir shells
     &                             npol,   ! in    : nr. polarizations
     &                             nbf,    ! in    : nr. basis functions
     &                             nmo,    ! in    : nr. MOs
     &                             debug)  ! in    : =.true. for debugging
c
c Purpose: Updating g_rhs with g_fock from shell_fock_build()
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)
c         

       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "prop.fh"
#include "bgj.fh"
#include "case.fh"
      integer g_rhs,
     &        g_d1,g_d2,g_fock,g_s10
      integer geom,basis,rtdb
      integer npol,disp,ndens,ispin
      integer vectors(npol),nocc(npol),nvirt(npol)
      integer ifld,ndir,nbf,nmo,shift
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3),
     &        clo(3), chi(3)
      logical debug
      double precision jfac(12),kfac(12),tol2e
      external shell_fock_build,
     &         shell_fock_build_cam,
     &         add_fock ! located in hnd_shift_zora.F
      ndir=3 ! nr directions (x,y,z)
c -------- Creating g_d2 ---------------START
c Note.- g_d2 =(g_d1 g_d1)
      ndens=ndir*npol
      clo(1) = ndens*2
      clo(2) = nbf
      clo(3) = nbf
      chi(1) =  1  
      chi(2) = -1 
      chi(3) = -1
      if (.not.nga_create(MT_DBL,3,clo,'g_d2 matrix',
     &                    chi,g_d2)) 
     &  call errquit('gprelim_fock: nga_create failed g_d2',
     &                0,GA_ERR)
       call ga_zero(g_d2)
       blo(1) = 1
       bhi(1) = ndens   
       blo(2) = 1
       bhi(2) = nbf
       blo(3) = 1
       bhi(3) = nbf
      do ispin=1,npol 
       disp=ndens*(ispin-1) 
       alo(1) = disp+1
       ahi(1) = disp+ndens   
       alo(2) = 1
       ahi(2) = nbf
       alo(3) = 1
       ahi(3) = nbf
       call nga_copy_patch('n',g_d1,blo,bhi,
     &                         g_d2,alo,ahi) 
      enddo ! end-loop-ispin
c -------- Creating g_d2 ---------------END
c     Build "first order fock matrix"
      if (use_theory.eq.'dft') then
         if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $     call errquit('hess_cphf: rtdb_put of xc_active failed',0,
     &       RTDB_ERR)
         if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 2))
     $     call errquit('hess_cphf: rtdb_put of calc_type failed',0,
     &       RTDB_ERR)
         if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $     call errquit('hess_cphf: rtdb_put of j_derfit failed',0,
     &       RTDB_ERR)
      endif
      clo(1) = ndir*npol*2
      clo(2) = nbf
      clo(3) = nbf
      chi(1) =  1  
      chi(2) = -1 
      chi(3) = -1
      if (.not.nga_create(MT_DBL,3,clo,'Fock matrix',chi,g_fock)) call 
     &    errquit('giao_b1: nga_create failed g_fock',0,GA_ERR)
      call ga_zero(g_fock)
c      if (ga_nodeid().eq.0)
c     &  write(*,*) 'cam_exch=',cam_exch
c from hnd_giaox.F
c Note: Just the exchange: jfac = 0.d0 (see above)
      if (.not.cam_exch) then

         call shell_fock_build(geom, basis,0,ndir*npol*2,
     $                         jfac,kfac,tol2e,
     &                         g_d2,  ! input
     &                         g_fock,! output
     &                         .false.)

      else

         call shell_fock_build_cam(geom, basis,0,ndir*npol*2,
     $                         jfac,kfac,tol2e,
     &                         g_d2,  ! input
     &                         g_fock,! output
     &                         .false.)

      end if

c      if (ga_nodeid().eq.0) then
c        write(*,70) nocc(1)   ,nocc(2),
c     &              nvirt(1)  ,nvirt(2),npol,nbf,nmo
c 70    format('BEF_add_fock:: nocc =(',i5,',',i5,') ',
c     &        'nvirt=(',i5,',',i5,') ',
c     &        '(npol,nbf,nmo)=(',i3,',',i3,',',i3,')')
c      endif


      if (debug) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '------- g_fock-nw ---- START'
       call ga_print(g_fock)
       if (ga_nodeid().eq.0)
     &  write(*,*) '------- g_fock-nw ---- END'
      endif ! end-if-debug
      if(use_theory.eq.'dft') then
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, 0))
     $      call errquit('giaox: rtdb_put failed',0,RTDB_ERR)
      endif
c Note.- add_fock() is defined in hnd_gshift_zora.F

      call add_fock(g_rhs, ! out: accumulated rhs expression
     &             g_fock, !  in: Fock-term
     &            vectors, !  in: MO  coeffs
     &                nbf, !  in: nr. basis functions
     &                nmo, !  in: nr. MOs (occ+virt)
     &               npol, !  in: nr. of polarizations
     &               nocc, !  in: nr. occ     MOs
     &              nvirt) !  in: nr. virtual MOs

      if (debug) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- g_rhs-AFT-shfock-inside-- START'
        call ga_print(g_rhs)
      if (ga_nodeid().eq.0)
     &  write(*,*) '---- g_rhs-AFT-shfock-inside---  END'
      endif ! end-if-debug
       if (.not.ga_destroy(g_d2)) call 
     &    errquit('giao_b1: ga_destroy failed g_fock',0,GA_ERR)
       if (.not.ga_destroy(g_fock)) call 
     &   errquit('giao_b1: ga_destroy failed g_fock',0,GA_ERR)
       if (.not.ga_destroy(g_d1)) call 
     &   errquit('giao_b1: ga_destroy failed g_d10',0,GA_ERR)
      return
      end

      subroutine get_d1_giao_b1(g_d1,   ! out:
     &                          g_u,    ! in :
     &                          vectors,! in : MO vectors
     &                          nocc,   ! in : nr. occ shells
     &                          npol,   ! in : nr. polarizations
     &                          nbf,    ! in : nr. basis functions
     &                          nmo,    ! in : nr. MOs
     &                          debug)  ! in : =.true. for debugging
c
c Purpose: get_d1 in giao_b1_movecs()
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)

       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "prop.fh"
      integer g_d1,g_s10
      integer ispin,npol
      integer vectors(npol),nocc(npol),
     &        g_u(npol)
      integer ifld,ndir,nbf,nmo,disp
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3),
     &        clo(3), chi(3),
     &        dlo(3), dhi(3)
      logical debug
      ndir=3 ! nr directions (x,y,z)
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 1
      if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',alo,g_s10)) call 
     &    errquit('giao_b1: nga_create failed g_s01',0,GA_ERR)
      clo(1) = ndir*npol
      clo(2) = nbf
      clo(3) = nbf
      chi(1) =  1  
      chi(2) = -1 
      chi(3) = -1
      if (.not.nga_create(MT_DBL,3,clo,'D10 matrix',chi,g_d1)) call 
     &    errquit('giao_b1: nga_create failed g_d1',0,GA_ERR)
      call ga_zero(g_d1)
      do ispin=1,npol
       call ga_zero(g_s10)
       blo(1) = 1
       bhi(1) = nbf
       blo(2) = 1
       clo(2) = 1
       chi(2) = nbf
       clo(3) = 1
       chi(3) = nbf
       dlo(1) = 1
       dhi(1) = nbf
       dlo(2) = 1
       dhi(2) = nocc(ispin)
       disp=ndir*(ispin-1) ! for (clo,chi)
c     Create "perturbed density matrix" for closed-closed g_u block
       do ifld=1,ndir
        alo(1) = 1
        ahi(1) = nmo
        alo(2) = 1
        ahi(2) = nocc(ispin)
        alo(3) = ifld
        ahi(3) = ifld
        dlo(3) = 1
        dhi(3) = 1
        bhi(2) = nmo 
        if (debug) then
         if (ga_nodeid().eq.0) then
          write(*,17) ifld,
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               dlo(1),dhi(1),dlo(2),dhi(2),
     &               dlo(3),dhi(3)
 17       format('1(',i3,')::alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'dlo-dhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ')
         endif
        endif ! end-if-debug
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                vectors(ispin),blo,bhi,  
     &                g_u(ispin)    ,alo,ahi,
     &                g_s10         ,dlo,dhi)  
        alo(1) = 1
        ahi(1) = nocc(ispin)
        alo(2) = 1
        ahi(2) = nbf ! nmo fix lindep 05-02-12
        alo(3) = 1
        ahi(3) = 1
        bhi(2) = nocc(ispin)
        clo(1) = disp+ifld
        chi(1) = disp+ifld
c     Minus sign as we subtract it from the RHS as we do not include 
c     it in the LHS
        if (debug) then
         if (ga_nodeid().eq.0) then
          write(*,16) ifld,
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               clo(1),chi(1),clo(2),chi(2),
     &               clo(3),chi(3)
 16       format('2(',i3,')::alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'clo-chi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ')
         endif
        endif ! end-if-debug
        call nga_matmul_patch('n','t',-1.0d0,0.0d0,
     &             vectors(ispin),blo,bhi,
     &                      g_s10,alo,ahi,
     &                       g_d1,clo,chi)  ! nbf x nbf  for (dir,ispin)
       enddo ! end-loop-ifld
      enddo ! end-loop-ispin
      if (.not.ga_destroy(g_s10)) call 
     &  errquit('giao_b1: ga_destroy failed g_s10',0,GA_ERR)
      return
      end

      subroutine update_rhs_threeAOints(
     &                       g_rhs,  ! in/out: RHS used for cphf2/3
     &                       vectors,! in    : MO vectors
     &                       rtdb,   ! in    : rtdb  handle 
     &                       basis,  ! in    : basis handle 
     &                       nocc,   ! in    : nr occ  shells
     &                       nvirt,  ! in    : nr virt shells
     &                       npol,   ! in    : nr. polarizations
     &                       nbf,    ! in    : nr. basis functions
     &                       nmo,    ! in    : nr. MOs
     &                       debug)  ! in    : logical for debugging
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)

       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "geom.fh"
#include "prop.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "apiP.fh"
#include "bgj.fh"
#include "case.fh"
      integer npol
      integer g_rhs,g_s10,g_s10_1,vectors(npol) 
      integer rtdb,basis
      integer nocc(npol),nvirt(npol),disp,shift,ispin
      integer ifld,ndir,nbf,nmo
      integer alo(3),ahi(3), 
     &        blo(3),bhi(3)
      integer nbq,nextbq,ncosbq ! for COSMO (adding solvent effects)
      logical oskel,debug
      integer nat
      parameter (nat=1)
      double precision origin(3)
      data origin/0d0,0d0,0d0/
      external giao_aotomo
c      external geom_extbq_ncenter
c     Current CPHF does not handle symmetry 
c     Making C1 geometry and store it on rtdb
      oskel = .false.
      ndir=3 ! nr directions (x,y,z)
c     Get H10 in GA, reusing g_s10 array
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = ndir
       if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',
     &                    alo,g_s10_1)) 
     &  call errquit('gprelim_fock: nga_create failed g_s10_1',
     &               0,GA_ERR)
       call ga_zero(g_s10_1)
      ahi(3) = ndir*npol
       if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',
     &                    alo,g_s10)) 
     &  call errquit('gprelim_fock: nga_create failed g_s10',
     &               0,GA_ERR)
      call ga_zero(g_s10)
      call int_giao_1ega(basis,basis,
     &                   g_s10_1,'l10',  ! out: g_s10 
     &                   origin,nat,oskel)
      call int_giao_1ega(basis,basis,
     &                   g_s10_1,'tv10', ! out: g_s10 updated 
     &                   origin,nat,oskel)
c
c     Get external and cosmo bq contribution
      nbq    = 0
      nextbq = 0
      ncosbq = 0
      if(geom_extbq_on()) nextbq = geom_extbq_ncenter()
      nbq = nextbq ! external bq's
      if (rtdb_get(rtdb,'cosmo:nefc',mt_int,1,ncosbq))
     &    nbq = ncosbq ! cosmo bq's
      if (nextbq.gt.0.and.ncosbq.gt.0)
     &    nbq = nextbq + ncosbq  ! tally up cosmo and external bqs
c
c     if (ga_nodeid().eq.0) write(6,*) "nbq: ", nbq
      if (nbq.gt.0) then
        call int_giao_1ega(basis,basis,
     &                     g_s10_1,'bq10', ! out
     &                     origin,nat,oskel)
      end if
c --------- g_s10_1 --> g_s10 --------- START
       blo(1) = 1
       bhi(1) = nbf ! nmo fix lindep 05-02-12
       blo(2) = 1
       bhi(2) = nbf ! nmo fix lindep 05-02-12
       blo(3) = 1
       bhi(3) = ndir
      do ispin=1,npol  
       disp=ndir*(ispin-1) 
       alo(1) = 1
       ahi(1) = nbf ! nmo fix lindep 05-02-12
       alo(2) = 1
       ahi(2) = nbf ! nmo fix lindep 05-02-12
       alo(3) = disp+1
       ahi(3) = disp+ndir  
       call nga_copy_patch('n',g_s10_1,blo,bhi,
     &                           g_s10,alo,ahi) 
      enddo ! end-loop-ispin
c --------- g_s10_1 --> g_s10 --------- END
      if (.not.ga_destroy(g_s10_1)) call 
     &    errquit('giao_b1: ga_destroy failed g_s10',0,GA_ERR)
c
c     ga_rhs(a,i) = ga_rhs(a,i) + H10(a,i)
c     Transform H10 to MO and add to g_rhs
      call giao_aotomo(g_s10,vectors,nocc,nvirt,npol,ndir,nbf)
      do ispin=1,npol
       disp=ndir*(ispin-1) 
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = disp+1
       ahi(3) = disp+ndir
       shift=nocc(1)*nvirt(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = ndir
       call nga_add_patch(1.0d0,g_rhs,blo,bhi,
     &                    1.0d0,g_s10,alo,ahi,
     &                          g_rhs,blo,bhi)
      enddo ! end-loop-ispin
c      if (ga_nodeid().eq.0)
c     & write(*,*) '------- g_s10-nw ---- START'
c      call ga_print(g_s10)
c      if (ga_nodeid().eq.0)
c     & write(*,*) '------- g_s10-nw ---- END'
c
c     Cleanup g_s10 as we do not need it right now
      if (.not.ga_destroy(g_s10)) call 
     &    errquit('giao_b1: ga_destroy failed g_s10',0,GA_ERR)
      return
      end

      subroutine update_rhs_eS10(
     &               g_rhs,  !in/out:
     &               g_u,    !out:
     &               g_sket1,!out:
     &               eval,   !in : energy values
     &               vectors,!in : MO vectors
     &               nocc,   !in : nr.   occupied MOs
     &               nvirt,  !in : nr. unoccupied MOs
     &               npol,   !in : nr. polarizations
     &               nbf,    !in : nr. basis functions
     &               nmo,    !in : nr. MOs
     &               basis,  !in : basis handle
     &               debug)  !in : logical var for debugging
c
c Author : Fredy W. Aquino
c Date   : 03-15-12
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach and appears on nwchem-devtrunk (date:03-02-12)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "geom.fh"
#include "prop.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "apiP.fh"
#include "bgj.fh"
#include "case.fh"
      integer basis
      integer npol,nbf,nmo     
      integer vectors(npol)
      integer nocc(npol),nvirt(npol),
     &        iocc,disp,disp1,ispin,
     &        ndir,shift,
     &        alo(3),ahi(3),
     &        blo(3),bhi(3),
     &        clo(3),chi(3)
      double precision eval(nmo*npol),toscl
      integer g_rhs,g_u(npol),g_s10,g_s10_1,g_sket1
      double precision origin(3)
      data origin/0d0,0d0,0d0/
      logical oskel,debug
      integer nat
      parameter (nat=1)
c     Current CPHF does not handle symmetry 
c     Making C1 geometry and store it on rtdb
      oskel=.false.
      ndir=3 ! nr. directions (x,y,z)
c     NGA dimension arrays for copying will be the same every time
c     Also third NGA dimension for any of the three dimensional
c     arrays will be the same everytime (running from 1 to 3)
c     So, lets define them once and for all in blo and bhi    
c     Get S10 in GA and transform to MO set (virt,occ)
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = ndir
      if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',alo,g_s10_1)) call 
     &    errquit('giao_b1: nga_create failed g_s01',0,GA_ERR)
      call ga_zero(g_s10_1)
      ahi(3) = ndir*npol
       if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',
     &                    alo,g_s10)) 
     &  call errquit('gprelim_fock: nga_create failed g_s10',
     &               0,GA_ERR)
      call ga_zero(g_s10)
      call int_giao_1ega(basis,basis,
     &                   g_s10_1,'s10', ! out: g_s10 FA-LBL-B
     &                   origin,nat,oskel)
c -------- create g_s10 --------------- START
       blo(1) = 1
       bhi(1) = nbf ! nmo fix-lindep 05-02-12
       blo(2) = 1
       bhi(2) = nbf ! nmo fix-lindep 05-02-12
       blo(3) = 1
       bhi(3) = ndir
      do ispin=1,npol  
       disp=ndir*(ispin-1) 
       alo(1) = 1
       ahi(1) = nbf ! nmo fix-lindep 05-02-12
       alo(2) = 1
       ahi(2) = nbf ! nmo fix-lindep 05-02-12
       alo(3) = disp+1
       ahi(3) = disp+ndir
       call nga_copy_patch('n',g_s10_1,blo,bhi,
     &                         g_s10  ,alo,ahi) 
      enddo ! end-loop-ispin
       if (.not.ga_destroy(g_s10_1)) call 
     &   errquit('giao_b1: ga_destroy failed vectors',0,GA_ERR)
c -------- create g_s10 --------------- END
c After giao_aotomo valid dim(g_s10): nmo x nmo
      call giao_aotomo(g_s10,vectors,nocc,nvirt,npol,ndir,nbf) 
      if (debug) write (luout,*) 'S10 done'
c     while we are calculating integrals, let's also determine
c     the 'half' overlap derivative, used later in the calling
c     routine
      call ga_zero(g_sket1)
      call int_giao_1ega(basis,basis,
     &                   g_sket1,'srxRb', ! out: g_sket1 FA-LBL-A
     &                   origin,nat,oskel)
      if (debug) write (luout,*) 'S1-ket done'
c     g_sket1 will not be used further here. It is one of the 
c     output results of this routine.
c     Broceed with the computation of the B-field perturbed
c     MO coefficients.
      do ispin=1,npol
c     ga_rhs(a,i) = ga_rhs(a,i) - e(i) * S10(a,i)
c     Scale (occ,virt) block g_s10 with - (minus) eigenvalues 
       disp1=ndir*(ispin-1)
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(3) = disp1+1
       ahi(3) = disp1+ndir
       clo(1) = 1
       chi(1) = nocc(ispin)
       clo(2) = 1
       chi(2) = nocc(ispin)
       clo(3) = 1
       chi(3) = ndir
c       disp   = nmo*(ispin-1)
       disp   = nbf*(ispin-1) ! fix-lindep 05-02-12
       do iocc=1,nocc(ispin)
        alo(2) = iocc
        ahi(2) = iocc
        toscl  =-eval(disp+iocc) 
        call nga_scale_patch(g_s10,alo,ahi,toscl) 
       enddo ! end-loop-iocc
c     Copy to ga_rhs 
c     alo(1) and ahi(1) the same as before
       alo(2) = 1
       ahi(2) = nocc(ispin)
       shift=nocc(1)*nvirt(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = ndir
       call nga_copy_patch('n',g_s10,alo,ahi,
     &                         g_rhs,blo,bhi)
      
c     Construct occ-occ part of the three U matrices
c     Occ-occ blocks for each field direction are defined as -1/2 S10
c     Scale (occ,occ) block g_s10 with -1/2 and add to g_u
c     alo(2) and ahi(2) will stay as 1 and nclosed(1) for a while
       alo(1) = 1
       ahi(1) = nocc(ispin)
       call nga_scale_patch(g_s10,alo,ahi,-0.5d0)
       call nga_copy_patch('n',g_s10     ,alo,ahi,
     &                         g_u(ispin),clo,chi)
      enddo ! end-loop-ispin
       if (debug) write (luout,*) 'S10 in occ-occ done'
       if (.not.ga_destroy(g_s10)) call 
     &   errquit('giao_b1: ga_destroy failed vectors',0,GA_ERR)
      return
      end
c $Id$
