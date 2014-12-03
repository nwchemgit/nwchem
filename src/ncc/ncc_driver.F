c===============================================================c
c                                                               c
c       NCC - a New Coupled-Cluster code for NWChem             c
c                                                               c
c       Developed by:                                           c
c                                                               c
c               Jeff R. Hammond                                 c
c               Leadership Computing Facility                   c
c               Argonne National Laboratory                     c
c               jhammond@mcs.anl.gov                            c
c                                                               c
c               Karol Kowalski                                  c
c               Environmental Molecular Sciences Laboratory     c
c               Pacific Northwest National Laboratory           c
c               karol.kowalski@pnl.gov                          c
c                                                               c
c               Marta Włoch                                     c
c               Department of Chemistry                         c
c               Michigan Technological University               c
c               wloch@mtu.edu                                   c
c                                                               c
c===============================================================c
      logical function ncc_driver(rtdb)
c
c $Id$
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "inp.fh"
#include "bas.fh"
#include "geom.fh"
#include "sym.fh"
#include "schwarz.fh"
c===============================================================c
c                                                               c
#include "ncc_input.fh"
c                                                               c
c      ...defines these common variables:                       c
c                                                               c
c      integer debugmode                                        c
c      double precision thresh                                  c
c      character*20 method                                      c
c                                                               c
c      I am using common variables to pass from this routine    c
c      into the NCC module so that I do not have to use the     c
c      RTDB since that would involve disk access.               c
c                                                               c
c===============================================================c
c
c     object handles
c
      integer rtdb             ! RTDB handle
      integer geom             ! GEOM handle
      integer basis            ! BASIS handle
c
c     low-level (GA, etc.) variables
c
      integer bytes                           ! Number of bytes in a double
      integer i,j,k,l
c
c     SCF code variables
c
      integer reference
      integer mult                            ! Ground state spin multiplicity
      integer ipol
      integer nclosed, nopen
      integer nbf,nshells                     ! Number of AO basis functions and shells
      integer nbfold
      integer nocc(2)                         ! Number of alpha, beta electrons
      integer nvir(2)                         ! Number of alpha, beta electrons
      integer nmo(2)                          ! Number of alpha, beta orbitals
      integer nfc(2),nfv(2)                   ! Number of frozen cores/virtuals
c
c     GA and MA handles
c
      integer g_moveca, k_moveca, l_moveca
      integer g_movecb, k_movecb, l_movecb
      integer k_eiga,l_eiga
      integer k_eigb,l_eigb
      integer k_occa,l_occa
      integer k_occb,l_occb
      integer k_irra,l_irra
      integer k_irrb,l_irrb
c
c     low-level (GA, etc.) variables
c
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
c
c     SCF code variables
c
      double precision scf_energy
c
c     four-index variables
c
      double precision blk_norm
      integer k_v2oooo,l_v2oooo
      integer k_v2ooov,l_v2ooov
      integer k_v2oovv,l_v2oovv
      integer k_v2ovvv,l_v2ovvv
      integer k_v2vvvv,l_v2vvvv
c
c     amplitudes
c
      integer t2info(8)                ! T2 amplitude "struct"
c
c     low-level (GA, etc.) variables
c
      logical nodezero         ! am I node 0?
      logical debug            ! debug mode?
      logical stat             ! status
c
c     SCF code variables
c
      logical scf_status       ! did SCF finish?
c
c     four-index variables
c
      logical use_symm
c
c     SCF code variables
c
      character*4 scftype
      character*8 group
      character*255 scf_movecs ! SCF MO vector file name
      character*255 title      ! Title
      character*255 basisname  ! Basis set name
c
c     external routines
c
      logical scf
      external scf
      logical int_normalize
      external int_normalize
      logical movecs_read
      external movecs_read
      logical movecs_read_header
      external movecs_read_header
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      double precision dnrm2
      external dnrm2
c
#ifdef DEBUG_PRINT
      print*,'top of ncc_driver'
      debug = (ga_nodeid().eq.0)
c      debug = .true. ! debug from all nodes
#else
      debug = .false.
#endif
c
      ncc_driver = .false.
c
      nodezero = (ga_nodeid().eq.0)
c
      bytes = ma_sizeof(mt_dbl,1,mt_byte)
c
c===============================================================c
c                                                               c
c     Run the SCF code, get GEOM, BASIS and MOVECS information  c
c                                                               c
c===============================================================c
c
      scf_status = scf(rtdb)
c
      if (.not.rtdb_get(rtdb,'scf:energy',mt_dbl,1,scf_energy)) then
        call errquit('ncc_driver: rtdb_get failed',0,RTDB_ERR)
      endif
c
c     ---------------------
c     Geometry and symmetry
c     ---------------------
c
      if (.not.geom_create(geom,'geometry')) then
        call errquit('ncc_driver: geom_create failed',0,GEOM_ERR)
      endif
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) then
        call errquit('ncc_driver: no geometry found',0,GEOM_ERR)
      endif
c
      call sym_group_name(geom,group)
      use_symm = (group(1:2) .ne. 'C1')
c
c     ---------
c     Basis set
c     ---------
c
      if (.not.bas_create(basis,'ao basis')) then
        call errquit('ncc_driver: bas_create failed',0,BASIS_ERR)
      endif
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis')) then
        call errquit('ncc_driver: no ao basis found',0,BASIS_ERR)
      endif
      if (.not.bas_numbf(basis,nbf)) then
        call errquit('ncc_driver: failed to read nbf',0,BASIS_ERR)
      endif
c
c     ----------------------
c     Integral code
c     ----------------------
c
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      if (.not.int_normalize(rtdb,basis)) then
        call errquit('ncc_driver: int_normalize failed',0,INT_ERR)
      endif
c
c     ------------------------------------------------------
c     Number of electrons and orbitals and spin multiplicity
c     ------------------------------------------------------
c
      if(.not.rtdb_cget(rtdb,'scf:scftype',1,scftype)) then
        call errquit('ncc_driver: failed to read scftype',0,RTDB_ERR)
      endif
      if (scftype.eq.'RHF') then
        reference = 0
      else if(scftype.eq.'ROHF')then
        reference = 1
      else if (scftype.eq.'UHF') then
        reference = 2
      else
        call errquit('ncc_driver: invalid reference',0,RTDB_ERR)
        reference = -1
      endif
c
      if(.not. geom_freeze(rtdb, geom,'ncc',nfc(1))) then
         call errquit('ncc_driver: geom_freeze failed',0,GEOM_ERR)
      endif
      nfc(2) = nfc(1) ! always freeze in pairs
c
      if (.not.rtdb_get(rtdb,'scf:nclosed',mt_int,1,nclosed)) then
        call errquit('ncc_driver: failed to read nclosed',0,RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'scf:nopen',mt_int,1,nopen)) then
        call errquit('ncc_driver: failed to read nopen',0,RTDB_ERR)
      endif
      nocc(1) = nclosed + nopen
      nocc(2) = nclosed
c
c     ----------
c     MO vectors
c     ----------
c
      if (.not.rtdb_cget(rtdb,'scf:output vectors',1,scf_movecs)) then
        call errquit('ncc_driver: SCF movecs not found',0,RTDB_ERR)
      endif

      g_moveca = ga_create_atom_blocked(geom,basis,'alpha movecs')
      g_movecb = ga_create_atom_blocked(geom,basis,'beta movecs')
c
      if (.not. movecs_read_header(scf_movecs,title,basisname,
     1          scftype,nbfold,ipol,nmo,2)) then
        call errquit('ncc_driver: cannot read MO vectors header',0,
     1               INPUT_ERR)
      endif
      nmo(2) = nmo(1)
c
      if (nbf.ne.nbfold) then
        call errquit('ncc_driver: nbf-nbfold = ',nbf-nbfold,CALC_ERR)
      endif
c
      if (.not.ma_push_get(mt_dbl,nbf,'eiga',l_eiga,k_eiga)) then
          call errquit ('ncc_driver: ma_push_get eiga',nbf,MA_ERR)
      endif
      if (.not.ma_push_get(mt_dbl,nbf,'eigb',l_eigb,k_eigb)) then
          call errquit ('ncc_driver: ma_push_get eigb',nbf,MA_ERR)
      endif
c
      if (.not.ma_push_get(mt_dbl,nbf,'occa',l_occa,k_occa)) then
          call errquit ('ncc_driver: ma_push_get occa',nbf,MA_ERR)
      endif
      if (.not.ma_push_get(mt_dbl,nbf,'occb',l_occb,k_occb)) then
          call errquit ('ncc_driver: ma_push_get occb',nbf,MA_ERR)
      endif
c
      if (.not.ma_push_get(mt_int,nbf,'irra',l_irra,k_irra)) then
          call errquit ('ncc_driver: ma_push_get irra',nbf,MA_ERR)
      endif
      if (.not.ma_push_get(mt_int,nbf,'irrb',l_irrb,k_irrb)) then
          call errquit ('ncc_driver: ma_push_get irrb',nbf,MA_ERR)
      endif
c
c     read alpha movecs
c
      stat = movecs_read(scf_movecs,1,dbl_mb(k_occa),
     1                   dbl_mb(k_eiga),g_moveca)
      if (.not.stat) then
        call errquit('ncc_driver: cannot read MO vectors',1,CALC_ERR)
      endif
      if (.not.rtdb_get(rtdb,'scf:alpha irreps',mt_int,nbf,
     1                  int_mb(k_irra))) then
        call errquit('ncc_driver: no alpha irreps found',1,RTDB_ERR)
      endif
c
c     read beta movecs
c     max(reference,1) means beta is different only for UHF
c
      stat = movecs_read(scf_movecs,max(reference,1),dbl_mb(k_occb),
     1                   dbl_mb(k_eigb),g_movecb)
      if (.not.stat) then
        call errquit('ncc_driver: cannot read MO vectors',2,CALC_ERR)
      endif
      if (reference .eq. 0) then ! RHF
        do i=1,nbf
          int_mb(k_irrb+i-1) = int_mb(k_irra+i-1)
        enddo
      else ! ROHF or UHF
        if (.not.rtdb_get(rtdb,'scf:beta irreps',mt_int,nbf,
     1                    int_mb(k_irrb))) then
          call errquit('ncc_driver: no beta irreps found',1,RTDB_ERR)
        endif
      endif
c
c===============================================================c
c                                                               c
c     Author information etc.                                   c
c                                                               c
c===============================================================c
c
      if (nodezero) then
        write(LuOut,1000)
        call util_flush(LuOut)
      endif
c
c===============================================================c
c                                                               c
c     Four-index transformation algorithms                      c
c                                                               c
c===============================================================c
c
      if (use_symm) then
!        call errquit('ncc_driver: disable symmetry',0,INPUT_ERR)
        write(LuOut,*) ' BE CAREFUL!  Using symmetry: ',group
        call util_flush(LuOut)
      endif
c
c     =======================
c     Get MO vectors
c     =======================
c
      if (nodezero) call bas_print_labels(basis)
c
      call scf_lindep(rtdb, geom, basis, g_moveca, nmo(1))
      if (nodezero) write(LuOut,2001) 'alpha',nbf,nmo(1)
c
      if (.not.ma_push_get(mt_dbl,nbf*nmo(1),'moveca',
     1                     l_moveca,k_moveca)) then
          call errquit ('ncc_driver: ma_push_get moveca',nbf,MA_ERR)
      endif
      call ga_get(g_moveca,1,nbf,1,nmo(1),dbl_mb(k_moveca),nbf)
!       call ga_print(g_moveca)
!       call ncc_print_movecs(nbf,nmo(1),dbl_mb(k_moveca))
c      call ncc_print_dvec(nbf,dbl_mb(k_eiga))
c      call ncc_print_dvec(nbf,dbl_mb(k_occa))
c      call ncc_print_ivec(nbf,int_mb(k_irra))
c
      call scf_lindep(rtdb, geom, basis, g_movecb, nmo(2))
      if (nodezero) write(LuOut,2001) 'beta ',nbf,nmo(2)
c
      if (.not.ma_push_get(mt_dbl,nbf*nmo(2),'movecb',
     1                     l_movecb,k_movecb)) then
          call errquit ('ncc_driver: ma_push_get movecb',nbf,MA_ERR)
      endif
      call ga_get(g_moveca,1,nbf,1,nmo(2),dbl_mb(k_movecb),nbf)
!       call ga_print(g_movecb)
!       call ncc_print_movecs(nbf,nmo(2),dbl_mb(k_movecb))
c      call ncc_print_dvec(nbf,dbl_mb(k_eigb))
c      call ncc_print_dvec(nbf,dbl_mb(k_occb))
c      call ncc_print_ivec(nbf,int_mb(k_irrb))
c
c     =======================
c     Computer 2-e integrals
c     =======================
c
      nvir(1) = nmo(1) - nocc(1)
      nvir(2) = nmo(2) - nocc(2)
! c
! c     allocate v2oooo
! c
!       if (.not.ma_push_get(mt_dbl,nocc(1)*nocc(1)*nocc(1)*nocc(1),
!      1                     'v2oooo',l_v2oooo,k_v2oooo)) then
!           call errquit ('ncc_driver: ma_push_get v2oooo',
!      1                  nocc(1)*nocc(1)*nocc(1)*nocc(1),MA_ERR)
!       endif
! c
! c     allocate v2ooov
! c
!       if (.not.ma_push_get(mt_dbl,nocc(1)*nocc(1)*nocc(1)*nvir(1),
!      1                     'v2ooov',l_v2ooov,k_v2ooov)) then
!           call errquit ('ncc_driver: ma_push_get v2ooov',
!      1                  nocc(1)*nocc(1)*nocc(1)*nvir(1),MA_ERR)
!       endif
! c
! c     allocate v2oovv
! c
!       if (.not.ma_push_get(mt_dbl,nocc(1)*nocc(1)*nvir(1)*nvir(1),
!      1                     'v2oovv',l_v2oovv,k_v2oovv)) then
!           call errquit ('ncc_driver: ma_push_get v2oovv',
!      1                  nocc(1)*nocc(1)*nvir(1)*nvir(1),MA_ERR)
!       endif
! c
! c     allocate v2ovvv
! c
!       if (.not.ma_push_get(mt_dbl,nocc(1)*nvir(1)*nvir(1)*nvir(1),
!      1                     'v2ovvv',l_v2ovvv,k_v2ovvv)) then
!           call errquit ('ncc_driver: ma_push_get v2ovvv',
!      1                  nocc(1)*nvir(1)*nvir(1)*nvir(1),MA_ERR)
!       endif
! c
! c     allocate v2vvvv
! c
!       if (.not.ma_push_get(mt_dbl,nvir(1)*nvir(1)*nvir(1)*nvir(1),
!      1                     'v2vvvv',l_v2vvvv,k_v2vvvv)) then
!           call errquit ('ncc_driver: ma_push_get v2vvvv',
!      1                  nvir(1)*nvir(1)*nvir(1)*nvir(1),MA_ERR)
!       endif
c
!       call ncc_fourindex(rtdb,basis,
!      1                   nocc(1),nocc(2),nvir(1),nvir(2),
!      2                   nmo(1),nmo(2),nbf,
!      3                   dbl_mb(k_v2oooo),dbl_mb(k_v2ooov),
!      4                   dbl_mb(k_v2oovv),dbl_mb(k_v2ovvv),
!      5                   dbl_mb(k_v2vvvv),
!      6                   dbl_mb(k_moveca),dbl_mb(k_movecb))
c
      call ncc_doubles_create_aa(t2info,      ! output array
     1                           "T2",        ! character name for this array
     2                           nocc(1),     ! number of occupied orbitals
     3                           nvir(1),     ! number of virtual orbitals
     4                           occtile,     ! tilesize of occupied orbitals
     5                           virtile)     ! tilesize of virtual orbitals
c
      call ncc_doubles_destroy(t2info)
c
c===============================================================c
c                                                               c
c     Close any open object references and stack blocks         c
c                                                               c
c===============================================================c
c
#ifdef DETAILED_FREE
!       if (.not.ma_pop_stack(l_v2vvvv)) then
!         call errquit('ncc_driver: MA problem v2vvvv ',0,MA_ERR)
!       endif
!       if (.not.ma_pop_stack(l_v2ovvv)) then
!         call errquit('ncc_driver: MA problem v2ovvv ',0,MA_ERR)
!       endif
!       if (.not.ma_pop_stack(l_v2oovv)) then
!         call errquit('ncc_driver: MA problem v2oovv ',0,MA_ERR)
!       endif
!       if (.not.ma_pop_stack(l_v2ooov)) then
!         call errquit('ncc_driver: MA problem v2ooov ',0,MA_ERR)
!       endif
!       if (.not.ma_pop_stack(l_v2oooo)) then
!         call errquit('ncc_driver: MA problem v2oooo ',0,MA_ERR)
!       endif
c
      if (.not.ma_pop_stack(l_movecb)) then
        call errquit('ncc_driver: MA problem movecb ',0,MA_ERR)
      endif
      if (.not.ma_pop_stack(l_moveca)) then
        call errquit('ncc_driver: MA problem moveca ',0,MA_ERR)
      endif
      if (.not.ma_pop_stack(l_irrb)) then
        call errquit('ncc_driver: MA problem irrb ',0,MA_ERR)
      endif
      if (.not.ma_pop_stack(l_irra)) then
        call errquit('ncc_driver: MA problem irra ',0,MA_ERR)
      endif
      if (.not.ma_pop_stack(l_occb)) then
        call errquit('ncc_driver: MA problem occb ',0,MA_ERR)
      endif
      if (.not.ma_pop_stack(l_occa)) then
        call errquit('ncc_driver: MA problem occa ',0,MA_ERR)
      endif
      if (.not.ma_pop_stack(l_eigb)) then
        call errquit('ncc_driver: MA problem eigb ',0,MA_ERR)
      endif
      if (.not.ma_pop_stack(l_eiga)) then
        call errquit('ncc_driver: MA problem eiga ',0,MA_ERR)
      endif
#else
      if (.not. ma_chop_stack(l_eiga)) then
        call errquit('ncc_driver: stack corrupt ',0, MA_ERR)
      endif
#endif
c
      call schwarz_tidy()
      call int_terminate()
c
      if (.not.ga_destroy(g_moveca)) then
        call errquit('ncc_driver: failed to destroy g_moveca',1,GA_ERR)
      endif
      if (.not.ga_destroy(g_movecb)) then
        call errquit('ncc_driver: failed to destroy g_movecb',2,GA_ERR)
      endif
c
      if (.not.bas_destroy(basis)) then
        call errquit('ncc_driver: bas_destroy failed',0,BASIS_ERR)
      endif
      if (.not.geom_destroy(geom)) then
        call errquit('tce_tidy: geom_destroy failed',0,GEOM_ERR)
      endif
c
c===============================================================c
c                                                               c
c                           THE END                             c
c                                                               c
c===============================================================c
c
      ncc_driver = .true.
c
#ifdef DEBUG_PRINT
      print*,'end of ncc_driver'
#endif
c
      return
c
 1000 format(/,' ==============================================',/,/,
     &         3x,'NCC - a New Coupled-Cluster code for NWChem',/,
     &         3x,'Developed by:',/,/,
     &         3x,'Jeff R. Hammond',/,
     &         3x,'Leadership Computing Facility',/,
     &         3x,'Argonne National Laboratory',/,
     &         3x,'jhammond@mcs.anl.gov',/,/,
     &         3x,'Karol Kowalski',/,
     &         3x,'Environmental Molecular Sciences Laboratory',/,
     &         3x,'Pacific Northwest National Laboratory',/,
     &         3x,'karol.kowalski@pnl.gov',/,/,
     &         3x,'Marta Włoch',/,
     &         3x,'Department of Chemistry',/,
     &         3x,'Michigan Technological University',/,
     &         3x,'wloch@mtu.edu',/,/,
     &         ' ==============================================',/)
 2001 format(/,3x,'for ',a5,' MO vectors:',/,
     &         3x,'number of basis functions    = ',i8,/,
     &         3x,'number of molecular orbitals = ',i8,/)

c
      end
