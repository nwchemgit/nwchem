c===============================================================c
c                                                               c
c       Socrates - an NWChem module for teaching                c
c       Copyright © 2009 Jeff Hammond                           c
c                                                               c
c       Developed by:                                           c
c                                                               c
c               Jeff R. Hammond                                 c
c               Leadership Computing Facility                   c
c               Argonne National Laboratory                     c
c               jhammond@mcs.anl.gov                            c
c                                                               c
c===============================================================c
      logical function socrates_scf_ga(rtdb,geom,num_alfa,num_beta)
c
c $Id$
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "bas.fh"
#include "geom.fh"
#include "sym.fh"
#include "schwarz.fh"
c
c     object handles
c
      integer rtdb             ! RTDB handle
      integer geom             ! GEOM handle
      integer basis            ! BASIS handle
c
c     GA variables
c
      integer bytes            ! Number of bytes in a double
c
c     INT variables
c
      integer i,j,k,l
      integer nbas             ! number of basis sets used
      parameter (nbas=1)
      integer nbf,nshells
      integer s_int2e,k_int2e,l_int2e
      integer s_scr2e,k_scr2e,l_scr2e
      integer ish,ilo,ihi,irng
      integer jsh,jlo,jhi,jrng
      integer ksh,klo,khi,krng
      integer lsh,llo,lhi,lrng
      integer ao_offset
c
c     SCF variables
c
      integer g_smat
      integer g_hmat
      integer g_dmata
      integer g_dmatb
      integer g_fmata
      integer g_fmatb
      integer num_elec,num_alfa, num_beta
c
c     timing variables
c
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
c
c     INT variables
c
      double precision int_thresh
c
c     SCF variables
c
      double precision tol2e
      double precision schwarz_ij,schwarz_kl
c
c     primitive variables
c
      logical nodezero         ! am I node 0?
      logical debug            ! debug mode?
      logical stat             ! status
      logical use_symm         ! use symmetry in fock building
c
c     SYM variables
c
      character*8 group
c
c     external routines
c
      logical int_normalize
      external int_normalize
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      double precision dnrm2
      external dnrm2
c
      parameter (use_symm = .false.)
c
#ifdef DEBUG_PRINT
      debug = (GA_NodeId().eq.0) ! debug print on nodezero only
c      debug = .true. ! debug print everywhere
#else
      parameter (debug = .false.)
#endif
c
      if (debug) then
        write(LuOut,*) 'top of socrates_scf_ga'
        call util_flush(LuOut)
      endif
c
      socrates_scf_ga = .false.
c
      nodezero=(ga_nodeid().eq.0)
c
c      bytes = ma_sizeof(mt_dbl,1,mt_byte)
c
c===============================================================c
c                                                               c
c     Initialize BASIS and INTEGRAL objects                     c
c                                                               c
c===============================================================c
c
c     ---------
c     Basis set
c     ---------
c
      if (.not.bas_create(basis,'ao basis')) then
        call errquit(__FILE__,__LINE__,BASIS_ERR)
      endif
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis')) then
        call errquit(__FILE__,__LINE__,BASIS_ERR)
      endif
      if (.not.bas_numbf(basis,nbf)) then
        call errquit(__FILE__,__LINE__,BASIS_ERR)
      endif
      if (.not.bas_numcont(basis,nshells)) then
        call errquit(__FILE__,__LINE__,BASIS_ERR)
      endif
c
      if (nodezero) then
        call bas_print_labels(basis)
        write(LuOut,*)
        call util_flush(LuOut)
      endif
c
c     ---------
c     Integrals
c     ---------
c
      if (.not.rtdb_get(rtdb,'socrates:tol2e',mt_dbl,1,tol2e)) then
        tol2e = 1.0d-10
      endif
      if (.not.rtdb_get(rtdb,'socrates:int_thresh',mt_dbl,1,
     1                  int_thresh)) then
        int_thresh = 1.0d-20
      endif
      call int_acc_set(int_thresh)
c
      call int_init(rtdb,nbas,basis)
      call schwarz_init(geom,basis)
      if (.not.int_normalize(rtdb,basis)) then
        if (nodezero) write(LuOut,*) 'int_normalize failed'
        call errquit(__FILE__,__LINE__,INT_ERR)
      endif
c
      call int_mem_2e4c(s_int2e,s_scr2e)
      if (.not.ma_push_get(mt_dbl,s_int2e,'int2e',l_int2e,k_int2e)) then
        call errquit(__FILE__,__LINE__,MA_ERR)
      endif
      if (.not.ma_push_get(mt_dbl,s_scr2e,'scr2e',l_scr2e,k_scr2e)) then
        call errquit(__FILE__,__LINE__,MA_ERR)
      endif
c
c===============================================================c
c                                                               c
c     Begin SCF code                                            c
c                                                               c
c===============================================================c
c
c     initialize overlap matrix
c
      g_smat  = ga_create_atom_blocked(geom, basis,'smat')
      call int_1e_ga(basis,basis,g_smat,'overlap',use_symm)
c
c     initialize Hcore matrix
c
      g_hmat  = ga_create_atom_blocked(geom, basis,'hmata')
      call int_1e_ga(basis,basis,g_hmat,'kinetic',use_symm)
      call int_1e_ga(basis,basis,g_hmat,'potential',use_symm)
c
c     initialize density matrices
c
      g_dmata = ga_create_atom_blocked(geom, basis,'dmata')
      g_dmatb = ga_create_atom_blocked(geom, basis,'dmatb')
!
!     initialize density matrix somehow here
!
c
c     create Fock matrices
c
      g_fmata = ga_create_atom_blocked(geom, basis,'fmata')
      g_fmatb = ga_create_atom_blocked(geom, basis,'fmatb')

c
c
c


c
#ifdef DEBUG_PRINT
      if (nodezero) then
        write(LuOut,*) 'Printing all non-negligible AOs'
        call util_flush(LuOut)
      endif
#endif
      do ish=1,nshells
        if (.not.bas_cn2bfr(basis,ish,ilo,ihi)) then
          call errquit(__FILE__,__LINE__,BASIS_ERR)
        endif
        irng = ihi - ilo + 1
        do jsh=1,nshells
          if (.not.bas_cn2bfr(basis,jsh,jlo,jhi)) then
            call errquit(__FILE__,__LINE__,BASIS_ERR)
          endif
          jrng = jhi - jlo + 1
          schwarz_ij = schwarz_shell(ish,jsh)
          do ksh=1,nshells
            if (.not.bas_cn2bfr(basis,ksh,klo,khi)) then
              call errquit(__FILE__,__LINE__,BASIS_ERR)
            endif
            krng = khi - klo + 1
            do lsh=1,nshells
              if (.not.bas_cn2bfr(basis,lsh,llo,lhi)) then
                call errquit(__FILE__,__LINE__,BASIS_ERR)
              endif
              lrng = lhi - llo + 1
              schwarz_kl = schwarz_shell(ksh,lsh)
              if ((schwarz_ij*schwarz_kl).ge.tol2e) then
                call int_2e4c(basis,ish,jsh,basis,ksh,lsh,
     1                        s_scr2e,dbl_mb(k_scr2e),
     2                        s_int2e,dbl_mb(k_int2e))
#ifdef DEBUG_PRINT
                call socrates_print_ao2e(ilo,ihi,jlo,jhi,
     1                                   klo,khi,llo,lhi,
     2                                   dbl_mb(k_int2e),tol2e)
#endif
              endif
            enddo
          enddo
        enddo
      enddo
#ifdef DEBUG_PRINT
      if (nodezero) then
        write(LuOut,*) 'End of integral printing'
        write(LuOut,*)
        call util_flush(LuOut)
      endif
#endif
c
c     ==========
c     Deallocate
c     ==========
c
      stat = ga_destroy(g_dmata)
      if (.not.stat) then
        call errquit(__FILE__,__LINE__,GEOM_ERR)
      endif
      stat = ga_destroy(g_dmatb)
      if (.not.stat) then
        call errquit(__FILE__,__LINE__,GEOM_ERR)
      endif
      stat = ga_destroy(g_fmata)
      if (.not.stat) then
        call errquit(__FILE__,__LINE__,GEOM_ERR)
      endif
      stat = ga_destroy(g_fmatb)
      if (.not.stat) then
        call errquit(__FILE__,__LINE__,GEOM_ERR)
      endif
      stat = ga_destroy(g_hmat)
      if (.not.stat) then
        call errquit(__FILE__,__LINE__,GEOM_ERR)
      endif
      stat = ga_destroy(g_smat)
      if (.not.stat) then
        call errquit(__FILE__,__LINE__,GEOM_ERR)
      endif
c
c===============================================================c
c                                                               c
c     Close any open object references and stack blocks         c
c                                                               c
c===============================================================c
c
#ifdef DETAILED_FREE
      if (.not.ma_pop_stack(l_scr2e)) then
        call errquit('ncc_driver: MA problem scr2e ',0,MA_ERR)
      endif
      if (.not.ma_pop_stack(l_int2e)) then
        call errquit('ncc_driver: MA problem int2e ',0,MA_ERR)
      endif
#else
      if (.not. ma_chop_stack(l_int2e)) then
        call errquit('ncc_driver: stack corrupt ',0, MA_ERR)
      endif
#endif
c
      call schwarz_tidy()
      call int_terminate()
c
      if (.not.bas_destroy(basis)) then
        call errquit(__FILE__,__LINE__,BASIS_ERR)
      endif
c
c===============================================================c
c                                                               c
c                           THE END                             c
c                                                               c
c===============================================================c
c
      socrates_scf_ga = .true.
c
      if (debug) then
        write(LuOut,*) 'end of socrates_scf_ga'
        call util_flush(LuOut)
      endif
c
      return
c
 2001 format(/,3x,'for ',a5,' MO vectors:',/,
     1         3x,'number of basis functions    = ',i8,/,
     2         3x,'number of molecular orbitals = ',i8,/)

c
      end
